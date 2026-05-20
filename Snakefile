# ================================================================================
# PacBio HiFi WGS 3차 분석 파이프라인 (Snakemake)
# ================================================================================

import os
import csv
import json
from pathlib import Path

# 설정 파일 로드
# --configfile 옵션으로 전달됨 (run_pipeline.sh에서 처리)
# 환경 변수 SNAKEMAKE_CONFIG는 하위 호환성을 위해 유지
CONFIG_FILE = os.environ.get("SNAKEMAKE_CONFIG", "config/config.yaml")
# configfile: CONFIG_FILE  # Snakemake 9.x에서는 --configfile 옵션 사용

# ── 다중 군(groups) 로딩 ──────────────────────────────────────────────────────
# samples.groups가 있으면 다중 군 모드, 없으면 control/experimental 2군으로 폴백
GROUPS = config["samples"].get("groups")
if GROUPS:
    ALL_SAMPLES = [s for grp in GROUPS.values() for s in grp]
else:
    GROUPS = {
        "control":      config["samples"]["control"],
        "experimental": config["samples"]["experimental"],
    }
    ALL_SAMPLES = config["samples"]["control"] + config["samples"]["experimental"]

GROUP_NAMES = list(GROUPS.keys())
GROUP_PAIRS = [(a, b) for i, a in enumerate(GROUP_NAMES) for b in GROUP_NAMES[i+1:]]

# 출력 디렉토리 설정
OUTPUT_DIR = config["paths"]["output_dir"]

# ================================================================================
# filelist.csv에서 파일 경로 매핑 로드
# ================================================================================

SAMPLE_FILES = {}

filelist_csv = config["paths"].get("filelist_csv", "filelist.csv")

if os.path.exists(filelist_csv):
    # filelist.csv가 있으면 사용
    with open(filelist_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_id = row['sample_id']
            SAMPLE_FILES[sample_id] = {
                'phased_small_variant_vcf': row['phased_small_variant_vcf'],
                'phased_sv_vcf': row['phased_sv_vcf'],
                'cpg_combined_bed': row['cpg_combined_bed'],
                'cpg_combined_bw': row['cpg_combined_bw'],
                # 선택적 입력 (없으면 빈 문자열)
                'phased_trgt_vcf': row.get('phased_trgt_vcf', ''),
                'cpg_hap1_bed': row.get('cpg_hap1_bed', ''),
                'cpg_hap2_bed': row.get('cpg_hap2_bed', ''),
            }
    print(f"✓ filelist.csv 로드 완료: {len(SAMPLE_FILES)}개 샘플")
else:
    # filelist.csv가 없으면 기존 방식 (wdl_out_dir 사용)
    print(f"⚠️  경고: {filelist_csv}를 찾을 수 없습니다.")
    print(f"   'python bin/prepare_inputs.py'를 먼저 실행하거나")
    print(f"   기존 방식을 사용하려면 paths.wdl_out_dir을 설정하세요.")
    
    # 기존 방식으로 폴백 (wdl_out_dir이 설정된 경우에만)
    wdl_out_dir = config["paths"].get("wdl_out_dir", "")
    if not wdl_out_dir:
        raise WorkflowError(
            f"{filelist_csv} 파일이 없고 paths.wdl_out_dir도 미설정입니다.\n"
            "  먼저 다음 명령어를 실행하세요:\n"
            "    python bin/prepare_inputs.py --config-file <config.yaml>\n"
            "  또는 config.yaml에 paths.wdl_out_dir을 설정하세요."
        )
    for sample in ALL_SAMPLES:
        SAMPLE_FILES[sample] = {
            'phased_small_variant_vcf': os.path.join(
                wdl_out_dir, f"phased_small_variant_{sample}.vcf.gz"
            ) if wdl_out_dir else "",
            'phased_sv_vcf': os.path.join(
                wdl_out_dir, f"phased_sv_{sample}.vcf.gz"
            ) if wdl_out_dir else "",
            'cpg_combined_bed': os.path.join(
                wdl_out_dir, f"cpg_combined_{sample}.bed"
            ) if wdl_out_dir else "",
            'cpg_combined_bw': os.path.join(
                wdl_out_dir, f"cpg_combined_{sample}.bw"
            ) if wdl_out_dir else "",
            'phased_trgt_vcf': "",
            'cpg_hap1_bed':    "",
            'cpg_hap2_bed':    "",
        }

# ================================================================================
# 최종 출력 파일 정의
# ================================================================================

rule all:
    input:
        # Small Variant 필터링 및 주석
        expand(
            os.path.join(OUTPUT_DIR, "small_variants", "{sample}.filtered.tsv"),
            sample=ALL_SAMPLES
        ),
        expand(
            os.path.join(OUTPUT_DIR, "small_variants", "{sample}.annotated.vcf.gz"),
            sample=ALL_SAMPLES
        ),
        expand(
            os.path.join(OUTPUT_DIR, "small_variants", "{sample}.annotated.tsv"),
            sample=ALL_SAMPLES
        ),

        # Structural Variant 필터링 + 유전자 영향 주석
        expand(
            os.path.join(OUTPUT_DIR, "structural_variants", "{sample}.filtered_sv.vcf.gz"),
            sample=ALL_SAMPLES
        ),
        expand(
            os.path.join(OUTPUT_DIR, "structural_variants", "{sample}.annotated_sv.tsv"),
            sample=ALL_SAMPLES
        ),

        # TRGT 반복서열 이상값 (파일이 있는 샘플만)
        [
            os.path.join(OUTPUT_DIR, "trgt", f"{s}.trgt_outliers.tsv")
            for s in ALL_SAMPLES
            if SAMPLE_FILES[s].get('phased_trgt_vcf')
        ],

        # Allele-specific methylation (hap1/hap2 파일이 있는 샘플만)
        [
            os.path.join(OUTPUT_DIR, "asm", f"{s}.asm_results.csv")
            for s in ALL_SAMPLES
            if SAMPLE_FILES[s].get('cpg_hap1_bed') and SAMPLE_FILES[s].get('cpg_hap2_bed')
        ],

        # 메틸화 파일 정리
        os.path.join(OUTPUT_DIR, "methylation", "methylation_summary.txt"),

        # ── Phase 4/5: 코호트 수준 통계 분석 ──────────────────────────────────
        # Small Variant burden test (per-sample canonical impact TSV 전처리 포함)
        expand(
            os.path.join(OUTPUT_DIR, "small_variants", "{sample}.canonical_impacts.tsv"),
            sample=ALL_SAMPLES
        ),
        os.path.join(OUTPUT_DIR, "cohort", "variant_burden_stats.csv"),

        # SV burden test
        os.path.join(OUTPUT_DIR, "cohort", "sv_burden_stats.csv"),

        # TRGT 반복서열 군 간 비교 (TRGT 데이터 있는 샘플만)
        *([os.path.join(OUTPUT_DIR, "cohort", "trgt_group_compare.csv")]
          if any(SAMPLE_FILES[s].get('phased_trgt_vcf') for s in ALL_SAMPLES)
          else []),

        # DMR pairwise (메틸화 데이터 있는 군 쌍에 대해)
        *[os.path.join(OUTPUT_DIR, "cohort", f"dmr_{g1}_vs_{g2}.csv")
          for g1, g2 in GROUP_PAIRS
          if all(any(SAMPLE_FILES[s].get('cpg_combined_bed') for s in GROUPS[g])
                 for g in [g1, g2])]

# ================================================================================
# Phase 1: Small Variant 필터링 (slivar)
# ================================================================================

rule filter_small_variants:
    input:
        vcf = lambda wildcards: SAMPLE_FILES[wildcards.sample]['phased_small_variant_vcf']
    output:
        filtered_vcf = os.path.join(OUTPUT_DIR, "small_variants", "{sample}.filtered.vcf.gz"),
        filtered_tsv = os.path.join(OUTPUT_DIR, "small_variants", "{sample}.filtered.tsv")
    params:
        species = config["parameters"]["general"]["species"],
        max_af = config["parameters"]["slivar"]["max_af"],
        min_gq = config["parameters"]["slivar"]["min_gq"],
        min_dp = config["parameters"]["slivar"]["min_dp"],
        hpo_terms = config["parameters"]["slivar"]["hpo_terms"],
        # Species-specific frequency database
        gnotate_db = lambda wildcards: config["parameters"]["slivar"]["frequency_db"].get(
            config["parameters"]["general"]["species"], ""
        ),
        # Species-specific frequency field name
        freq_field = lambda wildcards: config["parameters"]["slivar"]["frequency_field"].get(
            config["parameters"]["general"]["species"], "gnomad_popmax_af"
        )
    threads:
        config["resources"]["slivar"]["threads"]
    resources:
        mem_mb = config["resources"]["slivar"]["mem_mb"]
    log:
        os.path.join(OUTPUT_DIR, "logs", "slivar", "{sample}.log")
    conda:
        "environment.yaml"
    shell:
        """
        echo "Species: {params.species}" > {log}
        echo "Frequency database: {params.gnotate_db}" >> {log}
        echo "Frequency field: {params.freq_field}" >> {log}
        echo "Min GQ: {params.min_gq}, Min DP: {params.min_dp}" >> {log}
        
        # bcftools로 PASS 필터 + 품질 필터링
        bcftools view \
            -f PASS \
            -i 'FORMAT/GQ >= {params.min_gq} && FORMAT/DP >= {params.min_dp}' \
            -O z \
            -o {output.filtered_vcf} \
            {input.vcf} \
            2>> {log}
        
        # 인덱스 생성
        bcftools index -t {output.filtered_vcf} 2>> {log}
        
        # TSV 변환 (연구자용)
        bcftools query \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT\t%DP\t%GQ]\n' \
            -H \
            {output.filtered_vcf} \
            > {output.filtered_tsv} \
            2>> {log}
        
        echo "Filtering completed successfully" >> {log}
        """

# ================================================================================
# Phase 2: VEP 주석 추가
# ================================================================================

rule annotate_vep:
    input:
        vcf = os.path.join(OUTPUT_DIR, "small_variants", "{sample}.filtered.vcf.gz")
    output:
        annotated_vcf = os.path.join(OUTPUT_DIR, "small_variants", "{sample}.annotated.vcf.gz"),
        stats = os.path.join(OUTPUT_DIR, "small_variants", "{sample}.vep_stats.html")
    params:
        ref_genome = config["paths"]["ref_genome"],
        vep_cache = config["paths"]["vep_cache"],
        pick_order = config["parameters"]["vep"]["pick_order"],
        # Build plugin flags safely: pass each plugin as its own --plugin argument.
        # If a plugin requires parameters (e.g. file paths), specify them in the config
        # as a single string like "CADD,/path/to/file" so it becomes "--plugin CADD,/path/to/file".
        plugins_flags = " ".join(["--plugin %s" % p for p in config["parameters"]["vep"].get("plugins", [])])
    threads:
        config["resources"]["vep"]["threads"]
    resources:
        mem_mb = config["resources"]["vep"]["mem_mb"]
    log:
        os.path.join(OUTPUT_DIR, "logs", "vep", "{sample}.log")
    singularity:
        "docker://ensemblorg/ensembl-vep:release_110.1"
    shell:
        """
        # VEP 실행 (Singularity 컨테이너 내부)
        vep \
            --input_file {input.vcf} \
            --output_file {output.annotated_vcf} \
            --format vcf \
            --vcf \
            --compress_output bgzip \
            --fork {threads} \
            --species homo_sapiens \
            --assembly GRCh38 \
            --cache \
            --dir_cache {params.vep_cache} \
            --fasta {params.ref_genome} \
            --variant_class \
            --sift b \
            --polyphen b \
            --regulatory \
            --biotype \
            --canonical \
            --gene_phenotype \
            --pick_order {params.pick_order} \
            {params.plugins_flags} \
            --stats_file {output.stats} \
            --force_overwrite \
            2> {log}
        """

# VCF를 TSV로 변환 (bcftools는 호스트 conda 환경에서 실행)
rule vep_to_tsv:
    input:
        vcf = os.path.join(OUTPUT_DIR, "small_variants", "{sample}.annotated.vcf.gz")
    output:
        tsv = os.path.join(OUTPUT_DIR, "small_variants", "{sample}.annotated.tsv")
    log:
        os.path.join(OUTPUT_DIR, "logs", "vep", "{sample}_tsv.log")
    conda:
        "environment.yaml"
    shell:
        """
        bcftools query \
            -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO/CSQ\\n' \
            {input.vcf} \
            > {output.tsv} \
            2> {log}
        """

# ================================================================================
# Phase 3: Structural Variant 필터링 (svpack)
# ================================================================================

rule filter_sv:
    input:
        vcf = lambda wildcards: SAMPLE_FILES[wildcards.sample]['phased_sv_vcf']
    output:
        filtered_sv = os.path.join(OUTPUT_DIR, "structural_variants", "{sample}.filtered_sv.vcf.gz"),
        filtered_sv_tbi = os.path.join(OUTPUT_DIR, "structural_variants", "{sample}.filtered_sv.vcf.gz.tbi")
    params:
        min_sv_size = config["parameters"]["svpack"]["min_sv_size"],
        pass_only_flag = lambda wildcards: "--pass-only" if config["parameters"]["svpack"]["filter_pass_only"] else ""
    threads:
        config["resources"]["svpack"]["threads"]
    resources:
        mem_mb = config["resources"]["svpack"]["mem_mb"]
    log:
        os.path.join(OUTPUT_DIR, "logs", "svpack", "{sample}.log")
    conda:
        "environment.yaml"
    shell:
        """
        # svpack을 사용한 SV 필터링 (프로젝트 내 스크립트)
        # svpack filter는 표준 출력으로 결과를 출력하므로 리디렉션 필요
        python scripts/svpack filter \
            --min-svlen {params.min_sv_size} \
            {params.pass_only_flag} \
            {input.vcf} \
            2> {log} \
            | bgzip -c > {output.filtered_sv}
        
        # 인덱스 생성
        tabix -p vcf {output.filtered_sv} 2>> {log}
        """

# ================================================================================
# Phase 4: DMR 분석 준비 (CpG BED 파일 수집)
# ================================================================================

rule prepare_methylation_data:
    input:
        control_beds = lambda wildcards: [
            SAMPLE_FILES[sample]['cpg_combined_bed']
            for sample in config["samples"]["control"]
        ],
        experimental_beds = lambda wildcards: [
            SAMPLE_FILES[sample]['cpg_combined_bed']
            for sample in config["samples"]["experimental"]
        ]
    output:
        control_list = os.path.join(OUTPUT_DIR, "dmr_analysis", "control_samples.txt"),
        experimental_list = os.path.join(OUTPUT_DIR, "dmr_analysis", "experimental_samples.txt")
    run:
        # Control 샘플 목록 작성
        with open(output.control_list, 'w') as f:
            for bed_file in input.control_beds:
                f.write(f"{bed_file}\n")
        
        # Experimental 샘플 목록 작성
        with open(output.experimental_list, 'w') as f:
            for bed_file in input.experimental_beds:
                f.write(f"{bed_file}\n")

# ================================================================================
# Phase 5: DMR 분석 실행 (DSS in R)
# ================================================================================

rule run_dmr_analysis:
    input:
        control_list = os.path.join(OUTPUT_DIR, "dmr_analysis", "control_samples.txt"),
        experimental_list = os.path.join(OUTPUT_DIR, "dmr_analysis", "experimental_samples.txt"),
        r_script = "scripts/dmr_analysis.R"
    output:
        dmr_results = os.path.join(OUTPUT_DIR, "dmr_analysis", "dmr_results.csv"),
        dmr_plots = os.path.join(OUTPUT_DIR, "dmr_analysis", "dmr_plots.pdf")
    params:
        p_value = config["parameters"]["dmr"]["p_value_cutoff"],
        min_diff = config["parameters"]["dmr"]["min_methylation_diff"],
        min_cpg = config["parameters"]["dmr"]["min_cpg_sites"],
        smoothing = config["parameters"]["dmr"]["smoothing_span"]
    threads:
        config["resources"]["dmr"]["threads"]
    resources:
        mem_mb = config["resources"]["dmr"]["mem_mb"]
    log:
        os.path.join(OUTPUT_DIR, "logs", "dmr", "dmr_analysis.log")
    conda:
        "envs/dmr.yaml"
    shell:
        """
        Rscript {input.r_script} \
            --control {input.control_list} \
            --experimental {input.experimental_list} \
            --output-csv {output.dmr_results} \
            --output-pdf {output.dmr_plots} \
            --pvalue {params.p_value} \
            --min-diff {params.min_diff} \
            --min-cpg {params.min_cpg} \
            --smoothing {params.smoothing} \
            2> {log}
        """

# ================================================================================
# Phase 6: 메틸화 BigWig 파일 정리
# ================================================================================

rule merge_methylation_plots:
    input:
        bigwigs = lambda wildcards: [
            SAMPLE_FILES[sample]['cpg_combined_bw']
            for sample in ALL_SAMPLES
        ]
    output:
        summary = os.path.join(OUTPUT_DIR, "methylation", "methylation_summary.txt")
    run:
        # BigWig 파일 정리 및 요약
        import shutil
        
        os.makedirs(os.path.join(OUTPUT_DIR, "methylation"), exist_ok=True)
        
        with open(output.summary, 'w') as f:
            f.write("Sample\tBigWig_File\tFile_Size_MB\n")
            
            for sample in ALL_SAMPLES:
                bw_file = SAMPLE_FILES[sample]['cpg_combined_bw']
                dest_file = os.path.join(OUTPUT_DIR, "methylation", f"{sample}.bw")
                
                # BigWig 파일 복사
                if os.path.exists(bw_file):
                    shutil.copy2(bw_file, dest_file)
                    file_size = os.path.getsize(dest_file) / (1024 * 1024)  # MB
                    f.write(f"{sample}\t{dest_file}\t{file_size:.2f}\n")

# ================================================================================
# 유틸리티 규칙: 분석 요약 리포트 생성
# ================================================================================

# ================================================================================
# Phase 7: SV 유전자 영향 주석 (svpack consequence)
# ================================================================================

rule annotate_sv_consequence:
    input:
        filtered_sv = os.path.join(OUTPUT_DIR, "structural_variants", "{sample}.filtered_sv.vcf.gz"),
        gff3 = config["paths"]["gff3_file"]
    output:
        annotated_vcf = os.path.join(OUTPUT_DIR, "structural_variants", "{sample}.annotated_sv.vcf.gz"),
        annotated_tsv = os.path.join(OUTPUT_DIR, "structural_variants", "{sample}.annotated_sv.tsv")
    threads:
        config["resources"].get("sv_consequence", {}).get("threads", 2)
    resources:
        mem_mb = config["resources"].get("sv_consequence", {}).get("mem_mb", 4000)
    log:
        os.path.join(OUTPUT_DIR, "logs", "sv_consequence", "{sample}.log")
    conda:
        "environment.yaml"
    shell:
        """
        mkdir -p $(dirname {log})

        # svpack consequence로 유전자 영향 주석 추가
        bcftools view {input.filtered_sv} 2>> {log} \
            | python scripts/svpack consequence {input.gff3} 2>> {log} \
            | bgzip -c > {output.annotated_vcf}

        bcftools index -t {output.annotated_vcf} 2>> {log}

        # BCSQ 필드 파싱하여 TSV 생성
        bcftools query \
            -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/BCSQ\t[%GT]\n' \
            -H \
            {output.annotated_vcf} \
            > {output.annotated_tsv} \
            2>> {log}

        echo "SV consequence annotation completed" >> {log}
        """

# ================================================================================
# Phase 8: TRGT 반복서열 이상값 탐지
# ================================================================================

rule analyze_trgt:
    input:
        trgt_vcf = lambda wildcards: SAMPLE_FILES[wildcards.sample]['phased_trgt_vcf']
    output:
        outliers_tsv = os.path.join(OUTPUT_DIR, "trgt", "{sample}.trgt_outliers.tsv"),
        summary_tsv  = os.path.join(OUTPUT_DIR, "trgt", "{sample}.trgt_summary.tsv")
    params:
        min_sd = config["parameters"].get("trgt", {}).get("min_spanning_depth", 5),
        expansion_buffer = config["parameters"].get("trgt", {}).get("expansion_buffer", 0)
    threads:
        config["resources"].get("trgt", {}).get("threads", 2)
    resources:
        mem_mb = config["resources"].get("trgt", {}).get("mem_mb", 4000)
    log:
        os.path.join(OUTPUT_DIR, "logs", "trgt", "{sample}.log")
    conda:
        "environment.yaml"
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.outliers_tsv})

        python scripts/trgt_outlier.py \
            --input {input.trgt_vcf} \
            --output-outliers {output.outliers_tsv} \
            --output-summary {output.summary_tsv} \
            --min-sd {params.min_sd} \
            --expansion-buffer {params.expansion_buffer} \
            2> {log}

        echo "TRGT analysis completed" >> {log}
        """

# ================================================================================
# Phase 9: Allele-specific methylation (ASM) 분석
# ================================================================================

rule analyze_asm:
    input:
        hap1_bed = lambda wildcards: SAMPLE_FILES[wildcards.sample]['cpg_hap1_bed'],
        hap2_bed = lambda wildcards: SAMPLE_FILES[wildcards.sample]['cpg_hap2_bed'],
        r_script = "scripts/asm_analysis.R"
    output:
        asm_csv  = os.path.join(OUTPUT_DIR, "asm", "{sample}.asm_results.csv"),
        asm_pdf  = os.path.join(OUTPUT_DIR, "asm", "{sample}.asm_plots.pdf")
    params:
        p_value  = config["parameters"].get("asm", {}).get("p_value_cutoff", 0.05),
        min_diff = config["parameters"].get("asm", {}).get("min_methylation_diff", 0.2),
        min_cpg  = config["parameters"].get("asm", {}).get("min_cpg_sites", 3),
        smoothing = config["parameters"].get("asm", {}).get("smoothing_span", 500)
    threads:
        config["resources"].get("asm", {}).get("threads", 4)
    resources:
        mem_mb = config["resources"].get("asm", {}).get("mem_mb", 8000)
    log:
        os.path.join(OUTPUT_DIR, "logs", "asm", "{sample}.log")
    conda:
        "envs/dmr.yaml"
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.asm_csv})

        Rscript {input.r_script} \
            --hap1 {input.hap1_bed} \
            --hap2 {input.hap2_bed} \
            --sample-name {wildcards.sample} \
            --output-csv {output.asm_csv} \
            --output-pdf {output.asm_pdf} \
            --pvalue {params.p_value} \
            --min-diff {params.min_diff} \
            --min-cpg {params.min_cpg} \
            --smoothing {params.smoothing} \
            2> {log}

        echo "ASM analysis completed" >> {log}
        """

# ================================================================================
# 유틸리티 규칙: 분석 요약 리포트 생성
# ================================================================================

rule generate_summary_report:
    input:
        small_variants = expand(
            os.path.join(OUTPUT_DIR, "small_variants", "{sample}.annotated.tsv"),
            sample=ALL_SAMPLES
        ),
        sv_variants = expand(
            os.path.join(OUTPUT_DIR, "structural_variants", "{sample}.filtered_sv.vcf.gz"),
            sample=ALL_SAMPLES
        ),
        sv_annotated = expand(
            os.path.join(OUTPUT_DIR, "structural_variants", "{sample}.annotated_sv.tsv"),
            sample=ALL_SAMPLES
        ),
        dmr_results = os.path.join(OUTPUT_DIR, "dmr_analysis", "dmr_results.csv")
    output:
        report = os.path.join(OUTPUT_DIR, "analysis_summary_report.html")
    run:
        import datetime
        import csv as csv_module

        def count_lines(filepath):
            try:
                with open(filepath, 'r') as f:
                    return max(0, sum(1 for _ in f) - 1)  # 헤더 제외
            except Exception:
                return 0

        def read_dmr_summary(csv_path):
            try:
                with open(csv_path, 'r') as f:
                    reader = csv_module.DictReader(f)
                    rows = list(reader)
                    if not rows:
                        return 0, 0, 0
                    total = len(rows)
                    hyper = sum(1 for r in rows if float(r.get('diff_methy', 0) or 0) > 0)
                    hypo  = total - hyper
                    return total, hyper, hypo
            except Exception:
                return 0, 0, 0

        def collect_sv_genes(tsv_files):
            genes = set()
            for f in tsv_files:
                try:
                    with open(f, 'r') as fh:
                        for line in fh:
                            if line.startswith('#') or not line.strip():
                                continue
                            parts = line.strip().split('\t')
                            if len(parts) >= 6:
                                bcsq = parts[5]
                                for entry in bcsq.split(','):
                                    fields = entry.split('|')
                                    if len(fields) >= 2:
                                        genes.add(fields[1])
                except Exception:
                    pass
            genes.discard('.')
            genes.discard('')
            return sorted(genes)[:20]

        def collect_trgt_outliers():
            rows = []
            for s in ALL_SAMPLES:
                path = os.path.join(OUTPUT_DIR, "trgt", f"{s}.trgt_outliers.tsv")
                if not os.path.exists(path):
                    continue
                try:
                    with open(path, 'r') as fh:
                        reader = csv_module.DictReader(fh, delimiter='\t')
                        for r in reader:
                            r['sample'] = s
                            rows.append(r)
                except Exception:
                    pass
            return rows[:20]

        def collect_asm_counts():
            result = {}
            for s in ALL_SAMPLES:
                path = os.path.join(OUTPUT_DIR, "asm", f"{s}.asm_results.csv")
                result[s] = count_lines(path) if os.path.exists(path) else None
            return result

        dmr_total, dmr_hyper, dmr_hypo = read_dmr_summary(input.dmr_results)
        sv_genes = collect_sv_genes(input.sv_annotated)
        trgt_outliers = collect_trgt_outliers()
        asm_counts = collect_asm_counts()
        analysis_date = datetime.date.today().isoformat()

        trgt_rows_html = ""
        for r in trgt_outliers:
            trgt_rows_html += (
                f"<tr><td>{r.get('sample','')}</td><td>{r.get('trid','')}</td>"
                f"<td>{r.get('motif','')}</td><td>{r.get('allele1_len','')}</td>"
                f"<td>{r.get('allele2_len','')}</td><td>{r.get('expected_range','')}</td>"
                f"<td><b>{r.get('status','')}</b></td></tr>\n"
            )

        asm_rows_html = ""
        for s, cnt in asm_counts.items():
            display = str(cnt) if cnt is not None else "미분석"
            asm_rows_html += f"<tr><td>{s}</td><td>{display}</td></tr>\n"

        sv_genes_html = ", ".join(sv_genes) if sv_genes else "정보 없음"
        max_af  = config['parameters']['slivar']['max_af']
        min_gq  = config['parameters']['slivar']['min_gq']
        min_dp  = config['parameters']['slivar']['min_dp']
        n_small = len(input.small_variants)
        n_sv    = len(input.sv_variants)

        # TRGT 섹션 HTML (중첩 f-string 회피를 위해 사전 조합)
        if trgt_outliers:
            trgt_section = (
                "<table>\n"
                "<tr><th>샘플</th><th>Repeat ID</th><th>Motif</th>"
                "<th>Allele1 길이</th><th>Allele2 길이</th>"
                "<th>기대 범위</th><th>상태</th></tr>\n"
                + trgt_rows_html
                + "</table>\n"
                + f"<p>출력 위치: <code>{OUTPUT_DIR}/trgt/</code></p>"
            )
        else:
            trgt_section = "<p>TRGT 데이터 없음 (filelist.csv에 phased_trgt_vcf 미포함)</p>"

        # ASM 섹션 HTML
        if any(v is not None for v in asm_counts.values()):
            asm_section = (
                "<table>\n"
                "<tr><th>샘플</th><th>ASM 영역 수</th></tr>\n"
                + asm_rows_html
                + "</table>\n"
                + f"<p>출력 위치: <code>{OUTPUT_DIR}/asm/</code></p>"
            )
        else:
            asm_section = "<p>ASM 데이터 없음 (filelist.csv에 cpg_hap1_bed/cpg_hap2_bed 미포함)</p>"

        html = (
            "<!DOCTYPE html>\n"
            '<html lang="ko">\n'
            "<head>\n"
            '    <meta charset="UTF-8">\n'
            "    <title>PacBio HiFi WGS 3차 분석 리포트</title>\n"
            "    <style>\n"
            "        body { font-family: Arial, sans-serif; margin: 20px; background: #f9f9f9; }\n"
            "        h1 { color: #2c3e50; }\n"
            "        h2 { color: #34495e; border-bottom: 2px solid #3498db; padding-bottom: 4px; }\n"
            "        table { border-collapse: collapse; width: 100%; margin: 12px 0; background: white; }\n"
            "        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }\n"
            "        th { background-color: #3498db; color: white; }\n"
            "        tr:nth-child(even) { background-color: #f2f2f2; }\n"
            "        .stat { display: inline-block; margin: 8px 16px 8px 0; }\n"
            "        .stat-num { font-size: 2em; font-weight: bold; color: #3498db; }\n"
            "    </style>\n"
            "</head>\n"
            "<body>\n"
            f"    <h1>PacBio HiFi WGS 3차 분석 결과 요약</h1>\n"
            f"    <p>분석 일자: {analysis_date} &nbsp;|&nbsp; 총 샘플 수: <b>{len(ALL_SAMPLES)}</b></p>\n"
            "\n"
            "    <h2>1. Small Variant 분석</h2>\n"
            f'    <div class="stat"><div class="stat-num">{n_small}</div>샘플 주석 완료</div>\n'
            f"    <p>필터링 기준: AF &lt; {max_af}, GQ &ge; {min_gq}, DP &ge; {min_dp}</p>\n"
            f"    <p>출력 위치: <code>{OUTPUT_DIR}/small_variants/</code></p>\n"
            "\n"
            "    <h2>2. Structural Variant 분석</h2>\n"
            f'    <div class="stat"><div class="stat-num">{n_sv}</div>샘플 필터링 완료</div>\n'
            "    <p><b>유전자 영향 주석 (svpack consequence)</b> — 영향 유전자 상위 20개:</p>\n"
            f"    <p>{sv_genes_html}</p>\n"
            f"    <p>출력 위치: <code>{OUTPUT_DIR}/structural_variants/</code></p>\n"
            "\n"
            "    <h2>3. DMR 분석</h2>\n"
            f'    <div class="stat"><div class="stat-num">{dmr_total}</div>총 DMR</div>\n'
            f'    <div class="stat"><div class="stat-num" style="color:#e74c3c">{dmr_hyper}</div>Hyper-메틸화</div>\n'
            f'    <div class="stat"><div class="stat-num" style="color:#2980b9">{dmr_hypo}</div>Hypo-메틸화</div>\n'
            f"    <p>출력 위치: <code>{OUTPUT_DIR}/dmr_analysis/</code></p>\n"
            "\n"
            "    <h2>4. TRGT 반복서열 이상값</h2>\n"
            + trgt_section + "\n"
            "\n"
            "    <h2>5. Allele-specific Methylation (ASM)</h2>\n"
            + asm_section + "\n"
            "\n"
            "    <h2>파일 위치 요약</h2>\n"
            "    <ul>\n"
            f"        <li>Small Variants: <code>{OUTPUT_DIR}/small_variants/</code></li>\n"
            f"        <li>Structural Variants: <code>{OUTPUT_DIR}/structural_variants/</code></li>\n"
            f"        <li>DMR 분석: <code>{OUTPUT_DIR}/dmr_analysis/</code></li>\n"
            f"        <li>메틸화 데이터: <code>{OUTPUT_DIR}/methylation/</code></li>\n"
            f"        <li>TRGT 반복서열: <code>{OUTPUT_DIR}/trgt/</code></li>\n"
            f"        <li>ASM 분석: <code>{OUTPUT_DIR}/asm/</code></li>\n"
            "    </ul>\n"
            "</body>\n"
            "</html>"
        )
        with open(output.report, 'w') as f:
            f.write(html)

# ================================================================================
# Phase 10: VEP canonical impact 추출 (cohort burden 분석용 전처리)
# ================================================================================

rule extract_vep_canonical:
    input:
        vcf = os.path.join(OUTPUT_DIR, "small_variants", "{sample}.annotated.vcf.gz")
    output:
        tsv = os.path.join(OUTPUT_DIR, "small_variants", "{sample}.canonical_impacts.tsv")
    log:
        os.path.join(OUTPUT_DIR, "logs", "vep", "{sample}_canonical.log")
    conda:
        "environment.yaml"
    shell:
        """
        mkdir -p $(dirname {log})

        # bcftools +split-vep로 canonical 전사체 SYMBOL/IMPACT 추출 시도
        if bcftools +split-vep --version &>/dev/null 2>&1 || \
           bcftools +split-vep {input.vcf} --help &>/dev/null 2>&1; then
            bcftools +split-vep {input.vcf} \
                -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%IMPACT\t%Consequence\t[%GT]\n' \
                -d -A tab -s canonical \
                > {output.tsv} 2>> {log}
        else
            # fallback: Python 스크립트로 CSQ 파싱
            python scripts/extract_vep_fields.py \
                --input {input.vcf} \
                --output {output.tsv} \
                2>> {log}
        fi

        echo "VEP canonical extraction completed" >> {log}
        """

# ================================================================================
# Phase 11: Small Variant burden test (군 간 유전자별 비교)
# ================================================================================

rule cohort_variant_burden:
    input:
        tsvs = expand(
            os.path.join(OUTPUT_DIR, "small_variants", "{sample}.canonical_impacts.tsv"),
            sample=ALL_SAMPLES
        ),
        r_script = "scripts/variant_burden.R"
    output:
        stats_csv = os.path.join(OUTPUT_DIR, "cohort", "variant_burden_stats.csv"),
        plots_pdf = os.path.join(OUTPUT_DIR, "cohort", "variant_burden_plots.pdf")
    params:
        input_dir   = os.path.join(OUTPUT_DIR, "small_variants"),
        groups_json = json.dumps(GROUPS),
        min_impact  = lambda wildcards: config["parameters"].get("cohort", {}).get("min_impact", "HIGH|MODERATE"),
        fdr_cutoff  = lambda wildcards: config["parameters"].get("cohort", {}).get("fdr_cutoff", 0.05)
    threads:
        config["resources"].get("cohort_burden", {}).get("threads", 4)
    resources:
        mem_mb = config["resources"].get("cohort_burden", {}).get("mem_mb", 8000)
    log:
        os.path.join(OUTPUT_DIR, "logs", "cohort", "variant_burden.log")
    conda:
        "envs/dmr.yaml"
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.stats_csv})

        Rscript {input.r_script} \
            --input-dir '{params.input_dir}' \
            --groups '{params.groups_json}' \
            --min-impact '{params.min_impact}' \
            --fdr-cutoff {params.fdr_cutoff} \
            --output-csv {output.stats_csv} \
            --output-pdf {output.plots_pdf} \
            2> {log}

        echo "Variant burden analysis completed" >> {log}
        """

# ================================================================================
# Phase 12: SV burden test (군 간 유전자별 SV 비교)
# ================================================================================

rule cohort_sv_burden:
    input:
        tsvs = expand(
            os.path.join(OUTPUT_DIR, "structural_variants", "{sample}.annotated_sv.tsv"),
            sample=ALL_SAMPLES
        ),
        r_script = "scripts/sv_burden.R"
    output:
        stats_csv = os.path.join(OUTPUT_DIR, "cohort", "sv_burden_stats.csv"),
        plots_pdf = os.path.join(OUTPUT_DIR, "cohort", "sv_burden_plots.pdf")
    params:
        input_dir    = os.path.join(OUTPUT_DIR, "structural_variants"),
        groups_json  = json.dumps(GROUPS),
        sv_types     = lambda wildcards: json.dumps(
            config["parameters"].get("cohort", {}).get("sv_consequence_types", ["sv:cds", "sv:utr"])
        ),
        fdr_cutoff   = lambda wildcards: config["parameters"].get("cohort", {}).get("fdr_cutoff", 0.05)
    threads:
        config["resources"].get("cohort_burden", {}).get("threads", 4)
    resources:
        mem_mb = config["resources"].get("cohort_burden", {}).get("mem_mb", 8000)
    log:
        os.path.join(OUTPUT_DIR, "logs", "cohort", "sv_burden.log")
    conda:
        "envs/dmr.yaml"
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.stats_csv})

        Rscript {input.r_script} \
            --input-dir '{params.input_dir}' \
            --groups '{params.groups_json}' \
            --sv-types '{params.sv_types}' \
            --fdr-cutoff {params.fdr_cutoff} \
            --output-csv {output.stats_csv} \
            --output-pdf {output.plots_pdf} \
            2> {log}

        echo "SV burden analysis completed" >> {log}
        """

# ================================================================================
# Phase 13: TRGT 반복서열 군 간 비교
# ================================================================================

rule trgt_group_compare:
    input:
        summaries = [
            os.path.join(OUTPUT_DIR, "trgt", f"{s}.trgt_summary.tsv")
            for s in ALL_SAMPLES
            if SAMPLE_FILES.get(s, {}).get('phased_trgt_vcf')
        ],
        r_script = "scripts/trgt_group_compare.R"
    output:
        stats_csv = os.path.join(OUTPUT_DIR, "cohort", "trgt_group_compare.csv"),
        plots_pdf = os.path.join(OUTPUT_DIR, "cohort", "trgt_group_compare_plots.pdf")
    params:
        input_dir   = os.path.join(OUTPUT_DIR, "trgt"),
        groups_json = json.dumps(GROUPS),
        fdr_cutoff  = lambda wildcards: config["parameters"].get("cohort", {}).get("fdr_cutoff", 0.05)
    threads:
        config["resources"].get("trgt_compare", {}).get("threads", 2)
    resources:
        mem_mb = config["resources"].get("trgt_compare", {}).get("mem_mb", 4000)
    log:
        os.path.join(OUTPUT_DIR, "logs", "cohort", "trgt_group_compare.log")
    conda:
        "envs/dmr.yaml"
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.stats_csv})

        Rscript {input.r_script} \
            --input-dir '{params.input_dir}' \
            --groups '{params.groups_json}' \
            --fdr-cutoff {params.fdr_cutoff} \
            --output-csv {output.stats_csv} \
            --output-pdf {output.plots_pdf} \
            2> {log}

        echo "TRGT group comparison completed" >> {log}
        """

# ================================================================================
# Phase 14: DMR pairwise — 모든 군 쌍에 대해 자동 실행
# ================================================================================

rule prepare_group_methylation_list:
    input:
        beds = lambda wildcards: [
            SAMPLE_FILES[s]['cpg_combined_bed']
            for s in GROUPS[wildcards.group]
            if SAMPLE_FILES.get(s, {}).get('cpg_combined_bed')
        ]
    output:
        group_list = os.path.join(OUTPUT_DIR, "cohort", "methylation_{group}.txt")
    run:
        os.makedirs(os.path.dirname(output.group_list), exist_ok=True)
        with open(output.group_list, 'w') as f:
            for bed in input.beds:
                f.write(bed + '\n')

rule dmr_pairwise:
    input:
        g1_list  = os.path.join(OUTPUT_DIR, "cohort", "methylation_{group1}.txt"),
        g2_list  = os.path.join(OUTPUT_DIR, "cohort", "methylation_{group2}.txt"),
        r_script = "scripts/dmr_analysis.R"
    output:
        csv = os.path.join(OUTPUT_DIR, "cohort", "dmr_{group1}_vs_{group2}.csv"),
        pdf = os.path.join(OUTPUT_DIR, "cohort", "dmr_{group1}_vs_{group2}_plots.pdf")
    params:
        p_value   = config["parameters"]["dmr"]["p_value_cutoff"],
        min_diff  = config["parameters"]["dmr"]["min_methylation_diff"],
        min_cpg   = config["parameters"]["dmr"]["min_cpg_sites"],
        smoothing = config["parameters"]["dmr"]["smoothing_span"]
    threads:
        config["resources"]["dmr"]["threads"]
    resources:
        mem_mb = config["resources"]["dmr"]["mem_mb"]
    log:
        os.path.join(OUTPUT_DIR, "logs", "cohort", "dmr_{group1}_vs_{group2}.log")
    conda:
        "envs/dmr.yaml"
    shell:
        """
        mkdir -p $(dirname {log})

        Rscript {input.r_script} \
            --control {input.g1_list} \
            --experimental {input.g2_list} \
            --output-csv {output.csv} \
            --output-pdf {output.pdf} \
            --pvalue {params.p_value} \
            --min-diff {params.min_diff} \
            --min-cpg {params.min_cpg} \
            --smoothing {params.smoothing} \
            2> {log}

        echo "DMR pairwise ({wildcards.group1} vs {wildcards.group2}) completed" >> {log}
        """
