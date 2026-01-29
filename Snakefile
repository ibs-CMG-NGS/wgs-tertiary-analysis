# ================================================================================
# PacBio HiFi WGS 3차 분석 파이프라인 (Snakemake)
# ================================================================================

import os
import csv
from pathlib import Path

# 설정 파일 로드
# 환경 변수 SNAKEMAKE_CONFIG로 지정 가능, 기본값: config/config.yaml
CONFIG_FILE = os.environ.get("SNAKEMAKE_CONFIG", "config/config.yaml")
configfile: CONFIG_FILE

# 전체 샘플 목록 생성
ALL_SAMPLES = config["samples"]["control"] + config["samples"]["experimental"]

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
            }
    print(f"✓ filelist.csv 로드 완료: {len(SAMPLE_FILES)}개 샘플")
else:
    # filelist.csv가 없으면 기존 방식 (wdl_out_dir 사용)
    print(f"⚠️  경고: {filelist_csv}를 찾을 수 없습니다.")
    print(f"   'python bin/prepare_inputs.py'를 먼저 실행하거나")
    print(f"   기존 방식을 사용하려면 paths.wdl_out_dir을 설정하세요.")
    
    # 기존 방식으로 폴백
    for sample in ALL_SAMPLES:
        SAMPLE_FILES[sample] = {
            'phased_small_variant_vcf': os.path.join(
                config["paths"]["wdl_out_dir"],
                f"phased_small_variant_{sample}.vcf.gz"
            ),
            'phased_sv_vcf': os.path.join(
                config["paths"]["wdl_out_dir"],
                f"phased_sv_{sample}.vcf.gz"
            ),
            'cpg_combined_bed': os.path.join(
                config["paths"]["wdl_out_dir"],
                f"cpg_combined_{sample}.bed"
            ),
            'cpg_combined_bw': os.path.join(
                config["paths"]["wdl_out_dir"],
                f"cpg_combined_{sample}.bw"
            ),
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
        
        # Structural Variant 필터링
        expand(
            os.path.join(OUTPUT_DIR, "structural_variants", "{sample}.filtered_sv.vcf.gz"),
            sample=ALL_SAMPLES
        ),
        
        # DMR 분석 결과
        os.path.join(OUTPUT_DIR, "dmr_analysis", "dmr_results.csv"),
        os.path.join(OUTPUT_DIR, "dmr_analysis", "dmr_plots.pdf"),
        
        # 메틸화 파일 정리
        os.path.join(OUTPUT_DIR, "methylation", "methylation_summary.txt")

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
        max_af = config["parameters"]["slivar"]["max_af"],
        min_gq = config["parameters"]["slivar"]["min_gq"],
        min_dp = config["parameters"]["slivar"]["min_dp"],
        hpo_terms = config["parameters"]["slivar"]["hpo_terms"]
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
        # slivar를 사용한 필터링
        slivar expr \
            --vcf {input.vcf} \
            --pass-only \
            --out-vcf {output.filtered_vcf} \
            --info 'INFO.AF < {params.max_af}' \
            --sample-expr 'sample.GQ >= {params.min_gq} && sample.DP >= {params.min_dp}' \
            2> {log}
        
        # TSV 변환 (연구자용)
        slivar tsv \
            --vcf {output.filtered_vcf} \
            --csq-column CSQ \
            --info-field AF \
            --sample-field DP \
            --sample-field GQ \
            > {output.filtered_tsv} \
            2>> {log}
        """

# ================================================================================
# Phase 2: VEP 주석 추가
# ================================================================================

rule annotate_vep:
    input:
        vcf = os.path.join(OUTPUT_DIR, "small_variants", "{sample}.filtered.vcf.gz")
    output:
        annotated_vcf = os.path.join(OUTPUT_DIR, "small_variants", "{sample}.annotated.vcf.gz"),
        annotated_tsv = os.path.join(OUTPUT_DIR, "small_variants", "{sample}.annotated.tsv"),
        stats = os.path.join(OUTPUT_DIR, "small_variants", "{sample}.vep_stats.html")
    params:
        ref_genome = config["paths"]["ref_genome"],
        vep_cache = config["paths"]["vep_cache"],
        pick_order = config["parameters"]["vep"]["pick_order"],
        plugins = ",".join(config["parameters"]["vep"]["plugins"])
    threads:
        config["resources"]["vep"]["threads"]
    resources:
        mem_mb = config["resources"]["vep"]["mem_mb"]
    log:
        os.path.join(OUTPUT_DIR, "logs", "vep", "{sample}.log")
    conda:
        "environment.yaml"
    shell:
        """
        # VEP 실행
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
            --everything \
            --pick_order {params.pick_order} \
            --plugin {params.plugins} \
            --stats_file {output.stats} \
            --force_overwrite \
            2> {log}
        
        # VCF를 TSV로 변환 (연구자용)
        bcftools query \
            -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO/CSQ\\n' \
            {output.annotated_vcf} \
            > {output.annotated_tsv} \
            2>> {log}
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
        python scripts/svpack filter \
            --input {input.vcf} \
            --output {output.filtered_sv} \
            --min-svlen {params.min_sv_size} \
            {params.pass_only_flag} \
            2> {log}
        
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
        "environment.yaml"
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
        dmr_results = os.path.join(OUTPUT_DIR, "dmr_analysis", "dmr_results.csv")
    output:
        report = os.path.join(OUTPUT_DIR, "analysis_summary_report.html")
    run:
        # 간단한 HTML 리포트 생성
        with open(output.report, 'w') as f:
            f.write("""
<!DOCTYPE html>
<html>
<head>
    <title>PacBio HiFi WGS 3차 분석 리포트</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        h1 {{ color: #2c3e50; }}
        h2 {{ color: #34495e; border-bottom: 2px solid #3498db; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #3498db; color: white; }}
        tr:nth-child(even) {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
    <h1>PacBio HiFi WGS 3차 분석 결과 요약</h1>
    <p>분석 일자: {date}</p>
    
    <h2>1. Small Variant 분석</h2>
    <p>총 샘플 수: {total_samples}</p>
    <ul>
        <li>필터링된 변이: {small_variant_files}개 파일</li>
        <li>VEP 주석 완료: {annotated_files}개 파일</li>
    </ul>
    
    <h2>2. Structural Variant 분석</h2>
    <p>총 SV 파일 수: {sv_files}개</p>
    
    <h2>3. DMR 분석</h2>
    <p>DMR 결과: {dmr_file}</p>
    
    <h2>파일 위치</h2>
    <ul>
        <li>Small Variants: {output_dir}/small_variants/</li>
        <li>Structural Variants: {output_dir}/structural_variants/</li>
        <li>DMR 분석: {output_dir}/dmr_analysis/</li>
        <li>메틸화 데이터: {output_dir}/methylation/</li>
    </ul>
</body>
</html>
            """.format(
                date="2026-01-29",
                total_samples=len(ALL_SAMPLES),
                small_variant_files=len(input.small_variants),
                annotated_files=len(input.small_variants),
                sv_files=len(input.sv_variants),
                dmr_file=input.dmr_results,
                output_dir=OUTPUT_DIR
            ))
