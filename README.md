# PacBio HiFi WGS 3차 분석 파이프라인

PacBio WDL 파이프라인(singleton.wdl) 실행 결과를 입력받아 사용자 정의 필터링, 심화 주석, 후성유전학 분석, 반복서열 분석, 그리고 **다중 군 코호트 통계 비교**를 수행하는 Snakemake 기반의 자동화 파이프라인입니다.

## 📋 목차

- [개요](#개요)
- [프로젝트 구조](#프로젝트-구조)
- [시스템 요구사항](#시스템-요구사항)
- [설치](#설치)
- [설정](#설정)
- [사용법](#사용법)
- [출력 결과](#출력-결과)
- [다중 군 분석 (Phase 4/5)](#다중-군-분석-phase-45)
- [문제 해결](#문제-해결)

---

## 🎯 개요

본 파이프라인은 다음과 같은 분석 단계를 수행합니다:

### Phase 1: Small Variant 필터링 및 VEP 주석
- **질환 연관 변이 우선순위화**: `slivar`를 사용한 인구 집단 빈도 필터링 (gnomAD)
- **기능 예측 주석**: `VEP`를 통한 CADD, REVEL, 기능적 영향(IMPACT) 주석

### Phase 2: Structural Variant 필터링 및 유전자 영향 주석
- **인구 집단 빈도 필터링**: `svpack`을 사용한 SV 빈도 기반 필터링
- **유전자 영향 주석**: `svpack consequence` + Ensembl GFF3로 CDS/UTR 영향 주석 (`BCSQ` 필드)

### Phase 3: 후성유전학 심화 분석
- **DMR (Differential Methylation Region)**: `DSS`를 사용한 CpG 메틸화 차이 분석
- **ASM (Allele-Specific Methylation)**: Haplotype 1/2 간 메틸화 차이 분석 (각인 유전자 등)
- **메틸화 시각화**: IGV 호환 BigWig 파일 정리

### Phase 4: 반복서열 분석 (TRGT)
- **Expanded Repeat 탐지**: TRGT VCF에서 기대 범위를 벗어난 반복서열 이상값 탐지
- **전체 반복서열 요약**: 모든 PASS 반복서열의 allele 길이 요약 (군 간 비교 입력)

### Phase 5: 다중 군 코호트 통계 비교
- **Small Variant Burden Test**: 군별 유해 변이(HIGH/MODERATE IMPACT) 부담 유전자 비교
- **SV Burden Test**: 군별 유전자 파괴 SV 비교 (Fisher's exact / χ² + BH FDR)
- **TRGT 반복서열 군 간 비교**: 군별 allele 길이 차이 (Kruskal-Wallis + pairwise Wilcoxon)
- **DMR Pairwise 비교**: 모든 군 쌍에 대한 자동 DMR 비교

> **통계 주의사항**: 군당 2–3 샘플의 경우 검정력이 낮으므로, 결과는 **탐색적(exploratory)** 목적으로 해석하세요.

---

## 📁 프로젝트 구조

```
wgs-tertiary-pipeline/
│
├── Snakefile                    # 메인 워크플로우 (Phase 1–5)
├── README.md                    # 이 문서
├── environment.yaml             # Conda 환경 (기본 도구)
│
├── envs/
│   └── dmr.yaml                 # R/Bioconductor 환경 (DMR/ASM/Burden)
│
├── config/                      # 설정 파일
│   ├── config.yaml              # 실제 설정 (사용자 수정)
│   ├── config.example.yaml      # 설정 템플릿
│   ├── config.human.yaml        # Human GRCh38 설정 예시
│   └── cluster_config.yaml      # SLURM 리소스
│
├── bin/                         # 실행 스크립트
│   ├── prepare_inputs.py        # WDL 결과 파일 자동 탐색 → filelist.csv 생성
│   ├── run_pipeline.sh          # 메인 실행
│   ├── setup_check.sh           # 환경 검증
│   └── run_slurm.sh             # 클러스터 실행
│
├── scripts/                     # 분석 스크립트
│   ├── dmr_analysis.R           # DMR 분석 (DSS)
│   ├── asm_analysis.R           # Allele-specific methylation 분석
│   ├── trgt_outlier.py          # TRGT expanded repeat 탐지
│   ├── extract_vep_fields.py    # VEP canonical 전사체 IMPACT 추출
│   ├── variant_burden.R         # Small Variant burden test
│   ├── sv_burden.R              # SV burden test
│   ├── trgt_group_compare.R     # TRGT 반복서열 군 간 비교
│   └── svpack                   # PacBio svpack 도구
│
└── docs/                        # 추가 문서
    ├── QUICKSTART.md
    └── SUMMARY.md
```

---

## 💻 시스템 요구사항

### 필수 소프트웨어

| 도구 | 버전 | 용도 |
|------|------|------|
| Snakemake | ≥ 7.0 | 워크플로우 엔진 |
| Python | ≥ 3.8 | 전처리 스크립트 |
| bcftools | ≥ 1.17 | VCF 처리 |
| bgzip / tabix | htslib ≥ 1.17 | 인덱싱 |
| bedtools | ≥ 2.31 | 영역 교차 분석 |
| Docker 또는 Singularity | - | 컨테이너 실행 |

### 컨테이너 이미지

| 이미지 | 용도 |
|--------|------|
| `ensemblorg/ensembl-vep:release_110.1` | VEP 주석 |
| `bioconductor/bioconductor_docker:RELEASE_3_18` | R/DSS DMR·ASM·Burden 분석 |

> **참고**: `slivar`, `svpack`은 conda 환경 또는 로컬 설치 권장 (컨테이너 없이 동작).

### 하드웨어 권장사항

- **CPU**: 16 코어 이상
- **메모리**: 64 GB 이상 (DMR/VEP 동시 실행 시)
- **디스크**: 500 GB 이상

---

## 🔧 설치

### 1. Snakemake 설치

```bash
conda create -n snakemake_env -c conda-forge -c bioconda snakemake>=7.0
conda activate snakemake_env
```

### 2. 기본 Conda 환경 설치

```bash
conda env create -f environment.yaml
```

### 3. R/DMR 환경 설치

```bash
conda env create -f envs/dmr.yaml
```

### 4. 파이프라인 클론

```bash
git clone <repository-url>
cd wgs-tertiary-pipeline
```

---

## ⚙️ 설정

### 1. filelist.csv 생성 (필수)

`bin/prepare_inputs.py`를 실행하면 WDL 결과 디렉토리에서 필요한 파일을 자동으로 찾아 `filelist.csv`를 생성합니다:

```bash
python bin/prepare_inputs.py --config-file config/config.yaml
```

생성된 `filelist.csv`는 다음 컬럼을 포함합니다:

| 컬럼 | 설명 |
|------|------|
| `sample_id` | 샘플 ID |
| `phased_small_variant_vcf` | DeepVariant phased VCF |
| `phased_sv_vcf` | Sawfish phased SV VCF |
| `cpg_combined_bed` | 전체 haplotype CpG BED |
| `cpg_combined_bw` | CpG BigWig (IGV용) |
| `phased_trgt_vcf` | TRGT 반복서열 VCF (선택) |
| `cpg_hap1_bed` | Haplotype 1 CpG BED (선택) |
| `cpg_hap2_bed` | Haplotype 2 CpG BED (선택) |

> **선택 컬럼** (`phased_trgt_vcf`, `cpg_hap1/2_bed`)이 없으면 해당 분석(TRGT, ASM)은 자동으로 스킵됩니다.

### 2. config.yaml 편집

`config/config.example.yaml`을 복사하여 수정합니다:

```bash
cp config/config.example.yaml config/config.yaml
```

#### 기본 설정 (2군 비교)

```yaml
paths:
  batch_results_dir: "/data/hifi-wgs/batch_results"  # WDL 결과 최상위
  ref_genome:        "/data/references/GRCh38.fa"
  vep_cache:         "/data/vep_cache/homo_sapiens/110_GRCh38"
  gff3_file:         "/data/hifi-wdl-resources/GRCh38/ensembl.GRCh38.101.reformatted.gff3.gz"
  output_dir:        "./analysis_results"
  filelist_csv:      "filelist.csv"

samples:
  control:
    - ctrl_1
    - ctrl_2
  experimental:
    - treat_1
    - treat_2

parameters:
  slivar:
    max_af: 0.01
    min_gq: 20
    min_dp: 10
    hpo_terms: "HP:0001250,HP:0000707"
  dmr:
    p_value_cutoff: 0.05
    min_methylation_diff: 0.1
    min_cpg_sites: 3
```

#### 다중 군 비교 (3군 이상)

`samples.groups` 섹션을 추가하면 코호트 분석에서 모든 군 쌍을 자동 비교합니다:

```yaml
samples:
  control:       [ctrl_1, ctrl_2]       # 기존 DMR 2군 rule 호환용
  experimental:  [treat_1, treat_2]     # 기존 DMR 2군 rule 호환용
  groups:                               # 다중 군 정의 (있으면 이것 우선)
    control:
      - ctrl_1
      - ctrl_2
    treated_low:
      - treat_1
      - treat_2
    treated_high:
      - treat_3
      - treat_4

parameters:
  cohort:
    min_impact: "HIGH|MODERATE"          # Small Variant burden 필터
    sv_consequence_types:
      - "sv:cds"
      - "sv:utr"
    fdr_cutoff: 0.05
```

---

## 🚀 사용법

### 1. 환경 검증

```bash
bash bin/setup_check.sh --config config/config.yaml
```

### 2. WDL 결과 파일 탐색

```bash
python bin/prepare_inputs.py --config-file config/config.yaml
```

### 3. Dry-run

```bash
snakemake --dry-run --cores 1 --configfile config/config.yaml
```

### 4. 전체 파이프라인 실행

```bash
# Singularity 컨테이너 사용
snakemake --use-singularity --cores 16 --configfile config/config.yaml

# conda 환경 사용
snakemake --use-conda --cores 16 --configfile config/config.yaml
```

### 5. 특정 분석만 실행

```bash
# Small Variant burden test만
snakemake --use-conda --cores 8 --configfile config/config.yaml \
  analysis_results/cohort/variant_burden_stats.csv

# DMR pairwise (특정 군 쌍)
snakemake --use-conda --cores 4 --configfile config/config.yaml \
  analysis_results/cohort/dmr_control_vs_treated_low.csv

# TRGT 군 간 비교만
snakemake --use-conda --cores 4 --configfile config/config.yaml \
  analysis_results/cohort/trgt_group_compare.csv
```

### 6. 클러스터 환경 (SLURM)

```bash
snakemake \
  --use-singularity --use-conda \
  --cluster "sbatch -p normal -n {threads} --mem={resources.mem_mb}" \
  --jobs 20 --cores 128 \
  --configfile config/config.yaml
```

---

## 📊 출력 결과

### 전체 출력 구조

```
analysis_results/
│
├── small_variants/
│   ├── {sample}.filtered.vcf.gz          # slivar 필터링 결과
│   ├── {sample}.filtered.tsv
│   ├── {sample}.annotated.vcf.gz         # VEP 주석 추가
│   ├── {sample}.annotated.tsv
│   └── {sample}.canonical_impacts.tsv    # canonical 전사체 IMPACT 요약 (burden용)
│
├── structural_variants/
│   ├── {sample}.filtered_sv.vcf.gz       # svpack 필터링
│   ├── {sample}.annotated_sv.vcf.gz      # svpack consequence 주석
│   └── {sample}.annotated_sv.tsv         # BCSQ 유전자 영향 TSV
│
├── dmr_analysis/
│   ├── dmr_results.csv                   # Control vs Experimental DMR
│   └── dmr_plots.pdf                     # DMR 시각화 (5종)
│
├── trgt/
│   ├── {sample}.trgt_outliers.tsv        # Expanded repeat 이상값
│   └── {sample}.trgt_summary.tsv         # 전체 PASS repeat 요약 (군 비교용)
│
├── asm/
│   ├── {sample}.asm_results.csv          # Allele-specific methylation DMR
│   └── {sample}.asm_plots.pdf            # ASM 시각화 (4종)
│
├── methylation/
│   ├── {sample}.bw                       # CpG BigWig (IGV용)
│   └── methylation_summary.txt
│
├── cohort/                               # 다중 군 코호트 분석 결과
│   ├── variant_burden_stats.csv          # Small Variant burden 통계
│   ├── variant_burden_plots.pdf          # Burden 시각화 (4종)
│   ├── sv_burden_stats.csv               # SV burden 통계
│   ├── sv_burden_plots.pdf               # SV burden 시각화 (4종)
│   ├── trgt_group_compare.csv            # TRGT 반복서열 군 간 비교
│   ├── trgt_group_compare_plots.pdf      # TRGT 비교 시각화 (4종)
│   ├── dmr_{g1}_vs_{g2}.csv             # 군 쌍별 DMR (자동 생성)
│   ├── dmr_{g1}_vs_{g2}_plots.pdf
│   └── methylation_{group}.txt           # 군별 메틸화 파일 목록
│
├── logs/                                 # 각 규칙별 로그
│   ├── slivar/, vep/, svpack/, dmr/
│   ├── sv_consequence/, trgt/, asm/
│   └── cohort/
│
└── analysis_summary_report.html          # 전체 요약 리포트 (HTML)
```

---

### 주요 결과 파일 설명

#### Small Variant 결과

| 파일 | 설명 |
|------|------|
| `{sample}.filtered.tsv` | slivar 필터 통과 변이 (AF, GQ, DP 기준) |
| `{sample}.annotated.vcf.gz` | VEP 주석 포함 VCF (CADD, REVEL, IMPACT) |
| `{sample}.canonical_impacts.tsv` | canonical 전사체의 SYMBOL·IMPACT·Consequence |

#### Structural Variant 결과

| 파일 | 설명 |
|------|------|
| `{sample}.filtered_sv.vcf.gz` | svpack 필터 통과 SV |
| `{sample}.annotated_sv.tsv` | BCSQ 컬럼 포함: `sv:cds\|GENE\|...` 형태 |

#### TRGT 반복서열 결과

| 파일 | 설명 |
|------|------|
| `{sample}.trgt_outliers.tsv` | 기대 범위(ALLR) 초과 expanded repeat |
| `{sample}.trgt_summary.tsv` | 전체 PASS repeat (allele 길이, 상태 포함) |

컬럼: `trid`, `motif`, `chrom`, `pos`, `allele1_len`, `allele2_len`, `expected_range`, `spanning_depth`, `status`

#### ASM 결과

| 파일 | 설명 |
|------|------|
| `{sample}.asm_results.csv` | Hap1 vs Hap2 차이 메틸화 영역 |
| `{sample}.asm_plots.pdf` | 염색체 분포, delta-β histogram, Top 20 ASM, 크기 vs 유의성 |

#### 코호트 통계 결과 (Phase 5)

**variant_burden_stats.csv / sv_burden_stats.csv**

| 컬럼 | 설명 |
|------|------|
| `gene` | 유전자 이름 |
| `n_{group}` | 그룹별 변이 있는 샘플 수 |
| `pvalue` | Fisher's exact 또는 χ² p값 |
| `fdr` | BH 보정 FDR |
| `odds_ratio` | 2군일 때 Odds Ratio |

**trgt_group_compare.csv**

| 컬럼 | 설명 |
|------|------|
| `trid` | Tandem Repeat ID |
| `motif` | 반복 단위 서열 |
| `pvalue` | Kruskal-Wallis (3군+) 또는 Wilcoxon (2군) p값 |
| `fdr` | BH 보정 FDR |
| `median_len_{group}` | 그룹별 median allele 길이 |
| `pw_p_{g1}_vs_{g2}` | pairwise Wilcoxon p값 (유의 locus) |
| `pw_fdr_{g1}_vs_{g2}` | pairwise BH FDR |

---

## 🔬 다중 군 분석 (Phase 4/5)

### 설계 원칙

- **하위 호환성**: `samples.groups`가 없으면 기존 `control`/`experimental` 2군으로 동작
- **자동 군 쌍 생성**: k개 그룹이 있으면 k(k-1)/2개 비교 쌍 자동 생성
- **선택적 실행**: TRGT 파일이 없으면 TRGT 비교 자동 스킵; 메틸화 파일이 없으면 DMR 비교 스킵

### 통계 방법

| 분석 | 2군 | 3군 이상 |
|------|-----|---------|
| Variant/SV burden | Fisher's exact | χ² (simulate.p.value=TRUE if n<5) |
| TRGT allele 길이 | Wilcoxon rank-sum | Kruskal-Wallis + pairwise Wilcoxon |
| 다중 비교 보정 | BH (Benjamini-Hochberg) FDR | 동일 |
| DMR | DSS (DSS::callDMR) | 각 군 쌍별 독립 실행 |

### 시각화 출력

**variant_burden_plots.pdf** (4종):
1. Top 30 유전자 `-log10(FDR)` 막대 그래프
2. 군별 총 변이 부담 boxplot
3. 유의 유전자 bubble plot (군 × 유전자, 크기=샘플 수)
4. Heatmap (유전자 × 샘플, pheatmap)

**sv_burden_plots.pdf** (4종):
1. 군별 SVTYPE 분포 barplot
2. 군별 SV 크기 분포 violin plot
3. Top 20 유전자 `-log10(FDR)` 막대 그래프
4. 유의 유전자 bubble plot

**trgt_group_compare_plots.pdf** (4종):
1. Manhattan plot (`-log10(p)` vs 게놈 위치)
2. Top 20 유의 repeat violin/boxplot (군별 allele 길이)
3. 군별 expanded repeat 비율 비교
4. Effect size vs `-log10(p)` scatter plot

---

## 🔍 문제 해결

### 1. filelist.csv 파일 없음

```
WorkflowError: filelist.csv 파일이 없고 paths.wdl_out_dir도 미설정
```

**해결**:
```bash
python bin/prepare_inputs.py --config-file config/config.yaml
```

### 2. TRGT/ASM 분석이 실행되지 않음

TRGT 또는 ASM 분석은 filelist.csv의 해당 컬럼이 비어있으면 자동 스킵됩니다. `prepare_inputs.py`가 파일을 찾지 못한 경우 WDL 출력 디렉토리 구조를 확인하세요.

### 3. cohort 분석이 실행되지 않음

`parameters.cohort` 섹션이 config에 없으면 기본값을 사용합니다. `samples.groups`를 정의하지 않으면 기존 `control`/`experimental` 2군으로 코호트 분석이 실행됩니다.

### 4. VEP canonical 추출 오류

```bash
# bcftools +split-vep 플러그인 확인
bcftools +split-vep --help

# 없으면 Python fallback 사용 (자동)
python scripts/extract_vep_fields.py --input sample.annotated.vcf.gz --output /tmp/test.tsv
```

### 5. R 패키지 오류 (jsonlite/pheatmap)

```bash
conda env update -f envs/dmr.yaml --prune
```

### 6. VEP 캐시 오류

```bash
cd /path/to/vep_cache
wget http://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz
tar -xzf homo_sapiens_vep_110_GRCh38.tar.gz
```

### 7. 메모리 부족

```yaml
resources:
  vep:
    threads: 4
    mem_mb: 8000
  dmr:
    threads: 2
    mem_mb: 8000
  cohort_burden:
    threads: 2
    mem_mb: 4000
```

### 8. DMR이 발견되지 않음

```yaml
parameters:
  dmr:
    p_value_cutoff: 0.1
    min_methylation_diff: 0.05
    min_cpg_sites: 2
```

---

## 📚 참고 자료

- [PacBio HiFi Human WGS WDL](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL)
- [Snakemake 문서](https://snakemake.readthedocs.io/)
- [VEP 문서](https://www.ensembl.org/info/docs/tools/vep/index.html)
- [DSS Bioconductor](https://bioconductor.org/packages/release/bioc/html/DSS.html)
- [TRGT — PacBio](https://github.com/PacificBiosciences/trgt)
- [svpack — PacBio](https://github.com/PacificBiosciences/svpack)
- [slivar](https://github.com/brentp/slivar)

---

## 📝 변경 이력

### v3.0.0 (2026-05-20)
- **Phase 4/5 추가**: 다중 군 코호트 통계 분석
  - Small Variant burden test (`variant_burden.R`)
  - SV burden test (`sv_burden.R`)
  - TRGT 반복서열 군 간 비교 (`trgt_group_compare.R`)
  - DMR pairwise (모든 군 쌍 자동 비교, `dmr_pairwise` rule)
  - `samples.groups` config 지원 (하위 호환 유지)
  - VEP canonical IMPACT 추출 (`extract_vep_fields.py`)
- `envs/dmr.yaml`에 `r-pheatmap`, `r-jsonlite`, `r-dplyr`, `r-tidyr` 추가

### v2.0.0 (2026-04)
- **Phase 3 추가**: 신규 분석 모듈
  - SV 유전자 영향 주석 (`svpack consequence` + Ensembl GFF3)
  - TRGT expanded repeat 탐지 (`trgt_outlier.py`)
  - ASM 분석 (`asm_analysis.R`)
  - 분석 요약 리포트 (`analysis_summary_report.html`)
- `bin/prepare_inputs.py`: TRGT/hap1/hap2 파일 자동 탐색
- `dmr_analysis.R` 버그 수정 (조건 없는 중복 코드 제거)

### v1.0.0 (2026-01-29)
- 초기 릴리스
- Small/Structural Variant 필터링 및 VEP 주석
- DMR 분석, 메틸화 BigWig 정리
