# PacBio HiFi WGS 3차 분석 파이프라인

PacBio WDL 파이프라인(singleton.wdl) 실행 결과를 입력받아 사용자 정의 필터링, 심화 주석, 그리고 차이 메틸화 영역(DMR) 분석을 수행하는 Snakemake 기반의 자동화 파이프라인입니다.

## 📋 목차

- [개요](#개요)
- [프로젝트 구조](#프로젝트-구조)
- [시스템 요구사항](#시스템-요구사항)
- [설치](#설치)
- [설정](#설정)
- [사용법](#사용법)
- [출력 결과](#출력-결과)
- [문제 해결](#문제-해결)

---

## 🎯 개요

본 파이프라인은 다음과 같은 분석을 수행합니다:

### Phase 1: 변이 필터링 및 주석 보강
- **Small Variant Re-filtering**: `slivar`를 사용한 질환 연관 변이 우선순위화
- **Functional Annotation**: `VEP`를 통한 CADD, REVEL 등 기능 예측 점수 추가
- **Structural Variant Re-filtering**: `svpack`을 사용한 인구 집단 빈도 기반 필터링

### Phase 2: 후성유전학 심화 분석
- **Differential Methylation Analysis**: `DSS`를 사용한 차이 메틸화 영역(DMR) 분석
- **Methylation Visualization**: IGV 호환 BigWig 파일 정리

---

## 📁 프로젝트 구조

```
wgs-tertiary-pipeline/
│
├── Snakefile                    # 메인 워크플로우
├── README.md                    # 이 문서
├── environment.yaml             # Conda 환경
│
├── config/                      # 설정 파일
│   ├── config.yaml              # 실제 설정 (사용자 수정)
│   ├── config.example.yaml      # 설정 템플릿
│   └── cluster_config.yaml      # SLURM 리소스
│
├── bin/                         # 실행 스크립트
│   ├── run_pipeline.sh          # 메인 실행
│   ├── setup_check.sh           # 환경 검증
│   └── run_slurm.sh             # 클러스터 실행
│
├── scripts/                     # 분석 스크립트
│   └── dmr_analysis.R           # DMR 분석
│
└── docs/                        # 추가 문서
    ├── QUICKSTART.md            # 빠른 시작 가이드
    └── SUMMARY.md               # 프로젝트 요약
```

자세한 구조 설명은 [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)를 참고하세요.

---

## 💻 시스템 요구사항

### 필수 소프트웨어

- **Snakemake** (≥ 7.0)
- **Docker** 또는 **Singularity** (컨테이너 실행용)
- **Python** (≥ 3.8)

### 컨테이너 이미지

파이프라인은 다음 컨테이너를 자동으로 사용합니다:

- `quay.io/pacbio/slivar:latest` - Small variant 필터링
- `quay.io/pacbio/svpack:latest` - Structural variant 필터링
- `ensemblorg/ensembl-vep:latest` - VEP 주석
- `bioconductor/bioconductor_docker:latest` - R/DSS DMR 분석

### 하드웨어 권장사항

- **CPU**: 16 코어 이상
- **메모리**: 32GB 이상
- **디스크**: 500GB 이상 여유 공간

---

## 🔧 설치

### 1. Snakemake 설치

```bash
# Conda를 통한 설치 (권장)
conda install -c bioconda snakemake

# 또는 pip를 통한 설치
pip install snakemake
```

### 2. 컨테이너 런타임 설치

```bash
# Docker 설치 (Ubuntu 기준)
sudo apt-get update
sudo apt-get install docker.io

# 사용자를 docker 그룹에 추가
sudo usermod -aG docker $USER

# 또는 Singularity 설치
sudo apt-get install singularity-container
```

### 3. 파이프라인 클론

```bash
git clone <repository-url>
cd wgs-tertiary-pipeline
```

---

## ⚙️ 설정

### config/config.yaml 편집

분석 전 `config/config.yaml` 파일을 환경에 맞게 수정하세요:

```yaml
# 디렉토리 설정
paths:
  wdl_out_dir: "/path/to/wdl/out"        # WDL 파이프라인 결과 폴더
  ref_genome: "/path/to/GRCh38.fa"        # 참조 유전체
  vep_cache: "/path/to/vep_cache"         # VEP 캐시 디렉토리
  output_dir: "./3rd_analysis_results"    # 출력 디렉토리

# 샘플 및 그룹 설정
samples:
  control:
    - sample_A      # 대조군 샘플 ID
    - sample_B
  experimental:
    - sample_C      # 실험군 샘플 ID
    - sample_D

# 분석 파라미터 (필요시 조정)
parameters:
  slivar:
    max_af: 0.01          # 최대 인구 집단 빈도
    min_gq: 20            # 최소 Genotype Quality
    min_dp: 10            # 최소 Read Depth
    hpo_terms: "HP:0001250,HP:0000707"  # HPO 용어
  
  svpack:
    min_sv_size: 50       # 최소 SV 크기 (bp)
    filter_pass_only: true
  
  dmr:
    p_value_cutoff: 0.05
    min_methylation_diff: 0.1
```

### 주요 설정 항목 설명

#### 1. **paths.wdl_out_dir**
- WDL 파이프라인의 `out/` 디렉토리 경로
- 다음 파일들이 포함되어 있어야 합니다:
  - `phased_small_variant_{sample}.vcf.gz`
  - `phased_sv_{sample}.vcf.gz`
  - `cpg_combined_{sample}.bed`
  - `cpg_combined_{sample}.bw`

#### 2. **samples**
- `control`: 대조군 샘플 ID 목록
- `experimental`: 실험군 샘플 ID 목록
- 샘플 ID는 WDL 출력 파일명과 일치해야 합니다

#### 3. **parameters.slivar.hpo_terms**
- Human Phenotype Ontology 용어 (콤마로 구분)
- [HPO Browser](https://hpo.jax.org/app/)에서 검색 가능
- 예: `HP:0001250` (Seizures), `HP:0000707` (Nervous system abnormality)

---

## 🚀 사용법

### 0. WDL 결과 파일 자동 탐색 (권장)

WDL 파이프라인의 복잡한 디렉토리 구조에서 필요한 파일을 자동으로 찾습니다:

```bash
# config.yaml 설정 후 실행
python bin/prepare_inputs.py

# 실행 결과: filelist.csv 생성
# 이 파일은 Snakemake가 자동으로 사용합니다
```

자세한 내용은 [WDL 결과 파일 자동 탐색 가이드](docs/PREPARE_INPUTS_GUIDE.md)를 참조하세요.

### 1. Dry-run (테스트 실행)

실제 실행 전 파이프라인 구조를 확인:

```bash
snakemake --dry-run --cores 1
```

### 2. 전체 파이프라인 실행

```bash
# Docker 사용 시
snakemake --use-singularity --cores 16

# Singularity 사용 시 (HPC 환경)
snakemake --use-singularity --cores 16
```

### 3. 특정 규칙만 실행

```bash
# Small variant 필터링만 실행
snakemake --use-singularity --cores 8 \
  3rd_analysis_results/small_variants/sample_A.filtered.tsv

# DMR 분석만 실행
snakemake --use-singularity --cores 4 \
  3rd_analysis_results/dmr_analysis/dmr_results.csv
```

### 4. 클러스터 환경에서 실행

SLURM 클러스터 예시:

```bash
snakemake \
  --use-singularity \
  --cluster "sbatch -p normal -n {threads} --mem={resources.mem_mb}" \
  --jobs 10 \
  --cores 64
```

### 5. 실행 상태 확인

```bash
# DAG 시각화
snakemake --dag | dot -Tpdf > dag.pdf

# 규칙 그래프
snakemake --rulegraph | dot -Tpdf > rulegraph.pdf
```

---

## 📊 출력 결과

파이프라인 완료 후 다음과 같은 구조로 결과가 생성됩니다:

```
3rd_analysis_results/
├── small_variants/
│   ├── sample_A.filtered.vcf.gz       # 필터링된 Small Variant
│   ├── sample_A.filtered.tsv          # TSV 형식 (Excel 호환)
│   ├── sample_A.annotated.vcf.gz      # VEP 주석 추가
│   ├── sample_A.annotated.tsv         # TSV 형식
│   └── sample_A.vep_stats.html        # VEP 통계 리포트
│
├── structural_variants/
│   ├── sample_A.filtered_sv.vcf.gz    # 필터링된 SV
│   └── sample_A.filtered_sv.vcf.gz.tbi
│
├── dmr_analysis/
│   ├── dmr_results.csv                # DMR 결과 테이블
│   ├── dmr_plots.pdf                  # DMR 시각화 (5개 그래프)
│   ├── control_samples.txt
│   └── experimental_samples.txt
│
├── methylation/
│   ├── sample_A.bw                    # 메틸화 BigWig 파일 (IGV용)
│   ├── sample_B.bw
│   └── methylation_summary.txt
│
├── logs/
│   ├── slivar/
│   ├── vep/
│   ├── svpack/
│   └── dmr/
│
└── analysis_summary_report.html       # 전체 요약 리포트
```

### 주요 결과 파일 설명

#### 1. **Small Variant 결과**

- **`.filtered.tsv`**: 엑셀에서 열 수 있는 필터링된 변이 목록
  - 컬럼: CHROM, POS, REF, ALT, AF, GQ, DP 등
  
- **`.annotated.vcf.gz`**: VEP 주석이 추가된 VCF
  - CADD, REVEL 점수 포함
  - CSQ 필드에 기능 예측 정보

#### 2. **DMR 분석 결과**

- **`dmr_results.csv`**: DMR 목록
  - 컬럼: chr, start, end, length, nCG, meanMethy1, meanMethy2, diff_methy, pvalue, fdr
  
- **`dmr_plots.pdf`**: 5개의 시각화
  1. DMR 분포 (염색체별)
  2. 메틸화 차이 분포
  3. Volcano plot
  4. Top 20 DMR
  5. DMR 크기 vs 유의성

#### 3. **메틸화 BigWig 파일**

- IGV, UCSC Genome Browser 등에서 시각화 가능
- 각 샘플의 CpG 메틸화 수준 표시

---

## 🔍 문제 해결

### 1. 파일을 찾을 수 없음 (FileNotFoundError)

**증상**: `phased_small_variant_*.vcf.gz` 파일을 찾을 수 없다는 에러

**해결**:
- `config/config.yaml`의 `wdl_out_dir` 경로가 올바른지 확인
- WDL 파이프라인 출력 파일명이 예상과 일치하는지 확인
- 샘플 ID가 정확한지 확인

```bash
# WDL 출력 디렉토리 확인
ls -l /path/to/wdl/out/phased_*.vcf.gz
```

### 2. 컨테이너 이미지 다운로드 실패

**증상**: Docker/Singularity 이미지를 가져올 수 없음

**해결**:
```bash
# Docker 이미지 수동 다운로드
docker pull quay.io/pacbio/slivar:latest
docker pull quay.io/pacbio/svpack:latest
docker pull ensemblorg/ensembl-vep:latest
docker pull bioconductor/bioconductor_docker:latest
```

### 3. VEP 캐시 오류

**증상**: VEP 실행 시 캐시를 찾을 수 없음

**해결**:
```bash
# VEP 캐시 다운로드 (GRCh38)
cd /path/to/vep_cache
wget http://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz
tar -xzf homo_sapiens_vep_110_GRCh38.tar.gz
```

### 4. R 패키지 오류 (DSS)

**증상**: DMR 분석 시 R 패키지를 찾을 수 없음

**해결**:
```bash
# 컨테이너 내에서 R 패키지 설치 (필요시)
docker run -it bioconductor/bioconductor_docker:latest R

# R 콘솔에서
BiocManager::install(c("DSS", "bsseq"))
```

### 5. 메모리 부족

**증상**: `Killed` 또는 `Out of Memory` 에러

**해결**:
- `config/config.yaml`의 리소스 설정 조정:
```yaml
resources:
  vep:
    threads: 4  # 스레드 수 감소
    mem_mb: 8000  # 메모리 감소
```

### 6. DMR이 발견되지 않음

**증상**: "유의한 DMR이 발견되지 않았습니다" 메시지

**해결**:
- 파라미터 완화:
```yaml
parameters:
  dmr:
    p_value_cutoff: 0.1  # 완화
    min_methylation_diff: 0.05  # 완화
    min_cpg_sites: 2  # 감소
```

---

## 📚 추가 정보

### 파이프라인 구조

```
wgs-tertiary-pipeline/
├── Snakefile               # 메인 파이프라인
├── config/
│   ├── config.yaml         # 설정 파일
│   ├── config.example.yaml # 설정 예제
│   └── cluster_config.yaml # SLURM 리소스 설정
├── scripts/
│   └── dmr_analysis.R      # DMR 분석 R 스크립트
├── bin/
│   ├── run_pipeline.sh     # 파이프라인 실행 스크립트
│   ├── setup_check.sh      # 환경 검증 스크립트
│   └── run_slurm.sh        # SLURM 실행 스크립트
├── docs/
│   ├── QUICKSTART.md       # 빠른 시작 가이드
│   └── SUMMARY.md          # 프로젝트 요약
├── README.md               # 이 문서
└── environment.yaml        # Conda 환경 (선택사항)
```

### 참고 자료

- [PacBio WDL Workflows](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL)
- [Snakemake 문서](https://snakemake.readthedocs.io/)
- [slivar 문서](https://github.com/brentp/slivar)
- [VEP 문서](https://www.ensembl.org/info/docs/tools/vep/index.html)
- [DSS Bioconductor](https://bioconductor.org/packages/release/bioc/html/DSS.html)

### 라이센스

이 파이프라인은 MIT 라이센스 하에 배포됩니다.

### 문의

버그 리포트나 기능 제안은 이슈 트래커를 통해 제출해주세요.

---

## 📝 변경 이력

### v1.0.0 (2026-01-29)
- 초기 릴리스
- Small/Structural Variant 필터링 및 주석
- DMR 분석 기능 구현
- 자동화된 리포트 생성
