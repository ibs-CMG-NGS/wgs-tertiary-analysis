# PacBio HiFi WGS 3차 분석 파이프라인 - 프로젝트 요약

## 🎯 프로젝트 개요

PacBio WDL 파이프라인의 출력 결과를 입력으로 받아 다음 분석을 자동화하는 Snakemake 파이프라인입니다:

1. **Small Variant 재필터링** (slivar)
2. **기능 주석 추가** (VEP)
3. **Structural Variant 재필터링** (svpack)
4. **차이 메틸화 영역 분석** (DSS/R)
5. **통합 리포트 생성**

---

## 📁 프로젝트 구조

```
wgs-tertiary-pipeline/
│
├── Snakefile                   # 메인 파이프라인 (Snakemake 워크플로우)
├── environment.yaml            # Conda 환경 정의
│
├── config/
│   ├── config.yaml             # 실제 설정 파일 (사용자가 편집)
│   ├── config.example.yaml     # 설정 예제 (템플릿)
│   └── cluster_config.yaml     # SLURM 리소스 설정
│
├── scripts/
│   └── dmr_analysis.R          # DMR 분석 R 스크립트 (DSS)
│
├── bin/
│   ├── run_pipeline.sh         # 파이프라인 실행 스크립트
│   ├── setup_check.sh          # 환경 검증 스크립트
│   └── run_slurm.sh            # SLURM 클러스터 실행 스크립트
│
├── docs/
│   ├── QUICKSTART.md           # 빠른 시작 가이드
│   └── SUMMARY.md              # 이 문서 (프로젝트 요약)
│
├── README.md                   # 상세 사용 설명서
└── .gitignore                  # Git 제외 파일 목록
```

---

## 🚀 핵심 기능

### Phase 1: 변이 필터링 및 주석

#### 1.1 Small Variant 필터링 (slivar)
- **입력**: `phased_small_variant_{sample}.vcf.gz`
- **필터링 기준**:
  - 인구 집단 빈도 (AF < 0.01)
  - Genotype Quality (GQ ≥ 20)
  - Read Depth (DP ≥ 10)
  - HPO 기반 질환 연관성
- **출력**: 
  - `{sample}.filtered.vcf.gz`
  - `{sample}.filtered.tsv` (Excel 호환)

#### 1.2 VEP 주석
- **입력**: 필터링된 VCF
- **추가 정보**:
  - CADD 점수
  - REVEL 점수
  - 기능 예측 (missense, nonsense 등)
  - 전사체 정보
- **출력**:
  - `{sample}.annotated.vcf.gz`
  - `{sample}.annotated.tsv`
  - `{sample}.vep_stats.html`

#### 1.3 Structural Variant 필터링 (svpack)
- **입력**: `phased_sv_{sample}.vcf.gz`
- **필터링**: 
  - 최소 크기 (≥ 50bp)
  - PASS 필터만
  - 인구 집단 빈도
- **출력**: `{sample}.filtered_sv.vcf.gz`

### Phase 2: 후성유전학 분석

#### 2.1 DMR 분석 (DSS)
- **입력**: 
  - Control 그룹: `cpg_combined_{control}.bed`
  - Experimental 그룹: `cpg_combined_{experimental}.bed`
- **분석**:
  - 통계적 유의성 검정 (p < 0.05)
  - 메틸화 차이 계산 (≥ 10%)
  - FDR 보정
- **출력**:
  - `dmr_results.csv` (DMR 목록)
  - `dmr_plots.pdf` (5개 시각화)

#### 2.2 메틸화 시각화
- **입력**: `cpg_combined_{sample}.bw`
- **처리**: BigWig 파일 정리 및 IGV 호환 준비
- **출력**: `methylation/{sample}.bw`

---

## 📊 분석 결과

### 출력 디렉토리 구조

```
3rd_analysis_results/
│
├── small_variants/
│   ├── sample_A.filtered.vcf.gz       # 필터링된 변이
│   ├── sample_A.filtered.tsv          # TSV 형식
│   ├── sample_A.annotated.vcf.gz      # VEP 주석 추가
│   ├── sample_A.annotated.tsv         # TSV 형식
│   └── sample_A.vep_stats.html        # VEP 통계
│
├── structural_variants/
│   └── sample_A.filtered_sv.vcf.gz    # 필터링된 SV
│
├── dmr_analysis/
│   ├── dmr_results.csv                # DMR 테이블
│   ├── dmr_plots.pdf                  # 시각화 (5개)
│   ├── control_samples.txt
│   └── experimental_samples.txt
│
├── methylation/
│   ├── sample_A.bw                    # BigWig (IGV용)
│   └── methylation_summary.txt
│
├── logs/                              # 모든 로그
│   ├── slivar/
│   ├── vep/
│   ├── svpack/
│   └── dmr/
│
└── analysis_summary_report.html       # 전체 요약
```

---

## ⚙️ 설정 파라미터

### config.yaml 주요 섹션

```yaml
# 1. 경로 설정
paths:
  wdl_out_dir: "/path/to/wdl/out"      # WDL 출력
  ref_genome: "/path/to/GRCh38.fa"     # 참조 유전체
  vep_cache: "/path/to/vep_cache"      # VEP 캐시
  output_dir: "./3rd_analysis_results" # 출력

# 2. 샘플 그룹
samples:
  control: [sample_A, sample_B]
  experimental: [sample_C, sample_D]

# 3. 분석 파라미터
parameters:
  slivar:
    max_af: 0.01           # 인구 집단 빈도
    min_gq: 20             # 최소 GQ
    min_dp: 10             # 최소 DP
    hpo_terms: "HP:..."    # HPO 용어
  
  svpack:
    min_sv_size: 50        # 최소 SV 크기
  
  dmr:
    p_value_cutoff: 0.05   # p-value
    min_methylation_diff: 0.1  # 메틸화 차이

# 4. 컨테이너
containers:
  pacbio_slivar: "quay.io/pacbio/slivar:latest"
  pacbio_svpack: "quay.io/pacbio/svpack:latest"
  vep: "ensemblorg/ensembl-vep:latest"
  r_dss: "bioconductor/bioconductor_docker:latest"

# 5. 리소스
resources:
  vep:
    threads: 8
    mem_mb: 16000
```

---

## 🔧 사용 방법

### 1. 빠른 시작 (로컬)

```bash
# 1. 설정 파일 준비
cp config/config.example.yaml config/config.yaml
nano config/config.yaml  # 경로와 샘플 수정

# 2. 환경 확인
bash bin/setup_check.sh

# 3. Dry-run
./bin/run_pipeline.sh --dry-run

# 4. 실행
./bin/run_pipeline.sh --cores 16
```

### 2. SLURM 클러스터

```bash
# 클러스터 작업 제출
sbatch bin/run_slurm.sh

# 상태 확인
squeue -u $USER
```

### 3. 부분 실행

```bash
# Small variant만
./bin/run_pipeline.sh --target filter_small_variants

# DMR만
./bin/run_pipeline.sh --target run_dmr_analysis
./run_pipeline.sh --target filter_small_variants

# DMR만
./run_pipeline.sh --target run_dmr_analysis
```

---

## 🔬 Snakemake 규칙 (Rules)

### 1. `filter_small_variants`
- **입력**: `phased_small_variant_{sample}.vcf.gz`
- **도구**: slivar
- **출력**: 필터링된 VCF 및 TSV

### 2. `annotate_vep`
- **입력**: 필터링된 VCF
- **도구**: VEP (ensembl)
- **출력**: 주석된 VCF 및 TSV

### 3. `filter_sv`
- **입력**: `phased_sv_{sample}.vcf.gz`
- **도구**: svpack
- **출력**: 필터링된 SV VCF

### 4. `prepare_methylation_data`
- **입력**: 모든 `cpg_combined_{sample}.bed`
- **출력**: 그룹별 샘플 목록

### 5. `run_dmr_analysis`
- **입력**: 샘플 목록, R 스크립트
- **도구**: R/DSS
- **출력**: DMR 결과 CSV 및 PDF

### 6. `merge_methylation_plots`
- **입력**: 모든 `.bw` 파일
- **출력**: 정리된 BigWig 및 요약

### 7. `generate_summary_report`
- **입력**: 모든 결과 파일
- **출력**: HTML 리포트

---

## 📦 의존성

### 소프트웨어
- Snakemake ≥ 7.0
- Python ≥ 3.8
- Docker 또는 Singularity

### 컨테이너 이미지 (자동 다운로드)
- PacBio slivar
- PacBio svpack
- Ensembl VEP
- Bioconductor/DSS

### R 패키지 (컨테이너 내 포함)
- DSS
- bsseq
- ggplot2
- data.table

---

## 🎓 주요 도구 설명

### slivar
- **목적**: VCF 필터링 및 질환 연관 변이 우선순위화
- **특징**: 
  - 빠른 필터링
  - HPO 기반 Phrank 점수
  - 유연한 표현식
- **문서**: https://github.com/brentp/slivar

### VEP (Variant Effect Predictor)
- **목적**: 변이의 기능적 영향 예측
- **특징**:
  - CADD, REVEL 등 통합
  - 전사체 레벨 주석
  - 플러그인 시스템
- **문서**: https://www.ensembl.org/vep

### svpack
- **목적**: Structural variant 처리
- **특징**:
  - PacBio 최적화
  - 인구 집단 빈도 필터링
- **문서**: PacBio 공식 문서

### DSS
- **목적**: 차이 메틸화 분석
- **특징**:
  - Bisulfite sequencing 특화
  - 통계적 모델링
  - DMR 호출
- **문서**: https://bioconductor.org/packages/DSS

---

## 📈 성능 및 리소스

### 예상 실행 시간 (4 샘플 기준)

| 단계 | 시간 | 코어 | 메모리 |
|------|------|------|--------|
| Small Variant 필터링 | ~1h | 4 | 8GB |
| VEP 주석 | ~4h | 8 | 16GB |
| SV 필터링 | ~30m | 2 | 4GB |
| DMR 분석 | ~2h | 4 | 16GB |
| **전체** | **~8h** | 16 | 32GB |

### 디스크 사용량
- 입력 데이터: ~50GB/샘플
- 출력 결과: ~20GB/샘플
- 임시 파일: ~10GB/샘플
- **권장 여유 공간**: 500GB

---

## 🛠️ 문제 해결

### 일반적인 문제

1. **파일을 찾을 수 없음**
   - config/config.yaml의 경로 확인
   - 샘플 ID와 파일명 일치 확인
   - 샘플 ID와 파일명 일치 확인

2. **메모리 부족**
   - 리소스 설정 조정
   - 배치 크기 감소

3. **컨테이너 오류**
   - Docker/Singularity 권한 확인
   - 이미지 수동 다운로드

4. **DMR 없음**
   - 파라미터 완화
   - 샘플 수 확인
   - 그룹 간 차이 확인

자세한 내용은 README.md의 "문제 해결" 섹션 참조

---

## 📚 문서

- **README.md**: 상세 사용 설명서
- **QUICKSTART.md**: 5분 시작 가이드
- **이 파일 (SUMMARY.md)**: 프로젝트 전체 요약

---

## 🔄 워크플로우 DAG

```
WDL 출력
    ↓
┌───────────────────────────────────────┐
│  filter_small_variants (slivar)       │
└───────────────┬───────────────────────┘
                ↓
        ┌──────────────────┐
        │  annotate_vep    │
        └──────────────────┘
                
┌───────────────────────────────────────┐
│  filter_sv (svpack)                   │
└───────────────────────────────────────┘

┌───────────────────────────────────────┐
│  prepare_methylation_data             │
└───────────────┬───────────────────────┘
                ↓
        ┌──────────────────┐
        │  run_dmr_analysis│
        └──────────────────┘
                
┌───────────────────────────────────────┐
│  merge_methylation_plots              │
└───────────────────────────────────────┘
                ↓
        ┌──────────────────┐
        │  Summary Report  │
        └──────────────────┘
```

---

## ✨ 주요 특징

✅ **완전 자동화**: 설정 후 한 번의 명령으로 전체 분석  
✅ **재현 가능**: 컨테이너 기반으로 환경 일관성 보장  
✅ **확장 가능**: Snakemake 병렬 처리로 수십~수백 샘플 처리  
✅ **클러스터 지원**: SLURM 등 HPC 환경 지원  
✅ **사용자 친화적**: TSV 출력으로 Excel 분석 가능  
✅ **시각화**: 자동 그래프 및 리포트 생성  

---

## 📝 버전 정보

- **버전**: 1.0.0
- **날짜**: 2026-01-29
- **라이센스**: MIT
- **저자**: NGS core, center for memory and glioscience, institute for basic science 

---

## 🆘 지원

- **이슈 리포트**: GitHub Issues
- **문서**: README.md, QUICKSTART.md
- **예제**: config.example.yaml

---

**Happy Analysis! 🧬🔬**
