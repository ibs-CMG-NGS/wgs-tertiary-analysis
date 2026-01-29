# 빠른 시작 가이드 (Quick Start Guide)

이 가이드는 PacBio HiFi WGS 3차 분석 파이프라인을 빠르게 시작하기 위한 단계별 지침입니다.

## ⚡ 5분 안에 시작하기

### 1단계: 환경 설정 (처음 1회만)

```bash
# Conda 환경 생성 및 활성화
conda env create -f environment.yaml
conda activate pacbio-tertiary

# 또는 기존 환경에 Snakemake만 설치
conda install -c bioconda snakemake
```

### 2단계: 설정 파일 준비

```bash
# 예제 설정 파일 복사
cp config/config.example.yaml config/config.yaml

# 설정 파일 편집
nano config/config.yaml  # 또는 vim, emacs 등
```

**필수 수정 항목:**
- `paths.wdl_out_dir`: WDL 파이프라인 출력 디렉토리 경로
- `paths.ref_genome`: 참조 유전체 경로
- `paths.vep_cache`: VEP 캐시 경로
- `samples.control`: 대조군 샘플 ID 목록
- `samples.experimental`: 실험군 샘플 ID 목록

### 3단계: 설정 확인

```bash
# 설정 검증 스크립트 실행
bash bin/setup_check.sh
```

### 4단계: Dry-run 테스트

```bash
# 실제 실행 없이 파이프라인 구조 확인
./bin/run_pipeline.sh --dry-run
```

### 5단계: 파이프라인 실행

```bash
# 전체 파이프라인 실행 (16 코어)
./bin/run_pipeline.sh

# 또는 직접 Snakemake 실행
snakemake --use-singularity --cores 16
```

---

## 📝 상세 설정 가이드

### WDL 출력 디렉토리 구조 확인

WDL 파이프라인 실행 후 다음과 같은 파일들이 있어야 합니다:

```
/path/to/wdl/out/
├── phased_small_variant_sample_A.vcf.gz
├── phased_small_variant_sample_B.vcf.gz
├── phased_sv_sample_A.vcf.gz
├── phased_sv_sample_B.vcf.gz
├── cpg_combined_sample_A.bed
├── cpg_combined_sample_B.bed
├── cpg_combined_sample_A.bw
└── cpg_combined_sample_B.bw
```

### 샘플 ID 매칭

`config/config.yaml`의 샘플 ID는 **파일명의 접미사 부분**과 일치해야 합니다:

```yaml
samples:
  control:
    - sample_A    # ← phased_small_variant_sample_A.vcf.gz
    - sample_B    # ← phased_small_variant_sample_B.vcf.gz
```

### VEP 캐시 다운로드

```bash
# 캐시 디렉토리 생성
mkdir -p /path/to/vep_cache
cd /path/to/vep_cache

# GRCh38 캐시 다운로드 (약 20GB)
wget http://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz

# 압축 해제
tar -xzf homo_sapiens_vep_110_GRCh38.tar.gz
```

---

## 🎯 사용 사례별 실행 방법

### Case 1: 로컬 워크스테이션 (16 코어, 32GB RAM)

```bash
./bin/run_pipeline.sh --cores 16
```

### Case 2: 로컬 워크스테이션 (리소스 제한)

`config/config.yaml`에서 리소스 조정:

```yaml
resources:
  vep:
    threads: 4      # 8 → 4로 감소
    mem_mb: 8000    # 16000 → 8000으로 감소
```

그 후 실행:

```bash
./bin/run_pipeline.sh --cores 8
```

### Case 3: SLURM 클러스터

```bash
# SLURM 작업 제출
sbatch bin/run_slurm.sh

# 작업 상태 확인
squeue -u $USER

# 로그 확인
tail -f logs/slurm_*.out
```

### Case 4: 특정 분석만 실행

```bash
# Small variant 필터링만
./bin/run_pipeline.sh --target filter_small_variants

# DMR 분석만
./bin/run_pipeline.sh --target run_dmr_analysis

# 특정 샘플의 VEP 주석만
./bin/run_pipeline.sh --target 3rd_analysis_results/small_variants/sample_A.annotated.vcf.gz
```

---

## 🔧 파라미터 튜닝 가이드

### 희귀 질환 변이 분석 (엄격한 필터)

```yaml
parameters:
  slivar:
    max_af: 0.001      # 0.1% 이하
    min_gq: 30         # 높은 품질
    min_dp: 20         # 높은 coverage
  
  dmr:
    p_value_cutoff: 0.01    # 엄격한 유의성
    min_methylation_diff: 0.2  # 큰 차이만
```

### 일반 변이 분석 (완화된 필터)

```yaml
parameters:
  slivar:
    max_af: 0.05       # 5% 이하
    min_gq: 15
    min_dp: 8
  
  dmr:
    p_value_cutoff: 0.1
    min_methylation_diff: 0.05
```

### HPO 용어 설정 (질환별)

```yaml
# 신경발달장애
hpo_terms: "HP:0001263,HP:0000729,HP:0001250"

# 심혈관질환
hpo_terms: "HP:0001627,HP:0001638,HP:0001645"

# 암
hpo_terms: "HP:0002664,HP:0100013"
```

→ [HPO Browser](https://hpo.jax.org/app/)에서 검색

---

## 📊 결과 확인

### 디렉토리 구조

```
3rd_analysis_results/
├── small_variants/          # Small variant 결과
│   ├── *.filtered.tsv       ← Excel에서 열기
│   ├── *.annotated.tsv      ← Excel에서 열기
│   └── *.vep_stats.html     ← 브라우저에서 열기
│
├── dmr_analysis/
│   ├── dmr_results.csv      ← Excel에서 열기
│   └── dmr_plots.pdf        ← PDF 뷰어에서 열기
│
└── analysis_summary_report.html  ← 브라우저에서 열기
```

### 주요 결과 파일

1. **`small_variants/*.filtered.tsv`**
   - 필터링된 small variant 목록
   - 엑셀에서 정렬/필터링 가능

2. **`dmr_results.csv`**
   - 차이 메틸화 영역 목록
   - FDR로 정렬되어 있음
   - 상위 결과가 가장 유의함

3. **`dmr_plots.pdf`**
   - 5개의 시각화:
     - DMR 염색체별 분포
     - 메틸화 차이 히스토그램
     - Volcano plot
     - Top 20 DMR
     - DMR 크기 vs 유의성

---

## ❓ 자주 묻는 질문 (FAQ)

### Q1: "FileNotFoundError: phased_small_variant_*.vcf.gz"

**A:** 샘플 ID가 파일명과 일치하지 않습니다.

```bash
# 실제 파일명 확인
ls /path/to/wdl/out/phased_*.vcf.gz

# config/config.yaml의 샘플 ID를 파일명과 일치시키기
```

### Q2: VEP 실행이 너무 느림

**A:** 스레드 수를 늘리거나, 오프라인 캐시를 사용하세요.

```yaml
resources:
  vep:
    threads: 16  # 증가
```

### Q3: DMR이 하나도 발견되지 않음

**A:** 파라미터를 완화하거나 샘플 수를 늘리세요.

```yaml
parameters:
  dmr:
    p_value_cutoff: 0.1      # 완화
    min_methylation_diff: 0.05  # 완화
```

### Q4: 메모리 부족 에러

**A:** 리소스를 줄이거나, 큰 메모리의 노드를 사용하세요.

```yaml
resources:
  vep:
    mem_mb: 8000  # 감소
```

### Q5: 중간에 실패했는데 재시작하려면?

**A:** Snakemake는 자동으로 완료된 단계를 건너뜁니다.

```bash
# 그냥 다시 실행
./run_pipeline.sh

# 강제 재실행이 필요하면
snakemake --use-singularity --cores 16 --forceall
```

---

## 📚 다음 단계

- 자세한 내용은 [README.md](README.md) 참조
- 문제 발생 시 [Issues](../../issues) 등록
- 파이프라인 커스터마이징은 `Snakefile` 수정

---

**🎉 분석 성공을 기원합니다!**
