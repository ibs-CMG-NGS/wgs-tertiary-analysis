# 🎯 완전한 워크플로우 가이드

## 📋 전체 프로세스

```
WDL 파이프라인 실행
        ↓
batch_results/ 생성
        ↓
config.yaml 작성 ← 샘플 ID, batch_results 경로
        ↓
prepare_inputs.py 실행 ← 파일 자동 탐색
        ↓
filelist.csv 생성
        ↓
Snakemake 파이프라인 실행
        ↓
3rd_analysis_results/ 생성
```

---

## 🚀 단계별 실행 가이드

### Step 1: WDL 파이프라인 완료 확인

```bash
# WDL 결과 구조 확인
ls -la /data_4tb/hifi-human-wgs-wdl-custom/batch_results/

# 각 샘플의 out 디렉토리 확인
ls -la /data_4tb/hifi-human-wgs-wdl-custom/batch_results/KTY9537/*/out/
```

예상 출력:
```
batch_results/
├── KTY9537/
│   └── 20260120_101704_humanwgs_singleton/
│       └── out/
│           ├── phased_small_variant_vcf
│           ├── phased_sv_vcf
│           ├── cpg_combined_bed
│           └── cpg_combined_bw
├── LDK6217/
└── 15_JANG_SJ/
```

### Step 2: 3차 분석 파이프라인 설정

```bash
# 프로젝트 디렉토리로 이동
cd wgs-tertiary-pipeline

# 설정 파일 생성
cp config/config.example.yaml config/config.yaml

# 설정 편집
nano config/config.yaml
```

**config.yaml 수정 내용:**

```yaml
paths:
  # ✅ WDL batch_results 경로 (필수)
  batch_results_dir: "/data_4tb/hifi-human-wgs-wdl-custom/batch_results"
  
  # ✅ 참조 유전체
  ref_genome: "/reference/GRCh38/GRCh38.fa"
  
  # ✅ VEP 캐시
  vep_cache: "/data/vep_cache"

samples:
  # ✅ 샘플 ID = batch_results 하위 폴더명
  control:
    - KTY9537
    - LDK6217
  experimental:
    - 15_JANG_SJ
    - 16_BYUN_HY
```

### Step 3: 입력 파일 자동 탐색

```bash
# 파일 자동 탐색 실행
python bin/prepare_inputs.py

# 또는 dry-run으로 먼저 테스트
python bin/prepare_inputs.py --dry-run
```

**성공 시 출력:**

```
================================================================================
WDL 결과 파일 탐색 시작
================================================================================

[CONTROL] 그룹 처리 중...
  └─ KTY9537 탐색 중... ✓ 완료
  └─ LDK6217 탐색 중... ✓ 완료

[EXPERIMENTAL] 그룹 처리 중...
  └─ 15_JANG_SJ 탐색 중... ✓ 완료
  └─ 16_BYUN_HY 탐색 중... ✓ 완료

✓ filelist.csv 생성 완료: 4개 샘플
```

**생성된 filelist.csv 확인:**

```bash
head -n 3 filelist.csv
```

### Step 4: 환경 검증

```bash
# 설정 및 환경 확인
bash bin/setup_check.sh
```

**예상 출력:**

```
✓ snakemake 발견됨
✓ docker 발견됨
✓ config/config.yaml
✓ Snakefile
✓ scripts/dmr_analysis.R
✓ config.yaml 유효성 검사 통과
  - Control 샘플: 2개
  - Experimental 샘플: 2개
✓ WDL 출력 디렉토리: /data_4tb/.../batch_results
✓ phased VCF 파일: 4개 발견

모든 확인 통과!
```

### Step 5: Dry-run 테스트

```bash
# 실제 실행 없이 구조 확인
./bin/run_pipeline.sh --dry-run
```

**확인 사항:**
- 모든 입력 파일 경로가 올바른지
- 필요한 모든 규칙이 실행될지
- 리소스 할당이 적절한지

### Step 6: 파이프라인 실행

```bash
# 로컬 워크스테이션
./bin/run_pipeline.sh --cores 16

# 또는 SLURM 클러스터
sbatch bin/run_slurm.sh
```

### Step 7: 결과 확인

```bash
# 결과 디렉토리 구조
ls -R 3rd_analysis_results/

# Small variant 결과
ls -lh 3rd_analysis_results/small_variants/

# DMR 분석 결과
cat 3rd_analysis_results/dmr_analysis/dmr_results.csv | head

# 요약 리포트
open 3rd_analysis_results/analysis_summary_report.html
```

---

## 📊 예상 결과 구조

```
3rd_analysis_results/
├── small_variants/
│   ├── KTY9537.filtered.vcf.gz
│   ├── KTY9537.filtered.tsv          ← Excel에서 열기
│   ├── KTY9537.annotated.vcf.gz
│   ├── KTY9537.annotated.tsv         ← Excel에서 열기
│   ├── KTY9537.vep_stats.html
│   └── ... (다른 샘플들)
│
├── structural_variants/
│   ├── KTY9537.filtered_sv.vcf.gz
│   └── ...
│
├── dmr_analysis/
│   ├── dmr_results.csv                ← 주요 결과
│   ├── dmr_plots.pdf                  ← 시각화
│   ├── control_samples.txt
│   └── experimental_samples.txt
│
├── methylation/
│   ├── KTY9537.bw                     ← IGV 시각화용
│   ├── LDK6217.bw
│   └── methylation_summary.txt
│
├── logs/
│   ├── slivar/
│   ├── vep/
│   ├── svpack/
│   └── dmr/
│
└── analysis_summary_report.html       ← 전체 요약
```

---

## 🔧 트러블슈팅

### Issue 1: prepare_inputs.py 실행 오류

```bash
# 오류: "샘플 디렉토리를 찾을 수 없습니다"

# 해결: 실제 디렉토리 확인
ls /data_4tb/hifi-human-wgs-wdl-custom/batch_results/

# config.yaml의 샘플 ID와 폴더명 일치시키기
```

### Issue 2: filelist.csv에 누락된 파일

```bash
# 특정 샘플의 파일 직접 확인
ls -la /data_4tb/.../batch_results/KTY9537/*/out/

# 필수 파일:
# - phased_small_variant_vcf
# - phased_sv_vcf
# - cpg_combined_bed
# - cpg_combined_bw
```

### Issue 3: Snakemake 실행 오류

```bash
# 로그 확인
tail -f 3rd_analysis_results/logs/slivar/KTY9537.log

# 특정 규칙만 재실행
snakemake --use-singularity --cores 4 \
  3rd_analysis_results/small_variants/KTY9537.filtered.vcf.gz
```

---

## 💡 주요 팁

### 1. 특정 샘플만 먼저 테스트

```yaml
# config.yaml에서 하나만 설정
samples:
  control:
    - KTY9537
  experimental: []
```

```bash
python bin/prepare_inputs.py
./bin/run_pipeline.sh --dry-run
./bin/run_pipeline.sh --cores 8
```

### 2. 실행 시간 단축

```yaml
# 리소스 증가
resources:
  vep:
    threads: 16  # 기본 8 → 16
    mem_mb: 32000  # 기본 16000 → 32000
```

### 3. 중간 결과 보존

```bash
# Snakemake는 자동으로 중간 파일 삭제하지 않음
# 디스크 공간 확인
df -h
```

### 4. 병렬 실행

```bash
# 샘플별로 독립적으로 실행
# Snakemake가 자동으로 병렬화
./bin/run_pipeline.sh --cores 32  # 많은 코어 할당
```

---

## 📚 추가 문서

- [빠른 시작 가이드](docs/QUICKSTART.md)
- [WDL 결과 파일 자동 탐색](docs/PREPARE_INPUTS_GUIDE.md)
- [프로젝트 구조](PROJECT_STRUCTURE.md)
- [기술 요약](docs/SUMMARY.md)

---

## 🎓 실제 사용 예시

### 예시 1: 4개 샘플 분석

```bash
# 1. 설정
cat > config/config.yaml << 'EOF'
paths:
  batch_results_dir: "/data_4tb/hifi-human-wgs-wdl-custom/batch_results"
  ref_genome: "/reference/GRCh38.fa"
  vep_cache: "/data/vep_cache"
  
samples:
  control: [KTY9537, LDK6217]
  experimental: [15_JANG_SJ, 16_BYUN_HY]
EOF

# 2. 파일 탐색
python bin/prepare_inputs.py

# 3. 실행
./bin/run_pipeline.sh --cores 16

# 예상 시간: ~8시간 (16코어 기준)
```

### 예시 2: 대규모 배치 (20개 샘플)

```bash
# 1. 샘플 목록 준비
samples:
  control: [S01, S02, ..., S10]
  experimental: [S11, S12, ..., S20]

# 2. 클러스터 실행
sbatch bin/run_slurm.sh

# 예상 시간: ~12시간 (클러스터 병렬 처리)
```

---

**작성일**: 2026-01-29  
**버전**: 1.0
