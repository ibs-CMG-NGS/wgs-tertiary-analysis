# WDL 결과 파일 자동 탐색 가이드

## 📋 개요

WDL 파이프라인의 복잡한 디렉토리 구조에서 3차 분석에 필요한 파일들을 자동으로 찾아주는 `prepare_inputs.py` 스크립트 사용 가이드입니다.

## 🗂️ WDL 결과 디렉토리 구조

```
batch_results/
├── KTY9537/
│   └── 20260120_101704_humanwgs_singleton/
│       ├── out/
│       │   ├── phased_small_variant_vcf
│       │   ├── phased_sv_vcf
│       │   ├── cpg_combined_bed
│       │   ├── cpg_combined_bw
│       │   └── ... (기타 파일들)
│       ├── call-*/ (워크플로우 단계별 디렉토리)
│       └── workflow.log
├── LDK6217/
│   └── 20260125_093412_humanwgs_singleton/
│       └── out/
│           └── ...
└── 15_JANG_SJ/
    └── 20260128_144521_humanwgs_singleton/
        └── out/
            └── ...
```

## 🚀 사용 방법

### 1단계: config.yaml 설정

```yaml
paths:
  # WDL batch_results 최상위 디렉토리 지정
  batch_results_dir: "/data_4tb/hifi-human-wgs-wdl-custom/batch_results"
  
  # 기타 설정...
  ref_genome: "/path/to/GRCh38.fa"
  vep_cache: "/path/to/vep_cache"

samples:
  control:
    - KTY9537
    - LDK6217
  experimental:
    - 15_JANG_SJ
    - 16_BYUN_HY
```

### 2단계: 입력 파일 자동 탐색

```bash
# 기본 실행 (config.yaml 사용)
python bin/prepare_inputs.py

# 또는 경로 직접 지정
python bin/prepare_inputs.py -b /data_4tb/hifi-human-wgs-wdl-custom/batch_results

# Dry-run 모드 (파일 생성 없이 테스트)
python bin/prepare_inputs.py --dry-run
```

### 실행 결과

```
================================================================================
WDL 결과 파일 탐색 시작
================================================================================

Batch Results 디렉토리: /data_4tb/hifi-human-wgs-wdl-custom/batch_results

[CONTROL] 그룹 처리 중...
  └─ KTY9537 탐색 중... ✓ 완료
  └─ LDK6217 탐색 중... ✓ 완료

[EXPERIMENTAL] 그룹 처리 중...
  └─ 15_JANG_SJ 탐색 중... ✓ 완료
  └─ 16_BYUN_HY 탐색 중... ✓ 완료

파일 목록 저장 중: filelist.csv
✓ 완료: 4개 샘플의 파일 정보 저장

================================================================================
처리 요약
================================================================================
  총 샘플 수: 4
  성공: 4
  경고/오류: 0
  출력 파일: filelist.csv
================================================================================

✓ 모든 샘플 처리 완료!

다음 단계:
  1. filelist.csv 파일을 확인하세요
  2. Snakefile을 실행하여 3차 분석을 시작하세요:
     ./bin/run_pipeline.sh --cores 16
```

### 3단계: 생성된 filelist.csv 확인

```csv
sample_id,group,phased_small_variant_vcf,phased_small_variant_vcf_index,phased_sv_vcf,phased_sv_vcf_index,cpg_combined_bed,cpg_combined_bed_index,cpg_combined_bw
KTY9537,control,/data_4tb/.../KTY9537/20260120_101704_humanwgs_singleton/out/phased_small_variant_vcf,,/data_4tb/.../KTY9537/.../out/phased_sv_vcf,,/data_4tb/.../out/cpg_combined_bed,,/data_4tb/.../out/cpg_combined_bw
LDK6217,control,...,...,...,...,...,...,...
15_JANG_SJ,experimental,...,...,...,...,...,...,...
16_BYUN_HY,experimental,...,...,...,...,...,...,...
```

### 4단계: 파이프라인 실행

```bash
# filelist.csv가 있으면 Snakemake가 자동으로 사용합니다
./bin/run_pipeline.sh --cores 16
```

## 🔍 스크립트 동작 원리

### 1. 샘플 디렉토리 탐색
- `batch_results/{sample_id}/` 디렉토리 찾기

### 2. WDL 실행 디렉토리 탐색
- `*_humanwgs_singleton` 패턴의 디렉토리 찾기
- 여러 개 있을 경우 가장 최근 것 선택 (timestamp 기준)

### 3. 필수 파일 위치 확인
- `out/phased_small_variant_vcf`
- `out/phased_sv_vcf`
- `out/cpg_combined_bed`
- `out/cpg_combined_bw`

### 4. 파일 경로 CSV 저장
- 각 샘플의 실제 파일 경로 기록
- Snakemake가 이를 읽어서 분석 수행

## ⚙️ 고급 옵션

### 커스텀 출력 파일명

```bash
python bin/prepare_inputs.py -o my_samples.csv
```

그 후 `config.yaml`에서:

```yaml
paths:
  filelist_csv: "my_samples.csv"
```

### 특정 그룹만 처리

`config.yaml`에서 필요한 샘플만 나열:

```yaml
samples:
  control:
    - KTY9537  # 이 샘플만 처리
  experimental: []
```

### 여러 batch_results 디렉토리 처리

각 디렉토리에 대해 별도로 실행 후 CSV 병합:

```bash
python bin/prepare_inputs.py -b /path/to/batch1 -o batch1.csv
python bin/prepare_inputs.py -b /path/to/batch2 -o batch2.csv

# CSV 병합
cat batch1.csv > filelist.csv
tail -n +2 batch2.csv >> filelist.csv
```

## 🐛 문제 해결

### 오류: "샘플 디렉토리를 찾을 수 없습니다"

**원인**: `batch_results/{sample_id}` 디렉토리가 없음

**해결**:
```bash
# 실제 디렉토리 확인
ls /data_4tb/hifi-human-wgs-wdl-custom/batch_results/

# config.yaml의 샘플 ID와 디렉토리명 일치시키기
```

### 오류: "WDL 실행 디렉토리를 찾을 수 없습니다"

**원인**: `*_humanwgs_singleton` 디렉토리가 없음

**해결**:
```bash
# WDL 실행이 완료되었는지 확인
ls /data_4tb/hifi-human-wgs-wdl-custom/batch_results/KTY9537/

# 출력: 20260120_101704_humanwgs_singleton (이런 형식이어야 함)
```

### 경고: "필수 파일 누락"

**원인**: `out/` 디렉토리에 일부 파일이 없음

**해결**:
- WDL 파이프라인이 정상적으로 완료되었는지 확인
- `workflow.log` 확인하여 오류 체크
- 누락된 파일만 제외하고 나머지는 정상 처리됨

## 📊 워크플로우 다이어그램

```
config.yaml 작성
      ↓
python bin/prepare_inputs.py
      ↓
filelist.csv 생성
      ↓
Snakemake 실행
      ↓
3차 분석 결과
```

## 💡 팁

### 1. 파일 존재 여부만 빠르게 확인

```bash
python bin/prepare_inputs.py --dry-run
```

### 2. 특정 샘플만 테스트

`config.yaml`에서 하나의 샘플만 설정 후:

```bash
python bin/prepare_inputs.py
./bin/run_pipeline.sh --dry-run
```

### 3. 자동화 스크립트

```bash
#!/bin/bash
# 전체 파이프라인 자동 실행

# 1. 입력 파일 탐색
python bin/prepare_inputs.py || exit 1

# 2. 환경 확인
bash bin/setup_check.sh || exit 1

# 3. Dry-run
./bin/run_pipeline.sh --dry-run || exit 1

# 4. 실제 실행
./bin/run_pipeline.sh --cores 16
```

---

**작성일**: 2026-01-29  
**버전**: 1.0
