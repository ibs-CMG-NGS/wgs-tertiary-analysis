# DMR 분석 최적화 가이드

## 문제점
DMR 분석은 매번 30분~1시간 이상 소요되지만, 입력 데이터가 변하지 않으면 결과도 동일합니다. 
기존 구조에서는 매번 파이프라인을 실행할 때마다 DMR 분석이 재실행되었습니다.

## 해결 방법

### 🎯 **DMR 분석을 선택적 타겟으로 변경**

이제 DMR 분석은 **기본적으로 실행되지 않으며**, 필요할 때만 명시적으로 실행합니다.

---

## 실행 방법

### 1. **일반 파이프라인 실행 (DMR 제외)**
```bash
bash bin/run_pipeline.sh --config config/config_hifisolve.yaml --cores 16
```
- Small variant 필터링/주석
- Structural variant 필터링
- 메틸화 파일 정리

**✅ DMR 분석 건너뜀 → 빠른 실행**

---

### 2. **DMR 분석 포함하여 실행**
```bash
bash bin/run_pipeline.sh --config config/config_hifisolve.yaml --cores 16 --dmr
```
- 모든 분석 + DMR 분석

---

### 3. **DMR 분석만 실행** (이미 필터링 완료된 경우)
```bash
# 규칙 이름으로 실행
bash bin/run_pipeline.sh --config config/config_hifisolve.yaml --target run_dmr_analysis

# 또는 파일 경로로 실행
snakemake --cores 16 output/dmr_analysis/dmr_results.csv output/dmr_analysis/dmr_plots.pdf
```
- 기존 메틸화 데이터를 사용하여 DMR 분석만 재실행
- 매개변수 변경 후 재분석할 때 유용 (예: p-value, delta 임계값 조정)

---

### 4. **직접 Snakemake 명령 사용**
```bash
# DMR 분석만 실행
snakemake --cores 16 output/dmr_analysis/dmr_results.csv output/dmr_analysis/dmr_plots.pdf

# DMR 분석 포함 전체 실행
snakemake --cores 16 all output/dmr_analysis/dmr_results.csv output/dmr_analysis/dmr_plots.pdf
```

---

## Snakemake 캐싱 동작

Snakemake는 **파일 타임스탬프**를 기반으로 재실행 여부를 결정합니다:

### ✅ **재실행되지 않는 경우**
- 출력 파일(`output/dmr_analysis/dmr_results.csv`, `dmr_plots.pdf`)이 존재
- 입력 파일(메틸화 BED 파일)이 출력 파일보다 **오래됨**
- R 스크립트(`scripts/dmr_analysis.R`)가 변경되지 않음

### ⚠️ **재실행되는 경우**
- 출력 파일이 없거나 삭제됨
- 입력 파일이 업데이트됨 (더 최신)
- R 스크립트가 수정됨
- `--forcerun` 옵션 사용

---

## 매개변수 변경 시 재분석

DMR 분석 매개변수만 변경한 경우:

```bash
# 1. config.yaml에서 매개변수 수정
#    dmr_analysis:
#      p_value: 0.01  # 0.05에서 변경
#      min_delta: 0.2  # 0.1에서 변경

# 2. DMR 분석 강제 재실행
snakemake --cores 16 --forcerun run_dmr_analysis dmr_analysis
```

---

## 시간 절약 효과

| 시나리오 | 기존 (매번 실행) | 최적화 (선택적 실행) |
|---------|-----------------|---------------------|
| VEP 주석만 재실행 | ~1시간 (DMR 포함) | ~15분 (DMR 제외) |
| 설정 변경 후 재실행 | ~1시간 | ~15분 |
| DMR 매개변수 조정 | ~1시간 (전체) | ~30분 (DMR만) |

---

## 문제 해결

### Q: 이미 DMR 결과가 있는데 다시 실행됨
```bash
# 출력 파일 확인
ls -lh output/dmr_analysis/

# Snakemake가 최신 상태로 인식하는지 확인
snakemake --dry-run output/dmr_analysis/dmr_results.csv
```

### Q: 입력 파일은 그대로인데 강제로 재실행하고 싶음
```bash
# 출력 파일 삭제 후 재실행
rm output/dmr_analysis/dmr_results.csv output/dmr_analysis/dmr_plots.pdf
snakemake --cores 16 output/dmr_analysis/dmr_results.csv output/dmr_analysis/dmr_plots.pdf

# 또는 --forcerun 사용
snakemake --cores 16 --forcerun run_dmr_analysis output/dmr_analysis/dmr_results.csv
```

---

## 요약

**이전 방식:**
```bash
bash bin/run_pipeline.sh  # 항상 DMR 포함 → 느림
```

**새로운 방식:**
```bash
# 빠른 실행 (DMR 제외)
bash bin/run_pipeline.sh

# 필요할 때만 DMR 포함
bash bin/run_pipeline.sh --dmr

# DMR만 재실행 (규칙 이름)
bash bin/run_pipeline.sh --target run_dmr_analysis

# 또는 파일 경로 직접 지정
snakemake --cores 16 output/dmr_analysis/dmr_results.csv
```

**✅ 이제 DMR 분석 없이 파이프라인을 빠르게 실행할 수 있습니다!**
