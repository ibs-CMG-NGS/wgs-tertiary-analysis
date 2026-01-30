# Species-Specific Configuration Guide

## 지원하는 종 (Species)

본 파이프라인은 **Human**과 **Mouse** WGS 데이터 분석을 모두 지원합니다.

### 🧬 Human WGS 분석

**설정 파일:** `config/config.human.yaml`

**주요 특징:**
- gnomAD v3/v4 인구 집단 빈도 데이터베이스 사용
- 희귀 변이 필터링 (AF < 0.01)
- HPO 용어 기반 임상 우선순위화
- VEP 플러그인: CADD, REVEL

**빈도 데이터베이스 다운로드:**
```bash
# gnomAD v3.1.2 for GRCh38 (약 40GB)
wget https://github.com/brentp/slivar/releases/download/v0.3.1/gnomad.hg38.v3.1.2.zip
```

**설정 예시:**
```yaml
parameters:
  general:
    species: "human"
    reference_version: "GRCh38"
  
  slivar:
    max_af: 0.01
    frequency_db:
      human: "/data/references/gnomad.hg38.v3.1.2.zip"
    frequency_field:
      human: "gnomad_popmax_af"
```

---

### 🐭 Mouse WGS 분석

**설정 파일:** `config/config.mouse.yaml`

**주요 특징:**
- Mouse Genomes Project (MGP) 빈도 데이터베이스 (선택)
- 품질 중심 필터링 (GQ ≥ 30, DP ≥ 20)
- 기능 예측 중심 분석
- VEP 플러그인: dbNSFP

**빈도 데이터베이스 다운로드 (선택):**
```bash
# MGP v5 (Mouse 30+ strains)
wget ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP.vcf.gz
wget ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP.vcf.gz.tbi
```

**설정 예시:**
```yaml
parameters:
  general:
    species: "mouse"
    reference_version: "GRCm39"
  
  slivar:
    max_af: 0.05
    min_gq: 30        # 더 높은 품질 기준
    min_dp: 20
    frequency_db:
      mouse: "/data/references/mgp.v5.merged.snps_all.dbSNP.vcf.gz"
    frequency_field:
      mouse: "AF"
```

---

## 사용 방법

### Human 분석 실행
```bash
# 1. 설정 파일 복사 및 수정
cp config/config.human.yaml config/config.yaml
nano config/config.yaml  # 경로 수정

# 2. 파이프라인 실행
bash bin/run_pipeline.sh --config config/config.yaml --cores 16
```

### Mouse 분석 실행
```bash
# 1. 설정 파일 복사 및 수정
cp config/config.mouse.yaml config/config.yaml
nano config/config.yaml  # 경로 수정

# 2. 파이프라인 실행
bash bin/run_pipeline.sh --config config/config.yaml --cores 16
```

---

## 종별 차이점 요약

| 항목 | Human | Mouse |
|-----|-------|-------|
| **빈도 DB** | gnomAD v3 (필수 권장) | MGP v5 (선택) |
| **max_af** | 0.01 (1%) | 0.05 (5%) 또는 미사용 |
| **min_gq** | 20 | 30 |
| **min_dp** | 10 | 20 |
| **HPO terms** | 사용 | 미사용 |
| **VEP plugins** | CADD, REVEL | dbNSFP |
| **필터링 전략** | 빈도 + 기능 | 품질 + 기능 |

---

## 빈도 데이터베이스 없이 실행

두 종 모두 빈도 DB 없이 실행 가능합니다:

```yaml
slivar:
  frequency_db:
    human: ""
    mouse: ""
```

이 경우 품질 기반 필터링(GQ, DP)만 수행됩니다.
