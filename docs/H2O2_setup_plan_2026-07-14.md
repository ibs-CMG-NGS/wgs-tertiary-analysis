# H2O2 astrocyte WGS — wgs-tertiary-analysis 파이프라인 연동 계획

**작성일:** 2026-07-14
**목적:** 지금까지 수동/스크래치 스크립트로 진행한 H2O2 12샘플(A/B/C/D × 3배치) 군별 비교 분석을,
이미 구축돼 있던 `wgs-tertiary-analysis` Snakemake 파이프라인(3차 분석, 가변+시각화)으로 정식 이관.
`hifi-human-wgs-wdl-custom`(2차, 고정 WDL)과 역할을 분리하는 기존 base-mod 세션의 아키텍처 원칙을
WGS 쪽에도 동일하게 적용.

관련 문서: [`H2O2_variant_count_qc_analysis_2026-07-14.md`](../../hifi-human-wgs-wdl-custom/H2O2_variant_count_qc_analysis_2026-07-14.md), [`H2O2_key_findings_summary_2026-07-14.md`](../../hifi-human-wgs-wdl-custom/H2O2_key_findings_summary_2026-07-14.md)

---

## 1. 현재 상태 조사 결과

### 1.1 `wgs-tertiary-analysis`는 이미 상당히 구축돼 있으나 미실행 상태

- git 저장소, `Snakefile`(14 phase), `config/`, `envs/`, `bin/`, `docs/` 전부 존재.
- **`groups:` 다중 그룹 설정을 이미 지원** — control/experimental 2군 하드코딩이 아니라 임의 개수 그룹 가능. 우리 A/B/C/D 4그룹 구조에 그대로 맞음.
- **`filelist.csv`로 파일 경로를 직접 매핑**하는 구조 — 샘플이 물리적으로 한 디렉토리에 없어도(우리 경우 `JJ_dis_8tb`/`hdd1_1tb`/`hdd2_1tb` 3곳에 분산) 문제없이 연결 가능.
- 전용 conda env(`wgs-tertiary-pipeline`)가 이미 생성되어 있음(다른 목적으로 이미 사용 중, bcftools 포함 확인됨).
- **`3rd_analysis_results/`, `filelist.csv` 등 실제 실행 산출물이 전혀 없음** → 이 파이프라인은 스캐폴딩만 되어 있고 **한 번도 end-to-end로 돌아간 적이 없는 미검증 상태**.
- `config/config.mouse.yaml`은 범용 예시(`C57BL6_WT_1` 등 placeholder)로만 채워져 있어, H2O2 실제 경로로 새로 작성 필요.

### 1.2 Rule 목록과 의존성 (14 phase)

| Phase | Rule | VEP 필요? | 비고 |
|---|---|---|---|
| — | `filter_small_variants` | 아니오 | slivar 필터링 |
| — | `annotate_vep` / `vep_to_tsv` / `extract_vep_canonical` | **예** | mouse VEP 캐시 필요 (미보유) |
| — | `filter_sv` | 아니오 | |
| — | `prepare_methylation_data` / `run_dmr_analysis` / `merge_methylation_plots` | 아니오 | DSS 기반 DMR, R |
| — | `annotate_sv_consequence` | 아니오 (GFF3만 필요, VEP 아님) | svpack consequence, **Ensembl GFF3 포맷 필요** (기존에 base-mod에서 쓴 GENCODE GTF와 포맷이 다를 수 있어 확인 필요) |
| — | `analyze_trgt` | 아니오 | 반복서열(TRGT) |
| — | `analyze_asm` | 아니오 | assembly 관련 (우리는 `run_assembly=false`로 돌렸으므로 스킵 또는 조건부 처리 필요) |
| — | `generate_summary_report` | 아니오 | |
| 11 | `cohort_variant_burden` | **예** (VEP 주석 TSV 필요) | 소변이 유전자별 군간 burden test |
| 12 | `cohort_sv_burden` | 아니오 (annotate_sv_consequence만 필요) | SV 유전자별 군간 burden test |
| 13 | `trgt_group_compare` | 아니오 | 반복서열 군간 비교 |
| 14 | `prepare_group_methylation_list` / `dmr_pairwise` | 아니오 | **모든 군 쌍(6개 조합)에 대해 DMR 자동 실행** |

### 1.3 핵심 판단: VEP 없이도 상당 부분 진행 가능

**VEP가 실제로 필요한 건 `cohort_variant_burden`(소변이 유전자별 burden test) 하나뿐.** SV burden, TRGT 군간 비교, DMR pairwise는 전부 VEP 없이 돌아감 — 이번 세션에서 이미 발견한 이질성/LOH 신호(§3~5, 상세 문서 참고)를 **정식 통계 파이프라인으로 재현·확장**하는 데 필요한 건 대부분 VEP 비의존 경로로 커버됨.

---

## 2. 실행 계획 (VEP 제외, 1단계)

### 2.1 `filelist.csv` 작성

`bin/prepare_inputs.py`는 단일 `batch_results_dir` 가정이라 우리 상황(3개 드라이브 분산)에는 그대로 안 맞음 → `prepare_inputs.py`의 `find_wdl_output_dir`/`find_required_files` 함수를 재사용하되, 12개 샘플 각각의 실제 경로를 직접 지정하는 별도 스크립트로 `filelist.csv` 생성 (기존 `prepare_inputs.py` 파일 자체는 수정하지 않음).

12개 샘플 실행 디렉토리 (이미 확인된 경로):

| 샘플 | 그룹 | 경로 |
|---|---|---|
| H2O2-A01-ctrl1 | A | `/mnt/JJ_dis_8tb/h2o2-wgs/H2O2-A01-ctrl1/20260517_163205_humanwgs_singleton` |
| H2O2-B01 | B | `/mnt/JJ_dis_8tb/h2o2-wgs/H2O2-B01/20260518_163817_humanwgs_singleton` |
| H2O2-C01 | C | `/mnt/JJ_dis_8tb/h2o2-wgs/H2O2-C01/20260519_205414_humanwgs_singleton` |
| H2O2-D01 | D | `/mnt/JJ_dis_8tb/h2o2-wgs/H2O2-D01/20260520_184847_humanwgs_singleton` |
| H2O2-2-A01 | A | `/mnt/JJ_dis_8tb/h2o2-wgs/H2O2-2-A01/20260523_170354_humanwgs_singleton` |
| H2O2-2-B01 | B | `/mnt/JJ_dis_8tb/h2o2-wgs/H2O2-2-B01/20260524_171651_humanwgs_singleton` |
| H2O2-2-C01 | C | `/mnt/JJ_dis_8tb/h2o2-wgs/H2O2-2-C01/20260525_182243_humanwgs_singleton` |
| H2O2-2-D01 | D | `/mnt/JJ_dis_8tb/h2o2-wgs/H2O2-2-D01/20260526_181743_humanwgs_singleton` |
| H2O2-3-A01 | A | `/mnt/hdd1_1tb/h2o2-wgs/H2O2-3-A01/20260710_102952_humanwgs_singleton` |
| H2O2-3-B01 | B | `/mnt/hdd1_1tb/h2o2-wgs/H2O2-3-B01/20260711_070817_humanwgs_singleton` |
| H2O2-3-C01 | C | `/mnt/hdd2_1tb/h2o2-wgs/H2O2-3-C01/20260712_023952_humanwgs_singleton` |
| H2O2-3-D01 | D | `/mnt/hdd2_1tb/h2o2-wgs/H2O2-3-D01/20260713_134139_humanwgs_singleton` |

### 2.2 `config.yaml` 작성

`config/config.mouse.yaml`을 베이스로 복사, 다음을 실제 값으로 교체:
- `paths.batch_results_dir`: 사용 안 함(filelist.csv 모드이므로 무시됨, 명시만)
- `paths.ref_genome`: `/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/mouse_GRCm39.fasta`
- `paths.output_dir`: `./h2o2_analysis_results`
- `paths.filelist_csv`: `filelist_h2o2.csv`
- `samples.groups`: A/B/C/D 4그룹, 위 표의 12개 샘플 ID
- `parameters.general.species/reference_version`: mouse / GRCm39
- VEP 관련 섹션은 일단 그대로 두되(1단계에서는 미사용), `gff3_file` 경로는 §2.3에서 확인 후 채움

### 2.3 GFF3 파일 확인/준비 (`annotate_sv_consequence`용)

svpack consequence가 요구하는 정확한 GFF3 포맷 확인 필요. 기존 base-mod 작업에서 받은 GENCODE GTF(`gencode.vM36.annotation.sorted.gtf.bgz`)로 대체 가능한지, 아니면 Ensembl GFF3를 별도로 받아야 하는지 실행 전 확인.

### 2.4 Dry-run 검증

```bash
conda activate wgs-tertiary-pipeline
cd /data_4tb/shared/wgs-tertiary-analysis
snakemake --configfile config/config.yaml -n \
  cohort_sv_burden trgt_group_compare dmr_pairwise generate_summary_report
```

VEP 의존 타깃(`cohort_variant_burden` 등)은 이번 1단계 목표에서 명시적으로 제외.

### 2.5 실제 실행 + 결과를 상세 문서(§6 이후)에 반영

---

## 3. 2단계 (추후, 별도 승인 후 진행)

- Mouse VEP 캐시 다운로드 (`mus_musculus` GRCm39, 수 GB) 및 `cohort_variant_burden` 실행
- `analyze_asm` 관련 rule은 `run_assembly=false`로 돌렸으므로 스킵 여부 확인 (assembly 산출물 없음)

## 4. 3단계 — base-mod 분석 (사용자 지시대로 이 작업 이후 진행)

---

## 5. 리스크 / 확인 필요 사항

- `annotate_sv_consequence`가 요구하는 GFF3 포맷이 기존 GENCODE GTF와 호환되는지 사전 확인 필요 (§2.3)
- 파이프라인이 한 번도 실행된 적이 없어, dry-run 이후에도 실제 실행 시 예상 못 한 에러가 나올 가능성 있음 — 단계적으로(rule 하나씩) 검증 권장
- `envs/dmr.yaml`, `environment.yaml` 등 conda env가 실제로 필요한 R 패키지(DSS 등)를 다 포함하는지 미확인 — 첫 실행 시 확인 필요
