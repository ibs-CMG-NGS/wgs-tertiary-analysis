# H2O2 그룹비교 개선 계획 (joint calling + 샘플단위 요약통계 + 유전자셋 검정)

**작성일:** 2026-07-15
**배경:** 기존 `cohort_sv_burden`/`cohort_variant_burden`은 (1) single-sample VCF 기반이고 (2) 유전자 단위로
수백~천 개를 개별 검정해 n=3/group 다중검정 부담이 커 FDR 통과가 사실상 불가능했음(§ 대화 참고,
SV burden 최저 FDR 0.29, VEP burden 최저 FDR 1.0). 이를 보완하는 3가지 방법을 순서대로 적용한다.

---

## 1단계: 입력을 joint-called VCF로 교체

### 현황 확인 완료
- 4웰(A/B/C/D) 모두 joint calling(GLnexus+Sawfish, `joint.wdl`) 완료돼 있음:
  - `/mnt/hdd1_1tb/h2o2-wgs/joint-H2O2-A01/20260714_054055_joint`
  - `/mnt/hdd1_1tb/h2o2-wgs/joint-H2O2-B01/20260714_062818_joint`
  - `/mnt/hdd2_1tb/h2o2-wgs/joint-H2O2-C01/20260714_071318_joint`
  - `/mnt/hdd2_1tb/h2o2-wgs/joint-H2O2-D01/20260714_075510_joint`
- 각 웰 디렉토리의 `out/split_joint_small_variant_vcfs/{0,1,2}/`, `out/split_joint_structural_variant_vcfs/{0,1,2}/`
  안에 샘플별로 분리된 VCF 존재. 파일명 자체에 sample_id가 포함돼 있어(`{sample_id}.{well}-joint.joint.GRCm39.*.vcf.gz`)
  인덱스 번호(0/1/2) 순서에 의존하지 않고 `sample_id.*.vcf.gz` 글롭으로 안전하게 찾을 수 있음.
  인덱스 파일(`.tbi`)은 `_indices` 서픽스 디렉토리에 동일 이름으로 존재.
- methylation(cpg_combined_bed 등)과 TRGT VCF는 joint calling 대상이 아니므로(WDL이 small/SV variant만 joint
  genotyping) 기존 개별 파일 그대로 유지 — 이번 교체는 `phased_small_variant_vcf`/`phased_sv_vcf` 두 필드만 해당.

### 작업
1. `bin/make_h2o2_filelist.py`를 확장(또는 `--source joint` 플래그 추가)해서 joint VCF 경로로 채운
   `filelist_h2o2_joint.csv` 생성. 기존 개별-콜 파일리스트(`filelist_h2o2.csv`)는 그대로 두고 별도 파일로 관리
   (비교를 위해 두 버전 다 보존).
2. `config/config_joint.yaml`(config.yaml 복사 + `filelist_csv: filelist_h2o2_joint.csv`, `output_dir` 분리)
   생성.
3. `cohort_sv_burden`, `annotate_vep`→...→`cohort_variant_burden` 체인 재실행 (output_dir가 다르므로 전체
   파이프라인 재실행 필요 — filter_sv/filter_small_variants부터).
4. 결과 비교: 기존 single-sample 버전과 joint 버전의 top 유전자 목록·raw p-value·효과크기(Cramér's V 등) 대조.

---

## 2단계: 샘플 단위 요약통계 비교 (다중검정 부담 최소화)

기존 세션에서 이미 사용했던 방법(하나의 검정만 수행, n=3×4군 그대로 사용)을 정식 Snakemake rule로 편입.

### 지표
- VAF 분포 hom-peak 비율 (0.9~1.0 구간 비율/피크 높이)
- 샘플당 총 PASS 변이 수, 총 SV 개수
- Joint genotype 불일치율 (joint VCF 필요 → 1단계 완료 후 가능)
- (선택) LOH 후보 변이 개수 — 1단계의 joint VCF 필요

### 통계
- 그룹(A/B/C/D) 간 Kruskal-Wallis
- 그룹이 실제 처치강도 순서(A<C<B<D 등)를 가진다면 Kendall's W(순서 일치도, 이전 분석에서 이미 사용)
- 검정 수가 지표 개수(4~5개) 뿐이므로 BH 보정해도 크게 불리하지 않음

### 작업
1. 신규 스크립트 `scripts/sample_level_summary.R` (또는 python): 샘플별 VCF에서 위 지표 계산 → 1행/샘플 CSV.
2. 신규 스크립트 `scripts/group_summary_compare.R`: Kruskal-Wallis + Kendall's W 계산 + boxplot/순위 시각화.
3. Snakefile에 rule 추가: `compute_sample_metrics`(per-sample) → `cohort_summary_compare`(집계).
4. config에 지표 목록/파라미터 추가.

---

## 3단계: 유전자셋(pathway) 단위 집계 검정

개별 유전자 대신 생물학적으로 관련된 유전자 묶음(산화스트레스 반응, DNA 손상 복구 등) 전체를 하나의 단위로
검정 — 검정 개수를 수백~수천 개에서 소수(유전자셋 개수)로 줄여 다중검정 부담을 크게 낮춤.

### 유전자셋 확보
- 로컬에 mouse GO 주석 패키지(org.Mm.eg.db 등) 없음 확인됨 → MGI GO annotation 파일
  (`http://current.geneontology.org/annotations/mgi.gaf.gz`)을 다운로드해서 관련 GO term의 유전자 목록 추출:
  - GO:0006979 (response to oxidative stress)
  - GO:0034599 (cellular response to oxidative stress)
  - GO:0006281 (DNA repair)
  - GO:0006974 (DNA damage response)
  - (astrocyte 관련 있으면 추가 검토)
- gene_association 파일은 MGI ID 기반이므로 GENCODE vM36 gene_name과 매핑 필요 (MGI ID ↔ symbol은 GAF 파일에
  symbol 컬럼이 포함돼 있어 별도 변환 불요, symbol로 바로 매칭)

### 통계
- 유전자셋 내 "이 세트의 유전자 중 하나라도 SV/변이가 있는 샘플" 비율을 군별 비교 (Fisher/chi-square, 세트당 1회)
- 또는 세트 내 총 이벤트 수(유전자별 카운트 합)를 군별 비교 (Kruskal-Wallis, 연속값이라 검정력 더 좋음)
- 세트 개수(4~5개)만큼만 검정 → BH 보정 부담 거의 없음

### 작업
1. `scripts/fetch_go_genesets.sh` 또는 python: MGI GAF 다운로드 + 관련 GO term 필터 → `resources/genesets_oxidative_dna_repair.json`
2. `scripts/geneset_burden.R`: 기존 `sv_burden.R`/`variant_burden.R`의 gene×sample 행렬을 재사용해 유전자셋 단위로
   집계 후 검정 (신규 유전자별 파싱 불필요, 기존 산출물 재활용 가능하도록 설계).
3. Snakefile에 rule 추가: `cohort_geneset_burden` (SV/VEP 각각, 또는 통합).

---

## 실행 순서 및 확인 지점
1. 1단계 완료 → 결과(유전자 목록/효과크기 변화) 보고
2. 2단계 완료 → 결과(Kruskal-Wallis/Kendall's W) 보고
3. 3단계 완료 → 결과(유전자셋별 검정) 보고
각 단계 종료 시 다음 단계로 자동 진행(사용자가 이미 "순서대로 진행" 승인함), 단 예상 밖 이슈 발견 시 즉시 보고.
