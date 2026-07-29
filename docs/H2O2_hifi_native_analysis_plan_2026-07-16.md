# H2O2 astrocyte WGS — HiFi 롱리드 고유 기능 활용 분석 계획

**작성일:** 2026-07-16
**배경:** 지금까지(§0~§5, `H2O2_tertiary_analysis_results_2026-07-16.md`)의 분석은 bulk VAF 분포, 유전자별 presence/absence burden 등 **short-read WGS로도 동일하게 할 수 있는 방법**이었다. HiFi 롱리드 고유 기능(phasing, 분자 단위 read, allele-specific methylation)을 실제로 활용하지 않았다는 지적에 따라, 아래 4가지를 순서대로 검증한다. 이미 완료된 항목(카피수 안정성, §2.3)은 제외하고 나머지 3개를 진행한다.

**사전 조사 결과(오늘 확인):**
- phased VCF에 `PS`(phase set) FORMAT 태그 존재 — 최대 블록 4,012개 변이, 대부분의 이형접합 사이트가 유효한 PS 값을 가짐.
- **`MF`/`MD`/`MT`(변이별 methylation fraction/depth/type) FORMAT 필드가 VCF 헤더엔 정의돼 있으나 실제로는 313,942개 이형접합 레코드 전부 비어있음(`.`)** — 애초에 계획했던 "VCF 안에서 genotype-methylation 직접 연결" 경로는 **사용 불가로 확인, 폐기**. 대신 기존 `analyze_asm`(hap1/hap2 BED 비교) 경로로 진행.
- `merged_haplotagged_bam`에 `HP:i:` 태그(read별 haplotype 1/2 배정) 존재 확인. 단 **BAM 인덱스(.bai)가 없어 region 기반 조회 불가 — 사전에 `samtools index` 필요**(12개 BAM, 인덱스 파일 자체는 작아 디스크 부담 적음).

---

## 1단계: Haplotype 방향성 allelic imbalance (PS 블록 기반)

**목적:** 진짜 클론 선택이 있다면 특정 haplotype 쪽으로 **일관된 방향의** allelic imbalance가 phase block 내에서 여러 이형접합 사이트에 걸쳐 나타나야 한다. 개별 SNP를 독립적으로 보는 지금 방식(§2)으론 이걸 못 잡는다.

**방법:**
1. 각 샘플의 phased VCF에서 이형접합(`0|1`/`1|0`) 사이트를 `PS` 태그로 그룹핑.
2. 각 PS 블록 내에서, phase-0(`|`앞) allele의 VAF가 0.5에서 벗어나는 방향이 블록 내 여러 사이트에서 얼마나 일관되는지 검정(예: 블록 내 사이트별 "phase-0 allele VAF - 0.5"의 부호가 통계적으로 한쪽에 쏠리는지, sign test 또는 평균 편향).
3. 큰 PS 블록(사이트 수 많은 것) 위주로 분석 — 작은 블록은 노이즈가 큼.
4. 샘플별로 "게놈 전체에서 편향된 큰 PS 블록의 비율"을 지표화 → §2 방식대로 군간 비교(Kruskal-Wallis).

**입력:** 12샘플의 `phased_small_variant_vcf`
**출력:** `scripts/haplotype_imbalance.py`, `h2o2_analysis_results/cohort/hap_imbalance/{sample}.hap_imbalance.tsv`, 시각화: `scripts/plot_hifi_native_results.R` → `h2o2_analysis_results/cohort/hap_imbalance_plots.pdf`

### 결과 (완료, 2026-07-22)

큰 PS 블록(≥10 이형접합 사이트) 내 phase-0 allele VAF 편향의 이항검정(raw p<0.05) 비율(`frac_biased`)을 군별 평균 낸 결과:

| 군 | frac_biased 평균 | mean_abs_bias 평균 |
|---|---|---|
| A | 0.120 | 0.052 |
| C | 0.349 | 0.076 |
| B | 0.627 | 0.100 |
| D | 0.708 | 0.111 |

**A<C<B<D — 기존 이질성 그라디언트와 정확히 일치.** Kruskal-Wallis H=9.974(df=3), n=3×4군 완전분리와 같은 통계량(다른 지표들과 동일한 H값 — n=3/group 완전분리의 특성). Haplotype phasing 정보를 실제로 활용한 첫 검증이며, 개별 SNP 단위가 아니라 phase block(여러 SNP가 연결된 haplotype 구간) 단위로 일관된 편향을 봤다는 점에서 §2의 hom-peak%와는 원리가 다른 독립적 확인.

---

## 2단계: Allele-specific methylation (ASM) — 기존 파이프라인 rule 실행

**목적:** §2.3(카피수)과 같은 논리로, 메틸화 차원에서도 이질성 그라디언트가 재현되는지 확인. 정정: 처음 제안 때 "genotype과 methylation을 한 분자에서 직접 연결"한다고 했었는데, 실제로는 **한 샘플 내 hap1 vs hap2 메틸화 비교**(각인 스타일 ASM, DSS 기반)이며 특정 변이와 직접 연결되진 않음 — 그래도 안 써본 haplotype 기반 메틸화 지표이므로 진행 가치 있음.

**방법:**
1. `analyze_asm` rule을 12샘플 전체에 대해 실행 (이미 파이프라인에 있음, `scripts/asm_analysis.R`, 입력 `cpg_hap1_bed`/`cpg_hap2_bed` 12샘플 전부 확인됨).
2. 출력(`{sample}.asm_results.csv`)에서 샘플별 "유의 ASM 영역 수", "평균 |메틸화차이|" 등을 요약 지표로 추출.
3. §2 방식대로 군간 Kruskal-Wallis + 기존 이질성 지표(hom-peak%, discordance_rate, 카피수 안정성)와 순위 일치 여부 확인.

**입력:** 12샘플의 `cpg_hap1_bed`/`cpg_hap2_bed`
**출력:** `h2o2_analysis_results/asm/{sample}.asm_results.csv`, `{sample}.asm_plots.pdf` (파이프라인 표준 출력), 요약: `h2o2_analysis_results/cohort/asm_summary.tsv`, `asm_kruskal_wallis.csv`

### 결과 (완료, 2026-07-22)

샘플별 ASM 후보 영역 수(`n_asm_regions`)의 군별 평균:

| 군 | n_asm_regions 평균 |
|---|---|
| A | 97,478 |
| C | 61,978 |
| B | 41,780 |
| D | 20,955 |

**A>C>B>D — 기존 이질성 그라디언트와 정확히 일치.** Kruskal-Wallis H=9.9744(df=3, 완전분리). 다만 이 지표는 "A/B haplotype을 구분할 수 있는(=이형접합) 지점이 많을수록" 후보 영역도 자연히 많아지는 구조라, §2의 `n_het_snv` 감소와 부분적으로 같은 원인(이형접합 밀도)을 공유하는 지표임을 유의 — 완전히 독립적인 신호는 아니고 "메틸화 쪽에서도 같은 패턴이 보인다"는 보너스 확인 정도로 취급.

---

## 3단계: 분자(read) 단위 haplotype 일관성 확인

**목적:** 10-20kb 리드 하나가 여러 이형접합 SNP를 동시에 스팬할 때, 그 SNP들의 조합이 phasing 알고리즘이 배정한 2-haplotype 모델과 일관되는지 확인. 진짜 서브클론이 2개보다 많이 섞여 있다면 "unphased"/모호 배정 비율이 올라가거나, 국소적으로 haplotype 배정이 불안정해지는 패턴이 나타날 수 있다.

**방법(1차, 간단한 프록시부터):**
1. 12개 BAM(`merged_haplotagged_bam`) 인덱싱(`samtools index`) — 사전 필요.
2. 이형접합 SNP가 밀집된 영역(큰 PS 블록과 겹치는 곳)에서 read별 `HP` 태그 분포 확인: HP=1, HP=2, HP 없음(미배정)의 비율을 샘플별로 집계.
3. "미배정 read 비율"이 높을수록 그 영역의 haplotype 구조가 단순 2-클론 모델에 안 맞는다는 신호로 해석 — 샘플별 전체 미배정 비율을 지표화해 군간 비교.

**(2차, 여유 있으면) 더 엄밀한 버전:** pysam으로 특정 큰 PS 블록 내 2개 이상 SNP를 동시에 커버하는 read들을 뽑아, 각 read의 다중 SNP allele 조합이 정확히 2개의 haplotype 패턴으로만 분리되는지, 아니면 "제3의 조합"이 유의한 빈도로 나타나는지 확인(더 직접적인 서브클론 증거지만 구현 복잡도 높음 — 1차 결과 보고 진행 여부 결정).

**입력:** 12샘플의 `merged_haplotagged_bam`(상염색체 대표 8구간, 5Mb씩 샘플링 — 전체 게놈 스캔은 40GB×12개 기준 비현실적이라 축소)
**출력:** `scripts/read_level_haplotype_consistency.py`, `h2o2_analysis_results/cohort/hp_consistency_combined.tsv`, `hp_consistency_kruskal_wallis.csv`, 시각화: `scripts/plot_hifi_native_results.R` → `h2o2_analysis_results/cohort/hp_consistency_plots.pdf` (n_het_snv 상관 산점도로 재해석 근거 시각화, Spearman rho=-0.993)

### 결과 (완료, 2026-07-22) — 방향은 일치하나 해석 정정 필요

12개 BAM 전체를 `samtools index`(신버전 1.23.1 사용, 구버전 1.3.1은 스레딩 미지원)로 인덱싱 후, 대표 구간 8곳(각 5Mb)에서 HP 태그 미배정(unassigned) read 비율을 계산:

| 군 | unassigned_frac 평균 |
|---|---|
| A | 0.489 |
| C | 0.693 |
| B | 0.814 |
| D | 0.909 |

Kruskal-Wallis H=10.38(df=3). **A<C<B<D로 방향은 기존 그라디언트와 일치**하지만, 원래 세운 가설("서브클론이 여러 개면 phasing이 헷갈려 미배정↑")과는 **반대 방향**이다 — 가설대로면 가장 이질적인 A가 미배정률이 가장 높아야 하는데 실제로는 A가 가장 낮다.

**정정된 해석**: phasing은 read가 이형접합 사이트를 걸쳐야만 HP1/HP2를 배정할 수 있다. D는 이미 §2에서 확인했듯 이형접합 사이트 자체가 A의 1/3 수준(`n_het_snv`)이라, phasing에 쓸 "닻"이 애초에 부족해 미배정 read가 많아진 것 — 다클론 혼란이 아니라 이형접합 밀도 감소라는 **이미 아는 원인의 재관측**이다. 즉 이 지표는 §2단계(ASM 후보 영역 수)와 마찬가지로 독립적인 새 증거라기보다, 같은 근본 원인(이형접합 밀도)을 다른 측정 방식으로 다시 포착한 것으로 문서화한다.

2차(더 정밀한 read-level 다중 SNP 조합 분석)는 이 재해석에 따라 **진행하지 않음** — 1차 결과가 애초 가설(서브클론 혼란)을 지지하지 않아, 더 정교한 버전을 만들어도 같은 함정(이형접합 밀도 confound)에 빠질 가능성이 높다고 판단.

---

## 실행 순서 및 확인 지점

1단계 → 결과 보고 → 2단계 → 결과 보고 → 3단계(1차 프록시) → 결과 보고 → (필요시) 3단계 2차 정밀 버전. 각 단계 결과가 기존 이질성 그라디언트(A>C>B>D)와 같은 방향인지가 핵심 확인 포인트. 예상 밖 결과나 데이터 문제 발견 시 즉시 보고.

## 종합 결론 (2026-07-22 완료)

3단계 모두 A<C<B<D(또는 역방향 A>C>B>D) 그라디언트와 방향이 일치했다. 그중 **1단계(haplotype 방향성 allelic imbalance)만 원리상 완전히 독립적인 새 증거**이고, 2단계(ASM 후보 수)와 3단계(HP 미배정률)는 결과적으로 §2의 `n_het_snv`(이형접합 밀도) 감소라는 이미 아는 원인을 다른 측정 방식으로 재확인한 것에 가깝다 — 그래도 "메틸화 차원", "read 단위 phasing 차원"에서도 같은 방향이 일관되게 나온다는 점 자체는 그라디언트의 신뢰도를 보강한다.
