# H2O2 성상세포 WGS 연구 — 발표 자료용 소스 문서

**작성일:** 2026-07-25
**용도:** WGS 생물정보학을 전혀 모르는 청중에게 배경지식을 전달하면서, 우리 실험 결과(세포집단 이질성 그라디언트 A<C<B<D)를 효과적으로 전달하기 위한 슬라이드별 콘텐츠 소스. AI 다이어그램 생성 도구가 이 문서를 근거로 모식도를 그릴 수 있도록, 각 단계의 입력/출력/파일구조를 구체적으로 기술했다.

**전체 구성:**
- **Part 1**: HiFi 롱리드 WGS 분석 파이프라인 (unaligned BAM → 변이 탐색까지)
- **Part 2**: Small variant/SV 결과, VCF 포맷, 시각화 파일 목록
- **Part 3**: VAF 기반 세포 계통(클론) 분석, 생물학적 검증 계획, 용어집

---

# PART 1. HiFi 롱리드 WGS 분석 파이프라인 (unaligned BAM → 변이 탐색)

## Slide 1-1. HiFi 롱리드 시퀀싱이란 (배경)

**핵심 메시지**: 우리는 PacBio Revio 장비로 "HiFi(High-Fidelity) 롱리드"를 생성했다. 이는 기존 short-read(일루미나, 150bp 내외)와 대비되는 개념이다.

- **읽는 길이(read length)**: 약 15,000~20,000 bp (기존 short-read의 100배 이상)
- **정확도**: 약 99.9% (Q30 이상) — 길면서도 정확한 것이 HiFi의 핵심 차별점 (PacBio의 CCS: Circular Consensus Sequencing 기술로 같은 분자를 여러 번 반복 읽어 오차를 상쇄)
- **롱리드의 장점**: 하나의 read가 여러 개의 변이 지점(SNP)을 동시에 커버 → "이 두 변이가 같은 염색체(같은 부모 유래)에 있는지"를 직접 알 수 있음 (→ Part 1-6 phasing에서 활용)
- **모식도 제안**: short-read(짧은 막대 다수, 조각난 퍼즐) vs HiFi long-read(긴 막대, 하나로 이어진 퍼즐 조각) 비교 그림

## Slide 1-2. 전체 파이프라인 개요 (플로우차트)

**모식도 제안**: 아래 순서를 세로 플로우차트로. 각 박스 = 1개 분석 단계, 화살표 = 데이터 흐름, 박스 옆에 사용 도구명 표기.

```
[1] Unaligned BAM (원시 HiFi 리드)
        │  QC 리포트 생성
        ▼
[2] 참조 유전체 정렬 (pbmm2)  →  Aligned BAM
        │
        ├──▶ [3] Small variant 호출 (DeepVariant)  →  gVCF → VCF
        │         │
        │         └──▶ [3.5] Joint genotyping (GLnexus, 웰당 3반복 합동)  →  Joint VCF
        │
        └──▶ [4] 구조변이(SV) 호출 (Sawfish)  →  SV VCF + 카피수 신호
        │         └──▶ [3.5] Joint SV 호출 (Sawfish joint, 웰당 3반복)
        │
        ▼
[5] Phasing/haplotype 재구성 (HiPhase)  →  Phased VCF + Haplotagged BAM
        │
        ├──▶ [6] 메틸화(5mC) 호출 (pb-CpG-tools)  →  CpG bed/bigwig (hap1/hap2별)
        ├──▶ [7] 반복서열(STR) 유전형분석 (TRGT)  →  Repeat VCF
        └──▶ [8] 미토콘드리아 변이 분석 (MitoRSaw)  →  Mito VCF
```

이 파이프라인은 PacBio 공식 `HiFi-human-WGS-WDL` 워크플로우를 기반으로 하며, 본 프로젝트에서는 인간(GRCh38) 대신 마우스 참조유전체(GRCm39)로 실행하도록 구성을 변경해 사용했다.

## Slide 1-3. 1단계 — Unaligned BAM (원시 데이터)

| 항목 | 내용 |
|---|---|
| 입력 | PacBio Revio 시퀀서 raw output (분자 단위 신호) |
| 출력 | `{sample}.unaligned.bam` — 정렬 전 HiFi read 서열 + 품질값 + kinetics 태그(메틸화 호출에 재사용됨) |
| 파일 구조 | 표준 BAM(Binary Alignment Map) 포맷이지만 아직 참조 위치 정보 없음. 각 레코드=1개 read, 서열(SEQ)+품질(QUAL)+PacBio 전용 태그(예: `fi`/`fp`/`ri`/`rp` — kinetics 정보) 포함 |
| QC 산출물 | `read_length_plot`(read 길이 분포), `read_quality_plot`(품질 분포) — 시퀀싱이 잘 됐는지 1차 확인 |

## Slide 1-4. 2단계 — 참조 유전체 정렬 (Alignment)

| 항목 | 내용 |
|---|---|
| 도구 | **pbmm2** (PacBio가 minimap2를 HiFi 데이터에 맞게 wrapping한 정렬 도구) |
| 목적 | 각 read가 게놈의 어느 위치에서 유래했는지 결정 |
| 입력 | unaligned BAM + 참조 유전체(FASTA, 본 프로젝트는 GRCm39 마우스) |
| 출력 | 정렬된 BAM (좌표 정렬, 인덱스 `.bai` 포함) |
| QC 산출물 | `mapq_distribution_plot`(정렬 신뢰도 MAPQ 분포), `mosdepth_summary`/`mosdepth_region_bed`/`mosdepth_depth_distribution_plot`(커버리지 깊이 — 우리 데이터는 통상 30배 내외) |
| **주의 (참조 편향, reference bias)** | 정렬 시 "완벽히 일치하는 read만" 선택되는 게 아니라 read 전체 길이 기준 유사도로 매칭되므로, 1개 염기가 다른(ALT 대립유전자를 가진) read도 정상적으로 정렬됨 — 이 점이 Part 3의 VAF 해석에서 중요 |

## Slide 1-5. 3단계 — Small Variant 호출 (SNV/Indel)

| 항목 | 내용 |
|---|---|
| 도구 | **DeepVariant** (Google/PacBio, CNN 기반 딥러닝 변이 호출기) |
| 목적 | 단일염기변이(SNV)와 소규모 삽입/결실(indel, 보통 <50bp) 탐지 |
| 입력 | 정렬된 BAM + 참조 FASTA |
| 1차 출력 | **gVCF**(genomic VCF) — 변이 위치뿐 아니라 "변이 없음이 확인된 구간"도 END= 태그로 블록화해 기록하는 포맷 (실제 예시는 Slide 2-1 참고) |
| 최종 출력 | VCF (변이만 남긴 표준 포맷 — Part 2에서 상세) |
| QC 산출물 | `small_variant_stats`(Ti/Tv ratio 등 통계), `snv_distribution_plot`/`indel_distribution_plot`(변이 길이/타입 분포) |

## Slide 1-5b. 3.5단계 — Joint Variant Calling (그룹 내 3반복 합동 유전형)

**핵심 메시지**: 위 3단계까지는 "샘플 1개씩" 변이를 호출한 것. 하지만 우리 실험은 각 조건(웰)마다 3개 반복 샘플이 있으므로, **같은 조건의 3개 반복을 한꺼번에 합쳐 유전형을 다시 판정(joint genotyping)** 한다. 이렇게 하면 한 샘플에서만 애매하게 호출된 위치도 3개를 함께 보고 일관되게 판정할 수 있고, "같은 조건인데 반복 웰마다 유전형이 얼마나 다른가"(=집단 이질성)를 정량화할 수 있다. → Part 3 방법 2(불일치율)의 직접 입력.

| 항목 | 내용 |
|---|---|
| 도구 (small variant) | **GLnexus** (v1.4.3) — 여러 샘플의 gVCF를 병합·재유전형(joint genotyping)하는 표준 도구 |
| 도구 (SV) | **Sawfish (joint 모드)** — 3반복을 함께 SV 호출 |
| 입력 | 같은 웰 3개 반복 샘플의 gVCF (small) / 정렬 BAM (SV) |
| 구성 단위 | 웰(그룹)별로 별도 실행 → `joint-H2O2-A01`, `joint-H2O2-B01`, `joint-H2O2-C01`, `joint-H2O2-D01` 4개 |
| 출력 (small) | 합동 VCF를 다시 샘플별로 분리한 `split_joint_small_variant_vcfs/{0,1,2}/` (반복 3개) |
| 출력 (SV) | `split_joint_structural_variant_vcfs/`, 그리고 `sv_copynum_bedgraph`/`sv_copynum_summary`/`sv_depth_bw`/`sv_maf_bw` 등 카피수·depth 트랙 |
| **개별 호출과의 차이** | 개별 DeepVariant VCF는 `FILTER=PASS`를 세팅하지만, GLnexus joint VCF는 **FILTER를 "."(미평가)로 남기고 `MONOALLELIC`만 별도 표시** → 분석 시 `FILTER=PASS` 대신 "MONOALLELIC만 제외" 규칙을 써야 함(엔지니어링 주의점) |

## Slide 1-6. 4단계 — 구조변이(SV) 호출

| 항목 | 내용 |
|---|---|
| 도구 | **Sawfish** (PacBio, 최신 SV 호출기 — 이전 세대 pbsv를 대체) |
| 목적 | 큰 구조변이(결실/중복/역위/전좌, 보통 ≥50bp) 및 카피수(CN) 변화 탐지 |
| 입력 | 정렬된 BAM |
| 출력 1 | SV VCF (Slide 2-3에서 상세) |
| 출력 2 (본 연구 핵심 활용) | `sv_copynum_bedgraph`(게놈 전역 depth 기반 카피수 구간), `sv_copynum_summary.json` — SV 호출의 부산물이지만 별도 계산 없이 카피수 안정성 분석에 재사용 |
| 기타 산출물 | `sv_depth_bw`, `sv_gc_bias_corrected_depth_bw`(GC 편향 보정 depth), `sv_maf_bw`(minor allele frequency track), `sv_supporting_reads`(SV를 지지하는 read 목록) |

## Slide 1-7. 5단계 — Phasing (haplotype 재구성) ★ HiFi 고유 기능

**핵심 메시지**: 이 단계가 short-read WGS와 가장 차별화되는 지점. 부모로부터 물려받은 두 벌의 염색체(haplotype 1, haplotype 2) 중 어느 변이가 어느 쪽에 있는지 재구성한다.

| 항목 | 내용 |
|---|---|
| 도구 | **HiPhase** (PacBio) |
| 원리 | 긴 read 하나가 여러 이형접합(heterozygous) 변이를 동시에 커버하면, 그 변이들이 "같은 haplotype"에 있다고 직접 연결 가능 (short-read는 read가 짧아 이 연결이 자주 끊김) |
| 입력 | 정렬된 BAM + small variant VCF + SV VCF (함께 phasing) |
| 출력 1 | **Phased VCF** — `FORMAT`에 `PS`(phase set ID) 태그 추가, genotype이 `0/1`(순서無) 대신 `0\|1`/`1\|0`(순서有, `\|`가 phasing 표시) 형태로 기록됨 |
| 출력 2 | **Haplotagged BAM**(`merged_haplotagged_bam`) — 각 read에 `HP:i:1`/`HP:i:2`(어느 haplotype 유래인지) 태그 부여. 이 정보가 Part 3의 haplotype 방향성 분석과 메틸화 haplotype 분리에 재사용됨 |
| QC 산출물 | `phase_blocks`(연속적으로 phasing된 구간 목록), `phase_haplotags`, `phase_stats`(블록 크기 통계) |
| 검증 확인 사항 | phased VCF 헤더에 `MF`/`MD`/`MT`(변이별 메틸화 비율/깊이/타입) FORMAT 필드가 정의돼 있으나, 실제로는 전량 비어있음(`.`) — 이 경로는 사용 불가로 확인 후 폐기, 대신 hap1/hap2 BED 비교 방식으로 메틸화 분석 진행(Slide 1-8) |

## Slide 1-8. 6단계 — 메틸화(5mC) 호출

| 항목 | 내용 |
|---|---|
| 도구 | **pb-CpG-tools** (PacBio) — HiFi kinetics 신호에서 직접 5mC(CpG 메틸화) 추정 (별도의 bisulfite 처리 불필요, HiFi의 또 다른 고유 장점) |
| 입력 | Haplotagged BAM (kinetics 정보 + HP 태그) |
| 출력 | `cpg_combined_bed`(전체), `cpg_hap1_bed`/`cpg_hap2_bed`(haplotype별 분리 — allele-specific methylation 분석용), 각각의 bigwig(`cpg_combined_bw`, `cpg_hap1_bw`) — genome browser 시각화용 |
| 후속 요약 도구 | **MethBat**(`methbat_profile`) — 메틸화 프로파일 요약 |

## Slide 1-9. 7~8단계 — 반복서열/미토콘드리아 (부가 분석)

| 항목 | 도구 | 출력 |
|---|---|---|
| 반복서열(STR) 유전형분석 | **TRGT** (Tandem Repeat Genotyper) | `phased_trgt_vcf` — 반복 모티프 개수/길이를 haplotype별로 산출(Slide 2-2 상세), `trgt_coverage_dropouts`, `trgt_spanning_reads` |
| 미토콘드리아 변이 | **MitoRSaw** | `mitorsaw_vcf`, `mitorsaw_hap_stats` — 미토콘드리아는 모계 유전이라 별도 haplotype 로직 사용 |
| Runs of Homozygosity | **bcftools roh** | `bcftools_roh_bed`, `bcftools_roh_out` — 동형접합 구간 탐지(근친교배/계통 분석 참고 지표) |

## Slide 1-10. Part 1 요약 표 (전체 도구/입출력 한눈에)

| 단계 | 도구 | 입력 | 출력 |
|---|---|---|---|
| 정렬 | pbmm2 | unaligned BAM + 참조 FASTA | aligned BAM |
| Small variant | DeepVariant | aligned BAM | gVCF → VCF |
| Joint genotyping | GLnexus | 웰당 3반복 gVCF | joint VCF (split 3반복) |
| SV/카피수 | Sawfish (개별+joint) | aligned BAM | SV VCF, copynum bedgraph/summary |
| Phasing | HiPhase | aligned BAM + VCF(small+SV) | phased VCF, haplotagged BAM |
| 메틸화 | pb-CpG-tools | haplotagged BAM | CpG bed/bigwig (hap1/hap2/combined) |
| 반복서열 | TRGT | phased 정보 | repeat VCF |
| 미토콘드리아 | MitoRSaw | aligned BAM | mito VCF |
| 커버리지 QC | mosdepth | aligned BAM | depth summary/plot |

---

# PART 2. Small Variant / SV 분석 결과, VCF 포맷, 시각화

## Slide 2-1. VCF 포맷 상세 설명 (일반 청중용)

**핵심 메시지**: VCF(Variant Call Format)는 "참조 유전체와 다른 지점"을 표 형태로 정리한 표준 텍스트 포맷.

**실제 데이터 예시** (DeepVariant 출력, 마우스 chr1):
```
#CHROM  POS      ID  REF  ALT    QUAL  FILTER  INFO  FORMAT                  H2O2-A01-ctrl1
chr1    3050118  .   G    A,<*>  31    PASS    .     GT:GQ:DP:AD:VAF:MID:PL  0/1:31:38:25,13,0:0.342105,0:small_model:30,0,51,990,990,990
```

| 컬럼 | 의미 | 위 예시 값 |
|---|---|---|
| CHROM / POS | 염색체 / 위치 | chr1, 3050118번째 염기 |
| REF / ALT | 참조 대립유전자 / 변이 대립유전자 | G(참조) → A(변이) |
| QUAL | 변이 호출 신뢰도(Phred 스케일) | 31 |
| FILTER | 품질 필터 통과 여부 | PASS (통과) |
| FORMAT의 GT | Genotype(유전형) | 0/1 = 이형접합(한쪽만 변이) |
| FORMAT의 DP | 이 위치를 커버한 read 총 개수(depth) | 38개 |
| FORMAT의 AD | 대립유전자별 지지 read 수(Allelic Depth) | REF 25개, ALT 13개 |
| **FORMAT의 VAF** | **Variant Allele Frequency = ALT지지 read / 전체 read** | 13/38 ≈ 0.342 |
| FORMAT의 PL | 유전형별 Phred-scaled likelihood | (낮을수록 그 유전형일 가능성 높음) |

**gVCF와 VCF의 차이**: gVCF는 "변이 없음이 확인된 구간"도 `END=` 태그로 블록 압축해 기록(예: `chr1 3050001 . C <*> 0 . END=3050117 GT:GQ:MIN_DP:PL 0/0:50:34:...`) — "검사는 했고 정상이었다"는 정보까지 남기는 포맷. 최종 분석에 쓰는 VCF는 이런 정상 구간을 제외하고 실제 변이만 남긴 버전.

**Phasing 후 GT 표기 변화**: `0/1`(순서 모름) → `0|1` 또는 `1|0`(순서 앎, `|`가 phasing 표시) + `PS`(phase set ID) 태그 추가 — Part 1-7 참고.

## Slide 2-1b. Joint VCF 포맷 (GLnexus) — 개별 VCF와의 차이

**핵심 메시지**: 같은 위치라도 개별 호출 VCF와 joint 호출 VCF는 필드가 다르다. joint VCF는 3반복을 함께 판정한 결과라 집단 이질성(방법 2)의 근거가 된다.

**실제 데이터 예시** (GLnexus 출력, 같은 chr1:3050118):
```
#CHROM  POS      ID                REF  ALT  QUAL  FILTER  INFO                    FORMAT             H2O2-2-A01
chr1    3050118  chr1_3050118_G_A  G    A    37    .       AF=0.5;AQ=37;AC=1;AN=2  GT:DP:AD:GQ:PL:RNC  0/1:63:46,17:38:37,0,54:...
```

| 항목 | 개별(DeepVariant) VCF | Joint(GLnexus) VCF |
|---|---|---|
| ID 컬럼 | `.` (비어있음) | `chr1_3050118_G_A` (위치·염기 기반 고유 ID 부여) |
| FILTER | `PASS` | `.` (미평가) + `MONOALLELIC`만 별도 표시 |
| INFO | 비어있음 | `AF`(대립빈도), `AQ`(대립품질), `AC`(ALT 개수), `AN`(전체 대립 수) |
| VAF 필드 | `VAF`(직접 제공) | 없음 → `AD`(46,17)에서 17/63으로 직접 계산 |
| FORMAT 태그 | GT:GQ:DP:AD:VAF:MID:PL | GT:DP:AD:GQ:PL:**RNC**(Reason for No Call, 미콜 사유) |
| 도구/버전 | DeepVariant | GLnexus v1.4.3 |

**분석에서의 활용**: 웰당 3반복을 다시 샘플별로 분리한 3개 VCF(`split_joint_small_variant_vcfs/{0,1,2}/`)를 bcftools merge로 재병합해, 같은 위치에서 3반복의 GT가 일치하는지 비교 → **불일치율(discordance_rate)** 산출(Part 3 방법 2). 위상 무시하고 allele을 정렬한 튜플로 정규화해 비교하며, 3반복 모두 유효 GT이고 최소 1개는 non-ref인 사이트만 대상으로 삼는다.

## Slide 2-1c. Per-sample VCF vs Joint Calling — 군간 비교 두 방식

**핵심 메시지**: 두 방식의 진짜 차이는 **"빈 칸(없는 행)의 의미"** — 어떤 샘플에 그 위치의 행이 없을 때, 그게 *참조형(정상)* 인지 *데이터 없음* 인지를 구분할 수 있느냐다.

**핵심 모식도 (한 위치 X, 2샘플 유전형 행렬):**
```
                 Position X
  Sample A   ●  0/1  (변이 호출됨)
  Sample B   ?  ───

  개별 VCF   : (빈 칸) → 참조형? 커버리지 부족? 필터 제외?  → 모호(ambiguous)
  Joint      : 0/0 (DP=63, GQ=38) → 확실히 참조형        → 해소(resolved)
```
gVCF의 참조확신 블록(`END=`/`MIN_DP`)이 있어야 joint가 모든 칸을 채워 **빈 칸 없는 정사각 행렬**을 만든다.

**비교 표:**
| 항목 | 개별(per-sample) VCF | Joint calling (GLnexus) |
|---|---|---|
| "없음(absent)"의 의미 | 모호 (missing인지 참조형인지 불명) | 모든 샘플을 모든 위치에서 유전형 판정(정사각 행렬) |
| 대립유전자 표현 | 샘플마다 다를 수 있음 | 전 샘플 통일 |
| 샘플 추가/제거 | 자유(독립 계산) | 코호트 전체 재실행 필요 |
| 아티팩트 위험 | 샘플별로 격리됨 | 전 샘플로 증폭될 수 있음 |
| 적합한 분석 | 샘플 내부 지표(VAF, het-hom peak) | 샘플 간 유전형 비교(불일치율, 코호트 AF) |

**본 연구 실사례 (강조 박스):** SV burden test —
- 개별: FDR<0.05 유의 유전자 **0개** → 보수적·신뢰
- Joint: **112개** "유의" → 그러나 상위가 전부 다카피 유전자족(Scgb/Sirpb) = 아티팩트
- → **교훈: joint가 항상 낫지 않다.** 여기선 카피수 아티팩트를 증폭했다.

**한 줄 결론**: 질문 성격으로 고른다 — *한 샘플 내부* → 개별 VCF, *샘플 간 유전형 맞대기* → joint calling.

**정직성 주의**: 우리 GLnexus 설정은 `revise_genotypes=false`라, joint의 이득은 "정사각 행렬"이지 "샘플 간 유전형 보정"은 아니었다(과장 금지). 필터도 개별=`FILTER=PASS`, joint=`FILTER="."`+`MONOALLELIC`이라 후처리 규칙이 다름(`require_filter_pass=false`).

## Slide 2-1d. bcftools — VCF를 읽고·거르고·값을 뽑는 표준 도구

**핵심 메시지**: VCF는 변이가 수십만~수백만 줄이라 텍스트 에디터로 못 연다. **bcftools**는 이 VCF/BCF를 효율적으로 다루는 표준 커맨드라인 도구로, 보통 **"VCF에서 필요한 값을 꺼내는 앞단"** 역할을 하고, 꺼낸 값의 통계·시각화는 Python/R로 넘긴다. (BCF = VCF의 압축 이진 버전, `samtools`와 같은 htslib 계열.)

**핵심 개념 — 분석은 두 동작으로 나뉜다:**
```
(a) 행 선택 (filtering)  :  bcftools view -i 'FILTER="PASS"'      → 어떤 변이를 남길까
(b) 값 추출 (extraction) :  bcftools query -f '%CHROM\t[%VAF]\n'  → 어떤 태그 값을 꺼낼까
```

**주요 서브커맨드 (본 연구 사용):**
| 명령 | 하는 일 | 본 연구 용도 |
|---|---|---|
| `bcftools view` | 헤더/레코드 보기, 행 필터링(`-i`/`-e`, `-f PASS`) | VCF 헤더·레코드 확인, PASS 필터 |
| `bcftools query` | 원하는 태그 값만 뽑아 표로 출력 | VAF/GT 추출 → het_hom_peak_frac 계산 |
| `bcftools merge` | 여러 샘플 VCF를 같은 좌표로 합침 | 방법 2 불일치율(3반복 재병합) |
| `bcftools stats` | Ti/Tv 등 기술통계 요약 | QC |
| `bcftools norm` | 다중대립·indel 표기 정규화 | 전처리 |
| `bcftools roh` | Runs of Homozygosity 탐지 | 동형접합 구간(Slide 1-9) |

**흐름 모식도:**
```
VCF/BCF ──bcftools(필터+추출)──▶ 값 테이블 ──Python/R──▶ 파생·통계·시각화
        (앞단: 값 꺼내기)                    (뒷단: 분석)
```

**한 줄 결론**: bcftools는 VCF에서 값을 꺼내는 앞단 도구다 — 무거운 변이 호출(DeepVariant/Sawfish)이 끝난 뒤, 이 값 추출과 후속 통계는 초~분 단위로 가볍게 반복할 수 있다(Tier-2 분석).

## Slide 2-2. 반복서열(STR) VCF 포맷 — TRGT 출력

일반 SNV/indel VCF와 다른 전용 필드 사용 (실제 헤더 확인):

| 필드 | 의미 |
|---|---|
| INFO: TRID | 반복서열 고유 ID |
| INFO: MOTIFS | 이 반복 영역을 구성하는 반복 모티프 서열 |
| FORMAT: AL | 대립유전자별 반복 길이 |
| FORMAT: SD | 대립유전자별 지지 read 수 |
| FORMAT: MC | 대립유전자별 모티프 반복 횟수 |
| FORMAT: AM | 대립유전자별 평균 메틸화 수준 |
| FORMAT: PS | phase set ID (small variant와 동일 개념) |

## Slide 2-3. SV(구조변이) VCF와 결과 파일

| 파일 | 내용 |
|---|---|
| `phased_sv_vcf` | 결실/중복/역위 등 구조변이 목록, haplotype 배정 포함 |
| `sv_copynum_bedgraph` | 게놈 전역을 일정 구간으로 나눠 각 구간의 추정 카피수(depth 기반) 기록 — bedgraph는 `염색체 시작 끝 값` 4열 텍스트 포맷 |
| `sv_copynum_summary.json` | 샘플 전체 카피수 통계 요약 |

**주요 발견 & 수정**: 단일 9.8Mb 크기의 거대 중복(duplication) 변이 하나가 582개 유전자 중 490개를 "유의미하게 다르다"고 오염시킨 아티팩트를 발견 → `--max-svlen` 옵션으로 초대형 SV를 제외하도록 파이프라인 수정. 또한 A/D 웰=수컷(XY), B/C 웰=암컷(XX)이라는 성염색체 교란 요인을 mosdepth depth ratio로 확인 후 `--exclude-chroms chrX,chrY`로 제외.

## Slide 2-4. Small Variant 결과 파일과 분석 도구 재정리

| 파일/분석 | 사용 도구 | 목적 |
|---|---|---|
| `phased_small_variant_vcf` | DeepVariant + HiPhase | 최종 SNV/indel 목록(phasing 포함) |
| VEP 유전자 영향 주석 | Ensembl **VEP** | 각 변이가 어떤 유전자/기능에 영향(missense, nonsense 등)을 주는지 주석 |
| 유전자별 변이 부담(burden) 검정 | `scripts/variant_burden.R` | 그룹 간 특정 유전자에 변이가 쏠렸는지 통계 검정 (Fisher's exact + 효과크기: odds ratio, Cramér's V) |
| GO 유전자 세트(geneset) 부담 검정 | `scripts/geneset_burden.R` + `fetch_go_genesets.py` | 개별 유전자가 아닌 기능적 유전자 그룹(GO term) 단위로 집계 |
| DMR(차등 메틸화 영역) 분석 | `scripts/dmr_analysis.R` (methylKit/DSS 계열) | hap1/hap2 또는 그룹 간 메틸화 차이가 큰 영역 탐색 |

**[슬라이드 하단 주석 가이드]** — 이 슬라이드(및 Slide 2-5, 2-5b) 하단에 작은 footnote로 "burden" 용어 설명을 넣어줘. 청중이 burden testing을 처음 들을 때 뜻을 바로 알 수 있도록. 아래 문구를 각주 스타일(작은 글씨)로:

> *Burden = 한 유전자(또는 유전자 세트)가 짊어진 변이의 총량(있음/없음 또는 개수). 개별 변이는 너무 드물어 검정력이 없으므로, 유전자 단위로 변이를 뭉쳐(collapse) "그 총량이 군마다 다른지"를 검정하는 방식. 예) 유전자별 n_A~n_D = 각 군 3샘플 중 그 유전자에 변이를 가진 수.*

(영어 덱이면: *Burden = the total load of variants carried by a gene (or gene set) — presence/absence or count. Because individual rare variants have too little power, variants are collapsed per gene and tested for group differences. e.g., n_A–n_D = how many of each group's 3 samples carry a variant in that gene.*)

## Slide 2-5. SV 결과 파일과 분석 도구 재정리

| 파일/분석 | 사용 도구 | 목적 |
|---|---|---|
| `phased_sv_vcf` | Sawfish | SV 목록 |
| SV 기능 영향 주석 | **svpack** (consequence 서브커맨드, gff3 참조) | SV가 어떤 유전자에 걸치는지 주석 |
| 유전자별 SV 부담 검정 | `scripts/sv_burden.R` | 그룹 간 특정 유전자에 SV가 쏠렸는지 검정 |
| 카피수 안정성 검정 | `scripts/copynum_stability_check.py` | 알려진 다카피 유전자 자리(7곳) 기준 카피수 변동 안정성 비교 |

## Slide 2-5b. Burden Test 결과 요약 — "개별 유전자 범인 찾기"의 한계 ★

**핵심 메시지**: 각 변이 유형(SNV/indel/SV)별로 "특정 유전자에 변이가 군별로 쏠렸는가"를 검정(Fisher's exact)했다. **결론: 다중검정(FDR) 보정을 하면 개별 유전자 수준의 유의한 신호는 사실상 전멸한다.** 이것이 우리가 "소수 원인 유전자 찾기"가 아니라 "게놈 전역 판세(Part 3)"로 방향을 튼 **정량적 근거**다.

### 검정별 결과 (실제 수치)

| 검정 | 도구/파일 | 검정 유전자 수 | p<0.05 (보정 전) | **FDR<0.05 (보정 후)** | 최소 FDR | 최상위 유전자 |
|---|---|---|---|---|---|---|
| Small variant (SNV+indel) burden | `variant_burden.R` → `variant_burden_stats.csv` | 80 | **1개** (B4galt3, p=0.039) | **0개** | 1.00 | B4galt3 |
| SV burden (단일샘플) | `sv_burden.R` → `sv_burden_stats.csv` | 278 | **1개** (Cyp2c29, p=0.046) | **0개** | 1.00 | Cyp2c29 |
| SV burden (**joint**) | `sv_burden.R` (joint) → `..._joint/.../sv_burden_stats.csv` | 277 | **112개** | **112개** (겉보기 유의!) | 0.024 | Scgb1b10, Sirpb1b, Scgb2b11 |
| GO geneset burden (SNV+indel, 기존) | `geneset_burden.R` → `geneset_burden_variant.csv` | 2 세트 | 계산 불가(NA) | NA | NA | — |
| GO geneset burden (SV, 기존) | `geneset_burden.R` → `geneset_burden_sv.csv` | 2 세트 | 0개 (최소 p=0.108) | 0개 | 0.216 | dna_repair |
| **GO geneset 확대 (SNV+indel, 큐레이션)** | `--preset curated` → `geneset_burden_variant_curated.csv` | **20 세트** | 0개 | **0개** | 0.39 | apoptosis(p=0.31) |
| **GO geneset 확대 (SV, 큐레이션)** | `--preset curated` → `geneset_burden_sv_curated.csv` | **20 세트** | 1개(mitochondrion p=0.036) | **0개** | 0.33 | mitochondrion→FDR 0.33 소멸 |

### 이 표를 어떻게 읽는가 (발표 포인트)

1. **보정 전 → 후의 극적 소멸**: SNV+indel은 80개 중 딱 1개(B4galt3)가 명목상 유의(p=0.039)했지만, 80개를 동시에 검정한 다중검정 보정(BH-FDR)을 하면 **FDR=1.00으로 유의성이 완전히 사라진다**. SV 단일샘플도 동일(Cyp2c29 p=0.046 → FDR=1.00). → "우연히 하나쯤 p<0.05 나오는 것"과 구분되지 않는다.
2. **joint SV의 112개 "유의"는 함정(아티팩트)**: joint SV burden만 112개가 FDR<0.05로 대량 유의하게 보이지만, 최상위가 전부 **다카피 유전자족**(Scgb=secretoglobin, Sirpb, 이전엔 Vmn1r/Mup/Pramel/Ear/Or)이다. 제외 리스트를 늘려도 새 유전자족이 계속 튀어나오는 "두더지잡기"라, 이 112개는 생물학이 아니라 **다카피 영역의 정렬/카피수 아티팩트**로 판단하고 **단일샘플 SV burden(0개 유의)을 신뢰**하기로 결정.
3. **GO geneset도 유의 없음 — 2개→20개로 넓혀도 동일**: 처음엔 산화스트레스·DNA repair 2개 세트만 검정했고 유의하지 않았다(최소 FDR 0.22). 이후 **큐레이션 20개 세트로 확대**(글루타티온·DNA수선 하위경로 BER/NER/MMR/HR/NHEJ·손상체크포인트·세포자멸·페롭토시스·자가포식·노화·ER스트레스·미토콘드리아·염증·신경교분화·저산소·칼슘신호)해도 **두 모드 모두 FDR<0.05 유의 경로 0개**(SNV+indel 최소 FDR 0.39, SV 최소 FDR 0.33). SV의 mitochondrion만 보정 전 p=0.036이었으나 20세트 보정 후 FDR 0.33으로 소멸(우연 기대치 안). 즉 "H2O2니까 산화스트레스/DNA repair 등에 변이가 쏠렸을 것"이라는 가설은 **표적·확대 검정 모두에서** 지지되지 않는다.
4. **그래서 관점 전환**: 개별 유전자·유전자세트 어느 층위에서도 "군을 가르는 범인 유전자"는 안 나온다. 대신 **게놈 전역의 이질성/클론성 지표(Part 3)** 에서 A<C<B<D 그라디언트가 일관되게 나온다 — 이게 이 연구의 핵심 발견이 "범인 유전자"가 아니라 "판세(HRD 유사 게놈 불안정성 총량)"인 이유다.

### 결과 파일

| 파일 | 내용 |
|---|---|
| `variant_burden_stats.csv` | SNV+indel 유전자별: gene, pvalue, odds_ratio, max_prop_diff, group_max/min, odds_ratio_extreme, cramers_v, n_A~n_D, **fdr** |
| `sv_burden_stats.csv` (단일샘플) | SV 유전자별 (동일 컬럼 구조) |
| `h2o2_analysis_results_joint/cohort/sv_burden_stats.csv` | joint SV 유전자별 (다카피 아티팩트 확인용, 단일샘플과 대조) |
| `geneset_burden_variant.csv` / `geneset_burden_sv.csv` | GO 세트별(기존 2세트): geneset, n_genes, pvalue_presence, pvalue_count, max_prop_diff, cramers_v, n_presence_A~D, mean_count_A~D, fdr_presence, fdr_count |
| `geneset_burden_variant_curated.csv` / `geneset_burden_sv_curated.csv` | **큐레이션 20세트 확대판** (동일 컬럼, `--preset curated`) — 넓혀도 유의 0개 확인용 |

### 그래프

PDF에서 SVG(개별 그림 파일)로 전환됨 — 모든 그림이 `h2o2_analysis_results/cohort/cohort_report.html`(파라미터화된 HTML 리포트, 표+그림 통합)에도 함께 모여 있으니 슬라이드 제작 시 이쪽을 먼저 열어보는 걸 권장.

| 디렉토리 | 파일 | 내용 |
|---|---|---|
| `variant_burden_svg/` | `01_top30_genes.svg`, `02_burden_boxplot.svg`, `03_effect_size_vs_significance.svg`, `04_top30_group_bubble.svg`, `05_variant_heatmap.svg` | SNV+indel 유전자별 burden (p-value/효과크기 분포) |
| `sv_burden_svg/` | `01_svtype_distribution.svg`, `02_svlen_distribution.svg`, `03_top20_genes.svg`, `04_effect_size_vs_significance.svg`, `05_top20_group_bubble.svg` | 단일샘플 SV burden |
| `h2o2_analysis_results_joint/cohort/sv_burden_svg/` | (동일 파일명 구조) | joint SV burden (다카피 아티팩트가 상위를 차지하는 모습) |
| `geneset_burden_variant_svg/` / `geneset_burden_sv_svg/` | `01_<geneset명>.svg`, `02_<geneset명>.svg` (세트별 1장씩, 기존 2세트) | GO 세트 단위 burden |
| `geneset_burden_variant_curated_svg/` / `geneset_burden_sv_curated_svg/` | `01_<geneset명>.svg` … `20_<geneset명>.svg` (큐레이션 20세트) | 큐레이션 20세트 확대판 |

## Slide 2-5c. (보조) 다카피 정렬 아티팩트란 — joint SV "112개"의 정체

**핵심 메시지**: joint SV burden의 112개 "유의"는 생물학이 아니라 **다카피 유전자족의 정렬 아티팩트**다. 왜 생기고 어떻게 알아챘는지.

**메커니즘 (모식도용 흐름):**
```
서열 거의 동일한 사본들(Or/Vmn1r/Mup/Scgb…)
   → read가 사본끼리 오매핑(multi-mapping), MAPQ↓
   → 참조 카피수 ≠ 샘플 카피수 → depth 왜곡
   → Sawfish가 가짜 CN/SV 호출
   → n=3에서 확률적 변동이 "가짜 군차이"로
   → joint calling이 전 샘플에 걸쳐 증폭 → 112개
```

**문제의 유전자족:**
| 족 | 정체 | 왜 취약 |
|---|---|---|
| Or | 후각 수용체 | 게놈 최대 유전자족(1000개+) |
| Vmn1r | 서골비 수용체 | 수백 사본, 카피수 다형성 |
| Mup / Scgb / Sirpb / Ear / Pramel | urinary protein·secretoglobin 등 | 최근 중복, 개체마다 카피수 상이 |

**아티팩트로 판단한 5가지 신호:**
1. **생물학적 부적절성** — 최상위가 후각수용체·secretoglobin, H2O2/성상세포와 무관
2. **유전자족 군집** — 게놈에 퍼지지 않고 몇 개 족에 몰림
3. **두더지잡기** — 한 족 제외하면 다른 족이 튀어나옴(체계적 아티팩트의 서명)
4. **단일샘플 vs joint 불일치** — 단일 0개 / joint 112개(방법 의존 = 아티팩트)
5. **문헌적 악명** — 매핑·카피수 추정 어렵기로 알려진 영역

**대응**: 제외 리스트 무한 확장 대신 **단일샘플 SV burden(0개)을 신뢰**. 롱리드가 다중매핑을 크게 줄이지만 거대 탠덤 어레이는 완전 해결 못 함 → 방법 무관하게 조심.

**한 줄 결론**: "유의해 보이는 112개"의 함정 사례 — p값이 아니라 **생물학적 타당성 + 방법 간 일관성**으로 진짜/가짜를 가려야 한다.

## Slide 2-6. 시각화 파일 전체 목록

**PDF → SVG로 전환됨** (한글은 `svglite`+`Noto Sans CJK KR`로 렌더링, 다중패널은 patchwork로 구성은 동일). 모든 SVG는 `h2o2_analysis_results/cohort/` 하위 디렉토리(별도 표기 제외)에 있고, 통계 표까지 한 번에 보려면 같은 위치의 **`cohort_report.html`**(파라미터화된 R Markdown 리포트 — 아래 전체 목록과 동일 섹션 구성의 탭/표 포함)을 여는 게 가장 빠르다. 슬라이드에 낱장으로 넣을 그림만 골라 쓸 땐 아래 개별 SVG 경로를 참고.

| 디렉토리/파일 | 내용 | 대응 분석 |
|---|---|---|
| `variant_burden_svg/01_top30_genes.svg` 등 5장 | SNV+indel 유전자별 burden (p/효과크기 분포) | Slide 2-4/2-5b |
| `sv_burden_svg/03_top20_genes.svg` 등 5장 | 단일샘플 SV burden | Slide 2-5/2-5b |
| `../../h2o2_analysis_results_joint/cohort/sv_burden_svg/` | **joint SV burden** (다카피 아티팩트 상위 차지 확인용) | Slide 2-5b |
| `geneset_burden_variant_svg/` / `geneset_burden_sv_svg/` | GO 세트(기존 2세트) 단위 burden, 세트당 1장 | Slide 2-4/2-5b |
| `geneset_burden_variant_curated_svg/` / `geneset_burden_sv_curated_svg/` | GO 세트 큐레이션 20세트 확대판, 세트당 1장(01~20) | Slide 2-5b |
| `dmr_A_vs_B_svg/` … `dmr_C_vs_D_svg/` (6쌍 디렉토리, 각 `01_dmr_by_chrom.svg`/`02_methylation_diff_histogram.svg`) | 군쌍별 차등 메틸화 영역(DMR) 결과 | Slide 2-4 |
| `trgt_group_compare_svg/01_manhattan.svg`, `02_top20_boxplot.svg` | 반복서열(STR) 군간 비교 | Slide 2-2 |
| `asm/{sample}.asm_plots_svg/` (12샘플, 각 `01_asm_by_chrom.svg` 등) | 샘플별 hap1 vs hap2 ASM 결과 | Part 3 방법 5 |
| `summary_compare_svg/01_metrics_by_group.svg` | 샘플단위 4개 지표(변이수/SV수/hom-peak%/VAF중앙값) 군간 boxplot | Part 3 |
| `summary_compare_svg/02_discordance_rank.svg` | 웰간 불일치율 barplot + 지표별 순위 일관성 heatmap | Part 3 방법 2 |
| `hifi_native_svg/01_hap_imbalance.svg` | Haplotype 방향성 allelic imbalance 결과 | Part 3 방법 4 |
| `hifi_native_svg/02_hp_consistency.svg` | Read 단위 haplotype 배정 일관성 + n_het_snv 상관 산점도 | Part 3 (참고지표) |
| `heterogeneity_consensus_svg/01_heterogeneity_consensus_heatmap.svg` | **5개 독립 방법 종합 heatmap (핵심 결론 그림)** | Part 3 |

**fig-atlas 번들** (논문/발표용 핵심 그림은 png+pdf+svg 3종 + 메타데이터까지 갖춘 별도 번들로도 제공): `cohort/fig_atlas_bundles/{snv_indel_burden_top30, sv_burden_single_top20, sv_burden_joint_top_hits, dmr_direction_by_chrom, trgt_manhattan, heterogeneity_consensus_heatmap}/outputs/figure.{png,pdf,svg}` — 슬라이드에 고해상도 PNG나 벡터 PDF가 필요하면 여기서 가져다 쓰면 된다.

---

# PART 3. VAF 기반 세포 계통(클론) 분석

## Slide 3-1. 실험 설계 및 핵심 질문

- **모델**: 마우스 성상세포(astrocyte), H2O2 산화스트레스 처치
- **4개 그룹** (각 n=3 반복 웰):
  - **A** = 대조군 (control)
  - **B** = 급성 H2O2 처치
  - **C** = 급성 H2O2 + 약물
  - **D** = 만성 H2O2 + 약물
- **핵심 질문**: 처치 강도가 강할수록 세포집단 내 특정 계통(clone)이 수적으로 우세해지는(균질화되는) 경향이 있는가?
- **핵심 발견 (요약)**: 5개의 독립적인 분석 방법 모두에서 **A < C < B < D** 순서로 이질성이 감소(=균질화/클론 우세 증가)하는 일관된 그라디언트 확인 (Kruskal-Wallis H≈10, n=3×4군 완전 분리)

## Slide 3-2. VAF란 무엇인가 (개념 설명, 비유 포함)

- **VAF(Variant Allele Frequency)** = 특정 위치에서 "변이 있는 read 수 / 전체 read 수"
- **단일세포라면**: 이형접합(부모 중 한쪽에서만 변이 물려받음)일 때 이론적으로 VAF ≈ 0.5 (두 벌의 염색체 중 하나만 변이)
- **하지만 우리 WGS는 "bulk"(수많은 세포를 한꺼번에 갈아서) 시퀀싱** → 측정되는 VAF는 그 웰 안 모든 세포의 평균값
- **비유**: 빨간 구슬(변이 있음)과 파란 구슬(변이 없음)이 섞인 통에서 한 움큼 집었을 때 빨간 구슬 비율 — 통 안의 실제 구성 비율을 반영

## Slide 3-3. Bulk VAF → 세포 비율 역추론 (핵심 논리)

**핵심 메시지**: 우리는 세포 비율을 직접 모르지만, bulk VAF 값으로부터 "대략 이런 비율이었을 것"이라 역산할 수 있다.

- 이형접합 변이 = 세포 하나당 염색체 2벌 중 1벌만 변이 보유
- 만약 집단 전체가 유전적으로 균일(모든 세포가 동일 비율로 REF/ALT 염색체 보유)하다면 VAF ≈ 0.5 근처에 몰려야 함
- 실제로 특정 위치에서 VAF가 0.9 이상으로 튀는 경우가 많아짐 = 그 위치에서 "변이를 상실한(또는 우세해진) 특정 계통"이 집단 내 수적으로 지배적이 되었다는 신호
- **역추론 예시**: 특정 유전자좌에서 원래 이형접합이었는데, 어떤 세포 계통에서 우연히 그 자리의 변이 사본이 사라지는 사건(LOH, 아래 용어집 참고)이 일어났다고 가정. 이 계통이 처치에 의해(생존 유리 또는 증식 유리로) 집단의 90%를 차지하게 되면, bulk VAF는 (0.9×~0 + 0.1×0.5) ≈ 0.05 또는 반대 방향으로 쏠려 관측된다 — 즉 **VAF가 0.5에서 멀어질수록, 그 국소 유전형을 가진 계통이 집단에서 차지하는 비중이 크다**는 뜻

## Slide 3-4. het_hom_peak_frac — 이질성의 핵심 지표

- **정의**: 전체 이형접합(heterozygous) PASS SNV 중 VAF > 0.9인 것의 비율
- **threshold를 0.9로 정한 이유**: 이형접합 사이트가 진짜 균질 집단이라면 시퀀싱 노이즈로도 VAF가 0.9까지 튀는 경우는 드묾 — 0.9는 "명백히 한쪽으로 쏠린" 사이트만 보수적으로 골라내는 기준. 0.8/0.75로 낮추면 걸리는 loci 수는 늘지만 노이즈(우연한 변동) 혼입 위험도 커짐 — 그룹 간 비교의 신뢰도를 위해 보수적 기준(0.9) 채택
- **우리 데이터 값** (군별 평균, A<C<B<D 방향): A가 가장 낮고(=이형접합 유지, 이질적) D가 가장 높음(=치우침 多, 균질/클론 우세)
- **암유전체학과의 연결**: 이 지표는 암유전체학의 **HRD(Homologous Recombination Deficiency) score** 개념과 본질적으로 같은 범주 — LOH(Loss of Heterozygosity), LST(Large-scale State Transition), TAI(Telomeric Allelic Imbalance) 3개 성분으로 구성되며, 개별 원인 유전자가 아니라 "게놈 전체의 불안정성/재구성 총량"을 재는 지표라는 점에서 동일한 발상. 상업화 예: Myriad Genetics의 myChoice CDx. 참고문헌: Abkevich et al. 2012, Popova et al. 2012, Birkbak et al. 2012

## Slide 3-4b. 쏠린 지점은 어디서 겹치나 — 중첩 분석과 후보 loci

**핵심 메시지**: het_hom_peak(VAF>0.9 쏠림 지점)이 군마다 **얼마나 겹치는지**, 겹치는 지점이 어떤 유전자 근처인지 4군 교차검증했다. 결과는 그라디언트(A<C<B<D)와 일치하되, **"같은 loci가 반복"이 아니라 "경향성"** 이라는 점이 핵심.

**방법**: 각 웰의 3배치 모두 VAF>0.9인 사이트("3중겹침")를 구해 → 10kb 클러스터링 → GENCODE vM36으로 유전자 주석 → 4군 교집합/차집합 비교. (도구: `bcftools query` + python + `tabix`)

**군별 3중겹침 개수 (그라디언트와 일치):**
| 군 | 3배치 모두 VAF>0.9 사이트 | 고유 로커스 |
|---|---|---|
| A | 40 | 26 |
| B | 41 | 19 |
| C | 44 | 29 |
| **D** | **51** | **38** |
→ D가 겹침·고유 로커스 모두 최다 = "가장 균질화(우세 계통 뚜렷)"라는 결론과 일치.

**문턱값 0.9의 타당성 (부수 검증)**: 0.9→0.6으로 낮추면 겹침 개수는 늘지만(D 51→3,696) 우연 대비 배율은 27,906배→126배로 급락 → **0.9가 신호대잡음비 최강 지점**이라 유지.

**군간 교차비교 — 4범주:**
| 범주 | loci | 해석 |
|---|---|---|
| 전체 4군 공통 | 2 (Grin2a, Gm21149/43391 근처) | 처치 무관 **일반 노이즈**(콜링 불안정 위치) — 배제 |
| **B+D 공통** | 2 (**Hsf2bp**, **Tusc3**) | "약물 보호 실패" 조건 특이 후보 — Hsf2bp는 스트레스반응 연관, **가장 신뢰도 높음** |
| D 고유 | 38 (**Piwil1**=유전체 안정성, Speer4e1·chr8:55.8Mb 클러스터 등) | D 특이 후보(단, Gm-접두어/다카피 다수 → 아티팩트 가능성 병존) |
| A/B/C 고유 | 26/19/29 | 각 군 고유 |

**중요 한계 (정직성)**: 웰당 3배치뿐이라 개별 지점 재현성 낮음(배치간 12~18% 겹침). 여기 후보들은 **정식 통계검정 안 거친 탐색적 결과** — IGV 확인·더 큰 표본 재현 전까지 확정 발견 아님. **핵심은 개별 loci가 아니라 "쏠림 총량의 그라디언트"**(=Slide 3-4). 결과 파일: `loh_hotspot_cross_validation.tsv`.

## Slide 3-4c. Induction vs. Bottleneck-Drift — 구분의 한계와 통합 Working Model ★

**핵심 메시지**: 관측된 이질성 그라디언트(A<C<B<D)를 설명하는 두 메커니즘 — **①유도(Induction)**: 산화손상이 복구 과정에서 LOH 생성 자체를 늘림, **②병목(Bottleneck)-Drift**: 세포사멸이 유효집단크기(Ne)를 줄여 selection 없이도 클론 고정을 가속 — 은 서로 배타적이지 않다. 오히려 **같은 산화손상 사건이 갈라지는 두 갈래의 세포 운명**으로 봐야 하며, 두 메커니즘이 동시에 작동해 서로를 증폭시키는 것이 가장 그럴듯한 working model이다.

**두 가설 비교**:
| 가설 | 기전 | "처치가 만든다" vs "이미 그런 애들만 산다" |
|---|---|---|
| 유도(Induction) | H2O2 손상 → DSB/재조합 복구 중 LOH 신규 생성 → 처치가 세고 길수록 누적 | "처치가 만든다"에 가까움 |
| 병목(Bottleneck)-Drift | 세포사멸↑ → 유효집단크기(Ne)↓ → **selection 없이도** drift만으로 클론 고정 가속 | "이미 그런 애들만 산다"의 확장판(순수 우연도 포함, selection 불필요) |

**왜 지금 데이터로 구분이 안 되는가**: het_hom_peak_frac(VAF>0.9) 같은 bulk 지표는 한 로커스가 검출되려면 그 LOH를 가진 계통이 집단의 **~80% 이상**을 차지해야 한다는 수학적 제약이 있다(VAF≈0.5+0.5f, f=계통 비율). 즉 이 지표는 본질적으로 "**LOH 누적량 × 클론 우세도(f)**"가 곱해진 값이라, 우세도만 커져도(=병목/drift만으로도, LOH 생성률은 그대로여도) 같은 그라디언트가 나올 수 있다. 두 메커니즘이 원리상 완전히 다른데도 **관측 가능한 결과가 같다**는 게 핵심 난제.

**새로 확인한 근거 (2026-07-28 재분석 — VAF 분포/ROH 교차검증)**:
- VAF 분포에서 **아직 고정 안 된 subclonal 구간(VAF 0.6~0.85)도 정확히 같은 그라디언트**로 나타남(Kruskal-Wallis H=10.38, raw p=0.0156, FDR=0.039) — 우세도만 커진 게 아니라 "진행 중인" LOH 신호 자체도 늘었다는 쪽에 힘을 보탬(induction 쪽 정황 강화)
- 이형접합 SNV 총 개수(n_het) 자체가 A(~31만)→D(~10만)로 **3~4배 차이**(H=10.38, p=0.0156) — 이 데이터셋에서 가장 크고 깨끗한 그라디언트
- bcftools ROH 교차검증(잔여 이형접합 구간 길이 non_roh_bp, ROH 블록 수)은 **유의하지 않음**(p=0.41, p=0.14) — n=3 검정력 한계 또는 방법론적 민감도 차이로 추정, 정직하게 병기(교차검증이 실패한 것을 숨기지 않음)
- 결과 파일: `vaf_distribution_summary.csv`, `vaf_distribution_hist.csv`, `roh_summary.csv`, `loh_evidence_kruskal_wallis.csv`, 그림 `loh_evidence_svg/`

**문헌적 근거**: 유도 가설은 mitotic recombination이 H2O2 유발 변이 빈도 증가분의 87%를 차지한다는 연구(효모, Cox lab 계열) 등으로 뒷받침되고, 병목 가설은 유효집단크기 감소가 selection과 무관하게 유전적 다양성을 급격히 낮춘다는 집단유전학·종양진화 문헌으로 뒷받침된다. 두 원리 모두 각자 확립돼 있으나, **포유류 성상세포·산화스트레스 처치 맥락에서 둘이 함께 작동한다는 통합 모델은 이 연구가 처음 제안**하는 프레임이다.

**한계**: 지금 bulk WGS 데이터는 원천적으로 두 기여분을 분리하지 못한다. 최종 답은 clonal barcoding(계통 추적) 실험의 몫이다(Slide 3-7 참고).

**모식도 제안 — 그림 A (슬라이드 본문용, 간단)**: "구분 불가능성" 개념도. 상단에 관측 결과 박스("클론 우세화 그라디언트 A<C<B<D"), 그 아래에서 위로 화살표 2개가 올라오는 구조 — 왼쪽 화살표 출발점 "가설① 유도: LOH 생성률↑", 오른쪽 화살표 출발점 "가설② 병목-Drift: 우세도↑(selection 불필요)". 두 화살표 사이에 물음표 아이콘과 "❓ 현재 데이터로는 구분 불가 — 두 화살표 다 같은 결과를 만듦" 캡션, 그 옆에 "✓ 상호 배타적 아님, 둘 다 기여 가능" 태그.

**모식도 제안 — 그림 B (통합 Working Model, 상세 스펙은 부록 C 참고)**: 산화스트레스 → 단일세포 운명 분기(정확한 복구/재조합 복구=LOH 생성/복구실패=세포사멸) → 집단 수준 결과(세포수 감소 + 특정 계통 비중 증가)까지 잇는 2단 구성 다이어그램. AI 다이어그램 생성 도구에 그대로 넘길 수 있는 상세 스펙은 **부록 C: Working Model 다이어그램 생성 가이드** 참고.

## Slide 3-5. 5개 독립 검증 방법 총정리 (표)

| # | 방법 | 사용 도구/입력 | 무엇을 측정하는가 | HiFi 고유 기능 활용 |
|---|---|---|---|---|
| 1 | het_hom_peak_frac | `sample_level_summary.py`, phased VCF | 이형접합 SNV 중 VAF>0.9 비율 | 아니오 (short-read로도 가능) |
| 2 | Joint genotype 불일치율 | `joint_genotype_discordance.py`, **GLnexus joint VCF**(Slide 1-5b/2-1b) | 동일 군 3반복 간 유전형 불일치 정도(균질할수록 낮음) | 아니오 |
| 3 | 카피수 안정성 | `copynum_stability_check.py`, `sv_copynum_bedgraph` | 다카피 유전자좌의 카피수 변동폭 | 아니오 |
| 4 | **Haplotype 방향성 allelic imbalance** | `haplotype_imbalance.py`, PS 블록(phased VCF) | Phase block 내 여러 SNP가 같은 방향으로 일관되게 치우치는지 | **예 — 완전히 독립적인 신규 증거** |
| 5 | Allele-specific methylation(ASM) 후보 영역 수 | `asm_analysis.R`(DSS), `cpg_hap1_bed`/`cpg_hap2_bed` | hap1 vs hap2 메틸화 차이 영역 수 | 예 (단, n_het_snv와 부분적으로 연동돼 완전 독립은 아님) |
| (참고) | Read 단위 HP 미배정률 | `read_level_haplotype_consistency.py`, haplotagged BAM | 이형접합 밀도 부족 시 read의 haplotype 배정 실패율 | 예 (단, 독립 증거라기보다 이형접합 밀도 감소의 재확인) |

**정직한 caveat**: 5개 방법 모두 같은 방향(A<C<B<D)을 가리키지만, 방법 4(haplotype 방향성 allelic imbalance)만 원리상 완전히 독립적인 증거다. 방법 5와 참고 지표는 결과적으로 "이형접합 밀도 감소"라는 이미 알려진 원인을 다른 측정 방식으로 재확인한 것에 가깝다 — 그래도 메틸화·read phasing 차원에서 같은 패턴이 일관되게 재현된다는 것 자체는 그라디언트의 신뢰도를 보강한다.

## Slide 3-6. 결과 파일 구조

| 파일 | 구조 |
|---|---|
| `{sample}.hap_imbalance.tsv` | 샘플별 PS 블록 단위 이항검정 결과: `n_ps_blocks_tested`, `n_biased_blocks`, `frac_biased`, `mean_abs_bias` |
| `hp_consistency_combined.tsv` | 샘플별 대표 8구간(5Mb×8) HP 태그 배정 결과: `unassigned_frac` |
| `{sample}.asm_results.csv` | hap1 vs hap2 메틸화 비교 영역별 결과 (DSS 통계량 포함) |
| `asm_summary.tsv` | 샘플별 유의 ASM 영역 수(`n_asm_regions`) 요약 |
| `summary_compare.csv` / `_rank_matrix.csv` / `_friedman.csv` | Kruskal-Wallis 결과, 군별 순위행렬, Friedman/Kendall's W 검정 결과 |

## Slide 3-7. 생물학적 검증 계획 (n=3 WGS의 한계 보완)

**핵심 메시지**: 현재까지의 발견은 전부 **bulk WGS(군당 n=3)로부터의 계산적 추론**이며, 독립적인 생물학적 검증이 필요하다.

| # | 방법 | 답하는 질문 | 기존 gDNA 재사용 | 비용/부담 | 소요기간 |
|---|---|---|---|---|---|
| 1 | Flow cytometry (DNA 함량/배수성) | 세포집단이 실제로 균질화되는가? | 불가 | 낮음 | 1~2일 |
| 2 | Single-cell WGS/scRNA-seq + CNV 추론 | 몇 개 클론이 있고 비율은? | 불가 | 매우 높음 | 수주 |
| 3 | Locus-specific FISH | 후보 유전자 카피수가 세포마다 다른가? | 불가 | 중간~높음 | probe 수주+실험 수일 |
| 4 | **Digital PCR (ddPCR)** | 후보 유전자 쏠림이 더 큰 n에서 재현되는가? | **가능(최대 장점)** | 낮음~중간 | probe 1~2주+실험 1일 |
| 5 | Clonal barcoding (lineage tracing) | 처치가 기존 계통을 선택하는가, 새로 유도하는가? | 완전 불가 | 높음 | 수개월(신규실험) |
| 6 | Time-course 샘플링 | 균질화가 점진적인가 급격한가? | 부분적 | 중간~높음 | 기존과 유사+시점관리 |

**권장 순서**: ① Flow cytometry로 균질화 자체 스크리닝 → ② ddPCR(기존 gDNA 재사용 가능)로 후보 유전자(Hsf2bp, Piwil1 등) 재현성 확인 → ③ FISH로 세포 단위 직접 검증 → ④ (예산 허락 시) single-cell 분석 또는 clonal barcoding으로 "유도 vs 선택" 논쟁에 최종 답.

## Slide 3-8. 용어집 (Glossary — 별도 슬라이드)

| 용어 | 설명 |
|---|---|
| **VAF** (Variant Allele Frequency) | 특정 위치에서 변이 대립유전자를 지지하는 read 비율 (ALT read 수 / 전체 read 수) |
| **Bulk VAF** | 다수 세포를 한꺼번에 시퀀싱했을 때 관측되는 VAF — 개별 세포 유전형이 아니라 세포집단 전체의 "가중 평균" 신호 |
| **단일세포 VAF vs Bulk VAF** | 단일세포는 이형접합 시 VAF=0.5(또는 0/1의 이산값)이지만, bulk는 세포 비율에 따라 0~1 사이 연속값으로 나타남 |
| **LOH** (Loss of Heterozygosity, 이형접합소실) | 원래 이형접합(2개 다른 대립유전자 보유)이었던 자리가 어떤 세포 계통에서 한쪽 사본을 잃어(또는 다른 쪽으로 대체되어) 동형접합처럼 되는 현상. 본 연구에서는 LOH 자체가 균질화의 "원인"이 아니라, LOH를 가진 계통이 집단에서 수적으로 우세해질 때 bulk 신호에 "드러나는" 마커로 해석 |
| **het_hom_peak_frac** (hom-peak %) | 이형접합 SNV 중 VAF>0.9(사실상 동형접합처럼 보이는 지점)의 비율 — 세포집단 균질화/클론 우세 정도를 재는 핵심 지표 |
| **Phasing (위상 분석)** | 한 개체가 부모로부터 물려받은 두 벌의 염색체(haplotype) 중, 각 변이가 어느 쪽에 속하는지 재구성하는 작업 |
| **Haplotype** | 부모 한쪽으로부터 함께 물려받아 같은 염색체 사본에 있는 변이들의 조합(세트) |
| **PS (Phase Set)** | Phasing이 성공적으로 연결된 하나의 연속 구간(블록)을 식별하는 ID — 같은 PS를 가진 변이들은 서로 상대적 위상(어느 haplotype인지)을 알 수 있음 |
| **HP 태그** | 개별 read가 haplotype 1 또는 2 중 어디서 유래했는지 표시하는 BAM 파일 태그 |
| **ASM** (Allele-Specific Methylation) | 같은 샘플 내에서 haplotype 1과 haplotype 2 간 메틸화 수준 차이 |
| **DMR** (Differentially Methylated Region) | 그룹/조건 간(또는 haplotype 간) 메틸화 수준이 유의하게 다른 유전체 영역 |
| **GT** (Genotype) | 유전형 표기, 예: 0/0(동형 참조), 0/1(이형접합), 1/1(동형 변이) |
| **DP** (Depth) | 특정 위치를 커버한 총 read 수 |
| **AD** (Allelic Depth) | 대립유전자별(REF, ALT) 지지 read 수 |
| **참조 편향** (Reference bias) | 정렬 시 ALT 대립유전자를 가진 read도 전체 서열 유사도로 정상 정렬되는 현상 — "정렬이 REF만 편애해서 VAF가 왜곡된다"는 오해를 방지하기 위한 개념 |
| **SV** (Structural Variant, 구조변이) | 50bp 이상의 큰 규모 유전체 변화(결실, 중복, 역위, 전좌 등) |
| **CNV** (Copy Number Variant, 카피수 변이) | 특정 유전체 구간의 카피(사본) 개수가 정상(2개)과 다른 경우 |
| **STR** (Short Tandem Repeat, 반복서열) | 짧은 서열 단위가 여러 번 반복되는 유전체 영역 |
| **gVCF** | 변이 위치뿐 아니라 "정상 확인된 구간"도 블록으로 기록하는 VCF 확장 포맷 |
| **Joint calling** (합동 유전형 판정) | 여러 샘플의 gVCF를 합쳐 모든 샘플을 모든 변이 위치에서 다시 유전형 판정 → 빈 칸 없는 정사각 행렬. 샘플 간 비교(불일치율·코호트)에 필수(Slide 2-1c) |
| **bcftools** | VCF/BCF를 읽고·거르고·값을 뽑고·합치는 표준 커맨드라인 도구(htslib 계열). VCF 분석의 '값 추출 앞단' 역할(Slide 2-1d) |
| **Kruskal-Wallis 검정** | 3개 이상 그룹 간 분포 차이를 순위 기반으로 비교하는 비모수 통계 검정 (본 연구처럼 군당 n=3인 소규모 표본에 적합) |
| **Kendall's W (Friedman 검정)** | 여러 지표(측정 방법)에 걸쳐 그룹 순위가 얼마나 일관되는지 재는 지표 — 여러 독립 증거가 "같은 결론"을 가리키는지 정량화 |
| **HRD score** (참고: 암유전체학 유사 개념) | LOH+LST(대규모 상태 전환)+TAI(말단부 대립유전자 불균형) 3성분으로 게놈 전체 불안정성을 재는 임상 지표(Myriad myChoice CDx 등에 활용) |

---

# 부록: 슬라이드 구성 제안 순서 (전체 발표 흐름)

1. 표지
2. (Part1) HiFi 시퀀싱 배경 — Slide 1-1
3. (Part1) 파이프라인 개요 플로우차트 — Slide 1-2
4. (Part1) 정렬~변이호출 단계별 — Slide 1-3~1-9 (필요시 통합/축약)
5. (Part1) 도구 요약표 — Slide 1-10
6. (전환) 실험 설계 소개 — Slide 3-1
7. (Part2) VCF 포맷 설명 — Slide 2-1~2-1b, 2-2
8. (Part2) 개별 VCF vs Joint calling 비교 — Slide 2-1c
9. (Part2) bcftools 소개 — Slide 2-1d
10. (Part2) SV/small variant 결과 파일 및 도구 재정리 — Slide 2-3~2-5
11. (Part2) Burden test 결과(FDR 소멸) — Slide 2-5b
12. (Part2, 보조/선택) 다카피 정렬 아티팩트 — Slide 2-5c
13. (Part2) 시각화 파일 목록 — Slide 2-6
14. (Part3) VAF 개념 + 역추론 논리 — Slide 3-2~3-3
15. (Part3) het_hom_peak_frac + HRD 유사성 — Slide 3-4
16. (Part3) 쏠린 지점 중첩·후보 loci 교차검증 — Slide 3-4b
16.5. (Part3) **Induction vs. Bottleneck-Drift 통합 Working Model** — Slide 3-4c, 그림 A(구분 불가능성 개념도) + 그림 B(working model, 부록 C 스펙으로 생성) 삽입
17. (Part3) **5개 독립 검증 결과 종합 (핵심 결론 슬라이드)** — Slide 3-5, `heterogeneity_consensus_svg/01_heterogeneity_consensus_heatmap.svg` (또는 `fig_atlas_bundles/heterogeneity_consensus_heatmap/outputs/figure.png`) 삽입
18. (Part3) 결과 파일 구조 — Slide 3-6
19. (Part3) 향후 생물학적 검증 계획 — Slide 3-7
20. 용어집(Glossary) — Slide 3-8
21. 결론 및 향후 계획

---

# 부록 C: Working Model 다이어그램 생성 가이드 (AI 디자인 도구 전달용)

**용도**: Slide 3-4c의 그림 B(induction+bottleneck 통합 working model)를 AI 다이어그램/이미지 생성 도구에 그대로 붙여넣어 초안을 만들기 위한 상세 스펙. 아래 지시문은 도구에 구애받지 않도록(Mermaid, Claude 이미지 생성, 일러스트레이터용 프롬프트 등 무엇이든) 텍스트로만 구조·요소·색상·라벨을 명시했다.

## 전체 목표

"산화스트레스 → 두 메커니즘(유도 + 병목) → 클론 우세화 → 관측된 이질성 그라디언트"를 하나의 흐름도로 표현한다. **두 메커니즘이 같은 원인(세포 하나가 겪는 DNA 손상)에서 갈라지는 "세포 운명 분기"의 다른 결과라는 점**을 시각적으로 강조하는 것이 이 그림의 핵심 메시지다 (즉 "유도 상자"와 "병목 상자"를 나란히 병렬로 그리지 말고, 반드시 하나의 분기 트리에서 갈라져 나오는 모습으로 그릴 것).

## 레이아웃 — 세로 방향 2단 구성 (Panel A 위, Panel B 아래, 화살표로 연결)

### Panel A (상단) — 단일 세포 운명 분기 트리

1. 최상단 박스: **"정상 세포 (다수 이형접합 유지)"** — 색상: 중립 회색(#888888), 사각형 노드
2. 아래로 화살표 1개, 화살표 라벨: **"H2O2 산화손상 (DSB, 염기손상)"**
3. 도착 노드: **"손상 세포"** — 주황/붉은 테두리로 강조(경고를 뜻하는 색, 예: 테두리 #C0392B)
4. 이 "손상 세포" 노드에서 **3개의 화살표가 부채꼴로 갈라짐(fork)**, 각각 다음과 같은 라벨과 도착 박스를 가짐:
   - 화살표 라벨 **"정확한 복구 (BER/NHEJ)"** → 박스 **"이형접합 유지"** (중립 회색, "안전" 경로)
   - 화살표 라벨 **"재조합 기반 복구 (HR/mitotic recombination)"** → 박스 **"LOH 발생 (신규 마커 획득)"** — **주황색(#E67E22)으로 강조**, 박스 위에 작은 태그 **"[유도 Induction]"**
   - 화살표 라벨 **"복구 실패 / 손상 과다"** → 박스 **"세포사멸"** — **파란색(#3498DB)으로 강조**, 박스 위에 작은 태그 **"[병목 Bottleneck 기여]"**
5. Panel A 하단에 작은 캡션 텍스트: *"같은 손상 사건이 induction의 원료(②)와 bottleneck의 압력(③)을 동시에 만든다 — 두 메커니즘은 별개 현상이 아니라 한 세포 운명 분기의 다른 결과"*

### Panel B (하단) — 집단(웰) 수준 결과

1. 왼쪽에 원(집단) 하나: **"처치 전 집단"** — 원 안에 작은 점 20~30개, 대부분 회색 점(정상), 그중 1~2개만 주황 점(우연히 존재하는 LOH 보유 세포)
2. 가로 화살표, 라벨: **"처치 강도 증가: A → C → B → D"**
3. 오른쪽에 **작은 원 4개를 나란히** 배치해 그라디언트를 표현 (왼쪽부터 A, C, B, D 순서 — 그라디언트 순서와 동일하게):
   - **A**: 점 개수 최다(예: 25개), 회색 다수 + 주황 소수, 고르게 섞임 → 라벨 "이질적"
   - **C**: 점 개수 약간 감소, 주황 비중 약간 증가
   - **B**: 점 개수 더 감소, 주황 비중 더 증가
   - **D**: 점 개수 최소(예: 8개), 대부분 주황(=한 계통이 지배) → 라벨 "균질/클론 우세"
   - 각 원의 테두리 색은 기존 발표 자료의 그룹 색상과 통일: **A=#E41A1C, B=#377EB8, C=#4DAF4A, D=#984EA3**
   - 각 원 아래 작은 캡션: *"세포 수↓(병목 심화) + 주황 비중↑(LOH 계통 우세)"*
4. Panel B 하단에 결과 readout 박스 (점선 테두리): **"→ 관측 결과: VAF 쏠림 증가(subclonal 0.6~0.85 + 고정 >0.9), 이형접합 SNV 수(n_het) 감소, het_hom_peak_frac↑"**
5. Panel A → Panel B로 이어지는 화살표 2개 (Panel A의 ②③ 박스에서 각각 아래로 내려와 Panel B로 진입):
   - ②(LOH 발생)에서 내려오는 화살표 라벨: **"다수 세포에서 반복 발생 → 원료 누적 (주황 점 개수/색 강도)"**
   - ③(세포사멸)에서 내려오는 화살표 라벨: **"다수 세포에서 반복 발생 → 집단 크기 감소(Ne↓) → drift 가속"**

## 색상 가이드 (요약)

| 색상 | 용도 |
|---|---|
| 회색 #888888 | 중립/정상/이형접합 유지 |
| 주황 #E67E22 | 유도(Induction)/LOH 관련 요소 |
| 파랑 #3498DB | 병목(Bottleneck)/세포사멸 관련 요소 |
| A=#E41A1C, B=#377EB8, C=#4DAF4A, D=#984EA3 | Panel B의 4개 원 테두리(기존 슬라이드 그룹 색상과 통일) |

## 텍스트/스타일 지침

- 라벨 언어: 국문 우선, 필요시 영문 병기(예: "유도(Induction)")
- 전체 방향: 위→아래 흐름(Panel A→B), Panel B 내부는 좌→우 시간/강도 흐름
- 전체 캡션(그림 최하단): *"Working model: 산화스트레스에 의한 세포 손상은 (1) 재조합 복구 과정에서 LOH를 생성하는 induction 경로와 (2) 세포사멸을 통해 유효집단크기를 줄이는 bottleneck 경로로 동시에 분기하며, 두 경로 모두 관측된 클론 우세화(이질성 감소) 그라디언트에 기여한다. 두 경로의 상대적 기여도는 현재 bulk WGS 데이터로는 분리되지 않는다."*
- 톤: 학술 발표용 — 만화적 과장 없이 깔끔한 흐름도/개념도 스타일(예: Nature/Cell 리뷰 논문의 graphical abstract 톤)

## 그림 A(구분 불가능성 개념도) 미니 스펙 — 참고

Slide 3-4c 본문에 곁들일 더 단순한 보조 그림. Panel B 없이 Panel A 스타일의 초소형 버전:
- 상단 박스: "관측: 클론 우세화 그라디언트 (A<C<B<D)"
- 아래에서 위로 화살표 2개: 왼쪽 "가설① 유도 — LOH 생성률↑", 오른쪽 "가설② 병목-Drift — 우세도↑(selection 불필요)"
- 두 화살표 사이 물음표 아이콘 + 캡션 "❓ 현재 데이터로는 구분 불가 — 두 화살표 다 같은 결과를 만듦"
- 옆에 작은 태그: "✓ 상호 배타적 아님, 둘 다 기여 가능"
