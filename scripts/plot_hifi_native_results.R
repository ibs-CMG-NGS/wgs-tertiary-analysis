#!/usr/bin/env Rscript
# =============================================================================
# HiFi 롱리드 고유 기능 분석(§docs/H2O2_hifi_native_analysis_plan_2026-07-16.md)
# 결과 시각화 — haplotype 방향성 allelic imbalance, read-level HP 일관성
# 고품질(통계주석·일관된 스타일·다패널 조합) 버전
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

source("scripts/report_common.R")

option_list <- list(
  make_option("--input-dir",      type="character",
              help="cohort 출력 디렉토리 (hap_imbalance/, hp_consistency_combined.tsv, sample_metrics_combined.tsv 포함)"),
  make_option("--output-svg-dir", type="character", help="시각화 SVG 저장 디렉토리")
)
opt <- parse_args(OptionParser(option_list=option_list))
base <- opt$`input-dir`
out_svg_dir <- opt$`output-svg-dir`

theme_set(theme_bw(base_size=13, base_family=FONT_FAMILY))
update_geom_defaults("text", list(family=FONT_FAMILY))

# ── 데이터 로드 ──────────────────────────────────────────────────────────────
hap_imb <- do.call(rbind, lapply(
  list.files(file.path(base, "hap_imbalance"), pattern="\\.tsv$", full.names=TRUE),
  function(f) read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
))
hap_imb$group <- factor(hap_imb$group, levels=GROUP_LEVELS)

hp_cons <- read.table(file.path(base, "hp_consistency_combined.tsv"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
hp_cons$group <- factor(hp_cons$group, levels=GROUP_LEVELS)

make_boxplot <- function(df, yvar, ylab, title, subtitle_prefix, log_y=FALSE) {
  kw <- kw_pvalue(df[[yvar]], df$group)
  means <- df %>% group_by(group) %>% summarise(m=mean(.data[[yvar]]), .groups="drop")
  p <- ggplot(df, aes(x=group, y=.data[[yvar]], color=group)) +
    geom_boxplot(outlier.shape=NA, width=0.5, linewidth=0.7) +
    geom_jitter(width=0.12, size=2.8, alpha=0.85) +
    geom_text(data=means, aes(x=group, y=m, label=sprintf("%.4g", m)),
              inherit.aes=FALSE, vjust=-1.6, size=3.3, color="grey20", fontface="italic", family=FONT_FAMILY) +
    scale_color_manual(values=GROUP_COLORS) +
    scale_x_discrete(labels=GROUP_LABELS) +
    labs(title=title,
         subtitle=paste(strwrap(sprintf("%s | Kruskal-Wallis H=%.2f, raw p=%.4f (n=3/군)", subtitle_prefix, kw$H, kw$p),
                                 width=55), collapse="\n"),
         x=NULL, y=ylab) +
    base_theme + theme(plot.subtitle=element_text(size=9.5, family=FONT_FAMILY))
  if (log_y) p <- p + scale_y_log10()
  p
}

# ── 1. Haplotype 방향성 allelic imbalance ───────────────────────────────────
p1 <- make_boxplot(hap_imb, "frac_biased", "편향된 PS 블록 비율 (frac_biased)",
                    "Haplotype 방향성 Allelic Imbalance",
                    "PS 블록 내 phase-0 allele VAF 편향 (이항검정 p<0.05)")

p2 <- make_boxplot(hap_imb, "mean_abs_bias", "평균 |VAF 편향| (mean_abs_bias)",
                    "Haplotype 방향성 Allelic Imbalance — 편향 크기",
                    "블록별 평균 |phase-0 allele VAF - 0.5|")

p3 <- make_boxplot(hap_imb, "n_ps_blocks_tested", "검정 대상 PS 블록 수 (log10)",
                    "검정 가능했던 PS 블록 수",
                    "참고: 이형접합 밀도가 낮을수록(D군) 큰 PS 블록 자체가 적어짐", log_y=TRUE)

combined_hap <- (p1 | p2) / p3 +
  plot_annotation(
    title="1단계: Haplotype 방향성 Allelic Imbalance (PS 블록 기반)",
    subtitle="개별 SNP가 아니라 phase block(연결된 haplotype 구간) 단위로 일관된 VAF 편향을 검정 — §2(hom-peak%)와 원리가 다른 독립적 확인",
    theme=theme(plot.title=element_text(face="bold", size=16, family=FONT_FAMILY),
                plot.subtitle=element_text(size=11, color="grey30", family=FONT_FAMILY))
  )

# ── 2. Read-level HP 일관성 ──────────────────────────────────────────────────
p4 <- make_boxplot(hp_cons, "unassigned_frac", "HP 미배정 read 비율 (unassigned_frac)",
                    "Read-level Haplotype 배정 미배정률",
                    "대표 상염색체 8구간(5Mb×8) 샘플링 기준")

# 산점도: n_het_snv(이형접합 밀도, 재해석 근거)와의 관계 — sample_metrics_combined.tsv 참조
sample_metrics <- read.table(file.path(base, "sample_metrics_combined.tsv"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
merged <- merge(hp_cons, sample_metrics[, c("sample","n_het_snv")], by="sample")
merged$group <- factor(merged$group, levels=GROUP_LEVELS)

cor_test <- cor.test(merged$n_het_snv, merged$unassigned_frac, method="spearman")
p5 <- ggplot(merged, aes(x=n_het_snv, y=unassigned_frac, color=group)) +
  geom_point(size=3.5, alpha=0.9) +
  geom_smooth(inherit.aes=FALSE, aes(x=n_het_snv, y=unassigned_frac), method="lm", se=TRUE,
              color="grey40", linetype="dashed", linewidth=0.6) +
  scale_color_manual(values=GROUP_COLORS, labels=GROUP_LABELS, name="그룹") +
  scale_x_continuous(labels=scales::comma) +
  labs(title="HP 미배정률 vs 이형접합 밀도(n_het_snv)",
       subtitle=sprintf("Spearman rho=%.3f, p=%.4f — 강한 음의 상관관계는 '다클론 혼란' 가설이 아니라\n'이형접합 밀도 감소로 phasing 닻이 부족해짐' 재해석을 뒷받침", cor_test$estimate, cor_test$p.value),
       x="샘플당 이형접합 SNV 수 (n_het_snv)", y="HP 미배정 read 비율") +
  base_theme + theme(legend.position="right",
                      legend.text=element_text(family=FONT_FAMILY),
                      legend.title=element_text(family=FONT_FAMILY))

combined_hp <- p4 | p5
combined_hp <- combined_hp +
  plot_annotation(
    title="3단계: Read-level Haplotype 일관성 — 방향은 일치하나 해석은 정정됨",
    subtitle="원 가설(다클론 혼란→미배정↑)과 반대 방향 관측 → 이형접합 밀도 감소가 실제 원인으로 재해석(우측 패널 근거)",
    theme=theme(plot.title=element_text(face="bold", size=16, family=FONT_FAMILY),
                plot.subtitle=element_text(size=11, color="grey30", family=FONT_FAMILY))
  )

save_svg_plots(
  list("01_hap_imbalance"=list(plot=combined_hap, width=13, height=10),
       "02_hp_consistency"=list(plot=combined_hp, width=14, height=7)),
  out_svg_dir
)
message("시각화 저장: ", out_svg_dir)
