#!/usr/bin/env Rscript
# =============================================================================
# ROH(bcftools roh) + VAF 분포 전체 모양(subclonal mid-range) 군간 비교.
# het_hom_peak_frac(VAF>0.9 이진 문턱)이 "우세도 x LOH누적량"을 못 가르는 문제를
# 보완하기 위한 두 가지 참고 지표:
#   1) frac_subclonal_mid_0.6_0.85: 아직 완전히 고정 안 된(=우세도 낮아도 보이는)
#      "진행 중" LOH 신호 비율
#   2) ROH non_roh_bp: bcftools roh HMM 기반 잔여 이형접합 구간 길이 (inbred 마우스라
#      total_roh_bp는 게놈 크기에 포화돼 변별력이 없어, 그 여집합을 사용)
#   3) (부가) n_het: 샘플당 이형접합 SNV 총 개수 자체 — 이형접합 밀도 감소 정도
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

source("scripts/report_common.R")

vaf_summary <- read.csv("h2o2_analysis_results/cohort/vaf_distribution_summary.csv", stringsAsFactors=FALSE)
vaf_hist    <- read.csv("h2o2_analysis_results/cohort/vaf_distribution_hist.csv", stringsAsFactors=FALSE)
roh_summary <- read.csv("h2o2_analysis_results/cohort/roh_summary.csv", stringsAsFactors=FALSE)

vaf_summary$group <- factor(vaf_summary$group, levels=GROUP_LEVELS)
roh_summary$group <- factor(roh_summary$group, levels=GROUP_LEVELS)
vaf_hist$group    <- factor(vaf_hist$group, levels=GROUP_LEVELS)

# ── Kruskal-Wallis: 지표별 군간 비교 ────────────────────────────────────────
kw_test <- function(df, col) {
  kt <- kruskal.test(df[[col]], df$group)
  data.frame(metric=col, statistic=unname(kt$statistic), df=unname(kt$parameter), pvalue=kt$p.value)
}

kw_results <- bind_rows(
  kw_test(vaf_summary, "n_het"),
  kw_test(vaf_summary, "frac_subclonal_mid_0.6_0.85"),
  kw_test(vaf_summary, "frac_peak_gt_0.9"),
  kw_test(roh_summary, "n_roh_blocks"),
  kw_test(roh_summary, "non_roh_bp")
)
kw_results$fdr <- p.adjust(kw_results$pvalue, method="BH")
write.csv(kw_results, "h2o2_analysis_results/cohort/loh_evidence_kruskal_wallis.csv", row.names=FALSE)
message("Kruskal-Wallis 결과:")
print(kw_results)

# ── 시각화 ───────────────────────────────────────────────────────────────────
plots <- list()

make_box <- function(df, yvar, ylab, title) {
  kw_row <- kw_results[kw_results$metric == yvar, ]
  means <- df %>% group_by(group) %>% summarise(v=mean(.data[[yvar]], na.rm=TRUE), .groups="drop")
  ggplot(df, aes(x=group, y=.data[[yvar]], color=group)) +
    geom_boxplot(outlier.shape=NA, width=0.5, linewidth=0.7) +
    geom_jitter(width=0.12, size=2.8, alpha=0.85) +
    geom_text(data=means, aes(x=group, y=v, label=sprintf("%.4g", v)),
              inherit.aes=FALSE, vjust=-1.6, size=3.2, color="grey20", fontface="italic", family=FONT_FAMILY) +
    scale_color_manual(values=GROUP_COLORS) +
    scale_x_discrete(labels=GROUP_LABELS) +
    labs(title=title,
         subtitle=wrap_txt(sprintf("Kruskal-Wallis H=%.2f, raw p=%.4f (n=3/군)", kw_row$statistic, kw_row$pvalue)),
         x=NULL, y=ylab) +
    base_theme
}

plots[["01_n_het"]] <- make_box(vaf_summary, "n_het", "샘플당 이형접합 SNV 수",
                                  "이형접합 밀도")
plots[["02_frac_subclonal_mid"]] <- make_box(vaf_summary, "frac_subclonal_mid_0.6_0.85",
                                               "VAF 0.6~0.85 비율 (진행중 LOH, 미고정)",
                                               "Subclonal(미고정) LOH 신호")
plots[["03_frac_peak"]] <- make_box(vaf_summary, "frac_peak_gt_0.9",
                                      "VAF>0.9 비율 (het_hom_peak_frac)",
                                      "고정된(우세) LOH 신호 — 기존 지표")
plots[["04_n_roh_blocks"]] <- make_box(roh_summary, "n_roh_blocks", "ROH 블록 수",
                                         "bcftools ROH 블록 수")
plots[["05_non_roh_bp"]] <- make_box(roh_summary, "non_roh_bp", "잔여 이형접합 가능 구간 길이 (bp)",
                                       "bcftools ROH 여집합(비-ROH) 길이")

# VAF 분포 전체 모양 (밀도 곡선, 군별 색상, 샘플은 옅은 선으로)
vaf_hist_mean <- vaf_hist %>%
  group_by(group, bin_start) %>%
  summarise(density=mean(density), .groups="drop")

plots[["06_vaf_density_by_group"]] <- ggplot() +
  geom_line(data=vaf_hist, aes(x=bin_start, y=density, group=sample, color=group), alpha=0.25, linewidth=0.4) +
  geom_line(data=vaf_hist_mean, aes(x=bin_start, y=density, color=group), linewidth=1.3) +
  geom_vline(xintercept=c(0.6, 0.85, 0.9), linetype="dashed", color="grey50") +
  scale_color_manual(values=GROUP_COLORS, labels=GROUP_LABELS, name="그룹") +
  coord_cartesian(xlim=c(0.4, 1.0)) +
  labs(title="이형접합 SNV VAF 분포 (군 평균, 굵은 선 = 군 평균, 옅은 선 = 개별 샘플)",
       subtitle="점선: 0.6/0.85(subclonal 구간 경계), 0.9(기존 hom-peak 문턱)",
       x="VAF", y="밀도 (샘플별 이형접합 SNV 중 비율)") +
  base_theme + theme(legend.position="right")

save_svg_plots(plots, "h2o2_analysis_results/cohort/loh_evidence_svg", width=12, height=8)
message("완료: h2o2_analysis_results/cohort/loh_evidence_svg, loh_evidence_kruskal_wallis.csv")
