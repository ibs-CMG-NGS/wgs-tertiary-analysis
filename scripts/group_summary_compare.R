#!/usr/bin/env Rscript
# =============================================================================
# 샘플단위 요약통계 군간비교 — 유전자별 다중검정 대신 지표당 1회 검정
# (1) 각 연속 지표에 대해 Kruskal-Wallis (n=3/group, 유전자burden보다 검정력 좋음)
# (2) 여러 지표에 걸친 군 순위가 일관되는지 Friedman test(=Kendall's W)로 확인
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

source("scripts/report_common.R")

option_list <- list(
  make_option("--metrics-tsv",     type="character", help="샘플별 지표 TSV (sample_level_summary.py 출력 합친 것)"),
  make_option("--discordance-tsv", type="character", help="웰별 joint genotype 불일치율 TSV"),
  make_option("--output-csv",      type="character", help="결과 CSV 경로"),
  make_option("--output-svg-dir",  type="character", help="시각화 SVG 저장 디렉토리")
)
opt <- parse_args(OptionParser(option_list=option_list))
out_svg_dir <- opt$`output-svg-dir`

metrics <- read.table(opt$`metrics-tsv`, header=TRUE, sep="\t", stringsAsFactors=FALSE)
disc    <- read.table(opt$`discordance-tsv`, header=TRUE, sep="\t", stringsAsFactors=FALSE)

per_sample_metric_cols <- c("total_pass_variants", "total_pass_sv", "het_hom_peak_frac", "median_het_vaf")

# ── (1) 지표별 Kruskal-Wallis (n=3/group) ───────────────────────────────────
kw_results <- lapply(per_sample_metric_cols, function(m) {
  kt <- tryCatch(kruskal.test(metrics[[m]], factor(metrics$group)), error=function(e) NULL)
  data.frame(metric=m,
             statistic=if(!is.null(kt)) unname(kt$statistic) else NA_real_,
             df=if(!is.null(kt)) unname(kt$parameter) else NA_real_,
             pvalue=if(!is.null(kt)) kt$p.value else NA_real_)
})
kw_df <- bind_rows(kw_results)
kw_df$fdr <- p.adjust(kw_df$pvalue, method="BH")

# ── (2) 군별 평균(지표) + 불일치율을 합쳐 순위행렬 구성 ──────────────────────
# 주의(방향성): het_hom_peak_frac/discordance_rate는 둘 다 "세포집단 이질성/클론성"을
# 직접 재는 확립된 지표라 순위일관성(Kendall's W) 검정에 포함한다.
# 반면 total_pass_variants/total_pass_sv/median_het_vaf는 시퀀싱 배치·커버리지 등
# 기술적 요인에도 영향받을 수 있어(이전 QC 분석에서 확인된 confound) 이질성 지표로
# 확정할 근거가 약함 — 개별 Kruskal-Wallis 결과표에는 남기되 순위일관성 검정에서는 제외.
group_names <- sort(unique(metrics$group))
n_groups <- length(group_names)
heterogeneity_metric_cols <- c("het_hom_peak_frac")  # 방향: 값이 클수록 균질(클론성) ↑

group_means <- metrics %>%
  group_by(group) %>%
  summarise(across(all_of(heterogeneity_metric_cols), mean), .groups="drop")

# 불일치율(웰당 1개 값)도 같은 행렬에 추가 (방향: 값이 클수록 이질성 ↑ — het_hom_peak_frac과 반대 방향)
disc_row <- disc %>% select(group, discordance_rate) %>% rename(value=discordance_rate) %>% mutate(metric="discordance_rate")

long_means <- group_means %>%
  pivot_longer(-group, names_to="metric", values_to="value") %>%
  bind_rows(disc_row)

# 각 지표(행)마다 그룹 순위 매기기 (1=최소, n_groups=최대)
# het_hom_peak_frac은 "클론성"(균질) 방향이므로, "이질성" 방향으로 통일하기 위해 순위를 반전
rank_df <- long_means %>%
  group_by(metric) %>%
  mutate(rank=rank(value)) %>%
  mutate(rank=if_else(metric=="het_hom_peak_frac", n_groups+1-rank, rank)) %>%
  ungroup()

rank_matrix <- rank_df %>%
  select(metric, group, rank) %>%
  pivot_wider(names_from=group, values_from=rank) %>%
  as.data.frame()
rownames(rank_matrix) <- rank_matrix$metric
rank_mat_numeric <- as.matrix(rank_matrix[, group_names])

write.csv(rank_matrix, sub("\\.csv$", "_rank_matrix.csv", opt$`output-csv`), row.names=FALSE)

# ── Friedman test (= Kendall's W 검정, m개 지표를 "judges"로, 그룹을 "treatment"로) ──
# 주의: m(지표 개수)이 작아 카이제곱 근사가 정확하지 않을 수 있음 — 참고용 지표로 취급.
friedman_res <- tryCatch(friedman.test(rank_mat_numeric), error=function(e) NULL)
n_groups <- ncol(rank_mat_numeric)
n_metrics <- nrow(rank_mat_numeric)
kendalls_w <- if (!is.null(friedman_res)) {
  unname(friedman_res$statistic) / (n_metrics * (n_groups - 1))
} else NA_real_

friedman_df <- data.frame(
  n_metrics = n_metrics,
  n_groups  = n_groups,
  kendalls_w = kendalls_w,
  chisq_statistic = if(!is.null(friedman_res)) unname(friedman_res$statistic) else NA_real_,
  df = if(!is.null(friedman_res)) unname(friedman_res$parameter) else NA_real_,
  pvalue_asymptotic = if(!is.null(friedman_res)) friedman_res$p.value else NA_real_,
  note = "m(지표 수)이 작아 카이제곱 근사는 참고용 — 정확검정 아님"
)

write.csv(kw_df, opt$`output-csv`, row.names=FALSE)
write.csv(friedman_df, sub("\\.csv$", "_friedman.csv", opt$`output-csv`), row.names=FALSE)

message("Kruskal-Wallis 결과:")
print(kw_df)
message("\n순위행렬:")
print(rank_matrix)
message("\nFriedman/Kendall's W 결과:")
print(friedman_df)

# ── 시각화 ───────────────────────────────────────────────────────────────────
metrics$group <- factor(metrics$group, levels=group_names)
disc$group    <- factor(disc$group, levels=group_names)

METRIC_TITLES <- c(
  total_pass_variants = "샘플당 PASS 소변이 수",
  total_pass_sv       = "샘플당 PASS SV 수",
  het_hom_peak_frac   = "Hom-peak % (이형접합 중 VAF>0.9 비율)",
  median_het_vaf      = "이형접합 SNV VAF 중앙값"
)

make_metric_boxplot <- function(m) {
  kw_row <- kw_df[kw_df$metric == m, ]
  means <- metrics %>% group_by(group) %>% summarise(v=mean(.data[[m]]), .groups="drop")
  ggplot(metrics, aes(x=group, y=.data[[m]], color=group)) +
    geom_boxplot(outlier.shape=NA, width=0.5, linewidth=0.7) +
    geom_jitter(width=0.12, size=2.8, alpha=0.85) +
    geom_text(data=means, aes(x=group, y=v, label=sprintf("%.4g", v)),
              inherit.aes=FALSE, vjust=-1.6, size=3.2, color="grey20", fontface="italic", family=FONT_FAMILY) +
    scale_color_manual(values=GROUP_COLORS) +
    scale_x_discrete(labels=GROUP_LABELS) +
    labs(title=METRIC_TITLES[[m]],
         subtitle=wrap_txt(sprintf("Kruskal-Wallis H=%.2f, p=%.4f (FDR=%.4f)",
                                    kw_row$statistic, kw_row$pvalue, kw_row$fdr)),
         x=NULL, y=m) +
    base_theme
}

p_metrics <- lapply(per_sample_metric_cols, make_metric_boxplot)
combined_metrics <- (p_metrics[[1]] | p_metrics[[2]]) / (p_metrics[[3]] | p_metrics[[4]]) +
  plot_annotation(
    title="샘플단위 요약통계 — 4개 지표 군별 비교",
    subtitle="유전자별 다중검정 대신 지표당 1회 검정 (n=3/군) — hom-peak%가 이질성 그라디언트의 핵심 지표",
    theme=theme(plot.title=element_text(face="bold", size=15, family=FONT_FAMILY),
                plot.subtitle=element_text(size=10.5, color="grey30", family=FONT_FAMILY))
  )

# 불일치율 barplot (웰당 1개 값이라 검정 없이 관측치만 표시)
p_disc <- ggplot(disc, aes(x=group, y=discordance_rate, fill=group)) +
  geom_col(width=0.6) +
  geom_text(aes(label=sprintf("%.3f", discordance_rate)), vjust=-0.6, size=3.5, family=FONT_FAMILY) +
  scale_fill_manual(values=GROUP_COLORS) +
  scale_x_discrete(labels=GROUP_LABELS) +
  labs(title="웰별 Joint Genotype 불일치율",
       subtitle="웰당 1개 관측치(3배치 합쳐 계산) — 통계검정 아님, 순위 비교용",
       x=NULL, y="불일치율 (discordance_rate)") +
  base_theme + theme(plot.title=element_text(size=13))

# 순위행렬 히트맵 (지표 x 그룹, 색=순위) — 여러 지표가 같은 군 순서를 가리키는지 시각화
rank_long <- rank_df
rank_long$group <- factor(rank_long$group, levels=group_names)
METRIC_HEATMAP_LABELS <- c(het_hom_peak_frac="Hom-peak %\n(방향 반전)", discordance_rate="Joint genotype\n불일치율")
rank_long$metric_label <- METRIC_HEATMAP_LABELS[rank_long$metric]

p_heat <- ggplot(rank_long, aes(x=group, y=metric_label, fill=rank)) +
  geom_tile(color="white", linewidth=1) +
  geom_text(aes(label=rank), color="black", size=5, family=FONT_FAMILY) +
  scale_fill_gradient(low="steelblue", high="firebrick", name="이질성\n순위", limits=c(1, n_groups)) +
  scale_x_discrete(labels=GROUP_LABELS) +
  labs(title="지표별 군 순위 일관성",
       subtitle=wrap_txt(sprintf("Kendall's W=%.3f, Friedman 근사 p=%.4f (지표 수 %d개, m이 작아 근사는 참고용)",
                         kendalls_w, friedman_df$pvalue_asymptotic, n_metrics)),
       x=NULL, y=NULL) +
  base_theme + theme(plot.title=element_text(size=13))

combined_disc_heat <- p_disc | p_heat
combined_disc_heat <- combined_disc_heat +
  plot_annotation(
    title="이질성 지표 재현성 — 웰간 불일치율과 순위 일관성",
    theme=theme(plot.title=element_text(face="bold", size=15, family=FONT_FAMILY))
  )

save_svg_plots(
  list("01_metrics_by_group"=list(plot=combined_metrics, width=13, height=10),
       "02_discordance_rank"=list(plot=combined_disc_heat, width=13, height=6)),
  out_svg_dir
)

message("시각화 저장: ", out_svg_dir)
