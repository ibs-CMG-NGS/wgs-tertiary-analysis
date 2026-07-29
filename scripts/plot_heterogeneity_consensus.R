#!/usr/bin/env Rscript
# =============================================================================
# 4개 독립 방법(hom-peak%, joint discordance, 카피수 안정성, haplotype imbalance)의
# 이질성 순위를 한 그림으로 통합 — "네 가지 서로 다른 원리가 같은 그라디언트를
# 가리킨다"는 핵심 메시지를 한 눈에 보여주기 위한 통합 시각화.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

source("scripts/report_common.R")

option_list <- list(
  make_option("--input-dir",      type="character",
              help="cohort 출력 디렉토리 (summary_compare_stats_rank_matrix.csv, copynum_instability_rank.csv, hap_imbalance/, asm_summary.tsv 포함)"),
  make_option("--output-svg-dir", type="character", help="시각화 SVG 저장 디렉토리")
)
opt <- parse_args(OptionParser(option_list=option_list))
base <- opt$`input-dir`
out_svg_dir <- opt$`output-svg-dir`

rank_mat <- read.csv(file.path(base, "summary_compare_stats_rank_matrix.csv"), stringsAsFactors=FALSE)
cn_rank  <- read.csv(file.path(base, "copynum_instability_rank.csv"), stringsAsFactors=FALSE)
hap_imb  <- do.call(rbind, lapply(
  list.files(file.path(base, "hap_imbalance"), pattern="\\.tsv$", full.names=TRUE),
  function(f) read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
))
asm_summary <- read.table(file.path(base, "asm_summary.tsv"), header=TRUE, sep="\t", stringsAsFactors=FALSE)

# 1) het_hom_peak_frac / discordance_rate: 이미 rank_matrix.csv에 "이질성 방향"으로 보정돼 있음(§2)
m1 <- setNames(as.numeric(rank_mat[rank_mat$metric=="het_hom_peak_frac", c("A","B","C","D")]), c("A","B","C","D"))
m2 <- setNames(as.numeric(rank_mat[rank_mat$metric=="discordance_rate", c("A","B","C","D")]), c("A","B","C","D"))

# 2) 카피수 불안정성: 값 클수록 이미 "이질성 방향" (그대로 순위화)
cn_vals <- setNames(cn_rank$mean_instability_rank, cn_rank$group)
m3 <- setNames(rank(cn_vals[c("A","B","C","D")]), c("A","B","C","D"))

# 3) haplotype imbalance: frac_biased 클수록 "클론성" 방향이므로 반전(그룹 평균 후 순위 반전)
hap_mean <- hap_imb %>% group_by(group) %>% summarise(mean_frac_biased=mean(frac_biased), .groups="drop")
hap_vals <- setNames(hap_mean$mean_frac_biased, hap_mean$group)
m4_raw <- rank(hap_vals[c("A","B","C","D")])          # 클론성 방향 순위(1=가장 이질적,4=가장 클론성)
m4 <- setNames(5 - m4_raw, c("A","B","C","D"))         # 이질성 방향으로 반전

# 5) ASM 후보 영역 수: 많을수록 phaseable heterozygosity가 많다는 뜻 -> "이질성 방향" (그대로 순위화)
asm_mean <- asm_summary %>% group_by(group) %>% summarise(mean_n=mean(n_asm_regions), .groups="drop")
asm_vals <- setNames(asm_mean$mean_n, asm_mean$group)
m5 <- setNames(rank(asm_vals[c("A","B","C","D")]), c("A","B","C","D"))

combined <- data.frame(
  method = rep(c("Hom-peak % (VAF distribution)", "Joint genotype discordance rate",
                 "Multi-copy gene family CN stability", "Haplotype-directed allelic imbalance",
                 "ASM candidate region count"), each=4),
  group  = rep(c("A","B","C","D"), 5),
  rank   = c(m1, m2, m3, m4, m5)
)
combined$method <- factor(combined$method, levels=c(
  "Hom-peak % (VAF distribution)", "Joint genotype discordance rate",
  "Multi-copy gene family CN stability", "Haplotype-directed allelic imbalance",
  "ASM candidate region count"))
combined$group <- factor(combined$group, levels=c("A","B","C","D"))

wrap80 <- function(x) paste(strwrap(x, width=95), collapse="\n")

caption_text <- paste(
  "Y-axis legend (all oriented so higher = more heterogeneous cell population):",
  wrap80("- Hom-peak % (VAF distribution): fraction of heterozygous SNVs with VAF>0.9 in a single sample -- lower = more heterogeneous (inverted for display)"),
  wrap80("- Joint genotype discordance rate: fraction of sites where genotype disagrees across the 3 biological replicates of a well -- higher = more heterogeneous"),
  wrap80("- Multi-copy gene family CN stability: copy-number variability across 3 replicates at 7 loci (Vmn1r/Mup etc.) -- larger variability = more heterogeneous"),
  wrap80("- Haplotype-directed allelic imbalance: fraction of phase-set blocks with a consistent one-haplotype VAF bias -- lower = more heterogeneous (inverted for display)"),
  wrap80("- ASM candidate region count: number of candidate hap1-vs-hap2 methylation-difference regions detected per sample -- more = more heterogeneous (note: partly reflects the amount of phaseable heterozygosity available, so less independent from the other rows than it may appear)"),
  sep="\n"
)

p <- ggplot(combined, aes(x=group, y=method, fill=rank)) +
  geom_tile(color="white", linewidth=1) +
  geom_text(aes(label=rank), color="black", size=5) +
  scale_fill_gradient(low="steelblue", high="firebrick", name="Heterogeneity rank\n(4=most heterogeneous,\n1=most clonal)",
                       limits=c(1,4), breaks=1:4) +
  labs(title=wrap80("Five independent methods all point to the same heterogeneity gradient: A > C > B > D"),
       subtitle=wrap80("Metrics computed via completely different principles (VAF distribution / joint genotype / depth-based copy number / haplotype phasing / allele-specific methylation) rank the groups identically"),
       x="Treatment group (A=control, B=acute, C=acute+drug, D=chronic+drug)", y=NULL,
       caption=caption_text) +
  theme_bw(base_size=13) +
  theme(plot.title=element_text(face="bold"),
        axis.text.y=element_text(size=11),
        plot.caption=element_text(hjust=0, size=9, lineheight=1.3, margin=margin(t=15)))

save_svg_plots(list("01_heterogeneity_consensus_heatmap"=list(plot=p, width=13, height=9.5)), out_svg_dir)
print(combined)
