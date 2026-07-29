#!/usr/bin/env Rscript
# =============================================================================
# fig-atlas 번들 생성기
# https://github.com/parkgilbong/fig-atlas-template 의
# manifests/external_pipeline_requirements.md 계약(완화판: figure.py 대신
# figure.R — metadata.yaml의 notes에 명시)에 맞춰, 선별된 핵심 그림을
# scripts/inputs/outputs/metadata 구조의 독립 재현 가능한 번들로 내보낸다.
#
# 각 번들의 scripts/figure.R은 이 스크립트가 문자열로 작성한 뒤 그대로
# `Rscript figure.R`로 실행해서 outputs/를 만든다 — 즉 "번들에 적힌 스크립트"와
# "실제로 그림을 만든 스크립트"가 항상 동일함이 보장된다 (계약서 §8 검증 요건).
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

source("scripts/report_common.R")

option_list <- list(
  make_option("--slug",          type="character", help="번들 slug"),
  make_option("--output-dir",    type="character", help="파이프라인 OUTPUT_DIR (예: h2o2_analysis_results)"),
  make_option("--bundle-dir",    type="character", help="번들을 생성할 디렉토리"),
  make_option("--figure-number", type="integer",   default=1, help="metadata.yaml figure_number")
)
opt <- parse_args(OptionParser(option_list=option_list))

slug       <- opt$slug
output_dir <- opt$`output-dir`
bundle_dir <- opt$`bundle-dir`
cohort_dir <- file.path(output_dir, "cohort")

# sv_burden_single_top20 / sv_burden_joint_top_hits는 둘 다 같은 파일명
# (cohort/sv_burden_stats.csv)을 읽으므로, --output-dir이 실제로 그 트리(싱글샘플
# vs joint)와 일치하는지 확인한다 — 안 그러면 joint 콜 데이터에 "single" 라벨이
# 붙거나, 반대로 싱글샘플 콜 데이터에 "joint 아티팩트" 주장이 잘못 붙을 수 있다.
is_joint_tree <- grepl("joint", output_dir, ignore.case=TRUE)
if (slug == "sv_burden_joint_top_hits" && !is_joint_tree) {
  stop("sv_burden_joint_top_hits는 joint 콜 트리(--output-dir에 'joint' 포함)에서만 빌드할 수 있습니다. 받은 --output-dir: ", output_dir)
}
if (slug == "sv_burden_single_top20" && is_joint_tree) {
  stop("sv_burden_single_top20은 싱글샘플 콜 트리에서만 빌드할 수 있습니다 (받은 --output-dir이 joint 트리로 보임): ", output_dir)
}

for (d in c("scripts", "inputs", "outputs", "metadata")) {
  dir.create(file.path(bundle_dir, d), recursive=TRUE, showWarnings=FALSE)
}
data_csv <- file.path(bundle_dir, "inputs", "data.csv")
figure_r <- file.path(bundle_dir, "scripts", "figure.R")

# ── 공통: 각 figure.R이 쓰는 스타일 상수 (report_common.R을 참조하지 않고
#    번들 자체에 직접 박아넣어 fig-atlas 저장소로 옮겨도 외부 의존 없이 동작하게 함) ──
style_header <- '
suppressPackageStartupMessages({ library(ggplot2); library(svglite) })
GROUP_COLORS <- c(A="#E41A1C", B="#377EB8", C="#4DAF4A", D="#984EA3")
GROUP_LABELS <- c(A="A (control)", B="B (acute)", C="C (acute+drug)", D="D (chronic+drug)")
df <- read.csv("../inputs/data.csv", stringsAsFactors=FALSE)
'

save_footer <- function(w, h) sprintf('
dir.create("../outputs", showWarnings=FALSE)
ggsave("../outputs/figure.png", p, width=%s, height=%s, dpi=300)
ggsave("../outputs/figure.pdf", p, width=%s, height=%s)
svglite::svglite("../outputs/figure.svg", width=%s, height=%s); print(p); grDevices::dev.off()
message("figure.png / figure.pdf / figure.svg written to ../outputs")
', w, h, w, h, w, h)

# ── slug 정의: 데이터 준비 + figure.R plotting 본문 + 메타데이터 ─────────────

if (slug == "snv_indel_burden_top30") {
  src_csv <- file.path(cohort_dir, "variant_burden_stats.csv")
  file.copy(src_csv, data_csv, overwrite=TRUE)
  plot_body <- '
top30 <- head(df[order(df$pvalue, na.last=TRUE), ], 30)
top30$log_p <- -log10(pmax(top30$pvalue, 1e-300))
top30$gene <- factor(top30$gene, levels=rev(top30$gene))
p <- ggplot(top30, aes(x=gene, y=log_p, fill=cramers_v)) +
  geom_col() +
  scale_fill_gradient(low="steelblue", high="firebrick", name="Cramer\'s V", limits=c(0,1)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
  coord_flip() +
  labs(title="Top 30 genes -- Small Variant Burden (raw p-value rank, color = effect size)",
       x=NULL, y="-log10(raw p-value)",
       subtitle="Dashed line = raw p 0.05 (not FDR-corrected)") +
  theme_bw(base_size=11)
'
  w <- 12; h <- 9
  meta <- list(
    figure_title = "SNV/indel gene-level burden — top 30 genes by raw p-value",
    panel_notes = c("Single panel: bar chart of top 30 genes ranked by raw chi-square/Fisher p-value, bar fill = Cramer's V effect size, dashed line at raw p=0.05."),
    input_files = "inputs/data.csv (= cohort/variant_burden_stats.csv, 80 genes tested)",
    statistics_tables = "inputs/data.csv",
    claim_boundary = "Raw p<0.05 for only 1/80 genes (B4galt3); FDR<0.05 for 0/80. No single gene separates the A/B/C/D groups after multiple-testing correction — this panel shows the raw-p ranking for exploratory/QC purposes only, not a validated hit list.",
    extra_notes = character(0)
  )
} else if (slug %in% c("sv_burden_single_top20", "sv_burden_joint_top_hits")) {
  src_csv <- file.path(cohort_dir, "sv_burden_stats.csv")
  file.copy(src_csv, data_csv, overwrite=TRUE)
  plot_body <- '
top20 <- head(df[order(df$pvalue, na.last=TRUE), ], 20)
top20$log_p <- -log10(pmax(top20$pvalue, 1e-300))
top20$gene <- factor(top20$gene, levels=rev(top20$gene))
p <- ggplot(top20, aes(x=gene, y=log_p, fill=cramers_v)) +
  geom_col() +
  scale_fill_gradient(low="steelblue", high="firebrick", name="Cramer\'s V", limits=c(0,1)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
  coord_flip() +
  labs(title="Top 20 genes -- SV Burden (raw p-value rank, color = effect size)",
       x=NULL, y="-log10(raw p-value)",
       subtitle="Dashed line = raw p 0.05 (not FDR-corrected)") +
  theme_bw(base_size=11)
'
  w <- 12; h <- 9
  if (slug == "sv_burden_single_top20") {
    claim <- "Raw p<0.05 for only 1/278 genes (Cyp2c29); FDR<0.05 for 0/278. No gene survives multiple-testing correction in the single-sample-call SV burden test."
  } else {
    claim <- "112/277 genes reach FDR<0.05 in the joint-called SV burden test, but this is judged to be a multi-copy gene-family alignment artifact (ambiguous read mapping across paralogs inflates apparent SV calls), not a genuine treatment effect. Do not cite these FDR values as evidence of a real signal without independent validation (e.g. long-read assembly or targeted re-sequencing at the flagged loci)."
  }
  meta <- list(
    figure_title = sprintf("SV gene-level burden — top 20 genes by raw p-value (%s)",
                            if (slug == "sv_burden_single_top20") "single-sample calls" else "joint calls"),
    panel_notes = c("Single panel: bar chart of top 20 genes ranked by raw chi-square/Fisher p-value, bar fill = Cramer's V effect size, dashed line at raw p=0.05."),
    input_files = sprintf("inputs/data.csv (= %s/sv_burden_stats.csv)", cohort_dir),
    statistics_tables = "inputs/data.csv",
    claim_boundary = claim,
    extra_notes = character(0)
  )
} else if (slug == "dmr_direction_by_chrom") {
  group_pair_csvs <- list.files(cohort_dir, pattern="^dmr_.*_vs_.*\\.csv$", full.names=TRUE)
  if (length(group_pair_csvs) == 0) stop("dmr_*_vs_*.csv 파일이 없습니다: ", cohort_dir)
  combined <- do.call(rbind, lapply(group_pair_csvs, function(f) {
    d <- read.csv(f, stringsAsFactors=FALSE)
    if (nrow(d) == 0) return(NULL)
    d$comparison <- sub("^dmr_(.*)\\.csv$", "\\1", basename(f))
    d
  }))
  write.csv(combined, data_csv, row.names=FALSE)
  plot_body <- '
p <- ggplot(df, aes(x=chr, fill=diff_methy > 0)) +
  geom_bar() +
  facet_wrap(~comparison, ncol=3) +
  scale_fill_manual(values=c("TRUE"="firebrick","FALSE"="steelblue"),
                     labels=c("Hypo","Hyper"), name="Direction") +
  labs(title="DMR count by chromosome, across all group-pair comparisons",
       x="Chromosome", y="DMR count") +
  theme_bw(base_size=10) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
'
  w <- 14; h <- 10
  meta <- list(
    figure_title = "DMR count by chromosome, all A/B/C/D pairwise comparisons",
    panel_notes = c("One facet per group-pair comparison (6 panels: A_vs_B .. C_vs_D); each bar = DMR count per chromosome, colored by hyper/hypo methylation direction."),
    input_files = sprintf("inputs/data.csv (concatenation of %d cohort/dmr_*_vs_*.csv files, with an added `comparison` column)", length(group_pair_csvs)),
    statistics_tables = "inputs/data.csv",
    claim_boundary = "DSS::callDMR() does not provide region-level p-values (pvalue/fdr columns are always NA) -- DMR selection here reflects only callDMR()'s internal per-CpG p-threshold + delta + minCG filters, not a corrected regional significance test. Counts are descriptive, not hypothesis-tested.",
    extra_notes = character(0)
  )
} else if (slug == "trgt_manhattan") {
  src_csv <- file.path(cohort_dir, "trgt_group_compare.csv")
  file.copy(src_csv, data_csv, overwrite=TRUE)
  plot_body <- '
manhattan_df <- df[!is.na(df$pvalue), ]
manhattan_df$chrom_num <- suppressWarnings(as.numeric(sub("chr", "", manhattan_df$chrom, ignore.case=TRUE)))
manhattan_df$pos_num   <- suppressWarnings(as.numeric(manhattan_df$pos))
manhattan_df <- manhattan_df[!is.na(manhattan_df$chrom_num) & !is.na(manhattan_df$pos_num), ]
manhattan_df <- manhattan_df[order(manhattan_df$chrom_num, manhattan_df$pos_num), ]
manhattan_df$log_p <- -log10(pmax(manhattan_df$pvalue, 1e-300))
offsets <- aggregate(pos_num ~ chrom_num, manhattan_df, max)
offsets <- offsets[order(offsets$chrom_num), ]
offsets$offset <- cumsum(c(0, head(offsets$pos_num, -1)))
manhattan_df <- merge(manhattan_df, offsets[, c("chrom_num", "offset")], by="chrom_num")
manhattan_df$abs_pos <- manhattan_df$pos_num + manhattan_df$offset
p <- ggplot(manhattan_df, aes(x=abs_pos, y=log_p, color=factor(chrom_num %% 2))) +
  geom_point(size=0.8, alpha=0.7) +
  scale_color_manual(values=c("0"="steelblue", "1"="navy"), guide="none") +
  labs(title="TRGT Manhattan Plot", x="Genomic position", y="-log10(p-value)") +
  theme_bw(base_size=10)
'
  w <- 12; h <- 9
  meta <- list(
    figure_title = "TRGT tandem-repeat locus comparison — Manhattan plot",
    panel_notes = c("Single panel: genome-wide Manhattan plot of -log10(p) per TRGT locus (Kruskal-Wallis across A/B/C/D allele lengths)."),
    input_files = "inputs/data.csv (= cohort/trgt_group_compare.csv)",
    statistics_tables = "inputs/data.csv",
    claim_boundary = "No Bonferroni/FDR significance line is drawn in this bundle version; treat as a descriptive overview of the raw p-value landscape, not a list of validated hits.",
    extra_notes = character(0)
  )
} else if (slug == "heterogeneity_consensus_heatmap") {
  rank_mat_path <- file.path(cohort_dir, "summary_compare_stats_rank_matrix.csv")
  cn_rank_path  <- file.path(cohort_dir, "copynum_instability_rank.csv")
  hap_imb_dir   <- file.path(cohort_dir, "hap_imbalance")
  asm_summary_path <- file.path(cohort_dir, "asm_summary.tsv")
  missing <- c(rank_mat_path, cn_rank_path, asm_summary_path)[!file.exists(c(rank_mat_path, cn_rank_path, asm_summary_path))]
  if (!dir.exists(hap_imb_dir)) missing <- c(missing, hap_imb_dir)
  if (length(missing) > 0) {
    stop("heterogeneity_consensus_heatmap 번들에 필요한 업스트림 파일이 없습니다 (수동 실행 스크립트 산출물 필요): ",
         paste(missing, collapse=", "))
  }
  rank_mat <- read.csv(rank_mat_path, stringsAsFactors=FALSE)
  cn_rank  <- read.csv(cn_rank_path, stringsAsFactors=FALSE)
  hap_imb  <- do.call(rbind, lapply(
    list.files(hap_imb_dir, pattern="\\.tsv$", full.names=TRUE),
    function(f) read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  ))
  asm_summary <- read.table(asm_summary_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)

  m1 <- setNames(as.numeric(rank_mat[rank_mat$metric=="het_hom_peak_frac", c("A","B","C","D")]), c("A","B","C","D"))
  m2 <- setNames(as.numeric(rank_mat[rank_mat$metric=="discordance_rate", c("A","B","C","D")]), c("A","B","C","D"))
  cn_vals <- setNames(cn_rank$mean_instability_rank, cn_rank$group)
  m3 <- setNames(rank(cn_vals[c("A","B","C","D")]), c("A","B","C","D"))
  hap_mean <- aggregate(frac_biased ~ group, hap_imb, mean)
  hap_vals <- setNames(hap_mean$frac_biased, hap_mean$group)
  m4 <- setNames(5 - rank(hap_vals[c("A","B","C","D")]), c("A","B","C","D"))
  asm_mean <- aggregate(n_asm_regions ~ group, asm_summary, mean)
  asm_vals <- setNames(asm_mean$n_asm_regions, asm_mean$group)
  m5 <- setNames(rank(asm_vals[c("A","B","C","D")]), c("A","B","C","D"))

  combined <- data.frame(
    method = rep(c("Hom-peak % (VAF distribution)", "Joint genotype discordance rate",
                   "Multi-copy gene family CN stability", "Haplotype-directed allelic imbalance",
                   "ASM candidate region count"), each=4),
    group  = rep(c("A","B","C","D"), 5),
    rank   = c(m1, m2, m3, m4, m5)
  )
  write.csv(combined, data_csv, row.names=FALSE)

  plot_body <- '
df$method <- factor(df$method, levels=unique(df$method))
df$group  <- factor(df$group, levels=c("A","B","C","D"))
p <- ggplot(df, aes(x=group, y=method, fill=rank)) +
  geom_tile(color="white", linewidth=1) +
  geom_text(aes(label=rank), color="black", size=5) +
  scale_fill_gradient(low="steelblue", high="firebrick", name="Heterogeneity rank",
                       limits=c(1, 4), breaks=1:4) +
  labs(title="Five independent methods point to the same heterogeneity gradient: A > C > B > D",
       x="Treatment group", y=NULL) +
  theme_bw(base_size=13) +
  theme(axis.text.y=element_text(size=11))
'
  w <- 13; h <- 9.5
  meta <- list(
    figure_title = "Heterogeneity gradient consensus across 5 independent methods",
    panel_notes = c("Single panel: tile heatmap, rows = 5 independent heterogeneity metrics, columns = groups A-D, cell = heterogeneity rank (1=most clonal, 4=most heterogeneous)."),
    input_files = c(
      "inputs/data.csv (derived table: method x group x rank, computed from 4 upstream sources)",
      "source tables: cohort/summary_compare_stats_rank_matrix.csv, cohort/copynum_instability_rank.csv, cohort/hap_imbalance/*.tsv, cohort/asm_summary.tsv"
    ),
    statistics_tables = "inputs/data.csv",
    claim_boundary = "Ranks are directionally harmonized across 5 differently-principled metrics (VAF distribution / joint genotype / depth-based copy number / haplotype phasing / allele-specific methylation) so that higher always means more heterogeneous; the ASM row is noted to be only partially independent of the others since it partly reflects the amount of phaseable heterozygosity available.",
    extra_notes = character(0)
  )
} else {
  stop("알 수 없는 slug: ", slug, " — 지원: snv_indel_burden_top30, sv_burden_single_top20, sv_burden_joint_top_hits, dmr_direction_by_chrom, trgt_manhattan, heterogeneity_consensus_heatmap")
}

# ── figure.R 작성 + 실행 (작성한 스크립트를 그대로 Rscript로 돌려서 outputs/ 생성) ──
writeLines(c("#!/usr/bin/env Rscript",
             "# Standalone, reproducible from this bundle alone.",
             "# Run from the bundle's scripts/ directory: Rscript figure.R",
             style_header, plot_body, save_footer(w, h)),
           figure_r)

oldwd <- getwd()
setwd(file.path(bundle_dir, "scripts"))
result <- tryCatch(system2("Rscript", "figure.R", stdout=TRUE, stderr=TRUE),
                    finally = setwd(oldwd))
cat(paste(result, collapse="\n"), "\n")
if (!all(file.exists(file.path(bundle_dir, "outputs", c("figure.png", "figure.pdf", "figure.svg"))))) {
  stop("figure.R 실행 후 outputs/figure.{png,pdf,svg} 중 일부가 생성되지 않았습니다: ", bundle_dir)
}

# ── metadata.yaml ────────────────────────────────────────────────────────────
write_fig_atlas_metadata(
  list(
    figure_number  = opt$`figure-number`,
    figure_slug    = slug,
    figure_title   = meta$figure_title,
    source_stem    = file.path("outputs", "figure"),
    primary_code   = "scripts/figure.R",
    panel_notes    = meta$panel_notes,
    input_files    = meta$input_files,
    output_files   = file.path("outputs", c("figure.png", "figure.pdf", "figure.svg")),
    statistics_tables = meta$statistics_tables,
    claim_boundary = meta$claim_boundary,
    extra_notes    = meta$extra_notes
  ),
  file.path(bundle_dir, "metadata", "metadata.yaml")
)

message("번들 생성 완료: ", bundle_dir)
