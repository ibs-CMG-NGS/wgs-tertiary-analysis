#!/usr/bin/env Rscript
# =============================================================================
# Small Variant Gene-level Burden Test
# 군별 유해 변이 부담 비교 (Fisher's exact / Chi-square + BH FDR)
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(ggplot2)
  library(dplyr)
})

source("scripts/report_common.R")

# ── CLI ──────────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--input-dir",      type="character", help="canonical_impacts.tsv 파일들이 있는 디렉토리"),
  make_option("--groups",         type="character", help="그룹 정의 JSON (예: {\"ctrl\":[\"s1\",\"s2\"]})"),
  make_option("--min-impact",     type="character", default="HIGH|MODERATE", help="포함할 IMPACT 패턴 (regex)"),
  make_option("--fdr-cutoff",     type="double",    default=0.05,            help="FDR 임계값"),
  make_option("--output-csv",     type="character", help="결과 CSV 경로"),
  make_option("--output-svg-dir", type="character", help="시각화 SVG 저장 디렉토리")
)
opt <- parse_args(OptionParser(option_list=option_list))

groups      <- fromJSON(opt$`groups`)
input_dir   <- opt$`input-dir`
min_impact  <- opt$`min-impact`
fdr_cutoff  <- opt$`fdr-cutoff`
out_csv     <- opt$`output-csv`
out_svg_dir <- opt$`output-svg-dir`

# ── 데이터 로드 ──────────────────────────────────────────────────────────────
sample_group <- stack(groups)
colnames(sample_group) <- c("sample", "group")
sample_group$sample <- as.character(sample_group$sample)
sample_group$group  <- as.character(sample_group$group)

all_data <- list()
for (s in sample_group$sample) {
  path <- file.path(input_dir, paste0(s, ".canonical_impacts.tsv"))
  if (!file.exists(path)) {
    message("WARNING: 파일 없음, 스킵: ", path)
    next
  }
  df <- tryCatch(
    read.table(path, header=TRUE, sep="\t", stringsAsFactors=FALSE, quote=""),
    error = function(e) { message("읽기 오류: ", path); NULL }
  )
  if (is.null(df) || nrow(df) == 0) next
  df$sample <- s
  all_data[[s]] <- df
}

if (length(all_data) == 0) {
  message("ERROR: 로드된 데이터가 없습니다.")
  write.csv(data.frame(), out_csv, row.names=FALSE)
  write_placeholder_svg(out_svg_dir, "로드된 데이터 없음")
  quit(status=0)
}

combined <- bind_rows(all_data)

# IMPACT 필터
combined <- combined[grepl(min_impact, combined$IMPACT, ignore.case=FALSE), ]

if (nrow(combined) == 0) {
  message("WARNING: 필터 후 데이터 없음 (min_impact=", min_impact, ")")
  write.csv(data.frame(), out_csv, row.names=FALSE)
  write_placeholder_svg(out_svg_dir, paste("필터 후 데이터 없음 (min_impact=", min_impact, ")"))
  quit(status=0)
}

# 유전자 정보 컬럼 표준화 (SYMBOL 또는 GENE)
if (!"SYMBOL" %in% colnames(combined) && "GENE" %in% colnames(combined)) {
  combined$SYMBOL <- combined$GENE
}
combined <- combined[!is.na(combined$SYMBOL) & combined$SYMBOL != "." & combined$SYMBOL != "", ]

# ── gene × sample 이진 행렬 ──────────────────────────────────────────────────
gene_sample <- combined %>%
  distinct(SYMBOL, sample) %>%
  mutate(has_variant = 1L)

all_genes   <- unique(gene_sample$SYMBOL)
all_samples <- unique(sample_group$sample)

mat <- matrix(0L, nrow=length(all_genes), ncol=length(all_samples),
              dimnames=list(all_genes, all_samples))
for (i in seq_len(nrow(gene_sample))) {
  g <- gene_sample$SYMBOL[i]
  s <- gene_sample$sample[i]
  if (g %in% rownames(mat) && s %in% colnames(mat)) mat[g, s] <- 1L
}

group_names <- names(groups)

# 그룹별 이진 벡터 집계
group_mat <- sapply(group_names, function(grp) {
  samps <- intersect(groups[[grp]], colnames(mat))
  if (length(samps) == 0) return(rep(0L, nrow(mat)))
  rowSums(mat[, samps, drop=FALSE])
})
colnames(group_mat) <- group_names

# ── 통계 검정 ────────────────────────────────────────────────────────────────
n_per_group <- sapply(group_names, function(grp) length(intersect(groups[[grp]], colnames(mat))))

results <- lapply(all_genes, function(gene) {
  counts <- group_mat[gene, ]   # 그룹별 variant 있는 샘플 수
  n_grp  <- n_per_group

  test <- group_presence_test(counts, n_grp)
  eff  <- compute_group_effect_size(counts, n_grp)

  row <- data.frame(gene=gene, pvalue=test$pvalue, odds_ratio=test$odds_ratio,
                     max_prop_diff=eff$max_prop_diff, group_max=eff$group_max, group_min=eff$group_min,
                     odds_ratio_extreme=eff$odds_ratio_extreme, cramers_v=eff$cramers_v,
                     stringsAsFactors=FALSE)
  for (grp in group_names) {
    row[[paste0("n_", grp)]] <- counts[grp]
  }
  row
})

result_df <- bind_rows(results)
result_df$fdr <- p.adjust(result_df$pvalue, method="BH")
result_df <- result_df[order(result_df$fdr, na.last=TRUE), ]

write.csv(result_df, out_csv, row.names=FALSE)
message("결과 저장: ", out_csv, " (", nrow(result_df), "개 유전자)")

# ── 시각화 ───────────────────────────────────────────────────────────────────
sig_df <- result_df[!is.na(result_df$fdr) & result_df$fdr < fdr_cutoff, ]
plots <- list()

# 1. Top 30 유전자 (raw p-value 기준 — FDR은 n=3/group에서 대부분 편평해 순위에 못 씀)
#    막대 색으로 Cramer's V(효과크기) 표시
top30 <- head(result_df[order(result_df$pvalue, na.last=TRUE), ], 30)
top30$log_p <- -log10(pmax(top30$pvalue, 1e-300))
top30$gene <- factor(top30$gene, levels=rev(top30$gene))

plots[["01_top30_genes"]] <- ggplot(top30, aes(x=gene, y=log_p, fill=cramers_v)) +
  geom_col() +
  scale_fill_gradient(low="steelblue", high="firebrick", name="Cramer's V\n(효과크기)", limits=c(0,1)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
  coord_flip() +
  labs(title="Top 30 유전자 — Small Variant Burden (raw p-value 순위, 색=효과크기)",
       x=NULL, y="-log10(raw p-value)",
       subtitle=paste("점선 = raw p 0.05 (FDR 아님) | IMPACT:", min_impact)) +
  theme_bw(base_size=11)

# 2. 군별 총 burden boxplot (샘플별 HIGH/MODERATE variant 수)
burden_per_sample <- combined %>%
  group_by(sample) %>%
  summarise(n_variants=n(), .groups="drop") %>%
  left_join(sample_group, by="sample")

plots[["02_burden_boxplot"]] <- ggplot(burden_per_sample, aes(x=group, y=n_variants, color=group)) +
  geom_boxplot(outlier.shape=NA, width=0.5) +
  geom_jitter(width=0.15, size=2) +
  scale_color_brewer(palette="Set1") +
  labs(title="군별 유해 변이 부담 (샘플당)",
       x="그룹", y=paste("변이 수 (IMPACT:", min_impact, ")"),
       color="그룹") +
  theme_bw(base_size=11)

# 3. 효과크기 vs 유의성 ("volcano"): x=log2(OR_extreme), y=-log10(raw p), color=Cramer's V, size=max_prop_diff
volcano_df <- result_df[!is.na(result_df$pvalue) & !is.na(result_df$odds_ratio_extreme), ]
if (nrow(volcano_df) > 0) {
  volcano_df$log2_or <- log2(pmax(volcano_df$odds_ratio_extreme, 1e-10))
  volcano_df$log_p   <- -log10(pmax(volcano_df$pvalue, 1e-300))
  label_df <- head(volcano_df[order(volcano_df$pvalue), ], 10)

  plots[["03_effect_size_vs_significance"]] <- ggplot(volcano_df, aes(x=log2_or, y=log_p, color=cramers_v, size=max_prop_diff)) +
    geom_point(alpha=0.7) +
    scale_color_gradient(low="steelblue", high="firebrick", name="Cramer's V", limits=c(0,1)) +
    scale_size_continuous(range=c(1, 6), name="침투율 차이\n(max_prop_diff)") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
    geom_text(data=label_df, aes(label=gene), size=2.8, color="black",
              vjust=-0.8, check_overlap=TRUE, show.legend=FALSE) +
    labs(title="유전자별 효과크기 vs 유의성 (Small Variant Burden)",
         subtitle="x축: 최고/최저 군 오즈비(Haldane-Anscombe 보정, fold-change 성격) — 우측일수록 군간 차이 큼",
         x="log2(odds ratio, 최고군 vs 최저군)", y="-log10(raw p-value)") +
    theme_bw(base_size=11)
}

# 4. Bubble plot — Top 30 유전자(raw p 순위) 군별 분포 (FDR 게이트 없이, 참고용 원본 유지)
top_sig <- head(result_df[order(result_df$pvalue, na.last=TRUE), ], 30)
if (nrow(top_sig) > 0) {
  long_df <- do.call(rbind, lapply(group_names, function(grp) {
    data.frame(gene  = top_sig$gene,
               group = grp,
               count = top_sig[[paste0("n_", grp)]],
               stringsAsFactors=FALSE)
  }))
  long_df$gene <- factor(long_df$gene,
                         levels=rev(top_sig$gene[order(top_sig$pvalue)]))

  plots[["04_top30_group_bubble"]] <- ggplot(long_df, aes(x=group, y=gene, size=count, color=group)) +
    geom_point(alpha=0.7) +
    scale_size_continuous(range=c(2, 12), name="샘플 수") +
    scale_color_brewer(palette="Set1", name="그룹") +
    labs(title="Top 30 유전자(raw p 순위) — 군별 분포",
         subtitle="주의: FDR 유의 기준 아님 (n=3/group에서 FDR<0.05 도달 사실상 불가)",
         x="그룹", y="유전자") +
    theme_bw(base_size=11) +
    theme(axis.text.y=element_text(size=8))
}

# 5. Heatmap (pheatmap 있으면 사용, 없으면 ggplot 대체)
if (nrow(result_df) > 0 && ncol(mat) > 0) {
  top_genes <- head(result_df$gene[!is.na(result_df$fdr)], 40)
  hm_mat <- mat[top_genes[top_genes %in% rownames(mat)], , drop=FALSE]

  if (requireNamespace("pheatmap", quietly=TRUE) && nrow(hm_mat) > 1 && ncol(hm_mat) > 1) {
    annot_col <- data.frame(Group=sample_group$group[match(colnames(hm_mat), sample_group$sample)],
                             row.names=colnames(hm_mat))
    dir.create(out_svg_dir, showWarnings=FALSE, recursive=TRUE)
    svglite::svglite(file.path(out_svg_dir, "05_variant_heatmap.svg"), width=12, height=9)
    pheatmap::pheatmap(hm_mat,
      color=c("white", "firebrick3"),
      annotation_col=annot_col,
      cluster_rows=TRUE, cluster_cols=FALSE,
      main="Top 40 유전자 — Variant 유무 Heatmap",
      fontsize_row=7, fontsize_col=8)
    grDevices::dev.off()
  } else {
    hm_long <- as.data.frame(as.table(hm_mat))
    colnames(hm_long) <- c("gene", "sample", "has_variant")
    hm_long <- left_join(hm_long, sample_group, by="sample")
    hm_long$gene   <- factor(hm_long$gene, levels=rev(rownames(hm_mat)))
    hm_long$sample <- factor(hm_long$sample,
                             levels=sample_group$sample[order(sample_group$group)])

    plots[["05_variant_heatmap"]] <- ggplot(hm_long, aes(x=sample, y=gene, fill=factor(has_variant))) +
      geom_tile(color="grey90") +
      scale_fill_manual(values=c("0"="white", "1"="firebrick3"),
                        name="Variant", labels=c("없음","있음")) +
      facet_grid(. ~ group, scales="free_x", space="free_x") +
      labs(title="Top 40 유전자 — Variant 유무 Heatmap",
           x="샘플", y="유전자") +
      theme_bw(base_size=9) +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            axis.text.y=element_text(size=7))
  }
}

save_svg_plots(plots, out_svg_dir, width=12, height=9)
message("시각화 저장: ", out_svg_dir)
