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

# ── CLI ──────────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--input-dir",   type="character", help="canonical_impacts.tsv 파일들이 있는 디렉토리"),
  make_option("--groups",      type="character", help="그룹 정의 JSON (예: {\"ctrl\":[\"s1\",\"s2\"]})"),
  make_option("--min-impact",  type="character", default="HIGH|MODERATE", help="포함할 IMPACT 패턴 (regex)"),
  make_option("--fdr-cutoff",  type="double",    default=0.05,            help="FDR 임계값"),
  make_option("--output-csv",  type="character", help="결과 CSV 경로"),
  make_option("--output-pdf",  type="character", help="시각화 PDF 경로")
)
opt <- parse_args(OptionParser(option_list=option_list))

groups      <- fromJSON(opt$`groups`)
input_dir   <- opt$`input-dir`
min_impact  <- opt$`min-impact`
fdr_cutoff  <- opt$`fdr-cutoff`
out_csv     <- opt$`output-csv`
out_pdf     <- opt$`output-pdf`

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
  pdf(out_pdf); dev.off()
  quit(status=0)
}

combined <- bind_rows(all_data)

# IMPACT 필터
combined <- combined[grepl(min_impact, combined$IMPACT, ignore.case=FALSE), ]

if (nrow(combined) == 0) {
  message("WARNING: 필터 후 데이터 없음 (min_impact=", min_impact, ")")
  write.csv(data.frame(), out_csv, row.names=FALSE)
  pdf(out_pdf); dev.off()
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

  if (length(group_names) == 2) {
    # 2×2 Fisher's exact
    tbl <- matrix(c(counts[1], n_grp[1] - counts[1],
                    counts[2], n_grp[2] - counts[2]),
                  nrow=2)
    ft <- tryCatch(fisher.test(tbl), error=function(e) NULL)
    pval <- if (!is.null(ft)) ft$p.value else NA_real_
    or   <- if (!is.null(ft) && !is.null(ft$estimate)) unname(ft$estimate) else NA_real_
  } else {
    # χ² 검정 (3군+)
    tbl <- rbind(counts, n_grp - counts)
    ct <- tryCatch(chisq.test(tbl, simulate.p.value=(any(tbl < 5))), error=function(e) NULL)
    pval <- if (!is.null(ct)) ct$p.value else NA_real_
    or   <- NA_real_
  }

  row <- data.frame(gene=gene, pvalue=pval, odds_ratio=or, stringsAsFactors=FALSE)
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

pdf(out_pdf, width=12, height=9)

# 1. -log10(FDR) barplot (top 30)
top30 <- head(result_df[!is.na(result_df$fdr), ], 30)
top30$log_fdr <- -log10(pmax(top30$fdr, 1e-300))
top30$gene <- factor(top30$gene, levels=rev(top30$gene))

p1 <- ggplot(top30, aes(x=gene, y=log_fdr)) +
  geom_col(aes(fill=log_fdr), show.legend=FALSE) +
  scale_fill_gradient(low="steelblue", high="firebrick") +
  geom_hline(yintercept=-log10(fdr_cutoff), linetype="dashed", color="red") +
  coord_flip() +
  labs(title="Top 30 유전자 — Small Variant Burden",
       x=NULL, y="-log10(FDR)",
       subtitle=paste("FDR 임계값:", fdr_cutoff, "| IMPACT:", min_impact)) +
  theme_bw(base_size=11)
print(p1)

# 2. 군별 총 burden boxplot (샘플별 HIGH/MODERATE variant 수)
burden_per_sample <- combined %>%
  group_by(sample) %>%
  summarise(n_variants=n(), .groups="drop") %>%
  left_join(sample_group, by="sample")

p2 <- ggplot(burden_per_sample, aes(x=group, y=n_variants, color=group)) +
  geom_boxplot(outlier.shape=NA, width=0.5) +
  geom_jitter(width=0.15, size=2) +
  scale_color_brewer(palette="Set1") +
  labs(title="군별 유해 변이 부담 (샘플당)",
       x="그룹", y=paste("변이 수 (IMPACT:", min_impact, ")"),
       color="그룹") +
  theme_bw(base_size=11)
print(p2)

# 3. Bubble plot — 유의 유전자 (FDR < 임계값)
if (nrow(sig_df) > 0) {
  top_sig <- head(sig_df, 30)
  long_df <- do.call(rbind, lapply(group_names, function(grp) {
    data.frame(gene  = top_sig$gene,
               group = grp,
               count = top_sig[[paste0("n_", grp)]],
               stringsAsFactors=FALSE)
  }))
  long_df$gene <- factor(long_df$gene,
                         levels=rev(top_sig$gene[order(top_sig$fdr)]))

  p3 <- ggplot(long_df, aes(x=group, y=gene, size=count, color=group)) +
    geom_point(alpha=0.7) +
    scale_size_continuous(range=c(2, 12), name="샘플 수") +
    scale_color_brewer(palette="Set1", name="그룹") +
    labs(title=paste("유의 유전자 (FDR <", fdr_cutoff, ") — 군별 분포"),
         x="그룹", y="유전자") +
    theme_bw(base_size=11) +
    theme(axis.text.y=element_text(size=8))
  print(p3)
} else {
  plot.new()
  text(0.5, 0.5, paste("FDR <", fdr_cutoff, "인 유의 유전자 없음"), cex=1.4)
}

# 4. Heatmap (pheatmap 있으면 사용, 없으면 ggplot 대체)
if (nrow(result_df) > 0 && ncol(mat) > 0) {
  top_genes <- head(result_df$gene[!is.na(result_df$fdr)], 40)
  hm_mat <- mat[top_genes[top_genes %in% rownames(mat)], , drop=FALSE]

  if (requireNamespace("pheatmap", quietly=TRUE) && nrow(hm_mat) > 1 && ncol(hm_mat) > 1) {
    annot_col <- data.frame(Group=sample_group$group[match(colnames(hm_mat), sample_group$sample)],
                             row.names=colnames(hm_mat))
    pheatmap::pheatmap(hm_mat,
      color=c("white", "firebrick3"),
      annotation_col=annot_col,
      cluster_rows=TRUE, cluster_cols=FALSE,
      main="Top 40 유전자 — Variant 유무 Heatmap",
      fontsize_row=7, fontsize_col=8)
  } else {
    hm_long <- as.data.frame(as.table(hm_mat))
    colnames(hm_long) <- c("gene", "sample", "has_variant")
    hm_long <- left_join(hm_long, sample_group, by="sample")
    hm_long$gene   <- factor(hm_long$gene, levels=rev(rownames(hm_mat)))
    hm_long$sample <- factor(hm_long$sample,
                             levels=sample_group$sample[order(sample_group$group)])

    p4 <- ggplot(hm_long, aes(x=sample, y=gene, fill=factor(has_variant))) +
      geom_tile(color="grey90") +
      scale_fill_manual(values=c("0"="white", "1"="firebrick3"),
                        name="Variant", labels=c("없음","있음")) +
      facet_grid(. ~ group, scales="free_x", space="free_x") +
      labs(title="Top 40 유전자 — Variant 유무 Heatmap",
           x="샘플", y="유전자") +
      theme_bw(base_size=9) +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            axis.text.y=element_text(size=7))
    print(p4)
  }
}

dev.off()
message("시각화 저장: ", out_pdf)
