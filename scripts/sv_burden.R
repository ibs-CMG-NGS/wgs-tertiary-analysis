#!/usr/bin/env Rscript
# =============================================================================
# SV Gene-level Burden Test
# 군별 유전자 파괴 SV 부담 비교 (Fisher's exact / Chi-square + BH FDR)
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ── CLI ──────────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--input-dir",   type="character", help="annotated_sv.tsv 파일들이 있는 디렉토리"),
  make_option("--groups",      type="character", help="그룹 정의 JSON"),
  make_option("--sv-types",    type="character", default='["sv:cds","sv:utr"]',
              help="SV consequence 타입 JSON 배열"),
  make_option("--fdr-cutoff",  type="double",    default=0.05, help="FDR 임계값"),
  make_option("--output-csv",  type="character", help="결과 CSV 경로"),
  make_option("--output-pdf",  type="character", help="시각화 PDF 경로")
)
opt <- parse_args(OptionParser(option_list=option_list))

groups     <- fromJSON(opt$`groups`)
input_dir  <- opt$`input-dir`
sv_types   <- fromJSON(opt$`sv-types`)
fdr_cutoff <- opt$`fdr-cutoff`
out_csv    <- opt$`output-csv`
out_pdf    <- opt$`output-pdf`

group_names <- names(groups)

# ── 샘플 목록 ────────────────────────────────────────────────────────────────
sample_group <- stack(groups)
colnames(sample_group) <- c("sample", "group")
sample_group$sample <- as.character(sample_group$sample)
sample_group$group  <- as.character(sample_group$group)

# ── 데이터 로드 + BCSQ 파싱 ──────────────────────────────────────────────────
parse_bcsq <- function(bcsq_str) {
  # BCSQ 형식: "sv:cds|GENE|transcript|biotype|strand|aa_change" (콤마 구분 다중 entry)
  # → gene 이름과 consequence 타입 추출
  if (is.na(bcsq_str) || bcsq_str == "." || bcsq_str == "") return(data.frame())
  entries <- strsplit(bcsq_str, ",")[[1]]
  do.call(rbind, lapply(entries, function(e) {
    parts <- strsplit(e, "\\|")[[1]]
    if (length(parts) < 2) return(NULL)
    data.frame(sv_consequence=parts[1], gene=parts[2], stringsAsFactors=FALSE)
  }))
}

all_data <- list()
for (s in sample_group$sample) {
  path <- file.path(input_dir, paste0(s, ".annotated_sv.tsv"))
  if (!file.exists(path)) {
    message("WARNING: 파일 없음, 스킵: ", path)
    next
  }
  df <- tryCatch(
    read.table(path, header=TRUE, sep="\t", stringsAsFactors=FALSE, comment.char="",
               col.names=c("CHROM","POS","END","SVTYPE","SVLEN","BCSQ","GT")),
    error = function(e) { message("읽기 오류: ", path); NULL }
  )
  if (is.null(df) || nrow(df) == 0) next

  # BCSQ 파싱
  bcsq_rows <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    parsed <- parse_bcsq(df$BCSQ[i])
    if (nrow(parsed) == 0) return(NULL)
    cbind(df[i, c("CHROM","POS","SVTYPE","SVLEN")], parsed)
  }))

  if (is.null(bcsq_rows) || nrow(bcsq_rows) == 0) next
  bcsq_rows$sample <- s
  all_data[[s]] <- bcsq_rows
}

if (length(all_data) == 0) {
  message("ERROR: 로드된 SV 데이터가 없습니다.")
  write.csv(data.frame(), out_csv, row.names=FALSE)
  pdf(out_pdf); dev.off()
  quit(status=0)
}

combined <- bind_rows(all_data)

# sv_types 필터
combined <- combined[combined$sv_consequence %in% sv_types, ]
combined <- combined[!is.na(combined$gene) & combined$gene != "." & combined$gene != "", ]

if (nrow(combined) == 0) {
  message("WARNING: 필터 후 SV 데이터 없음 (sv_types: ", paste(sv_types, collapse=","), ")")
  write.csv(data.frame(), out_csv, row.names=FALSE)
  pdf(out_pdf); dev.off()
  quit(status=0)
}

# ── gene × sample 이진 행렬 ──────────────────────────────────────────────────
all_genes   <- unique(combined$gene)
all_samples <- unique(sample_group$sample)

mat <- matrix(0L, nrow=length(all_genes), ncol=length(all_samples),
              dimnames=list(all_genes, all_samples))
gene_sample_pairs <- combined %>% distinct(gene, sample)
for (i in seq_len(nrow(gene_sample_pairs))) {
  g <- gene_sample_pairs$gene[i]
  s <- gene_sample_pairs$sample[i]
  if (g %in% rownames(mat) && s %in% colnames(mat)) mat[g, s] <- 1L
}

n_per_group <- sapply(group_names, function(grp)
  length(intersect(groups[[grp]], colnames(mat))))

group_mat <- sapply(group_names, function(grp) {
  samps <- intersect(groups[[grp]], colnames(mat))
  if (length(samps) == 0) return(rep(0L, nrow(mat)))
  rowSums(mat[, samps, drop=FALSE])
})
colnames(group_mat) <- group_names

# ── 통계 검정 ────────────────────────────────────────────────────────────────
results <- lapply(all_genes, function(gene) {
  counts <- group_mat[gene, ]
  n_grp  <- n_per_group

  if (length(group_names) == 2) {
    tbl <- matrix(c(counts[1], n_grp[1] - counts[1],
                    counts[2], n_grp[2] - counts[2]),
                  nrow=2)
    ft <- tryCatch(fisher.test(tbl), error=function(e) NULL)
    pval <- if (!is.null(ft)) ft$p.value else NA_real_
    or   <- if (!is.null(ft) && !is.null(ft$estimate)) unname(ft$estimate) else NA_real_
  } else {
    tbl <- rbind(counts, n_grp - counts)
    ct <- tryCatch(chisq.test(tbl, simulate.p.value=(any(tbl < 5))), error=function(e) NULL)
    pval <- if (!is.null(ct)) ct$p.value else NA_real_
    or   <- NA_real_
  }

  row <- data.frame(gene=gene, pvalue=pval, odds_ratio=or, stringsAsFactors=FALSE)
  for (grp in group_names) row[[paste0("n_", grp)]] <- counts[grp]
  row
})

result_df <- bind_rows(results)
result_df$fdr <- p.adjust(result_df$pvalue, method="BH")
result_df <- result_df[order(result_df$fdr, na.last=TRUE), ]
write.csv(result_df, out_csv, row.names=FALSE)
message("결과 저장: ", out_csv)

# ── 시각화 ───────────────────────────────────────────────────────────────────
pdf(out_pdf, width=12, height=9)

# 1. SVTYPE 분포 barplot (군별)
svtype_grp <- combined %>%
  left_join(sample_group, by="sample") %>%
  count(group, SVTYPE, name="count")

p1 <- ggplot(svtype_grp, aes(x=SVTYPE, y=count, fill=group)) +
  geom_col(position="dodge") +
  scale_fill_brewer(palette="Set1") +
  labs(title="군별 SV 타입 분포",
       x="SV 타입", y="이벤트 수 (유전자 수준)", fill="그룹") +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle=30, hjust=1))
print(p1)

# 2. SV 크기 분포 violin (군별)
svlen_df <- combined %>%
  left_join(sample_group, by="sample") %>%
  mutate(SVLEN_abs=abs(as.numeric(SVLEN))) %>%
  filter(!is.na(SVLEN_abs) & SVLEN_abs > 0)

if (nrow(svlen_df) > 0) {
  p2 <- ggplot(svlen_df, aes(x=group, y=log10(SVLEN_abs), color=group)) +
    geom_violin(fill=NA) +
    geom_boxplot(width=0.15, outlier.shape=NA) +
    scale_color_brewer(palette="Set1") +
    labs(title="군별 SV 크기 분포",
         x="그룹", y="log10(SV 크기, bp)", color="그룹") +
    theme_bw(base_size=11)
  print(p2)
}

# 3. -log10(FDR) barplot (top 20)
top20 <- head(result_df[!is.na(result_df$fdr), ], 20)
top20$log_fdr <- -log10(pmax(top20$fdr, 1e-300))
top20$gene <- factor(top20$gene, levels=rev(top20$gene))

p3 <- ggplot(top20, aes(x=gene, y=log_fdr)) +
  geom_col(aes(fill=log_fdr), show.legend=FALSE) +
  scale_fill_gradient(low="steelblue", high="firebrick") +
  geom_hline(yintercept=-log10(fdr_cutoff), linetype="dashed", color="red") +
  coord_flip() +
  labs(title="Top 20 유전자 — SV Burden",
       x=NULL, y="-log10(FDR)") +
  theme_bw(base_size=11)
print(p3)

# 4. Top 20 유의 유전자 bubble plot
top_n <- min(20, nrow(result_df))
top_genes <- head(result_df$gene[!is.na(result_df$fdr)], top_n)
if (length(top_genes) > 0) {
  long_df <- do.call(rbind, lapply(group_names, function(grp) {
    data.frame(gene  = top_genes,
               group = grp,
               count = result_df[[paste0("n_", grp)]][match(top_genes, result_df$gene)],
               stringsAsFactors=FALSE)
  }))
  long_df$gene <- factor(long_df$gene, levels=rev(top_genes))

  p4 <- ggplot(long_df, aes(x=group, y=gene, size=count, color=group)) +
    geom_point(alpha=0.7) +
    scale_size_continuous(range=c(2, 12), name="샘플 수") +
    scale_color_brewer(palette="Set1", name="그룹") +
    labs(title="Top 유전자 — 군별 SV 분포",
         x="그룹", y="유전자") +
    theme_bw(base_size=11) +
    theme(axis.text.y=element_text(size=8))
  print(p4)
}

dev.off()
message("시각화 저장: ", out_pdf)
