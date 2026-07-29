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

source("scripts/report_common.R")

# ── CLI ──────────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--input-dir",   type="character", help="annotated_sv.tsv 파일들이 있는 디렉토리"),
  make_option("--groups",      type="character", help="그룹 정의 JSON"),
  make_option("--sv-types",    type="character", default='["sv:cds","sv:utr"]',
              help="SV consequence 타입 JSON 배열"),
  make_option("--max-svlen",   type="double",    default=1000000,
              help="최대 |SVLEN|(bp). 이보다 큰 SV는 burden test에서 제외 (염색체 단위 초대형 콜이 다수 유전자에 중복 카운트되는 것 방지) [default= %default]"),
  make_option("--exclude-chroms", type="character", default="",
              help="제외할 염색체(콤마구분, 예: chrX,chrY) — 군간 성별 구성이 다르면 성염색체 커버리지 차이가 SV로 오콜링되어 burden test를 교란함"),
  make_option("--exclude-gene-prefixes", type="character", default="",
              help="제외할 유전자명 접두어(콤마구분, 예: Vmn1r,Mup,Pramel,Ear,Gm) — 마우스에서 개체간 카피수 다형성이 원래 큰 것으로 알려진 다중카피 유전자군을 제외해 처치효과와 개체 유전배경차를 혼동하지 않기 위함"),
  make_option("--fdr-cutoff",     type="double",    default=0.05, help="FDR 임계값"),
  make_option("--output-csv",     type="character", help="결과 CSV 경로"),
  make_option("--output-svg-dir", type="character", help="시각화 SVG 저장 디렉토리")
)
opt <- parse_args(OptionParser(option_list=option_list))

groups     <- fromJSON(opt$`groups`)
input_dir  <- opt$`input-dir`
sv_types   <- fromJSON(opt$`sv-types`)
max_svlen  <- opt$`max-svlen`
exclude_chroms <- if (nchar(opt$`exclude-chroms`) > 0) strsplit(opt$`exclude-chroms`, ",")[[1]] else character(0)
exclude_gene_prefixes <- if (nchar(opt$`exclude-gene-prefixes`) > 0) strsplit(opt$`exclude-gene-prefixes`, ",")[[1]] else character(0)
fdr_cutoff  <- opt$`fdr-cutoff`
out_csv     <- opt$`output-csv`
out_svg_dir <- opt$`output-svg-dir`

group_names <- names(groups)

# ── 샘플 목록 ────────────────────────────────────────────────────────────────
sample_group <- stack(groups)
colnames(sample_group) <- c("sample", "group")
sample_group$sample <- as.character(sample_group$sample)
sample_group$group  <- as.character(sample_group$group)

# ── 데이터 로드 + BCSQ 파싱 (parse_bcsq()는 report_common.R 공유 함수) ─────────
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
  if (length(exclude_chroms) > 0) df <- df[!(df$CHROM %in% exclude_chroms), ]
  if (nrow(df) == 0) next

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
  write_placeholder_svg(out_svg_dir, "로드된 SV 데이터 없음")
  quit(status=0)
}

combined <- bind_rows(all_data)

# sv_types 필터
combined <- combined[combined$sv_consequence %in% sv_types, ]
combined <- combined[!is.na(combined$gene) & combined$gene != "." & combined$gene != "", ]

# 초대형 SV(염색체 단위 콜) 제외: 다수 유전자에 걸쳐 중복 카운트되어
# burden test를 왜곡하는 것을 방지 (예: 9.8Mb 단일 duplication이 490개 유전자로 집계된 사례)
n_before_svlen <- nrow(combined)
combined <- combined[abs(as.numeric(combined$SVLEN)) <= max_svlen, ]
n_removed <- n_before_svlen - nrow(combined)
if (n_removed > 0) {
  message(sprintf("최대 SVLEN(%.0fbp) 초과로 %d개 SV 레코드 제외됨", max_svlen, n_removed))
}

# 개체간 카피수 다형성이 원래 큰 것으로 알려진 다중카피 유전자군 제외
# (Vmn1r/Mup/Pramel/Ear/Gm-prefix 등 — 처치효과와 개체 유전배경차 혼동 방지)
if (length(exclude_gene_prefixes) > 0) {
  pattern <- paste0("^(", paste(exclude_gene_prefixes, collapse="|"), ")")
  n_before_genefilter <- nrow(combined)
  combined <- combined[!grepl(pattern, combined$gene), ]
  n_removed_gene <- n_before_genefilter - nrow(combined)
  message(sprintf("다중카피 유전자군 접두어(%s) 매칭으로 %d개 레코드 제외됨",
                  paste(exclude_gene_prefixes, collapse=","), n_removed_gene))
}

if (nrow(combined) == 0) {
  message("WARNING: 필터 후 SV 데이터 없음 (sv_types: ", paste(sv_types, collapse=","), ")")
  write.csv(data.frame(), out_csv, row.names=FALSE)
  write_placeholder_svg(out_svg_dir, paste("필터 후 SV 데이터 없음 (sv_types:", paste(sv_types, collapse=","), ")"))
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

  test <- group_presence_test(counts, n_grp)
  eff  <- compute_group_effect_size(counts, n_grp)

  row <- data.frame(gene=gene, pvalue=test$pvalue, odds_ratio=test$odds_ratio,
                     max_prop_diff=eff$max_prop_diff, group_max=eff$group_max, group_min=eff$group_min,
                     odds_ratio_extreme=eff$odds_ratio_extreme, cramers_v=eff$cramers_v,
                     stringsAsFactors=FALSE)
  for (grp in group_names) row[[paste0("n_", grp)]] <- counts[grp]
  row
})

result_df <- bind_rows(results)
result_df$fdr <- p.adjust(result_df$pvalue, method="BH")
result_df <- result_df[order(result_df$fdr, na.last=TRUE), ]
write.csv(result_df, out_csv, row.names=FALSE)
message("결과 저장: ", out_csv)

# ── 시각화 ───────────────────────────────────────────────────────────────────
plots <- list()

# 1. SVTYPE 분포 barplot (군별)
svtype_grp <- combined %>%
  left_join(sample_group, by="sample") %>%
  count(group, SVTYPE, name="count")

plots[["01_svtype_distribution"]] <- ggplot(svtype_grp, aes(x=SVTYPE, y=count, fill=group)) +
  geom_col(position="dodge") +
  scale_fill_brewer(palette="Set1") +
  labs(title="군별 SV 타입 분포",
       x="SV 타입", y="이벤트 수 (유전자 수준)", fill="그룹") +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle=30, hjust=1))

# 2. SV 크기 분포 violin (군별)
svlen_df <- combined %>%
  left_join(sample_group, by="sample") %>%
  mutate(SVLEN_abs=abs(as.numeric(SVLEN))) %>%
  filter(!is.na(SVLEN_abs) & SVLEN_abs > 0)

if (nrow(svlen_df) > 0) {
  plots[["02_svlen_distribution"]] <- ggplot(svlen_df, aes(x=group, y=log10(SVLEN_abs), color=group)) +
    geom_violin(fill=NA) +
    geom_boxplot(width=0.15, outlier.shape=NA) +
    scale_color_brewer(palette="Set1") +
    labs(title="군별 SV 크기 분포",
         x="그룹", y="log10(SV 크기, bp)", color="그룹") +
    theme_bw(base_size=11)
}

# 3. Top 20 유전자 (raw p-value 기준 — FDR이 n=3/group에서는 대부분 편평해서 raw p로 순위)
#    막대 색으로 Cramer's V(효과크기)도 함께 표시
top20 <- head(result_df[order(result_df$pvalue, na.last=TRUE), ], 20)
top20$log_p <- -log10(pmax(top20$pvalue, 1e-300))
top20$gene <- factor(top20$gene, levels=rev(top20$gene))

plots[["03_top20_genes"]] <- ggplot(top20, aes(x=gene, y=log_p, fill=cramers_v)) +
  geom_col() +
  scale_fill_gradient(low="steelblue", high="firebrick", name="Cramer's V\n(효과크기)", limits=c(0,1)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
  coord_flip() +
  labs(title="Top 20 유전자 — SV Burden (raw p-value 순위, 색=효과크기)",
       subtitle="점선 = raw p 0.05 (FDR 아님 — n=3/group에서는 FDR 임계값 도달이 사실상 불가능)",
       x=NULL, y="-log10(raw p-value)") +
  theme_bw(base_size=11)

# 4. 효과크기 vs 유의성 ("volcano"): x=log2(OR_extreme), y=-log10(raw p), color=Cramer's V, size=max_prop_diff
volcano_df <- result_df[!is.na(result_df$pvalue) & !is.na(result_df$odds_ratio_extreme), ]
if (nrow(volcano_df) > 0) {
  volcano_df$log2_or <- log2(pmax(volcano_df$odds_ratio_extreme, 1e-10))
  volcano_df$log_p   <- -log10(pmax(volcano_df$pvalue, 1e-300))
  label_df <- head(volcano_df[order(volcano_df$pvalue), ], 10)

  plots[["04_effect_size_vs_significance"]] <- ggplot(volcano_df, aes(x=log2_or, y=log_p, color=cramers_v, size=max_prop_diff)) +
    geom_point(alpha=0.7) +
    scale_color_gradient(low="steelblue", high="firebrick", name="Cramer's V", limits=c(0,1)) +
    scale_size_continuous(range=c(1, 6), name="침투율 차이\n(max_prop_diff)") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
    geom_text(data=label_df, aes(label=gene), size=2.8, color="black",
              vjust=-0.8, check_overlap=TRUE, show.legend=FALSE) +
    labs(title="유전자별 효과크기 vs 유의성 (SV burden)",
         subtitle="x축: 최고/최저 군 오즈비(Haldane-Anscombe 보정, fold-change 성격) — 우측일수록 군간 차이 큼",
         x="log2(odds ratio, 최고군 vs 최저군)", y="-log10(raw p-value)") +
    theme_bw(base_size=11)
}

# 5. Top 20 유전자 — 군별 SV 보유 샘플 수 bubble plot (raw p 기준 top20 유전자 대상)
top_genes <- as.character(top20$gene)
if (length(top_genes) > 0) {
  long_df <- do.call(rbind, lapply(group_names, function(grp) {
    data.frame(gene  = top_genes,
               group = grp,
               count = result_df[[paste0("n_", grp)]][match(top_genes, result_df$gene)],
               stringsAsFactors=FALSE)
  }))
  long_df$gene <- factor(long_df$gene, levels=rev(top_genes))

  plots[["05_top20_group_bubble"]] <- ggplot(long_df, aes(x=group, y=gene, size=count, color=group)) +
    geom_point(alpha=0.7) +
    scale_size_continuous(range=c(2, 12), name="샘플 수") +
    scale_color_brewer(palette="Set1", name="그룹") +
    labs(title="Top 유전자(raw p 순위) — 군별 SV 분포",
         x="그룹", y="유전자") +
    theme_bw(base_size=11) +
    theme(axis.text.y=element_text(size=8))
}

save_svg_plots(plots, out_svg_dir, width=12, height=9)
message("시각화 저장: ", out_svg_dir)
