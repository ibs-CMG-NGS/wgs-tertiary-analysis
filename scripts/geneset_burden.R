#!/usr/bin/env Rscript
# =============================================================================
# 유전자셋(GO pathway) 단위 집계 검정
# 개별 유전자 수백~천 개 대신, 생물학적으로 관련된 유전자 묶음(산화스트레스 반응,
# DNA손상복구 등) 전체를 하나의 단위로 검정 — 검정 개수를 유전자셋 개수(소수)로 줄여
# 다중검정 부담을 크게 낮춘다.
# 입력: --mode sv (annotated_sv.tsv, BCSQ 파싱) 또는 --mode variant (canonical_impacts.tsv)
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(ggplot2)
  library(dplyr)
})

source("scripts/report_common.R")

option_list <- list(
  make_option("--mode",           type="character", help="sv 또는 variant"),
  make_option("--input-dir",      type="character", help="입력 TSV들이 있는 디렉토리"),
  make_option("--groups",         type="character", help="그룹 정의 JSON"),
  make_option("--genesets-json",  type="character", help="유전자셋 JSON (fetch_go_genesets.py 출력)"),
  make_option("--sv-types",       type="character", default='["sv:cds","sv:utr"]',
              help="(mode=sv) SV consequence 타입 JSON 배열"),
  make_option("--max-svlen",      type="double", default=1000000, help="(mode=sv) 최대 |SVLEN|"),
  make_option("--min-impact",     type="character", default="HIGH|MODERATE", help="(mode=variant) IMPACT 정규식"),
  make_option("--fdr-cutoff",     type="double", default=0.05),
  make_option("--output-csv",     type="character"),
  make_option("--output-svg-dir", type="character", help="시각화 SVG 저장 디렉토리")
)
opt <- parse_args(OptionParser(option_list=option_list))

groups    <- fromJSON(opt$`groups`)
genesets  <- fromJSON(opt$`genesets-json`)
mode      <- opt$mode
out_svg_dir <- opt$`output-svg-dir`

sample_group <- stack(groups)
colnames(sample_group) <- c("sample", "group")
sample_group$sample <- as.character(sample_group$sample)
sample_group$group  <- as.character(sample_group$group)
group_names <- names(groups)
n_per_group <- sapply(group_names, function(grp) length(intersect(groups[[grp]], sample_group$sample)))

# ── 샘플별 유전자 이벤트 로드 (gene, sample) 쌍 + 이벤트 카운트 ────────────────
# parse_bcsq()는 report_common.R 공유 함수
load_sv_sample <- function(path, max_svlen, sv_types) {
  df <- tryCatch(
    read.table(path, header=TRUE, sep="\t", stringsAsFactors=FALSE, comment.char="",
               col.names=c("CHROM","POS","END","SVTYPE","SVLEN","BCSQ","GT")),
    error=function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(data.frame(gene=character(0)))
  df <- df[abs(as.numeric(df$SVLEN)) <= max_svlen, ]
  rows <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    parsed <- parse_bcsq(df$BCSQ[i])
    if (nrow(parsed) == 0) return(NULL)
    parsed
  }))
  if (is.null(rows)) return(data.frame(gene=character(0)))
  rows <- rows[rows$sv_consequence %in% sv_types, ]
  data.frame(gene=rows$gene[!is.na(rows$gene) & rows$gene != "."])
}

load_variant_sample <- function(path, min_impact) {
  df <- tryCatch(read.table(path, header=TRUE, sep="\t", stringsAsFactors=FALSE, comment.char=""),
                 error=function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(data.frame(gene=character(0)))
  df <- df[grepl(min_impact, df$IMPACT) & !is.na(df$SYMBOL) & df$SYMBOL != ".", ]
  data.frame(gene=df$SYMBOL)
}

sv_types <- fromJSON(opt$`sv-types`)
all_events <- list()
for (s in sample_group$sample) {
  fname <- if (mode == "sv") paste0(s, ".annotated_sv.tsv") else paste0(s, ".canonical_impacts.tsv")
  path <- file.path(opt$`input-dir`, fname)
  if (!file.exists(path)) { message("WARNING: 파일 없음, 스킵: ", path); next }
  ev <- if (mode == "sv") load_sv_sample(path, opt$`max-svlen`, sv_types) else load_variant_sample(path, opt$`min-impact`)
  if (nrow(ev) == 0) next
  ev$sample <- s
  all_events[[s]] <- ev
}
combined <- bind_rows(all_events)

# ── 유전자셋별 집계: 샘플당 (이벤트 있는 유전자 수>=1 여부, 총 이벤트 수) ───────
results <- lapply(names(genesets), function(gs_name) {
  gs_genes <- genesets[[gs_name]]
  gs_events <- combined[combined$gene %in% gs_genes, ]

  per_sample_count <- sapply(sample_group$sample, function(s) sum(gs_events$sample == s))
  names(per_sample_count) <- sample_group$sample
  per_sample_presence <- as.integer(per_sample_count > 0)
  names(per_sample_presence) <- sample_group$sample

  # (a) presence 기반 카이제곱/피셔 (기존 유전자 burden과 동일한 검정 방식)
  presence_by_group <- sapply(group_names, function(g) sum(per_sample_presence[groups[[g]]], na.rm=TRUE))
  p_presence <- group_presence_test(presence_by_group, n_per_group)$pvalue

  # (b) 총 이벤트 수 기반 Kruskal-Wallis (연속값, 검정력 더 좋음)
  count_group <- factor(sample_group$group[match(names(per_sample_count), sample_group$sample)])
  kt <- tryCatch(kruskal.test(per_sample_count, count_group), error=function(e) NULL)
  p_count <- if (!is.null(kt)) kt$p.value else NA_real_

  # 효과크기: 침투율 차이, Cramer's V (presence 기준)
  eff <- compute_group_effect_size(presence_by_group, n_per_group)

  row <- data.frame(geneset=gs_name, n_genes=length(gs_genes),
                     pvalue_presence=p_presence, pvalue_count=p_count,
                     max_prop_diff=eff$max_prop_diff, cramers_v=eff$cramers_v,
                     stringsAsFactors=FALSE)
  for (g in group_names) {
    row[[paste0("n_presence_", g)]] <- presence_by_group[g]
    row[[paste0("mean_count_", g)]] <- mean(per_sample_count[groups[[g]]], na.rm=TRUE)
  }
  attr(row, "per_sample_count") <- per_sample_count
  row
})

result_df <- bind_rows(results)
result_df$fdr_presence <- p.adjust(result_df$pvalue_presence, method="BH")
result_df$fdr_count    <- p.adjust(result_df$pvalue_count, method="BH")
write.csv(result_df, opt$`output-csv`, row.names=FALSE)
message("결과 저장: ", opt$`output-csv`)
print(result_df)

# ── 시각화 ───────────────────────────────────────────────────────────────────
plots <- list()
for (i in seq_along(results)) {
  gs_name <- names(genesets)[i]
  per_sample_count <- attr(results[[i]], "per_sample_count")
  plot_df <- data.frame(sample=names(per_sample_count), count=per_sample_count) %>%
    left_join(sample_group, by="sample")
  kw_p <- result_df$pvalue_count[result_df$geneset == gs_name]
  p <- ggplot(plot_df, aes(x=group, y=count, color=group)) +
    geom_boxplot(outlier.shape=NA, width=0.5) +
    geom_jitter(width=0.15, size=2) +
    scale_color_brewer(palette="Set1") +
    labs(title=paste0("유전자셋: ", gs_name, " (후보 ", length(genesets[[gs_name]]), "개 유전자 중 이벤트)"),
         subtitle=sprintf("샘플당 이벤트 수, Kruskal-Wallis p=%.4f", kw_p),
         x="그룹", y="유전자셋 내 이벤트 수 (샘플당)") +
    theme_bw(base_size=12)
  plots[[sprintf("%02d_%s", i, gsub("[^A-Za-z0-9_]+", "_", gs_name))]] <- p
}
save_svg_plots(plots, out_svg_dir, width=10, height=7)
message("시각화 저장: ", out_svg_dir)
