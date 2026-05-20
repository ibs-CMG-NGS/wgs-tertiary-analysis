#!/usr/bin/env Rscript
# =============================================================================
# TRGT 반복서열 군 간 비교
# trid별 allele 길이를 Kruskal-Wallis (3군+) 또는 Wilcoxon (2군) 검정
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
  make_option("--input-dir",  type="character", help="trgt_summary.tsv 파일들이 있는 디렉토리"),
  make_option("--groups",     type="character", help="그룹 정의 JSON"),
  make_option("--fdr-cutoff", type="double",    default=0.05, help="FDR 임계값"),
  make_option("--output-csv", type="character", help="결과 CSV 경로"),
  make_option("--output-pdf", type="character", help="시각화 PDF 경로")
)
opt <- parse_args(OptionParser(option_list=option_list))

groups     <- fromJSON(opt$`groups`)
input_dir  <- opt$`input-dir`
fdr_cutoff <- opt$`fdr-cutoff`
out_csv    <- opt$`output-csv`
out_pdf    <- opt$`output-pdf`

group_names <- names(groups)

sample_group <- stack(groups)
colnames(sample_group) <- c("sample", "group")
sample_group$sample <- as.character(sample_group$sample)
sample_group$group  <- as.character(sample_group$group)

# ── 데이터 로드 ──────────────────────────────────────────────────────────────
all_data <- list()
for (s in sample_group$sample) {
  path <- file.path(input_dir, paste0(s, ".trgt_summary.tsv"))
  if (!file.exists(path)) {
    message("WARNING: 파일 없음, 스킵: ", path)
    next
  }
  df <- tryCatch(
    read.table(path, header=TRUE, sep="\t", stringsAsFactors=FALSE),
    error = function(e) { message("읽기 오류: ", path); NULL }
  )
  if (is.null(df) || nrow(df) == 0) next
  df$sample <- s
  all_data[[s]] <- df
}

if (length(all_data) == 0) {
  message("ERROR: 로드된 TRGT 데이터가 없습니다.")
  write.csv(data.frame(), out_csv, row.names=FALSE)
  pdf(out_pdf); dev.off()
  quit(status=0)
}

combined <- bind_rows(all_data) %>% left_join(sample_group, by="sample")

# allele 길이: allele1_len, allele2_len 중 더 긴 것을 사용
if (!"allele1_len" %in% colnames(combined)) {
  message("ERROR: allele1_len 컬럼 없음")
  write.csv(data.frame(), out_csv, row.names=FALSE)
  pdf(out_pdf); dev.off()
  quit(status=0)
}

combined$max_allele_len <- pmax(
  suppressWarnings(as.numeric(combined$allele1_len)),
  suppressWarnings(as.numeric(combined$allele2_len)),
  na.rm=TRUE
)

# trid별로 최소 2개 그룹에 데이터 있는 것만 분석
trid_coverage <- combined %>%
  filter(!is.na(max_allele_len)) %>%
  group_by(trid) %>%
  summarise(n_groups=n_distinct(group), .groups="drop") %>%
  filter(n_groups >= 2)

combined <- combined %>% filter(trid %in% trid_coverage$trid)

if (nrow(combined) == 0) {
  message("WARNING: 2개 이상 그룹에서 공통 TRGT locus 없음")
  write.csv(data.frame(), out_csv, row.names=FALSE)
  pdf(out_pdf); dev.off()
  quit(status=0)
}

all_trids <- unique(combined$trid)
message("분석할 TRGT locus 수: ", length(all_trids))

# ── 통계 검정 ────────────────────────────────────────────────────────────────
# 그룹별 median allele 길이 계산용 헬퍼
group_medians <- combined %>%
  filter(!is.na(max_allele_len)) %>%
  group_by(trid, group) %>%
  summarise(median_len=median(max_allele_len, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from=group, values_from=median_len,
              names_prefix="median_len_")

results <- lapply(all_trids, function(tid) {
  sub <- combined %>% filter(trid == tid, !is.na(max_allele_len))
  grp_data <- split(sub$max_allele_len, sub$group)
  grp_data  <- grp_data[sapply(grp_data, length) > 0]

  if (length(grp_data) < 2) return(NULL)

  if (length(grp_data) >= 3) {
    kt <- tryCatch(kruskal.test(max_allele_len ~ group, data=sub), error=function(e) NULL)
    pval <- if (!is.null(kt)) kt$p.value else NA_real_
    test_type <- "kruskal"
  } else {
    # 2군
    wt <- tryCatch(
      wilcox.test(grp_data[[1]], grp_data[[2]], exact=FALSE),
      error=function(e) NULL
    )
    pval <- if (!is.null(wt)) wt$p.value else NA_real_
    test_type <- "wilcoxon"
  }

  # 기본 정보
  row_info <- sub[1, c("trid"), drop=FALSE]
  for (col in c("motif","chrom","pos","expected_range")) {
    if (col %in% colnames(sub)) row_info[[col]] <- sub[[col]][1]
  }
  row_info$pvalue    <- pval
  row_info$test_type <- test_type
  row_info
})

result_df <- bind_rows(Filter(Negate(is.null), results))
result_df$fdr <- p.adjust(result_df$pvalue, method="BH")
result_df <- left_join(result_df, group_medians, by="trid")
result_df <- result_df[order(result_df$fdr, na.last=TRUE), ]

# ── pairwise Wilcoxon (유의 locus에 대해) ────────────────────────────────────
if (length(group_names) >= 3) {
  pairs <- combn(group_names, 2, simplify=FALSE)
  nom_sig <- result_df$trid[!is.na(result_df$pvalue) & result_df$pvalue < 0.1]

  if (length(nom_sig) > 0) {
    pw_results <- lapply(nom_sig, function(tid) {
      sub <- combined %>% filter(trid == tid, !is.na(max_allele_len))
      pw_row <- data.frame(trid=tid)
      for (pr in pairs) {
        g1 <- sub$max_allele_len[sub$group == pr[1]]
        g2 <- sub$max_allele_len[sub$group == pr[2]]
        if (length(g1) == 0 || length(g2) == 0) next
        wt <- tryCatch(wilcox.test(g1, g2, exact=FALSE), error=function(e) NULL)
        col_name <- paste0("pw_p_", pr[1], "_vs_", pr[2])
        pw_row[[col_name]] <- if (!is.null(wt)) wt$p.value else NA_real_
      }
      pw_row
    })
    pw_df <- bind_rows(pw_results)
    # BH correction per pairwise column
    for (col in grep("^pw_p_", colnames(pw_df), value=TRUE)) {
      pw_df[[sub("^pw_p_", "pw_fdr_", col)]] <- p.adjust(pw_df[[col]], method="BH")
    }
    result_df <- left_join(result_df, pw_df, by="trid")
  }
}

write.csv(result_df, out_csv, row.names=FALSE)
message("결과 저장: ", out_csv, " (", nrow(result_df), "개 locus)")

# ── 시각화 ───────────────────────────────────────────────────────────────────
sig_df <- result_df[!is.na(result_df$fdr) & result_df$fdr < fdr_cutoff, ]

pdf(out_pdf, width=12, height=9)

# 1. Manhattan plot
if ("chrom" %in% colnames(result_df) && "pos" %in% colnames(result_df)) {
  manhattan_df <- result_df %>%
    filter(!is.na(pvalue)) %>%
    mutate(
      chrom_num = as.numeric(sub("chr", "", chrom, ignore.case=TRUE)),
      pos_num   = suppressWarnings(as.numeric(pos)),
      log_p     = -log10(pmax(pvalue, 1e-300))
    ) %>%
    filter(!is.na(chrom_num) & !is.na(pos_num)) %>%
    arrange(chrom_num, pos_num)

  if (nrow(manhattan_df) > 0) {
    chrom_offsets <- manhattan_df %>%
      group_by(chrom_num) %>%
      summarise(max_pos=max(pos_num), .groups="drop") %>%
      arrange(chrom_num) %>%
      mutate(offset=cumsum(lag(max_pos, default=0)))

    manhattan_df <- manhattan_df %>%
      left_join(chrom_offsets[, c("chrom_num","offset")], by="chrom_num") %>%
      mutate(abs_pos=pos_num + offset)

    sig_threshold <- -log10(fdr_cutoff / nrow(manhattan_df))

    p1 <- ggplot(manhattan_df, aes(x=abs_pos, y=log_p, color=factor(chrom_num %% 2))) +
      geom_point(size=0.8, alpha=0.7) +
      geom_hline(yintercept=sig_threshold, linetype="dashed", color="red") +
      scale_color_manual(values=c("0"="steelblue","1"="navy"), guide="none") +
      labs(title="TRGT Manhattan Plot",
           x="염색체 (게놈 위치)", y="-log10(p값)") +
      theme_bw(base_size=10)
    print(p1)
  }
}

# 2. Top 20 유의 repeat violin/boxplot
top_n <- min(20, nrow(result_df))
top_trids <- head(result_df$trid[!is.na(result_df$fdr)], top_n)
if (length(top_trids) > 0) {
  plot_df <- combined %>%
    filter(trid %in% top_trids, !is.na(max_allele_len)) %>%
    mutate(trid_label=if_else("motif" %in% colnames(combined),
                              paste0(trid, "\n(", motif, ")"),
                              trid))

  trid_order <- top_trids[top_trids %in% plot_df$trid]
  plot_df$trid_label <- factor(plot_df$trid_label,
                                levels=unique(plot_df$trid_label[match(trid_order, plot_df$trid)]))

  p2 <- ggplot(plot_df, aes(x=group, y=max_allele_len, color=group)) +
    geom_boxplot(outlier.shape=NA, width=0.5) +
    geom_jitter(width=0.15, size=1.5, alpha=0.7) +
    scale_color_brewer(palette="Set1") +
    facet_wrap(~trid_label, scales="free_y", ncol=4) +
    labs(title="Top TRGT Locus — 군별 Allele 길이 분포",
         x="그룹", y="최대 Allele 길이 (bp)", color="그룹") +
    theme_bw(base_size=9) +
    theme(axis.text.x=element_text(angle=30, hjust=1),
          strip.text=element_text(size=6))
  print(p2)
}

# 3. Expanded repeat 비율 비교 (expanded 컬럼이 있을 때)
if ("status" %in% colnames(combined)) {
  exp_rate <- combined %>%
    group_by(group, sample) %>%
    summarise(
      n_total    = n(),
      n_expanded = sum(status == "expanded", na.rm=TRUE),
      pct_expanded = 100 * n_expanded / n_total,
      .groups="drop"
    )

  p3 <- ggplot(exp_rate, aes(x=group, y=pct_expanded, color=group)) +
    geom_boxplot(outlier.shape=NA, width=0.5) +
    geom_jitter(width=0.15, size=2) +
    scale_color_brewer(palette="Set1") +
    labs(title="군별 Expanded Repeat 비율",
         x="그룹", y="Expanded repeat 비율 (%)", color="그룹") +
    theme_bw(base_size=11)
  print(p3)
}

# 4. Effect size vs -log10(p) scatter
if ("median_len_" %in% paste(colnames(result_df), collapse="")) {
  med_cols <- grep("^median_len_", colnames(result_df), value=TRUE)
  if (length(med_cols) >= 2) {
    result_df$effect_size <- apply(result_df[, med_cols, drop=FALSE], 1,
                                   function(x) diff(range(x, na.rm=TRUE)))
    scatter_df <- result_df %>%
      filter(!is.na(pvalue) & !is.na(effect_size)) %>%
      mutate(log_p=-log10(pmax(pvalue, 1e-300)),
             sig=(fdr < fdr_cutoff) & !is.na(fdr))

    p4 <- ggplot(scatter_df, aes(x=effect_size, y=log_p, color=sig)) +
      geom_point(alpha=0.6, size=1.5) +
      scale_color_manual(values=c("FALSE"="grey60","TRUE"="firebrick"),
                         name=paste("FDR <", fdr_cutoff)) +
      labs(title="Effect size vs 유의성",
           x="Allele 길이 차이 (median 범위, bp)",
           y="-log10(p값)") +
      theme_bw(base_size=11)
    print(p4)
  }
}

dev.off()
message("시각화 저장: ", out_pdf)
