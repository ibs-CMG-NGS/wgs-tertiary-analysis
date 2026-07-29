#!/usr/bin/env Rscript
# ================================================================================
# ASM (Allele-Specific Methylation) 분석 스크립트
# DSS를 사용하여 동일 샘플의 hap1 vs hap2 CpG 메틸화를 비교한다.
# 유의한 ASM 영역 = 부모 유래 메틸화 패턴 차이(각인 등)를 나타낸다.
# ================================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(DSS)
  library(bsseq)
  library(ggplot2)
  library(data.table)
})

source("scripts/report_common.R")

# ================================================================================
# 커맨드 라인 인자
# ================================================================================

option_list <- list(
  make_option(c("--hap1"), type="character", default=NULL,
              help="Hap1 CpG BED 파일 경로", metavar="character"),
  make_option(c("--hap2"), type="character", default=NULL,
              help="Hap2 CpG BED 파일 경로", metavar="character"),
  make_option(c("--sample-name"), type="character", default="sample",
              help="샘플 이름 (그래프 제목 등에 사용)", metavar="character"),
  make_option(c("--output-csv"), type="character", default="asm_results.csv",
              help="ASM 결과 CSV 파일 경로 [default= %default]", metavar="character"),
  make_option(c("--output-svg-dir"), type="character", default="asm_svg",
              help="ASM 시각화 SVG 저장 디렉토리 [default= %default]", metavar="character"),
  make_option(c("--pvalue"), type="double", default=0.05,
              help="p-value 임계값 [default= %default]", metavar="number"),
  make_option(c("--min-diff"), type="double", default=0.2,
              help="최소 메틸화 차이 (0-1) [default= %default]", metavar="number"),
  make_option(c("--min-cpg"), type="integer", default=3,
              help="ASM 영역 내 최소 CpG 사이트 수 [default= %default]", metavar="integer"),
  make_option(c("--smoothing"), type="integer", default=500,
              help="DSS 평활화 윈도우 크기 (bp) [default= %default]", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$hap1) || is.null(opt$hap2)) {
  print_help(opt_parser)
  stop("--hap1 및 --hap2 인자가 필요합니다.", call.=FALSE)
}

cat("================================================================================\n")
cat("ASM 분석 시작\n")
cat("================================================================================\n")
cat(sprintf("샘플: %s\n", opt$`sample-name`))
cat(sprintf("Hap1 BED: %s\n", opt$hap1))
cat(sprintf("Hap2 BED: %s\n", opt$hap2))
cat(sprintf("p-value 임계값: %.3f\n", opt$pvalue))
cat(sprintf("최소 메틸화 차이: %.2f\n", opt$`min-diff`))
cat(sprintf("최소 CpG 사이트: %d\n", opt$`min-cpg`))
cat("================================================================================\n\n")

# BED 파일 읽기: read_cpg_bed()는 report_common.R 공유 함수(dmr_analysis.R와 동일 형식).
# on_error="stop": hap1/hap2 둘 다 필수이므로 로드 실패 시 즉시 중단.
hap1_data <- read_cpg_bed(opt$hap1, "Hap1", on_error="stop")
hap2_data <- read_cpg_bed(opt$hap2, "Hap2", on_error="stop")

# ================================================================================
# DSS BSseq 객체 생성
# ================================================================================

cat("\nDSS BSseq 객체 생성 중...\n")
bs <- makeBSseqData(
  dat = list(hap1_data, hap2_data),
  sampleNames = c("Hap1", "Hap2")
)
cat(sprintf("총 CpG 사이트 수: %s\n", format(nrow(bs), big.mark=",")))

# ================================================================================
# DML 테스트 (Hap1 vs Hap2)
# ================================================================================

cat("\nDML 테스트 수행 중 (smoothing=TRUE)...\n")
dmlTest <- DMLtest(
  BSobj       = bs,
  group1      = "Hap1",
  group2      = "Hap2",
  smoothing   = TRUE,
  smoothing.span = opt$smoothing
)
cat("DML 테스트 완료\n")

# ================================================================================
# ASM 영역 호출
# ================================================================================

cat(sprintf("\nASM 영역 호출 중 (p < %.3f, delta > %.2f)...\n",
            opt$pvalue, opt$`min-diff`))

asms <- callDMR(
  dmlTest,
  p.threshold = opt$pvalue,
  delta       = opt$`min-diff`,
  minlen      = opt$`min-cpg`,
  minCG       = opt$`min-cpg`
)

if (nrow(asms) == 0) {
  cat("\n경고: 유의한 ASM 영역이 발견되지 않았습니다.\n")
  write.csv(data.frame(), opt$`output-csv`, row.names=FALSE)
  write_placeholder_svg(opt$`output-svg-dir`, sprintf("유의한 ASM 없음 (샘플: %s)", opt$`sample-name`))
  quit(save="no", status=0)
}

cat(sprintf("발견된 ASM 영역 수: %d\n", nrow(asms)))

# ================================================================================
# 결과 테이블 정리
# ================================================================================

asm_cols <- colnames(asms)

asm_results <- tryCatch({
  data.frame(
    chr        = if ("chr"       %in% asm_cols) asms$chr       else asms[[1]],
    start      = if ("start"     %in% asm_cols) asms$start     else asms[[2]],
    end        = if ("end"       %in% asm_cols) asms$end       else asms[[3]],
    length     = if ("length"    %in% asm_cols) asms$length    else (asms[[3]] - asms[[2]]),
    nCG        = if ("nCG"       %in% asm_cols) asms$nCG       else NA,
    methy_hap1 = if ("meanMethy1" %in% asm_cols) asms$meanMethy1 else NA,
    methy_hap2 = if ("meanMethy2" %in% asm_cols) asms$meanMethy2 else NA,
    diff_methy = if ("diff.Methy" %in% asm_cols) asms$diff.Methy else NA,
    areaStat   = if ("areaStat"  %in% asm_cols) asms$areaStat  else NA
  )
}, error = function(e) {
  cat(sprintf("결과 파싱 오류: %s — 원본 저장\n", e$message))
  as.data.frame(asms)
})

if ("pval" %in% asm_cols) {
  asm_results$pvalue <- asms$pval
  asm_results$fdr    <- p.adjust(asms$pval, method="BH")
} else {
  asm_results$pvalue <- NA
  asm_results$fdr    <- NA
}

if (!all(is.na(asm_results$fdr))) {
  asm_results <- asm_results[order(asm_results$fdr), ]
} else if (!all(is.na(asm_results$pvalue))) {
  asm_results <- asm_results[order(asm_results$pvalue), ]
}

write.csv(asm_results, opt$`output-csv`, row.names=FALSE)
cat(sprintf("결과 저장 완료: %s (%d 행)\n", opt$`output-csv`, nrow(asm_results)))

# ================================================================================
# 시각화 (4종)
# ================================================================================

cat(sprintf("\nASM 시각화 중: %s\n", opt$`output-svg-dir`))
plots <- list()

has_diff  <- "diff_methy" %in% colnames(asm_results) && !all(is.na(asm_results$diff_methy))
has_fdr   <- "fdr"        %in% colnames(asm_results) && !all(is.na(asm_results$fdr))
has_length <- "length"    %in% colnames(asm_results)
has_nCG   <- "nCG"        %in% colnames(asm_results) && !all(is.na(asm_results$nCG))

sample_label <- opt$`sample-name`

# 1. 염색체별 ASM 분포
if (has_diff) {
  p1 <- ggplot(asm_results, aes(x=chr, fill=diff_methy > 0)) +
    geom_bar() +
    scale_fill_manual(values=c("steelblue", "tomato"),
                      labels=c("Hap2 > Hap1 (Hyper-Hap2)", "Hap1 > Hap2 (Hyper-Hap1)"),
                      name="메틸화 방향") +
    theme_minimal() +
    labs(title=sprintf("[%s] ASM 분포 (염색체별)", sample_label),
         x="염색체", y="ASM 영역 수") +
    theme(axis.text.x=element_text(angle=45, hjust=1))
} else {
  p1 <- ggplot(asm_results, aes(x=chr)) +
    geom_bar(fill="steelblue") +
    theme_minimal() +
    labs(title=sprintf("[%s] ASM 분포 (염색체별)", sample_label),
         x="염색체", y="ASM 영역 수") +
    theme(axis.text.x=element_text(angle=45, hjust=1))
}
plots[["01_asm_by_chrom"]] <- p1

# 2. 메틸화 차이 분포 (Hap1 - Hap2)
if (has_diff) {
  plots[["02_methylation_diff_histogram"]] <- ggplot(asm_results, aes(x=diff_methy)) +
    geom_histogram(bins=40, fill="steelblue", color="white") +
    geom_vline(xintercept=0, linetype="dashed", color="red") +
    theme_minimal() +
    labs(title=sprintf("[%s] 메틸화 차이 분포 (Hap1 - Hap2)", sample_label),
         x="메틸화 차이 (Hap1 - Hap2)",
         y="빈도")
}

# 3. Top 20 ASM (FDR 기준)
if (has_fdr && has_diff) {
  top_asm <- head(asm_results, 20)
  top_asm$region <- paste0(top_asm$chr, ":",
                           format(top_asm$start, scientific=FALSE), "-",
                           format(top_asm$end,   scientific=FALSE))
  top_asm$log10_fdr <- -log10(top_asm$fdr + 1e-300)

  plots[["03_top20_asm"]] <- ggplot(top_asm, aes(x=reorder(region, log10_fdr), y=diff_methy)) +
    geom_bar(stat="identity", aes(fill=diff_methy > 0)) +
    coord_flip() +
    scale_fill_manual(values=c("steelblue", "tomato"),
                      labels=c("Hyper-Hap2", "Hyper-Hap1"), name="") +
    theme_minimal() +
    labs(title=sprintf("[%s] Top 20 ASM 영역 (FDR 기준)", sample_label),
         x="Genomic Region", y="메틸화 차이 (Hap1 - Hap2)") +
    theme(axis.text.y=element_text(size=8))
}

# 4. ASM 크기 vs 유의성
if (has_fdr && has_diff && has_length) {
  asm_results$log10_fdr <- -log10(asm_results$fdr + 1e-300)

  if (has_nCG) {
    p4 <- ggplot(asm_results, aes(x=length, y=log10_fdr)) +
      geom_point(aes(color=diff_methy, size=nCG), alpha=0.6) +
      scale_color_gradient2(low="steelblue", mid="white", high="tomato", midpoint=0,
                            name="메틸화 차이") +
      scale_size_continuous(name="CpG 수") +
      theme_minimal() +
      labs(title=sprintf("[%s] ASM 크기와 유의성", sample_label),
           x="ASM 길이 (bp)", y="-log10(FDR)")
  } else {
    p4 <- ggplot(asm_results, aes(x=length, y=log10_fdr)) +
      geom_point(aes(color=diff_methy), alpha=0.6) +
      scale_color_gradient2(low="steelblue", mid="white", high="tomato", midpoint=0,
                            name="메틸화 차이") +
      theme_minimal() +
      labs(title=sprintf("[%s] ASM 크기와 유의성", sample_label),
           x="ASM 길이 (bp)", y="-log10(FDR)")
  }
  plots[["04_size_vs_significance"]] <- p4
}

save_svg_plots(plots, opt$`output-svg-dir`, width=12, height=8)
cat("시각화 완료\n")

# ================================================================================
# 요약 통계
# ================================================================================

cat("\n================================================================================\n")
cat("ASM 분석 요약\n")
cat("================================================================================\n")
cat(sprintf("샘플: %s\n", sample_label))
cat(sprintf("총 ASM 영역 수: %d\n", nrow(asm_results)))

if (has_diff) {
  cat(sprintf("Hyper-Hap1 ASM: %d (%.1f%%)\n",
              sum(asm_results$diff_methy > 0, na.rm=TRUE),
              100 * sum(asm_results$diff_methy > 0, na.rm=TRUE) / nrow(asm_results)))
  cat(sprintf("Hyper-Hap2 ASM: %d (%.1f%%)\n",
              sum(asm_results$diff_methy < 0, na.rm=TRUE),
              100 * sum(asm_results$diff_methy < 0, na.rm=TRUE) / nrow(asm_results)))
  cat(sprintf("평균 메틸화 차이: %.3f\n", mean(abs(asm_results$diff_methy), na.rm=TRUE)))
}
if (has_length) {
  cat(sprintf("평균 ASM 길이: %.0f bp\n", mean(asm_results$length, na.rm=TRUE)))
}
if (has_nCG) {
  cat(sprintf("평균 CpG 사이트 수: %.1f\n", mean(asm_results$nCG, na.rm=TRUE)))
}
cat("================================================================================\n")
cat(sprintf("\n분석 완료!\n결과 파일: %s\n그래프 디렉토리: %s\n",
            opt$`output-csv`, opt$`output-svg-dir`))
