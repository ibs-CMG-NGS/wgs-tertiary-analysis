#!/usr/bin/env Rscript
# ================================================================================
# DMR (Differential Methylation Region) 분석 스크립트
# DSS (Dispersion Shrinkage for Sequencing data) 사용
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
# 커맨드 라인 인자 파싱
# ================================================================================

option_list <- list(
  make_option(c("--control"), type="character", default=NULL,
              help="Control 샘플 BED 파일 목록 (텍스트 파일)", metavar="character"),
  make_option(c("--experimental"), type="character", default=NULL,
              help="Experimental 샘플 BED 파일 목록 (텍스트 파일)", metavar="character"),
  make_option(c("--output-csv"), type="character", default="dmr_results.csv",
              help="DMR 결과 CSV 파일 경로 [default= %default]", metavar="character"),
  make_option(c("--output-svg-dir"), type="character", default="dmr_svg",
              help="DMR 시각화 SVG 저장 디렉토리 [default= %default]", metavar="character"),
  make_option(c("--pvalue"), type="double", default=0.05,
              help="p-value 임계값 [default= %default]", metavar="number"),
  make_option(c("--min-diff"), type="double", default=0.1,
              help="최소 메틸화 차이 (0-1) [default= %default]", metavar="number"),
  make_option(c("--min-cpg"), type="integer", default=3,
              help="DMR 내 최소 CpG 사이트 수 [default= %default]", metavar="integer"),
  make_option(c("--smoothing"), type="integer", default=500,
              help="평활화 윈도우 크기 (bp) [default= %default]", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 필수 인자 확인
if (is.null(opt$control) | is.null(opt$experimental)) {
  print_help(opt_parser)
  stop("--control 및 --experimental 인자가 필요합니다.", call.=FALSE)
}

cat("================================================================================\n")
cat("DMR 분석 시작\n")
cat("================================================================================\n")
cat(sprintf("Control 샘플: %s\n", opt$control))
cat(sprintf("Experimental 샘플: %s\n", opt$experimental))
cat(sprintf("p-value 임계값: %.3f\n", opt$pvalue))
cat(sprintf("최소 메틸화 차이: %.2f\n", opt$`min-diff`))
cat(sprintf("최소 CpG 사이트: %d\n", opt$`min-cpg`))
cat("================================================================================\n\n")

# ================================================================================
# Control 및 Experimental 샘플 로드 (read_cpg_bed()는 report_common.R 공유 함수)
# ================================================================================

# Control 샘플 목록 읽기
control_files <- readLines(opt$control)
control_files <- control_files[nchar(control_files) > 0]  # 빈 줄 제거

# Experimental 샘플 목록 읽기
experimental_files <- readLines(opt$experimental)
experimental_files <- experimental_files[nchar(experimental_files) > 0]  # 빈 줄 제거

cat(sprintf("\nControl 샘플 수: %d\n", length(control_files)))
cat(sprintf("Experimental 샘플 수: %d\n\n", length(experimental_files)))

# 데이터 로드
control_data_list <- lapply(control_files, read_cpg_bed, on_error="warn")
experimental_data_list <- lapply(experimental_files, read_cpg_bed, on_error="warn")

# NULL 제거 (로드 실패한 파일)
control_data_list <- control_data_list[!sapply(control_data_list, is.null)]
experimental_data_list <- experimental_data_list[!sapply(experimental_data_list, is.null)]

if (length(control_data_list) == 0 | length(experimental_data_list) == 0) {
  stop("유효한 데이터를 로드할 수 없습니다.", call.=FALSE)
}

# ================================================================================
# DSS 객체 생성
# ================================================================================

cat("\nDSS BSseq 객체 생성 중...\n")

# 모든 샘플 데이터를 하나의 리스트로 결합
all_data_list <- c(control_data_list, experimental_data_list)

# 샘플 이름 생성
control_names <- paste0("Control_", seq_along(control_data_list))
experimental_names <- paste0("Experimental_", seq_along(experimental_data_list))
all_sample_names <- c(control_names, experimental_names)

cat(sprintf("Control 샘플: %s\n", paste(control_names, collapse=", ")))
cat(sprintf("Experimental 샘플: %s\n", paste(experimental_names, collapse=", ")))

# BSseq 객체 생성
all_bsseq <- makeBSseqData(
  dat = all_data_list,
  sampleNames = all_sample_names
)

cat("BSseq 객체 생성 완료\n")
cat(sprintf("총 샘플 수: %d\n", length(all_sample_names)))
cat(sprintf("총 CpG 사이트 수: %s\n", format(nrow(all_bsseq), big.mark=",")))

# ================================================================================
# DML (Differential Methylation Loci) 테스트
# ================================================================================

cat("\nDML 테스트 수행 중...\n")

# DML 테스트
dmlTest <- DMLtest(
  BSobj = all_bsseq,
  group1 = control_names,
  group2 = experimental_names,
  smoothing = TRUE,
  smoothing.span = opt$smoothing
)

cat("DML 테스트 완료\n")

# ================================================================================
# DMR (Differential Methylation Region) 호출
# ================================================================================

cat(sprintf("\nDMR 호출 중 (p-value < %.3f, delta > %.2f)...\n",
            opt$pvalue, opt$`min-diff`))

dmrs <- callDMR(
  dmlTest,
  p.threshold = opt$pvalue,
  delta = opt$`min-diff`,
  minlen = opt$`min-cpg`,
  minCG = opt$`min-cpg`
)

if (nrow(dmrs) == 0) {
  cat("\n경고: 유의한 DMR이 발견되지 않았습니다.\n")
  cat("파라미터를 조정하거나 샘플 수를 늘려보세요.\n\n")

  # 빈 결과 저장
  write.csv(data.frame(), opt$`output-csv`, row.names=FALSE)
  write_placeholder_svg(opt$`output-svg-dir`, "유의한 DMR 없음")

  quit(save="no", status=0)
}

cat(sprintf("발견된 DMR 수: %d\n", nrow(dmrs)))

# DMR 객체 구조 확인 (디버깅용)
cat("\nDMR 객체 구조:\n")
cat(sprintf("  클래스: %s\n", class(dmrs)))
cat(sprintf("  컬럼명: %s\n", paste(colnames(dmrs), collapse=", ")))
if (nrow(dmrs) > 0) {
  cat("  첫 번째 DMR:\n")
  print(head(dmrs, 1))
}

# ================================================================================
# DMR 결과 저장
# ================================================================================

cat(sprintf("\nDMR 결과 저장 중: %s\n", opt$`output-csv`))

# DSS callDMR 결과는 data.frame이지만 컬럼명이 다를 수 있음
# 일반적인 컬럼: chr, start, end, length, nCG, meanMethy1, meanMethy2, diff.Methy, areaStat
# 컬럼명 확인 후 적절히 매핑
dmr_cols <- colnames(dmrs)

# 결과 테이블 작성
dmr_results <- tryCatch({
  data.frame(
    chr = if("chr" %in% dmr_cols) dmrs$chr else dmrs[[1]],
    start = if("start" %in% dmr_cols) dmrs$start else dmrs[[2]],
    end = if("end" %in% dmr_cols) dmrs$end else dmrs[[3]],
    length = if("length" %in% dmr_cols) dmrs$length else (dmrs[[3]] - dmrs[[2]]),
    nCG = if("nCG" %in% dmr_cols) dmrs$nCG else NA,
    meanMethy1 = if("meanMethy1" %in% dmr_cols) dmrs$meanMethy1 else NA,
    meanMethy2 = if("meanMethy2" %in% dmr_cols) dmrs$meanMethy2 else NA,
    diff_methy = if("diff.Methy" %in% dmr_cols) dmrs$diff.Methy else NA,
    areaStat = if("areaStat" %in% dmr_cols) dmrs$areaStat else NA
  )
}, error = function(e) {
  cat(sprintf("\n에러 발생: %s\n", e$message))
  cat("DMR 객체를 그대로 저장합니다.\n")
  as.data.frame(dmrs)
})

# p-value: DSS::callDMR()은 영역(region) 단위 p-value를 제공하지 않는다 (areaStat만 반환).
# 시도해본 areaStat/sqrt(nCG) 기반 결합검정(Stouffer 근사, 개별 CpG 독립 가정)은 인접 CpG의
# 공간적 상관관계(smoothing으로 인해 더 심함)를 무시해 심하게 반코저버티브였다 (전체 영역의
# 99%+ 가 FDR<0.05로 나오는 등 비현실적) — 실사용 불가로 판단해 폐기.
# 엄밀한 영역 단위 검정을 하려면 permutation test(그룹 라벨을 섞어 areaStat의 귀무분포를
# 영역 길이/nCG별로 추정)가 필요하며 이는 별도 작업으로 남긴다.
# 현재는 정직하게 NA로 두고, callDMR() 자체에 이미 적용된 개별 CpG p-threshold + delta +
# minCG 필터를 "선별 기준"으로, areaStat을 "효과크기 순위"로만 사용한다.
dmr_results$pvalue <- NA
dmr_results$fdr <- NA

# FDR로 정렬 (FDR이 존재하는 경우), 없으면 |areaStat| 내림차순(효과크기 순위)으로 정렬
if ("fdr" %in% colnames(dmr_results) && !all(is.na(dmr_results$fdr))) {
  dmr_results <- dmr_results[order(dmr_results$fdr), ]
} else if ("pvalue" %in% colnames(dmr_results) && !all(is.na(dmr_results$pvalue))) {
  dmr_results <- dmr_results[order(dmr_results$pvalue), ]
} else if ("areaStat" %in% colnames(dmr_results)) {
  dmr_results <- dmr_results[order(-abs(dmr_results$areaStat)), ]
}

write.csv(dmr_results, opt$`output-csv`, row.names=FALSE)

cat("DMR 결과 저장 완료\n")
cat(sprintf("  저장된 DMR 수: %d\n", nrow(dmr_results)))
cat(sprintf("  컬럼: %s\n", paste(colnames(dmr_results), collapse=", ")))

# ================================================================================
# 시각화
# ================================================================================

cat(sprintf("\nDMR 시각화 중: %s\n", opt$`output-svg-dir`))

plots <- list()

# diff_methy가 있는지 확인
has_diff_methy <- "diff_methy" %in% colnames(dmr_results) && !all(is.na(dmr_results$diff_methy))

# 1. DMR 분포 (염색체별)
if (has_diff_methy) {
  p1 <- ggplot(dmr_results, aes(x=chr, fill=diff_methy > 0)) +
    geom_bar() +
    theme_minimal() +
    labs(title="DMR 분포 (염색체별)",
         x="염색체", y="DMR 수",
         fill="메틸화 방향")
} else {
  p1 <- ggplot(dmr_results, aes(x=chr)) +
    geom_bar() +
    theme_minimal() +
    labs(title="DMR 분포 (염색체별)",
         x="염색체", y="DMR 수")
}
p1 <- p1 + theme(axis.text.x = element_text(angle=45, hjust=1))
plots[["01_dmr_by_chrom"]] <- p1

# 2. 메틸화 차이 분포
if (has_diff_methy) {
  plots[["02_methylation_diff_histogram"]] <- ggplot(dmr_results, aes(x=diff_methy)) +
    geom_histogram(bins=50, fill="steelblue", color="black") +
    geom_vline(xintercept=0, linetype="dashed", color="red") +
    theme_minimal() +
    labs(title="메틸화 차이 분포",
         x="메틸화 차이 (Experimental - Control)",
         y="빈도")
}

# 3. Volcano plot
has_fdr <- "fdr" %in% colnames(dmr_results) && !all(is.na(dmr_results$fdr))
if (has_diff_methy && has_fdr) {
  dmr_results$log10_fdr <- -log10(dmr_results$fdr + 1e-300)  # 0 방지
  plots[["03_volcano"]] <- ggplot(dmr_results, aes(x=diff_methy, y=log10_fdr)) +
    geom_point(aes(color=abs(diff_methy) > opt$`min-diff` & fdr < opt$pvalue),
               alpha=0.6) +
    geom_hline(yintercept=-log10(opt$pvalue), linetype="dashed", color="red") +
    geom_vline(xintercept=c(-opt$`min-diff`, opt$`min-diff`),
               linetype="dashed", color="blue") +
    theme_minimal() +
    labs(title="DMR Volcano Plot",
         x="메틸화 차이",
         y="-log10(FDR)") +
    scale_color_manual(values=c("grey", "red"),
                       labels=c("Not significant", "Significant"),
                       name="")
}

# 4. Top 20 DMR (FDR 기준)
if (has_fdr && has_diff_methy) {
  top_dmrs <- head(dmr_results, 20)
  top_dmrs$region <- paste(top_dmrs$chr, ":",
                           format(top_dmrs$start, scientific=FALSE), "-",
                           format(top_dmrs$end, scientific=FALSE), sep="")
  top_dmrs$log10_fdr <- -log10(top_dmrs$fdr + 1e-300)

  plots[["04_top20_dmr"]] <- ggplot(top_dmrs, aes(x=reorder(region, log10_fdr), y=diff_methy)) +
    geom_bar(stat="identity", aes(fill=diff_methy > 0)) +
    coord_flip() +
    theme_minimal() +
    labs(title="Top 20 DMR (FDR 기준)",
         x="Genomic Region",
         y="메틸화 차이",
         fill="메틸화 방향") +
    scale_fill_manual(values=c("blue", "red"),
                      labels=c("Hypo", "Hyper"),
                      name="") +
    theme(axis.text.y = element_text(size=8))
}

# 5. DMR 크기 vs 유의성
if (has_fdr && has_diff_methy && "length" %in% colnames(dmr_results)) {
  has_nCG <- "nCG" %in% colnames(dmr_results) && !all(is.na(dmr_results$nCG))

  dmr_results$log10_fdr <- -log10(dmr_results$fdr + 1e-300)

  if (has_nCG) {
    p5 <- ggplot(dmr_results, aes(x=length, y=log10_fdr)) +
      geom_point(aes(color=diff_methy, size=nCG), alpha=0.6) +
      scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0,
                            name="메틸화 차이") +
      scale_size_continuous(name="CpG 수") +
      theme_minimal() +
      labs(title="DMR 크기와 유의성",
           x="DMR 길이 (bp)",
           y="-log10(FDR)")
  } else {
    p5 <- ggplot(dmr_results, aes(x=length, y=log10_fdr)) +
      geom_point(aes(color=diff_methy), alpha=0.6) +
      scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0,
                            name="메틸화 차이") +
      theme_minimal() +
      labs(title="DMR 크기와 유의성",
           x="DMR 길이 (bp)",
           y="-log10(FDR)")
  }
  plots[["05_size_vs_significance"]] <- p5
}

save_svg_plots(plots, opt$`output-svg-dir`, width=12, height=8)

cat("시각화 완료\n")

# ================================================================================
# 요약 통계
# ================================================================================

cat("\n================================================================================\n")
cat("DMR 분석 요약\n")
cat("================================================================================\n")
cat(sprintf("총 DMR 수: %d\n", nrow(dmr_results)))

if (has_diff_methy) {
  cat(sprintf("Hyper-메틸화 DMR: %d (%.1f%%)\n",
              sum(dmr_results$diff_methy > 0, na.rm=TRUE),
              100 * sum(dmr_results$diff_methy > 0, na.rm=TRUE) / nrow(dmr_results)))
  cat(sprintf("Hypo-메틸화 DMR: %d (%.1f%%)\n",
              sum(dmr_results$diff_methy < 0, na.rm=TRUE),
              100 * sum(dmr_results$diff_methy < 0, na.rm=TRUE) / nrow(dmr_results)))
  cat(sprintf("평균 메틸화 차이: %.3f\n", mean(abs(dmr_results$diff_methy), na.rm=TRUE)))
}
if ("length" %in% colnames(dmr_results)) {
  cat(sprintf("평균 DMR 길이: %.0f bp\n", mean(dmr_results$length, na.rm=TRUE)))
}
if ("nCG" %in% colnames(dmr_results)) {
  cat(sprintf("평균 CpG 사이트 수: %.1f\n", mean(dmr_results$nCG, na.rm=TRUE)))
}
cat("================================================================================\n")

cat("\n분석 완료!\n")
cat(sprintf("결과 파일: %s\n", opt$`output-csv`))
cat(sprintf("그래프 디렉토리: %s\n", opt$`output-svg-dir`))
