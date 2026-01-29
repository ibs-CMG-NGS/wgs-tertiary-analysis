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
  make_option(c("--output-pdf"), type="character", default="dmr_plots.pdf",
              help="DMR 시각화 PDF 파일 경로 [default= %default]", metavar="character"),
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
# BED 파일 읽기 함수
# ================================================================================

read_cpg_bed <- function(bed_file) {
  cat(sprintf("파일 로딩 중: %s\n", bed_file))
  
  # PacBio CpG BED 형식: chr, start, end, methylation_level, coverage
  # 실제 형식에 맞게 조정 필요
  tryCatch({
    dt <- fread(bed_file, header=FALSE)
    
    # 컬럼명 설정 (PacBio 형식 가정)
    # 실제 BED 파일 형식에 따라 수정 필요
    if (ncol(dt) >= 5) {
      colnames(dt) <- c("chr", "start", "end", "methylation", "coverage")[1:ncol(dt)]
    } else if (ncol(dt) >= 4) {
      colnames(dt) <- c("chr", "pos", "coverage", "methylation")[1:ncol(dt)]
      dt$start <- dt$pos
      dt$end <- dt$pos + 1
    }
    
    # DSS 형식으로 변환
    dss_data <- data.frame(
      chr = dt$chr,
      pos = dt$start,
      N = as.integer(dt$coverage),
      X = as.integer(round(dt$coverage * dt$methylation))
    )
    
    # 필터링: coverage > 0
    dss_data <- dss_data[dss_data$N > 0, ]
    
    return(dss_data)
  }, error = function(e) {
    cat(sprintf("오류: %s 파일을 읽을 수 없습니다. %s\n", bed_file, e$message))
    return(NULL)
  })
}

# ================================================================================
# Control 및 Experimental 샘플 로드
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
control_data_list <- lapply(control_files, read_cpg_bed)
experimental_data_list <- lapply(experimental_files, read_cpg_bed)

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

# Control BSseq 객체
control_bsseq_list <- lapply(seq_along(control_data_list), function(i) {
  data <- control_data_list[[i]]
  makeBSseqData(
    dat = list(data),
    sampleNames = paste0("Control_", i)
  )
})

# Experimental BSseq 객체
experimental_bsseq_list <- lapply(seq_along(experimental_data_list), function(i) {
  data <- experimental_data_list[[i]]
  makeBSseqData(
    dat = list(data),
    sampleNames = paste0("Experimental_", i)
  )
})

# BSseq 객체 결합
all_bsseq <- combine(c(control_bsseq_list, experimental_bsseq_list))

cat("BSseq 객체 생성 완료\n")
cat(sprintf("총 CpG 사이트 수: %d\n", nrow(all_bsseq)))

# ================================================================================
# DML (Differential Methylation Loci) 테스트
# ================================================================================

cat("\nDML 테스트 수행 중...\n")

# 그룹 정의
design <- data.frame(
  sample = c(paste0("Control_", 1:length(control_data_list)),
             paste0("Experimental_", 1:length(experimental_data_list))),
  group = c(rep("Control", length(control_data_list)),
            rep("Experimental", length(experimental_data_list)))
)

# DML 테스트
dmlTest <- DMLtest(
  BSobj = all_bsseq,
  group1 = paste0("Control_", 1:length(control_data_list)),
  group2 = paste0("Experimental_", 1:length(experimental_data_list)),
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
  pdf(opt$`output-pdf`)
  plot.new()
  text(0.5, 0.5, "유의한 DMR 없음", cex=2)
  dev.off()
  
  quit(save="no", status=0)
}

cat(sprintf("발견된 DMR 수: %d\n", nrow(dmrs)))

# ================================================================================
# DMR 결과 저장
# ================================================================================

cat(sprintf("\nDMR 결과 저장 중: %s\n", opt$`output-csv`))

# 결과 테이블 작성
dmr_results <- data.frame(
  chr = dmrs$chr,
  start = dmrs$start,
  end = dmrs$end,
  length = dmrs$length,
  nCG = dmrs$nCG,
  meanMethy1 = dmrs$meanMethy1,  # Control 평균
  meanMethy2 = dmrs$meanMethy2,  # Experimental 평균
  diff_methy = dmrs$diff.Methy,
  areaStat = dmrs$areaStat,
  pvalue = dmrs$pval,
  fdr = p.adjust(dmrs$pval, method="BH")
)

# FDR로 정렬
dmr_results <- dmr_results[order(dmr_results$fdr), ]

write.csv(dmr_results, opt$`output-csv`, row.names=FALSE)

cat("DMR 결과 저장 완료\n")

# ================================================================================
# 시각화
# ================================================================================

cat(sprintf("\nDMR 시각화 중: %s\n", opt$`output-pdf`))

pdf(opt$`output-pdf`, width=12, height=8)

# 1. DMR 분포 (염색체별)
p1 <- ggplot(dmr_results, aes(x=chr, fill=diff_methy > 0)) +
  geom_bar() +
  theme_minimal() +
  labs(title="DMR 분포 (염색체별)",
       x="염색체", y="DMR 수",
       fill="메틸화 방향") +
  scale_fill_manual(values=c("blue", "red"),
                    labels=c("Hypo (감소)", "Hyper (증가)")) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

print(p1)

# 2. 메틸화 차이 분포
p2 <- ggplot(dmr_results, aes(x=diff_methy)) +
  geom_histogram(bins=50, fill="steelblue", color="black") +
  geom_vline(xintercept=0, linetype="dashed", color="red") +
  theme_minimal() +
  labs(title="메틸화 차이 분포",
       x="메틸화 차이 (Experimental - Control)",
       y="빈도")

print(p2)

# 3. Volcano plot
dmr_results$log10_fdr <- -log10(dmr_results$fdr)
p3 <- ggplot(dmr_results, aes(x=diff_methy, y=log10_fdr)) +
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

print(p3)

# 4. Top 20 DMR (FDR 기준)
top_dmrs <- head(dmr_results, 20)
top_dmrs$region <- paste(top_dmrs$chr, ":", 
                         format(top_dmrs$start, scientific=FALSE), "-",
                         format(top_dmrs$end, scientific=FALSE), sep="")

p4 <- ggplot(top_dmrs, aes(x=reorder(region, -log10_fdr), y=diff_methy)) +
  geom_bar(stat="identity", aes(fill=diff_methy > 0)) +
  coord_flip() +
  theme_minimal() +
  labs(title="Top 20 DMR (FDR 기준)",
       x="", y="메틸화 차이") +
  scale_fill_manual(values=c("blue", "red"),
                    labels=c("Hypo", "Hyper"),
                    name="") +
  theme(axis.text.y = element_text(size=8))

print(p4)

# 5. DMR 크기 vs 유의성
p5 <- ggplot(dmr_results, aes(x=length, y=log10_fdr)) +
  geom_point(aes(color=diff_methy, size=nCG), alpha=0.6) +
  scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0,
                        name="메틸화 차이") +
  scale_size_continuous(name="CpG 수") +
  theme_minimal() +
  labs(title="DMR 크기와 유의성",
       x="DMR 길이 (bp)",
       y="-log10(FDR)")

print(p5)

dev.off()

cat("시각화 완료\n")

# ================================================================================
# 요약 통계
# ================================================================================

cat("\n================================================================================\n")
cat("DMR 분석 요약\n")
cat("================================================================================\n")
cat(sprintf("총 DMR 수: %d\n", nrow(dmr_results)))
cat(sprintf("Hyper-메틸화 DMR: %d (%.1f%%)\n", 
            sum(dmr_results$diff_methy > 0),
            100 * sum(dmr_results$diff_methy > 0) / nrow(dmr_results)))
cat(sprintf("Hypo-메틸화 DMR: %d (%.1f%%)\n", 
            sum(dmr_results$diff_methy < 0),
            100 * sum(dmr_results$diff_methy < 0) / nrow(dmr_results)))
cat(sprintf("평균 DMR 길이: %.0f bp\n", mean(dmr_results$length)))
cat(sprintf("평균 CpG 사이트 수: %.1f\n", mean(dmr_results$nCG)))
cat(sprintf("평균 메틸화 차이: %.3f\n", mean(abs(dmr_results$diff_methy))))
cat("================================================================================\n")

cat("\n분석 완료!\n")
cat(sprintf("결과 파일: %s\n", opt$`output-csv`))
cat(sprintf("그래프 파일: %s\n", opt$`output-pdf`))
