#!/usr/bin/env Rscript
# =============================================================================
# cohort 리포트 공통 헬퍼 — 여러 cohort 분석 스크립트(variant_burden.R,
# sv_burden.R, geneset_burden.R, dmr_analysis.R, asm_analysis.R,
# trgt_group_compare.R, group_summary_compare.R, plot_hifi_native_results.R)에서
# 중복되던 상수/함수를 한 곳으로 모은 것. 각 스크립트는 이 파일을 source()해서 쓴다.
#
# Snakemake는 항상 저장소 루트에서 `Rscript scripts/xxx.R`를 실행하므로, 각
# 스크립트는 단순히 source("scripts/report_common.R")로 이 파일을 불러온다.
# report/cohort_report.Rmd도 knit_root_dir을 저장소 루트로 맞춰 동일하게 source한다.
# =============================================================================

suppressPackageStartupMessages({
  library(svglite)
})

# ── H2O2 4군 공통 스타일 (group_summary_compare.R / plot_hifi_native_results.R 중복 통합) ──
FONT_FAMILY  <- "Noto Sans CJK KR"
GROUP_LEVELS <- c("A", "B", "C", "D")
GROUP_COLORS <- c(A="#E41A1C", B="#377EB8", C="#4DAF4A", D="#984EA3")
GROUP_LABELS <- c(A="A\n(대조군)", B="B\n(급성)", C="C\n(급성+약물)", D="D\n(만성+약물)")

base_theme <- ggplot2::theme_bw(base_size=13, base_family=FONT_FAMILY) +
  ggplot2::theme(
    text = ggplot2::element_text(family=FONT_FAMILY),
    plot.title = ggplot2::element_text(face="bold", size=13, family=FONT_FAMILY),
    plot.subtitle = ggplot2::element_text(size=9.5, color="grey30", family=FONT_FAMILY),
    legend.position = "none",
    panel.grid.minor = ggplot2::element_blank()
  )

wrap_txt <- function(x, width=55) paste(strwrap(x, width=width), collapse="\n")

kw_pvalue <- function(values, groups) {
  kt <- kruskal.test(values, groups)
  list(H=unname(kt$statistic), p=kt$p.value)
}

# ── 그룹별 카운트 → 효과크기 (variant_burden.R / sv_burden.R 중복 통합) ─────────
# counts: named numeric, 그룹별 이벤트/carrier 수. n_grp: named numeric, 그룹별 표본 수.
# 반환: max_prop_diff(침투율 차이), group_max/group_min, odds_ratio_extreme
#       (Haldane-Anscombe 보정 최고/최저군 오즈비), cramers_v(카이제곱 기반 효과크기)
compute_group_effect_size <- function(counts, n_grp) {
  props <- counts / n_grp
  max_prop_diff <- max(props) - min(props)
  grp_max <- names(which.max(props))[1]
  grp_min <- names(which.min(props))[1]

  a  <- counts[grp_max] + 0.5; b <- n_grp[grp_max] - counts[grp_max] + 0.5
  cc <- counts[grp_min] + 0.5; d <- n_grp[grp_min] - counts[grp_min] + 0.5
  or_extreme <- unname((a * d) / (b * cc))

  ct_v <- tryCatch(
    suppressWarnings(chisq.test(rbind(counts, n_grp - counts), correct=FALSE)),
    error=function(e) NULL
  )
  cramers_v <- if (!is.null(ct_v)) unname(sqrt(ct_v$statistic / sum(n_grp))) else NA_real_

  list(max_prop_diff=max_prop_diff, group_max=grp_max, group_min=grp_min,
       odds_ratio_extreme=or_extreme, cramers_v=cramers_v)
}

# ── 그룹별 이진/카운트 벡터에 대한 Fisher(2군)/Chi-square(3군+) 검정 ────────────
# variant_burden.R / sv_burden.R의 "2군이면 Fisher, 3군+면 chi-square" 분기 통합.
group_presence_test <- function(counts, n_grp) {
  if (length(n_grp) == 2) {
    tbl <- matrix(c(counts[1], n_grp[1] - counts[1],
                    counts[2], n_grp[2] - counts[2]), nrow=2)
    ft <- tryCatch(fisher.test(tbl), error=function(e) NULL)
    list(pvalue = if (!is.null(ft)) ft$p.value else NA_real_,
         odds_ratio = if (!is.null(ft) && !is.null(ft$estimate)) unname(ft$estimate) else NA_real_)
  } else {
    tbl <- rbind(counts, n_grp - counts)
    ct <- tryCatch(chisq.test(tbl, simulate.p.value=(any(tbl < 5))), error=function(e) NULL)
    list(pvalue = if (!is.null(ct)) ct$p.value else NA_real_, odds_ratio = NA_real_)
  }
}

# ── BCSQ 필드 파서 (sv_burden.R / geneset_burden.R 중복 통합) ──────────────────
# BCSQ 형식: "sv:cds|GENE|transcript|biotype|strand|aa_change" (콤마 구분 다중 entry)
parse_bcsq <- function(bcsq_str) {
  if (is.na(bcsq_str) || bcsq_str == "." || bcsq_str == "") return(data.frame())
  entries <- strsplit(bcsq_str, ",")[[1]]
  do.call(rbind, lapply(entries, function(e) {
    parts <- strsplit(e, "\\|")[[1]]
    if (length(parts) < 2) return(NULL)
    data.frame(sv_consequence=parts[1], gene=parts[2], stringsAsFactors=FALSE)
  }))
}

# ── PacBio CpG BED → DSS 입력 형식 (dmr_analysis.R / asm_analysis.R 중복 통합) ──
# on_error="warn": 실패 시 NULL 반환하고 계속 진행 (dmr_analysis.R 방식, 여러 샘플 중 일부 스킵 허용)
# on_error="stop": 실패 시 즉시 stop() (asm_analysis.R 방식, hap1/hap2 둘 다 필수)
read_cpg_bed <- function(bed_file, label=NULL, on_error=c("warn", "stop")) {
  on_error <- match.arg(on_error)
  fail <- function(msg) {
    if (on_error == "stop") stop(msg, call.=FALSE)
    cat(sprintf("  ✗ %s\n", msg))
    NULL
  }
  tag <- if (!is.null(label)) sprintf("[%s] ", label) else ""
  cat(sprintf("%s파일 로딩 중: %s\n", tag, bed_file))

  if (!file.exists(bed_file)) return(fail(sprintf("파일이 존재하지 않습니다: %s", bed_file)))
  file_size <- file.info(bed_file)$size
  if (file_size == 0) return(fail(sprintf("파일이 비어있습니다: %s", bed_file)))
  cat(sprintf("  파일 크기: %.2f MB\n", file_size / (1024^2)))

  tryCatch({
    dt <- data.table::fread(bed_file, skip="#chrom", header=TRUE, sep="\t", showProgress=FALSE)
    if (nrow(dt) == 0) return(fail(sprintf("파일에 데이터가 없습니다: %s", bed_file)))
    if (ncol(dt) < 7) return(fail(sprintf("컬럼 수 부족 (최소 7개 필요, 현재 %d개): %s", ncol(dt), bed_file)))

    expected_cols <- c("chrom", "begin", "end", "mod_score", "type",
                        "cov", "est_mod_count", "est_unmod_count", "discretized_mod_score")
    if (ncol(dt) >= length(expected_cols)) colnames(dt)[1:length(expected_cols)] <- expected_cols

    dss_data <- data.frame(
      chr = as.character(dt$chrom),
      pos = as.integer(dt$begin),
      N   = as.integer(dt$cov),
      X   = as.integer(dt$est_mod_count)
    )
    dss_data <- na.omit(dss_data)
    dss_data <- dss_data[dss_data$N > 0, ]
    if (nrow(dss_data) == 0) return(fail("유효한 CpG 사이트가 없습니다 (coverage > 0)"))

    cat(sprintf("  ✓ %s CpG 사이트 로드됨\n", format(nrow(dss_data), big.mark=",")))
    dss_data
  }, error = function(e) fail(sprintf("%s: %s", bed_file, e$message)))
}

# ── SVG 저장 헬퍼 ───────────────────────────────────────────────────────────
# plots: named list. 이름이 곧 파일명 stem이 됨 (예: "01_top30_genes" ->
#        01_top30_genes.svg). 값은 ggplot/patchwork 객체이거나,
#        list(plot=..., width=, height=)로 개별 크기를 지정할 수 있음.
save_svg_plots <- function(plots, out_dir, width=12, height=9) {
  dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
  for (nm in names(plots)) {
    entry <- plots[[nm]]
    if (is.list(entry) && !is.null(entry$plot) && !inherits(entry, "ggplot") && !inherits(entry, "patchwork")) {
      p <- entry$plot
      w <- if (!is.null(entry$width)) entry$width else width
      h <- if (!is.null(entry$height)) entry$height else height
    } else {
      p <- entry
      w <- width
      h <- height
    }
    out_path <- file.path(out_dir, paste0(nm, ".svg"))
    svglite::svglite(out_path, width=w, height=h)
    print(p)
    grDevices::dev.off()
  }
  message("SVG 저장 (", length(plots), "개): ", out_dir)
}

# 데이터 없음 등으로 플롯을 생성할 수 없을 때 Snakemake directory() output을
# 비어있지 않게 만들기 위한 placeholder SVG.
write_placeholder_svg <- function(out_dir, message_text, width=8, height=6) {
  dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
  svglite::svglite(file.path(out_dir, "00_no_data.svg"), width=width, height=height)
  plot.new()
  text(0.5, 0.5, message_text, cex=1.3)
  grDevices::dev.off()
}

# ── fig-atlas 번들 metadata.yaml 작성 헬퍼 ──────────────────────────────────
# fields: figure_number, figure_slug, figure_title, source_stem, primary_code,
#         panel_notes(character vector), input_files/output_files/statistics_tables
#         (character vector, 경로), claim_boundary(character), extra_notes(character vector, 선택)
write_fig_atlas_metadata <- function(fields, path) {
  meta <- list(
    figure_number     = fields$figure_number,
    figure_slug       = fields$figure_slug,
    figure_title      = fields$figure_title,
    source_stem       = fields$source_stem,
    primary_code      = fields$primary_code,
    output_formats    = list("png", "pdf", "svg"),
    panel_notes       = as.list(fields$panel_notes),
    input_files       = as.list(fields$input_files),
    output_files      = as.list(fields$output_files),
    statistics_tables = as.list(fields$statistics_tables),
    claim_boundary    = fields$claim_boundary,
    notes = as.list(c(
      "primary_code is R (figure.R), not Python, per upstream pipeline language — deviation from external_pipeline_requirements.md §2.2, agreed with atlas maintainer.",
      fields$extra_notes
    )),
    environment = list(
      r_version    = paste(R.version$major, R.version$minor, sep="."),
      package_versions = list(
        ggplot2 = as.character(utils::packageVersion("ggplot2")),
        svglite = as.character(utils::packageVersion("svglite"))
      ),
      generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    )
  )
  dir.create(dirname(path), showWarnings=FALSE, recursive=TRUE)
  yaml::write_yaml(meta, path)
  message("메타데이터 저장: ", path)
}

# ── report/cohort_report.Rmd 전용 헬퍼 ──────────────────────────────────────
# (DT/knitr는 여기서 library()하지 않고 ::로만 사용 — 배치 스크립트에서 이 파일을
#  source()해도 불필요한 패키지가 로드되지 않도록)

# rel_dir(예: "variant_burden_svg")에 있는 SVG들을 순서대로 <img> 태그로 삽입.
# cohort_dir(예: "h2o2_analysis_results/cohort")는 렌더링될 HTML 파일이 위치할
# 디렉토리 — img src는 그 디렉토리 기준 상대경로로 써서 결과물을 옮겨도 안 깨지게 함.
# results='asis' 청크에서 호출해야 함.
svg_gallery <- function(rel_dir, cohort_dir) {
  abs_dir <- file.path(cohort_dir, rel_dir)
  if (!dir.exists(abs_dir)) { cat("_(그림 없음:", rel_dir, ")_\n\n"); return(invisible(NULL)) }
  files <- sort(list.files(abs_dir, pattern="\\.svg$"))
  if (length(files) == 0) { cat("_(그림 없음:", rel_dir, ")_\n\n"); return(invisible(NULL)) }
  for (f in files) {
    label <- tools::file_path_sans_ext(f)
    cat(sprintf('\n\n**%s**\n\n<img src="%s" style="max-width:100%%; height:auto;" />\n\n',
                label, file.path(rel_dir, f)))
  }
  invisible(NULL)
}

# csv_path의 결과 테이블을 정렬/검색 가능한 DT 테이블로 렌더. fdr_col이 존재하면
# fdr_cutoff 미만인 셀을 강조색으로 표시. 일반 R 청크(결과 자동 출력)에서 호출.
#
# 주의(max_rows): DT는 client-side 모드에서 pageLength와 무관하게 전체 데이터를
# 브라우저로 보내 임베드한다 — DMR 결과처럼 5~6만 행짜리 테이블을 그대로 넘기면
# 리포트 HTML이 테이블 하나당 수 MB씩 부풀어(총 40MB+ 관측) 열기 느려지고 공유가
# 어려워진다. max_rows 초과 시 fdr_col(없으면 |areaStat|) 기준 정렬 후 상위
# max_rows행만 표시하고, 잘렸다는 사실을 caption으로 명시한다 — 전체 결과는 원본
# CSV를 계속 참고하면 된다.
fdr_table <- function(csv_path, fdr_col="fdr", fdr_cutoff=0.05, page_length=15, max_rows=500) {
  if (!file.exists(csv_path)) return(cat("_(데이터 없음:", basename(csv_path), ")_\n\n"))
  df <- tryCatch(read.csv(csv_path, stringsAsFactors=FALSE), error=function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(cat("_(빈 결과:", basename(csv_path), ")_\n\n"))

  n_total <- nrow(df)
  caption <- NULL
  if (n_total > max_rows) {
    if (fdr_col %in% colnames(df) && !all(is.na(df[[fdr_col]]))) {
      df <- df[order(df[[fdr_col]], na.last=TRUE), ]
    } else if ("areaStat" %in% colnames(df)) {
      df <- df[order(-abs(df$areaStat)), ]
    }
    df <- head(df, max_rows)
    caption <- htmltools::tags$caption(
      style="caption-side: top; text-align: left; color: #a94442; font-weight: bold;",
      sprintf("상위 %d행만 표시 (전체 %s행) — 전체 결과는 %s 참고",
              max_rows, format(n_total, big.mark=","), basename(csv_path))
    )
  }

  dt <- DT::datatable(df, filter="top", options=list(pageLength=page_length, scrollX=TRUE),
                       rownames=FALSE, caption=caption)
  if (fdr_col %in% colnames(df)) {
    dt <- DT::formatStyle(dt, fdr_col,
      backgroundColor = DT::styleInterval(fdr_cutoff, c("#ffd6d6", "white")))
  }
  dt
}
