################################################################################
## make_summary_table_extended()
##
## Generates a summary table identical in structure to Table 1 in the manuscript
## but expanded to include ALL 5 replicates:
##   E_15_1, E_15_2, E_15_3A, E_15_3B, E_15_3C
## across all 3 pipelines: ProxiMeta, Bin3C, MetaTOR
##
## ── REQUIRES ────────────────────────────────────────────────────────────────
##   df_combined must already exist in your environment (built by the script
##   above). It must have these 6 columns:
##     cluster_id | Completeness | Contamination | Replicate | Software | GenomeSize
##
##   Software values  : "Proximeta", "Bin3C", "MetaTOR"
##   Replicate values : "E_15_1", "E_15_2", "E_15_3A", "E_15_3B", "E_15_3C"
##   GenomeSize units : base pairs (bp)
##
##   NOTE: genome size (100kb-10Mb) is already filtered before df_combined
##   is built, so no genome size re-filtering is needed here.
##
## ── PARAMETERS ──────────────────────────────────────────────────────────────
##   df                   — df_combined tibble
##   contamination_cutoff — max contamination % to include (default 10)
##   completeness_cutoff  — min completeness % to include  (default 0)
##   output_format        — "kable" | "latex" | "df"
##                          "kable"  -> HTML table, view in RStudio Viewer pane
##                          "latex"  -> LaTeX string (paste into manuscript)
##                          "df"     -> raw data frame (for custom formatting)
##   save_latex           — if TRUE, saves .tex file (only when output_format="latex")
##   latex_path           — file path for .tex output
##
## ── QUICK START ─────────────────────────────────────────────────────────────
##   source("make_summary_table_extended.R")
##   make_summary_table_extended(df_combined)               # HTML preview
##   make_summary_table_extended(df_combined, output_format = "latex",
##                               save_latex = TRUE,
##                               latex_path = "Table1_extended.tex")
################################################################################

make_summary_table_extended <- function(
    df,
    contamination_cutoff = 10,
    completeness_cutoff  = 0,
    output_format        = "kable",
    save_latex           = FALSE,
    latex_path           = "Table1_extended.tex"
) {

  ## ── 0. packages ─────────────────────────────────────────────────────────────
  suppressPackageStartupMessages({
    library(dplyr)
    library(knitr)
    library(kableExtra)
  })

  ## ── 1. exact levels matching your script ────────────────────────────────────
  ## Software stored as "Proximeta" (lowercase x) in df_combined
  soft_internal <- c("Proximeta", "Bin3C",  "MetaTOR")
  soft_display  <- c("ProxiMeta", "Bin3C",  "MetaTOR")   # display name in table

  rep_internal  <- c("E_15_1", "E_15_2", "E_15_3A", "E_15_3B", "E_15_3C")
  rep_display   <- c("E 15 1", "E 15 2", "E 15 3A",  "E 15 3B",  "E 15 3C")

  ## ── 2. apply contamination + completeness filter ────────────────────────────
  ## (genome size already filtered upstream — no need to re-filter)
  df_filt <- df %>%
    dplyr::filter(
      Contamination <= contamination_cutoff,
      Completeness  >= completeness_cutoff,
      Replicate     %in% rep_internal,
      Software      %in% soft_internal
    )

  if (nrow(df_filt) == 0)
    stop("No rows remain after filtering. Check contamination_cutoff / completeness_cutoff.")

  ## ── 3. genome size formatter (bp -> Mb or kb) ────────────────────────────────
  fmt_gs <- function(x_bp) {
    if (is.na(x_bp)) return("--")
    if (x_bp >= 1e6) {
      paste0(formatC(x_bp / 1e6, digits = 4, format = "fg", flag = "#"), " Mb")
    } else {
      paste0(formatC(x_bp / 1e3, digits = 4, format = "fg", flag = "#"), " kb")
    }
  }

  ## percentage formatter — trims trailing zeros, keeps up to 3 sig figs
  fmt_pct <- function(x) {
    if (is.na(x)) return("--")
    paste0(formatC(x, digits = 3, format = "g"), "%")
  }

  ## ── 4. compute statistics per Software x Replicate ──────────────────────────
  stats <- df_filt %>%
    dplyr::group_by(Software, Replicate) %>%
    dplyr::summarise(
      Bins_N       = dplyr::n(),
      Comp_mean    = mean(Completeness,   na.rm = TRUE),
      Comp_median  = median(Completeness,  na.rm = TRUE),
      Cont_mean    = mean(Contamination,  na.rm = TRUE),
      Cont_median  = median(Contamination, na.rm = TRUE),
      GS_mean_bp   = mean(GenomeSize,     na.rm = TRUE),
      GS_median_bp = median(GenomeSize,    na.rm = TRUE),
      .groups = "drop"
    )

  ## ── 5. lookup helper ────────────────────────────────────────────────────────
  get_cell <- function(sw, rr, stat) {
    row <- stats %>% dplyr::filter(Software == sw, Replicate == rr)
    if (nrow(row) == 0) return("--")
    switch(stat,
      "Bins_N"      = as.character(row$Bins_N),
      "Comp_mean"   = fmt_pct(row$Comp_mean),
      "Comp_median" = fmt_pct(row$Comp_median),
      "Cont_mean"   = fmt_pct(row$Cont_mean),
      "Cont_median" = fmt_pct(row$Cont_median),
      "GS_mean"     = fmt_gs(row$GS_mean_bp),
      "GS_median"   = fmt_gs(row$GS_median_bp),
      "--"
    )
  }

  ## ── 6. build the display data frame ─────────────────────────────────────────
  ##
  ## Row layout per pipeline (mirrors Table 1 exactly):
  ##   Variable             Stat    E151  E152  E153A  E153B  E153C
  ##   Bins                 Sigma
  ##   Bin Completeness     x-bar
  ##                        x-tilde
  ##   Bin Contamination    x-bar
  ##                        x-tilde
  ##   Bin Genome Size      x-bar
  ##                        x-tilde
  ##
  ## Unicode: Sigma=\u03A3  x-bar=x\u0304  x-tilde=x\u0303

  row_defs <- data.frame(
    Variable = c(
      "Bins",
      "Bin Completeness",  "Bin Completeness",
      "Bin Contamination", "Bin Contamination",
      "Bin Genome Size",   "Bin Genome Size"
    ),
    Stat = c(
      "\u03A3",
      "x\u0304", "x\u0303",
      "x\u0304", "x\u0303",
      "x\u0304", "x\u0303"
    ),
    stat_key = c(
      "Bins_N",
      "Comp_mean",  "Comp_median",
      "Cont_mean",  "Cont_median",
      "GS_mean",    "GS_median"
    ),
    stringsAsFactors = FALSE
  )
  n_rows_per_pipe <- nrow(row_defs)  # 7

  build_block <- function(sw_int, sw_disp) {
    block <- row_defs
    for (ri in seq_along(rep_internal)) {
      rr  <- rep_internal[ri]
      col <- rep_display[ri]
      block[[col]] <- vapply(
        block$stat_key,
        function(k) get_cell(sw_int, rr, k),
        character(1)
      )
    }
    ## Pipeline label only on first row of the block
    block$Pipeline <- c(sw_disp, rep("", n_rows_per_pipe - 1))

    ## Variable label only on first row of each variable group
    seen <- character(0)
    for (i in seq_len(nrow(block))) {
      v <- block$Variable[i]
      if (v %in% seen) {
        block$Variable[i] <- ""
      } else {
        seen <- c(seen, v)
      }
    }
    block[, c("Pipeline", "Variable", "Stat", rep_display)]
  }

  tbl_df <- do.call(
    rbind,
    mapply(build_block, soft_internal, soft_display, SIMPLIFY = FALSE)
  )
  rownames(tbl_df) <- NULL

  ## ── 7. return raw data frame ─────────────────────────────────────────────────
  if (output_format == "df") return(tbl_df)

  ## ── 8. row indices for pipeline separator lines ──────────────────────────────
  ## Rows 7 and 14 get a thick bottom border (between pipeline blocks)
  pipe_break_rows <- seq(n_rows_per_pipe,
                         n_rows_per_pipe * (length(soft_internal) - 1),
                         by = n_rows_per_pipe)

  ## ── 9a. HTML / kable output ──────────────────────────────────────────────────
  if (output_format %in% c("kable", "html")) {

    cap <- paste0(
      "Ranges of summary statistics from all Hi-C library replicates ",
      "(E 15 1, E 15 2, E 15 3A, E 15 3B, E 15 3C) of retrieved MAGs from wastewater, ",
      "by different binning pipelines. ",
      "Completeness >= ", completeness_cutoff, "%, contamination <= ",
      contamination_cutoff, "%, genome size 100 kb - 10 Mb. ",
      "E 15 3A is the primary replicate; E 15 3B and E 15 3C are QC reruns."
    )

    tbl <- knitr::kable(
      tbl_df,
      format   = "html",
      booktabs = TRUE,
      linesep  = "",
      align    = c("l", "l", "c", rep("c", length(rep_display))),
      caption  = cap,
      escape   = FALSE
    ) %>%
      kableExtra::kable_styling(
        bootstrap_options = c("striped", "hover", "condensed", "bordered"),
        full_width = FALSE,
        font_size  = 11,
        position   = "center"
      ) %>%
      kableExtra::column_spec(1, bold = TRUE, width = "2cm") %>%
      kableExtra::column_spec(2, width = "3.8cm") %>%
      kableExtra::column_spec(3, width = "0.8cm") %>%
      ## Shade E 15 3B (col 6) and E 15 3C (col 7) — QC reruns
      kableExtra::column_spec(6, background = "#FFFDE7") %>%
      kableExtra::column_spec(7, background = "#FFFDE7") %>%
      ## Thick line after ProxiMeta and Bin3C blocks
      kableExtra::row_spec(
        pipe_break_rows,
        extra_css = "border-bottom: 2px solid #333;"
      ) %>%
      kableExtra::add_header_above(
        c(" " = 3,
          "Independent Hi-C replicates" = 3,
          "QC reruns of E 15 3A" = 2),
        bold      = TRUE,
        font_size = 10
      ) %>%
      kableExtra::footnote(
        general = paste0(
          "x\u0304 = mean; x\u0303 = median; \u03A3 = count. ",
          "Genome sizes in Mb (megabases) or kb (kilobases). ",
          "Shaded columns (E 15 3B, E 15 3C) are quality-control reruns ",
          "of Hi-C library E 15 3A; main analyses use E 15 3A."
        ),
        general_title     = "Note: ",
        footnote_as_chunk = TRUE
      )

    return(tbl)
  }

  ## ── 9b. LaTeX output ─────────────────────────────────────────────────────────
  if (output_format == "latex") {

    cap <- paste0(
      "Ranges of summary statistics from all Hi-C library replicates ",
      "(E\\,15\\,1, E\\,15\\,2, E\\,15\\,3A, E\\,15\\,3B, E\\,15\\,3C) ",
      "of retrieved MAGs from wastewater, by different binning pipelines. ",
      "Completeness $\\geq$ ", completeness_cutoff, "\\% (CheckM), ",
      "contamination $\\leq$ ", contamination_cutoff,
      "\\%, genome size 100\\,kb--10\\,Mb. ",
      "E\\,15\\,3A is the primary replicate used in the main comparative analysis; ",
      "E\\,15\\,3B and E\\,15\\,3C are quality-control reruns of the same Hi-C library."
    )

    tbl <- knitr::kable(
      tbl_df,
      format   = "latex",
      booktabs = TRUE,
      linesep  = "",
      align    = c("l", "l", "c", rep("c", length(rep_display))),
      caption  = cap,
      escape   = FALSE
    ) %>%
      kableExtra::kable_styling(
        latex_options = c("hold_position", "scale_down"),
        font_size     = 9
      ) %>%
      kableExtra::column_spec(1, bold = TRUE) %>%
      ## Light yellow shading for QC columns
      kableExtra::column_spec(6, background = "yellow!15") %>%
      kableExtra::column_spec(7, background = "yellow!15") %>%
      kableExtra::row_spec(pipe_break_rows, hline_after = TRUE) %>%
      kableExtra::add_header_above(
        c(" " = 3,
          "Independent Hi-C replicates" = 3,
          "QC reruns of E\\,15\\,3A" = 2),
        bold   = TRUE,
        escape = FALSE
      ) %>%
      kableExtra::footnote(
        general = paste0(
          "$\\\\bar{x}$ = mean; $\\\\tilde{x}$ = median; $\\\\Sigma$ = count. ",
          "Genome sizes in Mb or kb. ",
          "Shaded columns (E\\\\,15\\\\,3B, E\\\\,15\\\\,3C) are QC reruns of E\\\\,15\\\\,3A."
        ),
        general_title     = "Note: ",
        footnote_as_chunk = TRUE,
        escape            = FALSE
      )

    if (save_latex) {
      writeLines(as.character(tbl), latex_path)
      message("LaTeX table written to: ", latex_path)
    }

    return(tbl)
  }

  stop("output_format must be one of: 'kable', 'html', 'latex', 'df'")
}


################################################################################
## VALIDATION HELPER — validate_vs_table1()
##
## Cross-checks ALL 9 cells from manuscript Table 1 (3 pipelines x 3 replicates)
## against df_combined, then prints a summary of the new E_15_3B / E_15_3C values.
##
## ALL 9 expected values from Table 1
## (Contamination <= 10%, GenomeSize 100kb-10Mb already filtered upstream):
##
##   ProxiMeta / E_15_1:  N=126  comp mean=40.4%  median=31.3%
##   ProxiMeta / E_15_2:  N=126  comp mean=42.6%  median=30.0%
##   ProxiMeta / E_15_3A: N=131  comp mean=41.2%  median=36.2%
##   Bin3C     / E_15_1:  N=116  comp mean=34.9%  median=9.56%
##   Bin3C     / E_15_2:  N=108  comp mean=38.2%  median=34.7%
##   Bin3C     / E_15_3A: N=111  comp mean=36.1%  median=22.5%
##   MetaTOR   / E_15_1:  N=101  comp mean=56.5%  median=61.1%
##   MetaTOR   / E_15_2:  N=89   comp mean=57.2%  median=61.3%
##   MetaTOR   / E_15_3A: N=88   comp mean=55.0%  median=59.4%
##
## Confirmed new replicates (from live session 2026-03-09):
##   ProxiMeta / E_15_3B: N=133  comp 40.8%/28.9%  cont 1.181%  GS 1.297 Mb
##   ProxiMeta / E_15_3C: N=131  comp 41.3%/28.9%  cont 1.290%  GS 1.349 Mb
##   Bin3C     / E_15_3B: N=108  comp 37.4%/30.2%  cont 0.849%  GS 1.222 Mb
##   Bin3C     / E_15_3C: N=108  comp 37.1%/27.2%  cont 0.816%  GS 1.211 Mb
##   MetaTOR   / E_15_3B: N=89   comp 54.6%/58.5%  cont 3.017%  GS 1.791 Mb
##   MetaTOR   / E_15_3C: N=88   comp 55.0%/59.4%  cont 2.987%  GS 1.801 Mb
################################################################################

validate_vs_table1 <- function(df) {

  cat("\n=== Validation: ALL 9 cells from manuscript Table 1 ===\n")
  cat("Filter: Contamination <= 10% (GenomeSize already filtered upstream)\n\n")

  df_filt <- df %>% dplyr::filter(Contamination <= 10)

  ## All 9 values from Table 1: 3 pipelines x 3 original replicates
  ## All contamination values read directly from manuscript Table 1
  expected <- list(
    ## ProxiMeta
    list(sw="Proximeta", rr="E_15_1",  N=126,
         comp_mean=40.4, comp_med=31.3, cont_mean=1.27,  cont_med=0.0200),
    list(sw="Proximeta", rr="E_15_2",  N=126,
         comp_mean=42.6, comp_med=30.0, cont_mean=1.24,  cont_med=0.260),
    list(sw="Proximeta", rr="E_15_3A", N=131,
         comp_mean=41.2, comp_med=36.2, cont_mean=1.17,  cont_med=0.0500),
    ## Bin3C
    list(sw="Bin3C",     rr="E_15_1",  N=116,
         comp_mean=34.9, comp_med=9.56, cont_mean=0.816, cont_med=0.0000),
    list(sw="Bin3C",     rr="E_15_2",  N=108,
         comp_mean=38.2, comp_med=34.7, cont_mean=0.789, cont_med=0.0000),
    list(sw="Bin3C",     rr="E_15_3A", N=111,
         comp_mean=36.1, comp_med=22.5, cont_mean=0.843, cont_med=0.0000),
    ## MetaTOR
    list(sw="MetaTOR",   rr="E_15_1",  N=101,
         comp_mean=56.5, comp_med=61.1, cont_mean=3.53,  cont_med=3.04),
    list(sw="MetaTOR",   rr="E_15_2",  N=89,
         comp_mean=57.2, comp_med=61.3, cont_mean=3.41,  cont_med=2.97),
    list(sw="MetaTOR",   rr="E_15_3A", N=88,
         comp_mean=55.0, comp_med=59.4, cont_mean=2.97,  cont_med=2.35)
  )

  all_pass <- TRUE
  for (e in expected) {
    d <- df_filt %>% dplyr::filter(Software == e$sw, Replicate == e$rr)

    n_ok     <- nrow(d) == e$N
    cm_ok    <- abs(mean(d$Completeness)    - e$comp_mean) < 0.15
    cmed_ok  <- abs(median(d$Completeness)  - e$comp_med)  < 0.15
    ktm_ok   <- abs(mean(d$Contamination)   - e$cont_mean) < 0.015
    ktmed_ok <- abs(median(d$Contamination) - e$cont_med)  < 0.005
    pass     <- n_ok & cm_ok & cmed_ok & ktm_ok & ktmed_ok
    all_pass <- all_pass & pass

    cat(sprintf(
      "  %-12s / %-10s  N=%d [%s]\n",
      e$sw, e$rr, nrow(d),
      ifelse(n_ok, "OK", paste0("EXPECTED=", e$N))
    ))
    cat(sprintf(
      "    comp  mean=%.2f%% [%s]  median=%.2f%% [%s]\n",
      mean(d$Completeness),
      ifelse(cm_ok,   "OK", paste0("EXPECTED=", e$comp_mean)),
      median(d$Completeness),
      ifelse(cmed_ok, "OK", paste0("EXPECTED=", e$comp_med))
    ))
    cat(sprintf(
      "    cont  mean=%.4f%% [%s]  median=%.4f%% [%s]\n\n",
      mean(d$Contamination),
      ifelse(ktm_ok,   "OK", paste0("EXPECTED=", e$cont_mean)),
      median(d$Contamination),
      ifelse(ktmed_ok, "OK", paste0("EXPECTED=", e$cont_med))
    ))
  }

  cat("\n")
  if (all_pass) {
    cat("All checks PASSED -- df_combined matches Table 1.\n")
  } else {
    cat("Some checks FAILED -- review df_combined construction.\n")
  }

  cat("\n--- New replicates (E_15_3B, E_15_3C) summary ---\n")
  for (sw in c("Proximeta", "Bin3C", "MetaTOR")) {
    for (rr in c("E_15_3B", "E_15_3C")) {
      d <- df_filt %>% dplyr::filter(Software == sw, Replicate == rr)
      if (nrow(d) == 0) {
        cat(sprintf("  %-12s / %s: NOT FOUND in df_combined\n", sw, rr))
      } else {
        cat(sprintf(
          "  %-12s / %s: N=%d  comp mean=%.1f%%  median=%.1f%%  cont mean=%.3f%%  GS mean=%.3f Mb\n",
          sw, rr, nrow(d),
          mean(d$Completeness),  median(d$Completeness),
          mean(d$Contamination), mean(d$GenomeSize) / 1e6
        ))
      }
    }
  }
  invisible(NULL)
}


################################################################################
## USAGE
## Paste these lines into your script AFTER df_combined is fully built.
################################################################################
##
## source("make_summary_table_extended.R")
##
## # Step 1: confirm df_combined matches existing Table 1
## validate_vs_table1(df_combined)
##
## # Step 2a: HTML preview in RStudio Viewer
## make_summary_table_extended(df_combined)
##
## # Step 2b: LaTeX for manuscript
## make_summary_table_extended(
##   df_combined,
##   output_format = "latex",
##   save_latex    = TRUE,
##   latex_path    = "~/Desktop/Table1_extended.tex"
## )
##
## # Step 2c: Raw data frame for Word/Excel
## tbl_df <- make_summary_table_extended(df_combined, output_format = "df")
## View(tbl_df)
##
## # Step 2d: Different completeness cutoff (e.g. >= 50%)
## make_summary_table_extended(df_combined, completeness_cutoff = 50)
##
################################################################################
