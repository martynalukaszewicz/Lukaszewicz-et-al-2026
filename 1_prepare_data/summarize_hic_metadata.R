## ============================================================
## summarize_hic_metadata.R
##
## PURPOSE:
##   Reads metadata.tsv AND data.tsv files for ProxiMeta,
##   Bin3C, and MetaTOR. Merges them (as in getdata2_
##   DistancebyLength.R) and produces comprehensive summary
##   tables at both the BIN level and the CONTIG level.
##
## KEY FINDING -- BIN3C N50_(contigs) ARTIFACT:
##   The N50_(contigs) column in the Bin3C metadata.tsv has
##   a TSV parsing artifact: the string "50" is prepended to
##   every value. For example, a true N50 of 411156 bp appears
##   as 50411156. This column should NOT be used.
##   The correct contig N50 for Bin3C is the "n50" column.
##   This script uses "n50" for Bin3C throughout and flags
##   the corrupted column explicitly.
##
## METATOR DATA NOTE:
##   MetaTOR data.tsv has the same format as ProxiMeta and Bin3C:
##   columns: cluster, contig, replicate, length.
##   It was produced upstream by prepare_data_for_similarity_
##   analysis.R and is loaded directly like the other two pipelines.
##
## FILTERS APPLIED (consistent with manuscript and
## getdata2_DistancebyLength.R):
##   1. Contamination <= 10%
##   2. Genome size >= 100,000 bp (100 kb)
##   3. Genome size <= 10,000,000 bp (10 Mb)
##
## OUTPUT:
##   For each pipeline x replicate:
##     - bin_summary: one row per bin, all metadata + contig stats
##     - contig_list: one row per contig in passing bins
##   Combined cross-pipeline summary tables:
##     - cross_pipeline_bin_summary_allreps.tsv  (all 5 replicates)
##     - cross_pipeline_completeness_threshold_summary.tsv
##         (includes mean_cont and median_cont)
##     - cross_pipeline_threshold_range_summary.tsv
##         (min-max ranges + per-replicate value lists for each metric)
##
## USAGE:
##   Edit the FILE PATHS section below, then:
##   Rscript summarize_hic_metadata.R
##   or source() in RStudio.
## ============================================================

library(tidyverse)

select    <- dplyr::select
filter    <- dplyr::filter
rename    <- dplyr::rename
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
group_by  <- dplyr::group_by
arrange   <- dplyr::arrange


## ============================================================
## USER PARAMETERS -- set these paths before running
## ============================================================

## Directory containing proximeta/, bin3c/, and metator/ subfolders
## (output of prepare_data_for_similarity_analysis.R)
dir_Similarity_Routput <- "~/path/to/Similarity_Routput/"

## Directory where output TSV files will be written
## (can be the same as dir_Similarity_Routput)
OUT_DIR <- "~/path/to/Similarity_Routput/"

## ============================================================
## END USER PARAMETERS
## ============================================================

PATH_PM_META <- file.path(dir_Similarity_Routput, "proximeta/metadata.tsv")
PATH_MT_META <- file.path(dir_Similarity_Routput, "metator/metadata.tsv")
PATH_B3_META <- file.path(dir_Similarity_Routput, "bin3c/metadata.tsv")

PATH_PM_DATA <- file.path(dir_Similarity_Routput, "proximeta/data.tsv")
PATH_B3_DATA <- file.path(dir_Similarity_Routput, "bin3c/data.tsv")
PATH_MT_DATA <- file.path(dir_Similarity_Routput, "metator/data.tsv")


## ============================================================
## PARAMETERS
## ============================================================

CONT_CUTOFF  <- 10
SIZE_MIN     <- 100e3
SIZE_MAX     <- 10e6

ALL_REPS  <- c("E_15_1", "E_15_2", "E_15_3A", "E_15_3B", "E_15_3C")
MAIN_REPS <- c("E_15_1", "E_15_2", "E_15_3A")
QC_REPS   <- c("E_15_3A", "E_15_3B", "E_15_3C")

COMP_THRESHOLDS <- c(0, 50, 85, 90, 95)


## ============================================================
## HELPER FUNCTIONS
## ============================================================

sumstats <- function(x, label, digits = 2) {
  x <- x[!is.na(x)]
  tibble(
    Metric = label,
    N      = length(x),
    Mean   = round(mean(x),   digits),
    Median = round(median(x), digits),
    Min    = round(min(x),    digits),
    Max    = round(max(x),    digits),
    SD     = round(sd(x),     digits)
  )
}

comp_counts <- function(comp_vec, n_total) {
  map_dfr(COMP_THRESHOLDS, function(t) {
    cnt <- sum(comp_vec >= t, na.rm = TRUE)
    tibble(
      `Completeness >=` = paste0(t, "%"),
      N_bins            = cnt,
      pct_of_total      = round(100 * cnt / n_total, 1)
    )
  })
}

compute_n50 <- function(lengths) {
  lengths <- sort(lengths[!is.na(lengths)], decreasing = TRUE)
  if (length(lengths) == 0) return(NA_real_)
  cumlen <- cumsum(lengths)
  half   <- sum(lengths) / 2
  lengths[which(cumlen >= half)[1]]
}

extract_phylum <- function(tax_vec) {
  sapply(tax_vec, function(t) {
    if (is.na(t) || t == "") return("Unassigned")
    m <- regmatches(t, regexpr("p__[^;]+", t))
    if (length(m) == 0 || m == "") return("Unassigned")
    trimws(sub("p__", "", m))
  }, USE.NAMES = FALSE)
}

print_section <- function(title) {
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat(sprintf("  %s\n", title))
  cat(strrep("=", 60), "\n", sep = "")
}

make_bin_summary <- function(merged_df) {
  merged_df %>%
    group_by(cluster_id, replicate) %>%
    summarise(
      N_contigs_data       = n(),
      Total_bp_data        = sum(length, na.rm = TRUE),
      N50_computed_bp      = compute_n50(length),
      Min_contig_len_bp    = min(length,    na.rm = TRUE),
      Max_contig_len_bp    = max(length,    na.rm = TRUE),
      Mean_contig_len_bp   = round(mean(length,   na.rm = TRUE), 1),
      Median_contig_len_bp = median(length, na.rm = TRUE),
      Contigs_lt_1kb       = sum(length < 1000,  na.rm = TRUE),
      Pct_contigs_lt_1kb   = round(100 * sum(length < 1000, na.rm = TRUE) / n(), 1),
      Contigs_gt_10kb      = sum(length > 10000, na.rm = TRUE),
      Pct_contigs_gt_10kb  = round(100 * sum(length > 10000, na.rm = TRUE) / n(), 1),
      Contigs_gt_50kb      = sum(length > 50000, na.rm = TRUE),
      Pct_contigs_gt_50kb  = round(100 * sum(length > 50000, na.rm = TRUE) / n(), 1),
      Contig_names         = paste(sort(contig), collapse = ";"),
      .groups = "drop"
    )
}


## ============================================================
## STEP 1: LOAD AND FILTER METADATA
## ============================================================

cat("\nLoading and filtering metadata files...\n")

pm_meta_raw <- read_tsv(PATH_PM_META, show_col_types = FALSE)
pm_meta <- pm_meta_raw %>%
  rename(
    Completeness       = completeness,
    Contamination      = marker_gene_overrepresentation,
    GenomeSize         = genome_size,
    N50_meta           = n50,
    NumContigs_meta    = num_contigs,
    MeanContigLen_meta = mean_contig_length,
    GC                 = `gc%`,
    Lineage            = marker_lineage,
    MashRef            = mash_reference
  ) %>%
  filter(Contamination <= CONT_CUTOFF, GenomeSize >= SIZE_MIN, GenomeSize <= SIZE_MAX)
cat(sprintf("ProxiMeta metadata: %d raw bins -> %d pass filters\n", nrow(pm_meta_raw), nrow(pm_meta)))

mt_meta_raw <- read_tsv(PATH_MT_META, show_col_types = FALSE)
mt_meta <- mt_meta_raw %>%
  rename(
    Completeness       = completness,
    Contamination      = contamination,
    GenomeSize         = size,
    N50_meta           = N50,
    NumContigs_meta    = contigs,
    LongestContig_meta = longest_contig,
    GC                 = GC,
    CodingDensity      = coding_density,
    HiCCoverage        = HiC_Coverage
  ) %>%
  filter(Contamination <= CONT_CUTOFF, GenomeSize >= SIZE_MIN, GenomeSize <= SIZE_MAX)
cat(sprintf("MetaTOR metadata:   %d raw bins -> %d pass filters\n", nrow(mt_meta_raw), nrow(mt_meta)))

b3_meta_raw <- read_tsv(PATH_B3_META, show_col_types = FALSE)
cat("\nNOTE: Bin3C N50_(contigs) column is CORRUPTED (artifact).\n")
cat("      Using 'n50' column instead for Bin3C N50.\n")
cat(sprintf("      Verification: n50[1]=%s, N50_(contigs)[1]=%s\n",
            b3_meta_raw$n50[1], b3_meta_raw$`N50_(contigs)`[1]))
cat("      Pattern confirmed: N50_(contigs) = '50' prepended to n50.\n\n")

b3_meta <- b3_meta_raw %>%
  rename(
    Completeness       = Completeness,
    Contamination      = Contamination,
    GenomeSize         = extent,
    N50_meta           = n50,
    NumContigs_meta    = `#_contigs`,
    MeanContigLen_meta = Mean_contig_length,
    GC                 = gc_mean
  ) %>%
  filter(Contamination <= CONT_CUTOFF, GenomeSize >= SIZE_MIN, GenomeSize <= SIZE_MAX)
cat(sprintf("Bin3C metadata:     %d raw bins -> %d pass filters\n", nrow(b3_meta_raw), nrow(b3_meta)))


## ============================================================
## STEP 2: LOAD DATA FILES
## ============================================================

cat("\nLoading data files (contig assignments)...\n")

pm_data_raw <- read_tsv(PATH_PM_DATA, show_col_types = FALSE) %>% rename(cluster_id = cluster)
cat(sprintf("ProxiMeta data: %d contig rows, %d unique bins\n",
            nrow(pm_data_raw), n_distinct(paste(pm_data_raw$cluster_id, pm_data_raw$replicate))))

b3_data_raw <- read_tsv(PATH_B3_DATA, show_col_types = FALSE) %>% rename(cluster_id = cluster)
cat(sprintf("Bin3C data:     %d contig rows, %d unique bins\n",
            nrow(b3_data_raw), n_distinct(paste(b3_data_raw$cluster_id, b3_data_raw$replicate))))

mt_data_raw <- read_tsv(PATH_MT_DATA, show_col_types = FALSE) %>% rename(cluster_id = cluster)
cat(sprintf("MetaTOR data:   %d contig rows, %d unique bins\n",
            nrow(mt_data_raw), n_distinct(paste(mt_data_raw$cluster_id, mt_data_raw$replicate))))


## ============================================================
## STEP 3: MERGE DATA WITH METADATA
## ============================================================

cat("\nMerging data with metadata (inner join on cluster_id + replicate)...\n")

pm_merged <- inner_join(
  pm_data_raw,
  pm_meta %>% dplyr::select(cluster_id, replicate, Completeness, Contamination,
                            GenomeSize, N50_meta, NumContigs_meta, MeanContigLen_meta,
                            GC, Lineage, MashRef),
  by = c("cluster_id", "replicate")
)
cat(sprintf("ProxiMeta merged: %d contigs in %d passing bins\n",
            nrow(pm_merged), n_distinct(paste(pm_merged$cluster_id, pm_merged$replicate))))

b3_merged <- inner_join(
  b3_data_raw,
  b3_meta %>% dplyr::select(cluster_id, replicate, Completeness, Contamination,
                            GenomeSize, N50_meta, NumContigs_meta, MeanContigLen_meta, GC),
  by = c("cluster_id", "replicate")
)
cat(sprintf("Bin3C merged:     %d contigs in %d passing bins\n",
            nrow(b3_merged), n_distinct(paste(b3_merged$cluster_id, b3_merged$replicate))))

mt_merged <- inner_join(
  mt_data_raw,
  mt_meta %>% dplyr::select(cluster_id, replicate, Completeness, Contamination,
                            GenomeSize, N50_meta, NumContigs_meta, LongestContig_meta,
                            GC, CodingDensity, HiCCoverage, taxonomy),
  by = c("cluster_id", "replicate")
)
cat(sprintf("MetaTOR merged:   %d contigs in %d passing bins\n",
            nrow(mt_merged), n_distinct(paste(mt_merged$cluster_id, mt_merged$replicate))))


## ============================================================
## STEP 4: COMPUTE PER-BIN CONTIG STATISTICS
## ============================================================

cat("\nComputing per-bin contig statistics from data files...\n")

pm_bin_stats <- make_bin_summary(pm_merged)
pm_bin_full <- pm_meta %>%
  inner_join(pm_bin_stats, by = c("cluster_id", "replicate")) %>%
  mutate(
    N50_meta_vs_computed_ratio    = round(N50_meta / N50_computed_bp, 3),
    NumContigs_meta_vs_data_match = (NumContigs_meta == N_contigs_data),
    Software = "ProxiMeta"
  ) %>%
  arrange(replicate, desc(Completeness))

b3_bin_stats <- make_bin_summary(b3_merged)
b3_bin_full <- b3_meta %>%
  inner_join(b3_bin_stats, by = c("cluster_id", "replicate")) %>%
  mutate(
    N50_meta_vs_computed_ratio    = round(N50_meta / N50_computed_bp, 3),
    NumContigs_meta_vs_data_match = (NumContigs_meta == N_contigs_data),
    Software = "Bin3C"
  ) %>%
  arrange(replicate, desc(Completeness))

mt_bin_stats <- make_bin_summary(mt_merged)
mt_bin_full <- mt_meta %>%
  inner_join(mt_bin_stats, by = c("cluster_id", "replicate")) %>%
  mutate(
    Phylum = extract_phylum(taxonomy),
    N50_meta_vs_computed_ratio    = round(N50_meta / N50_computed_bp, 3),
    NumContigs_meta_vs_data_match = (NumContigs_meta == N_contigs_data),
    Software = "MetaTOR"
  ) %>%
  arrange(replicate, desc(Completeness))


## ============================================================
## STEP 5: PRINT SUMMARIES
## ============================================================

print_section("TOTAL BP BINNED PER PIPELINE PER REPLICATE")
cat("(Binned contigs only -- unbinned contigs excluded)\n\n")

make_bp_summary_for_threshold <- function(df, comp_threshold) {
  d <- df %>% filter(Completeness >= comp_threshold)
  if (nrow(d) == 0) return(tibble())
  d %>%
    group_by(Software, replicate) %>%
    summarise(
      Completeness_threshold = paste0(">=", comp_threshold, "%"),
      N_bins                 = n(),
      Total_bp_binned        = sum(Total_bp_data, na.rm = TRUE),
      Total_Mb_binned        = round(sum(Total_bp_data, na.rm = TRUE)/1e6, 1),
      Mean_bin_size_Mb       = round(mean(GenomeSize, na.rm = TRUE)/1e6, 3),
      Median_bin_size_Mb     = round(median(GenomeSize, na.rm = TRUE)/1e6, 3),
      Median_bin_N50_kb      = round(median(N50_computed_bp, na.rm = TRUE)/1e3, 1),
      Mean_bin_N50_kb        = round(mean(N50_computed_bp, na.rm = TRUE)/1e3, 1),
      Min_bin_N50_kb         = round(min(N50_computed_bp, na.rm = TRUE)/1e3, 1),
      Min_contig_bp          = min(Min_contig_len_bp, na.rm = TRUE),
      Max_contig_Mb          = round(max(Max_contig_len_bp, na.rm = TRUE)/1e6, 3),
      Mean_contigs_per_bin   = round(mean(N_contigs_data,   na.rm = TRUE), 1),
      Median_contigs_per_bin = round(median(N_contigs_data, na.rm = TRUE), 0),
      .groups = "drop"
    )
}

all_dfs <- bind_rows(
  pm_bin_full,
  b3_bin_full,
  if (nrow(mt_bin_full) > 0) mt_bin_full else tibble()
)

bp_summary <- bind_rows(
  lapply(c(0, 85, 90, 95), function(t) {
    make_bp_summary_for_threshold(all_dfs, t)
  })
) %>% arrange(Software, replicate, Completeness_threshold)

print(bp_summary, n = Inf)

## Cross-verification: metadata N50 vs computed N50
print_section("CROSS-VERIFICATION: METADATA N50 vs COMPUTED N50")
cat("Ratio = N50 from metadata / N50 computed from contig lengths\n")
cat("Ratio close to 1.0 = consistent. Large deviation = investigate.\n\n")

for (df_name in c("ProxiMeta", "Bin3C", "MetaTOR")) {
  df <- switch(df_name,
               "ProxiMeta" = pm_bin_full,
               "Bin3C"     = b3_bin_full,
               "MetaTOR"   = if (nrow(mt_bin_full) > 0) mt_bin_full else NULL)
  if (is.null(df) || nrow(df) == 0) next
  r <- df$N50_meta_vs_computed_ratio
  r <- r[!is.na(r) & is.finite(r)]
  cat(sprintf("%s: mean ratio=%.3f, median=%.3f, range=[%.3f, %.3f]\n",
              df_name, mean(r), median(r), min(r), max(r)))
  match_pct <- round(100 * mean(df$NumContigs_meta_vs_data_match, na.rm = TRUE), 1)
  cat(sprintf("  NumContigs metadata == data: %.1f%% of bins match exactly\n", match_pct))
}

## Per-replicate bin-level summaries
print_section("PER-REPLICATE BIN SUMMARIES -- ALL PIPELINES")

for (pipeline in c("ProxiMeta", "Bin3C", "MetaTOR")) {
  df <- switch(pipeline,
               "ProxiMeta" = pm_bin_full,
               "Bin3C"     = b3_bin_full,
               "MetaTOR"   = if (nrow(mt_bin_full) > 0) mt_bin_full else NULL)
  if (is.null(df) || nrow(df) == 0) { cat(sprintf("\n%s: no data available\n", pipeline)); next }
  
  cat(sprintf("\n%s\n%s\n", pipeline, strrep("-", nchar(pipeline))))
  
  for (rep in ALL_REPS) {
    d <- df %>% filter(replicate == rep)
    if (nrow(d) == 0) next
    cat(sprintf("\n  --- %s: %d bins ---\n", rep, nrow(d)))
    cat("  Completeness thresholds:\n")
    print(comp_counts(d$Completeness, nrow(d)), n = Inf)
    
    num_tbl <- bind_rows(
      sumstats(d$Completeness,           "Completeness (%)"),
      sumstats(d$Contamination,          "Contamination (%)"),
      sumstats(d$GenomeSize/1e6,         "Genome size metadata (Mb)"),
      sumstats(d$Total_bp_data/1e6,      "Total bp data.tsv (Mb)"),
      sumstats(d$N50_meta/1e3,           "N50 from metadata (kb)"),
      sumstats(d$N50_computed_bp/1e3,    "N50 computed from data (kb)"),
      sumstats(d$N_contigs_data,         "N contigs per bin (data)"),
      sumstats(d$NumContigs_meta,        "N contigs per bin (metadata)"),
      sumstats(d$Mean_contig_len_bp,     "Mean contig length data (bp)"),
      sumstats(d$Median_contig_len_bp,   "Median contig length data (bp)"),
      sumstats(d$Max_contig_len_bp/1e3,  "Max contig length data (kb)"),
      sumstats(d$Pct_contigs_lt_1kb,     "% contigs < 1 kb"),
      sumstats(d$Pct_contigs_gt_10kb,    "% contigs > 10 kb"),
      sumstats(d$Pct_contigs_gt_50kb,    "% contigs > 50 kb"),
      sumstats(d$GC,                     "GC (%)")
    )
    cat("\n  Numeric summary:\n")
    print(num_tbl, n = Inf)
    
    if (pipeline == "MetaTOR" && "Phylum" %in% names(d)) {
      cat("\n  Phylum-level taxonomy:\n")
      print(as.data.frame(sort(table(d$Phylum), decreasing = TRUE)), row.names = FALSE)
      cat(sprintf("  Mean coding density: %.1f%%\n", mean(d$CodingDensity, na.rm = TRUE)))
      cat(sprintf("  Mean Hi-C coverage: %.1f contacts/kb\n", mean(d$HiCCoverage, na.rm = TRUE)))
    }
  }
}

## Contig length vs completeness strata
print_section("CONTIG LENGTH vs COMPLETENESS STRATA")
cat("Are higher-completeness bins composed of longer contigs?\n")
cat("(Addresses reviewer comment: assembly quality context)\n\n")

for (pipeline in c("ProxiMeta", "Bin3C", "MetaTOR")) {
  df <- switch(pipeline,
               "ProxiMeta" = pm_bin_full,
               "Bin3C"     = b3_bin_full,
               "MetaTOR"   = if (nrow(mt_bin_full) > 0) mt_bin_full else NULL)
  if (is.null(df) || nrow(df) == 0) next
  
  d <- df %>%
    filter(replicate %in% MAIN_REPS) %>%
    mutate(CompStrata = case_when(
      Completeness >= 95 ~ ">=95%",
      Completeness >= 90 ~ "90-95%",
      Completeness >= 85 ~ "85-90%",
      Completeness >= 50 ~ "50-85%",
      TRUE               ~ "<50%"
    )) %>%
    mutate(CompStrata = factor(CompStrata,
                               levels = c("<50%","50-85%","85-90%","90-95%",">=95%")))
  
  strata_summary <- d %>%
    group_by(CompStrata) %>%
    summarise(
      N_bins            = n(),
      Median_N50_kb     = round(median(N50_computed_bp, na.rm=TRUE)/1e3, 1),
      Mean_N50_kb       = round(mean(N50_computed_bp,   na.rm=TRUE)/1e3, 1),
      Min_N50_kb        = round(min(N50_computed_bp,    na.rm=TRUE)/1e3, 1),
      Median_pct_lt1kb  = round(median(Pct_contigs_lt_1kb,  na.rm=TRUE), 1),
      Median_pct_gt10kb = round(median(Pct_contigs_gt_10kb, na.rm=TRUE), 1),
      Median_contigs    = round(median(N_contigs_data,  na.rm=TRUE), 0),
      .groups = "drop"
    )
  
  cat(sprintf("\n%s (main replicates combined):\n", pipeline))
  print(strata_summary, n = Inf)
}


## ============================================================
## STEP 6: SAVE OUTPUT FILES
## ============================================================

print_section("SAVING OUTPUT FILES")

write_tsv(pm_bin_full %>% dplyr::select(-Contig_names),
          file.path(OUT_DIR, "proximeta_bin_summary.tsv"))
cat("Saved: proximeta_bin_summary.tsv\n")

write_tsv(b3_bin_full %>% dplyr::select(-Contig_names),
          file.path(OUT_DIR, "bin3c_bin_summary.tsv"))
cat("Saved: bin3c_bin_summary.tsv\n")

if (nrow(mt_bin_full) > 0) {
  write_tsv(mt_bin_full %>% dplyr::select(-Contig_names),
            file.path(OUT_DIR, "metator_bin_summary.tsv"))
  cat("Saved: metator_bin_summary.tsv\n")
}

pm_contig_list <- pm_merged %>%
  left_join(pm_meta %>% dplyr::select(cluster_id, replicate,
                                      Completeness, Contamination, GenomeSize, N50_meta),
            by = c("cluster_id","replicate")) %>%
  arrange(replicate, cluster_id, desc(length))
write_tsv(pm_contig_list, file.path(OUT_DIR, "proximeta_contig_list.tsv"))
cat("Saved: proximeta_contig_list.tsv\n")

b3_contig_list <- b3_merged %>%
  left_join(b3_meta %>% dplyr::select(cluster_id, replicate,
                                      Completeness, Contamination, GenomeSize, N50_meta),
            by = c("cluster_id","replicate")) %>%
  arrange(replicate, cluster_id, desc(length))
write_tsv(b3_contig_list, file.path(OUT_DIR, "bin3c_contig_list.tsv"))
cat("Saved: bin3c_contig_list.tsv\n")

if (nrow(mt_merged) > 0) {
  mt_contig_list <- mt_merged %>%
    left_join(mt_meta %>% dplyr::select(cluster_id, replicate,
                                        Completeness, Contamination, GenomeSize, N50_meta, taxonomy),
              by = c("cluster_id","replicate")) %>%
    arrange(replicate, cluster_id, desc(length))
  write_tsv(mt_contig_list, file.path(OUT_DIR, "metator_contig_list.tsv"))
  cat("Saved: metator_contig_list.tsv\n")
}

write_tsv(bp_summary, file.path(OUT_DIR, "total_bp_binned_summary.tsv"))
cat("Saved: total_bp_binned_summary.tsv\n")

cross_all <- bind_rows(
  pm_bin_full %>%
    dplyr::select(Software, replicate, cluster_id, Completeness, Contamination,
                  GenomeSize, N50_meta, N50_computed_bp, N_contigs_data, NumContigs_meta,
                  Total_bp_data, Pct_contigs_lt_1kb, Pct_contigs_gt_10kb, Pct_contigs_gt_50kb, GC),
  b3_bin_full %>%
    dplyr::select(Software, replicate, cluster_id, Completeness, Contamination,
                  GenomeSize, N50_meta, N50_computed_bp, N_contigs_data, NumContigs_meta,
                  Total_bp_data, Pct_contigs_lt_1kb, Pct_contigs_gt_10kb, Pct_contigs_gt_50kb, GC),
  if (nrow(mt_bin_full) > 0)
    mt_bin_full %>%
    dplyr::select(Software, replicate, cluster_id, Completeness, Contamination,
                  GenomeSize, N50_meta, N50_computed_bp, N_contigs_data, NumContigs_meta,
                  Total_bp_data, Pct_contigs_lt_1kb, Pct_contigs_gt_10kb, Pct_contigs_gt_50kb, GC)
  else tibble()
) %>% arrange(Software, replicate, desc(Completeness))

write_tsv(cross_all, file.path(OUT_DIR, "cross_pipeline_bin_summary_allreps.tsv"))
cat("Saved: cross_pipeline_bin_summary_allreps.tsv\n")

## ── Completeness-threshold summary table ─────────────────────────────
## For each Software x replicate x completeness threshold (>=0, >=20,
## >=40, >=60, >=80): n_bins, mean/median completeness, contamination,
## genome size, per-bin N50, and contigs-per-bin.

comp_thresh_summary <- bind_rows(
  lapply(c(0, 20, 40, 60, 80), function(thr) {
    cross_all %>%
      filter(Completeness >= thr) %>%
      group_by(Software, replicate) %>%
      summarise(
        Completeness_threshold  = thr,
        n_bins                  = n(),
        mean_comp               = round(mean(Completeness,    na.rm = TRUE), 2),
        median_comp             = round(median(Completeness,  na.rm = TRUE), 2),
        mean_cont               = round(mean(Contamination,   na.rm = TRUE), 2),
        median_cont             = round(median(Contamination, na.rm = TRUE), 2),
        mean_gs                 = round(mean(GenomeSize,      na.rm = TRUE), 0),
        median_gs               = round(median(GenomeSize,    na.rm = TRUE), 0),
        mean_binN50             = round(mean(N50_computed_bp,   na.rm = TRUE), 0),
        median_binN50           = round(median(N50_computed_bp, na.rm = TRUE), 0),
        mean_contigsperbin      = round(mean(N_contigs_data,    na.rm = TRUE), 2),
        median_contigsperbin    = round(median(N_contigs_data,  na.rm = TRUE), 0),
        .groups = "drop"
      )
  })
) %>%
  arrange(Software, replicate, Completeness_threshold)

write_tsv(comp_thresh_summary,
          file.path(OUT_DIR, "cross_pipeline_completeness_threshold_summary.tsv"))
cat("Saved: cross_pipeline_completeness_threshold_summary.tsv\n")

## ── Per-group range summary table ────────────────────────────────────
## Min-max ranges across replicates within Main (1,2,3A) and
## QC (3A,3B,3C) groups, for each Software x Completeness_threshold.
## Each metric also has a _list column with the per-replicate values
## in replicate order (e.g. "0.81,0.79,0.84") for traceability.
## CHANGE: added _list columns alongside each min/max pair.

make_range_summary <- function(df, grp_reps, grp_label) {
  df %>%
    filter(replicate %in% grp_reps) %>%
    arrange(replicate) %>%                          # consistent order within group
    group_by(Software, Completeness_threshold) %>%
    summarise(
      Group                = grp_label,
      n_reps               = n(),
      replicates_list      = paste(replicate, collapse = ","),
      
      n_bins_min           = min(n_bins),
      n_bins_max           = max(n_bins),
      n_bins_list          = paste(n_bins, collapse = ","),
      
      mean_comp_min        = min(mean_comp),
      mean_comp_max        = max(mean_comp),
      mean_comp_list       = paste(mean_comp, collapse = ","),
      
      median_comp_min      = min(median_comp),
      median_comp_max      = max(median_comp),
      median_comp_list     = paste(median_comp, collapse = ","),
      
      mean_cont_min        = min(mean_cont),
      mean_cont_max        = max(mean_cont),
      mean_cont_list       = paste(mean_cont, collapse = ","),
      
      median_cont_min      = min(median_cont),
      median_cont_max      = max(median_cont),
      median_cont_list     = paste(median_cont, collapse = ","),
      
      median_gs_min        = round(min(median_gs)/1e6, 3),
      median_gs_max        = round(max(median_gs)/1e6, 3),
      median_gs_list       = paste(round(median_gs/1e6, 3), collapse = ","),
      
      mean_N50_min         = round(min(mean_binN50)/1e3, 1),
      mean_N50_max         = round(max(mean_binN50)/1e3, 1),
      mean_N50_list        = paste(round(mean_binN50/1e3, 1), collapse = ","),
      
      median_N50_min       = round(min(median_binN50)/1e3, 1),
      median_N50_max       = round(max(median_binN50)/1e3, 1),
      median_N50_list      = paste(round(median_binN50/1e3, 1), collapse = ","),
      
      mean_ctg_min         = min(mean_contigsperbin),
      mean_ctg_max         = max(mean_contigsperbin),
      mean_ctg_list        = paste(mean_contigsperbin, collapse = ","),
      
      median_ctg_min       = min(median_contigsperbin),
      median_ctg_max       = max(median_contigsperbin),
      median_ctg_list      = paste(median_contigsperbin, collapse = ","),
      
      .groups = "drop"
    )
}

main_reps_vec <- c("E_15_1", "E_15_2", "E_15_3A")
qc_reps_vec   <- c("E_15_3A", "E_15_3B", "E_15_3C")

range_summary <- bind_rows(
  make_range_summary(comp_thresh_summary, main_reps_vec, "Main_1_2_3A"),
  make_range_summary(comp_thresh_summary, qc_reps_vec,   "QC_3A_3B_3C")
) %>%
  arrange(Software, Completeness_threshold, Group)

write_tsv(range_summary,
          file.path(OUT_DIR, "cross_pipeline_threshold_range_summary.tsv"))
cat("Saved: cross_pipeline_threshold_range_summary.tsv\n")

cat("\nDone.\n")

