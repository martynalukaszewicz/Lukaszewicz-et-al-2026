########################################################################
## metadata_clustering_completeness_contamination_manuscript.R
##
## Generates completeness/contamination and genome size plots for
## MAGs recovered from wastewater by ProxiMeta, Bin3C, and MetaTOR,
## across Hi-C library replicates E_15_1, E_15_2, E_15_3A, E_15_3B,
## E_15_3C. Also validates and generates the extended summary table
## (Table 1 in the manuscript).
##
## OUTPUTS (PNG files saved to dir_plots):
##   [replicate]_Completeness>=0%_Contamination<=10%_Genome_Size=100kb-10Mb
##     _notitle_variablewidth_horizontallegend_nolegend.png
##     _notitle_variablewidth_horizontallegend_nolegend_gs.png
##     _notitle_variablewidth_horizontallegend_nolegend_correlations_gs_comp_cont.png
##     _nolegend.png
##
## Authors: Martyna Lukaszewicz, Benjamin J. Ridenhour,  Eva M. Top, Thibault Stalder
## University of Idaho
########################################################################


# -----------------------------------------------------------------------
# USER PARAMETERS — set these paths before running
# -----------------------------------------------------------------------

# Path to processed data TSV files (output of prepare_data_for_similarity_analysis.R)
dir_Similarity_Routput <- "~/path/to/Similarity_Routput/"

# Path to DistancebyLength folder (for sourcing getdata2_DistancebyLength.R)
dir_out <- "~/path/to/DistancebyLength/"

# Path to folder where this script and make_summary_table_extended.R live
dir_cur <- "~/path/to/Top_Lab/"

# Path where PNG plots will be saved
dir_plots <- "~/path/to/metadata_clustering_completeness_contamination_plots/"

# -----------------------------------------------------------------------
# END USER PARAMETERS
# -----------------------------------------------------------------------


rm(list = ls())

library(ggplot2)
library(dplyr)
library(readr)
library(tidyverse)
library(scales)
library(ggtext)

## Pin dplyr verbs to avoid masking by MASS or other packages
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

source(file.path(dir_out, "getdata2_DistancebyLength.R"))


# -----------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------

genomesize_cutoff    <- c("100Kb" = 100e3, "10Mb" = 10e6)
completeness_cutoff  <- c(0)       ## any completeness
contamination_cutoff <- c(10)      ## contamination <= 10%


# -----------------------------------------------------------------------
# Load metadata for all three pipelines
# -----------------------------------------------------------------------

load1 <- function(x) { read_tsv(x, col_names = TRUE, show_col_types = FALSE) }

## Bin3C
bin3c <- load1(file.path(dir_Similarity_Routput, "bin3c/metadata.tsv")) %>%
  dplyr::rename(Replicate = replicate) %>%
  dplyr::filter(extent >= genomesize_cutoff["100Kb"] &
                extent <= genomesize_cutoff["10Mb"]) %>%
  dplyr::arrange(Replicate, desc(Completeness))
Bin3C <- bin3c[, c("cluster_id", "Completeness", "Contamination", "Replicate")]

## ProxiMeta
proximeta <- load1(file.path(dir_Similarity_Routput, "proximeta/metadata.tsv")) %>%
  dplyr::rename(Contamination = marker_gene_overrepresentation,
                Completeness  = completeness,
                Replicate     = replicate) %>%
  dplyr::filter(genome_size >= genomesize_cutoff["100Kb"] &
                genome_size <= genomesize_cutoff["10Mb"]) %>%
  dplyr::arrange(Replicate, desc(Completeness))
Proximeta <- proximeta[, c("cluster_id", "Completeness", "Contamination", "Replicate")]

## MetaTOR
metator <- load1(file.path(dir_Similarity_Routput, "metator/metadata.tsv")) %>%
  dplyr::rename(Completeness  = completness,
                Contamination = contamination,
                Replicate     = replicate) %>%
  dplyr::filter(size >= genomesize_cutoff["100Kb"] &
                size <= genomesize_cutoff["10Mb"]) %>%
  dplyr::arrange(Replicate, desc(Completeness))
MetaTOR <- metator[, c("cluster_id", "Completeness", "Contamination", "Replicate")]

MetaTOR$Software   <- "MetaTOR"
Proximeta$Software <- "Proximeta"
Bin3C$Software     <- "Bin3C"

df_combined <- rbind(Proximeta, Bin3C, MetaTOR)


# -----------------------------------------------------------------------
# Attach GenomeSize from getdata2 output
# -----------------------------------------------------------------------

df_combined$GenomeSize <- NA_real_

## getdata2 with no contamination filter so all bins are present
data2        <- getdata2(contamination_cutoff = 10000,
                         genomesize_cutoff,
                         dir_Similarity_Routput)
data2_tibble <- as_tibble(data2)

soft <- c("Proximeta", "Bin3C", "MetaTOR")
reps <- c("E_15_1", "E_15_2", "E_15_3A", "E_15_3B", "E_15_3C")

empty_df <- as_tibble(data.frame(matrix(ncol = ncol(df_combined), nrow = 0)))
colnames(empty_df) <- colnames(df_combined)

for (ss in seq_along(soft)) {
  for (rr in seq_along(reps)) {

    data2_tibble_temp <- data2_tibble %>%
      dplyr::filter(Software == soft[ss], replicate == reps[rr])

    ## one row per cluster_id (take first occurrence)
    xx <- data2_tibble_temp %>%
      mutate(new = row_number()) %>%
      group_by(cluster_id) %>%
      summarise(Index = toString(new), .groups = "drop")
    data2_tibble_temp <- data2_tibble_temp[
      as.numeric(sapply(strsplit(xx$Index, ","), `[[`, 1)), ]

    df_combined_temp <- df_combined %>%
      dplyr::filter(Software == soft[ss], Replicate == reps[rr])

    for (nn in seq_len(nrow(df_combined_temp))) {
      nn_temp <- data2_tibble_temp %>%
        dplyr::filter(Software   == df_combined_temp$Software[nn],
                      replicate  == df_combined_temp$Replicate[nn],
                      cluster_id == df_combined_temp$cluster_id[nn])
      df_combined_temp$GenomeSize[nn] <- nn_temp$GenomeSize
    }
    empty_df <- rbind(empty_df, df_combined_temp)
  }
}

df_combined <- empty_df


# -----------------------------------------------------------------------
# Summary table (Table 1 extended)
# -----------------------------------------------------------------------

source(file.path(dir_cur, "make_summary_table_extended.R"))

validate_vs_table1(df_combined)

## HTML preview in RStudio Viewer:
make_summary_table_extended(df_combined)

## LaTeX for manuscript (uncomment to use):
# make_summary_table_extended(
#   df_combined,
#   output_format = "latex",
#   save_latex    = TRUE,
#   latex_path    = file.path(dir_cur, "Table1_extended.tex")
# )


# -----------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------

## helper function: combine columns to equal-length matrix columns
add.col <- function(df, new.col1, new.col2, new.col3, new.col4) {
  n.row <- dim(df)[1]
  length(new.col1) <- n.row
  length(new.col2) <- n.row
  length(new.col3) <- n.row
  length(new.col4) <- n.row
  cbind(df, new.col1, new.col2, new.col3, new.col4)
}

## shape mapping
shape_vals <- c(
  "ProxiMeta Completeness"  = 17,  ## filled triangle  (ProxiMeta)
  "Bin3C Completeness"      = 16,  ## filled circle    (Bin3C)
  "MetaTOR Completeness"    = 15,  ## filled square    (MetaTOR)
  "ProxiMeta Contamination" = 17,  ## filled triangle  (ProxiMeta)
  "Bin3C Contamination"     = 16,  ## filled circle    (Bin3C)
  "MetaTOR Contamination"   = 15   ## filled square    (MetaTOR)
)

shape_vals_gs <- c(
  "ProxiMeta Genome Size" = 17,  ## filled triangle (ProxiMeta)
  "Bin3C Genome Size"     = 16,  ## filled circle   (Bin3C)
  "MetaTOR Genome Size"   = 15   ## filled square   (MetaTOR)
)

name <- unique(df_combined$Replicate)

## fixed plot dimensions
pt_size     <- 1.0
pt_legend   <- 0.175
pt_size_gs  <- 1.0
text_size   <- 9.3
width_save  <- 4.2
height_save <- 1.8

lwd_size        <- 1.5
lwd_thick       <- 0.1
lwd_size_legend <- 3
lwd_thick_legend <- 0.2


## ── Main plotting loop ────────────────────────────────────────────────

for (i in seq_along(name)) {
  for (j in seq_along(completeness_cutoff)) {
    for (k in seq_along(contamination_cutoff)) {

      df_cutoff <- df_combined %>%
        dplyr::filter(Replicate     == name[i],
                      Completeness  >= completeness_cutoff[j],
                      Contamination <= contamination_cutoff[k])

      if (completeness_cutoff[j] < 0)
        df_cutoff <- df_cutoff %>% dplyr::filter(Completeness > 0)

      df_max_num_clust <- df_cutoff %>%
        add_count(Software) %>%
        dplyr::filter(n %in% tail(sort(unique(n)), 1)) %>%
        arrange(desc(n))
      max_clust <- df_max_num_clust$n[1]

      ## ── reshape into wide per-pipeline columns ──────────────────────

      Proximeta_cutoff <- df_cutoff %>%
        dplyr::filter(Software == "Proximeta") %>%
        dplyr::rename(Proximeta_cluster_id   = cluster_id,
                      Proximeta_Completeness = Completeness,
                      Proximeta_Contamination = Contamination,
                      Proximeta_GenomeSize   = GenomeSize) %>%
        as.matrix()

      Bin3C_cutoff <- df_cutoff %>%
        dplyr::filter(Software == "Bin3C") %>%
        dplyr::rename(Bin3C_cluster_id    = cluster_id,
                      Bin3C_Completeness  = Completeness,
                      Bin3C_Contamination = Contamination,
                      Bin3C_GenomeSize    = GenomeSize) %>%
        as.matrix()

      MetaTOR_cutoff <- df_cutoff %>%
        dplyr::filter(Software == "MetaTOR") %>%
        dplyr::rename(MetaTOR_cluster_id    = cluster_id,
                      MetaTOR_Completeness  = Completeness,
                      MetaTOR_Contamination = Contamination,
                      MetaTOR_GenomeSize    = GenomeSize) %>%
        as.matrix()

      df <- cbind(id = 1:max_clust)
      df <- add.col(df,
                    Proximeta_cutoff[, "Proximeta_cluster_id"],
                    Proximeta_cutoff[, "Proximeta_Completeness"],
                    Proximeta_cutoff[, "Proximeta_Contamination"],
                    Proximeta_cutoff[, "Proximeta_GenomeSize"])
      df <- add.col(df,
                    Bin3C_cutoff[, "Bin3C_cluster_id"],
                    Bin3C_cutoff[, "Bin3C_Completeness"],
                    Bin3C_cutoff[, "Bin3C_Contamination"],
                    Bin3C_cutoff[, "Bin3C_GenomeSize"])
      df <- add.col(df,
                    MetaTOR_cutoff[, "MetaTOR_cluster_id"],
                    MetaTOR_cutoff[, "MetaTOR_Completeness"],
                    MetaTOR_cutoff[, "MetaTOR_Contamination"],
                    MetaTOR_cutoff[, "MetaTOR_GenomeSize"])
      colnames(df) <- c(
        "id",
        "Proximeta_cluster_id", "Proximeta_Completeness",
        "Proximeta_Contamination", "Proximeta_GenomeSize",
        "Bin3C_cluster_id", "Bin3C_Completeness",
        "Bin3C_Contamination", "Bin3C_GenomeSize",
        "MetaTOR_cluster_id", "MetaTOR_Completeness",
        "MetaTOR_Contamination", "MetaTOR_GenomeSize"
      )

      df <- as_tibble(as.data.frame(df)) %>%
        mutate_at(c("Proximeta_Completeness", "Proximeta_Contamination",
                    "Proximeta_GenomeSize",
                    "Bin3C_Completeness", "Bin3C_Contamination",
                    "Bin3C_GenomeSize",
                    "MetaTOR_Completeness", "MetaTOR_Contamination",
                    "MetaTOR_GenomeSize"), as.numeric) %>%
        mutate_at("id", as.integer)

      df$Proximeta_GenomeSize <- df$Proximeta_GenomeSize / 1e6
      df$Bin3C_GenomeSize     <- df$Bin3C_GenomeSize     / 1e6
      df$MetaTOR_GenomeSize   <- df$MetaTOR_GenomeSize   / 1e6

      ## ── title string ────────────────────────────────────────────────

      title <- paste0(name[i],
                      "_Completeness>=0%",
                      "_Contamination<=", contamination_cutoff[k], "%",
                      "_Genome_Size=100kb-10Mb")

      ## ── mean/median annotations ─────────────────────────────────────

      mean_proximeta_comp   <- as.numeric(summary(df$Proximeta_Completeness)["Mean"])
      median_proximeta_comp <- as.numeric(summary(df$Proximeta_Completeness)["Median"])
      mean_bin3c_comp       <- as.numeric(summary(df$Bin3C_Completeness)["Mean"])
      median_bin3c_comp     <- as.numeric(summary(df$Bin3C_Completeness)["Median"])
      mean_metator_comp     <- as.numeric(summary(df$MetaTOR_Completeness)["Mean"])
      median_metator_comp   <- as.numeric(summary(df$MetaTOR_Completeness)["Median"])

      table_mean_median  <- c(mean_proximeta_comp, median_proximeta_comp,
                               mean_bin3c_comp,    median_bin3c_comp,
                               mean_metator_comp,  median_metator_comp)
      table_vals         <- c(5, 20, 35, 85, 100, 115)
      rank_mean_med_vals <- table_vals[rank(table_mean_median)]

      ## ── completeness/contamination plot (no legend) ─────────────────

      p_nolegend <- ggplot(df, aes(x = id)) +

        geom_hline(yintercept = mean_proximeta_comp,   linetype = "dotted", col = "darkblue",    lwd = lwd_thick) +
        annotate("richtext", x = rank_mean_med_vals[1], y = mean_proximeta_comp,   label = "AVG", col = "darkblue",    fill = alpha("white", 0.5), label.size = NA, label.padding = unit(0.08, "lines"), size = lwd_size) +
        geom_hline(yintercept = median_proximeta_comp, linetype = "dashed",  col = "darkblue",    lwd = lwd_thick) +
        annotate("richtext", x = rank_mean_med_vals[2], y = median_proximeta_comp, label = "MED", col = "darkblue",    fill = alpha("white", 0.5), label.size = NA, label.padding = unit(0.08, "lines"), size = lwd_size) +

        geom_hline(yintercept = mean_bin3c_comp,       linetype = "dotted", col = "dodgerblue",  lwd = lwd_thick) +
        annotate("richtext", x = rank_mean_med_vals[3], y = mean_bin3c_comp,       label = "AVG", col = "dodgerblue",  fill = alpha("white", 0.5), label.size = NA, label.padding = unit(0.08, "lines"), size = lwd_size) +
        geom_hline(yintercept = median_bin3c_comp,     linetype = "dashed",  col = "dodgerblue",  lwd = lwd_thick) +
        annotate("richtext", x = rank_mean_med_vals[4], y = median_bin3c_comp,     label = "MED", col = "dodgerblue",  fill = alpha("white", 0.5), label.size = NA, label.padding = unit(0.08, "lines"), size = lwd_size) +

        geom_hline(yintercept = mean_metator_comp,     linetype = "dotted", col = "chartreuse3", lwd = lwd_thick) +
        annotate("richtext", x = rank_mean_med_vals[5], y = mean_metator_comp,     label = "AVG", col = "chartreuse3", fill = alpha("white", 0.5), label.size = NA, label.padding = unit(0.08, "lines"), size = lwd_size) +
        geom_hline(yintercept = median_metator_comp,   linetype = "dashed",  col = "chartreuse3", lwd = lwd_thick) +
        annotate("richtext", x = rank_mean_med_vals[6], y = median_metator_comp,   label = "MED", col = "chartreuse3", fill = alpha("white", 0.5), label.size = NA, label.padding = unit(0.08, "lines"), size = lwd_size) +

        geom_point(aes(y = Proximeta_Completeness,  colour = "ProxiMeta Completeness",  shape = "ProxiMeta Completeness"),  size = pt_size, alpha = 0.7) +
        geom_point(aes(y = Bin3C_Completeness,      colour = "Bin3C Completeness",      shape = "Bin3C Completeness"),      size = pt_size, alpha = 0.7) +
        geom_point(aes(y = MetaTOR_Completeness,    colour = "MetaTOR Completeness",    shape = "MetaTOR Completeness"),    size = pt_size, alpha = 0.7) +
        geom_point(aes(y = Proximeta_Contamination, colour = "ProxiMeta Contamination", shape = "ProxiMeta Contamination"), size = pt_size, alpha = 0.7) +
        geom_point(aes(y = Bin3C_Contamination,     colour = "Bin3C Contamination",     shape = "Bin3C Contamination"),     size = pt_size, alpha = 0.7) +
        geom_point(aes(y = MetaTOR_Contamination,   colour = "MetaTOR Contamination",   shape = "MetaTOR Contamination"),   size = pt_size, alpha = 0.7) +

        scale_colour_manual(
          "Completeness/ \nContamination",
          values = c("ProxiMeta Completeness"  = "darkblue",
                     "Bin3C Completeness"      = "dodgerblue",
                     "MetaTOR Completeness"    = "chartreuse3",
                     "ProxiMeta Contamination" = "red3",
                     "Bin3C Contamination"     = "orange2",
                     "MetaTOR Contamination"   = "yellow1")) +
        scale_shape_manual("Completeness/ \nContamination", values = shape_vals) +
        xlab("Completeness Bin Rank Index") +
        ylab(expression(paste("Completeness/ \nContamination (%)"))) +
        theme(
          legend.title      = element_text(size = 10),
          legend.text       = element_text(size = 8),
          panel.grid.major  = element_blank(),
          panel.grid.minor  = element_blank(),
          panel.background  = element_blank(),
          legend.key        = element_rect(colour = NA, fill = NA),
          plot.margin       = unit(c(-0.015, -0.029, -0.01, 0.09), "inches"),
          text              = element_text(size = text_size),
          axis.title.y      = element_text(vjust = -2.0),
          legend.position   = "none"
        ) +
        scale_x_continuous(breaks = c(0, 25, 50, 75, 100, 125)) +
        scale_y_continuous(breaks = seq(from = 0, to = 100, by = 10))

      ## ── genome size plot (no legend) ────────────────────────────────

      df_gs_proximeta <- df[c("Proximeta_cluster_id", "Proximeta_Completeness",
                              "Proximeta_Contamination", "Proximeta_GenomeSize")] %>%
        as.data.frame() %>% arrange(desc(Proximeta_GenomeSize))
      df_gs_bin3c     <- df[c("Bin3C_cluster_id", "Bin3C_Completeness",
                              "Bin3C_Contamination", "Bin3C_GenomeSize")] %>%
        as.data.frame() %>% arrange(desc(Bin3C_GenomeSize))
      df_gs_metator   <- df[c("MetaTOR_cluster_id", "MetaTOR_Completeness",
                              "MetaTOR_Contamination", "MetaTOR_GenomeSize")] %>%
        as.data.frame() %>% arrange(desc(MetaTOR_GenomeSize))

      df_gs <- cbind(id_gs = 1:max_clust)
      df_gs <- add.col(df_gs,
                       df_gs_proximeta[, "Proximeta_cluster_id"],
                       df_gs_proximeta[, "Proximeta_Completeness"],
                       df_gs_proximeta[, "Proximeta_Contamination"],
                       df_gs_proximeta[, "Proximeta_GenomeSize"])
      df_gs <- add.col(df_gs,
                       df_gs_bin3c[, "Bin3C_cluster_id"],
                       df_gs_bin3c[, "Bin3C_Completeness"],
                       df_gs_bin3c[, "Bin3C_Contamination"],
                       df_gs_bin3c[, "Bin3C_GenomeSize"])
      df_gs <- add.col(df_gs,
                       df_gs_metator[, "MetaTOR_cluster_id"],
                       df_gs_metator[, "MetaTOR_Completeness"],
                       df_gs_metator[, "MetaTOR_Contamination"],
                       df_gs_metator[, "MetaTOR_GenomeSize"])
      colnames(df_gs) <- c(
        "id_gs",
        "Proximeta_cluster_id", "Proximeta_Completeness",
        "Proximeta_Contamination", "Proximeta_GenomeSize",
        "Bin3C_cluster_id", "Bin3C_Completeness",
        "Bin3C_Contamination", "Bin3C_GenomeSize",
        "MetaTOR_cluster_id", "MetaTOR_Completeness",
        "MetaTOR_Contamination", "MetaTOR_GenomeSize"
      )
      df_gs <- as_tibble(df_gs) %>%
        mutate_at(c("Proximeta_Completeness", "Proximeta_Contamination",
                    "Proximeta_GenomeSize",
                    "Bin3C_Completeness", "Bin3C_Contamination",
                    "Bin3C_GenomeSize",
                    "MetaTOR_Completeness", "MetaTOR_Contamination",
                    "MetaTOR_GenomeSize"), as.numeric) %>%
        mutate_at("id_gs", as.integer)

      mean_proximeta_gs   <- as.numeric(summary(df_gs$Proximeta_GenomeSize)["Mean"])
      median_proximeta_gs <- as.numeric(summary(df_gs$Proximeta_GenomeSize)["Median"])
      mean_bin3c_gs       <- as.numeric(summary(df_gs$Bin3C_GenomeSize)["Mean"])
      median_bin3c_gs     <- as.numeric(summary(df_gs$Bin3C_GenomeSize)["Median"])
      mean_metator_gs     <- as.numeric(summary(df_gs$MetaTOR_GenomeSize)["Mean"])
      median_metator_gs   <- as.numeric(summary(df_gs$MetaTOR_GenomeSize)["Median"])

      table_mean_median_gs  <- c(mean_proximeta_gs, median_proximeta_gs,
                                  mean_bin3c_gs,    median_bin3c_gs,
                                  mean_metator_gs,  median_metator_gs)
      table_vals_gs          <- c(0.1, 15, 29, 85, 100, 115)
      rank_mean_med_vals_gs  <- table_vals_gs[rank(table_mean_median_gs)]

      p_nolegend_gs <- ggplot(df_gs, aes(x = id_gs)) +

        geom_hline(yintercept = mean_proximeta_gs,   linetype = "dotted", col = "darkblue",    lwd = lwd_thick) +
        annotate("richtext", x = rank_mean_med_vals_gs[1], y = mean_proximeta_gs,   label = "AVG", col = "darkblue",    fill = alpha("white", 0.5), label.size = NA, label.padding = unit(0.08, "lines"), size = lwd_size) +
        geom_hline(yintercept = median_proximeta_gs, linetype = "dashed",  col = "darkblue",    lwd = lwd_thick) +
        annotate("richtext", x = rank_mean_med_vals_gs[2], y = median_proximeta_gs, label = "MED", col = "darkblue",    fill = alpha("white", 0.5), label.size = NA, label.padding = unit(0.08, "lines"), size = lwd_size) +

        geom_hline(yintercept = mean_bin3c_gs,       linetype = "dotted", col = "dodgerblue",  lwd = lwd_thick) +
        annotate("richtext", x = rank_mean_med_vals_gs[3], y = mean_bin3c_gs,       label = "AVG", col = "dodgerblue",  fill = alpha("white", 0.5), label.size = NA, label.padding = unit(0.08, "lines"), size = lwd_size) +
        geom_hline(yintercept = median_bin3c_gs,     linetype = "dashed",  col = "dodgerblue",  lwd = lwd_thick) +
        annotate("richtext", x = rank_mean_med_vals_gs[4], y = median_bin3c_gs,     label = "MED", col = "dodgerblue",  fill = alpha("white", 0.5), label.size = NA, label.padding = unit(0.08, "lines"), size = lwd_size) +

        geom_hline(yintercept = mean_metator_gs,     linetype = "dotted", col = "chartreuse3", lwd = lwd_thick) +
        annotate("richtext", x = rank_mean_med_vals_gs[5], y = mean_metator_gs,     label = "AVG", col = "chartreuse3", fill = alpha("white", 0.5), label.size = NA, label.padding = unit(0.08, "lines"), size = lwd_size) +
        geom_hline(yintercept = median_metator_gs,   linetype = "dashed",  col = "chartreuse3", lwd = lwd_thick) +
        annotate("richtext", x = rank_mean_med_vals_gs[6], y = median_metator_gs,   label = "MED", col = "chartreuse3", fill = alpha("white", 0.5), label.size = NA, label.padding = unit(0.08, "lines"), size = lwd_size) +

        geom_hline(yintercept = 0.1, linetype = "blank", col = "#888888") +
        annotate("richtext", x = 0.01, y = 0.1, label = "100 kb", col = "black", fill = alpha("white", 0.5), label.size = NA, label.padding = unit(0.08, "lines"), size = lwd_size) +

        geom_point(aes(y = Proximeta_GenomeSize, colour = "ProxiMeta Genome Size", shape = "ProxiMeta Genome Size"), size = pt_size_gs, alpha = 0.7) +
        geom_point(aes(y = Bin3C_GenomeSize,     colour = "Bin3C Genome Size",     shape = "Bin3C Genome Size"),     size = pt_size_gs, alpha = 0.7) +
        geom_point(aes(y = MetaTOR_GenomeSize,   colour = "MetaTOR Genome Size",   shape = "MetaTOR Genome Size"),   size = pt_size_gs, alpha = 0.7) +

        scale_colour_manual("Genome Size",
                            values = c("ProxiMeta Genome Size" = "darkblue",
                                       "Bin3C Genome Size"     = "dodgerblue",
                                       "MetaTOR Genome Size"   = "chartreuse3")) +
        scale_shape_manual("Genome Size", values = shape_vals_gs) +
        xlab("Genome Size Bin Rank Index") +
        ylab(expression(paste("Genome Size (Mb)"))) +
        theme(
          legend.title     = element_text(size = 10),
          legend.text      = element_text(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.key       = element_rect(colour = NA, fill = NA),
          plot.margin      = unit(c(-0.015, -0.029, -0.01, 0.09), "inches"),
          text             = element_text(size = text_size),
          axis.title.y     = element_text(vjust = -0.0),
          legend.position  = "none"
        ) +
        scale_x_continuous(breaks = c(0, 25, 50, 75, 100, 125)) +
        scale_y_continuous(
          breaks = c(0, 0.1, 1, 2, 3, 4, 5, 6, 7),
          labels = c("0", "", "1", "2", "3", "4", "5", "6", "7"),
          limits = c(0, 7))

      ## ── correlation plot (genome size vs completeness/contamination) ─

      ## median genome size per 20% completeness interval
      comp_low_cutoff  <- seq(from = 0,  to = 80, by = 20)
      comp_high_cutoff <- c(seq(from = 20, to = 80, by = 20), 100.1)
      cont_low_cutoff  <- seq(from = 0,  to = 8,  by = 2)
      cont_high_cutoff <- c(seq(from = 2, to = 8, by = 2), 10.1)

      med_comp_mat <- matrix(NA, nrow = length(comp_low_cutoff), ncol = length(soft),
                             dimnames = list(
                               paste0("comp_", comp_low_cutoff, "_", comp_high_cutoff),
                               soft))
      med_cont_mat <- matrix(NA, nrow = length(cont_low_cutoff), ncol = length(soft),
                             dimnames = list(
                               paste0("cont_", cont_low_cutoff, "_", cont_high_cutoff),
                               soft))

      for (sss in seq_along(soft)) {
        comp_vals <- as.numeric(na.omit(unlist(
          df[, grepl(paste0(soft[sss], "_Completeness"), colnames(df))])))
        cont_vals <- as.numeric(na.omit(unlist(
          df[, grepl(paste0(soft[sss], "_Contamination"), colnames(df))])))
        gs_vals   <- as.numeric(na.omit(unlist(
          df[, grepl(paste0(soft[sss], "_GenomeSize"), colnames(df))])))

        for (sd_row in seq_along(comp_low_cutoff)) {
          idx <- comp_vals >= comp_low_cutoff[sd_row] &
                 comp_vals <  comp_high_cutoff[sd_row]
          med_comp_mat[sd_row, sss] <- median(gs_vals[idx], na.rm = TRUE)
        }
        for (sd_row in seq_along(cont_low_cutoff)) {
          idx <- cont_vals >= cont_low_cutoff[sd_row] &
                 cont_vals <  cont_high_cutoff[sd_row]
          med_cont_mat[sd_row, sss] <- median(gs_vals[idx], na.rm = TRUE)
        }
      }

      ## expand to per-unit index for smooth step lines
      expand_to_idx <- function(mat, low_cut, high_cut, n_out) {
        out <- as_tibble(matrix(0, nrow = n_out, ncol = length(soft)))
        colnames(out) <- soft
        out$idx <- seq(from = 0, to = n_out - 1, by = 1)
        for (ccc in seq_along(high_cut)) {
          rows <- out$idx >= low_cut[ccc] & out$idx < high_cut[ccc]
          out[rows, seq_along(soft)] <- as.list(mat[ccc, ])
        }
        out
      }
      med_comp_mat_idx <- expand_to_idx(med_comp_mat, comp_low_cutoff, comp_high_cutoff, 101)
      med_cont_mat_idx <- expand_to_idx(med_cont_mat, cont_low_cutoff, cont_high_cutoff, 11)

      p_correlations_nolegend <- ggplot() +

        geom_point(df, mapping = aes(y = Proximeta_GenomeSize, x = Proximeta_Completeness,
                                     colour = "ProxiMeta Completeness",
                                     shape  = "ProxiMeta Completeness"),
                   size = pt_size * .4, alpha = 0.7) +
        geom_point(df, mapping = aes(y = Bin3C_GenomeSize, x = Bin3C_Completeness,
                                     colour = "Bin3C Completeness",
                                     shape  = "Bin3C Completeness"),
                   size = pt_size * .4, alpha = 0.7) +
        geom_point(df, mapping = aes(y = MetaTOR_GenomeSize, x = MetaTOR_Completeness,
                                     colour = "MetaTOR Completeness",
                                     shape  = "MetaTOR Completeness"),
                   size = pt_size * .4, alpha = 0.7) +
        geom_point(df, mapping = aes(y = Proximeta_GenomeSize, x = Proximeta_Contamination,
                                     colour = "ProxiMeta Contamination",
                                     shape  = "ProxiMeta Contamination"),
                   size = pt_size * .075, alpha = 0.7) +
        geom_point(df, mapping = aes(y = Bin3C_GenomeSize, x = Bin3C_Contamination,
                                     colour = "Bin3C Contamination",
                                     shape  = "Bin3C Contamination"),
                   size = pt_size * .075, alpha = 0.7) +
        geom_point(df, mapping = aes(y = MetaTOR_GenomeSize, x = MetaTOR_Contamination,
                                     colour = "MetaTOR Contamination",
                                     shape  = "MetaTOR Contamination"),
                   size = pt_size * .075, alpha = 0.7) +

        scale_colour_manual(
          "Completeness/\nContamination",
          values = c("ProxiMeta Completeness"  = "darkblue",
                     "Bin3C Completeness"      = "dodgerblue",
                     "MetaTOR Completeness"    = "chartreuse3",
                     "ProxiMeta Contamination" = "red3",
                     "Bin3C Contamination"     = "orange2",
                     "MetaTOR Contamination"   = "yellow1")) +
        scale_shape_manual("Completeness/\nContamination", values = shape_vals) +
        ylab("Genome Size (Mb)") +
        xlab(expression(paste("Completeness/Contamination (%)"))) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.key       = element_rect(colour = NA, fill = NA),
          plot.margin      = unit(c(-0.015, 0.01, -0.01, 0.0), "inches"),
          text             = element_text(size = text_size),
          axis.title.y     = element_text(vjust = 0.05),
          legend.position  = "none"
        ) +
        scale_x_continuous(breaks = seq(from = 0, to = 100, by = 10)) +
        scale_y_continuous(
          breaks = c(0, 0.1, 1, 2, 3, 4, 5, 6, 7),
          labels = c("0", "", "1", "2", "3", "4", "5", "6", "7"),
          limits = c(0, 7)) +
        geom_hline(yintercept = 0.1, linetype = "1F", col = "#888888") +
        annotate("richtext", x = 99, y = 0.1, label = "100 kb", col = "black",
                 fill = alpha("white", 0.5), label.size = NA,
                 label.padding = unit(0.08, "lines"), size = lwd_size) +
        ## median genome size step lines
        geom_line(med_comp_mat_idx,
                  mapping = aes(y = Proximeta, x = idx),
                  col = "darkblue",    linetype = "solid", lwd = lwd_thick * 7) +
        geom_line(med_comp_mat_idx,
                  mapping = aes(y = Bin3C, x = idx),
                  col = "dodgerblue",  linetype = "solid", lwd = lwd_thick * 7) +
        geom_line(med_comp_mat_idx,
                  mapping = aes(y = MetaTOR, x = idx),
                  col = "chartreuse3", linetype = "solid", lwd = lwd_thick * 7) +
        geom_line(med_cont_mat_idx,
                  mapping = aes(y = Proximeta, x = idx),
                  col = "red3",    linetype = "solid", lwd = lwd_thick * 4) +
        geom_line(med_cont_mat_idx,
                  mapping = aes(y = Bin3C, x = idx),
                  col = "orange2", linetype = "solid", lwd = lwd_thick * 4) +
        geom_line(med_cont_mat_idx,
                  mapping = aes(y = MetaTOR, x = idx),
                  col = "yellow1", linetype = "solid", lwd = lwd_thick * 4)

      ## ── save nolegend outputs ────────────────────────────────────────

      base <- paste0(gsub("%", "", title),
                     "_notitle_variablewidth_horizontallegend_nolegend")

      ggplot2::ggsave(
        filename = paste0(base, ".png"),
        plot     = p_nolegend,
        width    = width_save, height = height_save,
        path     = dir_plots)

      ggplot2::ggsave(
        filename = paste0(base, "_gs.png"),
        plot     = p_nolegend_gs,
        width    = width_save, height = height_save,
        path     = dir_plots)

      ggplot2::ggsave(
        filename = paste0(base, "_correlations_gs_comp_cont.png"),
        plot     = p_correlations_nolegend,
        width    = width_save, height = height_save,
        path     = dir_plots)

      base_nolegend <- paste0(gsub("%", "", title), "_nolegend")
      ggplot2::ggsave(
        filename = paste0(base_nolegend, ".png"),
        plot     = p_nolegend,
        width    = width_save, height = height_save,
        path     = dir_plots)

      message("Saved plots for: ", title)

    }  ## end k (contamination)
  }  ## end j (completeness)
}  ## end i (replicate)

message("metadata_clustering_completeness_contamination_manuscript.R complete.")
