#######################################################################
## plot_MedianDissimilarity_figures_perbinningpipeline.R
##
## Generates four figures of median dissimilarity of NN pairs across
## completeness levels, one panel per software (ProxiMeta, Bin3C, MetaTOR).
##
## CAN BE RUN IN TWO WAYS:
##
##   (A) STANDALONE — run this script directly.
##       Set dir_out below, then run the whole script.
##       It reads NN_[Software]_summary.tsv files from dir_out itself.
##
##   (B) SOURCED from DistancebyLength.R — the all_software_r_data_summary
##       object is already in the environment and dir_out is already set,
##       so the data-reading block below is skipped automatically.
##
## REQUIRES (standalone mode):
##   dir_out — path to DistancebyLength output folder (set below)
##   NN_Proximeta_summary.tsv, NN_Bin3C_summary.tsv,
##   NN_MetaTOR_summary.tsv  — produced by DistancebyLength.R / saveNNsummary()
##
## OUTPUTS (saved to dir_out):
##   fig_all20_pairs_perbinningpipeline.png
##   fig_main_6pairs_perbinningpipeline.png
##   fig_qc_abc_6pairs_perbinningpipeline.png
##   fig_no_abc_14pairs_perbinningpipeline.png
##
## COLOR ENCODING STRATEGY:
##   Related pairs share a hue family so forward/reverse comparisons
##   between the same two replicates are visually grouped:
##     Replicate 1 <-> Replicate 2   : red/orange family
##     Replicate 1 <-> Replicate 3A  : blue/cyan family
##     Replicate 2 <-> Replicate 3A  : green/lime family
##     Replicate 1/2 <-> Replicate 3B: purple/brown family
##     Replicate 1/2 <-> Replicate 3C: pink/teal family
##     3A/3B/3C internal : dark/neutral family
########################################################################

library(ggplot2)
library(dplyr)
library(readr)
library(cowplot)

## remove cowplot grey background from the plot_grid canvas
theme_set(theme_cowplot(font_size = 10))


# -----------------------------------------------------------------------
# USER PARAMETER (standalone mode only)
# Set this path if running the script directly rather than sourcing it
# from DistancebyLength.R. When sourced, dir_out is already set.
# -----------------------------------------------------------------------
## dir_out <- "~/path/to/DistancebyLength/"
# -----------------------------------------------------------------------


## ── 0. Build all_software_r_data_summary if not already in environment ──
##
## When sourced from DistancebyLength.R, all_software_r_data_summary
## already exists and this block is skipped.
## When run standalone, this block reads the NN_[Software]_summary.tsv
## files produced by saveNNsummary() and assembles the object.

## Check dir_out is set before proceeding
if (!exists("dir_out")) {
  stop("dir_out is not set.",
       "
Set it at the top of this script before running:",
       "
  dir_out <- "~/path/to/DistancebyLength/"")
}

## rename_pairs must be defined before the data-reading block uses it
rename_pairs <- function(x) {
  gsub("E 15 ", "Replicate ", x)
}

if (!exists("all_software_r_data_summary")) {

  message("all_software_r_data_summary not found -- reading from NN summary files...")

  ## Mirrors exactly the id_tag=="Software" block in DistancebyLength.R
  ## NN_[Software]_summary.tsv files must already exist in dir_out subfolders

  folders <- c("Proximeta", "Bin3C", "MetaTOR")

  all_software_r_data_list         <- list()
  all_software_r_data_summary_list <- list()

  for (rr in seq_along(folders)) {

    path   <- paste0(dir_out, folders[rr], "/")
    fname  <- paste0(path, "NN_", folders[rr], "_summary.tsv")

    if (!file.exists(fname)) {
      stop("Required file not found: ", fname,
           "
Run DistancebyLength.R first to generate NN summary files.")
    }

    r_data <- read_tsv(fname, show_col_types = FALSE)

    ## normalise: underscores -> spaces (matches DistancebyLength.R behaviour)
    r_data$GenomeSize  <- gsub("K", "k", r_data$GenomeSize)
    r_data$`NN Pairs` <- gsub("_", " ", r_data$`NN Pairs`)
    r_data$Software    <- folders[rr]

    all_software_r_data_list[[rr]] <- r_data

    r_data_summary <- r_data %>%
      group_by(Software, NN_name, `NN Pairs`, GenomeSize,
               `Completeness>=`, `Contamination<=`) %>%
      summarise(
        `Mean Dissimilarity Value of NN Pairs`   = mean(`NN Value`,   na.rm = TRUE),
        `Median Dissimilarity Value of NN Pairs` = median(`NN Value`, na.rm = TRUE),
        .groups = "drop"
      )

    all_software_r_data_summary_list[[rr]] <- r_data_summary
  }

  ## bind_rows guarantees a tibble regardless of R version
  all_software_r_data         <- dplyr::bind_rows(all_software_r_data_list)
  all_software_r_data_summary <- dplyr::bind_rows(all_software_r_data_summary_list)

  message("all_software_r_data_summary assembled: ",
          nrow(all_software_r_data_summary), " rows.")
}

## ── 1. Color palette — grouped by replicate pair family ───────────────

pair_colors_20 <- c(
  ## rep 1 <-> rep 2 : red/orange family
  "Replicate 1 to Replicate 2"    = "#E6194B",   ## vivid red
  "Replicate 2 to Replicate 1"    = "#F58231",   ## orange
  
  ## rep 1 <-> rep 3A : blue/cyan family
  "Replicate 1 to Replicate 3A"   = "#4363D8",   ## vivid blue
  "Replicate 3A to Replicate 1"   = "#42D4F4",   ## cyan
  
  ## rep 2 <-> rep 3A : green/lime family
  "Replicate 2 to Replicate 3A"   = "#3CB44B",   ## vivid green
  "Replicate 3A to Replicate 2"   = "#BFEF45",   ## lime
  
  ## rep 1/2 <-> rep 3B : purple/brown family
  "Replicate 1 to Replicate 3B"   = "#911EB4",   ## purple
  "Replicate 3B to Replicate 1"   = "#DCBEFF",   ## lavender
  "Replicate 2 to Replicate 3B"   = "#9A6324",   ## brown
  "Replicate 3B to Replicate 2"   = "#FFD700",   ## gold
  
  ## rep 1/2 <-> rep 3C : pink/teal family
  "Replicate 1 to Replicate 3C"   = "#F032E6",   ## magenta
  "Replicate 3C to Replicate 1"   = "#FABED4",   ## pink
  "Replicate 2 to Replicate 3C"   = "#469990",   ## teal
  "Replicate 3C to Replicate 2"   = "#AAFFC3",   ## mint
  
  ## internal 3A/3B/3C QC pairs : dark/neutral family
  "Replicate 3A to Replicate 3B"  = "#000075",   ## navy
  "Replicate 3B to Replicate 3A"  = "#808000",   ## olive
  "Replicate 3A to Replicate 3C"  = "#800000",   ## maroon
  "Replicate 3C to Replicate 3A"  = "#A9A9A9",   ## gray
  "Replicate 3B to Replicate 3C"  = "#000000",   ## black
  "Replicate 3C to Replicate 3B"  = "#E6BEFF"    ## lilac
)

## ── 2. Add renamed column to assembled summary ───────────────────────────────

## Apply rename_pairs to produce short display labels ("Rep X" format)
## This works identically whether the object came from DistancebyLength.R
## or was built by the standalone reading block above, since both leave
## NN Pairs in "E 15 X to E 15 Y" format at this point.
all_software_r_data_summary$`NN Pairs renamed` <- rename_pairs(
  as.character(all_software_r_data_summary$`NN Pairs`)
)

## ── 3. Define the four pair sets ─────────────────────────────────────

set_all20 <- rename_pairs(c(
  "E 15 1 to E 15 2",   "E 15 2 to E 15 1",
  "E 15 1 to E 15 3A",  "E 15 3A to E 15 1",
  "E 15 2 to E 15 3A",  "E 15 3A to E 15 2",
  "E 15 1 to E 15 3B",  "E 15 3B to E 15 1",
  "E 15 1 to E 15 3C",  "E 15 3C to E 15 1",
  "E 15 2 to E 15 3B",  "E 15 3B to E 15 2",
  "E 15 2 to E 15 3C",  "E 15 3C to E 15 2",
  "E 15 3A to E 15 3B", "E 15 3B to E 15 3A",
  "E 15 3A to E 15 3C", "E 15 3C to E 15 3A",
  "E 15 3B to E 15 3C", "E 15 3C to E 15 3B"
))

set_main_6 <- rename_pairs(c(
  "E 15 1 to E 15 2",  "E 15 2 to E 15 1",
  "E 15 1 to E 15 3A", "E 15 3A to E 15 1",
  "E 15 2 to E 15 3A", "E 15 3A to E 15 2"
))

set_abc_6 <- rename_pairs(c(
  "E 15 3A to E 15 3B", "E 15 3B to E 15 3A",
  "E 15 3A to E 15 3C", "E 15 3C to E 15 3A",
  "E 15 3B to E 15 3C", "E 15 3C to E 15 3B"
))

set_no_abc_14 <- rename_pairs(c(
  "E 15 1 to E 15 2",  "E 15 2 to E 15 1",
  "E 15 1 to E 15 3A", "E 15 3A to E 15 1",
  "E 15 2 to E 15 3A", "E 15 3A to E 15 2",
  "E 15 1 to E 15 3B", "E 15 3B to E 15 1",
  "E 15 1 to E 15 3C", "E 15 3C to E 15 1",
  "E 15 2 to E 15 3B", "E 15 3B to E 15 2",
  "E 15 2 to E 15 3C", "E 15 3C to E 15 2"
))

figure_sets <- list(
  all20    = list(pairs = set_all20,     filename = "fig_all20_pairs_perbinningpipeline"),
  main     = list(pairs = set_main_6,    filename = "fig_main_6pairs_perbinningpipeline"),
  qc_abc   = list(pairs = set_abc_6,     filename = "fig_qc_abc_6pairs_perbinningpipeline"),
  no_abc   = list(pairs = set_no_abc_14, filename = "fig_no_abc_14pairs_perbinningpipeline")
)

## ── 4. Panel-building function (y_breaks passed in) ───────────────────

make_panel <- function(data_s, pipeline_name, border_col,
                       pair_colors, y_max, y_breaks,
                       show_legend = FALSE) {
  
  ggplot(data_s,
         aes(x     = `Completeness>=`,
             y     = `Median Dissimilarity Value of NN Pairs`,
             color = `NN Pairs renamed`,
             group = `NN Pairs renamed`)) +
    geom_line(linewidth = 1.0) +
    geom_hline(yintercept = 0.05, linetype = "dashed",
               color = "black", linewidth = 0.6) +
    scale_color_manual(values = pair_colors, name = "Hi-C Library Pair") +
    scale_x_continuous(breaks = seq(0, 100, by = 5),
                       limits = c(0, 100)) +
    scale_y_continuous(breaks = y_breaks,
                       limits = c(0, y_max)) +
    xlab(expression("(%) Completeness" >= "")) +
    ylab("Median Dissimilarity Values of NN Pairs") +
    ggtitle(pipeline_name) +
    theme_bw() +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_blank(),
      strip.background   = element_blank(),
      panel.border       = element_rect(colour = "grey", fill = NA),
      plot.background    = element_rect(colour = border_col,
                                        fill = "white", linewidth = 2),
      plot.margin        = margin(t = 0.25, r = 0.25,
                                  b = 0.25, l = 0.25, unit = "cm"),
      axis.text.x        = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y        = element_text(size = 7),
      axis.title         = element_text(size = 10),
      plot.title         = element_text(size = 11, face = "bold"),
      legend.background  = element_rect(fill = "white", colour = NA),
      legend.key.size    = unit(0.4, "cm"),
      legend.key.height  = unit(0.6, "cm"),
      legend.spacing.y   = unit(0.15, "cm"),
      legend.text        = element_text(size = 8),
      legend.title       = element_text(size = 9, face = "bold"),
      legend.margin      = margin(0, 0, 0, 4),
      legend.position    = if (show_legend) "right" else "none"
    )
}

## ── 5. Loop over figure sets and save ────────────────────────────────

for (fig_name in names(figure_sets)) {
  
  fig_info   <- figure_sets[[fig_name]]
  pairs_keep <- fig_info$pairs
  col_subset <- pair_colors_20[names(pair_colors_20) %in% pairs_keep]
  
  data_fig <- all_software_r_data_summary %>%
    filter(`NN Pairs renamed` %in% pairs_keep) %>%
    mutate(`NN Pairs renamed` = factor(`NN Pairs renamed`,
                                       levels = pairs_keep))
  
  ## shared y upper limit across all three pipelines, rounded to 0.025
  ## minimum of 0.075 so the 0.05 reference line is always visible
  y_max_raw <- max(data_fig$`Median Dissimilarity Value of NN Pairs`,
                   na.rm = TRUE)
  y_max     <- max(ceiling(y_max_raw / 0.025) * 0.025, 0.075)

  ## adaptive y breaks: fine (0.025) when y_max <= 0.075, coarse (0.05) otherwise
  y_breaks  <- if (y_max <= 0.075) seq(0, y_max, by = 0.025) else seq(0, 1, by = 0.05)

  ## pull contamination and genome size labels from data for the title
  cont <- unique(data_fig$`Contamination<=`)
  gs   <- unique(data_fig$GenomeSize)
  
  data_prox    <- data_fig %>% filter(Software == "Proximeta")
  data_bin3c   <- data_fig %>% filter(Software == "Bin3C")
  data_metator <- data_fig %>% filter(Software == "MetaTOR")
  
  p_prox    <- make_panel(data_prox,    "ProxiMeta", "darkblue",
                          col_subset, y_max, y_breaks, show_legend = FALSE)
  p_bin3c   <- make_panel(data_bin3c,   "Bin3C",     "dodgerblue",
                          col_subset, y_max, y_breaks, show_legend = FALSE)
  p_metator <- make_panel(data_metator, "MetaTOR",   "chartreuse3",
                          col_subset, y_max, y_breaks, show_legend = FALSE)
  
  legend_grob <- suppressWarnings(cowplot::get_legend(
    make_panel(data_metator, "MetaTOR", "chartreuse3",
               col_subset, y_max, y_breaks, show_legend = TRUE)
  ))
  
  panels <- plot_grid(
    p_prox, p_bin3c, p_metator,
    nrow  = 3, ncol = 1,
    align = "v", axis = "lr"
  )
  
  panels_with_legend <- plot_grid(panels, legend_grob,
                                  ncol = 2, rel_widths = c(1, 0.28))
  
  ## figure title — matches style of original plots
  title_label <- paste0(
    "Median Dissimilarity Value of Nearest Neighbors (NN) Pairs:\n",
    "Contamination\u226410% Genome Size=", gs
  )
  
  title_grob <- ggdraw() +
    draw_label(title_label,
               fontface  = "plain",
               fontfamily = "Times New Roman",
               size      = 12,
               hjust     = 0.5,
               vjust     = 0.5)
  
  p_final <- plot_grid(title_grob, panels_with_legend,
                       ncol        = 1,
                       rel_heights = c(0.05, 1)) +
    theme(plot.background = element_rect(fill = "white", colour = NA))
  
  print(p_final)
  
  ggplot2::ggsave(
    filename = paste0(fig_info$filename, ".png"),
    plot     = p_final,
    width    = 8,
    height   = 10.4,
    path     = dir_out,
    bg       = "white"
  )
  
  message("Saved: ", fig_info$filename, ".png  -->  ", dir_out)
}

message("plot_MedianDissimilarity_figures_perbinningpipeline.R complete.")
