########################################################################
## plot_MedianDissimilarity_figures_perreplicate.R
##
## Generates three figures of median dissimilarity of NN pairs across
## completeness levels, one panel per Hi-C library replicate.
## Replaces the boxplot Figure 4 with a line plot consistent with Fig 3.
##
## HOW TO SOURCE FROM YOUR MAIN SCRIPT:
##   source(file.path(dir_out, "plot_MedianDissimilarity_figures_perreplicate.R"))
##
## REQUIRES (must exist in environment before sourcing):
##   dir_out  — output directory path (string)
##
## OUTPUTS (saved to dir_out):
##   fig_main_3panels_perhic.png    — Replicate 1, Replicate 2, Replicate 3A
##   fig_qc_abc_3panels_perhic.png  — Replicate 3A, Replicate 3B, Replicate 3C
##   fig_all5panels_perhic.png      — Replicate 1, Replicate 2, Replicate 3A, Replicate 3B, Replicate 3C
##
## COLOR ENCODING STRATEGY:
##   Related pipeline pairs share a hue family:
##     Bin3C   <-> ProxiMeta : blue/cyan family
##     MetaTOR <-> ProxiMeta : red/orange family
##     MetaTOR <-> Bin3C     : green/olive family
########################################################################

library(ggplot2)
library(dplyr)
library(readr)
library(cowplot)

theme_set(theme_cowplot(font_size = 10))

TNR <- "Times New Roman"

## ── 1. Pipeline pair colors and order ────────────────────────────────

pairs_ord <- c(
  "Bin3C to ProxiMeta",
  "ProxiMeta to Bin3C",
  "MetaTOR to ProxiMeta",
  "ProxiMeta to MetaTOR",
  "MetaTOR to Bin3C",
  "Bin3C to MetaTOR"
)

pipeline_colors <- c(
  "Bin3C to ProxiMeta"    = "#4363D8",   ## blue
  "ProxiMeta to Bin3C"    = "#42D4F4",   ## cyan
  "MetaTOR to ProxiMeta"  = "#E6194B",   ## red
  "ProxiMeta to MetaTOR"  = "#F58231",   ## orange
  "MetaTOR to Bin3C"      = "#3CB44B",   ## green
  "Bin3C to MetaTOR"      = "#808000"    ## olive
)

## ── 2. Border colors per replicate ───────────────────────────────────

rep_border_colors <- c(
  "Replicate 1"  = "#3AAFA9",
  "Replicate 2"  = "#3AAFA9",
  "Replicate 3A" = "#3AAFA9",
  "Replicate 3B" = "#E8C84A",
  "Replicate 3C" = "#E8C84A"
)

## ── 3. Read and assemble replicate summary data ───────────────────────

folders_all <- c("E_15_1", "E_15_2", "E_15_3A", "E_15_3B", "E_15_3C")

all_reps_r_data_summary <- c()

for (rr in 1:length(folders_all)) {

  path   <- paste0(dir_out, folders_all[rr], "/")
  r_data <- read_tsv(paste0(path, "NN_", folders_all[rr], "_summary.tsv"),
                     show_col_types = FALSE)

  r_data <- r_data %>%
    rename(`Dissimilarity Value of NN Pairs` = `NN Value`)

  r_data$GenomeSize <- gsub("Kb", "kb", r_data$GenomeSize)
  r_data$`NN Pairs` <- gsub("Proximeta", "ProxiMeta", r_data$`NN Pairs`)

  ## keep only the six pipeline comparison pairs and set factor order
  r_data <- r_data %>%
    filter(`NN Pairs` %in% pairs_ord) %>%
    mutate(`NN Pairs` = factor(`NN Pairs`, levels = pairs_ord))

  ## display label for panel title
  r_data$Replicate <- gsub("E_15_", "Replicate ", folders_all[rr])
  r_data$Replicate <- gsub("_", "", r_data$Replicate)

  ## compute median per group
  r_data_summary <- r_data %>%
    group_by(Replicate, `NN Pairs`, GenomeSize,
             `Completeness>=`, `Contamination<=`) %>%
    summarise(
      `Median Dissimilarity Value of NN Pairs` =
        median(`Dissimilarity Value of NN Pairs`, na.rm = TRUE),
      .groups = "drop"
    )

  if (rr == 1) {
    all_reps_r_data_summary <- r_data_summary
  } else {
    all_reps_r_data_summary <- rbind(all_reps_r_data_summary,
                                      r_data_summary)
  }
}

message("all_reps_r_data_summary collected")

## ── 4. Figure set definitions ─────────────────────────────────────────

figure_sets <- list(
  main_3 = list(
    reps     = c("Replicate 1", "Replicate 2", "Replicate 3A"),
    filename = "fig_main_3panels_perhic"
  ),
  qc_abc = list(
    reps     = c("Replicate 3A", "Replicate 3B", "Replicate 3C"),
    filename = "fig_qc_abc_3panels_perhic"
  ),
  all_5 = list(
    reps     = c("Replicate 1", "Replicate 2", "Replicate 3A",
                 "Replicate 3B", "Replicate 3C"),
    filename = "fig_all5panels_perhic"
  )
)

## ── 5. Panel-building function ────────────────────────────────────────

make_rep_panel <- function(data_s, rep_name, border_col,
                           y_max, show_legend = FALSE) {

  ggplot(data_s,
         aes(x     = `Completeness>=`,
             y     = `Median Dissimilarity Value of NN Pairs`,
             color = `NN Pairs`,
             group = `NN Pairs`)) +
    geom_line(linewidth = 1.0) +
    geom_hline(yintercept = 0.05, linetype = "dashed",
               color = "black", linewidth = 0.6) +
    scale_color_manual(values = pipeline_colors,
                       name   = "Metagenomic Binning\nPipeline Pair") +
    scale_x_continuous(breaks = seq(0, 100, by = 5),
                       limits = c(0, 100)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.025),
                       limits = c(0, y_max)) +
    xlab(expression("(%) Completeness" >= "")) +
    ylab("Median dissimilarity of NN pairs") +
    ggtitle(rep_name) +
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
      legend.margin      = margin(0, 0, 0, 10),
      legend.position    = if (show_legend) "right" else "none"
    )
}

## ── 6. Loop over figure sets and save ────────────────────────────────

for (fig_name in names(figure_sets)) {

  fig_info  <- figure_sets[[fig_name]]
  reps_keep <- fig_info$reps

  data_fig <- all_reps_r_data_summary %>%
    filter(Replicate %in% reps_keep) %>%
    mutate(Replicate = factor(Replicate, levels = reps_keep))

  ## shared y upper limit across all replicates, rounded to 0.025
  y_max_raw <- max(data_fig$`Median Dissimilarity Value of NN Pairs`,
                   na.rm = TRUE)
  y_max     <- ceiling(y_max_raw / 0.025) * 0.025

  ## pull contamination and genome size for title
  cont <- unique(data_fig$`Contamination<=`)
  gs   <- unique(data_fig$GenomeSize)

  ## build panels without legend
  panel_list <- lapply(reps_keep, function(rep_name) {
    d <- data_fig %>% filter(Replicate == rep_name)
    make_rep_panel(d, rep_name,
                   border_col  = rep_border_colors[rep_name],
                   y_max       = y_max,
                   show_legend = FALSE)
  })

  ## extract legend from first panel built with legend on
  legend_grob <- suppressWarnings(
    cowplot::get_legend(
      make_rep_panel(
        data_fig %>% filter(Replicate == reps_keep[1]),
        reps_keep[1],
        border_col  = rep_border_colors[reps_keep[1]],
        y_max       = y_max,
        show_legend = TRUE
      )
    )
  )

  panels <- plot_grid(
    plotlist = panel_list,
    nrow     = length(reps_keep),
    ncol     = 1,
    align    = "v",
    axis     = "lr"
  )

  panels_with_legend <- plot_grid(panels, legend_grob,
                                  ncol = 2, rel_widths = c(1, 0.30))

  title_label <- paste0(
    "Median Dissimilarity Value of Nearest Neighbors (NN) Pairs:\n",
    "Contamination\u226410% Genome Size=", gs
  )

  title_grob <- ggdraw() +
    draw_label(title_label,
               fontface   = "plain",
               fontfamily = TNR,
               size       = 11,
               hjust      = 0.5,
               vjust      = 0.5)

  p_final <- plot_grid(title_grob, panels_with_legend,
                       ncol        = 1,
                       rel_heights = c(0.04, 1)) +
    theme(plot.background = element_rect(fill = "white", colour = NA))

  print(p_final)

  fig_height <- 3.5 * length(reps_keep) + 0.5

  ggplot2::ggsave(
    filename = paste0(fig_info$filename, ".png"),
    plot     = p_final,
    width    = 8,
    height   = fig_height,
    path     = dir_out,
    bg       = "white"
  )

  message("Saved: ", fig_info$filename, ".png  -->  ", dir_out)
}

message("plot_MedianDissimilarity_figures_perreplicate.R complete.")
