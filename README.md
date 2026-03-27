# Hi-C precision in bacterial genomes reconstruction

**Lukaszewicz et al. 2026**

Martyna Lukaszewicz, Benjamin J. Ridenhour, Eva M. Top, Thibault Stalder  
University of Idaho

---

## Overview

This repository contains all R scripts used for data processing, dissimilarity calculations, and figure generation in:

> Lukaszewicz M, Ridenhour BJ, Top EM, Stalder T. *Hi-C precision in bacterial genomes reconstruction.* (2026)

The study evaluates the reproducibility of Hi-C-based metagenomic binning pipelines (ProxiMeta, Bin3C v0.1.1, MetaTOR v1.1.5) using three Hi-C library replicates of a wastewater effluent sample, and quantifies the precision of metagenome-assembled genome (MAG) reconstruction using Jaccard distance.

---

## Data availability

Processed pipeline outputs underlying all analyses — including contig-to-bin assignment tables, per-bin CheckM quality metrics, pairwise Jaccard dissimilarity matrices, and summary statistics — are deposited at Zenodo:

> **DOI:** [10.5281/zenodo.19260186](https://doi.org/10.5281/zenodo.19260186)

Raw sequencing files (Hi-C and shotgun FASTQ reads) and the MEGAHIT shotgun assembly are no longer available, as these libraries were prepared and sequenced in 2020–2021 and files were not retained through subsequent storage transitions. The wastewater sample used in this study was collected from the same Moscow WWTP described in Stalder et al. 2019 (ISME J, BioProject PRJNA506462), but represents an independent sample collected on 15 November 2018 from the outflow of the facility. The shotgun and Hi-C libraries generated for the present study are separate sequencing runs not deposited under PRJNA506462.

---

## Repository structure

```
Lukaszewicz-et-al-2026/
├── 1_prepare_data/
│   ├── prepare_data_for_similarity_analysis.R
│   └── summarize_hic_metadata.R
├── 2_distance_calculations/
│   ├── DistancebyLength.R
│   ├── getdata2_DistancebyLength.R
│   ├── savesimMat_DistancebyLength.R
│   ├── savesummary_DistancebyLength.R
│   └── saveNNsummary_DistancebyLength.R
└── 3_figures_tables/
    ├── make_summary_table_extended.R
    ├── metadata_clustering_completeness_contamination_manuscript.R
    ├── plot_MedianDissimilarity_figures_perbinningpipeline.R
    └── plot_MedianDissimilarity_figures_perreplicate.R
```

---

## Script descriptions

### 1_prepare_data/

**`prepare_data_for_similarity_analysis.R`**  
Reads raw contig-to-bin assignment files and CheckM report files produced by ProxiMeta, Bin3C, and MetaTOR, and reformats them into standardised `data.tsv` (contig assignments with lengths) and `metadata.tsv` (per-bin CheckM quality metrics) files for each pipeline. These files are the primary inputs to all downstream scripts.

- Input: `contig_ids_per_clusters/` and `report/` directories
- Output: `Similarity_Routput/[pipeline]/data.tsv` and `metadata.tsv`

**`summarize_hic_metadata.R`**  
Reads `data.tsv` and `metadata.tsv` for all three pipelines, merges them, and produces comprehensive per-bin and per-contig summary statistics. Also generates cross-pipeline summary tables and the contig length stratification analysis used to contextualise assembly quality as a confounding factor.

- Input: `Similarity_Routput/[pipeline]/data.tsv` and `metadata.tsv`
- Output: `[pipeline]_bin_summary.tsv`, `[pipeline]_contig_list.tsv`, `cross_pipeline_bin_summary_allreps.tsv`, `total_bp_binned_summary.tsv`

### 2_distance_calculations/

**`DistancebyLength.R`**  
Main orchestration script. Sources the four helper functions below and runs the full Jaccard dissimilarity calculation pipeline. Run twice: once with `id_tag <- "replicate"` (cross-pipeline comparisons) and once with `id_tag <- "Software"` (cross-replicate comparisons).

**`getdata2_DistancebyLength.R`**  
Helper function `getdata2()`. Loads and merges `data.tsv` and `metadata.tsv` for all three pipelines, applies genome size and contamination filters, and returns the combined `data2` object used by `savesimMat()`.

**`savesimMat_DistancebyLength.R`**  
Helper function `savesimMat()`. Computes pairwise Jaccard dissimilarity matrices between MAG bins across pipeline or replicate pairs, and saves them as CSV files. Skips existing files to allow safe reruns.

**`savesummary_DistancebyLength.R`**  
Helper function `savesummary()`. Reads simMat CSV files, extracts nearest-neighbour dissimilarity values, and saves per-folder summary TSV files.

**`saveNNsummary_DistancebyLength.R`**  
Helper function `saveNNsummary()`. Removes NAs from summary files and saves cleaned `NN_[folder]_summary.tsv` files used for plotting.

### 3_figures_tables/

**`make_summary_table_extended.R`**  
Defines `make_summary_table_extended()` and `validate_vs_table1()`. Generates the extended summary table (Table 1 in the manuscript) across all five replicates and three pipelines, in HTML, LaTeX, or data frame format. Must be sourced before `metadata_clustering_completeness_contamination_manuscript.R`.

**`metadata_clustering_completeness_contamination_manuscript.R`**  
Generates completeness/contamination and genome size scatter plots for all five Hi-C library replicates. Also calls `validate_vs_table1()` and `make_summary_table_extended()` to produce Table 1.

- Sources: `getdata2_DistancebyLength.R`, `make_summary_table_extended.R`
- Input: `Similarity_Routput/[pipeline]/metadata.tsv`
- Output: PNG plots per replicate

**`plot_MedianDissimilarity_figures_perbinningpipeline.R`**  
Generates median nearest-neighbour dissimilarity line plots per binning pipeline across completeness thresholds. Can be run standalone or sourced from `DistancebyLength.R`.

- Input: `NN_[Software]_summary.tsv` files from `DistancebyLength/[Software]/`
- Output: `fig_main_6pairs_perbinningpipeline.png`, `fig_qc_abc_6pairs_perbinningpipeline.png`, `fig_all20_pairs_perbinningpipeline.png`, `fig_no_abc_14pairs_perbinningpipeline.png`

**`plot_MedianDissimilarity_figures_perreplicate.R`**  
Generates median nearest-neighbour dissimilarity line plots per Hi-C library replicate across completeness thresholds. Can be run standalone or sourced from `DistancebyLength.R`.

- Input: `NN_[replicate]_summary.tsv` files from `DistancebyLength/[replicate]/`
- Output: `fig_main_3panels_perhic.png`, `fig_qc_abc_3panels_perhic.png`, `fig_all5panels_perhic.png`

---

## How to run

### Step 1 — Prepare data

Set paths in `prepare_data_for_similarity_analysis.R` and run it to generate `data.tsv` and `metadata.tsv` per pipeline. If starting from the Zenodo deposit, this step is already done — download the processed files directly.

Optionally, run `summarize_hic_metadata.R` using the `data.tsv` and `metadata.tsv` files as input to produce comprehensive summary tables at both the bin level and the contig level, including per-bin contig statistics, cross-pipeline comparisons, and contig length stratification by completeness.

### Step 2 — Calculate dissimilarities

Set `dir_Similarity_Routput` and `dir_out` in `DistancebyLength.R`. Run twice:

```r
id_tag <- "replicate"   # cross-pipeline comparisons
source("2_distance_calculations/DistancebyLength.R")

id_tag <- "Software"    # cross-replicate comparisons
source("2_distance_calculations/DistancebyLength.R")
```

### Step 3 — Generate figures and tables

Run scripts in `3_figures_tables/` after setting path variables at the top of each script. `make_summary_table_extended.R` must be available before running `metadata_clustering_completeness_contamination_manuscript.R`.

---

## Dependencies

All scripts require R (≥ 4.0) and the following packages:

```r
install.packages(c("tidyverse", "stringi", "ggplot2", "cowplot",
                   "ggtext", "scales", "knitr", "kableExtra"))
```

---

## License

GPL-2.0 — see `LICENSE` file.
