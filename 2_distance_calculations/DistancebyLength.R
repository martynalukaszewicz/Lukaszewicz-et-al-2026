########################################################################
## DistancebyLength.R
##
## Calculates pairwise Jaccard dissimilarity values between
## metagenome-assembled genome (MAG) bins across:
##
##   (1) Cross-pipeline comparisons (same Hi-C library, different pipeline):
##       MetaTOR vs ProxiMeta, MetaTOR vs Bin3C, ProxiMeta vs Bin3C
##       — for each of: E_15_1, E_15_2, E_15_3A, E_15_3B, E_15_3C
##       Set: id_tag <- "replicate"
##
##   (2) Cross-replicate comparisons (same pipeline, different Hi-C library):
##       E_15_3A vs E_15_3B, E_15_3A vs E_15_3C, E_15_3B vs E_15_3C, etc.
##       — for each of: ProxiMeta, Bin3C, MetaTOR
##       Set: id_tag <- "Software"
##
## REQUIRES:
##   The four helper scripts below must be in dir_out:
##     getdata2_DistancebyLength.R
##     savesimMat_DistancebyLength.R
##     savesummary_DistancebyLength.R
##     saveNNsummary_DistancebyLength.R
##
##   The processed data files must exist in dir_Similarity_Routput:
##     proximeta/data.tsv, proximeta/metadata.tsv
##     bin3c/data.tsv,     bin3c/metadata.tsv
##     metator/data.tsv,   metator/metadata.tsv
##
## OUTPUTS (written to dir_out subfolders):
##   [replicate or Software]/
##     *.csv                        — pairwise Jaccard simMat matrices
##     [folder]_summary.tsv         — NN dissimilarity summary per group
##     NN_[folder]_summary.tsv      — same, with NAs removed
##
## Authors: Benjamin J. Ridenhour, Martyna Lukaszewicz, Eva M. Top,  Thibault Stalder
## University of Idaho
########################################################################


# -----------------------------------------------------------------------
# USER PARAMETERS — set these paths before running
# -----------------------------------------------------------------------

# Path to processed data TSV files (output of prepare_data_for_similarity_analysis.R)
dir_Similarity_Routput <- "~/path/to/Similarity_Routput/"

# Path to this script's folder (also where helper scripts live and outputs are written)
dir_out <- "~/path/to/DistancebyLength/"

# -----------------------------------------------------------------------
# END USER PARAMETERS
# -----------------------------------------------------------------------


rm(list = ls())

library(tidyverse)
library(stringr)
library(ggplot2)
library(dplyr)
library(readr)
library(cowplot)

## Pin dplyr verbs to avoid masking conflicts
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename


# -----------------------------------------------------------------------
# Source helper functions
# -----------------------------------------------------------------------

source(file.path(dir_out, "getdata2_DistancebyLength.R"))
source(file.path(dir_out, "savesimMat_DistancebyLength.R"))
source(file.path(dir_out, "savesummary_DistancebyLength.R"))
source(file.path(dir_out, "saveNNsummary_DistancebyLength.R"))


# -----------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------

contamination_cutoff <- 10
genomesize_cutoff    <- c("100Kb" = 100e3, "10Mb" = 10e6)
complete             <- c(0, 85, 90, 95)   ## completeness thresholds for simMat calculations

## Set id_tag to control which comparison is run:
##   "replicate" — cross-pipeline comparisons per Hi-C library replicate
##   "Software"  — cross-replicate comparisons per pipeline
##
## Run this script TWICE: once per id_tag value.

id_tag <- "replicate"   ## change to "Software" for cross-replicate analysis
## id_tag <- "Software"


# -----------------------------------------------------------------------
# Load and prepare data
# -----------------------------------------------------------------------

data2        <- getdata2(contamination_cutoff, genomesize_cutoff, dir_Similarity_Routput)
data2_tibble <- as_tibble(data2)


# -----------------------------------------------------------------------
# Run calculations
# -----------------------------------------------------------------------

## Step 1: compute and save pairwise Jaccard dissimilarity matrices
## savesimMat skips any CSV that already exists on disk -- safe to re-run
## after interruption or after adding new thresholds to complete above.
## Delete existing CSVs to force recompute.
savesimMat(data2, id_tag, dir_out, complete)

## Step 2: summarise matrices into NN dissimilarity tables
## savesummary processes ALL CSV files found in each subfolder, not just
## those matching the current complete vector -- so CSVs from previous
## runs with different thresholds are automatically included.
savesummary(id_tag, dir_out)

## Step 3: remove NAs from summary tables
saveNNsummary(id_tag, dir_out)

message("DistancebyLength.R complete for id_tag = '", id_tag, "'.")
message("Output written to: ", dir_out)
