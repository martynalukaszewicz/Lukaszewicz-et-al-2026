## prepare_data_for_similarity_analysis.R
##
## This script creates and saves data.tsv and metadata.tsv files per pipeline
## in (output_dir_main / software) directories, for use as input to
## DistanceByLength.R
##
## Pipelines handled: ProxiMeta, Bin3C, MetaTOR
## Replicates handled: E_15_1, E_15_2, E_15_3A, E_15_3B, E_15_3C
##
## Output files per pipeline:
##   data.tsv     -- contig-to-bin assignments with lengths
##   metadata.tsv -- per-bin CheckM quality metrics and Mash taxonomy
##
## Authors: Martyna Lukaszewicz, Benjamin J. Ridenhour, Eva M. Top,  Thibault Stalder
## University of Idaho
## -----------------------------------------------------------------------


# Clear environment first -- must be before everything else
rm(list = ls())

library(tidyverse)
library(stringi)

# Pin dplyr verbs explicitly to avoid masking conflicts with MASS or
# other packages that export select() or filter()
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename


# -----------------------------------------------------------------------
# USER PARAMETERS -- set these paths before running
# -----------------------------------------------------------------------

# Directory containing contig_ids_per_clusters subfolders per pipeline
dir_ids <- "~/path/to/contig_ids_per_clusters/"

# Subfolder names inside dir_ids, one per pipeline
software_folders <- c("proximeta/preprocess/combined", "bin3c", "metator")

# Directory containing CheckM report subfolders per pipeline
report_dir_main <- "~/path/to/report/"

# Output directory: data.tsv and metadata.tsv will be written here
# under subdirectories named by pipeline (proximeta/, bin3c/, metator/)
output_dir_main <- "~/path/to/Similarity_Routput/"

# -----------------------------------------------------------------------
# END USER PARAMETERS
# -----------------------------------------------------------------------


# -----------------------------------------------------------------------
# load_data()
#
# Reads raw contig-cluster assignment files for one pipeline, reformats
# them into a tidy data frame (cluster, contig, replicate, length),
# and writes data.tsv to the output directory.
#
# Arguments:
#   software_folder -- one element of software_folders (character)
#
# Returns:
#   tidy data frame for that pipeline
# -----------------------------------------------------------------------

load_data <- function(software_folder) {

  dir_path <- file.path(dir_ids, software_folder)

  # Inner loader: read one file and tag it with its filename
  load1 <- function(x, y) {
    read_delim(x, col_names = FALSE, delim = "\t", show_col_types = FALSE) %>%
      mutate(Sample = y)
  }

  # Read all .txt files in the pipeline directory
  temp <- list.files(path = dir_path, pattern = "*.txt", full.names = TRUE)
  if (length(temp) == 0) {
    stop(paste0("No .txt files found in: ", dir_path))
  }

  data <- data.frame()
  for (i in seq_along(temp)) {
    df1 <- load1(temp[i], temp[i])
    data <- rbind(data, df1)
  }
  data$Sample <- basename(data$Sample)


  # -------------------------------------------------------------------
  # Pipeline-specific parsing
  # -------------------------------------------------------------------

  if (software_folder == "proximeta/preprocess/combined") {

    message("Parsing ProxiMeta contig assignments...")

    data <- data %>%
      separate(X2, into = c("contig", "length"), sep = " ", remove = TRUE) %>%
      dplyr::select(X1, contig, length, Sample) %>%
      separate(X1, into = c("cluster", "junk"), sep = ".f", remove = TRUE) %>%
      dplyr::select(cluster, junk, contig, length, Sample) %>%
      separate(junk, into = c("junk", "contig2"), sep = ":>", remove = TRUE) %>%
      dplyr::select(cluster, contig2, contig, length, Sample) %>%
      separate(contig, into = c("junk2", "contig"), sep = ">", remove = TRUE) %>%
      dplyr::select(cluster, contig2, contig, length, Sample)

    if (identical(data$contig2, data$contig)) {
      message("ProxiMeta: contig2 and contig columns are identical -- dropping contig2.")
      data <- dplyr::select(data, -contig2)
    } else {
      warning("ProxiMeta: contig2 and contig columns are NOT identical -- check input files.")
    }

    data <- data %>%
      separate(Sample, into = c("replicate", "junk3"), sep = "_co", remove = TRUE) %>%
      dplyr::select(cluster, contig, replicate, length)

    # Update label for output directory
    software_folder <- "proximeta"

  } else if (software_folder == "bin3c") {

    message("Parsing Bin3C contig assignments...")

    data <- data %>%
      dplyr::rename(cluster = X1) %>%
      separate(cluster,
               into = c("cluster", "contig", "junk", "junk1", "length"),
               sep = " ", remove = TRUE) %>%
      dplyr::select(cluster, contig, length, Sample) %>%
      separate(contig, into = c("junk2", "contig"), sep = "g:", remove = TRUE) %>%
      dplyr::select(cluster, contig, length, Sample) %>%
      separate(Sample, into = c("junk3", "replicate"), sep = "3c_", remove = TRUE) %>%
      dplyr::select(cluster, contig, length, replicate) %>%
      separate(replicate, into = c("replicate", "junk4"), sep = "_clu", remove = TRUE) %>%
      dplyr::select(cluster, contig, length, replicate) %>%
      separate(cluster, into = c("junk", "cluster"), sep = ">", remove = TRUE) %>%
      dplyr::select(cluster, contig, length, replicate) %>%
      separate(cluster, into = c("cluster", "junk"), sep = "_", remove = TRUE) %>%
      dplyr::select(cluster, contig, replicate, length) %>%
      relocate(length, .after = replicate)

  } else if (software_folder == "metator") {

    message("Parsing MetaTOR contig assignments...")

    data <- data %>%
      dplyr::rename(cluster = X1) %>%
      separate(cluster, into = c("cluster", "contig"),
               sep = ".fa:>", remove = TRUE) %>%
      separate(contig,
               into = c("contig", "junk1", "junk2", "length"),
               sep = " ", remove = TRUE) %>%
      dplyr::select(cluster, contig, length, Sample) %>%
      separate(length, into = c("junk", "length"), sep = "len=", remove = TRUE) %>%
      dplyr::select(cluster, contig, length, Sample) %>%
      separate(Sample, into = c("junk", "replicate"), sep = "r_", remove = TRUE) %>%
      dplyr::select(cluster, contig, length, replicate) %>%
      separate(replicate, into = c("replicate", "junk"), sep = "_b", remove = TRUE) %>%
      dplyr::select(cluster, contig, length, replicate) %>%
      relocate(length, .after = replicate)

  } else {
    stop(paste0("Unrecognised software_folder: ", software_folder,
                "\nExpected one of: 'proximeta/preprocess/combined', 'bin3c', 'metator'"))
  }

  message("Replicates found: ", paste(unique(data$replicate), collapse = ", "))

  # Write data.tsv
  output_dir <- file.path(output_dir_main, software_folder)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  write_tsv(data, file.path(output_dir, "data.tsv"))
  message("Written: ", file.path(output_dir, "data.tsv"))

  return(data)
}


# -----------------------------------------------------------------------
# run_analysis()
#
# Orchestrates load_data() for all pipelines, then reads CheckM report
# files for each pipeline and writes metadata.tsv per pipeline.
# -----------------------------------------------------------------------

run_analysis <- function() {

  # --- Step 1: parse and write data.tsv for all pipelines ---
  results <- lapply(software_folders, load_data)

  # --- Step 2: parse and write metadata.tsv for all pipelines ---

  # Extract top-level folder names (proximeta, bin3c, metator) for report dirs
  report_software_folders <- sapply(
    stri_split_fixed(software_folders, pattern = "/", n = -1, simplify = FALSE),
    `[`, 1
  )
  report_dirs <- file.path(report_dir_main, report_software_folders)

  # Standardised replicate name mapping
  standardise_replicate <- function(x) {
    dplyr::case_when(
      grepl("E_15_1",  x) ~ "E_15_1",
      grepl("E_15_2",  x) ~ "E_15_2",
      grepl("E_15_3A", x) ~ "E_15_3A",
      grepl("E_15_3B", x) ~ "E_15_3B",
      grepl("E_15_3C", x) ~ "E_15_3C",
      TRUE ~ x
    )
  }

  for (report_dir in report_dirs) {

    pipeline <- basename(report_dir)
    message("\nProcessing metadata for: ", pipeline)

    # Choose the right reader for each pipeline's report format
    if (pipeline == "proximeta") {
      load3 <- function(x, y) read_tsv(x,  col_names = TRUE, show_col_types = FALSE) %>% mutate(Sample = y)
    } else if (pipeline == "bin3c") {
      load3 <- function(x, y) read_csv(x,  col_names = TRUE, show_col_types = FALSE) %>% mutate(Sample = y)
    } else if (pipeline == "metator") {
      load3 <- function(x, y) read_tsv(x,  col_names = TRUE, show_col_types = FALSE) %>% mutate(Sample = y)
    } else {
      stop(paste0("Unrecognised pipeline in report_dirs: ", pipeline))
    }

    # List per-replicate report files
    temp            <- list.files(path = report_dir, pattern = "*E", full.names = TRUE)
    replicate_names <- standardise_replicate(basename(temp))

    if (length(temp) == 0) {
      warning(paste0("No report files found in: ", report_dir, " -- skipping."))
      next
    }

    metadata <- data.frame()

    if (pipeline == "proximeta") {

      for (i in seq_along(temp)) {
        df1 <- load3(temp[i], temp[i])
        # Keep columns up to and including kmer_coherence
        df2 <- df1[, 1:grep("kmer_coherence", colnames(df1))]
        df2$software  <- "proximeta"
        df2$replicate <- replicate_names[i]
        metadata <- rbind(metadata, df2)
      }

    } else if (pipeline == "bin3c") {

      for (i in seq_along(temp)) {
        df1 <- load3(temp[i], temp[i]) %>%
          dplyr::rename(cluster_id = name) %>%
          dplyr::select(-id)
        df1$software  <- "bin3c"
        df1$replicate <- replicate_names[i]
        df1$Sample    <- basename(df1$Sample)
        metadata <- rbind(metadata, df1)
      }

    } else if (pipeline == "metator") {

      for (i in seq_along(temp)) {
        df1 <- load3(temp[i], temp[i])
        # MetaTOR report files have an unnamed first column which read_tsv
        # may load as "X1", "...1", or another auto-name depending on the
        # file and tidyverse version. Rename whatever the first column is
        # called to cluster_id, regardless of its current name.
        colnames(df1)[1] <- "cluster_id"
        df1$software  <- "metator"
        df1$replicate <- replicate_names[i]
        df1$Sample    <- basename(df1$Sample)
        metadata <- rbind(metadata, df1)
      }

    }

    # Write metadata.tsv
    output_dir <- file.path(output_dir_main, pipeline)
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    write_tsv(metadata, file.path(output_dir, "metadata.tsv"))
    message("Written: ", file.path(output_dir, "metadata.tsv"))
  }

  message("\nAll pipelines complete.")
  return(results)
}


# -----------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------

results <- run_analysis()
