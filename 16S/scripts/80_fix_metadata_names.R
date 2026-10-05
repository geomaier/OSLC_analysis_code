#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

source("config/config.R")

stopifnot(file.exists(METADATA_CSV))
meta <- read.csv(
  METADATA_CSV,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

meta <- meta[, names(meta) != "", drop = FALSE]
stopifnot(METADATA_SAMPLE_COL == "Sample_clean")
sample_ids <- rownames(meta)
stopifnot(!anyNA(sample_ids), all(nzchar(sample_ids)), !anyDuplicated(sample_ids))
stopifnot(METADATA_SAMPLE_COL %in% names(meta))
existing_ids <- as.character(meta[[METADATA_SAMPLE_COL]])
has_id <- !is.na(existing_ids) & nzchar(trimws(existing_ids))
stopifnot(all(existing_ids[has_id] == sample_ids[has_id]))
meta[[METADATA_SAMPLE_COL]] <- sample_ids

out_csv <- file.path("outputs", "metadata_fixed.csv")
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
fwrite(meta, out_csv)
cat("Saved: ", out_csv, "\n", sep = "")
