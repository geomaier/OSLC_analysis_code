#!/usr/bin/env Rscript
# scripts/80_fix_metadata_names.R

suppressPackageStartupMessages(library(data.table))

source("config/config.R")

stopifnot(file.exists(METADATA_CSV))
meta <- fread(METADATA_CSV)

fix_name <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub(" - ", "-", x)
  x <- gsub(" +", " ", x)
  x <- gsub(" ", "_", x)
  x
}

meta[[METADATA_SAMPLE_COL]] <- fix_name(meta[[METADATA_SAMPLE_COL]])

# Optional one-off corrections
if (exists("METADATA_PATCHES") && nrow(METADATA_PATCHES) > 0) {
  for (i in seq_len(nrow(METADATA_PATCHES))) {
    meta[[METADATA_SAMPLE_COL]][meta[[METADATA_SAMPLE_COL]] == METADATA_PATCHES$from[i]] <- METADATA_PATCHES$to[i]
  }
}

out_csv <- file.path("outputs", "metadata_fixed.csv")
fwrite(meta, out_csv)
cat("Saved: ", out_csv, "\n", sep = "")
