#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
})

source("config/config.R")

# Inputs
derep_path <- file.path(PATH_FILT, "derep")
stopifnot(dir.exists(derep_path))

err_path <- file.path("outputs", "qc", "error_model_PacBio.rds")
stopifnot(file.exists(err_path))
err <- readRDS(err_path)

# Output directory for per-sample dada objects
dada_dir <- file.path(PATH_FILT, "dada_per_sample")
dir.create(dada_dir, showWarnings = FALSE)

derep_files <- sort(list.files(derep_path, pattern = "_derep\\.rds$", full.names = TRUE))
stopifnot(length(derep_files) > 0)

sample_names <- sub("_derep\\.rds$", "", basename(derep_files))

cat("Found ", length(derep_files), " derep files.\n", sep = "")

set.seed(1)

for (i in seq_along(derep_files)) {
  sname    <- sample_names[i]
  out_file <- file.path(dada_dir, paste0(sname, "_dada.rds"))

  if (file.exists(out_file)) {
    cat("Skipping ", i, "/", length(derep_files), ": ", sname, " (already exists)\n", sep = "")
    next
  }

  cat("Running dada ", i, "/", length(derep_files), ": ", sname, "\n", sep = "")

  drp <- readRDS(derep_files[i])

  dds <- dada(
    drp,
    err         = err,
    BAND_SIZE   = PACBIO_BAND_SIZE,
    multithread = DADA_multithread
  )

  saveRDS(dds, out_file)

  rm(drp, dds)
  gc()
}

cat("DONE. Per-sample dada objects in: ", dada_dir, "\n", sep = "")
