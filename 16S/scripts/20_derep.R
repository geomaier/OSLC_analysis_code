#!/usr/bin/env Rscript
# scripts/20_derep.R

suppressPackageStartupMessages({
  library(dada2)
})

source("config/config.R")

# Uses filtered FASTQs produced by scripts/10_filter_trim.R
stopifnot(dir.exists(PATH_FILT))

derep_path <- file.path(PATH_FILT, "derep")
dir.create(derep_path, showWarnings = FALSE)

filt_files <- sort(list.files(PATH_FILT, pattern = "_filt\\.fastq\\.gz$", full.names = TRUE))
stopifnot(length(filt_files) > 0)

sample.names <- sub("_filt\\.fastq\\.gz$", "", basename(filt_files))

cat("Found ", length(filt_files), " filtered FASTQs.\n", sep = "")

for (i in seq_along(filt_files)) {
  cat("Dereplicating ", i, "/", length(filt_files), ": ", sample.names[i], "\n", sep = "")
  drp <- derepFastq(filt_files[i])
  saveRDS(drp, file.path(derep_path, paste0(sample.names[i], "_derep.rds")))
  rm(drp)
  gc()
}

cat("DONE. Dereplicated files saved to: ", derep_path, "\n", sep = "")
