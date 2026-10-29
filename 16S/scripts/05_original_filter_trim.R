#!/usr/bin/env Rscript

library(dada2)

source("config/config.R")

stopifnot(dir.exists(PATH_RAW_FASTQ_NOV2025))

fn <- sort(list.files(
  PATH_RAW_FASTQ_NOV2025,
  pattern = "\\.fastq\\.gz$",
  full.names = TRUE
))
stopifnot(length(fn) > 0)

sample.names <- sub("\\.fastq\\.gz$", "", basename(fn))
filt <- file.path(
  PATH_FILT_NOV2025,
  paste0(sample.names, "_filt.fastq.gz")
)

dir.create(PATH_FILT_NOV2025, recursive = TRUE, showWarnings = FALSE)

filtered <- filterAndTrim(
  fwd = fn,
  filt = filt,
  truncLen = FILTER_truncLen,
  minQ = FILTER_minQ,
  minLen = FILTER_minLen,
  maxLen = FILTER_maxLen,
  maxN = FILTER_maxN,
  maxEE = FILTER_maxEE,
  rm.phix = FILTER_rm_phix,
  multithread = FILTER_multithread
)

write.csv(filtered, file.path(PATH_FILT_NOV2025, "filtering_summary.csv"))
print(filtered)
