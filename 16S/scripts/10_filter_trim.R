#!/usr/bin/env Rscript

library(dada2)

source("config/config.R")

needed <- c(
  "PATH_RAW_FASTQ", "PATH_FILT",
  "FILTER_truncLen", "FILTER_minQ", "FILTER_minLen", "FILTER_maxLen",
  "FILTER_maxN", "FILTER_maxEE", "FILTER_rm_phix", "FILTER_multithread"
)

missing <- needed[!vapply(needed, exists, logical(1))]
if (length(missing) > 0) {
  stop("Missing configuration: ", paste(missing, collapse = ", "))
}

stopifnot(dir.exists(PATH_RAW_FASTQ))

fn <- sort(list.files(
  PATH_RAW_FASTQ,
  pattern = "\\.fastq\\.gz$",
  full.names = TRUE
))

stopifnot(length(fn) > 0)

sample.names <- sub("\\.fastq\\.gz$", "", basename(fn))
filt <- file.path(PATH_FILT, paste0(sample.names, "_filt.fastq.gz"))

dir.create(PATH_FILT, recursive = TRUE, showWarnings = FALSE)

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

write.csv(
  filtered,
  file.path(PATH_FILT, "filtering_summary.csv")
)