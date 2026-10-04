#!/usr/bin/env Rscript
# scripts/10_filter_trim.R

suppressPackageStartupMessages({
  library(dada2)
  library(data.table)
})

# Load local config (do NOT commit config/config.R if it has private paths)
source("config/config.R")

needed <- c(
  "PATH_RAW_FASTQ","PATH_FILT",
  "FILTER_truncLen","FILTER_minQ","FILTER_minLen","FILTER_maxLen",
  "FILTER_maxN","FILTER_maxEE","FILTER_rm_phix","FILTER_multithread"
)
missing <- needed[!vapply(needed, exists, logical(1))]
if (length(missing) > 0) stop("Missing in config/config.R: ", paste(missing, collapse = ", "))

dir.create(PATH_FILT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("outputs", "qc"), recursive = TRUE, showWarnings = FALSE)

# List all FASTQ files
fn <- sort(list.files(PATH_RAW_FASTQ, pattern = "\\.fastq\\.gz$", full.names = TRUE))
stopifnot(length(fn) > 0)

# Optional warning for expected PC file
pc_expected <- "BATCH_9_PC_trimmed_unmapped.fastq.gz"
if (!any(basename(fn) == pc_expected)) {
  warning("Expected PC file not found in PATH_RAW_FASTQ: ", pc_expected)
}

# Sample names
sample.names <- sub("\\.fastq\\.gz$", "", basename(fn))
filt <- file.path(PATH_FILT, paste0(sample.names, "_filt.fastq.gz"))

# Filter & trim
filtered <- filterAndTrim(
  fwd  = fn,
  filt = filt,
  truncLen = FILTER_truncLen,
  minQ     = FILTER_minQ,
  minLen   = FILTER_minLen,
  maxLen   = FILTER_maxLen,
  maxN     = FILTER_maxN,
  maxEE    = FILTER_maxEE,
  rm.phix  = FILTER_rm_phix,
  multithread = FILTER_multithread
)

# Save summary table
filter_summary <- data.frame(
  sample    = rownames(filtered),
  reads_in  = filtered[, "reads.in"],
  reads_out = filtered[, "reads.out"],
  pass_rate = round(100 * filtered[, "reads.out"] / filtered[, "reads.in"], 1),
  stringsAsFactors = FALSE
)

out_csv <- file.path("outputs", "qc", "filterAndTrim_summary.csv")
fwrite(filter_summary, out_csv)

cat("Saved filtering summary to: ", out_csv, "\n")
cat("Filtered FASTQs written to: ", PATH_FILT, "\n")
