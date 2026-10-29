#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
})

source("config/config.R")

derep_path <- file.path(PATH_FILT, "derep")
stopifnot(dir.exists(derep_path))

stopifnot(exists("ERROR_DEREP_FILES"))
stopifnot(length(ERROR_DEREP_FILES) >= 1)

derep_full <- file.path(derep_path, ERROR_DEREP_FILES)
missing <- ERROR_DEREP_FILES[!file.exists(derep_full)]
if (length(missing) > 0) stop("Missing derep files: ", paste(missing, collapse = ", "))

drp_list <- lapply(derep_full, readRDS)

set.seed(1)
err <- learnErrors(
  drp_list,
  errorEstimationFunction = PacBioErrfun,
  BAND_SIZE = PACBIO_BAND_SIZE,
  multithread = LEARNERR_multithread
)

dir.create(file.path("outputs", "qc"), recursive = TRUE, showWarnings = FALSE)

saveRDS(err, file.path("outputs", "qc", "error_model_PacBio.rds"))

png(file.path("outputs", "qc", "error_plot.png"), width = 1200, height = 800)
print(plotErrors(err, nominalQ = TRUE))
dev.off()

cat("Saved outputs/qc/error_model_PacBio.rds\n")
cat("Saved outputs/qc/error_plot.png\n")
