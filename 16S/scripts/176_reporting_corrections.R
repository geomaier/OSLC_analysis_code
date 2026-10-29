#!/usr/bin/env Rscript
source("config/config.R")
source("R/reporting_corrections.R")
run_reporting_corrections(PATH_WORK,
  deseq_run = if (exists("DESEQ2_REVIEWED_RUN")) DESEQ2_REVIEWED_RUN else NULL,
  matched_run = if (exists("MATCHED_REVIEWED_RUN")) MATCHED_REVIEWED_RUN else NULL)
