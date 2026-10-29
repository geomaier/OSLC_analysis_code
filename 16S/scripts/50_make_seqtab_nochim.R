#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
})

source("config/config.R")

dada_dir <- file.path(PATH_FILT, "dada_per_sample")
stopifnot(dir.exists(dada_dir))

dada_files <- sort(list.files(dada_dir, pattern = "_dada\\.rds$", full.names = TRUE))
stopifnot(length(dada_files) > 0)

sample_names <- sub("_dada\\.rds$", "", basename(dada_files))
names(dada_files) <- sample_names

cat("Loading ", length(dada_files), " per-sample dada objects...\n", sep = "")

dds_list <- lapply(dada_files, readRDS)

cat("Building sequence table...\n")
seqtab <- makeSequenceTable(dds_list)

dir.create(file.path("outputs", "dada2"), recursive = TRUE, showWarnings = FALSE)

seqtab_raw_out <- file.path("outputs", "dada2", "seqtab_PacBio_raw.rds")
saveRDS(seqtab, seqtab_raw_out)
cat("Saved raw seqtab: ", seqtab_raw_out, "\n", sep = "")

cat("Removing chimeras (", CHIMERA_method, ")...\n", sep = "")
seqtab_nochim <- removeBimeraDenovo(
  seqtab,
  method      = CHIMERA_method,
  multithread = CHIMERA_multithread
)

seqtab_nochim_out <- file.path("outputs", "dada2", "seqtab_PacBio_nochim.rds")
saveRDS(seqtab_nochim, seqtab_nochim_out)

cat("Saved no-chim seqtab: ", seqtab_nochim_out, "\n", sep = "")
cat("Raw dims: ", nrow(seqtab), " samples x ", ncol(seqtab), " ASVs\n", sep = "")
cat("Nochim dims: ", nrow(seqtab_nochim), " samples x ", ncol(seqtab_nochim), " ASVs\n", sep = "")
