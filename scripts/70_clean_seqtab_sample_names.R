#!/usr/bin/env Rscript
# scripts/70_clean_seqtab_sample_names.R

source("config/config.R")

seqtab_in  <- file.path("outputs", "dada2", "seqtab_PacBio_nochim.rds")
seqtab_out <- file.path("outputs", "dada2", "seqtab_PacBio_nochim_cleaned_names.rds")
stopifnot(file.exists(seqtab_in))

seqtab <- readRDS(seqtab_in)
rownames(seqtab) <- gsub("_trimmed_unmapped$", "", rownames(seqtab))

saveRDS(seqtab, seqtab_out)
cat("Saved: ", seqtab_out, "\n", sep = "")
