#!/usr/bin/env Rscript

source("config/config.R")

seqtab_in  <- file.path("outputs", "dada2", "seqtab_PacBio_nochim.rds")
seqtab_out <- file.path("outputs", "dada2", "seqtab_PacBio_nochim_cleaned_names.rds")
stopifnot(file.exists(seqtab_in))

seqtab <- readRDS(seqtab_in)
samples <- gsub("_trimmed_unmapped", "", rownames(seqtab))
samples <- gsub("_filt$", "", samples)

name_corrections <- c(
  "B0-SSL_SIG_7241" = "B0-HP_SIG_7241",
  "B1-SSL_DESC_7154" = "B1-HP_DESC_7154",
  "B1-SSL_DESC_7663" = "B1-HP_DESC_7663",
  "B1-SSL_RECT_8977" = "B1-HP_RECT_8977",
  "B2-SSL_DESC_1554" = "B2-HP_DESC_1554",
  "B2-SSL_SIG_6830" = "B2-HP_SIG_6830",
  "B2-SSL_SIG_7154" = "B2-HP_SIG_7154",
  "B2-SSL_SIG_7784" = "B2-HP_SIG_7784",
  "B3-HP_DESC_2102" = "B3-SSL_DESC_2102",
  "B3-HP_TRVS_7241" = "B3-SSL_TRVS_7241"
)

idx <- match(samples, names(name_corrections))
samples[!is.na(idx)] <- unname(name_corrections[idx[!is.na(idx)]])
stopifnot(!anyDuplicated(samples))
rownames(seqtab) <- samples

saveRDS(seqtab, seqtab_out)
cat("Saved: ", seqtab_out, "\n", sep = "")
