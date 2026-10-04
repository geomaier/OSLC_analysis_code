#!/usr/bin/env Rscript
# scripts/110_prevalence_filter_5pct.R

suppressPackageStartupMessages({
  library(phyloseq)
})

source("config/config.R")

# Input: output of scripts/105_rarefy_deep_samples_95pct.R
phy_f <- file.path("outputs", "phyloseq_noNC_dropLowRead_deepRarefied95.rds")
stopifnot(file.exists(phy_f))
phy <- readRDS(phy_f)

# Prevalence filter fraction (e.g. 0.05)
stopifnot(exists("PREV_FRAC"))

otu_mat <- as(otu_table(phy), "matrix")
if (taxa_are_rows(phy)) otu_mat <- t(otu_mat)

n_samp <- nsamples(phy)
min_samples <- ceiling(PREV_FRAC * n_samp)

prev <- colSums(otu_mat > 0)
keep_taxa <- names(prev)[prev >= min_samples]

phy2 <- prune_taxa(keep_taxa, phy)

# Final output name (matches your original naming)
out_f <- file.path("outputs", "phyloseq_after_reseq_noNC_5_percent.rds")
saveRDS(phy2, out_f)

cat("Prevalence threshold: ASVs present in ≥ ", min_samples, " of ", n_samp,
    " samples (", PREV_FRAC * 100, "%)\n", sep = "")
cat("Saved: ", out_f, "\n", sep = "")
cat("Samples: ", nsamples(phy2), " | ASVs: ", ntaxa(phy2), "\n", sep = "")
