#!/usr/bin/env Rscript
# scripts/100_remove_NC_asvs_drop_lowread.R

suppressPackageStartupMessages({
  library(phyloseq)
})

source("config/config.R")

phy_f <- file.path("outputs", "phyloseq_raw.rds")
stopifnot(file.exists(phy_f))
phy <- readRDS(phy_f)

# Identify NC samples by pattern in sample name
samps <- sample_names(phy)
nc_idx <- grepl(NC_PATTERN1, samps, fixed = TRUE) | grepl(NC_PATTERN2, samps, fixed = TRUE)
nc_samples <- samps[nc_idx]

cat("NC samples detected: ", length(nc_samples), "\n", sep = "")
if (length(nc_samples) > 0) print(nc_samples)

otu_mat <- as(otu_table(phy), "matrix")
if (taxa_are_rows(phy)) otu_mat <- t(otu_mat)

# ASVs present in any NC sample
nc_taxa <- character(0)
if (length(nc_samples) > 0) {
  nc_otu <- otu_mat[nc_samples, , drop = FALSE]
  nc_taxa <- colnames(otu_mat)[colSums(nc_otu > 0) > 0]
}

cat("ASVs present in NC samples: ", length(nc_taxa), "\n", sep = "")

# Remove those ASVs globally
if (length(nc_taxa) > 0) {
  phy <- prune_taxa(!(taxa_names(phy) %in% nc_taxa), phy)
}

# Drop samples with <= MIN_READS_KEEP reads
reads <- sample_sums(phy)
phy <- prune_samples(reads > MIN_READS_KEEP, phy)

out_f <- file.path("outputs", "phyloseq_noNC_dropLowRead.rds")
saveRDS(phy, out_f)

cat("Saved: ", out_f, "\n", sep = "")
cat("Samples: ", nsamples(phy), " | ASVs: ", ntaxa(phy), "\n", sep = "")
