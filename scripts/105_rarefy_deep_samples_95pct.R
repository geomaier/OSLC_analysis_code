#!/usr/bin/env Rscript
# scripts/105_rarefy_deep_samples_95pct.R

suppressPackageStartupMessages({
  library(phyloseq)
})

source("config/config.R")

in_f  <- file.path("outputs", "phyloseq_noNC_dropLowRead.rds")
out_f <- file.path("outputs", "phyloseq_noNC_dropLowRead_deepRarefied95.rds")
stopifnot(file.exists(in_f))

phy <- readRDS(in_f)

reads <- sample_sums(phy)
depth_95 <- as.numeric(quantile(reads, 0.95))

deep_samples <- names(reads[reads > depth_95])

cat("95th percentile depth: ", depth_95, "\n", sep = "")
cat("Deep samples to rarefy: ", length(deep_samples), "\n", sep = "")

if (length(deep_samples) > 0) {
  phy_deep    <- prune_samples(deep_samples, phy)
  phy_shallow <- prune_samples(!sample_names(phy) %in% deep_samples, phy)
  
  phy_deep_rar <- rarefy_even_depth(
    phy_deep,
    sample.size = depth_95,
    rngseed = RAREFY_seed,
    replace = FALSE,
    verbose = TRUE
  )
  
  phy <- merge_phyloseq(phy_shallow, phy_deep_rar)
}

saveRDS(phy, out_f)

cat("Saved: ", out_f, "\n", sep = "")
cat("Samples: ", nsamples(phy), " | ASVs: ", ntaxa(phy), "\n", sep = "")
