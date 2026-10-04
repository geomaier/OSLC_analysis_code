#!/usr/bin/env Rscript
# scripts/105_rarefy_deep_samples_95pct.R

suppressPackageStartupMessages({
  library(phyloseq)
})

source("config/config.R")

for (branch in c("keep7699", "without7699")) {
  in_f <- file.path(
    "outputs",
    paste0("phyloseq_March26_prevalence_unfiltered_biological_", branch, ".rds")
  )
  out_f <- file.path("outputs", paste0("phyloseq_March26_post_rarefaction_", branch, ".rds"))
  stopifnot(file.exists(in_f))
  
  phy <- readRDS(in_f)
  
  reads <- sample_sums(phy)
  depth_95 <- as.numeric(quantile(reads, 0.95))
  deep_samples <- names(reads[reads > depth_95])
  
  cat("Branch: ", branch, "\n", sep = "")
  cat("95th percentile depth: ", depth_95, "\n", sep = "")
  cat("Deep samples to rarefy: ", length(deep_samples), "\n", sep = "")
  
  if (length(deep_samples) > 0) {
    phy_deep    <- prune_samples(deep_samples, phy)
    phy_shallow <- prune_samples(!sample_names(phy) %in% deep_samples, phy)
    
    phy_deep_rar <- rarefy_even_depth(
      phy_deep,
      sample.size = floor(depth_95),
      rngseed = RAREFY_seed,
      replace = FALSE,
      verbose = TRUE
    )
    
    phy <- merge_phyloseq(phy_shallow, phy_deep_rar)
  }
  phy <- prune_taxa(taxa_sums(phy) > 0, phy)
  
  saveRDS(phy, out_f)
  reads_after <- sample_sums(phy)
  read_table <- data.frame(SampleID = names(reads_after), Reads = as.numeric(reads_after))
  read_table <- read_table[order(read_table$Reads, decreasing = TRUE), ]
  write.csv(
    read_table,
    file.path("outputs", paste0("sample_read_depths_post_rarefaction_", branch, ".csv")),
    row.names = FALSE
  )
  
  if (branch == "without7699") {
    saveRDS(phy, file.path("outputs", "phyloseq_noNC_dropLowRead_deepRarefied95.rds"))
  }
  
  cat("Saved: ", out_f, "\n", sep = "")
  cat("Samples: ", nsamples(phy), " | ASVs: ", ntaxa(phy), "\n", sep = "")
}
