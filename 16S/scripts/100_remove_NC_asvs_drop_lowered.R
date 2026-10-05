#!/usr/bin/env Rscript

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
pc_samples <- samps[grepl("PC", samps, fixed = TRUE)]
if (length(nc_samples) == 0L) stop("No negative controls matched the configured sample-name patterns.")

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
phy <- prune_taxa(taxa_sums(phy) > 0, phy)

# Drop samples with <= MIN_READS_KEEP reads
reads <- sample_sums(phy)
phy <- prune_samples(reads > MIN_READS_KEEP, phy)
phy <- prune_taxa(taxa_sums(phy) > 0, phy)

controls <- c(nc_samples, pc_samples)
phy <- prune_samples(!(sample_names(phy) %in% controls), phy)
phy <- prune_taxa(taxa_sums(phy) > 0, phy)

for (branch in c("keep7699", "without7699")) {
  samples_to_remove <- intersect("B3-H_ADJ_CAEC_6830", sample_names(phy))
  if (branch == "without7699") {
    samples_to_remove <- c(
      samples_to_remove,
      sample_names(phy)[grepl("7699", sample_names(phy), fixed = TRUE)]
    )
  }
  samples_to_remove <- unique(samples_to_remove)
  phy_branch <- prune_samples(!(sample_names(phy) %in% samples_to_remove), phy)
  phy_branch <- prune_taxa(taxa_sums(phy_branch) > 0, phy_branch)

  out_f <- file.path(
    "outputs",
    paste0("phyloseq_March26_prevalence_unfiltered_biological_", branch, ".rds")
  )
  saveRDS(phy_branch, out_f)
  write.csv(
    data.frame(SampleID = samples_to_remove),
    file.path("outputs", paste0("removed_specific_samples_", branch, ".csv")),
    row.names = FALSE
  )

  if (branch == "without7699") {
    saveRDS(phy_branch, file.path("outputs", "phyloseq_noNC_dropLowRead.rds"))
  }

  cat("Saved: ", out_f, "\n", sep = "")
  cat("Samples: ", nsamples(phy_branch), " | ASVs: ", ntaxa(phy_branch), "\n", sep = "")
}
