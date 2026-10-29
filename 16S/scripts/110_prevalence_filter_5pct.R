#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(phyloseq)
})

source("config/config.R")

stopifnot(exists("PREV_FRAC"), PREV_FRAC == 0.05)
fractions <- c("5pct" = PREV_FRAC)

for (branch in c("keep7699", "without7699")) {
  phy_f <- file.path("outputs", paste0("phyloseq_March26_post_rarefaction_", branch, ".rds"))
  stopifnot(file.exists(phy_f))
  phy <- readRDS(phy_f)

  otu_mat <- as(otu_table(phy), "matrix")
  if (taxa_are_rows(phy)) otu_mat <- t(otu_mat)

  n_samp <- nsamples(phy)
  prev <- colSums(otu_mat > 0)

  for (label in names(fractions)) {
    min_samples <- ceiling(fractions[[label]] * n_samp)
    keep_taxa <- names(prev)[prev >= min_samples]
    phy2 <- prune_taxa(keep_taxa, phy)
    phy2 <- prune_taxa(taxa_sums(phy2) > 0, phy2)

    out_f <- file.path("outputs", paste0("phyloseq_final_March26_", label, "_", branch, ".rds"))
    saveRDS(phy2, out_f)

    if (label == "5pct") {
      final_reads <- sample_sums(phy2)
      final_sample_df <- data.frame(SampleID = names(final_reads), Reads = as.numeric(final_reads))
      final_sample_df <- final_sample_df[order(final_sample_df$SampleID), ]
      write.csv(
        final_sample_df,
        file.path("outputs", paste0("final_samples_5pct_", branch, ".csv")),
        row.names = FALSE
      )
      if (branch == "without7699") {
        saveRDS(phy2, file.path("outputs", "phyloseq_after_reseq_noNC_5_percent.rds"))
      }
    }

    cat("Prevalence threshold: ASVs present in >= ", min_samples, " of ", n_samp,
        " samples (", fractions[[label]] * 100, "%)\n", sep = "")
    cat("Saved: ", out_f, "\n", sep = "")
    cat("Samples: ", nsamples(phy2), " | ASVs: ", ntaxa(phy2), "\n", sep = "")
  }
}
