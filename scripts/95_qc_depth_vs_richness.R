#!/usr/bin/env Rscript
# scripts/95_qc_depth_vs_richness.R

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(data.table)
})

phy_f <- file.path("outputs", "phyloseq_raw.rds")
stopifnot(file.exists(phy_f))
phy <- readRDS(phy_f)

otu_mat <- as(otu_table(phy), "matrix")
if (taxa_are_rows(phy)) otu_mat <- t(otu_mat)

df <- data.frame(
  Sample = rownames(otu_mat),
  ASVs   = rowSums(otu_mat > 0),
  Reads  = rowSums(otu_mat),
  stringsAsFactors = FALSE
)

r <- cor(df$ASVs, df$Reads, method = "pearson")

dir.create(file.path("outputs", "qc"), recursive = TRUE, showWarnings = FALSE)
fwrite(df, file.path("outputs", "qc", "depth_vs_asv_table.csv"))

p <- ggplot(df, aes(x = Reads, y = ASVs)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Read depth vs ASV richness",
    subtitle = paste0("Pearson r = ", round(r, 3)),
    x = "Total reads per sample",
    y = "ASVs detected"
  ) +
  theme_minimal(base_size = 14)

ggsave(file.path("outputs", "qc", "depth_vs_asv_richness.png"), p, width = 7, height = 5, dpi = 300)

cat("Saved outputs/qc/depth_vs_asv_richness.png\n")
