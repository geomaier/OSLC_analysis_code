#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
  library(data.table)
})

source("config/config.R")

seqtab_nochim_in <- file.path("outputs", "dada2", "seqtab_PacBio_nochim.rds")
stopifnot(file.exists(seqtab_nochim_in))

stopifnot(exists("SILVA_TRAIN_FA"))
stopifnot(exists("SILVA_SPEC_FA"))

seqtab_nochim <- readRDS(seqtab_nochim_in)

cat("Seqtab: ", nrow(seqtab_nochim), " samples x ", ncol(seqtab_nochim), " ASVs\n", sep = "")

cat("Assigning taxonomy...\n")
taxa <- assignTaxonomy(
  seqtab_nochim,
  SILVA_TRAIN_FA,
  multithread = TAX_multithread,
  tryRC       = TRUE
)

cat("Adding species...\n")
taxa <- addSpecies(taxa, SILVA_SPEC_FA)

dir.create(file.path("outputs", "dada2"), recursive = TRUE, showWarnings = FALSE)

tax_rds_out <- file.path("outputs", "dada2", "taxonomy_PacBio_SILVA_v138.1.rds")
saveRDS(taxa, tax_rds_out)

tax_df <- as.data.frame(taxa)
tax_df$ASV <- rownames(tax_df)
tax_df <- tax_df[, c("ASV", setdiff(colnames(tax_df), "ASV"))]

tax_csv_out <- file.path("outputs", "dada2", "taxonomy_PacBio_SILVA_v138.1.csv")
fwrite(tax_df, tax_csv_out)

cat("Saved taxonomy RDS: ", tax_rds_out, "\n", sep = "")
cat("Saved taxonomy CSV: ", tax_csv_out, "\n", sep = "")
