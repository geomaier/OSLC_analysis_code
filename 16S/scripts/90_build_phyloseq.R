#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(phyloseq)
  library(data.table)
})

source("config/config.R")

seqtab_f <- file.path("outputs", "dada2", "seqtab_PacBio_nochim_cleaned_names.rds")
tax_f    <- file.path("outputs", "dada2", "taxonomy_PacBio_SILVA_v138.1.rds")
meta_f   <- file.path("outputs", "metadata_fixed.csv")

stopifnot(file.exists(seqtab_f), file.exists(tax_f), file.exists(meta_f))

seqtab <- readRDS(seqtab_f)
taxa   <- readRDS(tax_f)
meta   <- as.data.frame(fread(meta_f))

stopifnot(!anyDuplicated(colnames(seqtab)), !anyDuplicated(rownames(taxa)))
stopifnot(all(colnames(seqtab) %in% rownames(taxa)))

# Align ASVs
common_asvs <- intersect(colnames(seqtab), rownames(taxa))
seqtab <- seqtab[, common_asvs, drop = FALSE]
taxa   <- taxa[common_asvs, , drop = FALSE]
stopifnot(all(colnames(seqtab) == rownames(taxa)))

# Align samples
stopifnot(METADATA_SAMPLE_COL %in% names(meta))
stopifnot(!anyNA(meta[[METADATA_SAMPLE_COL]]))
stopifnot(!anyDuplicated(meta[[METADATA_SAMPLE_COL]]), !anyDuplicated(rownames(seqtab)))
stopifnot(all(rownames(seqtab) %in% meta[[METADATA_SAMPLE_COL]]))
meta <- meta[meta[[METADATA_SAMPLE_COL]] %in% rownames(seqtab), , drop = FALSE]
meta <- meta[match(rownames(seqtab), meta[[METADATA_SAMPLE_COL]]), , drop = FALSE]
stopifnot(all(meta[[METADATA_SAMPLE_COL]] == rownames(seqtab)))
rownames(meta) <- meta[[METADATA_SAMPLE_COL]]

physeq <- phyloseq(
  otu_table(seqtab, taxa_are_rows = FALSE),
  tax_table(as.matrix(taxa)),
  sample_data(meta)
)
physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)

out_rds <- file.path("outputs", "phyloseq_raw.rds")
saveRDS(physeq, out_rds)

cat("Saved: ", out_rds, "\n", sep = "")
cat("Samples: ", nsamples(physeq), " | ASVs: ", ntaxa(physeq), "\n", sep = "")
