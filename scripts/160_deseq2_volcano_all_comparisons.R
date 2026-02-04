#!/usr/bin/env Rscript

## =============================================================================
## 160 — DESeq2 differential abundance (genus-level) + volcano plots
##
## Purpose
##   Run DESeq2 on genus-level counts for multiple pairwise contrasts between
##   pooled sample groups and output (i) results tables and (ii) volcano plots
##   (direction-controlled: “UP” = higher in the first group in the contrast).
##
## Key processing steps
##   - Drop zero-library samples
##   - Build Plot_Group by pooling sample types (SSL/HP; H_ADJ_SSL/H_ADJ_HP; etc.)
##   - Agglomerate to Genus (Unclassified for missing genus labels)
##   - DESeq2 with size factors (type="poscounts") and covariates if present
##   - For each requested comparison:
##       * results table CSV
##       * volcano PNG+PDF
##       * summary of significant taxa (FDR < 0.05)
##
## Output
##   - outputs/deseq2_volcano/DESeq2_<A>_vs_<B>.csv
##   - outputs/deseq2_volcano/Volcano_DESeq2_<A>_vs_<B>.png/.pdf
##   - outputs/deseq2_volcano/SUMMARY_nSig_per_comparison.csv
##
## Requirements
##   - phyloseq object with sample_data column: New_Sample_Type (TYPE_COL)
##   - tax_table contains "Genus" (or last available rank is used)
## =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(grid)   # unit()
})

## -------------------------------
## 0) INPUTS (GitHub-friendly)
## -------------------------------
physeq_path <- file.path("outputs", "phyloseq_after_reseq_noNC_5_percent.rds")
stopifnot(file.exists(physeq_path))

TYPE_COL <- "New_Sample_Type"

out_dir <- file.path("outputs", "deseq2_volcano")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

groups_universe <- c("H_COL_ND","H_ADJ_SSL_HP","SSL_HP","H_ADJ_CA","CA")

COL_UP   <- "red"     # UP = higher in FIRST group in contrast
COL_DOWN <- "green"
COL_NS   <- "grey70"

set.seed(1)

## -------------------------------
## 1) LOAD + DROP ZERO-LIBRARY
## -------------------------------
phy0 <- readRDS(physeq_path)

otu0 <- as(otu_table(phy0), "matrix")
if (taxa_are_rows(phy0)) otu0 <- t(otu0)
phy0 <- prune_samples(rowSums(otu0) > 0, phy0)

SD0 <- as.data.frame(sample_data(phy0))
if (!(TYPE_COL %in% colnames(SD0))) {
  stop(paste0("Column '", TYPE_COL, "' not found. Available: ",
              paste(colnames(SD0), collapse = ", ")))
}

## -------------------------------
## 2) BUILD Plot_Group (pooling logic)
## -------------------------------
raw <- trimws(as.character(sample_data(phy0)[[TYPE_COL]]))
canon <- toupper(raw)
canon <- gsub("[[:space:]]+", "_", canon)
canon <- gsub("-", "_", canon)

Plot_Group <- canon
Plot_Group[Plot_Group %in% c("SSL","HP","SSL_HP")] <- "SSL_HP"
Plot_Group[Plot_Group %in% c("H_ADJ_SSL","H_ADJ_HP","H_ADJ_SSL_HP","H_ADJ")] <- "H_ADJ_SSL_HP"
Plot_Group[Plot_Group %in% c("H_COL_ND","H_COL","H_COLND")] <- "H_COL_ND"
Plot_Group[Plot_Group %in% c("CA","CONVENTIONAL_ADENOMA","CONVENTIONAL_ADENOMAS")] <- "CA"
Plot_Group[grepl("H_ADJ_CA", Plot_Group)] <- "H_ADJ_CA"
Plot_Group[grepl("ADJ", Plot_Group) & grepl("CA", Plot_Group)] <- "H_ADJ_CA"
Plot_Group[grepl("ADEN", Plot_Group) & !grepl("ADJ", Plot_Group)] <- "CA"

sample_data(phy0)$Plot_Group <- Plot_Group

phy0 <- prune_samples(sample_data(phy0)$Plot_Group %in% groups_universe, phy0)
sample_data(phy0)$Plot_Group <- factor(as.character(sample_data(phy0)$Plot_Group),
                                       levels = groups_universe)

cat("\nN per Plot_Group:\n")
print(table(sample_data(phy0)$Plot_Group))

## -------------------------------
## 3) GENUS-LEVEL AGGLOMERATION
## -------------------------------
TT <- as(tax_table(phy0), "matrix")
if (is.null(TT) || ncol(TT) == 0) stop("tax_table is missing or empty.")
rank_use <- if ("Genus" %in% colnames(TT)) "Genus" else colnames(TT)[ncol(TT)]

phy_gen <- tax_glom(phy0, taxrank = rank_use, NArm = FALSE)

TTg <- as(tax_table(phy_gen), "matrix")
tax_label <- TTg[, rank_use]
tax_label[is.na(tax_label) | tax_label == ""] <- "Unclassified"
taxa_names(phy_gen) <- make.unique(as.character(tax_label))

## -------------------------------
## 4) DESeq2 model (poscounts) + safe covariates
## -------------------------------
metaD <- as.data.frame(sample_data(phy_gen), stringsAsFactors = FALSE)

covars <- c("Batch","Location_Colon","BRAF","GENDER")
covars <- covars[covars %in% colnames(metaD)]

keep_cov <- character(0)
for (v in covars) {
  metaD[[v]] <- droplevels(as.factor(metaD[[v]]))
  if (!any(is.na(metaD[[v]])) && nlevels(metaD[[v]]) >= 2) keep_cov <- c(keep_cov, v)
}

metaD$Plot_Group <- factor(as.character(metaD$Plot_Group), levels = groups_universe)

design_vars <- c("Plot_Group", keep_cov)
ok <- complete.cases(metaD[, design_vars, drop = FALSE])
if (sum(!ok) > 0) {
  cat("\nDropping samples with NA in design vars:", sum(!ok), "\n")
  phy_gen <- prune_samples(rownames(metaD)[ok], phy_gen)
  metaD <- metaD[ok, , drop = FALSE]
}
sample_data(phy_gen) <- sample_data(metaD)

design_fml <- if (length(keep_cov) > 0) {
  as.formula(paste("~", paste(keep_cov, collapse = " + "), "+ Plot_Group"))
} else {
  ~ Plot_Group
}
cat("\nDESeq2 design:", deparse(design_fml), "\n")

dds <- phyloseq_to_deseq2(phy_gen, design_fml)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds, sfType = "poscounts", quiet = TRUE)

## -------------------------------
## 5) Volcano plot function
## -------------------------------
make_volcano <- function(res_df, title, out_prefix,
                         padj_cut = 0.05,
                         label_top_n = 25,
                         lfc_cap = 30) {
  df <- res_df
  df$taxon <- rownames(df)
  
  df$padj2 <- df$padj
  df$mlog10 <- -log10(pmax(df$padj2, 1e-300))
  
  df$LFC <- df$log2FoldChange
  df$LFC[df$LFC >  lfc_cap] <-  lfc_cap
  df$LFC[df$LFC < -lfc_cap] <- -lfc_cap
  
  df$sig <- is.finite(df$padj2) & df$padj2 < padj_cut
  
  df$dir <- "NS"
  df$dir[df$sig & df$log2FoldChange > 0] <- "UP"
  df$dir[df$sig & df$log2FoldChange < 0] <- "DOWN"
  
  lab <- df[df$sig & is.finite(df$mlog10) & is.finite(df$LFC), ]
  lab <- lab[order(lab$padj2), , drop = FALSE]
  if (nrow(lab) > label_top_n) lab <- lab[1:label_top_n, , drop = FALSE]
  
  p <- ggplot(df, aes(x = LFC, y = mlog10)) +
    geom_point(aes(color = dir), size = 1.8, alpha = 0.85) +
    geom_vline(xintercept = 0, linewidth = 0.35) +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", linewidth = 0.35) +
    scale_color_manual(
      values = c("UP" = COL_UP, "DOWN" = COL_DOWN, "NS" = COL_NS),
      breaks = c("UP","DOWN"),
      labels = c("Up in first group", "Down in first group")
    ) +
    labs(
      title = title,
      x = "log2 fold-change (first / second)",
      y = expression(-log[10]("FDR (BH)")),
      color = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "right",
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      plot.title = element_text(face = "bold"),
      legend.key.height = unit(0.6, "lines")
    )
  
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = lab,
      aes(label = taxon),
      size = 3.0,
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0,
      segment.linewidth = 0.25,
      max.overlaps = Inf
    )
  } else {
    p <- p + geom_text(data = lab, aes(label = taxon),
                       size = 3.0, vjust = -0.4, check_overlap = TRUE)
  }
  
  ggsave(file.path(out_dir, paste0(out_prefix, ".png")),
         p, width = 6.5, height = 4.8, dpi = 600)
  ggsave(file.path(out_dir, paste0(out_prefix, ".pdf")),
         p, width = 6.5, height = 4.8, useDingbats = FALSE)
  
  invisible(p)
}

## -------------------------------
## 6) DEFINE COMPARISONS (direction-controlled)
## -------------------------------
comparisons <- list(
  c("H_ADJ_SSL_HP", "H_COL_ND"),
  c("H_COL_ND",     "SSL_HP"),
  c("H_COL_ND",     "CA"),
  c("H_ADJ_CA",     "CA"),
  c("CA",           "SSL_HP"),
  c("H_COL_ND",     "H_ADJ_CA"),
  c("H_COL_ND",     "H_ADJ_SSL_HP"),
  c("SSL_HP",       "H_ADJ_CA")
)

lvl <- levels(sample_data(phy_gen)$Plot_Group)
comparisons <- Filter(function(x) all(x %in% lvl), comparisons)

cat("\nComparisons to run:\n")
print(do.call(rbind, comparisons))

## -------------------------------
## 7) RUN COMPARISONS
## -------------------------------
sum_df <- data.frame()

for (cmp in comparisons) {
  gA <- cmp[1]; gB <- cmp[2]
  cat("\n--- Running contrast:", gA, "vs", gB, "---\n")
  
  res <- results(dds, contrast = c("Plot_Group", gA, gB))
  res_df <- as.data.frame(res)
  
  csv_name <- paste0("DESeq2_", gA, "_vs_", gB, ".csv")
  write.csv(data.frame(taxon = rownames(res_df), res_df),
            file.path(out_dir, csv_name),
            row.names = FALSE)
  
  title <- paste0("Differentially abundant genera (DESeq2)\n", gA, " vs ", gB)
  prefix <- paste0("Volcano_DESeq2_", gA, "_vs_", gB)
  
  make_volcano(res_df, title, prefix, padj_cut = 0.05, label_top_n = 25, lfc_cap = 30)
  
  n_sig <- sum(is.finite(res_df$padj) & res_df$padj < 0.05, na.rm = TRUE)
  sum_df <- rbind(sum_df, data.frame(
    Comparison = paste0(gA, " vs ", gB),
    n_sig_FDR05 = n_sig,
    stringsAsFactors = FALSE
  ))
}

write.csv(sum_df, file.path(out_dir, "SUMMARY_nSig_per_comparison.csv"), row.names = FALSE)

cat("\nDONE.\nVolcano plots + CSVs written to:\n", out_dir, "\n", sep = "")
print(sum_df)
