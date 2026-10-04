#!/usr/bin/env Rscript

## =============================================================================
## 140 — Beta diversity: pairwise PCoA (Aitchison / CLR Euclidean) + PERMANOVA
##
## Purpose
##   Generate pairwise PCoA plots between pre-defined groups using Aitchison
##   distance (Euclidean distance on CLR-transformed relative abundances),
##   and report PERMANOVA + betadisper p-values for each comparison.
##
## Key processing steps
##   - Drop zero-library samples
##   - Create Plot_Group by pooling sample types:
##       * SSL + HP -> SSL_HP
##       * H_ADJ_SSL + H_ADJ_HP -> H_ADJ_SSL_HP
##       * H_COL_ND kept as H_COL_ND (displayed as H_COL)
##       * CA and H_ADJ_CA inferred from labels
##   - For each requested pairwise comparison:
##       * compute CLR Euclidean distances
##       * PERMANOVA (adonis2) and betadisper/permutest
##       * PCoA (cmdscale) and plot with 95% ellipses
##
## Output
##   - One PNG + one PDF per comparison
##   - CSV summary of PERMANOVA R²/p and betadisper p
##
## Requirements
##   - phyloseq object with sample_data column: New_Sample_Type (or TYPE_COL)
## =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(dplyr)
})

## -------------------------------
## 0) INPUTS (GitHub-friendly)
## -------------------------------
physeq_path <- file.path("outputs", "phyloseq_after_reseq_noNC_5_percent.rds")
stopifnot(file.exists(physeq_path))

TYPE_COL <- "New_Sample_Type"

out_dir <- file.path("outputs", "beta_pairwise_PCoA_exports")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Analysis groups
groups_universe <- c("H_COL_ND","H_ADJ_SSL_HP","SSL_HP","H_ADJ_CA","CA")

## Colours (keys are analysis group names)
pal <- c(
  "H_COL_ND"     = "#00A65A",
  "H_ADJ_SSL_HP" = "#F39C12",
  "SSL_HP"       = "#E74C3C",
  "H_ADJ_CA"     = "#2980B9",
  "CA"           = "#8E44AD"
)

## Display labels (H_COL_ND shown as H_COL)
lab_map <- c(
  "H_COL_ND"     = "H_COL",
  "H_ADJ_SSL_HP" = "H_ADJ_SSL_HP",
  "SSL_HP"       = "SSL_HP",
  "H_ADJ_CA"     = "H_ADJ_CA",
  "CA"           = "CA"
)

## Requested pairwise comparisons
comparisons <- list(
  c("SSL_HP", "H_ADJ_SSL_HP"),
  c("SSL_HP", "H_COL_ND"),
  c("CA",     "H_COL_ND"),
  c("CA",     "H_ADJ_CA"),
  c("SSL_HP", "CA"),
  c("H_ADJ_SSL_HP", "H_COL_ND")
)

## -------------------------------
## 1) LOAD + DROP ZERO-LIBRARY
## -------------------------------
physeq <- readRDS(physeq_path)

SD <- as.data.frame(sample_data(physeq))
if (!(TYPE_COL %in% colnames(SD))) {
  stop(paste0("Column '", TYPE_COL, "' not found. Available: ",
              paste(colnames(SD), collapse = ", ")))
}

otu <- as(otu_table(physeq), "matrix")
if (taxa_are_rows(physeq)) otu <- t(otu)
libsize <- rowSums(otu)
cat("Zero-library samples removed:", sum(libsize == 0), "\n")
physeq <- prune_samples(libsize > 0, physeq)

## -------------------------------
## 2) BUILD POOLED GROUP COLUMN: Plot_Group
## -------------------------------
raw <- trimws(as.character(sample_data(physeq)[[TYPE_COL]]))

canon <- toupper(raw)
canon <- gsub("[[:space:]]+", "_", canon)
canon <- gsub("-", "_", canon)

Plot_Group <- canon

Plot_Group[Plot_Group %in% c("SSL","HP","SSL_HP")] <- "SSL_HP"
Plot_Group[Plot_Group %in% c("H_ADJ_SSL","H_ADJ_HP","H_ADJ_SSL_HP","H_ADJ")] <- "H_ADJ_SSL_HP"
Plot_Group[Plot_Group %in% c("H_COL_ND","H_COL","H_COLND","H_COL_ND_")] <- "H_COL_ND"

Plot_Group[Plot_Group %in% c("CA","CONVENTIONAL_ADENOMA","CONVENTIONAL_ADENOMAS")] <- "CA"
Plot_Group[grepl("H_ADJ_CA", Plot_Group)] <- "H_ADJ_CA"
Plot_Group[grepl("ADJ", Plot_Group) & grepl("CA", Plot_Group)] <- "H_ADJ_CA"
Plot_Group[grepl("ADEN", Plot_Group) & !grepl("ADJ", Plot_Group)] <- "CA"

sample_data(physeq)$Plot_Group <- Plot_Group

cat("\nCounts after pooling (Plot_Group):\n")
print(table(sample_data(physeq)$Plot_Group, useNA = "ifany"))

phy_all <- prune_samples(sample_data(physeq)$Plot_Group %in% groups_universe, physeq)
sample_data(phy_all)$Plot_Group <- factor(as.character(sample_data(phy_all)$Plot_Group),
                                          levels = groups_universe)

cat("\nN per universe group:\n")
print(table(sample_data(phy_all)$Plot_Group))

## -------------------------------
## 3) Helper: one pairwise plot
## -------------------------------
make_pair_plot <- function(phy_in, g1, g2, out_dir, pal, lab_map) {
  
  pair <- c(g1, g2)
  
  phy <- prune_samples(sample_data(phy_in)$Plot_Group %in% pair, phy_in)
  sample_data(phy)$Plot_Group <- factor(as.character(sample_data(phy)$Plot_Group), levels = pair)
  
  # Aitchison distance (CLR Euclidean)
  otu2 <- as(otu_table(phy), "matrix")
  if (taxa_are_rows(phy)) otu2 <- t(otu2)
  otu2 <- otu2[, colSums(otu2) > 0, drop = FALSE]
  
  X <- otu2 + 0.5
  X <- X / rowSums(X)
  clr <- log(X) - rowMeans(log(X))
  dist_obj <- dist(clr, method = "euclidean")
  
  meta <- data.frame(sample_data(phy), stringsAsFactors = FALSE)
  meta$Group <- factor(as.character(meta$Plot_Group), levels = pair)
  meta$SampleID <- rownames(meta)
  stopifnot(nrow(meta) == attr(dist_obj, "Size"))
  
  set.seed(1)
  perm <- adonis2(dist_obj ~ Group, data = meta, permutations = 9999)
  
  bd <- betadisper(dist_obj, group = meta$Group)
  bd_perm <- permutest(bd, permutations = 9999)
  
  ord <- cmdscale(dist_obj, eig = TRUE, k = 2)
  coords <- as.data.frame(ord$points)
  colnames(coords) <- c("PCoA1","PCoA2")
  coords$SampleID <- rownames(coords)
  
  df <- merge(coords, meta, by = "SampleID", all.x = TRUE)
  var_exp <- ord$eig / sum(ord$eig[ord$eig > 0])
  
  stat_txt <- paste0(
    "PERMANOVA: R²=", sprintf("%.3f", perm$R2[1]),
    ", p=", format.pval(perm$`Pr(>F)`[1], eps = 1e-3),
    "\n",
    "betadisper p=", format.pval(bd_perm$tab[1,"Pr(>F)"], eps = 1e-3)
  )
  
  p <- ggplot(df, aes(PCoA1, PCoA2)) +
    stat_ellipse(aes(group = Group, fill = Group),
                 geom = "polygon", type = "t", level = 0.95,
                 alpha = 0.22, color = NA) +
    geom_point(aes(color = Group), size = 2.1, alpha = 0.9) +
    scale_color_manual(values = pal[pair], breaks = pair, labels = lab_map[pair], drop = FALSE) +
    scale_fill_manual(values = pal[pair], breaks = pair, labels = lab_map[pair], drop = FALSE) +
    labs(
      title = paste0("PCoA on CLR Euclidean distances\n", lab_map[g1], " vs ", lab_map[g2]),
      x = paste0("PCoA1 (", round(100*var_exp[1], 1), "%)"),
      y = paste0("PCoA2 (", round(100*var_exp[2], 1), "%)")
    ) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.03, vjust = -0.4,
             label = stat_txt, size = 3.2) +
    theme_classic(base_size = 11) +
    theme(
      legend.title = element_blank(),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      plot.title = element_text(size = 11)
    )
  
  fn_stub <- paste0("PCoA_pair_", g1, "_vs_", g2)
  ggsave(file.path(out_dir, paste0(fn_stub, ".png")), p, width = 5.4, height = 4.2, dpi = 600)
  ggsave(file.path(out_dir, paste0(fn_stub, ".pdf")), p, width = 5.4, height = 4.2)
  
  list(perm = perm, bd_p = bd_perm$tab[1,"Pr(>F)"])
}

## -------------------------------
## 4) Run all comparisons + summary CSV
## -------------------------------
sum_df <- data.frame()

for (cmp in comparisons) {
  g1 <- cmp[1]; g2 <- cmp[2]
  cat("\n--- Running:", g1, "vs", g2, "---\n")
  
  out <- make_pair_plot(phy_all, g1, g2, out_dir, pal, lab_map)
  
  sum_df <- rbind(sum_df, data.frame(
    Comparison = paste0(lab_map[g1], " vs ", lab_map[g2]),
    R2 = unname(out$perm$R2[1]),
    p_PERMANOVA = unname(out$perm$`Pr(>F)`[1]),
    p_betadisper = unname(out$bd_p),
    stringsAsFactors = FALSE
  ))
}

write.csv(sum_df, file.path(out_dir, "Pairwise_PERMANOVA_betadisper_summary.csv"),
          row.names = FALSE)

cat("\nDONE.\nSaved plots + summary here:\n", out_dir, "\n")
print(sum_df)
