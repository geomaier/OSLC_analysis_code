#!/usr/bin/env Rscript

## =============================================================================
## 150 — Pairwise envfit driver plots (slide-ready)
##
## Purpose
##   For selected pairwise comparisons, compute Aitchison distance (CLR Euclidean)
##   at genus level and use envfit to identify genera whose CLR profiles align
##   with the ordination separation. Produces slide-ready PCoA plots with top
##   significant vectors and exports a drivers table per comparison.
##
## Comparisons
##   (1) SSL_HP vs H_COL_ND
##   (2) H_ADJ_SSL_HP vs H_COL_ND
##
## Key processing steps
##   - Drop zero-library samples
##   - Build Plot_Group by pooling sample types (SSL/HP; H_ADJ_SSL/H_ADJ_HP; etc.)
##   - Agglomerate to Genus
##   - Compute CLR (pseudocount 0.5) and Euclidean distance (Aitchison)
##   - PCoA (cmdscale) on CLR Euclidean distances
##   - envfit vectors on CLR features (9999 permutations)
##   - Keep only vectors with envfit p < 0.01; plot top 5 by (p, r²)
##
## Output
##   - PNG + PDF per comparison
##   - CSV drivers table per comparison (envfit stats + prevalence/median abundance)
##
## Requirements
##   - phyloseq object with sample_data column: New_Sample_Type (or TYPE_COL)
## =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(grid)   # unit()
})

## -------------------------
## 0) INPUTS (GitHub-friendly)
## -------------------------
physeq_path <- file.path("outputs", "phyloseq_after_reseq_noNC_5_percent.rds")
stopifnot(file.exists(physeq_path))

TYPE_COL <- "New_Sample_Type"

out_dir <- file.path("outputs", "beta_envfit_pairwise_drivers")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pairs <- list(
  c("SSL_HP",       "H_COL_ND"),
  c("H_ADJ_SSL_HP", "H_COL_ND")
)

pal <- c(
  "H_COL_ND"     = "#00A65A",
  "H_ADJ_SSL_HP" = "#F39C12",
  "SSL_HP"       = "#E74C3C",
  "H_ADJ_CA"     = "#2980B9",
  "CA"           = "#8E44AD"
)

lab_map <- c(
  "H_COL_ND"     = "H_COL",
  "H_ADJ_SSL_HP" = "H_ADJ_SSL_HP",
  "SSL_HP"       = "SSL_HP",
  "H_ADJ_CA"     = "H_ADJ_CA",
  "CA"           = "CA"
)

TOP_K_VECTORS <- 5
ENVFIT_P_CUT  <- 0.01
ELLIPSE_ALPHA <- 0.05

set.seed(1)

## -------------------------
## Helpers
## -------------------------
build_plot_group <- function(ps, type_col) {
  raw <- trimws(as.character(sample_data(ps)[[type_col]]))
  
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
  
  sample_data(ps)$Plot_Group <- Plot_Group
  ps
}

genus_glom <- function(ps) {
  TT <- as(tax_table(ps), "matrix")
  if (is.null(TT) || ncol(TT) == 0) stop("tax_table is missing or empty.")
  
  rank_use <- if ("Genus" %in% colnames(TT)) "Genus" else colnames(TT)[ncol(TT)]
  pg <- tax_glom(ps, taxrank = rank_use, NArm = FALSE)
  
  TTg <- as(tax_table(pg), "matrix")
  gen <- TTg[, rank_use]
  gen[is.na(gen) | gen == ""] <- "Unclassified"
  
  taxa_names(pg) <- make.unique(as.character(gen))
  pg
}

clr_from_phy <- function(phy_gen, pseudocount = 0.5) {
  otuG <- as(otu_table(phy_gen), "matrix")
  if (taxa_are_rows(phy_gen)) otuG <- t(otuG)  # samples x genera
  otuG <- otuG[, colSums(otuG) > 0, drop = FALSE]
  
  X <- otuG + pseudocount
  X <- X / rowSums(X)
  clr <- log(X) - rowMeans(log(X))
  clr <- as.matrix(clr)
  mode(clr) <- "numeric"
  
  keep_var <- apply(clr, 2, sd, na.rm = TRUE) > 0
  clr <- clr[, keep_var, drop = FALSE]
  
  list(clr = clr, otuG = otuG)
}

make_pair_envfit <- function(ps_in, g1, g2) {
  pair <- c(g1, g2)
  
  ps <- prune_samples(sample_data(ps_in)$Plot_Group %in% pair, ps_in)
  sample_data(ps)$Plot_Group <- factor(as.character(sample_data(ps)$Plot_Group), levels = pair)
  
  phy_gen <- genus_glom(ps)
  
  tmp <- clr_from_phy(phy_gen, pseudocount = 0.5)
  clr <- tmp$clr
  otuG <- tmp$otuG
  dist_obj <- dist(clr, method = "euclidean")
  
  meta <- as(sample_data(phy_gen), "data.frame")
  meta$Group <- factor(as.character(meta$Plot_Group), levels = pair)
  meta$SampleID <- rownames(meta)
  stopifnot(nrow(meta) == attr(dist_obj, "Size"))
  
  ord <- cmdscale(dist_obj, eig = TRUE, k = 2)
  coords <- as.data.frame(ord$points)
  colnames(coords) <- c("PCoA1","PCoA2")
  coords$SampleID <- rownames(coords)
  
  df_ord <- merge(coords, meta, by = "SampleID", all.x = TRUE)
  var_exp <- ord$eig / sum(ord$eig[ord$eig > 0])
  
  fit <- envfit(ord$points, clr, permutations = 9999)
  vec_mat <- scores(fit, display = "vectors")
  if (is.null(vec_mat) || nrow(vec_mat) == 0) stop("envfit returned no vectors.")
  
  vec <- as.data.frame(vec_mat)
  ax1 <- names(vec)[1]; ax2 <- names(vec)[2]
  vec$taxon <- rownames(vec)
  vec$r2 <- fit$vectors$r
  vec$p  <- fit$vectors$pvals
  
  vec <- vec[is.finite(vec$p) & vec$p < ENVFIT_P_CUT, , drop = FALSE]
  vec <- vec[order(vec$p, -vec$r2), , drop = FALSE]
  vec_top <- head(vec, TOP_K_VECTORS)
  
  rng1 <- diff(range(df_ord$PCoA1, na.rm = TRUE))
  rng2 <- diff(range(df_ord$PCoA2, na.rm = TRUE))
  arrow_mult <- 0.45 * min(rng1, rng2)
  
  vec_top$xend <- vec_top[[ax1]] * arrow_mult
  vec_top$yend <- vec_top[[ax2]] * arrow_mult
  
  rel <- otuG / rowSums(otuG)
  rel[!is.finite(rel)] <- 0
  # only keep taxa still present in CLR matrix
  rel <- rel[, colnames(clr), drop = FALSE]
  
  drivers <- vec_top %>%
    transmute(
      taxon = taxon,
      envfit_r2 = r2,
      envfit_p  = p,
      prev_g1 = mean(rel[meta$Group == g1, taxon] > 0, na.rm = TRUE),
      prev_g2 = mean(rel[meta$Group == g2, taxon] > 0, na.rm = TRUE),
      med_g1  = median(rel[meta$Group == g1, taxon], na.rm = TRUE),
      med_g2  = median(rel[meta$Group == g2, taxon], na.rm = TRUE)
    ) %>%
    arrange(envfit_p, desc(envfit_r2))
  
  p <- ggplot(df_ord, aes(PCoA1, PCoA2, color = Group)) +
    stat_ellipse(aes(group = Group), level = 0.95, linewidth = 0.3, alpha = ELLIPSE_ALPHA) +
    geom_point(size = 2.2, alpha = 0.92) +
    scale_color_manual(values = pal[pair], breaks = pair, labels = lab_map[pair]) +
    geom_segment(
      data = vec_top,
      aes(x = 0, y = 0, xend = xend, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.8,
      arrow = arrow(length = unit(2.2, "mm"))
    ) +
    ggrepel::geom_text_repel(
      data = vec_top,
      aes(x = xend, y = yend, label = taxon),
      inherit.aes = FALSE,
      size = 3.5,
      box.padding = 0.25,
      point.padding = 0.20,
      min.segment.length = 0,
      segment.linewidth = 0.25,
      max.overlaps = Inf
    ) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.2,
             label = paste0("→ ", lab_map[g1], "-aligned taxa"), size = 3.2) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2,
             label = paste0("← ", lab_map[g2], "-aligned taxa"), size = 3.2) +
    theme_classic(base_size = 11) +
    theme(
      legend.title = element_blank(),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      plot.title = element_text(face = "bold")
    ) +
    labs(
      title = "Drivers of beta-diversity separation (envfit)",
      subtitle = paste0(lab_map[g1], " vs ", lab_map[g2],
                        " | vectors: top ", TOP_K_VECTORS, " (envfit p < ", ENVFIT_P_CUT, ")"),
      x = paste0("PCoA1 (", round(100*var_exp[1], 1), "%)"),
      y = paste0("PCoA2 (", round(100*var_exp[2], 1), "%)")
    )
  
  list(plot = p, drivers = drivers)
}

## -------------------------
## Run
## -------------------------
phy0 <- readRDS(physeq_path)

otu0 <- as(otu_table(phy0), "matrix")
if (taxa_are_rows(phy0)) otu0 <- t(otu0)
phy0 <- prune_samples(rowSums(otu0) > 0, phy0)

stopifnot(TYPE_COL %in% colnames(as.data.frame(sample_data(phy0))))
phy0 <- build_plot_group(phy0, TYPE_COL)

for (cmp in pairs) {
  g1 <- cmp[1]; g2 <- cmp[2]
  cat("\n--- Running:", g1, "vs", g2, "---\n")
  
  out <- make_pair_envfit(phy0, g1, g2)
  stub <- paste0("PAIR_envfit_", g1, "_vs_", g2, "_top", TOP_K_VECTORS, "_p", ENVFIT_P_CUT)
  
  ggsave(file.path(out_dir, paste0(stub, ".png")),
         out$plot, width = 6.6, height = 4.8, dpi = 600)
  ggsave(file.path(out_dir, paste0(stub, ".pdf")),
         out$plot, width = 6.6, height = 4.8, useDingbats = FALSE)
  
  write.csv(out$drivers,
            file.path(out_dir, paste0(stub, "_drivers.csv")),
            row.names = FALSE)
  
  if (interactive()) {
    try(print(out$plot), silent = TRUE)
  }
}

cat("\nDONE. Saved outputs to:\n", out_dir, "\n")
