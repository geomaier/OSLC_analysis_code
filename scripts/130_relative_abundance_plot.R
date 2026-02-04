#!/usr/bin/env Rscript

## =============================================================================
## 130 — Figure 4B: Genus-level relative abundance (Top 12 overall)
##
## Purpose
##   Create a genus-level stacked bar plot of relative abundance using the Top 12
##   genera ranked by mean relative abundance across all retained samples.
##
## Key processing steps
##   - Remove positive controls (sample IDs ending in "_PC") before ranking
##   - Drop samples with missing/blank New_Sample_Type
##   - Keep groups: H_COL_ND, H_ADJ_SSL, SSL, HP, CA, H_ADJ_CA
##   - Pool SSL and HP into "SSL/HP" (exact label)
##   - Agglomerate to Genus; transform to relative abundance per sample
##   - Rank genera globally; plot Top 12 + "Other"
##   - Order samples within each facet by Fusobacterium abundance (descending)
##
## Output
##   - PDF and PNG saved to outputs/figures (configurable below)
##
## Requirements
##   - phyloseq object with sample_data column: New_Sample_Type
##   - tax_table contains: Genus
## =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(grid)   # unit()
})

## -------------------------------
## 0) INPUTS / OUTPUTS (GitHub-friendly: repo-relative)
## -------------------------------
phyloseq_in <- file.path("outputs", "phyloseq_after_reseq_noNC_5_percent.rds")
stopifnot(file.exists(phyloseq_in))

out_dir <- file.path("outputs", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_pdf <- file.path(out_dir, "Fig4B_REL_top12_OVERALL_PAIRED_noPC_noXlabels.pdf")
out_png <- file.path(out_dir, "Fig4B_REL_top12_OVERALL_PAIRED_noPC_noXlabels.png")

## -------------------------------
## 1) Load phyloseq
## -------------------------------
phy <- readRDS(phyloseq_in)

## -------------------------------
## 2) Remove positive controls (PC) FIRST
## -------------------------------
meta0 <- data.frame(sample_data(phy))
meta0$sample <- rownames(meta0)
pc_samples <- meta0$sample[grepl("_PC$", meta0$sample)]
cat("PC samples removed:", length(pc_samples), "\n")

phy_bio <- subset_samples(phy, !(sample_names(phy) %in% pc_samples))
phy_bio <- prune_taxa(taxa_sums(phy_bio) > 0, phy_bio)

## -------------------------------
## 3) Clean sample type + drop blanks
## -------------------------------
sample_data(phy_bio)$New_Sample_Type <- trimws(as.character(sample_data(phy_bio)$New_Sample_Type))
sample_data(phy_bio)$New_Sample_Type[is.na(sample_data(phy_bio)$New_Sample_Type)] <- ""

phy_bio <- subset_samples(phy_bio, New_Sample_Type != "")
phy_bio <- prune_taxa(taxa_sums(phy_bio) > 0, phy_bio)

## -------------------------------
## 4) Keep groups + pool SSL and HP into SSL/HP (exact)
## -------------------------------
groups_keep <- c("H_COL_ND", "H_ADJ_SSL", "SSL", "HP", "CA", "H_ADJ_CA")
phy_bio <- subset_samples(phy_bio, New_Sample_Type %in% groups_keep)
phy_bio <- prune_taxa(taxa_sums(phy_bio) > 0, phy_bio)

sd <- data.frame(sample_data(phy_bio))
sd$Group <- as.character(sd$New_Sample_Type)
sd$Group[sd$Group %in% c("SSL", "HP")] <- "SSL/HP"
sd$Group <- factor(sd$Group, levels = c("H_COL_ND","H_ADJ_SSL","SSL/HP","CA","H_ADJ_CA"))
sample_data(phy_bio)$Group <- sd$Group

cat("Group counts:\n")
print(table(sample_data(phy_bio)$Group))

## -------------------------------
## 5) Genus agglomeration + relative abundance
## -------------------------------
tax <- as.data.frame(tax_table(phy_bio))
if (!("Genus" %in% colnames(tax))) stop("No 'Genus' column in tax_table().")

tax$Genus <- as.character(tax$Genus)
tax$Genus[is.na(tax$Genus) | tax$Genus == "" | tax$Genus == "NA"] <- "Unclassified"
tax_table(phy_bio) <- tax_table(as.matrix(tax))

phy_gen <- tax_glom(phy_bio, taxrank = "Genus", NArm = FALSE)
phy_gen_rel <- transform_sample_counts(phy_gen, function(x) x / sum(x))

df <- psmelt(phy_gen_rel)
df$Genus <- as.character(df$Genus)
df$Group <- as.character(df$Group)

## -------------------------------
## 6) Top 12 genera overall (mean rel abundance across all samples)
## -------------------------------
top_n <- 12

rank_all <- df %>%
  group_by(Genus) %>%
  summarise(mean_all = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_all))

top_genera <- rank_all$Genus[seq_len(min(top_n, nrow(rank_all)))]
cat("Top genera used (overall, n=", length(top_genera), "):\n", paste(top_genera, collapse=", "), "\n")

df$Taxon_plot <- ifelse(df$Genus %in% top_genera, df$Genus, "Other")

df2 <- df %>%
  group_by(Group, Sample, Taxon_plot) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

## -------------------------------
## 7) Order samples within each facet by Fusobacterium abundance (desc)
##    (ordering only; Fusobacterium may be outside Top 12)
## -------------------------------
fus_raw <- df %>%
  filter(Genus == "Fusobacterium") %>%
  group_by(Group, Sample) %>%
  summarise(fuso = sum(Abundance, na.rm = TRUE), .groups = "drop")

all_samp <- df2 %>% distinct(Group, Sample)
fus_raw <- merge(all_samp, fus_raw, by = c("Group","Sample"), all.x = TRUE)
fus_raw$fuso[is.na(fus_raw$fuso)] <- 0

ord <- fus_raw[order(factor(fus_raw$Group, levels = levels(sd$Group)), -fus_raw$fuso), ]
df2$Sample <- factor(as.character(df2$Sample), levels = ord$Sample)

## -------------------------------
## 8) Facet labels with (n=)
## -------------------------------
n_df  <- df2 %>% distinct(Group, Sample) %>% count(Group, name = "n")
n_map <- setNames(n_df$n, as.character(n_df$Group))

df2$Group_lab <- paste0(as.character(df2$Group), "\n(n=", n_map[as.character(df2$Group)], ")")
df2$Group_lab <- factor(
  df2$Group_lab,
  levels = paste0(levels(sd$Group), "\n(n=", n_map[levels(sd$Group)], ")")
)

## -------------------------------
## 9) Colours: Brewer "Paired" for Top 12 + Other grey
## -------------------------------
paired_cols <- brewer.pal(12, "Paired")

df2$Taxon_plot <- factor(df2$Taxon_plot, levels = c(top_genera, "Other"))

fill_map <- setNames(paired_cols[seq_along(top_genera)], top_genera)
fill_map["Other"] <- "grey95"

## -------------------------------
## 10) Plot (relative abundance, no x labels)
## -------------------------------
pB <- ggplot(df2, aes(x = Sample, y = Abundance, fill = Taxon_plot)) +
  geom_col(width = 0.90) +
  facet_grid(~ Group_lab, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = fill_map, drop = FALSE) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "Relative abundance", fill = "Taxa") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "grey85", colour = "grey60", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 9, lineheight = 0.9),
    panel.spacing.x = unit(0.6, "lines"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(4, "mm")
  )

print(pB)

## -------------------------------
## 11) Save
## -------------------------------
ggsave(out_pdf, pB, width = 14.5, height = 3.6, units = "in", device = cairo_pdf)
ggsave(out_png, pB, width = 14.5, height = 3.6, units = "in", dpi = 600)

cat("Saved:\n", out_pdf, "\n", out_png, "\n")
