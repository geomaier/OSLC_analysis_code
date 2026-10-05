#!/usr/bin/env Rscript
# Run from the analysis directory; local paths are read from config/config.R.
(function() {
set.seed(1)
old_options <- options(stringsAsFactors = FALSE)
on.exit(options(old_options), add = TRUE)

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(RColorBrewer)
  library(grid)
})

if (capabilities("cairo")) {
  options(bitmapType = "cairo")
}

config <- new.env(parent = baseenv())
if (file.exists(file.path("config", "config.R"))) {
  sys.source(file.path("config", "config.R"), envir = config)
}
work_dir <- if (exists("PATH_WORK", envir = config, inherits = FALSE)) {
  get("PATH_WORK", envir = config, inherits = FALSE)
} else if (file.exists("phyloseq_final_March26_5pct_without7699.rds")) getwd() else "outputs"
if (!dir.exists(work_dir)) stop("Input directory does not exist: ", work_dir)
run_id <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
physeq_path <- file.path(work_dir, "phyloseq_final_March26_5pct_without7699.rds")
if (!file.exists(physeq_path)) stop("Input file does not exist: ", physeq_path)
out_dir <- file.path(work_dir, "Section4_5_1_Relative_abundance_FINAL", run_id)
if (dir.exists(out_dir) || !dir.create(out_dir, recursive = TRUE)) stop("Cannot create output directory.")
write.csv(data.frame(File = basename(physeq_path), MD5 = unname(tools::md5sum(physeq_path))),
          file.path(out_dir, "input_checksum.csv"), row.names = FALSE)

group_levels <- c( "SSL_HP", "H_ADJ_SSL_HP", "H_COL" )

group_display <- c( "SSL_HP" = "SSL_HP", "H_ADJ_SSL_HP" = "H_ADJ_SSL_HP", "H_COL" = "H_COL" )

expected_group_counts <- c( "SSL_HP" = 38L, "H_ADJ_SSL_HP" = 12L, "H_COL" = 14L )

top_n <- 12L

font_family <- "Arial"
font_pt <- 7
legend_pt <- 6
strip_pt <- 7
export_dpi <- 600

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "na", "N/A", "Unassigned", "unassigned")] <- NA_character_
  x
}

make_group <- function(x) {
  group <- clean_character(x)
  
  group[
    group %in% c( "SSL", "HP", "SSL_HP" )
  ] <- "SSL_HP"
  
  group[
    group %in% c( "H_ADJ_SSL", "H_ADJ_HP", "H_ADJ_SSL_HP" )
  ] <- "H_ADJ_SSL_HP"

  group[group == "H_COL"] <- "H_COL"
  
  group[!group %in% group_levels] <- NA_character_
  
  factor(group, levels = group_levels)
}

make_lineage_aware_genus <- function(tax_df) {
  genus <- clean_character(tax_df$Genus)
  genus <- gsub("^g__", "", genus)
  genus <- gsub("^Genus:", "", genus)
  genus <- trimws(genus)
  
  parent_ranks <- c("Family", "Order", "Class", "Phylum")
  parent_ranks <- parent_ranks[parent_ranks %in% colnames(tax_df)]
  
  parent <- rep(NA_character_, nrow(tax_df))
  
  for (rank_name in parent_ranks) {
    values <- clean_character(tax_df[[rank_name]])
    use <- is.na(parent) & !is.na(values)
    parent[use] <- values[use]
  }
  
  row_ids <- rownames(tax_df)
  if (is.null(row_ids)) {
    row_ids <- paste0("ASV_", seq_len(nrow(tax_df)))
  }
  
  missing_genus <- is.na(genus)
  
  genus[
    missing_genus & !is.na(parent)
  ] <- paste0( "Unclassified_", gsub("[^[:alnum:]_]+", "_", parent[missing_genus & !is.na(parent)]) )
  
  genus[
    missing_genus & is.na(parent)
  ] <- paste0( "Unclassified_", gsub("[^[:alnum:]_]+", "_", row_ids[missing_genus & is.na(parent)]) )
  
  lineage_ranks <- c("Phylum", "Class", "Order", "Family")
  lineage_ranks <- lineage_ranks[lineage_ranks %in% colnames(tax_df)]
  
  if (length(lineage_ranks) > 0L) {
    lineage_key <- apply(
      tax_df[, lineage_ranks, drop = FALSE],
      1,
      function(z) {
        z <- clean_character(z)
        z[is.na(z)] <- "Unclassified"
        paste(z, collapse = "|")
      }
    )
  } else {
    lineage_key <- rep("Unresolved_lineage", nrow(tax_df))
  }
  
  map <- data.frame(
    Original_genus = genus,
    Lineage_key = lineage_key,
    Parent = ifelse(is.na(parent), "Unresolved", parent),
    stringsAsFactors = FALSE
  ) %>%
    distinct() %>%
    arrange(Original_genus, Lineage_key) %>%
    group_by(Original_genus) %>%
    mutate(
      N_lineages = n(),
      Lineage_number = row_number(),
      Analysis_genus = ifelse(
        N_lineages > 1L,
        paste0(
          Original_genus,
          "__",
          gsub("[^[:alnum:]_]+", "_", Parent),
          "__lineage",
          Lineage_number
        ),
        Original_genus
      )
    ) %>%
    ungroup()
  
  key_map <- paste( map$Original_genus, map$Lineage_key, sep = "\r" )
  
  key_tax <- paste( genus, lineage_key, sep = "\r" )
  
  final_genus <- map$Analysis_genus[
    match(key_tax, key_map)
  ]
  
  if (anyNA(final_genus)) {
    stop("At least one ASV could not be assigned a lineage-aware genus label.")
  }
  
  list( genus = final_genus, mapping = map )
}

save_plot <- function(
    plot,
    stem,
    width_mm,
    height_mm
) {
  pdf_device <- if (capabilities("cairo")) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }
  
  ggsave(
    file.path(out_dir, paste0(stem, ".pdf")),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = pdf_device,
    bg = "white"
  )
  
  ggsave(
    file.path(out_dir, paste0(stem, ".png")),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = export_dpi,
    bg = "white"
  )
  
  ggsave(
    file.path(out_dir, paste0(stem, ".jpg")),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = export_dpi,
    quality = 100,
    bg = "white"
  )
}

if (!file.exists(physeq_path)) {
  stop( "The final phyloseq object was not found:\n", physeq_path )
}

phy <- readRDS(physeq_path)

if (!inherits(phy, "phyloseq")) {
  stop("The input RDS is not a phyloseq object.")
}

if (
  is.null(otu_table(phy, errorIfNULL = FALSE)) ||
  is.null(sample_data(phy, errorIfNULL = FALSE)) ||
  is.null(tax_table(phy, errorIfNULL = FALSE))
) {
  stop("The phyloseq object must contain counts, metadata and taxonomy.")
}

if (!"Sample_Type" %in% colnames(sample_data(phy))) {
  stop("Sample_Type was not found in sample_data(phy).")
}

if (any(grepl("7699", sample_names(phy), fixed = TRUE))) {
  stop(
    "The selected object still contains patient 7699. ",
    "Use phyloseq_final_March26_5pct_without7699.rds."
  )
}

otu_check <- as(otu_table(phy), "matrix")

if (
  any(!is.finite(otu_check)) ||
  any(otu_check < 0) ||
  any(abs(otu_check - round(otu_check)) > 1e-8)
) {
  stop("The input must contain non-negative integer ASV counts.")
}

sample_data(phy)$Group <- make_group( sample_data(phy)$Sample_Type )

phy <- subset_samples( phy, !is.na(Group) )

phy <- prune_taxa( taxa_sums(phy) > 0, phy )

observed_counts <- table( factor( sample_data(phy)$Group, levels = group_levels ) )

if (!all( as.integer(observed_counts) == expected_group_counts[group_levels] )) {
  stop(
    "Unexpected sample counts. Expected ",
    paste( names(expected_group_counts), expected_group_counts, sep = "=", collapse = ", " ),
    "; observed ",
    paste( names(observed_counts), as.integer(observed_counts), sep = "=", collapse = ", " ),
    "."
  )
}

tax_df <- as.data.frame( tax_table(phy), check.names = FALSE, stringsAsFactors = FALSE )

if (!"Genus" %in% colnames(tax_df)) {
  stop("The taxonomy table has no Genus column.")
}

genus_result <- make_lineage_aware_genus(tax_df)
tax_df$Genus <- genus_result$genus

tax_table(phy) <- tax_table( as.matrix(tax_df) )

write.csv(
  genus_result$mapping,
  file.path( out_dir, "Genus_lineage_mapping.csv" ),
  row.names = FALSE
)

write.csv(
  data.frame(
    Input_object = physeq_path,
    Samples = nsamples(phy),
    ASVs = ntaxa(phy),
    New_phyloseq_object_saved = FALSE,
    stringsAsFactors = FALSE
  ),
  file.path( out_dir, "Section4_5_1_input_summary.csv" ),
  row.names = FALSE
)

cat("\nFinal saved phyloseq object loaded:\n")
cat(physeq_path, "\n")
cat("Samples:", nsamples(phy), "\n")
cat("ASVs:", ntaxa(phy), "\n")
cat("No new phyloseq object was written.\n")

phy_genus <- tax_glom( phy, taxrank = "Genus", NArm = FALSE )

phy_genus <- prune_taxa( taxa_sums(phy_genus) > 0, phy_genus )

phy_genus_relative <- transform_sample_counts(
  phy_genus,
  function(x) {
    if (sum(x) == 0) {
      x
    } else {
      x / sum(x)
    }
  }
)

relative_long <- psmelt( phy_genus_relative ) %>%
  mutate(
    Genus = as.character(Genus),
    Group = factor( as.character(Group), levels = group_levels )
  ) %>%
  filter( !is.na(Group), !is.na(Genus), Genus != "" )

all_sample_genus <- expand_grid( Sample = sample_names(phy_genus_relative), Genus = unique(relative_long$Genus) )

sample_metadata <- data.frame( sample_data(phy_genus_relative), check.names = FALSE, stringsAsFactors = FALSE ) %>%
  rownames_to_column("Sample") %>%
  select( Sample, Group )

relative_complete <- all_sample_genus %>%
  left_join( relative_long %>% select( Sample, Genus, Abundance ), by = c("Sample", "Genus") ) %>%
  mutate( Abundance = replace_na(Abundance, 0) ) %>%
  left_join( sample_metadata, by = "Sample" ) %>%
  mutate( Group = factor( as.character(Group), levels = group_levels ) )

genus_summary <- relative_complete %>%
  group_by( Group, Genus ) %>%
  summarise(
    N_samples = n(),
    Mean_relative_abundance = mean(Abundance),
    Median_relative_abundance = median(Abundance),
    Prevalence = mean(Abundance > 0),
    Maximum_relative_abundance = max(Abundance),
    .groups = "drop"
  )

overall_rank <- relative_complete %>%
  group_by(Genus) %>%
  summarise(
    Mean_relative_abundance = mean(Abundance),
    Median_relative_abundance = median(Abundance),
    Prevalence = mean(Abundance > 0),
    Maximum_relative_abundance = max(Abundance),
    .groups = "drop"
  ) %>%
  arrange( desc(Mean_relative_abundance), desc(Prevalence), Genus )

top_genera <- head( overall_rank$Genus, top_n )

write.csv(
  genus_summary,
  file.path( out_dir, "Genus_relative_abundance_and_prevalence_by_group.csv" ),
  row.names = FALSE
)

write.csv(
  overall_rank,
  file.path( out_dir, "Genus_relative_abundance_rank_overall.csv" ),
  row.names = FALSE
)

write.csv(
  data.frame( Rank = seq_along(top_genera), Genus = top_genera, stringsAsFactors = FALSE ),
  file.path( out_dir, "Top12_genera_used_in_figure.csv" ),
  row.names = FALSE
)

plot_data <- relative_complete %>%
  mutate( Taxon_plot = ifelse( Genus %in% top_genera, Genus, "Other" ) ) %>%
  group_by( Group, Sample, Taxon_plot ) %>%
  summarise( Abundance = sum(Abundance), .groups = "drop" )

sample_order <- plot_data %>%
  distinct( Group, Sample ) %>%
  arrange( Group, Sample ) %>%
  pull(Sample)

plot_data$Sample <- factor( plot_data$Sample, levels = sample_order )

group_n <- plot_data %>%
  distinct( Group, Sample ) %>%
  count( Group, name = "N" )

group_n_map <- setNames( group_n$N, as.character(group_n$Group) )

plot_data <- plot_data %>%
  mutate(
    Group_label = paste0( group_display[as.character(Group)], "\n(n = ", group_n_map[as.character(Group)], ")" )
  )

group_label_levels <- paste0( group_display[group_levels], "\n(n = ", group_n_map[group_levels], ")" )

plot_data$Group_label <- factor( plot_data$Group_label, levels = group_label_levels )

plot_data$Taxon_plot <- factor( plot_data$Taxon_plot, levels = c( top_genera, "Other" ) )

taxon_colours <- RColorBrewer::brewer.pal( 12, "Paired" )[seq_along(top_genera)]

fill_map <- setNames( taxon_colours, top_genera )

fill_map["Other"] <- "grey92"

p_sample <- ggplot( plot_data, aes( x = Sample, y = Abundance, fill = Taxon_plot ) ) +
  geom_col( width = 0.92, colour = NA ) +
  facet_grid( ~ Group_label, scales = "free_x", space = "free_x" ) +
  scale_fill_manual( values = fill_map, breaks = c( top_genera, "Other" ), drop = FALSE ) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    breaks = c( 0, 0.25, 0.50, 0.75, 1 ),
    expand = c(0, 0)
  ) +
  labs( x = NULL, y = "Relative abundance", fill = "Top 12 genera overall" ) +
  theme_bw( base_family = font_family, base_size = font_pt ) +
  theme(
    panel.grid = element_blank(),
    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    axis.text.y = element_text( family = font_family, size = font_pt, colour = "black" ),
    
    axis.title.y = element_text( family = font_family, size = font_pt, colour = "black", margin = margin(r = 4) ),
    
    axis.line = element_line( linewidth = 0.35, colour = "black" ),
    
    axis.ticks = element_line( linewidth = 0.35, colour = "black" ),
    
    strip.background = element_rect( fill = "grey90", colour = "grey55", linewidth = 0.30 ),
    
    strip.text = element_text(
      family = font_family,
      face = "bold",
      size = strip_pt,
      colour = "black",
      lineheight = 0.9
    ),
    
    panel.spacing.x = unit( 0.6, "lines" ),
    
    legend.position = "right",
    
    legend.title = element_text( family = font_family, face = "bold", size = font_pt, colour = "black" ),
    
    legend.text = element_text( family = font_family, size = legend_pt, colour = "black" ),
    
    legend.key.size = unit( 4, "mm" ),
    
    legend.background = element_blank(),
    
    plot.background = element_rect( fill = "white", colour = NA ),
    
    panel.background = element_rect( fill = "white", colour = NA ),
    
    plot.margin = margin( 4, 4, 4, 4, "mm" )
  )

save_plot(
  p_sample,
  "Figure_4_5_1_Top12_genus_relative_abundance_by_sample",
  width_mm = 250,
  height_mm = 70
)

group_mean_plot_data <- plot_data %>%
  group_by( Group, Taxon_plot ) %>%
  summarise( Mean_abundance = mean(Abundance), .groups = "drop" ) %>%
  mutate(
    Group = factor( Group, levels = group_levels ),
    Taxon_plot = factor( Taxon_plot, levels = c( top_genera, "Other" ) )
  )

p_group_mean <- ggplot( group_mean_plot_data, aes( x = Group, y = Mean_abundance, fill = Taxon_plot ) ) +
  geom_col( width = 0.70, colour = "black", linewidth = 0.20 ) +
  scale_fill_manual( values = fill_map, breaks = c( top_genera, "Other" ), drop = FALSE ) +
  scale_x_discrete( labels = group_display[group_levels], drop = FALSE ) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    breaks = c( 0, 0.25, 0.50, 0.75, 1 ),
    expand = c(0, 0)
  ) +
  labs( x = NULL, y = "Mean relative abundance", fill = "Top 12 genera overall" ) +
  theme_classic( base_family = font_family, base_size = font_pt ) +
  theme(
    axis.text = element_text( family = font_family, size = font_pt, colour = "black" ),
    
    axis.title = element_text( family = font_family, size = font_pt, colour = "black" ),
    
    legend.title = element_text( family = font_family, face = "bold", size = font_pt ),
    
    legend.text = element_text( family = font_family, size = legend_pt, colour = "black" ),
    
    legend.position = "right",
    
    legend.key.height = unit( 4, "mm" ),
    
    legend.key.width = unit( 4, "mm" ),
    
    legend.box.margin = margin( 3, 3, 3, 3, "mm" ),
    
    plot.margin = margin( 7, 12, 7, 7, "mm" )
  )

save_plot(
  p_group_mean,
  "Supplementary_Top12_genus_mean_relative_abundance_by_group",
  width_mm = 190,
  height_mm = 115
)

writeLines(
  c(
    "SECTION 4.5.1 RELATIVE ABUNDANCE",
    "=================================",
    "",
    paste0("Input final phyloseq object: ", physeq_path),
    paste0("Samples: ", nsamples(phy)),
    paste0("ASVs in saved filtered object: ", ntaxa(phy)),
    "A new phyloseq object was not created or saved.",
    "Genus labels were prepared in memory for this script only.",
    paste0("Top genera shown: ", top_n, " selected by mean abundance across all samples"),
    "Relative-abundance plots are descriptive; no hypothesis test is attached to the stacked bars."
  ),
  file.path( out_dir, "README_Section4_5_1.txt" )
)

capture.output( sessionInfo(), file = file.path( out_dir, "sessionInfo_Section4_5_1.txt" ) )

cat("\nSECTION 4.5.1 COMPLETE\n")
cat("This script used the existing saved object and did not write a new phyloseq object.\n")

write.csv(relative_complete, file.path(out_dir, "sample_genus_relative_abundance.csv"), row.names = FALSE)
write.csv(plot_data, file.path(out_dir, "sample_plot_data.csv"), row.names = FALSE)
write.csv(group_mean_plot_data, file.path(out_dir, "group_mean_plot_data.csv"), row.names = FALSE)
cat("Outputs saved to: ", normalizePath(out_dir), "\n", sep = "")
})()
