#!/usr/bin/env Rscript
# Genus-level Bray-Curtis analysis. Run from the 16S directory.

set.seed(1)

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(vegan)
  library(tibble)
  library(ggplot2)
  library(grid)
})

if (capabilities("cairo")) {
  options(bitmapType = "cairo")
}

# Paths and settings

physeq_path <- file.path(
  "outputs",
  "phyloseq_final_March26_5pct_without7699.rds"
)

out_dir <- file.path(
  "outputs",
  "Beta_diversity_BrayCurtis_Genus_FINAL",
  format(Sys.time(), "%Y%m%d_%H%M%S")
)

if (dir.exists(out_dir)) stop("Output directory already exists: ", out_dir)
if (!dir.create(out_dir, recursive = TRUE)) stop("Could not create: ", out_dir)

patient_col_manual <- NULL

n_permutations <- 9999L
ellipse_level <- 0.68
export_dpi <- 600

font_family <- "Arial"
font_pt <- 7
pt_to_gg <- function(x) x / 2.845276

group_levels <- c(
  "H_COL",
  "H_ADJ_SSL_HP",
  "SSL_HP"
)

group_display <- c(
  "H_COL" = "H_COL",
  "H_ADJ_SSL_HP" = "H_ADJ_SSL_HP",
  "SSL_HP" = "SSL_HP"
)

expected_group_counts <- c(
  "H_COL" = 14L,
  "H_ADJ_SSL_HP" = 12L,
  "SSL_HP" = 38L
)

pal <- c(
  "H_COL" = "#1313EC",
  "H_ADJ_SSL_HP" = "#E9C61D",
  "SSL_HP" = "#FF0000"
)

pal_fill <- c(
  "H_COL" = "#AEB8FF",
  "H_ADJ_SSL_HP" = "#F6E79A",
  "SSL_HP" = "#FFB3B3"
)

pairwise_plot_order <- c(
  "H_ADJ_SSL_HP vs H_COL",
  "SSL_HP vs H_COL",
  "SSL_HP vs H_ADJ_SSL_HP"
)

# General helpers

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "na", "N/A", "Unassigned", "unassigned")] <-
    NA_character_
  x
}

first_existing_column <- function(column_names, candidates) {
  hit <- candidates[candidates %in% column_names]
  if (length(hit) == 0L) return(NULL)
  hit[1]
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("_+$|^_+", "", x)
}

format_p <- function(p, prefix = "P") {
  if (length(p) == 0L || is.na(p)) {
    return(paste0(prefix, " = NA"))
  }

  if (p < 0.001) {
    return(paste0(prefix, " < 0.001"))
  }

  paste0(
    prefix,
    " = ",
    formatC(p, format = "f", digits = 3)
  )
}

make_group <- function(x) {
  group <- clean_character(x)

  group[
    group %in% c(
      "SSL",
      "HP",
      "SSL_HP"
    )
  ] <- "SSL_HP"

  group[
    group %in% c(
      "H_ADJ_SSL",
      "H_ADJ_HP",
      "H_ADJ_SSL_HP"
    )
  ] <- "H_ADJ_SSL_HP"

  group[
    group %in% c(
      "H_COL",
      "H_COL_D",
      "H_COL_ND"
    )
  ] <- "H_COL"

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
  ] <- paste0(
    "Unclassified_",
    gsub(
      "[^[:alnum:]_]+",
      "_",
      parent[missing_genus & !is.na(parent)]
    )
  )

  genus[
    missing_genus & is.na(parent)
  ] <- paste0(
    "Unclassified_",
    gsub(
      "[^[:alnum:]_]+",
      "_",
      row_ids[missing_genus & is.na(parent)]
    )
  )

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

  key_map <- paste(
    map$Original_genus,
    map$Lineage_key,
    sep = "\r"
  )

  key_tax <- paste(
    genus,
    lineage_key,
    sep = "\r"
  )

  final_genus <- map$Analysis_genus[
    match(key_tax, key_map)
  ]

  if (anyNA(final_genus)) {
    stop(
      "At least one ASV could not be assigned a lineage-aware genus label."
    )
  }

  list(
    genus = final_genus,
    mapping = map
  )
}

# Metadata, distance and model helpers

get_meta <- function(phy) {
  meta <- data.frame(
    sample_data(phy),
    check.names = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    rownames_to_column("SampleID")

  meta$Group <- factor(
    as.character(meta$Group),
    levels = group_levels
  )

  patient_candidates <- c(
    "Patient",
    "patient",
    "Patient_ID",
    "patient_id",
    "PatientID",
    "Subject",
    "Subject_ID",
    "subject_id"
  )

  if (!is.null(patient_col_manual)) {
    if (!patient_col_manual %in% colnames(meta)) {
      stop(
        "patient_col_manual was not found in the metadata: ",
        patient_col_manual
      )
    }

    patient_column <- patient_col_manual
  } else {
    patient_column <- first_existing_column(
      colnames(meta),
      patient_candidates
    )
  }

  if (!is.null(patient_column)) {
    patient <- clean_character(meta[[patient_column]])
    patient_source <- patient_column
  } else {
    patient <- sub(
      ".*_([0-9]+)$",
      "\\1",
      meta$SampleID
    )

    patient[patient == meta$SampleID] <- NA_character_
    patient_source <- "Final numeric token in SampleID"
  }

  if (anyNA(patient)) {
    bad_samples <- meta$SampleID[is.na(patient)]

    stop(
      "Patient IDs could not be resolved for: ",
      paste(bad_samples, collapse = ", "),
      ". Set patient_col_manual to the correct metadata column."
    )
  }

  meta$Patient <- factor(patient)
  meta$Reads <- as.numeric(sample_sums(phy))
  meta$Patient_source <- patient_source

  rownames(meta) <- meta$SampleID
  meta
}

get_bray_curtis <- function(phy) {
  otu <- as(otu_table(phy), "matrix")

  if (taxa_are_rows(phy)) {
    otu <- t(otu)
  }

  otu <- otu[, colSums(otu) > 0, drop = FALSE]

  if (
    any(!is.finite(otu)) ||
    any(otu < 0)
  ) {
    stop("The genus count matrix contains negative or non-finite values.")
  }

  if (any(rowSums(otu) <= 0)) {
    stop("At least one sample has zero total abundance.")
  }

  relative_abundance <- sweep(
    otu,
    1,
    rowSums(otu),
    "/"
  )

  bray_distance <- vegan::vegdist(
    relative_abundance,
    method = "bray",
    binary = FALSE
  )

  list(
    counts = otu,
    relative_abundance = relative_abundance,
    distance = bray_distance
  )
}

safe_adonis2 <- function(
    distance_object,
    metadata,
    use_patient_strata = FALSE
) {
  set.seed(1)

  result <- tryCatch(
    {
      if (use_patient_strata) {
        vegan::adonis2(
          distance_object ~ Group,
          data = metadata,
          permutations = n_permutations,
          by = NULL,
          strata = metadata$Patient
        )
      } else {
        vegan::adonis2(
          distance_object ~ Group,
          data = metadata,
          permutations = n_permutations,
          by = NULL
        )
      }
    },
    error = function(e) e
  )

  if (inherits(result, "error")) {
    stop(
      "PERMANOVA failed",
      if (use_patient_strata) {
        " with permutations restricted within Patient"
      } else {
        " with free permutations"
      },
      ": ",
      conditionMessage(result)
    )
  }

  result
}

extract_adonis_stats <- function(adonis_object) {
  result <- as.data.frame(
    adonis_object,
    check.names = FALSE
  )

  term_rows <- setdiff(
    rownames(result),
    c("Residual", "Total")
  )

  if (length(term_rows) == 0L) {
    stop("No model term was found in the adonis2 result.")
  }

  term_row <- term_rows[1]
  p_column <- grep("^Pr\\(", colnames(result), value = TRUE)
  f_column <- grep("^F", colnames(result), value = TRUE)

  list(
    R2 = if ("R2" %in% colnames(result)) {
      as.numeric(result[term_row, "R2"])
    } else {
      NA_real_
    },

    F = if (length(f_column) > 0L) {
      as.numeric(result[term_row, f_column[1]])
    } else {
      NA_real_
    },

    P = if (length(p_column) > 0L) {
      as.numeric(result[term_row, p_column[1]])
    } else {
      NA_real_
    }
  )
}

run_betadisper <- function(distance_object, group) {
  set.seed(1)

  fit <- vegan::betadisper(
    distance_object,
    group,
    type = "median",
    bias.adjust = FALSE
  )

  permutation_result <- vegan::permutest(
    fit,
    permutations = n_permutations
  )

  list(
    fit = fit,
    P = as.numeric(
      permutation_result$tab[1, "Pr(>F)"]
    )
  )
}

get_pcoa <- function(distance_object, metadata) {
  ordination <- vegan::wcmdscale(
    distance_object,
    eig = TRUE,
    k = 2,
    add = "lingoes"
  )

  coordinates <- as.data.frame(
    ordination$points,
    check.names = FALSE
  )

  if (ncol(coordinates) < 2L) {
    stop("The PCoA returned fewer than two coordinate axes.")
  }

  coordinates <- coordinates[, 1:2, drop = FALSE]
  colnames(coordinates) <- c("PCoA1", "PCoA2")
  coordinates$SampleID <- rownames(coordinates)

  coordinates <- coordinates %>%
    left_join(
      metadata %>%
        select(
          SampleID,
          Patient,
          Group,
          Reads
        ),
      by = "SampleID"
    )

  eigenvalues <- as.numeric(ordination$eig)
  positive_eigenvalues <- eigenvalues[eigenvalues > 0]

  if (length(positive_eigenvalues) >= 2L) {
    variance_explained <- c(
      eigenvalues[1] / sum(positive_eigenvalues),
      eigenvalues[2] / sum(positive_eigenvalues)
    )
  } else if (length(positive_eigenvalues) == 1L) {
    variance_explained <- c(1, 0)
  } else {
    variance_explained <- c(NA_real_, NA_real_)
  }

  additive_constant <- if (!is.null(ordination$ac)) {
    as.numeric(ordination$ac)
  } else {
    NA_real_
  }

  list(
    coordinates = coordinates,
    variance_explained = variance_explained,
    eigenvalues = eigenvalues,
    additive_constant = additive_constant
  )
}

# Figure helpers

save_beta_plot <- function(
    plot,
    stem,
    width = 5.7,
    height = 4.4
) {
  pdf_device <- if (capabilities("cairo")) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }

  ggsave(
    file.path(out_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    device = pdf_device,
    bg = "white"
  )

  ggsave(
    file.path(out_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    dpi = export_dpi,
    bg = "white"
  )

  ggsave(
    file.path(out_dir, paste0(stem, ".jpg")),
    plot,
    width = width,
    height = height,
    dpi = export_dpi,
    quality = 100,
    bg = "white"
  )
}

make_beta_plot <- function(
    coordinates,
    variance_explained,
    plot_title,
    statistics_text,
    groups_present
) {
  coordinates$Group <- factor(
    as.character(coordinates$Group),
    levels = groups_present
  )

  group_counts <- table(coordinates$Group)
  ellipse_counts <- group_counts[groups_present]

  draw_ellipse <- (
    length(ellipse_counts) > 0L &&
      all(!is.na(ellipse_counts)) &&
      all(ellipse_counts >= 3L)
  )

  p <- ggplot(
    coordinates,
    aes(
      x = PCoA1,
      y = PCoA2
    )
  )

  if (draw_ellipse) {
    p <- p +
      stat_ellipse(
        aes(fill = Group),
        geom = "polygon",
        alpha = 0.26,
        colour = NA,
        level = ellipse_level,
        show.legend = FALSE
      ) +
      stat_ellipse(
        aes(colour = Group),
        geom = "path",
        linewidth = 0.45,
        level = ellipse_level,
        show.legend = FALSE
      )
  }

  p +
    geom_point(
      aes(colour = Group),
      shape = 16,
      size = 2.2,
      alpha = 0.95
    ) +
    scale_colour_manual(
      values = pal,
      breaks = groups_present,
      limits = groups_present,
      labels = group_display[groups_present],
      drop = FALSE
    ) +
    scale_fill_manual(
      values = pal_fill,
      breaks = groups_present,
      limits = groups_present,
      labels = group_display[groups_present],
      drop = FALSE
    ) +
    annotate(
      "text",
      x = Inf,
      y = -Inf,
      hjust = 1.03,
      vjust = -0.45,
      label = statistics_text,
      size = pt_to_gg(font_pt),
      family = font_family,
      colour = "black",
      lineheight = 0.9
    ) +
    labs(
      title = plot_title,
      x = paste0(
        "PCoA1 (",
        round(100 * variance_explained[1], 1),
        "%)"
      ),
      y = paste0(
        "PCoA2 (",
        round(100 * variance_explained[2], 1),
        "%)"
      )
    ) +
    theme_classic(
      base_family = font_family,
      base_size = font_pt
    ) +
    theme(
      legend.title = element_blank(),
      legend.position = "right",
      legend.key = element_blank(),
      legend.background = element_blank(),

      legend.text = element_text(
        family = font_family,
        size = font_pt,
        colour = "black"
      ),

      legend.key.size = unit(3.2, "mm"),
      legend.spacing.y = unit(0.5, "mm"),

      plot.title = element_text(
        family = font_family,
        face = "bold",
        size = font_pt,
        colour = "black",
        hjust = 0.5
      ),

      plot.tag = element_text(
        family = font_family,
        face = "bold",
        size = 8,
        colour = "black"
      ),

      axis.title = element_text(
        family = font_family,
        size = font_pt,
        colour = "black"
      ),

      axis.text = element_text(
        family = font_family,
        size = font_pt,
        colour = "black"
      ),

      axis.line = element_line(
        linewidth = 0.35,
        colour = "black"
      ),

      axis.ticks = element_line(
        linewidth = 0.35,
        colour = "black"
      ),

      axis.ticks.length = unit(1.2, "mm"),

      plot.background = element_rect(
        fill = "white",
        colour = NA
      ),

      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),

      plot.margin = margin(
        4,
        4,
        4,
        4,
        "mm"
      )
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(
          shape = 16,
          size = 2.6,
          alpha = 1
        )
      ),
      fill = "none"
    )
}

draw_plot_grid <- function(plot_list, nrow, ncol) {
  grid::grid.newpage()

  layout <- grid::grid.layout(
    nrow = nrow,
    ncol = ncol
  )

  grid::pushViewport(
    grid::viewport(layout = layout)
  )

  for (i in seq_along(plot_list)) {
    row_i <- ((i - 1L) %/% ncol) + 1L
    col_i <- ((i - 1L) %% ncol) + 1L

    print(
      plot_list[[i]],
      vp = grid::viewport(
        layout.pos.row = row_i,
        layout.pos.col = col_i
      )
    )
  }

  grid::popViewport()
  invisible(TRUE)
}

save_plot_grid <- function(
    plot_list,
    stem,
    nrow,
    ncol,
    width,
    height
) {
  pdf_file <- file.path(out_dir, paste0(stem, ".pdf"))
  png_file <- file.path(out_dir, paste0(stem, ".png"))
  jpg_file <- file.path(out_dir, paste0(stem, ".jpg"))

  pdf_device <- if (capabilities("cairo")) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }

  pdf_device(
    pdf_file,
    width = width,
    height = height,
    bg = "white"
  )
  draw_plot_grid(plot_list, nrow, ncol)
  grDevices::dev.off()

  grDevices::png(
    png_file,
    width = width,
    height = height,
    units = "in",
    res = export_dpi,
    bg = "white"
  )
  draw_plot_grid(plot_list, nrow, ncol)
  grDevices::dev.off()

  grDevices::jpeg(
    jpg_file,
    width = width,
    height = height,
    units = "in",
    res = export_dpi,
    quality = 100,
    bg = "white"
  )
  draw_plot_grid(plot_list, nrow, ncol)
  grDevices::dev.off()
}

# Load, validate and prepare the final phyloseq object

if (!file.exists(physeq_path)) {
  stop(
    "The final phyloseq object was not found:\n",
    physeq_path
  )
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
  stop(
    "The phyloseq object must contain counts, sample metadata and taxonomy."
  )
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

if (!taxa_are_rows(phy)) {
  otu_check <- t(otu_check)
}

if (
  any(!is.finite(otu_check)) ||
  any(otu_check < 0) ||
  any(abs(otu_check - round(otu_check)) > 1e-8)
) {
  stop("The input must contain non-negative integer ASV counts.")
}

sample_data(phy)$Group <- make_group(
  sample_data(phy)$Sample_Type
)

phy <- subset_samples(
  phy,
  !is.na(Group)
)

phy <- prune_taxa(
  taxa_sums(phy) > 0,
  phy
)

observed_group_counts <- table(
  factor(
    sample_data(phy)$Group,
    levels = group_levels
  )
)

if (!all(
  as.integer(observed_group_counts) ==
    expected_group_counts[group_levels]
)) {
  stop(
    "Unexpected sample counts. Expected ",
    paste(
      names(expected_group_counts),
      expected_group_counts,
      sep = "=",
      collapse = ", "
    ),
    "; observed ",
    paste(
      names(observed_group_counts),
      as.integer(observed_group_counts),
      sep = "=",
      collapse = ", "
    ),
    "."
  )
}

tax_df <- as.data.frame(
  tax_table(phy),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"Genus" %in% colnames(tax_df)) {
  stop("The taxonomy table has no Genus column.")
}

genus_result <- make_lineage_aware_genus(tax_df)
tax_df$Genus <- genus_result$genus

tax_table(phy) <- tax_table(
  as.matrix(tax_df)
)

write.csv(
  genus_result$mapping,
  file.path(
    out_dir,
    "Genus_lineage_mapping.csv"
  ),
  row.names = FALSE
)

phy_genus <- tax_glom(
  phy,
  taxrank = "Genus",
  NArm = FALSE
)

phy_genus <- prune_taxa(
  taxa_sums(phy_genus) > 0,
  phy_genus
)

genus_labels <- as.character(
  tax_table(phy_genus)[, "Genus"]
)

if (
  anyNA(genus_labels) ||
  any(genus_labels == "")
) {
  stop("Missing genus labels remain after genus-level agglomeration.")
}

taxa_names(phy_genus) <- make.unique(genus_labels)

metadata_genus <- get_meta(phy_genus)

sample_data(phy_genus)$Patient <- factor(
  metadata_genus[
    sample_names(phy_genus),
    "Patient"
  ]
)

sample_data(phy_genus)$Group <- factor(
  as.character(
    sample_data(phy_genus)$Group
  ),
  levels = group_levels
)

if (any(as.character(sample_data(phy_genus)$Patient) == "7699")) {
  stop("Patient 7699 is still present according to the resolved Patient IDs.")
}

all_bray_input <- get_bray_curtis(phy_genus)

write.csv(
  data.frame(
    SampleID = rownames(all_bray_input$counts),
    all_bray_input$counts,
    check.names = FALSE
  ),
  file.path(
    out_dir,
    "Genus_count_matrix_all_64_samples.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    SampleID = rownames(all_bray_input$relative_abundance),
    all_bray_input$relative_abundance,
    check.names = FALSE
  ),
  file.path(
    out_dir,
    "Genus_relative_abundance_matrix_all_64_samples.csv"
  ),
  row.names = FALSE
)

sample_manifest <- metadata_genus %>%
  select(
    SampleID,
    Patient,
    Group,
    Reads,
    Patient_source
  )

write.csv(
  sample_manifest,
  file.path(
    out_dir,
    "Sample_manifest_BrayCurtis_Genus.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    Input_object = physeq_path,
    Samples = nsamples(phy_genus),
    Patients = n_distinct(metadata_genus$Patient),
    ASVs_before_genus_agglomeration = ntaxa(phy),
    Genera_after_agglomeration = ntaxa(phy_genus),
    Additional_prevalence_filter_applied = FALSE,
    Distance = "Bray-Curtis",
    Taxonomic_level = "Genus",
    Abundance_scale = "Relative abundance",
    stringsAsFactors = FALSE
  ),
  file.path(
    out_dir,
    "BrayCurtis_Genus_input_summary.csv"
  ),
  row.names = FALSE
)

cat("\n===== FINAL GENUS-LEVEL INPUT =====\n")
cat("Input:", physeq_path, "\n")
cat("Samples:", nsamples(phy_genus), "\n")
cat("Patients:", n_distinct(metadata_genus$Patient), "\n")
cat("ASVs retained in 5% object:", ntaxa(phy), "\n")
cat("Genus-level categories:", ntaxa(phy_genus), "\n")
cat("No additional prevalence filtering was applied.\n")

# Run one Bray-Curtis analysis

run_beta_analysis <- function(
    phy,
    groups,
    comparison_label,
    use_patient_strata = FALSE
) {
  groups_present <- groups[groups %in% group_levels]

  metadata_full <- get_meta(phy)

  sample_ids <- metadata_full$SampleID[
    !is.na(metadata_full$Group) &
      as.character(metadata_full$Group) %in% groups_present
  ]

  phy_subset <- prune_samples(
    sample_ids,
    phy
  )

  phy_subset <- prune_taxa(
    taxa_sums(phy_subset) > 0,
    phy_subset
  )

  metadata <- get_meta(phy_subset)

  metadata$Group <- factor(
    as.character(metadata$Group),
    levels = groups_present
  )

  metadata <- metadata[
    !is.na(metadata$Group),
    ,
    drop = FALSE
  ]

  rownames(metadata) <- metadata$SampleID

  group_counts <- table(metadata$Group)

  if (
    nrow(metadata) == 0L ||
    length(group_counts) < 2L ||
    any(group_counts == 0L)
  ) {
    stop(
      "Insufficient data for ",
      comparison_label,
      "."
    )
  }

  bray_object <- get_bray_curtis(phy_subset)

  metadata <- metadata[
    rownames(bray_object$relative_abundance),
    ,
    drop = FALSE
  ]

  metadata$Group <- factor(
    as.character(metadata$Group),
    levels = groups_present
  )

  if (!identical(
    labels(bray_object$distance),
    rownames(metadata)
  )) {
    stop(
      "The Bray-Curtis distance labels and metadata rows are not aligned for ",
      comparison_label,
      "."
    )
  }

  patient_group_table <- table(
    metadata$Patient,
    metadata$Group
  )

  paired_patients <- rownames(patient_group_table)[
    rowSums(patient_group_table > 0) ==
      length(groups_present)
  ]

  n_paired_patients <- length(paired_patients)

  if (
    use_patient_strata &&
    n_paired_patients < 2L
  ) {
    stop(
      comparison_label,
      " was configured for Patient-restricted permutations, ",
      "but fewer than two patients contributed both tissue groups."
    )
  }

  permanova <- safe_adonis2(
    distance_object = bray_object$distance,
    metadata = metadata,
    use_patient_strata = use_patient_strata
  )

  permanova_stats <- extract_adonis_stats(
    permanova
  )

  dispersion <- run_betadisper(
    bray_object$distance,
    metadata$Group
  )

  pcoa <- get_pcoa(
    bray_object$distance,
    metadata
  )

  summary_row <- data.frame(
    Comparison = comparison_label,
    Taxonomic_level = "Genus",
    Abundance_scale = "Relative abundance",
    Distance = "Bray-Curtis",
    PCoA_correction = "Lingoes",
    PCoA_additive_constant = pcoa$additive_constant,
    N_samples = nrow(metadata),
    N_patients = n_distinct(metadata$Patient),
    N_genera_in_distance = ncol(bray_object$relative_abundance),
    Group_1 = groups_present[1],
    Group_1_n = sum(metadata$Group == groups_present[1]),
    Group_2 = if (length(groups_present) >= 2L) {
      groups_present[2]
    } else {
      NA_character_
    },
    Group_2_n = if (length(groups_present) >= 2L) {
      sum(metadata$Group == groups_present[2])
    } else {
      NA_integer_
    },
    Group_3 = if (length(groups_present) >= 3L) {
      groups_present[3]
    } else {
      NA_character_
    },
    Group_3_n = if (length(groups_present) >= 3L) {
      sum(metadata$Group == groups_present[3])
    } else {
      NA_integer_
    },
    Permutation_scheme = if (use_patient_strata) {
      "Permutations restricted within Patient"
    } else {
      "Free permutations"
    },
    N_permutations = n_permutations,
    N_paired_patients = n_paired_patients,
    R2_PERMANOVA = permanova_stats$R2,
    F_PERMANOVA = permanova_stats$F,
    P_PERMANOVA_raw = permanova_stats$P,
    P_betadisper = dispersion$P,
    PCoA1_percent = 100 * pcoa$variance_explained[1],
    PCoA2_percent = 100 * pcoa$variance_explained[2],
    stringsAsFactors = FALSE
  )

  list(
    summary = summary_row,
    coordinates = pcoa$coordinates,
    variance_explained = pcoa$variance_explained,
    eigenvalues = pcoa$eigenvalues,
    additive_constant = pcoa$additive_constant,
    groups_present = groups_present,
    metadata = metadata,
    count_matrix = bray_object$counts,
    relative_abundance = bray_object$relative_abundance,
    distance = bray_object$distance,
    permanova = permanova,
    betadisper = dispersion$fit
  )
}

# Overall three-group analysis

overall_result <- run_beta_analysis(
  phy = phy_genus,
  groups = group_levels,
  comparison_label = "Overall_3groups",
  use_patient_strata = FALSE
)

overall_summary <- overall_result$summary

overall_statistics_text <- paste0(
  "R² = ",
  sprintf("%.3f", overall_summary$R2_PERMANOVA),
  ", ",
  format_p(
    overall_summary$P_PERMANOVA_raw,
    "PERMANOVA P"
  ),
  "\n",
  format_p(
    overall_summary$P_betadisper,
    "betadisper P"
  )
)

overall_plot <- make_beta_plot(
  coordinates = overall_result$coordinates,
  variance_explained = overall_result$variance_explained,
  plot_title = "All groups",
  statistics_text = overall_statistics_text,
  groups_present = group_levels
)

save_beta_plot(
  overall_plot,
  "PCoA_BrayCurtis_Genus_overall_3groups"
)

write.csv(
  overall_summary,
  file.path(
    out_dir,
    "Overall_BrayCurtis_Genus_PERMANOVA_betadisper_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  overall_result$coordinates,
  file.path(
    out_dir,
    "PCoA_BrayCurtis_Genus_overall_3groups_coordinates.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    Axis = seq_along(overall_result$eigenvalues),
    Eigenvalue = overall_result$eigenvalues
  ),
  file.path(
    out_dir,
    "PCoA_BrayCurtis_Genus_overall_3groups_eigenvalues.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    SampleID = rownames(overall_result$relative_abundance),
    overall_result$relative_abundance,
    check.names = FALSE
  ),
  file.path(
    out_dir,
    "Overall_BrayCurtis_Genus_relative_abundance_input.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.matrix(overall_result$distance),
  file.path(
    out_dir,
    "Overall_BrayCurtis_Genus_distance_matrix.csv"
  ),
  row.names = TRUE
)

capture.output(
  overall_result$permanova,
  file = file.path(
    out_dir,
    "Overall_BrayCurtis_Genus_PERMANOVA_full_output.txt"
  )
)

# Three pairwise analyses

pairwise_specifications <- list(
  list(
    groups = c("H_ADJ_SSL_HP", "H_COL"),
    label = "H_ADJ_SSL_HP vs H_COL",
    use_patient_strata = FALSE
  ),
  list(
    groups = c("SSL_HP", "H_COL"),
    label = "SSL_HP vs H_COL",
    use_patient_strata = FALSE
  ),
  list(
    groups = c("SSL_HP", "H_ADJ_SSL_HP"),
    label = "SSL_HP vs H_ADJ_SSL_HP",
    use_patient_strata = TRUE
  )
)

pairwise_results <- lapply(
  pairwise_specifications,
  function(specification) {
    run_beta_analysis(
      phy = phy_genus,
      groups = specification$groups,
      comparison_label = specification$label,
      use_patient_strata = specification$use_patient_strata
    )
  }
)

names(pairwise_results) <- vapply(
  pairwise_specifications,
  function(x) x$label,
  character(1)
)

pairwise_summary <- bind_rows(
  lapply(
    pairwise_results,
    `[[`,
    "summary"
  )
)

pairwise_summary$P_PERMANOVA_Holm <- p.adjust(
  pairwise_summary$P_PERMANOVA_raw,
  method = "holm"
)

pairwise_summary$Significant_after_Holm <- (
  !is.na(pairwise_summary$P_PERMANOVA_Holm) &
    pairwise_summary$P_PERMANOVA_Holm < 0.05
)

pairwise_summary <- pairwise_summary %>%
  mutate(
    Comparison = factor(
      Comparison,
      levels = pairwise_plot_order
    )
  ) %>%
  arrange(Comparison) %>%
  mutate(
    Comparison = as.character(Comparison)
  )

write.csv(
  pairwise_summary,
  file.path(
    out_dir,
    "Pairwise_BrayCurtis_Genus_PERMANOVA_betadisper_raw_and_Holm.csv"
  ),
  row.names = FALSE
)

pairwise_plots <- list()

for (comparison_label in pairwise_plot_order) {
  result <- pairwise_results[[comparison_label]]

  result_row <- pairwise_summary[
    pairwise_summary$Comparison == comparison_label,
    ,
    drop = FALSE
  ]

  statistics_text <- paste0(
    "R² = ",
    sprintf("%.3f", result_row$R2_PERMANOVA),
    "\n",
    format_p(
      result_row$P_PERMANOVA_raw,
      "PERMANOVA raw P"
    ),
    "\n",
    format_p(
      result_row$P_PERMANOVA_Holm,
      "Holm P"
    ),
    "\n",
    format_p(
      result_row$P_betadisper,
      "betadisper P"
    )
  )

  plot_title <- paste(
    group_display[result$groups_present[1]],
    "vs",
    group_display[result$groups_present[2]]
  )

  p <- make_beta_plot(
    coordinates = result$coordinates,
    variance_explained = result$variance_explained,
    plot_title = plot_title,
    statistics_text = statistics_text,
    groups_present = result$groups_present
  )

  pairwise_plots[[comparison_label]] <- p

  stem <- paste0(
    "PCoA_BrayCurtis_Genus_",
    safe_name(comparison_label)
  )

  save_beta_plot(
    p,
    stem
  )

  write.csv(
    result$coordinates,
    file.path(
      out_dir,
      paste0(stem, "_coordinates.csv")
    ),
    row.names = FALSE
  )

  write.csv(
    data.frame(
      Axis = seq_along(result$eigenvalues),
      Eigenvalue = result$eigenvalues
    ),
    file.path(
      out_dir,
      paste0(stem, "_eigenvalues.csv")
    ),
    row.names = FALSE
  )

  write.csv(
    data.frame(
      SampleID = rownames(result$relative_abundance),
      result$relative_abundance,
      check.names = FALSE
    ),
    file.path(
      out_dir,
      paste0(stem, "_relative_abundance_input.csv")
    ),
    row.names = FALSE
  )

  write.csv(
    as.matrix(result$distance),
    file.path(
      out_dir,
      paste0(stem, "_distance_matrix.csv")
    ),
    row.names = TRUE
  )

  capture.output(
    result$permanova,
    file = file.path(
      out_dir,
      paste0(stem, "_PERMANOVA_full_output.txt")
    )
  )
}

# Combined three-panel pairwise figure

pairwise_panel_plots <- lapply(
  seq_along(pairwise_plot_order),
  function(i) {
    pairwise_plots[[pairwise_plot_order[i]]] +
      labs(tag = LETTERS[i]) +
      theme(
        plot.tag.position = c(0.01, 0.99)
      )
  }
)

save_plot_grid(
  plot_list = pairwise_panel_plots,
  stem = "Figure_BrayCurtis_Genus_pairwise_3panel",
  nrow = 1,
  ncol = 3,
  width = 17.1,
  height = 4.6
)

overall_plus_paired <- list(
  overall_plot +
    labs(tag = "A") +
    theme(plot.tag.position = c(0.01, 0.99)),

  pairwise_plots[["SSL_HP vs H_ADJ_SSL_HP"]] +
    labs(tag = "B") +
    theme(plot.tag.position = c(0.01, 0.99))
)

save_plot_grid(
  plot_list = overall_plus_paired,
  stem = "PCoA_BrayCurtis_Genus_main_overall_plus_SSLHP_vs_HADJSSLHP",
  nrow = 1,
  ncol = 2,
  width = 11.4,
  height = 4.6
)

# Combined result table

all_beta_results <- bind_rows(
  overall_summary %>%
    mutate(
      P_PERMANOVA_Holm = NA_real_,
      Significant_after_Holm = NA
    ),
  pairwise_summary
)

write.csv(
  all_beta_results,
  file.path(
    out_dir,
    "BrayCurtis_Genus_ALL_PERMANOVA_betadisper_results.csv"
  ),
  row.names = FALSE
)

# Reproducibility notes

writeLines(
  c(
    "FINAL GENUS-LEVEL BRAY-CURTIS BETA-DIVERSITY ANALYSIS",
    "======================================================",
    "",
    paste0("Input object: ", physeq_path),
    paste0("Output directory: ", out_dir),
    paste0("Samples: ", nsamples(phy_genus)),
    paste0("Patients: ", n_distinct(metadata_genus$Patient)),
    paste0("ASVs in existing 5%-filtered object: ", ntaxa(phy)),
    paste0("Lineage-aware genus categories: ", ntaxa(phy_genus)),
    "Patient 7699 is excluded by the selected input object.",
    "No additional prevalence filtering is applied.",
    "",
    "Distance:",
    "  Bray-Curtis dissimilarity calculated from sample-level genus relative abundance.",
    "  No pseudocount is used.",
    "  binary = FALSE; therefore abundance, not presence/absence alone, is used.",
    "",
    "Ordination:",
    "  PCoA generated with vegan::wcmdscale.",
    "  Lingoes correction is applied only to make the Bray-Curtis ordination Euclidean.",
    "  PERMANOVA and betadisper use the original uncorrected Bray-Curtis distances.",
    "",
    "Permutation schemes:",
    "  Overall three-group analysis: free permutations.",
    "  SSL_HP versus H_COL: free permutations.",
    "  H_ADJ_SSL_HP versus H_COL: free permutations.",
    "  SSL_HP versus H_ADJ_SSL_HP: permutations restricted within Patient.",
    paste0("  Number of permutations: ", n_permutations, "."),
    "",
    "Multiple testing:",
    "  Overall three-group PERMANOVA P value is not adjusted.",
    "  The three pairwise PERMANOVA P values are Holm-adjusted as one family.",
    "  betadisper P values are diagnostic and are not Holm-adjusted.",
    "",
    "Interpretation:",
    "  This is a genus-level Bray-Curtis sensitivity analysis.",
    "  The primary compositional beta-diversity analysis remains CLR/Aitchison at ASV level."
  ),
  file.path(
    out_dir,
    "README_BrayCurtis_Genus_FINAL.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "sessionInfo_BrayCurtis_Genus_FINAL.txt"
  )
)

# Print key results

cat("\n===== OVERALL BRAY-CURTIS SUMMARY =====\n")
print(overall_summary)

cat("\n===== PAIRWISE BRAY-CURTIS SUMMARY =====\n")
print(pairwise_summary)

cat("\nResults saved to:\n")
cat(out_dir, "\n")
cat("\nGENUS-LEVEL BRAY-CURTIS ANALYSIS COMPLETE\n")
