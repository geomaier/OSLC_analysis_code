#!/usr/bin/env Rscript
# ASV-level Aitchison beta diversity. Run from the 16S directory.

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

physeq_path <- file.path("outputs", "phyloseq_final_March26_5pct_without7699.rds")

out_dir <- file.path("outputs", "Beta_diversity_CLR_Aitchison_FINAL",
                     format(Sys.time(), "%Y%m%d_%H%M%S"))
if (dir.exists(out_dir)) stop("Output directory already exists: ", out_dir)
if (!dir.create(out_dir, recursive = TRUE)) stop("Could not create: ", out_dir)

patient_col_manual <- NULL

n_permutations <- 9999
clr_pseudocount <- 0.5
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
  "H_ADJ_SSL_HP"= "H_ADJ_SSL_HP",
  "SSL_HP" = "SSL_HP"
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

# Helpers

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

first_existing_column <- function(column_names, candidates) {
  hit <- candidates[candidates %in% column_names]
  if (length(hit) == 0L) return(NULL)
  hit[1]
}

format_p <- function(p, prefix = "p") {
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
    "Subject_ID"
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
  } else {
    patient <- sub(
      ".*_([0-9]+)$",
      "\\1",
      meta$SampleID
    )

    patient[patient == meta$SampleID] <- NA_character_
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

  rownames(meta) <- meta$SampleID
  meta
}

get_clr_dist <- function(phy) {
  otu <- as(otu_table(phy), "matrix")

  if (taxa_are_rows(phy)) {
    otu <- t(otu)
  }

  otu <- otu[, colSums(otu) > 0, drop = FALSE]

  if (any(otu < 0) || any(!is.finite(otu))) {
    stop("The ASV table contains negative or non-finite values.")
  }

  if (any(rowSums(otu) <= 0)) {
    stop("At least one sample has zero total abundance.")
  }

  composition <- otu + clr_pseudocount
  composition <- composition / rowSums(composition)

  clr <- log(composition) - rowMeans(log(composition))

  list(
    clr = clr,
    distance = stats::dist(clr, method = "euclidean")
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
  ordination <- stats::cmdscale(
    distance_object,
    eig = TRUE,
    k = 2
  )

  coordinates <- as.data.frame(
    ordination$points,
    check.names = FALSE
  )

  colnames(coordinates) <- c(
    "PCoA1",
    "PCoA2"
  )

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

  positive_eigenvalues <- ordination$eig[
    ordination$eig > 0
  ]

  if (length(positive_eigenvalues) >= 2L) {
    variance_explained <- positive_eigenvalues /
      sum(positive_eigenvalues)
  } else if (length(positive_eigenvalues) == 1L) {
    variance_explained <- c(1, 0)
  } else {
    variance_explained <- c(NA_real_, NA_real_)
  }

  list(
    coordinates = coordinates,
    variance_explained = variance_explained
  )
}

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
        round(
          100 * variance_explained[1],
          1
        ),
        "%)"
      ),
      y = paste0(
        "PCoA2 (",
        round(
          100 * variance_explained[2],
          1
        ),
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

run_beta_analysis <- function(
    phy,
    groups,
    comparison_label,
    use_patient_strata = FALSE
) {
  groups_present <- group_levels[
    group_levels %in% groups
  ]

  metadata_full <- get_meta(phy)

  sample_ids <- metadata_full$SampleID[
    !is.na(metadata_full$Group) &
      as.character(metadata_full$Group) %in%
      groups_present
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

  clr_object <- get_clr_dist(phy_subset)

  metadata <- metadata[
    rownames(clr_object$clr),
    ,
    drop = FALSE
  ]

  metadata$Group <- factor(
    as.character(metadata$Group),
    levels = groups_present
  )

  patient_group_table <- table(
    metadata$Patient,
    metadata$Group
  )

  paired_patients <- rownames(
    patient_group_table
  )[
    rowSums(patient_group_table > 0) ==
      length(groups_present)
  ]

  n_paired_patients <- length(
    paired_patients
  )

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
    distance_object = clr_object$distance,
    metadata = metadata,
    use_patient_strata = use_patient_strata
  )

  permanova_stats <- extract_adonis_stats(
    permanova
  )

  dispersion <- run_betadisper(
    clr_object$distance,
    metadata$Group
  )

  pcoa <- get_pcoa(
    clr_object$distance,
    metadata
  )

  summary_row <- data.frame(
    Comparison = comparison_label,
    N_samples = nrow(metadata),
    N_patients = n_distinct(metadata$Patient),
    Group_1 = groups_present[1],
    Group_1_n = sum(
      metadata$Group == groups_present[1]
    ),
    Group_2 = if (
      length(groups_present) >= 2L
    ) {
      groups_present[2]
    } else {
      NA_character_
    },
    Group_2_n = if (
      length(groups_present) >= 2L
    ) {
      sum(metadata$Group == groups_present[2])
    } else {
      NA_integer_
    },
    Group_3 = if (
      length(groups_present) >= 3L
    ) {
      groups_present[3]
    } else {
      NA_character_
    },
    Group_3_n = if (
      length(groups_present) >= 3L
    ) {
      sum(metadata$Group == groups_present[3])
    } else {
      NA_integer_
    },
    Permutation_scheme = if (
      use_patient_strata
    ) {
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
    stringsAsFactors = FALSE
  )

  list(
    summary = summary_row,
    coordinates = pcoa$coordinates,
    variance_explained = pcoa$variance_explained,
    groups_present = groups_present,
    metadata = metadata,
    distance = clr_object$distance
  )
}

# Load phyloseq and build the group variable

if (!file.exists(physeq_path)) {
  stop(
    "The phyloseq object was not found: ",
    physeq_path
  )
}

physeq <- readRDS(physeq_path)

if (!inherits(physeq, "phyloseq")) {
  stop("The loaded object is not a phyloseq object.")
}

if (!"Sample_Type" %in% colnames(sample_data(physeq))) {
  stop("Sample_Type was not found in sample_data(physeq).")
}

if (any(grepl("7699", sample_names(physeq), fixed = TRUE))) {
  stop(
    "The selected object still contains patient 7699. ",
    "Use phyloseq_final_March26_5pct_without7699.rds."
  )
}

group <- clean_character(
  sample_data(physeq)$Sample_Type
)

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

sample_data(physeq)$Group <- factor(
  group,
  levels = group_levels
)

physeq <- subset_samples(
  physeq,
  !is.na(Group)
)

physeq <- prune_taxa(
  taxa_sums(physeq) > 0,
  physeq
)

cat("Loaded phyloseq object:\n")
cat("Samples:", nsamples(physeq), "\n")
cat("ASVs:", ntaxa(physeq), "\n")

cat("\nSamples by group:\n")
print(
  table(
    sample_data(physeq)$Group,
    useNA = "ifany"
  )
)

# Overall three-group analysis

overall_result <- run_beta_analysis(
  phy = physeq,
  groups = group_levels,
  comparison_label = "Overall_3groups",
  use_patient_strata = FALSE
)

overall_summary <- overall_result$summary

overall_statistics_text <- paste0(
  "R² = ",
  sprintf(
    "%.3f",
    overall_summary$R2_PERMANOVA
  ),
  ", ",
  format_p(
    overall_summary$P_PERMANOVA_raw,
    "PERMANOVA p"
  ),
  "\n",
  format_p(
    overall_summary$P_betadisper,
    "betadisper p"
  )
)

overall_plot <- make_beta_plot(
  coordinates = overall_result$coordinates,
  variance_explained =
    overall_result$variance_explained,
  plot_title =
    "All groups",
  statistics_text = overall_statistics_text,
  groups_present = group_levels
)

save_beta_plot(
  overall_plot,
  "PCoA_CLR_overall_3groups"
)

write.csv(
  overall_summary,
  file.path(
    out_dir,
    "Overall_PERMANOVA_betadisper_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  overall_result$coordinates,
  file.path(
    out_dir,
    "PCoA_CLR_overall_3groups_coordinates.csv"
  ),
  row.names = FALSE
)

# Three pairwise analyses

pairwise_specifications <- list(
  list(
    groups = c("SSL_HP", "H_COL"),
    label = "SSL_HP vs H_COL",
    use_patient_strata = FALSE
  ),
  list(
    groups = c(
      "H_ADJ_SSL_HP",
      "H_COL"
    ),
    label = "H_ADJ_SSL_HP vs H_COL",
    use_patient_strata = FALSE
  ),
  list(
    groups = c(
      "SSL_HP",
      "H_ADJ_SSL_HP"
    ),
    label = "SSL_HP vs H_ADJ_SSL_HP",
    use_patient_strata = TRUE
  )
)

pairwise_results <- lapply(
  pairwise_specifications,
  function(specification) {
    run_beta_analysis(
      phy = physeq,
      groups = specification$groups,
      comparison_label = specification$label,
      use_patient_strata =
        specification$use_patient_strata
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
  pairwise_summary$P_PERMANOVA_Holm < 0.05
)

write.csv(
  pairwise_summary,
  file.path(
    out_dir,
    "Pairwise_PERMANOVA_betadisper_raw_and_Holm.csv"
  ),
  row.names = FALSE
)

for (
  comparison_label in names(
    pairwise_results
  )
) {
  result <- pairwise_results[[
    comparison_label
  ]]

  result_row <- pairwise_summary[
    pairwise_summary$Comparison ==
      comparison_label,
    ,
    drop = FALSE
  ]

  statistics_text <- paste0(
    "R² = ",
    sprintf(
      "%.3f",
      result_row$R2_PERMANOVA
    ),
    "\n",
    format_p(
      result_row$P_PERMANOVA_raw,
      "PERMANOVA raw p"
    ),
    "\n",
    format_p(
      result_row$P_PERMANOVA_Holm,
      "Holm p"
    ),
    "\n",
    format_p(
      result_row$P_betadisper,
      "betadisper p"
    )
  )

  plot_title <- paste(
    group_display[
      result$groups_present[1]
    ],
    "vs",
    group_display[
      result$groups_present[2]
    ]
  )

  p <- make_beta_plot(
    coordinates = result$coordinates,
    variance_explained =
      result$variance_explained,
    plot_title = plot_title,
    statistics_text = statistics_text,
    groups_present = result$groups_present
  )

  stem <- paste0(
    "PCoA_CLR_",
    gsub(
      "[^A-Za-z0-9]+",
      "_",
      comparison_label
    )
  )

  save_beta_plot(
    p,
    stem
  )

  write.csv(
    result$coordinates,
    file.path(
      out_dir,
      paste0(
        stem,
        "_coordinates.csv"
      )
    ),
    row.names = FALSE
  )
}

# Compact main figure

overall_coords <- overall_result$coordinates
overall_coords$Panel <- paste0(
  "Overall\nPCoA1 ",
  round(
    100 * overall_result$variance_explained[1],
    1
  ),
  "%; PCoA2 ",
  round(
    100 * overall_result$variance_explained[2],
    1
  ),
  "%"
)

ssl_adj_result <- pairwise_results[[
  "SSL_HP vs H_ADJ_SSL_HP"
]]

ssl_adj_summary <- pairwise_summary[
  pairwise_summary$Comparison ==
    "SSL_HP vs H_ADJ_SSL_HP",
  ,
  drop = FALSE
]

ssl_adj_coords <- ssl_adj_result$coordinates
ssl_adj_coords$Panel <- paste0(
  "SSL/HP vs adjacent\nPCoA1 ",
  round(
    100 * ssl_adj_result$variance_explained[1],
    1
  ),
  "%; PCoA2 ",
  round(
    100 * ssl_adj_result$variance_explained[2],
    1
  ),
  "%"
)

compact_coordinates <- bind_rows(
  overall_coords,
  ssl_adj_coords
)

panel_levels <- c(
  unique(overall_coords$Panel),
  unique(ssl_adj_coords$Panel)
)

compact_coordinates$Panel <- factor(
  compact_coordinates$Panel,
  levels = panel_levels
)

compact_statistics <- data.frame(
  Panel = factor(
    panel_levels,
    levels = panel_levels
  ),
  x = Inf,
  y = -Inf,
  label = c(
    paste0(
      "R² = ",
      sprintf(
        "%.3f",
        overall_summary$R2_PERMANOVA
      ),
      "\n",
      format_p(
        overall_summary$P_PERMANOVA_raw,
        "PERMANOVA p"
      ),
      "\n",
      format_p(
        overall_summary$P_betadisper,
        "betadisper p"
      )
    ),
    paste0(
      "R² = ",
      sprintf(
        "%.3f",
        ssl_adj_summary$R2_PERMANOVA
      ),
      "\n",
      format_p(
        ssl_adj_summary$P_PERMANOVA_raw,
        "PERMANOVA raw p"
      ),
      "\n",
      format_p(
        ssl_adj_summary$P_PERMANOVA_Holm,
        "Holm p"
      ),
      "\n",
      format_p(
        ssl_adj_summary$P_betadisper,
        "betadisper p"
      )
    )
  ),
  stringsAsFactors = FALSE
)

ellipse_data <- compact_coordinates %>%
  group_by(
    Panel,
    Group
  ) %>%
  mutate(
    N_group_panel = n()
  ) %>%
  ungroup() %>%
  group_by(Panel) %>%
  mutate(
    Draw_ellipse = all(
      N_group_panel >= 3
    )
  ) %>%
  ungroup() %>%
  filter(Draw_ellipse)

compact_plot <- ggplot(
  compact_coordinates,
  aes(
    x = PCoA1,
    y = PCoA2
  )
) +
  stat_ellipse(
    data = ellipse_data,
    aes(fill = Group),
    geom = "polygon",
    alpha = 0.26,
    colour = NA,
    level = ellipse_level,
    show.legend = FALSE
  ) +
  stat_ellipse(
    data = ellipse_data,
    aes(colour = Group),
    geom = "path",
    linewidth = 0.45,
    level = ellipse_level,
    show.legend = FALSE
  ) +
  geom_point(
    aes(colour = Group),
    shape = 16,
    size = 1.7,
    alpha = 0.95
  ) +
  geom_text(
    data = compact_statistics,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1.03,
    vjust = -0.45,
    family = font_family,
    size = pt_to_gg(font_pt),
    colour = "black",
    lineheight = 0.9
  ) +
  facet_wrap(
    ~ Panel,
    nrow = 1,
    scales = "free"
  ) +
  scale_colour_manual(
    values = pal,
    breaks = group_levels,
    limits = group_levels,
    labels = group_display[group_levels],
    drop = FALSE
  ) +
  scale_fill_manual(
    values = pal_fill,
    breaks = group_levels,
    limits = group_levels,
    labels = group_display[group_levels],
    drop = FALSE
  ) +
  labs(
    x = "PCoA1",
    y = "PCoA2"
  ) +
  theme_classic(
    base_family = font_family,
    base_size = font_pt
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.key = element_blank(),
    legend.background = element_blank(),

    legend.text = element_text(
      family = font_family,
      size = font_pt,
      colour = "black"
    ),

    legend.key.size = unit(3.2, "mm"),
    legend.spacing.x = unit(1.5, "mm"),

    strip.background = element_blank(),

    strip.text = element_text(
      family = font_family,
      size = font_pt,
      face = "bold",
      colour = "black",
      lineheight = 0.9
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

    axis.ticks.length = unit(
      1.2,
      "mm"
    ),

    panel.spacing.x = unit(
      4,
      "mm"
    ),

    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),

    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),

    plot.margin = margin(
      3,
      3,
      3,
      3,
      "mm"
    )
  ) +
  guides(
    colour = guide_legend(
      nrow = 1,
      override.aes = list(
        shape = 16,
        size = 2.4,
        alpha = 1
      )
    ),
    fill = "none"
  )

pdf_device <- if (capabilities("cairo")) {
  grDevices::cairo_pdf
} else {
  grDevices::pdf
}

ggsave(
  file.path(
    out_dir,
    "PCoA_CLR_main_overall_plus_SSLHP_vs_HADJSSLHP.pdf"
  ),
  compact_plot,
  width = 125,
  height = 66,
  units = "mm",
  device = pdf_device,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "PCoA_CLR_main_overall_plus_SSLHP_vs_HADJSSLHP.png"
  ),
  compact_plot,
  width = 125,
  height = 66,
  units = "mm",
  dpi = export_dpi,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "PCoA_CLR_main_overall_plus_SSLHP_vs_HADJSSLHP.jpg"
  ),
  compact_plot,
  width = 125,
  height = 66,
  units = "mm",
  dpi = export_dpi,
  quality = 100,
  bg = "white"
)

# Reproducibility notes

writeLines(
  c(
    "FINAL CLR/AITCHISON BETA-DIVERSITY ANALYSIS",
    "============================================",
    "",
    paste0("Input object: ", physeq_path),
    paste0("CLR pseudocount: ", clr_pseudocount),
    paste0("Permutations: ", n_permutations),
    paste0("Ellipse level: ", ellipse_level),
    "",
    "Distance:",
    "  Euclidean distance calculated between CLR-transformed ASV profiles.",
    "  This is the Aitchison distance.",
    "",
    "Permutation schemes:",
    "  Overall three-group analysis: free permutations.",
    "  SSL_HP versus H_COL: free permutations.",
    "  H_ADJ_SSL_HP versus H_COL: free permutations.",
    "  SSL_HP versus H_ADJ_SSL_HP: permutations restricted within Patient.",
    "",
    "Multiple testing:",
    "  Overall three-group PERMANOVA P value is not adjusted.",
    "  The three pairwise PERMANOVA P values are Holm-adjusted as one family.",
    "  Both raw and Holm-adjusted pairwise P values are exported and plotted.",
    "",
    "Dispersion:",
    "  betadisper uses the spatial median; permutest requests 9,999 free permutations.",
    "  betadisper P values are reported as diagnostic tests and are not Holm-adjusted."
  ),
  file.path(
    out_dir,
    "README_Beta_diversity_FINAL.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "sessionInfo_Beta_diversity_FINAL.txt"
  )
)

cat("\n===== OVERALL SUMMARY =====\n")
print(overall_summary)

cat("\n===== PAIRWISE SUMMARY =====\n")
print(pairwise_summary)

cat("\nResults saved to:\n", out_dir, "\n")
cat("\nBETA DIVERSITY COMPLETE\n")
