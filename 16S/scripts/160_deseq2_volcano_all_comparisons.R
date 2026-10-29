#!/usr/bin/env Rscript
# Genus DESeq2: patient-tissue summed counts and all available profiles.
rm(list = ls())
set.seed(20260820)

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(phyloseq)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
})

if (capabilities("cairo")) {
  options(bitmapType = "cairo")
}

config <- new.env(parent = baseenv())
if (file.exists("config/config.R")) sys.source("config/config.R", envir = config)
work_dir <- if (exists("PATH_WORK", config, inherits = FALSE)) config$PATH_WORK else "outputs"
if (!dir.exists(work_dir)) stop("Input directory does not exist: ", work_dir)

physeq_path <- file.path(work_dir, "phyloseq_final_March26_5pct_without7699.rds")

out_dir <- file.path(work_dir, "Section4_5_2_DESeq2_GENUS_ALL_PROFILES_THIRD_FINAL")
out_dir <- file.path(out_dir, paste0("corrected_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
if (dir.exists(out_dir)) stop("Output folder already exists: ", out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ancom_root <- file.path(work_dir, "Section4_5_3_ANCOMBC2_GENUS_ALL_PROFILES_THIRD_FINAL")
make_ancom_concordance <- FALSE

group_levels <- c(
  "H_COL",
  "H_ADJ_SSL_HP",
  "SSL_HP"
)

expected_group_counts <- c(
  "H_COL" = 14L,
  "H_ADJ_SSL_HP" = 12L,
  "SSL_HP" = 38L
)

alpha <- 0.05
export_dpi <- 600

font_family <- "Arial"
font_pt <- 7
legend_pt <- 6
label_pt <- 7

pt_to_gg <- function(x) {
  x / 2.845276
}

volcano_colours <- c(
  "Up" = "#A52A2A",
  "Down" = "green4",
  "NS" = "grey75"
)

volcano_width_mm <- 160
volcano_height_mm <- 110

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "na", "N/A")] <- NA_character_
  x
}

make_group <- function(x) {
  group <- clean_character(x)

  group[group %in% c("SSL", "HP", "SSL_HP")] <- "SSL_HP"

  group[
    group %in% c(
      "H_ADJ_SSL",
      "H_ADJ_HP",
      "H_ADJ_SSL_HP"
    )
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
  ] <- paste0(
    "Unclassified_",
    gsub("[^[:alnum:]_]+", "_", parent[missing_genus & !is.na(parent)])
  )

  genus[
    missing_genus & is.na(parent)
  ] <- paste0(
    "Unclassified_",
    gsub("[^[:alnum:]_]+", "_", row_ids[missing_genus & is.na(parent)])
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

  key_map <- paste(map$Original_genus, map$Lineage_key, sep = "\r")
  key_tax <- paste(genus, lineage_key, sep = "\r")

  final_genus <- map$Analysis_genus[match(key_tax, key_map)]

  if (anyNA(final_genus)) {
    stop("At least one ASV could not be assigned a lineage-aware genus label.")
  }

  list(genus = final_genus, mapping = map)
}

first_existing_column <- function(
    column_names,
    candidates
) {
  hit <- candidates[
    candidates %in% column_names
  ]

  if (length(hit) == 0L) {
    return(NULL)
  }

  hit[1]
}

save_ggplot <- function(
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
    file.path(
      out_dir,
      paste0(stem, ".pdf")
    ),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = pdf_device,
    bg = "white"
  )

  ggsave(
    file.path(
      out_dir,
      paste0(stem, ".png")
    ),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = export_dpi,
    bg = "white"
  )

  ggsave(
    file.path(
      out_dir,
      paste0(stem, ".jpg")
    ),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = export_dpi,
    quality = 100,
    bg = "white"
  )
}

safe_name <- function(x) {
  x <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    x
  )

  gsub(
    "_+$|^_+",
    "",
    x
  )
}

if (!file.exists(physeq_path)) {
  stop(
    "The final phyloseq object was not found:\n",
    physeq_path
  )
}

# Input validation and lineage-aware genus aggregation.
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
    "The phyloseq object must contain counts, metadata and taxonomy."
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
  stop("DESeq2 requires non-negative integer raw counts.")
}

sample_data(phy)$Group <- make_group(
  sample_data(phy)$Sample_Type
)

phy <- subset_samples(phy, !is.na(Group))
phy <- prune_taxa(taxa_sums(phy) > 0, phy)

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

tax_table(phy) <- tax_table(as.matrix(tax_df))

write.csv(
  genus_result$mapping,
  file.path(out_dir, "Genus_lineage_mapping.csv"),
  row.names = FALSE
)

metadata <- data.frame(
  sample_data(phy),
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  rownames_to_column("SampleID")

metadata$Group <- factor(
  as.character(metadata$Group),
  levels = group_levels
)

observed_group_counts <- table(
  factor(metadata$Group, levels = group_levels)
)

if (!all(
  as.integer(observed_group_counts) ==
  expected_group_counts[group_levels]
)) {
  stop(
    "Unexpected group counts. Expected ",
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

cat("\nThe existing final phyloseq object was loaded.\n")
cat("No new phyloseq object was created or saved.\n")

patient_column <- first_existing_column(
  colnames(metadata),
  c(
    "Patient",
    "patient",
    "Patient_ID",
    "patient_id",
    "PatientID"
  )
)

patient_from_sample_id <- sub(
  ".*_([0-9]+)$",
  "\\1",
  metadata$SampleID
)

patient_from_sample_id[
  patient_from_sample_id ==
    metadata$SampleID
] <- NA_character_

if (!is.null(patient_column)) {
  patient <- clean_character(
    metadata[[patient_column]]
  )

  patient[
    is.na(patient)
  ] <- patient_from_sample_id[
    is.na(patient)
  ]
} else {
  patient <- patient_from_sample_id
}

if (anyNA(patient)) {
  stop(
    "Patient IDs could not be resolved for all samples."
  )
}

metadata$Patient <- factor(
  patient
)

rownames(metadata) <- metadata$SampleID

metadata <- metadata[
  sample_names(phy),
  ,
  drop = FALSE
]

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
  tax_table(phy_genus)[
    ,
    "Genus"
  ]
)

if (
  anyNA(genus_labels) ||
  any(genus_labels == "")
) {
  stop(
    "Missing genus labels remain after agglomeration."
  )
}

taxa_names(phy_genus) <- make.unique(
  genus_labels
)

genus_counts <- as(
  otu_table(phy_genus),
  "matrix"
)

if (!taxa_are_rows(phy_genus)) {
  genus_counts <- t(
    genus_counts
  )
}

genus_counts <- round(
  genus_counts
)

genus_taxonomy <- as.data.frame(
  tax_table(phy_genus),
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  rownames_to_column(
    "Taxon"
  ) %>%
  mutate(
    Genus = taxa_names(
      phy_genus
    )
  )

metadata <- metadata[
  colnames(genus_counts),
  ,
  drop = FALSE
]

if (!identical(
  rownames(metadata),
  colnames(genus_counts)
)) {
  stop(
    "Metadata and genus-count columns are not aligned."
  )
}

# Patient-tissue count aggregation.
make_patient_group_pseudobulk <- function(
    group_1,
    group_2,
    paired_only = FALSE
) {
  comparison_groups <- c(
    group_1,
    group_2
  )

  meta_sub <- metadata %>%
    filter(
      Group %in%
        comparison_groups
    ) %>%
    mutate(
      Group = factor(
        as.character(Group),
        levels = c(
          group_2,
          group_1
        )
      )
    )

  if (isTRUE(paired_only)) {
    stop(
      "paired_only = TRUE is prohibited in the main genus-level DESeq2 analysis. ",
      "All three primary contrasts must use every available Patient x Group profile ",
      "with Patient adjustment when patients contribute both groups.",
      call. = FALSE
    )
  }

  meta_sub <- meta_sub %>%
    mutate(
      Patient = droplevels(
        factor(Patient)
      ),

      Unit_ID = paste(
        Patient,
        Group,
        sep = "__"
      )
    )

  unit_ids <- unique(
    meta_sub$Unit_ID
  )

  pseudobulk_counts <- vapply(
    unit_ids,
    function(unit_id) {
      sample_ids <- meta_sub$SampleID[
        meta_sub$Unit_ID ==
          unit_id
      ]

      rowSums(
        genus_counts[
          ,
          sample_ids,
          drop = FALSE
        ]
      )
    },
    FUN.VALUE = numeric(
      nrow(genus_counts)
    )
  )

  if (is.null(dim(pseudobulk_counts))) {
    pseudobulk_counts <- matrix(
      pseudobulk_counts,
      nrow = nrow(genus_counts),
      dimnames = list(
        rownames(genus_counts),
        unit_ids
      )
    )
  } else {
    rownames(pseudobulk_counts) <-
      rownames(genus_counts)

    colnames(pseudobulk_counts) <-
      unit_ids
  }

  unit_metadata <- lapply(
    unit_ids,
    function(unit_id) {
      one <- meta_sub[
        meta_sub$Unit_ID ==
          unit_id,
        ,
        drop = FALSE
      ]

      data.frame(
        Unit_ID = unit_id,
        Patient = as.character(
          one$Patient[1]
        ),
        Group = as.character(
          one$Group[1]
        ),
        N_biopsies = nrow(one),
        Biopsy_sample_IDs = paste(
          one$SampleID,
          collapse = ";"
        ),
        stringsAsFactors = FALSE
      )
    }
  ) %>%
    bind_rows()

  unit_metadata$Group <- factor(
    unit_metadata$Group,
    levels = c(
      group_2,
      group_1
    )
  )

  unit_metadata$Patient <- factor(
    unit_metadata$Patient
  )

  rownames(unit_metadata) <-
    unit_metadata$Unit_ID

  unit_metadata <- unit_metadata[
    colnames(pseudobulk_counts),
    ,
    drop = FALSE
  ]

  keep_taxa <- rowSums(
    pseudobulk_counts
  ) > 0

  pseudobulk_counts <-
    pseudobulk_counts[
      keep_taxa,
      ,
      drop = FALSE
    ]

  if (!identical(
    colnames(pseudobulk_counts),
    rownames(unit_metadata)
  )) {
    stop(
      "Pseudobulk counts and metadata are not aligned."
    )
  }

  list(
    counts = pseudobulk_counts,
    metadata = unit_metadata,
    paired_only = paired_only
  )
}

# Model fitting, complete result tables and volcano plots.
run_deseq2_comparison <- function(
    group_1,
    group_2,
    paired_only = FALSE
) {
  comparison <- paste(
    group_1,
    "-",
    group_2
  )

  filename_tag <- safe_name(
    comparison
  )

  pseudobulk <- make_patient_group_pseudobulk(
    group_1 = group_1,
    group_2 = group_2,
    paired_only = paired_only
  )

  col_data <- pseudobulk$metadata

  expected_unit_ids <- metadata %>%
    filter(as.character(Group) %in% c(group_1, group_2)) %>%
    transmute(Unit_ID = paste(Patient, Group, sep = "__")) %>%
    distinct() %>%
    pull(Unit_ID) %>%
    as.character()

  observed_unit_ids <- rownames(col_data)

  if (
    length(expected_unit_ids) != nrow(col_data) ||
    !setequal(expected_unit_ids, observed_unit_ids)
  ) {
    stop(
      "The ", comparison,
      " model did not contain every available Patient x Group profile. ",
      "The main analysis must not restrict to matched patients.",
      call. = FALSE
    )
  }

  profile_design_audit <- col_data %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    count(Group, name = "N_patient_group_profiles") %>%
    mutate(
      Comparison = comparison,
      Model = if (anyDuplicated(col_data$Patient)) "~ Patient + Group" else "~ Group",
      Uses_all_available_profiles = TRUE,
      Patient_fixed_effect = anyDuplicated(col_data$Patient) > 0L
    )

  write.csv(
    profile_design_audit,
    file.path(
      out_dir,
      paste0("DESeq2_", filename_tag, "_ALL_PROFILE_DESIGN_AUDIT.csv")
    ),
    row.names = FALSE
  )

  if (isTRUE(paired_only)) {
    stop(
      "A paired model was requested in the main genus-level DESeq2 analysis. ",
      "This analysis is locked to all available Patient x Group profiles and Patient adjustment when patients contribute both groups.",
      call. = FALSE
    )
  }

  patient_adjusted <- anyDuplicated(col_data$Patient) > 0L
  design_formula <- if (patient_adjusted) ~ Patient + Group else ~ Group
  model_description <- paste("All available Patient x Group profiles:", paste(deparse(design_formula), collapse = ""))

  design_matrix <- model.matrix(
    design_formula,
    data = col_data
  )

  if (
    qr(design_matrix)$rank <
    ncol(design_matrix)
  ) {
    stop(
      "The DESeq2 design matrix is not full rank for ",
      comparison,
      "."
    )
  }

  dds <- DESeqDataSetFromMatrix(
    countData = round(
      pseudobulk$counts
    ),
    colData = col_data,
    design = design_formula
  )

  dds <- dds[
    rowSums(
      counts(dds)
    ) > 0,
  ]

  dds <- estimateSizeFactors(
    dds,
    type = "poscounts"
  )

  dds <- tryCatch(
    DESeq(
      dds,
      test = "Wald",
      fitType = "parametric",
      quiet = FALSE
    ),
    error = function(e) {
      message(
        "Parametric dispersion fit failed for ",
        comparison,
        "; retrying with fitType = 'local'."
      )

      DESeq(
        dds,
        test = "Wald",
        fitType = "local",
        quiet = FALSE
      )
    }
  )

  contrast_vector <- c(
    "Group",
    group_1,
    group_2
  )

  result_mle <- results(
    dds,
    contrast = contrast_vector,
    alpha = alpha,
    independentFiltering = TRUE
  )

  result_df <- data.frame(
    Taxon = rownames(
      result_mle
    ),

    baseMean = result_mle$baseMean,

    log2FoldChange =
      result_mle$log2FoldChange,

    lfcSE =
      result_mle$lfcSE,

    Wald_statistic =
      result_mle$stat,

    P_value =
      result_mle$pvalue,

    FDR =
      result_mle$padj,

    stringsAsFactors = FALSE
  ) %>%
    left_join(
      genus_taxonomy,
      by = "Taxon"
    ) %>%
    mutate(
      Comparison = comparison,
      First_group = group_1,
      Second_group = group_2,
      Model = model_description,
      N_model_units = ncol(dds),
      N_patients = n_distinct(
        col_data$Patient
      ),

      Direction = case_when(
        is.na(log2FoldChange) ~
          NA_character_,

        log2FoldChange > 0 ~
          paste0(
            "Higher in ",
            group_1
          ),

        log2FoldChange < 0 ~
          paste0(
            "Higher in ",
            group_2
          ),

        TRUE ~
          "No direction"
      ),

      Significant_FDR_0_05 =
        !is.na(FDR) &
        FDR < alpha
    ) %>%
    arrange(
      FDR,
      P_value,
      desc(
        abs(
          log2FoldChange
        )
      )
    )

  significant_df <- result_df %>%
    filter(
      Significant_FDR_0_05
    )

  write.csv(
    result_df,
    file.path(
      out_dir,
      paste0(
        "DESeq2_",
        filename_tag,
        "_full_results.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )

  write.csv(
    significant_df,
    file.path(
      out_dir,
      paste0(
        "DESeq2_",
        filename_tag,
        "_FDR_significant.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )

  write.csv(
    col_data,
    file.path(
      out_dir,
      paste0(
        "DESeq2_",
        filename_tag,
        "_pseudobulk_manifest.csv"
      )
    ),
    row.names = FALSE
  )

  saveRDS(
    dds,
    file.path(
      out_dir,
      paste0(
        "DESeq2_",
        filename_tag,
        "_dds.rds"
      )
    )
  )

  # Unavailable adjusted P values remain in the results table, not at FDR = 1.
  volcano_df <- result_df %>%
    filter(is.finite(FDR), FDR >= 0, FDR <= 1, is.finite(log2FoldChange)) %>%
    mutate(
      FDR_plot = pmax(FDR, .Machine$double.xmin),

      Label = ifelse(
        is.na(Genus) | Genus == "",
        Taxon,
        Genus
      ),

      Plot_status = case_when(
        Significant_FDR_0_05 &
          is.finite(log2FoldChange) &
          log2FoldChange > 0 ~ "Up",

        Significant_FDR_0_05 &
          is.finite(log2FoldChange) &
          log2FoldChange < 0 ~ "Down",

        TRUE ~ "NS"
      ),

      Plot_status = factor(
        Plot_status,
        levels = c("Up", "Down", "NS")
      )
    )

  label_df <- volcano_df %>%
    filter(
      Plot_status %in% c("Up", "Down"),
      !is.na(FDR),
      FDR < alpha,
      !is.na(Label),
      Label != "",
      is.finite(log2FoldChange),
      is.finite(FDR_plot)
    ) %>%
    arrange(
      FDR,
      P_value,
      desc(abs(log2FoldChange))
    )

  write.csv(
    label_df %>%
      select(
        any_of(
          c(
            "Label",
            "Genus",
            "Taxon",
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "Wald_statistic",
            "P_value",
            "FDR",
            "Plot_status"
          )
        )
      ),
    file.path(
      out_dir,
      paste0(
        "Volcano_DESeq2_",
        filename_tag,
        "_labelled_genera.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )

  p_volcano <- ggplot(
    volcano_df,
    aes(
      x = log2FoldChange,
      y = -log10(FDR_plot)
    )
  ) +
    geom_point(
      aes(colour = Plot_status),
      size = 1.45,
      alpha = 0.85
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = Label),
      family = font_family,
      size = pt_to_gg(label_pt),
      colour = "black",
      max.overlaps = Inf,
      force = 5,
      force_pull = 0.08,
      box.padding = 0.30,
      point.padding = 0.20,
      segment.size = 0.20,
      min.segment.length = 0,
      seed = 1
    ) +
    scale_colour_manual(
      values = volcano_colours,
      breaks = c("Up", "Down", "NS"),
      limits = c("Up", "Down", "NS"),
      drop = FALSE
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0.10, 0.14))
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.14))
    ) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed",
      linewidth = 0.35,
      colour = "grey50"
    ) +
    geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed",
      linewidth = 0.35,
      colour = "grey50"
    ) +
    labs(
      title = paste0(
        group_1,
        " vs ",
        group_2
      ),
      x = "log2 fold change",
      y = "-log10(FDR)",
      colour = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(
      base_family = font_family,
      base_size = font_pt
    ) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",

      legend.text = element_text(
        family = font_family,
        size = legend_pt,
        colour = "black"
      ),

      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.key.size = grid::unit(3.2, "mm"),
      legend.spacing.x = grid::unit(2, "mm"),

      plot.title = element_text(
        family = font_family,
        face = "bold",
        size = font_pt,
        colour = "black",
        hjust = 0.5
      ),

      axis.title.x = element_text(
        family = font_family,
        size = font_pt,
        colour = "black",
        margin = margin(t = 4)
      ),

      axis.title.y = element_text(
        family = font_family,
        size = font_pt,
        colour = "black",
        margin = margin(r = 4)
      ),

      axis.text.x = element_text(
        family = font_family,
        size = font_pt,
        colour = "black"
      ),

      axis.text.y = element_text(
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

      axis.ticks.length = grid::unit(1.2, "mm"),

      plot.background = element_rect(
        fill = "white",
        colour = NA
      ),

      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),

      plot.margin = margin(
        5,
        10,
        5,
        6,
        "mm"
      )
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(
          size = 2.2,
          alpha = 1
        )
      )
    )

  save_ggplot(
    p_volcano,
    paste0(
      "Volcano_DESeq2_",
      filename_tag
    ),
    width_mm = volcano_width_mm,
    height_mm = volcano_height_mm
  )

  pdf_device <- if (capabilities("cairo")) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }

  pdf(
    file.path(
      out_dir,
      paste0(
        "MA_DESeq2_",
        filename_tag,
        ".pdf"
      )
    ),
    width = 6,
    height = 5
  )

  plotMA(
    result_mle,
    alpha = alpha,
    main = paste0(
      group_1,
      " vs ",
      group_2
    )
  )

  dev.off()

  png(
    file.path(
      out_dir,
      paste0(
        "MA_DESeq2_",
        filename_tag,
        ".png"
      )
    ),
    width = 6,
    height = 5,
    units = "in",
    res = export_dpi,
    type = if (
      capabilities("cairo")
    ) {
      "cairo"
    } else {
      getOption(
        "bitmapType"
      )
    }
  )

  plotMA(
    result_mle,
    alpha = alpha,
    main = paste0(
      group_1,
      " vs ",
      group_2
    )
  )

  dev.off()

  pdf(
    file.path(
      out_dir,
      paste0(
        "Dispersion_DESeq2_",
        filename_tag,
        ".pdf"
      )
    ),
    width = 6,
    height = 5
  )

  plotDispEsts(
    dds
  )

  dev.off()

  summary_row <- data.frame(
    Comparison = comparison,
    First_group = group_1,
    Second_group = group_2,
    Patient_fixed_effect = patient_adjusted,
    Model = model_description,
    Log2_fold_change_type = "Unshrunk DESeq2 maximum-likelihood estimate",
    N_model_units = ncol(dds),
    N_patients = n_distinct(
      col_data$Patient
    ),
    N_genera_tested = nrow(result_df),
    N_FDR_significant = nrow(
      significant_df
    ),
    Minimum_raw_P = ifelse(
      all(
        is.na(
          result_df$P_value
        )
      ),
      NA_real_,
      min(
        result_df$P_value,
        na.rm = TRUE
      )
    ),
    Minimum_FDR = ifelse(
      all(
        is.na(
          result_df$FDR
        )
      ),
      NA_real_,
      min(
        result_df$FDR,
        na.rm = TRUE
      )
    ),
    stringsAsFactors = FALSE
  )

  list(
    result = result_df,
    significant = significant_df,
    summary = summary_row,
    dds = dds,
    manifest = col_data
  )
}

comparison_specifications <- list(
  list(
    group_1 = "SSL_HP",
    group_2 = "H_COL",
    paired_only = FALSE
  ),

  list(
    group_1 = "H_ADJ_SSL_HP",
    group_2 = "H_COL",
    paired_only = FALSE
  ),

  list(
    group_1 = "SSL_HP",
    group_2 = "H_ADJ_SSL_HP",
    paired_only = FALSE
  )
)

if (any(vapply(
  comparison_specifications,
  function(x) isTRUE(x$paired_only),
  logical(1)
))) {
  stop(
    "At least one main DESeq2 contrast was configured as paired. ",
    "All three main contrasts must use every available Patient x Group profile ",
    "with Patient adjustment when patients contribute both groups.",
    call. = FALSE
  )
}

main_contrast_design_audit <- bind_rows(lapply(
  comparison_specifications,
  function(x) {
    data.frame(
      Comparison = paste(x$group_1, "-", x$group_2),
      First_group = x$group_1,
      Second_group = x$group_2,
      Uses_all_available_profiles = TRUE,
      Patient_fixed_effect = x$group_2 == "H_ADJ_SSL_HP",
      Model = if (x$group_2 == "H_ADJ_SSL_HP") "~ Patient + Group" else "~ Group",
      stringsAsFactors = FALSE
    )
  }
))

write.csv(
  main_contrast_design_audit,
  file.path(out_dir, "MAIN_CONTRAST_DESIGN_AUDIT.csv"),
  row.names = FALSE
)

deseq_results <- lapply(
  comparison_specifications,
  function(specification) {
    run_deseq2_comparison(
      group_1 =
        specification$group_1,

      group_2 =
        specification$group_2,

      paired_only =
        specification$paired_only
    )
  }
)

names(deseq_results) <- vapply(
  comparison_specifications,
  function(x) {
    paste(
      x$group_1,
      "-",
      x$group_2
    )
  },
  character(1)
)

summary_table <- bind_rows(
  lapply(
    deseq_results,
    function(x) x$summary
  )
)

all_results <- bind_rows(
  lapply(
    deseq_results,
    function(x) x$result
  )
)

all_significant <- bind_rows(
  lapply(
    deseq_results,
    function(x) x$significant
  )
)

write.csv(
  summary_table,
  file.path(
    out_dir,
    "DESeq2_primary_result_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  all_results,
  file.path(
    out_dir,
    "DESeq2_all_three_comparisons_full_results.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  all_significant,
  file.path(
    out_dir,
    "DESeq2_all_three_comparisons_FDR_significant.csv"
  ),
  row.names = FALSE,
  na = ""
)

if (isTRUE(make_ancom_concordance)) {
  latest_path_file <- file.path(
    ancom_root,
    "LATEST_RUN_PATH.txt"
  )

  if (file.exists(latest_path_file)) {
    ancom_latest_dir <- readLines(
      latest_path_file,
      warn = FALSE
    )[1]

    ancom_result_file <- file.path(
      ancom_latest_dir,
      "ANCOMBC2_genus_pairwise_tidy.csv"
    )

    if (file.exists(ancom_result_file)) {
      ancom_results <- read.csv(
        ancom_result_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      concordance_summary <- list()

      for (
        comparison_name in names(
          deseq_results
        )
      ) {
        deseq_one <- deseq_results[[
          comparison_name
        ]]$result

        ancom_one <- ancom_results %>%
          filter(
            .data$comparison ==
              comparison_name
          )

        concordance <- deseq_one %>%
          select(
            Genus,
            DESeq2_log2FC =
              log2FoldChange,
            DESeq2_P =
              P_value,
            DESeq2_FDR =
              FDR,
            DESeq2_significant =
              Significant_FDR_0_05
          ) %>%
          inner_join(
            ancom_one %>%
              transmute(
                Genus = taxon,
                ANCOM_log2FC =
                  lfc /
                  log(2),
                ANCOM_P =
                  p_value,
                ANCOM_q =
                  q_value,
                ANCOM_sensitivity_stable =
                  sensitivity_stable_difference
              ),
            by = "Genus"
          ) %>%
          filter(
            is.finite(
              DESeq2_log2FC
            ),
            is.finite(
              ANCOM_log2FC
            )
          ) %>%
          mutate(
            Label = ifelse(
              DESeq2_significant |
                ANCOM_sensitivity_stable,
              Genus,
              NA_character_
            )
          )

        filename_tag <- safe_name(
          comparison_name
        )

        write.csv(
          concordance,
          file.path(
            out_dir,
            paste0(
              "ANCOMBC2_DESeq2_concordance_",
              filename_tag,
              ".csv"
            )
          ),
          row.names = FALSE
        )

        correlation <- if (
          nrow(concordance) >= 3L
        ) {
          suppressWarnings(
            cor.test(
              concordance$ANCOM_log2FC,
              concordance$DESeq2_log2FC,
              method = "spearman",
              exact = FALSE
            )
          )
        } else {
          NULL
        }

        rho <- if (
          is.null(correlation)
        ) {
          NA_real_
        } else {
          unname(
            correlation$estimate
          )
        }

        correlation_p <- if (
          is.null(correlation)
        ) {
          NA_real_
        } else {
          correlation$p.value
        }

        concordance_summary[[
          comparison_name
        ]] <- data.frame(
          Comparison = comparison_name,
          N_common_genera =
            nrow(concordance),
          Spearman_rho =
            rho,
          P_value =
            correlation_p,
          stringsAsFactors = FALSE
        )

        if (nrow(concordance) > 0L) {
          p_concordance <- ggplot(
            concordance,
            aes(
              x = ANCOM_log2FC,
              y = DESeq2_log2FC
            )
          ) +
            geom_hline(
              yintercept = 0,
              colour = "grey65",
              linewidth = 0.30
            ) +
            geom_vline(
              xintercept = 0,
              colour = "grey65",
              linewidth = 0.30
            ) +
            geom_abline(
              intercept = 0,
              slope = 1,
              linetype = "dashed",
              linewidth = 0.35,
              colour = "grey45"
            ) +
            geom_point(
              size = 1.5,
              alpha = 0.75
            ) +
            ggrepel::geom_text_repel(
              aes(
                label = Label
              ),
              family = font_family,
              size = pt_to_gg(
                label_pt
              ),
              max.overlaps = Inf,
              seed = 1,
              na.rm = TRUE
            ) +
            labs(
              title = paste0(
                "ANCOM-BC2 vs DESeq2: ",
                comparison_name
              ),
              subtitle = paste0(
                "Spearman rho = ",
                ifelse(
                  is.na(rho),
                  "NA",
                  sprintf(
                    "%.2f",
                    rho
                  )
                ),
                "; p = ",
                ifelse(
                  is.na(correlation_p),
                  "NA",
                  formatC(
                    correlation_p,
                    format = "g",
                    digits = 3
                  )
                )
              ),
              x =
                "ANCOM-BC2 log2 fold change",
              y =
                "DESeq2 unshrunk log2 fold change"
            ) +
            theme_classic(
              base_family = font_family,
              base_size = font_pt
            ) +
            theme(
              plot.title = element_text(
                face = "bold",
                hjust = 0.5
              ),
              plot.subtitle = element_text(
                hjust = 0.5,
                size = 5.5
              )
            )

          save_ggplot(
            p_concordance,
            paste0(
              "Supplementary_ANCOMBC2_DESeq2_concordance_",
              filename_tag
            ),
            width_mm = 100,
            height_mm = 85
          )
        }
      }

      write.csv(
        bind_rows(
          concordance_summary
        ),
        file.path(
          out_dir,
          "ANCOMBC2_DESeq2_concordance_summary.csv"
        ),
        row.names = FALSE
      )
    } else {
      warning(
        "ANCOM-BC2 tidy result file was not found; concordance plots were skipped."
      )
    }
  } else {
    warning(
      "ANCOM-BC2 LATEST_RUN_PATH.txt was not found; concordance plots were skipped."
    )
  }
}

writeLines(
  c(
    "FINAL SECTION 4.5.2 DESEQ2 PRIMARY ANALYSIS",
    "===============================================",
    "",
    paste0("Input: ", physeq_path),
    paste0("DESeq2 version: ", packageVersion("DESeq2")),
    "",
    "Role:",
    "  Primary genus-level differential-abundance analysis; ANCOM-BC2 is the compositional sensitivity analysis.",
    "",
    "Input counts:",
    "  Raw integer genus counts derived in memory from the existing final 5%-filtered object.",
    "  Size factors estimated with type = poscounts.",
    "",
    "Repeated-biopsy handling:",
    "  Counts summed within Patient x Group to create patient-level pseudobulks.",
    "  No patient-matching restriction is applied in any of the three main contrasts.",
    "",
    "Models:",
    "  SSL_HP vs H_COL: ~ Group",
    "  H_ADJ_SSL_HP vs H_COL: ~ Group",
    "  SSL_HP vs H_ADJ_SSL_HP: all available Patient x Group profiles; ~ Patient + Group",
    "  Age is not included in any primary DESeq2 model.",
    "",
    "Inference:",
    "  Wald tests.",
    "  BH FDR within each contrast.",
    "  Significance is based on DESeq2 P values/FDR; LFC shrinkage is not used.",
    "  Original unshrunk maximum-likelihood log2 fold changes are used in plots and tables.",
    "  BH FDR is controlled separately within each prespecified contrast; no additional cross-contrast correction is applied.",
    "",
    "Figures:",
    "  Three volcano plots using unshrunk log2 fold changes.",
    "  Volcano colours: red Up, green Down and grey NS.",
    "  Every FDR-significant genus is labelled; no top-N restriction is used.",
    "  max.overlaps = Inf prevents ggrepel from silently dropping labels.",
    "  Three MA plots using unshrunk log2 fold changes.",
    "  Three dispersion-estimate plots.",
    "  Optional ANCOM-BC2/DESeq2 effect-direction concordance plots."
  ),
  file.path(
    out_dir,
    "README_DESeq2_PRIMARY_NO_AGE_FINAL.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "sessionInfo_DESeq2_PRIMARY_NO_AGE_FINAL.txt"
  )
)

cat("\n===== DESEQ2 PRIMARY SUMMARY =====\n")
print(summary_table)

cat("\nResults saved to:\n")
cat(out_dir, "\n")
