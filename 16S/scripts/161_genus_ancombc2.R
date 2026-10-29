#!/usr/bin/env Rscript
rm(list = ls())
set.seed(20260823)

options(stringsAsFactors = FALSE)
options(contrasts = c("contr.treatment", "contr.poly"))

suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(lme4)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# Version-specific missing-value alignment guard; retained for reproducibility.
patch_ancombc_bias_em_if_needed <- function() {
  namespace <- asNamespace("ANCOMBC")

  if (!exists(".bias_em", envir = namespace, inherits = FALSE)) {
    stop(
      "The installed ANCOMBC namespace does not contain .bias_em. ",
      "Do not run this script with an unknown package implementation."
    )
  }

  installed_version <- packageVersion("ANCOMBC")
  bias_em_function <- get(
    ".bias_em",
    envir = namespace,
    inherits = FALSE
  )

  source_before <- paste(
    deparse(body(bias_em_function)),
    collapse = "\n"
  )

  already_corrected <- grepl(
    "neither_na",
    source_before,
    fixed = TRUE
  )

  if (already_corrected) {
    return(
      data.frame(
        ANCOMBC_version = as.character(installed_version),
        Bias_EM_patch_required = FALSE,
        Bias_EM_patch_applied_in_memory = FALSE,
        Patch_status = "Installed ANCOMBC already contains the corrected NA alignment",
        stringsAsFactors = FALSE
      )
    )
  }

  if (installed_version != package_version("2.4.0")) {
    stop(
      "The installed ANCOMBC version is ",
      as.character(installed_version),
      " and its .bias_em function does not contain the corrected NA alignment. ",
      "This script will patch only the audited ANCOMBC 2.4.0 implementation. ",
      "Install a current official ANCOMBC release before running the analysis."
    )
  }

  body_parts <- as.list(
    body(bias_em_function)
  )

  if (length(body_parts) < 5L) {
    stop("Unexpected ANCOMBC 2.4.0 .bias_em function structure.")
  }

  normalise_expression <- function(expression_object) {
    output <- paste(
      deparse(expression_object),
      collapse = ""
    )

    output <- gsub(
      "\\s+",
      "",
      output
    )

    output <- gsub(
      "<-",
      "=",
      output,
      fixed = TRUE
    )

    output
  }

  observed_opening <- vapply(
    body_parts[2:4],
    normalise_expression,
    character(1)
  )

  expected_opening <- c(
    "beta=beta[!is.na(beta)]",
    "nu0=var_hat",
    "nu0=nu0[!is.na(nu0)]"
  )

  if (!identical(observed_opening, expected_opening)) {
    stop(
      "The ANCOMBC 2.4.0 .bias_em function did not match the audited buggy ",
      "opening statements. The namespace was not modified."
    )
  }

  corrected_body_parts <- c(
    body_parts[1],
    list(
      quote(nu0 <- var_hat),
      quote(neither_na <- !(is.na(beta) | is.na(nu0))),
      quote(beta <- beta[neither_na]),
      quote(nu0 <- nu0[neither_na])
    ),
    body_parts[5:length(body_parts)]
  )

  body(bias_em_function) <- as.call(
    corrected_body_parts
  )

  assignInNamespace(
    x = ".bias_em",
    value = bias_em_function,
    ns = "ANCOMBC"
  )

  source_after <- paste(
    deparse(
      body(
        get(
          ".bias_em",
          envir = asNamespace("ANCOMBC"),
          inherits = FALSE
        )
      )
    ),
    collapse = "\n"
  )

  if (!grepl("neither_na", source_after, fixed = TRUE)) {
    stop("The ANCOMBC 2.4.0 .bias_em in-memory patch could not be verified.")
  }

  data.frame(
    ANCOMBC_version = as.character(installed_version),
    Bias_EM_patch_required = TRUE,
    Bias_EM_patch_applied_in_memory = TRUE,
    Patch_status = paste0(
      "Applied the corrected joint NA mask to beta and var_hat in .bias_em; ",
      "installed package files were not modified"
    ),
    stringsAsFactors = FALSE
  )
}

config <- new.env(parent = baseenv())
if (file.exists("config/config.R")) sys.source("config/config.R", envir = config)
work_dir <- if (exists("PATH_WORK", config, inherits = FALSE)) config$PATH_WORK else "outputs"
if (!dir.exists(work_dir)) stop("Input directory does not exist: ", work_dir)

base_dir <- work_dir

physeq_path <- file.path(
  base_dir,
  "phyloseq_final_March26_5pct_without7699.rds"
)

out_root <- file.path(
  base_dir,
  "Section4_5_3_ANCOMBC2_GENUS_ALL_PROFILES_THIRD_FINAL"
)

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

run_id <- paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
out_dir <- file.path(out_root, run_id)
if (dir.exists(out_dir)) stop("Output folder already exists: ", out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(
  out_dir,
  file.path(out_root, "LATEST_RUN_PATH.txt")
)

ancombc_patch_status <- patch_ancombc_bias_em_if_needed()

write.csv(
  ancombc_patch_status,
  file.path(out_dir, "00_ANCOMBC_BIAS_EM_PATCH_STATUS.csv"),
  row.names = FALSE
)

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
p_adjust_method <- "BH"
n_cores <- 1L

ancom_prv_cut <- 0

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "na", "N/A", "Unassigned", "unassigned")] <- NA_character_
  x
}

make_group <- function(x) {
  group <- clean_character(x)

  group[group %in% c("SSL", "HP", "SSL_HP")] <- "SSL_HP"

  group[
    group %in% c("H_ADJ_SSL", "H_ADJ_HP", "H_ADJ_SSL_HP")
  ] <- "H_ADJ_SSL_HP"

  group[group %in% c("H_COL", "H_COL_D", "H_COL_ND")] <- "H_COL"
  group[!group %in% group_levels] <- NA_character_

  factor(group, levels = group_levels)
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

as_logical_flag <- function(x) {
  if (is.logical(x)) {
    return(!is.na(x) & x)
  }

  !is.na(x) & toupper(as.character(x)) %in% c("TRUE", "T", "1")
}

as_result_df <- function(x) {
  if (is.null(x)) return(data.frame())

  out <- as.data.frame(x, check.names = FALSE)

  if (!"taxon" %in% colnames(out)) {
    out <- rownames_to_column(out, var = "taxon")
  }

  out
}

companion_column <- function(df, prefix, term, default = NA) {
  column_name <- paste0(prefix, term)

  if (column_name %in% colnames(df)) {
    return(df[[column_name]])
  }

  rep(default, nrow(df))
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

  genus[missing_genus & !is.na(parent)] <- paste0(
    "Unclassified_",
    gsub(
      "[^[:alnum:]_]+",
      "_",
      parent[missing_genus & !is.na(parent)]
    )
  )

  genus[missing_genus & is.na(parent)] <- paste0(
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

  key_map <- paste(map$Original_genus, map$Lineage_key, sep = "\r")
  key_tax <- paste(genus, lineage_key, sep = "\r")
  final_genus <- map$Analysis_genus[match(key_tax, key_map)]

  if (anyNA(final_genus)) {
    stop("At least one ASV could not be assigned a lineage-aware genus label.")
  }

  list(genus = final_genus, mapping = map)
}

if (!file.exists(physeq_path)) {
  stop("The final phyloseq object was not found:\n", physeq_path)
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
if (!taxa_are_rows(phy)) otu_check <- t(otu_check)

if (
  any(!is.finite(otu_check)) ||
  any(otu_check < 0) ||
  any(abs(otu_check - round(otu_check)) > 1e-8)
) {
  stop("ANCOM-BC2 requires non-negative integer raw counts.")
}

sample_data(phy)$Group <- make_group(sample_data(phy)$Sample_Type)
phy <- subset_samples(phy, !is.na(Group))
phy <- prune_taxa(taxa_sums(phy) > 0, phy)

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

patient_column <- first_existing_column(
  colnames(metadata),
  c(
    "Patient",
    "patient",
    "Patient_ID",
    "patient_id",
    "PatientID",
    "Subject",
    "Subject_ID"
  )
)

patient_from_sample_id <- sub(
  ".*_([0-9]+)$",
  "\\1",
  metadata$SampleID
)

patient_from_sample_id[
  patient_from_sample_id == metadata$SampleID
] <- NA_character_

if (!is.null(patient_column)) {
  patient <- clean_character(metadata[[patient_column]])
  patient[is.na(patient)] <- patient_from_sample_id[is.na(patient)]
  patient_source <- patient_column
} else {
  patient <- patient_from_sample_id
  patient_source <- "Final numeric token in SampleID"
}

if (anyNA(patient)) {
  stop("Patient IDs could not be resolved for all samples.")
}

metadata$Patient <- factor(patient)
rownames(metadata) <- metadata$SampleID
metadata <- metadata[sample_names(phy), , drop = FALSE]

sample_data(phy)$Group <- factor(metadata$Group, levels = group_levels)
sample_data(phy)$Patient <- factor(metadata$Patient)

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
  file.path(out_dir, "00_genus_lineage_mapping.csv"),
  row.names = FALSE
)

phy_genus <- tax_glom(
  phy,
  taxrank = "Genus",
  NArm = FALSE
)

phy_genus <- prune_taxa(taxa_sums(phy_genus) > 0, phy_genus)

genus_labels <- as.character(tax_table(phy_genus)[, "Genus"])

if (anyNA(genus_labels) || any(genus_labels == "")) {
  stop("Missing genus labels remain after genus-level agglomeration.")
}

taxa_names(phy_genus) <- make.unique(genus_labels)

genus_counts <- as(otu_table(phy_genus), "matrix")
if (!taxa_are_rows(phy_genus)) genus_counts <- t(genus_counts)
genus_counts <- round(genus_counts)

genus_taxonomy <- as(tax_table(phy_genus), "matrix")

metadata <- metadata[colnames(genus_counts), , drop = FALSE]

if (!identical(rownames(metadata), colnames(genus_counts))) {
  stop("Metadata and genus-count columns are not aligned.")
}

write.csv(
  data.frame(
    Input_object = physeq_path,
    Samples = nsamples(phy),
    Patients = nlevels(metadata$Patient),
    ASVs_in_input_object = ntaxa(phy),
    Genus_level_units = nrow(genus_counts),
    ANCOMBC2_prv_cut = ancom_prv_cut,
    Additional_prevalence_filter_applied = FALSE,
    Patient_source = patient_source,
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "00_input_summary.csv"),
  row.names = FALSE
)

# Patient-tissue count aggregation.
make_patient_group_pseudobulk <- function(
    group_1,
    group_2,
    paired_only = FALSE
) {
  comparison_groups <- c(group_1, group_2)

  meta_sub <- metadata %>%
    filter(as.character(Group) %in% comparison_groups) %>%
    mutate(
      Group = factor(
        as.character(Group),
        levels = c(group_2, group_1)
      ),
      Patient = factor(Patient)
    )

  if (isTRUE(paired_only)) {
    stop(
      "paired_only = TRUE is prohibited in the main genus-level ANCOM-BC2 analysis. ",
      "All three primary contrasts must use every available Patient x Group profile ",
      "with Patient adjustment when patients contribute both groups.",
      call. = FALSE
    )
  }

  meta_sub <- meta_sub %>%
    mutate(
      Patient = droplevels(factor(Patient)),
      Group = factor(as.character(Group), levels = c(group_2, group_1)),
      Unit_ID = paste(Patient, Group, sep = "__")
    ) %>%
    arrange(Group, Patient, SampleID)

  unit_manifest <- meta_sub %>%
    group_by(Unit_ID, Patient, Group) %>%
    summarise(
      N_biopsies = n(),
      Biopsy_sample_IDs = paste(SampleID, collapse = ";"),
      .groups = "drop"
    ) %>%
    mutate(
      Patient = droplevels(factor(Patient)),
      Group = factor(as.character(Group), levels = c(group_2, group_1))
    )

  unit_ids <- unit_manifest$Unit_ID

  pseudobulk_counts <- vapply(
    unit_ids,
    function(unit_id) {
      sample_ids <- meta_sub$SampleID[meta_sub$Unit_ID == unit_id]

      rowSums(
        genus_counts[, sample_ids, drop = FALSE]
      )
    },
    FUN.VALUE = numeric(nrow(genus_counts))
  )

  if (is.null(dim(pseudobulk_counts))) {
    pseudobulk_counts <- matrix(
      pseudobulk_counts,
      nrow = nrow(genus_counts),
      dimnames = list(rownames(genus_counts), unit_ids)
    )
  } else {
    rownames(pseudobulk_counts) <- rownames(genus_counts)
    colnames(pseudobulk_counts) <- unit_ids
  }

  keep_taxa <- rowSums(pseudobulk_counts) > 0
  pseudobulk_counts <- pseudobulk_counts[keep_taxa, , drop = FALSE]

  unit_manifest <- as.data.frame(unit_manifest, stringsAsFactors = FALSE)
  rownames(unit_manifest) <- unit_manifest$Unit_ID
  unit_manifest <- unit_manifest[colnames(pseudobulk_counts), , drop = FALSE]

  unit_manifest$Patient <- droplevels(factor(unit_manifest$Patient))
  unit_manifest$Group <- factor(
    as.character(unit_manifest$Group),
    levels = c(group_2, group_1)
  )

  if (!identical(colnames(pseudobulk_counts), rownames(unit_manifest))) {
    stop("Pseudobulk counts and unit metadata are not aligned.")
  }

  tax_sub <- genus_taxonomy[rownames(pseudobulk_counts), , drop = FALSE]

  phy_pseudobulk <- phyloseq(
    otu_table(pseudobulk_counts, taxa_are_rows = TRUE),
    tax_table(tax_sub),
    sample_data(unit_manifest)
  )

  list(
    phyloseq = phy_pseudobulk,
    counts = pseudobulk_counts,
    metadata = unit_manifest,
    paired_only = paired_only
  )
}

make_prevalence_tables <- function(pb_counts, pb_metadata) {
  groups_present <- levels(pb_metadata$Group)

  prevalence_long <- bind_rows(
    lapply(
      groups_present,
      function(group_name) {
        unit_ids <- rownames(pb_metadata)[
          as.character(pb_metadata$Group) == group_name
        ]

        data.frame(
          taxon = rownames(pb_counts),
          Group = group_name,
          N_patient_tissue_units = length(unit_ids),
          N_detected = rowSums(pb_counts[, unit_ids, drop = FALSE] > 0),
          stringsAsFactors = FALSE
        ) %>%
          mutate(
            Prevalence = N_detected / N_patient_tissue_units
          )
      }
    )
  )

  prevalence_wide <- prevalence_long %>%
    pivot_wider(
      names_from = Group,
      values_from = c(
        N_patient_tissue_units,
        N_detected,
        Prevalence
      ),
      names_glue = "{.value}_{Group}"
    )

  list(long = prevalence_long, wide = prevalence_wide)
}

# Fit each contrast and retain test availability, structural zeros and warnings.
run_ancombc2_contrast <- function(
    group_1,
    group_2,
    paired_only = FALSE
) {
  comparison <- paste(group_1, "-", group_2)
  tag <- safe_name(comparison)
  comparison_dir <- file.path(out_dir, tag)
  dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)

  pb <- make_patient_group_pseudobulk(
    group_1 = group_1,
    group_2 = group_2,
    paired_only = paired_only
  )

  pb_phy <- pb$phyloseq
  pb_counts <- pb$counts
  pb_metadata <- pb$metadata

  expected_unit_ids <- metadata %>%
    filter(as.character(Group) %in% c(group_1, group_2)) %>%
    transmute(Unit_ID = paste(Patient, Group, sep = "__")) %>%
    distinct() %>%
    pull(Unit_ID) %>%
    as.character()

  observed_unit_ids <- rownames(pb_metadata)

  if (
    length(expected_unit_ids) != nrow(pb_metadata) ||
    !setequal(expected_unit_ids, observed_unit_ids)
  ) {
    stop(
      "The ", comparison,
      " ANCOM-BC2 model did not contain every available Patient x Group profile. ",
      "The main analysis must not restrict to matched patients.",
      call. = FALSE
    )
  }

  profile_design_audit <- pb_metadata %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    count(Group, name = "N_patient_group_profiles") %>%
    mutate(
      Comparison = comparison,
      Fixed_effect_formula = if (anyDuplicated(pb_metadata$Patient)) "Patient + Group" else "Group",
      Uses_all_available_profiles = TRUE,
      Patient_fixed_effect = anyDuplicated(pb_metadata$Patient) > 0L
    )

  write.csv(
    profile_design_audit,
    file.path(
      comparison_dir,
      paste0("ANCOMBC2_", tag, "_ALL_PROFILE_DESIGN_AUDIT.csv")
    ),
    row.names = FALSE
  )

  if (isTRUE(paired_only)) {
    stop(
      "A paired model was requested in the main genus-level ANCOM-BC2 analysis. ",
      "This analysis is locked to all available Patient x Group profiles and Patient adjustment when patients contribute both groups.",
      call. = FALSE
    )
  }

  fix_formula <- if (anyDuplicated(pb_metadata$Patient)) "Patient + Group" else "Group"

  design_matrix <- model.matrix(
    as.formula(paste0("~ ", fix_formula)),
    data = pb_metadata
  )

  if (qr(design_matrix)$rank < ncol(design_matrix)) {
    stop("The fixed-effect design matrix is not full rank for ", comparison, ".")
  }

  write.csv(
    pb_metadata,
    file.path(comparison_dir, paste0("ANCOMBC2_", tag, "_pseudobulk_manifest.csv")),
    row.names = FALSE
  )

  write.csv(
    data.frame(
      taxon = rownames(pb_counts),
      pb_counts,
      check.names = FALSE
    ),
    file.path(comparison_dir, paste0("ANCOMBC2_", tag, "_pseudobulk_counts.csv")),
    row.names = FALSE
  )

  prevalence <- make_prevalence_tables(pb_counts, pb_metadata)

  write.csv(
    prevalence$long,
    file.path(comparison_dir, paste0("ANCOMBC2_", tag, "_prevalence_long.csv")),
    row.names = FALSE
  )

  captured_warnings <- character(0)

  fit <- withCallingHandlers(
    ANCOMBC::ancombc2(
      data = pb_phy,
      assay_name = "counts",

      tax_level = NULL,

      fix_formula = fix_formula,
      rand_formula = NULL,

      p_adj_method = p_adjust_method,

      pseudo = 0,
      pseudo_sens = TRUE,

      prv_cut = ancom_prv_cut,
      lib_cut = 0,
      s0_perc = 0.05,

      group = "Group",
      struc_zero = TRUE,
      neg_lb = FALSE,

      alpha = alpha,
      n_cl = n_cores,
      verbose = TRUE,

      global = FALSE,
      pairwise = FALSE,
      dunnet = FALSE,
      trend = FALSE,

      iter_control = list(
        tol = 1e-2,
        max_iter = 20,
        verbose = TRUE
      ),

      em_control = list(
        tol = 1e-5,
        max_iter = 100
      ),

      lme_control = lme4::lmerControl()
    ),
    warning = function(w) {
      captured_warnings <<- c(captured_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  saveRDS(
    fit,
    file.path(comparison_dir, paste0("ANCOMBC2_", tag, "_full_output.rds"))
  )

  captured_warnings <- unique(captured_warnings)

  writeLines(
    if (length(captured_warnings) == 0L) {
      "No warnings were emitted."
    } else {
      captured_warnings
    },
    file.path(comparison_dir, paste0("ANCOMBC2_", tag, "_warnings.txt"))
  )

  res <- as_result_df(fit$res)

  if (nrow(res) == 0L) {
    stop("ANCOM-BC2 returned no primary model-based results for ", comparison, ".")
  }

  group_lfc_columns <- grep(
    "^lfc_Group",
    colnames(res),
    value = TRUE
  )

  if (length(group_lfc_columns) != 1L) {
    stop(
      "Expected exactly one Group coefficient in fit$res for ",
      comparison,
      ", but found: ",
      paste(group_lfc_columns, collapse = ", "),
      "."
    )
  }

  lfc_column <- group_lfc_columns[1]
  term <- sub("^lfc_", "", lfc_column)

  model_results <- data.frame(
    taxon = as.character(res$taxon),
    Comparison = comparison,
    First_group = group_1,
    Second_group = group_2,
    Patient_fixed_effect = anyDuplicated(pb_metadata$Patient) > 0L,
    Fixed_effect_formula = fix_formula,

    LFC_natural_log = as.numeric(res[[lfc_column]]),
    SE_natural_log = as.numeric(companion_column(res, "se_", term)),
    W = as.numeric(companion_column(res, "W_", term)),
    P_value = as.numeric(companion_column(res, "p_", term)),

    BH_FDR = as.numeric(companion_column(res, "q_", term)),

    ANCOMBC_diff_abn_flag = as_logical_flag(
      companion_column(res, "diff_", term, FALSE)
    ),

    passed_pseudocount_sensitivity = as_logical_flag(
      companion_column(res, "passed_ss_", term, FALSE)
    ),

    stringsAsFactors = FALSE
  ) %>%
    mutate(

      Holm_adjusted_P_sensitivity = p.adjust(
        P_value,
        method = "holm"
      ),

      LFC_log2 = LFC_natural_log / log(2),
      CI95_low_natural_log = LFC_natural_log - 1.96 * SE_natural_log,
      CI95_high_natural_log = LFC_natural_log + 1.96 * SE_natural_log,
      CI95_low_log2 = CI95_low_natural_log / log(2),
      CI95_high_log2 = CI95_high_natural_log / log(2),

      BH_significant = !is.na(BH_FDR) & BH_FDR < alpha,

      Holm_significant_sensitivity =
        !is.na(Holm_adjusted_P_sensitivity) &
        Holm_adjusted_P_sensitivity < alpha,

      model_based_BH_and_pseudocount_stable =
        ANCOMBC_diff_abn_flag &
        BH_significant &
        passed_pseudocount_sensitivity,

      model_based_Holm_and_pseudocount_stable =
        Holm_significant_sensitivity &
        passed_pseudocount_sensitivity,

      ANCOMBC_flag_matches_BH = case_when(
        is.na(BH_FDR) ~ TRUE,
        TRUE ~ ANCOMBC_diff_abn_flag == BH_significant
      ),

      Direction = case_when(
        is.na(LFC_natural_log) ~ NA_character_,
        LFC_natural_log > 0 ~ paste0("Higher in ", group_1),
        LFC_natural_log < 0 ~ paste0("Higher in ", group_2),
        TRUE ~ "No direction"
      )
    ) %>%
    left_join(prevalence$wide, by = "taxon") %>%
    arrange(BH_FDR, P_value, desc(abs(LFC_natural_log)))

  X <- model.matrix(as.formula(paste("~", fix_formula)), pb_metadata)
  ft <- as.matrix(fit$feature_table)
  ft <- ft[model_results$taxon, rownames(pb_metadata), drop = FALSE]
  contrast_row <- rep(0, ncol(X)); contrast_row[match(term, colnames(X))] <- 1
  test_info <- t(vapply(seq_len(nrow(ft)), function(i) {
    Xi <- X[ft[i, ] > 0, , drop = FALSE]
    rank <- qr(Xi)$rank
    c(Residual_df = nrow(Xi) - rank,
      Contrast_estimable = qr(rbind(Xi, contrast_row))$rank == rank)
  }, numeric(2)))
  model_results$Residual_df <- test_info[, "Residual_df"]
  model_results$Tissue_test_estimable <- test_info[, "Contrast_estimable"] == 1 &
    model_results$Residual_df > 0 & is.finite(model_results$LFC_natural_log) &
    is.finite(model_results$SE_natural_log) & model_results$SE_natural_log > 0 &
    is.finite(model_results$W) & is.finite(model_results$P_value)
  model_results$Reporting_P <- ifelse(model_results$Tissue_test_estimable, model_results$P_value, NA_real_)
  model_results$Reporting_BH_FDR <- ifelse(model_results$Tissue_test_estimable, model_results$BH_FDR, NA_real_)
  model_results$model_based_BH_and_pseudocount_stable <- model_results$model_based_BH_and_pseudocount_stable & model_results$Tissue_test_estimable
  model_results$model_based_Holm_and_pseudocount_stable <- model_results$model_based_Holm_and_pseudocount_stable & model_results$Tissue_test_estimable

  model_supported_BH <- model_results %>%
    filter(model_based_BH_and_pseudocount_stable)

  model_supported_Holm <- model_results %>%
    filter(model_based_Holm_and_pseudocount_stable)

  write.csv(
    model_results,
    file.path(
      comparison_dir,
      paste0("ANCOMBC2_", tag, "_model_based_all.csv")
    ),
    row.names = FALSE,
    na = ""
  )

  write.csv(
    model_supported_BH,
    file.path(
      comparison_dir,
      paste0("ANCOMBC2_", tag, "_model_based_BH_FDR_supported.csv")
    ),
    row.names = FALSE,
    na = ""
  )

  write.csv(
    model_supported_Holm,
    file.path(
      comparison_dir,
      paste0("ANCOMBC2_", tag, "_model_based_Holm_sensitivity.csv")
    ),
    row.names = FALSE,
    na = ""
  )

  res_zero <- as_result_df(fit$zero_ind)

  structural_zero_results <- data.frame(
    taxon = character(0),
    Comparison = character(0),
    First_group = character(0),
    Second_group = character(0),
    Patient_fixed_effect = logical(0),
    Evidence_type = character(0),
    Structural_zero_in_group = character(0),
    N_units_in_structural_zero_group = integer(0),
    N_detected_in_structural_zero_group = integer(0),
    Prevalence_in_structural_zero_group = numeric(0),
    Direction = character(0),
    Structural_zero_interpretation = character(0),
    Structural_zero_column = character(0),
    Structural_zero = logical(0),
    stringsAsFactors = FALSE
  )

  if (nrow(res_zero) > 0L) {
    zero_columns <- grep(
      "^structural_zero",
      colnames(res_zero),
      value = TRUE
    )

    if (length(zero_columns) > 0L) {
      structural_zero_results <- res_zero %>%
        select(taxon, all_of(zero_columns)) %>%
        pivot_longer(
          cols = all_of(zero_columns),
          names_to = "Structural_zero_column",
          values_to = "Structural_zero"
        ) %>%
        mutate(
          Structural_zero = as_logical_flag(Structural_zero),
          Structural_zero_in_group = sub(
            "^structural_zero \\(Group = (.*)\\)$",
            "\\1",
            Structural_zero_column
          )
        ) %>%
        filter(
          Structural_zero,
          Structural_zero_in_group %in% c(group_1, group_2)
        ) %>%
        left_join(
          prevalence$long %>%
            rename(
              Structural_zero_in_group = Group,
              N_units_in_structural_zero_group = N_patient_tissue_units,
              N_detected_in_structural_zero_group = N_detected,
              Prevalence_in_structural_zero_group = Prevalence
            ),
          by = c("taxon", "Structural_zero_in_group")
        ) %>%
        left_join(
          prevalence$wide,
          by = "taxon"
        ) %>%
        mutate(
          Comparison = comparison,
          First_group = group_1,
          Second_group = group_2,
          Patient_fixed_effect = anyDuplicated(pb_metadata$Patient) > 0L,

          Evidence_type = paste0(
            "Structural zero; group-restricted detection; ",
            "no finite model-based LFC or adjusted P value"
          ),

          Direction = case_when(
            Structural_zero_in_group == group_2 ~
              paste0("Detected in ", group_1, "; structural zero in ", group_2),

            Structural_zero_in_group == group_1 ~
              paste0("Detected in ", group_2, "; structural zero in ", group_1),

            TRUE ~ NA_character_
          ),

          Structural_zero_interpretation = paste0(
            "Group-level structural-zero screening; all available Patient x Group ",
            "profiles were retained and no patient matching was applied"
          )
        ) %>%
        select(
          taxon,
          Comparison,
          First_group,
          Second_group,
          Patient_fixed_effect,
          Evidence_type,
          Structural_zero_in_group,
          N_units_in_structural_zero_group,
          N_detected_in_structural_zero_group,
          Prevalence_in_structural_zero_group,
          Direction,
          Structural_zero_interpretation,
          everything()
        )
    }
  }

  write.csv(
    structural_zero_results,
    file.path(
      comparison_dir,
      paste0("ANCOMBC2_", tag, "_structural_zero_findings.csv")
    ),
    row.names = FALSE,
    na = ""
  )

  critical_pattern <- paste(
    c(
      "NaNs produced",
      "longer object length is not a multiple of shorter object length",
      "failed to converge",
      "convergence failure",
      "singular Hessian",
      "numerically singular",
      "near-zero eigenvalue",
      "not positive definite",
      "trust region",
      "bobyqa",
      "zero variances have been detected",
      "task [0-9]+ failed",
      "all scheduled cores encountered errors"
    ),
    collapse = "|"
  )

  critical_warnings <- captured_warnings[
    grepl(
      critical_pattern,
      captured_warnings,
      ignore.case = TRUE,
      perl = TRUE
    )
  ]

  noncritical_warnings <- setdiff(
    captured_warnings,
    critical_warnings
  )

  finite_test_rows <-
    is.finite(model_results$W) &
    is.finite(model_results$P_value)

  degenerate_all_p_one <-
    any(finite_test_rows) &&
    all(abs(model_results$P_value[finite_test_rows] - 1) < 1e-12) &&
    any(abs(model_results$W[finite_test_rows]) > 1e-8)

  degenerate_all_p_missing <- all(is.na(model_results$P_value))

  invalid_probability_values <-
    any(
      model_results$P_value < 0 |
        model_results$P_value > 1,
      na.rm = TRUE
    ) ||
    any(
      model_results$BH_FDR < 0 |
        model_results$BH_FDR > 1,
      na.rm = TRUE
    )

  n_flag_mismatches <- sum(
    !model_results$ANCOMBC_flag_matches_BH,
    na.rm = TRUE
  )

  supported_rows_finite <- all(
    is.finite(model_supported_BH$LFC_natural_log) &
      is.finite(model_supported_BH$SE_natural_log) &
      is.finite(model_supported_BH$P_value) &
      is.finite(model_supported_BH$BH_FDR)
  )

  critical_result_problem <-
    degenerate_all_p_one ||
    degenerate_all_p_missing ||
    invalid_probability_values ||
    n_flag_mismatches > 0L ||
    !supported_rows_finite

  run_status <- if (
    length(critical_warnings) > 0L ||
    critical_result_problem
  ) {
    "REVIEW_REQUIRED"
  } else if (length(noncritical_warnings) > 0L) {
    "COMPLETED_WITH_WARNINGS"
  } else {
    "PASS"
  }

  status_row <- data.frame(
    Comparison = comparison,
    Run_status = run_status,
    Fixed_effect_formula = fix_formula,
    Random_effect_formula = "NULL",
    Multiple_testing_primary = "BH FDR within contrast",
    Multiple_testing_sensitivity = "Holm within contrast",
    Pseudocount_sensitivity_required = TRUE,
    Additional_prevalence_filter_applied = FALSE,
    ANCOMBC2_prv_cut = ancom_prv_cut,
    N_patient_tissue_units = nsamples(pb_phy),
    N_patients = n_distinct(pb_metadata$Patient),
    N_genera_entering_contrast = ntaxa(pb_phy),
    N_model_output_genera = nrow(model_results),
    N_model_tested_genera = sum(model_results$Tissue_test_estimable),
    N_model_based_BH_supported = nrow(model_supported_BH),
    N_model_based_Holm_supported = nrow(model_supported_Holm),
    N_structural_zero_findings = nrow(structural_zero_results),

    Minimum_model_P = ifelse(
      all(is.na(model_results$P_value)),
      NA_real_,
      min(model_results$P_value, na.rm = TRUE)
    ),

    Minimum_model_BH_FDR = ifelse(
      all(is.na(model_results$BH_FDR)),
      NA_real_,
      min(model_results$BH_FDR, na.rm = TRUE)
    ),

    Minimum_model_Holm_adjusted_P = ifelse(
      all(is.na(model_results$Holm_adjusted_P_sensitivity)),
      NA_real_,
      min(model_results$Holm_adjusted_P_sensitivity, na.rm = TRUE)
    ),

    N_warnings = length(captured_warnings),
    N_critical_warnings = length(critical_warnings),
    N_noncritical_warnings = length(noncritical_warnings),
    N_ANCOMBC_flag_BH_mismatches = n_flag_mismatches,
    Invalid_probability_values = invalid_probability_values,
    Supported_rows_finite = supported_rows_finite,
    Degenerate_all_P_equal_1_with_nonzero_W = degenerate_all_p_one,
    Degenerate_all_P_missing = degenerate_all_p_missing,
    stringsAsFactors = FALSE
  )

  write.csv(
    status_row,
    file.path(
      comparison_dir,
      paste0("ANCOMBC2_", tag, "_run_status.csv")
    ),
    row.names = FALSE
  )

  writeLines(
    if (length(critical_warnings) == 0L) {
      "No critical warnings were emitted."
    } else {
      critical_warnings
    },
    file.path(
      comparison_dir,
      paste0("ANCOMBC2_", tag, "_critical_warnings.txt")
    )
  )

  writeLines(
    if (length(noncritical_warnings) == 0L) {
      "No non-critical warnings were emitted."
    } else {
      noncritical_warnings
    },
    file.path(
      comparison_dir,
      paste0("ANCOMBC2_", tag, "_noncritical_warnings.txt")
    )
  )

  if (run_status == "REVIEW_REQUIRED") {
    warning(
      comparison,
      " requires review. Do not report this contrast until its status and ",
      "warning files have been inspected."
    )
  }

  list(
    comparison = comparison,
    status = status_row,
    model_all = model_results,
    model_supported_BH = model_supported_BH,
    model_supported_Holm = model_supported_Holm,
    structural_zero_findings = structural_zero_results,
    manifest = pb_metadata
  )
}

# Primary comparison specifications.
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
    "At least one main ANCOM-BC2 contrast was configured as paired. ",
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
      Fixed_effect_formula = if (x$group_2 == "H_ADJ_SSL_HP") "Patient + Group" else "Group",
      Random_effect_formula = "NULL",
      stringsAsFactors = FALSE
    )
  }
))

write.csv(
  main_contrast_design_audit,
  file.path(out_dir, "MAIN_CONTRAST_DESIGN_AUDIT.csv"),
  row.names = FALSE
)

contrast_results <- lapply(
  comparison_specifications,
  function(specification) {
    run_ancombc2_contrast(
      group_1 = specification$group_1,
      group_2 = specification$group_2,
      paired_only = specification$paired_only
    )
  }
)

names(contrast_results) <- vapply(
  comparison_specifications,
  function(x) paste(x$group_1, "-", x$group_2),
  character(1)
)

run_summary <- bind_rows(
  lapply(contrast_results, function(x) x$status)
) %>%
  mutate(
    ANCOMBC_version = as.character(packageVersion("ANCOMBC")),
    Bias_EM_patch_applied_in_memory =
      ancombc_patch_status$Bias_EM_patch_applied_in_memory[1],
    Bias_EM_patch_status =
      ancombc_patch_status$Patch_status[1]
  )

all_model_results <- bind_rows(
  lapply(contrast_results, function(x) x$model_all)
)

all_model_supported_BH <- bind_rows(
  lapply(contrast_results, function(x) x$model_supported_BH)
)

all_model_supported_Holm <- bind_rows(
  lapply(contrast_results, function(x) x$model_supported_Holm)
)

all_structural_zero_findings <- bind_rows(
  lapply(contrast_results, function(x) x$structural_zero_findings)
)

model_based_thesis_table <- all_model_supported_BH %>%
  select(
    taxon,
    Comparison,
    First_group,
    Second_group,
    Patient_fixed_effect,
    Direction,
    LFC_natural_log,
    LFC_log2,
    CI95_low_natural_log,
    CI95_high_natural_log,
    CI95_low_log2,
    CI95_high_log2,
    SE_natural_log,
    W,
    P_value,
    BH_FDR,
    Holm_adjusted_P_sensitivity,
    passed_pseudocount_sensitivity,
    starts_with("N_patient_tissue_units_"),
    starts_with("N_detected_"),
    starts_with("Prevalence_")
  )

structural_zero_thesis_table <- all_structural_zero_findings %>%
  select(
    taxon,
    Comparison,
    First_group,
    Second_group,
    Patient_fixed_effect,
    Evidence_type,
    Direction,
    Structural_zero_in_group,
    N_units_in_structural_zero_group,
    N_detected_in_structural_zero_group,
    Prevalence_in_structural_zero_group,
    Structural_zero_interpretation,
    starts_with("N_patient_tissue_units_"),
    starts_with("N_detected_"),
    starts_with("Prevalence_")
  )

write.csv(
  run_summary,
  file.path(out_dir, "ANCOMBC2_ALL_CONTRASTS_RUN_SUMMARY.csv"),
  row.names = FALSE
)

write.csv(
  all_model_results,
  file.path(out_dir, "ANCOMBC2_ALL_CONTRASTS_MODEL_BASED_ALL.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  all_model_supported_BH,
  file.path(
    out_dir,
    "ANCOMBC2_ALL_CONTRASTS_MODEL_BASED_BH_FDR_SUPPORTED.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  all_model_supported_Holm,
  file.path(
    out_dir,
    "ANCOMBC2_ALL_CONTRASTS_MODEL_BASED_HOLM_SENSITIVITY.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  all_structural_zero_findings,
  file.path(
    out_dir,
    "ANCOMBC2_ALL_CONTRASTS_STRUCTURAL_ZERO_FINDINGS.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  model_based_thesis_table,
  file.path(
    out_dir,
    "Supplementary_Table_ANCOMBC2_MODEL_BASED_BH_FDR.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  structural_zero_thesis_table,
  file.path(
    out_dir,
    "Supplementary_Table_ANCOMBC2_STRUCTURAL_ZEROS.csv"
  ),
  row.names = FALSE,
  na = ""
)

readme_lines <- c(
  "FINAL CHAPTER 4.5.3 ANCOM-BC2 GENUS-LEVEL PSEUDOBULK ANALYSIS",
  "================================================================",
  "",
  paste0("Input: ", physeq_path),
  paste0("ANCOMBC version: ", packageVersion("ANCOMBC")),
  paste0(
    "ANCOMBC .bias_em patch applied in memory: ",
    ancombc_patch_status$Bias_EM_patch_applied_in_memory[1]
  ),
  paste0(
    "Patch status: ",
    ancombc_patch_status$Patch_status[1]
  ),
  paste0("Output directory: ", out_dir),
  "",
  "Taxonomic-level decision:",
  "  The final analysis is genus-level.",
  "  ASV-level ANCOM-BC2 is not used as the primary sensitivity analysis.",
  "  Genus level matches the primary DESeq2 analysis and reduces the number",
  "  of sparse, difficult-to-interpret tests.",
  "",
  "Input filtering:",
  "  The input phyloseq object was already filtered at 5% ASV prevalence.",
  "  No further prevalence or abundance threshold was applied.",
  "  ANCOM-BC2 prv_cut = 0.",
  "  Features with zero counts across every unit in a specific contrast were",
  "  removed because they contain no information for that contrast.",
  "",
  "Genus preparation:",
  "  Lineage-aware genus labels were created in memory.",
  "  Counts were agglomerated to genus before pseudobulk construction.",
  "  tax_level = NULL in ANCOM-BC2 because every input feature was already a",
  "  genus-level unit.",
  "",
  "Repeated-biopsy handling:",
  "  Counts were summed within Patient x Group.",
  "  No patient-matching restriction is applied in any of the three main contrasts.",
  "",
  "Models:",
  "  SSL_HP vs H_COL: Group; rand_formula = NULL.",
  "  H_ADJ_SSL_HP vs H_COL: Group; rand_formula = NULL.",
  "  SSL_HP vs H_ADJ_SSL_HP: Patient + Group using all available Patient x Group profiles;",
  "  rand_formula = NULL.",
  "",
  "ANCOM-BC2 testing:",
  "  Three separate prespecified two-group analyses.",
  "  The single Group coefficient was extracted from fit$res.",
  "  global = FALSE and pairwise = FALSE.",
  "  Primary multiplicity correction: BH FDR across genera within each",
  "  prespecified contrast.",
  "  Stricter sensitivity correction: Holm across genera within each contrast.",
  "  No additional correction was applied across the three prespecified",
  "  biological contrasts.",
  "  Pseudocount sensitivity was enabled and is required for a reported",
  "  model-based finding.",
  "",
  "Structural zeros:",
  "  Structural-zero findings are exported separately from finite model-based",
  "  estimates.",
  "  They do not receive a finite fold change, P value or adjusted P value.",
  "  Structural-zero classification is a group-level screen in all three",
  "  main contrasts; no patient-matching restriction is applied.",
  "",
  "Reporting rule:",
  "  Primary model-based result = BH FDR < 0.05, ANCOM-BC2 diff_abn = TRUE,",
  "  and passed_ss = TRUE.",
  "  The Holm table is a stricter sensitivity result.",
  "  Inspect warning files and coefficient-level availability before reporting.",
  "  Retain REVIEW_REQUIRED in the raw archive; record any subsequent review separately.",
  "  Report only available tissue tests supported by the documented review."
)

writeLines(
  readme_lines,
  file.path(out_dir, "README_ANCOMBC2_GENUS_PSEUDOBULK_FINAL.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "sessionInfo_ANCOMBC2_GENUS_PSEUDOBULK_FINAL.txt")
)

cat("\n===== ANCOM-BC2 FINAL RUN SUMMARY =====\n")
print(run_summary)

cat("\nFinal taxonomic level: genus.\n")
cat("No additional prevalence filter was applied: prv_cut = 0.\n")
cat("Primary multiple-testing correction: BH FDR within each contrast.\n")
cat("Holm-adjusted P values were exported as a stricter sensitivity analysis.\n")
cat("Results saved to:\n", out_dir, "\n")

if (any(run_summary$Run_status == "REVIEW_REQUIRED")) {
  stop(
    "At least one contrast requires review. Inspect ",
    "ANCOMBC2_ALL_CONTRASTS_RUN_SUMMARY.csv and the contrast warning files. ",
    "Inspect the warnings and record a review before reporting flagged results.",
    call. = FALSE
  )
}

cat("\nANCOM-BC2 ANALYSIS COMPLETE\n")
