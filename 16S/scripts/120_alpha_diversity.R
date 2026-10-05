#!/usr/bin/env Rscript
# Alpha diversity before and after prevalence filtering. Run from the 16S directory.

set.seed(1)

options(stringsAsFactors = FALSE)
options(contrasts = c("contr.treatment", "contr.poly"))

required_packages <- c(
  "phyloseq",
  "dplyr",
  "tibble",
  "ggplot2",
  "lme4",
  "lmerTest"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "The following required package(s) are missing: ",
    paste(missing_packages, collapse = ", "),
    "."
  )
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
})

if (capabilities("cairo")) {
  options(bitmapType = "cairo")
}

# Settings

base_dir <- "outputs"

physeq_unfiltered_path <- file.path(
  base_dir,
  "phyloseq_March26_post_rarefaction_without7699.rds"
)

physeq_5pct_path <- file.path(
  base_dir,
  "phyloseq_final_March26_5pct_without7699.rds"
)

out_dir <- file.path(
  base_dir,
  "Alpha_diversity_FINAL_BOTH_NO_EMMEANS_without7699",
  format(Sys.time(), "%Y%m%d_%H%M%S")
)

if (dir.exists(out_dir)) stop("Output directory already exists: ", out_dir)
if (!dir.create(out_dir, recursive = TRUE)) stop("Could not create: ", out_dir)

patient_col_manual <- NULL

groups_keep <- c(
  "SSL_HP",
  "H_ADJ_SSL_HP",
  "H_COL"
)

metrics <- c(
  "Observed",
  "Shannon",
  "Simpson"
)

dataset_levels <- c(
  "No prevalence filter",
  "5% prevalence filter"
)

dataset_stubs <- c(
  "No prevalence filter" = "NO_PREVALENCE_FILTER",
  "5% prevalence filter" = "5PCT_PREVALENCE_FILTER"
)

comparison_order <- c(
  "SSL_HP vs H_ADJ_SSL_HP",
  "SSL_HP vs H_COL",
  "H_ADJ_SSL_HP vs H_COL"
)

comparison_labels <- c(
  "SSL_HP vs H_ADJ_SSL_HP" = "SSL/HP versus adjacent",
  "SSL_HP vs H_COL" = "SSL/HP versus healthy colon",
  "H_ADJ_SSL_HP vs H_COL" = "Adjacent versus healthy colon"
)

group_labels <- c(
  "SSL_HP" = "SSL/HP",
  "H_ADJ_SSL_HP" = "Adjacent",
  "H_COL" = "Healthy colon"
)

metric_titles <- c(
  "Observed" = "Observed ASV richness",
  "Shannon" = "Shannon diversity",
  "Simpson" = "Simpson diversity"
)

metric_y_labels <- c(
  "Observed" = "Observed ASVs",
  "Shannon" = "Shannon index",
  "Simpson" = "Simpson index"
)

paper_colours <- c(
  "SSL_HP" = "#FF0000",
  "H_ADJ_SSL_HP" = "#E9C61D",
  "H_COL" = "#1313EC"
)

font_family <- "sans"
export_dpi <- 600

primary_plot_p_type <- "holm"

if (!primary_plot_p_type %in% c("holm", "raw")) {
  stop("primary_plot_p_type must be either 'holm' or 'raw'.")
}

# General helper functions

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

normalise_ids <- function(x) {
  x <- trimws(as.character(x))
  gsub("\\.", "-", x)
}

first_existing_column <- function(column_names, candidates) {
  hit <- candidates[candidates %in% column_names]
  if (length(hit) == 0L) return(NULL)
  hit[1]
}

extract_last_numeric_token <- function(x) {
  token_list <- regmatches(
    as.character(x),
    gregexpr("[0-9]+", as.character(x), perl = TRUE)
  )

  vapply(
    token_list,
    function(z) {
      if (length(z) == 0L) return(NA_character_)
      tail(z, 1L)
    },
    character(1)
  )
}

otu_samples_by_taxa <- function(ps) {
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)
  otu
}

check_count_phyloseq <- function(ps, object_name) {
  if (!inherits(ps, "phyloseq")) {
    stop(object_name, " is not a phyloseq object.")
  }

  otu <- otu_samples_by_taxa(ps)

  if (
    any(!is.finite(otu)) ||
    any(otu < 0) ||
    any(abs(otu - round(otu)) > 1e-8)
  ) {
    stop(
      object_name,
      " must contain non-negative integer ASV counts."
    )
  }

  if (any(sample_sums(ps) <= 0)) {
    stop(object_name, " contains at least one zero-depth sample.")
  }

  invisible(TRUE)
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

p_significance <- function(p) {
  if (length(p) == 0L || is.na(p)) return(NA_character_)
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

safe_lmer_ml <- function(formula_object, data) {
  tryCatch(
    lmerTest::lmer(
      formula_object,
      data = data,
      REML = FALSE,
      na.action = na.fail,
      control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 200000)
      )
    ),
    error = function(e) {
      structure(
        list(message = conditionMessage(e)),
        class = "fit_error"
      )
    }
  )
}

likelihood_ratio_p <- function(reduced_model, full_model) {
  if (
    inherits(reduced_model, "fit_error") ||
    inherits(full_model, "fit_error")
  ) {
    return(NA_real_)
  }

  comparison <- tryCatch(
    anova(reduced_model, full_model),
    error = function(e) NULL
  )

  if (is.null(comparison)) return(NA_real_)

  p_column <- grep(
    "^Pr\\(",
    colnames(comparison),
    value = TRUE
  )

  if (length(p_column) == 0L) return(NA_real_)

  as.numeric(comparison[nrow(comparison), p_column[1]])
}

extract_fixed_effect <- function(fit, coefficient_name) {
  coefficient_table <- as.data.frame(
    summary(fit)$coefficients,
    check.names = FALSE
  )

  if (!coefficient_name %in% rownames(coefficient_table)) {
    stop(
      "Coefficient '",
      coefficient_name,
      "' was not found. Available coefficients: ",
      paste(rownames(coefficient_table), collapse = ", ")
    )
  }

  p_column <- grep(
    "^Pr\\(",
    colnames(coefficient_table),
    value = TRUE
  )

  if (length(p_column) == 0L) {
    stop("No fixed-effect P-value column was found in the model summary.")
  }

  estimate <- as.numeric(
    coefficient_table[coefficient_name, "Estimate"]
  )

  standard_error <- as.numeric(
    coefficient_table[coefficient_name, "Std. Error"]
  )

  degrees_freedom <- if ("df" %in% colnames(coefficient_table)) {
    as.numeric(coefficient_table[coefficient_name, "df"])
  } else {
    NA_real_
  }

  critical_value <- if (is.finite(degrees_freedom)) {
    qt(0.975, df = degrees_freedom)
  } else {
    qnorm(0.975)
  }

  c(
    Estimate = estimate,
    Std_Error = standard_error,
    DF = degrees_freedom,
    CI_low = estimate - critical_value * standard_error,
    CI_high = estimate + critical_value * standard_error,
    P_raw = as.numeric(
      coefficient_table[coefficient_name, p_column[1]]
    )
  )
}

get_convergence_message <- function(fit) {
  messages <- fit@optinfo$conv$lme4$messages

  if (is.null(messages) || length(messages) == 0L) {
    return("")
  }

  paste(messages, collapse = " | ")
}

get_variance_components <- function(fit) {
  variance_table <- as.data.frame(lme4::VarCorr(fit))

  patient_variance <- variance_table$vcov[
    variance_table$grp == "Patient" & is.na(variance_table$var2)
  ]

  if (length(patient_variance) == 0L) {
    patient_variance <- NA_real_
  } else {
    patient_variance <- patient_variance[1]
  }

  c(
    Patient_random_intercept_variance = patient_variance,
    Residual_variance = sigma(fit)^2
  )
}

save_tiff_safely <- function(
    plot,
    filename,
    width,
    height,
    dpi = 600
) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(
      filename = filename,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      device = ragg::agg_tiff,
      compression = "lzw",
      background = "white"
    )

    return(invisible(TRUE))
  }

  if (capabilities("cairo")) {
    ggplot2::ggsave(
      filename = filename,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      device = "tiff",
      type = "cairo",
      compression = "lzw",
      bg = "white"
    )

    return(invisible(TRUE))
  }

  warning(
    "TIFF was not exported because neither ragg nor Cairo is available. ",
    "The PDF and 600-dpi PNG files were still created."
  )

  invisible(FALSE)
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

save_plot_grid_pdf <- function(
    plot_list,
    filename,
    nrow,
    ncol,
    width,
    height
) {
  grDevices::pdf(
    filename,
    width = width,
    height = height,
    onefile = TRUE,
    bg = "white"
  )

  on.exit(grDevices::dev.off(), add = TRUE)
  draw_plot_grid(plot_list, nrow, ncol)
  invisible(TRUE)
}

save_plot_grid_png <- function(
    plot_list,
    filename,
    nrow,
    ncol,
    width,
    height,
    dpi = 600
) {
  grDevices::png(
    filename,
    width = width,
    height = height,
    units = "in",
    res = dpi,
    bg = "white"
  )

  on.exit(grDevices::dev.off(), add = TRUE)
  draw_plot_grid(plot_list, nrow, ncol)
  invisible(TRUE)
}

# Load and validate both phyloseq objects

if (!file.exists(physeq_unfiltered_path)) {
  stop(
    "The non-prevalence-filtered phyloseq object was not found: ",
    physeq_unfiltered_path
  )
}

if (!file.exists(physeq_5pct_path)) {
  stop(
    "The 5% prevalence-filtered phyloseq object was not found: ",
    physeq_5pct_path
  )
}

physeq_unfiltered_full <- readRDS(physeq_unfiltered_path)
physeq_5pct_full <- readRDS(physeq_5pct_path)

check_count_phyloseq(
  physeq_unfiltered_full,
  "Non-prevalence-filtered object"
)

check_count_phyloseq(
  physeq_5pct_full,
  "5% prevalence-filtered object"
)

for (object_name in c("physeq_unfiltered_full", "physeq_5pct_full")) {
  ps <- get(object_name)

  if (any(grepl("7699", sample_names(ps), fixed = TRUE))) {
    stop(
      object_name,
      " still contains a sample whose name includes 7699."
    )
  }
}

if (!setequal(
  sample_names(physeq_unfiltered_full),
  sample_names(physeq_5pct_full)
)) {
  only_unfiltered <- setdiff(
    sample_names(physeq_unfiltered_full),
    sample_names(physeq_5pct_full)
  )

  only_filtered <- setdiff(
    sample_names(physeq_5pct_full),
    sample_names(physeq_unfiltered_full)
  )

  write.csv(
    data.frame(SampleID = only_unfiltered),
    file.path(out_dir, "Samples_only_in_unfiltered_object.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(SampleID = only_filtered),
    file.path(out_dir, "Samples_only_in_5pct_object.csv"),
    row.names = FALSE
  )

  stop(
    "The two phyloseq objects do not contain the same samples. ",
    "Inspect the two sample-difference CSV files."
  )
}

otu_unfiltered_full <- otu_samples_by_taxa(physeq_unfiltered_full)
otu_5pct_full <- otu_samples_by_taxa(physeq_5pct_full)

if (!all(colnames(otu_5pct_full) %in% colnames(otu_unfiltered_full))) {
  stop("The 5% object contains ASVs absent from the unfiltered object.")
}

otu_unfiltered_overlap <- otu_unfiltered_full[
  rownames(otu_5pct_full),
  colnames(otu_5pct_full),
  drop = FALSE
]

if (!isTRUE(all.equal(
  unname(otu_unfiltered_overlap),
  unname(otu_5pct_full),
  check.attributes = FALSE
))) {
  stop(
    "Counts for retained ASVs differ between the unfiltered and 5% objects."
  )
}

prevalence_threshold <- ceiling(
  0.05 * nsamples(physeq_unfiltered_full)
)

unfiltered_prevalence <- colSums(otu_unfiltered_full > 0)
expected_5pct_taxa <- names(
  unfiltered_prevalence[
    unfiltered_prevalence >= prevalence_threshold
  ]
)

if (!setequal(expected_5pct_taxa, taxa_names(physeq_5pct_full))) {
  missing_expected <- setdiff(
    expected_5pct_taxa,
    taxa_names(physeq_5pct_full)
  )

  unexpected_retained <- setdiff(
    taxa_names(physeq_5pct_full),
    expected_5pct_taxa
  )

  write.csv(
    data.frame(ASV = missing_expected),
    file.path(out_dir, "Expected_5pct_ASVs_missing_from_filtered_object.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(ASV = unexpected_retained),
    file.path(out_dir, "Unexpected_ASVs_in_5pct_filtered_object.csv"),
    row.names = FALSE
  )

  stop(
    "The filtered object is not the exact 5% prevalence derivative of the ",
    "unfiltered object. Inspect the ASV-difference CSV files."
  )
}

cat("\n===== FULL OBJECT VALIDATION PASSED =====\n")
cat("Unfiltered samples:", nsamples(physeq_unfiltered_full), "\n")
cat("Unfiltered ASVs:", ntaxa(physeq_unfiltered_full), "\n")
cat("5% filtered samples:", nsamples(physeq_5pct_full), "\n")
cat("5% filtered ASVs:", ntaxa(physeq_5pct_full), "\n")
cat("5% prevalence threshold:", prevalence_threshold, "samples\n")

# Create the three final analysis groups in both objects

add_analysis_group <- function(ps) {
  if (!"Sample_Type" %in% colnames(sample_data(ps))) {
    stop("Sample_Type was not found in sample_data().")
  }

  metadata <- data.frame(
    sample_data(ps),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  group <- clean_character(metadata$Sample_Type)

  group[group %in% c("SSL", "HP", "SSL_HP")] <- "SSL_HP"

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

  group[!group %in% groups_keep] <- NA_character_

  sample_data(ps)$Group <- factor(
    group,
    levels = groups_keep
  )

  ps <- subset_samples(ps, !is.na(Group))
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)

  sample_data(ps)$Group <- factor(
    as.character(sample_data(ps)$Group),
    levels = groups_keep
  )

  ps
}

physeq_unfiltered <- add_analysis_group(physeq_unfiltered_full)
physeq_5pct <- add_analysis_group(physeq_5pct_full)

if (!setequal(
  sample_names(physeq_unfiltered),
  sample_names(physeq_5pct)
)) {
  stop(
    "The two objects no longer contain the same samples after group subsetting."
  )
}

cat("\n===== FINAL THREE-GROUP SAMPLE COUNTS =====\n")
print(table(sample_data(physeq_unfiltered)$Group))

# Identify patient IDs from the unfiltered object

metadata_master <- data.frame(
  sample_data(physeq_unfiltered),
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  rownames_to_column("SampleID") %>%
  mutate(
    SampleID_norm = normalise_ids(SampleID)
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
  if (!patient_col_manual %in% colnames(metadata_master)) {
    stop(
      "patient_col_manual was not found in the metadata: ",
      patient_col_manual
    )
  }

  patient_column <- patient_col_manual
} else {
  patient_column <- first_existing_column(
    colnames(metadata_master),
    patient_candidates
  )
}

if (!is.null(patient_column)) {
  patient <- clean_character(
    metadata_master[[patient_column]]
  )

  cat(
    "\nPatient IDs taken from metadata column: ",
    patient_column,
    "\n",
    sep = ""
  )
} else {
  patient <- extract_last_numeric_token(
    metadata_master$SampleID
  )

  cat(
    "\nNo patient metadata column was detected. Patient IDs were extracted ",
    "from the final numeric token in each SampleID.\n",
    sep = ""
  )
}

if (any(is.na(patient))) {
  write.csv(
    metadata_master[is.na(patient), , drop = FALSE],
    file.path(out_dir, "Samples_with_missing_patient_ID.csv"),
    row.names = FALSE
  )

  stop(
    "At least one Patient ID is missing. Inspect ",
    "Samples_with_missing_patient_ID.csv and set patient_col_manual if needed."
  )
}

patient <- as.character(patient)

if (any(patient == "7699")) {
  stop("Patient 7699 is still present according to the Patient identifiers.")
}

metadata_model <- data.frame(
  SampleID = metadata_master$SampleID,
  SampleID_norm = metadata_master$SampleID_norm,
  Group = factor(
    metadata_master$Group,
    levels = groups_keep
  ),
  Patient = factor(patient),
  stringsAsFactors = FALSE
)

if (anyDuplicated(metadata_model$SampleID_norm)) {
  stop("Duplicated normalised sample IDs were found in the metadata.")
}

if (!any(table(metadata_model$Patient) > 1L)) {
  stop(
    "No repeated patients were detected. Check the Patient identifiers before ",
    "using a mixed-effects model."
  )
}

reads_common <- data.frame(
  SampleID = sample_names(physeq_unfiltered),
  SampleID_norm = normalise_ids(sample_names(physeq_unfiltered)),
  Reads = as.numeric(sample_sums(physeq_unfiltered)),
  stringsAsFactors = FALSE
)

if (anyDuplicated(reads_common$SampleID_norm)) {
  stop("Duplicated normalised sample IDs were found in the read-depth table.")
}

# Calculate alpha diversity for both objects

make_analysis_table <- function(ps, dataset_label) {
  alpha_table <- estimate_richness(
    ps,
    measures = metrics
  ) %>%
    rownames_to_column("SampleID_alpha") %>%
    mutate(
      SampleID_norm = normalise_ids(SampleID_alpha)
    )

  if (anyDuplicated(alpha_table$SampleID_norm)) {
    stop(
      "Duplicated normalised sample IDs were found in the alpha table for ",
      dataset_label,
      "."
    )
  }

  dataset_reads <- data.frame(
    SampleID_norm = normalise_ids(sample_names(ps)),
    Reads_after_preprocessing = as.numeric(sample_sums(ps)),
    stringsAsFactors = FALSE
  )

  analysis_table <- alpha_table %>%
    left_join(
      reads_common,
      by = "SampleID_norm"
    ) %>%
    left_join(
      dataset_reads,
      by = "SampleID_norm"
    ) %>%
    left_join(
      metadata_model,
      by = "SampleID_norm",
      suffix = c("", "_metadata")
    ) %>%
    transmute(
      Dataset = dataset_label,
      SampleID,
      SampleID_alpha,
      SampleID_norm,
      Patient,
      Group,
      Reads,
      Reads_after_preprocessing,
      Observed,
      Shannon,
      Simpson
    )

  required_columns <- c(
    "SampleID",
    "Patient",
    "Group",
    "Reads",
    "Reads_after_preprocessing",
    metrics
  )

  missing_required <- !complete.cases(
    analysis_table[, required_columns, drop = FALSE]
  )

  if (any(missing_required)) {
    missing_filename <- paste0(
      "Samples_with_missing_required_data_",
      dataset_stubs[[dataset_label]],
      ".csv"
    )

    write.csv(
      analysis_table[missing_required, , drop = FALSE],
      file.path(out_dir, missing_filename),
      row.names = FALSE
    )

    stop(
      "Missing required values were found in ",
      dataset_label,
      ". Inspect ",
      missing_filename,
      "."
    )
  }

  if (nrow(analysis_table) != nsamples(ps)) {
    stop(
      "The alpha table for ",
      dataset_label,
      " contains ",
      nrow(analysis_table),
      " rows but its phyloseq object contains ",
      nsamples(ps),
      " samples."
    )
  }

  analysis_table %>%
    mutate(
      Dataset = factor(Dataset, levels = dataset_levels),
      Patient = droplevels(factor(Patient)),
      Group = factor(Group, levels = groups_keep)
    )
}

analysis_unfiltered <- make_analysis_table(
  physeq_unfiltered,
  "No prevalence filter"
)

analysis_5pct <- make_analysis_table(
  physeq_5pct,
  "5% prevalence filter"
)

analysis_df <- bind_rows(
  analysis_unfiltered,
  analysis_5pct
) %>%
  mutate(
    Dataset = factor(Dataset, levels = dataset_levels),
    Group = factor(Group, levels = groups_keep)
  )

write.csv(
  analysis_df,
  file.path(out_dir, "Alpha_FINAL_analysis_input_BOTH.csv"),
  row.names = FALSE
)

write.csv(
  metadata_model,
  file.path(out_dir, "Alpha_FINAL_sample_patient_group_alignment.csv"),
  row.names = FALSE
)

preprocessing_comparison <- analysis_unfiltered %>%
  select(
    SampleID,
    Patient,
    Group,
    Reads,
    Observed_unfiltered = Observed,
    Shannon_unfiltered = Shannon,
    Simpson_unfiltered = Simpson
  ) %>%
  left_join(
    analysis_5pct %>%
      select(
        SampleID,
        Reads_after_5pct = Reads_after_preprocessing,
        Observed_5pct = Observed,
        Shannon_5pct = Shannon,
        Simpson_5pct = Simpson
      ),
    by = "SampleID"
  ) %>%
  mutate(
    Observed_change_5pct_minus_unfiltered =
      Observed_5pct - Observed_unfiltered,
    Shannon_change_5pct_minus_unfiltered =
      Shannon_5pct - Shannon_unfiltered,
    Simpson_change_5pct_minus_unfiltered =
      Simpson_5pct - Simpson_unfiltered
  )

write.csv(
  preprocessing_comparison,
  file.path(out_dir, "Alpha_FINAL_per_sample_preprocessing_comparison.csv"),
  row.names = FALSE
)

cat("\n===== FINAL ANALYSIS SAMPLE AND PATIENT COUNTS =====\n")
print(
  analysis_df %>%
    group_by(Dataset, Group) %>%
    summarise(
      N_samples = n(),
      N_patients = n_distinct(Patient),
      .groups = "drop"
    )
)

# Descriptive summaries

descriptive_summary <- analysis_df %>%
  group_by(Dataset, Group) %>%
  summarise(
    N_samples = n(),
    N_patients = n_distinct(Patient),

    Reads_median = median(Reads),
    Reads_Q1 = quantile(Reads, 0.25),
    Reads_Q3 = quantile(Reads, 0.75),

    Reads_after_preprocessing_median = median(Reads_after_preprocessing),
    Reads_after_preprocessing_Q1 = quantile(Reads_after_preprocessing, 0.25),
    Reads_after_preprocessing_Q3 = quantile(Reads_after_preprocessing, 0.75),

    Observed_median = median(Observed),
    Observed_Q1 = quantile(Observed, 0.25),
    Observed_Q3 = quantile(Observed, 0.75),

    Shannon_median = median(Shannon),
    Shannon_Q1 = quantile(Shannon, 0.25),
    Shannon_Q3 = quantile(Shannon, 0.75),

    Simpson_median = median(Simpson),
    Simpson_Q1 = quantile(Simpson, 0.25),
    Simpson_Q3 = quantile(Simpson, 0.75),

    .groups = "drop"
  )

write.csv(
  descriptive_summary,
  file.path(out_dir, "Alpha_FINAL_descriptive_summary_BOTH.csv"),
  row.names = FALSE
)

# Fit the mixed models and all three pairwise comparisons

fit_metric_reference <- function(
    model_data,
    metric,
    reference_group
) {
  model_data <- model_data %>%
    mutate(
      Group = relevel(
        factor(Group, levels = groups_keep),
        ref = reference_group
      ),
      Patient = droplevels(factor(Patient))
    )

  full_formula <- as.formula(
    paste0(metric, " ~ Group + Reads + (1|Patient)")
  )

  reduced_group_formula <- as.formula(
    paste0(metric, " ~ Reads + (1|Patient)")
  )

  reduced_reads_formula <- as.formula(
    paste0(metric, " ~ Group + (1|Patient)")
  )

  full_fit <- safe_lmer_ml(
    full_formula,
    model_data
  )

  if (inherits(full_fit, "fit_error")) {
    stop(
      metric,
      " model with reference '",
      reference_group,
      "' failed: ",
      full_fit$message
    )
  }

  reduced_group_fit <- safe_lmer_ml(
    reduced_group_formula,
    model_data
  )

  reduced_reads_fit <- safe_lmer_ml(
    reduced_reads_formula,
    model_data
  )

  list(
    fit = full_fit,
    data = model_data,
    formula = paste(deparse(full_formula), collapse = ""),
    group_p = likelihood_ratio_p(
      reduced_group_fit,
      full_fit
    ),
    reads_p = likelihood_ratio_p(
      reduced_reads_fit,
      full_fit
    )
  )
}

model_objects <- list()
overall_results_list <- list()
pairwise_results_list <- list()

for (dataset_label in dataset_levels) {
  dataset_data <- analysis_df %>%
    filter(Dataset == dataset_label) %>%
    droplevels()

  for (metric in metrics) {
    model_hcol <- fit_metric_reference(
      model_data = dataset_data,
      metric = metric,
      reference_group = "H_COL"
    )

    model_hadj <- fit_metric_reference(
      model_data = dataset_data,
      metric = metric,
      reference_group = "H_ADJ_SSL_HP"
    )

    model_key_stub <- paste0(
      dataset_stubs[[dataset_label]],
      "_",
      metric
    )

    model_objects[[
      paste0(model_key_stub, "_H_COL_reference")
    ]] <- model_hcol$fit

    model_objects[[
      paste0(model_key_stub, "_H_ADJ_reference")
    ]] <- model_hadj$fit

    singular_model <- lme4::isSingular(
      model_hcol$fit,
      tol = 1e-4
    )

    convergence_message <- get_convergence_message(
      model_hcol$fit
    )

    variance_components <- get_variance_components(
      model_hcol$fit
    )

    if (singular_model) {
      warning(
        dataset_label,
        " — ",
        metric,
        " model is singular. Inspect the variance components and diagnostics."
      )
    }

    if (nzchar(convergence_message)) {
      warning(
        dataset_label,
        " — ",
        metric,
        " convergence message: ",
        convergence_message
      )
    }

    ssl_vs_adj <- extract_fixed_effect(
      model_hadj$fit,
      "GroupSSL_HP"
    )

    ssl_vs_hcol <- extract_fixed_effect(
      model_hcol$fit,
      "GroupSSL_HP"
    )

    adj_vs_hcol <- extract_fixed_effect(
      model_hcol$fit,
      "GroupH_ADJ_SSL_HP"
    )

    pairwise_metric <- data.frame(
      Dataset = dataset_label,
      Metric = metric,
      Comparison = comparison_order,
      Comparison_label = unname(
        comparison_labels[comparison_order]
      ),

      Estimate_first_minus_second = c(
        ssl_vs_adj["Estimate"],
        ssl_vs_hcol["Estimate"],
        adj_vs_hcol["Estimate"]
      ),

      Std_Error = c(
        ssl_vs_adj["Std_Error"],
        ssl_vs_hcol["Std_Error"],
        adj_vs_hcol["Std_Error"]
      ),

      DF = c(
        ssl_vs_adj["DF"],
        ssl_vs_hcol["DF"],
        adj_vs_hcol["DF"]
      ),

      CI_low = c(
        ssl_vs_adj["CI_low"],
        ssl_vs_hcol["CI_low"],
        adj_vs_hcol["CI_low"]
      ),

      CI_high = c(
        ssl_vs_adj["CI_high"],
        ssl_vs_hcol["CI_high"],
        adj_vs_hcol["CI_high"]
      ),

      P_raw = c(
        ssl_vs_adj["P_raw"],
        ssl_vs_hcol["P_raw"],
        adj_vs_hcol["P_raw"]
      ),

      stringsAsFactors = FALSE
    )

    pairwise_metric$P_Holm <- p.adjust(
      pairwise_metric$P_raw,
      method = "holm"
    )

    pairwise_metric$Significance_Holm <- vapply(
      pairwise_metric$P_Holm,
      p_significance,
      character(1)
    )

    pairwise_metric$Significant_Holm_0_05 <-
      !is.na(pairwise_metric$P_Holm) &
      pairwise_metric$P_Holm < 0.05

    pairwise_results_list[[model_key_stub]] <- pairwise_metric

    overall_results_list[[model_key_stub]] <- data.frame(
      Dataset = dataset_label,
      Metric = metric,
      Model = model_hcol$formula,
      Read_depth_covariate = paste0(
        "Post-rarefaction read counts from ",
        "the non-prevalence-filtered object"
      ),
      Estimation = "Maximum likelihood",
      Pairwise_method = paste0(
        "Fixed-effect coefficient tests after re-levelling Group; ",
        "Satterthwaite degrees of freedom from lmerTest"
      ),
      N_samples = nrow(model_hcol$data),
      N_patients = nlevels(model_hcol$data$Patient),
      Group_omnibus_LRT_p = model_hcol$group_p,
      Reads_LRT_p = model_hcol$reads_p,
      Singular = singular_model,
      Convergence_message = convergence_message,
      Patient_random_intercept_variance = unname(
        variance_components["Patient_random_intercept_variance"]
      ),
      Residual_variance = unname(
        variance_components["Residual_variance"]
      ),
      stringsAsFactors = FALSE
    )
  }
}

overall_results <- bind_rows(overall_results_list) %>%
  mutate(
    Dataset = factor(Dataset, levels = dataset_levels),
    Metric = factor(Metric, levels = metrics)
  ) %>%
  arrange(Dataset, Metric) %>%
  mutate(
    Dataset = as.character(Dataset),
    Metric = as.character(Metric)
  )

pairwise_results <- bind_rows(pairwise_results_list) %>%
  mutate(
    Dataset = factor(Dataset, levels = dataset_levels),
    Metric = factor(Metric, levels = metrics),
    Comparison = factor(Comparison, levels = comparison_order)
  ) %>%
  arrange(Dataset, Metric, Comparison) %>%
  mutate(
    Dataset = as.character(Dataset),
    Metric = as.character(Metric),
    Comparison = as.character(Comparison)
  )

pairwise_results <- pairwise_results %>%
  group_by(Dataset) %>%
  mutate(
    P_BH_across_all_9_tests_in_dataset = p.adjust(
      P_raw,
      method = "BH"
    )
  ) %>%
  ungroup()

overall_results <- overall_results %>%
  group_by(Dataset) %>%
  mutate(
    Group_omnibus_p_BH_across_3_metrics = p.adjust(
      Group_omnibus_LRT_p,
      method = "BH"
    )
  ) %>%
  ungroup()

write.csv(
  overall_results,
  file.path(out_dir, "Alpha_FINAL_overall_results_BOTH_NO_EMMEANS.csv"),
  row.names = FALSE
)

write.csv(
  pairwise_results,
  file.path(out_dir, "Alpha_FINAL_pairwise_results_BOTH_NO_EMMEANS.csv"),
  row.names = FALSE
)

write.csv(
  pairwise_results %>%
    select(
      Dataset,
      Metric,
      Comparison_label,
      Estimate_first_minus_second,
      CI_low,
      CI_high,
      P_raw,
      P_Holm,
      Significance_Holm,
      Significant_Holm_0_05
    ),
  file.path(out_dir, "Alpha_FINAL_pairwise_report_table.csv"),
  row.names = FALSE
)

saveRDS(
  model_objects,
  file.path(out_dir, "Alpha_FINAL_model_objects_BOTH_NO_EMMEANS.rds")
)

# Figure functions

make_pairwise_annotations <- function(
    dataset_label,
    metric,
    p_type = c("holm", "raw")
) {
  p_type <- match.arg(p_type)

  result <- pairwise_results %>%
    filter(
      Dataset == dataset_label,
      Metric == metric
    )

  if (nrow(result) != 3L) {
    stop(
      "Expected three pairwise results for ",
      dataset_label,
      " — ",
      metric,
      ", but found ",
      nrow(result),
      "."
    )
  }

  p_column <- if (p_type == "holm") {
    "P_Holm"
  } else {
    "P_raw"
  }

  p_prefix <- if (p_type == "holm") {
    "Holm P"
  } else {
    "Raw P"
  }

  plot_data <- analysis_df %>%
    filter(Dataset == dataset_label)

  y_values <- plot_data[[metric]]
  y_min <- min(y_values, na.rm = TRUE)
  y_max <- max(y_values, na.rm = TRUE)
  y_range <- y_max - y_min

  if (!is.finite(y_range) || y_range == 0) {
    y_range <- 1
  }

  group_position <- setNames(
    seq_along(groups_keep),
    groups_keep
  )

  annotation <- data.frame(
    Comparison = comparison_order,
    Group1 = c(
      "SSL_HP",
      "SSL_HP",
      "H_ADJ_SSL_HP"
    ),
    Group2 = c(
      "H_ADJ_SSL_HP",
      "H_COL",
      "H_COL"
    ),
    y = y_max + c(
      0.14,
      0.34,
      0.54
    ) * y_range,
    stringsAsFactors = FALSE
  )

  annotation$x1 <- as.numeric(
    group_position[annotation$Group1]
  )

  annotation$x2 <- as.numeric(
    group_position[annotation$Group2]
  )

  annotation$P_value <- result[[p_column]][
    match(annotation$Comparison, result$Comparison)
  ]

  annotation$label <- vapply(
    annotation$P_value,
    format_p,
    character(1),
    prefix = p_prefix
  )

  annotation$y_tick <- annotation$y - 0.035 * y_range
  annotation$y_text <- annotation$y + 0.045 * y_range

  list(
    annotation = annotation,
    y_min = y_min,
    y_max = y_max,
    y_range = y_range
  )
}

make_alpha_plot <- function(
    dataset_label,
    metric,
    p_type = c("holm", "raw")
) {
  p_type <- match.arg(p_type)

  plot_data <- analysis_df %>%
    filter(Dataset == dataset_label) %>%
    droplevels()

  annotation_object <- make_pairwise_annotations(
    dataset_label,
    metric,
    p_type
  )

  annotation <- annotation_object$annotation

  overall_p <- overall_results$Group_omnibus_LRT_p[
    overall_results$Dataset == dataset_label &
      overall_results$Metric == metric
  ][1]

  group_counts <- table(
    factor(plot_data$Group, levels = groups_keep)
  )

  x_axis_labels <- paste0(
    group_labels[groups_keep],
    "\n(n = ",
    as.integer(group_counts[groups_keep]),
    ")"
  )

  subtitle_text <- paste0(
    dataset_label,
    " | Overall group LRT ",
    format_p(overall_p)
  )

  ggplot(
    plot_data,
    aes(
      x = Group,
      y = .data[[metric]]
    )
  ) +
    geom_boxplot(
      aes(fill = Group),
      width = 0.66,
      outlier.shape = NA,
      linewidth = 0.40,
      colour = "black",
      alpha = 0.90
    ) +
    geom_point(
      aes(colour = Group),
      position = position_jitter(
        width = 0.08,
        height = 0,
        seed = 1
      ),
      size = 1.45,
      alpha = 0.70
    ) +
    scale_fill_manual(
      values = paper_colours,
      limits = groups_keep,
      guide = "none"
    ) +
    scale_colour_manual(
      values = paper_colours,
      limits = groups_keep,
      guide = "none"
    ) +
    scale_x_discrete(
      labels = x_axis_labels,
      drop = FALSE
    ) +
    coord_cartesian(
      ylim = c(
        annotation_object$y_min,
        annotation_object$y_max +
          0.78 * annotation_object$y_range
      ),
      clip = "off"
    ) +
    geom_segment(
      data = annotation,
      aes(
        x = x1,
        xend = x2,
        y = y,
        yend = y
      ),
      inherit.aes = FALSE,
      linewidth = 0.35
    ) +
    geom_segment(
      data = annotation,
      aes(
        x = x1,
        xend = x1,
        y = y,
        yend = y_tick
      ),
      inherit.aes = FALSE,
      linewidth = 0.35
    ) +
    geom_segment(
      data = annotation,
      aes(
        x = x2,
        xend = x2,
        y = y,
        yend = y_tick
      ),
      inherit.aes = FALSE,
      linewidth = 0.35
    ) +
    geom_text(
      data = annotation,
      aes(
        x = (x1 + x2) / 2,
        y = y_text,
        label = label
      ),
      inherit.aes = FALSE,
      family = font_family,
      size = 7 / 2.845276
    ) +
    labs(
      x = NULL,
      y = metric_y_labels[[metric]],
      title = metric_titles[[metric]],
      subtitle = subtitle_text
    ) +
    theme_classic(
      base_family = font_family,
      base_size = 8
    ) +
    theme(
      text = element_text(
        family = font_family,
        colour = "black"
      ),
      axis.text = element_text(
        family = font_family,
        colour = "black",
        size = 7
      ),
      axis.text.x = element_text(
        lineheight = 0.95
      ),
      axis.title = element_text(
        family = font_family,
        colour = "black",
        size = 8
      ),
      plot.title = element_text(
        family = font_family,
        face = "bold",
        hjust = 0.5,
        size = 8
      ),
      plot.subtitle = element_text(
        family = font_family,
        hjust = 0.5,
        size = 7
      ),
      axis.line = element_line(
        linewidth = 0.35,
        colour = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.35,
        colour = "black"
      ),
      plot.margin = margin(
        5,
        5,
        5,
        5,
        unit = "mm"
      )
    )
}

# Export all individual and combined figures

primary_plots <- list()

for (dataset_label in dataset_levels) {
  dataset_stub <- dataset_stubs[[dataset_label]]

  for (metric in metrics) {
    primary_plot <- make_alpha_plot(
      dataset_label = dataset_label,
      metric = metric,
      p_type = primary_plot_p_type
    )

    raw_plot <- make_alpha_plot(
      dataset_label = dataset_label,
      metric = metric,
      p_type = "raw"
    )

    plot_key <- paste(dataset_label, metric, sep = "__")
    primary_plots[[plot_key]] <- primary_plot

    primary_stem <- paste0(
      "Alpha_",
      metric,
      "_",
      dataset_stub,
      "_HOLM"
    )

    raw_stem <- paste0(
      "Alpha_",
      metric,
      "_",
      dataset_stub,
      "_NOMINAL_RAW_P_VALUES"
    )

    ggsave(
      filename = file.path(
        out_dir,
        paste0(primary_stem, ".pdf")
      ),
      plot = primary_plot,
      width = 3.3,
      height = 5.5,
      units = "in",
      device = grDevices::pdf,
      bg = "white"
    )

    ggsave(
      filename = file.path(
        out_dir,
        paste0(primary_stem, ".png")
      ),
      plot = primary_plot,
      width = 3.3,
      height = 5.5,
      units = "in",
      dpi = export_dpi,
      bg = "white"
    )

    save_tiff_safely(
      plot = primary_plot,
      filename = file.path(
        out_dir,
        paste0(primary_stem, ".tif")
      ),
      width = 3.3,
      height = 5.5,
      dpi = export_dpi
    )

    ggsave(
      filename = file.path(
        out_dir,
        paste0(raw_stem, ".pdf")
      ),
      plot = raw_plot,
      width = 3.3,
      height = 5.5,
      units = "in",
      device = grDevices::pdf,
      bg = "white"
    )

    ggsave(
      filename = file.path(
        out_dir,
        paste0(raw_stem, ".png")
      ),
      plot = raw_plot,
      width = 3.3,
      height = 5.5,
      units = "in",
      dpi = export_dpi,
      bg = "white"
    )
  }
}

for (dataset_label in dataset_levels) {
  dataset_stub <- dataset_stubs[[dataset_label]]

  plot_list <- lapply(
    metrics,
    function(metric) {
      primary_plots[[
        paste(dataset_label, metric, sep = "__")
      ]]
    }
  )

  save_plot_grid_pdf(
    plot_list = plot_list,
    filename = file.path(
      out_dir,
      paste0(
        "Alpha_ALL_THREE_METRICS_",
        dataset_stub,
        "_HOLM.pdf"
      )
    ),
    nrow = 1,
    ncol = 3,
    width = 9.9,
    height = 5.5
  )

  save_plot_grid_png(
    plot_list = plot_list,
    filename = file.path(
      out_dir,
      paste0(
        "Alpha_ALL_THREE_METRICS_",
        dataset_stub,
        "_HOLM.png"
      )
    ),
    nrow = 1,
    ncol = 3,
    width = 9.9,
    height = 5.5,
    dpi = export_dpi
  )
}

all_six_plots <- list()

for (dataset_label in dataset_levels) {
  for (metric in metrics) {
    all_six_plots[[length(all_six_plots) + 1L]] <-
      primary_plots[[
        paste(dataset_label, metric, sep = "__")
      ]]
  }
}

save_plot_grid_pdf(
  plot_list = all_six_plots,
  filename = file.path(
    out_dir,
    "Alpha_ALL_SIX_PLOTS_BOTH_FILTERS_HOLM.pdf"
  ),
  nrow = 2,
  ncol = 3,
  width = 9.9,
  height = 11.0
)

save_plot_grid_png(
  plot_list = all_six_plots,
  filename = file.path(
    out_dir,
    "Alpha_ALL_SIX_PLOTS_BOTH_FILTERS_HOLM.png"
  ),
  nrow = 2,
  ncol = 3,
  width = 9.9,
  height = 11.0,
  dpi = export_dpi
)

# Model diagnostics

for (dataset_label in dataset_levels) {
  dataset_stub <- dataset_stubs[[dataset_label]]

  for (metric in metrics) {
    fit <- model_objects[[
      paste0(
        dataset_stub,
        "_",
        metric,
        "_H_COL_reference"
      )
    ]]

    diagnostic_file <- file.path(
      out_dir,
      paste0(
        "Alpha_",
        metric,
        "_",
        dataset_stub,
        "_diagnostics.pdf"
      )
    )

    grDevices::pdf(
      diagnostic_file,
      width = 6.5,
      height = 3.2,
      bg = "white"
    )

    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(1, 2))

    plot(
      fitted(fit),
      resid(fit),
      xlab = "Fitted values",
      ylab = "Residuals",
      main = paste(metric, "residuals")
    )

    abline(h = 0, lty = 2)

    qqnorm(
      resid(fit),
      main = paste(metric, "normal Q-Q plot")
    )

    qqline(resid(fit), lty = 2)

    par(old_par)
    grDevices::dev.off()
  }
}

# Reproducibility notes

readme_lines <- c(
  "FINAL ALPHA-DIVERSITY ANALYSIS — BOTH PREVALENCE VERSIONS",
  "=========================================================",
  "",
  paste0("Primary input: ", physeq_unfiltered_path),
  paste0("Sensitivity input: ", physeq_5pct_path),
  paste0("Output directory: ", out_dir),
  "",
  paste0(
    "5% prevalence threshold in the full without-7699 object: ",
    prevalence_threshold,
    " samples"
  ),
  paste0(
    "Unfiltered full object: ",
    nsamples(physeq_unfiltered_full),
    " samples; ",
    ntaxa(physeq_unfiltered_full),
    " ASVs"
  ),
  paste0(
    "5% filtered full object: ",
    nsamples(physeq_5pct_full),
    " samples; ",
    ntaxa(physeq_5pct_full),
    " ASVs"
  ),
  paste0(
    "Final three-group analysis: ",
    nsamples(physeq_unfiltered),
    " samples; ",
    nlevels(metadata_model$Patient),
    " patients"
  ),
  "",
  "Final models:",
  "  Observed ~ Group + Reads + (1|Patient)",
  "  Shannon ~ Group + Reads + (1|Patient)",
  "  Simpson ~ Group + Reads + (1|Patient)",
  "",
  "Read-depth covariate:",
  "  Post-rarefaction read counts from the non-prevalence-filtered object.",
  "  The same Reads value is used in both preprocessing analyses.",
  "  Reads is not log-transformed, centred or standardised.",
  "",
  "Inference:",
  "  Full and reduced models are fitted by maximum likelihood.",
  "  Overall Group and Reads effects are nested likelihood-ratio tests.",
  paste0(
    "  Pairwise comparisons are fixed-effect coefficient tests after ",
    "re-levelling Group."
  ),
  "  lmerTest supplies Satterthwaite degrees of freedom and P values.",
  "  Holm correction is applied to the three comparisons within each metric",
  "  and preprocessing version.",
  "  Main figures show Holm-adjusted pairwise P values.",
  "  Separate NOMINAL_RAW_P_VALUES figures are also exported.",
  "",
  "Interpretation:",
  "  The no-prevalence-filter Shannon analysis is the primary analysis.",
  "  The 5% prevalence-filtered Shannon analysis is the sensitivity analysis.",
  paste0(
    "  Observed in the 5% object means the number of ASVs retained after ",
    "the 5% filter, not total unfiltered richness."
  ),
  "  Both prevalence versions contain the same samples.",
  "",
  "Reporting:",
  "  Overall group likelihood-ratio-test P values are reported without adjustment.",
  "  Pairwise reporting uses the P_Holm column.",
  "  Additional BH-adjusted columns are retained in the full results tables."
)

writeLines(
  readme_lines,
  file.path(out_dir, "README_Alpha_FINAL_BOTH_NO_EMMEANS.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "sessionInfo_Alpha_FINAL_BOTH_NO_EMMEANS.txt"
  )
)

# Print key results

cat("\n===== FINAL OVERALL RESULTS =====\n")
print(overall_results)

cat("\n===== FINAL PAIRWISE RESULTS =====\n")
print(pairwise_results)

cat("\nPrimary analysis:\n")
cat("  Shannon — No prevalence filter\n")
cat("Sensitivity analysis:\n")
cat("  Shannon — 5% prevalence filter\n")

cat("\nAll files were saved to:\n")
cat(out_dir, "\n")
