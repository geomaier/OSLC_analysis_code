# MUC5AC qPCR: technical-replicate QC and longitudinal mixed models.
# Input requirements and output locations are described in ../README.md.

rm(list = ls())

required_path <- function(variable) {
  value <- Sys.getenv(variable)
  if (!nzchar(value)) stop("Set ", variable, " to the input directory before running this script.")
  normalizePath(value, winslash="/", mustWork=TRUE)
}


# Create a short, unused child directory; preserve earlier analysis outputs.
review_output_dir <- function(parent) {
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(9999L)) {
    candidate <- file.path(parent, paste0("r", i))
    if (!file.exists(candidate) && dir.create(candidate, showWarnings = FALSE)) return(candidate)
  }
  stop("No unused review output directory beneath: ", parent)
}

# Configuration and technical-replicate QC thresholds ----
set.seed(26082026)
options(stringsAsFactors = FALSE)
options(contrasts = c("contr.treatment", "contr.poly"))

SCRIPT_VERSION <- "M26F_MUC5AC_FIXED_V2_2026-08-17"
message("\nRunning corrected script: ", SCRIPT_VERSION)
message("This version removes Excel-export artefact rows before constructing the analysis table.\n")

base_dir <- required_path("OSLC_ORGANOID_PCR_DIR")

input_file <- file.path(base_dir, "Organoids_MUC5AC_ANALYSIS.csv")
output_dir <- review_output_dir(file.path(base_dir, "MUC5AC_Organoids_Final_Analysis_Aug26"))

overwrite_existing_outputs <- FALSE
install_missing_packages <- FALSE

outlier_min_distance_from_median <- 1.50
outlier_min_deviation_ratio <- 2.50
outlier_max_remaining_pair_range <- 1.25

high_dispersion_range_threshold <- 1.25
minimum_technical_replicates <- 2L

expected_clear_outlier_samples <- c(
  "GM_7699_SSL_D_P2",
  "GM_7699_H_ADJ_D_P1",
  "GM_11135_H_ADJ_A_P6",
  "GM_11725_SSL_A_P4",
  "GM_12147_SSL_T_P3",
  "GM_12147_SSL_T_P4",
  "GM_12147_SSL_T_P8",
  "GM_12147_H_ADJ_T_P3",
  "GM_12147_H_ADJ_T_P4",
  "GM_12147_H_ADJ_T_P7"
)

font_family <- "Arial"
base_font_pt <- 7
axis_font_pt <- 7
legend_font_pt <- 6
plot_title_pt <- 8
panel_tag_pt <- 8
export_dpi <- 600

paper_cols <- c(
  SSL_HP = "#FF0000",
  H_ADJ_SSL_HP = "#E9C61D"
)

facet_labels <- c(
  SSL_HP = "SSL/HP-derived organoids",
  H_ADJ_SSL_HP = "H_ADJ_SSL_HP organoids"
)

# Required packages ----
cran_packages <- c(
  "data.table", "dplyr", "tidyr", "stringr", "ggplot2",
  "tibble", "scales", "nlme", "patchwork"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!install_missing_packages) {
      stop(
        "Package '", pkg, "' is not installed. Install it or set ",
        "install_missing_packages <- TRUE."
      )
    }
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

for (pkg in cran_packages) install_if_missing(pkg)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(tibble)
  library(scales)
  library(nlme)
  library(patchwork)
  library(grid)
})

if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

if (!dir.exists(base_dir)) {
  stop("The input directory does not exist:\n", base_dir)
}

if (!file.exists(input_file)) {
  stop(
    "The input file was not found:\n", input_file, "\n\n",
    "Confirm that it is named exactly Organoids_MUC5AC_ANALYSIS.csv."
  )
}

if (overwrite_existing_outputs && dir.exists(output_dir)) {
  base_norm <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
  out_norm <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

  if (identical(base_norm, out_norm)) {
    stop("Refusing to delete output_dir because it equals base_dir.")
  }

  unlink(output_dir, recursive = TRUE, force = TRUE)
  if (dir.exists(output_dir)) {
    stop(
      "The previous output directory could not be removed. Close any open ",
      "output files and rerun:\n", output_dir
    )
  }
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

subdirs <- list(
  qc = file.path(output_dir, "01_QC"),
  tables = file.path(output_dir, "02_Tables"),
  figures = file.path(output_dir, "03_Figures"),
  statistics = file.path(output_dir, "04_Statistics")
)

for (d in subdirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Input, plotting and model-summary helpers ----
safe_filename <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

to_numeric <- function(x) {
  x <- clean_character(x)
  x[toupper(x) %in% c("N/A", "NA", "NAN", "NULL", ".")] <- NA_character_
  suppressWarnings(as.numeric(gsub(",", "", x, fixed = TRUE)))
}

column_key <- function(x) {
  toupper(gsub("[^A-Za-z0-9]", "", trimws(as.character(x))))
}

detect_column <- function(df, candidates, required = TRUE, label = "column") {
  nms <- names(df)
  idx <- match(column_key(candidates), column_key(nms))
  idx <- idx[!is.na(idx)]

  if (length(idx) > 0L) return(nms[idx[1]])

  if (required) {
    stop(
      "Could not identify the ", label, ". Tried: ",
      paste(candidates, collapse = ", "), "\nAvailable columns:\n",
      paste(nms, collapse = "\n")
    )
  }

  NA_character_
}

normalise_organoid_type <- function(x, sample_name) {
  raw <- toupper(gsub("[^A-Za-z0-9]", "", clean_character(x)))
  name_key <- toupper(as.character(sample_name))

  out <- rep(NA_character_, length(raw))
  out[raw %in% c("SSL", "HP", "SSLHP")] <- "SSL_HP"
  out[raw %in% c("HADJ", "HADJSSLHP", "HEALTHYADJACENT", "ADJACENT")] <-
    "H_ADJ_SSL_HP"

  missing <- is.na(out)
  out[missing & grepl("_H_ADJ_", name_key)] <- "H_ADJ_SSL_HP"
  out[missing & grepl("_(SSL|HP)_", name_key)] <- "SSL_HP"
  out
}

normalise_braf <- function(x) {
  key <- toupper(gsub("[^A-Za-z0-9]", "", clean_character(x)))
  out <- rep(NA_character_, length(key))
  out[grepl("V600E|MUTANT|MUT", key)] <- "V600E"
  out[key %in% c("WT", "WILDTYPE", "BRAFWT")] <- "WT"
  out
}

derive_trajectory_id <- function(sample_name) {
  x <- toupper(trimws(as.character(sample_name)))
  x <- stringr::str_replace(x, "_P[0-9]+$", "")
  x
}

derive_patient_id <- function(sample_name) {
  x <- toupper(trimws(as.character(sample_name)))
  patient <- stringr::str_extract(x, "^GM_[A-Z0-9]+")
  patient[is.na(patient)] <- x[is.na(patient)]
  patient
}

nature_theme <- function(base_size = base_font_pt) {
  theme_classic(base_family = font_family, base_size = base_size) +
    theme(
      text = element_text(family = font_family, colour = "black"),
      plot.title = element_text(
        family = font_family,
        face = "bold",
        size = plot_title_pt,
        hjust = 0
      ),
      plot.subtitle = element_text(
        family = font_family,
        size = 6,
        hjust = 0
      ),
      axis.title = element_text(
        family = font_family,
        size = axis_font_pt,
        colour = "black"
      ),
      axis.text = element_text(
        family = font_family,
        size = axis_font_pt,
        colour = "black"
      ),
      legend.title = element_text(
        family = font_family,
        size = legend_font_pt,
        colour = "black"
      ),
      legend.text = element_text(
        family = font_family,
        size = legend_font_pt,
        colour = "black"
      ),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = grid::unit(1.2, "mm"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(5, 7, 5, 5, "mm")
    )
}

theme_set(nature_theme())

save_plot_set <- function(plot, directory, stem, width_mm, height_mm) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  stem <- safe_filename(stem)
  prefix <- file.path(directory, stem)
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

  ggsave(
    paste0(prefix, ".pdf"), plot,
    width = width_mm, height = height_mm, units = "mm",
    device = pdf_device, bg = "white"
  )

  ggsave(
    paste0(prefix, ".png"), plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = export_dpi, bg = "white"
  )

  ggsave(
    paste0(prefix, ".tif"), plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = export_dpi, device = grDevices::tiff,
    compression = "lzw", bg = "white"
  )

  ggsave(
    paste0(prefix, ".jpg"), plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = export_dpi, device = grDevices::jpeg,
    quality = 100, bg = "white"
  )
}

extract_lme_fixed <- function(fit, model_name) {
  tt <- as.data.frame(summary(fit)$tTable)
  tt$Term <- rownames(tt)
  rownames(tt) <- NULL

  data.frame(
    Model = model_name,
    Term = tt$Term,
    Estimate = tt$Value,
    SE = tt$Std.Error,
    DF = tt$DF,
    t_value = tt$`t-value`,
    P_value = tt$`p-value`,
    stringsAsFactors = FALSE
  )
}

write_model_failure <- function(model_name, error_object) {
  data.frame(
    Model = model_name,
    Term = NA_character_,
    Estimate = NA_real_,
    SE = NA_real_,
    DF = NA_real_,
    t_value = NA_real_,
    P_value = NA_real_,
    Model_status = paste0("FAILED: ", conditionMessage(error_object)),
    stringsAsFactors = FALSE
  )
}

# CSV validation and metadata preparation ----
header_only <- data.table::fread(
  input_file,
  nrows = 0L,
  data.table = FALSE,
  check.names = FALSE
)

if (ncol(header_only) < 10L) {
  stop(
    "The input CSV contains fewer than 10 columns. Found ", ncol(header_only),
    ". Inspect the file before rerunning."
  )
}

raw <- data.table::fread(
  input_file,
  select = seq_len(10L),
  colClasses = "character",
  data.table = FALSE,
  check.names = FALSE,
  fill = TRUE,
  na.strings = c("", "NA", "N/A", "n/a", "NaN", "NULL")
)

names(raw) <- trimws(names(raw))
raw$.Input_row <- seq_len(nrow(raw))

sample_col <- detect_column(
  raw,
  c("Sample Name", "Sample_Name", "Sample", "SampleID", "Sample_ID", "Name"),
  label = "sample-name column"
)

sample_values_all <- clean_character(raw[[sample_col]])
valid_sample_row <- !is.na(sample_values_all) & grepl("^GM_", sample_values_all)

nonempty_anywhere <- rowSums(
  !is.na(raw[setdiff(names(raw), ".Input_row")])
) > 0L

malformed_nonempty <- raw[!valid_sample_row & nonempty_anywhere, , drop = FALSE]

import_audit <- data.frame(
  Input_file = input_file,
  Rows_read_from_CSV = nrow(raw),
  Valid_rows_with_Sample_Name_beginning_GM = sum(valid_sample_row),
  Discarded_rows_total = sum(!valid_sample_row),
  Discarded_nonempty_malformed_rows = nrow(malformed_nonempty),
  Columns_in_original_CSV = ncol(header_only),
  Columns_read_for_analysis = 10L,
  stringsAsFactors = FALSE
)

data.table::fwrite(
  import_audit,
  file.path(subdirs$qc, "MUC5AC_CSV_import_audit.csv")
)

if (nrow(malformed_nonempty) > 0L) {
  data.table::fwrite(
    utils::head(malformed_nonempty, 200L),
    file.path(subdirs$qc, "MUC5AC_discarded_malformed_row_examples_first200.csv")
  )
}

raw <- raw[valid_sample_row, , drop = FALSE]

if (nrow(raw) == 0L) {
  stop("No valid sample rows beginning with 'GM_' were retained from the CSV.")
}

message(
  "CSV import: ", import_audit$Rows_read_from_CSV,
  " rows read; ", import_audit$Valid_rows_with_Sample_Name_beginning_GM,
  " genuine sample rows retained; ", import_audit$Discarded_rows_total,
  " blank/malformed Excel-export rows ignored."
)

passage_col <- detect_column(raw, c("Passage", "Passage_Number", "P"), label = "passage column")
type_col <- detect_column(raw, c("Organoid_Type", "Organoid Type", "Type"), label = "organoid-type column")
topology_col <- detect_column(raw, c("Topology", "Location", "Site"), required = FALSE)
gender_col <- detect_column(raw, c("Gender", "Sex"), required = FALSE)
braf_col <- detect_column(raw, c("BRAF", "BRAF_Status", "BRAF Status"), required = FALSE)
infection_col <- detect_column(raw, c("Infection_Status", "Infection Status", "Infection"), required = FALSE)
dct1_col <- detect_column(raw, c("deltaCT1", "deltaCt1", "DeltaCT1", "dCt1"), label = "deltaCT1 column")
dct2_col <- detect_column(raw, c("deltaCT2", "deltaCt2", "DeltaCT2", "dCt2"), label = "deltaCT2 column")
dct3_col <- detect_column(raw, c("deltaCT3", "deltaCt3", "DeltaCT3", "dCt3"), label = "deltaCT3 column")

sample_name <- clean_character(raw[[sample_col]])
passage <- to_numeric(raw[[passage_col]])
organoid_group <- normalise_organoid_type(raw[[type_col]], sample_name)

input <- data.frame(
  Input_row = raw$.Input_row,
  Sample_Name = sample_name,
  Passage = passage,
  Organoid_Type_raw = clean_character(raw[[type_col]]),
  Organoid_Group = organoid_group,
  Patient_ID = derive_patient_id(sample_name),
  Trajectory_ID = derive_trajectory_id(sample_name),
  Topology = if (!is.na(topology_col)) clean_character(raw[[topology_col]]) else NA_character_,
  Gender = if (!is.na(gender_col)) clean_character(raw[[gender_col]]) else NA_character_,
  BRAF_raw = if (!is.na(braf_col)) clean_character(raw[[braf_col]]) else NA_character_,
  BRAF_Status = if (!is.na(braf_col)) normalise_braf(raw[[braf_col]]) else NA_character_,
  Infection_Status = if (!is.na(infection_col)) clean_character(raw[[infection_col]]) else NA_character_,
  deltaCT1 = to_numeric(raw[[dct1_col]]),
  deltaCT2 = to_numeric(raw[[dct2_col]]),
  deltaCT3 = to_numeric(raw[[dct3_col]]),
  stringsAsFactors = FALSE
)

if (any(is.na(input$Sample_Name) | input$Sample_Name == "")) {
  stop("A genuine retained row still has a missing Sample Name.")
}

if (anyDuplicated(input$Sample_Name)) {
  duplicate_names <- unique(input$Sample_Name[duplicated(input$Sample_Name)])
  stop(
    "Duplicated Sample Name values were found:\n",
    paste(duplicate_names, collapse = "\n")
  )
}

if (any(!is.finite(input$Passage))) {
  bad <- input %>% dplyr::filter(!is.finite(Passage))
  data.table::fwrite(
    bad,
    file.path(subdirs$qc, "MUC5AC_rows_with_unresolved_passage.csv")
  )
  stop("At least one retained sample has an unresolved passage.")
}

if (any(is.na(input$Organoid_Group))) {
  bad <- input %>% dplyr::filter(is.na(Organoid_Group))
  data.table::fwrite(
    bad,
    file.path(subdirs$qc, "MUC5AC_rows_with_unresolved_organoid_group.csv")
  )
  stop("At least one retained sample has an unresolved organoid group.")
}

input$Organoid_Group <- factor(
  input$Organoid_Group,
  levels = c("H_ADJ_SSL_HP", "SSL_HP")
)

input$BRAF_Status <- factor(
  input$BRAF_Status,
  levels = c("WT", "V600E")
)

input <- input %>%
  arrange(Organoid_Group, Patient_ID, Trajectory_ID, Passage, Sample_Name)

data.table::fwrite(
  input,
  file.path(subdirs$qc, "MUC5AC_input_cleaned_valid_sample_rows.csv")
)

data.table::fwrite(
  input %>%
    select(
      Sample_Name, Passage, Organoid_Type_raw, Organoid_Group,
      Topology, Gender, BRAF_raw, BRAF_Status, Infection_Status,
      deltaCT1, deltaCT2, deltaCT3
    ),
  file.path(output_dir, "Organoids_MUC5AC_ANALYSIS_CLEANED_VALID_ROWS.csv")
)

# Technical-replicate QC and sample means ----
qc_triplicate <- function(values) {
  stopifnot(length(values) == 3L)

  replicate_names <- paste0("deltaCT", seq_len(3L))
  valid_idx <- which(is.finite(values))
  valid_values <- values[valid_idx]
  n_valid <- length(valid_values)

  used_mask <- rep(FALSE, 3L)
  excluded_replicate <- NA_integer_
  outlier_removed <- FALSE
  high_dispersion <- FALSE
  status <- NA_character_
  comment <- NA_character_

  raw_range <- if (n_valid >= 2L) diff(range(valid_values)) else NA_real_
  median_raw <- if (n_valid >= 1L) stats::median(valid_values) else NA_real_

  if (n_valid < minimum_technical_replicates) {
    status <- "EXCLUDED_INSUFFICIENT_TECHNICAL_REPLICATES"
    comment <- paste0("Only ", n_valid, " finite technical replicate(s) available.")
  } else if (n_valid == 2L) {
    used_mask[valid_idx] <- TRUE
    high_dispersion <- diff(range(valid_values)) > high_dispersion_range_threshold
    status <- if (high_dispersion) {
      "INCLUDED_TWO_REPLICATES_HIGH_DISPERSION"
    } else {
      "INCLUDED_TWO_REPLICATES"
    }
    comment <- "One technical replicate was missing; the two available values were used."
  } else {
    median_value <- stats::median(values)
    deviations <- abs(values - median_value)
    ordered_deviations <- sort(deviations, decreasing = TRUE)
    max_deviation <- ordered_deviations[1]
    second_deviation <- ordered_deviations[2]
    max_idx <- which.max(deviations)
    unique_maximum <- sum(abs(deviations - max_deviation) < 1e-12) == 1L
    deviation_ratio <- max_deviation / max(second_deviation, 1e-12)
    remaining_values <- values[-max_idx]
    remaining_range <- diff(range(remaining_values))

    clear_outlier <-
      unique_maximum &&
      max_deviation >= outlier_min_distance_from_median &&
      deviation_ratio >= outlier_min_deviation_ratio &&
      remaining_range <= outlier_max_remaining_pair_range

    if (clear_outlier) {
      used_mask[] <- TRUE
      used_mask[max_idx] <- FALSE
      excluded_replicate <- max_idx
      outlier_removed <- TRUE
      high_dispersion <- remaining_range > high_dispersion_range_threshold
      status <- "INCLUDED_AFTER_CLEAR_OUTLIER_REMOVAL"
      comment <- paste0(
        replicate_names[max_idx], " removed as a clearly isolated technical outlier."
      )
    } else {
      used_mask[] <- TRUE
      high_dispersion <- raw_range > high_dispersion_range_threshold
      status <- if (high_dispersion) {
        "INCLUDED_ALL_THREE_HIGH_DISPERSION_REVIEW"
      } else {
        "INCLUDED_ALL_THREE"
      }
      comment <- if (high_dispersion) {
        paste0(
          "No single replicate met the clear-outlier rule; all three were retained, ",
          "but the triplicate range exceeded ", high_dispersion_range_threshold,
          " deltaCt."
        )
      } else {
        "All three technical replicates retained."
      }
    }
  }

  used_values <- values[used_mask & is.finite(values)]
  n_used <- length(used_values)
  used_range <- if (n_used >= 2L) diff(range(used_values)) else NA_real_
  mean_deltaCt <- if (n_used >= minimum_technical_replicates) mean(used_values) else NA_real_
  sd_deltaCt <- if (n_used >= 2L) stats::sd(used_values) else NA_real_
  sem_deltaCt <- if (n_used >= 2L) sd_deltaCt / sqrt(n_used) else NA_real_

  primary_included <- is.finite(mean_deltaCt)
  strict_included <- primary_included && !high_dispersion

  data.frame(
    N_valid_raw = n_valid,
    N_used = n_used,
    Raw_range_deltaCt = raw_range,
    Median_raw_deltaCt = median_raw,
    Outlier_removed = outlier_removed,
    Excluded_replicate_number = excluded_replicate,
    Excluded_replicate_name = ifelse(
      is.na(excluded_replicate),
      NA_character_,
      replicate_names[excluded_replicate]
    ),
    Excluded_replicate_value = ifelse(
      is.na(excluded_replicate),
      NA_real_,
      values[excluded_replicate]
    ),
    Included_replicates = paste(replicate_names[used_mask], collapse = ";"),
    Used_range_deltaCt = used_range,
    Mean_deltaCt = mean_deltaCt,
    SD_deltaCt = sd_deltaCt,
    SEM_deltaCt = sem_deltaCt,
    High_dispersion_flag = high_dispersion,
    Primary_analysis_included = primary_included,
    Strict_QC_sensitivity_included = strict_included,
    QC_status = status,
    QC_comment = comment,
    stringsAsFactors = FALSE
  )
}

qc_results <- lapply(seq_len(nrow(input)), function(i) {
  qc_triplicate(as.numeric(input[i, c("deltaCT1", "deltaCT2", "deltaCT3")]))
})

sample_qc <- dplyr::bind_cols(input, dplyr::bind_rows(qc_results))

# Replicate-level QC export; missing exclusion index means no well was removed ----
replicate_qc <- sample_qc %>%
  select(
    Sample_Name, Patient_ID, Trajectory_ID, Passage, Organoid_Group,
    deltaCT1, deltaCT2, deltaCT3,
    Excluded_replicate_number, Primary_analysis_included,
    Strict_QC_sensitivity_included, QC_status
  ) %>%
  pivot_longer(
    cols = c(deltaCT1, deltaCT2, deltaCT3),
    names_to = "Replicate",
    values_to = "deltaCt"
  ) %>%
  mutate(
    Replicate_number = as.integer(stringr::str_extract(Replicate, "[0-9]+")),
    Finite_measurement = is.finite(deltaCt),
    Automatically_excluded_as_outlier =
      Finite_measurement & !is.na(Excluded_replicate_number) &
      Replicate_number == Excluded_replicate_number,
    Used_in_primary_sample_mean =
      Finite_measurement & !Automatically_excluded_as_outlier &
      Primary_analysis_included
  )

data.table::fwrite(
  sample_qc,
  file.path(subdirs$qc, "MUC5AC_technical_replicate_QC_COMPLETE.csv")
)

data.table::fwrite(
  replicate_qc,
  file.path(subdirs$qc, "MUC5AC_technical_replicate_QC_LONG.csv")
)

qc_summary <- sample_qc %>%
  count(QC_status, name = "N_samples") %>%
  arrange(desc(N_samples), QC_status)

data.table::fwrite(
  qc_summary,
  file.path(subdirs$qc, "MUC5AC_technical_replicate_QC_SUMMARY.csv")
)

clear_outliers_observed <- sample_qc %>%
  filter(Outlier_removed) %>%
  pull(Sample_Name)

missing_expected <- setdiff(expected_clear_outlier_samples, clear_outliers_observed)
unexpected_observed <- setdiff(clear_outliers_observed, expected_clear_outlier_samples)

outlier_expectation_audit <- data.frame(
  Audit_item = c(
    "Expected clear-outlier samples",
    "Observed clear-outlier samples",
    "Expected but not observed",
    "Observed but not expected"
  ),
  Samples = c(
    paste(expected_clear_outlier_samples, collapse = ";"),
    paste(clear_outliers_observed, collapse = ";"),
    paste(missing_expected, collapse = ";"),
    paste(unexpected_observed, collapse = ";")
  ),
  stringsAsFactors = FALSE
)

data.table::fwrite(
  outlier_expectation_audit,
  file.path(subdirs$qc, "MUC5AC_expected_vs_observed_clear_outliers.csv")
)

if (length(missing_expected) > 0L || length(unexpected_observed) > 0L) {
  warning(
    "The observed automatic outlier list differs from the supplied-sheet audit. ",
    "Inspect MUC5AC_expected_vs_observed_clear_outliers.csv."
  )
}

# Primary and strict-QC datasets, with P1-relative values ----
sample_qc <- sample_qc %>%
  mutate(
    Minus_deltaCt = -Mean_deltaCt,
    Relative_expression_2_minus_deltaCt = 2^(-Mean_deltaCt),
    Passage_centered = Passage - 1,
    Passage_factor = factor(Passage, levels = sort(unique(Passage)))
  )

analysis_primary <- sample_qc %>%
  filter(Primary_analysis_included) %>%
  arrange(Organoid_Group, Trajectory_ID, Passage)

analysis_strict <- sample_qc %>%
  filter(Strict_QC_sensitivity_included) %>%
  arrange(Organoid_Group, Trajectory_ID, Passage)

p1_reference <- analysis_primary %>%
  filter(Passage == 1) %>%
  select(
    Trajectory_ID,
    P1_Mean_deltaCt = Mean_deltaCt,
    P1_Minus_deltaCt = Minus_deltaCt
  )

analysis_p1_relative <- analysis_primary %>%
  inner_join(p1_reference, by = "Trajectory_ID") %>%
  mutate(
    DeltaDeltaCt_vs_P1 = Mean_deltaCt - P1_Mean_deltaCt,
    Log2_fold_change_vs_P1 = -DeltaDeltaCt_vs_P1,
    Fold_change_vs_P1 = 2^(-DeltaDeltaCt_vs_P1)
  )

data.table::fwrite(
  sample_qc,
  file.path(subdirs$tables, "MUC5AC_all_samples_with_QC_and_derived_values.csv")
)

data.table::fwrite(
  analysis_primary,
  file.path(subdirs$tables, "MUC5AC_primary_sample_level_analysis.csv")
)

data.table::fwrite(
  analysis_strict,
  file.path(subdirs$tables, "MUC5AC_strict_QC_sensitivity_sample_level_analysis.csv")
)

data.table::fwrite(
  analysis_p1_relative,
  file.path(subdirs$tables, "MUC5AC_P1_relative_delta_delta_Ct_and_fold_change.csv")
)

trajectory_summary <- analysis_primary %>%
  group_by(Patient_ID, Trajectory_ID, Organoid_Group, BRAF_Status, Topology) %>%
  summarise(
    N_passages = n(),
    First_passage = min(Passage),
    Last_passage = max(Passage),
    Passages = paste(sort(unique(Passage)), collapse = ";"),
    N_clear_outlier_wells_removed = sum(Outlier_removed),
    N_high_dispersion_samples = sum(High_dispersion_flag),
    .groups = "drop"
  )

data.table::fwrite(
  trajectory_summary,
  file.path(subdirs$tables, "MUC5AC_organoid_trajectory_summary.csv")
)

summarise_by_group_passage <- function(dat, analysis_label) {
  dat %>%
    group_by(Organoid_Group, Passage) %>%
    summarise(
      Analysis = analysis_label,
      N_samples = n(),
      N_trajectories = n_distinct(Trajectory_ID),
      Mean_deltaCt = mean(Mean_deltaCt),
      SD_deltaCt = ifelse(n() > 1L, sd(Mean_deltaCt), NA_real_),
      SEM_deltaCt = ifelse(n() > 1L, sd(Mean_deltaCt) / sqrt(n()), NA_real_),
      Mean_minus_deltaCt = mean(Minus_deltaCt),
      SD_minus_deltaCt = ifelse(n() > 1L, sd(Minus_deltaCt), NA_real_),
      SEM_minus_deltaCt = ifelse(n() > 1L, sd(Minus_deltaCt) / sqrt(n()), NA_real_),
      Geometric_mean_relative_expression = 2^(Mean_minus_deltaCt),
      Geometric_mean_lower_SEM = ifelse(
        is.finite(SEM_minus_deltaCt),
        2^(Mean_minus_deltaCt - SEM_minus_deltaCt),
        NA_real_
      ),
      Geometric_mean_upper_SEM = ifelse(
        is.finite(SEM_minus_deltaCt),
        2^(Mean_minus_deltaCt + SEM_minus_deltaCt),
        NA_real_
      ),
      Median_relative_expression = median(Relative_expression_2_minus_deltaCt),
      .groups = "drop"
    )
}

group_summary_primary <- summarise_by_group_passage(
  analysis_primary,
  "Primary: clear outlier wells removed; high-dispersion samples retained"
)

group_summary_strict <- summarise_by_group_passage(
  analysis_strict,
  "Strict QC: high-dispersion samples excluded"
)

data.table::fwrite(
  bind_rows(group_summary_primary, group_summary_strict),
  file.path(subdirs$tables, "MUC5AC_group_passage_descriptive_statistics.csv")
)

# Trajectory and expression figures ----
make_relative_expression_trajectory_plot <- function(
    dat,
    group_summary,
    title_text,
    subtitle_text = NULL
) {
  ggplot(
    dat,
    aes(
      x = Passage,
      y = Relative_expression_2_minus_deltaCt,
      group = Trajectory_ID,
      colour = Organoid_Group
    )
  ) +
    geom_line(linewidth = 0.38, alpha = 0.42) +
    geom_point(
      aes(fill = Organoid_Group),
      shape = 21,
      size = 1.75,
      stroke = 0.30,
      colour = "black",
      alpha = 0.88
    ) +
    geom_errorbar(
      data = group_summary,
      aes(
        x = Passage,
        ymin = Geometric_mean_lower_SEM,
        ymax = Geometric_mean_upper_SEM,
        colour = Organoid_Group,
        group = Organoid_Group
      ),
      inherit.aes = FALSE,
      width = 0.10,
      linewidth = 0.40,
      na.rm = TRUE
    ) +
    geom_line(
      data = group_summary,
      aes(
        x = Passage,
        y = Geometric_mean_relative_expression,
        colour = Organoid_Group,
        group = Organoid_Group
      ),
      inherit.aes = FALSE,
      linewidth = 1.05
    ) +
    geom_point(
      data = group_summary,
      aes(
        x = Passage,
        y = Geometric_mean_relative_expression,
        fill = Organoid_Group
      ),
      inherit.aes = FALSE,
      shape = 21,
      size = 2.35,
      stroke = 0.42,
      colour = "black"
    ) +
    facet_wrap(
      ~Organoid_Group,
      nrow = 1,
      labeller = as_labeller(facet_labels)
    ) +
    scale_colour_manual(values = paper_cols, guide = "none", drop = FALSE) +
    scale_fill_manual(values = paper_cols, guide = "none", drop = FALSE) +
    scale_x_continuous(
      breaks = 1:9,
      limits = c(0.8, 9.2),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_log10(
      labels = scales::label_scientific(digits = 1),
      expand = expansion(mult = c(0.08, 0.14))
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Passage",
      y = expression("GAPDH-normalised " * italic(MUC5AC) * " expression (" * 2^{-Delta * Ct} * ")")
    ) +
    nature_theme() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(
        family = font_family,
        face = "bold",
        size = 7,
        colour = "black"
      ),
      panel.spacing.x = unit(8, "mm"),
      plot.margin = margin(6, 8, 6, 6, "mm")
    )
}

p_relative_primary <- make_relative_expression_trajectory_plot(
  analysis_primary,
  group_summary_primary,
  "MUC5AC expression changes during serial organoid passage",
  "Thin lines show individual organoid trajectories; thick lines show group geometric means +/- SEM"
)

save_plot_set(
  p_relative_primary,
  subdirs$figures,
  "MUC5AC_TaqMan_longitudinal_trajectories_FINAL",
  183,
  92
)

p_relative_strict <- make_relative_expression_trajectory_plot(
  analysis_strict,
  group_summary_strict,
  "MUC5AC expression during serial passage: strict technical-QC sensitivity",
  paste0(
    "Samples with a final technical-replicate range > ",
    high_dispersion_range_threshold,
    " deltaCt were excluded"
  )
)

save_plot_set(
  p_relative_strict,
  subdirs$figures,
  "MUC5AC_TaqMan_longitudinal_trajectories_STRICT_QC_SENSITIVITY",
  183,
  92
)

group_summary_minus <- group_summary_primary

p_minus_deltaCt <- ggplot(
  analysis_primary,
  aes(
    x = Passage,
    y = Minus_deltaCt,
    group = Trajectory_ID,
    colour = Organoid_Group
  )
) +
  geom_line(linewidth = 0.38, alpha = 0.42) +
  geom_point(
    aes(fill = Organoid_Group),
    shape = 21,
    size = 1.75,
    stroke = 0.30,
    colour = "black",
    alpha = 0.88
  ) +
  geom_errorbar(
    data = group_summary_minus,
    aes(
      x = Passage,
      ymin = Mean_minus_deltaCt - SEM_minus_deltaCt,
      ymax = Mean_minus_deltaCt + SEM_minus_deltaCt,
      colour = Organoid_Group,
      group = Organoid_Group
    ),
    inherit.aes = FALSE,
    width = 0.10,
    linewidth = 0.40,
    na.rm = TRUE
  ) +
  geom_line(
    data = group_summary_minus,
    aes(
      x = Passage,
      y = Mean_minus_deltaCt,
      colour = Organoid_Group,
      group = Organoid_Group
    ),
    inherit.aes = FALSE,
    linewidth = 1.05
  ) +
  geom_point(
    data = group_summary_minus,
    aes(x = Passage, y = Mean_minus_deltaCt, fill = Organoid_Group),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.35,
    stroke = 0.42,
    colour = "black"
  ) +
  facet_wrap(~Organoid_Group, nrow = 1, labeller = as_labeller(facet_labels)) +
  scale_colour_manual(values = paper_cols, guide = "none", drop = FALSE) +
  scale_fill_manual(values = paper_cols, guide = "none", drop = FALSE) +
  scale_x_continuous(breaks = 1:9, limits = c(0.8, 9.2)) +
  labs(
    title = "MUC5AC expression across serial organoid passage",
    subtitle = "Higher -deltaCt values indicate greater GAPDH-normalised expression",
    x = "Passage",
    y = expression(-Delta * Ct)
  ) +
  nature_theme() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(
      family = font_family,
      face = "bold",
      size = 7,
      colour = "black"
    ),
    panel.spacing.x = unit(8, "mm")
  )

save_plot_set(
  p_minus_deltaCt,
  subdirs$figures,
  "MUC5AC_TaqMan_longitudinal_minus_deltaCt",
  183,
  92
)

p_p1_relative <- ggplot(
  analysis_p1_relative,
  aes(
    x = Passage,
    y = Log2_fold_change_vs_P1,
    group = Trajectory_ID,
    colour = Organoid_Group
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.30, colour = "grey45") +
  geom_line(linewidth = 0.42, alpha = 0.65) +
  geom_point(
    aes(fill = Organoid_Group),
    shape = 21,
    size = 1.9,
    stroke = 0.32,
    colour = "black"
  ) +
  facet_wrap(~Organoid_Group, nrow = 1, labeller = as_labeller(facet_labels)) +
  scale_colour_manual(values = paper_cols, guide = "none", drop = FALSE) +
  scale_fill_manual(values = paper_cols, guide = "none", drop = FALSE) +
  scale_x_continuous(breaks = 1:9, limits = c(0.8, 9.2)) +
  labs(
    title = "MUC5AC expression relative to P1 within each organoid trajectory",
    subtitle = "Trajectories without a valid P1 measurement are excluded",
    x = "Passage",
    y = expression(log[2] * " fold change relative to P1")
  ) +
  nature_theme() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(family = font_family, face = "bold", size = 7),
    panel.spacing.x = unit(8, "mm")
  )

save_plot_set(
  p_p1_relative,
  subdirs$figures,
  "MUC5AC_TaqMan_P1_relative_log2_fold_change_SUPPLEMENTARY",
  183,
  92
)

heatmap_grid <- tidyr::expand_grid(
  Trajectory_ID = unique(analysis_primary$Trajectory_ID),
  Passage = sort(unique(input$Passage))
) %>%
  left_join(
    analysis_primary %>%
      select(
        Trajectory_ID, Passage, Organoid_Group, Patient_ID,
        BRAF_Status, Minus_deltaCt, Mean_deltaCt, QC_status
      ),
    by = c("Trajectory_ID", "Passage")
  ) %>%
  group_by(Trajectory_ID) %>%
  fill(Organoid_Group, Patient_ID, BRAF_Status, .direction = "downup") %>%
  ungroup()

trajectory_order <- analysis_primary %>%
  distinct(Organoid_Group, Patient_ID, Trajectory_ID) %>%
  arrange(Organoid_Group, Patient_ID, Trajectory_ID) %>%
  pull(Trajectory_ID)

heatmap_grid <- heatmap_grid %>%
  mutate(
    Trajectory_ID = factor(Trajectory_ID, levels = trajectory_order),
    Passage_factor = factor(Passage, levels = rev(sort(unique(Passage)))),
    Cell_label = ifelse(is.finite(Mean_deltaCt), sprintf("%.1f", Mean_deltaCt), "")
  )

heat_midpoint <- median(heatmap_grid$Minus_deltaCt, na.rm = TRUE)
heat_limit <- max(abs(heatmap_grid$Minus_deltaCt - heat_midpoint), na.rm = TRUE)

p_heatmap <- ggplot(
  heatmap_grid,
  aes(x = Trajectory_ID, y = Passage_factor, fill = Minus_deltaCt)
) +
  geom_tile(colour = "white", linewidth = 0.30) +
  geom_text(
    aes(label = Cell_label),
    family = font_family,
    size = 5.0 / 2.845276,
    colour = "black"
  ) +
  facet_grid(
    cols = vars(Organoid_Group),
    scales = "free_x",
    space = "free_x",
    labeller = as_labeller(facet_labels)
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = heat_midpoint,
    limits = c(heat_midpoint - heat_limit, heat_midpoint + heat_limit),
    oob = scales::squish,
    na.value = "grey92",
    name = expression(-Delta * Ct)
  ) +
  labs(
    title = "MUC5AC expression across individual organoid trajectories",
    subtitle = "Cell values show mean deltaCt after technical-replicate QC",
    x = NULL,
    y = "Passage"
  ) +
  theme_bw(base_family = font_family, base_size = base_font_pt) +
  theme(
    text = element_text(family = font_family, colour = "black"),
    plot.title = element_text(face = "bold", size = plot_title_pt, hjust = 0),
    plot.subtitle = element_text(size = 6, hjust = 0),
    panel.grid = element_blank(),
    panel.spacing.x = unit(1.5, "mm"),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.25),
    strip.text = element_text(size = 6.2, face = "bold"),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 5.2),
    axis.text.y = element_text(size = 6),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.margin = margin(5, 8, 8, 5, "mm")
  )

save_plot_set(
  p_heatmap,
  subdirs$figures,
  "MUC5AC_TaqMan_sample_by_passage_heatmap_SUPPLEMENTARY",
  183,
  115
)

# Mixed models and sensitivity analyses ----
fit_lme_safely <- function(fixed_formula, random_formula, data, model_name) {
  fit <- tryCatch(
    nlme::lme(
      fixed = fixed_formula,
      random = random_formula,
      data = data,
      method = "REML",
      na.action = stats::na.omit,
      control = nlme::lmeControl(
        returnObject = TRUE,
        opt = "optim",
        maxIter = 200,
        msMaxIter = 200
      )
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(
      fit = NULL,
      fixed = write_model_failure(model_name, fit),
      anova = data.frame(
        Model = model_name,
        Model_status = paste0("FAILED: ", conditionMessage(fit)),
        stringsAsFactors = FALSE
      )
    ))
  }

  fixed_table <- extract_lme_fixed(fit, model_name)
  fixed_table$Model_status <- "OK"

  anova_table <- as.data.frame(anova(fit))
  anova_table$Term <- rownames(anova_table)
  rownames(anova_table) <- NULL
  anova_table$Model <- model_name

  list(fit = fit, fixed = fixed_table, anova = anova_table)
}

primary_model <- fit_lme_safely(
  Minus_deltaCt ~ Passage_centered * Organoid_Group,
  ~1 | Trajectory_ID,
  analysis_primary,
  "Primary continuous-passage model: random intercept by trajectory"
)

strict_model <- fit_lme_safely(
  Minus_deltaCt ~ Passage_centered * Organoid_Group,
  ~1 | Trajectory_ID,
  analysis_strict,
  "Strict-QC continuous-passage model: random intercept by trajectory"
)

factor_passage_model <- fit_lme_safely(
  Minus_deltaCt ~ Passage_factor + Organoid_Group,
  ~1 | Trajectory_ID,
  analysis_primary,
  "Primary factor-passage sensitivity: random intercept by trajectory"
)

nested_patient_model <- fit_lme_safely(
  Minus_deltaCt ~ Passage_centered * Organoid_Group,
  ~1 | Patient_ID / Trajectory_ID,
  analysis_primary,
  "Patient/trajectory nested random-intercept sensitivity"
)

separate_group_models <- lapply(levels(analysis_primary$Organoid_Group), function(group_name) {
  dat_group <- analysis_primary %>%
    filter(Organoid_Group == group_name) %>%
    droplevels()

  fit_lme_safely(
    Minus_deltaCt ~ Passage_centered,
    ~1 | Trajectory_ID,
    dat_group,
    paste0("Separate passage slope: ", group_name)
  )
})

all_fixed_tables <- bind_rows(
  primary_model$fixed,
  strict_model$fixed,
  factor_passage_model$fixed,
  nested_patient_model$fixed,
  bind_rows(lapply(separate_group_models, `[[`, "fixed"))
)

all_anova_tables <- bind_rows(
  primary_model$anova,
  strict_model$anova,
  factor_passage_model$anova,
  nested_patient_model$anova,
  bind_rows(lapply(separate_group_models, `[[`, "anova"))
)

data.table::fwrite(
  all_fixed_tables,
  file.path(subdirs$statistics, "MUC5AC_mixed_effects_fixed_effects.csv")
)

data.table::fwrite(
  all_anova_tables,
  file.path(subdirs$statistics, "MUC5AC_mixed_effects_ANOVA_tables.csv")
)

model_summary_file <- file.path(subdirs$statistics, "MUC5AC_mixed_effects_model_summaries.txt")
summary_connection <- file(model_summary_file, open = "wt")

writeLines("MUC5AC MIXED-EFFECTS MODEL SUMMARIES\n", summary_connection)

model_objects <- list(
  Primary = primary_model$fit,
  Strict_QC = strict_model$fit,
  Factor_passage = factor_passage_model$fit,
  Nested_patient_trajectory = nested_patient_model$fit
)

for (model_name in names(model_objects)) {
  writeLines(paste0("\n\n===== ", model_name, " =====\n"), summary_connection)
  if (is.null(model_objects[[model_name]])) {
    writeLines("Model failed; see CSV audit.\n", summary_connection)
  } else {
    writeLines(capture.output(summary(model_objects[[model_name]])), summary_connection)
    writeLines("\nANOVA:\n", summary_connection)
    writeLines(capture.output(anova(model_objects[[model_name]])), summary_connection)
  }
}

for (i in seq_along(separate_group_models)) {
  group_name <- levels(analysis_primary$Organoid_Group)[i]
  writeLines(paste0("\n\n===== Separate slope: ", group_name, " =====\n"), summary_connection)
  if (is.null(separate_group_models[[i]]$fit)) {
    writeLines("Model failed; see CSV audit.\n", summary_connection)
  } else {
    writeLines(capture.output(summary(separate_group_models[[i]]$fit)), summary_connection)
  }
}

close(summary_connection)

# Methods, figure legend and session records ----
methods_text <- c(
  "MUC5AC expression was measured by TaqMan quantitative PCR and normalised to GAPDH using deltaCt = Ct(MUC5AC) - Ct(GAPDH).",
  "Technical reactions were performed in triplicate and were not treated as independent biological observations.",
  paste0(
    "A single technical replicate was removed only when its absolute distance from the triplicate median was at least ",
    outlier_min_distance_from_median,
    " deltaCt, this deviation was at least ",
    outlier_min_deviation_ratio,
    "-fold greater than the second-largest deviation, and the remaining pair differed by no more than ",
    outlier_max_remaining_pair_range,
    " deltaCt."
  ),
  "The mean of the retained technical replicates was calculated for each organoid line at each passage.",
  paste0(
    "Samples with a final technical-replicate range greater than ",
    high_dispersion_range_threshold,
    " deltaCt were retained but flagged in the primary analysis and excluded in a strict-QC sensitivity analysis."
  ),
  "GAPDH-normalised relative expression was represented as 2^(-deltaCt), and statistical analyses were performed using -deltaCt so that larger values indicated greater MUC5AC expression.",
  "Passage-associated changes were analysed using random-intercept linear mixed-effects models with repeated observations grouped by organoid trajectory."
)

writeLines(
  methods_text,
  file.path(output_dir, "THESIS_METHODS_MUC5AC_TAQMAN.txt")
)

figure_legend <- paste0(
  "TaqMan qPCR measurement of MUC5AC expression across serial passage of patient-derived colonic organoids. ",
  "Each thin line connects repeated samples from one organoid trajectory, and the thicker line shows the group geometric mean with SEM. ",
  "SSL and HP-derived organoids are shown together as SSL/HP-derived organoids, while matched adjacent non-lesional organoids are shown as H_ADJ_SSL_HP. ",
  "MUC5AC expression was normalised to GAPDH and is displayed as 2^(-deltaCt) on a logarithmic axis. ",
  "Technical triplicates were averaged after removal of clearly isolated technical outliers according to the pre-specified QC rule. ",
  "Technical replicates were not treated as independent biological observations."
)

writeLines(
  figure_legend,
  file.path(output_dir, "THESIS_FIGURE_LEGEND_MUC5AC_TAQMAN.txt")
)

run_notes <- c(
  "MUC5AC ORGANOID TAQMAN ANALYSIS: RUN NOTES",
  "===========================================",
  "",
  paste0("Run date: ", Sys.time()),
  paste0("Input file: ", input_file),
  paste0("Output directory: ", output_dir),
  "",
  "CSV import:",
  paste0("  Rows read: ", import_audit$Rows_read_from_CSV),
  paste0("  Genuine sample rows retained: ", import_audit$Valid_rows_with_Sample_Name_beginning_GM),
  paste0("  Blank/malformed rows discarded: ", import_audit$Discarded_rows_total),
  paste0("  Non-empty malformed rows among discarded rows: ", import_audit$Discarded_nonempty_malformed_rows),
  "",
  "Technical replicate QC:",
  paste0("  Clear outlier wells removed: ", sum(sample_qc$Outlier_removed)),
  paste0("  High-dispersion samples flagged: ", sum(sample_qc$High_dispersion_flag)),
  paste0("  Samples in primary analysis: ", nrow(analysis_primary)),
  paste0("  Samples in strict-QC sensitivity: ", nrow(analysis_strict)),
  paste0("  Organoid trajectories in primary analysis: ", n_distinct(analysis_primary$Trajectory_ID)),
  "",
  "Primary figure:",
  "  03_Figures/MUC5AC_TaqMan_longitudinal_trajectories_FINAL.*",
  "",
  "Before thesis submission:",
  "  Confirm every automatically removed well against the original amplification curve and raw MUC5AC/GAPDH Ct values."
)

writeLines(run_notes, file.path(output_dir, "RUN_NOTES.txt"))

writeLines(
  capture.output(sessionInfo()),
  file.path(subdirs$qc, "SESSION_INFO.txt")
)

message("\nDONE: final MUC5AC organoid TaqMan analysis complete.")
message("Output directory:\n", output_dir)
message("\nCSV rows read: ", import_audit$Rows_read_from_CSV)
message("Genuine sample rows retained: ", nrow(input))
message("Primary analysed samples: ", nrow(analysis_primary))
message("Strict-QC analysed samples: ", nrow(analysis_strict))
message("Clear technical outlier wells removed: ", sum(sample_qc$Outlier_removed))
message("High-dispersion samples flagged: ", sum(sample_qc$High_dispersion_flag))
message("\nInspect first:")
message("  01_QC/MUC5AC_CSV_import_audit.csv")
message("  01_QC/MUC5AC_technical_replicate_QC_COMPLETE.csv")
message("  03_Figures/MUC5AC_TaqMan_longitudinal_trajectories_FINAL.*")
message("  04_Statistics/MUC5AC_mixed_effects_fixed_effects.csv")
