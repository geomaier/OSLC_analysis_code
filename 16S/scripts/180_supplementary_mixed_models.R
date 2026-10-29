#!/usr/bin/env Rscript
# Supplementary alpha-diversity tables using a common pre-prevalence read-depth covariate.
(function() {
set.seed(1)

options(stringsAsFactors = FALSE)
options(contrasts = c("contr.treatment", "contr.poly"))

required_packages <- c(
  "phyloseq",
  "dplyr",
  "tidyr",
  "tibble",
  "lme4",
  "lmerTest"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following required package(s) before running: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(lme4)
  library(lmerTest)
})

config <- new.env(parent = baseenv())
if (file.exists("config/config.R")) sys.source("config/config.R", envir = config)
work_dir <- if (exists("PATH_WORK", config, inherits = FALSE)) config$PATH_WORK else "outputs"
if (!dir.exists(work_dir)) stop("Input directory does not exist: ", work_dir)
base_dir <- work_dir

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
  "Supplementary_alpha_COMMON_READS",
  paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
)

if (dir.exists(out_dir) || !dir.create(out_dir, recursive = TRUE)) stop("Cannot create output directory.")

patient_col  <- "Patient"
age_col      <- "AGE"
gender_col   <- "GENDER"
location_col <- "Location_Colon"

strict_expected_counts <- TRUE

expected_group_counts <- c(
  SSL_HP = 38L,
  H_ADJ_SSL_HP = 12L,
  H_COL = 14L
)

expected_group_patients <- c(
  SSL_HP = 23L,
  H_ADJ_SSL_HP = 10L,
  H_COL = 8L
)

group_levels <- c(
  "SSL_HP",
  "H_ADJ_SSL_HP",
  "H_COL"
)

metrics <- c(
  "Observed",
  "Shannon",
  "Simpson"
)

comparison_order <- c(
  "SSL_HP vs H_ADJ_SSL_HP",
  "SSL_HP vs H_COL",
  "H_ADJ_SSL_HP vs H_COL"
)

table_number_descriptive <- "ST4.2.1"
table_number_models      <- "ST4.2.2"
table_number_covariates  <- "ST4.2.3"

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "na", "N/A", "NULL", "None")] <- NA_character_
  x
}

normalise_ids <- function(x) {
  gsub("\\.", "-", clean_character(x))
}

make_group <- function(x) {
  x <- toupper(clean_character(x))
  x <- gsub("[^A-Z0-9]+", "_", x)

  out <- dplyr::case_when(
    x %in% c("SSL", "HP", "SSL_HP") ~ "SSL_HP",
    x %in% c("H_ADJ_SSL", "H_ADJ_HP", "H_ADJ_SSL_HP") ~ "H_ADJ_SSL_HP",
    x %in% c("H_COL", "H_COL_D", "H_COL_ND") ~ "H_COL",
    TRUE ~ NA_character_
  )

  factor(out, levels = group_levels)
}

extract_patient_from_sample_id <- function(sample_id) {
  patient <- sub(".*_([0-9]+)$", "\\1", as.character(sample_id))
  patient[patient == as.character(sample_id)] <- NA_character_
  patient
}

safe_p_adjust <- function(p, method = "holm") {
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  if (any(ok)) {
    out[ok] <- p.adjust(p[ok], method = method)
  }
  out
}

format_p <- function(x) {
  vapply(
    x,
    function(z) {
      if (!is.finite(z)) return("—")
      if (z < 0.001) return("<0.001")
      sprintf("%.3f", z)
    },
    character(1)
  )
}

format_num <- function(x, digits = 3L) {
  vapply(
    x,
    function(z) {
      if (!is.finite(z)) return("—")
      formatC(z, format = "f", digits = digits)
    },
    character(1)
  )
}

format_warning_text <- function(x) {
  x <- unique(clean_character(x))
  x <- x[!is.na(x)]
  if (length(x) == 0L) "" else paste(x, collapse = " | ")
}

prepare_physeq <- function(path, label) {

  if (!file.exists(path)) {
    stop("Missing phyloseq object: ", path, call. = FALSE)
  }

  phy <- readRDS(path)

  if (!inherits(phy, "phyloseq")) {
    stop(label, " is not a phyloseq object.", call. = FALSE)
  }

  if (
    is.null(otu_table(phy, errorIfNULL = FALSE)) ||
    is.null(sample_data(phy, errorIfNULL = FALSE))
  ) {
    stop(label, " must contain counts and sample metadata.", call. = FALSE)
  }

  if (!"Sample_Type" %in% colnames(sample_data(phy))) {
    stop("Sample_Type is missing from ", label, ".", call. = FALSE)
  }

  if (any(grepl("7699", sample_names(phy), fixed = TRUE))) {
    stop(label, " still contains patient 7699.", call. = FALSE)
  }

  otu_check <- as(otu_table(phy), "matrix")

  if (any(!is.finite(otu_check)) || any(otu_check < 0)) {
    stop(label, " contains invalid abundance values.", call. = FALSE)
  }

  sample_data(phy)$Group <- make_group(sample_data(phy)$Sample_Type)

  phy <- subset_samples(phy, !is.na(Group))
  phy <- prune_taxa(taxa_sums(phy) > 0, phy)

  observed_counts <- table(
    factor(
      as.character(sample_data(phy)$Group),
      levels = names(expected_group_counts)
    )
  )

  if (
    strict_expected_counts &&
    !identical(
      as.integer(observed_counts),
      as.integer(expected_group_counts)
    )
  ) {
    stop(
      "Unexpected group counts in ", label, ". Expected ",
      paste(
        names(expected_group_counts),
        expected_group_counts,
        sep = "=",
        collapse = ", "
      ),
      "; observed ",
      paste(
        names(observed_counts),
        as.integer(observed_counts),
        sep = "=",
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  list(
    phy = phy,
    label = label,
    path = path
  )
}

unfiltered_object <- prepare_physeq(
  physeq_unfiltered_path,
  "No_prevalence_filter"
)

filtered_object <- prepare_physeq(
  physeq_5pct_path,
  "Filtered_5pct"
)

unfiltered_ids <- normalise_ids(sample_names(unfiltered_object$phy))
filtered_ids   <- normalise_ids(sample_names(filtered_object$phy))

sample_set_audit <- tibble(
  SampleID_normalised = sort(unique(c(unfiltered_ids, filtered_ids))),
  Present_in_unfiltered = SampleID_normalised %in% unfiltered_ids,
  Present_in_filtered = SampleID_normalised %in% filtered_ids
)

write.csv(
  sample_set_audit,
  file.path(out_dir, "00_sample_set_audit.csv"),
  row.names = FALSE
)

if (!setequal(unfiltered_ids, filtered_ids)) {
  stop(
    "The prevalence-unfiltered and 5%-filtered objects do not contain the same samples.",
    call. = FALSE
  )
}

build_metadata <- function(phy) {

  metadata <- data.frame(
    sample_data(phy),
    check.names = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    rownames_to_column("SampleID") %>%
    mutate(
      SampleID_norm = normalise_ids(SampleID),
      Group = make_group(Sample_Type)
    )

  required_columns <- c(
    patient_col,
    age_col,
    gender_col,
    location_col
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(metadata)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Missing metadata columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  patient <- clean_character(metadata[[patient_col]])
  patient_from_id <- extract_patient_from_sample_id(metadata$SampleID)

  patient[is.na(patient)] <- patient_from_id[is.na(patient)]

  if (anyNA(patient)) {
    stop("Patient IDs could not be resolved.", call. = FALSE)
  }

  out <- tibble(
    SampleID = metadata$SampleID,
    SampleID_norm = metadata$SampleID_norm,
    Group = factor(metadata$Group, levels = group_levels),
    Patient = factor(patient),
    Age = suppressWarnings(
      as.numeric(clean_character(metadata[[age_col]]))
    ),
    Gender = factor(clean_character(metadata[[gender_col]])),
    Location = factor(clean_character(metadata[[location_col]]))
  )

  if (any(out$Age < 0 | out$Age > 120, na.rm = TRUE)) {
    stop("Age contains implausible values.", call. = FALSE)
  }

  out
}

metadata_unfiltered <- build_metadata(unfiltered_object$phy)
metadata_filtered   <- build_metadata(filtered_object$phy)

metadata_crosscheck <- metadata_unfiltered %>%
  select(
    SampleID_norm,
    Group_unfiltered = Group,
    Patient_unfiltered = Patient
  ) %>%
  left_join(
    metadata_filtered %>%
      select(
        SampleID_norm,
        Group_filtered = Group,
        Patient_filtered = Patient
      ),
    by = "SampleID_norm"
  ) %>%
  mutate(
    Group_matches =
      as.character(Group_unfiltered) ==
      as.character(Group_filtered),
    Patient_matches =
      as.character(Patient_unfiltered) ==
      as.character(Patient_filtered)
  )

write.csv(
  metadata_crosscheck,
  file.path(out_dir, "00_metadata_crosscheck.csv"),
  row.names = FALSE
)

if (
  any(!metadata_crosscheck$Group_matches) ||
  any(!metadata_crosscheck$Patient_matches)
) {
  stop(
    "Group or Patient metadata differ between the two objects.",
    call. = FALSE
  )
}

patient_counts <- metadata_unfiltered %>%
  distinct(Group, Patient) %>%
  count(Group, name = "N_patients")

patient_counts_named <- setNames(
  patient_counts$N_patients,
  as.character(patient_counts$Group)
)

if (strict_expected_counts) {
  observed_patient_counts <- as.integer(
    patient_counts_named[names(expected_group_patients)]
  )

  if (
    anyNA(observed_patient_counts) ||
    !identical(
      observed_patient_counts,
      as.integer(expected_group_patients)
    )
  ) {
    stop(
      "Unexpected patient counts. Expected ",
      paste(
        names(expected_group_patients),
        expected_group_patients,
        sep = "=",
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }
}

reads_common <- tibble(
  SampleID_norm = normalise_ids(
    sample_names(unfiltered_object$phy)
  ),
  Reads = as.numeric(
    sample_sums(unfiltered_object$phy)
  )
)

if (anyDuplicated(reads_common$SampleID_norm)) {
  stop(
    "Duplicated normalised sample IDs were found in the common read-depth table.",
    call. = FALSE
  )
}

if (any(!is.finite(reads_common$Reads) | reads_common$Reads <= 0)) {
  stop(
    "The common post-rarefaction read-depth table contains invalid values.",
    call. = FALSE
  )
}

write.csv(
  reads_common,
  file.path(out_dir, "00_COMMON_READ_DEPTH_COVARIATE.csv"),
  row.names = FALSE
)

build_alpha_table <- function(
    object,
    metadata,
    dataset_label,
    analysis_role
) {

  alpha_table <- estimate_richness(
    object$phy,
    measures = metrics
  ) %>%
    rownames_to_column("SampleID_alpha") %>%
    mutate(
      SampleID_norm = normalise_ids(SampleID_alpha)
    )

  out <- alpha_table %>%
    left_join(
      reads_common,
      by = "SampleID_norm"
    ) %>%
    left_join(
      metadata %>%
        select(
          SampleID_norm,
          SampleID,
          Group,
          Patient,
          Age,
          Gender,
          Location
        ),
      by = "SampleID_norm"
    ) %>%
    mutate(
      Dataset = dataset_label,
      Analysis_role = analysis_role,
      Group = factor(Group, levels = group_levels),
      Patient = droplevels(factor(Patient))
    )

  required <- c(
    "SampleID",
    "Group",
    "Patient",
    "Reads",
    metrics
  )

  bad <- !complete.cases(
    out[, required, drop = FALSE]
  )

  if (any(bad)) {
    write.csv(
      out[bad, , drop = FALSE],
      file.path(
        out_dir,
        paste0(
          "00_missing_alpha_data_",
          analysis_role,
          ".csv"
        )
      ),
      row.names = FALSE
    )

    stop(
      "Missing required alpha-diversity data in ",
      dataset_label,
      ".",
      call. = FALSE
    )
  }

  if (nrow(out) != nsamples(object$phy)) {
    stop(
      "Alpha-diversity sample alignment failed for ",
      dataset_label,
      ".",
      call. = FALSE
    )
  }

  out
}

alpha_unfiltered <- build_alpha_table(
  unfiltered_object,
  metadata_unfiltered,
  dataset_label = "No prevalence filter",
  analysis_role = "PRIMARY"
)

alpha_filtered <- build_alpha_table(
  filtered_object,
  metadata_filtered,
  dataset_label = "5% prevalence filter",
  analysis_role = "SENSITIVITY"
)

write.csv(
  alpha_unfiltered,
  file.path(
    out_dir,
    "Alpha_input_PRIMARY_no_prevalence_filter_COMMON_READS.csv"
  ),
  row.names = FALSE
)

write.csv(
  alpha_filtered,
  file.path(
    out_dir,
    "Alpha_input_SENSITIVITY_5pct_filter_COMMON_READS.csv"
  ),
  row.names = FALSE
)

reads_check <- alpha_unfiltered %>%
  select(
    SampleID_norm,
    Reads_primary = Reads
  ) %>%
  left_join(
    alpha_filtered %>%
      select(
        SampleID_norm,
        Reads_sensitivity = Reads
      ),
    by = "SampleID_norm"
  ) %>%
  mutate(
    Reads_identical =
      Reads_primary == Reads_sensitivity
  )

write.csv(
  reads_check,
  file.path(
    out_dir,
    "00_READ_DEPTH_IDENTICAL_BETWEEN_ANALYSES_CHECK.csv"
  ),
  row.names = FALSE
)

if (any(!reads_check$Reads_identical)) {
  stop(
    "The sequencing-depth covariate is not identical between analyses.",
    call. = FALSE
  )
}

fit_lmer_captured <- function(formula_object, data) {

  warnings <- character(0)

  fit <- tryCatch(
    withCallingHandlers(
      lmerTest::lmer(
        formula_object,
        data = data,
        REML = FALSE,
        control = lme4::lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 200000)
        )
      ),
      warning = function(w) {
        warnings <<- c(
          warnings,
          conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      structure(
        list(message = conditionMessage(e)),
        class = "fit_error"
      )
    }
  )

  list(
    fit = fit,
    warnings = unique(warnings)
  )
}

extract_lrt <- function(reduced, full) {

  if (
    inherits(reduced, "fit_error") ||
    inherits(full, "fit_error")
  ) {
    return(
      c(
        Chi_square = NA_real_,
        DF = NA_real_,
        P_value = NA_real_
      )
    )
  }

  tab <- tryCatch(
    anova(reduced, full),
    error = function(e) NULL
  )

  if (is.null(tab)) {
    return(
      c(
        Chi_square = NA_real_,
        DF = NA_real_,
        P_value = NA_real_
      )
    )
  }

  p_col <- grep(
    "^Pr\\(",
    colnames(tab),
    value = TRUE
  )

  chi_col <- grep(
    "^Chisq$|^L\\.Ratio$",
    colnames(tab),
    value = TRUE
  )

  df_col <- grep(
    "^Chi Df$|^Df$",
    colnames(tab),
    value = TRUE
  )

  c(
    Chi_square =
      if (length(chi_col)) {
        as.numeric(tab[nrow(tab), chi_col[1]])
      } else {
        NA_real_
      },
    DF =
      if (length(df_col)) {
        as.numeric(tab[nrow(tab), df_col[1]])
      } else {
        NA_real_
      },
    P_value =
      if (length(p_col)) {
        as.numeric(tab[nrow(tab), p_col[1]])
      } else {
        NA_real_
      }
  )
}

random_effect_summary <- function(fit) {

  if (inherits(fit, "fit_error")) {
    return(
      c(
        Patient_variance = NA_real_,
        Residual_variance = NA_real_,
        ICC = NA_real_
      )
    )
  }

  variance_table <- as.data.frame(
    lme4::VarCorr(fit)
  )

  patient_variance <- variance_table$vcov[
    variance_table$grp == "Patient" &
      is.na(variance_table$var2)
  ]

  residual_variance <- variance_table$vcov[
    variance_table$grp == "Residual"
  ]

  if (
    length(patient_variance) != 1L ||
    length(residual_variance) != 1L
  ) {
    return(
      c(
        Patient_variance = NA_real_,
        Residual_variance = NA_real_,
        ICC = NA_real_
      )
    )
  }

  c(
    Patient_variance = patient_variance,
    Residual_variance = residual_variance,
    ICC =
      patient_variance /
      (patient_variance + residual_variance)
  )
}

fit_alpha_base <- function(
    data,
    metric,
    analysis_role
) {

  needed <- c(
    metric,
    "Group",
    "Reads",
    "Patient"
  )

  model_data <- data[
    complete.cases(
      data[, needed, drop = FALSE]
    ),
    needed,
    drop = FALSE
  ]

  model_data$Group <- droplevels(
    factor(
      model_data$Group,
      levels = group_levels
    )
  )

  model_data$Patient <- droplevels(
    factor(model_data$Patient)
  )

  model_data$Reads <- as.numeric(
    model_data$Reads
  )

  full_formula <- as.formula(
    paste0(
      metric,
      " ~ Group + Reads + (1 | Patient)"
    )
  )

  no_group_formula <- as.formula(
    paste0(
      metric,
      " ~ Reads + (1 | Patient)"
    )
  )

  no_reads_formula <- as.formula(
    paste0(
      metric,
      " ~ Group + (1 | Patient)"
    )
  )

  full_object <- fit_lmer_captured(
    full_formula,
    model_data
  )

  no_group_object <- fit_lmer_captured(
    no_group_formula,
    model_data
  )

  no_reads_object <- fit_lmer_captured(
    no_reads_formula,
    model_data
  )

  group_test <- extract_lrt(
    no_group_object$fit,
    full_object$fit
  )

  reads_test <- extract_lrt(
    no_reads_object$fit,
    full_object$fit
  )

  random_summary <- random_effect_summary(
    full_object$fit
  )

  data.frame(
    Analysis_role = analysis_role,
    Dataset = unique(data$Dataset),
    Metric = metric,
    N_samples = nrow(model_data),
    N_patients = nlevels(model_data$Patient),
    Read_depth_term =
      "Raw post-rarefaction sample sum from the no-prevalence-filter object",
    Model = paste(
      deparse(full_formula),
      collapse = ""
    ),
    Group_LRT_Chi_square =
      group_test["Chi_square"],
    Group_LRT_DF =
      group_test["DF"],
    Group_LRT_P =
      group_test["P_value"],
    Reads_LRT_Chi_square =
      reads_test["Chi_square"],
    Reads_LRT_DF =
      reads_test["DF"],
    Reads_LRT_P =
      reads_test["P_value"],
    Patient_variance =
      random_summary["Patient_variance"],
    Residual_variance =
      random_summary["Residual_variance"],
    ICC =
      random_summary["ICC"],
    Singular =
      if (inherits(full_object$fit, "fit_error")) {
        NA
      } else {
        lme4::isSingular(
          full_object$fit,
          tol = 1e-4
        )
      },
    Warning_messages =
      format_warning_text(
        c(
          full_object$warnings,
          no_group_object$warnings,
          no_reads_object$warnings
        )
      ),
    stringsAsFactors = FALSE
  )
}

primary_results <- bind_rows(
  lapply(
    metrics,
    function(metric) {
      fit_alpha_base(
        alpha_unfiltered,
        metric,
        "PRIMARY_no_prevalence_filter"
      )
    }
  )
)

sensitivity_results <- bind_rows(
  lapply(
    metrics,
    function(metric) {
      fit_alpha_base(
        alpha_filtered,
        metric,
        "SENSITIVITY_5pct_prevalence_filter"
      )
    }
  )
)

alpha_model_results <- bind_rows(
  primary_results,
  sensitivity_results
)

if (any(!is.finite(alpha_model_results$Group_LRT_P))) {
  stop(
    "At least one alpha-diversity Group LRT P value is missing.",
    call. = FALSE
  )
}

write.csv(
  alpha_model_results,
  file.path(
    out_dir,
    "FINAL_alpha_overall_mixed_models_COMMON_READS.csv"
  ),
  row.names = FALSE
)

fit_reference_model <- function(
    data,
    metric,
    reference_group
) {

  model_data <- data %>%
    select(
      all_of(metric),
      Group,
      Reads,
      Patient
    ) %>%
    filter(
      complete.cases(.)
    ) %>%
    mutate(
      Group = relevel(
        factor(
          Group,
          levels = group_levels
        ),
        ref = reference_group
      ),
      Patient = droplevels(
        factor(Patient)
      ),
      Reads = as.numeric(Reads)
    )

  formula_object <- as.formula(
    paste0(
      metric,
      " ~ Group + Reads + (1 | Patient)"
    )
  )

  fit_object <- fit_lmer_captured(
    formula_object,
    model_data
  )

  list(
    fit = fit_object$fit,
    warnings = fit_object$warnings,
    data = model_data,
    formula = paste(
      deparse(formula_object),
      collapse = ""
    )
  )
}

extract_fixed_effect <- function(
    fit,
    coefficient
) {

  if (inherits(fit, "fit_error")) {
    return(
      c(
        Estimate = NA_real_,
        Std_Error = NA_real_,
        DF = NA_real_,
        T = NA_real_,
        P = NA_real_
      )
    )
  }

  coef_table <- coef(
    summary(fit)
  )

  if (!coefficient %in% rownames(coef_table)) {
    return(
      c(
        Estimate = NA_real_,
        Std_Error = NA_real_,
        DF = NA_real_,
        T = NA_real_,
        P = NA_real_
      )
    )
  }

  p_col <- grep(
    "^Pr\\(",
    colnames(coef_table),
    value = TRUE
  )

  df_col <- intersect(
    c("df", "DF"),
    colnames(coef_table)
  )

  t_col <- intersect(
    c("t value", "t.value"),
    colnames(coef_table)
  )

  c(
    Estimate =
      as.numeric(
        coef_table[
          coefficient,
          "Estimate"
        ]
      ),
    Std_Error =
      as.numeric(
        coef_table[
          coefficient,
          "Std. Error"
        ]
      ),
    DF =
      if (length(df_col)) {
        as.numeric(
          coef_table[
            coefficient,
            df_col[1]
          ]
        )
      } else {
        NA_real_
      },
    T =
      if (length(t_col)) {
        as.numeric(
          coef_table[
            coefficient,
            t_col[1]
          ]
        )
      } else {
        NA_real_
      },
    P =
      if (length(p_col)) {
        as.numeric(
          coef_table[
            coefficient,
            p_col[1]
          ]
        )
      } else {
        NA_real_
      }
  )
}

run_pairwise_for_dataset <- function(
    data,
    dataset_label
) {

  output <- list()

  for (metric in metrics) {

    model_hcol <- fit_reference_model(
      data,
      metric,
      "H_COL"
    )

    model_hadj <- fit_reference_model(
      data,
      metric,
      "H_ADJ_SSL_HP"
    )

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

    result <- data.frame(
      Dataset = dataset_label,
      Metric = metric,
      Comparison = comparison_order,
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
      T = c(
        ssl_vs_adj["T"],
        ssl_vs_hcol["T"],
        adj_vs_hcol["T"]
      ),
      P_raw = c(
        ssl_vs_adj["P"],
        ssl_vs_hcol["P"],
        adj_vs_hcol["P"]
      ),
      stringsAsFactors = FALSE
    )

    result$P_Holm <- p.adjust(
      result$P_raw,
      method = "holm"
    )

    output[[metric]] <- result
  }

  bind_rows(output)
}

pairwise_primary <- run_pairwise_for_dataset(
  alpha_unfiltered,
  "No prevalence filter (primary)"
)

pairwise_sensitivity <- run_pairwise_for_dataset(
  alpha_filtered,
  "5% prevalence filter (sensitivity)"
)

pairwise_results <- bind_rows(
  pairwise_primary,
  pairwise_sensitivity
)

write.csv(
  pairwise_results,
  file.path(
    out_dir,
    "FINAL_alpha_pairwise_comparisons_Holm_COMMON_READS.csv"
  ),
  row.names = FALSE
)

fit_alpha_covariate <- function(
    data,
    metric,
    covariate
) {

  needed <- c(
    metric,
    "Group",
    "Reads",
    "Patient",
    covariate
  )

  model_data <- data[
    complete.cases(
      data[, needed, drop = FALSE]
    ),
    needed,
    drop = FALSE
  ]

  model_data$Group <- droplevels(
    factor(
      model_data$Group,
      levels = group_levels
    )
  )

  model_data$Patient <- droplevels(
    factor(model_data$Patient)
  )

  model_data$Reads <- as.numeric(
    model_data$Reads
  )

  if (covariate %in% c("Gender", "Location")) {
    model_data[[covariate]] <- droplevels(
      factor(model_data[[covariate]])
    )
  }

  covariate_available <- if (covariate == "Age") {
    length(unique(model_data[[covariate]])) >= 2L
  } else {
    nlevels(model_data[[covariate]]) >= 2L
  }

  if (
    nrow(model_data) == 0L ||
    nlevels(model_data$Group) < 2L ||
    !covariate_available
  ) {
    return(
      data.frame(
        Metric = metric,
        Covariate = covariate,
        N_samples = nrow(model_data),
        N_patients = nlevels(model_data$Patient),
        Fixed_design_rank_deficient = NA,
        Group_P_base_same_rows = NA_real_,
        Group_P_adjusted = NA_real_,
        Reads_P_adjusted = NA_real_,
        Covariate_P = NA_real_,
        Singular_adjusted = NA,
        Warning_messages = "",
        Note = "Covariate unavailable or invariant",
        stringsAsFactors = FALSE
      )
    )
  }

  fixed_matrix <- model.matrix(
    stats::reformulate(
      c(
        "Group",
        "Reads",
        covariate
      )
    ),
    data = model_data
  )

  rank_deficient <-
    qr(fixed_matrix)$rank <
    ncol(fixed_matrix)

  if (rank_deficient) {
    return(
      data.frame(
        Metric = metric,
        Covariate = covariate,
        N_samples = nrow(model_data),
        N_patients = nlevels(model_data$Patient),
        Fixed_design_rank_deficient = TRUE,
        Group_P_base_same_rows = NA_real_,
        Group_P_adjusted = NA_real_,
        Reads_P_adjusted = NA_real_,
        Covariate_P = NA_real_,
        Singular_adjusted = NA,
        Warning_messages = "",
        Note =
          "Group and covariate are not jointly estimable in these complete-case data",
        stringsAsFactors = FALSE
      )
    )
  }

  base_formula <- as.formula(
    paste0(
      metric,
      " ~ Group + Reads + (1 | Patient)"
    )
  )

  base_no_group_formula <- as.formula(
    paste0(
      metric,
      " ~ Reads + (1 | Patient)"
    )
  )

  adjusted_formula <- as.formula(
    paste0(
      metric,
      " ~ Group + Reads + ",
      covariate,
      " + (1 | Patient)"
    )
  )

  adjusted_no_group_formula <- as.formula(
    paste0(
      metric,
      " ~ Reads + ",
      covariate,
      " + (1 | Patient)"
    )
  )

  adjusted_no_reads_formula <- as.formula(
    paste0(
      metric,
      " ~ Group + ",
      covariate,
      " + (1 | Patient)"
    )
  )

  base_fit <- fit_lmer_captured(
    base_formula,
    model_data
  )

  base_no_group_fit <- fit_lmer_captured(
    base_no_group_formula,
    model_data
  )

  adjusted_fit <- fit_lmer_captured(
    adjusted_formula,
    model_data
  )

  adjusted_no_group_fit <- fit_lmer_captured(
    adjusted_no_group_formula,
    model_data
  )

  adjusted_no_reads_fit <- fit_lmer_captured(
    adjusted_no_reads_formula,
    model_data
  )

  base_group_test <- extract_lrt(
    base_no_group_fit$fit,
    base_fit$fit
  )

  adjusted_group_test <- extract_lrt(
    adjusted_no_group_fit$fit,
    adjusted_fit$fit
  )

  reads_test <- extract_lrt(
    adjusted_no_reads_fit$fit,
    adjusted_fit$fit
  )

  covariate_test <- extract_lrt(
    base_fit$fit,
    adjusted_fit$fit
  )

  data.frame(
    Metric = metric,
    Covariate = covariate,
    N_samples = nrow(model_data),
    N_patients = nlevels(model_data$Patient),
    Fixed_design_rank_deficient = FALSE,
    Group_P_base_same_rows =
      base_group_test["P_value"],
    Group_P_adjusted =
      adjusted_group_test["P_value"],
    Reads_P_adjusted =
      reads_test["P_value"],
    Covariate_P =
      covariate_test["P_value"],
    Singular_adjusted =
      if (inherits(adjusted_fit$fit, "fit_error")) {
        NA
      } else {
        lme4::isSingular(
          adjusted_fit$fit,
          tol = 1e-4
        )
      },
    Warning_messages =
      format_warning_text(
        c(
          base_fit$warnings,
          base_no_group_fit$warnings,
          adjusted_fit$warnings,
          adjusted_no_group_fit$warnings,
          adjusted_no_reads_fit$warnings
        )
      ),
    Note =
      if (inherits(adjusted_fit$fit, "fit_error")) {
        paste0(
          "Model failed: ",
          adjusted_fit$fit$message
        )
      } else {
        ""
      },
    stringsAsFactors = FALSE
  )
}

alpha_covariate_results <- bind_rows(
  lapply(
    metrics,
    function(metric) {
      bind_rows(
        lapply(
          c(
            "Age",
            "Gender",
            "Location"
          ),
          function(covariate) {
            fit_alpha_covariate(
              alpha_unfiltered,
              metric,
              covariate
            )
          }
        )
      )
    }
  )
) %>%
  group_by(Metric) %>%
  mutate(
    Covariate_P_BH_within_metric =
      safe_p_adjust(
        Covariate_P,
        method = "BH"
      )
  ) %>%
  ungroup()

write.csv(
  alpha_covariate_results,
  file.path(
    out_dir,
    "FINAL_alpha_covariate_sensitivity_PRIMARY_unfiltered_COMMON_READS.csv"
  ),
  row.names = FALSE
)

alpha_descriptive <- bind_rows(
  alpha_unfiltered,
  alpha_filtered
) %>%
  group_by(
    Analysis_role,
    Dataset,
    Group
  ) %>%
  summarise(
    N_samples = n(),
    N_patients = n_distinct(Patient),

    Reads_median = median(Reads),
    Reads_Q1 = quantile(Reads, 0.25),
    Reads_Q3 = quantile(Reads, 0.75),

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
  alpha_descriptive,
  file.path(
    out_dir,
    "FINAL_alpha_descriptive_statistics_COMMON_READS.csv"
  ),
  row.names = FALSE
)

alpha_descriptive_publication <- alpha_descriptive %>%
  mutate(
    Dataset = case_when(
      Analysis_role == "PRIMARY" ~
        "Before 5% prevalence filter (primary)",
      Analysis_role == "SENSITIVITY" ~
        "5% prevalence-filtered (sensitivity)",
      TRUE ~ Dataset
    ),

    Group = as.character(Group),

    `Observed richness, median (IQR)` = sprintf(
      "%.1f (%.1f–%.1f)",
      Observed_median,
      Observed_Q1,
      Observed_Q3
    ),

    `Shannon, median (IQR)` = sprintf(
      "%.3f (%.3f–%.3f)",
      Shannon_median,
      Shannon_Q1,
      Shannon_Q3
    ),

    `Simpson, median (IQR)` = sprintf(
      "%.3f (%.3f–%.3f)",
      Simpson_median,
      Simpson_Q1,
      Simpson_Q3
    )
  ) %>%
  select(
    Dataset,
    Group,
    `n samples` = N_samples,
    `n patients` = N_patients,
    `Observed richness, median (IQR)`,
    `Shannon, median (IQR)`,
    `Simpson, median (IQR)`
  )

alpha_models_publication <- alpha_model_results %>%
  mutate(
    Dataset = case_when(
      Analysis_role ==
        "PRIMARY_no_prevalence_filter" ~
        "Before 5% prevalence filter (primary)",

      Analysis_role ==
        "SENSITIVITY_5pct_prevalence_filter" ~
        "5% prevalence-filtered (sensitivity)",

      TRUE ~ Dataset
    ),

    `Group χ²` =
      format_num(
        Group_LRT_Chi_square,
        3
      ),

    `Group df` =
      as.character(
        as.integer(
          Group_LRT_DF
        )
      ),

    `Group P` =
      format_p(
        Group_LRT_P
      ),

    `Reads χ²` =
      format_num(
        Reads_LRT_Chi_square,
        3
      ),

    `Reads df` =
      as.character(
        as.integer(
          Reads_LRT_DF
        )
      ),

    `Reads P` =
      format_p(
        Reads_LRT_P
      ),

    `ICC` =
      format_num(
        ICC,
        3
      )
  ) %>%
  select(
    Dataset,
    Metric,
    `n samples` = N_samples,
    `n patients` = N_patients,
    `Group χ²`,
    `Group df`,
    `Group P`,
    `Reads χ²`,
    `Reads df`,
    `Reads P`,
    `ICC`
  )

alpha_covariate_publication <- alpha_covariate_results %>%
  mutate(
    `Group P, base` =
      format_p(
        Group_P_base_same_rows
      ),

    `Group P, adjusted` =
      format_p(
        Group_P_adjusted
      ),

    `Reads P, adjusted` =
      format_p(
        Reads_P_adjusted
      ),

    `Covariate P` =
      format_p(
        Covariate_P
      ),

    `Covariate BH FDR` =
      format_p(
        Covariate_P_BH_within_metric
      )
  ) %>%
  select(
    Metric,
    Covariate,
    `n samples` = N_samples,
    `n patients` = N_patients,
    `Group P, base`,
    `Group P, adjusted`,
    `Reads P, adjusted`,
    `Covariate P`,
    `Covariate BH FDR`
  )

write.csv(
  alpha_descriptive_publication,
  file.path(
    out_dir,
    "PUBLICATION_TABLE_alpha_descriptive_FINAL.csv"
  ),
  row.names = FALSE
)

write.csv(
  alpha_models_publication,
  file.path(
    out_dir,
    "PUBLICATION_TABLE_alpha_mixed_models_FINAL.csv"
  ),
  row.names = FALSE
)

write.csv(
  alpha_covariate_publication,
  file.path(
    out_dir,
    "PUBLICATION_TABLE_alpha_covariates_FINAL.csv"
  ),
  row.names = FALSE
)

word_packages <- c(
  "officer",
  "flextable"
)

word_packages_available <- all(
  vapply(
    word_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
)

if (word_packages_available) {

  suppressPackageStartupMessages({
    library(officer)
    library(flextable)
  })

  publication_flextable <- function(
    data,
    first_columns_left = 2L
  ) {

    ft <- flextable::flextable(data)

    ft <- flextable::font(
      ft,
      fontname = "Arial",
      part = "all"
    )

    ft <- flextable::fontsize(
      ft,
      size = 8,
      part = "all"
    )

    ft <- flextable::bold(
      ft,
      part = "header",
      bold = TRUE
    )

    ft <- flextable::color(
      ft,
      part = "header",
      color = "white"
    )

    ft <- flextable::bg(
      ft,
      part = "header",
      bg = "#A6A6A6"
    )

    ft <- flextable::align(
      ft,
      part = "header",
      align = "center"
    )

    ft <- flextable::align(
      ft,
      part = "body",
      align = "center"
    )

    left_columns <- seq_len(
      min(
        first_columns_left,
        ncol(data)
      )
    )

    ft <- flextable::align(
      ft,
      j = left_columns,
      part = "body",
      align = "left"
    )

    ft <- flextable::bold(
      ft,
      j = 1,
      part = "body",
      bold = TRUE
    )

    if (nrow(data) >= 2L) {
      ft <- flextable::bg(
        ft,
        i = seq(
          2L,
          nrow(data),
          by = 2L
        ),
        part = "body",
        bg = "#E7E6E6"
      )
    }

    thin_border <- officer::fp_border(
      color = "#BFBFBF",
      width = 0.5
    )

    strong_border <- officer::fp_border(
      color = "black",
      width = 0.8
    )

    ft <- flextable::border_remove(ft)

    ft <- flextable::hline_top(
      ft,
      border = strong_border,
      part = "header"
    )

    ft <- flextable::hline_bottom(
      ft,
      border = strong_border,
      part = "header"
    )

    ft <- flextable::hline(
      ft,
      border = thin_border,
      part = "body"
    )

    ft <- flextable::hline_bottom(
      ft,
      border = strong_border,
      part = "body"
    )

    ft <- flextable::valign(
      ft,
      valign = "center",
      part = "all"
    )

    ft <- flextable::padding(
      ft,
      padding = 2.5,
      part = "all"
    )

    ft <- flextable::autofit(ft)

    ft <- flextable::set_table_properties(
      ft,
      layout = "autofit",
      width = 1
    )

    ft
  }

  make_landscape_doc <- function() {

    doc <- officer::read_docx()

    section <- officer::prop_section(
      page_size =
        officer::page_size(
          orient = "landscape"
        ),
      page_margins =
        officer::page_mar(
          top = 0.55,
          bottom = 0.55,
          left = 0.55,
          right = 0.55
        )
    )

    officer::body_set_default_section(
      doc,
      section
    )
  }

  add_caption <- function(
    doc,
    number,
    title
  ) {

    officer::body_add_fpar(
      doc,
      officer::fpar(
        officer::ftext(
          paste0(
            "Supplementary Table ",
            number,
            ". "
          ),
          officer::fp_text(
            font.family = "Arial",
            font.size = 9,
            bold = TRUE,
            italic = TRUE
          )
        ),
        officer::ftext(
          title,
          officer::fp_text(
            font.family = "Arial",
            font.size = 9,
            italic = TRUE
          )
        )
      )
    )
  }

  add_note <- function(
    doc,
    text
  ) {

    officer::body_add_fpar(
      doc,
      officer::fpar(
        officer::ftext(
          text,
          officer::fp_text(
            font.family = "Arial",
            font.size = 8
          )
        )
      )
    )
  }

  ft_descriptive <- publication_flextable(
    alpha_descriptive_publication,
    2L
  )

  ft_models <- publication_flextable(
    alpha_models_publication,
    2L
  )

  ft_covariates <- publication_flextable(
    alpha_covariate_publication,
    2L
  )

  combined_doc <- make_landscape_doc()

  combined_doc <- add_caption(
    combined_doc,
    table_number_descriptive,
    paste0(
      "Descriptive alpha-diversity statistics ",
      "before and after the 5% prevalence filter."
    )
  )

  combined_doc <- flextable::body_add_flextable(
    combined_doc,
    ft_descriptive
  )

  combined_doc <- add_note(
    combined_doc,
    paste0(
      "Values are medians with interquartile ranges. ",
      "The post-rarefaction dataset before prevalence filtering is the ",
      "primary alpha-diversity analysis; the 5% prevalence-filtered ",
      "dataset is a sensitivity analysis."
    )
  )

  combined_doc <- officer::body_add_par(
    combined_doc,
    ""
  )

  combined_doc <- add_caption(
    combined_doc,
    table_number_models,
    paste0(
      "Linear mixed-effects models of alpha diversity ",
      "and sensitivity to prevalence filtering."
    )
  )

  combined_doc <- flextable::body_add_flextable(
    combined_doc,
    ft_models
  )

  combined_doc <- add_note(
    combined_doc,
    paste0(
      "Model: diversity metric ~ Group + Reads + (1 | Patient), ",
      "fitted by maximum likelihood. Reads is the raw post-rarefaction ",
      "sample sum from the prevalence-unfiltered object and is identical ",
      "in both analyses. Group and Reads P values are likelihood-ratio tests. ",
      "ICC is the proportion of total model variance attributable to Patient."
    )
  )

  combined_doc <- officer::body_add_par(
    combined_doc,
    ""
  )

  combined_doc <- add_caption(
    combined_doc,
    table_number_covariates,
    paste0(
      "Sensitivity of alpha-diversity group effects ",
      "to age, gender and anatomical location."
    )
  )

  combined_doc <- flextable::body_add_flextable(
    combined_doc,
    ft_covariates
  )

  combined_doc <- add_note(
    combined_doc,
    paste0(
      "The primary no-prevalence-filter model was used. ",
      "Each covariate was added separately using the same complete-case ",
      "samples: diversity metric ~ Group + Reads + covariate + (1 | Patient). ",
      "Reads is the common raw post-rarefaction sequencing-depth covariate. ",
      "Covariate BH FDR values were calculated across Age, Gender and Location ",
      "within each diversity metric."
    )
  )

  combined_docx <- file.path(
    out_dir,
    "Chapter_4_Supplementary_Alpha_Tables_FINAL_CORRECTED.docx"
  )

  print(
    combined_doc,
    target = combined_docx
  )

  cat(
    "\nEditable Word tables written to:\n",
    combined_docx,
    "\n",
    sep = ""
  )

} else {

  writeLines(
    c(
      "The statistical analyses completed successfully.",
      "",
      "Editable Word tables were not created because officer and/or flextable",
      "are not installed.",
      "",
      "Install them with:",
      "install.packages(c('officer', 'flextable'))",
      "",
      "Then rerun this script."
    ),
    file.path(
      out_dir,
      "WORD_TABLES_NOT_CREATED_INSTALL_PACKAGES.txt"
    )
  )
}

readme <- c(
  "FINAL CHAPTER 4 SUPPLEMENTARY ALPHA-DIVERSITY ANALYSIS",
  "======================================================",
  "",
  paste0(
    "Primary object: ",
    physeq_unfiltered_path
  ),
  paste0(
    "Sensitivity object: ",
    physeq_5pct_path
  ),
  "",
  "FINAL MODEL:",
  "  Metric ~ Group + Reads + (1 | Patient)",
  "",
  "CRITICAL READ-DEPTH RULE:",
  "  Reads comes from the post-rarefaction object before prevalence filtering.",
  "  The identical Reads value is used for each sample in BOTH analyses.",
  "  The 5% filter therefore changes the ASVs used to calculate diversity,",
  "  but does not redefine sequencing depth.",
  "",
  "DATASET ROLES:",
  "  PRIMARY: no prevalence filter.",
  "  SENSITIVITY: 5% prevalence filter.",
  "",
  "COVARIATE SENSITIVITY:",
  "  Age, Gender and Location are added one at a time to the PRIMARY model.",
  "  Covariate BH FDR is calculated across the three covariates within each metric.",
  "",
  "PAIRWISE COMPARISONS:",
  "  Fixed-effect coefficient tests after re-levelling Group.",
  "  Holm correction is applied across the three pairwise comparisons",
  "  separately within each alpha-diversity metric and prevalence version.",
  "",
  "MAIN OUTPUTS:",
  "  FINAL_alpha_overall_mixed_models_COMMON_READS.csv",
  "  FINAL_alpha_pairwise_comparisons_Holm_COMMON_READS.csv",
  "  FINAL_alpha_covariate_sensitivity_PRIMARY_unfiltered_COMMON_READS.csv",
  "  FINAL_alpha_descriptive_statistics_COMMON_READS.csv",
  "  PUBLICATION_TABLE_alpha_descriptive_FINAL.csv",
  "  PUBLICATION_TABLE_alpha_mixed_models_FINAL.csv",
  "  PUBLICATION_TABLE_alpha_covariates_FINAL.csv",
  "  Chapter_4_Supplementary_Alpha_Tables_FINAL_CORRECTED.docx",
  "",
  "IMPORTANT:",
  "  This script does not reproduce the old patient-averaged beta-diversity tables.",
  "  Those tables should not be used for the current final Section 4.3 analysis."
)

writeLines(
  readme,
  file.path(
    out_dir,
    "README_FIRST.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "sessionInfo.txt"
  )
)

cat("\n============================================================\n")
cat("FINAL ALPHA-DIVERSITY ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat("Results saved to:\n")
cat(out_dir, "\n\n")

cat("Overall mixed models:\n")
print(
  alpha_model_results[
    ,
    c(
      "Analysis_role",
      "Dataset",
      "Metric",
      "N_samples",
      "N_patients",
      "Group_LRT_P",
      "Reads_LRT_P",
      "ICC",
      "Singular"
    )
  ],
  row.names = FALSE
)

cat("\nAge / Gender / Location sensitivity:\n")
print(
  alpha_covariate_results[
    ,
    c(
      "Metric",
      "Covariate",
      "N_samples",
      "N_patients",
      "Group_P_base_same_rows",
      "Group_P_adjusted",
      "Reads_P_adjusted",
      "Covariate_P",
      "Covariate_P_BH_within_metric",
      "Singular_adjusted"
    )
  ],
  row.names = FALSE
)

cat("\nPairwise comparisons:\n")
print(
  pairwise_results[
    ,
    c(
      "Dataset",
      "Metric",
      "Comparison",
      "Estimate_first_minus_second",
      "P_raw",
      "P_Holm"
    )
  ],
  row.names = FALSE
)

cat("\nDONE.\n")

})()
