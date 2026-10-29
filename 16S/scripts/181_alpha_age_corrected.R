#!/usr/bin/env Rscript
# Corrected-age sensitivity analysis for supplementary table ST4.2.3.
# Set ALPHA_AGE_INPUT_DIR and ALPHA_AGE_OUTPUT_DIR; see ../README.md.
# Refits three age models and checks the original age results. Gender/Location
# raw P values are retained; BH is reapplied across covariates within each metric.
# Requires lme4 and lmerTest.
alpha_input_dir <- Sys.getenv("ALPHA_AGE_INPUT_DIR")
if (!nzchar(alpha_input_dir) || !dir.exists(alpha_input_dir)) {
  stop("Set ALPHA_AGE_INPUT_DIR to the private folder containing the three age-correction CSVs.")
}
output_parent <- Sys.getenv("ALPHA_AGE_OUTPUT_DIR")
if (!nzchar(output_parent)) stop("Set ALPHA_AGE_OUTPUT_DIR to an existing output parent directory.")

run_alpha_age_correction <- function(output_parent) {
  needed_packages <- c("lme4", "lmerTest")
  missing_packages <- needed_packages[!vapply(
    needed_packages, requireNamespace, logical(1), quietly = TRUE
  )]
  if (length(missing_packages)) {
    stop("Install the required packages first: ",
         paste(missing_packages, collapse = ", "), call. = FALSE)
  }
  previous_options <- options(
    stringsAsFactors = FALSE,
    contrasts = c("contr.treatment", "contr.poly")
  )
  on.exit(options(previous_options), add = TRUE)
  group_levels <- c("SSL_HP", "H_ADJ_SSL_HP", "H_COL")
  metrics <- c("Observed", "Shannon", "Simpson")

  if (!dir.exists(output_parent) && interactive() && .Platform$OS.type == "windows") {
    selected <- utils::choose.dir(
      default = getwd(),
      caption = "Select your existing alpha-diversity results folder"
    )
    if (is.na(selected)) stop("No output folder was selected.", call. = FALSE)
    output_parent <- selected
  }
  if (!dir.exists(output_parent)) {
    stop("Set output_parent near the top of this script to your existing alpha results folder.",
         call. = FALSE)
  }

  # Original primary input, exactly as supplied in the saved alpha run.
  original <- read.csv(file.path(alpha_input_dir, "Alpha_input_PRIMARY_no_prevalence_filter_COMMON_READS.csv"), colClasses = "character",
                       check.names = FALSE, stringsAsFactors = FALSE)
  for (column in c("Observed", "Shannon", "Simpson", "Reads", "Age")) {
    original[[column]] <- as.numeric(original[[column]])
  }

  # Corrected ages are calculated at biopsy and match the submitted Methods.
  # Historical HP/SSL sample-name aliases are recorded only for the join;
  # the original analysis SampleID, Patient and Group values are retained.
  age_key <- read.csv(file.path(alpha_input_dir, "alpha_corrected_age_key.csv"), colClasses = "character",
                     check.names = FALSE, stringsAsFactors = FALSE)
  age_key$Age_corrected <- as.numeric(age_key$Age_corrected)

  # Previous sensitivity output. Six Gender/Location rows supply the unchanged
  # raw P values required by the original within-metric BH correction.
  previous <- read.csv(file.path(alpha_input_dir, "FINAL_alpha_covariate_sensitivity_PRIMARY_unfiltered_COMMON_READS.csv"), check.names = FALSE,
                      stringsAsFactors = FALSE)

  required <- c("SampleID", "Patient", "Group", "Age", "Gender", "Location",
                "Reads", metrics)
  stopifnot(all(required %in% names(original)), nrow(original) == 64L,
            nrow(age_key) == 64L, !anyDuplicated(original$SampleID),
            !anyDuplicated(age_key$SampleID),
            setequal(original$SampleID, age_key$SampleID),
            length(unique(original$Patient)) == 31L,
            !anyNA(original[, required, drop = FALSE]))
  stopifnot(all(vapply(original[, c(metrics, "Reads", "Age")],
                      function(x) all(is.finite(x)), logical(1))))
  stopifnot(identical(as.integer(table(factor(original$Group, levels = group_levels))),
                      c(38L, 12L, 14L)))
  idx <- match(original$SampleID, age_key$SampleID)
  stopifnot(identical(original$Patient, age_key$Patient[idx]),
            all(is.finite(age_key$Age_corrected)),
            all(age_key$Age_corrected >= 0 & age_key$Age_corrected <= 120),
            all(age_key$Age_corrected == floor(age_key$Age_corrected)))
  corrected <- original
  corrected$Age <- age_key$Age_corrected[idx]
  other_columns <- setdiff(names(original), "Age")
  stopifnot(identical(original[, other_columns], corrected[, other_columns]))
  changes <- data.frame(
    SampleID = original$SampleID, Patient = original$Patient,
    Table_sample_name = age_key$Table_sample_name[idx],
    Age_previous = original$Age, Age_corrected = corrected$Age,
    Age_changed = original$Age != corrected$Age,
    stringsAsFactors = FALSE
  )
  stopifnot(sum(changes$Age_changed) == 43L)
  expected_keys <- as.vector(t(outer(metrics, c("Age", "Gender", "Location"), paste, sep = "|")))
  result_keys <- paste(previous$Metric, previous$Covariate, sep = "|")
  stopifnot(nrow(previous) == 9L, !anyDuplicated(result_keys),
            setequal(result_keys, expected_keys),
            all(previous$N_samples == 64), all(previous$N_patients == 31),
            all(is.finite(previous$Covariate_P)))

  output_stem <- paste0("AGE_CORRECTED_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
  out_dir <- file.path(output_parent, output_stem)
  suffix <- 1L
  while (dir.exists(out_dir)) {
    out_dir <- file.path(output_parent, paste0(output_stem, "_", suffix))
    suffix <- suffix + 1L
  }
  if (!dir.create(out_dir, recursive = FALSE)) stop("Could not create output folder.")
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  writeLines("IN PROGRESS. Use results only when COMPLETED.txt is present.",
             file.path(out_dir, "RUN_STATUS.txt"))
  write.csv(original, file.path(out_dir, "00_ORIGINAL_PRIMARY_INPUT.csv"), row.names = FALSE)
  write.csv(changes, file.path(out_dir, "00_AGE_CORRECTIONS_ALL_64.csv"), row.names = FALSE)
  write.csv(corrected, file.path(out_dir, "Alpha_input_PRIMARY_AGE_CORRECTED.csv"), row.names = FALSE)
  message("64 biopsies; 31 patients; 43 corrected ages. Output: ", out_dir)

  # Model helpers copied from Chapter4_alpha_final_corrected.R.
  clean_character <- function(x) {
    x <- trimws(as.character(x))
    x[x %in% c("", "NA", "na", "N/A", "NULL", "None")] <- NA_character_
    x
  }

  format_warning_text <- function(x) {
    x <- unique(clean_character(x))
    x <- x[!is.na(x)]
    if (length(x) == 0L) "" else paste(x, collapse = " | ")
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

  # Record every fit without changing the fitting function or model settings.
  original_fitter <- fit_lmer_captured
  fit_records <- list()
  fitted_models <- list()
  fit_stage <- "previous_age_reproduction"
  fit_lmer_captured <- function(formula_object, data) {
    object <- original_fitter(formula_object, data)
    fit_id <- paste0(fit_stage, "_", length(fit_records) + 1L)
    fit <- object$fit
    failed <- inherits(fit, "fit_error")
    conv_messages <- if (failed) character(0) else unlist(fit@optinfo$conv$lme4$messages)
    opt_codes <- if (failed) numeric(0) else unlist(fit@optinfo$conv$opt)
    critical_messages <- conv_messages[!grepl("boundary (singular)", conv_messages, fixed = TRUE)]
    converged <- !failed && !length(critical_messages) &&
      (!length(opt_codes) || all(opt_codes == 0))
    fit_records[[fit_id]] <<- data.frame(
      Fit_ID = fit_id, Stage = fit_stage,
      Formula = paste(deparse(formula_object), collapse = " "),
      N_samples = nrow(data), N_patients = length(unique(data$Patient)),
      Failed = failed, Converged = converged,
      Singular = if (failed) NA else lme4::isSingular(fit, tol = 1e-4),
      Log_likelihood = if (failed) NA_real_ else as.numeric(logLik(fit)),
      Convergence_messages = paste(conv_messages, collapse = " | "),
      Warnings = paste(object$warnings, collapse = " | "),
      Error = if (failed) fit$message else "",
      stringsAsFactors = FALSE
    )
    fitted_models[[fit_id]] <<- fit
    object
  }
  save_diagnostics <- function() {
    diagnostics <- do.call(rbind, fit_records)
    rownames(diagnostics) <- NULL
    write.csv(diagnostics, file.path(out_dir, "FIT_DIAGNOSTICS.csv"), row.names = FALSE)
    saveRDS(fitted_models, file.path(out_dir, "MODELS_AGE_CORRECTION.rds"))
    writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
    diagnostics
  }

  # Verify that the supplied original input reproduces its previous AGE rows.
  reproduced <- do.call(rbind, lapply(metrics, function(metric) {
    fit_alpha_covariate(original, metric, "Age")
  }))
  previous_age <- previous[previous$Covariate == "Age", , drop = FALSE]
  previous_age <- previous_age[match(metrics, previous_age$Metric), , drop = FALSE]
  p_columns <- c("Group_P_base_same_rows", "Group_P_adjusted", "Reads_P_adjusted", "Covariate_P")
  # Absolute P-value tolerance is only a numerical reproduction check, not a
  # threshold used for biological inference or for choosing an analysis.
  reproduction_tolerance <- 1e-4
  reproduction_check <- do.call(rbind, lapply(p_columns, function(column) {
    difference <- abs(reproduced[[column]] - previous_age[[column]])
    data.frame(Metric = metrics, Test = column,
               Saved_P = previous_age[[column]], Reproduced_P = reproduced[[column]],
               Absolute_difference = difference,
               Pass = is.finite(difference) & difference <= reproduction_tolerance,
               stringsAsFactors = FALSE)
  }))
  write.csv(reproduction_check, file.path(out_dir, "00_PREVIOUS_AGE_MODEL_REPRODUCTION.csv"),
            row.names = FALSE)
  diagnostics <- save_diagnostics()
  if (!all(reproduction_check$Pass) || any(diagnostics$Failed | !diagnostics$Converged)) {
    stop("The previous AGE models were not reproduced reliably. See the reproduction check and fit diagnostics in ",
         out_dir, ". No corrected supplementary table has been released.", call. = FALSE)
  }

  fit_stage <- "corrected_age"
  new_age <- do.call(rbind, lapply(metrics, function(metric) {
    fit_alpha_covariate(corrected, metric, "Age")
  }))
  rownames(new_age) <- NULL
  diagnostics <- save_diagnostics()
  if (any(diagnostics$Failed | !diagnostics$Converged) ||
      any(!is.finite(as.matrix(new_age[, p_columns]))) ||
      any(new_age$Fixed_design_rank_deficient) ||
      any(new_age$N_samples != 64L | new_age$N_patients != 31L)) {
    stop("A corrected-age fit needs review. See FIT_DIAGNOSTICS.csv in ", out_dir,
         ". No corrected supplementary table has been released.", call. = FALSE)
  }

  updated <- previous
  shared_columns <- intersect(names(new_age), names(updated))
  age_rows <- match(paste(metrics, "Age", sep = "|"), result_keys)
  updated[age_rows, shared_columns] <- new_age[match(metrics, new_age$Metric), shared_columns]
  unchanged_rows <- previous$Covariate != "Age"
  stopifnot(identical(previous$Covariate_P[unchanged_rows],
                      updated$Covariate_P[unchanged_rows]))
  for (metric in metrics) {
    rows <- which(updated$Metric == metric)
    stopifnot(length(rows) == 3L, all(is.finite(updated$Covariate_P[rows])))
    updated$Covariate_P_BH_within_metric[rows] <-
      p.adjust(updated$Covariate_P[rows], method = "BH", n = 3L)
  }
  stopifnot(identical(original[, other_columns], corrected[, other_columns]))

  format_p <- function(x) vapply(x, function(p) {
    if (!is.finite(p)) return("NA")
    if (p < 0.001) return("<0.001")
    sprintf("%.3f", p)
  }, character(1))
  publication <- data.frame(
    Metric = updated$Metric, Covariate = updated$Covariate,
    "n samples" = updated$N_samples, "n patients" = updated$N_patients,
    "Group P, base" = format_p(updated$Group_P_base_same_rows),
    "Group P, adjusted" = format_p(updated$Group_P_adjusted),
    "Reads P, adjusted" = format_p(updated$Reads_P_adjusted),
    "Covariate P" = format_p(updated$Covariate_P),
    "Covariate BH FDR" = format_p(updated$Covariate_P_BH_within_metric),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  comparison <- data.frame(
    Metric = updated$Metric, Covariate = updated$Covariate,
    Previous_covariate_P = previous$Covariate_P,
    Corrected_covariate_P = updated$Covariate_P,
    Previous_covariate_BH_FDR = previous$Covariate_P_BH_within_metric,
    Corrected_covariate_BH_FDR = updated$Covariate_P_BH_within_metric,
    Previous_group_P_adjusted = previous$Group_P_adjusted,
    Corrected_group_P_adjusted = updated$Group_P_adjusted,
    Previous_reads_P_adjusted = previous$Reads_P_adjusted,
    Corrected_reads_P_adjusted = updated$Reads_P_adjusted,
    Previous_BH_below_005 = previous$Covariate_P_BH_within_metric < 0.05,
    Corrected_BH_below_005 = updated$Covariate_P_BH_within_metric < 0.05,
    stringsAsFactors = FALSE
  )
  write.csv(updated, file.path(out_dir, "FINAL_alpha_covariate_sensitivity_AGE_CORRECTED.csv"),
            row.names = FALSE)
  write.csv(publication, file.path(out_dir, "ST4_2_3_ALPHA_COVARIATES_AGE_CORRECTED.csv"),
            row.names = FALSE)
  write.csv(comparison, file.path(out_dir, "BEFORE_AFTER_ALL_COVARIATES.csv"), row.names = FALSE)
  write.csv(comparison[comparison$Covariate == "Age", , drop = FALSE],
            file.path(out_dir, "BEFORE_AFTER_AGE_ONLY.csv"), row.names = FALSE)

  readme <- c(
    "CHAPTER 4 SUPPLEMENTARY ALPHA DIVERSITY WITH CORRECTED AGES",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "Input: the embedded original primary unfiltered alpha table, 64 biopsies from 31 patients.",
    "Original input filename: Alpha_input_PRIMARY_no_prevalence_filter_COMMON_READS.csv",
    "Original input SHA256: 0731f6260c528a145da462a5eba94ba71288fe7ad2f9b2d3793f513df6a482cf",
    "Original sensitivity-result SHA256: 306609252186400001b046c2942c61a8787776211473c8d0d0e48f221d4d7a4d",
    "Age source: Table_2_Updated_Ages (1).docx, matching Methods FINAL(7).",
    "Corrected ages are supplied in alpha_corrected_age_key.csv and reflect age at biopsy.",
    "43 ages were replaced; all other input fields and the 64-sample order were retained.",
    "Original sample and patient identifiers were retained.",
    "This script updates the age values; it does not independently recalculate DOB-based ages.",
    "Model for each metric: Metric ~ Group + Reads + Age + (1 | Patient).",
    "ML fits, bobyqa optimizer, maxfun 200000; likelihood-ratio tests as in the original script.",
    "Reads is the unchanged raw post-rarefaction pre-prevalence-filter count.",
    "No diversity metrics are recalculated and no samples are added or removed.",
    "Six Gender/Location raw-P rows were retained from the original sensitivity output.",
    "BH was recalculated across Age, Gender and Location separately within each metric (three tests).",
    "Accordingly, Gender/Location adjusted P values may change although their raw P values do not.",
    "Main alpha-group comparisons, alpha plots, beta analyses and Chapter 6 models were not replaced.",
    "Use ST4_2_3_ALPHA_COVARIATES_AGE_CORRECTED.csv for the supplementary sensitivity table.",
    "Use BEFORE_AFTER_AGE_ONLY.csv to check whether the Chapter 4 age-sensitivity sentence changes.",
    "FIT_DIAGNOSTICS.csv records warnings, convergence and singularity for every fitted model.",
    "A singular fit is reported for review; it is not silently substituted with a different model.",
    "The age-corrected models were run only when this output folder was generated in R."
  )
  writeLines(readme, file.path(out_dir, "README_RESULTS.txt"))
  writeLines("COMPLETED", file.path(out_dir, "RUN_STATUS.txt"))
  writeLines("Completed corrected-age fits and within-metric BH adjustment.",
             file.path(out_dir, "COMPLETED.txt"))
  print(publication, row.names = FALSE)
  message("Finished. Results saved in: ", out_dir)
  invisible(out_dir)
}

run_alpha_age_correction(output_parent)
