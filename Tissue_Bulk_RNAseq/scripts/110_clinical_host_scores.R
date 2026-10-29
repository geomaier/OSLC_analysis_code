#!/usr/bin/env Rscript
# Exploratory clinical associations using the corrected age at biopsy.
# Score scaling uses 63 biopsies; clinical models use 37 lesions from 22 patients.
# Model: Score_z ~ Age10 + Gender + Region + BRAF + (1 | Patient).
# ML fits, Satterthwaite tests and BH correction across 16 primary tests.
# Command: Rscript scripts/110_clinical_host_scores.R corrected_master_csv output_parent
# Alternatively set OSLC_CLINICAL_INPUT and OSLC_CLINICAL_OUTPUT.
# The private age-validation CSV is required; see the folder README.
# Scientific source run: 11 October 2026. Packaging validation is recorded separately.

options(stringsAsFactors = FALSE, contrasts = c("contr.treatment", "contr.poly"))
set.seed(1)
THIS_SCRIPT <- local({
  all_args <- commandArgs(trailingOnly = FALSE)
  command_file <- sub("^--file=", "", all_args[grepl("^--file=", all_args)])
  source_files <- unlist(lapply(sys.frames(), function(f) {
    if (!is.null(f$ofile)) as.character(f$ofile) else NULL
  }))
  candidates <- unique(c(command_file, rev(source_files)))
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates)) normalizePath(candidates[1], winslash = "/", mustWork = TRUE) else NA_character_
})
args <- commandArgs(trailingOnly = TRUE)
input_csv <- if (length(args) >= 1L) args[1] else Sys.getenv("OSLC_CLINICAL_INPUT")
if (!nzchar(input_csv) || !file.exists(input_csv)) {
  stop("Provide the corrected clinical master CSV as argument 1 or OSLC_CLINICAL_INPUT.")
}
input_csv <- normalizePath(input_csv, winslash = "/", mustWork = TRUE)
bundle_dir <- dirname(input_csv)
output_parent <- if (length(args) >= 2L) args[2] else
  Sys.getenv("OSLC_CLINICAL_OUTPUT", file.path(getwd(), "outputs", "Clinical_Age_Corrected_Results"))
validation_file <- Sys.getenv("OSLC_CLINICAL_AGE_KEY",
  file.path(bundle_dir, "clinical_age_validation.csv"))
if (!file.exists(validation_file)) stop("Clinical age-validation CSV not found: ", validation_file)
out_dir <- file.path(output_parent, paste0("clinical_age_corrected_", format(Sys.time(), "%Y%m%d_%H%M%S")))
if (dir.exists(out_dir)) stop("Output folder already exists; run again with a different timestamp.")
required <- c("lme4", "lmerTest")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install packages first: install.packages(c(", paste(sprintf("'%s'", missing), collapse = ", "), "))")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) stop("Could not create output folder.")

# Retain complete copies of the input and rerun code with the new results.
snapshot_dir <- file.path(out_dir, "inputs_and_code")
dir.create(snapshot_dir, showWarnings = FALSE)
file.copy(input_csv, file.path(snapshot_dir, "MASTER_METADATA_CH4_AGE_CORRECTED.csv"), overwrite = FALSE)
if (!file.copy(validation_file, file.path(snapshot_dir, "clinical_age_validation.csv"), overwrite = FALSE)) {
  stop("Could not retain the private age-validation input with the run records.")
}
for (name in c("Metadata_March26_AGE_CORRECTED.csv", "Age_Correction_Audit.csv")) {
  path <- file.path(bundle_dir, name)
  if (file.exists(path)) file.copy(path, file.path(snapshot_dir, name), overwrite = FALSE)
}
if (!is.na(THIS_SCRIPT)) file.copy(THIS_SCRIPT, file.path(snapshot_dir, basename(THIS_SCRIPT)), overwrite = FALSE)

# Keep the supplied Batch column separate from the main biological model.
RUN_BATCH_CHECK <- TRUE
scores <- c("Gastric_Metaplasia_Score", "Interferon_Alpha_Score", "Inflammatory_Response_Score", "TNFA_NFKB_Score")
score_labels <- setNames(c("Gastric metaplasia", "Interferon-alpha", "Inflammatory response", "TNF-alpha / NF-kB"), scores)
predictors <- c("Age10", "Gender", "Region", "BRAF")
coef_names <- c(Age10 = "Age10", Gender = "GenderM", Region = "RegionProximal", BRAF = "BRAFV600E")
units <- c(Age10 = "SD score per 10 years older", Gender = "SD score: M minus F",
           Region = "SD score: proximal minus distal", BRAF = "SD score: V600E minus WT")

# 1. Metadata ---------------------------------------------------------------
d <- read.csv(input_csv, check.names = FALSE, na.strings = c("", "NA", "N/A"))
required_cols <- c("SampleID_phyloseq", "Patient", "Group", "Sample_Type_Original", "Age", "Gender", "Location", "BRAF", scores)
if (RUN_BATCH_CHECK) required_cols <- c(required_cols, "Batch")
if (length(setdiff(required_cols, names(d)))) stop("Missing columns: ", paste(setdiff(required_cols, names(d)), collapse = ", "))
for (v in c("SampleID_phyloseq", "Patient", "Group", "Sample_Type_Original", "Gender", "Location", "BRAF")) d[[v]] <- trimws(as.character(d[[v]]))
if (anyNA(d$SampleID_phyloseq) || any(!nzchar(d$SampleID_phyloseq)) || anyDuplicated(d$SampleID_phyloseq)) stop("Missing or duplicate sample IDs.")
if (anyNA(d$Patient) || any(!nzchar(d$Patient))) stop("Missing patient IDs.")
if (any(d$Patient == "7699")) stop("Patient 7699 must be excluded.")
for (v in c("Age", scores)) {
  old <- d[[v]]; d[[v]] <- suppressWarnings(as.numeric(old))
  if (any(!is.na(old) & is.na(d[[v]])) || any(is.infinite(d[[v]]))) stop("Invalid numeric values in ", v)
}
# Exact confirmed ages for this age-only rerun; retained original sample labels.
age_validation <- read.csv(validation_file, colClasses = "character",
                           check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("SampleID", "Age") %in% names(age_validation)) ||
    nrow(age_validation) != 63L || anyNA(age_validation$SampleID) ||
    any(!nzchar(age_validation$SampleID)) || anyDuplicated(age_validation$SampleID)) {
  stop("Clinical age-validation input must contain the 63 unique SampleID and Age records.")
}
age_validation$Age <- suppressWarnings(as.numeric(age_validation$Age))
if (any(!is.finite(age_validation$Age))) stop("Invalid ages in the clinical validation input.")
expected_ages <- setNames(age_validation$Age, age_validation$SampleID)
if (!setequal(d$SampleID_phyloseq, names(expected_ages)) || nrow(d) != 63L)
  stop("The input differs from the original 63-row clinical cohort.")
if (anyNA(d$Age) || any(d$Age != unname(expected_ages[d$SampleID_phyloseq])))
  stop("Age does not match the confirmed age-at-biopsy table. Use the supplied corrected master CSV.")
if (any(d$Age < 0 | d$Age > 120, na.rm = TRUE)) stop("Implausible age.")
if (any(!is.na(d$Gender) & !d$Gender %in% c("F", "M"))) stop("Gender must be coded F/M.")
if (any(!is.na(d$BRAF) & !d$BRAF %in% c("WT", "V600E"))) stop("BRAF must be coded WT/V600E.")
for (v in c("Age", "Gender")) {
  n <- tapply(d[[v]], d$Patient, function(x) length(unique(x[!is.na(x)])))
  if (any(n > 1L)) stop(v, " is inconsistent within a patient.")
}
region_map <- c("Caecum" = "Proximal", "Ascending Colon" = "Proximal", "Transverse Colon" = "Proximal",
                "Descending Colon" = "Distal", "Sigmoid Colon" = "Distal", "Rectum" = "Distal")
if (any(!is.na(d$Location) & !d$Location %in% names(region_map))) stop("Unknown Location; update region_map explicitly.")
d$Region <- factor(unname(region_map[d$Location]), levels = c("Distal", "Proximal"))
d$Gender <- factor(d$Gender, levels = c("F", "M"))
d$BRAF <- factor(d$BRAF, levels = c("WT", "V600E"))
d$Patient <- factor(d$Patient)
if (RUN_BATCH_CHECK) d$Batch <- factor(trimws(as.character(d$Batch)))
d$Age10 <- (d$Age - mean(d$Age, na.rm = TRUE)) / 10
# Preserve the earlier scaling: scores standardised over the full 63-row
# metadata cohort, then the same scaling used in every clinical subset.
scaling <- data.frame(Score = scores, Mean = NA_real_, SD = NA_real_)
for (i in seq_along(scores)) {
  s <- scores[i]; mu <- mean(d[[s]], na.rm = TRUE); sig <- sd(d[[s]], na.rm = TRUE)
  if (!is.finite(sig) || sig <= 0) stop("Constant or missing score: ", s)
  d[[paste0(s, "_z")]] <- (d[[s]] - mu) / sig
  scaling$Mean[i] <- mu; scaling$SD[i] <- sig
}
lesions <- d[!is.na(d$Group) & d$Group == "SSL_HP", ]
if (nrow(lesions) != 37L || length(unique(lesions$Patient)) != 22L)
  stop("Clinical cohort changed: expected 37 SSL/HP lesions from 22 patients.")
if (anyNA(lesions$Sample_Type_Original) || any(!lesions$Sample_Type_Original %in% c("SSL", "HP"))) stop("Unexpected lesion histology.")
write.csv(scaling, file.path(out_dir, "00_score_scaling.csv"), row.names = FALSE)
keep_cols <- unique(c(required_cols, "Region", "Age10", paste0(scores, "_z")))
write.csv(lesions[, keep_cols], file.path(out_dir, "00_clinical_input.csv"), row.names = FALSE)
write.csv(data.frame(Variable = keep_cols, Missing = vapply(lesions[, keep_cols], function(x) sum(is.na(x)), integer(1))),
          file.path(out_dir, "00_missingness.csv"), row.names = FALSE)
counts <- list()
for (v in c("Gender", "Region", "BRAF", if (RUN_BATCH_CHECK) "Batch")) {
  for (level in unique(as.character(lesions[[v]][!is.na(lesions[[v]])]))) {
    x <- lesions[!is.na(lesions[[v]]) & as.character(lesions[[v]]) == level, ]
    counts[[length(counts) + 1L]] <- data.frame(Variable = v, Level = level, N_samples = nrow(x), N_patients = length(unique(x$Patient)))
  }
}
write.csv(do.call(rbind, counts), file.path(out_dir, "00_predictor_counts.csv"), row.names = FALSE)
if (RUN_BATCH_CHECK) for (v in c("Region", "Gender", "BRAF"))
  write.csv(as.data.frame(table(lesions$Batch, lesions[[v]])), file.path(out_dir, paste0("00_Batch_by_", v, ".csv")), row.names = FALSE)

# 2. Fit complete-case ML models and retain warnings -------------------------
fits <- list(); diagnostics <- list(); fit_log <- list()
fit_model <- function(data, score, terms, id) {
  needed <- c("SampleID_phyloseq", "Patient", paste0(score, "_z"), terms)
  x <- data[complete.cases(data[, needed, drop = FALSE]), needed, drop = FALSE]
  for (v in names(x)) if (is.factor(x[[v]])) x[[v]] <- droplevels(x[[v]])
  x$Y <- x[[paste0(score, "_z")]]
  formula_text <- paste("Y ~", paste(terms, collapse = " + "), "+ (1 | Patient)")
  info <- data.frame(Model = id, Score = score, N_samples = nrow(x), N_patients = length(unique(x$Patient)),
    Formula = formula_text, Status = "not_estimable", Singular = NA, Convergence_OK = NA, ICC = NA_real_, Warnings = "", Note = "")
  record <- function(fit = NULL) {
    fit_log[[length(fit_log) + 1L]] <<- info
    list(fit = fit, info = info, data = x)
  }
  if (nrow(x) < 10L || length(unique(x$Patient)) < 5L || length(unique(x$Patient)) >= nrow(x)) {
    info$Note <- "Insufficient complete cases, patients or repeated biopsies."; return(record())
  }
  if (!is.finite(sd(x$Y)) || sd(x$Y) <= 0 || any(vapply(x[, terms, drop = FALSE], function(z) length(unique(z)) < 2L, logical(1)))) {
    info$Note <- "Invariant outcome or model term."; return(record())
  }
  mm <- model.matrix(reformulate(terms), data = x)
  if (qr(mm)$rank < ncol(mm) || nrow(x) - ncol(mm) < 8L) {
    info$Note <- "Rank-deficient or insufficient residual dimension; model not fitted."; return(record())
  }
  warns <- character(); msgs <- character()
  fit <- tryCatch(withCallingHandlers(lmerTest::lmer(as.formula(formula_text), data = x, REML = FALSE,
    control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }), error = function(e) e)
  if (inherits(fit, "error")) { info$Note <- conditionMessage(fit); return(record()) }
  info$Singular <- lme4::isSingular(fit, tol = 1e-4)
  conv <- fit@optinfo$conv$lme4$messages
  conv <- conv[!grepl("boundary.*singular", conv, ignore.case = TRUE)]
  opt <- fit@optinfo$conv$opt
  info$Convergence_OK <- length(conv) == 0L && (is.null(opt) || all(opt == 0))
  info$Status <- if (info$Convergence_OK) "ok" else "review_fit"
  info$Warnings <- paste(unique(c(warns, msgs, conv)), collapse = " | ")
  vc <- as.data.frame(lme4::VarCorr(fit)); vp <- vc$vcov[vc$grp == "Patient"]; ve <- vc$vcov[vc$grp == "Residual"]
  if (length(vp) == 1L && length(ve) == 1L) info$ICC <- vp / (vp + ve)
  if (info$Singular) info$Note <- "Patient variance estimated near zero; inspect diagnostics."
  if ("Batch" %in% terms) {
    small <- any(tapply(as.character(x$Patient), x$Batch, function(z) length(unique(z))) < 3L)
    if (small) info$Note <- paste(info$Note, "Some batch levels have fewer than three patients.")
  }
  fits[[id]] <<- fit
  diagnostics[[id]] <<- list(fit = fit, info = info)
  record(fit)
}

extract_effects <- function(object, terms, block) {
  info <- object$info; tab <- if (!is.null(object$fit)) coef(summary(object$fit)) else NULL
  do.call(rbind, lapply(terms, function(v) {
    r <- data.frame(Block = block, Score = info$Score, Predictor = v, Effect_unit = units[v],
      N_samples = info$N_samples, N_patients = info$N_patients, Estimate = NA_real_, SE = NA_real_, DF = NA_real_,
      CI_low = NA_real_, CI_high = NA_real_, P_raw = NA_real_, P_for_correction = NA_real_,
      Singular = info$Singular, Convergence_OK = info$Convergence_OK, Status = info$Status, Note = info$Note)
    cn <- unname(coef_names[v])
    if (!is.null(tab) && cn %in% rownames(tab)) {
      r$Estimate <- tab[cn, "Estimate"]; r$SE <- tab[cn, "Std. Error"]; r$DF <- tab[cn, "df"]
      r$P_raw <- tab[cn, "Pr(>|t|)"]
      if (!is.na(r$DF) && r$DF > 0) {
        delta <- qt(.975, r$DF) * r$SE; r$CI_low <- r$Estimate - delta; r$CI_high <- r$Estimate + delta
      }
      if (r$Status == "ok" && is.finite(r$P_raw)) r$P_for_correction <- r$P_raw
    }
    r
  }))
}
adjust_planned <- function(p, method = "BH") {
  out <- rep(NA_real_, length(p)); ok <- is.finite(p)
  if (any(ok)) out[ok] <- p.adjust(p[ok], method = method, n = length(p))
  out
}

# 3. Four main models and labelled sensitivities -----------------------------
main <- list(); batch_adjusted <- list(); batch_tests <- list()
for (score in scores) {
  cat("Fitting ", score, "\n", sep = "")
  fit <- fit_model(lesions, score, predictors, paste0("MAIN_", score))
  main[[score]] <- extract_effects(fit, predictors, "Main_SSL_HP")
  if (RUN_BATCH_CHECK) {
    # Full and reduced batch models MUST use exactly the same samples.
    vars <- c("Patient", paste0(score, "_z"), predictors, "Batch")
    x <- lesions[complete.cases(lesions[, vars]), ]
    full <- fit_model(x, score, c(predictors, "Batch"), paste0("BATCH_FULL_", score))
    reduced <- fit_model(x, score, predictors, paste0("BATCH_REDUCED_", score))
    batch_adjusted[[score]] <- extract_effects(full, predictors, "Sensitivity_batch_adjusted")
    r <- data.frame(Score = score, N_samples = nrow(x), N_patients = length(unique(x$Patient)),
      Batch_levels = length(unique(x$Batch)), Chi_square = NA_real_, DF = NA_real_, P_raw = NA_real_,
      P_for_correction = NA_real_, Status = "not_estimable", Singular_full = full$info$Singular,
      Note = "Supplied Batch definition is unverified. This is a technical sensitivity, not a biological effect.")
    if (!is.null(full$fit) && !is.null(reduced$fit)) {
      if (!identical(rownames(model.frame(full$fit)), rownames(model.frame(reduced$fit)))) stop("Batch model row alignment failed.")
      lf <- logLik(full$fit); lr <- logLik(reduced$fit)
      r$DF <- as.numeric(attr(lf, "df") - attr(lr, "df"))
      statistic <- 2 * (as.numeric(lf) - as.numeric(lr))
      if (r$DF > 0 && statistic >= -1e-7 && full$info$Status == "ok" && reduced$info$Status == "ok") {
        r$Chi_square <- max(0, statistic)
        r$P_raw <- pchisq(r$Chi_square, df = r$DF, lower.tail = FALSE)
        r$P_for_correction <- r$P_raw; r$Status <- "ok"
      } else r$Note <- paste(r$Note, "Model comparison or convergence requires review.")
    }
    batch_tests[[score]] <- r
  }
}
primary <- do.call(rbind, main); rownames(primary) <- NULL
stopifnot(nrow(primary) == 16L)
primary$P_BH_16 <- adjust_planned(primary$P_for_correction)
primary$Below_005 <- !is.na(primary$P_BH_16) & primary$P_BH_16 < .05
primary <- primary[order(primary$P_BH_16, primary$P_raw, na.last = TRUE), ]
write.csv(primary, file.path(out_dir, "01_PRIMARY_CLINICAL_16_TESTS.csv"), row.names = FALSE)
write.csv(primary[primary$Below_005, ], file.path(out_dir, "02_PRIMARY_BH_BELOW_005.csv"), row.names = FALSE)
if (RUN_BATCH_CHECK) {
  b <- do.call(rbind, batch_adjusted); b$P_BH_16 <- adjust_planned(b$P_for_correction)
  write.csv(b, file.path(out_dir, "03_BATCH_ADJUSTED_CLINICAL.csv"), row.names = FALSE)
  bt <- do.call(rbind, batch_tests); bt$P_BH_4 <- adjust_planned(bt$P_for_correction)
  write.csv(bt, file.path(out_dir, "04_BATCH_OMNIBUS.csv"), row.names = FALSE)
}
write.csv(do.call(rbind, fit_log), file.path(out_dir, "05_FIT_LOG.csv"), row.names = FALSE)
saveRDS(fits, file.path(out_dir, "models.rds"))

# 4. Diagnostics and main effect plots ---------------------------------------
pdf(file.path(out_dir, "model_diagnostics.pdf"), width = 10, height = 5.5)
for (id in names(diagnostics)) {
  f <- diagnostics[[id]]$fit
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
  plot(fitted(f), resid(f), pch = 16, col = "#42566B", xlab = "Fitted score (SD)", ylab = "Residual")
  abline(h = 0, lty = 2, col = "grey50"); lines(lowess(fitted(f), resid(f)), col = "#B24736", lwd = 2)
  qqnorm(resid(f), pch = 16, col = "#42566B", main = "Residual normal Q-Q plot"); qqline(resid(f), col = "#B24736", lwd = 2)
  mtext(id, outer = TRUE, cex = .8)
}
dev.off()
labels <- c(Age10 = "Age (+10 years)", Gender = "M vs F", Region = "Proximal vs distal", BRAF = "BRAF V600E vs WT")
pdf(file.path(out_dir, "clinical_forest_plots.pdf"), width = 9, height = 5.5)
for (score in scores) {
  q <- primary[primary$Score == score, ]; q <- q[match(predictors, q$Predictor), ]
  good <- is.finite(q$CI_low) & is.finite(q$CI_high); yy <- rev(seq_along(predictors))
  limits <- range(c(-.1, .1, q$CI_low[good], q$CI_high[good])); pad <- diff(limits) * .1
  par(mar = c(5, 12, 4, 2))
  plot(NA, xlim = limits + c(-pad, pad), ylim = c(.5, length(predictors) + .5), yaxt = "n", ylab = "",
       xlab = "Difference in signature score (SD); ordinary 95% CI", main = score_labels[score])
  axis(2, at = yy, labels = labels[predictors], las = 1, cex.axis = .85); abline(v = 0, lty = 2, col = "grey60")
  col <- ifelse(q$Below_005, "#B24736", "#42566B")
  segments(q$CI_low[good], yy[good], q$CI_high[good], yy[good], col = col[good], lwd = 2)
  points(q$Estimate[good], yy[good], pch = ifelse(q$Singular[good], 1, 16), col = col[good], cex = 1.2)
  if (any(!good)) text(rep(limits[1], sum(!good)), yy[!good], "Not estimable; see table", adj = 0, cex = .7)
  mtext("Red: BH across all 16 main clinical tests < 0.05. Open circle: singular fit.", side = 1, line = 3.5, cex = .7)
}
dev.off()

# 5. Methods and interpretation ---------------------------------------------
writeLines(c(
  "FOUR-SCORE CLINICAL ASSOCIATIONS ONLY",
  paste("Input:", input_csv), paste("Lesions:", nrow(lesions), "Patients:", length(unique(lesions$Patient))),
  "SSL and HP are pooled as SSL_HP throughout all analyses.",
  "Main model: Score_z ~ Age10 + Gender + Region + BRAF + (1 | Patient).",
  "Each score is fitted separately. Clinical predictors are mutually adjusted; associations are conditional, not causal effects.",
  "Scores are used as supplied, with SD scaling fixed across the full metadata cohort for comparability with the preceding analysis.",
  "Age remains continuous (per 10 years). Region: caecum/ascending/transverse versus descending/sigmoid/rectum.",
  "BRAF is a sample-level variable tested within pooled SSL_HP lesions. Normal tissue is not labelled as an additional WT control group.",
  "ML estimation; Satterthwaite coefficient P values and ordinary (not multiplicity-adjusted) 95% confidence intervals.",
  "Main multiplicity: BH across 16 tests (four clinical predictors x four scores).",
  "Batch-adjusted clinical effects: separate BH across 16 tests, using pooled SSL_HP lesions.",
  "Batch overall association: four likelihood-ratio tests, BH across four. All levels are tested together; there is no single beta for a multi-level batch factor.",
  "For batch LRTs, full and reduced models use exactly the same complete-case samples and are both ML fits.",
  "The supplied Batch label has not been established as RNA-sequencing batch. Confirm its meaning before including it as an expression-related technical covariate in the paper.",
  "Batch has sparse levels in this cohort. Inspect counts, confounding, residuals, convergence and uncertainty.",
  "No microbiome abundance, presence, alpha diversity or 16S read-depth terms are used.",
  "Lesion size and dysplasia are not in this metadata table and are not invented or automatically included.",
  "The score-specific complete-case sample/patient counts and all skipped-model reasons are exported.",
  "Nonconverged/failed tests remain in the planned correction denominator and are not counted as significant.",
  "The analysis scope was narrowed after inspection of the broader 40-test screen. Treat this as an exploratory clinical analysis; do not present it as prospectively prespecified or conceal the preceding screen.",
  "Report beta, ordinary 95% CI, exact correction scope, sample/patient counts and adjustment variables. Report null results as well as associations.",
  "Age-only rerun: confirmed completed years at the biopsy date; 42 of 63 metadata ages and 24 of 37 lesion ages changed.",
  "Sample labels, patient cluster labels, anatomy, cohort, scores, score scaling, model formula and multiplicity scopes are retained from the previous clinical run.",
  "The original sample IDs containing 8172 remain unchanged; their age uses the confirmed DOB for patient 3735.",
  "Previous results dated 26 September 2026 are included only for numerical comparison; inspect current fit diagnostics before interpreting changes.",
  "The code and data were prepared without executing R in the authoring environment; this folder is produced by the user's R execution.",
  "Methods: https://lme4.github.io/lme4/reference/lmer.html ; https://cran.r-project.org/package=lmerTest ; https://www.strobe-statement.org/"
), file.path(out_dir, "README_METHODS.txt"))
write.csv(data.frame(Input = input_csv, MD5 = unname(tools::md5sum(input_csv))), file.path(out_dir, "input_checksum.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
cat("\nDONE. Results saved to:\n", out_dir, "\n", sep = "")
print(primary[, c("Score", "Predictor", "Estimate", "CI_low", "CI_high", "P_raw", "P_BH_16", "Singular", "Status")], row.names = FALSE)


# 6. Compare new results with the actual previous completed clinical run --------
previous_dir <- file.path(bundle_dir, "previous_results")
compare_results <- function(filename, keys, p_column, output_name) {
  old_path <- file.path(previous_dir, filename)
  new_path <- file.path(out_dir, filename)
  if (!file.exists(old_path) || !file.exists(new_path)) {
    warning("Comparison skipped: missing ", filename)
    return(invisible(NULL))
  }
  old <- read.csv(old_path, check.names = FALSE)
  new <- read.csv(new_path, check.names = FALSE)
  values <- intersect(c("N_samples", "N_patients", "Estimate", "SE", "DF", "CI_low", "CI_high",
                        "Chi_square", "P_raw", p_column, "Singular", "Singular_full", "Convergence_OK", "Status"),
                      intersect(names(old), names(new)))
  ans <- merge(old[, c(keys, values), drop = FALSE], new[, c(keys, values), drop = FALSE],
               by = keys, all = TRUE, suffixes = c("_Previous", "_CorrectedAge"), sort = FALSE)
  if (nrow(ans) != nrow(old) || nrow(ans) != nrow(new)) stop("Result comparison has unmatched or duplicated tests.")
  for (v in intersect(c("Estimate", "P_raw", p_column), values))
    ans[[paste0(v, "_Change")]] <- ans[[paste0(v, "_CorrectedAge")]] - ans[[paste0(v, "_Previous")]]
  p_old <- ans[[paste0(p_column, "_Previous")]]
  p_new <- ans[[paste0(p_column, "_CorrectedAge")]]
  ans$BH_threshold_comparison <- ifelse(is.na(p_old) | is.na(p_new), "Review missing test",
    ifelse(p_old < .05 & p_new < .05, "Below 0.05 in both",
    ifelse(p_old >= .05 & p_new < .05, "Now below 0.05",
    ifelse(p_old < .05 & p_new >= .05, "No longer below 0.05", "At least 0.05 in both"))))
  write.csv(ans, file.path(out_dir, output_name), row.names = FALSE)
  invisible(ans)
}
compare_results("01_PRIMARY_CLINICAL_16_TESTS.csv", c("Score", "Predictor"), "P_BH_16", "06_PRIMARY_OLD_VS_CORRECTED_AGE.csv")
if (RUN_BATCH_CHECK) {
  compare_results("03_BATCH_ADJUSTED_CLINICAL.csv", c("Score", "Predictor"), "P_BH_16", "07_BATCH_ADJUSTED_OLD_VS_CORRECTED_AGE.csv")
  compare_results("04_BATCH_OMNIBUS.csv", "Score", "P_BH_4", "08_BATCH_OMNIBUS_OLD_VS_CORRECTED_AGE.csv")
}
old_scaling <- file.path(previous_dir, "00_score_scaling.csv")
if (file.exists(old_scaling)) {
  old <- read.csv(old_scaling)
  k <- match(scaling$Score, old$Score)
  if (anyNA(k) || any(abs(scaling$Mean - old$Mean[k]) > 1e-12) || any(abs(scaling$SD - old$SD[k]) > 1e-12))
    stop("Score scaling changed from the earlier clinical analysis; do not interpret this as an age-only comparison.")
}
cat("\nAge-corrected clinical rerun complete. Review 05_FIT_LOG.csv and 06_PRIMARY_OLD_VS_CORRECTED_AGE.csv.\n")
