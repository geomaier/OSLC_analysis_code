#!/usr/bin/env Rscript
# Beta diversity: patient-level control comparisons and exact-site paired inference.
# Run from 16S; set BETA_WORK_DIR, BETA_MATCHED_RUN and optional BETA_OUTPUT_DIR.
# Primary Aitchison comparisons: SSL_HP/H_COL (31 profiles), adjacent/H_COL
# (18 profiles), and lesion/adjacent (24 biopsies, 12 pairs, 10 patients).
# Bray-Curtis analyses use the same three comparisons as sensitivity analyses.
# Holm correction is within distance metric; dispersion P values are diagnostic.
# CLR and relative abundance are calculated per biopsy before patient averaging.
# Exact-site inference enumerates all 1,023 non-identity patient-level swaps.
# Reference statistics are validation targets only, never substituted for fits.
# Outputs are written to a new timestamped folder.

local({
PLOT_REVISION <- "UNIFIED_FIGURES_v7"
cat("\nRunning figure revision: ", PLOT_REVISION, "\n", sep = "")
PATH_WORK <- Sys.getenv("BETA_WORK_DIR")
if (!nzchar(PATH_WORK)) stop("Set BETA_WORK_DIR to the directory containing the final phyloseq input.")
OUTPUT_ROOT <- Sys.getenv("BETA_OUTPUT_DIR", unset = PATH_WORK)
MATCHED_RUN_DIR <- Sys.getenv(
  "BETA_MATCHED_RUN",
  unset = file.path(PATH_WORK, "Chapter4_MATCHED_SSLHP_HADJSSLHP_without7699",
                    "run_20261006_203539_746132")
)
# Exported-count mode supports validation against saved analysis inputs.
INPUT_MODE <- Sys.getenv("BETA_INPUT_MODE", unset = "phyloseq")
EXPORTED_INPUT_DIR <- Sys.getenv("BETA_INPUT_DIR", unset = file.path(PATH_WORK, "inputs"))
CHECK_PREVIOUS_RUN <- TRUE
stopifnot(INPUT_MODE %in% c("phyloseq", "exported_counts"))

required_packages <- c("phyloseq", "dplyr", "vegan", "permute", "tibble",
                       "ggplot2", "grid", "gridExtra")
missing_packages <- required_packages[!vapply(required_packages,
  function(p) requireNamespace(p, quietly = TRUE), logical(1))]
if (length(missing_packages)) stop("Required R packages are missing: ",
  paste(missing_packages, collapse = ", "), ". No packages are installed by this script.")
suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(vegan)
  library(tibble)
  library(ggplot2)
  library(grid)
})
options(stringsAsFactors = FALSE)
if (capabilities("cairo")) options(bitmapType = "cairo")
if (!dir.exists(PATH_WORK)) stop("Input folder not found: ", PATH_WORK)
PATH_WORK <- normalizePath(PATH_WORK, winslash = "/", mustWork = TRUE)
phyloseq_input <- file.path(PATH_WORK, "phyloseq_final_March26_5pct_without7699.rds")
if (INPUT_MODE == "phyloseq" && !file.exists(phyloseq_input))
  stop("Phyloseq input not found: ", phyloseq_input)
if (!dir.exists(MATCHED_RUN_DIR)) stop("Matched input folder not found: ", MATCHED_RUN_DIR)

matched_definitions <- data.frame(
  Analysis = c("Aitchison_ASV", "BrayCurtis_ASV", "BrayCurtis_Genus"),
  directory = c("02_AITCHISON", "03_BRAY_CURTIS", "03_BRAY_CURTIS/Genus"),
  stem = c("PCoA_Aitchison_EXACT_MATCHED", "PCoA_Bray_Curtis_EXACT_MATCHED",
           "PCoA_Genus_Bray_Curtis_EXACT_MATCHED"),
  label = c("Aitchison distance", "Bray-Curtis distance", "Genus Bray-Curtis distance"),
  stringsAsFactors = FALSE
)
matched_manifest_path <- file.path(MATCHED_RUN_DIR,
  "00_PAIRING_AND_QC/07_FINAL_MATCHED_SAMPLE_MANIFEST.csv")
matched_input_files <- c(matched_manifest_path,
  file.path(MATCHED_RUN_DIR, matched_definitions$directory,
            paste0(matched_definitions$stem, "_full_objects.rds")),
  file.path(MATCHED_RUN_DIR, matched_definitions$directory,
            paste0(matched_definitions$stem, "_statistics.csv")))
if (any(!file.exists(matched_input_files))) stop(
  "Missing original matched inputs:\n",
  paste(matched_input_files[!file.exists(matched_input_files)], collapse = "\n"),
  "\nKeep MATCHED_RUN_DIR pointing to the original complete matched run.")


storage_diagnostics <- function(directory) {
  lines <- c(paste0("Output folder: ", directory),
             paste0("Folder write access reported by R: ", file.access(directory, 2L) == 0L))
  if (.Platform$OS.type == "unix") {
    for (command in c("df", "quota")) {
      if (!nzchar(Sys.which(command))) next
      args <- if (command == "df") c("-h", shQuote(directory)) else "-s"
      answer <- tryCatch(suppressWarnings(system2(command, args, stdout = TRUE,
        stderr = TRUE, timeout = 5)), error = function(e) conditionMessage(e))
      lines <- c(lines, paste0(command, ":"), head(answer, 15L))
    }
  }
  lines
}

storage_stop <- function(destination, error) {
  stop(paste(c(
    "OUTPUT SAVE FAILED. This run is incomplete.",
    paste0("File: ", destination),
    paste0("R reported: ", conditionMessage(error)),
    storage_diagnostics(dirname(destination)),
    "Check free disk space and, on a university server, your storage quota.",
    "If necessary, set OUTPUT_ROOT near the top to another writable folder with space, then Source again.",
    "Keep PATH_WORK and MATCHED_RUN_DIR pointing to your existing input files."
  ), collapse = "\n"), call. = FALSE)
}

save_beta_rds <- function(object, file) {
  staging <- tempfile(pattern = ".beta-rds-", tmpdir = dirname(file), fileext = ".rds")
  on.exit(unlink(staging), add = TRUE)
  tryCatch({
    base::saveRDS(object, staging, compress = "gzip")
    # Prefer a rename within the destination filesystem; allow an existing file
    # to be refreshed within this new run (the overall and pairwise steps share it).
    installed <- suppressWarnings(file.rename(staging, file))
    if (!installed) {
      if (!file.copy(staging, file, overwrite = TRUE)) {
        stop("Could not copy the completed RDS into its output location.")
      }
    }
    invisible(base::readRDS(file))
    invisible(file)
  }, error = function(e) storage_stop(file, e))
}

check_output_storage <- function(directory) {
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE)) {
    storage_stop(directory, simpleError("Could not create the output folder."))
  }
  probe <- tempfile(pattern = ".beta-write-check-", tmpdir = directory, fileext = ".rds")
  on.exit(unlink(probe), add = TRUE)
  tryCatch({
    payload <- rep(as.raw(0:255), 4096L)
    base::saveRDS(payload, probe, compress = FALSE)
    if (!identical(payload, base::readRDS(probe))) stop("Output write/read check did not match.")
  }, error = function(e) storage_stop(probe, e))
  invisible(TRUE)
}

check_output_storage(OUTPUT_ROOT)
OUTPUT_ROOT <- normalizePath(OUTPUT_ROOT, winslash = "/", mustWork = TRUE)
cat("\nInput folder: ", PATH_WORK, "\nOutput folder: ", OUTPUT_ROOT, "\n", sep = "")

RUN_ID <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid(),
                 "_UNIFIED_FIGURES")


# Analysis and figure settings.
patient_col_manual <- NULL
n_permutations <- 9999L
clr_pseudocount <- 0.5
ellipse_level <- 0.68
export_dpi <- 600
font_family <- "Arial"
font_pt <- 7
pt_to_gg <- function(x) x / 2.845276
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
group_levels <- c("H_COL", "H_ADJ_SSL_HP", "SSL_HP")
group_display <- setNames(group_levels, group_levels)
pal <- c(H_COL = "#1313EC", H_ADJ_SSL_HP = "#E9C61D", SSL_HP = "#FF0000")
pal_fill <- c(H_COL = "#AEB8FF", H_ADJ_SSL_HP = "#F6E79A", SSL_HP = "#FFB3B3")
group_colours <- pal[c("H_ADJ_SSL_HP", "SSL_HP")]

CONTROL_LABELS <- c("SSL_HP vs H_COL", "H_ADJ_SSL_HP vs H_COL")
MATCHED_LABEL <- "SSL_HP vs H_ADJ_SSL_HP (exact-site matched)"
COMPARISONS <- c(CONTROL_LABELS, MATCHED_LABEL)
ANALYSES <- c("Aitchison_ASV", "BrayCurtis_ASV", "BrayCurtis_Genus")

source_args <- commandArgs(trailingOnly = FALSE)
source_files <- sub("^--file=", "", grep("^--file=", source_args, value = TRUE))
source_frames <- vapply(sys.frames(), function(e) {
  if (exists("ofile", envir = e, inherits = FALSE))
    as.character(get("ofile", envir = e))[1L] else ""
}, character(1))
file_flags <- which(source_args %in% c("-f", "--file"))
separate_files <- source_args[file_flags[file_flags < length(source_args)] + 1L]
script_candidates <- c(source_files, separate_files, rev(source_frames[nzchar(source_frames)]))
script_candidates <- script_candidates[file.exists(script_candidates)]
script_to_archive <- if (length(script_candidates))
  normalizePath(script_candidates[1L], winslash = "/", mustWork = TRUE) else NA_character_

run_dir <- file.path(OUTPUT_ROOT, "Beta_FINAL_three_comparisons", RUN_ID)
if (dir.exists(run_dir)) stop("Output folder already exists: ", run_dir)
if (!dir.create(run_dir, recursive = TRUE)) stop("Could not create: ", run_dir)
primary_dir <- file.path(run_dir, "01_PRIMARY_AITCHISON")
supplement_dir <- file.path(run_dir, "02_SUPPLEMENTARY_BRAY_CURTIS")
audit_dir <- file.path(run_dir, "03_AUDIT")
for (p in c(primary_dir, supplement_dir, audit_dir))
  if (!dir.create(p)) stop("Could not create: ", p)
out_dir <- primary_dir
warnings_seen <- character()
writeLines("Run started; RUN_COMPLETE.txt is written only after all checks and exports.",
           file.path(run_dir, "RUN_STATUS.txt"))

safe_write_csv <- function(x, path) {
  staging <- tempfile(pattern = ".beta-csv-", tmpdir = dirname(path), fileext = ".csv")
  on.exit(unlink(staging), add = TRUE)
  tryCatch({
    write.csv(as.data.frame(x, check.names = FALSE), staging, row.names = FALSE, na = "")
    if (file.info(staging)$size <= 0) stop("Empty CSV export.")
    if (!file.rename(staging, path)) stop("Could not install CSV: ", path)
  }, error = function(e) storage_stop(path, e))
}
write_matrix <- function(x, path) {
  safe_write_csv(data.frame(SampleID = rownames(x), x, check.names = FALSE), path)
}
counts_from_phy <- function(phy) {
  x <- as(phyloseq::otu_table(phy), "matrix")
  if (phyloseq::taxa_are_rows(phy)) x <- t(x)
  storage.mode(x) <- "double"
  if (any(!is.finite(x)) || any(x < 0) || any(rowSums(x) <= 0) ||
      any(abs(x - round(x)) > 1e-8)) stop("Expected non-negative integer counts and positive library sizes.")
  if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)))
    stop("Duplicate sample or feature identifiers in the count matrix.")
  x
}
analysis_method <- function(analysis) {
  if (analysis == "Aitchison_ASV") "Aitchison" else "Bray-Curtis"
}
comparison_annotation <- function(row) {
  paste0("R² = ", sprintf("%.4f", row$R2_PERMANOVA),
    "\n", format_final_p(row$P_PERMANOVA_raw, "PERMANOVA raw p"),
    "\n", format_final_p(row$P_PERMANOVA_Holm, "Holm p"),
    "\n", format_final_p(row$P_betadisper, "betadisper p"))
}
summary_row <- function(analysis, label, metadata, test_metadata, n_features,
                        fit, dispersion_test, permutations, unit, pairs = 0L) {
  st <- extract_adonis_stats(fit)
  group <- as.character(metadata$Group)
  dp <- as.numeric(dispersion_test$tab[1L, "Pr(>F)"])
  df <- as.numeric(dispersion_test$tab[1L, "F"])
  if (any(!is.finite(c(st$R2, st$F, st$P, dp, df)))) stop("Non-finite test results: ", label)
  data.frame(
    Analysis = analysis, Role = if (analysis == "Aitchison_ASV") "Primary" else "Sensitivity",
    Comparison = label, Distance = analysis_method(analysis),
    Feature_level = if (analysis == "BrayCurtis_Genus") "Genus" else "ASV",
    Profile_scale = if (analysis == "Aitchison_ASV") "CLR" else "Relative abundance",
    N_biopsies = nrow(metadata), N_patients = length(unique(metadata$Patient)),
    N_test_profiles = nrow(test_metadata), N_site_pairs = as.integer(pairs),
    N_features = n_features, H_COL_biopsies = sum(group == "H_COL"),
    SSL_HP_biopsies = sum(group == "SSL_HP"),
    H_ADJ_SSL_HP_biopsies = sum(group == "H_ADJ_SSL_HP"),
    Statistical_unit = unit,
    Permutations = nrow(permutations), Seed = 1L,
    R2_PERMANOVA = st$R2, F_PERMANOVA = st$F, P_PERMANOVA_raw = st$P,
    P_PERMANOVA_Holm = NA_real_, F_betadisper = df, P_betadisper = dp,
    R2_definition = "Tissue-group sum of squares / total sum of squares",
    stringsAsFactors = FALSE
  )
}
write_test_audit <- function(result, directory, counts = NULL) {
  safe_write_csv(result$metadata, file.path(directory, "test_metadata.csv"))
  safe_write_csv(result$biopsy_metadata, file.path(directory, "contributing_biopsies.csv"))
  if (!is.null(counts)) write_matrix(counts, file.path(directory, "biopsy_counts.csv"))
  if (!is.null(result$profiles))
    write_matrix(result$profiles, file.path(directory, "test_profiles.csv"))
  safe_write_csv(data.frame(Term = rownames(result$permanova),
    as.data.frame(result$permanova), check.names = FALSE),
    file.path(directory, "PERMANOVA_full_table.csv"))
  safe_write_csv(data.frame(Term = rownames(result$dispersion_test$tab),
    as.data.frame(result$dispersion_test$tab), check.names = FALSE),
    file.path(directory, "dispersion_full_table.csv"))
  safe_write_csv(as.data.frame(result$permutations),
                 file.path(directory, "permutation_indices.csv"))
  write_matrix(as.matrix(result$distance), file.path(directory, "distance_matrix.csv"))
  save_beta_rds(result, file.path(directory, "analysis_objects.rds"))
}
run_control <- function(counts, metadata, analysis, tissue) {
  label <- paste(tissue, "vs H_COL")
  keep <- as.character(metadata$Group) %in% c(tissue, "H_COL")
  md <- metadata[keep, , drop = FALSE]
  x <- counts[md$SampleID, , drop = FALSE]
  md$Group <- droplevels(md$Group)
  stopifnot(identical(rownames(x), as.character(md$SampleID)))
  p <- make_patient_profiles(x, md, analysis_method(analysis))
  pm <- p$metadata
  stopifnot(!anyDuplicated(pm$Patient))
  expected <- if (tissue == "SSL_HP") c(52L, 31L) else c(26L, 18L)
  stopifnot(nrow(md) == expected[1L], nrow(pm) == expected[2L])
  if (analysis != "Aitchison_ASV")
    stopifnot(max(abs(rowSums(p$profiles) - 1)) < 1e-10)
  set.seed(1)
  permutations <- permute::shuffleSet(nrow(pm), nset = n_permutations,
    control = permute::how(within = permute::Within(type = "free")))
  stopifnot(nrow(permutations) == 9999L)
  fit <- vegan::adonis2(p$distance ~ Group, data = pm,
                       permutations = permutations, by = NULL)
  dispersion <- vegan::betadisper(p$distance, pm$Group, type = "median", bias.adjust = TRUE)
  dispersion_test <- vegan::permutest(dispersion, permutations = permutations)
  result <- c(p, list(
    permanova = fit, dispersion = dispersion, dispersion_test = dispersion_test,
    permutations = permutations,
    summary = summary_row(analysis, label, md, pm, ncol(x), fit,
      dispersion_test, permutations, "Independent patient-tissue mean profiles")
  ))
  target <- file.path(audit_dir, analysis, gsub("[^A-Za-z0-9]+", "_", label))
  dir.create(target, recursive = TRUE)
  result$audit_directory <- target
  result$counts <- x
  result
}
run_exact_site <- function(definition, matched_metadata, main_metadata) {
  input <- file.path(MATCHED_RUN_DIR, definition$directory,
                    paste0(definition$stem, "_full_objects.rds"))
  summary_input <- file.path(MATCHED_RUN_DIR, definition$directory,
                            paste0(definition$stem, "_statistics.csv"))
  old <- readRDS(input)
  distance <- old$distance
  if (!inherits(distance, "dist")) stop("Missing saved matched distance: ", input)
  labels <- attr(distance, "Labels")
  md <- matched_metadata[match(labels, matched_metadata$SampleID), , drop = FALSE]
  stopifnot(!is.null(labels), !anyNA(md$SampleID), !anyDuplicated(labels),
            identical(as.character(md$SampleID), as.character(labels)))
  required <- c("SampleID", "Patient", "Group", "PairID", "Site")
  stopifnot(all(required %in% names(md)), !anyNA(md[, required]),
            nrow(md) == 24L, length(unique(md$PairID)) == 12L,
            length(unique(md$Patient)) == 10L)
  for (v in required) md[[v]] <- as.character(md[[v]])
  md$Group <- factor(md$Group, levels = c("H_ADJ_SSL_HP", "SSL_HP"))
  stopifnot(!anyNA(md$Group), all(table(md$Group) == 12L))
  main_index <- match(md$SampleID, main_metadata$SampleID)
  stopifnot(!anyNA(main_index),
    identical(md$Patient, as.character(main_metadata$Patient[main_index])),
    identical(as.character(md$Group), as.character(main_metadata$Group[main_index])))
  for (ii in split(seq_len(nrow(md)), md$PairID)) {
    stopifnot(length(ii) == 2L, length(unique(md$Site[ii])) == 1L,
              length(unique(md$Patient[ii])) == 1L, length(unique(md$Group[ii])) == 2L)
  }
  if (!is.null(old$metadata)) {
    om <- as.data.frame(old$metadata, stringsAsFactors = FALSE)
    stopifnot(all(required %in% names(om)))
    oi <- match(md$SampleID, om$SampleID)
    stopifnot(!anyNA(oi))
    for (v in required) stopifnot(identical(as.character(md[[v]]), as.character(om[[v]][oi])))
  }
  d <- as.matrix(distance)
  stopifnot(all(is.finite(d)), all(d >= -1e-12),
            max(abs(d - t(d))) < 1e-10, max(abs(diag(d))) < 1e-10)
  permutations <- patient_swap_matrix(md, n_permutations)
  stopifnot(nrow(permutations) == 1023L,
            nrow(unique(as.data.frame(permutations))) == 1023L,
            all(apply(permutations, 1L, function(i) !identical(i, seq_len(24L)))))
  set.seed(1)
  fit <- vegan::adonis2(distance ~ Group, data = md,
                       permutations = permutations, by = "terms")
  dispersion <- vegan::betadisper(distance, md$Group, type = "median", bias.adjust = FALSE)
  set.seed(1)
  dispersion_test <- vegan::permutest(dispersion, permutations = permutations)
  # Independent enumeration verifies the exact tissue-test P value.
  observed_f <- extract_adonis_stats(fit)$F
  perm_f <- apply(permutations, 1L, function(i) descriptive_partition(distance, md$Group[i])$F)
  enumerated_p <- (1 + sum(perm_f >= observed_f - sqrt(.Machine$double.eps))) / 1024
  stopifnot(abs(enumerated_p - extract_adonis_stats(fit)$P) < 1e-10)
  saved_summary <- read.csv(summary_input, check.names = FALSE)
  stopifnot(nrow(saved_summary) == 1L, "N_features" %in% names(saved_summary))
  if ("PERMANOVA_R2" %in% names(saved_summary))
    stopifnot(abs(extract_adonis_stats(fit)$R2 - saved_summary$PERMANOVA_R2[1L]) < 1e-10)
  expected_features <- if (definition$Analysis == "BrayCurtis_Genus") 94L else 623L
  stopifnot(saved_summary$N_features[1L] == expected_features)
  result <- list(metadata = md, biopsy_metadata = md, profiles = NULL, distance = distance,
    permanova = fit, dispersion = dispersion, dispersion_test = dispersion_test,
    permutations = permutations,
    summary = summary_row(definition$Analysis, MATCHED_LABEL, md, md,
      expected_features, fit, dispersion_test, permutations,
      "Individual exact-site biopsies; synchronous pair swaps within patient", pairs = 12L))
  target <- file.path(audit_dir, definition$Analysis, "EXACT_SITE_MATCHED")
  dir.create(target, recursive = TRUE)
  result$audit_directory <- target
  safe_write_csv(data.frame(Distance_source = normalizePath(input, winslash = "/"),
    Saved_distances_reused = TRUE, Saved_P_values_reused = FALSE,
    Metadata_source = matched_manifest_path,
    Matched_profiles_averaged = FALSE, Observed_included = TRUE,
    Alternative_arrangements = 1023L, Enumerated_P = enumerated_p),
    file.path(target, "matched_provenance.csv"))
  result
}


# Close each graphics device before validating and writing the final export.
# Keep plot evaluation lazy: arranged grobs must be built on the export device.
save_beta_image <- function(filename, plot, ...) {
  extension <- tolower(tools::file_ext(filename))
  temporary <- tempfile(pattern = ".beta-plot-", tmpdir = dirname(filename),
                        fileext = paste0(".", extension))
  on.exit(unlink(temporary), add = TRUE)
  args <- c(list(filename = temporary), list(...))
  if (extension %in% c("png", "jpg", "jpeg") && capabilities("cairo")) {
    args$type <- if (extension == "png") "cairo-png" else "cairo"
  }
  do.call(function(...) ggplot2::ggsave(plot = plot, ...), args)
  bytes <- readBin(temporary, what = "raw", n = file.info(temporary)$size)
  if (length(bytes) < 100L) stop("Empty or incomplete plot export: ", filename)
  complete <- switch(extension,
    pdf = grepl("%%EOF", rawToChar(tail(bytes, 20L))),
    png = identical(tail(bytes, 12L), as.raw(c(0, 0, 0, 0, 73, 69, 78, 68, 174, 66, 96, 130))),
    jpg = identical(tail(bytes, 2L), as.raw(c(255, 217))),
    jpeg = identical(tail(bytes, 2L), as.raw(c(255, 217))), TRUE)
  if (!isTRUE(complete)) stop("Incomplete plot export: ", filename)
  writeBin(bytes, filename)
  invisible(filename)
}

make_patient_profiles <- function(counts, metadata, method) {
  if (method == "Aitchison") {
    z <- log(counts + 0.5)
    biopsy_profiles <- z - rowMeans(z)
  } else {
    biopsy_profiles <- counts / rowSums(counts)
  }
  keys <- paste(metadata$Group, metadata$Patient, sep = "::")
  key_levels <- unique(keys)
  n_biopsies <- tabulate(match(keys, key_levels), length(key_levels))
  profiles <- rowsum(biopsy_profiles, keys, reorder = FALSE) / n_biopsies
  pm <- metadata[match(key_levels, keys), c("Patient", "Group"), drop = FALSE]
  pm$SampleID <- key_levels
  pm$PairID <- as.character(pm$Patient)
  pm$N_biopsies <- n_biopsies
  pm$Group <- droplevels(factor(pm$Group))
  rownames(pm) <- key_levels
  metadata$N_biopsies <- 1L
  distance_for <- function(x) if (method == "Aitchison") stats::dist(x) else
    vegan::vegdist(x, method = "bray")
  stopifnot(identical(rownames(profiles), pm$SampleID),
            sum(pm$N_biopsies) == nrow(metadata))
  list(metadata = pm, profiles = profiles, distance = distance_for(profiles),
       biopsy_metadata = metadata, biopsy_profiles = biopsy_profiles,
       biopsy_distance = distance_for(biopsy_profiles))
}

descriptive_partition <- function(distance, group) {
  d2 <- as.matrix(distance)^2
  n <- nrow(d2)
  group <- droplevels(factor(group))
  total <- sum(d2) / (2 * n)
  within <- sum(vapply(split(seq_len(n), group), function(i)
    sum(d2[i, i, drop = FALSE]) / (2 * length(i)), numeric(1)))
  between <- total - within
  if (!is.finite(total) || total <= 0) stop("Distance matrix has no variation.")
  data.frame(Df = nlevels(group) - 1L, SumOfSqs = between,
    R2 = between / total,
    F = (between / (nlevels(group) - 1L)) / (within / (n - nlevels(group))),
    "Pr(>F)" = NA_real_, check.names = FALSE, row.names = "Model")
}

format_final_p <- function(p, label = "P") {
  if (!is.finite(p)) return(paste0(label, ": not tested"))
  paste0(label, if (p < 0.0001) " < 0.0001" else
    paste0(" = ", sprintf("%.4f", p)))
}

read_exported_counts <- function(filename) {
  x <- read.csv(file.path(EXPORTED_INPUT_DIR, filename), check.names = FALSE)
  stopifnot(names(x)[1L] == "SampleID", !anyDuplicated(x$SampleID))
  ans <- as.matrix(x[, -1L, drop = FALSE])
  storage.mode(ans) <- "double"
  rownames(ans) <- as.character(x$SampleID)
  stopifnot(all(is.finite(ans)), all(ans >= 0), all(rowSums(ans) > 0),
            !anyDuplicated(colnames(ans)))
  ans
}

read_main_phyloseq <- function(path) {
  if (INPUT_MODE == "phyloseq") return(readRDS(path))
  x <- read_exported_counts("ASV_counts.csv")
  md <- read.csv(file.path(EXPORTED_INPUT_DIR, "metadata.csv"),
                 check.names = FALSE, colClasses = "character")
  md <- md[match(rownames(x), md$SampleID), , drop = FALSE]
  stopifnot(identical(rownames(x), md$SampleID), !anyNA(md$Patient))
  rownames(md) <- md$SampleID
  md$SampleID <- NULL
  md$Sample_Type <- md$Group
  phyloseq::phyloseq(phyloseq::otu_table(x, taxa_are_rows = FALSE),
                    phyloseq::sample_data(md))
}

patient_swap_matrix <- function(metadata, n_permutations = 9999L) {
  needed <- c("SampleID", "Patient", "Group", "PairID")
  stopifnot(all(needed %in% names(metadata)), !anyNA(metadata[, needed]),
            !anyDuplicated(metadata$SampleID), length(unique(metadata$Group)) == 2L)
  pair_rows <- split(seq_len(nrow(metadata)), metadata$PairID)
  stopifnot(all(lengths(pair_rows) == 2L))
  for (ii in pair_rows) {
    stopifnot(length(unique(metadata$Patient[ii])) == 1L,
              length(unique(metadata$Group[ii])) == 2L)
  }
  patients <- sort(unique(as.character(metadata$Patient)))
  n_patients <- length(patients)
  if (n_patients < 2L) stop("At least two paired patients are needed.")
  n_possible <- 2^n_patients - 1
  if (n_possible <= n_permutations) {
    swaps <- outer(seq_len(n_possible), seq_len(n_patients) - 1L,
                   function(a, b) bitwAnd(bitwShiftR(as.integer(a), b), 1L))
  } else {
    set.seed(1)
    swaps <- matrix(0L, nrow = n_permutations, ncol = n_patients)
    for (i in seq_len(n_permutations)) {
      repeat {
        z <- sample.int(2L, n_patients, replace = TRUE) - 1L
        if (any(z)) break
      }
      swaps[i, ] <- z
    }
  }
  permutations <- matrix(rep(seq_len(nrow(metadata)), each = nrow(swaps)),
                         nrow = nrow(swaps))
  for (ii in pair_rows) {
    j <- match(as.character(metadata$Patient[ii[1L]]), patients)
    selected <- which(swaps[, j] == 1L)
    permutations[selected, ii[1L]] <- ii[2L]
    permutations[selected, ii[2L]] <- ii[1L]
  }
  stopifnot(all(apply(permutations, 1L, function(z) identical(sort(z), seq_len(nrow(metadata))))))
  attr(permutations, "scheme") <- "Synchronous lesion-adjacent swaps within patient; exact-site pairs preserved"
  permutations
}

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
          N_biopsies
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

  save_beta_image(
    file.path(out_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    device = pdf_device,
    bg = "white"
  )

  save_beta_image(
    file.path(out_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    dpi = export_dpi,
    bg = "white"
  )

  save_beta_image(
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

make_lineage_aware_genus <- function(tax_df) {
  # Preserve the genus branch's original missing-taxonomy handling.
  clean_character <- function(x) {
    x <- trimws(as.character(x))
    x[x %in% c("", "NA", "na", "N/A")] <- NA_character_
    x
  }
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

export_primary_figures <- function(results) {
  # Every panel goes through the SAME plotting function and export function.
  # All ordinations retain their own observations, coordinates and axis scales.
  panel_order <- c("H_ADJ_SSL_HP vs H_COL", "SSL_HP vs H_COL", MATCHED_LABEL)
  stems <- c("PCoA_CLR_H_ADJ_SSL_HP_vs_H_COL", "PCoA_CLR_SSL_HP_vs_H_COL",
             "PCoA_Aitchison_EXACT_MATCHED")
  plots <- setNames(vector("list", length(panel_order)), panel_order)
  for (i in seq_along(panel_order)) {
    result <- results$Aitchison_ASV[[panel_order[i]]]
    if (identical(panel_order[i], MATCHED_LABEL)) {
      ord <- vegan::wcmdscale(result$distance, k = 2, eig = TRUE,
                              add = FALSE, x.ret = TRUE)
      coordinates <- as.data.frame(ord$points, check.names = FALSE)
      names(coordinates) <- c("PCoA1", "PCoA2")
      coordinates$SampleID <- rownames(coordinates)
      coordinates <- dplyr::left_join(coordinates, result$metadata, by = "SampleID")
      variance_explained <- ord$eig[1:2] / sum(ord$eig[ord$eig > 0])
      groups_present <- c("H_ADJ_SSL_HP", "SSL_HP")
    } else {
      ord <- get_pcoa(result$biopsy_distance, result$biopsy_metadata)
      coordinates <- ord$coordinates
      variance_explained <- ord$variance_explained
      groups_present <- c("H_COL", if (i == 1L) "H_ADJ_SSL_HP" else "SSL_HP")
    }
    plot <- make_beta_plot(coordinates, variance_explained,
      paste(groups_present, collapse = " vs "),
      comparison_annotation(result$summary), groups_present)
    save_beta_plot(plot, stems[i], width = 5.7, height = 4.4)
    safe_write_csv(coordinates,
      file.path(result$audit_directory, "plotted_biopsy_coordinates.csv"))
    plots[[i]] <- plot
    for (extension in c("pdf", "png", "jpg")) {
      source_file <- file.path(primary_dir, paste0(stems[i], ".", extension))
      final_file <- file.path(primary_dir,
        paste0("Figure_4_3", LETTERS[i], ".", extension))
      if (!file.copy(source_file, final_file, overwrite = FALSE))
        stop("Could not create final figure file: ", final_file)
      if (!identical(unname(tools::md5sum(source_file)),
                     unname(tools::md5sum(final_file))))
        stop("Final figure copy did not match: ", final_file)
    }
  }
  panels <- lapply(seq_along(plots), function(i)
    plots[[i]] + ggplot2::labs(tag = LETTERS[i]) +
      ggplot2::theme(plot.tag.position = c(0.01, 0.99)))
  save_beta_plot(gridExtra::arrangeGrob(grobs = panels, nrow = 1),
    "Figure_Aitchison_ASV_pairwise_3panel", width = 17.1, height = 4.6)
  writeLines(c(PLOT_REVISION,
    "Use these three files together from this run's 01_PRIMARY_AITCHISON folder:",
    "Figure_4_3A.jpg : H_COL vs H_ADJ_SSL_HP",
    "Figure_4_3B.jpg : H_COL vs SSL_HP",
    "Figure_4_3C.jpg : H_ADJ_SSL_HP vs SSL_HP; 24 exact-site biopsies",
    "The corresponding PNG and PDF files are also supplied.",
    "All three use solid points, a right-hand legend, the same 68% ellipse style,",
    "and R-squared, raw PERMANOVA P, Holm P and betadisper P inside the plot.",
    "All three standalone exports are 5.7 x 4.4 inches; PNG/JPG use 600 dpi."
  ), file.path(primary_dir, "USE_THESE_FINAL_FIGURES.txt"))
  invisible(plots)
}


# ---------------- Main workflow ----------------
run_analysis <- function() {
  phy <- read_main_phyloseq(phyloseq_input)
  if (!inherits(phy, "phyloseq")) stop("The input is not a phyloseq object.")
  if (!"Sample_Type" %in% phyloseq::sample_variables(phy))
    stop("The metadata have no Sample_Type column.")
  group <- clean_character(phyloseq::sample_data(phy)$Sample_Type)
  group[group %in% c("SSL", "HP", "SSL_HP")] <- "SSL_HP"
  group[group %in% c("H_ADJ_SSL", "H_ADJ_HP", "H_ADJ_SSL_HP")] <- "H_ADJ_SSL_HP"
  group[group %in% c("H_COL", "H_COL_D", "H_COL_ND")] <- "H_COL"
  group[!group %in% group_levels] <- NA_character_
  phyloseq::sample_data(phy)$Group <- factor(group, levels = group_levels)
  phy <- phyloseq::prune_samples(!is.na(group), phy)
  phy <- phyloseq::prune_taxa(phyloseq::taxa_sums(phy) > 0, phy)
  asv_counts <- counts_from_phy(phy)
  metadata <- get_meta(phy)
  metadata <- metadata[rownames(asv_counts), , drop = FALSE]
  stopifnot(identical(rownames(asv_counts), as.character(metadata$SampleID)),
            nrow(metadata) == 64L, ncol(asv_counts) == 797L,
            length(unique(metadata$Patient)) == 31L)
  if (any(as.character(metadata$Patient) == "7699") ||
      any(grepl("7699", metadata$SampleID)))
    stop("Patient 7699 is still present. Use the original without7699 input.")
  # Patient 7669 is NOT excluded.
  expected_groups <- c(H_COL = 14L, H_ADJ_SSL_HP = 12L, SSL_HP = 38L)
  stopifnot(identical(as.integer(table(metadata$Group)), unname(expected_groups)))
  patients <- lapply(group_levels, function(g)
    unique(as.character(metadata$Patient[metadata$Group == g])))
  names(patients) <- group_levels
  stopifnot(length(patients$H_COL) == 8L, length(patients$H_ADJ_SSL_HP) == 10L,
            length(patients$SSL_HP) == 23L,
            !length(intersect(patients$H_COL, patients$SSL_HP)),
            !length(intersect(patients$H_COL, patients$H_ADJ_SSL_HP)),
            all(patients$H_ADJ_SSL_HP %in% patients$SSL_HP))
  cat("\nLoaded input: 64 biopsies, 31 patients, 797 ASVs.\n")
  print(table(metadata$Group))
  write_matrix(asv_counts, file.path(audit_dir, "ASV_counts_all_64.csv"))
  safe_write_csv(metadata, file.path(audit_dir, "sample_manifest_all_64.csv"))

  if (INPUT_MODE == "phyloseq") {
    tax_df <- as.data.frame(phyloseq::tax_table(phy),
                           check.names = FALSE, stringsAsFactors = FALSE)
    if (!"Genus" %in% names(tax_df)) stop("The taxonomy table has no Genus column.")
    genus_result <- make_lineage_aware_genus(tax_df)
    tax_df$Genus <- genus_result$genus
    phyloseq::tax_table(phy) <- phyloseq::tax_table(as.matrix(tax_df))
    safe_write_csv(genus_result$mapping, file.path(audit_dir, "Genus_lineage_mapping.csv"))
    phy_genus <- phyloseq::tax_glom(phy, taxrank = "Genus", NArm = FALSE)
    phy_genus <- phyloseq::prune_taxa(phyloseq::taxa_sums(phy_genus) > 0, phy_genus)
    genus_labels <- as.character(phyloseq::tax_table(phy_genus)[, "Genus"])
    stopifnot(!anyNA(genus_labels), all(nzchar(genus_labels)))
    phyloseq::taxa_names(phy_genus) <- make.unique(genus_labels)
    genus_counts <- counts_from_phy(phy_genus)
  } else {
    # Optional validation mode: use the original already-agglomerated count matrix.
    genus_counts <- read_exported_counts("Genus_counts.csv")
  }
  stopifnot(setequal(rownames(genus_counts), metadata$SampleID))
  genus_counts <- genus_counts[metadata$SampleID, , drop = FALSE]
  stopifnot(ncol(genus_counts) == 106L,
            all(rowSums(genus_counts) == rowSums(asv_counts)))
  write_matrix(genus_counts, file.path(audit_dir, "Genus_counts_all_64.csv"))
  matched_metadata <- read.csv(matched_manifest_path, check.names = FALSE,
                               colClasses = "character")
  stopifnot(all(c("SampleID", "Patient", "Group", "PairID", "Site") %in%
                names(matched_metadata)), !anyDuplicated(matched_metadata$SampleID))
  safe_write_csv(matched_metadata, file.path(audit_dir, "exact_site_sample_manifest.csv"))

  results <- setNames(vector("list", length(ANALYSES)), ANALYSES)
  for (analysis in ANALYSES) {
    cat("\nRunning ", analysis, ": two independent-control tests and exact-site test.\n", sep = "")
    x <- if (analysis == "BrayCurtis_Genus") genus_counts else asv_counts
    tests <- list(
      run_control(x, metadata, analysis, "SSL_HP"),
      run_control(x, metadata, analysis, "H_ADJ_SSL_HP"),
      run_exact_site(matched_definitions[matched_definitions$Analysis == analysis, ],
                     matched_metadata, metadata)
    )
    names(tests) <- COMPARISONS
    raw_p <- vapply(tests, function(z) z$summary$P_PERMANOVA_raw, numeric(1))
    stopifnot(length(raw_p) == 3L, all(is.finite(raw_p)),
              all(raw_p >= 0 & raw_p <= 1))
    adjusted_p <- stats::p.adjust(raw_p, method = "holm")
    for (j in seq_along(tests)) {
      tests[[j]]$summary$P_PERMANOVA_Holm <- unname(adjusted_p[j])
      tests[[j]]$summary$Holm_family <- paste(analysis, "three tissue comparisons")
      write_test_audit(tests[[j]], tests[[j]]$audit_directory, tests[[j]]$counts)
      safe_write_csv(tests[[j]]$summary,
                     file.path(tests[[j]]$audit_directory, "final_statistics.csv"))
    }
    results[[analysis]] <- tests
  }
  all_results <- dplyr::bind_rows(lapply(results, function(tests)
    dplyr::bind_rows(lapply(tests, function(z) z$summary))))
  stopifnot(nrow(all_results) == 9L, !anyNA(all_results$P_PERMANOVA_Holm))
  safe_write_csv(all_results, file.path(audit_dir, "ALL_NINE_RESULTS.csv"))
  primary <- all_results[all_results$Analysis == "Aitchison_ASV", , drop = FALSE]
  supplementary <- all_results[all_results$Analysis != "Aitchison_ASV", , drop = FALSE]
  safe_write_csv(primary, file.path(primary_dir, "PRIMARY_AITCHISON_THREE_COMPARISONS.csv"))
  safe_write_csv(supplementary,
    file.path(supplement_dir, "SUPPLEMENTARY_BRAY_CURTIS_SIX_TESTS.csv"))
  publication_table <- function(x) {
    data.frame(
      Metric = x$Analysis, Comparison = x$Comparison, Biopsies = x$N_biopsies,
      Patients = x$N_patients, Test_profiles = x$N_test_profiles,
      R_squared = sprintf("%.4f", x$R2_PERMANOVA),
      Pseudo_F = sprintf("%.3f", x$F_PERMANOVA),
      P_raw = sprintf("%.4f", x$P_PERMANOVA_raw),
      P_Holm_3 = sprintf("%.4f", x$P_PERMANOVA_Holm),
      Dispersion_P_unadjusted = sprintf("%.4f", x$P_betadisper),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }
  safe_write_csv(publication_table(primary),
    file.path(primary_dir, "PRIMARY_TABLE_FORMATTED.csv"))
  safe_write_csv(publication_table(supplementary),
    file.path(supplement_dir, "SUPPLEMENTARY_TABLE_FORMATTED.csv"))

  # Numerical references from the 16 October completed run. Checking only:
  # these numbers NEVER replace a fitted statistic or a calculated P value.
  reference <- data.frame(
    Analysis = rep(ANALYSES, each = 3L), Comparison = rep(COMPARISONS, 3L),
    Previous_R2 = c(.0433500552076952, .0699562848183679, .0157514075356653,
                    .0367476732438766, .0647855016421790, .0114493434067659,
                    .0806347226102084, .1208711762628840, .0133053356734541),
    Previous_raw_P = c(.0183, .0035, .349609375,
                       .1333, .1324, .4453125,
                       .0034, .0084, .625),
    Previous_dispersion_P = c(.2259, .5812, .5400390625,
                              .6471, .6798, .8486328125,
                              .7578, .8272, .8564453125)
  )
  stopifnot(identical(reference$Analysis, all_results$Analysis),
            identical(reference$Comparison, all_results$Comparison))
  reference$Current_R2 <- all_results$R2_PERMANOVA
  reference$Current_raw_P <- all_results$P_PERMANOVA_raw
  reference$Current_dispersion_P <- all_results$P_betadisper
  reference$R2_matches <- abs(reference$Previous_R2 - reference$Current_R2) < 1e-9
  reference$Raw_P_matches <- abs(reference$Previous_raw_P - reference$Current_raw_P) < 1e-12
  reference$Dispersion_P_matches <-
    abs(reference$Previous_dispersion_P - reference$Current_dispersion_P) < 1e-12
  safe_write_csv(reference, file.path(audit_dir, "CHECK_AGAINST_PREVIOUS_RUN.csv"))
  if (CHECK_PREVIOUS_RUN && !all(reference$R2_matches))
    stop("An R-squared differs from the previous run. Inspect CHECK_AGAINST_PREVIOUS_RUN.csv; ",
         "check input, taxon set and profile construction before using this run.")
  if (CHECK_PREVIOUS_RUN && (!all(reference$Raw_P_matches) ||
                             !all(reference$Dispersion_P_matches)))
    warning("Some permutation P values differ from the previous run. ",
      "The newly computed results have been retained, not replaced. ",
      "Inspect CHECK_AGAINST_PREVIOUS_RUN.csv, saved permutations, sample order, ",
      "RNG kind and package versions before interpreting any difference.")

  # PRIMARY FIGURES ONLY. Both control plots show their individual biopsies,
  # with inference performed on patient-tissue means, as stated in the captions.
  plots <- export_primary_figures(results)
  exact <- results$Aitchison_ASV[[MATCHED_LABEL]]

  # Explicit audit of inclusion: every eligible biopsy enters a control test.
  included <- unique(unlist(lapply(results$Aitchison_ASV[CONTROL_LABELS],
                                  function(z) as.character(z$biopsy_metadata$SampleID))))
  stopifnot(setequal(included, metadata$SampleID))
  inclusion <- metadata[, c("SampleID", "Patient", "Group"), drop = FALSE]
  inclusion$Lesion_control <- inclusion$Group %in% c("SSL_HP", "H_COL")
  inclusion$Adjacent_control <- inclusion$Group %in% c("H_ADJ_SSL_HP", "H_COL")
  inclusion$Exact_site <- inclusion$SampleID %in% exact$metadata$SampleID
  safe_write_csv(inclusion, file.path(audit_dir, "SAMPLE_INCLUSION_BY_COMPARISON.csv"))
  input_paths <- if (INPUT_MODE == "phyloseq") phyloseq_input else
    file.path(EXPORTED_INPUT_DIR, c("ASV_counts.csv", "Genus_counts.csv", "metadata.csv"))
  input_paths <- c(input_paths, matched_input_files)
  safe_write_csv(data.frame(Path = input_paths, MD5 = unname(tools::md5sum(input_paths))),
                 file.path(audit_dir, "INPUT_FILE_CHECKSUMS.csv"))
  save_beta_rds(results, file.path(audit_dir, "ALL_ANALYSIS_OBJECTS.rds"))
  writeLines(capture.output(sessionInfo()), file.path(audit_dir, "SESSION_INFO.txt"))
  writeLines(capture.output(RNGkind()), file.path(audit_dir, "RNG_KIND.txt"))
  writeLines(c(
    "FINAL THREE-COMPARISON BETA-DIVERSITY OUTPUT",
    "",
    "Use 01_PRIMARY_AITCHISON/Figure_4_3A, Figure_4_3B and Figure_4_3C (PNG/PDF/JPG).",
    "Combined version: 01_PRIMARY_AITCHISON/Figure_Aitchison_ASV_pairwise_3panel.pdf.",
    "A: adjacent versus healthy colon; B: lesion versus healthy colon;",
    "C: EXACT-SITE MATCHED lesion versus adjacent (not the earlier all-50 comparison).",
    "Every standalone panel uses the original pairwise size: 5.7 x 4.4 inches.",
    "No patient/biopsy counts were added to figures. The counts are in the tables.",
    "",
    "PRIMARY_AITCHISON_THREE_COMPARISONS.csv has complete numerical results.",
    "PRIMARY_TABLE_FORMATTED.csv is the presentation version (four-decimal P values).",
    "02_SUPPLEMENTARY_BRAY_CURTIS/SUPPLEMENTARY_BRAY_CURTIS_SIX_TESTS.csv",
    "has all six sensitivity tests: three ASV and three genus Bray-Curtis.",
    "SUPPLEMENTARY_TABLE_FORMATTED.csv is the presentation version.",
    "No Bray-Curtis figures, overall comparison, or full-50 lesion-adjacent test are exported.",
    "",
    "PATIENT HANDLING",
    "All 64 eligible biopsies contribute to at least one control comparison.",
    "Lesion/control: 52 biopsies, 31 independent patient-tissue mean profiles.",
    "Adjacent/control: 26 biopsies, 18 independent patient-tissue mean profiles.",
    "Control tests use 9,999 unrestricted permutations BETWEEN independent patient means.",
    "Control plots retain individual biopsies for visual continuity; they are not",
    "the independent observational units used in the PERMANOVA tests.",
    "Exact-site: 24 biopsies, 12 pairs, 10 patients; NO averaging across sites.",
    "Every patient's site pairs swap synchronously: all 1,023 alternatives plus",
    "the observed arrangement. P relies on joint within-patient exchangeability.",
    "This is a restricted permutation test, not a linear mixed-effects model.",
    "",
    "TRANSFORMATIONS",
    "Aitchison: per-biopsy log(count+0.5), subtract its row mean, then average",
    "CLR profiles within patient and tissue and use Euclidean distances.",
    "Bray-Curtis: divide each biopsy by its own library size, average relative",
    "abundance profiles within patient and tissue, then use Bray-Curtis distances.",
    "Genus Bray-Curtis first uses the unchanged lineage-aware genus aggregation.",
    "The control analyses retain the full-cohort feature set: 797 ASVs or 106 genera.",
    "Exact-site analyses reuse their original saved distances: 623 ASVs or 94 genera.",
    "This preserves the original matched preprocessing; those distances are not",
    "silently replaced by distances from the full-cohort feature set.",
    "",
    "MULTIPLE TESTING AND DIAGNOSTICS",
    "Holm correction across the THREE tissue tests separately within each metric.",
    "There is no cross-metric multiplicity adjustment; Bray-Curtis is sensitivity",
    "analysis, not additional independent confirmation. Report discordance.",
    "Matched raw P is preserved, but its Holm P is now part of the three-test family.",
    "Dispersion P values are unadjusted diagnostics, not part of that Holm family.",
    "Controls: betadisper median, bias.adjust=TRUE; same free patient permutations.",
    "Exact-site: betadisper median, bias.adjust=FALSE; synchronous patient swaps.",
    "The matched dispersion check is labelled a restricted diagnostic.",
    "R-squared uses the full distance matrix, never only the two plotted axes.",
    "No fixed Patient covariate or patient random intercept is fitted in PERMANOVA.",
    "",
    "REVIEW",
    "Check 03_AUDIT/CHECK_AGAINST_PREVIOUS_RUN.csv and WARNINGS.txt.",
    "Every P value is computed anew. Reference values are never inserted as results.",
    "Warnings are captured for review; RUN_COMPLETE.txt means exports finished,",
    "not that all warnings are necessarily ignorable.",
    "Archive the entire folder with its exact script, input checksums, session",
    "information, permutation matrices, counts, profiles and test objects.",
    "Metric choice was not made prospectively in this rerun; do not call it preregistered."
  ), file.path(run_dir, "READ_ME_FIRST.txt"))
  writeLines(c(
    "DRAFT METHODS -- reconcile with the final chapter wording",
    "",
    "ASV-level Aitchison distance was the primary beta-diversity measure.",
    "For each biopsy, counts were augmented by a pseudocount of 0.5 and",
    "centred-log-ratio transformed. In comparisons of lesion or adjacent mucosa",
    "with healthy colon, transformed profiles were averaged within patient",
    "and tissue group before calculating Euclidean distances. Each of these",
    "comparisons contained one mean profile per independent patient.",
    "PERMANOVA used 9,999 unrestricted permutations of these patient profiles.",
    "Individual biopsies were retained in ordinations for visualisation.",
    "",
    "The exact-site lesion-adjacent comparison retained 24 biopsies forming",
    "12 site pairs from ten patients, with no averaging across sites.",
    "The original matched distance matrices were retained. All site pairs",
    "from the same patient were swapped synchronously, enumerating 1,023",
    "alternative arrangements plus the observed arrangement. Validity assumes",
    "joint exchangeability of paired tissue labels within each patient.",
    "",
    "ASV- and genus-level Bray-Curtis analyses were reported as sensitivity",
    "analyses. Biopsies were converted to relative abundances before averaging",
    "within patient and tissue for the two control comparisons. The matched",
    "analyses retained individual biopsy profiles and synchronous swaps.",
    "Holm adjustment was applied across the three tissue comparisons within",
    "each metric. Multivariate dispersion was assessed as a diagnostic;",
    "dispersion P values were unadjusted."
  ), file.path(run_dir, "METHODS_DRAFT.txt"))
  writeLines(c(
    "FIGURE CAPTION DETAILS TO RETAIN",
    "Three ASV-level Aitchison ordinations: adjacent/control, lesion/control,",
    "and exact-site matched lesion/adjacent. Points represent individual biopsies.",
    "In the control comparisons, PERMANOVA is based on patient-tissue mean",
    "CLR profiles, not independent-biopsy permutations.",
    "In the exact-site comparison, 12 site pairs from ten patients are retained",
    "and all pairs from the same patient are swapped synchronously.",
    "The shown R-squared and P values refer to these explicitly stated tests;",
    "R-squared is not the percentage explained by the two ordination axes.",
    "Holm P values are adjusted across these three comparisons."
  ), file.path(run_dir, "FIGURE_CAPTION_NOTES.txt"))
  if (!is.na(script_to_archive)) {
    if (!file.copy(script_to_archive,
                   file.path(run_dir, "RUN_BETA_FINAL_PATIENT_LEVEL.R"), overwrite = FALSE))
      warning("The script could not be copied into the output folder; include it manually.")
  } else warning("Source-script path was not detected; include this R script manually in the archive.")
  invisible(all_results)
}

completed_results <- tryCatch(
  withCallingHandlers(run_analysis(), warning = function(w) {
    warnings_seen <<- c(warnings_seen, conditionMessage(w))
  }),
  error = function(e) {
    try(writeLines(c("RUN INCOMPLETE", conditionMessage(e)),
                   file.path(run_dir, "RUN_FAILED.txt")), silent = TRUE)
    try(writeLines(capture.output(sessionInfo()),
                   file.path(audit_dir, "SESSION_INFO.txt")), silent = TRUE)
    stop(conditionMessage(e), call. = FALSE)
  }
)
writeLines(if (length(warnings_seen)) unique(warnings_seen) else "No warnings were recorded.",
           file.path(run_dir, "WARNINGS.txt"))
writeLines("All three comparisons per metric and requested exports completed.",
           file.path(run_dir, "RUN_COMPLETE.txt"))
writeLines("COMPLETE. Review WARNINGS.txt and CHECK_AGAINST_PREVIOUS_RUN.csv.",
           file.path(run_dir, "RUN_STATUS.txt"))
cat("\nCOMPLETE OUTPUT FOLDER:\n", normalizePath(run_dir, winslash = "/"), "\n")
cat("\nFINAL FIGURES (", PLOT_REVISION, "):\n",
    file.path(primary_dir, c("Figure_4_3A.jpg", "Figure_4_3B.jpg", "Figure_4_3C.jpg")),
    sep = "\n")
print(completed_results[, c("Analysis", "Comparison", "R2_PERMANOVA",
                            "P_PERMANOVA_raw", "P_PERMANOVA_Holm")], row.names = FALSE)
if (requireNamespace("zip", quietly = TRUE)) {
  zip_path <- paste0(run_dir, ".zip")
  zip_ok <- tryCatch({
    zip::zipr(zipfile = zip_path,
      files = list.files(run_dir, all.files = TRUE, no.. = TRUE),
      root = run_dir, include_directories = TRUE)
    TRUE
  }, error = function(e) {
    cat("\nAutomatic ZIP failed: ", conditionMessage(e),
        "\nThe completed output folder can be zipped manually.\n", sep = "")
    FALSE
  })
  if (zip_ok) cat("\nZIP TO UPLOAD:\n", normalizePath(zip_path, winslash = "/"), "\n")
} else cat("\nZip this complete output folder manually and upload it for chapter review.\n")
invisible(completed_results)
})
