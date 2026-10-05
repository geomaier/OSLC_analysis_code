#!/usr/bin/env Rscript
# Patient-level sensitivity checks for the two H_COL comparisons.

run_patient_check <- function(
    work_dir,
    asv_dir = file.path(work_dir, "Beta_diversity_BrayCurtis_ASV_FINAL"),
    genus_dir = file.path(work_dir, "Beta_diversity_BrayCurtis_Genus_FINAL")
) {
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("The vegan package is required.")
  }

  work_dir <- normalizePath(work_dir, mustWork = TRUE)
  inputs <- c(
    asv_manifest = file.path(asv_dir, "00_sample_manifest.csv"),
    asv_counts = file.path(asv_dir, "ASV_count_matrix_all_64_samples.csv"),
    genus_manifest = file.path(genus_dir, "Sample_manifest_BrayCurtis_Genus.csv"),
    genus_counts = file.path(genus_dir, "Genus_count_matrix_all_64_samples.csv")
  )
  if (!all(file.exists(inputs))) {
    stop("Missing input files:\n", paste(inputs[!file.exists(inputs)], collapse = "\n"))
  }

  read_manifest <- function(path) {
    m <- read.csv(path, check.names = FALSE, colClasses = "character")
    needed <- c("SampleID", "Patient", "Group")
    stopifnot(all(needed %in% names(m)))
    m <- m[, needed, drop = FALSE]
    stopifnot(!anyNA(m), !anyDuplicated(m$SampleID))
    stopifnot(all(vapply(m, function(x) all(nzchar(trimws(x))), logical(1))))
    stopifnot(all(m$Group %in% c("H_COL", "SSL_HP", "H_ADJ_SSL_HP")))
    rownames(m) <- m$SampleID
    m
  }

  metadata <- read_manifest(inputs[["asv_manifest"]])
  genus_metadata <- read_manifest(inputs[["genus_manifest"]])
  stopifnot(setequal(metadata$SampleID, genus_metadata$SampleID))
  genus_metadata <- genus_metadata[metadata$SampleID, , drop = FALSE]
  stopifnot(identical(metadata$Patient, genus_metadata$Patient))
  stopifnot(identical(metadata$Group, genus_metadata$Group))

  read_counts <- function(path) {
    x <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
    stopifnot(is.numeric(x), !anyDuplicated(rownames(x)), !anyDuplicated(colnames(x)))
    stopifnot(all(is.finite(x)), all(x >= 0), all(x == floor(x)))
    stopifnot(setequal(rownames(x), metadata$SampleID))
    x <- x[metadata$SampleID, , drop = FALSE]
    stopifnot(all(rowSums(x) > 0))
    x
  }

  asv_counts <- read_counts(inputs[["asv_counts"]])
  genus_counts <- read_counts(inputs[["genus_counts"]])
  stopifnot(identical(as.numeric(rowSums(asv_counts)), as.numeric(rowSums(genus_counts))))

  out_dir <- file.path(work_dir, "Beta_patient_check", format(Sys.time(), "%Y%m%d_%H%M%S"))
  if (dir.exists(out_dir)) stop("Output directory already exists: ", out_dir)
  if (!dir.create(out_dir, recursive = TRUE)) stop("Could not create: ", out_dir)

  write.csv(metadata, file.path(out_dir, "input_samples.csv"), row.names = FALSE)
  write.csv(
    data.frame(Input = names(inputs), Path = unname(inputs), MD5 = unname(tools::md5sum(inputs))),
    file.path(out_dir, "input_files.csv"), row.names = FALSE
  )

  n_permutations <- 9999L
  seed <- 1L
  clr_pseudocount <- 0.5
  analyses <- c("Aitchison_ASV", "BrayCurtis_ASV", "BrayCurtis_Genus")
  comparisons <- c("SSL_HP", "H_ADJ_SSL_HP")
  results <- list()

  for (analysis in analyses) {
    all_counts <- if (analysis == "BrayCurtis_Genus") genus_counts else asv_counts

    for (tissue in comparisons) {
      m <- metadata[metadata$Group %in% c("H_COL", tissue), , drop = FALSE]
      groups_per_patient <- tapply(m$Group, m$Patient, function(z) length(unique(z)))
      if (any(groups_per_patient != 1L)) {
        stop("A patient contributes to both groups in ", tissue, " vs H_COL.")
      }

      x <- all_counts[m$SampleID, , drop = FALSE]
      x <- x[, colSums(x) > 0, drop = FALSE]

      if (analysis == "Aitchison_ASV") {
        log_counts <- log(x + clr_pseudocount)
        biopsy_profiles <- sweep(log_counts, 1, rowMeans(log_counts), "-")
        aggregation <- "Arithmetic mean of biopsy CLR profiles within each patient"
      } else {
        biopsy_profiles <- sweep(x, 1, rowSums(x), "/")
        aggregation <- "Arithmetic mean of biopsy relative abundances within each patient"
      }

      keys <- paste(m$Group, m$Patient, sep = "::")
      key_levels <- unique(keys)
      n_biopsies <- tabulate(match(keys, key_levels), nbins = length(key_levels))
      patient_profiles <- rowsum(biopsy_profiles, group = keys, reorder = FALSE)
      stopifnot(identical(rownames(patient_profiles), key_levels))
      patient_profiles <- sweep(patient_profiles, 1, n_biopsies, "/")

      pm <- m[match(key_levels, keys), c("Patient", "Group"), drop = FALSE]
      pm$ProfileID <- key_levels
      pm$N_biopsies <- n_biopsies
      pm$Group <- factor(pm$Group, levels = c("H_COL", tissue))
      rownames(pm) <- pm$ProfileID
      stopifnot(!anyDuplicated(pm$Patient), all(table(pm$Group) >= 2L))

      distance <- if (analysis == "Aitchison_ASV") {
        stats::dist(patient_profiles, method = "euclidean")
      } else {
        vegan::vegdist(patient_profiles, method = "bray", binary = FALSE)
      }
      stopifnot(identical(attr(distance, "Labels"), rownames(pm)))

      set.seed(seed)
      permanova <- vegan::adonis2(distance ~ Group, data = pm, permutations = n_permutations, by = NULL)
      stopifnot("Model" %in% rownames(permanova))
      stopifnot(all(is.finite(as.numeric(unlist(permanova["Model", c("R2", "F", "Pr(>F)")])))))
      dispersion <- vegan::betadisper(distance, group = pm$Group, type = "median", bias.adjust = FALSE)
      set.seed(seed)
      dispersion_test <- vegan::permutest(dispersion, permutations = n_permutations)

      label <- paste0(analysis, "_", tissue, "_vs_H_COL")
      results[[label]] <- data.frame(
        Analysis = analysis,
        Comparison = paste(tissue, "vs H_COL"),
        Aggregation = aggregation,
        N_biopsies = nrow(m),
        N_patients = nrow(pm),
        N_tissue_patients = sum(pm$Group == tissue),
        N_H_COL_patients = sum(pm$Group == "H_COL"),
        N_features = ncol(patient_profiles),
        R2 = permanova["Model", "R2"],
        F = permanova["Model", "F"],
        P_raw = permanova["Model", "Pr(>F)"],
        P_betadisper_raw = dispersion_test$tab[1, "Pr(>F)"],
        Permutations_requested = n_permutations,
        Seed = seed,
        stringsAsFactors = FALSE
      )

      write.csv(pm, file.path(out_dir, paste0(label, "_patients.csv")), row.names = FALSE)
      saveRDS(
        list(patient_metadata = pm, patient_profiles = patient_profiles,
             distance = distance, permanova = permanova,
             dispersion = dispersion, dispersion_test = dispersion_test),
        file.path(out_dir, paste0(label, ".rds"))
      )
    }
  }

  summary <- do.call(rbind, results)
  rownames(summary) <- NULL
  summary$P_Holm_conservative_family3 <- NA_real_
  for (analysis in analyses) {
    selected <- summary$Analysis == analysis
    summary$P_Holm_conservative_family3[selected] <- stats::p.adjust(
      summary$P_raw[selected], method = "holm", n = 3L
    )
  }
  write.csv(summary, file.path(out_dir, "patient_results.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  writeLines(c(
    "Patient-level sensitivity analysis of SSL_HP vs H_COL and H_ADJ_SSL_HP vs H_COL.",
    "One profile per patient is used separately within each two-group comparison.",
    "Every biopsy has equal weight when calculating its patient's profile.",
    "Aitchison: remove ASVs absent from the two-group subset; add 0.5 to biopsy counts;",
    "calculate biopsy CLR profiles; average CLR profiles within patient; Euclidean distance.",
    "Bray-Curtis: calculate biopsy relative abundances; average within patient; Bray-Curtis distance.",
    "Existing 5% prevalence filtering is retained. No new prevalence or read-depth cutoff is applied.",
    "PERMANOVA and dispersion tests use patient profiles and request 9,999 free permutations.",
    "The seed is 1. Dispersion uses the spatial median and the original default bias.adjust = FALSE.",
    "P_Holm_conservative_family3 uses Holm with n = 3 separately for each distance/feature analysis.",
    "Only two comparisons are run; the untested third comparison is conservatively reserved in the family.",
    "Dispersion P values are unadjusted diagnostics, as in the original analysis.",
    "This analysis estimates differences between patients' average sampled microbial profiles.",
    "The overall three-group test and lesion-adjacent comparison are not tested by this script."
  ), file.path(out_dir, "analysis_notes.txt"))

  print(summary[, c("Analysis", "Comparison", "N_patients", "R2", "P_raw",
                    "P_Holm_conservative_family3", "P_betadisper_raw")], row.names = FALSE)
  cat("\nSaved patient-level check to:\n", out_dir, "\n", sep = "")
  invisible(list(results = summary, output_directory = out_dir))
}

source("config/config.R")
stopifnot(exists("PATH_WORK"))
if (!exists("BETA_ASV_EXPORT_DIR")) {
  BETA_ASV_EXPORT_DIR <- file.path(PATH_WORK, "Beta_diversity_BrayCurtis_ASV_FINAL")
}
if (!exists("BETA_GENUS_EXPORT_DIR")) {
  BETA_GENUS_EXPORT_DIR <- file.path(PATH_WORK, "Beta_diversity_BrayCurtis_Genus_FINAL")
}
run_patient_check(PATH_WORK, BETA_ASV_EXPORT_DIR, BETA_GENUS_EXPORT_DIR)
