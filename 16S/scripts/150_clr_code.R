#!/usr/bin/env Rscript
# Targeted genus CLR and detection analyses, with prevalence-filter sensitivity.
suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(stringr)
  library(grid)
})

if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Package 'gridExtra' is required. Install it with install.packages('gridExtra').")
}

options(stringsAsFactors = FALSE, contrasts = c("contr.treatment", "contr.poly"))
for (pkg in c("lme4", "lmerTest")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg)
set.seed(20260821)

config <- new.env(parent = baseenv())
if (file.exists("config/config.R")) sys.source("config/config.R", envir = config)
work_dir <- if (exists("PATH_WORK", config, inherits = FALSE)) config$PATH_WORK else "outputs"
if (!dir.exists(work_dir)) stop("Input directory does not exist: ", work_dir)

out_root <-
  file.path(work_dir, "Section4_6_candidate_genera_ALL3_FULL_GROUP_COMPARISONS_FINAL")

out_root <- file.path(out_root, paste0("corrected_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
if (dir.exists(out_root)) stop("Output folder already exists: ", out_root)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

analysis_specs <- list(
  list(
    id = "01_PRIMARY_5pct_filtered_without7699",
    label = "Primary 5% ASV-prevalence-filtered analysis",
    prevalence_status = "5% ASV-prevalence filtered",
    physeq_path =
      file.path(work_dir, "phyloseq_final_March26_5pct_without7699.rds")
  ),
  list(
    id = "02_SENSITIVITY_prevalence_unfiltered_without7699",
    label = "Sensitivity analysis without the 5% ASV-prevalence filter",
    prevalence_status = "No ASV-prevalence filter after selective rarefaction",
    physeq_path =
      file.path(work_dir, "phyloseq_March26_post_rarefaction_without7699.rds")
  )
)

GROUP_VAR <- "Sample_Type"
PATIENT_VAR <- "Patient"

STRICT_EXPECTED_COUNTS <- TRUE

expected_sample_counts <- c(
  SSL_HP = 38L,
  H_ADJ_SSL_HP = 12L,
  H_COL = 14L
)

expected_patient_counts <- c(
  SSL_HP = 23L,
  H_ADJ_SSL_HP = 10L,
  H_COL = 8L
)

PLOT_GROUPS <- c("SSL_HP", "H_ADJ_SSL_HP", "H_COL")

group_labels <- c(
  SSL_HP = "SSL_HP",
  H_ADJ_SSL_HP = "H_ADJ_SSL_HP",
  H_COL = "H_COL"
)

group_cols <- c(
  SSL_HP = "#FF0000",
  H_ADJ_SSL_HP = "#E9C61D",
  H_COL = "#1313EC"
)

target_genera <- c(
  "Escherichia-Shigella",
  "Serratia",
  "Fusobacterium"
)

CLR_PSEUDOCOUNT <- 0.5

comparison_df <- tibble(
  Contrast = c(
    "SSL_HP_vs_H_ADJ_SSL_HP",
    "SSL_HP_vs_H_COL",
    "H_ADJ_SSL_HP_vs_H_COL"
  ),
  Group1 = c("SSL_HP", "SSL_HP", "H_ADJ_SSL_HP"),
  Group2 = c("H_ADJ_SSL_HP", "H_COL", "H_COL"),
  Design = c("partially_paired_full_group", "independent_full_group", "independent_full_group"),
  BracketOrder = c(1L, 3L, 2L)
)

ABUNDANCE_PLOT_PVALUE_COLUMN <- "FDR_BH_all"
ABUNDANCE_SIGNIFICANCE_THRESHOLD <- 0.05

font_family <- "Arial"
font_pt <- 7
axis_x_pt <- 6
label_pt <- 6
figure_width_mm <- 183
figure_height_mm <- 148

pt_to_gg <- function(x) x / 2.845276

msg <- function(...) cat(paste0(...), "\n")

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "NULL", "None")] <- NA_character_
  x
}

first_existing <- function(available, candidates) {
  hit <- candidates[candidates %in% available]
  if (length(hit) == 0L) NULL else hit[[1L]]
}

normalise_group <- function(x) {
  x <- trimws(as.character(x))
  key <- toupper(gsub("[^A-Za-z0-9]", "", x))

  dplyr::case_when(
    key %in% c(
      "HCOL", "HCOLD", "HCOLND", "HEALTHYCOL", "HEALTHYCOLON"
    ) ~ "H_COL",

    key %in% c(
      "HADJSSLHP", "HADJSSL", "HADJHP",
      "HADJACENTSSLHP", "HADJACENTSSL", "HADJACENTHP"
    ) ~ "H_ADJ_SSL_HP",

    key %in% c(
      "SSLHP", "SSL", "HP", "LESIONSSLHP", "LESIONSSL", "LESIONHP"
    ) ~ "SSL_HP",

    TRUE ~ NA_character_
  )
}

clean_tax_name <- function(x) {
  x <- clean_character(x)
  x <- gsub("^g__", "", x, ignore.case = TRUE)
  x <- gsub("^Genus:", "", x, ignore.case = TRUE)
  x <- trimws(x)

  unresolved <- is.na(x) |
    toupper(x) %in% c(
      "", "NA", "N/A", "UNCLASSIFIED", "UNCULTURED",
      "UNKNOWN", "UNASSIGNED"
    )

  x[unresolved] <- NA_character_
  x
}

taxon_key <- function(x) {
  x <- clean_tax_name(x)
  tolower(gsub("[^A-Za-z0-9]", "", x))
}

canonicalise_genus <- function(x) {
  cleaned <- clean_tax_name(x)
  key <- taxon_key(cleaned)

  cleaned[!is.na(key) & key == "escherichiashigella"] <-
    "Escherichia-Shigella"
  cleaned[!is.na(key) & key == "serratia"] <- "Serratia"
  cleaned[!is.na(key) & key == "fusobacterium"] <- "Fusobacterium"

  cleaned
}

p_to_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ ""
  )
}

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

check_pvalue_column <- function(column_name) {
  permitted <- c(
    "P_raw",
    "FDR_BH_all",
    "Holm_all",
    "FDR_BH_within_genus",
    "Holm_within_genus"
  )

  if (!(column_name %in% permitted)) {
    stopf(
      "P-value column '%s' is invalid. Choose one of: %s",
      column_name,
      paste(permitted, collapse = ", ")
    )
  }
}

check_pvalue_column(ABUNDANCE_PLOT_PVALUE_COLUMN)

capture_model <- function(expr) {
  warnings <- character()
  error <- NULL
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) { error <<- conditionMessage(e); NULL })
  list(value = value, warnings = unique(warnings), error = error)
}

full_group_model <- function(data, response, out_dir, genus_name) {
  d <- as.data.frame(data[data$Group %in% c("H_ADJ_SSL_HP", "SSL_HP"), ])
  d$Group <- factor(as.character(d$Group), levels = c("H_ADJ_SSL_HP", "SSL_HP"))
  d$Patient <- factor(d$Patient)
  d$Response <- d[[response]]
  stopifnot(nrow(d) == 33L, nlevels(d$Patient) == 23L,
            !anyDuplicated(paste(d$Patient, d$Group)), all(is.finite(d$Response)))
  model_dir <- file.path(out_dir, "models", gsub("[^A-Za-z0-9_-]", "_", genus_name), response)
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(d, file.path(model_dir, "input.csv"), row.names = FALSE)
  binary <- identical(response, "Detected")
  attempt <- capture_model({
    if (length(unique(d$Response)) < 2L) stop("The response is constant; a mixed model cannot be estimated.")
    if (binary && any(vapply(split(d$Response, d$Group), function(z) length(unique(z)) < 2L, logical(1)))) stop("A tissue group has only one detection outcome; the unpenalised logistic group effect is separated and is not reported.")
    if (binary) {
      lme4::glmer(Response ~ Group + (1 | Patient), data = d,
        family = stats::binomial(), nAGQ = 1L,
        control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
    } else {
      lmerTest::lmer(Response ~ Group + (1 | Patient), data = d, REML = TRUE,
        control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
    }
  })
  method <- if (binary) "Logistic mixed model; patient random intercept; Wald z test" else
    "Linear mixed model; patient random intercept; Satterthwaite t test"
  answer <- list(P = NA_real_, W = NA_real_, RankBiserial = NA_real_,
    Estimate = NA_real_, SE = NA_real_, Df = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
    OddsRatio = NA_real_, Status = "FAILED", Test = method)
  messages <- c(attempt$warnings, attempt$error)
  if (binary) {
    paired_detection <- merge(
      d[d$Group == "H_ADJ_SSL_HP", c("Patient", "Response")],
      d[d$Group == "SSL_HP", c("Patient", "Response")],
      by = "Patient", suffixes = c("_adjacent", "_lesion"))
    lesion_only <- sum(paired_detection$Response_adjacent == 0 & paired_detection$Response_lesion == 1)
    adjacent_only <- sum(paired_detection$Response_adjacent == 1 & paired_detection$Response_lesion == 0)
    if (lesion_only + adjacent_only > 0L && min(lesion_only, adjacent_only) == 0L) {
      messages <- c(messages, "Detection discordance among shared patients occurs in only one direction; review the likelihood boundary before reporting unpenalised logistic Wald inference.")
    }
  }
  if (!is.null(attempt$value)) {
    fit <- attempt$value
    saveRDS(fit, file.path(model_dir, "model.rds"))
    ct <- coef(summary(fit))
    write.csv(data.frame(Term = rownames(ct), ct, check.names = FALSE),
              file.path(model_dir, "coefficients.csv"), row.names = FALSE)
    if (!"GroupSSL_HP" %in% rownames(ct)) stop("The mixed-model tissue coefficient is missing.")
    coefficient <- ct["GroupSSL_HP", ]
    answer$Estimate <- unname(coefficient["Estimate"])
    answer$SE <- unname(coefficient["Std. Error"])
    answer$Df <- if (binary) Inf else unname(coefficient["df"])
    answer$W <- unname(coefficient[if (binary) "z value" else "t value"])
    answer$P <- unname(coefficient[if (binary) "Pr(>|z|)" else "Pr(>|t|)"])
    critical <- if (binary) qnorm(0.975) else qt(0.975, df = answer$Df)
    answer$CI_low <- answer$Estimate - critical * answer$SE
    answer$CI_high <- answer$Estimate + critical * answer$SE
    if (binary) {
      answer$OddsRatio <- exp(answer$Estimate)
      answer$CI_low <- exp(answer$CI_low)
      answer$CI_high <- exp(answer$CI_high)
    }
    singular <- lme4::isSingular(fit)
    convergence <- fit@optinfo$conv$lme4$messages
    optimizer_code <- fit@optinfo$conv$opt
    messages <- unique(c(messages, convergence,
      if (singular) "Random-intercept variance is on the boundary; review the fit."))
    finite <- all(is.finite(c(answer$Estimate, answer$SE, answer$W, answer$P,
                             answer$CI_low, answer$CI_high))) && answer$SE > 0 &&
              answer$P >= 0 && answer$P <= 1
    answer$Status <- if (!finite || singular || length(messages) ||
                         any(optimizer_code != 0)) "REVIEW_REQUIRED" else "COMPLETED"
    diagnostics <- data.frame(ProfileID = d$ProfileID, Group = d$Group, Patient = d$Patient,
      Response = d$Response, Fitted = fitted(fit), Residual = residuals(fit))
    write.csv(diagnostics, file.path(model_dir, "diagnostics.csv"), row.names = FALSE)
    if (!binary) {
      grDevices::pdf(file.path(model_dir, "residual_diagnostics.pdf"), width = 8, height = 4)
      old_par <- par(mfrow = c(1, 2))
      plot(diagnostics$Fitted, diagnostics$Residual, xlab = "Fitted mean CLR", ylab = "Residual")
      abline(h = 0, lty = 2)
      qqnorm(diagnostics$Residual); qqline(diagnostics$Residual)
      par(old_par)
      grDevices::dev.off()
    }
  }
  writeLines(if (length(messages)) messages else "None", file.path(model_dir, "warnings.txt"))
  write.csv(as.data.frame(answer), file.path(model_dir, "status.csv"), row.names = FALSE)
  answer
}

safe_independent_wilcox <- function(x, y) {
  x <- as.numeric(x[is.finite(x)])
  y <- as.numeric(y[is.finite(y)])

  if (length(x) < 2L || length(y) < 2L) {
    return(list(P = NA_real_, W = NA_real_, RankBiserial = NA_real_))
  }

  if (length(unique(c(x, y))) == 1L) {
    return(list(
      P = 1,
      W = length(x) * length(y) / 2,
      RankBiserial = 0
    ))
  }

  fit <- suppressWarnings(
    tryCatch(
      wilcox.test(
        x,
        y,
        paired = FALSE,
        exact = FALSE,
        correct = FALSE
      ),
      error = function(e) NULL
    )
  )

  if (is.null(fit)) {
    return(list(P = NA_real_, W = NA_real_, RankBiserial = NA_real_))
  }

  w <- unname(fit$statistic)
  rank_biserial <- 2 * w / (length(x) * length(y)) - 1

  list(
    P = as.numeric(fit$p.value),
    W = as.numeric(w),
    RankBiserial = as.numeric(rank_biserial)
  )
}

safe_fisher <- function(detected1, detected2) {
  detected1 <- as.integer(detected1)
  detected2 <- as.integer(detected2)

  table_2x2 <- matrix(
    c(
      sum(detected1 == 1L),
      sum(detected1 == 0L),
      sum(detected2 == 1L),
      sum(detected2 == 0L)
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(
      Group = c("Group1", "Group2"),
      Detection = c("Detected", "NotDetected")
    )
  )

  fit <- tryCatch(
    fisher.test(table_2x2, alternative = "two.sided"),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(list(
      P = NA_real_,
      OddsRatio = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      Table = table_2x2
    ))
  }

  odds_ratio <- if (length(fit$estimate) == 0L) NA_real_ else unname(fit$estimate)
  conf_int <- if (is.null(fit$conf.int)) c(NA_real_, NA_real_) else fit$conf.int

  list(
    P = as.numeric(fit$p.value),
    OddsRatio = as.numeric(odds_ratio),
    CI_low = as.numeric(conf_int[1L]),
    CI_high = as.numeric(conf_int[2L]),
    Table = table_2x2
  )
}

add_multiple_testing_columns <- function(df) {
  n_tests <- nrow(df)
  df %>%
    mutate(Family_complete = all(Status == "COMPLETED"),
      P_for_adjustment = ifelse(Status == "COMPLETED", P_raw, NA_real_),
      FDR_BH_all = p.adjust(P_for_adjustment, method = "BH", n = n_tests),
      Holm_all = p.adjust(P_for_adjustment, method = "holm", n = n_tests)) %>%
    group_by(Genus) %>%
    mutate(FDR_BH_within_genus = p.adjust(P_for_adjustment, method = "BH", n = 3L),
           Holm_within_genus = p.adjust(P_for_adjustment, method = "holm", n = 3L)) %>%
    ungroup()
}

theme_candidate_paper <- function() {
  theme_classic(base_family = font_family, base_size = font_pt) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(
        family = font_family,
        size = font_pt,
        face = "bold",
        colour = "black",
        hjust = 0
      ),
      axis.title.x = element_text(
        family = font_family,
        size = font_pt,
        colour = "black"
      ),
      axis.title.y = element_text(
        family = font_family,
        size = font_pt,
        colour = "black",
        margin = margin(r = 4)
      ),
      axis.text.x = element_text(
        family = font_family,
        size = axis_x_pt,
        colour = "black",
        angle = 35,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(
        family = font_family,
        size = font_pt,
        colour = "black"
      ),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = unit(1.2, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(
        family = font_family,
        size = font_pt,
        face = "italic",
        colour = "black"
      ),
      legend.position = "none",
      plot.margin = margin(3, 3, 3, 3, "mm")
    )
}

save_figure_formats <- function(plot, stem, width_mm, height_mm) {
  pdf_device <- if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }

  ggsave(
    paste0(stem, ".pdf"),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = pdf_device,
    bg = "white"
  )

  bitmap_device <- if (
    requireNamespace("ragg", quietly = TRUE)
  ) {
    list(
      png = ragg::agg_png,
      tiff = ragg::agg_tiff,
      jpg = ragg::agg_jpeg
    )
  } else {
    NULL
  }

  if (!is.null(bitmap_device)) {
    ggsave(
      paste0(stem, ".png"),
      plot,
      width = width_mm,
      height = height_mm,
      units = "mm",
      dpi = 600,
      device = bitmap_device$png,
      background = "white"
    )

    ggsave(
      paste0(stem, ".tiff"),
      plot,
      width = width_mm,
      height = height_mm,
      units = "mm",
      dpi = 600,
      device = bitmap_device$tiff,
      compression = "lzw",
      background = "white"
    )

    ggsave(
      paste0(stem, ".jpg"),
      plot,
      width = width_mm,
      height = height_mm,
      units = "mm",
      dpi = 600,
      device = bitmap_device$jpg,
      quality = 100,
      background = "white"
    )
  } else {
    if (isTRUE(capabilities("cairo"))) options(bitmapType = "cairo")

    ggsave(
      paste0(stem, ".png"),
      plot,
      width = width_mm,
      height = height_mm,
      units = "mm",
      dpi = 600,
      bg = "white"
    )

    try(
      ggsave(
        paste0(stem, ".tiff"),
        plot,
        width = width_mm,
        height = height_mm,
        units = "mm",
        dpi = 600,
        device = "tiff",
        compression = "lzw",
        bg = "white"
      ),
      silent = TRUE
    )

    ggsave(
      paste0(stem, ".jpg"),
      plot,
      width = width_mm,
      height = height_mm,
      units = "mm",
      dpi = 600,
      bg = "white"
    )
  }
}

run_candidate_analysis <- function(spec) {
  analysis_id <- spec$id
  analysis_label <- spec$label
  prevalence_status <- spec$prevalence_status
  physeq_path <- spec$physeq_path
  out_dir <- file.path(out_root, analysis_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  msg("\n============================================================")
  msg("RUNNING: ", analysis_label)
  msg("INPUT: ", physeq_path)
  msg("OUTPUT: ", out_dir)
  msg("============================================================")

  msg("Reading phyloseq object: ", physeq_path)

  if (!file.exists(physeq_path)) {
    stopf("Input phyloseq file does not exist: %s", physeq_path)
  }

  phy0 <- readRDS(physeq_path)

  if (!inherits(phy0, "phyloseq")) {
    stop("The loaded RDS object is not a phyloseq object.", call. = FALSE)
  }

  otu0 <- as(otu_table(phy0), "matrix")
  if (taxa_are_rows(phy0)) {
    sample_reads <- colSums(otu0)
  } else {
    sample_reads <- rowSums(otu0)
  }

  phy0 <- prune_samples(sample_reads > 0, phy0)
  phy0 <- prune_taxa(taxa_sums(phy0) > 0, phy0)

  metadata <- data.frame(
    sample_data(phy0),
    check.names = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    rownames_to_column("SampleID")

  if (!(GROUP_VAR %in% colnames(metadata))) {
    stopf(
      "Group column '%s' was not found. Available columns: %s",
      GROUP_VAR,
      paste(colnames(metadata), collapse = ", ")
    )
  }

  metadata$Group <- normalise_group(metadata[[GROUP_VAR]])

  resolved_patient_var <- if (PATIENT_VAR %in% colnames(metadata)) {
    PATIENT_VAR
  } else {
    first_existing(
      colnames(metadata),
      c("Patient", "Patient_ID", "PatientID", "patient", "patient_id")
    )
  }

  patient_from_id <- sub(".*[_-]([0-9]+)$", "\\1", metadata$SampleID)
  patient_from_id[patient_from_id == metadata$SampleID] <- NA_character_

  if (!is.null(resolved_patient_var)) {
    metadata$Patient <- clean_character(metadata[[resolved_patient_var]])

    comparable <- !is.na(metadata$Patient) & !is.na(patient_from_id)
    mismatch <- comparable & metadata$Patient != patient_from_id

    if (any(mismatch)) {
      stopf(
        "Patient metadata disagrees with SampleID for: %s",
        paste(metadata$SampleID[mismatch], collapse = ", ")
      )
    }

    metadata$Patient[is.na(metadata$Patient)] <-
      patient_from_id[is.na(metadata$Patient)]
    patient_source <- paste0("validated metadata column: ", resolved_patient_var)
  } else {
    metadata$Patient <- patient_from_id
    patient_source <- "terminal numeric token parsed from SampleID"
  }

  metadata <- metadata %>%
    filter(!is.na(Group), Group %in% PLOT_GROUPS)

  if (nrow(metadata) == 0L || anyNA(metadata$Patient)) {
    stop("No analysable samples remain, or at least one patient ID is unresolved.", call. = FALSE)
  }

  if (
    any(grepl("(^|[_-])7699$", metadata$SampleID)) ||
    "7699" %in% metadata$Patient
  ) {
    stop(
      "Sample/patient 7699 is present, but this script requires the without7699 object.",
      call. = FALSE
    )
  }

  phy <- prune_samples(metadata$SampleID, phy0)
  phy <- prune_taxa(taxa_sums(phy) > 0, phy)

  metadata <- metadata[
    match(sample_names(phy), metadata$SampleID),
    ,
    drop = FALSE
  ]

  if (!identical(metadata$SampleID, sample_names(phy))) {
    stop("Metadata and phyloseq sample order failed alignment.", call. = FALSE)
  }

  metadata$Group <- factor(metadata$Group, levels = PLOT_GROUPS)

  observed_sample_counts <- table(metadata$Group)
  observed_patient_counts <- metadata %>%
    distinct(Group, Patient) %>%
    count(Group, name = "N")
  observed_patient_counts <- setNames(
    observed_patient_counts$N,
    as.character(observed_patient_counts$Group)
  )

  msg("\nBiopsy counts:")
  print(observed_sample_counts)
  msg("\nPatient x group counts:")
  print(observed_patient_counts)

  if (STRICT_EXPECTED_COUNTS) {
    sample_check <- as.integer(observed_sample_counts[names(expected_sample_counts)])
    patient_check <- as.integer(observed_patient_counts[names(expected_patient_counts)])

    if (anyNA(sample_check) ||
        !identical(sample_check, as.integer(expected_sample_counts))) {
      stopf(
        "Unexpected sample counts. Expected %s; observed %s.",
        paste(
          names(expected_sample_counts),
          expected_sample_counts,
          sep = "=",
          collapse = ", "
        ),
        paste(
          names(observed_sample_counts),
          observed_sample_counts,
          sep = "=",
          collapse = ", "
        )
      )
    }

    if (anyNA(patient_check) ||
        !identical(patient_check, as.integer(expected_patient_counts))) {
      stopf(
        "Unexpected patient x group counts. Expected %s; observed %s.",
        paste(
          names(expected_patient_counts),
          expected_patient_counts,
          sep = "=",
          collapse = ", "
        ),
        paste(
          names(observed_patient_counts),
          observed_patient_counts,
          sep = "=",
          collapse = ", "
        )
      )
    }
  }

  cross_class_patients <- metadata %>%
    distinct(Patient, Group) %>%
    group_by(Patient) %>%
    summarise(
      Has_H_COL = any(Group == "H_COL"),
      Has_SSL_or_adjacent = any(Group %in% c("SSL_HP", "H_ADJ_SSL_HP")),
      .groups = "drop"
    ) %>%
    filter(Has_H_COL & Has_SSL_or_adjacent)

  if (nrow(cross_class_patients) > 0L) {
    stopf(
      "H_COL is not independent: these patients also occur in lesion/adjacent groups: %s",
      paste(cross_class_patients$Patient, collapse = ", ")
    )
  }

  lesion_adjacent_overlap <- metadata %>%
    filter(Group %in% c("SSL_HP", "H_ADJ_SSL_HP")) %>%
    distinct(Patient, Group) %>%
    count(Patient, name = "N_groups") %>%
    filter(N_groups == 2L) %>%
    pull(Patient)

  msg(
    "\nPatients represented in both SSL_HP and H_ADJ_SSL_HP: ",
    length(lesion_adjacent_overlap),
    " (AUDIT ONLY; no full-group comparison is subset)"
  )

  write.csv(
    tibble(
      Patient = as.character(lesion_adjacent_overlap),
      Purpose = paste0(
        "Overlap audit only. All 23 SSL_HP and all 10 H_ADJ_SSL_HP patient-level ",
        "profiles remain in the full-group comparison."
      )
    ),
    file.path(out_dir, "00_SSL_HP_H_ADJ_SSL_HP_overlap_AUDIT_ONLY.csv"),
    row.names = FALSE
  )

  metadata <- metadata %>%
    mutate(
      ProfileID = paste(Patient, as.character(Group), sep = "__")
    )

  write.csv(
    metadata,
    file.path(out_dir, "00_sample_manifest.csv"),
    row.names = FALSE
  )

  profile_manifest <- metadata %>%
    group_by(ProfileID, Patient, Group) %>%
    summarise(
      N_biopsies = n(),
      SampleIDs = paste(SampleID, collapse = ";"),
      .groups = "drop"
    ) %>%
    arrange(Group, Patient)

  write.csv(
    profile_manifest,
    file.path(out_dir, "00_patient_group_profile_manifest.csv"),
    row.names = FALSE
  )

  asv_counts <- as(otu_table(phy), "matrix")
  if (taxa_are_rows(phy)) asv_counts <- t(asv_counts)
  storage.mode(asv_counts) <- "double"

  if (any(!is.finite(asv_counts)) || any(asv_counts < 0)) {
    stop("The ASV count table contains invalid values.", call. = FALSE)
  }

  if (max(abs(asv_counts - round(asv_counts))) > 1e-7) {
    stop("The phyloseq abundance table is not an integer count table.", call. = FALSE)
  }

  asv_counts <- asv_counts[metadata$SampleID, , drop = FALSE]

  if (!identical(rownames(asv_counts), metadata$SampleID)) {
    stop("ASV counts and metadata failed row alignment.", call. = FALSE)
  }

  taxonomy <- data.frame(
    tax_table(phy),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  taxonomy <- taxonomy[colnames(asv_counts), , drop = FALSE]

  if (!identical(rownames(taxonomy), colnames(asv_counts))) {
    stop("Taxonomy and ASV counts failed alignment.", call. = FALSE)
  }

  possible_genus_cols <- c("Genus", "genus", "GENUS", "Rank6", "rank6", "g")
  GENUS_COL <- first_existing(colnames(taxonomy), possible_genus_cols)

  if (is.null(GENUS_COL)) {
    stopf(
      "No genus column was found. Available taxonomy columns: %s",
      paste(colnames(taxonomy), collapse = ", ")
    )
  }

  msg("\nUsing taxonomy genus column: ", GENUS_COL)

  genus <- canonicalise_genus(taxonomy[[GENUS_COL]])
  taxonomy_text <- apply(taxonomy, 1L, function(z) paste(z, collapse = ";"))

  keep_asv <- !is.na(genus) &
    !grepl("CHLOROPLAST|MITOCHONDR", toupper(taxonomy_text))

  if (sum(keep_asv) < 2L) {
    stop("Fewer than two classified bacterial genera remain.", call. = FALSE)
  }

  genus_counts <- t(
    rowsum(
      t(asv_counts[, keep_asv, drop = FALSE]),
      group = genus[keep_asv],
      reorder = TRUE
    )
  )

  genus_counts <- genus_counts[
    metadata$SampleID,
    colSums(genus_counts) > 0,
    drop = FALSE
  ]

  if (any(rowSums(genus_counts) <= 0)) {
    stop("At least one sample has no retained classified genus counts.", call. = FALSE)
  }

  missing_target_genera <- setdiff(target_genera, colnames(genus_counts))

  if (length(missing_target_genera) > 0L) {
    stopf(
      paste0(
        "The following target genera were not found after taxonomy cleaning: %s. ",
        "Inspect the taxonomy labels rather than silently adding zero-only rows."
      ),
      paste(missing_target_genera, collapse = ", ")
    )
  }

  genus_input_summary <- tibble(
    Genus = colnames(genus_counts),
    Total_reads = colSums(genus_counts),
    Biopsy_prevalence = colMeans(genus_counts > 0)
  ) %>%
    arrange(desc(Total_reads))

  write.csv(
    genus_input_summary,
    file.path(out_dir, "00_genus_input_summary.csv"),
    row.names = FALSE
  )

  write.csv(
    genus_input_summary %>% filter(Genus %in% target_genera),
    file.path(out_dir, "00_target_genera_input_summary.csv"),
    row.names = FALSE
  )

  log_counts <- log(genus_counts + CLR_PSEUDOCOUNT)
  clr_all <- sweep(log_counts, 1L, rowMeans(log_counts), FUN = "-")

  relative_abundance_all <- sweep(
    genus_counts,
    1L,
    rowSums(genus_counts),
    FUN = "/"
  )

  target_clr <- clr_all[, target_genera, drop = FALSE]
  target_ra <- relative_abundance_all[, target_genera, drop = FALSE]
  target_reads <- genus_counts[, target_genera, drop = FALSE]

  biopsy_long <- tibble(
    SampleID = rep(rownames(target_clr), times = length(target_genera)),
    Genus = rep(target_genera, each = nrow(target_clr)),
    CLR_abundance = as.vector(target_clr),
    Relative_abundance = as.vector(target_ra),
    Relative_abundance_percent = 100 * as.vector(target_ra),
    Reads = as.vector(target_reads)
  ) %>%
    left_join(
      metadata %>% select(SampleID, ProfileID, Patient, Group),
      by = "SampleID"
    ) %>%
    mutate(
      Detected = as.integer(Reads > 0),
      Genus = factor(Genus, levels = target_genera),
      Group = factor(Group, levels = PLOT_GROUPS)
    )

  if (anyNA(biopsy_long$Patient) || anyNA(biopsy_long$Group)) {
    stop("Biopsy-level target-genus data failed metadata alignment.", call. = FALSE)
  }

  write.csv(
    biopsy_long,
    file.path(out_dir, "Target_genera_biopsy_level_long.csv"),
    row.names = FALSE
  )

  patient_profile_long <- biopsy_long %>%
    group_by(ProfileID, Patient, Group, Genus) %>%
    summarise(
      N_biopsies = n(),
      Mean_CLR_abundance = mean(CLR_abundance),
      Median_CLR_abundance = median(CLR_abundance),
      Mean_relative_abundance = mean(Relative_abundance),
      Mean_relative_abundance_percent = mean(Relative_abundance_percent),
      Total_reads_across_biopsies = sum(Reads),
      N_detected_biopsies = sum(Detected),
      Detected_fraction_of_biopsies = mean(Detected),
      Detected = as.integer(any(Detected == 1L)),
      .groups = "drop"
    ) %>%
    mutate(
      Group = factor(Group, levels = PLOT_GROUPS),
      Group_num = as.integer(Group),
      Genus = factor(Genus, levels = target_genera)
    ) %>%
    arrange(Genus, Group, Patient)

  expected_profile_rows <- nrow(profile_manifest) * length(target_genera)
  if (nrow(patient_profile_long) != expected_profile_rows) {
    stop("Patient x group aggregation produced an unexpected number of rows.", call. = FALSE)
  }

  write.csv(
    patient_profile_long,
    file.path(out_dir, "Target_genera_patient_group_profiles_long.csv"),
    row.names = FALSE
  )

  biopsy_detection_summary <- biopsy_long %>%
    group_by(Genus, Group) %>%
    summarise(
      N_biopsies = n(),
      N_patients_represented = n_distinct(Patient),
      N_detected_biopsies = sum(Detected),
      N_undetected_biopsies = sum(Detected == 0L),
      Detection_percent_biopsies = 100 * mean(Detected),
      .groups = "drop"
    ) %>%
    mutate(
      Genus = factor(Genus, levels = target_genera),
      Group = factor(Group, levels = PLOT_GROUPS),
      Group_num = as.integer(Group),
      Detection_label = paste0(N_detected_biopsies, "/", N_biopsies),
      Label_y = Detection_percent_biopsies + 4
    ) %>%
    arrange(Genus, Group)

  write.csv(
    biopsy_detection_summary,
    file.path(out_dir, "Target_genera_biopsy_detection_summary_DESCRIPTIVE_ONLY.csv"),
    row.names = FALSE
  )

  within_patient_biopsy_detection <- biopsy_long %>%
    group_by(Genus, Group, Patient) %>%
    summarise(
      N_biopsies = n(),
      N_detected_biopsies = sum(Detected),
      N_undetected_biopsies = sum(Detected == 0L),
      Detection_pattern = case_when(
        N_detected_biopsies == 0L ~ "No biopsies detected",
        N_undetected_biopsies == 0L ~ "All biopsies detected",
        TRUE ~ "Mixed positive and negative biopsies"
      ),
      Detected_sample_IDs = paste(SampleID[Detected == 1L], collapse = ";"),
      Undetected_sample_IDs = paste(SampleID[Detected == 0L], collapse = ";"),
      .groups = "drop"
    ) %>%
    mutate(
      Genus = factor(Genus, levels = target_genera),
      Group = factor(Group, levels = PLOT_GROUPS)
    ) %>%
    arrange(Genus, Group, Patient)

  write.csv(
    within_patient_biopsy_detection,
    file.path(
      out_dir,
      "Target_genera_within_patient_biopsy_detection_AUDIT.csv"
    ),
    row.names = FALSE
  )

  within_patient_mixed_detection <- within_patient_biopsy_detection %>%
    filter(Detection_pattern == "Mixed positive and negative biopsies")

  write.csv(
    within_patient_mixed_detection,
    file.path(
      out_dir,
      "Target_genera_within_patient_MIXED_DETECTION_ONLY.csv"
    ),
    row.names = FALSE
  )

  within_patient_discordance_summary <- within_patient_biopsy_detection %>%
    group_by(Genus, Group) %>%
    summarise(
      N_patient_group_profiles = n(),
      N_profiles_with_multiple_biopsies = sum(N_biopsies > 1L),
      N_profiles_with_mixed_detection = sum(
        Detection_pattern == "Mixed positive and negative biopsies"
      ),
      .groups = "drop"
    ) %>%
    mutate(
      Genus = factor(Genus, levels = target_genera),
      Group = factor(Group, levels = PLOT_GROUPS)
    ) %>%
    arrange(Genus, Group)

  write.csv(
    within_patient_discordance_summary,
    file.path(
      out_dir,
      "Target_genera_within_patient_mixed_detection_SUMMARY.csv"
    ),
    row.names = FALSE
  )

  msg("\nWithin-patient mixed-detection summary:")
  print(within_patient_discordance_summary)
  msg("\nExact mixed-detection Patient x tissue-group profiles:")
  print(within_patient_mixed_detection, n = Inf)

  abundance_stats_list <- list()
  row_counter <- 0L

  for (genus_name in target_genera) {
    genus_data <- patient_profile_long %>%
      filter(as.character(Genus) == genus_name)

    observed_n <- genus_data %>%
      distinct(Patient, Group) %>%
      count(Group, name = "N")

    for (i in seq_len(nrow(comparison_df))) {
      row_counter <- row_counter + 1L

      contrast <- comparison_df$Contrast[i]
      group1 <- comparison_df$Group1[i]
      group2 <- comparison_df$Group2[i]
      design <- comparison_df$Design[i]

      x <- genus_data %>%
        filter(Group == group1) %>%
        pull(Mean_CLR_abundance)

      y <- genus_data %>%
        filter(Group == group2) %>%
        pull(Mean_CLR_abundance)

      expected_n1 <- unname(expected_patient_counts[group1])
      expected_n2 <- unname(expected_patient_counts[group2])

      if (STRICT_EXPECTED_COUNTS &&
          (length(x) != expected_n1 || length(y) != expected_n2)) {
        stopf(
          paste0(
            "Full-group abundance comparison %s for %s has unexpected sizes. ",
            "Expected %s=%d and %s=%d; observed %d and %d."
          ),
          contrast,
          genus_name,
          group1,
          expected_n1,
          group2,
          expected_n2,
          length(x),
          length(y)
        )
      }

      shared <- intersect(as.character(genus_data$Patient[genus_data$Group == group1]),
                          as.character(genus_data$Patient[genus_data$Group == group2]))
      if (contrast == "SSL_HP_vs_H_ADJ_SSL_HP") {
        if (length(shared) != 10L) stop("Expected 10 shared patients.")
        test_result <- full_group_model(genus_data, "Mean_CLR_abundance", out_dir, genus_name)
      } else {
        if (length(shared)) stop("An independent comparison contains shared patients.")
        test_result <- safe_independent_wilcox(x, y)
        test_result$Status <- if (is.finite(test_result$P)) "COMPLETED" else "FAILED"
        test_result$Test <- "Wilcoxon rank-sum test using all patient-level profiles in both groups"
        test_result$Estimate <- test_result$SE <- test_result$Df <- NA_real_
      }

      abundance_stats_list[[row_counter]] <- tibble(
        Genus = genus_name,
        Contrast = contrast,
        Group1 = group1,
        Group2 = group2,
        Design = design,
        Test = test_result$Test,
        Status = test_result$Status,
        N_shared_patients = length(shared),
        Model_estimate = test_result$Estimate,
        Model_SE = test_result$SE,
        Model_df = test_result$Df,
        N_group1 = length(x),
        N_group2 = length(y),
        N_pairs = NA_integer_,
        Mean_group1 = mean(x),
        Mean_group2 = mean(y),
        Median_group1 = median(x),
        Median_group2 = median(y),
        Median_difference_group1_minus_group2 = median(x) - median(y),
        Mean_difference_group1_minus_group2 = mean(x) - mean(y),
        Rank_biserial_correlation = test_result$RankBiserial,
        Statistic = test_result$W,
        P_raw = test_result$P
      )
    }
  }

  abundance_stats <- bind_rows(abundance_stats_list) %>%
    add_multiple_testing_columns() %>%
    left_join(
      comparison_df %>% select(Contrast, BracketOrder),
      by = "Contrast"
    ) %>%
    mutate(
      Genus = factor(Genus, levels = target_genera)
    ) %>%
    arrange(Genus, BracketOrder)

  if (nrow(abundance_stats) != length(target_genera) * 3L) {
    stop("Expected exactly nine abundance tests: 3 genera x 3 comparisons.", call. = FALSE)
  }

  write.csv(
    abundance_stats,
    file.path(
      out_dir,
      "Target_genera_patient_level_CLR_abundance_statistics_ALL3_FULL_GROUP_COMPARISONS.csv"
    ),
    row.names = FALSE
  )

  detection_stats_list <- list()
  row_counter <- 0L

  for (genus_name in target_genera) {
    genus_data <- patient_profile_long %>%
      filter(as.character(Genus) == genus_name)

    for (i in seq_len(nrow(comparison_df))) {
      row_counter <- row_counter + 1L

      contrast <- comparison_df$Contrast[i]
      group1 <- comparison_df$Group1[i]
      group2 <- comparison_df$Group2[i]
      design <- comparison_df$Design[i]

      x <- genus_data %>%
        filter(Group == group1) %>%
        pull(Detected)

      y <- genus_data %>%
        filter(Group == group2) %>%
        pull(Detected)

      expected_n1 <- unname(expected_patient_counts[group1])
      expected_n2 <- unname(expected_patient_counts[group2])

      if (STRICT_EXPECTED_COUNTS &&
          (length(x) != expected_n1 || length(y) != expected_n2)) {
        stopf(
          paste0(
            "Full-group detection comparison %s for %s has unexpected sizes. ",
            "Expected %s=%d and %s=%d; observed %d and %d."
          ),
          contrast,
          genus_name,
          group1,
          expected_n1,
          group2,
          expected_n2,
          length(x),
          length(y)
        )
      }

      shared <- intersect(as.character(genus_data$Patient[genus_data$Group == group1]),
                          as.character(genus_data$Patient[genus_data$Group == group2]))
      if (contrast == "SSL_HP_vs_H_ADJ_SSL_HP") {
        if (length(shared) != 10L) stop("Expected 10 shared patients.")
        test_result <- full_group_model(genus_data, "Detected", out_dir, genus_name)
      } else {
        if (length(shared)) stop("An independent comparison contains shared patients.")
        test_result <- safe_fisher(x, y)
        test_result$Status <- if (is.finite(test_result$P)) "COMPLETED" else "FAILED"
        test_result$Test <- "Fisher exact test using all patient-level profiles in both groups"
        test_result$Estimate <- test_result$SE <- test_result$Df <- NA_real_
      }

      detection_stats_list[[row_counter]] <- tibble(
        Genus = genus_name,
        Contrast = contrast,
        Group1 = group1,
        Group2 = group2,
        Design = design,
        Test = test_result$Test,
        Status = test_result$Status,
        N_shared_patients = length(shared),
        N_group1 = length(x),
        N_group2 = length(y),
        N_pairs = NA_integer_,
        N_detected_group1 = sum(x == 1L),
        N_detected_group2 = sum(y == 1L),
        Detection_percent_group1 = 100 * mean(x == 1L),
        Detection_percent_group2 = 100 * mean(y == 1L),
        Detection_difference_percentage_points =
          100 * (mean(x == 1L) - mean(y == 1L)),
        Group1_only_detected_pairs = NA_integer_,
        Group2_only_detected_pairs = NA_integer_,
        Discordant_pairs = NA_integer_,
        Odds_ratio = test_result$OddsRatio,
        Odds_ratio_CI_low = test_result$CI_low,
        Odds_ratio_CI_high = test_result$CI_high,
        P_raw = test_result$P
      )
    }
  }

  detection_stats <- bind_rows(detection_stats_list) %>%
    add_multiple_testing_columns() %>%
    left_join(
      comparison_df %>% select(Contrast, BracketOrder),
      by = "Contrast"
    ) %>%
    mutate(
      Genus = factor(Genus, levels = target_genera)
    ) %>%
    arrange(Genus, BracketOrder)

  if (nrow(detection_stats) != length(target_genera) * 3L) {
    stop("Expected exactly nine detection tests: 3 genera x 3 comparisons.", call. = FALSE)
  }

  write.csv(
    detection_stats,
    file.path(
      out_dir,
      "Target_genera_patient_level_detection_statistics_ALL3_FULL_GROUP_COMPARISONS.csv"
    ),
    row.names = FALSE
  )

  detection_summary <- patient_profile_long %>%
    group_by(Genus, Group) %>%
    summarise(
      N_patients = n(),
      N_detected_patients = sum(Detected),
      Detection_percent_patients = 100 * mean(Detected),
      .groups = "drop"
    ) %>%
    mutate(
      Genus = factor(Genus, levels = target_genera),
      Group = factor(Group, levels = PLOT_GROUPS),
      Group_num = as.integer(Group),
      Detection_label = paste0(N_detected_patients, "/", N_patients),
      Label_y = Detection_percent_patients + 4
    ) %>%
    arrange(Genus, Group)

  write.csv(
    detection_summary,
    file.path(out_dir, "Target_genera_patient_level_detection_summary.csv"),
    row.names = FALSE
  )

  abundance_y_lookup <- patient_profile_long %>%
    group_by(Genus) %>%
    summarise(
      Y_min = min(Mean_CLR_abundance, na.rm = TRUE),
      Y_max = max(Mean_CLR_abundance, na.rm = TRUE),
      Y_range = Y_max - Y_min,
      .groups = "drop"
    ) %>%
    mutate(
      Y_range = ifelse(!is.finite(Y_range) | Y_range <= 0, 1, Y_range)
    )

  abundance_brackets <- abundance_stats %>%
    mutate(
      P_for_plot = .data[[ABUNDANCE_PLOT_PVALUE_COLUMN]],
      Significance = p_to_stars(P_for_plot),
      X1 = match(Group1, PLOT_GROUPS),
      X2 = match(Group2, PLOT_GROUPS),
      X_mid = (X1 + X2) / 2
    ) %>%
    filter(
      !is.na(P_for_plot),
      P_for_plot <= ABUNDANCE_SIGNIFICANCE_THRESHOLD,
      Status == "COMPLETED", Family_complete,
      Significance != ""
    ) %>%
    group_by(Genus) %>%
    arrange(BracketOrder, .by_group = TRUE) %>%
    mutate(Plot_level = row_number()) %>%
    ungroup() %>%
    left_join(abundance_y_lookup, by = "Genus") %>%
    mutate(
      Y_bracket = Y_max + (0.10 + 0.16 * (Plot_level - 1L)) * Y_range,
      Y_tick = Y_bracket - 0.035 * Y_range,
      Y_text = Y_bracket + 0.035 * Y_range
    )

  write.csv(
    abundance_brackets,
    file.path(out_dir, "Target_genera_abundance_SIGNIFICANT_plot_brackets.csv"),
    row.names = FALSE
  )

  p_abundance <- ggplot(
    patient_profile_long,
    aes(x = Group_num, y = Mean_CLR_abundance)
  ) +
    geom_boxplot(
      aes(group = Group),
      width = 0.44,
      outlier.shape = NA,
      fill = NA,
      colour = "black",
      linewidth = 0.35
    ) +
    geom_point(
      aes(colour = Group),
      position = position_jitter(width = 0.08, height = 0, seed = 20260821),
      size = 1.45,
      alpha = 0.9
    ) +
    geom_segment(
      data = abundance_brackets,
      aes(x = X1, xend = X2, y = Y_bracket, yend = Y_bracket),
      inherit.aes = FALSE,
      linewidth = 0.35,
      colour = "grey35"
    ) +
    geom_segment(
      data = abundance_brackets,
      aes(x = X1, xend = X1, y = Y_bracket, yend = Y_tick),
      inherit.aes = FALSE,
      linewidth = 0.35,
      colour = "grey35"
    ) +
    geom_segment(
      data = abundance_brackets,
      aes(x = X2, xend = X2, y = Y_bracket, yend = Y_tick),
      inherit.aes = FALSE,
      linewidth = 0.35,
      colour = "grey35"
    ) +
    geom_text(
      data = abundance_brackets,
      aes(x = X_mid, y = Y_text, label = Significance),
      inherit.aes = FALSE,
      family = font_family,
      size = pt_to_gg(label_pt),
      colour = "black"
    ) +
    facet_wrap(~ Genus, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = group_cols, drop = FALSE) +
    scale_x_continuous(
      breaks = seq_along(PLOT_GROUPS),
      labels = group_labels[PLOT_GROUPS],
      limits = c(0.5, length(PLOT_GROUPS) + 0.5)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.10))) +
    labs(
      title = "A",
      x = NULL,
      y = "Patient-level mean CLR abundance"
    ) +
    theme_candidate_paper()

  p_detection <- ggplot(
    detection_summary,
    aes(x = Group_num, y = Detection_percent_patients, fill = Group)
  ) +
    geom_col(width = 0.62, colour = NA) +
    geom_text(
      aes(y = Label_y, label = Detection_label),
      family = font_family,
      size = pt_to_gg(label_pt),
      colour = "black"
    ) +
    facet_wrap(~ Genus, nrow = 1) +
    scale_fill_manual(values = group_cols, drop = FALSE) +
    scale_x_continuous(
      breaks = seq_along(PLOT_GROUPS),
      labels = group_labels[PLOT_GROUPS],
      limits = c(0.5, length(PLOT_GROUPS) + 0.5)
    ) +
    scale_y_continuous(
      limits = c(0, 110),
      breaks = c(0, 25, 50, 75, 100),
      labels = function(x) paste0(x, "%"),
      expand = c(0, 0)
    ) +
    labs(
      title = "B",
      x = NULL,
      y = "Patients positive"
    ) +
    theme_candidate_paper() +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank()
    )

  p_biopsy_detection <- ggplot(
    biopsy_detection_summary,
    aes(x = Group_num, y = Detection_percent_biopsies, fill = Group)
  ) +
    geom_col(width = 0.62, colour = NA) +
    geom_text(
      aes(y = Label_y, label = Detection_label),
      family = font_family,
      size = pt_to_gg(label_pt),
      colour = "black"
    ) +
    facet_wrap(~ Genus, nrow = 1) +
    scale_fill_manual(values = group_cols, drop = FALSE) +
    scale_x_continuous(
      breaks = seq_along(PLOT_GROUPS),
      labels = group_labels[PLOT_GROUPS],
      limits = c(0.5, length(PLOT_GROUPS) + 0.5)
    ) +
    scale_y_continuous(
      limits = c(0, 110),
      breaks = c(0, 25, 50, 75, 100),
      labels = function(x) paste0(x, "%"),
      expand = c(0, 0)
    ) +
    labs(
      title = "C",
      x = NULL,
      y = "Biopsies positive"
    ) +
    theme_candidate_paper() +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank()
    )

  save_figure_formats(
    p_abundance,
    file.path(out_dir, "Figure_4_6_panel_A_patient_CLR_abundance_ALL3_FULL_GROUP_COMPARISONS"),
    figure_width_mm,
    62
  )

  save_figure_formats(
    p_detection,
    file.path(out_dir, "Figure_4_6_panel_B_patient_detection_ALL3_FULL_GROUP_COMPARISONS"),
    figure_width_mm,
    44
  )

  save_figure_formats(
    p_biopsy_detection,
    file.path(out_dir, "Figure_4_6_panel_C_biopsy_detection_ALL3_FULL_GROUP_COMPARISONS"),
    figure_width_mm,
    44
  )

  g_abundance <- ggplotGrob(p_abundance)
  g_detection <- ggplotGrob(p_detection)
  g_biopsy_detection <- ggplotGrob(p_biopsy_detection)

  max_widths <- grid::unit.pmax(
    g_abundance$widths,
    g_detection$widths,
    g_biopsy_detection$widths
  )
  g_abundance$widths <- max_widths
  g_detection$widths <- max_widths
  g_biopsy_detection$widths <- max_widths

  combined_grob <- gridExtra::arrangeGrob(
    g_abundance,
    g_detection,
    g_biopsy_detection,
    ncol = 1,
    heights = c(1.08, 0.72, 0.72)
  )

  save_figure_formats(
    combined_grob,
    file.path(out_dir, "Figure_4_6_candidate_genera_ALL3_FULL_GROUP_COMPARISONS_FINAL"),
    figure_width_mm,
    figure_height_mm
  )

  abundance_compact <- abundance_stats %>%
    transmute(
      Analysis = "CLR abundance",
      Genus = as.character(Genus),
      Contrast,
      Design,
      Test,
      Status,
      N_group1,
      N_group2,
      N_pairs,
      Effect = ifelse(Design == "partially_paired_full_group", Model_estimate, Median_difference_group1_minus_group2),
      Effect_definition = ifelse(Design == "partially_paired_full_group", "Mixed-model mean CLR difference: SSL_HP minus H_ADJ_SSL_HP", "Median CLR difference: Group1 minus Group2"),
      P_raw,
      FDR_BH_all,
      Holm_all
    )

  detection_compact <- detection_stats %>%
    transmute(
      Analysis = "Detection",
      Genus = as.character(Genus),
      Contrast,
      Design,
      Test,
      Status,
      N_group1,
      N_group2,
      N_pairs,
      Effect = Detection_difference_percentage_points,
      Effect_definition = "Detection percentage-point difference: Group1 minus Group2",
      P_raw,
      FDR_BH_all,
      Holm_all
    )

  write.csv(
    bind_rows(abundance_compact, detection_compact),
    file.path(out_dir, "Target_genera_compact_statistical_summary.csv"),
    row.names = FALSE
  )

  settings_table <- tibble(
    Setting = c(
      "Analysis ID",
      "Analysis label",
      "Prevalence-filter status",
      "Input phyloseq object",
      "Input MD5",
      "Output directory",
      "Patient identifier source",
      "Target genera",
      "CLR pseudocount",
      "Abundance statistical unit",
      "Patient-detection statistical unit",
      "Biopsy-prevalence display unit",
      "Full patient-level group sizes",
      "Full biopsy-level group sizes",
      "Pairwise comparisons",
      "Abundance test",
      "Patient-detection test",
      "Biopsy-prevalence test",
      "Matched-patient subsetting in this script",
      "Patients represented in both SSL_HP and H_ADJ_SSL_HP",
      "Abundance plotted adjustment",
      "Abundance annotation rule",
      "Patient/biopsy prevalence plot annotations",
      "Figure panels"
    ),
    Value = c(
      analysis_id,
      analysis_label,
      prevalence_status,
      physeq_path,
      unname(tools::md5sum(physeq_path)),
      out_dir,
      patient_source,
      paste(target_genera, collapse = "; "),
      as.character(CLR_PSEUDOCOUNT),
      "Mean biopsy CLR within each Patient x tissue group",
      "Any biopsy positive within each Patient x tissue group",
      "Individual microbiome biopsies",
      "SSL_HP=23; H_ADJ_SSL_HP=10; H_COL=8",
      "SSL_HP=38; H_ADJ_SSL_HP=12; H_COL=14",
      paste(
        "SSL_HP vs H_ADJ_SSL_HP; SSL_HP vs H_COL;",
        "H_ADJ_SSL_HP vs H_COL"
      ),
      "H_COL: independent Wilcoxon; SSL_HP-adjacent: patient random-intercept linear model",
      "H_COL: Fisher exact; SSL_HP-adjacent: patient random-intercept logistic model",
      "None; biopsy prevalence is descriptive because biopsies are clustered within patients",
      "None. The matched analysis is separate.",
      as.character(length(lesion_adjacent_overlap)),
      ABUNDANCE_PLOT_PVALUE_COLUMN,
      paste0(
        "Stars/brackets only when ",
        ABUNDANCE_PLOT_PVALUE_COLUMN,
        " <= ",
        ABUNDANCE_SIGNIFICANCE_THRESHOLD,
        "; no ns annotations"
      ),
      "None",
      paste(
        "A = patient-level mean CLR abundance;",
        "B = patient-level prevalence; C = biopsy-level prevalence"
      )
    )
  )

  write.csv(
    settings_table,
    file.path(out_dir, "Analysis_settings.csv"),
    row.names = FALSE
  )

  readme <- c(
    "FINAL SECTION 4.6 TARGETED CANDIDATE-GENERA ANALYSIS WITH BIOPSY PREVALENCE",
    "",
    paste0("Analysis ID: ", analysis_id),
    paste0("Analysis label: ", analysis_label),
    paste0("Prevalence-filter status: ", prevalence_status),
    paste0("Run date: ", Sys.time()),
    paste0("Input: ", physeq_path),
    paste0("Input MD5: ", unname(tools::md5sum(physeq_path))),
    paste0("Output: ", out_dir),
    paste0("Patient identifier source: ", patient_source),
    "",
    "TARGETS",
    paste(target_genera, collapse = ", "),
    "These are candidate colorectal-neoplasia-associated genera, not organisms",
    "demonstrated to be pathogenic by genus-level 16S sequencing.",
    "",
    "ANALYTICAL UNIT",
    paste0("Biopsy-level CLR pseudocount: ", CLR_PSEUDOCOUNT, "."),
    "CLR was calculated using all retained classified bacterial genera.",
    "Biopsy CLR values were averaged within each Patient x tissue group.",
    "The full patient-level groups are SSL_HP=23, H_ADJ_SSL_HP=10 and H_COL=8.",
    "",
    "FULL-GROUP COMPARISONS",
    "Every genus is tested for all three pairwise group comparisons:",
    "SSL_HP vs H_ADJ_SSL_HP, SSL_HP vs H_COL, and H_ADJ_SSL_HP vs H_COL.",
    "Every comparison uses every available patient-level profile in each group.",
    "No comparison is restricted to the ten matched lesion-adjacent patients.",
    "The matched-patient analysis remains a separate analysis.",
    "",
    "ABUNDANCE",
    "The H_COL comparisons use independent Wilcoxon rank-sum tests.",
    "SSL_HP versus H_ADJ_SSL_HP uses a linear mixed model with a patient random intercept.",
    "",
    "DETECTION AND PREVALENCE",
    "A Patient x tissue-group profile is positive when at least one contributing",
    "biopsy contains the genus. The H_COL comparisons use Fisher's exact tests.",
    "SSL_HP versus H_ADJ_SSL_HP uses a logistic mixed model with a patient random intercept.",
    "Panel C separately shows the proportion of individual biopsies containing",
    "each genus, with detected/total biopsy counts printed above the bars.",
    "Biopsy-level prevalence is descriptive only because multiple biopsies can",
    "come from the same patient.",
    "",
    "WITHIN-PATIENT DETECTION AUDIT",
    "The output also identifies Patient x tissue-group profiles containing both",
    "genus-positive and genus-negative biopsies. This is the direct audit needed",
    "to establish whether within-patient biopsy-level contrasts are available.",
    "",
    "MULTIPLE TESTING",
    "Abundance and detection are corrected as separate families.",
    "Each family contains nine tests: three genera x three group comparisons.",
    "Abundance stars use BH-FDR across all nine abundance tests by default.",
    "Only significant abundance comparisons are annotated; no ns labels appear.",
    "Detection statistics are exported but are not annotated on the prevalence plots.",
    "Raw P values, BH-FDR and Holm values are retained in the CSV outputs.",
    "",
    "FIGURE",
    "Panel A retains all patient-level profiles: SSL_HP=23, H_ADJ_SSL_HP=10 and H_COL=8.",
    "Panel A points represent Patient x tissue-group profiles; boxplots show the",
    "median and interquartile range of patient-level mean CLR abundance.",
    "Panel B shows patient-level prevalence and detected/total patient counts.",
    "Panel C shows biopsy-level prevalence and detected/total biopsy counts:",
    "SSL_HP=38, H_ADJ_SSL_HP=12 and H_COL=14 biopsies.",
    "All three rows use the exact labels SSL_HP, H_ADJ_SSL_HP and H_COL.",
    "",
    "IMPORTANT INTERPRETATION NOTE",
    "Panel C establishes biopsy-level detection frequencies, but it does not by",
    "itself prove that positive and negative biopsies came from the same patient.",
    "Use Target_genera_within_patient_MIXED_DETECTION_ONLY.csv to identify the",
    "exact within-patient discordant profiles and their positive/negative samples.",
    "Ten patients occur in both SSL_HP and H_ADJ_SSL_HP. This script deliberately",
    "implements the requested full-group exploratory comparison without matched",
    "subsetting. The separate matched analysis should be reported alongside the",
    "direct SSL_HP-versus-H_ADJ_SSL_HP result.",
    "",
    "This is a targeted exploratory analysis. Genus-level 16S data do not establish",
    "pathogenic species, strains or virulence traits."
  )

  writeLines(
    readme,
    file.path(out_dir, "READ_ME_FIRST.txt")
  )

  writeLines(
    capture.output(sessionInfo()),
    file.path(out_dir, "sessionInfo.txt")
  )

  msg("\nAnalysis completed: ", analysis_label)
  msg("Output directory: ", out_dir)
  msg(
    "Main figure: ",
    file.path(out_dir, "Figure_4_6_candidate_genera_ALL3_FULL_GROUP_COMPARISONS_FINAL.pdf")
  )

  list(
    analysis_id = analysis_id,
    analysis_label = analysis_label,
    prevalence_status = prevalence_status,
    physeq_path = physeq_path,
    out_dir = out_dir,
    abundance_stats = abundance_stats,
    detection_stats = detection_stats,
    detection_summary = detection_summary,
    biopsy_detection_summary = biopsy_detection_summary,
    within_patient_biopsy_detection = within_patient_biopsy_detection,
    within_patient_mixed_detection = within_patient_mixed_detection,
    within_patient_discordance_summary = within_patient_discordance_summary,
    patient_profiles = patient_profile_long,
    sample_manifest = metadata %>%
      transmute(
        SampleID,
        Patient = as.character(Patient),
        Group = as.character(Group)
      ) %>%
      arrange(SampleID),
    genus_input_summary = genus_input_summary,
    abundance_plot = p_abundance,
    detection_plot = p_detection,
    biopsy_detection_plot = p_biopsy_detection,
    combined_grob = combined_grob
  )
}

results <- lapply(analysis_specs, run_candidate_analysis)
names(results) <- vapply(results, function(x) x$analysis_id, character(1))

primary_manifest <- results[[1]]$sample_manifest
sensitivity_manifest <- results[[2]]$sample_manifest

if (!identical(primary_manifest, sensitivity_manifest)) {
  write.csv(
    primary_manifest,
    file.path(out_root, "AUDIT_primary_sample_manifest.csv"),
    row.names = FALSE
  )
  write.csv(
    sensitivity_manifest,
    file.path(out_root, "AUDIT_sensitivity_sample_manifest.csv"),
    row.names = FALSE
  )
  stop(
    paste0(
      "Primary and sensitivity sample manifests are not identical. ",
      "The two audit manifests were written to the root output directory."
    ),
    call. = FALSE
  )
}

write.csv(
  primary_manifest,
  file.path(out_root, "AUDIT_identical_sample_manifest_primary_and_sensitivity.csv"),
  row.names = FALSE
)

input_feature_audit <- tibble(
  AnalysisID = vapply(results, function(x) x$analysis_id, character(1)),
  AnalysisLabel = vapply(results, function(x) x$analysis_label, character(1)),
  InputPhyloseq = vapply(results, function(x) x$physeq_path, character(1)),
  N_classified_genera_used_for_CLR = vapply(
    results,
    function(x) nrow(x$genus_input_summary),
    integer(1)
  ),
  N_target_genera = length(target_genera),
  N_samples = nrow(primary_manifest),
  Sample_manifest_identical = TRUE
)

write.csv(
  input_feature_audit,
  file.path(out_root, "PRIMARY_vs_UNFILTERED_input_feature_audit.csv"),
  row.names = FALSE
)

all_abundance_stats <- bind_rows(
  lapply(results, function(x) {
    x$abundance_stats %>%
      mutate(
        AnalysisID = x$analysis_id,
        AnalysisLabel = x$analysis_label,
        PrevalenceStatus = x$prevalence_status,
        InputPhyloseq = x$physeq_path,
        .before = 1
      )
  })
)

all_detection_stats <- bind_rows(
  lapply(results, function(x) {
    x$detection_stats %>%
      mutate(
        AnalysisID = x$analysis_id,
        AnalysisLabel = x$analysis_label,
        PrevalenceStatus = x$prevalence_status,
        InputPhyloseq = x$physeq_path,
        .before = 1
      )
  })
)

write.csv(
  all_abundance_stats,
  file.path(out_root, "PRIMARY_vs_UNFILTERED_CLR_abundance_statistics_LONG.csv"),
  row.names = FALSE
)

write.csv(
  all_detection_stats,
  file.path(out_root, "PRIMARY_vs_UNFILTERED_detection_statistics_LONG.csv"),
  row.names = FALSE
)

all_biopsy_detection_summary <- bind_rows(
  lapply(results, function(x) {
    x$biopsy_detection_summary %>%
      mutate(
        AnalysisID = x$analysis_id,
        AnalysisLabel = x$analysis_label,
        PrevalenceStatus = x$prevalence_status,
        InputPhyloseq = x$physeq_path,
        .before = 1
      )
  })
)

write.csv(
  all_biopsy_detection_summary,
  file.path(
    out_root,
    "PRIMARY_vs_UNFILTERED_biopsy_prevalence_DESCRIPTIVE_LONG.csv"
  ),
  row.names = FALSE
)

all_within_patient_discordance <- bind_rows(
  lapply(results, function(x) {
    x$within_patient_discordance_summary %>%
      mutate(
        AnalysisID = x$analysis_id,
        AnalysisLabel = x$analysis_label,
        PrevalenceStatus = x$prevalence_status,
        InputPhyloseq = x$physeq_path,
        .before = 1
      )
  })
)

write.csv(
  all_within_patient_discordance,
  file.path(
    out_root,
    "PRIMARY_vs_UNFILTERED_within_patient_mixed_detection_SUMMARY.csv"
  ),
  row.names = FALSE
)

abundance_comparison_wide <- all_abundance_stats %>%
  select(
    AnalysisID,
    Genus,
    Contrast,
    Design,
    N_group1,
    N_group2,
    N_pairs,
    Median_difference_group1_minus_group2,
    P_raw,
    FDR_BH_all,
    Holm_all
  ) %>%
  mutate(Genus = as.character(Genus)) %>%
  pivot_wider(
    names_from = AnalysisID,
    values_from = c(
      N_group1,
      N_group2,
      N_pairs,
      Median_difference_group1_minus_group2,
      P_raw,
      FDR_BH_all,
      Holm_all
    ),
    names_glue = "{.value}__{AnalysisID}"
  )

write.csv(
  abundance_comparison_wide,
  file.path(out_root, "PRIMARY_vs_UNFILTERED_CLR_abundance_statistics_WIDE.csv"),
  row.names = FALSE
)

detection_comparison_wide <- all_detection_stats %>%
  select(
    AnalysisID,
    Genus,
    Contrast,
    Design,
    N_group1,
    N_group2,
    N_pairs,
    Detection_difference_percentage_points,
    P_raw,
    FDR_BH_all,
    Holm_all
  ) %>%
  mutate(Genus = as.character(Genus)) %>%
  pivot_wider(
    names_from = AnalysisID,
    values_from = c(
      N_group1,
      N_group2,
      N_pairs,
      Detection_difference_percentage_points,
      P_raw,
      FDR_BH_all,
      Holm_all
    ),
    names_glue = "{.value}__{AnalysisID}"
  )

write.csv(
  detection_comparison_wide,
  file.path(out_root, "PRIMARY_vs_UNFILTERED_detection_statistics_WIDE.csv"),
  row.names = FALSE
)

root_readme <- c(
  "FINAL SECTION 4.6: TARGETED GENERA WITH PATIENT- AND BIOPSY-LEVEL PREVALENCE",
  "",
  paste0("Run date: ", Sys.time()),
  "",
  "PRIMARY",
  paste0("Input: ", analysis_specs[[1]]$physeq_path),
  paste0("Output: ", results[[1]]$out_dir),
  "",
  "SENSITIVITY",
  paste0("Input: ", analysis_specs[[2]]$physeq_path),
  paste0("Output: ", results[[2]]$out_dir),
  "The sensitivity input is post-selective-rarefaction and immediately before",
  "the 5% ASV-prevalence filter, so the sensitivity isolates that filter.",
  "",
  "DISPLAY RULES",
  paste0(
    "Abundance brackets/stars are shown only when ",
    ABUNDANCE_PLOT_PVALUE_COLUMN,
    " <= ",
    ABUNDANCE_SIGNIFICANCE_THRESHOLD,
    "."
  ),
  "No ns labels are shown.",
  "No significance annotations are shown on patient- or biopsy-level prevalence plots.",
  "Panel C is descriptive and displays detected/total biopsy counts.",
  "Both primary and sensitivity figures use the exact x-axis labels SSL_HP, H_ADJ_SSL_HP and H_COL.",
  "",
  "STATISTICS",
  "Each dataset is analysed separately using the identical patient-level pipeline.",
  "Every comparison retains all patient-level profiles: SSL_HP=23,",
  "H_ADJ_SSL_HP=10 and H_COL=8.",
  "All three pairwise comparisons are run for every genus. There is no matched",
  "subsetting in this script; the matched analysis remains separate.",
  "Abundance and patient-level detection are corrected as separate 9-test families",
  "within each dataset. Biopsy-level prevalence is not tested inferentially because",
  "biopsies are clustered within patients.",
  "The within-patient audit identifies profiles with both positive and negative",
  "biopsies; this is required before describing a comparison as within-patient.",
  "The unfiltered analysis is a sensitivity analysis, not a replacement for the",
  "prespecified 5% filtered primary analysis."
)

writeLines(root_readme, file.path(out_root, "READ_ME_FIRST.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))

msg("\nBoth Section 4.6 analyses completed.")
msg("Root output directory: ", out_root)
msg("Primary output: ", results[[1]]$out_dir)
msg("Sensitivity output: ", results[[2]]$out_dir)

write.csv(bind_rows(
  all_abundance_stats %>% transmute(AnalysisID, Outcome = "CLR", Genus, Contrast, Test, Status, Family_complete),
  all_detection_stats %>% transmute(AnalysisID, Outcome = "Detection", Genus, Contrast, Test, Status, Family_complete)),
  file.path(out_root, "RUN_STATUS.csv"), row.names = FALSE)
