#!/usr/bin/env Rscript
# Exact-site matched analyses and saved distance objects.
# Run script 144 for final beta-diversity inference from these distances.
rm(list = ls())
set.seed(20260823)

options(
  stringsAsFactors = FALSE,
  contrasts = c("contr.treatment", "contr.poly")
)

config <- new.env(parent = baseenv())
if (file.exists("config/config.R")) sys.source("config/config.R", envir = config)
work_dir <- if (exists("PATH_WORK", config, inherits = FALSE)) config$PATH_WORK else "outputs"
if (!dir.exists(work_dir)) stop("Input directory does not exist: ", work_dir)

base_dir <- work_dir

use_keep7699_object <- FALSE

physeq_without7699_path <- file.path(
  base_dir,
  "phyloseq_final_March26_5pct_without7699.rds"
)

physeq_keep7699_path <- file.path(
  base_dir,
  "phyloseq_final_March26_5pct_keep7699.rds"
)

physeq_path <- if (use_keep7699_object) {
  physeq_keep7699_path
} else {
  physeq_without7699_path
}

alpha_reads_without7699_path <- file.path(
  base_dir,
  "phyloseq_March26_post_rarefaction_without7699.rds"
)

alpha_reads_keep7699_path <- file.path(
  base_dir,
  "phyloseq_March26_post_rarefaction_keep7699.rds"
)

alpha_reads_path <- if (use_keep7699_object) {
  alpha_reads_keep7699_path
} else {
  alpha_reads_without7699_path
}

use_prevalence_unfiltered_reads_for_alpha <- TRUE

out_root <- file.path(
  base_dir,
  if (use_keep7699_object) {
    "Chapter4_MATCHED_SSLHP_HADJSSLHP_keep7699"
  } else {
    "Chapter4_MATCHED_SSLHP_HADJSSLHP_without7699"
  }
)

run_id <- paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
out_dir <- file.path(out_root, run_id)
if (dir.exists(out_dir)) stop("Output folder already exists: ", out_dir)

ambiguous_pair_action <- "exclude"

n_permutations <- 9999L
clr_pseudocount <- 0.5
alpha_level <- 0.05
ancombc_n_cores <- 1L

target_genera <- c(
  "Escherichia-Shigella",
  "Serratia",
  "Fusobacterium"
)

connect_pairs_in_alpha_plots <- TRUE
connect_pairs_in_beta_plots <- FALSE
ellipse_level <- 0.68

font_family <- "Arial"
export_dpi <- 600

group_levels <- c("H_ADJ_SSL_HP", "SSL_HP")

group_colours <- c(
  "H_ADJ_SSL_HP" = "#E9C61D",
  "SSL_HP" = "#FF0000"
)

group_fill_colours <- c(
  "H_ADJ_SSL_HP" = "#F6E79A",
  "SSL_HP" = "#FFB3B3"
)

required_packages <- c(
  "phyloseq",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "lme4",
  "lmerTest",
  "vegan",
  "permute",
  "DESeq2",
  "ggrepel",
  "gridExtra",
  "ANCOMBC"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following required package(s) before running the script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(vegan)
  library(permute)
  library(DESeq2)
  library(ggrepel)
  library(gridExtra)
  library(ANCOMBC)
})

# Version-specific missing-value alignment guard; retained for reproducibility.
patch_ancombc_bias_em_if_needed <- function() {
  namespace <- asNamespace("ANCOMBC")

  if (!exists(".bias_em", envir = namespace, inherits = FALSE)) {
    stop(
      "The installed ANCOMBC namespace does not contain .bias_em. ",
      "Do not run the analysis with an unknown implementation.",
      call. = FALSE
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

  if (grepl("neither_na", source_before, fixed = TRUE)) {
    return(
      data.frame(
        ANCOMBC_version = as.character(installed_version),
        Bias_EM_patch_required = FALSE,
        Bias_EM_patch_applied_in_memory = FALSE,
        Patch_status = "Installed ANCOMBC already contains the corrected joint NA mask",
        stringsAsFactors = FALSE
      )
    )
  }

  if (installed_version != package_version("2.4.0")) {
    stop(
      "ANCOMBC ",
      as.character(installed_version),
      " does not contain the corrected joint NA mask. ",
      "This script patches only the audited ANCOMBC 2.4.0 implementation. ",
      "Install a current official ANCOMBC release before running.",
      call. = FALSE
    )
  }

  body_parts <- as.list(body(bias_em_function))

  normalise_expression <- function(expression_object) {
    output <- paste(deparse(expression_object), collapse = "")
    output <- gsub("\\s+", "", output)
    gsub("<-", "=", output, fixed = TRUE)
  }

  if (length(body_parts) < 5L) {
    stop("Unexpected ANCOMBC 2.4.0 .bias_em structure.", call. = FALSE)
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
      "The ANCOMBC 2.4.0 .bias_em function does not match the audited ",
      "implementation. The namespace was not modified.",
      call. = FALSE
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

  body(bias_em_function) <- as.call(corrected_body_parts)

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
    stop("The in-memory ANCOMBC patch could not be verified.", call. = FALSE)
  }

  data.frame(
    ANCOMBC_version = as.character(installed_version),
    Bias_EM_patch_required = TRUE,
    Bias_EM_patch_applied_in_memory = TRUE,
    Patch_status = paste0(
      "Applied the corrected joint NA mask in memory; ",
      "installed package files were not modified"
    ),
    stringsAsFactors = FALSE
  )
}

if (capabilities("cairo")) {
  options(bitmapType = "cairo")
}

if (!ambiguous_pair_action %in% c("exclude", "stop")) {
  stop("ambiguous_pair_action must be either 'exclude' or 'stop'.", call. = FALSE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(out_dir, file.path(out_root, "LATEST_RUN_PATH.txt"))

subdirs <- list(
  audit = file.path(out_dir, "00_PAIRING_AND_QC"),
  alpha = file.path(out_dir, "01_ALPHA_DIVERSITY"),
  aitchison = file.path(out_dir, "02_AITCHISON"),
  bray = file.path(out_dir, "03_BRAY_CURTIS"),
  target_clr = file.path(out_dir, "04_TARGETED_CLR_AND_PREVALENCE"),
  deseq2_genus = file.path(out_dir, "05_DESEQ2_GENUS"),
  deseq2_species = file.path(out_dir, "06_DESEQ2_RESOLVED_SPECIES"),
  ancombc2_genus = file.path(out_dir, "07_ANCOMBC2_GENUS"),
  ancombc2_species = file.path(out_dir, "08_ANCOMBC2_RESOLVED_SPECIES"),
  combined = file.path(out_dir, "09_COMBINED_FIGURES")
)

invisible(lapply(subdirs, dir.create, recursive = TRUE, showWarnings = FALSE))

ancombc_patch_status <- patch_ancombc_bias_em_if_needed()

safe_patch_status_path <- file.path(
  subdirs$audit,
  "00_ANCOMBC_BIAS_EM_PATCH_STATUS.csv"
)

write.csv(
  ancombc_patch_status,
  safe_patch_status_path,
  row.names = FALSE
)

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

sample_id_key <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\.", "-", x)
  toupper(x)
}

clean_tax_name <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "NULL", "None")] <- NA_character_
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

clean_plain_genus <- function(x) canonicalise_genus(x)

format_adjusted_p <- function(p, prefix = "BH FDR") {
  if (length(p) == 0L || is.na(p) || !is.finite(p)) {
    return(paste0(prefix, " = NA"))
  }

  if (p < 0.001) {
    return(paste0(prefix, " < 0.001"))
  }

  paste0(prefix, " = ", formatC(p, format = "f", digits = 3))
}

exact_paired_binary_test <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  x <- as.integer(x[keep])
  y <- as.integer(y[keep])

  group1_only <- sum(x == 1L & y == 0L)
  group2_only <- sum(x == 0L & y == 1L)
  discordant <- group1_only + group2_only

  p_value <- if (discordant == 0L) {
    1
  } else {
    stats::binom.test(
      x = group1_only,
      n = discordant,
      p = 0.5,
      alternative = "two.sided"
    )$p.value
  }

  list(
    p_value = as.numeric(p_value),
    group1_only = group1_only,
    group2_only = group2_only,
    discordant = discordant
  )
}

spread_label_positions <- function(
    y,
    lower,
    upper,
    preferred_gap
) {
  y <- as.numeric(y)

  if (length(y) == 0L) return(numeric(0))
  if (length(y) == 1L) return(min(max(y, lower), upper))

  available <- upper - lower

  if (!is.finite(available) || available <= 0) {
    return(rep(lower, length(y)))
  }

  gap <- min(
    preferred_gap,
    available / (length(y) - 1L)
  )

  order_index <- order(y)
  positioned <- pmin(pmax(y[order_index], lower), upper)

  for (i in 2:length(positioned)) {
    positioned[i] <- max(positioned[i], positioned[i - 1L] + gap)
  }

  if (tail(positioned, 1L) > upper) {
    positioned <- positioned - (tail(positioned, 1L) - upper)
  }

  for (i in seq(length(positioned) - 1L, 1L, by = -1L)) {
    positioned[i] <- min(positioned[i], positioned[i + 1L] - gap)
  }

  if (positioned[1L] < lower) {
    positioned <- positioned + (lower - positioned[1L])
  }

  output <- numeric(length(y))
  output[order_index] <- positioned
  output
}

first_existing_column <- function(column_names, candidates) {
  hit <- candidates[candidates %in% column_names]
  if (length(hit) == 0L) return(NULL)
  hit[1L]
}

stop_with <- function(...) {
  stop(paste0(...), call. = FALSE)
}

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

safe_write_csv <- function(x, path) {
  write.csv(
    as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE),
    path,
    row.names = FALSE,
    na = ""
  )
}

pdf_device <- if (capabilities("cairo")) {
  grDevices::cairo_pdf
} else {
  grDevices::pdf
}

save_plot_set <- function(plot, directory, stem, width_mm, height_mm) {
  ggsave(
    file.path(directory, paste0(stem, ".pdf")),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = pdf_device,
    bg = "white"
  )

  ggsave(
    file.path(directory, paste0(stem, ".png")),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = export_dpi,
    bg = "white"
  )

  ggsave(
    file.path(directory, paste0(stem, ".jpg")),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = export_dpi,
    quality = 100,
    bg = "white"
  )
}

format_p_value <- function(p, prefix = "p") {
  if (length(p) == 0L || is.na(p) || !is.finite(p)) {
    return(paste0(prefix, " = NA"))
  }

  if (p < 0.0001) {
    paste0(prefix, " < 0.0001")
  } else {
    paste0(prefix, " = ", formatC(p, format = "f", digits = 4))
  }
}

p_to_label <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.0001) return("p < 0.0001")
  paste0("p = ", formatC(p, format = "f", digits = 4))
}

nature_theme <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      text = element_text(family = font_family, colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(colour = "black"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

make_group <- function(x) {
  x <- clean_character(x)

  x[x %in% c("SSL", "HP", "SSL_HP")] <- "SSL_HP"

  x[x %in% c(
    "H_ADJ_SSL",
    "H_ADJ_HP",
    "H_ADJ_SSL_HP"
  )] <- "H_ADJ_SSL_HP"

  x[!x %in% group_levels] <- NA_character_
  factor(x, levels = group_levels)
}

canonical_site <- function(x) {
  x <- clean_character(x)
  if (length(x) == 0L) return(character(0))

  z <- toupper(x)
  z <- gsub("[^A-Z0-9]", "", z)

  out <- rep(NA_character_, length(z))

  out[z %in% c("CAEC", "CAECUM", "CECUM", "CECAL")] <- "CAEC"
  out[z %in% c("ASC", "ASCENDING", "ASCENDINGCOLON")] <- "ASC"
  out[z %in% c("TRVS", "TRANSVERSE", "TRANSVERSECOLON")] <- "TRVS"
  out[z %in% c("DESC", "DESCENDING", "DESCENDINGCOLON")] <- "DESC"
  out[z %in% c("SIG", "SIGMOID", "SIGMOIDCOLON")] <- "SIG"
  out[z %in% c("RECT", "RECTUM", "RECTAL")] <- "RECT"

  out
}

extract_patient_from_sample_id <- function(sample_id) {
  out <- sub(".*_([0-9]+)$", "\\1", sample_id)
  out[out == sample_id] <- NA_character_
  clean_character(out)
}

extract_site_from_sample_id <- function(sample_id) {
  site_codes <- c("CAEC", "ASC", "TRVS", "DESC", "SIG", "RECT")

  vapply(
    sample_id,
    function(id) {
      hits <- site_codes[
        vapply(
          site_codes,
          function(site) grepl(paste0("_", site, "_"), id, fixed = TRUE),
          logical(1)
        )
      ]

      if (length(hits) == 1L) hits else NA_character_
    },
    character(1)
  )
}

as_taxa_rows_matrix <- function(phy) {
  x <- as(otu_table(phy), "matrix")
  if (!taxa_are_rows(phy)) x <- t(x)
  x
}

as_samples_rows_matrix <- function(phy) {
  t(as_taxa_rows_matrix(phy))
}

check_raw_integer_counts <- function(phy, analysis_name = "analysis") {
  x <- as_taxa_rows_matrix(phy)

  if (
    any(!is.finite(x)) ||
    any(x < 0) ||
    any(abs(x - round(x)) > 1e-8)
  ) {
    stop_with(
      analysis_name,
      " requires non-negative integer raw counts, but the abundance table does not satisfy this."
    )
  }

  invisible(TRUE)
}

make_blocked_permutation_control <- function(blocks, nperm) {
  control <- permute::how(nperm = nperm)
  setBlocks(control) <- factor(blocks)
  control
}

extract_lrt_p <- function(reduced_fit, full_fit) {
  tab <- anova(reduced_fit, full_fit)
  p_col <- grep("Pr\\(>Chisq\\)", colnames(tab), value = TRUE)

  if (length(p_col) != 1L || nrow(tab) < 2L) return(NA_real_)
  as.numeric(tab[2L, p_col])
}

get_convergence_message <- function(fit) {
  messages <- fit@optinfo$conv$lme4$messages
  if (is.null(messages)) return("")
  paste(messages, collapse = " | ")
}

clean_taxonomy_value <- function(x) {
  x <- clean_character(x)

  bad <- is.na(x) |
    grepl(
      "^(NA|N/A|UNKNOWN|UNCLASSIFIED|UNCULTURED|UNIDENTIFIED|UNASSIGNED)$",
      toupper(x)
    )

  x[bad] <- NA_character_
  x
}

make_lineage_aware_genus <- function(tax_df) {
  tax_df <- as.data.frame(
    tax_df,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  required_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

  for (rank in required_ranks) {
    if (!rank %in% colnames(tax_df)) tax_df[[rank]] <- NA_character_
    tax_df[[rank]] <- clean_taxonomy_value(tax_df[[rank]])
  }

  fallback_parent <- rep("Bacteria", nrow(tax_df))

  for (rank in c("Family", "Order", "Class", "Phylum", "Kingdom")) {
    use <- fallback_parent == "Bacteria" & !is.na(tax_df[[rank]])
    fallback_parent[use] <- tax_df[[rank]][use]
  }

  genus_base <- tax_df$Genus
  missing_genus <- is.na(genus_base)
  genus_base[missing_genus] <- paste0(
    "Unclassified_",
    gsub("[^A-Za-z0-9]+", "_", fallback_parent[missing_genus])
  )

  lineage_key <- apply(
    tax_df[, c("Phylum", "Class", "Order", "Family"), drop = FALSE],
    1,
    function(z) paste(ifelse(is.na(z), "NA", z), collapse = "|")
  )

  genus_label <- genus_base

  for (this_genus in unique(genus_base)) {
    idx <- which(genus_base == this_genus)
    distinct_lineages <- sort(unique(lineage_key[idx]))

    if (length(distinct_lineages) > 1L) {
      lineage_number <- match(lineage_key[idx], distinct_lineages)
      genus_label[idx] <- paste0(this_genus, "__lineage", lineage_number)
    }
  }

  mapping <- data.frame(
    Original_taxon = rownames(tax_df),
    Original_genus = tax_df$Genus,
    Fallback_parent = fallback_parent,
    Lineage_key = lineage_key,
    Genus_for_analysis = genus_label,
    stringsAsFactors = FALSE
  )

  list(genus = genus_label, mapping = mapping)
}

prepare_lineage_aware_asv_phyloseq <- function(phy, mapping_path = NULL) {
  tax_df <- as.data.frame(
    tax_table(phy),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!"Genus" %in% colnames(tax_df)) {
    stop("The taxonomy table has no Genus column.", call. = FALSE)
  }

  genus_result <- make_lineage_aware_genus(tax_df)
  tax_df$Genus <- genus_result$genus
  tax_table(phy) <- tax_table(as.matrix(tax_df))

  if (!is.null(mapping_path)) {
    safe_write_csv(genus_result$mapping, mapping_path)
  }

  phy
}

prepare_genus_phyloseq <- function(phy, mapping_path = NULL) {
  phy_labeled <- prepare_lineage_aware_asv_phyloseq(phy, mapping_path)

  phy_genus <- tax_glom(
    phy_labeled,
    taxrank = "Genus",
    NArm = FALSE
  )

  phy_genus <- prune_taxa(taxa_sums(phy_genus) > 0, phy_genus)

  genus_names <- as.character(tax_table(phy_genus)[, "Genus"])

  if (anyNA(genus_names) || any(genus_names == "")) {
    stop("Missing genus labels remain after lineage-aware genus aggregation.", call. = FALSE)
  }

  if (anyDuplicated(genus_names)) {
    stop("Duplicated genus labels remain after tax_glom().", call. = FALSE)
  }

  taxa_names(phy_genus) <- genus_names
  phy_genus
}

if (!file.exists(physeq_path)) {
  stop_with("The selected phyloseq object was not found: ", physeq_path)
}

phy_all <- readRDS(physeq_path)

if (!inherits(phy_all, "phyloseq")) {
  stop("The input RDS is not a phyloseq object.", call. = FALSE)
}

if (
  is.null(otu_table(phy_all, errorIfNULL = FALSE)) ||
  is.null(sample_data(phy_all, errorIfNULL = FALSE)) ||
  is.null(tax_table(phy_all, errorIfNULL = FALSE))
) {
  stop(
    "The phyloseq object must contain raw counts, sample metadata and taxonomy.",
    call. = FALSE
  )
}

check_raw_integer_counts(phy_all, "The complete matched workflow")

metadata <- data.frame(
  sample_data(phy_all),
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  rownames_to_column("SampleID")

sample_type_col <- first_existing_column(
  colnames(metadata),
  c("Sample_Type", "SampleType", "sample_type", "Group")
)

if (is.null(sample_type_col)) {
  stop("No Sample_Type/Group metadata column was found.", call. = FALSE)
}

patient_col <- first_existing_column(
  colnames(metadata),
  c("Patient", "patient", "Patient_ID", "patient_id", "PatientID")
)

location_col <- first_existing_column(
  colnames(metadata),
  c(
    "Location_Colon",
    "LocationColon",
    "Colon_Location",
    "Location",
    "location"
  )
)

metadata$Group <- make_group(metadata[[sample_type_col]])

metadata$Patient_metadata <- if (!is.null(patient_col)) {
  clean_character(metadata[[patient_col]])
} else {
  NA_character_
}

metadata$Patient_from_SampleID <- extract_patient_from_sample_id(metadata$SampleID)
metadata$Patient <- dplyr::coalesce(
  metadata$Patient_metadata,
  metadata$Patient_from_SampleID
)

metadata$Site_metadata <- if (!is.null(location_col)) {
  canonical_site(metadata[[location_col]])
} else {
  NA_character_
}

metadata$Site_from_SampleID <- extract_site_from_sample_id(metadata$SampleID)
metadata$Site <- dplyr::coalesce(
  metadata$Site_metadata,
  metadata$Site_from_SampleID
)

metadata$PairID <- ifelse(
  !is.na(metadata$Patient) & !is.na(metadata$Site),
  paste(metadata$Patient, metadata$Site, sep = "__"),
  NA_character_
)

metadata$Reads_in_5pct_object <- as.numeric(
  sample_sums(phy_all)[metadata$SampleID]
)

patient_conflicts <- metadata %>%
  filter(
    !is.na(Patient_metadata),
    !is.na(Patient_from_SampleID),
    Patient_metadata != Patient_from_SampleID
  )

site_conflicts <- metadata %>%
  filter(
    !is.na(Site_metadata),
    !is.na(Site_from_SampleID),
    Site_metadata != Site_from_SampleID
  )

safe_write_csv(
  patient_conflicts,
  file.path(subdirs$audit, "00A_patient_metadata_vs_sampleID_conflicts.csv")
)

safe_write_csv(
  site_conflicts,
  file.path(subdirs$audit, "00B_site_metadata_vs_sampleID_conflicts.csv")
)

if (nrow(patient_conflicts) > 0L || nrow(site_conflicts) > 0L) {
  stop_with(
    "Patient or exact colonic-site metadata disagree with the sample identifiers. ",
    "Inspect the 00A/00B audit files before matching."
  )
}

contains_7699 <- any(metadata$Patient == "7699", na.rm = TRUE)

if (use_keep7699_object && !contains_7699) {
  warning(
    "use_keep7699_object is TRUE, but patient 7699 was not detected in this input.",
    call. = FALSE
  )
}

if (!use_keep7699_object && contains_7699) {
  stop_with(
    "The without7699 branch was selected, but patient 7699 is present. ",
    "Check physeq_path or set use_keep7699_object <- TRUE."
  )
}

safe_write_csv(
  metadata,
  file.path(subdirs$audit, "00_all_samples_metadata_mapping.csv")
)

candidate_samples <- metadata %>%
  filter(!is.na(Group))

missing_pair_fields <- candidate_samples %>%
  filter(is.na(Patient) | is.na(Site) | is.na(PairID))

safe_write_csv(
  missing_pair_fields,
  file.path(subdirs$audit, "01_candidate_samples_missing_patient_or_exact_site.csv")
)

pairable_samples <- candidate_samples %>%
  filter(!is.na(PairID)) %>%
  mutate(Group = factor(Group, levels = group_levels))

pair_counts_long <- pairable_samples %>%
  count(PairID, Patient, Site, Group, name = "N_samples")

pair_counts_wide <- pair_counts_long %>%
  tidyr::pivot_wider(
    names_from = Group,
    values_from = N_samples,
    values_fill = 0
  )

for (group_name in group_levels) {
  if (!group_name %in% colnames(pair_counts_wide)) {
    pair_counts_wide[[group_name]] <- 0L
  }
}

pair_counts_wide <- pair_counts_wide %>%
  mutate(
    Pair_status = case_when(
      .data[["H_ADJ_SSL_HP"]] == 1L & .data[["SSL_HP"]] == 1L ~
        "VALID_EXACT_ONE_TO_ONE_PAIR",
      .data[["H_ADJ_SSL_HP"]] == 0L | .data[["SSL_HP"]] == 0L ~
        "UNMATCHED_MISSING_ONE_GROUP",
      TRUE ~
        "AMBIGUOUS_MULTIPLE_SAMPLES_IN_AT_LEAST_ONE_GROUP"
    )
  ) %>%
  arrange(Patient, Site)

safe_write_csv(
  pair_counts_long,
  file.path(subdirs$audit, "02_pair_counts_long.csv")
)

safe_write_csv(
  pair_counts_wide,
  file.path(subdirs$audit, "03_pair_status_all_patient_site_keys.csv")
)

ambiguous_pairs <- pair_counts_wide %>%
  filter(Pair_status == "AMBIGUOUS_MULTIPLE_SAMPLES_IN_AT_LEAST_ONE_GROUP")

unmatched_pairs <- pair_counts_wide %>%
  filter(Pair_status == "UNMATCHED_MISSING_ONE_GROUP")

valid_pairs <- pair_counts_wide %>%
  filter(Pair_status == "VALID_EXACT_ONE_TO_ONE_PAIR")

safe_write_csv(
  ambiguous_pairs,
  file.path(subdirs$audit, "04_ambiguous_patient_site_keys.csv")
)

safe_write_csv(
  unmatched_pairs,
  file.path(subdirs$audit, "05_unmatched_patient_site_keys.csv")
)

safe_write_csv(
  valid_pairs,
  file.path(subdirs$audit, "06_valid_exact_patient_site_pairs.csv")
)

if (nrow(ambiguous_pairs) > 0L && ambiguous_pair_action == "stop") {
  stop_with(
    nrow(ambiguous_pairs),
    " Patient__Site key(s) contain multiple samples in at least one group. ",
    "Inspect 04_ambiguous_patient_site_keys.csv."
  )
}

if (nrow(ambiguous_pairs) > 0L && ambiguous_pair_action == "exclude") {
  warning(
    nrow(ambiguous_pairs),
    " ambiguous Patient__Site key(s) were excluded rather than paired arbitrarily. ",
    "See 04_ambiguous_patient_site_keys.csv.",
    call. = FALSE
  )
}

if (nrow(valid_pairs) < 2L) {
  stop_with(
    "Fewer than two exact one-to-one Patient__Site pairs were found. ",
    "Inspect the pairing audit files."
  )
}

matched_manifest <- pairable_samples %>%
  filter(PairID %in% valid_pairs$PairID) %>%
  arrange(PairID, Group)

pair_validation <- matched_manifest %>%
  count(PairID, Group, name = "N") %>%
  tidyr::pivot_wider(names_from = Group, values_from = N, values_fill = 0)

if (
  nrow(pair_validation) != nrow(valid_pairs) ||
  any(pair_validation$H_ADJ_SSL_HP != 1L) ||
  any(pair_validation$SSL_HP != 1L)
) {
  stop("Internal exact-pair validation failed.", call. = FALSE)
}

matched_sample_ids <- matched_manifest$SampleID

phy_matched <- prune_samples(matched_sample_ids, phy_all)
phy_matched <- prune_taxa(taxa_sums(phy_matched) > 0, phy_matched)

matched_manifest <- matched_manifest %>%
  mutate(
    Patient = factor(Patient),
    Site = factor(Site),
    PairID = factor(PairID),
    Group = factor(Group, levels = group_levels)
  ) %>%
  as.data.frame(
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

if (anyDuplicated(matched_manifest$SampleID)) {
  duplicated_manifest <- matched_manifest[
    duplicated(matched_manifest$SampleID) |
      duplicated(matched_manifest$SampleID, fromLast = TRUE),
    ,
    drop = FALSE
  ]

  safe_write_csv(
    duplicated_manifest,
    file.path(subdirs$audit, "07A_duplicated_SampleID_in_matched_manifest.csv")
  )

  stop(
    "Duplicated SampleID values were found in the matched manifest. ",
    "Inspect 07A_duplicated_SampleID_in_matched_manifest.csv.",
    call. = FALSE
  )
}

phy_sample_order <- sample_names(phy_matched)
manifest_index <- match(phy_sample_order, matched_manifest$SampleID)
manifest_alignment_method <- "exact SampleID"

if (anyNA(manifest_index)) {
  manifest_keys <- sample_id_key(matched_manifest$SampleID)
  phy_keys <- sample_id_key(phy_sample_order)

  if (anyDuplicated(manifest_keys) || anyDuplicated(phy_keys)) {
    manifest_index <- rep(NA_integer_, length(phy_sample_order))
  } else {
    manifest_index <- match(phy_keys, manifest_keys)
    manifest_alignment_method <- "normalised SampleID key"
  }
}

if (anyNA(manifest_index)) {
  manifest_alignment_audit <- data.frame(
    Phyloseq_SampleID = phy_sample_order,
    Phyloseq_SampleID_key = sample_id_key(phy_sample_order),
    Matched_manifest_SampleID = ifelse(
      is.na(manifest_index),
      NA_character_,
      matched_manifest$SampleID[manifest_index]
    ),
    stringsAsFactors = FALSE
  )

  safe_write_csv(
    manifest_alignment_audit,
    file.path(subdirs$audit, "07B_matched_manifest_alignment_failure.csv")
  )

  stop(
    "The matched manifest could not be aligned to the matched phyloseq object. ",
    "Inspect 07B_matched_manifest_alignment_failure.csv.",
    call. = FALSE
  )
}

matched_manifest <- matched_manifest[manifest_index, , drop = FALSE]
rownames(matched_manifest) <- phy_sample_order

if (!identical(as.character(matched_manifest$SampleID), phy_sample_order)) {
  stop("Internal matched-manifest ordering validation failed.", call. = FALSE)
}

phy_sample_metadata <- matched_manifest[
  ,
  setdiff(colnames(matched_manifest), "SampleID"),
  drop = FALSE
]
rownames(phy_sample_metadata) <- phy_sample_order
sample_data(phy_matched) <- sample_data(phy_sample_metadata)

sample_metadata_check <- data.frame(
  sample_data(phy_matched),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!identical(rownames(sample_metadata_check), phy_sample_order)) {
  stop("sample_data alignment failed after constructing the matched object.", call. = FALSE)
}

safe_write_csv(
  data.frame(
    SampleID = phy_sample_order,
    Alignment_method = manifest_alignment_method,
    stringsAsFactors = FALSE
  ),
  file.path(subdirs$audit, "07C_matched_manifest_alignment_method.csv")
)

if (nsamples(phy_matched) != 2L * nrow(valid_pairs)) {
  stop("The matched phyloseq object does not contain exactly two samples per pair.", call. = FALSE)
}

check_raw_integer_counts(phy_matched, "The exact matched subset")

safe_write_csv(
  matched_manifest,
  file.path(subdirs$audit, "07_FINAL_MATCHED_SAMPLE_MANIFEST.csv")
)

saveRDS(
  phy_matched,
  file.path(subdirs$audit, "phyloseq_EXACT_MATCHED_5pct_subset.rds")
)

pairing_summary <- data.frame(
  Input_object = physeq_path,
  Input_samples = nsamples(phy_all),
  Input_ASVs = ntaxa(phy_all),
  Candidate_SSL_HP_samples = sum(candidate_samples$Group == "SSL_HP", na.rm = TRUE),
  Candidate_H_ADJ_SSL_HP_samples = sum(
    candidate_samples$Group == "H_ADJ_SSL_HP",
    na.rm = TRUE
  ),
  Exact_pairs = nrow(valid_pairs),
  Matched_samples = nsamples(phy_matched),
  Matched_patients = n_distinct(matched_manifest$Patient),
  Matched_ASVs_nonzero = ntaxa(phy_matched),
  Unmatched_patient_site_keys = nrow(unmatched_pairs),
  Ambiguous_patient_site_keys = nrow(ambiguous_pairs),
  stringsAsFactors = FALSE
)

safe_write_csv(
  pairing_summary,
  file.path(subdirs$audit, "08_pairing_summary.csv")
)

if (use_keep7699_object) {
  has_7699_trvs_pair <- "7699__TRVS" %in% as.character(valid_pairs$PairID)

  if (!has_7699_trvs_pair) {
    warning(
      "The keep-7699 object did not yield a valid exact 7699__TRVS pair. ",
      "The pairing audit shows whether one sample was absent or the key was ambiguous.",
      call. = FALSE
    )
  }
}

cat("\n===== EXACT MATCHING SUMMARY =====\n")
print(pairing_summary)
cat("\nFinal exact pairs:\n")
print(valid_pairs[, c("PairID", "Patient", "Site")])

target_sample_ids <- sample_names(phy_matched)
target_sample_keys <- sample_id_key(target_sample_ids)

reads_5pct <- as.numeric(sample_sums(phy_matched))
names(reads_5pct) <- target_sample_ids

alpha_reads_source <- "5% matched object sample_sums"
alpha_reads_alignment_method <- "same matched phyloseq object"
alpha_read_totals <- reads_5pct
alpha_source_sample_ids <- target_sample_ids

if (
  isTRUE(use_prevalence_unfiltered_reads_for_alpha) &&
  file.exists(alpha_reads_path)
) {
  phy_reads <- readRDS(alpha_reads_path)

  if (!inherits(phy_reads, "phyloseq")) {
    stop("The alpha read-depth source is not a phyloseq object.", call. = FALSE)
  }

  source_sample_ids <- sample_names(phy_reads)
  source_sample_keys <- sample_id_key(source_sample_ids)

  read_index <- match(target_sample_ids, source_sample_ids)
  read_alignment_method <- "exact SampleID"

  if (anyNA(read_index)) {
    if (
      anyDuplicated(target_sample_keys) ||
      anyDuplicated(source_sample_keys)
    ) {
      read_index <- rep(NA_integer_, length(target_sample_ids))
    } else {
      read_index <- match(target_sample_keys, source_sample_keys)
      read_alignment_method <- "normalised SampleID key"
    }
  }

  if (anyNA(read_index)) {
    read_alignment_failure <- data.frame(
      SampleID = target_sample_ids,
      SampleID_key = target_sample_keys,
      Matched_source_SampleID = ifelse(
        is.na(read_index),
        NA_character_,
        source_sample_ids[read_index]
      ),
      stringsAsFactors = FALSE
    )

    safe_write_csv(
      read_alignment_failure,
      file.path(subdirs$alpha, "Alpha_prevalence_unfiltered_read_alignment_failure.csv")
    )

    warning(
      "The preferred pre-prevalence read-depth object could not be aligned to ",
      "every matched sample. The script will stop before modelling. ",
      "See Alpha_prevalence_unfiltered_read_alignment_failure.csv.",
      call. = FALSE
    )
  } else {
    source_reads <- as.numeric(sample_sums(phy_reads))
    names(source_reads) <- source_sample_ids

    aligned_source_ids <- source_sample_ids[read_index]
    preferred_reads <- as.numeric(source_reads[aligned_source_ids])
    names(preferred_reads) <- target_sample_ids

    if (
      any(!is.finite(preferred_reads)) ||
      any(preferred_reads <= 0)
    ) {
      bad_preferred_reads <- data.frame(
        SampleID = target_sample_ids,
        Source_SampleID = aligned_source_ids,
        Preferred_Reads = preferred_reads,
        stringsAsFactors = FALSE
      )

      safe_write_csv(
        bad_preferred_reads[
          !is.finite(bad_preferred_reads$Preferred_Reads) |
            bad_preferred_reads$Preferred_Reads <= 0,
          ,
          drop = FALSE
        ],
        file.path(subdirs$alpha, "Alpha_invalid_prevalence_unfiltered_read_totals.csv")
      )

      warning(
        "The preferred pre-prevalence read totals contained missing, non-finite ",
        "or zero values. The script will stop before modelling. ",
        "See Alpha_invalid_prevalence_unfiltered_read_totals.csv.",
        call. = FALSE
      )
    } else {
      alpha_read_totals <- preferred_reads
      alpha_reads_source <- alpha_reads_path
      alpha_reads_alignment_method <- read_alignment_method
      alpha_source_sample_ids <- aligned_source_ids
    }
  }
} else if (isTRUE(use_prevalence_unfiltered_reads_for_alpha)) {
  warning(
    "The preferred pre-prevalence alpha read-depth object was not found: ",
    alpha_reads_path,
    ". The script will stop before modelling.",
    call. = FALSE
  )
}

if (!identical(alpha_reads_source, alpha_reads_path)) {
  stop("The common pre-prevalence read-depth source could not be used; the filtered-read fallback is not used for this analysis.")
}

alpha_read_totals <- as.numeric(alpha_read_totals[target_sample_ids])
names(alpha_read_totals) <- target_sample_ids

if (
  any(!is.finite(alpha_read_totals)) ||
  any(alpha_read_totals <= 0)
) {
  invalid_alpha_reads <- data.frame(
    SampleID = target_sample_ids,
    Reads = as.numeric(alpha_read_totals),
    stringsAsFactors = FALSE
  )

  safe_write_csv(
    invalid_alpha_reads[
      !is.finite(invalid_alpha_reads$Reads) |
        invalid_alpha_reads$Reads <= 0,
      ,
      drop = FALSE
    ],
    file.path(subdirs$alpha, "Alpha_invalid_final_read_totals.csv")
  )

  stop(
    "The final alpha-diversity read totals contain missing, non-finite or zero ",
    "values. Inspect Alpha_invalid_final_read_totals.csv.",
    call. = FALSE
  )
}

alpha_read_audit <- data.frame(
  SampleID = target_sample_ids,
  SampleID_key = target_sample_keys,
  Source_SampleID = alpha_source_sample_ids,
  Reads_used_for_alpha = as.numeric(alpha_read_totals[target_sample_ids]),
  Reads_in_5pct_object = as.numeric(reads_5pct[target_sample_ids]),
  Read_source = alpha_reads_source,
  Alignment_method = alpha_reads_alignment_method,
  stringsAsFactors = FALSE
)

safe_write_csv(
  alpha_read_audit,
  file.path(subdirs$alpha, "Alpha_read_depth_source_audit.csv")
)

alpha_directory <- subdirs$alpha
alpha_figure_directory <- subdirs$combined
phy_prevalence_unfiltered <- prune_samples(alpha_source_sample_ids, phy_reads)
sample_names(phy_prevalence_unfiltered) <- target_sample_ids[match(sample_names(phy_prevalence_unfiltered), alpha_source_sample_ids)]
stopifnot(setequal(sample_names(phy_prevalence_unfiltered), target_sample_ids))
for (alpha_dataset in c("No_prevalence_filter", "Filtered_5pct")) {
  phy_alpha <- if (alpha_dataset == "No_prevalence_filter") phy_prevalence_unfiltered else phy_matched
  subdirs$alpha <- if (alpha_dataset == "No_prevalence_filter") file.path(alpha_directory, alpha_dataset) else alpha_directory
  subdirs$combined <- if (alpha_dataset == "No_prevalence_filter") file.path(alpha_figure_directory, alpha_dataset) else alpha_figure_directory
  dir.create(subdirs$alpha, recursive = TRUE, showWarnings = FALSE)
  dir.create(subdirs$combined, recursive = TRUE, showWarnings = FALSE)
alpha_sample_order <- sample_names(phy_alpha)
alpha_counts <- as_samples_rows_matrix(phy_alpha)
alpha_counts <- alpha_counts[alpha_sample_order, , drop = FALSE]

zero_count_samples <- rowSums(alpha_counts) <= 0

if (any(zero_count_samples)) {
  safe_write_csv(
    data.frame(
      SampleID = alpha_sample_order[zero_count_samples],
      Counts_in_alpha_dataset = rowSums(alpha_counts)[zero_count_samples],
      stringsAsFactors = FALSE
    ),
    file.path(subdirs$alpha, "Alpha_zero_count_samples_in_matched_object.csv")
  )

  stop(
    "At least one matched sample has zero counts in the selected alpha dataset. Inspect ",
    "Alpha_zero_count_samples_in_matched_object.csv.",
    call. = FALSE
  )
}

alpha_estimates <- as.data.frame(
  estimate_richness(
    phy_alpha,
    measures = c("Observed", "Shannon", "Simpson")
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

alpha_estimates$SampleID_alpha <- rownames(alpha_estimates)
alpha_estimates$SampleID_key <- sample_id_key(alpha_estimates$SampleID_alpha)

if (anyDuplicated(alpha_estimates$SampleID_key)) {
  stop(
    "Duplicated normalised SampleID keys were returned by estimate_richness().",
    call. = FALSE
  )
}

alpha_index <- match(
  sample_id_key(alpha_sample_order),
  alpha_estimates$SampleID_key
)

alpha_alignment_audit <- data.frame(
  SampleID = alpha_sample_order,
  SampleID_key = sample_id_key(alpha_sample_order),
  Alpha_SampleID = ifelse(
    is.na(alpha_index),
    NA_character_,
    alpha_estimates$SampleID_alpha[alpha_index]
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  alpha_alignment_audit,
  file.path(subdirs$alpha, "Alpha_sample_alignment_audit.csv")
)

if (anyNA(alpha_index)) {
  stop(
    "estimate_richness() output could not be aligned to every matched sample. ",
    "Inspect Alpha_sample_alignment_audit.csv.",
    call. = FALSE
  )
}

manifest_index_alpha <- match(alpha_sample_order, matched_manifest$SampleID)

if (anyNA(manifest_index_alpha)) {
  stop(
    "The matched manifest could not be realigned for alpha diversity.",
    call. = FALSE
  )
}

alpha_df <- matched_manifest[
  manifest_index_alpha,
  c("SampleID", "PairID", "Patient", "Site", "Group"),
  drop = FALSE
]

alpha_df$Reads <- as.numeric(alpha_read_totals[alpha_sample_order])
alpha_df$Observed <- as.numeric(alpha_estimates$Observed[alpha_index])
alpha_df$Shannon <- as.numeric(alpha_estimates$Shannon[alpha_index])
alpha_df$Simpson <- as.numeric(alpha_estimates$Simpson[alpha_index])

alpha_df$PairID <- factor(alpha_df$PairID)
alpha_df$Patient <- factor(alpha_df$Patient)
alpha_df$Site <- factor(alpha_df$Site)
alpha_df$Group <- factor(alpha_df$Group, levels = group_levels)

required_alpha_columns <- c(
  "SampleID",
  "PairID",
  "Patient",
  "Site",
  "Group",
  "Reads",
  "Observed",
  "Shannon",
  "Simpson"
)

alpha_missing_matrix <- do.call(
  cbind,
  lapply(
    required_alpha_columns,
    function(column_name) {
      value <- alpha_df[[column_name]]

      if (is.numeric(value)) {
        is.na(value) | !is.finite(value)
      } else {
        is.na(value) | trimws(as.character(value)) == ""
      }
    }
  )
)
colnames(alpha_missing_matrix) <- required_alpha_columns

alpha_bad_rows <- apply(alpha_missing_matrix, 1, any)

if (any(alpha_bad_rows)) {
  incomplete <- alpha_df[alpha_bad_rows, , drop = FALSE]
  incomplete$Missing_fields <- apply(
    alpha_missing_matrix[alpha_bad_rows, , drop = FALSE],
    1,
    function(z) paste(colnames(alpha_missing_matrix)[z], collapse = "; ")
  )

  safe_write_csv(
    incomplete,
    file.path(subdirs$alpha, "Alpha_samples_with_missing_required_data.csv")
  )

  stop(
    "Missing or non-finite values remain in the matched alpha-diversity input. ",
    "Affected field(s): ",
    paste(sort(unique(unlist(strsplit(incomplete$Missing_fields, "; ", fixed = TRUE)))), collapse = ", "),
    ". Inspect Alpha_samples_with_missing_required_data.csv.",
    call. = FALSE
  )
}

if (any(alpha_df$Reads <= 0)) {
  stop("At least one matched sample has zero reads.", call. = FALSE)
}

if (!identical(as.character(alpha_df$SampleID), alpha_sample_order)) {
  stop("Final alpha-diversity sample order validation failed.", call. = FALSE)
}

safe_write_csv(
  alpha_df,
  file.path(subdirs$alpha, "Alpha_EXACT_MATCHED_analysis_input.csv")
)

alpha_metrics <- c("Observed", "Shannon", "Simpson")

alpha_titles <- c(
  "Observed" = "Observed ASV richness",
  "Shannon" = "Shannon diversity",
  "Simpson" = "Simpson diversity"
)

alpha_y_labels <- c(
  "Observed" = "Observed ASVs",
  "Shannon" = "Shannon index",
  "Simpson" = "Simpson index"
)

fit_alpha_metric <- function(metric) {
  full_formula <- as.formula(
    paste0(metric, " ~ Group + Reads + (1 | Patient)")
  )

  reduced_group_formula <- as.formula(
    paste0(metric, " ~ Reads + (1 | Patient)")
  )

  reduced_reads_formula <- as.formula(
    paste0(metric, " ~ Group + (1 | Patient)")
  )

  full_fit <- lmerTest::lmer(
    full_formula,
    data = alpha_df,
    REML = FALSE,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)
    )
  )

  reduced_group_fit <- lmerTest::lmer(
    reduced_group_formula,
    data = alpha_df,
    REML = FALSE,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)
    )
  )

  reduced_reads_fit <- lmerTest::lmer(
    reduced_reads_formula,
    data = alpha_df,
    REML = FALSE,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)
    )
  )

  coefficient_table <- as.data.frame(coef(summary(full_fit))) %>%
    rownames_to_column("Term")

  group_row <- coefficient_table %>%
    filter(Term == "GroupSSL_HP")

  if (nrow(group_row) != 1L) {
    stop_with("Could not identify the SSL_HP group coefficient for ", metric, ".")
  }

  overall_row <- data.frame(
    Metric = metric,
    Formula = paste(deparse(full_formula), collapse = ""),
    Dataset = alpha_dataset,
    N_patients = n_distinct(alpha_df$Patient),
    N_pairs = n_distinct(alpha_df$PairID),
    N_samples = nrow(alpha_df),
    Group_LRT_P = extract_lrt_p(reduced_group_fit, full_fit),
    Reads_LRT_P = extract_lrt_p(reduced_reads_fit, full_fit),
    Singular = lme4::isSingular(full_fit, tol = 1e-4),
    Convergence_message = get_convergence_message(full_fit),
    stringsAsFactors = FALSE
  )

  contrast_row <- data.frame(
    Metric = metric,
    Comparison = "SSL_HP vs H_ADJ_SSL_HP",
    Estimate_SSL_HP_minus_H_ADJ_SSL_HP = group_row$Estimate,
    SE = group_row$`Std. Error`,
    Df = group_row$df,
    T_value = group_row$`t value`,
    P_raw = group_row$`Pr(>|t|)`,
    P_Holm = group_row$`Pr(>|t|)`,
    Adjustment_note = "Only one prespecified contrast exists within this metric; Holm equals raw P.",
    stringsAsFactors = FALSE
  )

  list(
    fit = full_fit,
    overall = overall_row,
    contrast = contrast_row,
    coefficients = coefficient_table
  )
}

alpha_fits <- setNames(lapply(alpha_metrics, fit_alpha_metric), alpha_metrics)
saveRDS(alpha_fits, file.path(subdirs$alpha, "Alpha_models.rds"))

alpha_overall_results <- bind_rows(lapply(alpha_fits, `[[`, "overall"))
alpha_contrast_results <- bind_rows(lapply(alpha_fits, `[[`, "contrast"))
alpha_coefficient_results <- bind_rows(
  lapply(
    names(alpha_fits),
    function(metric) {
      alpha_fits[[metric]]$coefficients %>% mutate(Metric = metric, .before = 1)
    }
  )
)

alpha_descriptive <- alpha_df %>%
  group_by(Group) %>%
  summarise(
    N_samples = n(),
    N_pairs = n_distinct(PairID),
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

safe_write_csv(
  alpha_overall_results,
  file.path(subdirs$alpha, "Alpha_EXACT_MATCHED_overall_results.csv")
)

safe_write_csv(
  alpha_contrast_results,
  file.path(subdirs$alpha, "Alpha_EXACT_MATCHED_SSL_HP_vs_H_ADJ_SSL_HP.csv")
)

safe_write_csv(
  alpha_coefficient_results,
  file.path(subdirs$alpha, "Alpha_EXACT_MATCHED_full_coefficient_tables.csv")
)

safe_write_csv(
  alpha_descriptive,
  file.path(subdirs$alpha, "Alpha_EXACT_MATCHED_descriptive_summary.csv")
)

alpha_plots <- list()

for (metric in alpha_metrics) {
  metric_result <- alpha_contrast_results %>% filter(Metric == metric)
  p_text <- paste0("Matched model: ", p_to_label(metric_result$P_raw))

  p <- ggplot(
    alpha_df,
    aes(x = Group, y = .data[[metric]])
  )

  if (isTRUE(connect_pairs_in_alpha_plots)) {
    p <- p +
      geom_line(
        aes(group = PairID),
        colour = "grey70",
        linewidth = 0.35,
        alpha = 0.8
      )
  }

  p <- p +
    geom_boxplot(
      aes(fill = Group),
      width = 0.48,
      outlier.shape = NA,
      alpha = 0.34,
      colour = "black",
      linewidth = 0.35
    ) +
    geom_point(
      aes(fill = Group),
      shape = 21,
      size = 2.1,
      stroke = 0.35,
      colour = "black"
    ) +
    scale_fill_manual(values = group_colours, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    labs(
      title = alpha_titles[[metric]],
      subtitle = p_text,
      x = NULL,
      y = alpha_y_labels[[metric]]
    ) +
    nature_theme(7) +
    theme(
      legend.position = "none",
      plot.subtitle = element_text(hjust = 0.5, size = 6.5)
    )

  alpha_plots[[metric]] <- p

  save_plot_set(
    p,
    subdirs$alpha,
    paste0("Alpha_", metric, "_EXACT_MATCHED"),
    width_mm = 62,
    height_mm = 78
  )

  diagnostic_path <- file.path(
    subdirs$alpha,
    paste0("Alpha_", metric, "_EXACT_MATCHED_diagnostics.pdf")
  )

  grDevices::pdf(diagnostic_path, width = 6.5, height = 3.2)
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2))

  plot(
    fitted(alpha_fits[[metric]]$fit),
    resid(alpha_fits[[metric]]$fit),
    xlab = "Fitted values",
    ylab = "Residuals",
    main = paste(metric, "residuals")
  )
  abline(h = 0, lty = 2)

  qqnorm(
    resid(alpha_fits[[metric]]$fit),
    main = paste(metric, "normal Q-Q plot")
  )
  qqline(resid(alpha_fits[[metric]]$fit), lty = 2)

  par(old_par)
  grDevices::dev.off()
}

alpha_combined <- gridExtra::arrangeGrob(
  grobs = alpha_plots[alpha_metrics],
  nrow = 1
)

save_plot_set(
  alpha_combined,
  subdirs$combined,
  "Figure_EXACT_MATCHED_alpha_diversity_panel",
  width_mm = 183,
  height_mm = 78
)

}
subdirs$alpha <- alpha_directory
subdirs$combined <- alpha_figure_directory

make_beta_plot <- function(
    coordinates,
    variance_explained,
    plot_title,
    statistics_text,
    sample_text,
    connect_pairs = FALSE
) {
  coordinates$Group <- factor(coordinates$Group, levels = group_levels)

  p <- ggplot(
    coordinates,
    aes(x = PCoA1, y = PCoA2)
  )

  if (isTRUE(connect_pairs)) {
    p <- p +
      geom_line(
        aes(group = PairID),
        colour = "grey75",
        linewidth = 0.35,
        alpha = 0.8
      )
  }

  group_sizes <- table(coordinates$Group)

  if (length(group_sizes) == 2L && all(group_sizes >= 3L)) {
    p <- p +
      stat_ellipse(
        aes(colour = Group, fill = Group),
        geom = "polygon",
        type = "norm",
        level = ellipse_level,
        alpha = 0.14,
        linewidth = 0.55,
        show.legend = FALSE
      )
  }

  p +
    geom_point(
      aes(fill = Group),
      shape = 21,
      size = 2.35,
      stroke = 0.35,
      colour = "black",
      alpha = 0.95
    ) +
    scale_fill_manual(values = group_colours, drop = FALSE) +
    scale_colour_manual(values = group_colours, drop = FALSE) +
    labs(
      title = plot_title,
      subtitle = sample_text,
      caption = statistics_text,
      x = paste0("PCoA1 (", sprintf("%.1f", variance_explained[1]), "%)"),
      y = paste0("PCoA2 (", sprintf("%.1f", variance_explained[2]), "%)"),
      fill = NULL
    ) +
    coord_equal(clip = "off") +
    nature_theme(7) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 6.2,
        margin = margin(b = 2, unit = "mm")
      ),
      plot.caption = element_text(
        hjust = 0,
        size = 6.2,
        colour = "black",
        lineheight = 1.0,
        margin = margin(t = 2, unit = "mm")
      ),
      plot.margin = margin(7, 7, 7, 7, "mm")
    )
}

run_beta_analysis <- function(
    distance_object,
    metadata_beta,
    method_label,
    output_directory,
    file_stem,
    pcoa_add = FALSE,
    n_features = NA_integer_
) {
  if (!inherits(distance_object, "dist")) {
    stop("distance_object must be a dist object.", call. = FALSE)
  }

  distance_labels <- attr(distance_object, "Labels")

  if (is.null(distance_labels)) {
    stop("The distance object has no sample labels.", call. = FALSE)
  }

  metadata_beta <- metadata_beta[
    match(distance_labels, metadata_beta$SampleID),
    ,
    drop = FALSE
  ]

  if (anyNA(metadata_beta$SampleID)) {
    stop("Distance and metadata sample alignment failed.", call. = FALSE)
  }

  if (!identical(metadata_beta$SampleID, distance_labels)) {
    stop("Distance labels are not in the same order as beta metadata.", call. = FALSE)
  }

  metadata_beta$Group <- factor(metadata_beta$Group, levels = group_levels)
  metadata_beta$PairID <- factor(metadata_beta$PairID)

  permanova_control <- make_blocked_permutation_control(
    metadata_beta$Patient,
    n_permutations
  )

  set.seed(1)
  permanova <- vegan::adonis2(
    distance_object ~ Group,
    data = metadata_beta,
    permutations = permanova_control,
    by = "terms"
  )

  dispersion <- vegan::betadisper(
    distance_object,
    group = metadata_beta$Group,
    type = "median",
    bias.adjust = FALSE
  )

  dispersion_control <- n_permutations
  set.seed(1)

  dispersion_test <- vegan::permutest(
    dispersion,
    permutations = dispersion_control
  )

  ordination <- vegan::wcmdscale(
    distance_object,
    k = 2,
    eig = TRUE,
    add = pcoa_add,
    x.ret = TRUE
  )

  coordinates <- as.data.frame(ordination$points)
  colnames(coordinates)[1:2] <- c("PCoA1", "PCoA2")
  coordinates$SampleID <- rownames(coordinates)

  coordinates <- coordinates %>%
    left_join(
      metadata_beta %>%
        select(SampleID, PairID, Patient, Site, Group),
      by = "SampleID"
    )

  positive_eigenvalues <- ordination$eig[ordination$eig > 0]

  variance_explained <- if (length(positive_eigenvalues) >= 2L) {
    100 * ordination$eig[1:2] / sum(positive_eigenvalues)
  } else {
    c(NA_real_, NA_real_)
  }

  permanova_tab <- as.data.frame(
    permanova,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  permanova_term_name <- if ("Group" %in% rownames(permanova_tab)) {
    "Group"
  } else if ("Model" %in% rownames(permanova_tab)) {
    "Model"
  } else {
    candidate_rows <- setdiff(
      rownames(permanova_tab),
      c("Residual", "Total")
    )

    if (length(candidate_rows) != 1L) {
      safe_write_csv(
        permanova_tab %>% rownames_to_column("Term"),
        file.path(
          output_directory,
          paste0(file_stem, "_PERMANOVA_EXTRACTION_FAILURE.csv")
        )
      )

      stop(
        "Could not identify one PERMANOVA model row for ",
        method_label,
        ". Inspect ",
        paste0(file_stem, "_PERMANOVA_EXTRACTION_FAILURE.csv"),
        ".",
        call. = FALSE
      )
    }

    candidate_rows[1L]
  }

  permanova_row <- permanova_tab[
    permanova_term_name,
    ,
    drop = FALSE
  ]

  permanova_f_col <- grep(
    "^F$|^F value$|^F\\.value$",
    colnames(permanova_tab),
    value = TRUE
  )

  permanova_p_col <- grep(
    "^Pr\\(>F\\)$",
    colnames(permanova_tab),
    value = TRUE
  )

  permanova_ss_col <- grep(
    "^SumOfSqs$|^SumOfSqs\\.$",
    colnames(permanova_tab),
    value = TRUE
  )

  if (
    !"Df" %in% colnames(permanova_tab) ||
    !"R2" %in% colnames(permanova_tab) ||
    length(permanova_f_col) == 0L ||
    length(permanova_p_col) == 0L ||
    length(permanova_ss_col) == 0L
  ) {
    safe_write_csv(
      permanova_tab %>% rownames_to_column("Term"),
      file.path(
        output_directory,
        paste0(file_stem, "_PERMANOVA_COLUMN_FAILURE.csv")
      )
    )

    stop(
      "The expected PERMANOVA columns were not found for ",
      method_label,
      ". Inspect ",
      paste0(file_stem, "_PERMANOVA_COLUMN_FAILURE.csv"),
      ".",
      call. = FALSE
    )
  }

  permanova_df <- as.numeric(permanova_row[1L, "Df"])
  permanova_sumsq <- as.numeric(
    permanova_row[1L, permanova_ss_col[1L]]
  )
  permanova_r2 <- as.numeric(permanova_row[1L, "R2"])
  permanova_f <- as.numeric(
    permanova_row[1L, permanova_f_col[1L]]
  )
  permanova_p <- as.numeric(
    permanova_row[1L, permanova_p_col[1L]]
  )

  if (
    any(!is.finite(c(
      permanova_df,
      permanova_sumsq,
      permanova_r2,
      permanova_f,
      permanova_p
    )))
  ) {
    safe_write_csv(
      permanova_tab %>% rownames_to_column("Term"),
      file.path(
        output_directory,
        paste0(file_stem, "_PERMANOVA_NONFINITE_FAILURE.csv")
      )
    )

    stop(
      "PERMANOVA returned a missing or non-finite statistic for ",
      method_label,
      ". Inspect ",
      paste0(file_stem, "_PERMANOVA_NONFINITE_FAILURE.csv"),
      ".",
      call. = FALSE
    )
  }

  dispersion_tab <- as.data.frame(
    dispersion_test$tab,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  dispersion_f_col <- grep(
    "^F$|^F value$|^F\\.value$",
    colnames(dispersion_tab),
    value = TRUE
  )

  dispersion_p_col <- grep(
    "^Pr\\(>F\\)$",
    colnames(dispersion_tab),
    value = TRUE
  )

  dispersion_f <- if (length(dispersion_f_col) > 0L) {
    as.numeric(dispersion_tab[1L, dispersion_f_col[1L]])
  } else {
    NA_real_
  }

  dispersion_p <- if (length(dispersion_p_col) > 0L) {
    as.numeric(dispersion_tab[1L, dispersion_p_col[1L]])
  } else {
    NA_real_
  }

  summary_df <- data.frame(
    Method = method_label,
    Comparison = "SSL_HP vs H_ADJ_SSL_HP — exact Patient__Site pairs",
    N_pairs = n_distinct(metadata_beta$PairID),
    N_samples = nrow(metadata_beta),
    N_features = as.integer(n_features),
    PERMANOVA_term = permanova_term_name,
    PERMANOVA_Df = permanova_df,
    PERMANOVA_SumOfSqs = permanova_sumsq,
    PERMANOVA_R2 = permanova_r2,
    PERMANOVA_F = permanova_f,
    PERMANOVA_P = permanova_p,
    Betadisper_F = dispersion_f,
    Betadisper_P = dispersion_p,
    Requested_permutations = n_permutations,
    Permutation_blocks = "Patient",
    Dispersion_permutations = "Unrestricted",
    Seed = 1L,
    N_patients = n_distinct(metadata_beta$Patient),
    stringsAsFactors = FALSE
  )

  statistics_text <- paste0(
    "R² = ", sprintf("%.3f", summary_df$PERMANOVA_R2),
    "; pseudo-F = ", sprintf("%.3f", summary_df$PERMANOVA_F),
    "; ", format_p_value(summary_df$PERMANOVA_P, "PERMANOVA P"), "\n",
    "betadisper F = ",
    ifelse(
      is.finite(summary_df$Betadisper_F),
      sprintf("%.3f", summary_df$Betadisper_F),
      "NA"
    ),
    "; ", format_p_value(summary_df$Betadisper_P, "P")
  )

  sample_text <- paste0(
    "Exact Patient__Site pairs: n = ",
    n_distinct(metadata_beta$PairID),
    "; patients: ", n_distinct(metadata_beta$Patient), " (", nrow(metadata_beta), " biopsies)"
  )

  beta_plot <- make_beta_plot(
    coordinates = coordinates,
    variance_explained = variance_explained,
    plot_title = method_label,
    statistics_text = statistics_text,
    sample_text = sample_text,
    connect_pairs = connect_pairs_in_beta_plots
  )

  safe_write_csv(
    summary_df,
    file.path(output_directory, paste0(file_stem, "_statistics.csv"))
  )

  safe_write_csv(
    permanova_tab %>% rownames_to_column("Term"),
    file.path(output_directory, paste0(file_stem, "_PERMANOVA_full_table.csv"))
  )

  safe_write_csv(
    dispersion_tab %>% rownames_to_column("Term"),
    file.path(output_directory, paste0(file_stem, "_betadisper_full_table.csv"))
  )

  safe_write_csv(
    coordinates,
    file.path(output_directory, paste0(file_stem, "_PCoA_coordinates.csv"))
  )

  safe_write_csv(
    data.frame(
      Axis = seq_along(ordination$eig),
      Eigenvalue = as.numeric(ordination$eig),
      stringsAsFactors = FALSE
    ),
    file.path(output_directory, paste0(file_stem, "_PCoA_eigenvalues.csv"))
  )

  saveRDS(
    list(
      distance = distance_object,
      permanova = permanova,
      betadisper = dispersion,
      betadisper_permutation_test = dispersion_test,
      pcoa = ordination
    ),
    file.path(output_directory, paste0(file_stem, "_full_objects.rds"))
  )

  save_plot_set(
    beta_plot,
    output_directory,
    file_stem,
    width_mm = 112,
    height_mm = 101
  )

  list(
    summary = summary_df,
    coordinates = coordinates,
    plot = beta_plot,
    permanova = permanova,
    dispersion = dispersion,
    dispersion_test = dispersion_test,
    ordination = ordination
  )
}

asv_counts_samples_rows <- as_samples_rows_matrix(phy_matched)

if (any(rowSums(asv_counts_samples_rows) <= 0)) {
  stop("At least one matched sample has zero ASV counts.", call. = FALSE)
}

aitchison_relative <- sweep(
  asv_counts_samples_rows + clr_pseudocount,
  1,
  rowSums(asv_counts_samples_rows + clr_pseudocount),
  "/"
)

aitchison_log <- log(aitchison_relative)
aitchison_clr <- aitchison_log - rowMeans(aitchison_log)
aitchison_distance <- dist(aitchison_clr, method = "euclidean")

metadata_beta <- matched_manifest %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  select(SampleID, PairID, Patient, Site, Group)

aitchison_result <- run_beta_analysis(
  distance_object = aitchison_distance,
  metadata_beta = metadata_beta,
  method_label = "Aitchison distance",
  output_directory = subdirs$aitchison,
  file_stem = "PCoA_Aitchison_EXACT_MATCHED",
  pcoa_add = FALSE,
  n_features = ncol(aitchison_clr)
)

safe_write_csv(
  cbind(
    SampleID = rownames(aitchison_clr),
    as.data.frame(aitchison_clr, check.names = FALSE)
  ),
  file.path(subdirs$aitchison, "Aitchison_CLR_matrix_EXACT_MATCHED.csv")
)

bray_relative <- sweep(
  asv_counts_samples_rows,
  1,
  rowSums(asv_counts_samples_rows),
  "/"
)

bray_distance <- vegan::vegdist(
  bray_relative,
  method = "bray",
  binary = FALSE
)

bray_result <- run_beta_analysis(
  distance_object = bray_distance,
  metadata_beta = metadata_beta,
  method_label = "Bray-Curtis distance",
  output_directory = subdirs$bray,
  file_stem = "PCoA_Bray_Curtis_EXACT_MATCHED",
  pcoa_add = "lingoes",
  n_features = ncol(bray_relative)
)

safe_write_csv(
  cbind(
    SampleID = rownames(bray_relative),
    as.data.frame(bray_relative, check.names = FALSE)
  ),
  file.path(subdirs$bray, "Bray_relative_abundance_matrix_EXACT_MATCHED.csv")
)

beta_summary <- bind_rows(
  aitchison_result$summary,
  bray_result$summary
)

safe_write_csv(
  beta_summary,
  file.path(subdirs$combined, "Beta_diversity_EXACT_MATCHED_summary.csv")
)

beta_combined <- gridExtra::arrangeGrob(
  aitchison_result$plot,
  bray_result$plot,
  nrow = 1
)

save_plot_set(
  beta_combined,
  subdirs$combined,
  "Figure_EXACT_MATCHED_Aitchison_and_Bray_Curtis",
  width_mm = 225,
  height_mm = 103
)

clean_species_value <- function(x) {
  x <- clean_character(x)
  x <- gsub("^s__", "", x, ignore.case = TRUE)
  x <- gsub("^Species:", "", x, ignore.case = TRUE)
  x <- gsub("_", " ", x, fixed = TRUE)
  x <- gsub("\\s+", " ", x, perl = TRUE)
  x <- trimws(x)

  bad_exact <- is.na(x) |
    toupper(x) %in% c(
      "NA", "N/A", "UNKNOWN", "UNCLASSIFIED", "UNCULTURED",
      "UNIDENTIFIED", "UNASSIGNED", "BACTERIUM", "METAGENOME",
      "SP", "SP."
    )

  bad_pattern <- grepl(
    "uncultured|unclassified|unidentified|unknown|unassigned|metagenome|environmental sample|bacterium$|\\bsp\\.?$|\\bspp\\.?$",
    x,
    ignore.case = TRUE,
    perl = TRUE
  )

  x[bad_exact | bad_pattern] <- NA_character_
  x
}

make_lineage_aware_species <- function(tax_df) {
  tax_df <- as.data.frame(
    tax_df,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  required_ranks <- c(
    "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
  )

  for (rank in required_ranks) {
    if (!rank %in% colnames(tax_df)) tax_df[[rank]] <- NA_character_
  }

  for (rank in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) {
    tax_df[[rank]] <- clean_taxonomy_value(tax_df[[rank]])
  }

  tax_df$Species <- clean_species_value(tax_df$Species)
  valid_species <- !is.na(tax_df$Species)

  species_base <- rep(NA_character_, nrow(tax_df))

  for (i in which(valid_species)) {
    genus_i <- tax_df$Genus[i]
    species_i <- tax_df$Species[i]

    already_binomial <- !is.na(genus_i) &&
      startsWith(
        tolower(species_i),
        paste0(tolower(genus_i), " ")
      )

    species_base[i] <- if (!is.na(genus_i) && !already_binomial) {
      paste(genus_i, species_i)
    } else {
      species_i
    }
  }

  lineage_key <- apply(
    tax_df[, c("Phylum", "Class", "Order", "Family", "Genus"), drop = FALSE],
    1,
    function(z) paste(ifelse(is.na(z), "NA", z), collapse = "|")
  )

  species_label <- species_base

  for (this_species in unique(species_base[!is.na(species_base)])) {
    idx <- which(species_base == this_species)
    distinct_lineages <- sort(unique(lineage_key[idx]))

    if (length(distinct_lineages) > 1L) {
      lineage_number <- match(lineage_key[idx], distinct_lineages)
      species_label[idx] <- paste0(this_species, "__lineage", lineage_number)
    }
  }

  mapping <- data.frame(
    Original_taxon = rownames(tax_df),
    Original_genus = tax_df$Genus,
    Original_species = tax_df$Species,
    Lineage_key = lineage_key,
    Species_for_analysis = species_label,
    Species_annotation_resolved = valid_species,
    stringsAsFactors = FALSE
  )

  list(
    species = species_label,
    valid = valid_species,
    mapping = mapping
  )
}

build_taxonomic_phyloseq <- function(
    phy,
    rank = c("Genus", "Species"),
    mapping_path
) {
  rank <- match.arg(rank)

  counts <- as_taxa_rows_matrix(phy)
  tax_df <- as.data.frame(
    tax_table(phy),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  tax_df <- tax_df[rownames(counts), , drop = FALSE]

  if (rank == "Genus") {
    if (!"Genus" %in% colnames(tax_df)) {
      stop("The taxonomy table has no Genus column.", call. = FALSE)
    }

    label_result <- make_lineage_aware_genus(tax_df)
    labels <- label_result$genus
    keep <- !is.na(labels) & labels != ""
    mapping <- label_result$mapping
  } else {
    species_col <- first_existing_column(
      colnames(tax_df),
      c("Species", "species", "SPECIES", "Rank7", "rank7", "s")
    )

    if (is.null(species_col)) {
      safe_write_csv(
        data.frame(
          Rank = "Species",
          Available = FALSE,
          Reason = "No Species column was present in the taxonomy table",
          stringsAsFactors = FALSE
        ),
        mapping_path
      )

      return(list(
        available = FALSE,
        reason = "No Species column was present in the taxonomy table",
        rank = rank,
        phy = NULL,
        n_features = 0L,
        n_input_asvs = nrow(counts),
        n_retained_asvs = 0L
      ))
    }

    if (species_col != "Species") {
      tax_df$Species <- tax_df[[species_col]]
    }

    label_result <- make_lineage_aware_species(tax_df)
    labels <- label_result$species
    keep <- label_result$valid & !is.na(labels) & labels != ""
    mapping <- label_result$mapping
  }

  safe_write_csv(mapping, mapping_path)

  if (sum(keep) < 2L) {
    return(list(
      available = FALSE,
      reason = paste0(
        "Fewer than two ASVs had usable ", tolower(rank), " labels"
      ),
      rank = rank,
      phy = NULL,
      n_features = 0L,
      n_input_asvs = nrow(counts),
      n_retained_asvs = sum(keep)
    ))
  }

  aggregated_counts <- rowsum(
    counts[keep, , drop = FALSE],
    group = labels[keep],
    reorder = TRUE
  )

  aggregated_counts <- aggregated_counts[
    rowSums(aggregated_counts) > 0,
    ,
    drop = FALSE
  ]

  if (nrow(aggregated_counts) < 2L) {
    return(list(
      available = FALSE,
      reason = paste0(
        "Fewer than two non-zero ", tolower(rank), " features remained after aggregation"
      ),
      rank = rank,
      phy = NULL,
      n_features = nrow(aggregated_counts),
      n_input_asvs = nrow(counts),
      n_retained_asvs = sum(keep)
    ))
  }

  aggregated_counts <- round(aggregated_counts)
  storage.mode(aggregated_counts) <- "integer"

  sample_meta <- data.frame(
    sample_data(phy),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  sample_meta <- sample_meta[colnames(aggregated_counts), , drop = FALSE]

  rank_taxonomy <- matrix(
    rownames(aggregated_counts),
    ncol = 1L,
    dimnames = list(
      rownames(aggregated_counts),
      rank
    )
  )

  rank_phy <- phyloseq::phyloseq(
    phyloseq::otu_table(aggregated_counts, taxa_are_rows = TRUE),
    phyloseq::tax_table(rank_taxonomy),
    phyloseq::sample_data(sample_meta)
  )

  list(
    available = TRUE,
    reason = "Available",
    rank = rank,
    phy = rank_phy,
    n_features = ntaxa(rank_phy),
    n_input_asvs = nrow(counts),
    n_retained_asvs = sum(keep)
  )
}

genus_bundle <- build_taxonomic_phyloseq(
  phy_matched,
  rank = "Genus",
  mapping_path = file.path(
    subdirs$audit,
    "09_lineage_aware_genus_mapping.csv"
  )
)

species_bundle <- list(available = FALSE, rank = "Species", phy = NULL,
  reason = "Species are analysed by the separate checked species workflow; this matched script runs genus models only.",
  n_features = 0L, n_input_asvs = ntaxa(phy_matched), n_retained_asvs = 0L)

if (!genus_bundle$available) {
  stop("The genus-level count object could not be built: ", genus_bundle$reason)
}

saveRDS(
  genus_bundle$phy,
  file.path(subdirs$audit, "phyloseq_EXACT_MATCHED_genus_raw_counts.rds")
)

if (species_bundle$available) {
  saveRDS(
    species_bundle$phy,
    file.path(subdirs$audit, "phyloseq_EXACT_MATCHED_resolved_species_raw_counts.rds")
  )
} else {
  writeLines(
    species_bundle$reason,
    file.path(subdirs$audit, "SPECIES_LEVEL_ANALYSIS_NOT_AVAILABLE.txt")
  )
  message(
    "Species-level DESeq2 and ANCOM-BC2 are handled separately: ",
    species_bundle$reason
  )
}

rank_input_summary <- bind_rows(
  data.frame(
    Rank = "Genus",
    Available = genus_bundle$available,
    N_input_ASVs = genus_bundle$n_input_asvs,
    N_ASVs_with_usable_rank_label = genus_bundle$n_retained_asvs,
    N_aggregated_features = genus_bundle$n_features,
    Reason = genus_bundle$reason,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Rank = "Resolved species",
    Available = species_bundle$available,
    N_input_ASVs = species_bundle$n_input_asvs,
    N_ASVs_with_usable_rank_label = species_bundle$n_retained_asvs,
    N_aggregated_features = species_bundle$n_features,
    Reason = species_bundle$reason,
    stringsAsFactors = FALSE
  )
)

safe_write_csv(
  rank_input_summary,
  file.path(subdirs$audit, "11_taxonomic_rank_input_summary.csv")
)

genus_beta_counts <- as_samples_rows_matrix(genus_bundle$phy)
genus_beta_relative <- sweep(genus_beta_counts, 1, rowSums(genus_beta_counts), "/")
genus_bray_directory <- file.path(subdirs$bray, "Genus")
dir.create(genus_bray_directory, recursive = TRUE, showWarnings = FALSE)
genus_bray_result <- run_beta_analysis(
  vegan::vegdist(genus_beta_relative, method = "bray", binary = FALSE), metadata_beta,
  "Genus Bray-Curtis distance", genus_bray_directory, "PCoA_Genus_Bray_Curtis_EXACT_MATCHED",
  pcoa_add = "lingoes", n_features = ncol(genus_beta_relative))
safe_write_csv(genus_bray_result$summary, file.path(subdirs$combined, "Genus_Bray_Curtis_EXACT_MATCHED_summary.csv"))

target_asv_counts <- as_samples_rows_matrix(phy_all)
target_asv_counts <- target_asv_counts[
  sample_names(phy_all),
  ,
  drop = FALSE
]

target_taxonomy <- as.data.frame(
  tax_table(phy_all),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

target_taxonomy <- target_taxonomy[
  colnames(target_asv_counts),
  ,
  drop = FALSE
]

if (!identical(rownames(target_taxonomy), colnames(target_asv_counts))) {
  stop("Targeted CLR taxonomy/count alignment failed.", call. = FALSE)
}

if (!"Genus" %in% colnames(target_taxonomy)) {
  stop("The taxonomy table has no Genus column for targeted CLR.", call. = FALSE)
}

plain_genus <- clean_plain_genus(target_taxonomy$Genus)
taxonomy_text <- apply(
  target_taxonomy,
  1,
  function(z) paste(z, collapse = ";")
)

keep_target_asv <- !is.na(plain_genus) &
  !grepl("CHLOROPLAST|MITOCHONDR", toupper(taxonomy_text))

if (sum(keep_target_asv) < 2L) {
  stop("Fewer than two classified bacterial genera remain for targeted CLR.")
}

plain_genus_counts <- t(
  rowsum(
    t(target_asv_counts[, keep_target_asv, drop = FALSE]),
    group = plain_genus[keep_target_asv],
    reorder = TRUE
  )
)

plain_genus_counts <- plain_genus_counts[
  sample_names(phy_all),
  colSums(plain_genus_counts) > 0,
  drop = FALSE
]

missing_target_genera <- setdiff(target_genera, colnames(plain_genus_counts))

if (length(missing_target_genera) > 0L) {
  stop(
    "The following targeted genera were not found after genus aggregation: ",
    paste(missing_target_genera, collapse = ", "),
    call. = FALSE
  )
}

target_log_counts <- log(plain_genus_counts + clr_pseudocount)
target_clr_all <- sweep(
  target_log_counts,
  1,
  rowMeans(target_log_counts),
  FUN = "-"
)

target_clr_all <- target_clr_all[sample_names(phy_matched), , drop = FALSE]
plain_genus_counts <- plain_genus_counts[sample_names(phy_matched), , drop = FALSE]

target_meta <- matched_manifest %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  select(SampleID, PairID, Patient, Site, Group)

target_clr_long <- as.data.frame(
  target_clr_all[, target_genera, drop = FALSE],
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(
    cols = all_of(target_genera),
    names_to = "Genus",
    values_to = "CLR_abundance"
  ) %>%
  left_join(target_meta, by = "SampleID") %>%
  mutate(
    PairID = factor(PairID),
    Group = factor(Group, levels = group_levels),
    Group_num = as.integer(Group),
    Genus = factor(Genus, levels = target_genera)
  ) %>%
  arrange(Genus, PairID, Group)

target_raw_long <- as.data.frame(
  plain_genus_counts[, target_genera, drop = FALSE],
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(
    cols = all_of(target_genera),
    names_to = "Genus",
    values_to = "Raw_count"
  ) %>%
  left_join(target_meta, by = "SampleID") %>%
  mutate(
    Detected = as.integer(Raw_count > 0),
    PairID = factor(PairID),
    Group = factor(Group, levels = group_levels),
    Group_num = as.integer(Group),
    Genus = factor(Genus, levels = target_genera)
  ) %>%
  arrange(Genus, PairID, Group)

if (anyNA(target_clr_long$PairID) || anyNA(target_raw_long$PairID)) {
  stop("Targeted CLR/prevalence metadata alignment failed.", call. = FALSE)
}

target_pair_check <- target_clr_long %>%
  count(Genus, PairID, Group, name = "N")

if (any(target_pair_check$N != 1L)) {
  safe_write_csv(
    target_pair_check,
    file.path(subdirs$target_clr, "Target_exact_pair_validation_failure.csv")
  )
  stop("Targeted data do not contain exactly one value per Genus x PairID x Group.")
}

safe_write_csv(target_clr_long, file.path(subdirs$target_clr, "Target_genera_biopsy_CLR_values.csv"))
target_clr_long <- target_clr_long %>%
  group_by(Genus, Patient, Group) %>%
  summarise(CLR_abundance = mean(CLR_abundance), N_biopsies = n(), .groups = "drop") %>%
  mutate(PairID = factor(Patient), ProfileID = paste(Patient, Group, sep = "__"),
         Group = factor(Group, levels = group_levels), Group_num = as.integer(Group))
target_patient_detection <- target_raw_long %>%
  group_by(Genus, Patient, Group) %>%
  summarise(Detected = as.integer(any(Detected > 0)), .groups = "drop") %>%
  mutate(ProfileID = paste(Patient, Group, sep = "__"), Group = factor(Group, levels = group_levels))

safe_write_csv(
  target_clr_long,
  file.path(subdirs$target_clr, "Target_genera_EXACT_MATCHED_CLR_values.csv")
)

safe_write_csv(
  target_raw_long,
  file.path(
    subdirs$target_clr,
    "Target_genera_EXACT_MATCHED_raw_counts_and_detection.csv"
  )
)

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
  stopifnot(nrow(d) == 20L, nlevels(d$Patient) == 10L,
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

target_clr_stats <- bind_rows(lapply(target_genera, function(genus_name) {
  d <- target_clr_long[target_clr_long$Genus == genus_name, , drop = FALSE]
  fit <- full_group_model(d, "CLR_abundance", subdirs$target_clr, genus_name)
  data.frame(Genus = genus_name, N_patients = n_distinct(d$Patient), N_profiles = nrow(d),
    Test = fit$Test, Status = fit$Status, Estimate = fit$Estimate, SE = fit$SE, Df = fit$Df, P_raw = fit$P)
})) %>% mutate(
  FDR_BH_across_3_genera = p.adjust(ifelse(Status == "COMPLETED", P_raw, NA_real_), method = "BH", n = 3L),
  P_Holm_across_3_genera = p.adjust(ifelse(Status == "COMPLETED", P_raw, NA_real_), method = "holm", n = 3L))

target_detection_stats <- bind_rows(lapply(target_genera, function(genus_name) {
  d <- target_patient_detection[target_patient_detection$Genus == genus_name, , drop = FALSE]
  fit <- full_group_model(d, "Detected", subdirs$target_clr, genus_name)
  data.frame(Genus = genus_name, N_patients = n_distinct(d$Patient), N_profiles = nrow(d),
    N_detected_SSL_HP = sum(d$Detected[d$Group == "SSL_HP"]),
    N_detected_H_ADJ_SSL_HP = sum(d$Detected[d$Group == "H_ADJ_SSL_HP"]),
    Test = fit$Test, Status = fit$Status, Odds_ratio = fit$OddsRatio, P_raw = fit$P)
})) %>% mutate(
  FDR_BH_across_3_genera = p.adjust(ifelse(Status == "COMPLETED", P_raw, NA_real_), method = "BH", n = 3L),
  P_Holm_across_3_genera = p.adjust(ifelse(Status == "COMPLETED", P_raw, NA_real_), method = "holm", n = 3L))
safe_write_csv(target_patient_detection, file.path(subdirs$target_clr, "Target_patient_detection_values.csv"))

safe_write_csv(
  target_clr_stats,
  file.path(subdirs$target_clr, "Target_genera_EXACT_MATCHED_CLR_paired_statistics.csv")
)

safe_write_csv(
  target_detection_stats,
  file.path(
    subdirs$target_clr,
    "Target_genera_EXACT_MATCHED_prevalence_paired_statistics.csv"
  )
)

target_clr_y <- target_clr_long %>%
  group_by(Genus) %>%
  summarise(
    Y_min = min(CLR_abundance, na.rm = TRUE),
    Y_max = max(CLR_abundance, na.rm = TRUE),
    Y_range = Y_max - Y_min,
    .groups = "drop"
  ) %>%
  mutate(
    Y_range = ifelse(!is.finite(Y_range) | Y_range <= 0, 1, Y_range),
    Y_bracket = Y_max + 0.12 * Y_range,
    Y_tick = Y_bracket - 0.035 * Y_range,
    Y_text = Y_bracket + 0.045 * Y_range
  )

target_clr_annotations <- target_clr_stats %>%
  mutate(
    Genus = factor(Genus, levels = target_genera),
    Label = vapply(
      FDR_BH_across_3_genera,
      format_adjusted_p,
      character(1),
      prefix = "BH FDR"
    )
  ) %>%
  left_join(target_clr_y, by = "Genus")

target_clr_plot <- ggplot(
  target_clr_long,
  aes(x = Group_num, y = CLR_abundance)
) +
  geom_line(
    aes(group = PairID),
    colour = "grey72",
    linewidth = 0.35,
    alpha = 0.85
  ) +
  geom_boxplot(
    aes(group = Group, fill = Group),
    width = 0.48,
    outlier.shape = NA,
    alpha = 0.30,
    colour = "black",
    linewidth = 0.35
  ) +
  geom_point(
    aes(fill = Group),
    shape = 21,
    size = 2.0,
    stroke = 0.35,
    colour = "black"
  ) +
  geom_segment(
    data = target_clr_annotations,
    aes(x = 1, xend = 2, y = Y_bracket, yend = Y_bracket),
    inherit.aes = FALSE,
    linewidth = 0.35
  ) +
  geom_segment(
    data = target_clr_annotations,
    aes(x = 1, xend = 1, y = Y_bracket, yend = Y_tick),
    inherit.aes = FALSE,
    linewidth = 0.35
  ) +
  geom_segment(
    data = target_clr_annotations,
    aes(x = 2, xend = 2, y = Y_bracket, yend = Y_tick),
    inherit.aes = FALSE,
    linewidth = 0.35
  ) +
  geom_text(
    data = target_clr_annotations,
    aes(x = 1.5, y = Y_text, label = Label),
    inherit.aes = FALSE,
    family = font_family,
    size = 2.35
  ) +
  facet_wrap(~ Genus, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = group_colours, limits = group_levels, drop = FALSE) +
  scale_x_continuous(
    breaks = seq_along(group_levels),
    labels = group_levels,
    limits = c(0.5, length(group_levels) + 0.5)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.22))) +
  labs(
    title = "A",
    x = NULL,
    y = "CLR abundance"
  ) +
  nature_theme(7) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "italic", size = 7),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
  )

target_detection_summary <- target_raw_long %>%
  group_by(Genus, Group, Group_num) %>%
  summarise(
    N_pairs = n(),
    N_detected = sum(Detected),
    Prevalence_percent = 100 * mean(Detected),
    .groups = "drop"
  ) %>%
  mutate(
    Detection_label = paste0(N_detected, "/", N_pairs),
    Label_y = pmin(Prevalence_percent + 4, 102)
  )

safe_write_csv(
  target_detection_summary,
  file.path(
    subdirs$target_clr,
    "Target_genera_EXACT_MATCHED_prevalence_summary.csv"
  )
)

target_prevalence_plot <- ggplot(
  target_detection_summary,
  aes(x = Group_num, y = Prevalence_percent, fill = Group)
) +
  geom_col(width = 0.62, colour = "black", linewidth = 0.25) +
  geom_text(
    aes(y = Label_y, label = Detection_label),
    family = font_family,
    size = 2.35,
    colour = "black"
  ) +
  facet_wrap(~ Genus, nrow = 1) +
  scale_fill_manual(values = group_colours, limits = group_levels, drop = FALSE) +
  scale_x_continuous(
    breaks = seq_along(group_levels),
    labels = group_levels,
    limits = c(0.5, length(group_levels) + 0.5)
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = c(0, 25, 50, 75, 100),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +
  labs(
    title = "B",
    x = NULL,
    y = "Detection in matched biopsies"
  ) +
  nature_theme(7) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
  )

save_plot_set(
  target_clr_plot,
  subdirs$target_clr,
  "Target_genera_EXACT_MATCHED_CLR_paired",
  width_mm = 210,
  height_mm = 88
)

save_plot_set(
  target_prevalence_plot,
  subdirs$target_clr,
  "Target_genera_EXACT_MATCHED_prevalence_paired",
  width_mm = 210,
  height_mm = 72
)

g_target_clr <- ggplotGrob(target_clr_plot)
g_target_prev <- ggplotGrob(target_prevalence_plot)
max_target_widths <- grid::unit.pmax(g_target_clr$widths, g_target_prev$widths)
g_target_clr$widths <- max_target_widths
g_target_prev$widths <- max_target_widths

target_combined <- gridExtra::arrangeGrob(
  g_target_clr,
  g_target_prev,
  ncol = 1,
  heights = c(1.1, 0.88)
)

save_plot_set(
  target_combined,
  subdirs$combined,
  "Figure_EXACT_MATCHED_target_CLR_and_prevalence",
  width_mm = 210,
  height_mm = 155
)

genus_counts_biopsy <- as_taxa_rows_matrix(genus_bundle$phy)
genus_meta_biopsy <- data.frame(sample_data(genus_bundle$phy), check.names = FALSE)
genus_meta_biopsy <- genus_meta_biopsy[colnames(genus_counts_biopsy), , drop = FALSE]
genus_meta_biopsy$ProfileID <- paste(genus_meta_biopsy$Patient, genus_meta_biopsy$Group, sep = "__")
genus_profiles <- genus_meta_biopsy %>% group_by(ProfileID, Patient, Group) %>%
  summarise(N_biopsies = n(), PairID = paste(as.character(PairID), collapse = ";"),
            Site = paste(as.character(Site), collapse = ";"), .groups = "drop") %>%
  arrange(Group, Patient) %>% as.data.frame()
rownames(genus_profiles) <- genus_profiles$ProfileID
genus_profiles$Patient <- factor(genus_profiles$Patient)
genus_profiles$Group <- factor(genus_profiles$Group, levels = group_levels)
genus_counts_patient <- t(rowsum(t(genus_counts_biopsy), genus_meta_biopsy$ProfileID, reorder = FALSE))
genus_counts_patient <- genus_counts_patient[, rownames(genus_profiles), drop = FALSE]
stopifnot(nrow(genus_profiles) == 20L, nlevels(genus_profiles$Patient) == 10L,
          sum(genus_counts_patient) == sum(genus_counts_biopsy),
          all(table(genus_profiles$Patient, genus_profiles$Group) == 1L))
genus_bundle$phy <- phyloseq(otu_table(genus_counts_patient, taxa_are_rows = TRUE),
  tax_table(genus_bundle$phy), sample_data(genus_profiles))
safe_write_csv(genus_profiles, file.path(subdirs$audit, "Genus_patient_group_profiles.csv"))
safe_write_csv(data.frame(Taxon = rownames(genus_counts_patient), genus_counts_patient, check.names = FALSE),
               file.path(subdirs$audit, "Genus_patient_group_counts.csv"))

wrap_taxon_label <- function(x, width = 34L) {
  vapply(
    as.character(x),
    function(z) paste(strwrap(z, width = width), collapse = "\n"),
    character(1)
  )
}

make_aligned_volcano <- function(
    result_df,
    rank_label,
    output_directory,
    filename_stem
) {
  # Unavailable adjusted P values remain in the results table, not at FDR = 1.
  plot_df <- result_df %>%
    filter(is.finite(FDR), FDR >= 0, FDR <= 1, is.finite(log2FoldChange)) %>%
    mutate(
      FDR_plot = pmax(FDR, .Machine$double.xmin),
      minus_log10_FDR = -log10(FDR_plot),
      Plot_category = case_when(
        Significant_FDR_0_05 & log2FoldChange > 0 ~ "Higher in SSL_HP",
        Significant_FDR_0_05 & log2FoldChange < 0 ~ "Higher in H_ADJ_SSL_HP",
        TRUE ~ "Not significant"
      )
    )

  label_df <- plot_df %>%
    filter(
      Significant_FDR_0_05,
      is.finite(log2FoldChange),
      is.finite(minus_log10_FDR)
    ) %>%
    mutate(
      Side = ifelse(log2FoldChange < 0, "Left", "Right")
    )

  finite_x <- plot_df$log2FoldChange[is.finite(plot_df$log2FoldChange)]
  max_abs_x <- if (length(finite_x) == 0L) 1 else max(abs(finite_x), na.rm = TRUE)
  if (!is.finite(max_abs_x) || max_abs_x <= 0) max_abs_x <- 1

  finite_y <- plot_df$minus_log10_FDR[is.finite(plot_df$minus_log10_FDR)]
  max_y <- if (length(finite_y) == 0L) 1 else max(finite_y, na.rm = TRUE)
  if (!is.finite(max_y) || max_y <= 0) max_y <- 1

  x_limit <- 2.00 * max_abs_x
  y_upper <- max(1.18 * max_y, -log10(alpha_level) + 1)
  preferred_gap <- max(0.40, y_upper / max(12, nrow(label_df) + 2L))

  if (nrow(label_df) > 0L) {
    label_df <- label_df %>%
      group_by(Side) %>%
      arrange(minus_log10_FDR, .by_group = TRUE) %>%
      mutate(
        Label_y = spread_label_positions(
          minus_log10_FDR,
          lower = 0.45,
          upper = y_upper - 0.35,
          preferred_gap = preferred_gap
        ),
        Label_x = ifelse(Side == "Left", -0.95 * x_limit, 0.95 * x_limit),
        Label_hjust = ifelse(Side == "Left", 0, 1),
        Segment_x = ifelse(Side == "Left", -1.12 * max_abs_x, 1.12 * max_abs_x),
        Display_label = wrap_taxon_label(Taxon, width = ifelse(rank_label == "Genus", 31L, 38L))
      ) %>%
      ungroup()
  }

  safe_write_csv(
    label_df,
    file.path(output_directory, paste0(filename_stem, "_label_layout.csv"))
  )

  colours <- c(
    "Higher in SSL_HP" = "#D73027",
    "Higher in H_ADJ_SSL_HP" = "#1A9850",
    "Not significant" = "grey72"
  )

  p <- ggplot(
    plot_df,
    aes(x = log2FoldChange, y = minus_log10_FDR)
  ) +
    geom_point(
      aes(colour = Plot_category),
      size = 1.70,
      alpha = 0.88,
      na.rm = TRUE
    ) +
    geom_hline(
      yintercept = -log10(alpha_level),
      linetype = 2,
      linewidth = 0.35
    ) +
    geom_vline(
      xintercept = 0,
      linetype = 2,
      linewidth = 0.35
    )

  if (nrow(label_df) > 0L) {
    p <- p +
      geom_segment(
        data = label_df,
        aes(
          x = log2FoldChange,
          y = minus_log10_FDR,
          xend = Segment_x,
          yend = Label_y
        ),
        inherit.aes = FALSE,
        linewidth = 0.28,
        colour = "grey45"
      ) +
      geom_text(
        data = label_df,
        aes(
          x = Label_x,
          y = Label_y,
          label = Display_label,
          hjust = Label_hjust
        ),
        inherit.aes = FALSE,
        family = font_family,
        size = 2.35,
        lineheight = 0.92,
        colour = "black"
      )
  }

  p <- p +
    scale_colour_manual(
      values = colours,
      breaks = c(
        "Higher in H_ADJ_SSL_HP",
        "Higher in SSL_HP",
        "Not significant"
      ),
      drop = FALSE
    ) +
    coord_cartesian(
      xlim = c(-x_limit, x_limit),
      ylim = c(0, y_upper),
      clip = "on"
    ) +
    labs(
      title = paste0("DESeq2 ", tolower(rank_label), ": SSL_HP vs H_ADJ_SSL_HP"),
      subtitle = "Exact Patient__Site pairs; positive log2 fold change indicates higher abundance in SSL_HP",
      x = expression(log[2]~fold~change),
      y = expression(-log[10]~FDR),
      colour = NULL
    ) +
    nature_theme(7) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.subtitle = element_text(hjust = 0.5, size = 6.2),
      plot.margin = margin(7, 7, 7, 7, "mm")
    )

  save_plot_set(
    p,
    output_directory,
    filename_stem,
    width_mm = if (rank_label == "Genus") 205 else 215,
    height_mm = 116
  )

  p
}

# Differential abundance on patient-group counts from exact site pairs.
run_deseq2_rank <- function(rank_bundle, output_directory) {
  if (!rank_bundle$available) return(NULL)

  rank_label <- rank_bundle$rank
  phy_rank <- rank_bundle$phy
  counts_rank <- as_taxa_rows_matrix(phy_rank)
  counts_rank <- round(counts_rank)
  storage.mode(counts_rank) <- "integer"

  col_data <- data.frame(
    sample_data(phy_rank),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  col_data <- col_data[colnames(counts_rank), , drop = FALSE]
  col_data$Patient <- factor(col_data$Patient)
  col_data$Group <- factor(col_data$Group, levels = group_levels)

  if (!identical(rownames(col_data), colnames(counts_rank))) {
    stop("DESeq2 count/metadata alignment failed for ", rank_label, ".")
  }

  if (any(table(col_data$Patient) != 2L)) {
    stop("DESeq2 requires exactly two samples per Patient for ", rank_label, ".")
  }

  design_formula <- ~ Patient + Group
  design_matrix <- model.matrix(design_formula, data = col_data)

  if (qr(design_matrix)$rank < ncol(design_matrix)) {
    stop("The DESeq2 design matrix is not full rank for ", rank_label, ".")
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_rank,
    colData = col_data,
    design = design_formula
  )

  dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]
  dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")

  dds <- tryCatch(
    DESeq2::DESeq(
      dds,
      test = "Wald",
      fitType = "parametric",
      quiet = FALSE
    ),
    error = function(e) {
      message(
        "Parametric dispersion fit failed for ", rank_label,
        "; retrying with fitType = 'local'."
      )
      DESeq2::DESeq(
        dds,
        test = "Wald",
        fitType = "local",
        quiet = FALSE
      )
    }
  )

  result <- DESeq2::results(
    dds,
    contrast = c("Group", "SSL_HP", "H_ADJ_SSL_HP"),
    alpha = alpha_level,
    independentFiltering = TRUE,
    cooksCutoff = TRUE
  )

  result_df <- as.data.frame(result) %>%
    rownames_to_column("Taxon") %>%
    mutate(
      Rank = rank_label,
      Comparison = "SSL_HP vs H_ADJ_SSL_HP",
      Direction = case_when(
        is.na(log2FoldChange) ~ NA_character_,
        log2FoldChange > 0 ~ "Higher in SSL_HP",
        log2FoldChange < 0 ~ "Higher in H_ADJ_SSL_HP",
        TRUE ~ "No direction"
      ),
      Significant_FDR_0_05 = !is.na(padj) & padj < alpha_level,
      FDR = padj,
      P_value = pvalue
    ) %>%
    arrange(FDR, P_value, desc(abs(log2FoldChange)))

  safe_write_csv(
    result_df,
    file.path(
      output_directory,
      paste0("DESeq2_EXACT_MATCHED_", tolower(rank_label), "_ALL_results.csv")
    )
  )

  safe_write_csv(
    result_df %>% filter(Significant_FDR_0_05),
    file.path(
      output_directory,
      paste0("DESeq2_EXACT_MATCHED_", tolower(rank_label), "_FDR_LT_0.05.csv")
    )
  )

  safe_write_csv(
    data.frame(
      SampleID = colnames(dds),
      Patient = as.character(SummarizedExperiment::colData(dds)$Patient),
      Group = as.character(SummarizedExperiment::colData(dds)$Group),
      Size_factor = sizeFactors(dds),
      stringsAsFactors = FALSE
    ),
    file.path(output_directory, "DESeq2_sample_manifest_and_size_factors.csv")
  )

  safe_write_csv(
    cbind(
      Taxon = rownames(dds),
      as.data.frame(DESeq2::counts(dds, normalized = TRUE), check.names = FALSE)
    ),
    file.path(output_directory, "DESeq2_normalized_counts.csv")
  )

  saveRDS(
    dds,
    file.path(output_directory, paste0("DESeq2_EXACT_MATCHED_", tolower(rank_label), "_dds.rds"))
  )

  volcano <- make_aligned_volcano(
    result_df = result_df,
    rank_label = rank_label,
    output_directory = output_directory,
    filename_stem = paste0(
      "DESeq2_EXACT_MATCHED_",
      tolower(rank_label),
      "_volcano_WIDE_FULL_LABELS"
    )
  )

  grDevices::pdf(
    file.path(output_directory, paste0("DESeq2_", tolower(rank_label), "_MA_plot.pdf")),
    width = 6,
    height = 5
  )
  DESeq2::plotMA(
    result,
    alpha = alpha_level,
    main = paste0("DESeq2 exact matched ", tolower(rank_label), " analysis")
  )
  grDevices::dev.off()

  grDevices::pdf(
    file.path(output_directory, paste0("DESeq2_", tolower(rank_label), "_dispersion_plot.pdf")),
    width = 6,
    height = 5
  )
  DESeq2::plotDispEsts(dds)
  grDevices::dev.off()

  summary_df <- data.frame(
    Analysis = paste0("DESeq2 ", rank_label, " — patient-group counts from exact site pairs"),
    Rank = rank_label,
    Formula = "~ Patient + Group",
    N_patients = nlevels(col_data$Patient),
    N_samples = nrow(col_data),
    N_features_tested = nrow(result_df),
    N_FDR_significant = sum(result_df$Significant_FDR_0_05, na.rm = TRUE),
    N_higher_SSL_HP = sum(
      result_df$Significant_FDR_0_05 & result_df$log2FoldChange > 0,
      na.rm = TRUE
    ),
    N_higher_H_ADJ_SSL_HP = sum(
      result_df$Significant_FDR_0_05 & result_df$log2FoldChange < 0,
      na.rm = TRUE
    ),
    Size_factor_method = "poscounts",
    Test = "Wald",
    Multiple_testing = "Benjamini-Hochberg FDR",
    LFC_shrinkage = FALSE,
    stringsAsFactors = FALSE
  )

  safe_write_csv(summary_df, file.path(output_directory, "DESeq2_summary.csv"))

  list(
    summary = summary_df,
    result = result_df,
    significant = result_df %>% filter(Significant_FDR_0_05),
    dds = dds,
    volcano = volcano
  )
}

deseq2_genus <- run_deseq2_rank(genus_bundle, subdirs$deseq2_genus)

deseq2_species <- if (species_bundle$available) {
  tryCatch(
    run_deseq2_rank(species_bundle, subdirs$deseq2_species),
    error = function(e) {
      writeLines(
        conditionMessage(e),
        file.path(subdirs$deseq2_species, "DESEQ2_SPECIES_ANALYSIS_FAILED.txt")
      )
      warning("Species-level DESeq2 failed: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
} else {
  NULL
}

as_ancom_result_df <- function(x) {
  if (is.null(x)) return(data.frame())

  out <- as.data.frame(
    x,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!"taxon" %in% colnames(out)) {
    out <- out %>% rownames_to_column("taxon")
  }

  out
}

as_ancom_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  !is.na(x) & toupper(as.character(x)) %in% c("TRUE", "T", "1")
}

find_ancom_column <- function(df, prefix, preferred_suffix = "GroupSSL_HP") {
  exact <- paste0(prefix, "_", preferred_suffix)
  if (exact %in% colnames(df)) return(exact)

  candidates <- grep(
    paste0("^", prefix, "_.*SSL_HP$"),
    colnames(df),
    value = TRUE
  )

  candidates <- candidates[
    !grepl("H_ADJ_SSL_HP", candidates, fixed = TRUE)
  ]

  if (length(candidates) == 0L) return(NA_character_)
  candidates[1L]
}

make_ancom_effect_plot <- function(
    tidy_df,
    rank_label,
    output_directory
) {
  significant_df <- tidy_df %>%
    filter(
      Sensitivity_stable_significant,
      is.finite(log2_abundance_difference),
      is.finite(log2_CI_low),
      is.finite(log2_CI_high)
    )

  colours <- c(
    "Higher in SSL_HP" = "#D73027",
    "Higher in H_ADJ_SSL_HP" = "#1A9850"
  )

  if (nrow(significant_df) > 0L) {
    significant_df <- significant_df %>%
      arrange(log2_abundance_difference) %>%
      mutate(
        Taxon = factor(Taxon, levels = Taxon),
        Display_taxon = wrap_taxon_label(as.character(Taxon), width = 42L)
      )

    label_map <- setNames(significant_df$Display_taxon, as.character(significant_df$Taxon))

    p <- ggplot(
      significant_df,
      aes(x = log2_abundance_difference, y = Taxon)
    ) +
      geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35) +
      geom_errorbar(
        aes(xmin = log2_CI_low, xmax = log2_CI_high),
        orientation = "y",
        width = 0,
        linewidth = 0.40
      ) +
      geom_point(aes(colour = Effect_direction), size = 2.2) +
      scale_colour_manual(values = colours, drop = FALSE) +
      scale_y_discrete(labels = label_map) +
      labs(
        title = paste0("ANCOM-BC2 ", tolower(rank_label), ": SSL_HP vs H_ADJ_SSL_HP"),
        subtitle = "BH-significant and pseudocount-sensitivity-stable features",
        x = "Bias-corrected log2 abundance difference",
        y = NULL,
        colour = NULL
      ) +
      nature_theme(7) +
      theme(
        legend.position = "bottom",
        plot.subtitle = element_text(hjust = 0.5, size = 6.2)
      )

    save_plot_set(
      p,
      output_directory,
      paste0("ANCOMBC2_EXACT_MATCHED_", tolower(rank_label), "_significant_effects"),
      width_mm = if (rank_label == "Genus") 190 else 235,
      height_mm = max(85, 8 * nrow(significant_df) + 38)
    )
  }

  top_df <- tidy_df %>%
    filter(is.finite(log2_abundance_difference), is.finite(se)) %>%
    arrange(q_value_BH, p_value, desc(abs(log2_abundance_difference))) %>%
    slice_head(n = 20L) %>%
    arrange(log2_abundance_difference) %>%
    mutate(
      Taxon = factor(Taxon, levels = Taxon),
      Display_taxon = wrap_taxon_label(as.character(Taxon), width = 42L),
      Result = ifelse(
        Sensitivity_stable_significant,
        "BH significant and sensitivity-stable",
        "Not confirmatory"
      )
    )

  if (nrow(top_df) > 0L) {
    label_map <- setNames(top_df$Display_taxon, as.character(top_df$Taxon))
    top_colours <- c(
      "BH significant and sensitivity-stable" = "#D55E00",
      "Not confirmatory" = "grey45"
    )

    p_top <- ggplot(
      top_df,
      aes(x = log2_abundance_difference, y = Taxon)
    ) +
      geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.35) +
      geom_errorbar(
        aes(xmin = log2_CI_low, xmax = log2_CI_high),
        orientation = "y",
        width = 0,
        linewidth = 0.35
      ) +
      geom_point(aes(colour = Result), size = 2.0) +
      scale_colour_manual(values = top_colours, drop = FALSE) +
      scale_y_discrete(labels = label_map) +
      labs(
        title = paste0("ANCOM-BC2 ", tolower(rank_label), " top model-based effects"),
        subtitle = "Exploratory display; inference is based on BH FDR plus pseudocount sensitivity",
        x = "Bias-corrected log2 abundance difference",
        y = NULL,
        colour = NULL
      ) +
      nature_theme(7) +
      theme(
        legend.position = "bottom",
        plot.subtitle = element_text(hjust = 0.5, size = 6.2)
      )

    save_plot_set(
      p_top,
      output_directory,
      paste0("Supplementary_ANCOMBC2_", tolower(rank_label), "_top20_effects_EXPLORATORY"),
      width_mm = if (rank_label == "Genus") 200 else 245,
      height_mm = max(120, 7 * nrow(top_df) + 38)
    )
  }
}

# ANCOM-BC2 sensitivity analysis and availability diagnostics.
run_ancombc2_rank <- function(rank_bundle, output_directory) {
  if (!rank_bundle$available) return(NULL)

  rank_label <- rank_bundle$rank
  phy_rank <- rank_bundle$phy

  sample_data(phy_rank)$Patient <- factor(sample_data(phy_rank)$Patient)
  sample_data(phy_rank)$Group <- factor(
    as.character(sample_data(phy_rank)$Group),
    levels = group_levels
  )

  check_raw_integer_counts(
    phy_rank,
    paste0("ANCOM-BC2 exact matched ", tolower(rank_label), " analysis")
  )

  manifest <- data.frame(
    sample_data(phy_rank),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  manifest$SampleID <- rownames(manifest)
  manifest <- manifest[
    match(sample_names(phy_rank), manifest$SampleID),
    c("SampleID", "Patient", "Site", "Group"),
    drop = FALSE
  ]
  manifest$Patient <- factor(manifest$Patient)
  manifest$Group <- factor(manifest$Group, levels = group_levels)

  design <- model.matrix(~ Patient + Group, data = manifest)
  if (qr(design)$rank < ncol(design)) {
    stop("The ANCOM-BC2 fixed-effect design is not full rank for ", rank_label, ".")
  }

  safe_write_csv(manifest, file.path(output_directory, "ANCOMBC2_sample_manifest.csv"))

  captured_warnings <- character(0)

  fit <- withCallingHandlers(
    ANCOMBC::ancombc2(
      data = phy_rank,
      assay_name = "counts",
      tax_level = NULL,
      fix_formula = "Patient + Group",
      rand_formula = NULL,
      p_adj_method = "BH",
      pseudo = 0,
      pseudo_sens = TRUE,
      prv_cut = 0,
      lib_cut = 0,
      s0_perc = 0.05,
      group = "Group",
      struc_zero = TRUE,
      neg_lb = FALSE,
      alpha = alpha_level,
      n_cl = ancombc_n_cores,
      verbose = TRUE,
      global = FALSE,
      pairwise = FALSE,
      dunnet = FALSE,
      trend = FALSE,
      iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
      em_control = list(tol = 1e-5, max_iter = 100)
    ),
    warning = function(w) {
      captured_warnings <<- c(captured_warnings, conditionMessage(w))
    }
  )

  saveRDS(
    fit,
    file.path(output_directory, paste0("ANCOMBC2_EXACT_MATCHED_", tolower(rank_label), "_full_output.rds"))
  )

  captured_warnings <- unique(captured_warnings)
  writeLines(
    if (length(captured_warnings) == 0L) "No warnings were emitted." else captured_warnings,
    file.path(output_directory, "ANCOMBC2_warnings.txt")
  )

  primary_raw <- as_ancom_result_df(fit$res)
  zero_raw <- as_ancom_result_df(fit$zero_ind)
  safe_write_csv(primary_raw, file.path(output_directory, "ANCOMBC2_primary_fixed_effect_results_raw.csv"))
  safe_write_csv(zero_raw, file.path(output_directory, "ANCOMBC2_structural_zeros_raw.csv"))

  lfc_col <- find_ancom_column(primary_raw, "lfc")
  se_col <- find_ancom_column(primary_raw, "se")
  w_col <- find_ancom_column(primary_raw, "W")
  p_col <- find_ancom_column(primary_raw, "p")
  q_col <- find_ancom_column(primary_raw, "q")
  diff_col <- find_ancom_column(primary_raw, "diff")
  passed_ss_col <- find_ancom_column(primary_raw, "passed_ss")

  required_cols <- c(lfc_col, se_col, w_col, p_col, q_col, diff_col, passed_ss_col)
  if (anyNA(required_cols)) {
    safe_write_csv(
      data.frame(
        Requested = c("lfc", "se", "W", "p", "q", "diff", "passed_ss"),
        Resolved_column = required_cols,
        stringsAsFactors = FALSE
      ),
      file.path(output_directory, "ANCOMBC2_column_resolution_failure.csv")
    )
    stop("Could not resolve the GroupSSL_HP coefficient columns for ", rank_label, ".")
  }

  tidy <- data.frame(
    Taxon = as.character(primary_raw$taxon),
    Rank = rank_label,
    Comparison = "SSL_HP vs H_ADJ_SSL_HP",
    lfc_natural_log = as.numeric(primary_raw[[lfc_col]]),
    se = as.numeric(primary_raw[[se_col]]),
    W = as.numeric(primary_raw[[w_col]]),
    p_value = as.numeric(primary_raw[[p_col]]),
    q_value_BH = as.numeric(primary_raw[[q_col]]),
    diff = as_ancom_flag(primary_raw[[diff_col]]),
    passed_ss = as_ancom_flag(primary_raw[[passed_ss_col]]),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      q_value_Holm = p.adjust(p_value, method = "holm"),
      Wald_CI_low = lfc_natural_log - 1.96 * se,
      Wald_CI_high = lfc_natural_log + 1.96 * se,
      log2_abundance_difference = lfc_natural_log / log(2),
      log2_CI_low = Wald_CI_low / log(2),
      log2_CI_high = Wald_CI_high / log(2),
      Sensitivity_stable_significant =
        !is.na(q_value_BH) & q_value_BH < alpha_level & passed_ss,
      Holm_sensitivity_stable =
        !is.na(q_value_Holm) & q_value_Holm < alpha_level & passed_ss,
      Effect_direction = case_when(
        is.na(lfc_natural_log) ~ NA_character_,
        lfc_natural_log > 0 ~ "Higher in SSL_HP",
        lfc_natural_log < 0 ~ "Higher in H_ADJ_SSL_HP",
        TRUE ~ "No direction"
      )
    ) %>%
    arrange(q_value_BH, p_value, desc(abs(lfc_natural_log)))

  X <- model.matrix(~ Patient + Group, manifest)
  ft <- as.matrix(fit$feature_table)[tidy$Taxon, rownames(manifest), drop = FALSE]
  e <- rep(0, ncol(X)); e[match("GroupSSL_HP", colnames(X))] <- 1
  test_info <- t(vapply(seq_len(nrow(ft)), function(i) {
    Xi <- X[ft[i, ] > 0, , drop = FALSE]
    rank <- qr(Xi)$rank
    c(df = nrow(Xi) - rank, estimable = qr(rbind(Xi, e))$rank == rank)
  }, numeric(2)))
  tidy$Residual_df <- test_info[, "df"]
  tidy$Contrast_estimable <- test_info[, "estimable"] == 1
  tidy$Test_available <- tidy$Residual_df > 0 & tidy$Contrast_estimable &
    is.finite(tidy$lfc_natural_log) & is.finite(tidy$se) & tidy$se > 0 &
    is.finite(tidy$W) & is.finite(tidy$p_value) & is.finite(tidy$q_value_BH)
  tidy$Reporting_P <- ifelse(tidy$Test_available, tidy$p_value, NA_real_)
  tidy$Reporting_BH_FDR <- ifelse(tidy$Test_available, tidy$q_value_BH, NA_real_)
  tidy$Reporting_Holm_P <- ifelse(tidy$Test_available, tidy$q_value_Holm, NA_real_)
  tidy$Reporting_log2_difference <- ifelse(tidy$Test_available, tidy$log2_abundance_difference, NA_real_)
  tidy$Sensitivity_stable_significant <- tidy$Sensitivity_stable_significant & tidy$Test_available
  tidy$Holm_sensitivity_stable <- tidy$Holm_sensitivity_stable & tidy$Test_available

  safe_write_csv(tidy, file.path(output_directory, "ANCOMBC2_tidy_all_model_based_results.csv"))
  safe_write_csv(
    tidy %>% filter(Sensitivity_stable_significant),
    file.path(output_directory, "ANCOMBC2_BH_SENSITIVITY_STABLE.csv")
  )
  safe_write_csv(
    tidy %>% filter(Holm_sensitivity_stable),
    file.path(output_directory, "ANCOMBC2_HOLM_SENSITIVITY_STABLE.csv")
  )

  structural_zero_long <- data.frame()
  if (nrow(zero_raw) > 0L) {
    zero_columns <- setdiff(colnames(zero_raw), "taxon")
    if (length(zero_columns) > 0L) {
      structural_zero_long <- zero_raw %>%
        pivot_longer(
          cols = all_of(zero_columns),
          names_to = "Group_column",
          values_to = "Structural_zero"
        ) %>%
        mutate(
          Structural_zero = as_ancom_flag(Structural_zero),
          Group = sub(
            "^structural_zero \\(Group = (.*)\\)$",
            "\\1",
            Group_column
          )
        ) %>%
        filter(Structural_zero)
    }
  }
  safe_write_csv(structural_zero_long, file.path(output_directory, "ANCOMBC2_structural_zeros.csv"))

  critical_patterns <- c(
    "NaNs produced",
    "longer object length is not a multiple",
    "failed to converge",
    "singular Hessian",
    "not positive definite",
    "bobyqa --",
    "trust region"
  )

  critical_warnings <- captured_warnings[
    vapply(
      captured_warnings,
      function(z) any(grepl(paste(critical_patterns, collapse = "|"), z, ignore.case = TRUE)),
      logical(1)
    )
  ]

  finite_p <- tidy$p_value[is.finite(tidy$p_value)]
  all_p_one <- length(finite_p) > 0L &&
    all(abs(finite_p - 1) < 1e-12) &&
    any(abs(tidy$W) > 1e-8, na.rm = TRUE)
  all_p_missing <- length(finite_p) == 0L

  run_status <- if (
    length(critical_warnings) > 0L || all_p_one || all_p_missing
  ) {
    "REVIEW_REQUIRED"
  } else if (length(captured_warnings) > 0L) {
    "COMPLETED_WITH_WARNINGS"
  } else {
    "PASS"
  }

  summary_df <- data.frame(
    Analysis = paste0("ANCOM-BC2 ", rank_label, " — patient-group counts from exact site pairs"),
    Rank = rank_label,
    Run_status = run_status,
    ANCOMBC_version = as.character(packageVersion("ANCOMBC")),
    Bias_EM_patch_applied_in_memory =
      ancombc_patch_status$Bias_EM_patch_applied_in_memory[1],
    Fixed_effect_formula = "Patient + Group",
    Random_effect_formula = "NULL",
    N_patients = n_distinct(manifest$Patient),
    N_samples = nrow(manifest),
    N_features_reported = nrow(tidy),
    N_available_tests = sum(tidy$Test_available),
    N_unavailable_tests = sum(!tidy$Test_available),
    N_BH_sensitivity_stable = sum(tidy$Sensitivity_stable_significant, na.rm = TRUE),
    N_Holm_sensitivity_stable = sum(tidy$Holm_sensitivity_stable, na.rm = TRUE),
    N_structural_zero_rows = nrow(structural_zero_long),
    Additional_prevalence_filter = FALSE,
    ANCOMBC_prv_cut = 0,
    N_warnings = length(captured_warnings),
    N_critical_warnings = length(critical_warnings),
    Degenerate_all_P_equal_1_with_nonzero_W = all_p_one,
    Degenerate_all_P_missing = all_p_missing,
    stringsAsFactors = FALSE
  )

  safe_write_csv(summary_df, file.path(output_directory, "ANCOMBC2_summary.csv"))
  writeLines(
    if (length(critical_warnings) == 0L) "No critical warnings were detected." else critical_warnings,
    file.path(output_directory, "ANCOMBC2_CRITICAL_WARNINGS.txt")
  )

  make_ancom_effect_plot(tidy, rank_label, output_directory)

  list(
    summary = summary_df,
    tidy = tidy,
    significant = tidy %>% filter(Sensitivity_stable_significant),
    structural_zero = structural_zero_long,
    fit = fit
  )
}

ancombc2_genus <- run_ancombc2_rank(genus_bundle, subdirs$ancombc2_genus)

ancombc2_species <- if (species_bundle$available) {
  tryCatch(
    run_ancombc2_rank(species_bundle, subdirs$ancombc2_species),
    error = function(e) {
      writeLines(
        conditionMessage(e),
        file.path(subdirs$ancombc2_species, "ANCOMBC2_SPECIES_ANALYSIS_FAILED.txt")
      )
      warning("Species-level ANCOM-BC2 failed: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
} else {
  NULL
}

# Cross-method reporting uses available tissue tests only.
make_cross_method_table <- function(deseq_object, ancom_object, rank_label) {
  if (is.null(deseq_object) || is.null(ancom_object)) return(NULL)

  out <- full_join(
    deseq_object$result %>%
      select(
        Taxon,
        DESeq2_log2FoldChange = log2FoldChange,
        DESeq2_P = P_value,
        DESeq2_FDR = FDR,
        DESeq2_significant = Significant_FDR_0_05
      ),
    ancom_object$tidy %>%
      select(
        Taxon,
        ANCOMBC2_test_available = Test_available,
        ANCOMBC2_log2_difference = Reporting_log2_difference,
        ANCOMBC2_P = Reporting_P,
        ANCOMBC2_BH_q = Reporting_BH_FDR,
        ANCOMBC2_Holm_q = Reporting_Holm_P,
        ANCOMBC2_passed_ss = passed_ss,
        ANCOMBC2_BH_sensitivity_stable = Sensitivity_stable_significant
      ),
    by = "Taxon"
  ) %>%
    mutate(
      Rank = rank_label,
      Same_effect_direction = case_when(
        is.na(DESeq2_FDR) | is.na(ANCOMBC2_test_available) |
          !ANCOMBC2_test_available | !is.finite(DESeq2_log2FoldChange) |
          !is.finite(ANCOMBC2_log2_difference) ~ NA,
        sign(DESeq2_log2FoldChange) == sign(ANCOMBC2_log2_difference) ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    arrange(
      desc(ANCOMBC2_BH_sensitivity_stable),
      desc(DESeq2_significant),
      ANCOMBC2_BH_q,
      DESeq2_FDR
    )

  safe_write_csv(
    out,
    file.path(
      subdirs$combined,
      paste0("DESeq2_ANCOMBC2_EXACT_MATCHED_", tolower(rank_label), "_cross_method.csv")
    )
  )

  out
}

cross_genus <- make_cross_method_table(deseq2_genus, ancombc2_genus, "Genus")
cross_species <- make_cross_method_table(deseq2_species, ancombc2_species, "Species")

analysis_status <- bind_rows(
  data.frame(
    Analysis = c(
      "Aitchison beta diversity",
      "Bray-Curtis beta diversity"
    ),
    Rank = "ASV",
    Status = "COMPLETED",
    N_features = beta_summary$N_features,
    Main_result = paste0(
      "R2=", sprintf("%.4f", beta_summary$PERMANOVA_R2),
      "; PERMANOVA p=", formatC(beta_summary$PERMANOVA_P, format = "g", digits = 4),
      "; betadisper p=", formatC(beta_summary$Betadisper_P, format = "g", digits = 4)
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Analysis = "Targeted CLR and prevalence",
    Rank = "Genus",
    Status = if (all(target_clr_stats$Status == "COMPLETED") && all(target_detection_stats$Status == "COMPLETED")) "COMPLETED" else "REVIEW_REQUIRED",
    N_features = length(target_genera),
    Main_result = "Patient-adjusted CLR/detection tests and descriptive biopsy detection exported",
    stringsAsFactors = FALSE
  ),
  data.frame(
    Analysis = "DESeq2",
    Rank = "Genus",
    Status = if (is.null(deseq2_genus)) "FAILED" else "COMPLETED",
    N_features = if (is.null(deseq2_genus)) NA_integer_ else deseq2_genus$summary$N_features_tested,
    Main_result = if (is.null(deseq2_genus)) NA_character_ else paste0(
      deseq2_genus$summary$N_FDR_significant, " FDR-significant genera"
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Analysis = "DESeq2",
    Rank = "Resolved species",
    Status = if (!species_bundle$available) {
      "NOT_AVAILABLE"
    } else if (is.null(deseq2_species)) {
      "FAILED"
    } else {
      "COMPLETED"
    },
    N_features = if (is.null(deseq2_species)) species_bundle$n_features else deseq2_species$summary$N_features_tested,
    Main_result = if (is.null(deseq2_species)) species_bundle$reason else paste0(
      deseq2_species$summary$N_FDR_significant, " FDR-significant resolved species"
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Analysis = "ANCOM-BC2",
    Rank = "Genus",
    Status = if (is.null(ancombc2_genus)) "FAILED" else ancombc2_genus$summary$Run_status,
    N_features = if (is.null(ancombc2_genus)) NA_integer_ else ancombc2_genus$summary$N_features_reported,
    Main_result = if (is.null(ancombc2_genus)) NA_character_ else paste0(
      ancombc2_genus$summary$N_BH_sensitivity_stable,
      " BH-significant and sensitivity-stable genera"
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Analysis = "ANCOM-BC2",
    Rank = "Resolved species",
    Status = if (!species_bundle$available) {
      "NOT_AVAILABLE"
    } else if (is.null(ancombc2_species)) {
      "FAILED"
    } else {
      ancombc2_species$summary$Run_status
    },
    N_features = if (is.null(ancombc2_species)) species_bundle$n_features else ancombc2_species$summary$N_features_reported,
    Main_result = if (is.null(ancombc2_species)) species_bundle$reason else paste0(
      ancombc2_species$summary$N_BH_sensitivity_stable,
      " BH-significant and sensitivity-stable resolved species"
    ),
    stringsAsFactors = FALSE
  )
)

safe_write_csv(
  analysis_status,
  file.path(out_dir, "MASTER_EXACT_MATCHED_ANALYSIS_STATUS.csv")
)

readme_lines <- c(
  "CHAPTER 4 — EXACTLY MATCHED SSL_HP VERSUS H_ADJ_SSL_HP",
  "======================================================",
  "",
  paste0("Run time: ", Sys.time()),
  paste0("Input object: ", physeq_path),
  paste0("Input MD5: ", unname(tools::md5sum(physeq_path))),
  paste0("Output directory: ", out_dir),
  paste0("Keep-7699 branch selected: ", use_keep7699_object),
  paste0("Exact Patient__Site pairs retained: ", nrow(valid_pairs)),
  "",
  "BETA DIVERSITY",
  "  adonis2 uses by = terms and the Group row is checked before statistics are exported.",
  "  Aitchison and Bray-Curtis use ASVs; genus Bray-Curtis is also exported separately.",
  "  PERMANOVA permutations are blocked within Patient; median-dispersion permutations are unrestricted.",
  "  Each plot contains R2, pseudo-F, PERMANOVA P, betadisper F and betadisper P.",
  "  The exact statistics are also exported as CSV tables.",
  "",
  "TARGETED CLR AND PREVALENCE",
  "  Targets: Escherichia-Shigella, Serratia and Fusobacterium.",
  "  Panel A shows patient-group mean CLR with the same mixed model as the main comparison.",
  "  Panel B shows detected/total biopsies, with 12 biopsies per group; patient detection is tested separately.",
  "  Prevalence test results remain in CSV outputs but are not printed on the prevalence figure.",
  "",
  "DESEQ2",
  "  Genus counts are summed within Patient and Group and use ~ Patient + Group; species are analysed separately.",
  "  H_ADJ_SSL_HP is the reference group.",
  "  Size factors: poscounts; Wald tests; BH FDR; no LFC shrinkage.",
  "  Volcano plots are exported at wide page size with labels kept inside the panel.",
  "",
  "ANCOM-BC2",
  paste0("  ANCOMBC version: ", packageVersion("ANCOMBC")),
  "  Genus counts are summed within Patient and Group and use Patient + Group; species are analysed separately.",
  "  No random effect and no multi-group pairwise module are used.",
  "  BH FDR is primary; Holm is exported as a stricter sensitivity check.",
  "  Pseudocount sensitivity and structural-zero detection are enabled.",
  "  prv_cut = 0: no second prevalence filter is imposed.",
  "",
  "SPECIES-LEVEL INTERPRETATION",
  "  Only ASVs with a resolved species annotation are included.",
  "  Unclassified or uncultured ASVs are not treated as species.",
  "  A 16S species label remains a taxonomic assignment and does not guarantee strain resolution.",
  "",
  "FILES TO OPEN FIRST",
  "  MASTER_EXACT_MATCHED_ANALYSIS_STATUS.csv",
  "  00_PAIRING_AND_QC/06_valid_exact_patient_site_pairs.csv",
  "  02_AITCHISON/PCoA_Aitchison_EXACT_MATCHED_statistics.csv",
  "  03_BRAY_CURTIS/PCoA_Bray_Curtis_EXACT_MATCHED_statistics.csv",
  "  09_COMBINED_FIGURES/Figure_EXACT_MATCHED_Aitchison_and_Bray_Curtis.pdf",
  "  09_COMBINED_FIGURES/Figure_EXACT_MATCHED_target_CLR_and_prevalence.pdf",
  "  05_DESEQ2_GENUS/DESeq2_EXACT_MATCHED_genus_ALL_results.csv",
  "  07_ANCOMBC2_GENUS/ANCOMBC2_tidy_all_model_based_results.csv"
)

writeLines(readme_lines, file.path(out_dir, "READ_ME_FIRST.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

cat("\n===== BETA-DIVERSITY RESULTS =====\n")
print(beta_summary)
cat("\n===== TARGETED CLR RESULTS =====\n")
print(target_clr_stats)
cat("\n===== TARGETED PREVALENCE RESULTS =====\n")
print(target_detection_stats)
cat("\n===== TAXONOMIC-RANK INPUT SUMMARY =====\n")
print(rank_input_summary)
cat("\n===== MASTER STATUS =====\n")
print(analysis_status)
cat("\nAll exact matched analyses completed.\n")
cat("Output directory:\n", out_dir, "\n")
cat("Open first:\n", file.path(out_dir, "READ_ME_FIRST.txt"), "\n")
