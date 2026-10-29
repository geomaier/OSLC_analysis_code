# Gastric-signature GSEA across passages and tissue-of-origin contrasts.
# Inputs, output locations and the submitted figure override are described in ../README.md.

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

# Required packages ----
cran_packages <- c(
  "data.table", "dplyr", "purrr", "tidyr", "stringr",
  "ggplot2", "readr", "scales"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0L) {
  stop("Install required packages before running: ", paste(missing_cran, collapse = ", "))
}
if (!requireNamespace("fgsea", quietly = TRUE)) {
  stop("Install the required Bioconductor package fgsea before running.")
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(fgsea)
})

# Input files and output directory ----
FILE_PATTERN <- paste0(
  "^.+_Passage_.+_P[0-9]+_v_P[0-9]+_DE_results\\.txt$"
)

SEARCH_ROOT <- Sys.getenv("OSLC_ORGANOID_ROOT", getwd())

DATA_DIR_OVERRIDE <- required_path("OSLC_ORGANOID_PASSAGE_DIR")

DATA_DIR <- DATA_DIR_OVERRIDE

if (is.na(DATA_DIR) || !nzchar(DATA_DIR)) {
  DATA_DIR <- NA_character_

  if (dir.exists(SEARCH_ROOT)) {
    candidate_files <- list.files(
      SEARCH_ROOT,
      pattern = FILE_PATTERN,
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    )

    if (length(candidate_files) > 0L) {
      candidate_counts <- sort(table(dirname(candidate_files)), decreasing = TRUE)
      DATA_DIR <- unname(names(candidate_counts)[1])

      if (length(candidate_counts) > 1L) {
        message(
          "Passage DE tables were found in more than one directory. ",
          "Using the directory containing the largest number of matching files:\n",
          DATA_DIR
        )
      }
    }
  }
}

if ((is.na(DATA_DIR) || !dir.exists(DATA_DIR)) &&
    .Platform$OS.type == "windows" && interactive()) {
  selected_dir <- utils::choose.dir(
    default = if (dir.exists(SEARCH_ROOT)) SEARCH_ROOT else getwd(),
    caption = paste(
      "Select the folder containing the CGG028_Passage_*_DE_results.txt files"
    )
  )

  if (!is.na(selected_dir) && nzchar(selected_dir)) {
    DATA_DIR <- selected_dir
  }
}

if (is.na(DATA_DIR) || !dir.exists(DATA_DIR)) {
  stop(
    "The passage DE-table directory could not be found.\n",
    "Select or enter the folder that directly contains files such as:\n",
    "CGG028_Passage_H_ADJ_P2_v_P1_DE_results.txt\n\n",
    "Search root checked:\n", SEARCH_ROOT
  )
}

DATA_DIR <- normalizePath(DATA_DIR, winslash = "/", mustWork = TRUE)

matching_de_files <- list.files(
  DATA_DIR,
  pattern = FILE_PATTERN,
  full.names = TRUE,
  recursive = FALSE,
  ignore.case = TRUE
)

if (length(matching_de_files) == 0L) {
  stop(
    "The selected directory exists, but it does not directly contain any ",
    "matching passage DE files:\n", DATA_DIR, "\n\n",
    "Select the folder containing files such as ",
    "CGG028_Passage_H_ADJ_P2_v_P1_DE_results.txt."
  )
}

message("Using DATA_DIR:\n", DATA_DIR)
message("Found ", length(matching_de_files), " matching passage DE tables.")

OUTPUT_DIR <- review_output_dir(file.path(
  DATA_DIR,
  "Final_24_gene_gastric_metaplasia_GSEA_THESIS_FINAL"
))

PRESERVE_EXISTING_OUTPUTS <- TRUE

if (PRESERVE_EXISTING_OUTPUTS && dir.exists(OUTPUT_DIR)) {
  existing_output_files <- list.files(
    OUTPUT_DIR,
    full.names = TRUE,
    recursive = FALSE,
    all.files = FALSE
  )
  existing_output_files <- existing_output_files[
    file.exists(existing_output_files) & !dir.exists(existing_output_files)
  ]

  if (length(existing_output_files) > 0L) {
    archive_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    archive_dir <- file.path(
      OUTPUT_DIR,
      paste0("PREVIOUS_OUTPUTS_", archive_timestamp)
    )

    archive_suffix <- 1L
    while (dir.exists(archive_dir)) {
      archive_suffix <- archive_suffix + 1L
      archive_dir <- file.path(
        OUTPUT_DIR,
        paste0(
          "PREVIOUS_OUTPUTS_",
          archive_timestamp,
          "_",
          archive_suffix
        )
      )
    }

    dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)

    copied_ok <- file.copy(
      from = existing_output_files,
      to = archive_dir,
      overwrite = FALSE,
      copy.mode = TRUE,
      copy.date = TRUE
    )

    if (!all(copied_ok)) {
      warning(
        "Some existing output files could not be archived before the new run. ",
        "Check: ", archive_dir,
        call. = FALSE
      )
    } else {
      message(
        "Archived ", length(existing_output_files),
        " existing output file(s) to:\n", archive_dir
      )
    }
  }
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# GSEA settings and fixed 24-gene signature ----
MIN_LOGCPM <- -Inf

USE_COMMON_GENE_UNIVERSE <- TRUE
MIN_RANKED_GENES <- 1000L

MIN_SIGNATURE_GENES_TO_RUN <- 10L
PREFERRED_SIGNATURE_GENES <- 15L

STAR_BY <- "q_global_BH"

GSEA_SEED <- 26082026L

EXPECTED_REFERENCE_PASSAGE <- "P1"

CROSS_SECTIONAL_FILE_NAMES <- c(
  P1 = "CGG028_First_Passage_SSL_HP_DE_results.txt",
  P4 = "CGG028_Fourth_Passage_SSL_HP_DE_results.txt",
  P8 = "CGG028_Eighth_Passage_SSL_HP_DE_results.txt"
)
CROSS_SECTIONAL_PASSAGE_ORDER <- names(CROSS_SECTIONAL_FILE_NAMES)
CROSS_SECTIONAL_CONTRAST <- "SSL_HP - H_ADJ_SSL_HP"
CROSS_SECTIONAL_POSITIVE_NES <- "enriched in SSL_HP relative to H_ADJ_SSL_HP"
CROSS_SECTIONAL_GENE_LOGFC_LIMIT <- 2

CROSS_SECTIONAL_STAR_BY <- "q_cross_sectional_BH"

RUN_PVALUE_RANK_SENSITIVITY <- TRUE
RUN_PER_TABLE_UNIVERSE_SENSITIVITY <- TRUE

NEGATIVE_LR_MIN_NULL_PVALUE <- 0.999
LR_NUMERICAL_ABS_TOL <- 1e-8
LR_NUMERICAL_REL_TOL <- 100 * .Machine$double.eps

METAPLASIA_SIGNATURE_TITLE <-
  "Serrated Lesion Gastric Metaplasia Gene Signature"
METAPLASIA_SIGNATURE_FILE_STEM <-
  "Serrated_Lesion_Gastric_Metaplasia_Gene_Signature"

gastric_metaplasia_24 <- c(
  "AQP5", "ANXA10", "MUC5AC", "VSIG1", "TCN1", "MUC6",
  "PSAPL1", "CLDN18", "SFTA2", "MUCL3", "GATA4", "KLK11",
  "DUOXA1", "SULT1C2", "IRX3", "AADAC", "FOXQ1", "IGFALS",
  "NR0B2", "C13ORF46", "FOLR1", "LYPD6B", "TFR2", "SPTSSB"
)

gastric_metaplasia_24 <- unique(toupper(trimws(gastric_metaplasia_24)))

if (length(gastric_metaplasia_24) != 24L) {
  stop(
    METAPLASIA_SIGNATURE_TITLE,
    " must contain exactly 24 unique genes."
  )
}

readr::write_lines(
  gastric_metaplasia_24,
  file.path(OUTPUT_DIR, "final_24_gene_gastric_metaplasia_signature_used.txt")
)
readr::write_lines(
  gastric_metaplasia_24,
  file.path(
    OUTPUT_DIR,
    paste0(METAPLASIA_SIGNATURE_FILE_STEM, "_genes_used.txt")
  )
)

if (!STAR_BY %in% c("q_global_BH", "pval")) {
  stop("STAR_BY must be either 'q_global_BH' or 'pval'.")
}

# DE-table import and ranking helpers ----
normalise_column_name <- function(x) {
  gsub("[^a-z0-9]", "", tolower(x))
}

pick_column <- function(df, candidates, required = TRUE, description = NULL) {
  available_normalised <- normalise_column_name(names(df))
  candidate_normalised <- normalise_column_name(candidates)
  hits <- match(candidate_normalised, available_normalised)
  hits <- hits[!is.na(hits)]

  if (length(hits) > 0L) {
    return(names(df)[hits[1]])
  }

  if (required) {
    stop(
      "Could not find ", description %||% "the required column", ". Tried: ",
      paste(candidates, collapse = ", "), ". Available columns: ",
      paste(names(df), collapse = ", ")
    )
  }

  NA_character_
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x[1])) y else x
}

to_numeric <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", as.character(x), fixed = TRUE)))
}

parse_de_filename <- function(file) {
  file_name <- basename(file)
  stem <- stringr::str_remove(
    file_name,
    stringr::regex("_DE_results\\.txt$", ignore_case = TRUE)
  )

  parsed <- stringr::str_match(
    stem,
    stringr::regex(
      "^(.*?)_Passage_(.+)_(P[0-9]+)_v_(P[0-9]+)$",
      ignore_case = TRUE
    )
  )

  if (is.na(parsed[1, 1])) {
    stop(
      "Filename does not match the required pattern: ", file_name,
      ". Expected, for example: ",
      "CGG028_Passage_H_ADJ_P2_v_P1_DE_results.txt"
    )
  }

  tibble::tibble(
    file = file,
    file_name = file_name,
    line_id = parsed[1, 2],
    condition = toupper(parsed[1, 3]),
    numerator_passage = toupper(parsed[1, 4]),
    reference_passage = toupper(parsed[1, 5]),
    comparison = paste0(
      toupper(parsed[1, 4]), " vs ", toupper(parsed[1, 5])
    ),
    contrast_definition_assumed = paste0(
      toupper(parsed[1, 4]), " - ", toupper(parsed[1, 5])
    ),
    positive_NES_interpretation = paste0(
      "enriched in ", toupper(parsed[1, 4]),
      " relative to ", toupper(parsed[1, 5])
    )
  )
}

make_cross_sectional_meta <- function(file, passage_label) {
  passage_label <- toupper(trimws(as.character(passage_label)))

  if (!passage_label %in% CROSS_SECTIONAL_PASSAGE_ORDER) {
    stop(
      "Unexpected cross-sectional passage label: ", passage_label,
      ". Expected one of: ",
      paste(CROSS_SECTIONAL_PASSAGE_ORDER, collapse = ", ")
    )
  }

  tibble::tibble(
    file = file,
    file_name = basename(file),
    line_id = "CGG028",
    condition = "SSL_HP_vs_H_ADJ_SSL_HP",
    numerator_passage = passage_label,
    reference_passage = "H_ADJ_SSL_HP",
    comparison = passage_label,
    contrast_definition_assumed = paste0(
      CROSS_SECTIONAL_CONTRAST, " at ", passage_label
    ),
    positive_NES_interpretation = CROSS_SECTIONAL_POSITIVE_NES
  )
}

read_and_rank_de_table <- function(file, meta_override = NULL) {
  meta <- if (is.null(meta_override)) {
    parse_de_filename(file)
  } else {
    meta_override
  }

  required_meta_columns <- c(
    "file", "file_name", "line_id", "condition",
    "numerator_passage", "reference_passage", "comparison",
    "contrast_definition_assumed", "positive_NES_interpretation"
  )

  missing_meta_columns <- setdiff(required_meta_columns, names(meta))
  if (length(missing_meta_columns) > 0L || nrow(meta) != 1L) {
    stop(
      "The DE-table metadata must contain exactly one row and these columns: ",
      paste(required_meta_columns, collapse = ", "),
      if (length(missing_meta_columns) > 0L) {
        paste0(". Missing: ", paste(missing_meta_columns, collapse = ", "))
      } else {
        "."
      }
    )
  }

  de <- data.table::fread(
    file,
    sep = "auto",
    header = TRUE,
    data.table = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN", "NULL")
  )

  if (nrow(de) < 2L) {
    stop("The DE table is empty or too short: ", basename(file))
  }

  gene_col <- pick_column(
    de,
    c(
      "SYMBOL", "gene_symbol", "GeneSymbol", "hgnc_symbol",
      "external_gene_name", "Gene", "genes", "gene"
    ),
    description = "a gene-symbol column"
  )

  logfc_col <- pick_column(
    de,
    c("logFC", "log2FC", "log2FoldChange", "log2_fold_change"),
    description = "a log2 fold-change column"
  )

  lr_col <- pick_column(
    de,
    c("LR", "LRT", "likelihood_ratio", "likelihood ratio"),
    required = FALSE
  )

  pvalue_col <- pick_column(
    de,
    c("PValue", "P.Value", "pvalue", "p_val", "pval", "p"),
    required = FALSE
  )

  fdr_col <- pick_column(
    de,
    c("FDR", "padj", "adj.P.Val", "adjusted_pvalue"),
    required = FALSE
  )

  logcpm_col <- pick_column(
    de,
    c("logCPM", "AveExpr", "baseMean"),
    required = FALSE
  )

  gene <- toupper(trimws(as.character(de[[gene_col]])))
  gene[gene %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  logfc <- to_numeric(de[[logfc_col]])

  pvalue <- if (!is.na(pvalue_col)) {
    to_numeric(de[[pvalue_col]])
  } else {
    rep(NA_real_, nrow(de))
  }

  if (
    !is.na(pvalue_col) &&
      any(
        !is.na(pvalue) &
          (!is.finite(pvalue) | pvalue < 0 | pvalue > 1)
      )
  ) {
    stop("Non-finite P values or P values outside [0, 1] were found in ",
         basename(file), ".")
  }

  fdr <- if (!is.na(fdr_col)) {
    to_numeric(de[[fdr_col]])
  } else {
    rep(NA_real_, nrow(de))
  }

  logcpm <- if (!is.na(logcpm_col)) {
    to_numeric(de[[logcpm_col]])
  } else {
    rep(NA_real_, nrow(de))
  }

  if (!is.na(lr_col)) {
    raw_lr <- de[[lr_col]]
    lr <- to_numeric(raw_lr)

    finite_lr <- lr[is.finite(lr)]
    n_finite_lr <- length(finite_lr)
    n_negative_lr <- sum(finite_lr < 0)
    n_missing_lr <- sum(is.na(lr))
    n_infinite_lr <- sum(is.infinite(lr))
    min_lr <- if (n_finite_lr > 0L) min(finite_lr) else NA_real_
    max_lr <- if (n_finite_lr > 0L) max(finite_lr) else NA_real_
    negative_lr_fraction <- if (n_finite_lr > 0L) {
      n_negative_lr / n_finite_lr
    } else {
      NA_real_
    }

    if (n_finite_lr == 0L) {
      stop("No finite LR values were found in ", basename(file), ".")
    }

    lr_negative_tolerance <- max(
      LR_NUMERICAL_ABS_TOL,
      LR_NUMERICAL_REL_TOL * max(1, max(abs(finite_lr)))
    )

    negative_idx <- which(is.finite(lr) & lr < 0)
    tiny_negative_idx <- which(
      is.finite(lr) & lr < 0 & lr >= -lr_negative_tolerance
    )
    material_negative_idx <- which(
      is.finite(lr) & lr < -lr_negative_tolerance
    )
    infinite_idx <- which(is.infinite(lr))

    n_tiny_negative_lr <- length(tiny_negative_idx)
    n_material_negative_lr <- length(material_negative_idx)
    negative_pvalues <- pvalue[negative_idx]
    n_negative_lr_missing_pvalue <- sum(!is.finite(negative_pvalues))
    n_negative_lr_nonnull_pvalue <- sum(
      is.finite(negative_pvalues) &
        negative_pvalues < NEGATIVE_LR_MIN_NULL_PVALUE
    )
    min_pvalue_among_negative_lr <- if (
      any(is.finite(negative_pvalues))
    ) {
      min(negative_pvalues[is.finite(negative_pvalues)])
    } else {
      NA_real_
    }

    diagnostic_idx <- sort(unique(c(negative_idx, infinite_idx)))
    lr_diagnostic_file <- NA_character_

    if (length(diagnostic_idx) > 0L) {
      input_stem <- stringr::str_remove(
        basename(file),
        stringr::regex("\\.txt$", ignore_case = TRUE)
      )
      lr_diagnostic_file <- file.path(
        OUTPUT_DIR,
        paste0(input_stem, "_LR_diagnostic.csv")
      )

      lr_diagnostic <- tibble::tibble(
        input_row = diagnostic_idx,
        Gene = gene[diagnostic_idx],
        logFC = logfc[diagnostic_idx],
        LR_raw = as.character(raw_lr[diagnostic_idx]),
        LR_numeric = lr[diagnostic_idx],
        PValue = pvalue[diagnostic_idx],
        FDR = fdr[diagnostic_idx],
        LR_numerical_tolerance = lr_negative_tolerance,
        classification = dplyr::case_when(
          is.infinite(lr[diagnostic_idx]) ~ "infinite_invalid",
          lr[diagnostic_idx] < -lr_negative_tolerance ~
            "negative_beyond_numerical_tolerance",
          TRUE ~ "tiny_negative_numerical"
        )
      )

      readr::write_csv(lr_diagnostic, lr_diagnostic_file, na = "")
    }

    if (n_infinite_lr > 0L) {
      stop(
        basename(file), " contains ", n_infinite_lr,
        " infinite LR value(s). Offending rows were written to:\n",
        lr_diagnostic_file
      )
    }

    if (n_negative_lr > 0L) {
      if (is.na(pvalue_col)) {
        stop(
          "Negative LR values were found in ", basename(file),
          ", but no PValue column is available to verify that they represent ",
          "null/numerical rows. Details were written to:\n",
          lr_diagnostic_file
        )
      }

      if (
        n_negative_lr_missing_pvalue > 0L ||
          n_negative_lr_nonnull_pvalue > 0L
      ) {
        stop(
          basename(file), " contains negative LR rows whose P values are ",
          "missing or below ", NEGATIVE_LR_MIN_NULL_PVALUE,
          ". This is inconsistent with a standard edgeR LR table. ",
          "Do not use abs(LR). Inspect:\n", lr_diagnostic_file
        )
      }

      warning(
        basename(file), ": setting ", n_negative_lr,
        " negative LR value", ifelse(n_negative_lr == 1L, "", "s"),
        " to 0 before ranking after confirming null P values. ",
        "Tiny negatives: ", n_tiny_negative_lr,
        "; negatives beyond the numerical tolerance: ",
        n_material_negative_lr,
        "; minimum LR = ",
        format(min_lr, scientific = TRUE, digits = 6),
        "; minimum P value among negative LR rows = ",
        format(min_pvalue_among_negative_lr, scientific = TRUE, digits = 6),
        ". These rows contribute rank 0; abs(LR) is not used. Details: ",
        lr_diagnostic_file,
        call. = FALSE
      )
    }

    lr_for_rank <- lr
    lr_for_rank[negative_idx] <- 0
    rank_metric <- sign(logfc) * sqrt(lr_for_rank)
    rank_method <- "sign(logFC) * sqrt(LR)"
  } else {
    n_finite_lr <- NA_integer_
    n_negative_lr <- NA_integer_
    n_tiny_negative_lr <- NA_integer_
    n_material_negative_lr <- NA_integer_
    n_missing_lr <- NA_integer_
    n_infinite_lr <- NA_integer_
    min_lr <- NA_real_
    max_lr <- NA_real_
    negative_lr_fraction <- NA_real_
    lr_negative_tolerance <- NA_real_
    n_negative_lr_missing_pvalue <- NA_integer_
    n_negative_lr_nonnull_pvalue <- NA_integer_
    min_pvalue_among_negative_lr <- NA_real_
    lr_diagnostic_file <- NA_character_

    if (is.na(pvalue_col)) {
      stop(
        basename(file),
        " has neither an LR column nor a nominal P-value column."
      )
    }

    if (any(pvalue < 0 | pvalue > 1, na.rm = TRUE)) {
      stop("P values outside [0, 1] were found in ", basename(file), ".")
    }

    p_safe <- pmax(pvalue, 1e-300)
    rank_metric <- sign(logfc) * -log10(p_safe)
    rank_method <- "sign(logFC) * -log10(PValue) [fallback]"
  }

  ranked <- tibble::tibble(
    Gene = gene,
    logFC = logfc,
    logCPM = logcpm,
    PValue = pvalue,
    FDR = fdr,
    rank_metric = rank_metric
  ) %>%
    dplyr::filter(
      !is.na(Gene),
      nzchar(Gene),
      is.finite(logFC),
      is.finite(rank_metric)
    )

  if (is.finite(MIN_LOGCPM)) {
    if (is.na(logcpm_col)) {
      stop(
        "MIN_LOGCPM is finite, but no logCPM-like column was found in ",
        basename(file), "."
      )
    }
    ranked <- ranked %>% dplyr::filter(is.finite(logCPM), logCPM >= MIN_LOGCPM)
  }

  ranked_uncollapsed <- ranked

  ranked <- ranked %>%
    dplyr::group_by(Gene) %>%
    dplyr::arrange(
      dplyr::desc(abs(rank_metric)),
      PValue,
      FDR,
      dplyr::desc(abs(logFC)),
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(n = 1L) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(rank_metric), Gene)

  if (nrow(ranked) < MIN_RANKED_GENES) {
    stop(
      basename(file), " contains only ", nrow(ranked),
      " unique ranked genes. GSEA requires the full tested-gene table, ",
      "not a table already filtered to significant DE genes."
    )
  }

  audit <- meta %>%
    dplyr::mutate(
      n_input_rows = nrow(de),
      n_unique_ranked_genes = nrow(ranked),
      gene_column = gene_col,
      logfc_column = logfc_col,
      lr_column = lr_col,
      pvalue_column = pvalue_col,
      fdr_column = fdr_col,
      logcpm_column = logcpm_col,
      rank_method = rank_method,
      n_finite_LR = n_finite_lr,
      n_negative_LR_raw = n_negative_lr,
      n_tiny_negative_LR = n_tiny_negative_lr,
      n_negative_LR_beyond_numerical_tolerance = n_material_negative_lr,
      n_negative_LR_floored_to_zero = n_negative_lr,
      n_missing_LR = n_missing_lr,
      n_infinite_LR = n_infinite_lr,
      negative_LR_fraction = negative_lr_fraction,
      LR_numerical_tolerance = lr_negative_tolerance,
      minimum_LR_before_flooring = min_lr,
      maximum_LR = max_lr,
      n_negative_LR_missing_PValue = n_negative_lr_missing_pvalue,
      n_negative_LR_with_PValue_below_null_threshold =
        n_negative_lr_nonnull_pvalue,
      minimum_PValue_among_negative_LR = min_pvalue_among_negative_lr,
      LR_diagnostic_file = lr_diagnostic_file
    )

  list(meta = meta, ranked = ranked, ranked_uncollapsed = ranked_uncollapsed, audit = audit)
}

make_rank_vector <- function(ranked, universe) {
  ranked_use <- ranked %>%
    dplyr::filter(Gene %in% universe) %>%
    dplyr::arrange(dplyr::desc(rank_metric), Gene)

  ranks <- ranked_use$rank_metric
  names(ranks) <- ranked_use$Gene

  if (length(ranks) < 2L || length(unique(ranks)) < 2L) {
    stop("The ranked gene vector has fewer than two distinct values.")
  }

  tied <- duplicated(ranks) | duplicated(ranks, fromLast = TRUE)
  if (any(tied)) {
    scale_value <- max(1, max(abs(ranks), na.rm = TRUE))
    within_tie <- ave(seq_along(ranks)[tied], ranks[tied], FUN = seq_along)
    ranks[tied] <- ranks[tied] -
      (within_tie - 1L) * (.Machine$double.eps * scale_value)
  }

  sort(ranks, decreasing = TRUE)
}

empty_gsea_row <- function(meta, n_ranked, overlap, missing, status) {
  meta %>%
    dplyr::mutate(
      pathway = METAPLASIA_SIGNATURE_TITLE,
      n_ranked_genes_used = n_ranked,
      intended_signature_size = length(gastric_metaplasia_24),
      represented_signature_size = length(overlap),
      represented_signature_genes = paste(overlap, collapse = ";"),
      missing_signature_genes = paste(missing, collapse = ";"),
      ES = NA_real_,
      NES = NA_real_,
      pval = NA_real_,
      fgsea_padj_within_file = NA_real_,
      log2err = NA_real_,
      leading_edge = NA_character_,
      status = status
    )
}

# Single-signature GSEA with recorded gene coverage ----
run_one_gsea <- function(rank_object, universe, seed) {
  meta <- rank_object$meta
  ranks <- make_rank_vector(rank_object$ranked, universe)
  overlap <- intersect(gastric_metaplasia_24, names(ranks))
  missing <- setdiff(gastric_metaplasia_24, names(ranks))

  if (length(overlap) < MIN_SIGNATURE_GENES_TO_RUN) {
    return(
      empty_gsea_row(
        meta = meta,
        n_ranked = length(ranks),
        overlap = overlap,
        missing = missing,
        status = paste0(
          "NOT_RUN: only ", length(overlap),
          " signature genes represented (<", MIN_SIGNATURE_GENES_TO_RUN, ")"
        )
      )
    )
  }

  set.seed(seed)

  fg <- fgsea::fgseaMultilevel(
    pathways = stats::setNames(
      list(gastric_metaplasia_24),
      METAPLASIA_SIGNATURE_TITLE
    ),
    stats = ranks,
    minSize = MIN_SIGNATURE_GENES_TO_RUN,
    maxSize = length(gastric_metaplasia_24),
    eps = 0,
    scoreType = "std",
    gseaParam = 1,
    sampleSize = 1001,
    nPermSimple = 10000,
    nproc = 1
  )

  if (nrow(fg) != 1L) {
    return(
      empty_gsea_row(
        meta = meta,
        n_ranked = length(ranks),
        overlap = overlap,
        missing = missing,
        status = "NOT_RUN: fgsea did not return exactly one pathway row"
      )
    )
  }

  if (!is.finite(fg$NES[1]) || !is.finite(fg$pval[1])) {
    return(
      empty_gsea_row(
        meta = meta,
        n_ranked = length(ranks),
        overlap = overlap,
        missing = missing,
        status = "NOT_RUN: fgsea returned a non-finite NES or P value"
      )
    )
  }

  overlap_status <- if (length(overlap) >= PREFERRED_SIGNATURE_GENES) {
    "OK"
  } else {
    paste0("LOW_OVERLAP: ", length(overlap), "/24 genes")
  }

  meta %>%
    dplyr::mutate(
      pathway = as.character(fg$pathway[1]),
      n_ranked_genes_used = length(ranks),
      intended_signature_size = length(gastric_metaplasia_24),
      represented_signature_size = length(overlap),
      represented_signature_genes = paste(overlap, collapse = ";"),
      missing_signature_genes = paste(missing, collapse = ";"),
      ES = as.numeric(fg$ES[1]),
      NES = as.numeric(fg$NES[1]),
      pval = as.numeric(fg$pval[1]),
      
      fgsea_padj_within_file = as.numeric(fg$padj[1]),
      log2err = as.numeric(fg$log2err[1]),
      leading_edge = paste(as.character(fg$leadingEdge[[1]]), collapse = ";"),
      status = overlap_status
    )
}

significance_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.0001 ~ "****",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

# Longitudinal passage contrasts ----
de_files <- list.files(
  DATA_DIR,
  pattern = FILE_PATTERN,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(de_files) == 0L) {
  stop(
    "No passage DE files were found in:\n", DATA_DIR,
    "\nExpected filenames such as ",
    "CGG028_Passage_H_ADJ_P2_v_P1_DE_results.txt"
  )
}

de_files <- sort(de_files)
message("Found ", length(de_files), " passage DE tables.")

contrast_manifest <- purrr::map_dfr(de_files, parse_de_filename)
readr::write_csv(
  contrast_manifest,
  file.path(OUTPUT_DIR, "passage_GSEA_contrast_manifest.csv"),
  na = ""
)

message(
  "Direction convention retained from the original script: Pn_v_P1 is ",
  "interpreted as Pn - P1. The contrast manifest was written to: ",
  file.path(OUTPUT_DIR, "passage_GSEA_contrast_manifest.csv")
)

unexpected_reference <- contrast_manifest %>%
  dplyr::filter(reference_passage != EXPECTED_REFERENCE_PASSAGE)

if (nrow(unexpected_reference) > 0L) {
  print(unexpected_reference)
  stop(
    "At least one file does not use the expected reference passage ",
    EXPECTED_REFERENCE_PASSAGE,
    ". Inspect passage_GSEA_contrast_manifest.csv before continuing."
  )
}

message(
  "Positive NES is interpreted as enrichment in the numerator passage Pn ",
  "relative to the P1 reference, following the filename convention."
)

rank_objects <- purrr::map(de_files, read_and_rank_de_table)

column_audit <- purrr::map_dfr(rank_objects, "audit")
readr::write_csv(
  column_audit,
  file.path(OUTPUT_DIR, "passage_GSEA_input_column_audit.csv"),
  na = ""
)

ranking_methods <- unique(column_audit$rank_method)
if (length(ranking_methods) != 1L) {
  stop(
    "Different ranking metrics would be used across files: ",
    paste(ranking_methods, collapse = " | "),
    ". All heatmap comparisons must use the same method."
  )
}

all_gene_universes <- purrr::map(rank_objects, ~ .x$ranked$Gene)

if (USE_COMMON_GENE_UNIVERSE) {
  gsea_universe <- Reduce(intersect, all_gene_universes)
} else {
  
  gsea_universe <- NULL
}

if (USE_COMMON_GENE_UNIVERSE && length(gsea_universe) < MIN_RANKED_GENES) {
  stop(
    "The common gene universe contains only ", length(gsea_universe),
    " genes. Check gene-symbol mapping and confirm that every input file ",
    "contains all tested genes."
  )
}

gsea_results <- purrr::imap_dfr(rank_objects, function(rank_object, i) {
  universe_i <- if (USE_COMMON_GENE_UNIVERSE) {
    gsea_universe
  } else {
    rank_object$ranked$Gene
  }

  run_one_gsea(
    rank_object = rank_object,
    universe = universe_i,
    seed = GSEA_SEED + as.integer(i)
  )
}) %>%
  dplyr::mutate(
    q_global_BH = stats::p.adjust(pval, method = "BH")
  )

star_p <- gsea_results[[STAR_BY]]
gsea_results$stars <- significance_stars(star_p)
gsea_results$star_basis <- if (STAR_BY == "q_global_BH") {
  "Global Benjamini-Hochberg q-value across all passage comparisons"
} else {
  "Nominal fgsea P value (exploratory)"
}

duplicate_cells <- gsea_results %>%
  dplyr::count(line_id, condition, comparison, name = "n_files") %>%
  dplyr::filter(n_files > 1L)

if (nrow(duplicate_cells) > 0L) {
  print(duplicate_cells)
  stop(
    "At least two DE files map to the same line/condition/comparison cell. ",
    "Remove or rename duplicate files before plotting."
  )
}

readr::write_csv(
  gsea_results,
  file.path(OUTPUT_DIR, "final_24_gene_passage_GSEA_results.csv"),
  na = ""
)

signature_presence <- purrr::map_dfr(rank_objects, function(rank_object) {
  table_gene_universe <- rank_object$ranked$Gene
  common_gene_universe <- if (USE_COMMON_GENE_UNIVERSE) {
    gsea_universe
  } else {
    Reduce(intersect, all_gene_universes)
  }

  rank_object$meta %>%
    dplyr::slice(rep(1L, length(gastric_metaplasia_24))) %>%
    dplyr::mutate(
      Gene = gastric_metaplasia_24,
      present_in_this_DE_table = Gene %in% table_gene_universe,
      present_in_common_universe = Gene %in% common_gene_universe,
      used_in_primary_GSEA = if (USE_COMMON_GENE_UNIVERSE) {
        present_in_common_universe
      } else {
        present_in_this_DE_table
      }
    )
})

readr::write_csv(
  signature_presence,
  file.path(OUTPUT_DIR, "final_24_gene_presence_by_DE_table_AND_common_universe.csv"),
  na = ""
)

signature_presence_summary <- signature_presence %>%
  dplyr::group_by(file_name, line_id, condition, comparison) %>%
  dplyr::summarise(
    intended_signature_size = length(gastric_metaplasia_24),
    n_present_in_this_DE_table = sum(present_in_this_DE_table),
    n_present_in_common_universe = sum(present_in_common_universe),
    n_used_in_primary_GSEA = sum(used_in_primary_GSEA),
    missing_from_this_DE_table = paste(Gene[!present_in_this_DE_table], collapse = ";"),
    missing_from_common_universe = paste(Gene[!present_in_common_universe], collapse = ";"),
    .groups = "drop"
  )

readr::write_csv(
  signature_presence_summary,
  file.path(OUTPUT_DIR, "final_24_gene_presence_summary.csv"),
  na = ""
)

represented_sizes_primary <- unique(
  gsea_results$represented_signature_size[is.finite(gsea_results$represented_signature_size)]
)

if (length(represented_sizes_primary) == 0L) {
  represented_signature_size_primary <- NA_integer_
} else if (length(represented_sizes_primary) == 1L) {
  represented_signature_size_primary <- as.integer(represented_sizes_primary)
} else {
  represented_signature_size_primary <- min(represented_sizes_primary)
  warning(
    "The represented signature size varies across primary comparisons. ",
    "The heatmap subtitle reports the minimum; inspect the result table."
  )
}

missing_signature_primary <- sort(unique(unlist(
  strsplit(
    paste(gsea_results$missing_signature_genes, collapse = ";"),
    split = ";",
    fixed = TRUE
  )
)))
missing_signature_primary <- missing_signature_primary[
  !is.na(missing_signature_primary) & missing_signature_primary != ""
]

comparison_order <- gsea_results %>%
  dplyr::distinct(comparison, numerator_passage, reference_passage) %>%
  dplyr::mutate(
    numerator_number = as.integer(stringr::str_remove(numerator_passage, "^P")),
    reference_number = as.integer(stringr::str_remove(reference_passage, "^P"))
  ) %>%
  dplyr::arrange(reference_number, numerator_number, comparison) %>%
  dplyr::pull(comparison)

comparison_order <- unique(comparison_order)

condition_priority <- c("H_ADJ", "H_ADJ_SSL_HP", "SSL_HP")

row_metadata <- gsea_results %>%
  dplyr::distinct(line_id, condition) %>%
  dplyr::mutate(
    condition_order = match(condition, condition_priority),
    condition_order = ifelse(
      is.na(condition_order),
      length(condition_priority) + dplyr::row_number(),
      condition_order
    ),
    row_label = if (dplyr::n_distinct(gsea_results$line_id) == 1L) {
      condition
    } else {
      paste(line_id, condition, sep = " | ")
    }
  ) %>%
  dplyr::arrange(line_id, condition_order, condition)

row_order <- row_metadata$row_label

plot_results <- gsea_results %>%
  dplyr::left_join(row_metadata, by = c("line_id", "condition"))

heatmap_data <- tidyr::expand_grid(
  row_label = row_order,
  comparison = comparison_order
) %>%
  dplyr::left_join(
    plot_results,
    by = c("row_label", "comparison")
  ) %>%
  dplyr::mutate(
    comparison = factor(comparison, levels = comparison_order),
    row_label = factor(row_label, levels = rev(row_order)),
    cell_label = dplyr::case_when(
      is.na(NES) & !is.na(status) ~ "NA",
      is.na(NES) ~ "",
      stars == "" ~ sprintf("%.2f", NES),
      TRUE ~ paste0(sprintf("%.2f", NES), "\n", stars)
    )
  )

finite_nes <- heatmap_data$NES[is.finite(heatmap_data$NES)]

if (length(finite_nes) == 0L) {
  stop(
    "No finite NES values were produced. Inspect ",
    file.path(OUTPUT_DIR, "final_24_gene_passage_GSEA_results.csv")
  )
}

colour_limit <- max(abs(finite_nes))
if (colour_limit == 0) {
  colour_limit <- 1
}

heatmap_data <- heatmap_data %>%
  dplyr::mutate(
    text_colour = ifelse(
      is.finite(NES) & abs(NES) >= 0.58 * colour_limit,
      "white",
      "black"
    )
  )

star_caption <- if (STAR_BY == "q_global_BH") {
  "Stars: global BH q-values (**** <0.0001; *** <0.001; ** <0.01; * <0.05)."
} else {
  "Stars: nominal fgsea P values (**** <0.0001; *** <0.001; ** <0.01; * <0.05)."
}

nes_heatmap <- ggplot(
  heatmap_data,
  aes(x = comparison, y = row_label, fill = NES)
) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(
    aes(label = cell_label, colour = text_colour),
    family = "Arial",
    fontface = "bold",
    size = 2.55,
    lineheight = 0.82,
    na.rm = TRUE
  ) +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-colour_limit, colour_limit),
    breaks = scales::breaks_pretty(n = 5),
    oob = scales::squish,
    na.value = "grey92",
    name = "NES"
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(24, "mm"),
      ticks = TRUE
    )
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  coord_fixed() +
  labs(
    title = METAPLASIA_SIGNATURE_TITLE,
    subtitle = paste0(
      "Preranked GSEA across organoid passage comparisons (",
      represented_signature_size_primary, "/",
      length(gastric_metaplasia_24), " genes represented)"
    ),
    x = NULL,
    y = NULL,
    caption = paste0(
      "Positive NES indicates enrichment in the numerator passage (for example, P2 in P2 vs P1); ",
      "negative NES indicates enrichment in P1. ",
      "Missing signature genes: ",
      ifelse(length(missing_signature_primary) == 0L, "none", paste(missing_signature_primary, collapse = ", ")),
      ".\n",
      star_caption
    )
  ) +
  theme_minimal(base_family = "Arial", base_size = 7) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      family = "Arial",
      size = 7,
      face = "bold",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      family = "Arial",
      size = 7,
      face = "bold",
      colour = "black"
    ),
    plot.title = element_text(
      family = "Arial",
      size = 9,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_text(
      family = "Arial",
      size = 7.5,
      hjust = 0
    ),
    plot.caption = element_text(
      family = "Arial",
      size = 6,
      hjust = 0,
      lineheight = 1.05,
      colour = "grey25",
      margin = margin(t = 2, unit = "mm")
    ),
    plot.caption.position = "plot",
    legend.position = "right",
    legend.direction = "vertical",
    legend.justification = "center",
    legend.box.just = "center",
    legend.title = element_text(
      family = "Arial",
      size = 7,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 1, unit = "mm")
    ),
    legend.text = element_text(
      family = "Arial",
      size = 6.5,
      margin = margin(l = 1, unit = "mm")
    ),
    legend.margin = margin(0, 0, 0, 1.5, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.spacing.y = grid::unit(0.5, "mm"),
    plot.margin = margin(5, 8, 6, 5, unit = "mm"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

figure_width_mm <- 183
figure_height_mm <- max(88, 50 + 13 * length(row_order))

primary_pdf_device <- if (capabilities("cairo")) {
  grDevices::cairo_pdf
} else {
  grDevices::pdf
}

ggsave(
  file.path(OUTPUT_DIR, "final_24_gene_passage_NES_heatmap.pdf"),
  nes_heatmap,
  device = primary_pdf_device,
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  bg = "white"
)

ggsave(
  file.path(OUTPUT_DIR, "final_24_gene_passage_NES_heatmap.png"),
  nes_heatmap,
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(OUTPUT_DIR, "final_24_gene_passage_NES_heatmap.tiff"),
  nes_heatmap,
  device = "tiff",
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

for (jpeg_extension in c("jpg", "jpeg")) {
  ggsave(
    file.path(
      OUTPUT_DIR,
      paste0("final_24_gene_passage_NES_heatmap.", jpeg_extension)
    ),
    nes_heatmap,
    device = "jpeg",
    width = figure_width_mm,
    height = figure_height_mm,
    units = "mm",
    dpi = 600,
    quality = 100,
    bg = "white"
  )
}

readr::write_csv(
  heatmap_data %>% dplyr::select(-text_colour),
  file.path(OUTPUT_DIR, "final_24_gene_passage_NES_heatmap_data.csv"),
  na = ""
)

# Cross-sectional SSL_HP versus H_ADJ_SSL_HP contrasts ----
cross_sectional_files <- file.path(
  DATA_DIR,
  unname(CROSS_SECTIONAL_FILE_NAMES)
)
names(cross_sectional_files) <- CROSS_SECTIONAL_PASSAGE_ORDER

missing_cross_sectional_files <- cross_sectional_files[
  !file.exists(cross_sectional_files)
]

if (length(missing_cross_sectional_files) > 0L) {
  stop(
    "The following required SSL_HP-versus-H_ADJ_SSL_HP DE tables are missing:\n",
    paste(missing_cross_sectional_files, collapse = "\n"),
    "\n\nThey must be located directly inside DATA_DIR:\n",
    DATA_DIR
  )
}

cross_sectional_manifest <- purrr::imap_dfr(
  cross_sectional_files,
  function(file, passage_label) {
    make_cross_sectional_meta(file, passage_label)
  }
)

readr::write_csv(
  cross_sectional_manifest,
  file.path(
    OUTPUT_DIR,
    "SSL_HP_vs_H_ADJ_SSL_HP_P1_P4_P8_contrast_manifest.csv"
  ),
  na = ""
)

message(
  "\nCross-sectional direction assumption: ",
  CROSS_SECTIONAL_CONTRAST,
  ". Positive NES/logFC is interpreted as higher in SSL_HP organoids."
)

cross_sectional_rank_objects <- purrr::imap(
  cross_sectional_files,
  function(file, passage_label) {
    read_and_rank_de_table(
      file = file,
      meta_override = make_cross_sectional_meta(file, passage_label)
    )
  }
)

cross_sectional_column_audit <- purrr::map_dfr(
  cross_sectional_rank_objects,
  "audit"
)

readr::write_csv(
  cross_sectional_column_audit,
  file.path(
    OUTPUT_DIR,
    "SSL_HP_vs_H_ADJ_SSL_HP_P1_P4_P8_input_column_audit.csv"
  ),
  na = ""
)

cross_sectional_ranking_methods <- unique(
  cross_sectional_column_audit$rank_method
)

if (length(cross_sectional_ranking_methods) != 1L) {
  stop(
    "Different ranking metrics would be used across the P1/P4/P8 tables: ",
    paste(cross_sectional_ranking_methods, collapse = " | "),
    ". All three comparisons must use the same ranking method."
  )
}

cross_sectional_gene_universes <- purrr::map(
  cross_sectional_rank_objects,
  ~ .x$ranked$Gene
)

cross_sectional_common_universe <- Reduce(
  intersect,
  cross_sectional_gene_universes
)

if (length(cross_sectional_common_universe) < MIN_RANKED_GENES) {
  stop(
    "The P1/P4/P8 common gene universe contains only ",
    length(cross_sectional_common_universe),
    " genes. Confirm that each input is a complete DE table rather than a ",
    "significance-filtered export."
  )
}

cross_sectional_gsea_results <- purrr::imap_dfr(
  cross_sectional_rank_objects,
  function(rank_object, passage_label) {
    passage_index <- match(
      as.character(passage_label),
      CROSS_SECTIONAL_PASSAGE_ORDER
    )

    if (is.na(passage_index)) {
      stop("Could not assign a deterministic seed to passage ", passage_label, ".")
    }

    run_one_gsea(
      rank_object = rank_object,
      universe = cross_sectional_common_universe,
      seed = GSEA_SEED + 5000L + passage_index
    )
  }
) %>%
  dplyr::mutate(
    comparison = factor(
      comparison,
      levels = CROSS_SECTIONAL_PASSAGE_ORDER
    ),
    q_cross_sectional_BH = stats::p.adjust(pval, method = "BH")
  ) %>%
  dplyr::arrange(comparison)

cross_sectional_star_p <- cross_sectional_gsea_results[[CROSS_SECTIONAL_STAR_BY]]

cross_sectional_gsea_results$stars <- significance_stars(
  cross_sectional_star_p
)
cross_sectional_gsea_results$star_basis <- paste0(
  "Benjamini-Hochberg q-value across the three cross-sectional ",
  "comparisons P1, P4 and P8"
)

readr::write_csv(
  cross_sectional_gsea_results,
  file.path(
    OUTPUT_DIR,
    "final_24_gene_SSL_HP_vs_H_ADJ_SSL_HP_P1_P4_P8_GSEA_results.csv"
  ),
  na = ""
)

cross_sectional_signature_presence <- purrr::map_dfr(
  cross_sectional_rank_objects,
  function(rank_object) {
    rank_object$meta %>%
      dplyr::slice(rep(1L, length(gastric_metaplasia_24))) %>%
      dplyr::mutate(
        Gene = gastric_metaplasia_24,
        present_in_this_DE_table = Gene %in% rank_object$ranked$Gene,
        present_in_cross_sectional_common_universe =
          Gene %in% cross_sectional_common_universe,
        used_in_cross_sectional_GSEA =
          present_in_cross_sectional_common_universe
      )
  }
)

readr::write_csv(
  cross_sectional_signature_presence,
  file.path(
    OUTPUT_DIR,
    "final_24_gene_SSL_HP_vs_H_ADJ_SSL_HP_P1_P4_P8_gene_presence.csv"
  ),
  na = ""
)

cross_sectional_represented_sizes <- unique(
  cross_sectional_gsea_results$represented_signature_size[
    is.finite(cross_sectional_gsea_results$represented_signature_size)
  ]
)

if (length(cross_sectional_represented_sizes) == 0L) {
  cross_sectional_represented_size <- NA_integer_
} else if (length(cross_sectional_represented_sizes) == 1L) {
  cross_sectional_represented_size <- as.integer(
    cross_sectional_represented_sizes
  )
} else {
  cross_sectional_represented_size <- min(
    cross_sectional_represented_sizes
  )
  warning(
    "The represented signature size varies across P1/P4/P8. ",
    "The subtitle reports the minimum; inspect the GSEA result table."
  )
}

cross_sectional_missing_genes <- sort(unique(unlist(
  strsplit(
    paste(
      cross_sectional_gsea_results$missing_signature_genes,
      collapse = ";"
    ),
    split = ";",
    fixed = TRUE
  )
)))
cross_sectional_missing_genes <- cross_sectional_missing_genes[
  !is.na(cross_sectional_missing_genes) &
    cross_sectional_missing_genes != ""
]

cross_sectional_heatmap_data <- cross_sectional_gsea_results %>%
  dplyr::mutate(
    Passage = factor(
      as.character(comparison),
      levels = CROSS_SECTIONAL_PASSAGE_ORDER
    ),
    Analysis = factor(
      "SSL_HP vs H_ADJ_SSL_HP",
      levels = "SSL_HP vs H_ADJ_SSL_HP"
    ),
    cell_label = dplyr::case_when(
      is.na(NES) ~ "NA",
      stars == "" ~ sprintf("%.2f", NES),
      TRUE ~ paste0(sprintf("%.2f", NES), "\n", stars)
    )
  )

cross_sectional_finite_nes <- cross_sectional_heatmap_data$NES[
  is.finite(cross_sectional_heatmap_data$NES)
]

if (length(cross_sectional_finite_nes) == 0L) {
  stop(
    "No finite cross-sectional NES values were produced. Inspect ",
    "final_24_gene_SSL_HP_vs_H_ADJ_SSL_HP_P1_P4_P8_GSEA_results.csv."
  )
}

cross_sectional_colour_limit <- max(abs(cross_sectional_finite_nes))
if (!is.finite(cross_sectional_colour_limit) ||
    cross_sectional_colour_limit == 0) {
  cross_sectional_colour_limit <- 1
}

cross_sectional_heatmap_data <- cross_sectional_heatmap_data %>%
  dplyr::mutate(
    text_colour = ifelse(
      is.finite(NES) &
        abs(NES) >= 0.58 * cross_sectional_colour_limit,
      "white",
      "black"
    )
  )

cross_sectional_subtitle_text <- stringr::str_wrap(
  paste0(
    "Preranked GSEA: SSL_HP vs H_ADJ_SSL_HP organoids at P1, P4 and P8 (",
    cross_sectional_represented_size,
    "/",
    length(gastric_metaplasia_24),
    " genes represented)"
  ),
  width = 100
)

cross_sectional_direction_caption <- stringr::str_wrap(
  paste0(
    "Positive NES indicates enrichment in SSL_HP relative to H_ADJ_SSL_HP ",
    "at the indicated passage; negative NES indicates enrichment in ",
    "H_ADJ_SSL_HP. Missing signature genes: ",
    ifelse(
      length(cross_sectional_missing_genes) == 0L,
      "none",
      paste(cross_sectional_missing_genes, collapse = ", ")
    ),
    "."
  ),
  width = 128
)

cross_sectional_significance_caption <- paste0(
  "Stars: BH q-values across P1, P4 and P8 ",
  "(**** <0.0001; *** <0.001; ** <0.01; * <0.05)."
)

cross_sectional_caption_text <- paste(
  cross_sectional_direction_caption,
  cross_sectional_significance_caption,
  sep = "\n"
)

cross_sectional_nes_heatmap <- ggplot2::ggplot(
  cross_sectional_heatmap_data,
  ggplot2::aes(x = Passage, y = Analysis, fill = NES)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.9) +
  ggplot2::geom_text(
    ggplot2::aes(label = cell_label, colour = text_colour),
    family = "Arial",
    fontface = "bold",
    size = 3.0,
    lineheight = 0.82,
    na.rm = TRUE
  ) +
  ggplot2::scale_colour_identity() +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -cross_sectional_colour_limit,
      cross_sectional_colour_limit
    ),
    breaks = scales::breaks_pretty(n = 5),
    oob = scales::squish,
    na.value = "grey92",
    name = "NES"
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(24, "mm"),
      ticks = TRUE
    )
  ) +
  ggplot2::scale_x_discrete(
    limits = CROSS_SECTIONAL_PASSAGE_ORDER,
    drop = FALSE
  ) +
  ggplot2::scale_y_discrete(drop = FALSE) +
  ggplot2::coord_fixed(ratio = 1.05, clip = "off") +
  ggplot2::labs(
    title = METAPLASIA_SIGNATURE_TITLE,
    subtitle = cross_sectional_subtitle_text,
    x = NULL,
    y = NULL,
    caption = cross_sectional_caption_text
  ) +
  ggplot2::theme_minimal(base_family = "Arial", base_size = 7) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      family = "Arial",
      size = 7,
      face = "bold",
      colour = "black"
    ),
    axis.text.y = ggplot2::element_text(
      family = "Arial",
      size = 7,
      face = "bold",
      colour = "black",
      margin = ggplot2::margin(r = 3, unit = "mm")
    ),
    plot.title = ggplot2::element_text(
      family = "Arial",
      size = 9,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = ggplot2::element_text(
      family = "Arial",
      size = 7.5,
      hjust = 0,
      lineheight = 1.05,
      margin = ggplot2::margin(b = 3, unit = "mm")
    ),
    plot.caption = ggplot2::element_text(
      family = "Arial",
      size = 6,
      hjust = 0,
      lineheight = 1.08,
      colour = "grey25",
      margin = ggplot2::margin(t = 3, unit = "mm")
    ),
    plot.caption.position = "plot",
    legend.position = "right",
    legend.direction = "vertical",
    legend.justification = "center",
    legend.box.just = "center",
    legend.title = ggplot2::element_text(
      family = "Arial",
      size = 7,
      face = "bold",
      hjust = 0.5,
      margin = ggplot2::margin(b = 1, unit = "mm")
    ),
    legend.text = ggplot2::element_text(
      family = "Arial",
      size = 6.5,
      margin = ggplot2::margin(l = 1, unit = "mm")
    ),
    legend.margin = ggplot2::margin(0, 0, 0, 2, unit = "mm"),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
    legend.spacing.y = grid::unit(0.5, "mm"),
    plot.margin = ggplot2::margin(6, 14, 10, 8, unit = "mm"),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA)
  )

cross_sectional_width_mm <- 183
cross_sectional_height_mm <- 92
cross_sectional_stem <- file.path(
  OUTPUT_DIR,
  "final_24_gene_SSL_HP_vs_H_ADJ_SSL_HP_P1_P4_P8_NES_heatmap"
)

ggplot2::ggsave(
  paste0(cross_sectional_stem, ".pdf"),
  cross_sectional_nes_heatmap,
  device = primary_pdf_device,
  width = cross_sectional_width_mm,
  height = cross_sectional_height_mm,
  units = "mm",
  bg = "white",
  limitsize = FALSE
)

ggplot2::ggsave(
  paste0(cross_sectional_stem, ".png"),
  cross_sectional_nes_heatmap,
  width = cross_sectional_width_mm,
  height = cross_sectional_height_mm,
  units = "mm",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

ggplot2::ggsave(
  paste0(cross_sectional_stem, ".tiff"),
  cross_sectional_nes_heatmap,
  device = "tiff",
  width = cross_sectional_width_mm,
  height = cross_sectional_height_mm,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white",
  limitsize = FALSE
)

for (jpeg_extension in c("jpg", "jpeg")) {
  ggplot2::ggsave(
    paste0(cross_sectional_stem, ".", jpeg_extension),
    cross_sectional_nes_heatmap,
    device = "jpeg",
    width = cross_sectional_width_mm,
    height = cross_sectional_height_mm,
    units = "mm",
    dpi = 600,
    quality = 100,
    bg = "white",
    limitsize = FALSE
  )
}

readr::write_csv(
  cross_sectional_heatmap_data %>%
    dplyr::select(-text_colour),
  file.path(
    OUTPUT_DIR,
    "final_24_gene_SSL_HP_vs_H_ADJ_SSL_HP_P1_P4_P8_NES_heatmap_data.csv"
  ),
  na = ""
)

cross_sectional_gene_logfc <- purrr::imap_dfr(
  cross_sectional_rank_objects,
  function(rank_object, passage_label) {
    tibble::tibble(
      Gene = gastric_metaplasia_24,
      Signature_order = seq_along(gastric_metaplasia_24)
    ) %>%
      dplyr::left_join(
        rank_object$ranked %>%
          dplyr::select(Gene, logFC, logCPM, PValue, FDR),
        by = "Gene"
      ) %>%
      dplyr::mutate(
        Passage = factor(
          passage_label,
          levels = CROSS_SECTIONAL_PASSAGE_ORDER
        ),
        Present_in_DE_table = is.finite(logFC),
        logFC_plot = pmax(
          pmin(logFC, CROSS_SECTIONAL_GENE_LOGFC_LIMIT),
          -CROSS_SECTIONAL_GENE_LOGFC_LIMIT
        )
      )
  }
)

readr::write_csv(
  cross_sectional_gene_logfc,
  file.path(
    OUTPUT_DIR,
    "final_24_gene_SSL_HP_vs_H_ADJ_SSL_HP_P1_P4_P8_gene_log2FC_data.csv"
  ),
  na = ""
)

cross_sectional_gene_heatmap <- ggplot2::ggplot(
  cross_sectional_gene_logfc %>%
    dplyr::mutate(
      Gene = factor(
        Gene,
        levels = rev(gastric_metaplasia_24)
      )
    ),
  ggplot2::aes(x = Passage, y = Gene, fill = logFC_plot)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.28) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -CROSS_SECTIONAL_GENE_LOGFC_LIMIT,
      CROSS_SECTIONAL_GENE_LOGFC_LIMIT
    ),
    oob = scales::squish,
    na.value = "grey90",
    name = "log2FC"
  ) +
  ggplot2::scale_x_discrete(
    limits = CROSS_SECTIONAL_PASSAGE_ORDER,
    drop = FALSE
  ) +
  ggplot2::labs(
    title = METAPLASIA_SIGNATURE_TITLE,
    subtitle = "Gene-level differential expression: SSL_HP vs H_ADJ_SSL_HP organoids",
    x = NULL,
    y = NULL,
    caption = paste0(
      "Positive log2FC indicates higher expression in SSL_HP; negative log2FC ",
      "indicates higher expression in H_ADJ_SSL_HP. Values are truncated at ",
      "±", CROSS_SECTIONAL_GENE_LOGFC_LIMIT, " for visualisation; missing genes are grey."
    )
  ) +
  ggplot2::theme_minimal(base_family = "Arial", base_size = 7) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      family = "Arial",
      size = 7,
      face = "bold",
      colour = "black"
    ),
    axis.text.y = ggplot2::element_text(
      family = "Arial",
      size = 6.5,
      colour = "black"
    ),
    plot.title = ggplot2::element_text(
      family = "Arial",
      size = 9,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = ggplot2::element_text(
      family = "Arial",
      size = 7.5,
      hjust = 0
    ),
    plot.caption = ggplot2::element_text(
      family = "Arial",
      size = 6,
      hjust = 0,
      lineheight = 1.05,
      colour = "grey25"
    ),
    plot.caption.position = "plot",
    legend.position = "right",
    legend.title = ggplot2::element_text(
      family = "Arial",
      size = 7,
      face = "bold"
    ),
    legend.text = ggplot2::element_text(
      family = "Arial",
      size = 6.5
    ),
    plot.margin = ggplot2::margin(5, 8, 6, 5, unit = "mm"),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA)
  )

cross_sectional_gene_stem <- file.path(
  OUTPUT_DIR,
  "final_24_gene_SSL_HP_vs_H_ADJ_SSL_HP_P1_P4_P8_gene_log2FC_heatmap"
)

ggplot2::ggsave(
  paste0(cross_sectional_gene_stem, ".pdf"),
  cross_sectional_gene_heatmap,
  device = primary_pdf_device,
  width = 120,
  height = 145,
  units = "mm",
  bg = "white"
)

ggplot2::ggsave(
  paste0(cross_sectional_gene_stem, ".png"),
  cross_sectional_gene_heatmap,
  width = 120,
  height = 145,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggplot2::ggsave(
  paste0(cross_sectional_gene_stem, ".tiff"),
  cross_sectional_gene_heatmap,
  device = "tiff",
  width = 120,
  height = 145,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

for (jpeg_extension in c("jpg", "jpeg")) {
  ggplot2::ggsave(
    paste0(cross_sectional_gene_stem, ".", jpeg_extension),
    cross_sectional_gene_heatmap,
    device = "jpeg",
    width = 120,
    height = 145,
    units = "mm",
    dpi = 600,
    quality = 100,
    bg = "white"
  )
}

cross_sectional_methods_text <- c(
  paste0(
    "Preranked GSEA was performed separately at P1, P4 and P8 using the prespecified ",
    METAPLASIA_SIGNATURE_TITLE,
    " (24 genes)."
  ),
  "Each edgeR contrast was interpreted as SSL_HP minus H_ADJ_SSL_HP at the corresponding passage.",
  paste0(
    "Genes were ranked by ",
    cross_sectional_ranking_methods,
    "; no logCPM, nominal P-value or FDR threshold was applied before ranking, and the primary analysis used the intersection of ranked gene symbols represented across all three DE tables."
  ),
  paste0(
    cross_sectional_represented_size,
    " of 24 signature genes were represented in the common ranked gene universe."
  ),
  "Nominal enrichment P values were adjusted across the three passage-specific comparisons using the Benjamini-Hochberg method.",
  "Positive NES values indicate enrichment in SSL_HP relative to H_ADJ_SSL_HP."
)

writeLines(
  cross_sectional_methods_text,
  file.path(
    OUTPUT_DIR,
    "THESIS_METHODS_SSL_HP_vs_H_ADJ_SSL_HP_P1_P4_P8.txt"
  )
)

cross_sectional_legend_text <- paste0(
  "Preranked GSEA of the ",
  METAPLASIA_SIGNATURE_TITLE,
  " in SSL_HP relative to H_ADJ_SSL_HP organoids at P1, P4 and P8. ",
  cross_sectional_represented_size,
  " of the 24 genes were represented in the common ranked gene universe",
  ifelse(
    length(cross_sectional_missing_genes) > 0L,
    paste0(
      " (missing: ",
      paste(cross_sectional_missing_genes, collapse = ", "),
      ")"
    ),
    ""
  ),
  ". Positive normalised enrichment scores (NES) indicate enrichment in SSL_HP relative to H_ADJ_SSL_HP at the indicated passage. Values within cells show NES; asterisks denote Benjamini-Hochberg-adjusted q values across the three comparisons (****q<0.0001, ***q<0.001, **q<0.01, *q<0.05)."
)

writeLines(
  cross_sectional_legend_text,
  file.path(
    OUTPUT_DIR,
    "THESIS_FIGURE_LEGEND_SSL_HP_vs_H_ADJ_SSL_HP_P1_P4_P8.txt"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "passage_GSEA_sessionInfo.txt")
)

summary_to_print <- gsea_results %>%
  dplyr::select(
    line_id,
    condition,
    comparison,
    represented_signature_size,
    NES,
    pval,
    q_global_BH,
    stars,
    status
  ) %>%
  dplyr::arrange(line_id, condition, comparison)

print(summary_to_print, n = Inf)

cross_sectional_summary_to_print <- cross_sectional_gsea_results %>%
  dplyr::select(
    comparison,
    represented_signature_size,
    NES,
    pval,
    q_cross_sectional_BH,
    stars,
    status
  ) %>%
  dplyr::arrange(comparison)

message("\nSSL_HP versus H_ADJ_SSL_HP at P1, P4 and P8:")
print(cross_sectional_summary_to_print, n = Inf)

message("\nLongitudinal ranking method: ", ranking_methods)
message(
  "Longitudinal gene universe: ",
  if (USE_COMMON_GENE_UNIVERSE) "common intersection across all longitudinal DE tables" else "per table"
)
message(
  "Cross-sectional ranking method: ",
  cross_sectional_ranking_methods
)
message(
  "Cross-sectional gene universe: common intersection across P1, P4 and P8 DE tables"
)
message("Outputs written to:\n", OUTPUT_DIR)

# Gene-universe and ranking sensitivities ----
make_pvalue_rank_object <- function(rank_object) {
  dat <- rank_object$ranked_uncollapsed

  if (is.null(dat) || nrow(dat) == 0L) {
    stop("No uncollapsed ranked rows were retained for P-value sensitivity.")
  }

  positive_p <- dat$PValue[is.finite(dat$PValue) & dat$PValue > 0]
  min_positive_p <- if (length(positive_p) > 0L) min(positive_p) else 1e-300

  dat <- dat %>%
    dplyr::mutate(
      PValue_safe = dplyr::case_when(
        !is.finite(PValue) ~ NA_real_,
        PValue <= 0 ~ min_positive_p * 0.1,
        TRUE ~ PValue
      ),
      rank_metric = sign(logFC) * -log10(PValue_safe)
    ) %>%
    dplyr::filter(is.finite(rank_metric), !is.na(Gene), Gene != "") %>%
    dplyr::group_by(Gene) %>%
    dplyr::arrange(
      dplyr::desc(abs(rank_metric)),
      PValue_safe,
      dplyr::desc(abs(logFC)),
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(n = 1L) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(rank_metric), Gene)

  list(
    meta = rank_object$meta,
    ranked = dat,
    audit = rank_object$audit
  )
}

sensitivity_results_list <- list(
  PRIMARY_SIGNED_SQRT_LR_COMMON_UNIVERSE = gsea_results %>%
    dplyr::mutate(
      analysis_variant = "Primary: signed sqrt(LR), common universe",
      ranking_metric_sensitivity = "sign(logFC) * sqrt(LR)",
      universe_sensitivity = "Common intersection"
    )
)

if (RUN_PVALUE_RANK_SENSITIVITY) {
  pvalue_rank_objects <- purrr::map(rank_objects, make_pvalue_rank_object)
  pvalue_gene_universes <- purrr::map(pvalue_rank_objects, ~ .x$ranked$Gene)
  pvalue_common_universe <- Reduce(intersect, pvalue_gene_universes)

  if (length(pvalue_common_universe) < MIN_RANKED_GENES) {
    stop(
      "The P-value sensitivity common universe contains only ",
      length(pvalue_common_universe), " genes."
    )
  }

  pvalue_results <- purrr::imap_dfr(pvalue_rank_objects, function(rank_object, i) {
    run_one_gsea(
      rank_object = rank_object,
      universe = pvalue_common_universe,
      seed = GSEA_SEED + 1000L + as.integer(i)
    )
  }) %>%
    dplyr::mutate(
      q_global_BH = stats::p.adjust(pval, method = "BH"),
      analysis_variant = "Sensitivity: signed -log10(P), common universe",
      ranking_metric_sensitivity = "sign(logFC) * -log10(PValue)",
      universe_sensitivity = "Common intersection"
    )

  sensitivity_results_list$PVALUE_RANK_COMMON_UNIVERSE <- pvalue_results
}

if (RUN_PER_TABLE_UNIVERSE_SENSITIVITY) {
  per_table_results <- purrr::imap_dfr(rank_objects, function(rank_object, i) {
    run_one_gsea(
      rank_object = rank_object,
      universe = rank_object$ranked$Gene,
      seed = GSEA_SEED + 2000L + as.integer(i)
    )
  }) %>%
    dplyr::mutate(
      q_global_BH = stats::p.adjust(pval, method = "BH"),
      analysis_variant = "Sensitivity: signed sqrt(LR), per-table universe",
      ranking_metric_sensitivity = "sign(logFC) * sqrt(LR)",
      universe_sensitivity = "Each DE table"
    )

  sensitivity_results_list$SIGNED_SQRT_LR_PER_TABLE_UNIVERSE <- per_table_results
}

sensitivity_results <- dplyr::bind_rows(sensitivity_results_list) %>%
  dplyr::mutate(
    NES_direction = dplyr::case_when(
      is.na(NES) ~ NA_character_,
      NES > 0 ~ "Positive",
      NES < 0 ~ "Negative",
      TRUE ~ "Zero"
    )
  )

readr::write_csv(
  sensitivity_results,
  file.path(OUTPUT_DIR, "final_24_gene_passage_GSEA_ALL_PRIMARY_AND_SENSITIVITY_RESULTS.csv"),
  na = ""
)

sensitivity_wide <- sensitivity_results %>%
  dplyr::select(
    line_id, condition, comparison, numerator_passage, reference_passage,
    analysis_variant, NES, pval, q_global_BH,
    represented_signature_size, missing_signature_genes, NES_direction
  ) %>%
  tidyr::pivot_wider(
    names_from = analysis_variant,
    values_from = c(
      NES, pval, q_global_BH, represented_signature_size,
      missing_signature_genes, NES_direction
    ),
    names_glue = "{.value}__{analysis_variant}"
  )

variant_names <- unique(sensitivity_results$analysis_variant)

primary_direction_col <- paste0(
  "NES_direction__",
  "Primary: signed sqrt(LR), common universe"
)

if (primary_direction_col %in% names(sensitivity_wide)) {
  sensitivity_direction_cols <- grep(
    "^NES_direction__",
    names(sensitivity_wide),
    value = TRUE
  )

  sensitivity_wide <- sensitivity_wide %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      all_available_variants_same_direction = {
        dirs <- c_across(dplyr::all_of(sensitivity_direction_cols))
        dirs <- dirs[!is.na(dirs) & dirs != ""]
        length(unique(dirs)) <= 1L
      }
    ) %>%
    dplyr::ungroup()
}

readr::write_csv(
  sensitivity_wide,
  file.path(OUTPUT_DIR, "final_24_gene_passage_GSEA_sensitivity_comparison.csv"),
  na = ""
)

sensitivity_summary <- sensitivity_results %>%
  dplyr::group_by(analysis_variant, ranking_metric_sensitivity, universe_sensitivity) %>%
  dplyr::summarise(
    n_comparisons = dplyr::n(),
    n_positive_NES = sum(NES > 0, na.rm = TRUE),
    n_negative_NES = sum(NES < 0, na.rm = TRUE),
    n_global_BH_q_lt_0_05 = sum(q_global_BH < 0.05, na.rm = TRUE),
    minimum_NES = min(NES, na.rm = TRUE),
    maximum_NES = max(NES, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  sensitivity_summary,
  file.path(OUTPUT_DIR, "final_24_gene_passage_GSEA_sensitivity_summary.csv"),
  na = ""
)

sensitivity_plot_data <- sensitivity_results %>%
  dplyr::mutate(
    comparison = factor(comparison, levels = comparison_order),
    condition = factor(condition, levels = condition_priority),
    analysis_variant = factor(analysis_variant, levels = variant_names),
    sensitivity_row = interaction(
      condition, analysis_variant,
      sep = " | ", drop = TRUE, lex.order = TRUE
    ),
    cell_label = ifelse(is.finite(NES), sprintf("%.2f", NES), "")
  )

sensitivity_limit <- max(abs(sensitivity_plot_data$NES), na.rm = TRUE)
if (!is.finite(sensitivity_limit) || sensitivity_limit == 0) sensitivity_limit <- 1

sensitivity_heatmap <- ggplot(
  sensitivity_plot_data,
  aes(x = comparison, y = sensitivity_row, fill = NES)
) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(
    aes(label = cell_label),
    family = "Arial", fontface = "bold", size = 2.2
  ) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0,
    limits = c(-sensitivity_limit, sensitivity_limit),
    oob = scales::squish,
    na.value = "grey92",
    name = "NES"
  ) +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = METAPLASIA_SIGNATURE_TITLE,
    subtitle = paste0(
      "Passage GSEA robustness analysis: primary signed-LR/common-universe ",
      "result and pre-specified sensitivity analyses"
    ),
    x = NULL, y = NULL,
    caption = "This robustness heatmap is supplementary/QC; retain the primary heatmap as the thesis figure."
  ) +
  theme_minimal(base_family = "Arial", base_size = 7) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      family = "Arial", size = 7, face = "bold",
      angle = 45, hjust = 1, vjust = 1
    ),
    axis.text.y = element_text(
      family = "Arial", size = 6.4, colour = "black"
    ),
    plot.title = element_text(
      family = "Arial", size = 9, face = "bold", hjust = 0
    ),
    plot.subtitle = element_text(family = "Arial", size = 7.5, hjust = 0),
    plot.caption = element_text(
      family = "Arial", size = 6, hjust = 0, colour = "grey25"
    ),
    legend.position = "right",
    plot.margin = margin(5, 8, 6, 5, unit = "mm")
  )

sensitivity_height_mm <- max(105, 58 + 9 * length(unique(sensitivity_plot_data$sensitivity_row)))
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

ggsave(
  file.path(OUTPUT_DIR, "final_24_gene_passage_NES_sensitivity_heatmap.pdf"),
  sensitivity_heatmap,
  device = pdf_device,
  width = 183,
  height = sensitivity_height_mm,
  units = "mm",
  bg = "white"
)

ggsave(
  file.path(OUTPUT_DIR, "final_24_gene_passage_NES_sensitivity_heatmap.png"),
  sensitivity_heatmap,
  width = 183,
  height = sensitivity_height_mm,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(OUTPUT_DIR, "final_24_gene_passage_NES_sensitivity_heatmap.tiff"),
  sensitivity_heatmap,
  device = "tiff",
  width = 183,
  height = sensitivity_height_mm,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# Methods, figure legend and session records ----
primary_n_sig <- sum(gsea_results$q_global_BH < 0.05, na.rm = TRUE)
primary_n_positive <- sum(gsea_results$NES > 0, na.rm = TRUE)
primary_n_total <- nrow(gsea_results)

methods_text <- c(
  paste0(
    "Preranked gene-set enrichment analysis was performed using fgseaMultilevel and the prespecified ",
    METAPLASIA_SIGNATURE_TITLE,
    " (24 genes)."
  ),
  "For each edgeR passage contrast, genes were ranked by sign(logFC) × sqrt(LR), which combines the direction of differential expression with the magnitude of the likelihood-ratio statistic.",
  "No logCPM, nominal P-value or FDR threshold was applied before ranking.",
  if (USE_COMMON_GENE_UNIVERSE) {
    "The primary analysis used the intersection of ranked gene symbols represented across all passage comparisons to maintain a common tested-gene universe."
  } else {
    "The primary analysis used the ranked gene universe available in each passage comparison."
  },
  paste0(
    represented_signature_size_primary, " of the 24 signature genes were represented in the primary ranked universe."
  ),
  "Nominal enrichment P values were adjusted across all passage comparisons using the Benjamini-Hochberg method.",
  "Positive NES values indicate enrichment in the numerator passage relative to P1.",
  "As robustness checks, the analysis was repeated using a signed -log10(PValue) ranking and using each comparison's own ranked-gene universe."
)

writeLines(
  methods_text,
  file.path(OUTPUT_DIR, "THESIS_METHODS_TEXT.txt")
)

figure_legend_text <- paste0(
  "Preranked GSEA of the ",
  METAPLASIA_SIGNATURE_TITLE,
  " across organoid passage comparisons. ",
  represented_signature_size_primary, " of the 24 genes were represented in the common ranked gene universe",
  ifelse(
    length(missing_signature_primary) > 0L,
    paste0(" (missing: ", paste(missing_signature_primary, collapse = ", "), ")"),
    ""
  ),
  ". Positive normalised enrichment scores (NES) indicate enrichment in the later numerator passage relative to P1. ",
  "Values within cells show NES; asterisks denote Benjamini-Hochberg-adjusted q values across all passage comparisons ",
  "(****q<0.0001, ***q<0.001, **q<0.01, *q<0.05)."
)

writeLines(
  figure_legend_text,
  file.path(OUTPUT_DIR, "THESIS_FIGURE_LEGEND_TEMPLATE.txt")
)

analysis_summary_text <- c(
  toupper(paste0(METAPLASIA_SIGNATURE_TITLE, ": PASSAGE GSEA SUMMARY")),
  "=======================================",
  "",
  paste0("Primary comparisons analysed: ", primary_n_total),
  paste0("Positive NES values: ", primary_n_positive, "/", primary_n_total),
  paste0("Global BH q < 0.05: ", primary_n_sig, "/", primary_n_total),
  paste0(
    "Signature genes represented in primary universe: ",
    represented_signature_size_primary, "/", length(gastric_metaplasia_24)
  ),
  paste0(
    "Missing signature genes: ",
    ifelse(length(missing_signature_primary) == 0L, "none", paste(missing_signature_primary, collapse = ", "))
  ),
  "",
  "Primary analysis:",
  "  Ranking = sign(logFC) * sqrt(LR)",
  paste0(
    "  Universe = ",
    ifelse(USE_COMMON_GENE_UNIVERSE, "common intersection across DE tables", "per DE table")
  ),
  "  No logCPM/P/FDR prefiltering",
  "  Multiplicity = BH across all passage comparisons",
  "",
  "Interpretation limitation:",
  "  These are separate Pn-vs-P1 enrichment contrasts and do not constitute a formal linear passage-trend or passage-by-organoid-type interaction test.",
  "  The direction is valid only if the original edgeR contrasts were truly Pn - P1."
)

writeLines(
  analysis_summary_text,
  file.path(OUTPUT_DIR, "ANALYSIS_SUMMARY_AND_INTERPRETATION.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "passage_GSEA_sessionInfo.txt")
)

message("\nThesis robustness analyses complete.")
message("Primary result file: final_24_gene_passage_GSEA_results.csv")
message("Sensitivity comparison: final_24_gene_passage_GSEA_sensitivity_comparison.csv")
message("Corrected gene-presence audit: final_24_gene_presence_by_DE_table_AND_common_universe.csv")
message("Thesis methods: THESIS_METHODS_TEXT.txt")
message("Figure legend: THESIS_FIGURE_LEGEND_TEMPLATE.txt")
message(
  "Cross-sectional NES heatmap: ",
  "final_24_gene_SSL_HP_vs_H_ADJ_SSL_HP_P1_P4_P8_NES_heatmap.*"
)
message(
  "Cross-sectional gene-level heatmap: ",
  "final_24_gene_SSL_HP_vs_H_ADJ_SSL_HP_P1_P4_P8_gene_log2FC_heatmap.*"
)
