# Sensitivity of gastric-signature membership to the recorded selection thresholds.
# Run from Tissue_Bulk_RNAseq after configuring config/config.R.
# The baseline must reconstruct the fixed 24-gene signature before sensitivity settings are evaluated.

source("config/config.R")

required_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
  "scales",
  "AnnotationDbi",
  "org.Hs.eg.db"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

set.seed(1)

base_dir <- PATH_BULK

params <- list(
  base_dir = base_dir,

  gtex_file = file.path(
    base_dir,
    "GM LIST",
    "GTEx_v10_TPM_by_tissue.csv"
  ),

  hpa_file = file.path(
    base_dir,
    "GM LIST",
    "rna_tissue_hpa.tsv"
  ),

  combined_de_file = file.path(
    base_dir,
    "Individual DE Tables",
    "COMBINED_DE_TABLE_ALL_COMPARISONS.csv"
  ),

  antrum_prefix = "ANTRUM_v_H_COL_DE_results",

  ssl_hp_hcol_de_file = file.path(
    base_dir,
    "Re-run GSEA",
    "CGG014_Maier_SSL_HP_v_H_COL_DE_results.txt"
  ),

  primary_signature_file_candidates = c(
    file.path(
      base_dir,
      "GM_publication_ready",
      "tables",
      "FINAL_24_GENE_GASTRIC_METAPLASIA_SIGNATURE.txt"
    ),
    file.path(
      base_dir,
      "GM_publication_ready",
      "tables",
      "FINAL_SSL_METAPLASIA_LIST.txt"
    ),
    file.path(
      base_dir,
      "GM_publication_ready",
      "CPM_TPM_signature_scoring_final",
      "tables",
      "FINAL_SSL_METAPLASIA_LIST.txt"
    )
  ),

  out_dir = file.path(
    base_dir,
    "GM_publication_ready",
    "FULL_SIGNATURE_THRESHOLD_SENSITIVITY_AUDIT_V3"
  ),

  pseudocount = 0.1,

  TAU_MIN              = 0.85,
  GTEX_STOMACH_TPM_MIN = 2,
  GTEX_COLON_TPM_MAX   = 1,
  GTEX_LOG2FC_MIN      = 2.5,

  HPA_STOMACH_NTPM_MIN = 2,
  HPA_COLON_NTPM_MAX   = 10,
  HPA_LOG2FC_MIN       = 3,
  HPA_MISSING_RULE     = "original_missing_as_zero",

  ANTRUM_LOG2FC_MIN = 1.5,
  ANTRUM_FDR_MAX    = 0.05,

  SSL_HP_LOG2FC_MIN = 0,
  SSL_HP_FDR_MAX    = 0.05,

  use_raw_p_if_adjusted_p_missing = FALSE,

  stop_if_baseline_reconstruction_mismatch = TRUE,
  stop_if_primary_file_differs_from_embedded_set = TRUE,

  write_individual_txt_for_every_setting = FALSE,

  established_primary_signature = c(
    "AQP5", "ANXA10", "MUC5AC", "VSIG1", "TCN1", "MUC6",
    "PSAPL1", "CLDN18", "SFTA2", "MUCL3", "GATA4", "KLK11",
    "DUOXA1", "SULT1C2", "IRX3", "AADAC", "FOXQ1", "IGFALS",
    "NR0B2", "C13ORF46", "FOLR1", "LYPD6B", "TFR2", "SPTSSB"
  )
)
params$primary_signature_file_candidates <- PATH_GASTRIC_SIGNATURE
params$out_dir <- OUTPUT_SENSITIVITY
ensure_new_output(params$out_dir)

sensitivity_values <- list(
  gtex_tau_min          = c(0.80, 0.85, 0.90),
  gtex_stomach_min      = c(1, 2, 5),
  gtex_colon_max        = c(0.5, 1, 2),
  gtex_log2fc_min       = c(2, log2(5), 2.5, 3),

  hpa_stomach_min       = c(1, 2, 5),
  hpa_colon_max         = c(1, 5, 10, 20, Inf),
  hpa_log2fc_min        = c(2, 2.5, 3),
  hpa_missing_rule      = c("original_missing_as_zero", "complete_case"),

  antrum_log2fc_min     = c(0, 1, 1.5, 2, log2(5), 2.5, 3),
  antrum_fdr_max        = c(0.01, 0.05, 0.10),

  ssl_hp_log2fc_min     = c(0, 0.5, 1, 1.5, 2),
  ssl_hp_fdr_max        = c(0.01, 0.05, 0.10)
)

dir.create(params$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(params$out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(params$out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(params$out_dir, "qc"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(params$out_dir, "selected_gene_lists"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(params$out_dir, "all_setting_gene_lists"), showWarnings = FALSE, recursive = TRUE)

to_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

safe_max <- function(x) {
  x <- to_num(x)
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

valid_gene_symbol <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x)))
}

clean_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[!valid_gene_symbol(x)] <- NA_character_
  x
}

tau_row <- function(v) {
  v <- to_num(v)
  n_nonmissing <- sum(!is.na(v))
  if (n_nonmissing < 2L) return(NA_real_)

  maximum <- max(v, na.rm = TRUE)
  if (!is.finite(maximum) || maximum <= 0) return(NA_real_)

  scaled <- v / maximum
  sum(1 - scaled, na.rm = TRUE) / (n_nonmissing - 1)
}

pick_first_existing <- function(candidates, object) {
  hits <- candidates[candidates %in% names(object)]
  if (length(hits) == 0L) return(NA_character_)
  hits[1]
}

require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " was not found:\n  ", path)
  }
}

require_columns <- function(df, required, label) {
  missing_columns <- setdiff(required, names(df))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required column(s) in ", label, ": ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

read_table_flexible <- function(path) {
  require_file(path, "Input file")
  data.table::fread(
    path,
    data.table = FALSE,
    fill = TRUE,
    check.names = FALSE
  )
}

write_csv_overwrite <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  temp_path <- paste0(path, ".tmp")

  if (file.exists(temp_path)) unlink(temp_path, force = TRUE)
  data.table::fwrite(df, temp_path)

  if (file.exists(path)) {
    removed <- suppressWarnings(file.remove(path))
    if (!removed) {
      stop(
        "Could not overwrite:\n  ", path,
        "\nClose the file and allow OneDrive to finish syncing, then rerun."
      )
    }
  }

  moved <- file.rename(temp_path, path)
  if (!moved) stop("Could not move temporary output into place:\n  ", path)
}

write_text_overwrite <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  temp_path <- paste0(path, ".tmp")

  if (file.exists(temp_path)) unlink(temp_path, force = TRUE)
  writeLines(as.character(x), con = temp_path, useBytes = TRUE)

  if (file.exists(path)) {
    removed <- suppressWarnings(file.remove(path))
    if (!removed) {
      stop(
        "Could not overwrite:\n  ", path,
        "\nClose the file and allow OneDrive to finish syncing, then rerun."
      )
    }
  }

  moved <- file.rename(temp_path, path)
  if (!moved) stop("Could not move temporary output into place:\n  ", path)
}

write_gene_list <- function(x, path) {
  x <- unique(clean_symbol(x))
  x <- x[valid_gene_symbol(x)]
  write_text_overwrite(x, path)
}

format_number_tag <- function(x) {
  if (length(x) == 0L || is.na(x)) return("NA")
  if (is.numeric(x) && is.infinite(x)) return("INF")
  out <- as.character(signif(as.numeric(x), 8))
  out <- gsub("-", "m", out, fixed = TRUE)
  out <- gsub("\\.", "p", out)
  out
}

safe_file_tag <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

collapse_gene_string <- function(x) {
  x <- unique(clean_symbol(x))
  x <- x[valid_gene_symbol(x)]
  if (length(x) == 0L) return("")
  paste(x, collapse = ";")
}

order_gene_set <- function(genes, primary_order) {
  genes <- unique(clean_symbol(genes))
  genes <- genes[valid_gene_symbol(genes)]
  c(
    primary_order[primary_order %in% genes],
    sort(setdiff(genes, primary_order))
  )
}

standardize_gene_vector <- function(x) {
  raw <- trimws(as.character(x))
  no_version <- sub("\\..*$", "", raw)

  nonmissing <- no_version[valid_gene_symbol(no_version)]
  ensembl_fraction <- if (length(nonmissing) == 0L) {
    0
  } else {
    mean(grepl("^ENSG[0-9]+$", toupper(nonmissing)))
  }

  if (is.finite(ensembl_fraction) && ensembl_fraction > 0.5) {
    mapping_keys <- unique(no_version[valid_gene_symbol(no_version)])

    mapped <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = mapping_keys,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )

    symbol <- unname(mapped[no_version])
    symbol <- clean_symbol(symbol)

    missing_map <- !valid_gene_symbol(symbol) & valid_gene_symbol(no_version)
    symbol[missing_map] <- clean_symbol(no_version[missing_map])
    symbol
  } else {
    clean_symbol(raw)
  }
}

standardize_de_table <- function(path) {
  df <- read_table_flexible(path)

  gene_col <- pick_first_existing(
    c(
      "gene_symbol", "GeneSymbol", "SYMBOL", "symbol", "gene", "Gene",
      "external_gene_name", "gene_name", "GeneName", "id", "ID",
      "gene_id", "Geneid", "ENSEMBL", "ensembl_gene_id"
    ),
    df
  )
  if (is.na(gene_col)) gene_col <- names(df)[1]

  log2fc_col <- pick_first_existing(
    c(
      "log2FoldChange", "log2FC", "logFC", "avg_log2FC", "LFC",
      "log2_fold_change"
    ),
    df
  )

  adjusted_p_col <- pick_first_existing(
    c(
      "padj", "FDR", "fdr", "qvalue", "q_value", "adj.P.Val",
      "p.adjust", "p_adj", "P_adj"
    ),
    df
  )

  raw_p_col <- pick_first_existing(
    c("pvalue", "p_value", "P.Value", "pval", "p", "P"),
    df
  )

  if (is.na(log2fc_col)) {
    stop("Could not identify the log2 fold-change column in:\n  ", path)
  }

  out <- data.frame(
    Gene = standardize_gene_vector(df[[gene_col]]),
    log2FC = to_num(df[[log2fc_col]]),
    adjusted_p = if (!is.na(adjusted_p_col)) {
      to_num(df[[adjusted_p_col]])
    } else {
      rep(NA_real_, nrow(df))
    },
    raw_p = if (!is.na(raw_p_col)) {
      to_num(df[[raw_p_col]])
    } else {
      rep(NA_real_, nrow(df))
    },
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(valid_gene_symbol(Gene), is.finite(log2FC)) %>%
    dplyr::group_by(Gene) %>%
    dplyr::arrange(
      dplyr::desc(abs(log2FC)),
      ifelse(is.finite(adjusted_p), adjusted_p, Inf),
      .by_group = TRUE
    ) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  attr(out, "gene_column_used") <- gene_col
  attr(out, "log2fc_column_used") <- log2fc_col
  attr(out, "adjusted_p_column_used") <- adjusted_p_col
  attr(out, "raw_p_column_used") <- raw_p_col
  out
}

load_primary_signature <- function(params) {
  embedded <- unique(clean_symbol(params$established_primary_signature))
  embedded <- embedded[valid_gene_symbol(embedded)]

  existing <- params$primary_signature_file_candidates[
    file.exists(params$primary_signature_file_candidates)
  ]

  if (length(existing) > 0L) {
    source_file <- existing[1]
    primary <- unique(clean_symbol(readLines(source_file, warn = FALSE)))
    primary <- primary[valid_gene_symbol(primary)]
    source_type <- "existing_primary_signature_file"
  } else {
    source_file <- NA_character_
    primary <- embedded
    source_type <- "embedded_established_24_gene_set"
  }

  list(
    genes = primary,
    embedded = embedded,
    source_file = source_file,
    source_type = source_type,
    exact_embedded_match = setequal(primary, embedded)
  )
}

primary_info <- load_primary_signature(params)
PRIMARY_SIGNATURE <- primary_info$genes
EMBEDDED_SIGNATURE <- primary_info$embedded

primary_check <- data.frame(
  metric = c(
    "primary_signature_source_type",
    "primary_signature_source_file",
    "primary_signature_size",
    "embedded_signature_size",
    "primary_matches_embedded_set_exactly"
  ),
  value = c(
    primary_info$source_type,
    ifelse(is.na(primary_info$source_file), "not_used", primary_info$source_file),
    length(PRIMARY_SIGNATURE),
    length(EMBEDDED_SIGNATURE),
    primary_info$exact_embedded_match
  ),
  stringsAsFactors = FALSE
)

write_csv_overwrite(
  primary_check,
  file.path(params$out_dir, "qc", "PRIMARY_SIGNATURE_FILE_CHECK.csv")
)

write_gene_list(
  PRIMARY_SIGNATURE,
  file.path(params$out_dir, "tables", "00_ESTABLISHED_PRIMARY_24_GENE_SIGNATURE.txt")
)

write_csv_overwrite(
  data.frame(
    primary_order = seq_along(PRIMARY_SIGNATURE),
    gene_symbol = PRIMARY_SIGNATURE,
    stringsAsFactors = FALSE
  ),
  file.path(params$out_dir, "tables", "00_ESTABLISHED_PRIMARY_24_GENE_SIGNATURE.csv")
)

if (length(PRIMARY_SIGNATURE) != 24L) {
  stop(
    "The established primary signature contains ",
    length(PRIMARY_SIGNATURE),
    " genes rather than 24."
  )
}

if (
  isTRUE(params$stop_if_primary_file_differs_from_embedded_set) &&
    !primary_info$exact_embedded_match
) {
  write_gene_list(
    setdiff(EMBEDDED_SIGNATURE, PRIMARY_SIGNATURE),
    file.path(params$out_dir, "qc", "GENES_IN_EMBEDDED_SET_NOT_PRIMARY_FILE.txt")
  )
  write_gene_list(
    setdiff(PRIMARY_SIGNATURE, EMBEDDED_SIGNATURE),
    file.path(params$out_dir, "qc", "GENES_IN_PRIMARY_FILE_NOT_EMBEDDED_SET.txt")
  )
  stop(
    "The existing primary signature file does not match the embedded established ",
    "24-gene set. See the QC files before proceeding."
  )
}

message("STEP 1/9: Loading GTEx and computing tissue-specificity metrics...")
require_file(params$gtex_file, "GTEx input")
GTEx <- data.table::fread(params$gtex_file, data.table = FALSE, check.names = FALSE)
require_columns(GTEx, "gene_symbol", "GTEx table")

GTEx$gene_uc <- clean_symbol(GTEx$gene_symbol)
GTEx$valid_gene_symbol <- valid_gene_symbol(GTEx$gene_uc)

tissue_cols <- setdiff(
  names(GTEx),
  c("gene_symbol", "gene_id", "Name", "Description", "gene_uc", "valid_gene_symbol")
)

GTEx[tissue_cols] <- lapply(GTEx[tissue_cols], to_num)

stomach_cols <- grep("^Stomach", names(GTEx), value = TRUE)
colon_cols <- grep("^Colon", names(GTEx), value = TRUE)
colon_cols <- colon_cols[!grepl("Muscularis", colon_cols, ignore.case = TRUE)]

if (length(stomach_cols) == 0L) stop("No GTEx stomach columns were detected.")
if (length(colon_cols) == 0L) stop("No GTEx colon columns were detected.")

write_text_overwrite(
  stomach_cols,
  file.path(params$out_dir, "qc", "GTEx_STOMACH_COLUMNS_USED.txt")
)
write_text_overwrite(
  colon_cols,
  file.path(params$out_dir, "qc", "GTEx_COLON_COLUMNS_USED.txt")
)

GTEx_expression_matrix <- as.matrix(GTEx[, tissue_cols, drop = FALSE])
storage.mode(GTEx_expression_matrix) <- "double"

GTEx$tau <- apply(GTEx_expression_matrix, 1, tau_row)

GTEx$GTEx_stomach_mean_TPM <- rowMeans(
  as.matrix(GTEx[, stomach_cols, drop = FALSE]),
  na.rm = TRUE
)
GTEx$GTEx_stomach_mean_TPM[is.nan(GTEx$GTEx_stomach_mean_TPM)] <- NA_real_

GTEx$GTEx_colon_max_TPM <- apply(
  as.matrix(GTEx[, colon_cols, drop = FALSE]),
  1,
  safe_max
)

GTEx$GTEx_stomach_vs_colon_log2FC <-
  log2(GTEx$GTEx_stomach_mean_TPM + params$pseudocount) -
  log2(GTEx$GTEx_colon_max_TPM + params$pseudocount)

GTEX_INPUT_GENES <- unique(GTEx$gene_uc[GTEx$valid_gene_symbol])

message("STEP 2/9: Loading HPA and computing validation metrics...")
require_file(params$hpa_file, "HPA input")
HPA <- read_table_flexible(params$hpa_file)
names(HPA) <- tolower(names(HPA))

hpa_gene_col <- if ("gene name" %in% names(HPA)) {
  "gene name"
} else if ("gene" %in% names(HPA)) {
  "gene"
} else {
  NA_character_
}

hpa_tissue_col <- if ("tissue" %in% names(HPA)) "tissue" else NA_character_
hpa_ntpm_col <- if ("ntpm" %in% names(HPA)) "ntpm" else NA_character_

if (any(is.na(c(hpa_gene_col, hpa_tissue_col, hpa_ntpm_col)))) {
  stop("Could not identify the HPA gene, tissue and nTPM columns.")
}

HPA$gene_uc <- clean_symbol(HPA[[hpa_gene_col]])
HPA$tissue_lc <- tolower(trimws(as.character(HPA[[hpa_tissue_col]])))
HPA$ntpm_numeric <- to_num(HPA[[hpa_ntpm_col]])

HPA <- HPA %>%
  dplyr::filter(valid_gene_symbol(gene_uc), nzchar(tissue_lc), is.finite(ntpm_numeric))

HPA_stomach <- HPA %>%
  dplyr::filter(tissue_lc == "stomach") %>%
  dplyr::group_by(gene_uc) %>%
  dplyr::summarise(
    HPA_stomach_nTPM = safe_max(ntpm_numeric),
    HPA_stomach_rows = dplyr::n(),
    .groups = "drop"
  )

HPA_colon <- HPA %>%
  dplyr::filter(grepl("^colon", tissue_lc)) %>%
  dplyr::group_by(gene_uc) %>%
  dplyr::summarise(
    HPA_colon_max_nTPM = safe_max(ntpm_numeric),
    HPA_colon_rows = dplyr::n(),
    .groups = "drop"
  )

HPA_highest_other <- HPA %>%
  dplyr::filter(tissue_lc != "stomach") %>%
  dplyr::arrange(gene_uc, dplyr::desc(ntpm_numeric), tissue_lc) %>%
  dplyr::group_by(gene_uc) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    gene_uc,
    HPA_highest_other_tissue = tissue_lc,
    HPA_highest_other_nTPM = ntpm_numeric
  )

EVIDENCE <- GTEx %>%
  dplyr::left_join(HPA_stomach, by = "gene_uc") %>%
  dplyr::left_join(HPA_colon, by = "gene_uc") %>%
  dplyr::left_join(HPA_highest_other, by = "gene_uc")

EVIDENCE$HPA_stomach_available <- is.finite(EVIDENCE$HPA_stomach_nTPM)
EVIDENCE$HPA_colon_available <- is.finite(EVIDENCE$HPA_colon_max_nTPM)
EVIDENCE$HPA_highest_other_available <- is.finite(EVIDENCE$HPA_highest_other_nTPM)

EVIDENCE$HPA_colon_for_original_log2FC <- ifelse(
  EVIDENCE$HPA_colon_available,
  EVIDENCE$HPA_colon_max_nTPM,
  0
)

EVIDENCE$HPA_stomach_vs_colon_log2FC_original <-
  log2(EVIDENCE$HPA_stomach_nTPM + params$pseudocount) -
  log2(EVIDENCE$HPA_colon_for_original_log2FC + params$pseudocount)

EVIDENCE$HPA_stomach_vs_colon_log2FC_complete_case <- ifelse(
  EVIDENCE$HPA_stomach_available & EVIDENCE$HPA_colon_available,
  log2(EVIDENCE$HPA_stomach_nTPM + params$pseudocount) -
    log2(EVIDENCE$HPA_colon_max_nTPM + params$pseudocount),
  NA_real_
)

EVIDENCE$HPA_stomach_vs_highest_other_log2FC <- ifelse(
  EVIDENCE$HPA_stomach_available & EVIDENCE$HPA_highest_other_available,
  log2(EVIDENCE$HPA_stomach_nTPM + params$pseudocount) -
    log2(EVIDENCE$HPA_highest_other_nTPM + params$pseudocount),
  NA_real_
)

message("STEP 3/9: Loading OSLC ANTRUM-versus-H_COL results...")
DE_combined <- read_table_flexible(params$combined_de_file)

antrum_log2fc_col <- paste0(params$antrum_prefix, ".logFC")
antrum_fdr_col <- paste0(params$antrum_prefix, ".FDR")
antrum_symbol_col <- paste0(params$antrum_prefix, ".SYMBOL")

if (!(antrum_symbol_col %in% names(DE_combined))) {
  alternative_symbol_col <- paste0(params$antrum_prefix, ".symbol")
  if (alternative_symbol_col %in% names(DE_combined)) {
    antrum_symbol_col <- alternative_symbol_col
  }
}

require_columns(
  DE_combined,
  c(antrum_symbol_col, antrum_log2fc_col, antrum_fdr_col),
  "combined OSLC DE table"
)

ANTRUM_DE <- data.frame(
  gene_symbol = standardize_gene_vector(DE_combined[[antrum_symbol_col]]),
  ANTRUM_vs_H_COL_log2FC = to_num(DE_combined[[antrum_log2fc_col]]),
  ANTRUM_vs_H_COL_FDR = to_num(DE_combined[[antrum_fdr_col]]),
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(valid_gene_symbol(gene_symbol)) %>%
  dplyr::arrange(
    dplyr::desc(abs(ANTRUM_vs_H_COL_log2FC)),
    ifelse(is.finite(ANTRUM_vs_H_COL_FDR), ANTRUM_vs_H_COL_FDR, Inf)
  ) %>%
  dplyr::distinct(gene_symbol, .keep_all = TRUE)

EVIDENCE <- EVIDENCE %>%
  dplyr::left_join(
    ANTRUM_DE,
    by = c("gene_uc" = "gene_symbol")
  )

message("STEP 4/9: Loading OSLC SSL_HP-versus-H_COL results...")
SSL_DE <- standardize_de_table(params$ssl_hp_hcol_de_file)

if (any(is.finite(SSL_DE$adjusted_p))) {
  SSL_DE$SSL_significance_value <- SSL_DE$adjusted_p
  SSL_SIGNIFICANCE_SOURCE <- "adjusted_p_or_FDR"
} else if (
  isTRUE(params$use_raw_p_if_adjusted_p_missing) &&
    any(is.finite(SSL_DE$raw_p))
) {
  SSL_DE$SSL_significance_value <- SSL_DE$raw_p
  SSL_SIGNIFICANCE_SOURCE <- "raw_p_used_because_adjusted_p_missing"
} else {
  stop(
    "The SSL_HP-versus-H_COL DE table contains no usable adjusted P/FDR ",
    "column. Raw P values were not substituted."
  )
}

SSL_FOR_JOIN <- SSL_DE %>%
  dplyr::transmute(
    gene_symbol = Gene,
    SSL_HP_vs_H_COL_log2FC = log2FC,
    SSL_HP_vs_H_COL_adjusted_p = adjusted_p,
    SSL_HP_vs_H_COL_raw_p = raw_p,
    SSL_HP_vs_H_COL_significance_value = SSL_significance_value
  )

EVIDENCE <- EVIDENCE %>%
  dplyr::left_join(
    SSL_FOR_JOIN,
    by = c("gene_uc" = "gene_symbol")
  )

write_csv_overwrite(
  EVIDENCE,
  file.path(params$out_dir, "tables", "ALL_GENES_ALL_FILTER_METRICS.csv")
)

baseline_setting_values <- list(
  gtex_tau_min = params$TAU_MIN,
  gtex_stomach_min = params$GTEX_STOMACH_TPM_MIN,
  gtex_colon_max = params$GTEX_COLON_TPM_MAX,
  gtex_log2fc_min = params$GTEX_LOG2FC_MIN,

  hpa_rule = "colon_focused",
  hpa_stomach_min = params$HPA_STOMACH_NTPM_MIN,
  hpa_colon_max = params$HPA_COLON_NTPM_MAX,
  hpa_log2fc_min = params$HPA_LOG2FC_MIN,
  hpa_missing_rule = params$HPA_MISSING_RULE,

  antrum_log2fc_min = params$ANTRUM_LOG2FC_MIN,
  antrum_fdr_max = params$ANTRUM_FDR_MAX,

  ssl_hp_log2fc_min = params$SSL_HP_LOG2FC_MIN,
  ssl_hp_fdr_max = params$SSL_HP_FDR_MAX
)

make_setting <- function(
    setting_id,
    analysis_family,
    setting_label,
    overrides = list(),
    varied_parameter = NA_character_,
    varied_value_label = NA_character_
) {
  values <- utils::modifyList(baseline_setting_values, overrides)

  data.frame(
    setting_id = as.character(setting_id),
    analysis_family = as.character(analysis_family),
    setting_label = as.character(setting_label),
    varied_parameter = as.character(varied_parameter),
    varied_value_label = as.character(varied_value_label),

    gtex_tau_min = as.numeric(values$gtex_tau_min),
    gtex_stomach_min = as.numeric(values$gtex_stomach_min),
    gtex_colon_max = as.numeric(values$gtex_colon_max),
    gtex_log2fc_min = as.numeric(values$gtex_log2fc_min),

    hpa_rule = as.character(values$hpa_rule),
    hpa_stomach_min = as.numeric(values$hpa_stomach_min),
    hpa_colon_max = as.numeric(values$hpa_colon_max),
    hpa_log2fc_min = as.numeric(values$hpa_log2fc_min),
    hpa_missing_rule = as.character(values$hpa_missing_rule),

    antrum_log2fc_min = as.numeric(values$antrum_log2fc_min),
    antrum_fdr_max = as.numeric(values$antrum_fdr_max),

    ssl_hp_log2fc_min = as.numeric(values$ssl_hp_log2fc_min),
    ssl_hp_fdr_max = as.numeric(values$ssl_hp_fdr_max),
    stringsAsFactors = FALSE
  )
}

compute_stage_passes <- function(evidence, setting) {
  pass_gtex <-
    evidence$valid_gene_symbol &
    is.finite(evidence$tau) &
    evidence$tau >= setting$gtex_tau_min &
    is.finite(evidence$GTEx_stomach_mean_TPM) &
    evidence$GTEx_stomach_mean_TPM >= setting$gtex_stomach_min &
    is.finite(evidence$GTEx_colon_max_TPM) &
    evidence$GTEx_colon_max_TPM <= setting$gtex_colon_max &
    is.finite(evidence$GTEx_stomach_vs_colon_log2FC) &
    evidence$GTEx_stomach_vs_colon_log2FC >= setting$gtex_log2fc_min

  if (identical(setting$hpa_rule, "none")) {
    pass_hpa <- evidence$valid_gene_symbol
  } else if (identical(setting$hpa_rule, "colon_focused")) {
    pass_hpa_stomach <-
      evidence$HPA_stomach_available &
      evidence$HPA_stomach_nTPM >= setting$hpa_stomach_min

    if (identical(setting$hpa_missing_rule, "original_missing_as_zero")) {
      pass_hpa_colon <-
        !evidence$HPA_colon_available |
        evidence$HPA_colon_max_nTPM <= setting$hpa_colon_max

      hpa_log2fc_use <- evidence$HPA_stomach_vs_colon_log2FC_original
    } else if (identical(setting$hpa_missing_rule, "complete_case")) {
      pass_hpa_colon <-
        evidence$HPA_colon_available &
        evidence$HPA_colon_max_nTPM <= setting$hpa_colon_max

      hpa_log2fc_use <- evidence$HPA_stomach_vs_colon_log2FC_complete_case
    } else {
      stop("Unknown HPA missing-value rule: ", setting$hpa_missing_rule)
    }

    pass_hpa_log2fc <-
      is.finite(hpa_log2fc_use) &
      hpa_log2fc_use >= setting$hpa_log2fc_min

    pass_hpa <- pass_hpa_stomach & pass_hpa_colon & pass_hpa_log2fc
  } else if (identical(setting$hpa_rule, "all_tissue_enriched")) {

    pass_hpa <-
      evidence$HPA_stomach_available &
      evidence$HPA_stomach_nTPM >= setting$hpa_stomach_min &
      evidence$HPA_highest_other_available &
      is.finite(evidence$HPA_stomach_vs_highest_other_log2FC) &
      evidence$HPA_stomach_vs_highest_other_log2FC >= setting$hpa_log2fc_min
  } else {
    stop("Unknown HPA rule: ", setting$hpa_rule)
  }

  pass_antrum <-
    is.finite(evidence$ANTRUM_vs_H_COL_log2FC) &
    is.finite(evidence$ANTRUM_vs_H_COL_FDR) &
    evidence$ANTRUM_vs_H_COL_log2FC >= setting$antrum_log2fc_min &
    evidence$ANTRUM_vs_H_COL_FDR <= setting$antrum_fdr_max

  pass_ssl <-
    is.finite(evidence$SSL_HP_vs_H_COL_log2FC) &
    is.finite(evidence$SSL_HP_vs_H_COL_significance_value) &
    evidence$SSL_HP_vs_H_COL_log2FC > setting$ssl_hp_log2fc_min &
    evidence$SSL_HP_vs_H_COL_significance_value <= setting$ssl_hp_fdr_max

  pass_after_gtex <- pass_gtex
  pass_after_hpa <- pass_after_gtex & pass_hpa
  pass_after_antrum <- pass_after_hpa & pass_antrum
  pass_final <- pass_after_antrum & pass_ssl

  list(
    pass_gtex = pass_gtex,
    pass_hpa = pass_hpa,
    pass_antrum = pass_antrum,
    pass_ssl = pass_ssl,
    pass_after_gtex = pass_after_gtex,
    pass_after_hpa = pass_after_hpa,
    pass_after_antrum = pass_after_antrum,
    pass_final = pass_final
  )
}

unique_genes_passing <- function(evidence, pass_vector) {
  genes <- unique(evidence$gene_uc[pass_vector & evidence$valid_gene_symbol])
  genes[valid_gene_symbol(genes)]
}

evaluate_setting <- function(evidence, setting, primary_signature) {
  passes <- compute_stage_passes(evidence, setting)

  genes_after_gtex <- unique_genes_passing(evidence, passes$pass_after_gtex)
  genes_after_hpa <- unique_genes_passing(evidence, passes$pass_after_hpa)
  genes_after_antrum <- unique_genes_passing(evidence, passes$pass_after_antrum)
  final_genes <- unique_genes_passing(evidence, passes$pass_final)
  final_genes <- order_gene_set(final_genes, primary_signature)

  primary_retained <- primary_signature[primary_signature %in% final_genes]
  primary_removed <- primary_signature[!(primary_signature %in% final_genes)]
  additional_genes <- setdiff(final_genes, primary_signature)

  union_size <- length(union(final_genes, primary_signature))
  jaccard <- if (union_size == 0L) {
    NA_real_
  } else {
    length(intersect(final_genes, primary_signature)) / union_size
  }

  stage_summary <- data.frame(
    n_input_GTEx_genes = length(GTEX_INPUT_GENES),
    n_after_GTEx = length(genes_after_gtex),
    n_after_HPA = length(genes_after_hpa),
    n_after_ANTRUM = length(genes_after_antrum),
    n_final = length(final_genes),
    n_primary_24_retained = length(primary_retained),
    n_primary_24_removed = length(primary_removed),
    n_additional_genes_not_in_primary = length(additional_genes),
    exact_match_to_primary_24 = setequal(final_genes, primary_signature),
    jaccard_with_primary_24 = jaccard,
    primary_genes_removed = collapse_gene_string(primary_removed),
    additional_genes_not_in_primary = collapse_gene_string(additional_genes),
    final_gene_list = collapse_gene_string(final_genes),
    stringsAsFactors = FALSE
  )

  gtex_set <- genes_after_gtex
  hpa_set <- genes_after_hpa
  antrum_set <- genes_after_antrum
  final_set <- final_genes

  primary_membership <- data.frame(
    gene_symbol = primary_signature,
    pass_after_GTEx = primary_signature %in% gtex_set,
    pass_after_HPA = primary_signature %in% hpa_set,
    pass_after_ANTRUM = primary_signature %in% antrum_set,
    pass_final = primary_signature %in% final_set,
    stringsAsFactors = FALSE
  )

  primary_membership$first_stage_failed <- dplyr::case_when(
    !primary_membership$pass_after_GTEx ~ "GTEx",
    primary_membership$pass_after_GTEx & !primary_membership$pass_after_HPA ~ "HPA",
    primary_membership$pass_after_HPA & !primary_membership$pass_after_ANTRUM ~ "ANTRUM",
    primary_membership$pass_after_ANTRUM & !primary_membership$pass_final ~ "SSL_HP",
    TRUE ~ "RETAINED"
  )

  final_membership <- data.frame(
    gene_symbol = final_genes,
    in_established_primary_24 = final_genes %in% primary_signature,
    primary_order = match(final_genes, primary_signature),
    stringsAsFactors = FALSE
  )

  list(
    summary = stage_summary,
    primary_membership = primary_membership,
    final_membership = final_membership,
    final_genes = final_genes,
    passes = passes
  )
}

run_settings_table <- function(settings, evidence, primary_signature) {
  if (nrow(settings) == 0L) {
    return(list(
      summary = data.frame(),
      primary_membership = data.frame(),
      final_membership = data.frame(),
      final_gene_lists = list()
    ))
  }

  summary_list <- vector("list", nrow(settings))
  primary_list <- vector("list", nrow(settings))
  final_list <- vector("list", nrow(settings))
  final_gene_lists <- vector("list", nrow(settings))
  names(final_gene_lists) <- settings$setting_id

  for (i in seq_len(nrow(settings))) {
    setting <- settings[i, , drop = FALSE]
    result <- evaluate_setting(evidence, setting, primary_signature)

    summary_list[[i]] <- dplyr::bind_cols(
      setting,
      result$summary
    )

    primary_meta <- data.frame(
      setting_id = rep(as.character(setting$setting_id), nrow(result$primary_membership)),
      analysis_family = rep(as.character(setting$analysis_family), nrow(result$primary_membership)),
      setting_label = rep(as.character(setting$setting_label), nrow(result$primary_membership)),
      stringsAsFactors = FALSE
    )

    primary_list[[i]] <- dplyr::bind_cols(
      primary_meta,
      result$primary_membership
    )

    if (nrow(result$final_membership) > 0L) {
      final_meta <- data.frame(
        setting_id = rep(as.character(setting$setting_id), nrow(result$final_membership)),
        analysis_family = rep(as.character(setting$analysis_family), nrow(result$final_membership)),
        setting_label = rep(as.character(setting$setting_label), nrow(result$final_membership)),
        stringsAsFactors = FALSE
      )

      final_list[[i]] <- dplyr::bind_cols(
        final_meta,
        result$final_membership
      )
    } else {
      final_list[[i]] <- data.frame(
        setting_id = character(0),
        analysis_family = character(0),
        setting_label = character(0),
        gene_symbol = character(0),
        in_established_primary_24 = logical(0),
        primary_order = integer(0),
        stringsAsFactors = FALSE
      )
    }

    final_gene_lists[[setting$setting_id]] <- result$final_genes

    if (isTRUE(params$write_individual_txt_for_every_setting)) {
      write_gene_list(
        result$final_genes,
        file.path(
          params$out_dir,
          "all_setting_gene_lists",
          paste0(safe_file_tag(setting$setting_id), ".txt")
        )
      )
    }
  }

  list(
    summary = dplyr::bind_rows(summary_list),
    primary_membership = dplyr::bind_rows(primary_list),
    final_membership = dplyr::bind_rows(final_list),
    final_gene_lists = final_gene_lists
  )
}

save_run_outputs <- function(run_object, prefix) {
  write_csv_overwrite(
    run_object$summary,
    file.path(params$out_dir, "tables", paste0(prefix, "_SUMMARY.csv"))
  )
  write_csv_overwrite(
    run_object$primary_membership,
    file.path(params$out_dir, "tables", paste0(prefix, "_PRIMARY_24_MEMBERSHIP.csv"))
  )
  write_csv_overwrite(
    run_object$final_membership,
    file.path(params$out_dir, "tables", paste0(prefix, "_FINAL_GENE_MEMBERSHIP.csv"))
  )
}

message("STEP 5/9: Building sensitivity settings...")

BASELINE_SETTINGS <- make_setting(
  setting_id = "BASELINE_ORIGINAL_PIPELINE",
  analysis_family = "BASELINE",
  setting_label = "Original thresholds"
)

ofat_specs <- list(
  list(
    parameter = "gtex_tau_min",
    values = sensitivity_values$gtex_tau_min,
    family_label = "GTEx tau minimum"
  ),
  list(
    parameter = "gtex_stomach_min",
    values = sensitivity_values$gtex_stomach_min,
    family_label = "GTEx stomach TPM minimum"
  ),
  list(
    parameter = "gtex_colon_max",
    values = sensitivity_values$gtex_colon_max,
    family_label = "GTEx colon TPM maximum"
  ),
  list(
    parameter = "gtex_log2fc_min",
    values = sensitivity_values$gtex_log2fc_min,
    family_label = "GTEx stomach-versus-colon log2FC minimum"
  ),
  list(
    parameter = "hpa_stomach_min",
    values = sensitivity_values$hpa_stomach_min,
    family_label = "HPA stomach nTPM minimum"
  ),
  list(
    parameter = "hpa_colon_max",
    values = sensitivity_values$hpa_colon_max,
    family_label = "HPA colon nTPM maximum"
  ),
  list(
    parameter = "hpa_log2fc_min",
    values = sensitivity_values$hpa_log2fc_min,
    family_label = "HPA stomach-versus-colon log2FC minimum"
  ),
  list(
    parameter = "hpa_missing_rule",
    values = sensitivity_values$hpa_missing_rule,
    family_label = "HPA missing-colon rule"
  ),
  list(
    parameter = "antrum_log2fc_min",
    values = sensitivity_values$antrum_log2fc_min,
    family_label = "ANTRUM-versus-H_COL log2FC minimum"
  ),
  list(
    parameter = "antrum_fdr_max",
    values = sensitivity_values$antrum_fdr_max,
    family_label = "ANTRUM-versus-H_COL FDR maximum"
  ),
  list(
    parameter = "ssl_hp_log2fc_min",
    values = sensitivity_values$ssl_hp_log2fc_min,
    family_label = "SSL_HP-versus-H_COL log2FC minimum"
  ),
  list(
    parameter = "ssl_hp_fdr_max",
    values = sensitivity_values$ssl_hp_fdr_max,
    family_label = "SSL_HP-versus-H_COL FDR maximum"
  )
)

ofat_rows <- list()
ofat_counter <- 1L

for (spec in ofat_specs) {
  for (value in spec$values) {
    value_label <- if (is.character(value)) {
      as.character(value)
    } else if (is.numeric(value) && is.infinite(value)) {
      "no ceiling"
    } else {
      as.character(signif(as.numeric(value), 8))
    }

    overrides <- list()
    overrides[[spec$parameter]] <- value

    ofat_rows[[ofat_counter]] <- make_setting(
      setting_id = paste0(
        "OFAT_",
        safe_file_tag(spec$parameter),
        "_",
        safe_file_tag(value_label)
      ),
      analysis_family = "ONE_FACTOR_AT_A_TIME",
      setting_label = paste0(spec$family_label, " = ", value_label),
      overrides = overrides,
      varied_parameter = spec$parameter,
      varied_value_label = value_label
    )
    ofat_counter <- ofat_counter + 1L
  }
}

ofat_rows[[ofat_counter]] <- make_setting(
  setting_id = "OFAT_HPA_RULE_NO_COLON_CEILING",
  analysis_family = "ONE_FACTOR_AT_A_TIME",
  setting_label = "HPA original rule without a separate colon ceiling",
  overrides = list(hpa_colon_max = Inf),
  varied_parameter = "hpa_rule_variant",
  varied_value_label = "no_colon_ceiling"
)
ofat_counter <- ofat_counter + 1L

ofat_rows[[ofat_counter]] <- make_setting(
  setting_id = "OFAT_HPA_RULE_FOURFOLD_VS_COLON",
  analysis_family = "ONE_FACTOR_AT_A_TIME",
  setting_label = "HPA stomach >=1 nTPM and >=4-fold versus colon; no colon ceiling",
  overrides = list(
    hpa_rule = "colon_focused",
    hpa_stomach_min = 1,
    hpa_colon_max = Inf,
    hpa_log2fc_min = 2,
    hpa_missing_rule = "complete_case"
  ),
  varied_parameter = "hpa_rule_variant",
  varied_value_label = "fourfold_vs_colon"
)
ofat_counter <- ofat_counter + 1L

ofat_rows[[ofat_counter]] <- make_setting(
  setting_id = "OFAT_HPA_RULE_FORMAL_ALL_TISSUE_FOURFOLD",
  analysis_family = "ONE_FACTOR_AT_A_TIME",
  setting_label = "HPA stomach >=1 nTPM and >=4-fold above every other tissue",
  overrides = list(
    hpa_rule = "all_tissue_enriched",
    hpa_stomach_min = 1,
    hpa_log2fc_min = 2,
    hpa_missing_rule = "complete_case"
  ),
  varied_parameter = "hpa_rule_variant",
  varied_value_label = "formal_all_tissue_fourfold"
)
ofat_counter <- ofat_counter + 1L

ofat_rows[[ofat_counter]] <- make_setting(
  setting_id = "OFAT_NO_HPA_FILTER",
  analysis_family = "ONE_FACTOR_AT_A_TIME",
  setting_label = "No HPA filter",
  overrides = list(hpa_rule = "none"),
  varied_parameter = "hpa_rule_variant",
  varied_value_label = "no_HPA_filter"
)

OFAT_SETTINGS <- dplyr::bind_rows(ofat_rows)

GTEX_GRID_RAW <- expand.grid(
  gtex_tau_min = sensitivity_values$gtex_tau_min,
  gtex_stomach_min = sensitivity_values$gtex_stomach_min,
  gtex_colon_max = sensitivity_values$gtex_colon_max,
  gtex_log2fc_min = sensitivity_values$gtex_log2fc_min,
  stringsAsFactors = FALSE
)

GTEX_GRID_SETTINGS <- dplyr::bind_rows(lapply(
  seq_len(nrow(GTEX_GRID_RAW)),
  function(i) {
    r <- GTEX_GRID_RAW[i, , drop = FALSE]
    make_setting(
      setting_id = paste0("GTEX_GRID_", sprintf("%03d", i)),
      analysis_family = "GTEX_WITHIN_STAGE_GRID",
      setting_label = paste0(
        "tau>=", signif(r$gtex_tau_min, 4),
        "; stomach>=", signif(r$gtex_stomach_min, 4),
        "; colon<=", signif(r$gtex_colon_max, 4),
        "; log2FC>=", signif(r$gtex_log2fc_min, 4)
      ),
      overrides = list(
        gtex_tau_min = r$gtex_tau_min,
        gtex_stomach_min = r$gtex_stomach_min,
        gtex_colon_max = r$gtex_colon_max,
        gtex_log2fc_min = r$gtex_log2fc_min
      )
    )
  }
))

HPA_GRID_RAW <- expand.grid(
  hpa_stomach_min = sensitivity_values$hpa_stomach_min,
  hpa_colon_max = sensitivity_values$hpa_colon_max,
  hpa_log2fc_min = sensitivity_values$hpa_log2fc_min,
  hpa_missing_rule = sensitivity_values$hpa_missing_rule,
  stringsAsFactors = FALSE
)

HPA_GRID_SETTINGS <- dplyr::bind_rows(lapply(
  seq_len(nrow(HPA_GRID_RAW)),
  function(i) {
    r <- HPA_GRID_RAW[i, , drop = FALSE]
    colon_label <- if (is.infinite(r$hpa_colon_max)) {
      "no ceiling"
    } else {
      as.character(signif(r$hpa_colon_max, 4))
    }

    make_setting(
      setting_id = paste0("HPA_GRID_", sprintf("%03d", i)),
      analysis_family = "HPA_WITHIN_STAGE_GRID",
      setting_label = paste0(
        "stomach>=", signif(r$hpa_stomach_min, 4),
        "; colon<=", colon_label,
        "; log2FC>=", signif(r$hpa_log2fc_min, 4),
        "; ", r$hpa_missing_rule
      ),
      overrides = list(
        hpa_rule = "colon_focused",
        hpa_stomach_min = r$hpa_stomach_min,
        hpa_colon_max = r$hpa_colon_max,
        hpa_log2fc_min = r$hpa_log2fc_min,
        hpa_missing_rule = r$hpa_missing_rule
      )
    )
  }
))

HPA_RULE_SETTINGS <- dplyr::bind_rows(
  make_setting(
    "HPA_RULE_ORIGINAL",
    "HPA_RULE_COMPARISON",
    "Original colon-focused HPA rule"
  ),
  make_setting(
    "HPA_RULE_ORIGINAL_COMPLETE_CASE",
    "HPA_RULE_COMPARISON",
    "Original HPA thresholds with complete-case colon measurement",
    overrides = list(hpa_missing_rule = "complete_case")
  ),
  make_setting(
    "HPA_RULE_NO_COLON_CEILING",
    "HPA_RULE_COMPARISON",
    "Original HPA stomach and eightfold criteria; no colon ceiling",
    overrides = list(hpa_colon_max = Inf)
  ),
  make_setting(
    "HPA_RULE_FOURFOLD_VS_COLON_MIN1",
    "HPA_RULE_COMPARISON",
    "Stomach >=1 nTPM and >=4-fold versus colon; no colon ceiling",
    overrides = list(
      hpa_rule = "colon_focused",
      hpa_stomach_min = 1,
      hpa_colon_max = Inf,
      hpa_log2fc_min = 2,
      hpa_missing_rule = "complete_case"
    )
  ),
  make_setting(
    "HPA_RULE_FOURFOLD_VS_COLON_MIN2",
    "HPA_RULE_COMPARISON",
    "Stomach >=2 nTPM and >=4-fold versus colon; no colon ceiling",
    overrides = list(
      hpa_rule = "colon_focused",
      hpa_stomach_min = 2,
      hpa_colon_max = Inf,
      hpa_log2fc_min = 2,
      hpa_missing_rule = "complete_case"
    )
  ),
  make_setting(
    "HPA_RULE_FORMAL_ALL_TISSUE_FOURFOLD_MIN1",
    "HPA_RULE_COMPARISON",
    "Stomach >=1 nTPM and >=4-fold above every other HPA tissue",
    overrides = list(
      hpa_rule = "all_tissue_enriched",
      hpa_stomach_min = 1,
      hpa_log2fc_min = 2,
      hpa_missing_rule = "complete_case"
    )
  ),
  make_setting(
    "HPA_RULE_FORMAL_ALL_TISSUE_FOURFOLD_MIN2",
    "HPA_RULE_COMPARISON",
    "Stomach >=2 nTPM and >=4-fold above every other HPA tissue",
    overrides = list(
      hpa_rule = "all_tissue_enriched",
      hpa_stomach_min = 2,
      hpa_log2fc_min = 2,
      hpa_missing_rule = "complete_case"
    )
  ),
  make_setting(
    "HPA_RULE_NONE",
    "HPA_RULE_COMPARISON",
    "No HPA filter",
    overrides = list(hpa_rule = "none")
  )
)

ANTRUM_GRID_RAW <- expand.grid(
  antrum_log2fc_min = sensitivity_values$antrum_log2fc_min,
  antrum_fdr_max = sensitivity_values$antrum_fdr_max,
  stringsAsFactors = FALSE
)

ANTRUM_GRID_SETTINGS <- dplyr::bind_rows(lapply(
  seq_len(nrow(ANTRUM_GRID_RAW)),
  function(i) {
    r <- ANTRUM_GRID_RAW[i, , drop = FALSE]
    make_setting(
      setting_id = paste0("ANTRUM_GRID_", sprintf("%03d", i)),
      analysis_family = "ANTRUM_WITHIN_STAGE_GRID",
      setting_label = paste0(
        "ANTRUM log2FC>=", signif(r$antrum_log2fc_min, 4),
        "; FDR<=", signif(r$antrum_fdr_max, 4)
      ),
      overrides = list(
        antrum_log2fc_min = r$antrum_log2fc_min,
        antrum_fdr_max = r$antrum_fdr_max
      )
    )
  }
))

SSL_GRID_RAW <- expand.grid(
  ssl_hp_log2fc_min = sensitivity_values$ssl_hp_log2fc_min,
  ssl_hp_fdr_max = sensitivity_values$ssl_hp_fdr_max,
  stringsAsFactors = FALSE
)

SSL_GRID_SETTINGS <- dplyr::bind_rows(lapply(
  seq_len(nrow(SSL_GRID_RAW)),
  function(i) {
    r <- SSL_GRID_RAW[i, , drop = FALSE]
    make_setting(
      setting_id = paste0("SSL_GRID_", sprintf("%03d", i)),
      analysis_family = "SSL_HP_WITHIN_STAGE_GRID",
      setting_label = paste0(
        "SSL_HP log2FC>", signif(r$ssl_hp_log2fc_min, 4),
        "; FDR<=", signif(r$ssl_hp_fdr_max, 4)
      ),
      overrides = list(
        ssl_hp_log2fc_min = r$ssl_hp_log2fc_min,
        ssl_hp_fdr_max = r$ssl_hp_fdr_max
      )
    )
  }
))

GTEX_SCENARIOS <- data.frame(
  gtex_scenario = c(
    "GTEX_ORIGINAL",
    "GTEX_MIN1_FIVEFOLD",
    "GTEX_TAU_RELAXED",
    "GTEX_TAU_STRICT",
    "GTEX_COLON_RELAXED",
    "GTEX_STRICT_COMBINED"
  ),
  gtex_tau_min = c(0.85, 0.85, 0.80, 0.90, 0.85, 0.90),
  gtex_stomach_min = c(2, 1, 2, 2, 2, 2),
  gtex_colon_max = c(1, 1, 1, 1, 2, 0.5),
  gtex_log2fc_min = c(2.5, log2(5), 2.5, 2.5, 2.5, 3),
  stringsAsFactors = FALSE
)

HPA_SCENARIOS <- data.frame(
  hpa_scenario = c(
    "HPA_ORIGINAL",
    "HPA_ORIGINAL_COMPLETE_CASE",
    "HPA_NO_COLON_CEILING",
    "HPA_FOURFOLD_VS_COLON",
    "HPA_FORMAL_ALL_TISSUE_FOURFOLD",
    "HPA_NONE"
  ),
  hpa_rule = c(
    "colon_focused",
    "colon_focused",
    "colon_focused",
    "colon_focused",
    "all_tissue_enriched",
    "none"
  ),
  hpa_stomach_min = c(2, 2, 2, 1, 1, 2),
  hpa_colon_max = c(10, 10, Inf, Inf, Inf, 10),
  hpa_log2fc_min = c(3, 3, 3, 2, 2, 3),
  hpa_missing_rule = c(
    "original_missing_as_zero",
    "complete_case",
    "complete_case",
    "complete_case",
    "complete_case",
    "complete_case"
  ),
  stringsAsFactors = FALSE
)

ANTRUM_SCENARIOS <- data.frame(
  antrum_scenario = c(
    "ANTRUM_DIRECTION_ONLY",
    "ANTRUM_ORIGINAL_1P5",
    "ANTRUM_LOG2FC_2",
    "ANTRUM_LOG2FC_FIVEFOLD",
    "ANTRUM_LOG2FC_2P5",
    "ANTRUM_LOG2FC_3",
    "ANTRUM_1P5_FDR_0P01"
  ),
  antrum_log2fc_min = c(0, 1.5, 2, log2(5), 2.5, 3, 1.5),
  antrum_fdr_max = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.01),
  stringsAsFactors = FALSE
)

SSL_SCENARIOS <- data.frame(
  ssl_scenario = c(
    "SSL_ORIGINAL_DIRECTION_ONLY",
    "SSL_LOG2FC_0P5",
    "SSL_LOG2FC_1",
    "SSL_LOG2FC_1P5",
    "SSL_DIRECTION_FDR_0P01"
  ),
  ssl_hp_log2fc_min = c(0, 0.5, 1, 1.5, 0),
  ssl_hp_fdr_max = c(0.05, 0.05, 0.05, 0.05, 0.01),
  stringsAsFactors = FALSE
)

COMBINED_INDEX <- expand.grid(
  gtex_i = seq_len(nrow(GTEX_SCENARIOS)),
  hpa_i = seq_len(nrow(HPA_SCENARIOS)),
  antrum_i = seq_len(nrow(ANTRUM_SCENARIOS)),
  ssl_i = seq_len(nrow(SSL_SCENARIOS)),
  stringsAsFactors = FALSE
)

COMBINED_SETTINGS <- dplyr::bind_rows(lapply(
  seq_len(nrow(COMBINED_INDEX)),
  function(i) {
    idx <- COMBINED_INDEX[i, , drop = FALSE]
    g <- GTEX_SCENARIOS[idx$gtex_i, , drop = FALSE]
    h <- HPA_SCENARIOS[idx$hpa_i, , drop = FALSE]
    a <- ANTRUM_SCENARIOS[idx$antrum_i, , drop = FALSE]
    s <- SSL_SCENARIOS[idx$ssl_i, , drop = FALSE]

    make_setting(
      setting_id = paste0("COMBINED_", sprintf("%04d", i)),
      analysis_family = "COMBINED_REPRESENTATIVE_GRID",
      setting_label = paste(
        g$gtex_scenario,
        h$hpa_scenario,
        a$antrum_scenario,
        s$ssl_scenario,
        sep = " | "
      ),
      overrides = list(
        gtex_tau_min = g$gtex_tau_min,
        gtex_stomach_min = g$gtex_stomach_min,
        gtex_colon_max = g$gtex_colon_max,
        gtex_log2fc_min = g$gtex_log2fc_min,

        hpa_rule = h$hpa_rule,
        hpa_stomach_min = h$hpa_stomach_min,
        hpa_colon_max = h$hpa_colon_max,
        hpa_log2fc_min = h$hpa_log2fc_min,
        hpa_missing_rule = h$hpa_missing_rule,

        antrum_log2fc_min = a$antrum_log2fc_min,
        antrum_fdr_max = a$antrum_fdr_max,

        ssl_hp_log2fc_min = s$ssl_hp_log2fc_min,
        ssl_hp_fdr_max = s$ssl_hp_fdr_max
      )
    ) %>%
      dplyr::mutate(
        gtex_scenario = g$gtex_scenario,
        hpa_scenario = h$hpa_scenario,
        antrum_scenario = a$antrum_scenario,
        ssl_scenario = s$ssl_scenario
      )
  }
))

write_csv_overwrite(BASELINE_SETTINGS, file.path(params$out_dir, "qc", "SETTINGS_BASELINE.csv"))
write_csv_overwrite(OFAT_SETTINGS, file.path(params$out_dir, "qc", "SETTINGS_ONE_FACTOR_AT_A_TIME.csv"))
write_csv_overwrite(GTEX_GRID_SETTINGS, file.path(params$out_dir, "qc", "SETTINGS_GTEX_GRID.csv"))
write_csv_overwrite(HPA_GRID_SETTINGS, file.path(params$out_dir, "qc", "SETTINGS_HPA_GRID.csv"))
write_csv_overwrite(HPA_RULE_SETTINGS, file.path(params$out_dir, "qc", "SETTINGS_HPA_RULE_COMPARISON.csv"))
write_csv_overwrite(ANTRUM_GRID_SETTINGS, file.path(params$out_dir, "qc", "SETTINGS_ANTRUM_GRID.csv"))
write_csv_overwrite(SSL_GRID_SETTINGS, file.path(params$out_dir, "qc", "SETTINGS_SSL_GRID.csv"))
write_csv_overwrite(COMBINED_SETTINGS, file.path(params$out_dir, "qc", "SETTINGS_COMBINED_REPRESENTATIVE_GRID.csv"))

message("STEP 6/9: Reconstructing the original 24-gene signature...")
BASELINE_RUN <- run_settings_table(BASELINE_SETTINGS, EVIDENCE, PRIMARY_SIGNATURE)
save_run_outputs(BASELINE_RUN, "01_BASELINE")

BASELINE_FINAL_GENES <- BASELINE_RUN$final_gene_lists[["BASELINE_ORIGINAL_PIPELINE"]]
BASELINE_EXACT_MATCH <- setequal(BASELINE_FINAL_GENES, PRIMARY_SIGNATURE)

BASELINE_DETAIL <- evaluate_setting(
  EVIDENCE,
  BASELINE_SETTINGS[1, , drop = FALSE],
  PRIMARY_SIGNATURE
)

write_gene_list(
  unique_genes_passing(EVIDENCE, BASELINE_DETAIL$passes$pass_after_gtex),
  file.path(params$out_dir, "selected_gene_lists", "BASELINE_01_AFTER_GTEX.txt")
)
write_gene_list(
  unique_genes_passing(EVIDENCE, BASELINE_DETAIL$passes$pass_after_hpa),
  file.path(params$out_dir, "selected_gene_lists", "BASELINE_02_AFTER_HPA.txt")
)
write_gene_list(
  unique_genes_passing(EVIDENCE, BASELINE_DETAIL$passes$pass_after_antrum),
  file.path(params$out_dir, "selected_gene_lists", "BASELINE_03_AFTER_ANTRUM_VS_H_COL.txt")
)
write_gene_list(
  BASELINE_FINAL_GENES,
  file.path(params$out_dir, "selected_gene_lists", "BASELINE_04_AFTER_SSL_HP_VS_H_COL_FINAL.txt")
)

write_gene_list(
  BASELINE_FINAL_GENES,
  file.path(params$out_dir, "selected_gene_lists", "BASELINE_RECONSTRUCTED_FINAL_SIGNATURE.txt")
)
write_gene_list(
  setdiff(PRIMARY_SIGNATURE, BASELINE_FINAL_GENES),
  file.path(params$out_dir, "qc", "BASELINE_PRIMARY_GENES_MISSING_FROM_RECONSTRUCTION.txt")
)
write_gene_list(
  setdiff(BASELINE_FINAL_GENES, PRIMARY_SIGNATURE),
  file.path(params$out_dir, "qc", "BASELINE_UNEXPECTED_GENES_IN_RECONSTRUCTION.txt")
)

if (
  isTRUE(params$stop_if_baseline_reconstruction_mismatch) &&
    !BASELINE_EXACT_MATCH
) {
  stop(
    "The original thresholds did not reconstruct the established 24-gene set. ",
    "The sensitivity analysis has stopped to prevent testing a different ",
    "pipeline. Inspect 01_BASELINE_SUMMARY.csv and the baseline QC gene lists."
  )
}

message("STEP 7/9: Running all threshold-sensitivity analyses...")

OFAT_RUN <- run_settings_table(OFAT_SETTINGS, EVIDENCE, PRIMARY_SIGNATURE)
GTEX_GRID_RUN <- run_settings_table(GTEX_GRID_SETTINGS, EVIDENCE, PRIMARY_SIGNATURE)
HPA_GRID_RUN <- run_settings_table(HPA_GRID_SETTINGS, EVIDENCE, PRIMARY_SIGNATURE)
HPA_RULE_RUN <- run_settings_table(HPA_RULE_SETTINGS, EVIDENCE, PRIMARY_SIGNATURE)
ANTRUM_GRID_RUN <- run_settings_table(ANTRUM_GRID_SETTINGS, EVIDENCE, PRIMARY_SIGNATURE)
SSL_GRID_RUN <- run_settings_table(SSL_GRID_SETTINGS, EVIDENCE, PRIMARY_SIGNATURE)
COMBINED_RUN <- run_settings_table(COMBINED_SETTINGS, EVIDENCE, PRIMARY_SIGNATURE)

save_run_outputs(OFAT_RUN, "02_ONE_FACTOR_AT_A_TIME")
save_run_outputs(GTEX_GRID_RUN, "03_GTEX_WITHIN_STAGE_GRID")
save_run_outputs(HPA_GRID_RUN, "04_HPA_WITHIN_STAGE_GRID")
save_run_outputs(HPA_RULE_RUN, "05_HPA_RULE_COMPARISON")
save_run_outputs(ANTRUM_GRID_RUN, "06_ANTRUM_WITHIN_STAGE_GRID")
save_run_outputs(SSL_GRID_RUN, "07_SSL_HP_WITHIN_STAGE_GRID")
save_run_outputs(COMBINED_RUN, "08_COMBINED_REPRESENTATIVE_GRID")

for (setting_id in names(HPA_RULE_RUN$final_gene_lists)) {
  write_gene_list(
    HPA_RULE_RUN$final_gene_lists[[setting_id]],
    file.path(
      params$out_dir,
      "selected_gene_lists",
      paste0(safe_file_tag(setting_id), "_FINAL_SIGNATURE.txt")
    )
  )
}

for (setting_id in names(OFAT_RUN$final_gene_lists)) {
  write_gene_list(
    OFAT_RUN$final_gene_lists[[setting_id]],
    file.path(
      params$out_dir,
      "selected_gene_lists",
      paste0(safe_file_tag(setting_id), "_FINAL_SIGNATURE.txt")
    )
  )
}

message("STEP 8/9: Summarising gene robustness and identifying sensitive genes...")

ALL_SUMMARIES <- dplyr::bind_rows(
  BASELINE_RUN$summary,
  OFAT_RUN$summary,
  GTEX_GRID_RUN$summary,
  HPA_GRID_RUN$summary,
  HPA_RULE_RUN$summary,
  ANTRUM_GRID_RUN$summary,
  SSL_GRID_RUN$summary,
  COMBINED_RUN$summary
)

ALL_PRIMARY_MEMBERSHIP <- dplyr::bind_rows(
  BASELINE_RUN$primary_membership,
  OFAT_RUN$primary_membership,
  GTEX_GRID_RUN$primary_membership,
  HPA_GRID_RUN$primary_membership,
  HPA_RULE_RUN$primary_membership,
  ANTRUM_GRID_RUN$primary_membership,
  SSL_GRID_RUN$primary_membership,
  COMBINED_RUN$primary_membership
)

ALL_FINAL_MEMBERSHIP <- dplyr::bind_rows(
  BASELINE_RUN$final_membership,
  OFAT_RUN$final_membership,
  GTEX_GRID_RUN$final_membership,
  HPA_GRID_RUN$final_membership,
  HPA_RULE_RUN$final_membership,
  ANTRUM_GRID_RUN$final_membership,
  SSL_GRID_RUN$final_membership,
  COMBINED_RUN$final_membership
)

write_csv_overwrite(
  ALL_SUMMARIES,
  file.path(params$out_dir, "tables", "09_ALL_SETTINGS_SUMMARY.csv")
)
write_csv_overwrite(
  ALL_PRIMARY_MEMBERSHIP,
  file.path(params$out_dir, "tables", "10_ALL_SETTINGS_PRIMARY_24_MEMBERSHIP.csv")
)
write_csv_overwrite(
  ALL_FINAL_MEMBERSHIP,
  file.path(params$out_dir, "tables", "11_ALL_SETTINGS_FINAL_GENE_MEMBERSHIP.csv")
)

FAMILY_IMPACT <- ALL_SUMMARIES %>%
  dplyr::group_by(analysis_family) %>%
  dplyr::summarise(
    n_settings_tested = dplyr::n(),
    n_settings_exactly_reproducing_primary_24 = sum(exact_match_to_primary_24),
    percentage_settings_exactly_reproducing_primary_24 =
      100 * mean(exact_match_to_primary_24),
    minimum_primary_genes_retained = min(n_primary_24_retained),
    maximum_primary_genes_retained = max(n_primary_24_retained),
    minimum_primary_genes_removed = min(n_primary_24_removed),
    maximum_primary_genes_removed = max(n_primary_24_removed),
    minimum_additional_genes = min(n_additional_genes_not_in_primary),
    maximum_additional_genes = max(n_additional_genes_not_in_primary),
    minimum_final_list_size = min(n_final),
    maximum_final_list_size = max(n_final),
    minimum_genes_after_GTEx = min(n_after_GTEx),
    maximum_genes_after_GTEx = max(n_after_GTEx),
    minimum_genes_after_HPA = min(n_after_HPA),
    maximum_genes_after_HPA = max(n_after_HPA),
    minimum_genes_after_ANTRUM = min(n_after_ANTRUM),
    maximum_genes_after_ANTRUM = max(n_after_ANTRUM),
    .groups = "drop"
  )

write_csv_overwrite(
  FAMILY_IMPACT,
  file.path(params$out_dir, "tables", "12_SENSITIVITY_FAMILY_IMPACT_SUMMARY.csv")
)

OFAT_PARAMETER_IMPACT <- OFAT_RUN$summary %>%
  dplyr::group_by(varied_parameter) %>%
  dplyr::summarise(
    n_values_tested = dplyr::n(),
    n_exact_primary_24 = sum(exact_match_to_primary_24),
    minimum_primary_genes_retained = min(n_primary_24_retained),
    maximum_primary_genes_retained = max(n_primary_24_retained),
    maximum_primary_genes_removed = max(n_primary_24_removed),
    maximum_additional_genes = max(n_additional_genes_not_in_primary),
    minimum_final_list_size = min(n_final),
    maximum_final_list_size = max(n_final),
    settings_changing_primary_24 = sum(!exact_match_to_primary_24),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    minimum_primary_genes_retained,
    dplyr::desc(settings_changing_primary_24),
    varied_parameter
  )

write_csv_overwrite(
  OFAT_PARAMETER_IMPACT,
  file.path(params$out_dir, "tables", "13_ONE_FACTOR_PARAMETER_IMPACT_RANKING.csv")
)

PRIMARY_RETENTION_BY_FAMILY <- ALL_PRIMARY_MEMBERSHIP %>%
  dplyr::filter(analysis_family != "BASELINE") %>%
  dplyr::group_by(gene_symbol, analysis_family) %>%
  dplyr::summarise(
    n_settings_tested = dplyr::n(),
    n_settings_retained = sum(pass_final),
    retention_percentage = 100 * mean(pass_final),
    n_failed_at_GTEx = sum(first_stage_failed == "GTEx"),
    n_failed_at_HPA = sum(first_stage_failed == "HPA"),
    n_failed_at_ANTRUM = sum(first_stage_failed == "ANTRUM"),
    n_failed_at_SSL_HP = sum(first_stage_failed == "SSL_HP"),
    .groups = "drop"
  )

write_csv_overwrite(
  PRIMARY_RETENTION_BY_FAMILY,
  file.path(params$out_dir, "tables", "14_PRIMARY_24_RETENTION_BY_SENSITIVITY_FAMILY.csv")
)

PRIMARY_RETENTION_WIDE <- PRIMARY_RETENTION_BY_FAMILY %>%
  dplyr::select(gene_symbol, analysis_family, retention_percentage) %>%
  tidyr::pivot_wider(
    names_from = analysis_family,
    values_from = retention_percentage,
    names_prefix = "retention_pct_"
  )

PRIMARY_OVERALL_RETENTION <- ALL_PRIMARY_MEMBERSHIP %>%
  dplyr::filter(analysis_family != "BASELINE") %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::summarise(
    n_all_settings_tested = dplyr::n(),
    n_all_settings_retained = sum(pass_final),
    overall_descriptive_retention_percentage = 100 * mean(pass_final),
    .groups = "drop"
  )

baseline_passes <- compute_stage_passes(EVIDENCE, BASELINE_SETTINGS[1, , drop = FALSE])

BASELINE_PRIMARY_EVIDENCE <- EVIDENCE %>%
  dplyr::mutate(pass_baseline_final = baseline_passes$pass_final) %>%
  dplyr::filter(pass_baseline_final, gene_uc %in% PRIMARY_SIGNATURE) %>%
  dplyr::arrange(
    gene_uc,
    dplyr::desc(ANTRUM_vs_H_COL_log2FC),
    dplyr::desc(GTEx_stomach_vs_colon_log2FC),
    dplyr::desc(tau),
    dplyr::desc(HPA_stomach_vs_colon_log2FC_original)
  ) %>%
  dplyr::distinct(gene_uc, .keep_all = TRUE) %>%
  dplyr::mutate(primary_order = match(gene_uc, PRIMARY_SIGNATURE)) %>%
  dplyr::arrange(primary_order)

if (!setequal(BASELINE_PRIMARY_EVIDENCE$gene_uc, PRIMARY_SIGNATURE)) {
  stop("Could not obtain one baseline-passing evidence row for every primary gene.")
}

BASELINE_MARGIN_TABLE <- BASELINE_PRIMARY_EVIDENCE %>%
  dplyr::transmute(
    primary_order,
    gene_symbol = gene_uc,

    tau,
    GTEx_stomach_mean_TPM,
    GTEx_colon_max_TPM,
    GTEx_stomach_vs_colon_log2FC,
    HPA_stomach_nTPM,
    HPA_colon_max_nTPM,
    HPA_stomach_vs_colon_log2FC = HPA_stomach_vs_colon_log2FC_original,
    HPA_highest_other_tissue,
    HPA_highest_other_nTPM,
    HPA_stomach_vs_highest_other_log2FC,
    ANTRUM_vs_H_COL_log2FC,
    ANTRUM_vs_H_COL_FDR,
    SSL_HP_vs_H_COL_log2FC,
    SSL_HP_vs_H_COL_adjusted_p,

    margin_GTEx_tau = tau - params$TAU_MIN,
    margin_GTEx_stomach_TPM = GTEx_stomach_mean_TPM - params$GTEX_STOMACH_TPM_MIN,
    margin_GTEx_colon_TPM = params$GTEX_COLON_TPM_MAX - GTEx_colon_max_TPM,
    margin_GTEx_log2FC = GTEx_stomach_vs_colon_log2FC - params$GTEX_LOG2FC_MIN,

    margin_HPA_stomach_nTPM = HPA_stomach_nTPM - params$HPA_STOMACH_NTPM_MIN,
    margin_HPA_colon_nTPM = params$HPA_COLON_NTPM_MAX - HPA_colon_max_nTPM,
    margin_HPA_log2FC = HPA_stomach_vs_colon_log2FC_original - params$HPA_LOG2FC_MIN,

    margin_ANTRUM_log2FC = ANTRUM_vs_H_COL_log2FC - params$ANTRUM_LOG2FC_MIN,
    margin_ANTRUM_FDR = params$ANTRUM_FDR_MAX - ANTRUM_vs_H_COL_FDR,

    margin_SSL_HP_log2FC = SSL_HP_vs_H_COL_log2FC - params$SSL_HP_LOG2FC_MIN,
    margin_SSL_HP_FDR = params$SSL_HP_FDR_MAX - SSL_HP_vs_H_COL_significance_value
  )

BASELINE_MARGIN_TABLE <- BASELINE_MARGIN_TABLE %>%
  dplyr::mutate(
    relative_margin_GTEx_tau = tau / params$TAU_MIN - 1,
    relative_margin_GTEx_stomach = GTEx_stomach_mean_TPM / params$GTEX_STOMACH_TPM_MIN - 1,
    relative_margin_GTEx_colon = 1 - GTEx_colon_max_TPM / params$GTEX_COLON_TPM_MAX,
    relative_margin_GTEx_log2FC = GTEx_stomach_vs_colon_log2FC / params$GTEX_LOG2FC_MIN - 1,

    relative_margin_HPA_stomach = HPA_stomach_nTPM / params$HPA_STOMACH_NTPM_MIN - 1,
    relative_margin_HPA_colon = 1 - HPA_colon_max_nTPM / params$HPA_COLON_NTPM_MAX,
    relative_margin_HPA_log2FC = HPA_stomach_vs_colon_log2FC / params$HPA_LOG2FC_MIN - 1,

    relative_margin_ANTRUM_log2FC = ANTRUM_vs_H_COL_log2FC / params$ANTRUM_LOG2FC_MIN - 1,
    relative_margin_ANTRUM_FDR = 1 - ANTRUM_vs_H_COL_FDR / params$ANTRUM_FDR_MAX,

    relative_margin_SSL_HP_log2FC = SSL_HP_vs_H_COL_log2FC,
    relative_margin_SSL_HP_FDR = 1 - SSL_HP_vs_H_COL_adjusted_p / params$SSL_HP_FDR_MAX
  )

relative_margin_cols <- c(
  "relative_margin_GTEx_tau",
  "relative_margin_GTEx_stomach",
  "relative_margin_GTEx_colon",
  "relative_margin_GTEx_log2FC",
  "relative_margin_HPA_stomach",
  "relative_margin_HPA_colon",
  "relative_margin_HPA_log2FC",
  "relative_margin_ANTRUM_log2FC",
  "relative_margin_ANTRUM_FDR",
  "relative_margin_SSL_HP_log2FC",
  "relative_margin_SSL_HP_FDR"
)

relative_margin_matrix <- as.matrix(
  BASELINE_MARGIN_TABLE[, relative_margin_cols, drop = FALSE]
)
storage.mode(relative_margin_matrix) <- "double"

BASELINE_MARGIN_TABLE$minimum_relative_margin <- apply(
  relative_margin_matrix,
  1,
  function(x) min(x[is.finite(x)], na.rm = TRUE)
)

BASELINE_MARGIN_TABLE$closest_original_filter <- apply(
  relative_margin_matrix,
  1,
  function(x) {
    if (!any(is.finite(x))) return(NA_character_)
    relative_margin_cols[which.min(ifelse(is.finite(x), x, Inf))]
  }
)

write_csv_overwrite(
  BASELINE_MARGIN_TABLE,
  file.path(params$out_dir, "tables", "15_PRIMARY_24_BASELINE_FILTER_MARGINS.csv")
)

GENE_ROBUSTNESS_RANKING <- PRIMARY_RETENTION_WIDE %>%
  dplyr::left_join(PRIMARY_OVERALL_RETENTION, by = "gene_symbol") %>%
  dplyr::left_join(
    BASELINE_MARGIN_TABLE %>%
      dplyr::select(
        gene_symbol,
        primary_order,
        minimum_relative_margin,
        closest_original_filter,
        ANTRUM_vs_H_COL_log2FC,
        ANTRUM_vs_H_COL_FDR,
        SSL_HP_vs_H_COL_log2FC,
        SSL_HP_vs_H_COL_adjusted_p,
        HPA_stomach_nTPM,
        HPA_colon_max_nTPM,
        HPA_stomach_vs_colon_log2FC,
        HPA_highest_other_tissue,
        HPA_stomach_vs_highest_other_log2FC
      ),
    by = "gene_symbol"
  )

retention_columns <- grep("^retention_pct_", names(GENE_ROBUSTNESS_RANKING), value = TRUE)

GENE_ROBUSTNESS_RANKING$minimum_family_retention_percentage <- apply(
  as.matrix(GENE_ROBUSTNESS_RANKING[, retention_columns, drop = FALSE]),
  1,
  function(x) min(x[is.finite(x)], na.rm = TRUE)
)

GENE_ROBUSTNESS_RANKING <- GENE_ROBUSTNESS_RANKING %>%
  dplyr::arrange(
    minimum_family_retention_percentage,
    overall_descriptive_retention_percentage,
    minimum_relative_margin,
    primary_order
  ) %>%
  dplyr::mutate(sensitivity_rank_1_most_sensitive = dplyr::row_number())

write_csv_overwrite(
  GENE_ROBUSTNESS_RANKING,
  file.path(params$out_dir, "tables", "16_PRIMARY_24_GENE_ROBUSTNESS_RANKING.csv")
)

MOST_SENSITIVE_10 <- GENE_ROBUSTNESS_RANKING %>%
  dplyr::slice_head(n = min(10L, nrow(GENE_ROBUSTNESS_RANKING)))

write_csv_overwrite(
  MOST_SENSITIVE_10,
  file.path(params$out_dir, "tables", "17_TEN_MOST_THRESHOLD_SENSITIVE_PRIMARY_GENES.csv")
)

CHANGED_SETTINGS <- ALL_SUMMARIES %>%
  dplyr::filter(!exact_match_to_primary_24) %>%
  dplyr::arrange(
    n_primary_24_retained,
    n_additional_genes_not_in_primary,
    analysis_family,
    setting_id
  )

write_csv_overwrite(
  CHANGED_SETTINGS,
  file.path(params$out_dir, "tables", "18_SETTINGS_THAT_CHANGED_THE_24_GENE_LIST.csv")
)

EXACT_SETTINGS <- ALL_SUMMARIES %>%
  dplyr::filter(exact_match_to_primary_24) %>%
  dplyr::arrange(analysis_family, setting_id)

write_csv_overwrite(
  EXACT_SETTINGS,
  file.path(params$out_dir, "tables", "19_SETTINGS_THAT_REPRODUCED_THE_EXACT_24_GENES.csv")
)

retention_plot_df <- PRIMARY_RETENTION_BY_FAMILY %>%
  dplyr::mutate(
    gene_symbol = factor(
      gene_symbol,
      levels = rev(GENE_ROBUSTNESS_RANKING$gene_symbol)
    ),
    analysis_family = factor(
      analysis_family,
      levels = c(
        "ONE_FACTOR_AT_A_TIME",
        "GTEX_WITHIN_STAGE_GRID",
        "HPA_WITHIN_STAGE_GRID",
        "HPA_RULE_COMPARISON",
        "ANTRUM_WITHIN_STAGE_GRID",
        "SSL_HP_WITHIN_STAGE_GRID",
        "COMBINED_REPRESENTATIVE_GRID"
      )
    )
  )

p_retention <- ggplot2::ggplot(
  retention_plot_df,
  ggplot2::aes(
    x = analysis_family,
    y = gene_symbol,
    fill = retention_percentage
  )
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.0f", retention_percentage)),
    size = 2.5
  ) +
  ggplot2::scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(0, 100),
    name = "Retained (%)"
  ) +
  ggplot2::labs(
    title = "Retention of the established 24 genes across sensitivity analyses",
    subtitle = "Numbers show the percentage of tested settings retaining each gene",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    panel.grid = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(params$out_dir, "plots", "PRIMARY_24_RETENTION_BY_FAMILY.pdf"),
  p_retention,
  width = 12,
  height = 8.5
)
ggplot2::ggsave(
  file.path(params$out_dir, "plots", "PRIMARY_24_RETENTION_BY_FAMILY.png"),
  p_retention,
  width = 12,
  height = 8.5,
  dpi = 320
)

ofat_plot_df <- OFAT_RUN$summary %>%
  dplyr::mutate(
    varied_value_plot = factor(
      paste(varied_parameter, varied_value_label, sep = "__"),
      levels = unique(paste(varied_parameter, varied_value_label, sep = "__")),
      labels = varied_value_label
    )
  )

p_ofat <- ggplot2::ggplot(
  ofat_plot_df,
  ggplot2::aes(
    x = varied_value_plot,
    y = n_primary_24_retained,
    group = 1
  )
) +
  ggplot2::geom_hline(yintercept = 24, linetype = 2) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~varied_parameter, scales = "free_x", ncol = 3) +
  ggplot2::scale_y_continuous(limits = c(0, 24), breaks = seq(0, 24, by = 4)) +
  ggplot2::labs(
    title = "One-factor-at-a-time threshold sensitivity",
    subtitle = "All other filters were fixed at their original values",
    x = "Tested value",
    y = "Established genes retained"
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    panel.grid = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(params$out_dir, "plots", "ONE_FACTOR_AT_A_TIME_PRIMARY_GENE_RETENTION.pdf"),
  p_ofat,
  width = 12,
  height = 12
)
ggplot2::ggsave(
  file.path(params$out_dir, "plots", "ONE_FACTOR_AT_A_TIME_PRIMARY_GENE_RETENTION.png"),
  p_ofat,
  width = 12,
  height = 12,
  dpi = 320
)

KEY_OFAT_PARAMETERS <- c(
  "gtex_tau_min",
  "gtex_log2fc_min",
  "hpa_log2fc_min",
  "antrum_log2fc_min",
  "ssl_hp_log2fc_min"
)

KEY_OFAT_LABELS <- c(
  gtex_tau_min = "GTEx tau minimum",
  gtex_log2fc_min = "GTEx stomach-versus-colon log2FC minimum",
  hpa_log2fc_min = "HPA stomach-versus-colon log2FC minimum",
  antrum_log2fc_min = "OSLC antrum-versus-H_COL log2FC minimum",
  ssl_hp_log2fc_min = "SSL_HP-versus-H_COL log2FC minimum"
)

KEY_OFAT_COMPOSITION <- OFAT_RUN$summary %>%
  dplyr::filter(varied_parameter %in% KEY_OFAT_PARAMETERS) %>%
  dplyr::mutate(
    parameter_label = unname(KEY_OFAT_LABELS[varied_parameter]),
    parameter_label = factor(
      parameter_label,
      levels = unname(KEY_OFAT_LABELS[KEY_OFAT_PARAMETERS])
    ),
    tested_value = varied_value_label,
    exact_same_24 = exact_match_to_primary_24,
    arithmetic_check =
      n_final == n_primary_24_retained + n_additional_genes_not_in_primary &
      24L == n_primary_24_retained + n_primary_24_removed
  ) %>%
  dplyr::select(
    setting_id,
    parameter_label,
    varied_parameter,
    tested_value,
    n_final,
    n_primary_24_retained,
    n_primary_24_removed,
    n_additional_genes_not_in_primary,
    exact_same_24,
    jaccard_with_primary_24,
    primary_genes_removed,
    additional_genes_not_in_primary,
    final_gene_list,
    arithmetic_check
  )

if (!all(KEY_OFAT_COMPOSITION$arithmetic_check)) {
  stop("Arithmetic check failed in key one-factor final-list composition table.")
}

write_csv_overwrite(
  KEY_OFAT_COMPOSITION,
  file.path(
    params$out_dir,
    "tables",
    "22_KEY_ONE_FACTOR_FINAL_LIST_COMPOSITION_AND_GENE_CHANGES.csv"
  )
)

OFAT_SETTING_METADATA <- OFAT_RUN$summary %>%
  dplyr::select(
    setting_id,
    varied_parameter,
    varied_value_label
  ) %>%
  dplyr::distinct()

KEY_OFAT_GENE_CHANGES_LONG <- dplyr::bind_rows(
  OFAT_RUN$final_membership %>%
    dplyr::left_join(OFAT_SETTING_METADATA, by = "setting_id") %>%
    dplyr::filter(
      varied_parameter %in% KEY_OFAT_PARAMETERS,
      !in_established_primary_24
    ) %>%
    dplyr::transmute(
      setting_id,
      varied_parameter,
      tested_value = varied_value_label,
      gene_symbol,
      change_type = "Additional gene included"
    ),
  OFAT_RUN$primary_membership %>%
    dplyr::left_join(OFAT_SETTING_METADATA, by = "setting_id") %>%
    dplyr::filter(
      varied_parameter %in% KEY_OFAT_PARAMETERS,
      !pass_final
    ) %>%
    dplyr::transmute(
      setting_id,
      varied_parameter,
      tested_value = varied_value_label,
      gene_symbol,
      change_type = "Established gene removed"
    )
) %>%
  dplyr::mutate(
    parameter_label = unname(KEY_OFAT_LABELS[varied_parameter]),
    parameter_label = factor(
      parameter_label,
      levels = unname(KEY_OFAT_LABELS[KEY_OFAT_PARAMETERS])
    )
  ) %>%
  dplyr::arrange(parameter_label, tested_value, change_type, gene_symbol)

write_csv_overwrite(
  KEY_OFAT_GENE_CHANGES_LONG,
  file.path(
    params$out_dir,
    "tables",
    "23_KEY_ONE_FACTOR_ADDED_AND_REMOVED_GENES_LONG.csv"
  )
)

key_ofat_order <- OFAT_RUN$summary %>%
  dplyr::filter(varied_parameter %in% KEY_OFAT_PARAMETERS) %>%
  dplyr::transmute(
    varied_parameter,
    varied_value_label,
    value_key = paste(varied_parameter, varied_value_label, sep = "__")
  ) %>%
  dplyr::distinct() %>%
  dplyr::pull(value_key)

key_ofat_comp_long <- KEY_OFAT_COMPOSITION %>%
  dplyr::mutate(
    value_key = factor(
      paste(varied_parameter, tested_value, sep = "__"),
      levels = key_ofat_order
    )
  ) %>%
  tidyr::pivot_longer(
    cols = c(n_primary_24_retained, n_additional_genes_not_in_primary),
    names_to = "component",
    values_to = "n_genes"
  ) %>%
  dplyr::mutate(
    component = factor(
      component,
      levels = c(
        "n_primary_24_retained",
        "n_additional_genes_not_in_primary"
      ),
      labels = c(
        "Established genes retained",
        "Additional genes included"
      )
    )
  )

key_ofat_label_df <- KEY_OFAT_COMPOSITION %>%
  dplyr::mutate(
    value_key = factor(
      paste(varied_parameter, tested_value, sep = "__"),
      levels = key_ofat_order
    ),
    result_label = paste0(
      "Final n=", n_final,
      "\nretained ", n_primary_24_retained,
      ", added ", n_additional_genes_not_in_primary,
      ", removed ", n_primary_24_removed
    )
  )

p_key_ofat_composition <- ggplot2::ggplot(
  key_ofat_comp_long,
  ggplot2::aes(x = value_key, y = n_genes, fill = component)
) +
  ggplot2::geom_hline(yintercept = 24, linetype = 2) +
  ggplot2::geom_col(colour = "black", linewidth = 0.3) +
  ggplot2::geom_text(
    data = key_ofat_label_df,
    ggplot2::aes(
      x = value_key,
      y = n_final,
      label = result_label
    ),
    inherit.aes = FALSE,
    vjust = -0.25,
    size = 2.7
  ) +
  ggplot2::facet_wrap(~parameter_label, scales = "free_x", ncol = 2) +
  ggplot2::scale_fill_manual(
    values = c(
      "Established genes retained" = "grey45",
      "Additional genes included" = "grey80"
    )
  ) +
  ggplot2::scale_x_discrete(
    labels = function(x) sub("^.*__", "", x)
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.22))
  ) +
  ggplot2::labs(
    title = "Sensitivity of final signature composition to key thresholds",
    subtitle = paste0(
      "Each threshold was varied separately; all other criteria were held at their original values. ",
      "Bars include both established genes retained and additional genes admitted."
    ),
    x = "Tested threshold",
    y = "Number of genes in final list",
    fill = NULL
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
    legend.position = "top",
    panel.grid = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(
    params$out_dir,
    "plots",
    "KEY_ONE_FACTOR_FINAL_SIGNATURE_COMPOSITION.pdf"
  ),
  p_key_ofat_composition,
  width = 12.5,
  height = 10
)

ggplot2::ggsave(
  file.path(
    params$out_dir,
    "plots",
    "KEY_ONE_FACTOR_FINAL_SIGNATURE_COMPOSITION.png"
  ),
  p_key_ofat_composition,
  width = 12.5,
  height = 10,
  dpi = 320
)

key_ofat_changes_long <- dplyr::bind_rows(
  KEY_OFAT_COMPOSITION %>%
    dplyr::transmute(
      setting_id,
      parameter_label,
      varied_parameter,
      tested_value,
      change_type = "Established genes removed",
      n_changed = -n_primary_24_removed
    ),
  KEY_OFAT_COMPOSITION %>%
    dplyr::transmute(
      setting_id,
      parameter_label,
      varied_parameter,
      tested_value,
      change_type = "Additional genes included",
      n_changed = n_additional_genes_not_in_primary
    )
) %>%
  dplyr::mutate(
    value_key = factor(
      paste(varied_parameter, tested_value, sep = "__"),
      levels = key_ofat_order
    )
  )

p_key_ofat_changes <- ggplot2::ggplot(
  key_ofat_changes_long,
  ggplot2::aes(x = value_key, y = n_changed, fill = change_type)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
  ggplot2::geom_col(colour = "black", linewidth = 0.25) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = ifelse(n_changed == 0, "", abs(n_changed)),
      vjust = ifelse(n_changed >= 0, -0.25, 1.15)
    ),
    size = 2.8
  ) +
  ggplot2::facet_wrap(~parameter_label, scales = "free_x", ncol = 2) +
  ggplot2::scale_fill_manual(
    values = c(
      "Established genes removed" = "grey40",
      "Additional genes included" = "grey80"
    )
  ) +
  ggplot2::scale_x_discrete(labels = function(x) sub("^.*__", "", x)) +
  ggplot2::labs(
    title = "Genes added to or removed from the established signature",
    subtitle = "Negative values indicate established genes removed; positive values indicate additional genes included",
    x = "Tested threshold",
    y = "Number of genes changed",
    fill = NULL
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
    legend.position = "top",
    panel.grid = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(
    params$out_dir,
    "plots",
    "KEY_ONE_FACTOR_ADDITIONS_AND_REMOVALS.pdf"
  ),
  p_key_ofat_changes,
  width = 12.5,
  height = 10
)

ggplot2::ggsave(
  file.path(
    params$out_dir,
    "plots",
    "KEY_ONE_FACTOR_ADDITIONS_AND_REMOVALS.png"
  ),
  p_key_ofat_changes,
  width = 12.5,
  height = 10,
  dpi = 320
)

HPA_RULE_COMPOSITION <- HPA_RULE_RUN$summary %>%
  dplyr::transmute(
    setting_id,
    setting_label,
    n_final,
    n_established_24_retained = n_primary_24_retained,
    n_established_24_removed = n_primary_24_removed,
    n_new_genes_added = n_additional_genes_not_in_primary,
    exact_same_24_gene_list = exact_match_to_primary_24,
    jaccard_with_established_24 = jaccard_with_primary_24,
    established_genes_removed = primary_genes_removed,
    new_genes_added = additional_genes_not_in_primary,
    final_gene_list
  ) %>%
  dplyr::mutate(
    arithmetic_check =
      n_final == n_established_24_retained + n_new_genes_added
  )

if (!all(HPA_RULE_COMPOSITION$arithmetic_check)) {
  stop(
    "Internal HPA comparison error: n_final does not equal retained established genes plus newly included genes."
  )
}

write_csv_overwrite(
  HPA_RULE_COMPOSITION,
  file.path(
    params$out_dir,
    "tables",
    "20_HPA_RULE_FINAL_LIST_COMPOSITION_AND_GENE_CHANGES.csv"
  )
)

HPA_RULE_GENE_CHANGES_LONG <- dplyr::bind_rows(
  HPA_RULE_RUN$final_membership %>%
    dplyr::filter(!in_established_primary_24) %>%
    dplyr::transmute(
      setting_id,
      setting_label,
      gene_symbol,
      change_type = "Newly included under alternative rule"
    ),
  HPA_RULE_RUN$primary_membership %>%
    dplyr::filter(!pass_final) %>%
    dplyr::transmute(
      setting_id,
      setting_label,
      gene_symbol,
      change_type = "Established gene removed under alternative rule"
    )
) %>%
  dplyr::arrange(setting_id, change_type, gene_symbol)

write_csv_overwrite(
  HPA_RULE_GENE_CHANGES_LONG,
  file.path(
    params$out_dir,
    "tables",
    "21_HPA_RULE_ADDED_AND_REMOVED_GENES_LONG.csv"
  )
)

hpa_plot_order <- rev(HPA_RULE_COMPOSITION$setting_label)

hpa_composition_long <- HPA_RULE_COMPOSITION %>%
  dplyr::select(
    setting_id,
    setting_label,
    n_established_24_retained,
    n_new_genes_added
  ) %>%
  tidyr::pivot_longer(
    cols = c(n_established_24_retained, n_new_genes_added),
    names_to = "gene_category",
    values_to = "n_genes"
  ) %>%
  dplyr::mutate(
    setting_label = factor(setting_label, levels = hpa_plot_order),
    gene_category = factor(
      gene_category,
      levels = c("n_established_24_retained", "n_new_genes_added"),
      labels = c("Established 24 retained", "Additional genes included")
    ),
    segment_label = dplyr::if_else(n_genes > 0, as.character(n_genes), "")
  )

hpa_label_df <- HPA_RULE_COMPOSITION %>%
  dplyr::mutate(
    setting_label = factor(setting_label, levels = hpa_plot_order),
    end_label = paste0(
      "Final n = ", n_final,
      "  |  retained ", n_established_24_retained,
      ", added ", n_new_genes_added,
      ", removed ", n_established_24_removed
    )
  )

hpa_x_max <- max(hpa_label_df$n_final, 24, na.rm = TRUE)
hpa_x_padding <- max(4, 0.38 * hpa_x_max)

p_hpa_composition <- ggplot2::ggplot(
  hpa_composition_long,
  ggplot2::aes(
    x = n_genes,
    y = setting_label,
    fill = gene_category
  )
) +
  ggplot2::geom_vline(xintercept = 24, linetype = 2) +
  ggplot2::geom_col(colour = "black", linewidth = 0.35) +
  ggplot2::geom_text(
    ggplot2::aes(label = segment_label),
    position = ggplot2::position_stack(vjust = 0.5),
    size = 3.2
  ) +
  ggplot2::geom_text(
    data = hpa_label_df,
    ggplot2::aes(
      x = n_final,
      y = setting_label,
      label = end_label
    ),
    inherit.aes = FALSE,
    hjust = -0.04,
    size = 3.1
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Established 24 retained" = "grey45",
      "Additional genes included" = "grey82"
    ),
    drop = FALSE
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, hpa_x_max + hpa_x_padding),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::labs(
    title = "Effect of alternative HPA validation rules on final-list composition",
    subtitle = paste0(
      "GTEx, ANTRUM and SSL_HP criteria were held at their original values. ",
      "The dashed line marks the original signature size (24 genes)."
    ),
    x = "Genes in the final list",
    y = NULL,
    fill = NULL
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    legend.position = "top",
    panel.grid = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(5.5, 130, 5.5, 5.5)
  ) +
  ggplot2::coord_cartesian(clip = "off")

ggplot2::ggsave(
  file.path(
    params$out_dir,
    "plots",
    "HPA_RULE_COMPARISON_FINAL_SIGNATURE_COMPOSITION.pdf"
  ),
  p_hpa_composition,
  width = 13.5,
  height = 7
)
ggplot2::ggsave(
  file.path(
    params$out_dir,
    "plots",
    "HPA_RULE_COMPARISON_FINAL_SIGNATURE_COMPOSITION.png"
  ),
  p_hpa_composition,
  width = 13.5,
  height = 7,
  dpi = 320
)

hpa_change_long <- dplyr::bind_rows(
  HPA_RULE_COMPOSITION %>%
    dplyr::transmute(
      setting_id,
      setting_label,
      change_type = "Additional genes included",
      n_change = n_new_genes_added
    ),
  HPA_RULE_COMPOSITION %>%
    dplyr::transmute(
      setting_id,
      setting_label,
      change_type = "Established genes removed",
      n_change = -n_established_24_removed
    )
) %>%
  dplyr::mutate(
    setting_label = factor(setting_label, levels = hpa_plot_order),
    change_type = factor(
      change_type,
      levels = c("Established genes removed", "Additional genes included")
    ),
    change_label = dplyr::if_else(
      n_change == 0,
      "",
      as.character(abs(n_change))
    )
  )

hpa_change_limit <- max(abs(hpa_change_long$n_change), 1, na.rm = TRUE)

p_hpa_changes <- ggplot2::ggplot(
  hpa_change_long,
  ggplot2::aes(
    x = n_change,
    y = setting_label,
    fill = change_type
  )
) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.45) +
  ggplot2::geom_col(colour = "black", linewidth = 0.3) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = change_label,
      hjust = dplyr::if_else(n_change < 0, 1.15, -0.15)
    ),
    size = 3.2
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Established genes removed" = "grey35",
      "Additional genes included" = "grey80"
    ),
    drop = FALSE
  ) +
  ggplot2::scale_x_continuous(
    limits = c(-1.30 * hpa_change_limit, 1.30 * hpa_change_limit),
    breaks = scales::pretty_breaks(n = 7),
    labels = function(x) abs(x)
  ) +
  ggplot2::labs(
    title = "Genes added to or removed from the established 24-gene list",
    subtitle = "Left: established genes removed; right: additional genes included",
    x = "Number of genes changed",
    y = NULL,
    fill = NULL
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    legend.position = "top",
    panel.grid = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(
    params$out_dir,
    "plots",
    "HPA_RULE_COMPARISON_ADDITIONS_AND_REMOVALS.pdf"
  ),
  p_hpa_changes,
  width = 11.5,
  height = 7
)
ggplot2::ggsave(
  file.path(
    params$out_dir,
    "plots",
    "HPA_RULE_COMPARISON_ADDITIONS_AND_REMOVALS.png"
  ),
  p_hpa_changes,
  width = 11.5,
  height = 7,
  dpi = 320
)

margin_plot_df <- BASELINE_MARGIN_TABLE %>%
  dplyr::mutate(
    gene_symbol = factor(
      gene_symbol,
      levels = BASELINE_MARGIN_TABLE$gene_symbol[
        order(BASELINE_MARGIN_TABLE$minimum_relative_margin)
      ]
    )
  )

p_margin <- ggplot2::ggplot(
  margin_plot_df,
  ggplot2::aes(x = gene_symbol, y = minimum_relative_margin)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::geom_col(fill = "grey70", colour = "black") +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Distance of each primary gene from its closest original filter",
    subtitle = "Smaller values indicate a gene closer to at least one original cut-off",
    x = NULL,
    y = "Minimum relative margin across original filters"
  ) +
  ggplot2::theme_classic(base_size = 11)

ggplot2::ggsave(
  file.path(params$out_dir, "plots", "PRIMARY_24_CLOSEST_ORIGINAL_FILTER_MARGIN.pdf"),
  p_margin,
  width = 8,
  height = 7
)
ggplot2::ggsave(
  file.path(params$out_dir, "plots", "PRIMARY_24_CLOSEST_ORIGINAL_FILTER_MARGIN.png"),
  p_margin,
  width = 8,
  height = 7,
  dpi = 320
)

family_lines <- character(nrow(FAMILY_IMPACT))
for (i in seq_len(nrow(FAMILY_IMPACT))) {
  family_lines[i] <- paste0(
    "- ", FAMILY_IMPACT$analysis_family[i], ": ",
    FAMILY_IMPACT$n_settings_tested[i], " settings; ",
    FAMILY_IMPACT$n_settings_exactly_reproducing_primary_24[i],
    " reproduced the exact 24 genes; minimum retained = ",
    FAMILY_IMPACT$minimum_primary_genes_retained[i], "/24; final-list range = ",
    FAMILY_IMPACT$minimum_final_list_size[i], "-",
    FAMILY_IMPACT$maximum_final_list_size[i], "."
  )
}

most_sensitive_gene_lines <- character(nrow(MOST_SENSITIVE_10))
for (i in seq_len(nrow(MOST_SENSITIVE_10))) {
  most_sensitive_gene_lines[i] <- paste0(
    i, ". ", MOST_SENSITIVE_10$gene_symbol[i],
    " -- minimum family retention ",
    sprintf("%.1f", MOST_SENSITIVE_10$minimum_family_retention_percentage[i]),
    "%; overall descriptive retention ",
    sprintf("%.1f", MOST_SENSITIVE_10$overall_descriptive_retention_percentage[i]),
    "%; closest original filter: ",
    MOST_SENSITIVE_10$closest_original_filter[i], "."
  )
}

most_influential_parameter_lines <- character(nrow(OFAT_PARAMETER_IMPACT))
for (i in seq_len(nrow(OFAT_PARAMETER_IMPACT))) {
  most_influential_parameter_lines[i] <- paste0(
    "- ", OFAT_PARAMETER_IMPACT$varied_parameter[i],
    ": minimum ", OFAT_PARAMETER_IMPACT$minimum_primary_genes_retained[i],
    "/24 retained; ", OFAT_PARAMETER_IMPACT$settings_changing_primary_24[i],
    " tested values changed the exact list."
  )
}

interpretation_lines <- c(
  "FULL SIGNATURE THRESHOLD SENSITIVITY AUDIT",
  "==========================================",
  "",
  paste0("Original pipeline reconstructs the established 24 genes exactly: ", BASELINE_EXACT_MATCH),
  paste0("SSL_HP significance column used: ", SSL_SIGNIFICANCE_SOURCE),
  "",
  "WHAT WAS TESTED",
  "---------------",
  "Every listed numerical threshold was varied one at a time. All selected combinations were also tested within the GTEx, HPA, ANTRUM and SSL_HP stages. Alternative HPA definitions were compared, including no separate colon ceiling and an all-tissue fourfold stomach-enrichment rule. A representative cross-stage grid then combined original, relaxed and stricter rules.",
  "",
  "RESULT BY ANALYSIS FAMILY",
  "-------------------------",
  family_lines,
  "",
  "ONE-FACTOR PARAMETERS ORDERED FROM MOST TO LEAST INFLUENTIAL",
  "------------------------------------------------------------",
  most_influential_parameter_lines,
  "",
  "TEN MOST THRESHOLD-SENSITIVE PRIMARY GENES",
  "------------------------------------------",
  most_sensitive_gene_lines,
  "",
  "HOW TO INTERPRET THIS",
  "---------------------",
  "A gene should not be called biologically invalid simply because it is removed by a stricter setting. The ranking identifies genes whose inclusion depends most strongly on the exact operational thresholds. Review their GTEx, HPA, ANTRUM and SSL_HP evidence in 15_PRIMARY_24_BASELINE_FILTER_MARGINS.csv and 16_PRIMARY_24_GENE_ROBUSTNESS_RANKING.csv.",
  "",
  "Retention of all 24 established genes does not by itself mean that the final list is unchanged. An alternative rule can retain all 24 and also admit additional genes. Exact reproduction requires n_primary_24_removed = 0, n_additional_genes_not_in_primary = 0 and exact_match_to_primary_24 = TRUE.",
  "",
  "The percentage retained across settings is descriptive. It depends on the range of values deliberately tested here and is not a statistical confidence level.",
  "",
  "KEY FILES",
  "---------",
  "01_BASELINE_SUMMARY.csv",
  "02_ONE_FACTOR_AT_A_TIME_SUMMARY.csv",
  "05_HPA_RULE_COMPARISON_SUMMARY.csv",
  "09_ALL_SETTINGS_SUMMARY.csv",
  "13_ONE_FACTOR_PARAMETER_IMPACT_RANKING.csv",
  "15_PRIMARY_24_BASELINE_FILTER_MARGINS.csv",
  "16_PRIMARY_24_GENE_ROBUSTNESS_RANKING.csv",
  "17_TEN_MOST_THRESHOLD_SENSITIVE_PRIMARY_GENES.csv",
  "18_SETTINGS_THAT_CHANGED_THE_24_GENE_LIST.csv",
  "20_HPA_RULE_FINAL_LIST_COMPOSITION_AND_GENE_CHANGES.csv",
  "21_HPA_RULE_ADDED_AND_REMOVED_GENES_LONG.csv",
  "",
  "No established signature or downstream thesis result was overwritten."
)

write_text_overwrite(
  interpretation_lines,
  file.path(params$out_dir, "00_READ_ME_FIRST_RESULT_INTERPRETATION.txt")
)

RUN_SUMMARY <- data.frame(
  metric = c(
    "baseline_exactly_matches_established_24",
    "established_signature_size",
    "baseline_reconstructed_signature_size",
    "number_one_factor_settings",
    "number_GTEx_grid_settings",
    "number_HPA_grid_settings",
    "number_HPA_rule_settings",
    "number_ANTRUM_grid_settings",
    "number_SSL_HP_grid_settings",
    "number_combined_representative_settings",
    "number_all_settings_including_baseline",
    "number_all_settings_reproducing_exact_24",
    "number_all_settings_changing_the_24_gene_list",
    "output_directory"
  ),
  value = as.character(c(
    BASELINE_EXACT_MATCH,
    length(PRIMARY_SIGNATURE),
    length(BASELINE_FINAL_GENES),
    nrow(OFAT_SETTINGS),
    nrow(GTEX_GRID_SETTINGS),
    nrow(HPA_GRID_SETTINGS),
    nrow(HPA_RULE_SETTINGS),
    nrow(ANTRUM_GRID_SETTINGS),
    nrow(SSL_GRID_SETTINGS),
    nrow(COMBINED_SETTINGS),
    nrow(ALL_SUMMARIES),
    sum(ALL_SUMMARIES$exact_match_to_primary_24),
    sum(!ALL_SUMMARIES$exact_match_to_primary_24),
    params$out_dir
  )),
  stringsAsFactors = FALSE
)

write_csv_overwrite(
  RUN_SUMMARY,
  file.path(params$out_dir, "qc", "RUN_SUMMARY.csv")
)

message("============================================================")
message("FULL THRESHOLD SENSITIVITY AUDIT COMPLETE")
message("============================================================")
message("Original reconstruction matches established 24 exactly: ", BASELINE_EXACT_MATCH)
message("Total settings evaluated, including baseline: ", nrow(ALL_SUMMARIES))
message("Settings reproducing the exact 24 genes: ", sum(ALL_SUMMARIES$exact_match_to_primary_24))
message("Settings changing the exact list: ", sum(!ALL_SUMMARIES$exact_match_to_primary_24))
message("Output folder: ", params$out_dir)
message("Read: 00_READ_ME_FIRST_RESULT_INTERPRETATION.txt")
message("Most important ranking: tables/16_PRIMARY_24_GENE_ROBUSTNESS_RANKING.csv")
