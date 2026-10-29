#!/usr/bin/env Rscript
# Preranked pathway GSEA from the final tissue differential-expression tables.
# Run from Tissue_Bulk_RNAseq after configuring config/config.R.
# Primary Hallmark ranks use sign(logFC) * sqrt(LR); all-gene and logCPM >= 2 results are exported.


source("config/config.R")

input_dir <- PATH_DE_TABLES

out_dir <- OUTPUT_GSEA

install_missing_packages <- FALSE

run_all_genes <- TRUE
run_logcpm_filtered <- TRUE
logcpm_cutoff <- 2

primary_filter_mode <- "LOGCPM_GE_2"

rank_metric_mode <- "signed_sqrt_LR"

gene_set_families_to_run <- c(
  "HALLMARK",
  "REACTOME",
  "KEGG",
  "BIOCARTA",
  "GO_BP",
  "GO_CC",
  "GO_MF"
)

fgsea_min_size <- 15
fgsea_max_size <- 500
fgsea_eps <- 0
fgsea_score_type <- "std"
fgsea_nproc <- 0
fgsea_nperm_simple <- 10000
fdr_alpha <- 0.05
set_seed <- 1

make_plots <- TRUE
make_top_hallmark_enrichment_plots <- TRUE
top_hallmark_enrichment_plots_per_comparison <- 3

use_fixed_final_file_list <- TRUE

cran_packages <- c(
  "data.table", "dplyr", "tidyr", "stringr", "tibble",
  "ggplot2", "readr", "purrr", "msigdbr"
)

bioc_packages <- c("fgsea")

install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!install_missing_packages) {
      stop(
        "Package '", pkg, "' is not installed.\n",
        "Install the required package before rerunning."
      )
    }

    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
}

for (pkg in cran_packages) install_if_missing(pkg, bioc = FALSE)
for (pkg in bioc_packages) install_if_missing(pkg, bioc = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(msigdbr)
  library(fgsea)
})

if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir)
}

ensure_new_output(out_dir)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tables_dir <- file.path(out_dir, "tables")
fgsea_dir <- file.path(tables_dir, "fgsea_by_comparison")
ranked_dir <- file.path(out_dir, "ranked_gene_tables")
plot_dir <- file.path(out_dir, "plots")
enrichment_plot_dir <- file.path(out_dir, "enrichment_plots")
metadata_dir <- file.path(out_dir, "metadata")

for (d in c(tables_dir, fgsea_dir, ranked_dir, plot_dir, enrichment_plot_dir, metadata_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

safe_filename <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

clean_gene_symbols <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x[x == ""] <- NA_character_
  x
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

detect_col <- function(df, candidates, required = FALSE, label = "column") {
  cn <- colnames(df)

  exact <- candidates[candidates %in% cn]
  if (length(exact) > 0) return(exact[1])

  lower_cn <- tolower(cn)
  lower_candidates <- tolower(candidates)
  idx <- match(lower_candidates, lower_cn)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0) return(cn[idx[1]])

  if (required) {
    stop(
      "Could not detect ", label, ". Available columns are:\n",
      paste(cn, collapse = "\n")
    )
  }

  NA_character_
}

detect_gene_col <- function(df) {
  detect_col(
    df,
    candidates = c("SYMBOL", "Symbol", "symbol", "gene_symbol", "Gene", "gene", "genes"),
    required = TRUE,
    label = "gene symbol column"
  )
}

detect_logfc_col <- function(df) {
  detect_col(
    df,
    candidates = c("logFC", "LogFC", "logfc", "log2FoldChange", "log2FC", "coef", "estimate"),
    required = TRUE,
    label = "logFC/log2 fold-change column"
  )
}

detect_logcpm_col <- function(df) {
  detect_col(
    df,
    candidates = c(
      "logCPM", "LogCPM", "logcpm",
      "log2CPM", "Log2CPM", "log2cpm", "log2_CPM", "log2_cpm",
      "mean_logCPM", "meanLogCPM", "mean_log2CPM", "meanLog2CPM",
      "AveExpr", "AveExpression", "average_logCPM", "averageLogCPM"
    ),
    required = FALSE,
    label = "logCPM/log2CPM/AveExpr column"
  )
}

detect_p_col <- function(df) {
  detect_col(
    df,
    candidates = c("PValue", "P.Value", "Pvalue", "pvalue", "p.value", "P", "p", "p_val", "pval"),
    required = TRUE,
    label = "nominal P-value column"
  )
}

detect_fdr_col <- function(df) {
  detect_col(
    df,
    candidates = c("FDR", "fdr", "padj", "adj.P.Val", "adj.PValue", "adj_p_value", "qvalue", "q_value", "p_val_adj"),
    required = FALSE,
    label = "FDR/adjusted P-value column"
  )
}

detect_lr_col <- function(df) {
  detect_col(
    df,
    candidates = c("LR", "LRT", "likelihood_ratio", "LikelihoodRatio"),
    required = FALSE,
    label = "edgeR likelihood-ratio statistic column"
  )
}

detect_signed_stat_col <- function(df) {
  detect_col(
    df,
    candidates = c("t", "T", "stat", "STAT", "Statistic", "WaldStatistic", "moderated_t", "tstat", "z", "Z", "score"),
    required = FALSE,
    label = "signed test statistic column"
  )
}

empty_fgsea_result <- function() {
  tibble::tibble(
    comparison = character(),
    filter_mode = character(),
    gene_set_family = character(),
    pathway = character(),
    pval = numeric(),
    padj = numeric(),
    log2err = numeric(),
    ES = numeric(),
    NES = numeric(),
    size = integer(),
    leadingEdge = character()
  )
}

if (use_fixed_final_file_list) {
  de_file_table <- tibble::tribble(
    ~comparison,                         ~file,

    "CA_v_H_COL",                        "CGG014_Maier_CA_v_H_COL_DE_results.txt",
    "ES_H_ADJ_SSL_HP",                   "CGG014_Maier_ES_H_ADJ_SSL_HP_DE_results.txt",
    "ES_SSL_HP_and_H_ADJ_SSL_HP",        "CGG014_Maier_ES_SSL_HP_and_H_ADJ_SSL_HP_DE_results.txt",
    "ES_SSL_HP",                         "CGG014_Maier_ES_SSL_HP_DE_results.txt",

    "Fuso_H_ADJ_SSL_HP",                 "CGG014_Maier_Fuso_H_ADJ_SSL_HP_DE_results.txt",
    "Fuso_SSL_HP_and_H_ADJ_SSL_HP",      "CGG014_Maier_Fuso_SSL_HP_and_H_ADJ_SSL_HP_DE_results.txt",
    "Fuso_SSL_HP",                       "CGG014_Maier_Fuso_SSL_HP_DE_results.txt",

    "H_ADJ_SSL_HP_v_H_COL",              "CGG014_Maier_H_ADJ_SSL_HP_v_H_COL_DE_results.txt",

    "Serr_H_ADJ_SSL_HP",                 "CGG014_Maier_Serr_H_ADJ_SSL_HP_DE_results.txt",
    "Serr_SSL_HP_and_H_ADJ_SSL_HP",      "CGG014_Maier_Serr_SSL_HP_and_H_ADJ_SSL_HP_DE_results.txt",
    "Serr_SSL_HP",                       "CGG014_Maier_Serr_SSL_HP_DE_results.txt",

    "SSL_HP_v_CA",                       "CGG014_Maier_SSL_HP_v_CA_DE_results.txt",
    "SSL_HP_v_H_ADJ_SSL_HP",             "CGG014_Maier_SSL_HP_v_H_ADJ_SSL_HP_DE_results.txt",
    "SSL_HP_v_H_COL",                    "CGG014_Maier_SSL_HP_v_H_COL_DE_results.txt",

    "V600E_v_WT_in_SSL_HP",              "CGG014_Maier_V600E_v_WT_in_SSL_HP_DE_results.txt",

    "Alt_Infected",                      "CGG028_Alt_Infected_DE_results.txt",
    "Initial_BRAF_SSL_HP",               "CGG028_Initial_BRAF_SSL_HP_DE_results.txt"
  ) %>%
    dplyr::mutate(
      path = file.path(input_dir, file),
      comparison = safe_filename(comparison)
    )
} else {
  de_files <- list.files(input_dir, pattern = "_DE_results\\.txt$", full.names = TRUE)

  if (length(de_files) == 0) {
    stop("No *_DE_results.txt files found in input_dir: ", input_dir)
  }

  de_file_table <- tibble::tibble(
    path = de_files,
    file = basename(de_files),
    comparison = file %>%
      stringr::str_replace("_DE_results\\.txt$", "") %>%
      stringr::str_replace("^CGG014_Maier_", "") %>%
      stringr::str_replace("^CGG028_", "") %>%
      safe_filename()
  )
}

missing_de_files <- de_file_table %>% dplyr::filter(!file.exists(path))

if (nrow(missing_de_files) > 0) {
  stop(
    "These DE files are missing from input_dir:\n",
    paste(missing_de_files$file, collapse = "\n"),
    "\n\nInput directory checked:\n", input_dir
  )
}

if (any(duplicated(de_file_table$comparison))) {
  stop("Duplicate comparison names found. Make each comparison name unique.")
}

readr::write_csv(de_file_table, file.path(metadata_dir, "DE_FILES_ANALYSED.csv"))

message("Found ", nrow(de_file_table), " DE result files:")
for (i in seq_len(nrow(de_file_table))) {
  message("  ", de_file_table$comparison[i], " -> ", de_file_table$path[i])
}

get_msig_collection <- function(collection_code) {
  out <- tryCatch(
    suppressWarnings(msigdbr::msigdbr(species = "Homo sapiens", collection = collection_code)),
    error = function(e) NULL
  )

  if (is.null(out) || nrow(out) == 0) {
    out <- tryCatch(
      suppressWarnings(msigdbr::msigdbr(species = "Homo sapiens", category = collection_code)),
      error = function(e) NULL
    )
  }

  if (is.null(out) || nrow(out) == 0) {
    stop("Could not load MSigDB collection: ", collection_code)
  }

  as.data.frame(out)
}

standardise_msig <- function(msig_df, collection_label) {
  gene_col <- detect_col(
    msig_df,
    candidates = c("gene_symbol", "human_gene_symbol", "db_gene_symbol", "symbol"),
    required = TRUE,
    label = "MSigDB gene symbol column"
  )

  pathway_col <- detect_col(
    msig_df,
    candidates = c("gs_name", "geneset", "gene_set_name", "gs_exact_source"),
    required = TRUE,
    label = "MSigDB pathway/gene-set name column"
  )

  subcollection_col <- detect_col(
    msig_df,
    candidates = c("gs_subcollection", "gs_subcat", "subcategory"),
    required = FALSE,
    label = "MSigDB subcollection column"
  )

  version_col <- detect_col(
    msig_df,
    candidates = c("db_version", "msigdb_version"),
    required = FALSE,
    label = "MSigDB database version column"
  )

  tibble::tibble(
    msigdb_collection = collection_label,
    msigdb_subcollection = if (!is.na(subcollection_col)) as.character(msig_df[[subcollection_col]]) else NA_character_,
    pathway = as.character(msig_df[[pathway_col]]),
    gene = clean_gene_symbols(msig_df[[gene_col]]),
    msigdb_version = if (!is.na(version_col)) as.character(msig_df[[version_col]]) else NA_character_
  ) %>%
    dplyr::filter(!is.na(pathway), !is.na(gene), gene != "") %>%
    dplyr::distinct()
}

pathway_list_from_df <- function(df) {
  split(df$gene, df$pathway) %>%
    lapply(function(x) sort(unique(clean_gene_symbols(x))))
}

message("\nLoading MSigDB collections with msigdbr...")

msig_h <- standardise_msig(get_msig_collection("H"), "H")
msig_c2 <- standardise_msig(get_msig_collection("C2"), "C2")
msig_c5 <- standardise_msig(get_msig_collection("C5"), "C5")

msig_all_metadata <- dplyr::bind_rows(msig_h, msig_c2, msig_c5) %>%
  dplyr::distinct()

readr::write_csv(msig_all_metadata, file.path(metadata_dir, "MSIGDB_GENE_SET_METADATA_USED_LONG.csv"))

msigdb_versions <- msig_all_metadata %>%
  dplyr::filter(!is.na(msigdb_version), msigdb_version != "") %>%
  dplyr::distinct(msigdb_version) %>%
  dplyr::arrange(msigdb_version)

readr::write_csv(msigdb_versions, file.path(metadata_dir, "MSIGDB_VERSION_USED.csv"))

c2_annot <- msig_c2 %>%
  dplyr::mutate(
    pathway_upper = toupper(pathway),
    subcollection_upper = toupper(dplyr::coalesce(msigdb_subcollection, ""))
  )

c5_annot <- msig_c5 %>%
  dplyr::mutate(
    pathway_upper = toupper(pathway),
    subcollection_upper = toupper(dplyr::coalesce(msigdb_subcollection, ""))
  )

hallmark_df <- msig_h

reactome_df <- c2_annot %>%
  dplyr::filter(stringr::str_detect(pathway_upper, "^REACTOME_") | stringr::str_detect(subcollection_upper, "REACTOME")) %>%
  dplyr::select(-pathway_upper, -subcollection_upper)

kegg_df <- c2_annot %>%
  dplyr::filter(stringr::str_detect(pathway_upper, "^KEGG_") | stringr::str_detect(subcollection_upper, "KEGG")) %>%
  dplyr::select(-pathway_upper, -subcollection_upper)

biocarta_df <- c2_annot %>%
  dplyr::filter(stringr::str_detect(pathway_upper, "^BIOCARTA_") | stringr::str_detect(subcollection_upper, "BIOCARTA")) %>%
  dplyr::select(-pathway_upper, -subcollection_upper)

go_bp_df <- c5_annot %>%
  dplyr::filter(
    stringr::str_detect(pathway_upper, "^GOBP_") |
      stringr::str_detect(subcollection_upper, "GO:BP|GOBP|BIOLOGICAL_PROCESS")
  ) %>%
  dplyr::select(-pathway_upper, -subcollection_upper)

go_cc_df <- c5_annot %>%
  dplyr::filter(
    stringr::str_detect(pathway_upper, "^GOCC_") |
      stringr::str_detect(subcollection_upper, "GO:CC|GOCC|CELLULAR_COMPONENT")
  ) %>%
  dplyr::select(-pathway_upper, -subcollection_upper)

go_mf_df <- c5_annot %>%
  dplyr::filter(
    stringr::str_detect(pathway_upper, "^GOMF_") |
      stringr::str_detect(subcollection_upper, "GO:MF|GOMF|MOLECULAR_FUNCTION")
  ) %>%
  dplyr::select(-pathway_upper, -subcollection_upper)

gene_set_families <- list(
  HALLMARK = pathway_list_from_df(hallmark_df),
  REACTOME = pathway_list_from_df(reactome_df),
  KEGG = pathway_list_from_df(kegg_df),
  BIOCARTA = pathway_list_from_df(biocarta_df),
  GO_BP = pathway_list_from_df(go_bp_df),
  GO_CC = pathway_list_from_df(go_cc_df),
  GO_MF = pathway_list_from_df(go_mf_df)
)

gene_set_families <- gene_set_families[names(gene_set_families) %in% gene_set_families_to_run]

gene_set_family_summary <- purrr::map_dfr(names(gene_set_families), function(fam) {
  pathways <- gene_set_families[[fam]]

  tibble::tibble(
    gene_set_family = fam,
    n_pathways = length(pathways),
    min_database_size = safe_min(lengths(pathways)),
    median_database_size = safe_median(lengths(pathways)),
    max_database_size = safe_max(lengths(pathways))
  )
})

readr::write_csv(gene_set_family_summary, file.path(metadata_dir, "GENE_SET_FAMILY_SUMMARY.csv"))

message("\nGene-set families loaded:")
print(gene_set_family_summary)

make_ranked_object <- function(df, comparison, filter_mode) {
  gene_col <- detect_gene_col(df)
  logfc_col <- detect_logfc_col(df)
  logcpm_col <- detect_logcpm_col(df)
  p_col <- detect_p_col(df)
  fdr_col <- detect_fdr_col(df)
  lr_col <- detect_lr_col(df)
  signed_stat_col <- detect_signed_stat_col(df)

  if (filter_mode == "LOGCPM_GE_2" && is.na(logcpm_col)) {
    stop("LOGCPM_GE_2 requested but no logCPM/log2CPM/AveExpr column was detected for ", comparison)
  }

  dat <- df %>%
    as.data.frame() %>%
    dplyr::mutate(
      Gene = clean_gene_symbols(.data[[gene_col]]),
      logFC = suppressWarnings(as.numeric(.data[[logfc_col]])),
      logCPM = if (!is.na(logcpm_col)) suppressWarnings(as.numeric(.data[[logcpm_col]])) else NA_real_,
      PValue = suppressWarnings(as.numeric(.data[[p_col]])),
      FDR = if (!is.na(fdr_col)) suppressWarnings(as.numeric(.data[[fdr_col]])) else NA_real_,
      LR = if (!is.na(lr_col)) suppressWarnings(as.numeric(.data[[lr_col]])) else NA_real_,
      signed_stat = if (!is.na(signed_stat_col)) suppressWarnings(as.numeric(.data[[signed_stat_col]])) else NA_real_
    )

  n_input_rows <- nrow(dat)

  dat <- dat %>%
    dplyr::filter(!is.na(Gene), Gene != "") %>%
    dplyr::filter(is.finite(logFC), is.finite(PValue))

  n_valid_rank_rows_before_expression_filter <- nrow(dat)

  if (filter_mode == "LOGCPM_GE_2") {
    dat <- dat %>%
      dplyr::filter(is.finite(logCPM), logCPM >= logcpm_cutoff)
  }

  n_after_expression_filter <- nrow(dat)

  if (n_after_expression_filter < 1000) {
    warning(
      "Only ", n_after_expression_filter,
      " rows retained before duplicate collapse for ", comparison,
      " / ", filter_mode, ". Check filtering and gene IDs."
    )
  }

  min_nonzero_p <- suppressWarnings(min(dat$PValue[dat$PValue > 0], na.rm = TRUE))
  if (!is.finite(min_nonzero_p)) min_nonzero_p <- 1e-300

  dat <- dat %>%
    dplyr::mutate(
      PValue_safe = dplyr::if_else(PValue <= 0 | is.na(PValue), min_nonzero_p * 0.1, PValue)
    )

  if (rank_metric_mode == "signed_sqrt_LR" && !is.na(lr_col) && sum(is.finite(dat$LR)) >= 100) {
    dat <- dat %>%
      dplyr::mutate(rank_metric = sign(logFC) * sqrt(pmax(LR, 0)))

    rank_metric_used <- paste0("sign(logFC) * sqrt(", lr_col, ")")
  } else if (rank_metric_mode == "signed_stat" && !is.na(signed_stat_col) && sum(is.finite(dat$signed_stat)) >= 100) {
    dat <- dat %>%
      dplyr::mutate(rank_metric = signed_stat)

    rank_metric_used <- paste0("signed statistic column: ", signed_stat_col)
  } else {
    dat <- dat %>%
      dplyr::mutate(rank_metric = sign(logFC) * -log10(PValue_safe))

    rank_metric_used <- "sign(logFC) * -log10(PValue)"
  }

  dat <- dat %>%
    dplyr::filter(is.finite(rank_metric), !is.na(rank_metric))

  duplicate_gene_summary <- dat %>%
    dplyr::count(Gene, name = "n_rows_for_symbol") %>%
    dplyr::summarise(
      n_symbols_before_collapse = dplyr::n(),
      n_symbols_with_duplicates = sum(n_rows_for_symbol > 1),
      max_rows_per_symbol = max(n_rows_for_symbol, na.rm = TRUE),
      .groups = "drop"
    )

  dat_unique <- dat %>%
    dplyr::group_by(Gene) %>%
    dplyr::arrange(dplyr::desc(abs(rank_metric)), PValue_safe, dplyr::desc(abs(logFC)), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(rank_metric), Gene)

  ranks <- dat_unique$rank_metric
  names(ranks) <- dat_unique$Gene

  if (any(duplicated(ranks))) {
    ranks <- ranks - seq_along(ranks) * 1e-12
    dat_unique$rank_metric <- ranks
  }

  ranks <- sort(ranks, decreasing = TRUE)

  metadata <- tibble::tibble(
    comparison = comparison,
    filter_mode = filter_mode,
    n_input_rows = n_input_rows,
    n_valid_rank_rows_before_expression_filter = n_valid_rank_rows_before_expression_filter,
    n_after_expression_filter = n_after_expression_filter,
    n_ranked_genes_after_symbol_collapse = length(ranks),
    rank_metric_mode = rank_metric_mode,
    rank_metric_used = rank_metric_used,
    gene_col = gene_col,
    logfc_col = logfc_col,
    logcpm_col = logcpm_col,
    p_col = p_col,
    fdr_col = fdr_col,
    lr_col = lr_col,
    signed_stat_col = signed_stat_col,
    n_symbols_with_duplicates = duplicate_gene_summary$n_symbols_with_duplicates,
    max_rows_per_symbol = duplicate_gene_summary$max_rows_per_symbol
  )

  list(
    ranks = ranks,
    cleaned_table = dat_unique,
    metadata = metadata
  )
}

make_pathway_overlap_table <- function(pathways, ranks, comparison, filter_mode, gene_set_family) {
  if (length(pathways) == 0) {
    return(tibble::tibble())
  }

  tibble::tibble(
    comparison = comparison,
    filter_mode = filter_mode,
    gene_set_family = gene_set_family,
    pathway = names(pathways),
    n_genes_database = lengths(pathways),
    n_genes_overlap = vapply(pathways, function(g) sum(g %in% names(ranks)), integer(1)),
    overlap_fraction = n_genes_overlap / n_genes_database,
    passes_size_filter = n_genes_overlap >= fgsea_min_size & n_genes_overlap <= fgsea_max_size
  )
}

prepare_pathways_for_ranks <- function(pathways, ranks) {
  pathways <- pathways[!is.na(names(pathways))]
  pathways <- lapply(pathways, function(x) sort(unique(clean_gene_symbols(x))))
  pathways <- lapply(pathways, function(x) x[!is.na(x) & x != ""])

  overlap_sizes <- vapply(pathways, function(g) sum(g %in% names(ranks)), integer(1))
  keep <- overlap_sizes >= fgsea_min_size & overlap_sizes <= fgsea_max_size

  pathways[keep]
}

run_fgsea_family <- function(pathways, ranks, comparison, filter_mode, gene_set_family) {
  pathways_prepared <- prepare_pathways_for_ranks(pathways, ranks)

  if (length(pathways_prepared) == 0) {
    warning(
      "No pathways passed overlap filters for ", comparison,
      " / ", filter_mode,
      " / ", gene_set_family
    )

    return(empty_fgsea_result())
  }

  set.seed(set_seed)

  fg <- fgsea::fgseaMultilevel(
    pathways = pathways_prepared,
    stats = ranks,
    minSize = fgsea_min_size,
    maxSize = fgsea_max_size,
    eps = fgsea_eps,
    scoreType = fgsea_score_type,
    nproc = fgsea_nproc,
    nPermSimple = fgsea_nperm_simple
  )

  fg <- as.data.frame(fg)

  if (!"padj" %in% colnames(fg)) {
    fg$padj <- p.adjust(fg$pval, method = "BH")
  }

  if (!"log2err" %in% colnames(fg)) {
    fg$log2err <- NA_real_
  }

  if (!"leadingEdge" %in% colnames(fg)) {
    fg$leadingEdge <- vector("list", nrow(fg))
  }

  fg %>%
    dplyr::arrange(padj, dplyr::desc(abs(NES))) %>%
    dplyr::mutate(
      comparison = comparison,
      filter_mode = filter_mode,
      gene_set_family = gene_set_family,
      leadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
    ) %>%
    dplyr::select(
      comparison, filter_mode, gene_set_family,
      pathway, pval, padj, log2err, ES, NES, size, leadingEdge
    )
}

plot_hallmark_nes_heatmap <- function(hallmark_results, filter_mode, outfile_prefix) {
  plot_df <- hallmark_results %>%
    dplyr::filter(.data$filter_mode == .env$filter_mode) %>%
    dplyr::mutate(
      comparison = factor(comparison, levels = unique(comparison)),
      pathway = factor(pathway, levels = rev(sort(unique(pathway)))),
      sig_label = dplyr::case_when(
        !is.na(padj) & padj < 0.001 ~ "***",
        !is.na(padj) & padj < 0.01 ~ "**",
        !is.na(padj) & padj < 0.05 ~ "*",
        !is.na(pval) & pval < 0.05 ~ ".",
        TRUE ~ ""
      )
    )

  if (nrow(plot_df) == 0) {
    message("No Hallmark data to plot for ", filter_mode)
    return(NULL)
  }

  p <- ggplot(plot_df, aes(x = comparison, y = pathway, fill = NES)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = sig_label), size = 2.8) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "NES"
    ) +
    labs(
      title = paste0("Hallmark preranked FGSEA NES | ", filter_mode),
      subtitle = "* FDR < 0.05; . nominal P < 0.05. FDR is within Hallmark family.",
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 6.8),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      panel.grid = element_blank()
    )

  ggsave(paste0(outfile_prefix, ".png"), p, width = 11, height = 10, dpi = 300)
  ggsave(paste0(outfile_prefix, ".pdf"), p, width = 11, height = 10)

  p
}

plot_top_pathway_dotplot <- function(results, filter_mode, gene_set_family, outfile_prefix, top_n_per_comparison = 10) {
  plot_df <- results %>%
    dplyr::filter(.data$filter_mode == .env$filter_mode, .data$gene_set_family == .env$gene_set_family) %>%
    dplyr::filter(is.finite(NES), is.finite(padj)) %>%
    dplyr::group_by(comparison) %>%
    dplyr::arrange(padj, dplyr::desc(abs(NES)), .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n_per_comparison) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      comparison = factor(comparison, levels = unique(comparison)),
      pathway = stringr::str_replace(pathway, paste0("^", gene_set_family, "_"), ""),
      pathway = factor(pathway, levels = rev(unique(pathway))),
      neg_log10_padj = -log10(pmax(padj, 1e-300))
    )

  if (nrow(plot_df) == 0) {
    message("No top pathway dotplot data for ", gene_set_family, " / ", filter_mode)
    return(NULL)
  }

  p <- ggplot(plot_df, aes(x = comparison, y = pathway, size = neg_log10_padj, fill = NES)) +
    geom_point(shape = 21, color = "black", alpha = 0.85) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "NES"
    ) +
    labs(
      title = paste0(gene_set_family, " top pathways | ", filter_mode),
      subtitle = paste0("Top ", top_n_per_comparison, " per comparison by FDR"),
      x = NULL,
      y = NULL,
      size = "-log10(FDR)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      panel.grid.major.y = element_blank()
    )

  ggsave(paste0(outfile_prefix, ".png"), p, width = 11, height = 8, dpi = 300)
  ggsave(paste0(outfile_prefix, ".pdf"), p, width = 11, height = 8)

  p
}

plot_enrichment_safely <- function(pathway_genes, ranks, comparison, filter_mode, pathway_name, outfile_prefix) {
  overlap <- intersect(pathway_genes, names(ranks))

  if (length(overlap) < fgsea_min_size) {
    message(
      "Skipping enrichment plot for ", comparison,
      " / ", filter_mode,
      " / ", pathway_name,
      ": only ", length(overlap), " overlapping genes."
    )

    return(NULL)
  }

  p <- fgsea::plotEnrichment(pathway_genes, ranks) +
    labs(
      title = paste0(pathway_name, "\n", comparison, " | ", filter_mode),
      x = "Ranked genes",
      y = "Running enrichment score"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave(paste0(outfile_prefix, ".png"), p, width = 7.5, height = 5.5, dpi = 300)
  ggsave(paste0(outfile_prefix, ".pdf"), p, width = 7.5, height = 5.5)

  p
}

all_fgsea_results <- list()
all_overlap_tables <- list()
all_run_metadata <- list()
rank_objects_for_plots <- list()

for (i in seq_len(nrow(de_file_table))) {
  comparison <- de_file_table$comparison[i]
  file_path <- de_file_table$path[i]

  message("\n\n############################################################")
  message("Reading: ", file_path)
  message("Comparison: ", comparison)
  message("############################################################")

  de <- data.table::fread(file_path, data.table = FALSE)

  filter_modes <- character()

  if (run_all_genes) {
    filter_modes <- c(filter_modes, "ALL_GENES")
  }

  if (run_logcpm_filtered) {
    filter_modes <- c(filter_modes, "LOGCPM_GE_2")
  }

  for (filter_mode in filter_modes) {
    message("\n--- Filter mode: ", filter_mode, " ---")

    rank_obj <- make_ranked_object(
      df = de,
      comparison = comparison,
      filter_mode = filter_mode
    )

    ranks <- rank_obj$ranks
    cleaned_table <- rank_obj$cleaned_table

    rank_objects_for_plots[[paste(comparison, filter_mode, sep = "__")]] <- list(
      comparison = comparison,
      filter_mode = filter_mode,
      ranks = ranks
    )

    ranked_outfile <- file.path(
      ranked_dir,
      paste0(safe_filename(comparison), "_", filter_mode, "_ranked_gene_table.csv")
    )

    cleaned_table %>%
      dplyr::arrange(dplyr::desc(rank_metric)) %>%
      readr::write_csv(ranked_outfile)

    all_run_metadata[[paste(comparison, filter_mode, sep = "__")]] <- rank_obj$metadata %>%
      dplyr::mutate(input_file = file_path, .after = "filter_mode")

    for (gene_set_family in names(gene_set_families)) {
      message("Running fgsea: ", comparison, " / ", filter_mode, " / ", gene_set_family)

      pathways <- gene_set_families[[gene_set_family]]

      overlap_one <- make_pathway_overlap_table(
        pathways = pathways,
        ranks = ranks,
        comparison = comparison,
        filter_mode = filter_mode,
        gene_set_family = gene_set_family
      )

      all_overlap_tables[[paste(comparison, filter_mode, gene_set_family, sep = "__")]] <- overlap_one

      fgsea_one <- run_fgsea_family(
        pathways = pathways,
        ranks = ranks,
        comparison = comparison,
        filter_mode = filter_mode,
        gene_set_family = gene_set_family
      )

      out_file_one <- file.path(
        fgsea_dir,
        paste0(safe_filename(comparison), "_", filter_mode, "_", gene_set_family, "_FGSEA.csv")
      )

      readr::write_csv(fgsea_one, out_file_one)

      all_fgsea_results[[paste(comparison, filter_mode, gene_set_family, sep = "__")]] <- fgsea_one
    }
  }
}

combined_fgsea <- dplyr::bind_rows(all_fgsea_results)
pathway_overlap_diagnostics <- dplyr::bind_rows(all_overlap_tables)
run_metadata <- dplyr::bind_rows(all_run_metadata)

readr::write_csv(combined_fgsea, file.path(tables_dir, "COMBINED_FGSEA_RESULTS_ALL_FAMILIES.csv"))
readr::write_csv(pathway_overlap_diagnostics, file.path(metadata_dir, "PATHWAY_OVERLAP_DIAGNOSTICS.csv"))
readr::write_csv(run_metadata, file.path(metadata_dir, "RUN_METADATA_COLUMNS_FILTERING_AND_RANKING.csv"))

for (fam in unique(combined_fgsea$gene_set_family)) {
  fam_df <- combined_fgsea %>%
    dplyr::filter(gene_set_family == fam) %>%
    dplyr::arrange(filter_mode, comparison, padj, dplyr::desc(abs(NES)))

  readr::write_csv(
    fam_df,
    file.path(tables_dir, paste0(fam, "_ALL_COMPARISONS_ALL_FILTERS.csv"))
  )

  fam_primary <- fam_df %>%
    dplyr::filter(filter_mode == primary_filter_mode)

  readr::write_csv(
    fam_primary,
    file.path(tables_dir, paste0(fam, "_PRIMARY_", primary_filter_mode, ".csv"))
  )
}

hallmark_primary <- combined_fgsea %>%
  dplyr::filter(gene_set_family == "HALLMARK", filter_mode == primary_filter_mode) %>%
  dplyr::arrange(comparison, padj, dplyr::desc(abs(NES)))

readr::write_csv(
  hallmark_primary,
  file.path(tables_dir, paste0("MAIN_HALLMARK_RESULTS_FOR_FIGURE_", primary_filter_mode, ".csv"))
)

significant_results <- combined_fgsea %>%
  dplyr::filter(!is.na(padj), padj < fdr_alpha) %>%
  dplyr::arrange(gene_set_family, filter_mode, comparison, padj, dplyr::desc(abs(NES)))

readr::write_csv(
  significant_results,
  file.path(tables_dir, paste0("SIGNIFICANT_FGSEA_RESULTS_FAMILY_SPECIFIC_FDR_LT_", fdr_alpha, ".csv"))
)

top10_each_comparison_family <- combined_fgsea %>%
  dplyr::filter(is.finite(padj), is.finite(NES)) %>%
  dplyr::group_by(gene_set_family, filter_mode, comparison) %>%
  dplyr::arrange(padj, dplyr::desc(abs(NES)), .by_group = TRUE) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

readr::write_csv(
  top10_each_comparison_family,
  file.path(tables_dir, "TOP10_PATHWAYS_PER_COMPARISON_FILTER_AND_FAMILY.csv")
)

summary_by_comparison_family <- combined_fgsea %>%
  dplyr::group_by(comparison, filter_mode, gene_set_family) %>%
  dplyr::summarise(
    n_pathways_tested = dplyr::n(),
    n_fdr_sig_positive = sum(padj < fdr_alpha & NES > 0, na.rm = TRUE),
    n_fdr_sig_negative = sum(padj < fdr_alpha & NES < 0, na.rm = TRUE),
    n_nominal_positive = sum(pval < 0.05 & NES > 0, na.rm = TRUE),
    n_nominal_negative = sum(pval < 0.05 & NES < 0, na.rm = TRUE),
    median_NES = safe_median(NES),
    min_padj = safe_min(padj),
    top_pathway_by_FDR = {
      idx <- which(is.finite(padj))
      if (length(idx) == 0) NA_character_ else pathway[idx[which.min(padj[idx])]]
    },
    top_pathway_NES = {
      idx <- which(is.finite(padj))
      if (length(idx) == 0) NA_real_ else NES[idx[which.min(padj[idx])]]
    },
    .groups = "drop"
  ) %>%
  dplyr::arrange(filter_mode, gene_set_family, comparison)

readr::write_csv(
  summary_by_comparison_family,
  file.path(tables_dir, "SUMMARY_BY_COMPARISON_FILTER_AND_GENE_SET_FAMILY.csv")
)

if (make_plots) {
  hallmark_results <- combined_fgsea %>%
    dplyr::filter(gene_set_family == "HALLMARK")

  if ("ALL_GENES" %in% hallmark_results$filter_mode) {
    plot_hallmark_nes_heatmap(
      hallmark_results = hallmark_results,
      filter_mode = "ALL_GENES",
      outfile_prefix = file.path(plot_dir, "HALLMARK_NES_HEATMAP_ALL_GENES")
    )
  }

  if ("LOGCPM_GE_2" %in% hallmark_results$filter_mode) {
    plot_hallmark_nes_heatmap(
      hallmark_results = hallmark_results,
      filter_mode = "LOGCPM_GE_2",
      outfile_prefix = file.path(plot_dir, "HALLMARK_NES_HEATMAP_LOGCPM_GE_2")
    )
  }

  for (fam in intersect(c("REACTOME", "KEGG", "BIOCARTA", "GO_BP", "GO_CC", "GO_MF"), unique(combined_fgsea$gene_set_family))) {
    if (primary_filter_mode %in% combined_fgsea$filter_mode) {
      plot_top_pathway_dotplot(
        results = combined_fgsea,
        filter_mode = primary_filter_mode,
        gene_set_family = fam,
        outfile_prefix = file.path(plot_dir, paste0(fam, "_TOP10_DOTPLOT_", primary_filter_mode)),
        top_n_per_comparison = 10
      )
    }
  }
}

if (make_top_hallmark_enrichment_plots && "HALLMARK" %in% names(gene_set_families)) {
  hallmark_top_for_plots <- combined_fgsea %>%
    dplyr::filter(gene_set_family == "HALLMARK", filter_mode == primary_filter_mode) %>%
    dplyr::filter(is.finite(padj), is.finite(NES)) %>%
    dplyr::group_by(comparison) %>%
    dplyr::arrange(padj, dplyr::desc(abs(NES)), .by_group = TRUE) %>%
    dplyr::slice_head(n = top_hallmark_enrichment_plots_per_comparison) %>%
    dplyr::ungroup()

  for (j in seq_len(nrow(hallmark_top_for_plots))) {
    comparison <- hallmark_top_for_plots$comparison[j]
    pathway <- hallmark_top_for_plots$pathway[j]
    key <- paste(comparison, primary_filter_mode, sep = "__")

    if (key %in% names(rank_objects_for_plots) && pathway %in% names(gene_set_families[["HALLMARK"]])) {
      out_prefix <- file.path(
        enrichment_plot_dir,
        paste0(safe_filename(comparison), "_", primary_filter_mode, "_", safe_filename(pathway))
      )

      plot_enrichment_safely(
        pathway_genes = gene_set_families[["HALLMARK"]][[pathway]],
        ranks = rank_objects_for_plots[[key]]$ranks,
        comparison = comparison,
        filter_mode = primary_filter_mode,
        pathway_name = pathway,
        outfile_prefix = out_prefix
      )
    }
  }
}

methods_text <- c(
  "Preranked gene-set enrichment analysis was performed with fgseaMultilevel using MSigDB gene sets obtained through msigdbr.",
  paste0("For each differential-expression contrast, analyses were performed using all valid genes and after filtering to genes with logCPM >= ", logcpm_cutoff, "."),
  "Gene identifiers were matched to MSigDB using HGNC gene symbols from the SYMBOL column.",
  "Duplicate gene symbols were collapsed by retaining the row with the largest absolute ranking statistic.",
  if (rank_metric_mode == "signed_sqrt_LR") {
    "For edgeR likelihood-ratio test tables, genes were ranked by sign(logFC) * sqrt(LR), giving a directional statistic from the non-directional LR statistic."
  } else {
    "Genes were ranked by sign(logFC) * -log10(nominal P value)."
  },
  paste0("Gene sets with fewer than ", fgsea_min_size, " or more than ", fgsea_max_size, " overlapping genes were excluded."),
  "Benjamini-Hochberg FDR correction was performed separately within each comparison, filter mode and gene-set family.",
  "Hallmark gene sets were treated as the primary compact pathway analysis; Reactome, KEGG, BioCarta and GO analyses were treated as secondary or exploratory analyses."
)

writeLines(methods_text, file.path(out_dir, "METHODS_TEXT_TEMPLATE.txt"))

package_versions <- tibble::tibble(
  package = c("R", "fgsea", "msigdbr", "data.table", "dplyr", "ggplot2"),
  version = c(
    paste0(R.version$major, ".", R.version$minor),
    as.character(utils::packageVersion("fgsea")),
    as.character(utils::packageVersion("msigdbr")),
    as.character(utils::packageVersion("data.table")),
    as.character(utils::packageVersion("dplyr")),
    as.character(utils::packageVersion("ggplot2"))
  )
)

readr::write_csv(package_versions, file.path(metadata_dir, "PACKAGE_VERSIONS.csv"))
writeLines(capture.output(sessionInfo()), file.path(metadata_dir, "SESSION_INFO.txt"))

message("\n\nDONE.")
message("Output directory:")
message(out_dir)

message("\nMain files to inspect first:")
message("1. tables/MAIN_HALLMARK_RESULTS_FOR_FIGURE_", primary_filter_mode, ".csv")
message("2. tables/HALLMARK_ALL_COMPARISONS_ALL_FILTERS.csv")
message("3. tables/COMBINED_FGSEA_RESULTS_ALL_FAMILIES.csv")
message("4. tables/SIGNIFICANT_FGSEA_RESULTS_FAMILY_SPECIFIC_FDR_LT_", fdr_alpha, ".csv")
message("5. tables/SUMMARY_BY_COMPARISON_FILTER_AND_GENE_SET_FAMILY.csv")
message("6. metadata/RUN_METADATA_COLUMNS_FILTERING_AND_RANKING.csv")
message("7. metadata/PATHWAY_OVERLAP_DIAGNOSTICS.csv")
message("8. plots/HALLMARK_NES_HEATMAP_", primary_filter_mode, ".png")
message("9. METHODS_TEXT_TEMPLATE.txt")

message("\nPrimary Hallmark result preview:")
print(hallmark_primary)
