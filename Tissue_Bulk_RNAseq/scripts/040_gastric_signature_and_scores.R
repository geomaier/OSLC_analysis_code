#!/usr/bin/env Rscript
# Derive and verify the 24-gene gastric signature, then export its heatmap and GSEA.
# Run from Tissue_Bulk_RNAseq with config/config.R configured.
# Sample-level scores and their statistical tests are implemented in script 090.

source("config/config.R")

options(stringsAsFactors = FALSE)
set.seed(1)

cran_packages <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "patchwork"
)

bioc_packages <- c(
  "AnnotationDbi", "org.Hs.eg.db", "fgsea"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0 || length(missing_bioc) > 0) {
  msg <- c("Missing R packages were detected.")
  if (length(missing_cran) > 0) {
    msg <- c(
      msg,
      paste0(
        "Install CRAN packages with:\ninstall.packages(c(",
        paste(sprintf("'%s'", missing_cran), collapse = ", "),
        "))"
      )
    )
  }
  if (length(missing_bioc) > 0) {
    msg <- c(
      msg,
      paste0(
        "Install Bioconductor packages with:\n",
        "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
        "BiocManager::install(c(",
        paste(sprintf("'%s'", missing_bioc), collapse = ", "),
        "))"
      )
    )
  }
  stop(paste(msg, collapse = "\n\n"), call. = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(fgsea)
  library(grid)
})

base_dir <- PATH_BULK

params <- list(

  gtex_file = file.path(
    base_dir, "GM LIST", "GTEx_v10_TPM_by_tissue.csv"
  ),
  hpa_file = file.path(
    base_dir, "GM LIST", "rna_tissue_hpa.tsv"
  ),
  combined_de_file = file.path(
    base_dir, "Individual DE Tables", "COMBINED_DE_TABLE_ALL_COMPARISONS.csv"
  ),
  de_files = c(
    SSL_HP_v_H_COL = file.path(
      base_dir, "Re-run GSEA", "CGG014_Maier_SSL_HP_v_H_COL_DE_results.txt"
    ),
    SSL_HP_v_H_ADJ_SSL_HP = file.path(
      base_dir, "Re-run GSEA", "CGG014_Maier_SSL_HP_v_H_ADJ_SSL_HP_DE_results.txt"
    ),
    SSL_HP_v_CA = file.path(
      base_dir, "Re-run GSEA", "CGG014_Maier_SSL_HP_v_CA_DE_results.txt"
    )
  ),

  out_dir = OUTPUT_GASTRIC,

  tau_min = 0.85,
  gtex_stomach_tpm_min = 2,
  gtex_colon_tpm_max = 1,
  gtex_stomach_colon_log2fc_min = 2.5,
  gtex_pseudocount = 0.1,

  hpa_stomach_ntpm_min = 2,
  hpa_colon_ntpm_max = 10,
  hpa_stomach_colon_log2fc_min = 3,
  hpa_pseudocount = 0.1,

  antrum_de_prefix = "ANTRUM_v_H_COL_DE_results",
  antrum_log2fc_min = 1.5,
  antrum_fdr_max = 0.05,

  final_signature_fdr_max = 0.05,
  final_signature_log2fc_min_exclusive = 0,
  expected_final_signature_size = 24L,
  strict_expected_signature_check = TRUE,
  expected_final_signature_genes = c(
    "AQP5", "ANXA10", "MUC5AC", "VSIG1", "TCN1", "MUC6",
    "PSAPL1", "CLDN18", "SFTA2", "MUCL3", "GATA4", "KLK11",
    "DUOXA1", "SULT1C2", "IRX3", "AADAC", "FOXQ1", "IGFALS",
    "NR0B2", "C13ORF46", "FOLR1", "LYPD6B", "TFR2", "SPTSSB"
  ),

  gsea_comparison_order = c(
    "CA_v_H_COL",
    "SSL_HP_v_H_COL",
    "H_ADJ_SSL_HP_v_H_COL"
  ),
  gsea_prefix_candidates = list(
    CA_v_H_COL = c("CA_v_H_COL_DE_results", "CA_v_H_COL"),
    SSL_HP_v_H_COL = c("SSL_HP_v_H_COL_DE_results", "SSL_HP_v_H_COL"),
    H_ADJ_SSL_HP_v_H_COL = c(
      "H_ADJ_SSL_HP_v_H_COL_DE_results",
      "H_ADJ_SSL_HP_v_H_COL"
    )
  ),
  gsea_comparison_labels = c(
    CA_v_H_COL = "CA vs H_COL",
    SSL_HP_v_H_COL = "SSL_HP vs H_COL",
    H_ADJ_SSL_HP_v_H_COL = "H_ADJ_SSL_HP vs H_COL"
  ),
  gsea_comparison_colours = c(
    CA_v_H_COL = "#825CA6",
    SSL_HP_v_H_COL = "#FF0000",
    H_ADJ_SSL_HP_v_H_COL = "#E9C61D"
  ),
  gsea_min_size = 1L,
  gsea_max_size = Inf,
  gsea_eps = 1e-10,
  gsea_nperm_simple = 1000L,
  gsea_score_type = "std",
  gsea_gsea_param = 1,

  require_all_final_signature_genes_in_gsea_rank = FALSE,

  gsea_plot_title = "Gastric metaplasia signature enrichment",
  gsea_panel_letter = "f",
  gsea_font_base_pt = 6,
  gsea_axis_title_pt = 6,
  gsea_axis_text_pt = 5,
  gsea_annotation_pt = 5,
  gsea_title_pt = 7,
  gsea_panel_label_pt = 8,
  gsea_curve_linewidth = 1.05,
  gsea_tick_linewidth = 0.35,

  font_family = "Arial",
  font_pt = 7,
  x_font_pt = 6,
  label_pt = 6,
  show_workflow_counts = FALSE,
  heatmap_colour_limit = NULL,
  figure_width_mm = 183,
  figure_dpi = 600
)

ensure_new_output(params$out_dir)

output_dirs <- c(
  params$out_dir,
  file.path(params$out_dir, "tables"),
  file.path(params$out_dir, "qc"),
  file.path(params$out_dir, "figures"),
  file.path(params$out_dir, "figures", "individual_panels"),
  file.path(params$out_dir, "figures", "thesis")
)

invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

if (.Platform$OS.type == "windows") {
  try(windowsFonts(Arial = windowsFont("Arial")), silent = TRUE)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

safe_max <- function(x) {
  x <- to_num(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

pick_first_existing <- function(candidates, x) {
  nms <- if (is.data.frame(x) || is.list(x)) names(x) else as.character(x)
  exact <- candidates[candidates %in% nms]
  if (length(exact) > 0) return(exact[1])

  idx <- match(tolower(candidates), tolower(nms))
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0) return(nms[idx[1]])

  NA_character_
}

assert_files_exist <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(
      "The following required input files were not found:\n",
      paste(missing, collapse = "\n"),
      call. = FALSE
    )
  }
}

read_table_flexible <- function(path) {
  data.table::fread(
    path,
    data.table = FALSE,
    fill = TRUE,
    check.names = FALSE
  )
}

write_csv_safe <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(df, path)
}

write_lines_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(as.character(x), con = path)
}

pt_to_mm <- function(x) x / 2.845276

save_plot_set <- function(plot, stem, width_mm, height_mm, dpi = params$figure_dpi) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

  ggsave(
    paste0(stem, ".pdf"), plot,
    width = width_mm, height = height_mm, units = "mm",
    device = pdf_device, bg = "white"
  )
  ggsave(
    paste0(stem, ".png"), plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = dpi, bg = "white"
  )
  ggsave(
    paste0(stem, ".jpg"), plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = dpi, device = grDevices::jpeg, quality = 100, bg = "white"
  )
  ggsave(
    paste0(stem, ".jpeg"), plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = dpi, device = grDevices::jpeg, quality = 100, bg = "white"
  )
  ggsave(
    paste0(stem, ".tiff"), plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = dpi, device = grDevices::tiff, compression = "lzw", bg = "white"
  )
}

prepare_gene_symbols <- function(ids) {
  ids <- trimws(as.character(ids))
  ids_novers <- sub("\\..*$", "", ids)
  ensembl_fraction <- mean(grepl("^ENSG", toupper(ids_novers)), na.rm = TRUE)
  ensembl_like <- is.finite(ensembl_fraction) && ensembl_fraction > 0.5

  if (ensembl_like) {
    mapped <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = unique(ids_novers),
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
    symbols <- unname(mapped[ids_novers])
    symbols <- toupper(trimws(symbols))

    missing_symbol <- is.na(symbols) | !nzchar(symbols)
    symbols[missing_symbol] <- toupper(ids_novers[missing_symbol])
  } else {
    symbols <- toupper(trimws(ids))
  }

  data.frame(
    gene_id = ids,
    gene_id_novers = ids_novers,
    gene_symbol = symbols,
    stringsAsFactors = FALSE
  )
}

standardise_gene_vector <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(character(0))
  ann <- prepare_gene_symbols(x)
  out <- unique(toupper(trimws(ann$gene_symbol)))
  out[!is.na(out) & nzchar(out)]
}

standardise_de_table <- function(path, require_fdr = FALSE) {
  df <- read_table_flexible(path)

  gene_col <- pick_first_existing(
    c(
      "SYMBOL", "symbol", "gene_symbol", "GeneSymbol", "Gene", "gene",
      "external_gene_name", "gene_name", "genes", "gene_id", "Geneid",
      "ENSEMBL", "ensembl_gene_id"
    ),
    df
  )
  if (is.na(gene_col)) gene_col <- names(df)[1]

  logfc_col <- pick_first_existing(
    c("logFC", "LogFC", "log2FC", "log2FoldChange", "LFC", "log2_fold_change"),
    df
  )
  fdr_col <- pick_first_existing(
    c("FDR", "fdr", "padj", "adj.P.Val", "qvalue", "q_value", "p_adj"),
    df
  )
  p_col <- pick_first_existing(
    c("PValue", "P.Value", "pvalue", "p.value", "pval", "p"),
    df
  )
  logcpm_col <- pick_first_existing(
    c("logCPM", "LogCPM", "log2CPM", "AveExpr", "average_logCPM"),
    df
  )

  if (is.na(logfc_col)) stop("No log2 fold-change column was found in: ", path)
  if (require_fdr && is.na(fdr_col)) {
    stop("No adjusted P-value/FDR column was found in: ", path)
  }

  gene_ann <- prepare_gene_symbols(df[[gene_col]])

  out <- data.frame(
    Gene = gene_ann$gene_symbol,
    log2FC = to_num(df[[logfc_col]]),
    FDR = if (!is.na(fdr_col)) to_num(df[[fdr_col]]) else NA_real_,
    PValue = if (!is.na(p_col)) to_num(df[[p_col]]) else NA_real_,
    logCPM = if (!is.na(logcpm_col)) to_num(df[[logcpm_col]]) else NA_real_,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(Gene), nzchar(Gene), is.finite(log2FC)) %>%
    dplyr::group_by(Gene) %>%
    dplyr::arrange(dplyr::desc(abs(log2FC)), FDR, PValue, .by_group = TRUE) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()

  attr(out, "columns_used") <- data.frame(
    file = path,
    gene_col = gene_col,
    logfc_col = logfc_col,
    fdr_col = fdr_col,
    p_col = p_col,
    logcpm_col = logcpm_col,
    stringsAsFactors = FALSE
  )

  out
}

extract_combined_de_contrast <- function(df, prefix) {
  symbol_col <- pick_first_existing(
    c(
      paste0(prefix, ".SYMBOL"),
      paste0(prefix, ".symbol"),
      paste0(prefix, ".Gene"),
      paste0(prefix, ".gene")
    ),
    df
  )
  logfc_col <- pick_first_existing(
    c(paste0(prefix, ".logFC"), paste0(prefix, ".log2FC")),
    df
  )
  fdr_col <- pick_first_existing(
    c(
      paste0(prefix, ".FDR"),
      paste0(prefix, ".padj"),
      paste0(prefix, ".adj.P.Val")
    ),
    df
  )

  if (any(is.na(c(symbol_col, logfc_col, fdr_col)))) {
    stop(
      "Could not identify SYMBOL, logFC and FDR columns for prefix: ", prefix
    )
  }

  out <- data.frame(
    Gene = toupper(trimws(as.character(df[[symbol_col]]))),
    log2FC = to_num(df[[logfc_col]]),
    FDR = to_num(df[[fdr_col]]),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(Gene), nzchar(Gene)) %>%
    dplyr::group_by(Gene) %>%
    dplyr::arrange(dplyr::desc(abs(log2FC)), FDR, .by_group = TRUE) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()

  attr(out, "columns_used") <- data.frame(
    prefix = prefix,
    symbol_col = symbol_col,
    logfc_col = logfc_col,
    fdr_col = fdr_col,
    stringsAsFactors = FALSE
  )

  out
}

tau_score <- function(v) {
  v <- to_num(v)
  n <- sum(is.finite(v))
  if (n < 2) return(NA_real_)
  vmax <- max(v, na.rm = TRUE)
  if (!is.finite(vmax) || vmax <= 0) return(NA_real_)
  sum(1 - v / vmax, na.rm = TRUE) / (n - 1)
}

generate_stomach_candidates <- function() {
  message("[1/9] Generating stomach-enriched candidate genes...")

  gtex <- fread(params$gtex_file, data.table = FALSE, check.names = FALSE)
  if (!"gene_symbol" %in% names(gtex)) {
    stop("GTEx file must contain a 'gene_symbol' column.")
  }

  tissue_cols <- base::setdiff(
    names(gtex), c("gene_symbol", "gene_id", "Name", "Description")
  )
  gtex[tissue_cols] <- lapply(gtex[tissue_cols], to_num)

  stomach_cols <- grep("^Stomach", names(gtex), value = TRUE)
  colon_cols <- grep("^Colon", names(gtex), value = TRUE)
  colon_cols <- colon_cols[!grepl("Muscularis", colon_cols, ignore.case = TRUE)]

  if (length(stomach_cols) == 0) stop("No GTEx stomach columns were detected.")
  if (length(colon_cols) == 0) stop("No GTEx colon columns remained after muscularis exclusion.")

  write_lines_safe(stomach_cols, file.path(params$out_dir, "qc", "GTEx_stomach_columns.txt"))
  write_lines_safe(colon_cols, file.path(params$out_dir, "qc", "GTEx_colon_columns.txt"))

  tissue_mat <- as.matrix(gtex[, tissue_cols, drop = FALSE])
  storage.mode(tissue_mat) <- "double"

  gtex$tau <- apply(tissue_mat, 1, tau_score)
  gtex$Stomach_TPM <- rowMeans(
    as.matrix(gtex[, stomach_cols, drop = FALSE]), na.rm = TRUE
  )
  gtex$Stomach_TPM[is.nan(gtex$Stomach_TPM)] <- NA_real_
  gtex$Colon_TPM <- apply(
    as.matrix(gtex[, colon_cols, drop = FALSE]), 1, safe_max
  )
  gtex$GTEx_stomach_colon_log2FC <-
    log2(gtex$Stomach_TPM + params$gtex_pseudocount) -
    log2(gtex$Colon_TPM + params$gtex_pseudocount)
  gtex$Gene <- toupper(trimws(as.character(gtex$gene_symbol)))

  gtex <- gtex %>%
    dplyr::mutate(
      pass_tau = is.finite(tau) & tau >= params$tau_min,
      pass_stomach = is.finite(Stomach_TPM) &
        Stomach_TPM >= params$gtex_stomach_tpm_min,
      pass_colon = is.finite(Colon_TPM) &
        Colon_TPM <= params$gtex_colon_tpm_max,
      pass_gtex_log2fc = is.finite(GTEx_stomach_colon_log2FC) &
        GTEx_stomach_colon_log2FC >= params$gtex_stomach_colon_log2fc_min,
      pass_all_gtex = pass_tau & pass_stomach & pass_colon & pass_gtex_log2fc
    )

  hpa <- fread(params$hpa_file, data.table = FALSE, check.names = FALSE)
  names(hpa) <- tolower(names(hpa))

  hpa_gene_col <- pick_first_existing(c("gene name", "gene"), hpa)
  hpa_tissue_col <- pick_first_existing(c("tissue"), hpa)
  hpa_ntpm_col <- pick_first_existing(c("ntpm"), hpa)

  if (any(is.na(c(hpa_gene_col, hpa_tissue_col, hpa_ntpm_col)))) {
    stop("HPA file must contain gene, tissue and nTPM columns.")
  }

  hpa$Gene <- toupper(trimws(as.character(hpa[[hpa_gene_col]])))
  hpa$tissue_lc <- tolower(trimws(as.character(hpa[[hpa_tissue_col]])))
  hpa$nTPM <- to_num(hpa[[hpa_ntpm_col]])

  hpa_stomach <- hpa %>%
    dplyr::filter(tissue_lc == "stomach") %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(HPA_stomach_nTPM = safe_max(nTPM), .groups = "drop")

  hpa_colon <- hpa %>%
    dplyr::filter(grepl("^colon", tissue_lc)) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(HPA_colon_nTPM = safe_max(nTPM), .groups = "drop")

  evidence <- gtex %>%
    dplyr::left_join(hpa_stomach, by = "Gene") %>%
    dplyr::left_join(hpa_colon, by = "Gene") %>%
    dplyr::mutate(
      HPA_colon_for_log2FC = ifelse(is.na(HPA_colon_nTPM), 0, HPA_colon_nTPM),
      HPA_stomach_colon_log2FC = log2(
        (HPA_stomach_nTPM + params$hpa_pseudocount) /
          (HPA_colon_for_log2FC + params$hpa_pseudocount)
      ),
      pass_hpa_stomach = is.finite(HPA_stomach_nTPM) &
        HPA_stomach_nTPM >= params$hpa_stomach_ntpm_min,
      pass_hpa_colon = is.na(HPA_colon_nTPM) |
        HPA_colon_nTPM <= params$hpa_colon_ntpm_max,
      pass_hpa_log2fc = is.finite(HPA_stomach_colon_log2FC) &
        HPA_stomach_colon_log2FC >= params$hpa_stomach_colon_log2fc_min,
      pass_all_hpa = pass_hpa_stomach & pass_hpa_colon & pass_hpa_log2fc,
      pass_external = pass_all_gtex & pass_all_hpa
    )

  external_candidates <- evidence %>%
    dplyr::filter(pass_external) %>%
    dplyr::arrange(dplyr::desc(GTEx_stomach_colon_log2FC), dplyr::desc(tau), dplyr::desc(HPA_stomach_colon_log2FC)) %>%
    dplyr::distinct(Gene, .keep_all = TRUE)

  combined_de <- fread(
    params$combined_de_file,
    data.table = FALSE,
    check.names = FALSE
  )
  antrum_de <- extract_combined_de_contrast(
    combined_de, params$antrum_de_prefix
  )

  evidence_full <- evidence %>%
    dplyr::left_join(
      antrum_de %>%
        dplyr::rename(
          ANTRUM_vs_H_COL_log2FC = log2FC,
          ANTRUM_vs_H_COL_FDR = FDR
        ),
      by = "Gene"
    ) %>%
    dplyr::mutate(
      pass_antrum = is.finite(ANTRUM_vs_H_COL_log2FC) &
        is.finite(ANTRUM_vs_H_COL_FDR) &
        ANTRUM_vs_H_COL_log2FC >= params$antrum_log2fc_min &
        ANTRUM_vs_H_COL_FDR <= params$antrum_fdr_max,
      pass_stomach_candidate = pass_external & pass_antrum
    )

  stomach_candidates <- evidence_full %>%
    dplyr::filter(pass_stomach_candidate) %>%
    dplyr::arrange(
      dplyr::desc(ANTRUM_vs_H_COL_log2FC),
      dplyr::desc(GTEx_stomach_colon_log2FC),
      dplyr::desc(tau)
    ) %>%
    dplyr::distinct(Gene, .keep_all = TRUE)

  if (nrow(stomach_candidates) == 0) {
    stop("No stomach-enriched genes passed all GTEx, HPA and antrum criteria.")
  }

  write_csv_safe(
    evidence_full,
    file.path(params$out_dir, "tables", "STOMACH_IDENTITY_all_gene_evidence.csv")
  )
  write_csv_safe(
    external_candidates,
    file.path(params$out_dir, "tables", "STOMACH_IDENTITY_external_candidates.csv")
  )
  write_csv_safe(
    stomach_candidates,
    file.path(params$out_dir, "tables", "STOMACH_IDENTITY_final_candidates.csv")
  )
  write_lines_safe(
    stomach_candidates$Gene,
    file.path(params$out_dir, "tables", "STOMACH_IDENTITY_final_candidate_genes.txt")
  )
  write_lines_safe(
    stomach_candidates$Gene,
    file.path(params$out_dir, "tables", "SIM_final_gene_list.txt")
  )

  write_csv_safe(
    attr(antrum_de, "columns_used"),
    file.path(params$out_dir, "qc", "ANTRUM_DE_columns_used.csv")
  )

  count_table <- data.frame(
    step = c(
      "GTEx rows with non-empty gene symbol",
      "GTEx genes passing tissue-specificity and stomach-versus-colon filters",
      "Genes also passing HPA validation",
      "Genes also passing gastric-antrum validation"
    ),
    n = c(
      sum(!is.na(evidence$Gene) & nzchar(evidence$Gene)),
      dplyr::n_distinct(evidence$Gene[evidence$pass_all_gtex %in% TRUE]),
      dplyr::n_distinct(evidence$Gene[evidence$pass_external %in% TRUE]),
      nrow(stomach_candidates)
    ),
    stringsAsFactors = FALSE
  )
  write_csv_safe(
    count_table,
    file.path(params$out_dir, "qc", "STOMACH_IDENTITY_workflow_counts.csv")
  )

  list(
    evidence = evidence_full,
    external_candidates = external_candidates,
    stomach_candidates = stomach_candidates,
    workflow_counts = count_table
  )
}

generate_final_signature <- function(stomach_candidates) {
  message("[2/9] Intersecting stomach-enriched candidates with SSL_HP-upregulated genes...")

  ssl_de <- standardise_de_table(
    params$de_files[["SSL_HP_v_H_COL"]],
    require_fdr = TRUE
  )

  if (!any(is.finite(ssl_de$FDR))) {
    stop("The SSL_HP versus H_COL DE table contains no finite adjusted P values.")
  }

  ssl_up <- ssl_de %>%
    dplyr::filter(
      is.finite(FDR),
      FDR <= params$final_signature_fdr_max,
      is.finite(log2FC),
      log2FC > params$final_signature_log2fc_min_exclusive
    ) %>%
    dplyr::pull(Gene) %>%
    unique()

  candidate_genes <- unique(stomach_candidates$Gene)
  final_genes <- base::intersect(candidate_genes, ssl_up)

  selection_audit <- stomach_candidates %>%
    dplyr::select(
      Gene,
      tau,
      Stomach_TPM,
      Colon_TPM,
      GTEx_stomach_colon_log2FC,
      HPA_stomach_nTPM,
      HPA_colon_nTPM,
      HPA_stomach_colon_log2FC,
      ANTRUM_vs_H_COL_log2FC,
      ANTRUM_vs_H_COL_FDR
    ) %>%
    dplyr::left_join(
      ssl_de %>%
        dplyr::select(
          Gene,
          SSL_HP_vs_H_COL_log2FC = log2FC,
          SSL_HP_vs_H_COL_FDR = FDR,
          SSL_HP_vs_H_COL_PValue = PValue,
          SSL_HP_vs_H_COL_logCPM = logCPM
        ),
      by = "Gene"
    ) %>%
    dplyr::mutate(
      passes_SSL_HP_vs_H_COL = Gene %in% ssl_up,
      selected_for_final_signature = Gene %in% final_genes,
      final_signature_logCPM_filter_applied = FALSE
    ) %>%
    dplyr::arrange(
      dplyr::desc(selected_for_final_signature),
      dplyr::desc(SSL_HP_vs_H_COL_log2FC),
      Gene
    )

  expected <- sort(unique(toupper(params$expected_final_signature_genes)))
  observed <- sort(unique(toupper(final_genes)))
  missing_expected <- base::setdiff(expected, observed)
  unexpected_observed <- base::setdiff(observed, expected)

  signature_check <- data.frame(
    metric = c(
      "expected_signature_size",
      "observed_signature_size",
      "number_missing_expected_genes",
      "number_unexpected_genes"
    ),
    value = c(
      params$expected_final_signature_size,
      length(observed),
      length(missing_expected),
      length(unexpected_observed)
    ),
    stringsAsFactors = FALSE
  )

  write_csv_safe(
    selection_audit,
    file.path(params$out_dir, "tables", "FINAL_METAPLASIA_SIGNATURE_selection_audit.csv")
  )
  write_csv_safe(
    signature_check,
    file.path(params$out_dir, "qc", "FINAL_METAPLASIA_SIGNATURE_expected_check.csv")
  )
  write_lines_safe(
    missing_expected,
    file.path(params$out_dir, "qc", "FINAL_METAPLASIA_SIGNATURE_missing_expected_genes.txt")
  )
  write_lines_safe(
    unexpected_observed,
    file.path(params$out_dir, "qc", "FINAL_METAPLASIA_SIGNATURE_unexpected_genes.txt")
  )
  write_csv_safe(
    attr(ssl_de, "columns_used"),
    file.path(params$out_dir, "qc", "SSL_HP_vs_H_COL_DE_columns_used.csv")
  )

  if (params$strict_expected_signature_check) {
    if (
      length(observed) != params$expected_final_signature_size ||
        length(missing_expected) > 0 ||
        length(unexpected_observed) > 0
    ) {
      stop(
        "The derived final signature does not match the expected 24-gene set. ",
        "Inspect the QC files in: ", file.path(params$out_dir, "qc")
      )
    }
  }

  final_ordered <- ssl_de %>%
    dplyr::filter(Gene %in% final_genes) %>%
    dplyr::arrange(dplyr::desc(log2FC), FDR, Gene) %>%
    dplyr::pull(Gene)

  write_lines_safe(
    ssl_up,
    file.path(params$out_dir, "tables", "SSL_HP_vs_H_COL_upregulated_FDR_le_0p05_log2FC_gt_0.txt")
  )
  write_lines_safe(
    final_ordered,
    file.path(params$out_dir, "tables", "FINAL_24_GENE_GASTRIC_METAPLASIA_SIGNATURE.txt")
  )
  write_lines_safe(
    final_ordered,
    file.path(params$out_dir, "tables", "FINAL_SSL_METAPLASIA_LIST.txt")
  )
  write_csv_safe(
    data.frame(
      order = seq_along(final_ordered),
      gene_symbol = final_ordered,
      stringsAsFactors = FALSE
    ),
    file.path(params$out_dir, "tables", "FINAL_24_GENE_GASTRIC_METAPLASIA_SIGNATURE.csv")
  )

  list(
    final_genes = final_ordered,
    ssl_up = ssl_up,
    ssl_de = ssl_de,
    selection_audit = selection_audit
  )
}

make_workflow_panel <- function(workflow_counts, final_n) {
  labels <- c(
    "GTEx v10\nnormal-tissue\nexpression",
    "Tissue-\nspecificity\nscoring",
    "Stomach-versus-\ncolon filtering",
    "HPA RNA\nvalidation",
    "Gastric antrum\nvalidation",
    "Stomach-enriched\ncandidate genes",
    paste0(
      "Intersect with genes\nupregulated in SSL_HP\nversus H_COL\n",
      "FDR ≤ 0.05; log₂FC > 0"
    ),
    paste0("Final ", final_n, "-gene\ngastric metaplasia\nsignature")
  )

  x <- c(1.0, 2.4, 3.8, 5.4, 6.8, 8.2, 9.9, 11.7)
  box_width <- c(1.05, 1.15, 1.25, 1.05, 1.15, 1.25, 1.55, 1.35)
  box_height <- 1.15

  boxes <- data.frame(
    x = x,
    label = labels,
    xmin = x - box_width / 2,
    xmax = x + box_width / 2,
    ymin = -box_height / 2,
    ymax = box_height / 2,
    stringsAsFactors = FALSE
  )

  arrows <- data.frame(
    x = boxes$xmax[-nrow(boxes)] + 0.04,
    xend = boxes$xmin[-1] - 0.04,
    y = 0,
    yend = 0
  )

  headers <- data.frame(
    xmin = c(0.35, 4.70, 9.00),
    xmax = c(4.50, 8.95, 12.45),
    ymin = 0.94,
    ymax = 1.27,
    label = c("Discovery", "Cross-dataset validation", "SSL-specific refinement"),
    stringsAsFactors = FALSE
  )

  p <- ggplot() +
    geom_rect(
      data = headers,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "white", colour = "black", linewidth = 0.35
    ) +
    geom_text(
      data = headers,
      aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
      family = params$font_family,
      fontface = "bold",
      size = pt_to_mm(params$font_pt)
    ) +
    geom_segment(
      data = arrows,
      aes(x = x, xend = xend, y = y, yend = yend),
      colour = "grey65",
      linewidth = 0.8,
      arrow = grid::arrow(type = "closed", length = unit(2.5, "mm"))
    ) +
    geom_rect(
      data = boxes,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "white", colour = "black", linewidth = 0.35
    ) +
    geom_text(
      data = boxes,
      aes(x = x, y = 0, label = label),
      family = params$font_family,
      size = pt_to_mm(params$font_pt),
      lineheight = 0.92
    ) +
    coord_cartesian(xlim = c(0.30, 12.50), ylim = c(-0.80, 1.36), clip = "off") +
    theme_void(base_family = params$font_family) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(2, 2, 2, 2, "mm")
    )

  if (isTRUE(params$show_workflow_counts)) {
    count_lookup <- c(
      workflow_counts$n[1],
      workflow_counts$n[1],
      workflow_counts$n[2],
      workflow_counts$n[3],
      workflow_counts$n[4],
      workflow_counts$n[4],
      NA,
      final_n
    )
    count_df <- data.frame(
      x = x,
      y = -0.70,
      label = ifelse(is.na(count_lookup), "", paste0("n = ", count_lookup))
    )
    p <- p + geom_text(
      data = count_df,
      aes(x = x, y = y, label = label),
      family = params$font_family,
      size = pt_to_mm(params$font_pt - 0.5)
    )
  }

  p
}

make_signature_heatmap <- function(final_genes) {
  message("[3/9] Building the final-signature differential-expression heatmap...")

  contrast_order <- c(
    "SSL_HP_v_H_COL",
    "SSL_HP_v_H_ADJ_SSL_HP",
    "SSL_HP_v_CA"
  )

  contrast_labels <- c(
    SSL_HP_v_H_COL = "SSL_HP vs H_COL",
    SSL_HP_v_H_ADJ_SSL_HP = "SSL_HP vs H_ADJ_SSL_HP",
    SSL_HP_v_CA = "SSL_HP vs CA"
  )

  tables <- lapply(contrast_order, function(comp) {
    tab <- standardise_de_table(params$de_files[[comp]], require_fdr = FALSE)
    columns_used <- attr(tab, "columns_used")
    tab$comparison <- comp
    write_csv_safe(
      columns_used,
      file.path(params$out_dir, "qc", paste0(comp, "_DE_columns_used.csv"))
    )
    tab
  })
  names(tables) <- contrast_order

  heatmap_long <- dplyr::bind_rows(tables) %>%
    dplyr::filter(Gene %in% final_genes) %>%
    dplyr::select(comparison, Gene, log2FC, FDR, PValue, logCPM) %>%
    tidyr::complete(comparison = contrast_order, Gene = final_genes) %>%
    dplyr::mutate(
      comparison = factor(comparison, levels = contrast_order),
      comparison_label = factor(
        contrast_labels[as.character(comparison)],

        levels = rev(contrast_labels[contrast_order])
      ),
      Gene = factor(Gene, levels = final_genes)
    )

  missing_cells <- heatmap_long %>% dplyr::filter(!is.finite(log2FC))
  if (nrow(missing_cells) > 0) {
    write_csv_safe(
      missing_cells,
      file.path(params$out_dir, "qc", "FINAL_SIGNATURE_heatmap_missing_cells.csv")
    )
    stop("At least one final-signature gene is missing from a heatmap DE contrast.")
  }

  observed_limit <- ceiling(max(abs(heatmap_long$log2FC), na.rm = TRUE))
  colour_limit <- params$heatmap_colour_limit %||% observed_limit
  if (!is.finite(colour_limit) || colour_limit <= 0) colour_limit <- 1

  heatmap_long <- heatmap_long %>%
    dplyr::mutate(log2FC_display = pmax(pmin(log2FC, colour_limit), -colour_limit))

  write_csv_safe(
    heatmap_long,
    file.path(params$out_dir, "tables", "FINAL_24_GENE_SIGNATURE_heatmap_data.csv")
  )

  p <- ggplot(
    heatmap_long,
    aes(x = Gene, y = comparison_label, fill = log2FC_display)
  ) +
    geom_tile(colour = "white", linewidth = 0.30) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-colour_limit, colour_limit),
      oob = scales::squish,
      name = expression(log[2] * "FC")
    ) +
    scale_x_discrete(drop = FALSE, expand = expansion(add = 0.03)) +
    scale_y_discrete(drop = FALSE, expand = expansion(add = 0.03)) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_family = params$font_family, base_size = params$font_pt) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        family = params$font_family,
        size = params$x_font_pt,
        angle = 60,
        hjust = 1,
        vjust = 1,
        colour = "black"
      ),
      axis.text.y = element_text(
        family = params$font_family,
        size = params$font_pt,
        colour = "black"
      ),
      legend.title = element_text(
        family = params$font_family,
        size = params$font_pt
      ),
      legend.text = element_text(
        family = params$font_family,
        size = params$font_pt
      ),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(2, 4, 3, 4, "mm")
    )

  list(plot = p, data = heatmap_long, colour_limit = colour_limit)
}

format_gsea_p <- function(p) {
  p <- to_num(p)
  if (length(p) != 1 || !is.finite(p)) return("NA")
  if (p < 0.001) return("<0.001")
  as.character(signif(p, 3))
}

find_combined_gsea_prefix <- function(df, candidates, comparison_name) {
  for (candidate in candidates) {
    if (paste0(candidate, ".logFC") %in% names(df) ||
        paste0(candidate, ".log2FC") %in% names(df)) {
      return(candidate)
    }
  }

  stop(
    "No matching combined-DE prefix was found for ", comparison_name,
    ". Tried: ", paste(candidates, collapse = ", ")
  )
}

build_m26f_prerank <- function(df, prefix, comparison_name) {
  gene_col <- pick_first_existing(
    c(
      paste0(prefix, ".SYMBOL"),
      paste0(prefix, ".symbol"),
      paste0(prefix, ".Gene"),
      paste0(prefix, ".gene")
    ),
    df
  )

  logfc_col <- pick_first_existing(
    c(paste0(prefix, ".logFC"), paste0(prefix, ".log2FC")),
    df
  )

  signed_stat_col <- pick_first_existing(
    c(
      paste0(prefix, ".t"),
      paste0(prefix, ".stat"),
      paste0(prefix, ".Wald"),
      paste0(prefix, ".z")
    ),
    df
  )

  p_col <- pick_first_existing(
    c(
      paste0(prefix, ".PValue"),
      paste0(prefix, ".pvalue"),
      paste0(prefix, ".P.Value"),
      paste0(prefix, ".p_val")
    ),
    df
  )

  fdr_col <- pick_first_existing(
    c(
      paste0(prefix, ".FDR"),
      paste0(prefix, ".adj.P.Val"),
      paste0(prefix, ".padj")
    ),
    df
  )

  if (is.na(gene_col) || is.na(logfc_col)) {
    stop(
      "The combined DE table lacks a gene-symbol or logFC column for ",
      comparison_name, " (prefix: ", prefix, ")."
    )
  }

  gene_ann <- prepare_gene_symbols(df[[gene_col]])
  log2fc <- to_num(df[[logfc_col]])
  pvalue <- if (!is.na(p_col)) to_num(df[[p_col]]) else rep(NA_real_, nrow(df))
  fdr <- if (!is.na(fdr_col)) to_num(df[[fdr_col]]) else rep(NA_real_, nrow(df))

  if (!is.na(signed_stat_col)) {
    rank_metric <- to_num(df[[signed_stat_col]])
    rank_metric_used <- "signed statistic column (t/stat/Wald/z)"
  } else if (!is.na(p_col)) {
    p_safe <- pmax(pvalue, 1e-300)
    rank_metric <- sign(log2fc) * -log10(p_safe)
    rank_metric_used <- "sign(log2FC) * -log10(nominal P value)"
  } else if (!is.na(fdr_col)) {
    fdr_safe <- pmax(fdr, 1e-300)
    rank_metric <- sign(log2fc) * -log10(fdr_safe)
    rank_metric_used <- "sign(log2FC) * -log10(FDR)"
  } else {
    rank_metric <- log2fc
    rank_metric_used <- "log2FC"
  }

  ranked <- data.frame(
    Gene = gene_ann$gene_symbol,
    log2FC = log2fc,
    PValue = pvalue,
    FDR = fdr,
    rank_metric = rank_metric,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(
      !is.na(Gene),
      nzchar(Gene),
      is.finite(log2FC),
      is.finite(rank_metric)
    ) %>%
    dplyr::group_by(Gene) %>%
    dplyr::arrange(
      dplyr::desc(abs(rank_metric)),
      PValue,
      FDR,
      dplyr::desc(abs(log2FC)),
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(rank_metric), Gene)

  if (nrow(ranked) < 2) {
    stop("Too few genes remained to run GSEA for ", comparison_name, ".")
  }

  ranks <- ranked$rank_metric
  names(ranks) <- ranked$Gene

  if (anyDuplicated(ranks) > 0) {
    scale_value <- max(1, max(abs(ranks), na.rm = TRUE))
    ranks <- ranks - seq_along(ranks) * (.Machine$double.eps * scale_value)
    ranked$rank_metric <- ranks
  }

  ranks <- sort(ranks, decreasing = TRUE)
  ranked <- ranked[match(names(ranks), ranked$Gene), , drop = FALSE]

  columns_used <- data.frame(
    comparison = comparison_name,
    prefix = prefix,
    gene_col = gene_col,
    logfc_col = logfc_col,
    signed_stat_col = signed_stat_col,
    p_col = p_col,
    fdr_col = fdr_col,
    rank_metric_used = rank_metric_used,
    n_ranked_genes = length(ranks),
    stringsAsFactors = FALSE
  )

  list(
    table = ranked,
    ranks = ranks,
    columns_used = columns_used,
    rank_metric_used = rank_metric_used
  )
}

calculate_running_enrichment <- function(ranks, pathway_genes, gsea_param = 1) {
  ranks <- sort(ranks, decreasing = TRUE)
  pathway_genes <- standardise_gene_vector(pathway_genes)
  pathway_genes <- base::intersect(pathway_genes, names(ranks))

  n_total <- length(ranks)
  n_hits <- length(pathway_genes)

  if (n_hits < 1 || n_hits >= n_total) {
    stop("The final signature has an invalid overlap with a GSEA ranked list.")
  }

  is_hit <- names(ranks) %in% pathway_genes
  weights <- abs(ranks)^gsea_param
  normaliser <- sum(weights[is_hit])

  if (!is.finite(normaliser) || normaliser <= 0) {
    weights <- rep(1, length(ranks))
    normaliser <- sum(weights[is_hit])
  }

  running_es <- cumsum(ifelse(is_hit, weights / normaliser, 0)) -
    cumsum(ifelse(!is_hit, 1 / (n_total - n_hits), 0))

  curve <- data.frame(
    rank_index = seq_along(ranks),
    rank_position_percent = 100 * seq_along(ranks) / n_total,
    running_enrichment_score = running_es,
    Gene = names(ranks),
    rank_metric = as.numeric(ranks),
    signature_hit = is_hit,
    stringsAsFactors = FALSE
  )

  list(
    curve = curve,
    hits = curve[curve$signature_hit, , drop = FALSE]
  )
}

run_final_signature_gsea <- function(final_genes) {
  message(
    "[4/9] Running three-contrast GSEA overlay for the final 24-gene signature..."
  )

  final_genes <- standardise_gene_vector(final_genes)
  combined_de <- data.table::fread(
    params$combined_de_file,
    data.table = FALSE,
    check.names = FALSE
  )

  comparison_meta <- data.frame(
    comparison = params$gsea_comparison_order,
    label = unname(params$gsea_comparison_labels[params$gsea_comparison_order]),
    colour = unname(params$gsea_comparison_colours[params$gsea_comparison_order]),
    stringsAsFactors = FALSE
  )

  rank_objects <- list()
  prefix_rows <- list()
  overlap_list <- list()
  missing_list <- list()

  for (comparison_name in params$gsea_comparison_order) {
    prefix <- find_combined_gsea_prefix(
      combined_de,
      params$gsea_prefix_candidates[[comparison_name]],
      comparison_name
    )

    rank_obj <- build_m26f_prerank(
      combined_de,
      prefix = prefix,
      comparison_name = comparison_name
    )

    overlap <- base::intersect(final_genes, names(rank_obj$ranks))
    missing_from_rank <- base::setdiff(final_genes, overlap)

    message(
      "    ", params$gsea_comparison_labels[[comparison_name]], ": ",
      length(overlap), "/", length(final_genes),
      " final-signature genes present in the ranked list."
    )
    if (length(missing_from_rank) > 0) {
      message(
        "      Absent from this ranked universe: ",
        paste(missing_from_rank, collapse = ", ")
      )
    }

    if (params$require_all_final_signature_genes_in_gsea_rank &&
        length(missing_from_rank) > 0) {
      stop(
        "The following final-signature genes were absent from the ",
        comparison_name, " GSEA ranked list: ",
        paste(missing_from_rank, collapse = ", ")
      )
    }

    if (length(overlap) < params$gsea_min_size) {
      stop(
        "Only ", length(overlap), " signature genes overlapped the ",
        comparison_name, " ranked list; gsea_min_size = ",
        params$gsea_min_size, "."
      )
    }

    rank_objects[[comparison_name]] <- rank_obj
    overlap_list[[comparison_name]] <- overlap
    missing_list[[comparison_name]] <- missing_from_rank
    prefix_rows[[comparison_name]] <- rank_obj$columns_used
  }

  ranking_methods <- unique(vapply(
    rank_objects,
    function(x) x$rank_metric_used,
    character(1)
  ))

  if (length(ranking_methods) != 1) {
    stop(
      "Different ranking metrics would be used across the three GSEA contrasts: ",
      paste(ranking_methods, collapse = " | "),
      ". A publication figure must use one consistent ranking method."
    )
  }

  set.seed(1)
  fgsea_rows <- list()
  curve_rows <- list()
  hit_rows <- list()
  position_rows <- list()

  for (i in seq_along(params$gsea_comparison_order)) {
    comparison_name <- params$gsea_comparison_order[i]
    ranks <- rank_objects[[comparison_name]]$ranks
    overlap <- overlap_list[[comparison_name]]

    pathway_name <- "FINAL_24_GENE_GASTRIC_METAPLASIA_SIGNATURE"

    fg <- fgsea::fgseaMultilevel(
      pathways = stats::setNames(list(overlap), pathway_name),
      stats = ranks,
      minSize = params$gsea_min_size,
      maxSize = params$gsea_max_size,
      eps = params$gsea_eps,
      scoreType = params$gsea_score_type,
      nproc = 0,
      nPermSimple = params$gsea_nperm_simple
    )

    fg <- as.data.frame(fg)
    if (nrow(fg) != 1) {
      stop("FGSEA did not return exactly one result for ", comparison_name, ".")
    }

    if ("leadingEdge" %in% names(fg)) {
      fg$leadingEdge <- vapply(
        fg$leadingEdge,
        paste,
        collapse = ";",
        FUN.VALUE = character(1)
      )
    }

    fg$comparison <- comparison_name
    fg$comparison_label <- params$gsea_comparison_labels[[comparison_name]]
    fg$rank_metric <- ranking_methods[1]
    fg$ranked_gene_count <- length(ranks)
    fg$signature_gene_count_intended <- length(final_genes)
    fg$signature_gene_count_in_rank <- length(overlap)
    fgsea_rows[[comparison_name]] <- fg

    enrichment <- calculate_running_enrichment(
      ranks,
      overlap,
      gsea_param = params$gsea_gsea_param
    )

    curve_rows[[comparison_name]] <- enrichment$curve %>%
      dplyr::mutate(
        comparison = comparison_name,
        comparison_label = params$gsea_comparison_labels[[comparison_name]],
        colour = params$gsea_comparison_colours[[comparison_name]]
      )

    hit_rows[[comparison_name]] <- enrichment$hits %>%
      dplyr::mutate(
        comparison = comparison_name,
        comparison_label = params$gsea_comparison_labels[[comparison_name]],
        colour = params$gsea_comparison_colours[[comparison_name]]
      )

    position_rows[[comparison_name]] <- data.frame(
      comparison = comparison_name,
      comparison_label = params$gsea_comparison_labels[[comparison_name]],
      Gene = final_genes,
      present_in_ranked_list = final_genes %in% names(ranks),
      rank_position = match(final_genes, names(ranks)),
      rank_position_percent = 100 * match(final_genes, names(ranks)) / length(ranks),
      rank_metric = unname(ranks[final_genes]),
      stringsAsFactors = FALSE
    )

    write_csv_safe(
      rank_objects[[comparison_name]]$table,
      file.path(
        params$out_dir,
        "tables",
        paste0("FINAL_24_GENE_GSEA_", comparison_name, "_ranked_genes.csv")
      )
    )
  }

  results_df <- dplyr::bind_rows(fgsea_rows) %>%
    dplyr::mutate(
      comparison = factor(comparison, levels = params$gsea_comparison_order)
    ) %>%
    dplyr::arrange(comparison) %>%
    dplyr::mutate(comparison = as.character(comparison))

  curve_df <- dplyr::bind_rows(curve_rows)
  hit_df <- dplyr::bind_rows(hit_rows)
  positions_df <- dplyr::bind_rows(position_rows)
  columns_df <- dplyr::bind_rows(prefix_rows)

  write_csv_safe(
    results_df,
    file.path(params$out_dir, "tables", "FINAL_24_GENE_GSEA_multicomparison_results.csv")
  )
  write_csv_safe(
    curve_df,
    file.path(params$out_dir, "tables", "FINAL_24_GENE_GSEA_multicomparison_curves.csv")
  )
  write_csv_safe(
    hit_df,
    file.path(params$out_dir, "tables", "FINAL_24_GENE_GSEA_multicomparison_hits.csv")
  )
  write_csv_safe(
    positions_df,
    file.path(params$out_dir, "tables", "FINAL_24_GENE_GSEA_signature_gene_positions_all_comparisons.csv")
  )
  write_csv_safe(
    columns_df,
    file.path(params$out_dir, "qc", "FINAL_24_GENE_GSEA_columns_and_ranking_all_comparisons.csv")
  )

  missing_long <- dplyr::bind_rows(lapply(
    names(missing_list),
    function(x) {
      missing_genes <- missing_list[[x]]
      if (length(missing_genes) == 0) {
        return(data.frame(
          comparison = character(0),
          missing_gene = character(0),
          stringsAsFactors = FALSE
        ))
      }
      data.frame(
        comparison = rep(x, length(missing_genes)),
        missing_gene = missing_genes,
        stringsAsFactors = FALSE
      )
    }
  ))
  write_csv_safe(
    missing_long,
    file.path(params$out_dir, "qc", "FINAL_24_GENE_GSEA_missing_signature_genes_all_comparisons.csv")
  )

  curve_df$comparison <- factor(
    curve_df$comparison,
    levels = params$gsea_comparison_order
  )
  hit_df$comparison <- factor(
    hit_df$comparison,
    levels = params$gsea_comparison_order
  )

  curve_min <- min(curve_df$running_enrichment_score, na.rm = TRUE)
  curve_max <- max(curve_df$running_enrichment_score, na.rm = TRUE)
  y_lower <- min(-0.5, floor((curve_min - 0.05) * 10) / 10)
  y_upper <- max(1.05, ceiling((curve_max + 0.05) * 10) / 10)
  y_span <- y_upper - y_lower

  hit_df$tick_y0 <- y_lower + 0.035 * y_span
  hit_df$tick_y1 <- y_lower + 0.078 * y_span

  annotation_layout <- data.frame(
    comparison = params$gsea_comparison_order,
    x_line_start = c(46, 46, 76),
    x_line_end = c(51, 51, 81),
    x_text = c(53, 53, 83),
    vertical_slot = c("top", "second", "top"),
    stringsAsFactors = FALSE
  )

  annotation_df <- comparison_meta %>%
    dplyr::left_join(
      results_df %>%
        dplyr::select(comparison, NES, pval, padj),
      by = "comparison"
    ) %>%
    dplyr::left_join(annotation_layout, by = "comparison") %>%
    dplyr::mutate(
      y = dplyr::if_else(
        vertical_slot == "top",
        y_upper - 0.045 * y_span,
        y_upper - 0.205 * y_span
      ),
      p_text = dplyr::if_else(
        is.finite(pval) & pval < 0.001,
        "p < 0.001",
        paste0("p = ", vapply(pval, format_gsea_p, character(1)))
      ),
      annotation = paste0(
        label,
        "\nNES = ", sprintf("%.2f", NES), ", ", p_text
      )
    )

  write_csv_safe(
    annotation_df,
    file.path(params$out_dir, "tables", "FINAL_24_GENE_GSEA_plot_annotations.csv")
  )

  colour_map <- stats::setNames(
    params$gsea_comparison_colours[params$gsea_comparison_order],
    params$gsea_comparison_order
  )

  p <- ggplot() +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.35,
      colour = "grey50"
    ) +
    geom_line(
      data = curve_df,
      aes(
        x = rank_position_percent,
        y = running_enrichment_score,
        colour = comparison
      ),
      linewidth = params$gsea_curve_linewidth,
      lineend = "round"
    ) +
    geom_segment(
      data = hit_df,
      aes(
        x = rank_position_percent,
        xend = rank_position_percent,
        y = tick_y0,
        yend = tick_y1,
        colour = comparison
      ),
      linewidth = params$gsea_tick_linewidth,
      alpha = 0.95
    ) +
    geom_segment(
      data = annotation_df,
      aes(
        x = x_line_start,
        xend = x_line_end,
        y = y,
        yend = y,
        colour = comparison
      ),
      linewidth = params$gsea_curve_linewidth,
      lineend = "butt"
    ) +
    geom_text(
      data = annotation_df,
      aes(x = x_text, y = y, label = annotation),
      hjust = 0,
      vjust = 0.5,
      family = params$font_family,
      size = pt_to_mm(params$gsea_annotation_pt),
      lineheight = 0.88,
      colour = "black"
    ) +
    scale_colour_manual(
      values = colour_map,
      breaks = params$gsea_comparison_order,
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = c(0, 25, 50, 75, 100),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      breaks = pretty(c(y_lower, y_upper), n = 4),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    coord_cartesian(
      xlim = c(0, 100),
      ylim = c(y_lower, y_upper),
      clip = "off"
    ) +
    labs(
      title = params$gsea_plot_title,
      x = "Ranked gene list position (%)",
      y = "Running enrichment score"
    ) +
    theme_classic(
      base_family = params$font_family,
      base_size = params$gsea_font_base_pt
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        family = params$font_family,
        size = params$gsea_title_pt,
        face = "bold",
        hjust = 0.5,
        colour = "black",
        margin = margin(b = 5)
      ),
      axis.title = element_text(
        family = params$font_family,
        size = params$gsea_axis_title_pt,
        colour = "black"
      ),
      axis.text = element_text(
        family = params$font_family,
        size = params$gsea_axis_text_pt,
        colour = "black"
      ),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = unit(1.2, "mm"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(4, 16, 5, 5, "mm")
    )

  p_with_panel_letter <- patchwork::wrap_plots(p) +
    patchwork::plot_annotation(
      tag_levels = list(c(params$gsea_panel_letter)),
      theme = theme(
        plot.tag = element_text(
          family = params$font_family,
          face = "bold",
          size = params$gsea_panel_label_pt,
          colour = "black"
        ),
        plot.tag.position = c(0.005, 0.995)
      )
    )

  list(
    plot = p,
    plot_with_panel_letter = p_with_panel_letter,
    results = results_df,
    curves = curve_df,
    hits = hit_df,
    ranks = lapply(rank_objects, `[[`, "ranks"),
    overlaps = overlap_list,
    missing_from_rank = missing_list,
    rank_metric_used = ranking_methods[1],
    annotation_data = annotation_df
  )
}

run_pipeline <- function() {
  required_files <- c(
    params$gtex_file,
    params$hpa_file,
    params$combined_de_file,
    unname(params$de_files)
  )
  assert_files_exist(required_files)

  candidate_obj <- generate_stomach_candidates()
  final_obj <- generate_final_signature(candidate_obj$stomach_candidates)

  workflow_plot <- make_workflow_panel(
    candidate_obj$workflow_counts,
    length(final_obj$final_genes)
  )
  heatmap_obj <- make_signature_heatmap(final_obj$final_genes)
  gsea_obj <- run_final_signature_gsea(final_obj$final_genes)

  save_plot_set(
    workflow_plot,
    file.path(params$out_dir, "figures", "individual_panels", "Panel_A_signature_workflow"),
    width_mm = params$figure_width_mm,
    height_mm = 50
  )
  save_plot_set(
    heatmap_obj$plot,
    file.path(params$out_dir, "figures", "individual_panels", "Panel_B_24_gene_log2FC_heatmap"),
    width_mm = params$figure_width_mm,
    height_mm = 48
  )
  save_plot_set(
    gsea_obj$plot,
    file.path(params$out_dir, "figures", "individual_panels", "Panel_C_FINAL_24_gene_GSEA_multicomparison"),
    width_mm = 105,
    height_mm = 72
  )
  save_plot_set(
    gsea_obj$plot_with_panel_letter,
    file.path(params$out_dir, "figures", "individual_panels", "Panel_F_Gastric_metaplasia_signature_enrichment"),
    width_mm = params$figure_width_mm,
    height_mm = 120
  )
  save_plot_set(
    gsea_obj$plot,
    file.path(params$out_dir, "figures", "thesis", "Figure_3_7C_gastric_metaplasia_signature_enrichment"),
    width_mm = 120,
    height_mm = 82
  )

  figure_3_7 <- wrap_plots(
    workflow_plot,
    heatmap_obj$plot,
    gsea_obj$plot,
    ncol = 1,
    heights = c(0.75, 0.75, 1.0)
  ) +
    plot_annotation(
      tag_levels = "A",
      theme = theme(
        plot.tag = element_text(
          family = params$font_family,
          face = "bold",
          size = 8
        )
      )
    )

  save_plot_set(
    figure_3_7,
    file.path(params$out_dir, "figures", "thesis", "Figure_3_7_signature_development"),
    width_mm = params$figure_width_mm,
    height_mm = 170
  )

  capture.output(sessionInfo(), file=file.path(params$out_dir, "SESSION_INFO.txt"))
  message("Signature discovery and gastric GSEA complete. Use 090_complete_score_analysis.R for score tests and score figures.")
  invisible(list(candidate_obj=candidate_obj, final_obj=final_obj, gsea_obj=gsea_obj))
}
results <- run_pipeline()
