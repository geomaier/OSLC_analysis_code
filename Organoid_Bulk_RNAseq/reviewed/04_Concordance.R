# Tissue-organoid concordance of differential-expression and pathway results.
# Reads existing DE tables; input requirements and output locations are in ../README.md.

rm(list = ls())

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

# Configuration and comparison manifest ----
set.seed(1)
options(stringsAsFactors = FALSE)
options(contrasts = c("contr.treatment", "contr.poly"))

input_dir <- required_path("OSLC_CONCORDANCE_DIR")

metaplasia_signature_candidates <- c(
  file.path(input_dir, "FINAL_SSL_METAPLASIA_LIST.txt"),
  file.path(input_dir, "FINAL_SSL_METAPLASIA_LIST.csv"),
  file.path(
    dirname(input_dir),
    "CPM_TPM_signature_scoring_final",
    "tables",
    "FINAL_SSL_METAPLASIA_LIST.txt"
  ),
  file.path(
    dirname(input_dir),
    "CPM_TPM_signature_scoring_final",
    "tables",
    "FINAL_SSL_METAPLASIA_LIST.csv"
  )
)

save_outputs_in_temp_first <- FALSE
if (save_outputs_in_temp_first) {
  out_dir <- file.path(tempdir(), "M26F_DE_PATHWAY_GENE_CONCORDANCE_FINAL_ALL_DE_TABLES")
} else {
  out_dir <- file.path(
    input_dir,
    "M26F_DE_PATHWAY_GENE_CONCORDANCE_FINAL_ALL_DE_TABLES"
  )
}

out_dir <- review_output_dir(out_dir)

overwrite_existing_outputs <- FALSE
install_missing_packages <- FALSE

require_all_manifest_files <- TRUE

hallmark_logcpm_cutoff <- 2
fdr_cutoff <- 0.05
strong_logfc_cutoff <- 1
fgsea_min_size <- 15L
fgsea_max_size <- 500L
fgsea_eps <- 0
fgsea_nperm_simple <- 10000L
fgsea_nproc <- 0L

font_family <- "Arial"
export_dpi <- 600
heatmap_logfc_limit <- 2.5
max_driver_genes_per_module <- 40L
max_scatter_labels_per_module <- 6L

write_alt_vs_all_individual_outputs <- TRUE

write_all_pairwise_individual_outputs <- TRUE

write_focused_ssl_infection_braf_outputs <- TRUE
focused_max_genes_per_pathway <- 20L
focused_minimum_significant_leading_edges <- 2L
focused_logfc_colour_limit <- 2

primary_comparison_order <- c(
  "Alt_Infected",
  "SSL_HP_v_H_COL",
  "SSL_HP_v_H_ADJ_SSL_HP",
  "ES_SSL_HP",
  "ES_H_ADJ_SSL_HP"
)

all_comparison_order <- c(
  "Alt_Infected",
  "SSL_HP_v_H_COL",
  "SSL_HP_v_H_ADJ_SSL_HP",
  "H_ADJ_SSL_HP_v_H_COL",
  "Serr_SSL_HP",
  "Serr_H_ADJ_SSL_HP",
  "Fuso_SSL_HP",
  "Fuso_H_ADJ_SSL_HP",
  "ES_SSL_HP",
  "ES_H_ADJ_SSL_HP",
  "V600E_v_WT_in_SSL_HP",
  "Initial_BRAF_SSL_HP_organoid"
)

de_file_table <- tibble::tribble(
  ~Comparison, ~Display, ~Short_label, ~Positive_logFC_means, ~Primary, ~Required, ~File,
  "Alt_Infected",
  "Infected vs Non-Infected Organoids",
  "Infected vs Non-Infected Organoids",
  "Higher in infected organoids",
  TRUE, TRUE,
  "CGG028_Alt_Infected_DE_results.txt",

  "SSL_HP_v_H_COL",
  "SSL_HP vs H_COL",
  "SSL_HP vs H_COL",
  "Higher in SSL_HP than H_COL",
  TRUE, TRUE,
  "CGG014_Maier_SSL_HP_v_H_COL_DE_results.txt",

  "SSL_HP_v_H_ADJ_SSL_HP",
  "SSL_HP vs H_ADJ_SSL_HP",
  "SSL_HP vs H_ADJ",
  "Higher in SSL_HP than H_ADJ_SSL_HP",
  TRUE, TRUE,
  "CGG014_Maier_SSL_HP_v_H_ADJ_SSL_HP_DE_results.txt",

  "ES_SSL_HP",
  "Escherichia-Shigella + vs - in SSL_HP",
  "ES+ vs ES- SSL_HP",
  "Higher in Escherichia-Shigella-positive SSL_HP",
  TRUE, TRUE,
  "CGG014_Maier_ES_SSL_HP_DE_results.txt",

  "ES_H_ADJ_SSL_HP",
  "Escherichia-Shigella + vs - in H_ADJ_SSL_HP",
  "ES+ vs ES- H_ADJ",
  "Higher in Escherichia-Shigella-positive H_ADJ_SSL_HP",
  TRUE, TRUE,
  "CGG014_Maier_ES_H_ADJ_SSL_HP_DE_results.txt",

  "H_ADJ_SSL_HP_v_H_COL",
  "H_ADJ_SSL_HP vs H_COL",
  "H_ADJ vs H_COL",
  "Higher in H_ADJ_SSL_HP than H_COL",
  FALSE, FALSE,
  "CGG014_Maier_H_ADJ_SSL_HP_v_H_COL_DE_results.txt",

  "Fuso_SSL_HP",
  "Fusobacterium + vs - in SSL_HP",
  "Fuso+ vs - SSL_HP",
  "Higher in Fusobacterium-positive SSL_HP",
  FALSE, FALSE,
  "CGG014_Maier_Fuso_SSL_HP_DE_results.txt",

  "Fuso_H_ADJ_SSL_HP",
  "Fusobacterium + vs - in H_ADJ_SSL_HP",
  "Fuso+ vs - H_ADJ",
  "Higher in Fusobacterium-positive H_ADJ_SSL_HP",
  FALSE, FALSE,
  "CGG014_Maier_Fuso_H_ADJ_SSL_HP_DE_results.txt",

  "Serr_SSL_HP",
  "Serratia + vs - in SSL_HP",
  "Serr+ vs - SSL_HP",
  "Higher in Serratia-positive SSL_HP",
  FALSE, FALSE,
  "CGG014_Maier_Serr_SSL_HP_DE_results.txt",

  "Serr_H_ADJ_SSL_HP",
  "Serratia + vs - in H_ADJ_SSL_HP",
  "Serr+ vs - H_ADJ",
  "Higher in Serratia-positive H_ADJ_SSL_HP",
  FALSE, FALSE,
  "CGG014_Maier_Serr_H_ADJ_SSL_HP_DE_results.txt",

  "V600E_v_WT_in_SSL_HP",
  "BRAF V600E vs WT in SSL_HP tissue",
  "Tissue V600E vs WT",
  "Higher in BRAF V600E SSL_HP tissue",
  FALSE, FALSE,
  "CGG014_Maier_V600E_v_WT_in_SSL_HP_DE_results.txt",

  "Initial_BRAF_SSL_HP_organoid",
  "BRAF V600E vs WT Organoids",
  "BRAF V600E vs WT Organoids",
  "Higher in BRAF V600E SSL_HP organoids",
  FALSE, FALSE,
  "CGG028_First_Passage_SSL_HP_V600E_WT.txt"
)

# Required packages ----
cran_packages <- c(
  "data.table", "dplyr", "tidyr", "purrr", "stringr", "tibble",
  "ggplot2", "ggrepel", "patchwork", "scales", "msigdbr"
)
bioc_packages <- c("fgsea", "AnnotationDbi", "org.Hs.eg.db")

install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!install_missing_packages) {
      stop(
        "Package '", pkg, "' is not installed. Install it or set ",
        "install_missing_packages <- TRUE."
      )
    }
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

for (pkg in cran_packages) install_if_missing(pkg, FALSE)
for (pkg in bioc_packages) install_if_missing(pkg, TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(msigdbr)
  library(fgsea)
  library(grid)
})

if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)
if (overwrite_existing_outputs && dir.exists(out_dir)) {
  unlink(out_dir, recursive = TRUE, force = TRUE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

subdirs <- list(
  qc = file.path(out_dir, "01_QC"),
  gsea = file.path(out_dir, "02_GSEA_and_leading_edges"),
  primary_plots = file.path(out_dir, "03_Primary_pairwise_scatterplots"),
  alt_plots = file.path(out_dir, "03_Infection_vs_each_available_scatterplots"),
  all_pairwise_plots = file.path(out_dir, "03_ALL_pairwise_scatterplots"),
  all_pairwise_single_module_plots = file.path(out_dir, "03_ALL_pairwise_single_module_scatterplots"),
  primary_tables = file.path(out_dir, "04_Primary_pairwise_gene_tables"),
  alt_tables = file.path(out_dir, "04_Infection_vs_each_available_gene_tables"),
  all_pairwise_tables = file.path(out_dir, "04_ALL_pairwise_gene_tables"),
  heatmaps = file.path(out_dir, "05_Driver_gene_heatmaps"),
  summary_plots = file.path(out_dir, "06_Pairwise_summary_heatmaps"),
  tables = file.path(out_dir, "07_Tables"),
  focused = file.path(out_dir, "08_Focused_SSL_Infection_BRAF")
)
for (d in subdirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Input and plotting helpers ----
clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

column_key <- function(x) {
  toupper(gsub("[^A-Za-z0-9]", "", trimws(as.character(x))))
}

safe_filename <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
}

detect_col <- function(df, candidates, required = FALSE, label = "column") {
  nms <- names(df)
  exact <- candidates[candidates %in% nms]
  if (length(exact) > 0L) return(exact[1])
  idx <- match(column_key(candidates), column_key(nms))
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0L) return(nms[idx[1]])
  if (required) {
    stop(
      "Could not detect ", label, ". Available columns:\n",
      paste(nms, collapse = "\n")
    )
  }
  NA_character_
}

first_existing_file <- function(paths) {
  hit <- unique(paths)[file.exists(unique(paths))]
  if (length(hit) == 0L) NA_character_ else hit[1]
}

read_gene_list <- function(file) {
  if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    dat <- fread(file, data.table = FALSE, check.names = FALSE)
    gene_col <- detect_col(
      dat,
      c("gene_symbol", "Gene", "gene", "SYMBOL", "Symbol", "symbol"),
      FALSE
    )
    if (is.na(gene_col)) gene_col <- names(dat)[1]
    genes <- dat[[gene_col]]
  } else {
    lines <- readLines(file, warn = FALSE)
    genes <- if (length(lines) == 1L && grepl("[,;\\t]", lines[1])) {
      unlist(strsplit(lines, "[,;\\t]"))
    } else {
      lines
    }
  }
  genes <- unique(toupper(trimws(as.character(genes))))
  genes[!is.na(genes) & genes != ""]
}

nature_theme <- function(base_size = 7) {
  theme_classic(base_family = font_family, base_size = base_size) +
    theme(
      text = element_text(family = font_family, colour = "black"),
      plot.title = element_text(face = "bold", size = 8, hjust = 0),
      plot.subtitle = element_text(size = 6, hjust = 0),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6, colour = "black"),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5.5),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.35),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(4, 4, 4, 4, "mm")
    )
}

save_plot_set <- function(plot, directory, stem, width_mm, height_mm) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(directory, safe_filename(stem))
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

  ggsave(
    paste0(prefix, ".pdf"), plot,
    width = width_mm, height = height_mm, units = "mm",
    device = pdf_device, bg = "white"
  )
  ggsave(
    paste0(prefix, ".png"), plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = export_dpi, bg = "white"
  )
  ggsave(
    paste0(prefix, ".tif"), plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = export_dpi, device = grDevices::tiff, compression = "lzw", bg = "white"
  )
}

make_empty_plot <- function(title_text, message_text) {
  ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5, label = message_text,
      family = font_family, size = 2.2
    ) +
    xlim(0, 1) + ylim(0, 1) +
    labs(title = title_text) +
    theme_void(base_family = font_family) +
    theme(plot.title = element_text(face = "bold", size = 8, hjust = 0))
}

# Required DE-table checks ----
resolve_direct_de_file <- function(file_name) {
  exact_path <- file.path(input_dir, file_name)
  if (file.exists(exact_path)) return(exact_path)

  candidates <- list.files(
    input_dir,
    full.names = TRUE,
    recursive = FALSE,
    include.dirs = FALSE
  )
  hits <- candidates[
    tolower(basename(candidates)) == tolower(file_name)
  ]

  if (length(hits) > 1L) {
    stop(
      "More than one case-insensitive match was found for ",
      file_name,
      " in input_dir:\n",
      paste(hits, collapse = "\n")
    )
  }

  if (length(hits) == 1L) return(hits[1])
  exact_path
}

de_file_table <- de_file_table %>%
  mutate(
    Required = if (require_all_manifest_files) TRUE else Required,
    Path = vapply(File, resolve_direct_de_file, character(1)),
    File_found = file.exists(Path),
    Included_in_full_analysis = File_found,
    Included_in_ALT_vs_each_outputs =
      File_found & Comparison != "Alt_Infected",
    Primary_column_meaning = ifelse(
      Primary,
      "Original five-comparison core subset",
      "Additional comparison; still fully analysed when File_found is TRUE"
    )
  )

missing_required <- de_file_table %>% filter(Required, !File_found)
if (nrow(missing_required) > 0L) {
  stop(
    "Required DE file(s) are missing directly from input_dir:\n",
    paste0(
      missing_required$Comparison,
      " -> ",
      missing_required$Path,
      collapse = "\n"
    ),
    "\n\nAll concordance DE inputs must already exist as DE-result tables in ",
    "the Concordance analysis folder. In particular, use ",
    "CGG028_Alt_Infected_DE_results.txt for infection and ",
    "CGG028_First_Passage_SSL_HP_V600E_WT.txt for BRAF V600E vs WT organoids. ",
    "No DE table is generated or replaced by this script."
  )
}

missing_optional <- de_file_table %>% filter(!Required, !File_found)
if (nrow(missing_optional) > 0L) {
  warning(
    "Optional comparison(s) skipped because their DE files were absent from input_dir:\n",
    paste0(
      missing_optional$Comparison,
      " -> ",
      missing_optional$Path,
      collapse = "\n"
    )
  )
}

active_manifest <- de_file_table %>%
  filter(File_found) %>%
  arrange(factor(Comparison, levels = all_comparison_order))
active_comparisons <- active_manifest$Comparison
primary_comparisons <- primary_comparison_order[
  primary_comparison_order %in% active_comparisons
]
if (length(primary_comparisons) != length(primary_comparison_order)) {
  stop("Not all five primary comparisons were found.")
}

fwrite(
  de_file_table,
  file.path(subdirs$qc, "COMPARISON_FILE_AND_ORIENTATION_MANIFEST.csv")
)

metaplasia_signature_file <- first_existing_file(
  metaplasia_signature_candidates
)
if (is.na(metaplasia_signature_file)) {
  stop(
    "FINAL_SSL_METAPLASIA_LIST was not found. Checked:\n",
    paste(metaplasia_signature_candidates, collapse = "\n")
  )
}
metaplasia_genes <- read_gene_list(metaplasia_signature_file)
if (length(metaplasia_genes) != 24L) {
  stop(
    "Expected exactly 24 unique final gastric-metaplasia genes, but found ",
    length(metaplasia_genes), ". File: ", metaplasia_signature_file
  )
}

# Pathway definitions ----
get_hallmark_collection <- function() {
  out <- tryCatch(
    suppressWarnings(msigdbr(
      species = "Homo sapiens",
      collection = "H"
    )),
    error = function(e) NULL
  )
  if (is.null(out) || nrow(out) == 0L) {
    out <- tryCatch(
      suppressWarnings(msigdbr(
        species = "Homo sapiens",
        category = "H"
      )),
      error = function(e) NULL
    )
  }
  if (is.null(out) || nrow(out) == 0L) {
    stop("Could not load the MSigDB Hallmark collection.")
  }
  as.data.frame(out)
}

hallmark_raw <- get_hallmark_collection()
hallmark_gene_col <- detect_col(
  hallmark_raw,
  c("gene_symbol", "human_gene_symbol", "db_gene_symbol", "symbol"),
  TRUE,
  "MSigDB Hallmark gene-symbol column"
)
hallmark_pathway_col <- detect_col(
  hallmark_raw,
  c("gs_name", "geneset", "gene_set_name", "gs_exact_source"),
  TRUE,
  "MSigDB Hallmark pathway column"
)

hallmark_df <- hallmark_raw %>%
  transmute(
    pathway = as.character(.data[[hallmark_pathway_col]]),
    gene = toupper(trimws(as.character(.data[[hallmark_gene_col]])))
  ) %>%
  filter(!is.na(pathway), pathway != "", !is.na(gene), gene != "") %>%
  distinct()

hallmark_pathways <- split(hallmark_df$gene, hallmark_df$pathway) %>%
  lapply(function(x) sort(unique(toupper(x))))

target_hallmark_names <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)
missing_hallmarks <- setdiff(target_hallmark_names, names(hallmark_pathways))
if (length(missing_hallmarks) > 0L) {
  stop("Missing required Hallmark set(s): ", paste(missing_hallmarks, collapse = ", "))
}

module_pathway_map <- c(
  "Interferon alpha response" = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "TNFA/NFKB signalling" = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "Inflammatory response" = "HALLMARK_INFLAMMATORY_RESPONSE",
  "Gastric metaplasia" = "FINAL_SSL_METAPLASIA_LIST"
)
module_order <- names(module_pathway_map)
module_gene_sets <- list(
  "Interferon alpha response" = hallmark_pathways[[
    "HALLMARK_INTERFERON_ALPHA_RESPONSE"
  ]],
  "TNFA/NFKB signalling" = hallmark_pathways[[
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
  ]],
  "Inflammatory response" = hallmark_pathways[[
    "HALLMARK_INFLAMMATORY_RESPONSE"
  ]],
  "Gastric metaplasia" = metaplasia_genes
)

signature_manifest <- imap_dfr(
  module_gene_sets,
  function(genes, module_name) {
    tibble(
      Module = module_name,
      Pathway = module_pathway_map[[module_name]],
      N_genes = length(unique(genes)),
      Source = ifelse(
        module_name == "Gastric metaplasia",
        metaplasia_signature_file,
        "MSigDB Hallmark collection via msigdbr"
      ),
      GSEA_filter = ifelse(
        module_name == "Gastric metaplasia",
        "No expression filter",
        paste0("DE-table logCPM >= ", hallmark_logcpm_cutoff)
      )
    )
  }
)
fwrite(signature_manifest, file.path(subdirs$qc, "SIGNATURE_MANIFEST.csv"))
writeLines(
  metaplasia_genes,
  file.path(subdirs$qc, "FINAL_SSL_METAPLASIA_LIST_EXACT_24_USED.txt")
)

# DE-table import and gene-symbol mapping ----
map_gene_symbols <- function(gene_ids, supplied_symbols = NULL) {
  ids <- clean_character(gene_ids)
  ids_no_version <- toupper(str_remove(ids, "\\.[0-9]+$"))
  symbols <- if (is.null(supplied_symbols)) {
    rep(NA_character_, length(ids))
  } else {
    clean_character(supplied_symbols)
  }
  symbols[toupper(symbols) %in% c("NA", "N/A", ".", "NULL")] <- NA_character_

  is_ensembl <- !is.na(ids_no_version) & grepl("^ENSG[0-9]+$", ids_no_version)
  need_mapping <- is.na(symbols) & is_ensembl
  if (any(need_mapping)) {
    mapped <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = unique(ids_no_version[need_mapping]),
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )
    symbols[need_mapping] <- unname(mapped[ids_no_version[need_mapping]])
  }
  fallback <- is.na(symbols) & !is.na(ids_no_version) & !is_ensembl
  symbols[fallback] <- ids_no_version[fallback]
  symbols <- toupper(trimws(symbols))
  symbols[symbols == "" | symbols == "NA"] <- NA_character_
  symbols
}

as_numeric_strict <- function(x, col, role, file) {
  raw <- clean_character(x)
  out <- suppressWarnings(as.numeric(raw))
  bad <- !is.na(raw) & is.na(out)
  if (any(bad)) {
    stop(
      "Non-numeric ", role, " values in ", basename(file),
      " column '", col, "': ",
      paste(head(unique(raw[bad]), 8), collapse = "; ")
    )
  }
  out
}

standardise_de_table <- function(file, comparison) {
  de <- fread(
    file,
    data.table = FALSE,
    check.names = FALSE,
    fill = TRUE,
    na.strings = c("", "NA", "N/A", "na")
  )
  names(de) <- trimws(names(de))
  keep_cols <- vapply(
    de,
    function(x) any(!is.na(x) & trimws(as.character(x)) != ""),
    logical(1)
  )
  de <- de[, keep_cols, drop = FALSE]

  symbol_col <- detect_col(
    de,
    c("SYMBOL", "Symbol", "symbol", "gene_symbol", "GeneSymbol", "Gene"),
    FALSE
  )
  id_col <- detect_col(
    de,
    c("genes", "Gene_ID", "gene_id", "ENSEMBL", "Ensembl_ID", "Geneid", "ID"),
    FALSE
  )
  if (is.na(symbol_col) && is.na(id_col)) id_col <- names(de)[1]

  logfc_col <- detect_col(
    de,
    c("logFC", "LogFC", "logfc", "log2FC", "log2FoldChange", "LFC"),
    TRUE,
    paste0(comparison, " logFC")
  )
  logcpm_col <- detect_col(
    de,
    c("logCPM", "LogCPM", "logcpm", "log2CPM", "AveExpr", "average_expression"),
    FALSE
  )
  lr_col <- detect_col(
    de,
    c("LR", "LRT", "likelihood_ratio", "LikelihoodRatio"),
    FALSE
  )
  t_col <- detect_col(
    de,
    c("t", "T", "stat", "STAT", "Statistic", "moderated_t"),
    FALSE
  )
  p_col <- detect_col(
    de,
    c("PValue", "P.Value", "Pvalue", "pvalue", "p.value", "p", "P"),
    TRUE,
    paste0(comparison, " P-value")
  )
  fdr_col <- detect_col(
    de,
    c("FDR", "fdr", "padj", "adj.P.Val", "adj_p", "qvalue", "q.value"),
    FALSE
  )

  gene_ids <- if (!is.na(id_col)) de[[id_col]] else de[[symbol_col]]
  supplied_symbols <- if (!is.na(symbol_col)) de[[symbol_col]] else NULL

  out <- tibble(
    Gene = map_gene_symbols(gene_ids, supplied_symbols),
    logFC = as_numeric_strict(de[[logfc_col]], logfc_col, "logFC", file),
    logCPM = if (!is.na(logcpm_col)) {
      as_numeric_strict(de[[logcpm_col]], logcpm_col, "logCPM", file)
    } else NA_real_,
    LR = if (!is.na(lr_col)) {
      as_numeric_strict(de[[lr_col]], lr_col, "LR", file)
    } else NA_real_,
    t = if (!is.na(t_col)) {
      as_numeric_strict(de[[t_col]], t_col, "t statistic", file)
    } else NA_real_,
    PValue = as_numeric_strict(de[[p_col]], p_col, "PValue", file),
    FDR = if (!is.na(fdr_col)) {
      as_numeric_strict(de[[fdr_col]], fdr_col, "FDR", file)
    } else NA_real_,
    Comparison = comparison,
    Source_file = normalizePath(file, winslash = "/", mustWork = TRUE)
  )

  if (any(is.finite(out$LR) & out$LR < -1e-8)) {
    stop(comparison, " contains materially negative LR values. Do not use abs(LR).")
  }
  out$LR[is.finite(out$LR) & out$LR < 0] <- 0

  out <- out %>%
    filter(!is.na(Gene), Gene != "", is.finite(logFC), is.finite(PValue))
  if (nrow(out) == 0L) stop("No valid rows remained for ", comparison)
  if (all(is.na(out$FDR))) out$FDR <- p.adjust(out$PValue, method = "BH")

  positive_p <- out$PValue[is.finite(out$PValue) & out$PValue > 0]
  min_positive_p <- if (length(positive_p) > 0L) min(positive_p) else 1e-300
  p_safe <- ifelse(
    !is.finite(out$PValue) | out$PValue <= 0,
    min_positive_p * 0.1,
    out$PValue
  )

  out <- out %>%
    mutate(
      collapse_rank = case_when(
        is.finite(LR) ~ sign(logFC) * sqrt(LR),
        is.finite(t) ~ t,
        TRUE ~ sign(logFC) * -log10(p_safe)
      )
    ) %>%
    group_by(Gene) %>%
    arrange(desc(abs(collapse_rank)), PValue, desc(abs(logFC)), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(-collapse_rank)

  qc <- tibble(
    Comparison = comparison,
    Source_file = unique(out$Source_file)[1],
    N_genes = nrow(out),
    N_logCPM = sum(is.finite(out$logCPM)),
    N_LR = sum(is.finite(out$LR)),
    N_t = sum(is.finite(out$t)),
    Symbol_column = ifelse(is.na(symbol_col), NA_character_, symbol_col),
    ID_column = ifelse(is.na(id_col), NA_character_, id_col),
    logFC_column = logfc_col,
    logCPM_column = ifelse(is.na(logcpm_col), NA_character_, logcpm_col),
    LR_column = ifelse(is.na(lr_col), NA_character_, lr_col),
    PValue_column = p_col,
    FDR_column = ifelse(is.na(fdr_col), "BH calculated", fdr_col)
  )
  fwrite(qc, file.path(subdirs$qc, paste0(safe_filename(comparison), "_DE_QC.csv")))
  out
}

de_list <- list()
for (i in seq_len(nrow(active_manifest))) {
  comparison <- active_manifest$Comparison[i]
  message("Reading DE table: ", comparison)
  de_list[[comparison]] <- standardise_de_table(
    active_manifest$Path[i],
    comparison
  )
  fwrite(
    de_list[[comparison]],
    file.path(subdirs$tables, paste0(safe_filename(comparison), "_standardised_DE.csv"))
  )
}

# Ranked GSEA: Hallmarks filtered, gastric signature unfiltered ----
make_ranked_vector <- function(de, comparison, apply_logcpm_filter) {
  dat <- de
  if (apply_logcpm_filter) {
    if (!any(is.finite(dat$logCPM))) {
      stop(comparison, " has no usable logCPM column for matched Hallmark GSEA.")
    }
    dat <- dat %>%
      filter(is.finite(logCPM), logCPM >= hallmark_logcpm_cutoff)
  }
  if (nrow(dat) < 100L) stop("Too few ranked genes for ", comparison)

  if (sum(is.finite(dat$LR)) >= 100L) {
    dat <- dat %>% filter(is.finite(LR))
    dat$rank_metric <- sign(dat$logFC) * sqrt(pmax(dat$LR, 0))
    method <- "sign(logFC) * sqrt(LR)"
  } else if (sum(is.finite(dat$t)) >= 100L) {
    dat <- dat %>% filter(is.finite(t))
    dat$rank_metric <- dat$t
    method <- "signed t statistic"
  } else {
    positive_p <- dat$PValue[is.finite(dat$PValue) & dat$PValue > 0]
    min_positive_p <- if (length(positive_p) > 0L) min(positive_p) else 1e-300
    safe_p <- ifelse(
      !is.finite(dat$PValue) | dat$PValue <= 0,
      min_positive_p * 0.1,
      dat$PValue
    )
    dat$rank_metric <- sign(dat$logFC) * -log10(safe_p)
    method <- "sign(logFC) * -log10(PValue)"
  }

  dat <- dat %>%
    filter(is.finite(rank_metric), !is.na(Gene), Gene != "") %>%
    arrange(desc(rank_metric), Gene)
  ranks <- dat$rank_metric
  names(ranks) <- dat$Gene
  if (anyDuplicated(ranks)) {
    ranks <- ranks - seq_along(ranks) * 1e-12
    dat$rank_metric <- ranks
  }
  ranks <- sort(ranks, decreasing = TRUE)

  list(
    ranks = ranks,
    table = dat,
    method = method,
    filter_mode = ifelse(
      apply_logcpm_filter,
      paste0("LOGCPM_GE_", hallmark_logcpm_cutoff),
      "ALL_VALID_GENES"
    )
  )
}

normalise_fgsea <- function(fg, comparison, family, rank_obj) {
  fg <- as.data.frame(fg)
  if (!"padj" %in% names(fg)) fg$padj <- p.adjust(fg$pval, method = "BH")
  if (!"log2err" %in% names(fg)) fg$log2err <- NA_real_
  if (!"leadingEdge" %in% names(fg)) fg$leadingEdge <- vector("list", nrow(fg))
  fg %>%
    mutate(
      Comparison = comparison,
      Gene_set_family = family,
      Filter_mode = rank_obj$filter_mode,
      Ranking_method = rank_obj$method,
      N_ranked_genes = length(rank_obj$ranks),
      leadingEdge = vapply(
        leadingEdge,
        function(x) paste(unique(toupper(x)), collapse = ";"),
        character(1)
      )
    )
}

run_fgsea_for_comparison <- function(de, comparison) {
  hallmark_rank <- make_ranked_vector(de, comparison, TRUE)
  metaplasia_rank <- make_ranked_vector(de, comparison, FALSE)

  set.seed(1)
  hallmark_fg <- fgseaMultilevel(
    pathways = hallmark_pathways,
    stats = hallmark_rank$ranks,
    minSize = fgsea_min_size,
    maxSize = fgsea_max_size,
    eps = fgsea_eps,
    nproc = fgsea_nproc,
    nPermSimple = fgsea_nperm_simple,
    scoreType = "std"
  ) %>%
    normalise_fgsea(
      comparison,
      "HALLMARK_FULL_COLLECTION",
      hallmark_rank
    )

  target_hallmark <- hallmark_fg %>%
    filter(pathway %in% target_hallmark_names) %>%
    mutate(
      P_adjustment_scope = "BH FDR across the complete Hallmark collection"
    )
  if (nrow(target_hallmark) != 3L) {
    stop("Not all three target Hallmark sets were returned for ", comparison)
  }

  set.seed(1)
  custom_fg <- fgseaMultilevel(
    pathways = list(FINAL_SSL_METAPLASIA_LIST = metaplasia_genes),
    stats = metaplasia_rank$ranks,
    minSize = 3L,
    maxSize = fgsea_max_size,
    eps = fgsea_eps,
    nproc = fgsea_nproc,
    nPermSimple = fgsea_nperm_simple,
    scoreType = "std"
  ) %>%
    normalise_fgsea(
      comparison,
      "CUSTOM_FINAL_24_GASTRIC_METAPLASIA",
      metaplasia_rank
    ) %>%
    mutate(
      P_adjustment_scope = paste0(
        "One custom pathway; padj equals nominal p and is not ",
        "Hallmark-family FDR"
      )
    )

  list(
    hallmark_rank = hallmark_rank,
    metaplasia_rank = metaplasia_rank,
    hallmark_all = hallmark_fg,
    target_results = bind_rows(target_hallmark, custom_fg)
  )
}

fgsea_objects <- list()
for (comparison in active_comparisons) {
  message("Running matched GSEA: ", comparison)
  fgsea_objects[[comparison]] <- run_fgsea_for_comparison(
    de_list[[comparison]],
    comparison
  )
  fwrite(
    fgsea_objects[[comparison]]$hallmark_rank$table,
    file.path(subdirs$gsea, paste0(comparison, "_Hallmark_ranked_genes.csv"))
  )
  fwrite(
    fgsea_objects[[comparison]]$metaplasia_rank$table,
    file.path(subdirs$gsea, paste0(comparison, "_Metaplasia_ALL_VALID_GENES_ranked_genes.csv"))
  )
  fwrite(
    fgsea_objects[[comparison]]$hallmark_all,
    file.path(subdirs$gsea, paste0(comparison, "_FULL_HALLMARK_FGSEA.csv"))
  )
}

comparison_processing_status <- active_manifest %>%
  transmute(
    Comparison,
    Display,
    Short_label,
    Primary,
    File,
    Path,
    File_found,
    DE_table_loaded = Comparison %in% names(de_list),
    Four_pathway_GSEA_completed = Comparison %in% names(fgsea_objects),
    Included_in_all_pair_numerical_analysis = TRUE,
    Included_in_ALT_vs_each_individual_outputs =
      Comparison != "Alt_Infected" & write_alt_vs_all_individual_outputs,
    Included_in_every_pair_individual_outputs =
      write_all_pairwise_individual_outputs
  )

fwrite(
  comparison_processing_status,
  file.path(subdirs$qc, "COMPARISON_PROCESSING_STATUS.csv")
)

target_fgsea_results <- bind_rows(
  lapply(fgsea_objects, function(x) x[["target_results"]])
) %>%
  mutate(
    Module = names(module_pathway_map)[
      match(pathway, unname(module_pathway_map))
    ],
    Module = factor(Module, levels = module_order),
    Positive_enrichment = is.finite(NES) & NES > 0
  ) %>%
  left_join(
    active_manifest %>% select(Comparison, Display, Short_label),
    by = "Comparison"
  ) %>%
  arrange(factor(Comparison, levels = active_comparisons), Module)

fwrite(
  target_fgsea_results,
  file.path(subdirs$tables, "TARGET_FOUR_PATHWAYS_FGSEA_RESULTS.csv")
)

parse_leading_edge <- function(x) {
  if (length(x) == 0L || is.na(x) || x == "") return(character(0))
  unique(toupper(strsplit(x, ";", fixed = TRUE)[[1]]))
}

module_gene_long_list <- list()
for (comparison in active_comparisons) {
  de <- de_list[[comparison]]
  display <- active_manifest$Display[
    match(comparison, active_manifest$Comparison)
  ]
  short <- active_manifest$Short_label[
    match(comparison, active_manifest$Comparison)
  ]

  for (module_name in module_order) {
    genes <- unique(toupper(module_gene_sets[[module_name]]))
    gsea_row <- target_fgsea_results %>%
      filter(
        Comparison == comparison,
        as.character(Module) == module_name
      ) %>%
      slice(1)
    if (nrow(gsea_row) != 1L) {
      stop("Missing GSEA result for ", comparison, " / ", module_name)
    }

    rank_obj <- if (module_name == "Gastric metaplasia") {
      fgsea_objects[[comparison]]$metaplasia_rank
    } else {
      fgsea_objects[[comparison]]$hallmark_rank
    }
    leading_edge <- parse_leading_edge(gsea_row$leadingEdge[1])

    one <- tibble(
      Module = module_name,
      Pathway = module_pathway_map[[module_name]],
      Signature_order = seq_along(genes),
      Gene = genes
    ) %>%
      left_join(
        de %>% select(Gene, logFC, logCPM, LR, t, PValue, FDR),
        by = "Gene"
      ) %>%
      left_join(
        rank_obj$table %>% select(Gene, rank_metric),
        by = "Gene"
      ) %>%
      mutate(
        Comparison = comparison,
        Comparison_display = display,
        Comparison_short = short,
        Present_in_DE = is.finite(logFC),
        Eligible_for_GSEA = is.finite(rank_metric),
        Upregulated = Eligible_for_GSEA & is.finite(logFC) & logFC > 0,
        Significant_up =
          Eligible_for_GSEA & is.finite(FDR) & FDR <= fdr_cutoff & logFC > 0,
        Strong_significant_up = Significant_up & logFC >= strong_logfc_cutoff,
        In_leading_edge = Gene %in% leading_edge,
        Positive_leading_edge =
          In_leading_edge & is.finite(gsea_row$NES[1]) & gsea_row$NES[1] > 0 &
          Eligible_for_GSEA & is.finite(logFC) & logFC > 0,
        GSEA_ES = gsea_row$ES[1],
        GSEA_NES = gsea_row$NES[1],
        GSEA_pval = gsea_row$pval[1],
        GSEA_padj = gsea_row$padj[1],
        GSEA_filter_mode = gsea_row$Filter_mode[1],
        GSEA_ranking_method = gsea_row$Ranking_method[1]
      )

    module_gene_long_list[[paste(comparison, module_name, sep = "__")]] <- one
  }
}

module_gene_long <- bind_rows(module_gene_long_list) %>%
  mutate(
    Module = factor(Module, levels = module_order),
    Comparison = factor(Comparison, levels = active_comparisons)
  )

fwrite(
  module_gene_long,
  file.path(subdirs$tables, "GENE_LEVEL_FOUR_SIGNATURES_ACROSS_ALL_DE_COMPARISONS.csv")
)

leading_edge_long <- module_gene_long %>%
  filter(In_leading_edge) %>%
  select(
    Comparison, Comparison_display, Module, Pathway, Gene,
    logFC, logCPM, PValue, FDR,
    Eligible_for_GSEA, Positive_leading_edge,
    GSEA_NES, GSEA_pval, GSEA_padj
  )
fwrite(leading_edge_long, file.path(subdirs$tables, "LEADING_EDGE_GENES_LONG.csv"))

leading_edge_membership <- module_gene_long %>%
  filter(Positive_leading_edge) %>%
  transmute(Module, Gene, Comparison = as.character(Comparison), Value = 1L) %>%
  distinct() %>%
  pivot_wider(names_from = Comparison, values_from = Value, values_fill = 0L)
fwrite(
  leading_edge_membership,
  file.path(subdirs$tables, "POSITIVE_LEADING_EDGE_MEMBERSHIP_MATRIX.csv")
)

make_consensus_table <- function(data, comparisons, scope) {
  data %>%
    filter(as.character(Comparison) %in% comparisons) %>%
    mutate(Comparison = as.character(Comparison)) %>%
    group_by(Module, Pathway, Gene) %>%
    summarise(
      Scope = scope,
      N_comparisons_total = length(comparisons),
      N_present = sum(Present_in_DE, na.rm = TRUE),
      N_GSEA_eligible = sum(Eligible_for_GSEA, na.rm = TRUE),
      N_positive_logFC = sum(Upregulated, na.rm = TRUE),
      N_significant_up = sum(Significant_up, na.rm = TRUE),
      N_strong_significant_up = sum(Strong_significant_up, na.rm = TRUE),
      N_positive_leading_edge = sum(Positive_leading_edge, na.rm = TRUE),
      Median_logFC = ifelse(
        any(is.finite(logFC)),
        median(logFC, na.rm = TRUE),
        NA_real_
      ),
      Min_logFC = ifelse(any(is.finite(logFC)), min(logFC, na.rm = TRUE), NA_real_),
      Max_logFC = ifelse(any(is.finite(logFC)), max(logFC, na.rm = TRUE), NA_real_),
      Positive_in_every_eligible_comparison =
        sum(Eligible_for_GSEA, na.rm = TRUE) > 0L &&
        sum(Upregulated, na.rm = TRUE) == sum(Eligible_for_GSEA, na.rm = TRUE),
      Positive_leading_edge_comparisons = paste(
        Comparison[Positive_leading_edge],
        collapse = ";"
      ),
      Significant_up_comparisons = paste(
        Comparison[Significant_up],
        collapse = ";"
      ),
      Positive_logFC_comparisons = paste(
        Comparison[Upregulated],
        collapse = ";"
      ),
      .groups = "drop"
    ) %>%
    arrange(
      factor(Module, levels = module_order),
      desc(N_positive_leading_edge),
      desc(N_significant_up),
      desc(N_positive_logFC),
      desc(Median_logFC),
      Gene
    )
}

consensus_primary <- make_consensus_table(
  module_gene_long,
  primary_comparisons,
  "PRIMARY_FIVE_COMPARISONS"
)
consensus_all <- make_consensus_table(
  module_gene_long,
  active_comparisons,
  "ALL_AVAILABLE_COMPARISONS"
)
fwrite(
  consensus_primary,
  file.path(subdirs$tables, "CONSENSUS_DRIVER_GENES_PRIMARY_FIVE_COMPARISONS.csv")
)
fwrite(
  consensus_all,
  file.path(subdirs$tables, "CONSENSUS_DRIVER_GENES_ALL_AVAILABLE_COMPARISONS.csv")
)

# Pairwise gene concordance and leading-edge overlap ----
safe_cor_test <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3L || sd(x) == 0 || sd(y) == 0) {
    return(list(estimate = NA_real_, p.value = NA_real_))
  }
  out <- suppressWarnings(
    tryCatch(cor.test(x, y, method = method, exact = FALSE), error = function(e) NULL)
  )
  if (is.null(out)) return(list(estimate = NA_real_, p.value = NA_real_))
  list(estimate = unname(out$estimate), p.value = out$p.value)
}

jaccard <- function(a, b) {
  union_n <- length(union(unique(a), unique(b)))
  if (union_n == 0L) NA_real_ else length(intersect(unique(a), unique(b))) / union_n
}

overlap_hypergeometric_p <- function(universe_n, set_a_n, set_b_n, overlap_n) {
  values <- c(universe_n, set_a_n, set_b_n, overlap_n)
  if (any(!is.finite(values)) || universe_n <= 0L ||
      set_a_n < 0L || set_b_n < 0L || overlap_n < 0L ||
      set_a_n > universe_n || set_b_n > universe_n ||
      overlap_n > min(set_a_n, set_b_n)) {
    return(NA_real_)
  }
  if (set_a_n == 0L || set_b_n == 0L) return(NA_real_)
  stats::phyper(
    q = overlap_n - 1L,
    m = set_a_n,
    n = universe_n - set_a_n,
    k = set_b_n,
    lower.tail = FALSE
  )
}

analyse_pair_module <- function(comp_a, comp_b, module_name) {
  a <- module_gene_long %>%
    filter(
      as.character(Comparison) == comp_a,
      as.character(Module) == module_name
    ) %>%
    transmute(
      Module = as.character(Module), Pathway, Gene, Signature_order,
      Present_A = Present_in_DE, Eligible_A = Eligible_for_GSEA,
      logFC_A = logFC, logCPM_A = logCPM, PValue_A = PValue, FDR_A = FDR,
      Rank_A = rank_metric, Up_A = Upregulated,
      Significant_up_A = Significant_up,
      Strong_significant_up_A = Strong_significant_up,
      Leading_edge_A = In_leading_edge,
      Positive_leading_edge_A = Positive_leading_edge,
      NES_A = GSEA_NES, GSEA_p_A = GSEA_pval, GSEA_padj_A = GSEA_padj
    )

  b <- module_gene_long %>%
    filter(
      as.character(Comparison) == comp_b,
      as.character(Module) == module_name
    ) %>%
    transmute(
      Module = as.character(Module), Pathway, Gene, Signature_order,
      Present_B = Present_in_DE, Eligible_B = Eligible_for_GSEA,
      logFC_B = logFC, logCPM_B = logCPM, PValue_B = PValue, FDR_B = FDR,
      Rank_B = rank_metric, Up_B = Upregulated,
      Significant_up_B = Significant_up,
      Strong_significant_up_B = Strong_significant_up,
      Leading_edge_B = In_leading_edge,
      Positive_leading_edge_B = Positive_leading_edge,
      NES_B = GSEA_NES, GSEA_p_B = GSEA_pval, GSEA_padj_B = GSEA_padj
    )

  pair <- inner_join(
    a,
    b,
    by = c("Module", "Pathway", "Gene", "Signature_order")
  ) %>%
    mutate(
      Comparison_A = comp_a,
      Comparison_B = comp_b,
      Primary_pair = comp_a %in% primary_comparisons && comp_b %in% primary_comparisons,
      Present_both = Present_A & Present_B,
      Eligible_both = Eligible_A & Eligible_B,
      Positive_both = Eligible_both & Up_A & Up_B,
      Significant_up_both = Eligible_both & Significant_up_A & Significant_up_B,
      Strong_significant_up_both =
        Eligible_both & Strong_significant_up_A & Strong_significant_up_B,
      Shared_positive_leading_edge =
        Positive_leading_edge_A & Positive_leading_edge_B,
      Shared_directional_leading_edge =
        Eligible_both & Leading_edge_A & Leading_edge_B,
      Same_direction =
        Eligible_both & is.finite(logFC_A) & is.finite(logFC_B) &
        sign(logFC_A) == sign(logFC_B),
      Opposite_direction =
        Eligible_both & is.finite(logFC_A) & is.finite(logFC_B) &
        sign(logFC_A) != sign(logFC_B),
      Concordance_status = case_when(
        Shared_positive_leading_edge ~ "Shared positive leading edge",
        Shared_directional_leading_edge ~ "Shared directional leading edge",
        Significant_up_both ~ "Significant up in both",
        Positive_both ~ "Positive in both",
        Opposite_direction ~ "Opposite direction",
        Positive_leading_edge_A | Positive_leading_edge_B ~
          "Positive leading edge in one",
        TRUE ~ "Other"
      )
    )

  eligible <- pair %>%
    filter(Eligible_both, is.finite(logFC_A), is.finite(logFC_B))
  spearman <- safe_cor_test(eligible$logFC_A, eligible$logFC_B, "spearman")
  pearson <- safe_cor_test(eligible$logFC_A, eligible$logFC_B, "pearson")

  common_universe <- pair$Gene[pair$Eligible_both]
  le_a <- pair$Gene[pair$Positive_leading_edge_A & pair$Eligible_both]
  le_b <- pair$Gene[pair$Positive_leading_edge_B & pair$Eligible_both]
  directional_le_a <- pair$Gene[pair$Leading_edge_A & pair$Eligible_both]
  directional_le_b <- pair$Gene[pair$Leading_edge_B & pair$Eligible_both]
  sig_a <- pair$Gene[pair$Significant_up_A & pair$Eligible_both]
  sig_b <- pair$Gene[pair$Significant_up_B & pair$Eligible_both]
  shared_le <- intersect(unique(le_a), unique(le_b))
  shared_directional_le <- intersect(
    unique(directional_le_a),
    unique(directional_le_b)
  )
  shared_sig <- intersect(unique(sig_a), unique(sig_b))

  summary <- tibble(
    Comparison_A = comp_a,
    Comparison_B = comp_b,
    Primary_pair = comp_a %in% primary_comparisons && comp_b %in% primary_comparisons,
    Module = module_name,
    Pathway = unique(pair$Pathway)[1],
    Fixed_module_size = nrow(pair),
    Genes_present_A = sum(pair$Present_A, na.rm = TRUE),
    Genes_present_B = sum(pair$Present_B, na.rm = TRUE),
    Genes_present_both = sum(pair$Present_both, na.rm = TRUE),
    Genes_GSEA_eligible_both = length(unique(common_universe)),
    Positive_in_both_n = sum(pair$Positive_both, na.rm = TRUE),
    Significant_up_in_both_n = sum(pair$Significant_up_both, na.rm = TRUE),
    Strong_significant_up_in_both_n = sum(
      pair$Strong_significant_up_both,
      na.rm = TRUE
    ),
    Positive_leading_edge_A_n = length(unique(le_a)),
    Positive_leading_edge_B_n = length(unique(le_b)),
    Shared_positive_leading_edge_n = length(shared_le),
    Positive_leading_edge_Jaccard = jaccard(le_a, le_b),
    Positive_leading_edge_overlap_P_value = overlap_hypergeometric_p(
      length(unique(common_universe)),
      length(unique(le_a)),
      length(unique(le_b)),
      length(shared_le)
    ),
    Directional_leading_edge_A_n = length(unique(directional_le_a)),
    Directional_leading_edge_B_n = length(unique(directional_le_b)),
    Shared_directional_leading_edge_n = length(shared_directional_le),
    Directional_leading_edge_Jaccard = jaccard(
      directional_le_a,
      directional_le_b
    ),
    Directional_leading_edge_overlap_P_value = overlap_hypergeometric_p(
      length(unique(common_universe)),
      length(unique(directional_le_a)),
      length(unique(directional_le_b)),
      length(shared_directional_le)
    ),
    Significant_up_A_n = length(unique(sig_a)),
    Significant_up_B_n = length(unique(sig_b)),
    Significant_up_Jaccard = jaccard(sig_a, sig_b),
    Significant_up_overlap_P_value = overlap_hypergeometric_p(
      length(unique(common_universe)),
      length(unique(sig_a)),
      length(unique(sig_b)),
      length(shared_sig)
    ),
    Same_direction_fraction = ifelse(
      nrow(eligible) > 0L,
      mean(sign(eligible$logFC_A) == sign(eligible$logFC_B)),
      NA_real_
    ),
    Both_positive_fraction = ifelse(
      nrow(eligible) > 0L,
      mean(eligible$logFC_A > 0 & eligible$logFC_B > 0),
      NA_real_
    ),
    Spearman_rho_logFC = spearman$estimate,
    Spearman_P_value = spearman$p.value,
    Pearson_r_logFC = pearson$estimate,
    Pearson_P_value = pearson$p.value,
    NES_A = unique(pair$NES_A)[1],
    NES_B = unique(pair$NES_B)[1],
    GSEA_P_A = unique(pair$GSEA_p_A)[1],
    GSEA_P_B = unique(pair$GSEA_p_B)[1],
    GSEA_padj_A = unique(pair$GSEA_padj_A)[1],
    GSEA_padj_B = unique(pair$GSEA_padj_B)[1],
    Shared_positive_leading_edge_genes = paste(shared_le, collapse = ";"),
    Shared_directional_leading_edge_genes = paste(
      shared_directional_le,
      collapse = ";"
    ),
    Shared_significant_up_genes = paste(shared_sig, collapse = ";")
  )

  list(gene_table = pair, summary = summary)
}

# All comparison-module tests and multiple-testing correction ----
all_pairs <- combn(active_comparisons, 2, simplify = FALSE)
pair_gene_tables <- list()
pair_summaries <- list()
for (pair in all_pairs) {
  for (module_name in module_order) {
    key <- paste(pair[1], pair[2], module_name, sep = "__")
    result <- analyse_pair_module(pair[1], pair[2], module_name)
    pair_gene_tables[[key]] <- result$gene_table
    pair_summaries[[key]] <- result$summary
  }
}

pairwise_gene_long <- bind_rows(pair_gene_tables) %>%
  mutate(Module = factor(Module, levels = module_order))
pairwise_summary <- bind_rows(pair_summaries) %>%
  mutate(
    Module = factor(Module, levels = module_order),
    Spearman_FDR_across_all_pair_module_tests = stats::p.adjust(
      Spearman_P_value,
      method = "BH"
    ),
    Positive_leading_edge_overlap_FDR_across_all_pair_module_tests = stats::p.adjust(
      Positive_leading_edge_overlap_P_value,
      method = "BH"
    ),
    Directional_leading_edge_overlap_FDR_across_all_pair_module_tests = stats::p.adjust(
      Directional_leading_edge_overlap_P_value,
      method = "BH"
    ),
    Significant_up_overlap_FDR_across_all_pair_module_tests = stats::p.adjust(
      Significant_up_overlap_P_value,
      method = "BH"
    )
  ) %>%
  group_by(Module) %>%
  mutate(
    Spearman_FDR_within_module = stats::p.adjust(Spearman_P_value, method = "BH"),
    Positive_leading_edge_overlap_FDR_within_module = stats::p.adjust(
      Positive_leading_edge_overlap_P_value,
      method = "BH"
    ),
    Directional_leading_edge_overlap_FDR_within_module = stats::p.adjust(
      Directional_leading_edge_overlap_P_value,
      method = "BH"
    ),
    Significant_up_overlap_FDR_within_module = stats::p.adjust(
      Significant_up_overlap_P_value,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  left_join(
    active_manifest %>%
      select(Comparison_A = Comparison, Display_A = Display, Short_A = Short_label),
    by = "Comparison_A"
  ) %>%
  left_join(
    active_manifest %>%
      select(Comparison_B = Comparison, Display_B = Display, Short_B = Short_label),
    by = "Comparison_B"
  ) %>%
  arrange(
    desc(Primary_pair),
    factor(Comparison_A, levels = active_comparisons),
    factor(Comparison_B, levels = active_comparisons),
    Module
  )

fwrite(
  pairwise_gene_long,
  file.path(subdirs$tables, "PAIRWISE_GENE_LEVEL_CONCORDANCE_ALL_PAIRS.csv")
)
fwrite(
  pairwise_summary,
  file.path(subdirs$tables, "PAIRWISE_CONCORDANCE_SUMMARY_ALL_PAIRS.csv")
)

primary_pairs <- combn(primary_comparisons, 2, simplify = FALSE)
for (pair in primary_pairs) {
  pair_dir <- file.path(
    subdirs$primary_tables,
    paste0(safe_filename(pair[1]), "__VS__", safe_filename(pair[2]))
  )
  dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
  for (module_name in module_order) {
    one <- pairwise_gene_long %>%
      filter(
        Comparison_A == pair[1],
        Comparison_B == pair[2],
        as.character(Module) == module_name
      ) %>%
      arrange(
        desc(Shared_positive_leading_edge),
        desc(Shared_directional_leading_edge),
        desc(Significant_up_both),
        desc(Positive_both),
        desc(pmin(abs(logFC_A), abs(logFC_B), na.rm = TRUE)),
        Gene
      )
    fwrite(
      one,
      file.path(pair_dir, paste0(safe_filename(module_name), "_gene_concordance.csv"))
    )
  }
  fwrite(
    pairwise_summary %>%
      filter(Comparison_A == pair[1], Comparison_B == pair[2]),
    file.path(pair_dir, "PAIR_SUMMARY_ALL_FOUR_MODULES.csv")
  )
}

write_pair_gene_tables <- function(pair, root_directory) {
  pair_dir <- file.path(
    root_directory,
    paste0(safe_filename(pair[1]), "__VS__", safe_filename(pair[2]))
  )
  dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)

  for (module_name in module_order) {
    one <- pairwise_gene_long %>%
      filter(
        Comparison_A == pair[1],
        Comparison_B == pair[2],
        as.character(Module) == module_name
      ) %>%
      arrange(
        desc(Shared_positive_leading_edge),
        desc(Shared_directional_leading_edge),
        desc(Significant_up_both),
        desc(Positive_both),
        desc(pmin(abs(logFC_A), abs(logFC_B), na.rm = TRUE)),
        Gene
      )

    fwrite(
      one,
      file.path(
        pair_dir,
        paste0(safe_filename(module_name), "_gene_concordance.csv")
      )
    )
  }

  fwrite(
    pairwise_summary %>%
      filter(Comparison_A == pair[1], Comparison_B == pair[2]),
    file.path(pair_dir, "PAIR_SUMMARY_ALL_FOUR_MODULES.csv")
  )
}

alt_other_comparisons <- setdiff(active_comparisons, "Alt_Infected")
alt_pairs <- lapply(
  alt_other_comparisons,
  function(other_comparison) c("Alt_Infected", other_comparison)
)

if (write_alt_vs_all_individual_outputs) {
  for (pair in alt_pairs) {
    write_pair_gene_tables(pair, subdirs$alt_tables)
  }
}

if (write_all_pairwise_individual_outputs) {
  for (pair in all_pairs) {
    write_pair_gene_tables(pair, subdirs$all_pairwise_tables)
  }
}

# Concordance figures and per-comparison exports ----
concordance_colours <- c(
  "Shared positive leading edge" = "#B2182B",
  "Shared directional leading edge only" = "#762A83",
  "Other pathway genes" = "grey75"
)

format_metric_number <- function(x, digits = 2L) {
  ifelse(is.na(x), "NA", sprintf(paste0("%.", digits, "f"), x))
}

format_metric_p <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(
      x < 0.001,
      formatC(x, format = "e", digits = 1),
      sprintf("%.3f", x)
    )
  )
}

make_pair_scatter_data <- function(comp_a, comp_b) {
  pairwise_gene_long %>%
    filter(
      Comparison_A == comp_a,
      Comparison_B == comp_b,
      Eligible_both,
      is.finite(logFC_A),
      is.finite(logFC_B)
    ) %>%
    mutate(
      Module = factor(as.character(Module), levels = module_order),
      Concordance_status = case_when(
        Shared_positive_leading_edge ~
          "Shared positive leading edge",
        Shared_directional_leading_edge ~
          "Shared directional leading edge only",
        TRUE ~ "Other pathway genes"
      ),
      Concordance_status = factor(
        Concordance_status,
        levels = names(concordance_colours)
      ),
      Label_priority = case_when(
        Shared_positive_leading_edge ~ 2L,
        Shared_directional_leading_edge ~ 1L,
        TRUE ~ 0L
      ),
      Label_score = pmin(abs(logFC_A), abs(logFC_B))
    )
}

make_pair_annotation <- function(comp_a, comp_b) {
  pairwise_summary %>%
    filter(Comparison_A == comp_a, Comparison_B == comp_b) %>%
    mutate(
      Module = factor(as.character(Module), levels = module_order),
      Same_NES_direction =
        is.finite(NES_A) & is.finite(NES_B) &
        sign(NES_A) == sign(NES_B),
      Both_positive_NES =
        is.finite(NES_A) & is.finite(NES_B) &
        NES_A > 0 & NES_B > 0,
      Leading_edge_line = case_when(
        Both_positive_NES ~ paste0(
          "Shared positive LE: ",
          Shared_positive_leading_edge_n,
          " genes; Jaccard = ",
          format_metric_number(Positive_leading_edge_Jaccard)
        ),
        Same_NES_direction ~ paste0(
          "Shared directional LE: ",
          Shared_directional_leading_edge_n,
          " genes; Jaccard = ",
          format_metric_number(Directional_leading_edge_Jaccard)
        ),
        TRUE ~ paste0(
          "Opposite NES directions; shared LE: ",
          Shared_directional_leading_edge_n,
          " genes; Jaccard = ",
          format_metric_number(Directional_leading_edge_Jaccard)
        )
      ),
      Annotation = paste0(
        Leading_edge_line,
        "\nSpearman \u03c1 = ",
        format_metric_number(Spearman_rho_logFC),
        "; BH-FDR = ",
        format_metric_p(
          Spearman_FDR_across_all_pair_module_tests
        )
      )
    )
}

make_pair_scatter <- function(comp_a, comp_b) {
  dat <- make_pair_scatter_data(comp_a, comp_b)
  if (nrow(dat) == 0L) {
    return(make_empty_plot(
      "Gene-level pathway concordance",
      "No eligible genes"
    ))
  }

  labels <- dat %>%
    filter(Label_priority > 0L) %>%
    group_by(Module) %>%
    arrange(
      desc(Label_priority),
      desc(Label_score),
      FDR_A,
      FDR_B,
      .by_group = TRUE
    ) %>%
    slice_head(n = max_scatter_labels_per_module) %>%
    ungroup()

  annotations <- make_pair_annotation(comp_a, comp_b)

  display_a <- active_manifest$Display[
    match(comp_a, active_manifest$Comparison)
  ]
  display_b <- active_manifest$Display[
    match(comp_b, active_manifest$Comparison)
  ]
  comparison_subtitle <- stringr::str_wrap(
    paste0(display_a, " versus ", display_b),
    width = 88
  )

  ggplot(dat, aes(logFC_A, logFC_B)) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.3,
      colour = "grey45"
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.3,
      colour = "grey45"
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.3,
      colour = "grey55"
    ) +
    geom_point(
      aes(colour = Concordance_status),
      size = 1.45,
      alpha = 0.88
    ) +
    geom_text_repel(
      data = labels,
      aes(label = Gene),
      family = font_family,
      size = 5.1 / 2.845276,
      colour = "black",
      box.padding = 0.25,
      point.padding = 0.12,
      min.segment.length = 0,
      segment.colour = "grey55",
      segment.linewidth = 0.25,
      max.overlaps = Inf,
      seed = 1,
      show.legend = FALSE
    ) +
    geom_text(
      data = annotations,
      aes(x = -Inf, y = Inf, label = Annotation),
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = 1.10,
      family = font_family,
      size = 5.2 / 2.845276,
      lineheight = 0.94
    ) +
    facet_wrap(~Module, ncol = 2, scales = "free") +
    scale_colour_manual(
      values = concordance_colours,
      drop = TRUE
    ) +
    labs(
      title = "Gene-level pathway concordance",
      subtitle = comparison_subtitle,
      x = paste0(display_a, " log2FC"),
      y = paste0(display_b, " log2FC"),
      colour = NULL
    ) +
    nature_theme() +
    theme(
      plot.title.position = "plot",
      plot.subtitle = element_text(
        family = font_family,
        size = 6.2,
        hjust = 0,
        lineheight = 0.95,
        margin = margin(b = 2, unit = "mm")
      ),
      legend.position = "bottom",
      strip.background = element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.3
      ),
      strip.text = element_text(
        face = "bold",
        size = 6.2
      ),
      panel.spacing = unit(4, "mm"),
      plot.margin = margin(5, 5, 5, 5, unit = "mm")
    ) +
    guides(
      colour = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(size = 2, alpha = 1)
      )
    )
}

make_single_module_scatter <- function(
    comp_a,
    comp_b,
    module_name
) {
  dat <- make_pair_scatter_data(comp_a, comp_b) %>%
    filter(as.character(Module) == module_name)

  if (nrow(dat) == 0L) {
    return(make_empty_plot(
      module_name,
      "No eligible genes"
    ))
  }

  labels <- dat %>%
    filter(Label_priority > 0L) %>%
    arrange(
      desc(Label_priority),
      desc(Label_score),
      FDR_A,
      FDR_B
    ) %>%
    slice_head(n = max_scatter_labels_per_module)

  annotation <- make_pair_annotation(comp_a, comp_b) %>%
    filter(as.character(Module) == module_name)

  display_a <- active_manifest$Display[
    match(comp_a, active_manifest$Comparison)
  ]
  display_b <- active_manifest$Display[
    match(comp_b, active_manifest$Comparison)
  ]
  comparison_subtitle <- stringr::str_wrap(
    paste0(display_a, " versus ", display_b),
    width = 58
  )

  ggplot(dat, aes(logFC_A, logFC_B)) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.3,
      colour = "grey45"
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.3,
      colour = "grey45"
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.3,
      colour = "grey55"
    ) +
    geom_point(
      aes(colour = Concordance_status),
      size = 1.65,
      alpha = 0.90
    ) +
    geom_text_repel(
      data = labels,
      aes(label = Gene),
      family = font_family,
      size = 5.6 / 2.845276,
      colour = "black",
      box.padding = 0.28,
      point.padding = 0.14,
      min.segment.length = 0,
      segment.colour = "grey55",
      segment.linewidth = 0.25,
      max.overlaps = Inf,
      seed = 1,
      show.legend = FALSE
    ) +
    geom_text(
      data = annotation,
      aes(x = -Inf, y = Inf, label = Annotation),
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = 1.05,
      family = font_family,
      size = 5.6 / 2.845276,
      lineheight = 0.94
    ) +
    scale_colour_manual(
      values = concordance_colours,
      drop = TRUE
    ) +
    labs(
      title = module_name,
      subtitle = comparison_subtitle,
      x = paste0(display_a, " log2FC"),
      y = paste0(display_b, " log2FC"),
      colour = NULL
    ) +
    nature_theme() +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(
        family = font_family,
        face = "bold",
        size = 8,
        hjust = 0
      ),
      plot.subtitle = element_text(
        family = font_family,
        size = 6.2,
        hjust = 0,
        lineheight = 0.95,
        margin = margin(b = 2, unit = "mm")
      ),
      legend.position = "bottom",
      plot.margin = margin(5, 7, 5, 5, unit = "mm")
    ) +
    guides(
      colour = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(size = 2.2, alpha = 1)
      )
    )
}

single_module_stub <- function(module_name) {
  gsub(
    "_+",
    "_",
    gsub("[^A-Za-z0-9]+", "_", module_name)
  )
}

save_all_single_module_plots <- function(comp_a, comp_b) {
  for (module_name in module_order) {
    save_plot_set(
      make_single_module_scatter(
        comp_a,
        comp_b,
        module_name
      ),
      subdirs$all_pairwise_single_module_plots,
      paste0(
        comp_a,
        "__VS__",
        comp_b,
        "__",
        single_module_stub(module_name),
        "__SCATTER"
      ),
      135,
      110
    )
  }
}

for (pair in primary_pairs) {
  save_plot_set(
    make_pair_scatter(pair[1], pair[2]),
    subdirs$primary_plots,
    paste0(
      pair[1],
      "__VS__",
      pair[2],
      "__FOUR_MODULE_CONCORDANCE"
    ),
    183,
    150
  )
}

if (write_alt_vs_all_individual_outputs) {
  for (pair in alt_pairs) {
    save_plot_set(
      make_pair_scatter(pair[1], pair[2]),
      subdirs$alt_plots,
      paste0(
        pair[1],
        "__VS__",
        pair[2],
        "__FOUR_MODULE_CONCORDANCE"
      ),
      183,
      150
    )
  }
}

if (write_all_pairwise_individual_outputs) {
  for (pair in all_pairs) {
    save_plot_set(
      make_pair_scatter(pair[1], pair[2]),
      subdirs$all_pairwise_plots,
      paste0(
        pair[1],
        "__VS__",
        pair[2],
        "__FOUR_MODULE_CONCORDANCE"
      ),
      183,
      150
    )
    save_all_single_module_plots(pair[1], pair[2])
  }
}

write_chapter5_es_organoid_concordance_figure <- TRUE

if (write_chapter5_es_organoid_concordance_figure) {
  chapter5_concordance_dir <- file.path(
    out_dir,
    "09_Chapter5_ES_Organoid_Concordance"
  )
  dir.create(
    chapter5_concordance_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  chapter5_modules <- c(
    "Interferon alpha response",
    "TNFA/NFKB signalling",
    "Inflammatory response"
  )

  chapter5_module_labels <- c(
    "Interferon alpha response" = "Interferon-alpha response",
    "TNFA/NFKB signalling" = "TNF-alpha/NF-kappaB signalling",
    "Inflammatory response" = "Inflammatory response"
  )

  chapter5_tissue_labels <- c(
    "ES_SSL_HP" = "SSL_HP",
    "ES_H_ADJ_SSL_HP" = "H_ADJ_SSL_HP"
  )

  chapter5_required_comparisons <- c(
    "Alt_Infected",
    names(chapter5_tissue_labels)
  )

  chapter5_missing_comparisons <- setdiff(
    chapter5_required_comparisons,
    active_comparisons
  )

  if (length(chapter5_missing_comparisons) > 0L) {
    warning(
      "Chapter 5 Escherichia-Shigella concordance figure skipped because ",
      "comparison(s) were missing: ",
      paste(chapter5_missing_comparisons, collapse = ", ")
    )
  } else {

    get_chapter5_pair <- function(tissue_comparison, module_name) {
      forward_data <- pairwise_gene_long %>%
        filter(
          Comparison_A == "Alt_Infected",
          Comparison_B == tissue_comparison,
          as.character(Module) == module_name,
          Eligible_both,
          is.finite(logFC_A),
          is.finite(logFC_B)
        )

      if (nrow(forward_data) > 0L) {
        plot_data <- forward_data %>%
          transmute(
            Gene,
            Organoid_logFC = logFC_A,
            Tissue_logFC = logFC_B,
            Shared_positive_leading_edge,
            Label_score = pmin(abs(logFC_A), abs(logFC_B))
          )

        summary_data <- pairwise_summary %>%
          filter(
            Comparison_A == "Alt_Infected",
            Comparison_B == tissue_comparison,
            as.character(Module) == module_name
          ) %>%
          slice(1)
      } else {
        reverse_data <- pairwise_gene_long %>%
          filter(
            Comparison_A == tissue_comparison,
            Comparison_B == "Alt_Infected",
            as.character(Module) == module_name,
            Eligible_both,
            is.finite(logFC_A),
            is.finite(logFC_B)
          )

        plot_data <- reverse_data %>%
          transmute(
            Gene,
            Organoid_logFC = logFC_B,
            Tissue_logFC = logFC_A,
            Shared_positive_leading_edge,
            Label_score = pmin(abs(logFC_A), abs(logFC_B))
          )

        summary_data <- pairwise_summary %>%
          filter(
            Comparison_A == tissue_comparison,
            Comparison_B == "Alt_Infected",
            as.character(Module) == module_name
          ) %>%
          slice(1)
      }

      if (nrow(plot_data) == 0L || nrow(summary_data) != 1L) {
        stop(
          "Could not obtain the required pairwise data for Alt_Infected versus ",
          tissue_comparison,
          " / ",
          module_name,
          "."
        )
      }

      list(data = plot_data, summary = summary_data)
    }

    chapter5_point_colours <- c(
      "Shared positive leading-edge genes" = "#B2182B",
      "Other pathway genes" = "grey75"
    )

    make_chapter5_concordance_panel <- function(
        tissue_comparison,
        module_name,
        show_x_title = FALSE,
        show_y_title = FALSE
    ) {
      pair_object <- get_chapter5_pair(
        tissue_comparison,
        module_name
      )
      dat <- pair_object$data %>%
        mutate(
          Point_group = ifelse(
            Shared_positive_leading_edge,
            "Shared positive leading-edge genes",
            "Other pathway genes"
          ),
          Point_group = factor(
            Point_group,
            levels = names(chapter5_point_colours)
          )
        )

      summary_row <- pair_object$summary

      labels <- dat %>%
        filter(Shared_positive_leading_edge) %>%
        arrange(desc(Label_score), Gene) %>%
        slice_head(n = 4L)

      statistics_subtitle <- paste0(
        "Shared positive LE: ",
        summary_row$Shared_positive_leading_edge_n,
        " genes; Jaccard = ",
        format_metric_number(
          summary_row$Positive_leading_edge_Jaccard,
          digits = 2L
        ),
        "\nSpearman rho = ",
        format_metric_number(
          summary_row$Spearman_rho_logFC,
          digits = 2L
        ),
        "; global BH-FDR = ",
        format_metric_p(
          summary_row$Spearman_FDR_across_all_pair_module_tests
        )
      )

      tissue_label <- unname(
        chapter5_tissue_labels[tissue_comparison]
      )
      panel_title <- paste0(
        unname(chapter5_module_labels[module_name]),
        " - ",
        tissue_label
      )

      x_title <- if (show_x_title) {
        "IHE3034-infected versus non-infected organoids log2FC"
      } else {
        NULL
      }

      y_title <- if (show_y_title) {
        paste0(
          "Escherichia-Shigella-positive versus negative\n",
          tissue_label,
          " log2FC"
        )
      } else {
        NULL
      }

      ggplot(dat, aes(Organoid_logFC, Tissue_logFC)) +
        geom_hline(
          yintercept = 0,
          linewidth = 0.30,
          colour = "grey45"
        ) +
        geom_vline(
          xintercept = 0,
          linewidth = 0.30,
          colour = "grey45"
        ) +
        geom_abline(
          slope = 1,
          intercept = 0,
          linetype = "dashed",
          linewidth = 0.30,
          colour = "grey55"
        ) +
        geom_point(
          aes(colour = Point_group),
          size = 1.55,
          alpha = 0.90
        ) +
        ggrepel::geom_text_repel(
          data = labels,
          aes(label = Gene),
          family = font_family,
          size = 4.8 / 2.845276,
          colour = "black",
          box.padding = 0.24,
          point.padding = 0.10,
          min.segment.length = 0,
          segment.colour = "grey55",
          segment.linewidth = 0.22,
          max.overlaps = Inf,
          seed = 1,
          show.legend = FALSE
        ) +
        scale_colour_manual(
          values = chapter5_point_colours,
          drop = FALSE
        ) +
        scale_x_continuous(
          expand = expansion(mult = c(0.06, 0.12))
        ) +
        scale_y_continuous(
          expand = expansion(mult = c(0.08, 0.14))
        ) +
        labs(
          title = panel_title,
          subtitle = statistics_subtitle,
          x = x_title,
          y = y_title,
          colour = NULL
        ) +
        nature_theme(base_size = 7) +
        theme(
          plot.title = element_text(
            face = "bold",
            size = 7.1,
            hjust = 0
          ),
          plot.subtitle = element_text(
            size = 5.2,
            hjust = 0,
            lineheight = 0.96,
            margin = margin(b = 1.5, unit = "mm")
          ),
          axis.title.x = if (show_x_title) {
            element_text(size = 6.0)
          } else {
            element_blank()
          },
          axis.title.y = if (show_y_title) {
            element_text(size = 6.0)
          } else {
            element_blank()
          },
          legend.position = "bottom",
          plot.margin = margin(3.5, 4.5, 3.5, 4.5, unit = "mm")
        ) +
        guides(
          colour = guide_legend(
            nrow = 1,
            byrow = TRUE,
            override.aes = list(size = 2.2, alpha = 1)
          )
        )
    }

    p_ch5_a <- make_chapter5_concordance_panel(
      "ES_SSL_HP",
      "Interferon alpha response",
      show_x_title = FALSE,
      show_y_title = TRUE
    )
    p_ch5_b <- make_chapter5_concordance_panel(
      "ES_H_ADJ_SSL_HP",
      "Interferon alpha response",
      show_x_title = FALSE,
      show_y_title = FALSE
    )
    p_ch5_c <- make_chapter5_concordance_panel(
      "ES_SSL_HP",
      "TNFA/NFKB signalling",
      show_x_title = FALSE,
      show_y_title = TRUE
    )
    p_ch5_d <- make_chapter5_concordance_panel(
      "ES_H_ADJ_SSL_HP",
      "TNFA/NFKB signalling",
      show_x_title = FALSE,
      show_y_title = FALSE
    )
    p_ch5_e <- make_chapter5_concordance_panel(
      "ES_SSL_HP",
      "Inflammatory response",
      show_x_title = TRUE,
      show_y_title = TRUE
    )
    p_ch5_f <- make_chapter5_concordance_panel(
      "ES_H_ADJ_SSL_HP",
      "Inflammatory response",
      show_x_title = TRUE,
      show_y_title = FALSE
    )

    chapter5_concordance_figure <-
      (p_ch5_a | p_ch5_b) /
      (p_ch5_c | p_ch5_d) /
      (p_ch5_e | p_ch5_f) +
      patchwork::plot_layout(
        guides = "collect",
        heights = c(1, 1, 1)
      ) +
      patchwork::plot_annotation(
        tag_levels = "A",
        title = paste0(
          "Concordance between the IHE3034-induced epithelial response and ",
          "the Escherichia-Shigella-associated tissue response"
        ),
        subtitle = paste0(
          "Left column: SSL_HP. Right column: H_ADJ_SSL_HP. ",
          "Gastric metaplasia is analysed separately in Section 5.5."
        ),
        theme = theme(
          plot.title = element_text(
            family = font_family,
            face = "bold",
            size = 9,
            hjust = 0
          ),
          plot.subtitle = element_text(
            family = font_family,
            size = 6.2,
            hjust = 0,
            margin = margin(b = 2, unit = "mm")
          ),
          plot.tag = element_text(
            family = font_family,
            face = "bold",
            size = 7.5
          )
        )
      ) &
      theme(legend.position = "bottom")

    save_plot_set(
      chapter5_concordance_figure,
      chapter5_concordance_dir,
      "FIGURE_5_4_4_IHE3034_ES_IMMUNE_PATHWAY_CONCORDANCE",
      180,
      225
    )

    chapter5_concordance_metrics <- pairwise_summary %>%
      filter(
        as.character(Module) %in% chapter5_modules,
        (
          Comparison_A == "Alt_Infected" &
            Comparison_B %in% names(chapter5_tissue_labels)
        ) |
          (
            Comparison_B == "Alt_Infected" &
              Comparison_A %in% names(chapter5_tissue_labels)
          )
      ) %>%
      mutate(
        Tissue_comparison = ifelse(
          Comparison_A == "Alt_Infected",
          Comparison_B,
          Comparison_A
        ),
        Tissue_group = unname(
          chapter5_tissue_labels[Tissue_comparison]
        )
      ) %>%
      arrange(
        factor(Module, levels = chapter5_modules),
        factor(
          Tissue_comparison,
          levels = names(chapter5_tissue_labels)
        )
      )

    fwrite(
      chapter5_concordance_metrics,
      file.path(
        chapter5_concordance_dir,
        "FIGURE_5_4_4_CONCORDANCE_METRICS.csv"
      )
    )
  }
}

# Driver-gene heatmaps and pairwise summaries ----
select_driver_genes <- function(consensus, max_per_module) {
  selected <- consensus %>%
    filter(N_positive_leading_edge >= 1L | N_significant_up >= 2L) %>%
    group_by(Module) %>%
    arrange(
      desc(N_positive_leading_edge),
      desc(N_significant_up),
      desc(N_positive_logFC),
      desc(Median_logFC),
      Gene,
      .by_group = TRUE
    ) %>%
    slice_head(n = max_per_module) %>%
    ungroup()

  missing_modules <- setdiff(module_order, as.character(unique(selected$Module)))
  if (length(missing_modules) > 0L) {
    fallback <- consensus %>%
      filter(as.character(Module) %in% missing_modules) %>%
      group_by(Module) %>%
      arrange(desc(N_positive_logFC), desc(Median_logFC), Gene, .by_group = TRUE) %>%
      slice_head(n = min(15L, max_per_module)) %>%
      ungroup()
    selected <- bind_rows(selected, fallback)
  }

  selected %>%
    mutate(Module = factor(as.character(Module), levels = module_order)) %>%
    arrange(
      Module,
      desc(N_positive_leading_edge),
      desc(N_significant_up),
      desc(N_positive_logFC),
      desc(Median_logFC),
      Gene
    )
}

make_driver_heatmap <- function(comparisons, consensus, title_text) {
  selected <- select_driver_genes(consensus, max_driver_genes_per_module)
  plot_df <- module_gene_long %>%
    filter(as.character(Comparison) %in% comparisons) %>%
    mutate(Comparison = as.character(Comparison)) %>%
    inner_join(selected %>% select(Module, Gene), by = c("Module", "Gene"))

  if (nrow(plot_df) == 0L) {
    return(list(plot = make_empty_plot(title_text, "No selected genes"), selected = selected, data = plot_df))
  }

  labels <- active_manifest$Short_label[match(comparisons, active_manifest$Comparison)]
  names(labels) <- comparisons

  gene_levels <- selected %>%
    mutate(
      Module = factor(as.character(Module), levels = module_order),
      Gene_key = paste(as.character(Module), Gene, sep = "|||")
    ) %>%
    arrange(
      Module,
      desc(N_positive_leading_edge),
      desc(N_significant_up),
      desc(N_positive_logFC),
      desc(Median_logFC),
      Gene
    ) %>%
    pull(Gene_key)

  plot_df <- plot_df %>%
    mutate(
      Module = factor(as.character(Module), levels = module_order),
      Comparison_display = factor(
        unname(labels[Comparison]),
        levels = unname(labels[comparisons])
      ),
      Gene_key = factor(
        paste(as.character(Module), Gene, sep = "|||"),
        levels = rev(gene_levels)
      ),
      logFC_plot = pmax(pmin(logFC, heatmap_logfc_limit), -heatmap_logfc_limit)
    )

  p <- ggplot(plot_df, aes(Comparison_display, Gene_key, fill = logFC_plot)) +
    geom_tile(colour = "white", linewidth = 0.20) +
    geom_point(
      data = plot_df %>% filter(Positive_leading_edge),
      shape = 21,
      size = 1.2,
      stroke = 0.2,
      fill = "black",
      colour = "black",
      show.legend = FALSE
    ) +
    facet_grid(
      rows = vars(Module),
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-heatmap_logfc_limit, heatmap_logfc_limit),
      oob = squish,
      na.value = "grey92",
      name = "log2FC"
    ) +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|\\|", "", x)) +
    labs(
      title = title_text,
      subtitle = paste0(
        "Union of positive leading-edge genes and recurrent FDR-up genes; ",
        "black dot = positive leading edge"
      ),
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_family = font_family, base_size = 6) +
    theme(
      text = element_text(family = font_family, colour = "black"),
      plot.title = element_text(face = "bold", size = 8, hjust = 0),
      plot.subtitle = element_text(size = 5.5, hjust = 0),
      panel.grid = element_blank(),
      panel.spacing.y = unit(1.4, "mm"),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.25),
      strip.text.y.left = element_text(face = "bold", size = 6, angle = 0),
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 5.2),
      axis.text.y = element_text(size = 5.0),
      axis.ticks = element_blank(),
      legend.position = "right",
      plot.margin = margin(4, 4, 4, 4, "mm")
    )

  list(plot = p, selected = selected, data = plot_df)
}

primary_heatmap <- make_driver_heatmap(
  primary_comparisons,
  consensus_primary,
  "Concordant pathway-driver genes across the five primary comparisons"
)
all_heatmap <- make_driver_heatmap(
  active_comparisons,
  consensus_all,
  "Concordant pathway-driver genes across all available comparisons"
)

save_plot_set(
  primary_heatmap$plot,
  subdirs$heatmaps,
  "PRIMARY_FIVE_COMPARISONS_driver_gene_log2FC_heatmap",
  183,
  max(150, min(300, 75 + 3.0 * nrow(primary_heatmap$selected)))
)
save_plot_set(
  all_heatmap$plot,
  subdirs$heatmaps,
  "ALL_AVAILABLE_COMPARISONS_driver_gene_log2FC_heatmap",
  max(183, 24 * length(active_comparisons)),
  max(170, min(330, 75 + 2.6 * nrow(all_heatmap$selected)))
)
fwrite(
  primary_heatmap$selected,
  file.path(subdirs$tables, "PRIMARY_driver_genes_selected_for_heatmap.csv")
)
fwrite(
  primary_heatmap$data,
  file.path(subdirs$tables, "PRIMARY_driver_gene_heatmap_long_table.csv")
)
fwrite(
  all_heatmap$selected,
  file.path(subdirs$tables, "ALL_AVAILABLE_driver_genes_selected_for_heatmap.csv")
)
fwrite(
  all_heatmap$data,
  file.path(subdirs$tables, "ALL_AVAILABLE_driver_gene_heatmap_long_table.csv")
)

make_symmetric_metric <- function(metric, comparisons = primary_comparisons) {
  dat <- pairwise_summary %>%
    filter(
      Comparison_A %in% comparisons,
      Comparison_B %in% comparisons
    ) %>%
    transmute(
      Module = as.character(Module),
      Comparison_A,
      Comparison_B,
      Value = .data[[metric]]
    )
  swapped <- dat %>%
    transmute(
      Module,
      Comparison_A = Comparison_B,
      Comparison_B = Comparison_A,
      Value
    )
  diagonal <- expand_grid(
    Module = module_order,
    Comparison_A = comparisons
  ) %>%
    mutate(Comparison_B = Comparison_A, Value = 1)
  bind_rows(dat, swapped, diagonal) %>%
    distinct(Module, Comparison_A, Comparison_B, .keep_all = TRUE)
}

make_metric_heatmap <- function(
    metric,
    title_text,
    mode = c("rho", "jaccard"),
    comparisons = primary_comparisons
) {
  mode <- match.arg(mode)
  dat <- make_symmetric_metric(metric, comparisons)
  short_map <- active_manifest$Short_label
  names(short_map) <- active_manifest$Comparison

  dat <- dat %>%
    mutate(
      Module = factor(Module, levels = module_order),
      A_label = factor(
        unname(short_map[Comparison_A]),
        levels = unname(short_map[comparisons])
      ),
      B_label = factor(
        unname(short_map[Comparison_B]),
        levels = rev(unname(short_map[comparisons]))
      ),
      Text = ifelse(is.na(Value), "", sprintf("%.2f", Value))
    )

  p <- ggplot(dat, aes(A_label, B_label, fill = Value)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label = Text), family = font_family, size = 4.8 / 2.845276) +
    facet_wrap(~Module, ncol = 2) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_bw(base_family = font_family, base_size = 6) +
    theme(
      text = element_text(family = font_family, colour = "black"),
      plot.title = element_text(face = "bold", size = 8, hjust = 0),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.25),
      strip.text = element_text(face = "bold", size = 6),
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 5),
      axis.text.y = element_text(size = 5),
      axis.ticks = element_blank(),
      legend.position = "right",
      panel.spacing = unit(3, "mm")
    )

  if (mode == "rho") {
    p + scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      oob = squish,
      na.value = "grey92",
      name = "Spearman rho"
    )
  } else {
    p + scale_fill_gradient(
      low = "white",
      high = "#B2182B",
      limits = c(0, 1),
      oob = squish,
      na.value = "grey92",
      name = "Jaccard"
    )
  }
}

save_plot_set(
  make_metric_heatmap(
    "Spearman_rho_logFC",
    "Pairwise log2FC concordance within each fixed pathway",
    "rho"
  ),
  subdirs$summary_plots,
  "PRIMARY_pairwise_Spearman_log2FC_concordance_heatmaps",
  183,
  155
)
save_plot_set(
  make_metric_heatmap(
    "Positive_leading_edge_Jaccard",
    "Pairwise overlap of positive GSEA leading-edge genes",
    "jaccard"
  ),
  subdirs$summary_plots,
  "PRIMARY_pairwise_positive_leading_edge_Jaccard_heatmaps",
  183,
  155
)
save_plot_set(
  make_metric_heatmap(
    "Directional_leading_edge_Jaccard",
    "Pairwise overlap of directional GSEA leading-edge genes",
    "jaccard"
  ),
  subdirs$summary_plots,
  "PRIMARY_pairwise_directional_leading_edge_Jaccard_heatmaps",
  183,
  155
)

all_pairwise_heatmap_width <- max(250, 22 * length(active_comparisons))
all_pairwise_heatmap_height <- max(220, 17 * length(active_comparisons))

save_plot_set(
  make_metric_heatmap(
    "Spearman_rho_logFC",
    "Pairwise log2FC concordance across all available DE comparisons",
    "rho",
    comparisons = active_comparisons
  ),
  subdirs$summary_plots,
  "ALL_AVAILABLE_pairwise_Spearman_log2FC_concordance_heatmaps",
  all_pairwise_heatmap_width,
  all_pairwise_heatmap_height
)

save_plot_set(
  make_metric_heatmap(
    "Positive_leading_edge_Jaccard",
    "Pairwise positive leading-edge overlap across all available DE comparisons",
    "jaccard",
    comparisons = active_comparisons
  ),
  subdirs$summary_plots,
  "ALL_AVAILABLE_pairwise_positive_leading_edge_Jaccard_heatmaps",
  all_pairwise_heatmap_width,
  all_pairwise_heatmap_height
)

save_plot_set(
  make_metric_heatmap(
    "Directional_leading_edge_Jaccard",
    "Pairwise directional leading-edge overlap across all available DE comparisons",
    "jaccard",
    comparisons = active_comparisons
  ),
  subdirs$summary_plots,
  "ALL_AVAILABLE_pairwise_directional_leading_edge_Jaccard_heatmaps",
  all_pairwise_heatmap_width,
  all_pairwise_heatmap_height
)

alt_vs_each <- pairwise_summary %>%
  filter(
    (Comparison_A == "Alt_Infected" &
       Comparison_B %in% setdiff(active_comparisons, "Alt_Infected")) |
      (Comparison_B == "Alt_Infected" &
         Comparison_A %in% setdiff(active_comparisons, "Alt_Infected"))
  ) %>%
  mutate(
    Other_comparison = ifelse(
      Comparison_A == "Alt_Infected",
      Comparison_B,
      Comparison_A
    ),
    Other_display = active_manifest$Display[
      match(Other_comparison, active_manifest$Comparison)
    ]
  ) %>%
  arrange(
    factor(Other_comparison, levels = setdiff(active_comparisons, "Alt_Infected")),
    Module
  )
fwrite(
  alt_vs_each,
  file.path(subdirs$tables, "INFECTION_VS_EACH_AVAILABLE_COMPARISON_SUMMARY.csv")
)

alt_short_map <- active_manifest$Short_label
names(alt_short_map) <- active_manifest$Comparison
alt_plot_data <- alt_vs_each %>%
  mutate(
    Module = factor(as.character(Module), levels = module_order),
    Other_label = factor(
      unname(alt_short_map[Other_comparison]),
      levels = unname(
        alt_short_map[setdiff(active_comparisons, "Alt_Infected")]
      )
    ),
    Shared_LE_label = ifelse(
      is.na(Positive_leading_edge_Jaccard),
      paste0("n = ", Shared_positive_leading_edge_n, "\nJ = NA"),
      paste0(
        "n = ", Shared_positive_leading_edge_n,
        "\nJ = ", sprintf("%.2f", Positive_leading_edge_Jaccard)
      )
    ),
    Directional_LE_label = ifelse(
      is.na(Directional_leading_edge_Jaccard),
      paste0("n = ", Shared_directional_leading_edge_n, "\nJ = NA"),
      paste0(
        "n = ", Shared_directional_leading_edge_n,
        "\nJ = ", sprintf("%.2f", Directional_leading_edge_Jaccard)
      )
    ),
    Rho_label = ifelse(
      is.na(Spearman_rho_logFC),
      "rho = NA",
      paste0("rho = ", sprintf("%.2f", Spearman_rho_logFC))
    )
  )

p_alt_directional_leading_edge <- ggplot(
  alt_plot_data,
  aes(x = Other_label, y = Module, fill = Directional_leading_edge_Jaccard)
) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(
    aes(label = Directional_LE_label),
    family = font_family,
    size = 4.7 / 2.845276,
    lineheight = 0.90
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#7B3294",
    limits = c(0, 1),
    oob = squish,
    na.value = "grey92",
    name = "Directional LE\nJaccard"
  ) +
  labs(
    title = "Directional leading-edge overlap with Infected vs Non-Infected Organoids",
    subtitle = "Same GSEA driver genes are counted even when pathway enrichment is in the opposite direction",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_family = font_family, base_size = 6) +
  theme(
    text = element_text(family = font_family, colour = "black"),
    plot.title = element_text(face = "bold", size = 8, hjust = 0),
    plot.subtitle = element_text(size = 5.5, hjust = 0),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 5.2),
    axis.text.y = element_text(size = 5.2),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

p_alt_leading_edge <- ggplot(
  alt_plot_data,
  aes(x = Other_label, y = Module, fill = Positive_leading_edge_Jaccard)
) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(
    aes(label = Shared_LE_label),
    family = font_family,
    size = 4.7 / 2.845276,
    lineheight = 0.90
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#B2182B",
    limits = c(0, 1),
    oob = squish,
    na.value = "grey92",
    name = "Leading-edge\nJaccard"
  ) +
  labs(
    title = "Infected vs Non-Infected Organoids versus every available DE comparison",
    subtitle = "Cell text: number of shared positive leading-edge genes and Jaccard overlap",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_family = font_family, base_size = 6) +
  theme(
    text = element_text(family = font_family, colour = "black"),
    plot.title = element_text(face = "bold", size = 8, hjust = 0),
    plot.subtitle = element_text(size = 5.5, hjust = 0),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 5.2),
    axis.text.y = element_text(size = 5.2),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

p_alt_logfc <- ggplot(
  alt_plot_data,
  aes(x = Other_label, y = Module, fill = Spearman_rho_logFC)
) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(
    aes(label = Rho_label),
    family = font_family,
    size = 4.7 / 2.845276
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    oob = squish,
    na.value = "grey92",
    name = "Spearman rho"
  ) +
  labs(
    title = "Pathway-restricted log2FC concordance with Infected vs Non-Infected Organoids",
    subtitle = "Spearman correlation across genes represented in both matched rank universes",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_family = font_family, base_size = 6) +
  theme(
    text = element_text(family = font_family, colour = "black"),
    plot.title = element_text(face = "bold", size = 8, hjust = 0),
    plot.subtitle = element_text(size = 5.5, hjust = 0),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 5.2),
    axis.text.y = element_text(size = 5.2),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

alt_summary_width <- max(183, 28 * length(alt_other_comparisons))

save_plot_set(
  p_alt_directional_leading_edge,
  subdirs$summary_plots,
  "INFECTION_vs_each_directional_leading_edge_overlap",
  alt_summary_width,
  100
)
save_plot_set(
  p_alt_leading_edge,
  subdirs$summary_plots,
  "INFECTION_vs_each_positive_leading_edge_overlap",
  alt_summary_width,
  100
)
save_plot_set(
  p_alt_logfc,
  subdirs$summary_plots,
  "INFECTION_vs_each_pathway_log2FC_Spearman_concordance",
  alt_summary_width,
  100
)
save_plot_set(
  p_alt_directional_leading_edge /
    p_alt_leading_edge /
    p_alt_logfc +
    patchwork::plot_annotation(tag_levels = "a"),
  subdirs$summary_plots,
  "INFECTION_vs_each_COMBINED_CONCORDANCE_SUMMARY",
  alt_summary_width,
  285
)

if (write_focused_ssl_infection_braf_outputs) {
  focused_comparisons <- c(
    "SSL_HP_v_H_COL",
    "Alt_Infected",
    "Initial_BRAF_SSL_HP_organoid"
  )

  focused_comparison_labels <- c(
    "SSL_HP_v_H_COL" = "SSL_HP vs H_COL",
    "Alt_Infected" = "Infected vs Non-Infected Organoids",
    "Initial_BRAF_SSL_HP_organoid" = "BRAF V600E vs WT Organoids"
  )

  focused_comparison_axis_labels <- c(
    "SSL_HP_v_H_COL" = "SSL_HP vs H_COL",
    "Alt_Infected" = "Infected vs Non-Infected\nOrganoids",
    "Initial_BRAF_SSL_HP_organoid" = "BRAF V600E vs WT\nOrganoids"
  )

  focused_modules <- c(
    "Interferon alpha response",
    "TNFA/NFKB signalling",
    "Inflammatory response"
  )

  focused_module_labels <- c(
    "Interferon alpha response" = "Interferon-α response",
    "TNFA/NFKB signalling" = "TNF-α/NF-κB signalling",
    "Inflammatory response" = "Inflammatory response"
  )

  focused_fdr_stars <- function(x) {
    dplyr::case_when(
      is.na(x) ~ "",
      x < 0.001 ~ "***",
      x < 0.01 ~ "**",
      x < 0.05 ~ "*",
      TRUE ~ ""
    )
  }

  focused_safe_mean_abs <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    mean(abs(x))
  }

  focused_safe_max_abs <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    max(abs(x))
  }

  focused_theme <- function(base_size = 6) {
    ggplot2::theme_bw(
      base_family = font_family,
      base_size = base_size
    ) +
      ggplot2::theme(
        text = ggplot2::element_text(
          family = font_family,
          colour = "black"
        ),
        plot.title = ggplot2::element_text(
          family = font_family,
          face = "bold",
          size = 8,
          hjust = 0
        ),
        plot.subtitle = ggplot2::element_text(
          family = font_family,
          size = 5.5,
          hjust = 0
        ),
        plot.caption = ggplot2::element_text(
          family = font_family,
          size = 5,
          hjust = 0,
          lineheight = 0.95
        ),
        axis.title = ggplot2::element_text(
          family = font_family,
          size = 6
        ),
        axis.text = ggplot2::element_text(
          family = font_family,
          size = 5.2,
          colour = "black"
        ),
        panel.grid = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        legend.title = ggplot2::element_text(
          family = font_family,
          size = 5.5
        ),
        legend.text = ggplot2::element_text(
          family = font_family,
          size = 5
        ),
        strip.background = ggplot2::element_rect(
          fill = "white",
          colour = "black",
          linewidth = 0.25
        ),
        strip.text = ggplot2::element_text(
          family = font_family,
          face = "bold",
          size = 6
        ),
        plot.margin = ggplot2::margin(
          3, 3, 3, 3,
          unit = "mm"
        )
      )
  }

  missing_focused_comparisons <- setdiff(
    focused_comparisons,
    active_comparisons
  )

  if (length(missing_focused_comparisons) > 0L) {
    warning(
      "Focused SSL/infection/BRAF figure skipped because comparison(s) were missing: ",
      paste(missing_focused_comparisons, collapse = ", ")
    )
  } else {
    focused_gsea <- target_fgsea_results %>%
      dplyr::filter(
        Comparison %in% focused_comparisons,
        as.character(Module) %in% focused_modules
      ) %>%
      dplyr::mutate(
        Comparison = as.character(Comparison),
        Comparison_axis = factor(
          unname(focused_comparison_axis_labels[Comparison]),
          levels = unname(
            focused_comparison_axis_labels[focused_comparisons]
          )
        ),
        Module_text = as.character(Module),
        Module_axis = factor(
          unname(focused_module_labels[Module_text]),
          levels = rev(
            unname(focused_module_labels[focused_modules])
          )
        ),
        Hallmark_FDR_stars = focused_fdr_stars(padj),
        Cell_label = paste0(
          sprintf("%.2f", NES),
          Hallmark_FDR_stars
        )
      )

    if (nrow(focused_gsea) != 9L) {
      stop(
        "Expected nine focused Hallmark results, but found ",
        nrow(focused_gsea),
        "."
      )
    }

    p_focused_gsea <- ggplot2::ggplot(
      focused_gsea,
      ggplot2::aes(
        x = Comparison_axis,
        y = Module_axis,
        fill = NES
      )
    ) +
      ggplot2::geom_tile(
        colour = "white",
        linewidth = 0.35
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = Cell_label),
        family = font_family,
        size = 5.2 / 2.845276
      ) +
      ggplot2::scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        limits = c(-3, 3),
        oob = scales::squish,
        na.value = "grey92",
        name = "NES"
      ) +
      ggplot2::labs(
        title = "Hallmark enrichment direction",
        subtitle = paste0(
          "Cell values are NES; stars denote Hallmark-family FDR ",
          "(* < 0.05, ** < 0.01, *** < 0.001)"
        ),
        x = NULL,
        y = NULL
      ) +
      focused_theme() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 0,
          hjust = 0.5,
          size = 5.5
        ),
        axis.text.y = ggplot2::element_text(size = 5.5),
        legend.position = "right"
      )

    focused_gene_data <- module_gene_long %>%
      dplyr::filter(
        as.character(Comparison) %in% focused_comparisons,
        as.character(Module) %in% focused_modules
      ) %>%
      dplyr::mutate(
        Comparison = as.character(Comparison),
        Module = as.character(Module),
        Significant_directional_leading_edge =
          In_leading_edge &
          is.finite(GSEA_padj) &
          GSEA_padj < 0.05
      )

    focused_gene_selection <- focused_gene_data %>%
      dplyr::group_by(Module, Gene) %>%
      dplyr::summarise(
        N_significant_directional_leading_edges = sum(
          Significant_directional_leading_edge,
          na.rm = TRUE
        ),
        N_directional_leading_edges = sum(
          In_leading_edge,
          na.rm = TRUE
        ),
        Mean_absolute_log2FC = focused_safe_mean_abs(logFC),
        Maximum_absolute_log2FC = focused_safe_max_abs(logFC),
        .groups = "drop"
      ) %>%
      dplyr::filter(
        N_significant_directional_leading_edges >=
          focused_minimum_significant_leading_edges
      ) %>%
      dplyr::mutate(
        Module = factor(Module, levels = focused_modules)
      ) %>%
      dplyr::arrange(
        Module,
        dplyr::desc(N_significant_directional_leading_edges),
        dplyr::desc(Mean_absolute_log2FC),
        dplyr::desc(Maximum_absolute_log2FC),
        Gene
      ) %>%
      dplyr::group_by(Module) %>%
      dplyr::slice_head(n = focused_max_genes_per_pathway) %>%
      dplyr::mutate(Selection_rank = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Module = as.character(Module))

    focused_selection_counts <- focused_gene_selection %>%
      dplyr::count(Module, name = "N_selected_genes")

    missing_focused_modules <- setdiff(
      focused_modules,
      focused_selection_counts$Module[
        focused_selection_counts$N_selected_genes > 0L
      ]
    )

    if (length(missing_focused_modules) > 0L) {
      warning(
        "Focused recurrent-gene heatmap omitted module(s) with no qualifying genes: ",
        paste(missing_focused_modules, collapse = ", ")
      )
    }

    focused_gene_heatmap_long <- focused_gene_selection %>%
      dplyr::select(Module, Gene, Selection_rank) %>%
      dplyr::mutate(.Cross_join_key = 1L) %>%
      dplyr::left_join(
        tibble::tibble(
          Comparison = focused_comparisons,
          .Cross_join_key = 1L
        ),
        by = ".Cross_join_key"
      ) %>%
      dplyr::select(-.Cross_join_key) %>%
      dplyr::left_join(
        focused_gene_data %>%
          dplyr::select(
            Module,
            Gene,
            Comparison,
            logFC,
            FDR,
            In_leading_edge,
            Significant_directional_leading_edge,
            GSEA_NES,
            GSEA_padj
          ),
        by = c("Module", "Gene", "Comparison")
      ) %>%
      dplyr::mutate(
        Comparison_axis = factor(
          unname(focused_comparison_axis_labels[Comparison]),
          levels = unname(
            focused_comparison_axis_labels[focused_comparisons]
          )
        ),
        logFC_plot = pmax(
          pmin(logFC, focused_logfc_colour_limit),
          -focused_logfc_colour_limit
        )
      )

    make_focused_pathway_gene_heatmap <- function(module_name) {
      dat <- focused_gene_heatmap_long %>%
        dplyr::filter(Module == module_name)

      selected_genes <- focused_gene_selection %>%
        dplyr::filter(Module == module_name) %>%
        dplyr::arrange(Selection_rank) %>%
        dplyr::pull(Gene)

      if (length(selected_genes) == 0L) {
        return(make_empty_plot(
          unname(focused_module_labels[module_name]),
          "No recurrent FDR-significant leading-edge genes"
        ))
      }

      dat <- dat %>%
        dplyr::mutate(
          Gene = factor(Gene, levels = rev(selected_genes))
        )

      ggplot2::ggplot(
        dat,
        ggplot2::aes(
          x = Comparison_axis,
          y = Gene,
          fill = logFC_plot
        )
      ) +
        ggplot2::geom_tile(
          colour = "white",
          linewidth = 0.25,
          width = 0.96,
          height = 0.92
        ) +
        ggplot2::geom_point(
          data = dat %>%
            dplyr::filter(Significant_directional_leading_edge),
          shape = 21,
          size = 1.35,
          stroke = 0.25,
          colour = "white",
          fill = "black"
        ) +
        ggplot2::scale_fill_gradient2(
          low = "#2166AC",
          mid = "white",
          high = "#B2182B",
          midpoint = 0,
          limits = c(
            -focused_logfc_colour_limit,
            focused_logfc_colour_limit
          ),
          oob = scales::squish,
          na.value = "grey92",
          name = "log2FC"
        ) +
        ggplot2::labs(
          title = paste0(
            unname(focused_module_labels[module_name]),
            "\n",
            length(selected_genes),
            " recurrent leading-edge genes"
          ),
          x = NULL,
          y = NULL
        ) +
        focused_theme() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            face = "bold",
            size = 6.5,
            hjust = 0
          ),
          axis.text.x = ggplot2::element_text(
            angle = 55,
            hjust = 1,
            vjust = 1,
            size = 4.7
          ),
          axis.text.y = ggplot2::element_text(size = 4.8),
          legend.position = "right",
          plot.margin = ggplot2::margin(
            2, 2, 2, 2,
            unit = "mm"
          )
        )
    }

    focused_gene_plots <- lapply(
      focused_modules,
      make_focused_pathway_gene_heatmap
    )
    names(focused_gene_plots) <- focused_modules

    focused_gene_row <- patchwork::wrap_plots(
      focused_gene_plots,
      ncol = 3,
      guides = "collect"
    ) &
      ggplot2::theme(legend.position = "right")

    focused_main_figure <- p_focused_gsea /
      focused_gene_row +
      patchwork::plot_layout(heights = c(0.55, 2.25)) +
      patchwork::plot_annotation(
        tag_levels = "a",
        title = paste0(
          "Inflammatory pathway and gene-level concordance across ",
          "tissue and organoid contrasts"
        ),
        caption = paste0(
          "Black dots mark genes in the directional leading edge of an ",
          "FDR-significant Hallmark pathway. Positive log2FC denotes higher ",
          "expression in SSL_HP, infected organoids or BRAF V600E organoids."
        ),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(
            family = font_family,
            face = "bold",
            size = 9,
            hjust = 0
          ),
          plot.caption = ggplot2::element_text(
            family = font_family,
            size = 5,
            hjust = 0,
            lineheight = 0.95
          ),
          plot.tag = ggplot2::element_text(
            family = font_family,
            face = "bold",
            size = 8
          )
        )
      )

    save_plot_set(
      focused_main_figure,
      subdirs$focused,
      "FOCUSED_MAIN_GSEA_AND_RECURRENT_GENE_HEATMAPS",
      183,
      180
    )

    focused_position <- stats::setNames(
      seq_along(focused_comparisons),
      focused_comparisons
    )

    focused_pairwise <- pairwise_summary %>%
      dplyr::filter(
        Comparison_A %in% focused_comparisons,
        Comparison_B %in% focused_comparisons,
        as.character(Module) %in% focused_modules
      ) %>%
      dplyr::mutate(
        Module = as.character(Module),
        Comparison_1 = dplyr::if_else(
          unname(focused_position[Comparison_A]) <
            unname(focused_position[Comparison_B]),
          Comparison_A,
          Comparison_B
        ),
        Comparison_2 = dplyr::if_else(
          unname(focused_position[Comparison_A]) <
            unname(focused_position[Comparison_B]),
          Comparison_B,
          Comparison_A
        )
      )

    focused_pair_levels <- c(
      paste0(
        "SSL_HP vs H_COL\nvs\n",
        "Infected vs Non-Infected Organoids"
      ),
      paste0(
        "SSL_HP vs H_COL\nvs\n",
        "BRAF V600E vs WT Organoids"
      ),
      paste0(
        "Infected vs Non-Infected Organoids\nvs\n",
        "BRAF V600E vs WT Organoids"
      )
    )

    focused_pairwise <- focused_pairwise %>%
      dplyr::mutate(
        Pair_display_text = dplyr::case_when(
          Comparison_1 == "SSL_HP_v_H_COL" &
            Comparison_2 == "Alt_Infected" ~
            focused_pair_levels[1],
          Comparison_1 == "SSL_HP_v_H_COL" &
            Comparison_2 == "Initial_BRAF_SSL_HP_organoid" ~
            focused_pair_levels[2],
          Comparison_1 == "Alt_Infected" &
            Comparison_2 == "Initial_BRAF_SSL_HP_organoid" ~
            focused_pair_levels[3],
          TRUE ~ NA_character_
        ),
        Pair_display = factor(
          Pair_display_text,
          levels = focused_pair_levels
        ),
        Module_display = factor(
          unname(focused_module_labels[Module]),
          levels = rev(
            unname(focused_module_labels[focused_modules])
          )
        ),
        Spearman_FDR_stars = focused_fdr_stars(
          Spearman_FDR_across_all_pair_module_tests
        ),
        Cell_label = paste0(
          "rho = ",
          ifelse(
            is.finite(Spearman_rho_logFC),
            sprintf("%.2f", Spearman_rho_logFC),
            "NA"
          ),
          Spearman_FDR_stars,
          "\nLE J = ",
          ifelse(
            is.finite(Directional_leading_edge_Jaccard),
            sprintf("%.2f", Directional_leading_edge_Jaccard),
            "NA"
          ),
          " (n = ",
          Shared_directional_leading_edge_n,
          ")"
        )
      )

    if (nrow(focused_pairwise) != 9L) {
      stop(
        "Expected nine focused pair-by-pathway rows, but found ",
        nrow(focused_pairwise),
        "."
      )
    }

    p_focused_pairwise <- ggplot2::ggplot(
      focused_pairwise,
      ggplot2::aes(
        x = Pair_display,
        y = Module_display,
        fill = Spearman_rho_logFC
      )
    ) +
      ggplot2::geom_tile(
        colour = "white",
        linewidth = 0.35
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = Cell_label),
        family = font_family,
        size = 4.8 / 2.845276,
        lineheight = 0.92
      ) +
      ggplot2::scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        limits = c(-1, 1),
        oob = scales::squish,
        na.value = "grey92",
        name = "Spearman rho"
      ) +
      ggplot2::labs(
        title = "Focused gene-level pathway concordance",
        subtitle = paste0(
          "rho uses all eligible pathway genes and stars use global BH FDR. ",
          "LE J is directional leading-edge Jaccard; n is the number of ",
          "shared directional leading-edge genes."
        ),
        x = NULL,
        y = NULL,
        caption = "* FDR < 0.05; ** FDR < 0.01; *** FDR < 0.001."
      ) +
      focused_theme() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 0,
          hjust = 0.5,
          size = 5.1,
          lineheight = 0.9
        ),
        axis.text.y = ggplot2::element_text(size = 5.5),
        legend.position = "right"
      )

    save_plot_set(
      p_focused_pairwise,
      subdirs$focused,
      "FOCUSED_PAIRWISE_SPEARMAN_AND_DIRECTIONAL_LEADING_EDGE_JACCARD",
      190,
      100
    )

    data.table::fwrite(
      focused_gsea,
      file.path(
        subdirs$focused,
        "FOCUSED_GSEA_NES_AND_HALLMARK_FDR.csv"
      )
    )

    data.table::fwrite(
      focused_gene_selection,
      file.path(
        subdirs$focused,
        "FOCUSED_SELECTED_RECURRENT_LEADING_EDGE_GENES.csv"
      )
    )

    data.table::fwrite(
      focused_gene_heatmap_long,
      file.path(
        subdirs$focused,
        "FOCUSED_GENE_HEATMAP_LONG_TABLE.csv"
      )
    )

    data.table::fwrite(
      focused_pairwise,
      file.path(
        subdirs$focused,
        "FOCUSED_PAIRWISE_CONCORDANCE_METRICS.csv"
      )
    )

    writeLines(
      c(
        "FOCUSED FIGURE INTERPRETATION NOTES",
        "===================================",
        "",
        "Comparison labels:",
        "  SSL_HP vs H_COL",
        "  Infected vs Non-Infected Organoids",
        "  BRAF V600E vs WT Organoids",
        "",
        "Main figure:",
        "  NES and Hallmark-family FDR show pathway direction and significance.",
        "  Recurrent genes must enter the directional leading edge of at least two FDR-significant pathway results.",
        "  Black dots identify significant directional leading-edge membership.",
        "",
        "Concordance panel:",
        "  Spearman rho measures whole-pathway log2FC concordance.",
        "  Spearman stars use BH FDR across every pair-module test.",
        "  Directional leading-edge Jaccard identifies shared GSEA driver genes even when enrichment directions oppose.",
        "",
        "The full analysis of every comparison and every pair remains in folders 01-07."
      ),
      file.path(
        subdirs$focused,
        "FOCUSED_FIGURE_METHODS_AND_INTERPRETATION.txt"
      )
    )
  }
}

writeLines(
  c(
    "DE PATHWAY GENE CONCORDANCE ANALYSIS",
    "====================================",
    "",
    paste0("Run date: ", Sys.time()),
    paste0("Input directory: ", input_dir),
    paste0("Output directory: ", out_dir),
    "",
    "Primary answer:",
    "  Shared_positive_leading_edge_n and Shared_positive_leading_edge_genes",
    "  identify whether the same genes drive concordant positive enrichment.",
    "  Shared_directional_leading_edge_n and Directional_leading_edge_Jaccard",
    "  identify the same GSEA driver genes even when enrichment directions differ.",
    "  Spearman_rho_logFC measures broader effect-size concordance.",
    "  Scatterplots show shared leading-edge n/Jaccard plus Spearman rho and global BH-FDR only.",
    "  Raw Spearman p values and all additional overlap statistics remain in the CSV tables.",
    "  Significant_up_in_both_n counts genes with positive logFC and FDR <= 0.05",
    "  in both comparisons.",
    "  Pairwise overlap statistics use the genes eligible in both matched rank universes.",
    "  Hypergeometric overlap p values and BH-adjusted values are exported as support.",
    "",
    "Method:",
    paste0("  Hallmark rank universe: DE-table logCPM >= ", hallmark_logcpm_cutoff, "."),
    "  Hallmark fgsea: complete Hallmark collection per comparison.",
    "  Final 24-gene gastric-metaplasia rank universe: ALL_VALID_GENES; no logCPM filter.",
    "  Primary edgeR rank: sign(logFC) * sqrt(LR); abs(LR) is never used.",
    "",
    "A positive leading-edge driver requires positive NES, leading-edge membership,",
    "GSEA eligibility and positive gene logFC.",
    "",
    "Comparison inclusion:",
    "  Every comparison with File_found == TRUE is fully analysed.",
    "  Primary == FALSE means only that the comparison is outside the original core five;",
    "  it does not mean that Fusobacterium, Serratia or BRAF was skipped.",
    "  This single script writes every possible pair as an individual table folder and scatterplot.",
    "  Every DE input is read directly from the GM_publication_ready/Concordance analysis folder; no DE result is generated from CPM.",
    "",
    "Circularity note:",
    "  The final 24-gene metaplasia list was partly defined by positive DE in",
    "  SSL_HP vs H_COL. Positive behaviour in that defining contrast is expected",
    "  by construction. Concordance with infection and ES contrasts is the",
    "  informative independent result.",
    "",
    "Inspect first:",
    "  07_Tables/INFECTION_VS_EACH_AVAILABLE_COMPARISON_SUMMARY.csv",
    "  07_Tables/PAIRWISE_CONCORDANCE_SUMMARY_ALL_PAIRS.csv",
    "  07_Tables/CONSENSUS_DRIVER_GENES_PRIMARY_FIVE_COMPARISONS.csv",
    "  07_Tables/LEADING_EDGE_GENES_LONG.csv",
    "  03_ALL_pairwise_scatterplots/",
    "  03_ALL_pairwise_single_module_scatterplots/",
    "  04_ALL_pairwise_gene_tables/",
    "  03_Primary_pairwise_scatterplots/",
    "  05_Driver_gene_heatmaps/",
    "  06_Pairwise_summary_heatmaps/",
    "  06_Pairwise_summary_heatmaps/INFECTION_vs_each_COMBINED_CONCORDANCE_SUMMARY.*",
    "  08_Focused_SSL_Infection_BRAF/"
  ),
  file.path(out_dir, "RUN_NOTES.txt")
)

# Package and session records ----
versions <- data.frame(
  Package = c(
    "R", "data.table", "dplyr", "tidyr", "ggplot2", "fgsea",
    "msigdbr", "AnnotationDbi", "org.Hs.eg.db"
  ),
  Version = c(
    paste0(R.version$major, ".", R.version$minor),
    as.character(packageVersion("data.table")),
    as.character(packageVersion("dplyr")),
    as.character(packageVersion("tidyr")),
    as.character(packageVersion("ggplot2")),
    as.character(packageVersion("fgsea")),
    as.character(packageVersion("msigdbr")),
    as.character(packageVersion("AnnotationDbi")),
    as.character(packageVersion("org.Hs.eg.db"))
  ),
  stringsAsFactors = FALSE
)
fwrite(versions, file.path(subdirs$qc, "PACKAGE_VERSIONS.csv"))
capture.output(sessionInfo(), file = file.path(subdirs$qc, "SESSION_INFO.txt"))

message("\nDONE.")
message("Output directory: ", out_dir)
message("Inspect first:")
message("  07_Tables/INFECTION_VS_EACH_AVAILABLE_COMPARISON_SUMMARY.csv")
message("  07_Tables/CONSENSUS_DRIVER_GENES_PRIMARY_FIVE_COMPARISONS.csv")
message("All DE tables were read directly from: ", input_dir)
message("  03_ALL_pairwise_scatterplots/")
message("  03_ALL_pairwise_single_module_scatterplots/")
message("  04_ALL_pairwise_gene_tables/")
message("  03_Primary_pairwise_scatterplots/")
message("  05_Driver_gene_heatmaps/")
message("  06_Pairwise_summary_heatmaps/")
message("  08_Focused_SSL_Infection_BRAF/")
