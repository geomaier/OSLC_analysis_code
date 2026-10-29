# Tissue GSEA and figures for BRAF V600E versus WT within SSL_HP.
# Run from Tissue_Bulk_RNAseq after configuring config/config.R.
# Hallmark analyses retain their expression filter; gastric-signature ranks have no logCPM filter.

source("config/config.R")

gc()

set.seed(1)
options(stringsAsFactors = FALSE)
options(contrasts = c("contr.treatment", "contr.poly"))

input_dir <- PATH_BRAF_DE

de_file <- "CGG014_Maier_V600E_v_WT_in_SSL_HP_DE_results.txt"
de_path <- file.path(input_dir, de_file)

out_dir <- OUTPUT_BRAF_OVERALL
ensure_new_output(out_dir)

install_missing_packages <- FALSE

require_exactly_24_metaplasia_genes <- TRUE

hallmark_logcpm_cutoff <- 2
minimum_finite_lr_values <- 100L
minimum_finite_signed_stat_values <- 100L
minimum_ranked_genes_warning <- 1000L
warn_if_finite_lr_fraction_below <- 0.95

fgsea_min_size <- 15L
fgsea_max_size <- 500L
fgsea_custom_min_size <- 3L
fgsea_eps <- 0
fgsea_score_type <- "std"
fgsea_gsea_param <- 1
fgsea_nproc <- 0L
fgsea_nperm_simple <- 10000L
set_seed <- 1L

module_reference_logcpm_cutoff <- 2
module_reference_logfc_cutoff <- 1
module_max_genes <- 12L
full_hallmark_include_genes_absent_from_de_table <- FALSE

logfc_colour_limit <- 2
nes_colour_limit <- 3

font_family <- "Arial"
font_base_size <- 6
axis_title_size <- 6
axis_text_size <- 5.0
legend_title_size <- 5.0
legend_text_size <- 5.0
plot_title_size <- 7
panel_tag_size <- 8
thin_line <- 0.25
export_dpi <- 600

positive_colour <- "#B2182B"
negative_colour <- "#2166AC"
mid_colour <- "white"
tile_border_colour <- "white"
tile_border_linewidth <- thin_line
tile_fraction <- 0.96
shared_label_pt <- 5.0
x_label_angle <- 55
curve_colour <- positive_colour

metaplasia_display_title <-
  "Serrated Lesion Gastric Metaplasia Gene Signature"
metaplasia_display_title_wrapped <- paste(
  "Serrated Lesion Gastric",
  "Metaplasia Gene Signature",
  sep = "\n"
)

contrast_key <- "BRAF_V600E_vs_WT_SSL_HP"
contrast_group <- "SSL_HP"
contrast_display_label <- "BRAF V600E vs WT\nSSL_HP"
contrast_long_label <- "BRAF V600E versus WT in SSL_HP"

metaplasia_signature_candidates <- PATH_GASTRIC_SIGNATURE

de_file_table <- data.frame(
  comparison = contrast_key,
  group = contrast_group,
  group_order = 1L,
  organism = "BRAF_V600E_vs_WT",
  file = de_file,
  display_order = 1L,
  stringsAsFactors = FALSE
)

de_file_table$path <- de_path
comparison_keys <- contrast_key
comparison_display_labels <- stats::setNames(
  contrast_display_label,
  contrast_key
)
comparison_long_labels <- stats::setNames(
  contrast_long_label,
  contrast_key
)
curve_colours <- stats::setNames(curve_colour, contrast_key)

cran_packages <- c(
  "data.table", "dplyr", "tidyr", "stringr", "purrr", "ggplot2",
  "patchwork", "scales", "msigdbr", "tibble", "readr"
)

bioc_packages <- c("fgsea")

install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!install_missing_packages) {
      stop(
        "Package '", pkg, "' is not installed. Install it or set ",
        "install_missing_packages <- TRUE (explicit opt-in)."
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

for (pkg in cran_packages) install_if_missing(pkg, bioc = FALSE)
for (pkg in bioc_packages) install_if_missing(pkg, bioc = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(msigdbr)
  library(tibble)
  library(readr)
  library(fgsea)
  library(grid)
})

if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir)
}

if (!file.exists(de_path)) {
  stop(
    "The required BRAF DE file is missing:\n",
    de_path
  )
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

subdirs <- list(
  qc = file.path(out_dir, "01_QC"),
  tables = file.path(out_dir, "02_Tables"),
  hallmark_curves = file.path(out_dir, "03_Hallmark_GSEA_curves"),
  metaplasia_curves = file.path(out_dir, "04_Gastric_metaplasia_GSEA_curve"),
  hallmark_heatmaps = file.path(out_dir, "05_Pathway_NES_heatmaps"),
  full_modules = file.path(out_dir, "06_Full_signature_gene_log2FC_heatmaps"),
  selected_modules = file.path(out_dir, "07_Selected_gene_log2FC_heatmaps"),
  combined = file.path(out_dir, "08_Combined_figures"),
  ranks = file.path(out_dir, "09_Ranked_gene_tables"),
  gene_summary = file.path(out_dir, "10_Gene_level_summary"),
  four_module = file.path(out_dir, "11_Four_module_heatmap_panel")
)

for (d in subdirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

write_test_file <- file.path(out_dir, "_write_test.tmp")
tryCatch(
  {
    writeLines("write test", write_test_file)
    unlink(write_test_file, force = TRUE)
  },
  error = function(e) {
    stop(
      "Cannot write to output directory:\n", out_dir,
      "\n\nR error:\n", conditionMessage(e)
    )
  }
)

data.table::fwrite(
  de_file_table,
  file.path(subdirs$qc, "BRAF_DE_FILE_USED.csv")
)

safe_filename <- function(x) {
  x |>
    as.character() |>
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_") |>
    stringr::str_replace_all("_+", "_") |>
    stringr::str_replace_all("^_|_$", "")
}

column_key <- function(x) {
  toupper(gsub("[^A-Za-z0-9]", "", trimws(as.character(x))))
}

detect_col <- function(df, candidates, required = FALSE, label = "column") {
  nms <- names(df)

  if (length(nms) == 0L) {
    if (required) stop("No columns are available while detecting ", label, ".")
    return(NA_character_)
  }

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

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

clean_gene_symbols <- function(x) {
  x <- clean_character(x)
  x <- toupper(x)
  x[x %in% c("NA", "N/A", ".", "NULL")] <- NA_character_
  x
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  min(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  stats::median(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  max(x)
}

first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0L) return(NA_character_)
  hit[1]
}

format_numeric_p <- function(p) {
  if (length(p) == 0L || is.na(p)) return("NA")
  if (p < 1e-3) return(formatC(p, format = "e", digits = 1))
  sprintf("%.3f", p)
}

significance_stars <- function(padj) {
  dplyr::case_when(
    is.na(padj) ~ "",
    padj < 0.001 ~ "***",
    padj < 0.01 ~ "**",
    padj < 0.05 ~ "*",
    TRUE ~ ""
  )
}

geom_text_size_from_pt <- function(pt) {
  pt / 2.845276
}

nature_theme <- function(base_size = font_base_size) {
  ggplot2::theme_classic(
    base_size = base_size,
    base_family = font_family
  ) +
    ggplot2::theme(
      text = ggplot2::element_text(family = font_family, colour = "black"),
      plot.title = ggplot2::element_text(
        family = font_family,
        face = "bold",
        hjust = 0,
        size = plot_title_size
      ),
      axis.title = ggplot2::element_text(
        family = font_family,
        size = axis_title_size,
        colour = "black"
      ),
      axis.text = ggplot2::element_text(
        family = font_family,
        size = axis_text_size,
        colour = "black"
      ),
      legend.title = ggplot2::element_text(
        family = font_family,
        size = legend_title_size,
        colour = "black"
      ),
      legend.text = ggplot2::element_text(
        family = font_family,
        size = legend_text_size,
        colour = "black"
      ),
      axis.line = ggplot2::element_line(colour = "black", linewidth = thin_line),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = thin_line),
      panel.grid = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(4, 4, 4, 4, "mm")
    )
}

ggplot2::theme_set(nature_theme())

vertical_heatmap_guide <- function(barheight_mm = 18, barwidth_mm = 2.0) {
  ggplot2::guide_colourbar(
    title.position = "right",
    title.hjust = 0.5,
    title.vjust = 0.5,
    title.theme = ggplot2::element_text(
      family = font_family,
      size = legend_title_size,
      angle = 90,
      hjust = 0.5,
      vjust = 0.5,
      colour = "black",
      margin = ggplot2::margin(l = 2, r = 0, unit = "mm")
    ),
    label.position = "right",
    label.hjust = 0.5,
    label.vjust = 0.5,
    label.theme = ggplot2::element_text(
      family = font_family,
      size = legend_text_size,
      colour = "black"
    ),
    barwidth = grid::unit(barwidth_mm, "mm"),
    barheight = grid::unit(barheight_mm, "mm"),
    ticks = TRUE
  )
}

save_pdf_safely <- function(plot, path, width_mm, height_mm) {
  cairo_success <- FALSE

  if (capabilities("cairo")) {
    cairo_success <- tryCatch(
      {
        ggplot2::ggsave(
          path,
          plot,
          width = width_mm,
          height = height_mm,
          units = "mm",
          device = grDevices::cairo_pdf,
          bg = "white"
        )
        TRUE
      },
      error = function(e) {
        warning(
          "Cairo PDF export failed for:\n", path,
          "\nTrying the standard PDF device.\n",
          conditionMessage(e),
          call. = FALSE
        )
        FALSE
      }
    )
  }

  if (!cairo_success) {
    base_pdf_success <- tryCatch(
      {
        ggplot2::ggsave(
          path,
          plot,
          width = width_mm,
          height = height_mm,
          units = "mm",
          device = grDevices::pdf,
          bg = "white"
        )
        TRUE
      },
      error = function(e) {
        warning(
          "Standard PDF export also failed for:\n", path,
          "\nThe run will continue and retain PNG/TIFF/JPG outputs.\n",
          conditionMessage(e),
          call. = FALSE
        )
        FALSE
      }
    )

    invisible(base_pdf_success)
  }
}

save_plot_set <- function(
    plot,
    directory,
    stem,
    width_mm,
    height_mm,
    dpi = export_dpi
) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  stem <- safe_filename(stem)
  prefix <- file.path(directory, stem)

  save_pdf_safely(
    plot,
    paste0(prefix, ".pdf"),
    width_mm,
    height_mm
  )

  ggplot2::ggsave(
    paste0(prefix, ".png"),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = dpi,
    bg = "white"
  )

  ggplot2::ggsave(
    paste0(prefix, ".tif"),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = dpi,
    device = grDevices::tiff,
    compression = "lzw",
    bg = "white"
  )

  ggplot2::ggsave(
    paste0(prefix, ".jpg"),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = dpi,
    device = grDevices::jpeg,
    quality = 100,
    bg = "white"
  )

  gc()
}

make_empty_plot <- function(title_text, message_text) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = message_text,
      family = font_family,
      size = 2.2
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(title = title_text, x = NULL, y = NULL) +
    ggplot2::theme_void(base_family = font_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = font_family,
        face = "bold",
        size = plot_title_size,
        hjust = 0
      )
    )
}

read_gene_list <- function(file) {
  if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    dat <- data.table::fread(
      file,
      data.table = FALSE,
      check.names = FALSE,
      na.strings = c("", "NA", "N/A", "na")
    )

    gene_col <- detect_col(
      dat,
      c("gene_symbol", "Gene", "gene", "SYMBOL", "Symbol", "symbol"),
      required = FALSE
    )

    if (is.na(gene_col)) gene_col <- names(dat)[1]
    genes <- dat[[gene_col]]
  } else {
    lines <- readLines(file, warn = FALSE)

    if (length(lines) == 1L && grepl("[,;\\t]", lines[1])) {
      genes <- unlist(strsplit(lines, "[,;\\t]"))
    } else {
      genes <- lines
    }
  }

  genes <- clean_gene_symbols(genes)
  genes <- genes[!is.na(genes) & genes != ""]
  genes <- genes[!stringr::str_detect(genes, "^#")]
  genes <- genes[!genes %in% c("GENE", "GENES", "SYMBOL", "GENE_SYMBOL")]
  unique(genes)
}

metaplasia_signature_file <- first_existing_file(
  metaplasia_signature_candidates
)

if (is.na(metaplasia_signature_file)) {
  stop(
    "The final gastric-metaplasia signature file was not found. Checked:\n",
    paste(metaplasia_signature_candidates, collapse = "\n")
  )
}

metaplasia_genes <- read_gene_list(metaplasia_signature_file)

if (length(metaplasia_genes) < 3L) {
  stop(
    "The final gastric-metaplasia signature contains fewer than three genes: ",
    metaplasia_signature_file
  )
}

if (length(metaplasia_genes) != 24L) {
  message_text <- paste0(
    "The final gastric-metaplasia file contains ",
    length(metaplasia_genes),
    " unique genes rather than the required 24: ",
    metaplasia_signature_file
  )

  if (require_exactly_24_metaplasia_genes) {
    stop(message_text)
  } else {
    warning(message_text, call. = FALSE)
  }
}

writeLines(
  c(
    paste0("Signature file: ", metaplasia_signature_file),
    paste0("Unique genes: ", length(metaplasia_genes)),
    "",
    metaplasia_genes
  ),
  file.path(subdirs$qc, "FINAL_GASTRIC_METAPLASIA_SIGNATURE_USED.txt")
)

read_de_table_auto <- function(file) {
  dat <- tryCatch(
    data.table::fread(
      file,
      data.table = FALSE,
      check.names = FALSE,
      showProgress = FALSE,
      fill = TRUE,
      na.strings = c("", "NA", "N/A", "na")
    ),
    error = function(e) NULL
  )

  if (!is.null(dat) && ncol(dat) > 1L) return(dat)

  dat <- tryCatch(
    readr::read_delim(
      file,
      delim = "\t",
      show_col_types = FALSE,
      progress = FALSE,
      name_repair = "unique"
    ),
    error = function(e) NULL
  )

  if (!is.null(dat) && ncol(dat) > 1L) return(as.data.frame(dat))

  stop("Could not read DE table: ", file)
}

convert_numeric_strict <- function(x, column_name, comparison) {
  original <- clean_character(x)
  converted <- suppressWarnings(as.numeric(original))
  invalid <- !is.na(original) & is.na(converted)

  if (any(invalid)) {
    examples <- utils::head(unique(original[invalid]), 8)
    stop(
      "Non-numeric values were found in '", column_name,
      "' for ", comparison, ": ",
      paste(examples, collapse = "; ")
    )
  }

  converted
}

standardise_de_table <- function(file, comparison) {
  raw <- read_de_table_auto(file)
  names(raw) <- trimws(names(raw))

  nonempty_columns <- vapply(
    raw,
    function(x) any(!is.na(x) & trimws(as.character(x)) != ""),
    logical(1)
  )
  raw <- raw[, nonempty_columns, drop = FALSE]

  gene_col <- detect_col(
    raw,
    c(
      "SYMBOL", "Symbol", "symbol", "gene_symbol", "GeneSymbol",
      "geneSymbol", "HGNC.symbol", "Gene", "gene"
    ),
    required = TRUE,
    label = paste0(comparison, " gene-symbol column")
  )

  logfc_col <- detect_col(
    raw,
    c("logFC", "LogFC", "logfc", "log2FoldChange", "log2FC", "LFC"),
    required = TRUE,
    label = paste0(comparison, " logFC column")
  )

  logcpm_col <- detect_col(
    raw,
    c(
      "logCPM", "LogCPM", "logcpm", "log2CPM", "Log2CPM",
      "log2cpm", "AveExpr", "AveExpression", "mean_logCPM",
      "meanLogCPM", "average_logCPM", "averageLogCPM"
    ),
    required = TRUE,
    label = paste0(comparison, " logCPM column")
  )

  p_col <- detect_col(
    raw,
    c(
      "PValue", "P.Value", "Pvalue", "pvalue", "p.value",
      "P", "p", "p_val", "pval"
    ),
    required = TRUE,
    label = paste0(comparison, " nominal P-value column")
  )

  fdr_col <- detect_col(
    raw,
    c(
      "FDR", "fdr", "padj", "adj.P.Val", "adj.PValue",
      "adj_p_value", "qvalue", "q_value", "p_val_adj"
    ),
    required = FALSE,
    label = paste0(comparison, " adjusted P-value column")
  )

  lr_col <- detect_col(
    raw,
    c("LR", "LRT", "likelihood_ratio", "LikelihoodRatio"),
    required = FALSE,
    label = paste0(comparison, " LR column")
  )

  signed_stat_col <- detect_col(
    raw,
    c(
      "t", "T", "stat", "STAT", "Statistic", "moderated_t",
      "tstat", "z", "Z", "score"
    ),
    required = FALSE,
    label = paste0(comparison, " signed-statistic column")
  )

  dat <- data.frame(
    original_row = seq_len(nrow(raw)),
    Gene = clean_gene_symbols(raw[[gene_col]]),
    logFC = convert_numeric_strict(raw[[logfc_col]], logfc_col, comparison),
    logCPM = convert_numeric_strict(raw[[logcpm_col]], logcpm_col, comparison),
    PValue = convert_numeric_strict(raw[[p_col]], p_col, comparison),
    FDR = if (!is.na(fdr_col)) {
      convert_numeric_strict(raw[[fdr_col]], fdr_col, comparison)
    } else {
      rep(NA_real_, nrow(raw))
    },
    LR = if (!is.na(lr_col)) {
      convert_numeric_strict(raw[[lr_col]], lr_col, comparison)
    } else {
      rep(NA_real_, nrow(raw))
    },
    signed_stat = if (!is.na(signed_stat_col)) {
      convert_numeric_strict(raw[[signed_stat_col]], signed_stat_col, comparison)
    } else {
      rep(NA_real_, nrow(raw))
    },
    comparison = comparison,
    source_file = basename(file),
    stringsAsFactors = FALSE
  )

  input_rows <- nrow(dat)

  invalid_p <- is.finite(dat$PValue) & (dat$PValue < 0 | dat$PValue > 1)
  if (any(invalid_p)) {
    stop(
      comparison, " contains nominal P values outside [0,1]."
    )
  }

  invalid_fdr <- is.finite(dat$FDR) & (dat$FDR < 0 | dat$FDR > 1)
  if (any(invalid_fdr)) {
    stop(
      comparison, " contains adjusted P values outside [0,1]."
    )
  }

  materially_negative_lr <- is.finite(dat$LR) & dat$LR < -1e-8
  if (any(materially_negative_lr)) {
    data.table::fwrite(
      dat[materially_negative_lr, , drop = FALSE],
      file.path(
        subdirs$qc,
        paste0(safe_filename(comparison), "_MATERIALLY_NEGATIVE_LR_ROWS.csv")
      )
    )

    stop(
      comparison, " contains ", sum(materially_negative_lr),
      " materially negative LR value(s). Do not use abs(LR)."
    )
  }

  dat$LR[is.finite(dat$LR) & dat$LR < 0] <- 0

  dat <- dat |>
    dplyr::filter(
      !is.na(Gene),
      Gene != "",
      is.finite(logFC),
      is.finite(logCPM),
      is.finite(PValue)
    )

  if (nrow(dat) == 0L) {
    stop("No valid DE rows remained for ", comparison, ".")
  }

  bh_fallback <- stats::p.adjust(dat$PValue, method = "BH")
  dat$FDR <- ifelse(is.finite(dat$FDR), dat$FDR, bh_fallback)

  column_mapping <- data.frame(
    comparison = comparison,
    role = c(
      "Gene", "logFC", "logCPM", "PValue", "FDR", "LR", "signed_stat"
    ),
    input_column = c(
      gene_col, logfc_col, logcpm_col, p_col, fdr_col, lr_col,
      signed_stat_col
    ),
    stringsAsFactors = FALSE
  )

  data.table::fwrite(
    column_mapping,
    file.path(
      subdirs$qc,
      paste0(safe_filename(comparison), "_DE_COLUMN_MAPPING.csv")
    )
  )

  structure_qc <- data.frame(
    comparison = comparison,
    source_file = basename(file),
    input_rows = input_rows,
    valid_rows = nrow(dat),
    unique_symbols_before_collapse = dplyr::n_distinct(dat$Gene),
    finite_logFC = sum(is.finite(dat$logFC)),
    finite_logCPM = sum(is.finite(dat$logCPM)),
    finite_LR = sum(is.finite(dat$LR)),
    finite_signed_stat = sum(is.finite(dat$signed_stat)),
    finite_PValue = sum(is.finite(dat$PValue)),
    finite_FDR = sum(is.finite(dat$FDR)),
    stringsAsFactors = FALSE
  )

  data.table::fwrite(
    structure_qc,
    file.path(
      subdirs$qc,
      paste0(safe_filename(comparison), "_DE_STRUCTURE_QC.csv")
    )
  )

  dat
}

standardised_de_list <- list()

for (i in seq_len(nrow(de_file_table))) {
  comp <- de_file_table$comparison[i]
  path <- de_file_table$path[i]

  message("Reading DE table: ", comp)
  message("  ", path)

  de_one <- standardise_de_table(path, comp)
  standardised_de_list[[comp]] <- de_one

  data.table::fwrite(
    de_one,
    file.path(
      subdirs$tables,
      paste0(safe_filename(comp), "_DE_STANDARDISED_ALL_VALID_ROWS.csv")
    )
  )
}

make_ranked_object <- function(
    de,
    comparison,
    apply_logcpm_filter
) {
  dat <- de

  filter_mode <- if (apply_logcpm_filter) {
    paste0("LOGCPM_GE_", hallmark_logcpm_cutoff)
  } else {
    "ALL_VALID_GENES"
  }

  if (apply_logcpm_filter) {
    dat <- dat |>
      dplyr::filter(logCPM >= hallmark_logcpm_cutoff)
  }

  n_before_rank <- nrow(dat)

  if (n_before_rank == 0L) {
    stop("No genes remained for ", comparison, " / ", filter_mode, ".")
  }

  positive_p <- dat$PValue[
    is.finite(dat$PValue) & dat$PValue > 0 & dat$PValue <= 1
  ]

  min_positive_p <- if (length(positive_p) > 0L) {
    min(positive_p)
  } else {
    1e-300
  }

  zero_p_replacement <- max(
    min_positive_p * 0.1,
    .Machine$double.xmin
  )

  dat$PValue_safe <- ifelse(
    dat$PValue == 0,
    zero_p_replacement,
    dat$PValue
  )

  n_finite_lr <- sum(is.finite(dat$LR))
  n_finite_signed_stat <- sum(is.finite(dat$signed_stat))
  finite_lr_fraction <- n_finite_lr / nrow(dat)

  if (n_finite_lr >= minimum_finite_lr_values) {
    if (finite_lr_fraction < warn_if_finite_lr_fraction_below) {
      warning(
        comparison, " / ", filter_mode, ": only ",
        round(100 * finite_lr_fraction, 1),
        "% of retained rows have finite LR values. The established workflow ",
        "still uses LR because >=", minimum_finite_lr_values,
        " values are available, but inspect the DE table before interpretation.",
        call. = FALSE
      )
    }

    dat$rank_metric <- sign(dat$logFC) * sqrt(dat$LR)
    rank_method <- "sign(logFC) * sqrt(LR)"
  } else if (n_finite_signed_stat >= minimum_finite_signed_stat_values) {
    dat$rank_metric <- dat$signed_stat
    rank_method <- "signed test statistic"
  } else {
    dat$rank_metric <- sign(dat$logFC) * -log10(dat$PValue_safe)
    rank_method <- "sign(logFC) * -log10(PValue)"
  }

  dat <- dat |>
    dplyr::filter(is.finite(rank_metric))

  if (nrow(dat) == 0L) {
    stop("No finite ranking values remained for ", comparison, ".")
  }

  duplicate_rows <- dat |>
    dplyr::group_by(Gene) |>
    dplyr::filter(dplyr::n() > 1L) |>
    dplyr::ungroup()

  if (nrow(duplicate_rows) > 0L) {
    data.table::fwrite(
      duplicate_rows,
      file.path(
        subdirs$qc,
        paste0(
          safe_filename(comparison),
          "_",
          safe_filename(filter_mode),
          "_DUPLICATE_SYMBOL_ROWS_BEFORE_COLLAPSE.csv"
        )
      )
    )
  }

  dat_unique <- dat |>
    dplyr::group_by(Gene) |>
    dplyr::arrange(
      dplyr::desc(abs(rank_metric)),
      PValue_safe,
      dplyr::desc(abs(logFC)),
      original_row,
      .by_group = TRUE
    ) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(rank_metric), Gene)

  ranks <- dat_unique$rank_metric
  names(ranks) <- dat_unique$Gene

  if (anyDuplicated(names(ranks)) > 0L) {
    stop("Duplicated gene symbols remain after collapse for ", comparison, ".")
  }

  dat_unique$rank_metric_before_tie_break <- dat_unique$rank_metric

  if (any(duplicated(ranks))) {
    ranks <- ranks - seq_along(ranks) * 1e-12
  }

  dat_unique$rank_metric <- ranks
  dat_unique$rank_metric_fgsea <- ranks
  dat_unique <- dat_unique |>
    dplyr::arrange(dplyr::desc(rank_metric), Gene)

  ranks <- sort(ranks, decreasing = TRUE)

  if (length(ranks) < minimum_ranked_genes_warning) {
    warning(
      comparison, " / ", filter_mode, " contains only ", length(ranks),
      " ranked genes after symbol collapse.",
      call. = FALSE
    )
  }

  metadata <- data.frame(
    comparison = comparison,
    filter_mode = filter_mode,
    rank_method = rank_method,
    n_rows_before_rank = n_before_rank,
    n_ranked_genes = length(ranks),
    n_finite_lr = n_finite_lr,
    finite_lr_fraction = finite_lr_fraction,
    n_finite_signed_stat = n_finite_signed_stat,
    n_duplicate_rows_before_collapse = nrow(duplicate_rows),
    min_logCPM = safe_min(dat_unique$logCPM),
    median_logCPM = safe_median(dat_unique$logCPM),
    min_logFC = safe_min(dat_unique$logFC),
    max_logFC = safe_max(dat_unique$logFC),
    stringsAsFactors = FALSE
  )

  list(
    ranks = ranks,
    table = dat_unique,
    metadata = metadata,
    method = rank_method,
    filter_mode = filter_mode
  )
}

hallmark_rank_objects <- list()
metaplasia_rank_objects <- list()

for (comp in comparison_keys) {
  hallmark_rank_objects[[comp]] <- make_ranked_object(
    standardised_de_list[[comp]],
    comp,
    apply_logcpm_filter = TRUE
  )

  metaplasia_rank_objects[[comp]] <- make_ranked_object(
    standardised_de_list[[comp]],
    comp,
    apply_logcpm_filter = FALSE
  )

  data.table::fwrite(
    hallmark_rank_objects[[comp]]$table,
    file.path(
      subdirs$ranks,
      paste0(safe_filename(comp), "_LOGCPM_GE_2_RANKED_GENES.csv")
    )
  )

  data.table::fwrite(
    metaplasia_rank_objects[[comp]]$table,
    file.path(
      subdirs$ranks,
      paste0(safe_filename(comp), "_ALL_VALID_GENES_FOR_GM_GSEA.csv")
    )
  )
}

rank_metadata <- dplyr::bind_rows(
  c(
    lapply(hallmark_rank_objects, `[[`, "metadata"),
    lapply(metaplasia_rank_objects, `[[`, "metadata")
  )
) |>
  dplyr::left_join(
    de_file_table |>
      dplyr::select(
        comparison,
        group,
        group_order,
        organism,
        display_order,
        file
      ),
    by = "comparison"
  ) |>
  dplyr::arrange(group_order, display_order, filter_mode)

data.table::fwrite(
  rank_metadata,
  file.path(subdirs$qc, "ALL_RANKED_GENE_METADATA.csv")
)

get_msig_hallmark <- function() {
  out <- tryCatch(
    suppressWarnings(
      msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
    ),
    error = function(e) NULL
  )

  if (is.null(out) || nrow(out) == 0L) {
    out <- tryCatch(
      suppressWarnings(
        msigdbr::msigdbr(species = "Homo sapiens", category = "H")
      ),
      error = function(e) NULL
    )
  }

  if (is.null(out) || nrow(out) == 0L) {
    stop("Could not load the MSigDB Hallmark collection through msigdbr.")
  }

  as.data.frame(out)
}

message("Loading MSigDB Hallmark collection...")
msig_h_raw <- get_msig_hallmark()

msig_gene_col <- detect_col(
  msig_h_raw,
  c("gene_symbol", "human_gene_symbol", "db_gene_symbol", "symbol"),
  required = TRUE,
  label = "MSigDB gene-symbol column"
)

msig_pathway_col <- detect_col(
  msig_h_raw,
  c("gs_name", "geneset", "gene_set_name", "gs_exact_source"),
  required = TRUE,
  label = "MSigDB pathway-name column"
)

msig_h <- data.frame(
  pathway = as.character(msig_h_raw[[msig_pathway_col]]),
  gene = clean_gene_symbols(msig_h_raw[[msig_gene_col]]),
  stringsAsFactors = FALSE
) |>
  dplyr::filter(
    !is.na(pathway), pathway != "",
    !is.na(gene), gene != ""
  ) |>
  dplyr::distinct(pathway, gene)

hallmark_pathways <- split(msig_h$gene, msig_h$pathway)
hallmark_pathways <- lapply(
  hallmark_pathways,
  function(x) sort(unique(clean_gene_symbols(x)))
)

get_required_pathway <- function(pathways, requested) {
  available <- names(pathways)
  hit <- available[toupper(available) == toupper(requested)]

  if (length(hit) == 0L) {
    stop("Required Hallmark pathway was not found: ", requested)
  }

  list(
    name = hit[1],
    genes = sort(unique(clean_gene_symbols(pathways[[hit[1]]])))
  )
}

ifna_info <- get_required_pathway(
  hallmark_pathways,
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"
)

inflammatory_info <- get_required_pathway(
  hallmark_pathways,
  "HALLMARK_INFLAMMATORY_RESPONSE"
)

tnfa_info <- get_required_pathway(
  hallmark_pathways,
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

msig_version_col <- detect_col(
  msig_h_raw,
  c("db_version", "msigdb_version", "version"),
  required = FALSE
)

msig_version_text <- if (!is.na(msig_version_col)) {
  paste(unique(na.omit(msig_h_raw[[msig_version_col]])), collapse = "; ")
} else {
  NA_character_
}

pathway_qc <- data.frame(
  pathway = c(ifna_info$name, inflammatory_info$name, tnfa_info$name),
  number_of_genes = c(
    length(ifna_info$genes),
    length(inflammatory_info$genes),
    length(tnfa_info$genes)
  ),
  stringsAsFactors = FALSE
)

data.table::fwrite(
  pathway_qc,
  file.path(subdirs$qc, "SELECTED_HALLMARK_PATHWAYS_USED.csv")
)

prepare_pathways_for_ranks <- function(pathways, ranks, min_size, max_size) {
  pathways <- lapply(pathways, function(x) {
    x <- clean_gene_symbols(x)
    sort(unique(x[!is.na(x) & x != ""]))
  })

  overlap <- vapply(
    pathways,
    function(g) sum(g %in% names(ranks)),
    integer(1)
  )

  pathways[overlap >= min_size & overlap <= max_size]
}

run_fgsea_one <- function(
    rank_object,
    pathways,
    comparison,
    family,
    min_size,
    max_size = fgsea_max_size
) {
  prepared <- prepare_pathways_for_ranks(
    pathways,
    rank_object$ranks,
    min_size,
    max_size
  )

  if (length(prepared) == 0L) {
    stop("No pathways passed overlap-size filtering for ", comparison, ".")
  }

  set.seed(set_seed)

  fg <- fgsea::fgseaMultilevel(
    pathways = prepared,
    stats = rank_object$ranks,
    minSize = min_size,
    maxSize = max_size,
    eps = fgsea_eps,
    scoreType = fgsea_score_type,
    gseaParam = fgsea_gsea_param,
    nproc = fgsea_nproc,
    nPermSimple = fgsea_nperm_simple
  )

  fg <- as.data.frame(fg)

  if (!"padj" %in% names(fg)) {
    fg$padj <- stats::p.adjust(fg$pval, method = "BH")
  }
  if (!"log2err" %in% names(fg)) fg$log2err <- NA_real_
  if (!"leadingEdge" %in% names(fg)) {
    fg$leadingEdge <- vector("list", nrow(fg))
  }

  fg |>
    dplyr::mutate(
      comparison = comparison,
      gene_set_family = family,
      filter_mode = rank_object$filter_mode,
      ranking_method = rank_object$method,
      n_ranked_genes = length(rank_object$ranks),
      leadingEdge = vapply(
        leadingEdge,
        paste,
        collapse = ";",
        FUN.VALUE = character(1)
      )
    ) |>
    dplyr::arrange(padj, dplyr::desc(abs(NES))) |>
    dplyr::select(
      comparison,
      filter_mode,
      ranking_method,
      n_ranked_genes,
      gene_set_family,
      pathway,
      pval,
      padj,
      log2err,
      ES,
      NES,
      size,
      leadingEdge
    )
}

custom_metaplasia_pathway <- list(
  FINAL_SSL_METAPLASIA_LIST = metaplasia_genes
)

hallmark_fgsea_list <- list()
metaplasia_fgsea_list <- list()

for (comp in comparison_keys) {
  message("Running Hallmark FGSEA: ", comp)

  hallmark_fgsea_list[[comp]] <- run_fgsea_one(
    hallmark_rank_objects[[comp]],
    hallmark_pathways,
    comp,
    "HALLMARK",
    min_size = fgsea_min_size
  )

  message("Running final gastric-metaplasia FGSEA: ", comp)

  metaplasia_fgsea_list[[comp]] <- run_fgsea_one(
    metaplasia_rank_objects[[comp]],
    custom_metaplasia_pathway,
    comp,
    "CUSTOM_GASTRIC_METAPLASIA_NO_LOGCPM_FILTER",
    min_size = fgsea_custom_min_size
  ) |>
    dplyr::mutate(
      p_value_interpretation = paste0(
        "Nominal custom-set P value. Only one custom pathway was tested; ",
        "padj equals pval and is not a Hallmark-family FDR."
      )
    )
}

hallmark_fgsea_all <- dplyr::bind_rows(hallmark_fgsea_list) |>
  dplyr::left_join(
    de_file_table |>
      dplyr::select(comparison, group, group_order, organism, display_order),
    by = "comparison"
  ) |>
  dplyr::arrange(group_order, display_order, padj, dplyr::desc(abs(NES)))

metaplasia_fgsea_all <- dplyr::bind_rows(metaplasia_fgsea_list) |>
  dplyr::left_join(
    de_file_table |>
      dplyr::select(comparison, group, group_order, organism, display_order),
    by = "comparison"
  ) |>
  dplyr::arrange(group_order, display_order)

combined_fgsea <- dplyr::bind_rows(
  hallmark_fgsea_all,
  metaplasia_fgsea_all
)

data.table::fwrite(
  hallmark_fgsea_all,
  file.path(subdirs$tables, "BRAF_V600E_VS_WT_SSLHP_HALLMARK_FGSEA_ALL_RESULTS.csv")
)

data.table::fwrite(
  metaplasia_fgsea_all,
  file.path(
    subdirs$tables,
    "BRAF_V600E_VS_WT_SSLHP_FINAL_GASTRIC_METAPLASIA_FGSEA_ALL_VALID_GENES.csv"
  )
)

data.table::fwrite(
  combined_fgsea,
  file.path(subdirs$tables, "BRAF_V600E_VS_WT_SSLHP_COMBINED_FGSEA_RESULTS.csv")
)

selected_hallmark_results <- hallmark_fgsea_all |>
  dplyr::filter(
    pathway %in% c(ifna_info$name, inflammatory_info$name, tnfa_info$name)
  ) |>
  dplyr::mutate(
    pathway_display = dplyr::case_when(
      pathway == ifna_info$name ~ "Interferon alpha response",
      pathway == inflammatory_info$name ~ "Inflammatory response",
      pathway == tnfa_info$name ~ "TNFA/NFKB signaling",
      TRUE ~ pathway
    )
  ) |>
  dplyr::arrange(group_order, display_order, pathway_display)

data.table::fwrite(
  selected_hallmark_results,
  file.path(
    subdirs$tables,
    "BRAF_V600E_VS_WT_SSLHP_SELECTED_IFNA_INFLAMMATORY_TNFA_FGSEA.csv"
  )
)

metaplasia_gene_audit <- purrr::map_dfr(
  comparison_keys,
  function(comp) {
    result <- metaplasia_fgsea_all |>
      dplyr::filter(comparison == comp) |>
      dplyr::slice(1)

    leading_edge <- if (nrow(result) == 1L && !is.na(result$leadingEdge[1])) {
      strsplit(result$leadingEdge[1], ";", fixed = TRUE)[[1]]
    } else {
      character(0)
    }

    data.frame(
      Gene = metaplasia_genes,
      signature_order = seq_along(metaplasia_genes),
      stringsAsFactors = FALSE
    ) |>
      dplyr::left_join(
        metaplasia_rank_objects[[comp]]$table |>
          dplyr::select(Gene, logFC, logCPM, LR, PValue, FDR),
        by = "Gene"
      ) |>
      dplyr::mutate(
        comparison = comp,
        present_in_DE_table = is.finite(logFC),
        present_in_ALL_VALID_GENES_rank =
          Gene %in% names(metaplasia_rank_objects[[comp]]$ranks),
        rank_metric = unname(metaplasia_rank_objects[[comp]]$ranks[Gene]),
        leading_edge = Gene %in% leading_edge,
        logCPM_filter_applied = FALSE
      )
  }
)

data.table::fwrite(
  metaplasia_gene_audit,
  file.path(subdirs$qc, "FINAL_GASTRIC_METAPLASIA_GSEA_GENE_AUDIT.csv")
)

collapse_de_for_heatmap <- function(de, comparison) {
  dat <- de

  positive_p <- dat$PValue[
    is.finite(dat$PValue) & dat$PValue > 0
  ]
  min_positive_p <- if (length(positive_p) > 0L) {
    min(positive_p)
  } else {
    1e-300
  }

  safe_p <- ifelse(
    dat$PValue == 0,
    max(min_positive_p * 0.1, .Machine$double.xmin),
    dat$PValue
  )

  n_finite_lr <- sum(is.finite(dat$LR))
  n_finite_signed_stat <- sum(is.finite(dat$signed_stat))

  dat$collapse_rank <- if (
    n_finite_lr >= minimum_finite_lr_values
  ) {
    sign(dat$logFC) * sqrt(dat$LR)
  } else if (
    n_finite_signed_stat >= minimum_finite_signed_stat_values
  ) {
    dat$signed_stat
  } else {
    sign(dat$logFC) * -log10(safe_p)
  }

  dat |>
    dplyr::filter(is.finite(collapse_rank)) |>
    dplyr::group_by(Gene) |>
    dplyr::arrange(
      dplyr::desc(abs(collapse_rank)),
      PValue,
      dplyr::desc(abs(logFC)),
      original_row,
      .by_group = TRUE
    ) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::mutate(comparison = comparison) |>
    dplyr::select(
      comparison,
      Gene,
      logFC,
      logCPM,
      PValue,
      FDR,
      LR,
      signed_stat
    )
}

heatmap_de <- collapse_de_for_heatmap(
  standardised_de_list[[contrast_key]],
  contrast_key
)

data.table::fwrite(
  heatmap_de,
  file.path(
    subdirs$tables,
    "BRAF_V600E_VS_WT_SSLHP_COLLAPSED_DE_FOR_GENE_HEATMAPS.csv"
  )
)

module_gene_sets <- list(
  "Interferon alpha response" = ifna_info$genes,
  "Inflammatory response" = inflammatory_info$genes,
  "TNFA/NFKB signaling" = tnfa_info$genes,
  "Gastric metaplasia" = metaplasia_genes
)

module_pathway_names <- c(
  "Interferon alpha response" = ifna_info$name,
  "Inflammatory response" = inflammatory_info$name,
  "TNFA/NFKB signaling" = tnfa_info$name,
  "Gastric metaplasia" = "FINAL_SSL_METAPLASIA_LIST"
)

module_order <- names(module_gene_sets)

module_display_title <- function(module_name) {
  if (identical(module_name, "Gastric metaplasia")) {
    return(metaplasia_display_title_wrapped)
  }
  module_name
}

parse_leading_edge <- function(result_df, pathway_name) {
  hit <- result_df |>
    dplyr::filter(toupper(pathway) == toupper(pathway_name)) |>
    dplyr::slice(1)

  if (
    nrow(hit) == 0L ||
    !"leadingEdge" %in% names(hit) ||
    is.na(hit$leadingEdge[1]) ||
    hit$leadingEdge[1] == ""
  ) {
    return(character(0))
  }

  unique(strsplit(hit$leadingEdge[1], ";", fixed = TRUE)[[1]])
}

module_gene_audit <- purrr::imap_dfr(
  module_gene_sets,
  function(genes, module_name) {
    genes <- unique(clean_gene_symbols(genes))
    genes <- genes[!is.na(genes) & genes != ""]

    pathway_name <- unname(module_pathway_names[[module_name]])

    if (module_name == "Gastric metaplasia") {
      leading_edge <- parse_leading_edge(
        metaplasia_fgsea_all,
        pathway_name
      )
    } else {
      leading_edge <- parse_leading_edge(
        hallmark_fgsea_all,
        pathway_name
      )
    }

    tibble::tibble(
      module = module_name,
      pathway = pathway_name,
      Gene = genes,
      signature_order = seq_along(genes)
    ) |>
      dplyr::left_join(
        heatmap_de |>
          dplyr::select(
            Gene,
            logFC,
            logCPM,
            PValue,
            FDR,
            LR,
            signed_stat
          ),
        by = "Gene"
      ) |>
      dplyr::mutate(
        present_in_DE_table = is.finite(logFC),
        expression_eligible = dplyr::if_else(
          module == "Gastric metaplasia",
          present_in_DE_table,
          present_in_DE_table &
            is.finite(logCPM) &
            logCPM >= hallmark_logcpm_cutoff
        ),
        direction = dplyr::case_when(
          !present_in_DE_table ~ "Absent",
          logFC > 0 ~ "Higher in BRAF V600E",
          logFC < 0 ~ "Higher in BRAF WT",
          TRUE ~ "No change"
        ),
        significant_up_in_V600E =
          expression_eligible &
          is.finite(FDR) &
          FDR < 0.05 &
          logFC > 0,
        significant_up_in_WT =
          expression_eligible &
          is.finite(FDR) &
          FDR < 0.05 &
          logFC < 0,
        leading_edge = Gene %in% leading_edge,
        contrast = contrast_long_label
      )
  }
) |>
  dplyr::arrange(
    factor(module, levels = module_order),
    signature_order
  )

module_gene_summary <- module_gene_audit |>
  dplyr::group_by(module, pathway) |>
  dplyr::summarise(
    n_signature_genes = dplyr::n(),
    n_present_in_DE_table = sum(present_in_DE_table, na.rm = TRUE),
    n_expression_eligible = sum(expression_eligible, na.rm = TRUE),
    n_higher_in_V600E = sum(
      expression_eligible & is.finite(logFC) & logFC > 0,
      na.rm = TRUE
    ),
    n_higher_in_WT = sum(
      expression_eligible & is.finite(logFC) & logFC < 0,
      na.rm = TRUE
    ),
    n_FDR_significant_higher_in_V600E = sum(
      significant_up_in_V600E,
      na.rm = TRUE
    ),
    n_FDR_significant_higher_in_WT = sum(
      significant_up_in_WT,
      na.rm = TRUE
    ),
    n_leading_edge = sum(leading_edge, na.rm = TRUE),
    mean_logFC_expression_eligible = if (
      any(expression_eligible & is.finite(logFC))
    ) {
      mean(logFC[expression_eligible & is.finite(logFC)], na.rm = TRUE)
    } else {
      NA_real_
    },
    median_logFC_expression_eligible = if (
      any(expression_eligible & is.finite(logFC))
    ) {
      stats::median(
        logFC[expression_eligible & is.finite(logFC)],
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    .groups = "drop"
  )

data.table::fwrite(
  module_gene_audit,
  file.path(
    subdirs$gene_summary,
    "ALL_SELECTED_MODULE_GENES_BRAF_V600E_VS_WT_SSLHP_LONG.csv"
  )
)

data.table::fwrite(
  module_gene_audit |>
    dplyr::filter(leading_edge),
  file.path(
    subdirs$gene_summary,
    "SELECTED_PATHWAY_LEADING_EDGE_GENES_BRAF_V600E_VS_WT_SSLHP.csv"
  )
)

data.table::fwrite(
  module_gene_summary,
  file.path(
    subdirs$gene_summary,
    "SELECTED_MODULE_GENE_DIRECTION_SUMMARY_BRAF_V600E_VS_WT_SSLHP.csv"
  )
)

calc_enrichment_curve <- function(ranks, pathway_genes) {
  ranks <- sort(ranks, decreasing = TRUE)

  pathway_genes <- unique(clean_gene_symbols(pathway_genes))
  pathway_genes <- pathway_genes[
    !is.na(pathway_genes) & pathway_genes != ""
  ]
  pathway_genes <- intersect(pathway_genes, names(ranks))

  n_all <- length(ranks)
  n_hits <- length(pathway_genes)

  if (n_hits < 1L || n_all <= n_hits) {
    return(tibble::tibble())
  }

  hits <- names(ranks) %in% pathway_genes
  weights <- abs(ranks)^fgsea_gsea_param
  norm_hit <- sum(weights[hits])

  if (!is.finite(norm_hit) || norm_hit == 0) {
    weights <- rep(1, length(ranks))
    norm_hit <- sum(weights[hits])
  }

  running_es <- cumsum(ifelse(hits, weights / norm_hit, 0)) -
    cumsum(ifelse(!hits, 1 / (n_all - n_hits), 0))

  tibble::tibble(
    position = seq_along(ranks),
    position_percent = 100 * position / n_all,
    running_es = running_es,
    hit = hits,
    gene = names(ranks)
  )
}

get_fgsea_result <- function(result_df, pathway_name) {
  hit <- result_df |>
    dplyr::filter(toupper(pathway) == toupper(pathway_name)) |>
    dplyr::slice(1)

  if (nrow(hit) == 0L) {
    stop("FGSEA result was not found for pathway: ", pathway_name)
  }

  hit
}

make_single_gsea_curve <- function(
    rank_object,
    pathway_genes,
    result_row,
    title_text,
    p_label
) {
  curve <- calc_enrichment_curve(
    rank_object$ranks,
    pathway_genes
  )

  if (nrow(curve) == 0L) {
    return(make_empty_plot(
      title_text,
      "No overlapping pathway genes"
    ))
  }

  p_value <- if (identical(p_label, "Hallmark FDR")) {
    result_row$padj[1]
  } else {
    result_row$pval[1]
  }

  annotation <- paste0(
    "NES = ",
    sprintf("%.2f", result_row$NES[1]),
    "\n",
    p_label,
    " = ",
    format_numeric_p(p_value),
    "\nGenes = ",
    result_row$size[1]
  )

  ggplot2::ggplot(
    curve,
    ggplot2::aes(
      x = position_percent,
      y = running_es
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = thin_line,
      colour = "grey45"
    ) +
    ggplot2::geom_line(
      colour = curve_colour,
      linewidth = 0.78
    ) +
    ggplot2::geom_rug(
      data = curve |>
        dplyr::filter(hit),
      ggplot2::aes(x = position_percent),
      sides = "b",
      inherit.aes = FALSE,
      colour = curve_colour,
      linewidth = 0.18
    ) +
    ggplot2::annotate(
      "text",
      x = 98,
      y = Inf,
      label = annotation,
      hjust = 1,
      vjust = 1.25,
      family = font_family,
      size = geom_text_size_from_pt(5.3)
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 25, 50, 75, 100),
      limits = c(0, 100),
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.06, 0.18))
    ) +
    ggplot2::labs(
      title = title_text,
      x = "Ranked gene list position (%)",
      y = "Running enrichment score"
    ) +
    nature_theme() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(
        family = font_family,
        face = "bold",
        size = plot_title_size,
        hjust = 0
      )
    )
}

curve_specifications <- list(
  IFNA = list(
    pathway = ifna_info$name,
    genes = ifna_info$genes,
    title = paste0(
      "Interferon alpha response | ",
      contrast_long_label
    ),
    family = "HALLMARK"
  ),
  INFLAMMATORY = list(
    pathway = inflammatory_info$name,
    genes = inflammatory_info$genes,
    title = paste0(
      "Inflammatory response | ",
      contrast_long_label
    ),
    family = "HALLMARK"
  ),
  TNFA_NFKB = list(
    pathway = tnfa_info$name,
    genes = tnfa_info$genes,
    title = paste0(
      "TNFA/NFKB signaling | ",
      contrast_long_label
    ),
    family = "HALLMARK"
  ),
  GASTRIC_METAPLASIA = list(
    pathway = "FINAL_SSL_METAPLASIA_LIST",
    genes = metaplasia_genes,
    title = paste0(
      metaplasia_display_title_wrapped,
      "\n",
      contrast_long_label
    ),
    family = "CUSTOM"
  )
)

gsea_curve_plots <- list()

for (curve_key in names(curve_specifications)) {
  spec <- curve_specifications[[curve_key]]

  if (spec$family == "HALLMARK") {
    result_row <- get_fgsea_result(
      hallmark_fgsea_all,
      spec$pathway
    )
    rank_object <- hallmark_rank_objects[[contrast_key]]
    p_label <- "Hallmark FDR"
    output_directory <- subdirs$hallmark_curves
  } else {
    result_row <- get_fgsea_result(
      metaplasia_fgsea_all,
      spec$pathway
    )
    rank_object <- metaplasia_rank_objects[[contrast_key]]
    p_label <- "Nominal p"
    output_directory <- subdirs$metaplasia_curves
  }

  p <- make_single_gsea_curve(
    rank_object = rank_object,
    pathway_genes = spec$genes,
    result_row = result_row,
    title_text = spec$title,
    p_label = p_label
  )

  gsea_curve_plots[[curve_key]] <- p

  save_plot_set(
    p,
    output_directory,
    paste0(
      "BRAF_V600E_VS_WT_SSLHP_",
      curve_key,
      "_GSEA_CURVE"
    ),
    92,
    72
  )
}

p_gsea_curves_2x2 <- patchwork::wrap_plots(
  gsea_curve_plots[c(
    "IFNA",
    "INFLAMMATORY",
    "TNFA_NFKB",
    "GASTRIC_METAPLASIA"
  )],
  ncol = 2,
  guides = "keep"
) +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = paste0(
      "Selected GSEA results: ",
      contrast_long_label
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = font_family,
        face = "bold",
        size = 8,
        hjust = 0
      ),
      plot.tag = ggplot2::element_text(
        family = font_family,
        face = "bold",
        size = panel_tag_size
      )
    )
  )

save_plot_set(
  p_gsea_curves_2x2,
  subdirs$combined,
  "FINAL_BRAF_V600E_VS_WT_SSLHP_SELECTED_GSEA_CURVES_2X2",
  183,
  150
)

short_hallmark_label <- function(x) {
  x_upper <- toupper(as.character(x))

  fallback <- x |>
    as.character() |>
    stringr::str_remove("^HALLMARK_") |>
    stringr::str_replace_all("_", " ") |>
    stringr::str_to_lower() |>
    stringr::str_to_sentence()

  dplyr::case_when(
    x_upper == "HALLMARK_INTERFERON_ALPHA_RESPONSE" ~
      "Interferon alpha response",
    x_upper == "HALLMARK_INTERFERON_GAMMA_RESPONSE" ~
      "Interferon gamma response",
    x_upper == "HALLMARK_TNFA_SIGNALING_VIA_NFKB" ~
      "TNFA/NFKB signaling",
    x_upper == "HALLMARK_INFLAMMATORY_RESPONSE" ~
      "Inflammatory response",
    x_upper == "HALLMARK_IL6_JAK_STAT3_SIGNALING" ~
      "IL-6/JAK/STAT3 signaling",
    x_upper == "HALLMARK_IL2_STAT5_SIGNALING" ~
      "IL-2/STAT5 signaling",
    x_upper == "HALLMARK_WNT_BETA_CATENIN_SIGNALING" ~
      "Wnt/beta-catenin signaling",
    x_upper == "HALLMARK_OXIDATIVE_PHOSPHORYLATION" ~
      "Oxidative phosphorylation",
    x_upper == "HALLMARK_UNFOLDED_PROTEIN_RESPONSE" ~
      "Unfolded protein response",
    x_upper == "HALLMARK_GLYCOLYSIS" ~ "Glycolysis",
    x_upper == "HALLMARK_COMPLEMENT" ~ "Complement",
    TRUE ~ fallback
  )
}

make_hallmark_nes_heatmap <- function(
    pathways = NULL,
    top_n = NULL,
    title_text = NULL
) {
  dat <- hallmark_fgsea_all |>
    dplyr::filter(is.finite(NES)) |>
    dplyr::group_by(pathway) |>
    dplyr::arrange(
      padj,
      pval,
      dplyr::desc(abs(NES)),
      .by_group = TRUE
    ) |>
    dplyr::slice(1) |>
    dplyr::ungroup()

  if (!is.null(pathways)) {
    dat <- dat |>
      dplyr::filter(pathway %in% pathways)
  } else if (!is.null(top_n)) {
    dat <- dat |>
      dplyr::arrange(
        padj,
        dplyr::desc(abs(NES)),
        pathway
      ) |>
      dplyr::slice_head(n = top_n)
  }

  if (nrow(dat) == 0L) {
    return(make_empty_plot(
      "Hallmark GSEA",
      "No requested Hallmark pathways were found"
    ))
  }

  pathway_order <- dat |>
    dplyr::arrange(
      dplyr::desc(abs(NES)),
      padj,
      pathway
    ) |>
    dplyr::pull(pathway)

  plot_df <- dat |>
    dplyr::mutate(
      comparison_display = factor(
        contrast_display_label,
        levels = contrast_display_label
      ),
      pathway_label_text = short_hallmark_label(pathway),
      pathway_label = factor(
        pathway_label_text,
        levels = rev(short_hallmark_label(pathway_order))
      ),
      significance = significance_stars(padj)
    )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = comparison_display,
      y = pathway_label,
      fill = NES
    )
  ) +
    ggplot2::geom_tile(
      colour = tile_border_colour,
      linewidth = tile_border_linewidth,
      width = tile_fraction,
      height = tile_fraction
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = significance),
      family = font_family,
      fontface = "bold",
      colour = "white",
      size = geom_text_size_from_pt(shared_label_pt)
    ) +
    ggplot2::scale_x_discrete(
      drop = FALSE,
      expand = ggplot2::expansion(add = 0.05)
    ) +
    ggplot2::scale_y_discrete(
      drop = FALSE,
      labels = function(x) stringr::str_wrap(x, width = 24),
      expand = ggplot2::expansion(add = 0.08)
    ) +
    ggplot2::scale_fill_gradient2(
      low = negative_colour,
      mid = mid_colour,
      high = positive_colour,
      midpoint = 0,
      limits = c(-nes_colour_limit, nes_colour_limit),
      oob = scales::squish,
      na.value = "grey95",
      name = "NES",
      guide = vertical_heatmap_guide()
    ) +
    ggplot2::coord_fixed(ratio = 1, clip = "off") +
    ggplot2::labs(title = title_text, x = NULL, y = NULL) +
    nature_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        family = font_family,
        size = shared_label_pt,
        angle = x_label_angle,
        hjust = 1,
        vjust = 1,
        colour = "black"
      ),
      axis.text.y = ggplot2::element_text(
        family = font_family,
        size = shared_label_pt,
        lineheight = 0.9,
        colour = "black"
      ),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      legend.direction = "vertical",
      plot.margin = ggplot2::margin(2, 2, 2, 2, unit = "mm")
    )
}

key_hallmark_pathways <- c(
  ifna_info$name,
  inflammatory_info$name,
  tnfa_info$name,
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING"
)

selected_three_hallmark_pathways <- c(
  ifna_info$name,
  inflammatory_info$name,
  tnfa_info$name
)

p_hallmark_all <- make_hallmark_nes_heatmap(
  pathways = sort(unique(hallmark_fgsea_all$pathway)),
  title_text = NULL
)

p_hallmark_top20 <- make_hallmark_nes_heatmap(
  top_n = 20L,
  title_text = NULL
)

p_hallmark_key <- make_hallmark_nes_heatmap(
  pathways = key_hallmark_pathways,
  title_text = NULL
)

p_hallmark_selected_three <- make_hallmark_nes_heatmap(
  pathways = selected_three_hallmark_pathways,
  title_text = NULL
)

save_plot_set(
  p_hallmark_all,
  subdirs$hallmark_heatmaps,
  "FINAL_ALL_HALLMARK_NES_HEATMAP_BRAF_V600E_VS_WT_SSLHP",
  95,
  220
)

save_plot_set(
  p_hallmark_top20,
  subdirs$hallmark_heatmaps,
  "FINAL_TOP20_HALLMARK_NES_HEATMAP_BRAF_V600E_VS_WT_SSLHP",
  95,
  135
)

save_plot_set(
  p_hallmark_key,
  subdirs$hallmark_heatmaps,
  "FINAL_KEY_HALLMARK_NES_HEATMAP_BRAF_V600E_VS_WT_SSLHP",
  95,
  92
)

save_plot_set(
  p_hallmark_selected_three,
  subdirs$hallmark_heatmaps,
  "FINAL_IFNA_INFLAMMATORY_TNFA_NES_HEATMAP_BRAF_V600E_VS_WT_SSLHP",
  95,
  70
)

make_metaplasia_nes_heatmap <- function() {
  plot_df <- metaplasia_fgsea_all |>
    dplyr::slice(1) |>
    dplyr::mutate(
      comparison_display = factor(
        contrast_display_label,
        levels = contrast_display_label
      ),
      pathway_display = factor(
        metaplasia_display_title_wrapped,
        levels = metaplasia_display_title_wrapped
      ),
      nominal_significance = significance_stars(pval)
    )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = comparison_display,
      y = pathway_display,
      fill = NES
    )
  ) +
    ggplot2::geom_tile(
      colour = tile_border_colour,
      linewidth = tile_border_linewidth,
      width = tile_fraction,
      height = tile_fraction
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = nominal_significance),
      family = font_family,
      fontface = "bold",
      colour = "white",
      size = geom_text_size_from_pt(shared_label_pt)
    ) +
    ggplot2::scale_fill_gradient2(
      low = negative_colour,
      mid = mid_colour,
      high = positive_colour,
      midpoint = 0,
      limits = c(-nes_colour_limit, nes_colour_limit),
      oob = scales::squish,
      na.value = "grey95",
      name = "NES",
      guide = vertical_heatmap_guide()
    ) +
    ggplot2::coord_fixed(ratio = 1, clip = "off") +
    ggplot2::labs(
      title = metaplasia_display_title_wrapped,
      subtitle = "GSEA; stars use nominal p",
      x = NULL,
      y = NULL
    ) +
    nature_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        family = font_family,
        size = shared_label_pt,
        angle = x_label_angle,
        hjust = 1,
        vjust = 1,
        colour = "black"
      ),
      axis.text.y = ggplot2::element_text(
        family = font_family,
        size = shared_label_pt,
        colour = "black"
      ),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      legend.position = "right",
      legend.direction = "vertical"
    )
}

p_metaplasia_nes <- make_metaplasia_nes_heatmap()

save_plot_set(
  p_metaplasia_nes,
  subdirs$hallmark_heatmaps,
  "FINAL_GASTRIC_METAPLASIA_NES_HEATMAP_BRAF_V600E_VS_WT_SSLHP_NOMINAL_P",
  95,
  58
)

p_selected_pathway_nes <- p_hallmark_selected_three /
  p_metaplasia_nes +
  patchwork::plot_layout(heights = c(1.2, 0.8), guides = "collect") +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = paste0(
      "Selected pathway NES: ",
      contrast_long_label
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = font_family,
        face = "bold",
        size = 8,
        hjust = 0
      ),
      plot.tag = ggplot2::element_text(
        family = font_family,
        face = "bold",
        size = panel_tag_size
      )
    )
  ) &
  ggplot2::theme(legend.position = "right")

save_plot_set(
  p_selected_pathway_nes,
  subdirs$combined,
  "FINAL_SELECTED_PATHWAY_NES_BRAF_V600E_VS_WT_SSLHP",
  110,
  125
)

make_gene_heatmap_df <- function(module_name, genes) {
  genes <- unique(clean_gene_symbols(genes))
  genes <- genes[!is.na(genes) & genes != ""]

  if (length(genes) == 0L) {
    return(tibble::tibble())
  }

  tibble::tibble(
    module = module_name,
    Gene = genes,
    gene_order = seq_along(genes)
  ) |>
    dplyr::left_join(
      heatmap_de |>
        dplyr::select(
          Gene,
          logFC,
          logCPM,
          PValue,
          FDR
        ),
      by = "Gene"
    ) |>
    dplyr::mutate(
      comparison = contrast_key,
      comparison_display = factor(
        contrast_display_label,
        levels = contrast_display_label
      ),
      Gene_display = factor(
        Gene,
        levels = rev(genes)
      ),
      logFC_colour = pmax(
        pmin(logFC, logfc_colour_limit),
        -logfc_colour_limit
      )
    )
}

make_gene_heatmap <- function(
    module_name,
    genes,
    title_text = module_display_title(module_name),
    show_legend = TRUE
) {
  genes <- unique(clean_gene_symbols(genes))
  genes <- genes[!is.na(genes) & genes != ""]

  if (length(genes) == 0L) {
    return(make_empty_plot(
      title_text,
      paste0(
        "No genes passed logCPM >= ",
        module_reference_logcpm_cutoff,
        " and logFC > ",
        module_reference_logfc_cutoff
      )
    ))
  }

  plot_df <- make_gene_heatmap_df(module_name, genes)
  legend_position <- if (show_legend) "right" else "none"

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = comparison_display,
      y = Gene_display,
      fill = logFC_colour
    )
  ) +
    ggplot2::geom_tile(
      colour = tile_border_colour,
      linewidth = tile_border_linewidth,
      width = tile_fraction,
      height = tile_fraction
    ) +
    ggplot2::scale_fill_gradient2(
      low = negative_colour,
      mid = mid_colour,
      high = positive_colour,
      midpoint = 0,
      limits = c(-logfc_colour_limit, logfc_colour_limit),
      breaks = c(-logfc_colour_limit, 0, logfc_colour_limit),
      labels = c(
        paste0("<= -", logfc_colour_limit),
        "0",
        paste0(">= +", logfc_colour_limit)
      ),
      oob = scales::squish,
      na.value = "grey90",
      name = "log2FC",
      guide = vertical_heatmap_guide()
    ) +
    ggplot2::scale_x_discrete(
      drop = FALSE,
      expand = ggplot2::expansion(add = 0.05)
    ) +
    ggplot2::scale_y_discrete(
      drop = FALSE,
      expand = ggplot2::expansion(add = 0.05)
    ) +
    ggplot2::coord_fixed(ratio = 1, clip = "off") +
    ggplot2::labs(
      title = title_text,
      x = NULL,
      y = NULL
    ) +
    nature_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        family = font_family,
        size = shared_label_pt,
        angle = x_label_angle,
        hjust = 1,
        vjust = 1,
        colour = "black"
      ),
      axis.text.y = ggplot2::element_text(
        family = font_family,
        size = shared_label_pt,
        colour = "black"
      ),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = legend_position,
      legend.direction = "vertical",
      plot.margin = ggplot2::margin(2, 2, 2, 2, unit = "mm")
    )
}

full_module_gene_lists <- list()
full_module_gene_qc <- list()

for (module_name in module_order) {
  module_qc <- module_gene_audit |>
    dplyr::filter(module == module_name)

  if (
    module_name != "Gastric metaplasia" &&
    !full_hallmark_include_genes_absent_from_de_table
  ) {
    module_qc <- module_qc |>
      dplyr::filter(present_in_DE_table)
  }

  if (module_name == "Gastric metaplasia") {
    module_qc <- module_qc |>
      dplyr::arrange(signature_order)
  } else {
    module_qc <- module_qc |>
      dplyr::arrange(
        dplyr::desc(abs(logFC)),
        signature_order,
        Gene
      )
  }

  full_module_gene_lists[[module_name]] <- module_qc$Gene
  full_module_gene_qc[[module_name]] <- module_qc |>
    dplyr::mutate(full_heatmap_order = dplyr::row_number())
}

full_module_gene_qc_all <- dplyr::bind_rows(full_module_gene_qc)

data.table::fwrite(
  full_module_gene_qc_all,
  file.path(
    subdirs$qc,
    "FULL_SIGNATURE_GENE_COVERAGE_AND_ORDER_QC_BRAF_V600E_VS_WT_SSLHP.csv"
  )
)

full_module_heatmap_long <- purrr::imap_dfr(
  full_module_gene_lists,
  function(genes, module_name) {
    make_gene_heatmap_df(module_name, genes)
  }
)

data.table::fwrite(
  full_module_heatmap_long,
  file.path(
    subdirs$tables,
    "FULL_SIGNATURE_LOG2FC_HEATMAP_LONG_BRAF_V600E_VS_WT_SSLHP.csv"
  )
)

full_module_plots <- list()

for (module_name in module_order) {
  genes <- full_module_gene_lists[[module_name]]

  p_one <- make_gene_heatmap(
    module_name = module_name,
    genes = genes,
    title_text = module_display_title(module_name),
    show_legend = TRUE
  )

  full_module_plots[[module_name]] <- p_one

  n_genes <- length(genes)
  height_mm <- max(
    85,
    min(340, 35 + 2.0 * n_genes)
  )

  save_plot_set(
    p_one,
    subdirs$full_modules,
    paste0(
      "FULL_SIGNATURE_",
      safe_filename(module_name),
      "_LOG2FC_HEATMAP_BRAF_V600E_VS_WT_SSLHP"
    ),
    85,
    height_mm
  )
}

if ("Gastric metaplasia" %in% names(full_module_plots)) {
  save_plot_set(
    full_module_plots[["Gastric metaplasia"]],
    subdirs$full_modules,
    "FINAL_24_GENE_GASTRIC_METAPLASIA_LOG2FC_HEATMAP_BRAF_V600E_VS_WT_SSLHP",
    85,
    100
  )
}

selected_candidate_qc <- module_gene_audit |>
  dplyr::mutate(
    selection_rule = dplyr::if_else(
      module == "Gastric metaplasia",
      "All final 24 gastric-metaplasia genes; no logCPM or logFC threshold",
      paste0(
        "logCPM >= ",
        module_reference_logcpm_cutoff,
        " and logFC > ",
        module_reference_logfc_cutoff,
        "; top ",
        module_max_genes,
        " by decreasing logFC"
      )
    ),
    passes_selection_filter = dplyr::case_when(
      module == "Gastric metaplasia" ~ TRUE,
      TRUE ~ present_in_DE_table &
        is.finite(logCPM) &
        logCPM >= module_reference_logcpm_cutoff &
        is.finite(logFC) &
        logFC > module_reference_logfc_cutoff
    )
  )

selected_module_genes <- purrr::map_dfr(
  module_order,
  function(module_name) {
    dat <- selected_candidate_qc |>
      dplyr::filter(module == module_name)

    if (module_name == "Gastric metaplasia") {
      dat |>
        dplyr::arrange(signature_order) |>
        dplyr::mutate(selection_rank = dplyr::row_number())
    } else {
      dat |>
        dplyr::filter(passes_selection_filter) |>
        dplyr::arrange(
          dplyr::desc(logFC),
          FDR,
          Gene
        ) |>
        dplyr::slice_head(n = module_max_genes) |>
        dplyr::mutate(selection_rank = dplyr::row_number())
    }
  }
)

missing_selected_modules <- setdiff(
  setdiff(module_order, "Gastric metaplasia"),
  unique(selected_module_genes$module)
)

if (length(missing_selected_modules) > 0L) {
  warning(
    "No genes passed the strict BRAF V600E-up selection rule for: ",
    paste(missing_selected_modules, collapse = ", "),
    ". Full-signature heatmaps are still written; the threshold was not ",
    "silently relaxed.",
    call. = FALSE
  )
}

data.table::fwrite(
  selected_candidate_qc,
  file.path(
    subdirs$qc,
    "PUBLICATION_SELECTED_MODULE_CANDIDATE_QC_BRAF_V600E_VS_WT_SSLHP.csv"
  )
)

data.table::fwrite(
  selected_module_genes,
  file.path(
    subdirs$tables,
    "PUBLICATION_SELECTED_GENES_BRAF_V600E_VS_WT_SSLHP.csv"
  )
)

selected_gene_lists <- split(
  selected_module_genes$Gene,
  selected_module_genes$module
)

for (module_name in module_order) {
  if (!module_name %in% names(selected_gene_lists)) {
    selected_gene_lists[[module_name]] <- character(0)
  }
}
selected_gene_lists <- selected_gene_lists[module_order]

selected_module_heatmap_long <- purrr::imap_dfr(
  selected_gene_lists,
  function(genes, module_name) {
    if (length(genes) == 0L) return(tibble::tibble())
    make_gene_heatmap_df(module_name, genes)
  }
)

data.table::fwrite(
  selected_module_heatmap_long,
  file.path(
    subdirs$tables,
    "PUBLICATION_SELECTED_LOG2FC_HEATMAP_LONG_BRAF_V600E_VS_WT_SSLHP.csv"
  )
)

selected_module_plots <- list()

for (module_name in module_order) {
  genes <- selected_gene_lists[[module_name]]

  p_one <- make_gene_heatmap(
    module_name = module_name,
    genes = genes,
    title_text = module_display_title(module_name),
    show_legend = TRUE
  )

  selected_module_plots[[module_name]] <- p_one

  selected_height <- if (module_name == "Gastric metaplasia") {
    100
  } else {
    82
  }

  save_plot_set(
    p_one,
    subdirs$selected_modules,
    paste0(
      "SELECTED_",
      safe_filename(module_name),
      "_LOG2FC_HEATMAP_BRAF_V600E_VS_WT_SSLHP"
    ),
    82,
    selected_height
  )
}

selected_module_plots_combined <- list()

for (i in seq_along(module_order)) {
  module_name <- module_order[i]
  selected_module_plots_combined[[module_name]] <- make_gene_heatmap(
    module_name = module_name,
    genes = selected_gene_lists[[module_name]],
    title_text = module_display_title(module_name),
    show_legend = i == length(module_order)
  )
}

p_selected_modules_combined <- patchwork::wrap_plots(
  selected_module_plots_combined[module_order],
  ncol = 4
) +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = paste0(
      "Selected gene-level log2FC heatmaps: ",
      contrast_long_label
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = font_family,
        face = "bold",
        size = 8,
        hjust = 0
      ),
      plot.tag = ggplot2::element_text(
        family = font_family,
        face = "bold",
        size = panel_tag_size
      )
    )
  )

save_plot_set(
  p_selected_modules_combined,
  subdirs$combined,
  "FINAL_SELECTED_GENE_MODULE_HEATMAPS_BRAF_V600E_VS_WT_SSLHP",
  220,
  145
)

four_module_get_msig_collection <- function(
    collection,
    subcollection = NULL
) {
  if (is.null(subcollection)) {
    attempts <- list(
      function() {
        msigdbr::msigdbr(
          species = "Homo sapiens",
          collection = collection
        )
      },
      function() {
        msigdbr::msigdbr(
          species = "Homo sapiens",
          category = collection
        )
      }
    )
  } else {
    attempts <- list(
      function() {
        msigdbr::msigdbr(
          species = "Homo sapiens",
          collection = collection,
          subcollection = subcollection
        )
      },
      function() {
        msigdbr::msigdbr(
          species = "Homo sapiens",
          category = collection,
          subcategory = subcollection
        )
      },
      function() {
        msigdbr::msigdbr(
          species = "Homo sapiens",
          collection = collection
        )
      },
      function() {
        msigdbr::msigdbr(
          species = "Homo sapiens",
          category = collection
        )
      }
    )
  }

  for (attempt in attempts) {
    result <- tryCatch(
      suppressWarnings(attempt()),
      error = function(e) NULL
    )

    if (!is.null(result) && nrow(result) > 0L) {
      return(as.data.frame(result))
    }
  }

  stop(
    "Could not load MSigDB collection ", collection,
    if (!is.null(subcollection)) {
      paste0(" / ", subcollection)
    } else {
      ""
    },
    " through msigdbr."
  )
}

four_module_split_msig_pathways <- function(
    msig_raw,
    collection_label
) {
  gene_col <- detect_col(
    msig_raw,
    c(
      "gene_symbol", "human_gene_symbol", "db_gene_symbol",
      "symbol"
    ),
    required = TRUE,
    label = paste0(collection_label, " gene-symbol column")
  )

  pathway_col <- detect_col(
    msig_raw,
    c(
      "gs_name", "geneset", "gene_set_name", "gs_exact_source"
    ),
    required = TRUE,
    label = paste0(collection_label, " pathway-name column")
  )

  clean_table <- data.frame(
    pathway = as.character(msig_raw[[pathway_col]]),
    gene = clean_gene_symbols(msig_raw[[gene_col]]),
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(
      !is.na(pathway),
      pathway != "",
      !is.na(gene),
      gene != ""
    ) |>
    dplyr::distinct(pathway, gene)

  pathways <- split(clean_table$gene, clean_table$pathway)
  lapply(
    pathways,
    function(x) sort(unique(clean_gene_symbols(x)))
  )
}

four_module_resolve_pathway <- function(
    named_pathways,
    requested_pathway,
    loose_pattern
) {
  available <- names(named_pathways)
  exact <- available[
    toupper(available) == toupper(requested_pathway)
  ]

  if (length(exact) > 0L) {
    resolved <- exact[1]
  } else {
    candidates <- available[
      stringr::str_detect(
        toupper(available),
        toupper(loose_pattern)
      )
    ]

    if (length(candidates) == 0L) {
      stop(
        "Could not resolve required pathway '",
        requested_pathway,
        "'."
      )
    }

    candidates <- candidates[
      order(nchar(candidates), candidates)
    ]
    resolved <- candidates[1]

    warning(
      "Exact pathway '", requested_pathway,
      "' was not found. Using '", resolved, "'.",
      call. = FALSE
    )
  }

  genes <- clean_gene_symbols(named_pathways[[resolved]])
  genes <- sort(unique(genes[!is.na(genes) & genes != ""]))

  list(
    requested_pathway = requested_pathway,
    resolved_pathway = resolved,
    genes = genes
  )
}

message(
  "Loading MSigDB C2 Reactome and C5 GO molecular-function sets ",
  "for the four-module heatmap panel..."
)

msig_c2_four_module <- four_module_get_msig_collection(
  collection = "C2",
  subcollection = "CP:REACTOME"
)

msig_c5_four_module <- four_module_get_msig_collection(
  collection = "C5",
  subcollection = "GO:MF"
)

reactome_pathways_four_module <- four_module_split_msig_pathways(
  msig_c2_four_module,
  "MSigDB C2 Reactome"
)

go_mf_pathways_four_module <- four_module_split_msig_pathways(
  msig_c5_four_module,
  "MSigDB C5 GO molecular function"
)

chemokine_activity_info_four_module <- four_module_resolve_pathway(
  named_pathways = go_mf_pathways_four_module,
  requested_pathway = "GOMF_CHEMOKINE_ACTIVITY",
  loose_pattern = "CHEMOKINE_ACTIVITY"
)

antimicrobial_info_four_module <- four_module_resolve_pathway(
  named_pathways = reactome_pathways_four_module,
  requested_pathway = "REACTOME_ANTIMICROBIAL_PEPTIDES",
  loose_pattern = "ANTIMICROBIAL_PEPTIDES"
)

canonical_chemokine_regex_four_module <- paste0(
  "^(CCL[0-9]+[A-Z0-9]*|",
  "CXCL[0-9]+[A-Z0-9]*|",
  "CX3CL[0-9]+[A-Z0-9]*|",
  "XCL[0-9]+[A-Z0-9]*|",
  "PF4|PF4V1|PPBP)$"
)

chemokine_gene_audit_four_module <- tibble::tibble(
  Gene = chemokine_activity_info_four_module$genes,
  source_gene_set =
    chemokine_activity_info_four_module$resolved_pathway,
  retained_for_heatmap = stringr::str_detect(
    Gene,
    canonical_chemokine_regex_four_module
  )
)

chemokine_genes_four_module <- chemokine_gene_audit_four_module |>
  dplyr::filter(retained_for_heatmap) |>
  dplyr::arrange(Gene) |>
  dplyr::pull(Gene)

if (length(chemokine_genes_four_module) == 0L) {
  stop(
    "GOMF_CHEMOKINE_ACTIVITY was resolved, but no canonical ",
    "chemokine ligands remained after filtering."
  )
}

four_module_gene_sets <- list(
  "Interferon alpha response" = ifna_info$genes,
  "Chemokines" = chemokine_genes_four_module,
  "Antimicrobial peptides" = antimicrobial_info_four_module$genes,
  "TNFA/NFKB signaling" = tnfa_info$genes
)

four_module_order <- names(four_module_gene_sets)

four_module_source_table <- tibble::tibble(
  module = four_module_order,
  requested_pathway = c(
    ifna_info$name,
    chemokine_activity_info_four_module$requested_pathway,
    antimicrobial_info_four_module$requested_pathway,
    tnfa_info$name
  ),
  resolved_pathway = c(
    ifna_info$name,
    chemokine_activity_info_four_module$resolved_pathway,
    antimicrobial_info_four_module$resolved_pathway,
    tnfa_info$name
  ),
  source = c(
    "MSigDB Hallmark",
    paste0(
      "MSigDB C5 GOMF_CHEMOKINE_ACTIVITY; canonical ",
      "CCL/CXCL/CX3CL/XCL/PF4-family ligands only"
    ),
    "MSigDB C2 Reactome",
    "MSigDB Hallmark"
  ),
  n_candidate_genes = vapply(
    four_module_gene_sets,
    length,
    integer(1)
  )
)

data.table::fwrite(
  chemokine_gene_audit_four_module,
  file.path(
    subdirs$qc,
    paste0(
      "FOUR_MODULE_GOMF_CHEMOKINE_ACTIVITY_RAW_AND_",
      "CANONICAL_RETAINED_GENES.csv"
    )
  )
)

data.table::fwrite(
  tibble::tibble(
    Gene = antimicrobial_info_four_module$genes,
    source_gene_set =
      antimicrobial_info_four_module$resolved_pathway
  ),
  file.path(
    subdirs$qc,
    "FOUR_MODULE_REACTOME_ANTIMICROBIAL_PEPTIDE_GENES.csv"
  )
)

data.table::fwrite(
  four_module_source_table,
  file.path(
    subdirs$qc,
    "FOUR_MODULE_PANEL_GENE_SET_SOURCES.csv"
  )
)

four_module_candidate_qc <- purrr::imap_dfr(
  four_module_gene_sets,
  function(genes, module_name) {
    genes <- clean_gene_symbols(genes)
    genes <- unique(genes[!is.na(genes) & genes != ""])

    tibble::tibble(
      module = module_name,
      Gene = genes,
      source_order = seq_along(genes)
    ) |>
      dplyr::left_join(
        heatmap_de |>
          dplyr::select(
            Gene,
            logFC,
            logCPM,
            PValue,
            FDR
          ),
        by = "Gene"
      ) |>
      dplyr::mutate(
        present_in_DE_table = is.finite(logFC),
        passes_selection_filter =
          present_in_DE_table &
          is.finite(logCPM) &
          logCPM >= module_reference_logcpm_cutoff &
          is.finite(logFC) &
          abs(logFC) > module_reference_logfc_cutoff,
        selection_rule = paste0(
          "BRAF organoid-matched rule: logCPM >= ",
          module_reference_logcpm_cutoff,
          " and abs(logFC) > ",
          module_reference_logfc_cutoff,
          "; top ",
          module_max_genes,
          " by decreasing abs(logFC)"
        )
      )
  }
)

four_module_selected_genes <- purrr::map_dfr(
  four_module_order,
  function(module_name) {
    four_module_candidate_qc |>
      dplyr::filter(
        module == module_name,
        passes_selection_filter
      ) |>
      dplyr::arrange(
        dplyr::desc(abs(logFC)),
        FDR,
        Gene
      ) |>
      dplyr::slice_head(n = module_max_genes) |>
      dplyr::mutate(
        selection_rank = dplyr::row_number()
      )
  }
)

missing_four_module_blocks <- setdiff(
  four_module_order,
  unique(four_module_selected_genes$module)
)

if (length(missing_four_module_blocks) > 0L) {
  warning(
    "No genes passed the organoid-matched BRAF four-module rule for: ",
    paste(missing_four_module_blocks, collapse = ", "),
    ". The corresponding module remains visible as a grey placeholder; ",
    "the threshold was not relaxed.",
    call. = FALSE
  )
}

data.table::fwrite(
  four_module_candidate_qc,
  file.path(
    subdirs$qc,
    paste0(
      "FOUR_MODULE_PANEL_BRAF_CANDIDATE_QC_",
      "LOGCPM_GE_2_ABS_LOGFC_GT_1.csv"
    )
  )
)

data.table::fwrite(
  four_module_selected_genes,
  file.path(
    subdirs$tables,
    paste0(
      "FOUR_MODULE_PANEL_BRAF_SELECTED_GENES_",
      "LOGCPM_GE_2_ABS_LOGFC_GT_1_TOP12.csv"
    )
  )
)

four_module_selected_gene_lists <- split(
  four_module_selected_genes$Gene,
  four_module_selected_genes$module
)

for (module_name in four_module_order) {
  if (!module_name %in% names(four_module_selected_gene_lists)) {
    four_module_selected_gene_lists[[module_name]] <- character(0)
  }
}

four_module_selected_gene_lists <-
  four_module_selected_gene_lists[four_module_order]

make_four_module_braf_heatmap_df <- function(
    module_name,
    genes
) {
  genes <- clean_gene_symbols(genes)
  genes <- unique(genes[!is.na(genes) & genes != ""])

  is_placeholder <- length(genes) == 0L
  if (is_placeholder) genes <- "No genes passed"

  gene_index <- tibble::tibble(
    Gene = genes,
    gene_order = seq_along(genes)
  )

  gene_index |>
    dplyr::left_join(
      heatmap_de |>
        dplyr::select(
          Gene,
          logFC,
          logCPM,
          PValue,
          FDR
        ),
      by = "Gene"
    ) |>
    dplyr::mutate(
      module = module_name,
      comparison = contrast_key,
      comparison_display = factor(
        contrast_display_label,
        levels = contrast_display_label
      ),
      placeholder = is_placeholder,
      logFC_colour = pmax(
        pmin(logFC, logfc_colour_limit),
        -logfc_colour_limit
      )
    )
}

build_four_module_braf_one_row <- function(
    df_list,
    selected_modules = NULL,
    force_block_height = NULL,
    show_legend = TRUE,
    legend_title = "log2FC"
) {
  module_names <- four_module_order[
    four_module_order %in% names(df_list)
  ]

  if (!is.null(selected_modules)) {
    module_names <- module_names[
      module_names %in% selected_modules
    ]
  }

  module_names <- module_names[
    vapply(
      module_names,
      function(module_name) {
        !is.null(df_list[[module_name]]) &&
          nrow(df_list[[module_name]]) > 0L
      },
      logical(1)
    )
  ]

  if (length(module_names) == 0L) {
    return(make_empty_plot(
      "Gene-level heatmaps",
      "No modules were available"
    ))
  }

  comparisons <- contrast_key
  comparison_labels <- contrast_display_label
  n_comparisons <- length(comparisons)
  max_label_length <- max(nchar(comparison_labels))

  module_n <- purrr::map_dfr(
    module_names,
    function(module_name) {
      tibble::tibble(
        module = module_name,
        n_genes = dplyr::n_distinct(
          df_list[[module_name]]$Gene
        )
      )
    }
  )

  if (!is.null(force_block_height) &&
      is.finite(force_block_height)) {
    block_height <- force_block_height
  } else {
    block_height <- max(module_n$n_genes, na.rm = TRUE)
  }
  block_height <- max(block_height, 1)

  x_start <- 3.55
  x_gap <- if (n_comparisons >= 3L) 6.75 else 6.25
  if (max_label_length > 18L) x_gap <- x_gap + 0.55
  if (max_label_length > 24L) x_gap <- x_gap + 0.65

  layout_table <- tibble::tibble(
    module = module_names,
    layout_index = seq_along(module_names)
  ) |>
    dplyr::left_join(module_n, by = "module") |>
    dplyr::mutate(
      x_offset = x_start +
        (layout_index - 1) * (n_comparisons + x_gap),
      y_top = block_height,
      y_bottom = 0
    )

  tile_half_width <- tile_fraction / 2

  tile_data <- purrr::map_dfr(
    module_names,
    function(module_name) {
      df <- df_list[[module_name]]
      layout_one <- layout_table |>
        dplyr::filter(module == module_name)

      genes_in_order <- df |>
        dplyr::distinct(Gene, gene_order) |>
        dplyr::arrange(gene_order) |>
        dplyr::mutate(
          gene_order_plot = dplyr::row_number()
        )

      comparison_positions <- tibble::tibble(
        comparison = comparisons,
        x_local = seq_along(comparisons)
      )

      row_height <- block_height / max(layout_one$n_genes, 1)

      df |>
        dplyr::left_join(
          genes_in_order,
          by = c("Gene", "gene_order")
        ) |>
        dplyr::left_join(
          comparison_positions,
          by = "comparison"
        ) |>
        dplyr::mutate(
          x_plot = layout_one$x_offset + x_local - 1,
          y_plot = layout_one$y_top -
            (gene_order_plot - 0.5) * row_height,
          tile_height = row_height * tile_fraction,
          xmin = x_plot - tile_half_width,
          xmax = x_plot + tile_half_width,
          ymin = y_plot - tile_height / 2,
          ymax = y_plot + tile_height / 2,
          module = module_name
        )
    }
  )

  gene_label_data <- tile_data |>
    dplyr::distinct(module, Gene, y_plot) |>
    dplyr::left_join(layout_table, by = "module") |>
    dplyr::mutate(x_plot = x_offset - 0.62)

  module_title_data <- layout_table |>
    dplyr::mutate(
      x_plot = x_offset + (n_comparisons - 1) / 2,
      y_plot = y_top + 0.82,
      title_label = stringr::str_wrap(module, width = 18)
    )

  x_label_data <- purrr::map_dfr(
    module_names,
    function(module_name) {
      layout_one <- layout_table |>
        dplyr::filter(module == module_name)

      tibble::tibble(
        module = module_name,
        comparison_label = comparison_labels,
        x_plot = layout_one$x_offset +
          seq_along(comparisons) - 1,
        y_plot = layout_one$y_bottom - 0.52
      )
    }
  )

  x_min <- min(layout_table$x_offset) - 3.05
  x_max <- max(layout_table$x_offset) +
    n_comparisons - 0.05

  bottom_extra <- if (max_label_length > 24L) 6.6 else 5.8
  if (max_label_length > 35L) bottom_extra <- 7.4

  y_min <- min(x_label_data$y_plot) - bottom_extra
  y_max <- max(layout_table$y_top) + 1.55

  legend_position <- if (show_legend) "right" else "none"

  ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = tile_data,
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = logFC_colour
      ),
      colour = tile_border_colour,
      linewidth = tile_border_linewidth
    ) +
    ggplot2::geom_text(
      data = gene_label_data,
      ggplot2::aes(
        x = x_plot,
        y = y_plot,
        label = Gene
      ),
      hjust = 1,
      family = font_family,
      size = geom_text_size_from_pt(shared_label_pt),
      colour = "black"
    ) +
    ggplot2::geom_text(
      data = module_title_data,
      ggplot2::aes(
        x = x_plot,
        y = y_plot,
        label = title_label
      ),
      hjust = 0.5,
      vjust = 0,
      family = font_family,
      fontface = "bold",
      size = geom_text_size_from_pt(plot_title_size),
      lineheight = 0.90,
      colour = "black"
    ) +
    ggplot2::geom_text(
      data = x_label_data,
      ggplot2::aes(
        x = x_plot,
        y = y_plot,
        label = comparison_label
      ),
      hjust = 1,
      vjust = 1,
      angle = x_label_angle,
      family = font_family,
      size = geom_text_size_from_pt(shared_label_pt),
      colour = "black"
    ) +
    ggplot2::scale_fill_gradient2(
      low = negative_colour,
      mid = mid_colour,
      high = positive_colour,
      midpoint = 0,
      limits = c(-logfc_colour_limit, logfc_colour_limit),
      breaks = c(-logfc_colour_limit, 0, logfc_colour_limit),
      labels = c(
        paste0("<= -", logfc_colour_limit),
        "0",
        paste0(">= +", logfc_colour_limit)
      ),
      oob = scales::squish,
      na.value = "grey90",
      name = legend_title,
      guide = vertical_heatmap_guide()
    ) +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = c(x_min, x_max),
      ylim = c(y_min, y_max),
      clip = "off"
    ) +
    ggplot2::theme_void(base_family = font_family) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = font_family,
        colour = "black"
      ),
      legend.position = legend_position,
      legend.direction = "vertical",
      legend.title = ggplot2::element_text(
        family = font_family,
        size = legend_title_size,
        angle = 90,
        hjust = 0.5,
        vjust = 0.5,
        colour = "black"
      ),
      legend.text = ggplot2::element_text(
        family = font_family,
        size = legend_text_size,
        colour = "black"
      ),
      plot.margin = ggplot2::margin(
        2,
        3,
        2,
        3,
        unit = "mm"
      )
    )
}

four_module_heatmap_df_list <- lapply(
  four_module_order,
  function(module_name) {
    make_four_module_braf_heatmap_df(
      module_name,
      four_module_selected_gene_lists[[module_name]]
    )
  }
)
names(four_module_heatmap_df_list) <- four_module_order

four_module_heatmap_long <- dplyr::bind_rows(
  four_module_heatmap_df_list,
  .id = "module_list_name"
)

data.table::fwrite(
  four_module_heatmap_long,
  file.path(
    subdirs$tables,
    paste0(
      "FOUR_MODULE_PANEL_DISPLAYED_LOG2FC_VALUES_",
      "BRAF_V600E_VS_WT_SSLHP.csv"
    )
  )
)

p_four_module_braf <- build_four_module_braf_one_row(
  df_list = four_module_heatmap_df_list,
  show_legend = TRUE,
  legend_title = "log2FC"
)

four_module_main_stem <- paste0(
  "FINAL_ALL4_IFNA_CHEMOKINES_ANTIMICROBIAL_TNFA_",
  "LOG2FC_BRAF_V600E_VS_WT_SSLHP_NATURE_EXACT"
)

save_plot_set(
  p_four_module_braf,
  subdirs$four_module,
  four_module_main_stem,
  183,
  82
)

save_plot_set(
  p_four_module_braf,
  subdirs$combined,
  four_module_main_stem,
  183,
  82
)

common_four_module_height <- max(
  vapply(
    four_module_heatmap_df_list,
    function(x) dplyr::n_distinct(x$Gene),
    integer(1)
  ),
  1L
)

for (module_name in four_module_order) {
  p_one_four_module <- build_four_module_braf_one_row(
    df_list = four_module_heatmap_df_list[module_name],
    selected_modules = module_name,
    force_block_height = common_four_module_height,
    show_legend = TRUE,
    legend_title = "log2FC"
  )

  save_plot_set(
    p_one_four_module,
    subdirs$four_module,
    paste0(
      "FINAL_SINGLE_MODULE_BRAF_V600E_VS_WT_SSLHP_",
      safe_filename(module_name),
      "_NATURE_EXACT"
    ),
    60,
    105
  )
}

writeLines(
  c(
    "ORGANOID-MATCHED FOUR-MODULE BRAF HEATMAP PANEL",
    "================================================",
    "",
    "Modules, left to right:",
    "  Interferon alpha response",
    "  Chemokines",
    "  Antimicrobial peptides",
    "  TNFA/NFKB signaling",
    "",
    "Gene-set sources:",
    paste0("  IFNA: ", ifna_info$name),
    paste0(
      "  Chemokines: ",
      chemokine_activity_info_four_module$resolved_pathway,
      "; canonical chemokine ligands only"
    ),
    paste0(
      "  Antimicrobial peptides: ",
      antimicrobial_info_four_module$resolved_pathway
    ),
    paste0("  TNFA/NFKB: ", tnfa_info$name),
    "",
    paste0(
      "Selection: logCPM >= ",
      module_reference_logcpm_cutoff,
      " and abs(logFC) > ",
      module_reference_logfc_cutoff,
      "; top ",
      module_max_genes,
      " per module by decreasing abs(logFC)."
    ),
    "This is the organoid BRAF-specific two-direction selection rule.",
    "It is additional to, and does not replace, the existing positive-direction heatmaps.",
    "Empty modules are shown as grey placeholders; thresholds are not relaxed.",
    "",
    "Appearance:",
    paste0("  Font: ", font_family),
    paste0("  Labels: ", shared_label_pt, " pt"),
    paste0("  X-label angle: ", x_label_angle, " degrees"),
    paste0(
      "  Colour scale: ", negative_colour,
      " / white / ", positive_colour
    ),
    paste0(
      "  Saturation: -", logfc_colour_limit,
      " to +", logfc_colour_limit, " log2FC"
    ),
    "  White tile borders and right-side vertical log2FC legend.",
    "  Main export size: 183 x 82 mm."
  ),
  file.path(
    subdirs$four_module,
    "FOUR_MODULE_PANEL_METHODS_AND_SELECTION.txt"
  )
)

rank_method_summary <- rank_metadata |>
  dplyr::select(
    comparison,
    filter_mode,
    rank_method,
    n_rows_before_rank,
    n_ranked_genes,
    n_finite_lr,
    finite_lr_fraction
  ) |>
  dplyr::arrange(filter_mode)

data.table::fwrite(
  rank_method_summary,
  file.path(
    subdirs$qc,
    "FINAL_RANK_METHOD_SUMMARY_BRAF_V600E_VS_WT_SSLHP.csv"
  )
)

writeLines(
  c(
    "FINAL BRAF V600E VERSUS WT SSL_HP GSEA AND HEATMAP RUN",
    "======================================================",
    "",
    paste0("Run date: ", Sys.time()),
    paste0("Input directory: ", input_dir),
    paste0("Input DE file: ", de_path),
    paste0("Output directory: ", out_dir),
    "",
    "Contrast orientation:",
    "  This script assumes the DE table is BRAF V600E SSL_HP minus BRAF WT SSL_HP.",
    "  Positive logFC/NES therefore means higher or enriched in BRAF V600E SSL_HP.",
    "  Negative logFC/NES means higher or enriched in BRAF WT SSL_HP.",
    "  Confirm this against the original edgeR design and contrast before interpretation.",
    "",
    "Ranking:",
    "  Primary edgeR rank: sign(logFC) * sqrt(LR), when >=100 finite LR values remain.",
    "  Fallback: signed statistic, then sign(logFC) * -log10(PValue).",
    "  Duplicate symbols are collapsed by largest absolute rank, then P value, then absolute logFC.",
    "  Deterministic rank ties are broken by subtracting seq_along(ranks) * 1e-12.",
    "",
    "Hallmark GSEA:",
    paste0("  Gene-universe filter: logCPM >= ", hallmark_logcpm_cutoff, "."),
    paste0("  minSize = ", fgsea_min_size, "."),
    paste0("  maxSize = ", fgsea_max_size, "."),
    paste0("  eps = ", fgsea_eps, "."),
    paste0("  gseaParam = ", fgsea_gsea_param, "."),
    "  The complete Hallmark collection is tested for this one BRAF contrast.",
    "  BH FDR is calculated across the Hallmark collection.",
    "",
    "Final gastric-metaplasia GSEA:",
    paste0("  Signature source: ", metaplasia_signature_file),
    paste0("  Required and used genes: ", length(metaplasia_genes), "."),
    "  Gene-universe filter: none; ALL_VALID_GENES is used.",
    "  One custom pathway is tested, so padj equals pval; report it as nominal P, not Hallmark FDR.",
    "",
    "Gene-level heatmaps:",
    "  Full-signature heatmaps are written for IFNA, inflammatory response, TNFA/NFKB and gastric metaplasia.",
    "  All 24 gastric-metaplasia genes are displayed in supplied signature order, including grey tiles where absent.",
    "  Hallmark genes are ordered by decreasing absolute log2FC in the BRAF contrast.",
    "",
    "Publication-selected heatmaps:",
    paste0(
      "  Non-GM rule: logCPM >= ",
      module_reference_logcpm_cutoff,
      " and logFC > ",
      module_reference_logfc_cutoff,
      "; top ",
      module_max_genes,
      " by decreasing logFC."
    ),
    "  This selects genes increased in BRAF V600E, not large changes in either direction.",
    "  The threshold is never silently relaxed; full-signature heatmaps remain available if a module is empty.",
    "  Gastric metaplasia always uses all 24 genes without an expression or effect-size threshold.",
    "",
    "Additional organoid-matched four-module panel:",
    "  Modules are IFNA, chemokines, antimicrobial peptides and TNFA/NFKB in one horizontal row.",
    paste0(
      "  Selection uses logCPM >= ",
      module_reference_logcpm_cutoff,
      " and abs(logFC) > ",
      module_reference_logfc_cutoff,
      "; top ",
      module_max_genes,
      " by decreasing abs(logFC), matching the organoid BRAF panel."
    ),
    "  Canonical chemokines come from GOMF_CHEMOKINE_ACTIVITY; antimicrobial genes come from REACTOME_ANTIMICROBIAL_PEPTIDES.",
    "  This additional two-direction panel does not replace the positive-direction heatmaps above.",
    "",
    "Pathway-level heatmaps:",
    "  Stars on Hallmark heatmaps use Hallmark-family BH FDR.",
    "  Stars on the custom gastric-metaplasia NES heatmap use nominal P.",
    "",
    "Nature-style appearance:",
    paste0("  Font: ", font_family, "."),
    paste0("  Gene/comparison labels: ", shared_label_pt, " pt."),
    paste0("  Comparison-label angle: ", x_label_angle, " degrees."),
    paste0(
      "  Colours: ",
      negative_colour,
      " / white / ",
      positive_colour,
      "."
    ),
    paste0(
      "  log2FC colour saturation: -",
      logfc_colour_limit,
      " to +",
      logfc_colour_limit,
      "."
    ),
    "  White tile borders and right-side vertical legends are used.",
    "",
    "Interpretation:",
    "  This analysis describes transcriptomic differences associated with BRAF status within SSL_HP.",
    "  It does not by itself establish that BRAF V600E directly causes a pathway or gene-expression change.",
    "  There is only one contrast, so cross-comparison concordance statistics are not applicable."
  ),
  file.path(out_dir, "RUN_NOTES.txt")
)

package_names <- c(
  "R",
  "data.table",
  "dplyr",
  "tidyr",
  "stringr",
  "purrr",
  "ggplot2",
  "patchwork",
  "scales",
  "msigdbr",
  "fgsea"
)

package_versions <- data.frame(
  Package = package_names,
  Version = c(
    paste0(R.version$major, ".", R.version$minor),
    vapply(
      package_names[-1],
      function(pkg) {
        as.character(utils::packageVersion(pkg))
      },
      character(1)
    )
  ),
  stringsAsFactors = FALSE
)

package_versions$MSigDB_version <- NA_character_
package_versions$MSigDB_version[
  package_versions$Package == "msigdbr"
] <- msig_version_text

data.table::fwrite(
  package_versions,
  file.path(
    subdirs$qc,
    "PACKAGE_AND_MSIGDB_VERSIONS.csv"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    subdirs$qc,
    "SESSION_INFO.txt"
  )
)

message("\nDONE.")
message("Output directory: ", out_dir)
message("Inspect first:")
message("  02_Tables/BRAF_V600E_VS_WT_SSLHP_HALLMARK_FGSEA_ALL_RESULTS.csv")
message("  02_Tables/BRAF_V600E_VS_WT_SSLHP_FINAL_GASTRIC_METAPLASIA_FGSEA_ALL_VALID_GENES.csv")
message("  03_Hallmark_GSEA_curves/")
message("  04_Gastric_metaplasia_GSEA_curve/")
message("  05_Pathway_NES_heatmaps/")
message("  06_Full_signature_gene_log2FC_heatmaps/")
message("  07_Selected_gene_log2FC_heatmaps/")
message("  08_Combined_figures/")
message("  10_Gene_level_summary/")
message("  11_Four_module_heatmap_panel/")
