# Tissue GSEA and figures stratified by bacterial genus detection.
# Run from Tissue_Bulk_RNAseq after configuring config/config.R.
# The logCPM >= 2 filter applies to Hallmark analyses; the gastric signature uses unfiltered valid ranks.

source("config/config.R")

gc()

set.seed(1)
options(stringsAsFactors = FALSE)
options(contrasts = c("contr.treatment", "contr.poly"))

input_dir <- PATH_DE_TABLES

out_dir <- OUTPUT_BACTERIAL
ensure_new_output(out_dir)

install_missing_packages <- FALSE

require_exactly_24_metaplasia_genes <- TRUE

organism_display_order <- c("Serr", "Fuso", "ES")

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

full_hallmark_include_genes_absent_from_all_six_de_tables <- FALSE

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

metaplasia_display_title <-
  "Serrated Lesion Gastric Metaplasia Gene Signature"
metaplasia_display_title_wrapped <- paste(
  "Serrated Lesion Gastric",
  "Metaplasia Gene Signature",
  sep = "\n"
)

nes_x_label_angle <- 0
nes_x_label_hjust <- 0.5
nes_x_label_vjust <- 1
module_x_label_angle <- 90
module_x_label_hjust <- 1
module_x_label_vjust <- 0.5
six_comparison_heatmap_width_mm <- 183

curve_colours_all <- c(
  ES_SSL_HP = "#FF0000",
  Fuso_SSL_HP = "#CC0000",
  Serr_SSL_HP = "#800000",
  ES_H_ADJ_SSL_HP = "#E9C61D",
  Fuso_H_ADJ_SSL_HP = "#C99F00",
  Serr_H_ADJ_SSL_HP = "#8A6F00"
)

metaplasia_signature_candidates <- PATH_GASTRIC_SIGNATURE

de_file_table <- data.frame(
  comparison = c(
    "Serr_SSL_HP",
    "Fuso_SSL_HP",
    "ES_SSL_HP",
    "Serr_H_ADJ_SSL_HP",
    "Fuso_H_ADJ_SSL_HP",
    "ES_H_ADJ_SSL_HP"
  ),
  group = c(
    "SSL_HP", "SSL_HP", "SSL_HP",
    "H_ADJ_SSL_HP", "H_ADJ_SSL_HP", "H_ADJ_SSL_HP"
  ),
  group_order = c(1L, 1L, 1L, 2L, 2L, 2L),
  organism = c("Serr", "Fuso", "ES", "Serr", "Fuso", "ES"),
  file = c(
    "CGG014_Maier_Serr_SSL_HP_DE_results.txt",
    "CGG014_Maier_Fuso_SSL_HP_DE_results.txt",
    "CGG014_Maier_ES_SSL_HP_DE_results.txt",
    "CGG014_Maier_Serr_H_ADJ_SSL_HP_DE_results.txt",
    "CGG014_Maier_Fuso_H_ADJ_SSL_HP_DE_results.txt",
    "CGG014_Maier_ES_H_ADJ_SSL_HP_DE_results.txt"
  ),
  stringsAsFactors = FALSE
)

de_file_table$display_order <- match(
  de_file_table$organism,
  organism_display_order
)

if (anyNA(de_file_table$display_order)) {
  stop(
    "At least one organism in de_file_table is absent from ",
    "organism_display_order: ",
    paste(
      unique(de_file_table$organism[is.na(de_file_table$display_order)]),
      collapse = ", "
    )
  )
}

de_file_table <- de_file_table[
  order(de_file_table$group_order, de_file_table$display_order),
  ,
  drop = FALSE
]
row.names(de_file_table) <- NULL

de_file_table$path <- file.path(input_dir, de_file_table$file)

comparison_keys <- de_file_table$comparison

sslhp_idx <- de_file_table$group == "SSL_HP"
sslhp_comparisons <- de_file_table$comparison[sslhp_idx][
  order(de_file_table$display_order[sslhp_idx])
]

hadj_idx <- de_file_table$group == "H_ADJ_SSL_HP"
hadj_comparisons <- de_file_table$comparison[hadj_idx][
  order(de_file_table$display_order[hadj_idx])
]

comparisons_by_group <- list(
  SSL_HP = sslhp_comparisons,
  H_ADJ_SSL_HP = hadj_comparisons
)

reference_comparisons <- c(
  SSL_HP = "ES_SSL_HP",
  H_ADJ_SSL_HP = "ES_H_ADJ_SSL_HP"
)

group_display_labels <- c(
  SSL_HP = "SSL_HP",
  H_ADJ_SSL_HP = "H_ADJ_SSL_HP"
)

organism_long_labels <- c(
  Serr = "Serratia",
  Fuso = "Fusobacterium",
  ES = "Escherichia-Shigella"
)

comparison_short_labels <- stats::setNames(
  de_file_table$organism,
  de_file_table$comparison
)

comparison_display_labels <- stats::setNames(
  paste0(
    de_file_table$organism,
    " + vs -\n",
    group_display_labels[de_file_table$group]
  ),
  de_file_table$comparison
)

comparison_long_labels <- stats::setNames(
  paste0(
    unname(organism_long_labels[de_file_table$organism]),
    " positive versus negative in ",
    group_display_labels[de_file_table$group]
  ),
  de_file_table$comparison
)

curve_colours <- curve_colours_all[comparison_keys]

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

missing_de_files <- de_file_table |>
  dplyr::filter(!file.exists(path))

if (nrow(missing_de_files) > 0L) {
  stop(
    "These required DE files are missing:\n",
    paste(missing_de_files$file, collapse = "\n"),
    "\n\nDirectory checked:\n",
    input_dir
  )
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

subdirs <- list(
  qc = file.path(out_dir, "01_QC"),
  tables = file.path(out_dir, "02_Tables"),
  hallmark_curves = file.path(out_dir, "03_Hallmark_GSEA_overlays"),
  metaplasia_curves = file.path(out_dir, "04_Gastric_metaplasia_GSEA_curves"),
  hallmark_heatmaps = file.path(out_dir, "05_Hallmark_NES_heatmaps"),
  full_modules = file.path(out_dir, "06_Full_signature_gene_log2FC_heatmaps"),
  selected_modules = file.path(out_dir, "07_Selected_gene_log2FC_heatmaps"),
  combined = file.path(out_dir, "08_Combined_figures"),
  ranks = file.path(out_dir, "09_Ranked_gene_tables"),
  concordance = file.path(out_dir, "10_Gene_concordance"),
  four_module = file.path(out_dir, "11_Four_module_heatmap_panels")
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
  file.path(subdirs$qc, "SIX_DE_FILES_USED.csv")
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
  file.path(subdirs$tables, "SIX_COMPARISONS_HALLMARK_FGSEA_ALL_RESULTS.csv")
)

data.table::fwrite(
  metaplasia_fgsea_all,
  file.path(
    subdirs$tables,
    "SIX_COMPARISONS_FINAL_GASTRIC_METAPLASIA_FGSEA_ALL_VALID_GENES.csv"
  )
)

data.table::fwrite(
  combined_fgsea,
  file.path(subdirs$tables, "SIX_COMPARISONS_COMBINED_FGSEA_RESULTS.csv")
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
    "SIX_COMPARISONS_SELECTED_IFNA_INFLAMMATORY_TNFA_FGSEA.csv"
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

find_fgsea_row <- function(fgsea_df, comparison, pathway_name) {
  hit <- fgsea_df |>
    dplyr::filter(
      .data$comparison == .env$comparison,
      toupper(.data$pathway) == toupper(.env$pathway_name)
    ) |>
    dplyr::slice(1)

  if (nrow(hit) == 0L) return(NULL)
  hit
}

selected_hallmark_curve_info <- list(
  IFNA = list(
    name = ifna_info$name,
    genes = ifna_info$genes,
    label = "Interferon alpha response"
  ),
  INFLAMMATORY = list(
    name = inflammatory_info$name,
    genes = inflammatory_info$genes,
    label = "Inflammatory response"
  ),
  TNFA_NFKB = list(
    name = tnfa_info$name,
    genes = tnfa_info$genes,
    label = "TNFA/NFKB signaling"
  )
)

make_hallmark_overlay_plot <- function(
    comparisons,
    pathway_name,
    pathway_genes,
    title_text
) {
  curve_data <- purrr::map_dfr(
    comparisons,
    function(comp) {
      curve <- calc_enrichment_curve(
        hallmark_rank_objects[[comp]]$ranks,
        pathway_genes
      )

      if (nrow(curve) == 0L) {
        stop(
          "No enrichment curve could be calculated for ",
          comp,
          " / ",
          pathway_name,
          "."
        )
      }

      result <- find_fgsea_row(
        hallmark_fgsea_all,
        comp,
        pathway_name
      )

      if (is.null(result)) {
        stop(
          "Hallmark FGSEA result was not found for ",
          comp,
          " / ",
          pathway_name,
          "."
        )
      }

      curve |>
        dplyr::mutate(
          comparison = comp,
          legend_label = paste0(
            comparison_short_labels[[comp]],
            "\nNES=",
            sprintf("%.2f", result$NES[1]),
            ", FDR=",
            format_numeric_p(result$padj[1])
          )
        )
    }
  )

  hit_data <- curve_data |>
    dplyr::filter(hit)

  legend_table <- curve_data |>
    dplyr::distinct(comparison, legend_label)

  legend_labels <- stats::setNames(
    legend_table$legend_label,
    legend_table$comparison
  )

  ggplot2::ggplot(
    curve_data,
    ggplot2::aes(
      x = position_percent,
      y = running_es,
      colour = comparison
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = thin_line,
      colour = "grey45"
    ) +
    ggplot2::geom_line(linewidth = 0.78) +
    ggplot2::geom_rug(
      data = hit_data,
      ggplot2::aes(
        x = position_percent,
        colour = comparison
      ),
      sides = "b",
      inherit.aes = FALSE,
      linewidth = 0.18
    ) +
    ggplot2::scale_colour_manual(
      values = curve_colours[comparisons],
      breaks = comparisons,
      labels = legend_labels[comparisons],
      guide = ggplot2::guide_legend(
        ncol = 1,
        byrow = TRUE,
        keyheight = grid::unit(3.5, "mm"),
        keywidth = grid::unit(7, "mm")
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 25, 50, 75, 100),
      limits = c(0, 100),
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.08))
    ) +
    ggplot2::labs(
      title = title_text,
      x = "Ranked gene list position (%)",
      y = "Running enrichment score",
      colour = NULL
    ) +
    nature_theme() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(
        family = font_family,
        size = 4.6,
        lineheight = 0.85
      ),
      legend.title = ggplot2::element_blank(),
      legend.spacing.y = grid::unit(0.4, "mm"),
      legend.box.spacing = grid::unit(1, "mm"),
      plot.margin = ggplot2::margin(2, 4, 7, 2, unit = "mm")
    )
}

hallmark_overlay_plots <- list()
hallmark_overlay_rows <- list()

for (group_name in names(comparisons_by_group)) {
  group_comparisons <- comparisons_by_group[[group_name]]
  group_plots <- list()

  for (pathway_key in names(selected_hallmark_curve_info)) {
    pathway_info <- selected_hallmark_curve_info[[pathway_key]]

    p <- make_hallmark_overlay_plot(
      comparisons = group_comparisons,
      pathway_name = pathway_info$name,
      pathway_genes = pathway_info$genes,
      title_text = paste0(
        group_display_labels[[group_name]],
        " | ",
        pathway_info$label
      )
    )

    storage_key <- paste(group_name, pathway_key, sep = "__")
    hallmark_overlay_plots[[storage_key]] <- p
    group_plots[[pathway_key]] <- p

    save_plot_set(
      p,
      subdirs$hallmark_curves,
      paste0(
        "HALLMARK_OVERLAY_",
        group_name,
        "_",
        pathway_key,
        "_SERR_FUSO_ES"
      ),
      90,
      78
    )
  }

  p_row <- patchwork::wrap_plots(
    group_plots[c("IFNA", "INFLAMMATORY", "TNFA_NFKB")],
    ncol = 3,
    guides = "keep"
  ) +
    patchwork::plot_annotation(
      title = paste0(
        group_display_labels[[group_name]],
        ": Hallmark GSEA overlays"
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          family = font_family,
          face = "bold",
          size = 8,
          hjust = 0
        )
      )
    )

  hallmark_overlay_rows[[group_name]] <- p_row

  save_plot_set(
    p_row,
    subdirs$combined,
    paste0(
      "HALLMARK_OVERLAY_ROW_",
      group_name,
      "_IFNA_INFLAMMATORY_TNFA"
    ),
    183,
    88
  )
}

p_hallmark_overlays_final <- patchwork::wrap_plots(
  hallmark_overlay_plots[["SSL_HP__IFNA"]],
  hallmark_overlay_plots[["SSL_HP__INFLAMMATORY"]],
  hallmark_overlay_plots[["SSL_HP__TNFA_NFKB"]],
  hallmark_overlay_plots[["H_ADJ_SSL_HP__IFNA"]],
  hallmark_overlay_plots[["H_ADJ_SSL_HP__INFLAMMATORY"]],
  hallmark_overlay_plots[["H_ADJ_SSL_HP__TNFA_NFKB"]],
  ncol = 3,
  guides = "keep"
) +
  patchwork::plot_annotation(
    tag_levels = "a",
    theme = ggplot2::theme(
      text = ggplot2::element_text(
        family = font_family,
        colour = "black"
      ),
      plot.tag = ggplot2::element_text(
        family = font_family,
        face = "bold",
        size = panel_tag_size
      )
    )
  ) &
  ggplot2::theme(
    text = ggplot2::element_text(family = font_family)
  )

save_plot_set(
  p_hallmark_overlays_final,
  subdirs$combined,
  "FINAL_2X3_HALLMARK_OVERLAYS_SSLHP_AND_HADJ_IFNA_INFLAMMATORY_TNFA_183MM",
  183,
  155
)

save_plot_set(
  p_hallmark_overlays_final,
  subdirs$combined,
  "FINAL_2X3_HALLMARK_OVERLAYS_SSLHP_AND_HADJ_IFNA_INFLAMMATORY_TNFA_220MM",
  220,
  165
)

make_metaplasia_curve_plot <- function(comp) {
  rank_object <- metaplasia_rank_objects[[comp]]
  result <- metaplasia_fgsea_all |>
    dplyr::filter(comparison == comp) |>
    dplyr::slice(1)

  curve <- calc_enrichment_curve(
    rank_object$ranks,
    metaplasia_genes
  )

  if (nrow(curve) == 0L || nrow(result) == 0L) {
    return(make_empty_plot(
      comparison_short_labels[[comp]],
      "No overlapping gastric-metaplasia genes"
    ))
  }

  annotation <- paste0(
    "NES = ",
    sprintf("%.2f", result$NES[1]),
    "\nNominal p = ",
    format_numeric_p(result$pval[1]),
    "\nGenes = ",
    result$size[1]
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
      colour = curve_colours[[comp]],
      linewidth = 0.75
    ) +
    ggplot2::geom_rug(
      data = curve |> dplyr::filter(hit),
      ggplot2::aes(x = position_percent),
      sides = "b",
      inherit.aes = FALSE,
      colour = curve_colours[[comp]],
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
      title = comparison_long_labels[[comp]],
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

metaplasia_curve_plots <- lapply(
  comparison_keys,
  make_metaplasia_curve_plot
)
names(metaplasia_curve_plots) <- comparison_keys

for (comp in comparison_keys) {
  save_plot_set(
    metaplasia_curve_plots[[comp]],
    subdirs$metaplasia_curves,
    paste0(comp, "_FINAL_GASTRIC_METAPLASIA_GSEA_CURVE"),
    90,
    72
  )
}

metaplasia_curve_rows <- list()

for (group_name in names(comparisons_by_group)) {
  group_comparisons <- comparisons_by_group[[group_name]]

  p_row <- patchwork::wrap_plots(
    metaplasia_curve_plots[group_comparisons],
    ncol = 3,
    guides = "keep"
  ) +
    patchwork::plot_annotation(
      title = paste0(
        metaplasia_display_title_wrapped,
        " | ",
        group_display_labels[[group_name]]
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          family = font_family,
          face = "bold",
          size = 8,
          hjust = 0
        )
      )
    )

  metaplasia_curve_rows[[group_name]] <- p_row

  save_plot_set(
    p_row,
    subdirs$combined,
    paste0(
      "FINAL_GASTRIC_METAPLASIA_GSEA_ROW_",
      group_name,
      "_SERR_FUSO_ES"
    ),
    183,
    84
  )
}

p_metaplasia_curves_2x3 <- patchwork::wrap_plots(
  metaplasia_curve_plots[comparison_keys],
  ncol = 3,
  guides = "keep"
) +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = paste0(
      metaplasia_display_title_wrapped,
      " | SSL_HP and H_ADJ_SSL_HP"
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
  p_metaplasia_curves_2x3,
  subdirs$combined,
  "FINAL_GASTRIC_METAPLASIA_GSEA_2X3_ALL_SIX_COMPARISONS",
  183,
  158
)

make_metaplasia_overlay_plot <- function(
    comparisons,
    title_text,
    legend_ncol = 1L
) {
  curve_data <- purrr::map_dfr(
    comparisons,
    function(comp) {
      result <- metaplasia_fgsea_all |>
        dplyr::filter(comparison == comp) |>
        dplyr::slice(1)

      curve <- calc_enrichment_curve(
        metaplasia_rank_objects[[comp]]$ranks,
        metaplasia_genes
      )

      if (nrow(curve) == 0L || nrow(result) == 0L) {
        stop(
          "No gastric-metaplasia curve/result was available for ",
          comp,
          "."
        )
      }

      curve |>
        dplyr::mutate(
          comparison = comp,
          legend_label = paste0(
            comparison_short_labels[[comp]],
            " | ",
            group_display_labels[[
              de_file_table$group[
                match(comp, de_file_table$comparison)
              ]
            ]],
            "\nNES=",
            sprintf("%.2f", result$NES[1]),
            ", p=",
            format_numeric_p(result$pval[1])
          )
        )
    }
  )

  hit_data <- curve_data |>
    dplyr::filter(hit)

  legend_table <- curve_data |>
    dplyr::distinct(comparison, legend_label)

  legend_labels <- stats::setNames(
    legend_table$legend_label,
    legend_table$comparison
  )

  ggplot2::ggplot(
    curve_data,
    ggplot2::aes(
      x = position_percent,
      y = running_es,
      colour = comparison
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = thin_line,
      colour = "grey45"
    ) +
    ggplot2::geom_line(linewidth = 0.78) +
    ggplot2::geom_rug(
      data = hit_data,
      ggplot2::aes(
        x = position_percent,
        colour = comparison
      ),
      sides = "b",
      inherit.aes = FALSE,
      linewidth = 0.18
    ) +
    ggplot2::scale_colour_manual(
      values = curve_colours[comparisons],
      breaks = comparisons,
      labels = legend_labels[comparisons],
      guide = ggplot2::guide_legend(
        ncol = legend_ncol,
        byrow = TRUE,
        keywidth = grid::unit(7, "mm")
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 25, 50, 75, 100),
      limits = c(0, 100),
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.08))
    ) +
    ggplot2::labs(
      title = title_text,
      x = "Ranked gene list position (%)",
      y = "Running enrichment score",
      colour = NULL
    ) +
    nature_theme() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(
        family = font_family,
        size = 4.7,
        lineheight = 0.88
      )
    )
}

p_metaplasia_overlay_sslhp <- make_metaplasia_overlay_plot(
  comparisons = sslhp_comparisons,
  title_text = paste0(metaplasia_display_title_wrapped, " | SSL_HP"),
  legend_ncol = 1L
)

p_metaplasia_overlay_hadj <- make_metaplasia_overlay_plot(
  comparisons = hadj_comparisons,
  title_text = paste0(metaplasia_display_title_wrapped, " | H_ADJ_SSL_HP"),
  legend_ncol = 1L
)

p_metaplasia_overlay_all6 <- make_metaplasia_overlay_plot(
  comparisons = comparison_keys,
  title_text = paste0(metaplasia_display_title_wrapped, " | all six comparisons"),
  legend_ncol = 2L
)

save_plot_set(
  p_metaplasia_overlay_sslhp,
  subdirs$metaplasia_curves,
  "FINAL_GASTRIC_METAPLASIA_GSEA_OVERLAY_SSL_HP_SERR_FUSO_ES",
  110,
  84
)

save_plot_set(
  p_metaplasia_overlay_hadj,
  subdirs$metaplasia_curves,
  "FINAL_GASTRIC_METAPLASIA_GSEA_OVERLAY_H_ADJ_SSL_HP_SERR_FUSO_ES",
  110,
  84
)

save_plot_set(
  p_metaplasia_overlay_all6,
  subdirs$metaplasia_curves,
  "FINAL_GASTRIC_METAPLASIA_GSEA_OVERLAY_ALL_SIX_COMPARISONS",
  140,
  100
)

make_metaplasia_nes_heatmap <- function() {
  plot_df <- metaplasia_fgsea_all |>
    dplyr::mutate(
      comparison_display = factor(
        unname(comparison_display_labels[comparison]),
        levels = unname(comparison_display_labels[comparison_keys])
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
        angle = nes_x_label_angle,
        hjust = nes_x_label_hjust,
        vjust = nes_x_label_vjust,
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

p_metaplasia_nes_heatmap <- make_metaplasia_nes_heatmap()

save_plot_set(
  p_metaplasia_nes_heatmap,
  subdirs$hallmark_heatmaps,
  "FINAL_GASTRIC_METAPLASIA_NES_HEATMAP_ALL_SIX_COMPARISONS_NOMINAL_P",
  six_comparison_heatmap_width_mm,
  70
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
    comparisons = comparison_keys,
    pathways = NULL,
    top_n_per_comparison = NULL,
    title_text = NULL,
    legend_title = "NES"
) {
  comparisons <- comparisons[comparisons %in% comparison_keys]

  if (length(comparisons) == 0L) {
    return(make_empty_plot(
      "Hallmark GSEA",
      "No valid comparisons were supplied"
    ))
  }

  dat <- hallmark_fgsea_all |>
    dplyr::filter(
      comparison %in% comparisons,
      is.finite(NES)
    ) |>
    dplyr::group_by(comparison, pathway) |>
    dplyr::arrange(
      padj,
      pval,
      dplyr::desc(abs(NES)),
      .by_group = TRUE
    ) |>
    dplyr::slice(1) |>
    dplyr::ungroup()

  if (nrow(dat) == 0L) {
    return(make_empty_plot("Hallmark GSEA", "No Hallmark results"))
  }

  if (!is.null(pathways)) {
    selected_pathways <- unique(pathways[pathways %in% dat$pathway])
  } else if (!is.null(top_n_per_comparison)) {
    selected_pathways <- dat |>
      dplyr::group_by(comparison) |>
      dplyr::arrange(
        padj,
        dplyr::desc(abs(NES)),
        .by_group = TRUE
      ) |>
      dplyr::slice_head(n = top_n_per_comparison) |>
      dplyr::ungroup() |>
      dplyr::distinct(pathway) |>
      dplyr::pull(pathway)
  } else {
    selected_pathways <- sort(unique(dat$pathway))
  }

  if (length(selected_pathways) == 0L) {
    return(make_empty_plot(
      "Hallmark GSEA",
      "No requested pathways were found"
    ))
  }

  pathway_order <- dat |>
    dplyr::filter(pathway %in% selected_pathways) |>
    dplyr::group_by(pathway) |>
    dplyr::summarise(
      max_abs_NES = max(abs(NES), na.rm = TRUE),
      best_FDR = safe_min(padj),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::desc(max_abs_NES),
      best_FDR,
      pathway
    ) |>
    dplyr::pull(pathway)

  comparison_display_levels <- unname(
    comparison_display_labels[comparisons]
  )

  plot_df <- tidyr::expand_grid(
    comparison = comparisons,
    pathway = pathway_order
  ) |>
    dplyr::left_join(
      dat |>
        dplyr::select(
          comparison,
          pathway,
          NES,
          pval,
          padj
        ),
      by = c("comparison", "pathway")
    ) |>
    dplyr::mutate(
      comparison_display = factor(
        unname(comparison_display_labels[comparison]),
        levels = comparison_display_levels
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
      name = legend_title,
      guide = vertical_heatmap_guide()
    ) +
    ggplot2::coord_fixed(ratio = 1, clip = "off") +
    ggplot2::labs(title = title_text, x = NULL, y = NULL) +
    nature_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        family = font_family,
        size = shared_label_pt,
        angle = nes_x_label_angle,
        hjust = nes_x_label_hjust,
        vjust = nes_x_label_vjust,
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

p_hallmark_all6_all <- make_hallmark_nes_heatmap(
  comparisons = comparison_keys,
  pathways = sort(unique(hallmark_fgsea_all$pathway)),
  title_text = NULL
)

p_hallmark_all6_top <- make_hallmark_nes_heatmap(
  comparisons = comparison_keys,
  top_n_per_comparison = 20L,
  title_text = NULL
)

p_hallmark_all6_key <- make_hallmark_nes_heatmap(
  comparisons = comparison_keys,
  pathways = key_hallmark_pathways,
  title_text = NULL
)

p_hallmark_all6_selected_three <- make_hallmark_nes_heatmap(
  comparisons = comparison_keys,
  pathways = selected_three_hallmark_pathways,
  title_text = NULL
)

p_hallmark_sslhp_selected_three <- make_hallmark_nes_heatmap(
  comparisons = sslhp_comparisons,
  pathways = selected_three_hallmark_pathways,
  title_text = "SSL_HP"
)

p_hallmark_hadj_selected_three <- make_hallmark_nes_heatmap(
  comparisons = hadj_comparisons,
  pathways = selected_three_hallmark_pathways,
  title_text = "H_ADJ_SSL_HP"
)

save_plot_set(
  p_hallmark_all6_all,
  subdirs$hallmark_heatmaps,
  "FINAL_ALL_HALLMARK_NES_HEATMAP_ALL_SIX_COMPARISONS_NATURE_EXACT",
  six_comparison_heatmap_width_mm,
  220
)

save_plot_set(
  p_hallmark_all6_top,
  subdirs$hallmark_heatmaps,
  "FINAL_TOP20_PER_COMPARISON_HALLMARK_NES_HEATMAP_ALL_SIX_NATURE_EXACT",
  six_comparison_heatmap_width_mm,
  165
)

save_plot_set(
  p_hallmark_all6_key,
  subdirs$hallmark_heatmaps,
  "FINAL_KEY_HALLMARK_NES_HEATMAP_ALL_SIX_NATURE_EXACT",
  six_comparison_heatmap_width_mm,
  110
)

save_plot_set(
  p_hallmark_all6_selected_three,
  subdirs$hallmark_heatmaps,
  "FINAL_IFNA_INFLAMMATORY_TNFA_NES_HEATMAP_ALL_SIX_COMPARISONS",
  six_comparison_heatmap_width_mm,
  82
)

save_plot_set(
  p_hallmark_sslhp_selected_three,
  subdirs$hallmark_heatmaps,
  "FINAL_IFNA_INFLAMMATORY_TNFA_NES_HEATMAP_SSLHP_ONLY",
  110,
  78
)

save_plot_set(
  p_hallmark_hadj_selected_three,
  subdirs$hallmark_heatmaps,
  "FINAL_IFNA_INFLAMMATORY_TNFA_NES_HEATMAP_HADJ_ONLY",
  110,
  78
)

p_selected_hallmark_by_group <- p_hallmark_sslhp_selected_three /
  p_hallmark_hadj_selected_three +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = "Selected Hallmark pathways by tissue compartment",
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
  p_selected_hallmark_by_group,
  subdirs$combined,
  "FINAL_SELECTED_HALLMARK_NES_HEATMAPS_SSLHP_AND_HADJ_TWO_ROWS",
  135,
  155
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

heatmap_de_list <- lapply(
  comparison_keys,
  function(comp) {
    collapse_de_for_heatmap(
      standardised_de_list[[comp]],
      comp
    )
  }
)
names(heatmap_de_list) <- comparison_keys

heatmap_de_all <- dplyr::bind_rows(heatmap_de_list) |>
  dplyr::left_join(
    de_file_table |>
      dplyr::select(
        comparison,
        group,
        group_order,
        organism,
        display_order
      ),
    by = "comparison"
  ) |>
  dplyr::arrange(group_order, display_order, Gene)

data.table::fwrite(
  heatmap_de_all,
  file.path(
    subdirs$tables,
    "SIX_COMPARISONS_COLLAPSED_DE_FOR_GENE_HEATMAPS.csv"
  )
)

module_gene_sets <- list(
  "Interferon alpha response" = ifna_info$genes,
  "Inflammatory response" = inflammatory_info$genes,
  "TNFA/NFKB signaling" = tnfa_info$genes,
  "Gastric metaplasia" = metaplasia_genes
)

module_order <- names(module_gene_sets)

module_display_title <- function(module_name, wrap_width = 18L) {
  if (identical(module_name, "Gastric metaplasia")) {
    return(stringr::str_wrap(metaplasia_display_title, width = 26L))
  }
  stringr::str_wrap(module_name, width = wrap_width)
}

build_module_concordance <- function(module_name, genes) {
  genes <- unique(clean_gene_symbols(genes))
  genes <- genes[!is.na(genes) & genes != ""]

  gene_index <- tibble::tibble(
    Gene = genes,
    signature_order = seq_along(genes)
  )

  long <- tidyr::expand_grid(
    comparison = comparison_keys,
    Gene = genes
  ) |>
    dplyr::left_join(gene_index, by = "Gene") |>
    dplyr::left_join(
      heatmap_de_all |>
        dplyr::select(
          comparison,
          Gene,
          logFC,
          logCPM,
          PValue,
          FDR,
          group,
          group_order,
          organism,
          display_order
        ),
      by = c("comparison", "Gene")
    ) |>
    dplyr::mutate(
      module = module_name,
      present = is.finite(logFC)
    )

  if (module_name == "Gastric metaplasia") {
    long <- long |>
      dplyr::mutate(expression_eligible = present)
  } else {
    long <- long |>
      dplyr::mutate(
        expression_eligible =
          present &
          is.finite(logCPM) &
          logCPM >= hallmark_logcpm_cutoff
      )
  }

  long <- long |>
    dplyr::mutate(
      up = present & expression_eligible & logFC > 0,
      down = present & expression_eligible & logFC < 0,
      significant_up = up & is.finite(FDR) & FDR < 0.05,
      significant_down = down & is.finite(FDR) & FDR < 0.05
    ) |>
    dplyr::arrange(signature_order, group_order, display_order)

  summary <- long |>
    dplyr::group_by(module, Gene, signature_order) |>
    dplyr::summarise(
      n_present = sum(present, na.rm = TRUE),
      n_expression_eligible = sum(expression_eligible, na.rm = TRUE),
      n_up = sum(up, na.rm = TRUE),
      n_down = sum(down, na.rm = TRUE),
      n_significant_up = sum(significant_up, na.rm = TRUE),
      n_significant_down = sum(significant_down, na.rm = TRUE),
      n_up_SSL_HP = sum(up & group == "SSL_HP", na.rm = TRUE),
      n_up_H_ADJ_SSL_HP = sum(
        up & group == "H_ADJ_SSL_HP",
        na.rm = TRUE
      ),
      n_significant_up_SSL_HP = sum(
        significant_up & group == "SSL_HP",
        na.rm = TRUE
      ),
      n_significant_up_H_ADJ_SSL_HP = sum(
        significant_up & group == "H_ADJ_SSL_HP",
        na.rm = TRUE
      ),
      mean_logFC_all_present = if (
        any(present)
      ) {
        mean(logFC[present], na.rm = TRUE)
      } else {
        NA_real_
      },
      mean_logFC_expression_eligible = if (
        any(expression_eligible & present)
      ) {
        mean(logFC[expression_eligible & present], na.rm = TRUE)
      } else {
        NA_real_
      },
      maximum_logFC = safe_max(logFC[present]),
      minimum_logFC = safe_min(logFC[present]),
      maximum_absolute_logFC = safe_max(abs(logFC[present])),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      all_six_present = n_present == length(comparison_keys),
      all_six_up = n_up == length(comparison_keys),
      all_six_down = n_down == length(comparison_keys),
      all_three_SSL_HP_up = n_up_SSL_HP == length(sslhp_comparisons),
      all_three_H_ADJ_SSL_HP_up =
        n_up_H_ADJ_SSL_HP == length(hadj_comparisons),
      up_in_both_tissue_groups =
        n_up_SSL_HP > 0 & n_up_H_ADJ_SSL_HP > 0,
      direction_class = dplyr::case_when(
        all_six_up ~ "Up in all six comparisons",
        all_six_down ~ "Down in all six comparisons",
        n_up > 0 & n_down == 0 ~ "Up where direction is non-zero",
        n_down > 0 & n_up == 0 ~ "Down where direction is non-zero",
        n_up > 0 & n_down > 0 ~ "Mixed direction",
        TRUE ~ "No non-zero eligible direction"
      )
    )

  list(long = long, summary = summary)
}

module_concordance_objects <- purrr::imap(
  module_gene_sets,
  function(genes, module_name) {
    build_module_concordance(module_name, genes)
  }
)

module_concordance_long <- dplyr::bind_rows(
  lapply(module_concordance_objects, `[[`, "long")
)

module_concordance_summary <- dplyr::bind_rows(
  lapply(module_concordance_objects, `[[`, "summary")
) |>
  dplyr::arrange(
    factor(module, levels = module_order),
    dplyr::desc(n_significant_up),
    dplyr::desc(n_up),
    dplyr::desc(mean_logFC_expression_eligible),
    Gene
  )

module_concordance_wide <- module_concordance_long |>
  dplyr::select(
    module,
    Gene,
    signature_order,
    comparison,
    logFC,
    logCPM,
    PValue,
    FDR,
    up,
    significant_up
  ) |>
  tidyr::pivot_wider(
    names_from = comparison,
    values_from = c(
      logFC,
      logCPM,
      PValue,
      FDR,
      up,
      significant_up
    ),
    names_glue = "{.value}__{comparison}"
  ) |>
  dplyr::left_join(
    module_concordance_summary,
    by = c("module", "Gene", "signature_order")
  )

data.table::fwrite(
  module_concordance_long,
  file.path(
    subdirs$concordance,
    "ALL_MODULE_GENES_ALL_SIX_COMPARISONS_LONG.csv"
  )
)

data.table::fwrite(
  module_concordance_summary,
  file.path(
    subdirs$concordance,
    "ALL_MODULE_GENES_CONCORDANCE_SUMMARY.csv"
  )
)

data.table::fwrite(
  module_concordance_wide,
  file.path(
    subdirs$concordance,
    "ALL_MODULE_GENES_CONCORDANCE_WIDE.csv"
  )
)

compute_pairwise_module_concordance <- function(module_name, long_df) {
  comparison_pairs <- utils::combn(
    comparison_keys,
    2,
    simplify = FALSE
  )

  purrr::map_dfr(
    comparison_pairs,
    function(pair) {
      comp_a <- pair[1]
      comp_b <- pair[2]

      a <- long_df |>
        dplyr::filter(
          module == module_name,
          comparison == comp_a
        ) |>
        dplyr::select(
          Gene,
          logFC_A = logFC,
          eligible_A = expression_eligible,
          up_A = up,
          significant_up_A = significant_up
        )

      b <- long_df |>
        dplyr::filter(
          module == module_name,
          comparison == comp_b
        ) |>
        dplyr::select(
          Gene,
          logFC_B = logFC,
          eligible_B = expression_eligible,
          up_B = up,
          significant_up_B = significant_up
        )

      joined <- dplyr::inner_join(a, b, by = "Gene") |>
        dplyr::filter(
          is.finite(logFC_A),
          is.finite(logFC_B),
          eligible_A,
          eligible_B
        )

      n_genes <- nrow(joined)

      rho <- if (n_genes >= 3L) {
        suppressWarnings(
          stats::cor(
            joined$logFC_A,
            joined$logFC_B,
            method = "spearman",
            use = "complete.obs"
          )
        )
      } else {
        NA_real_
      }

      up_union <- sum(joined$up_A | joined$up_B, na.rm = TRUE)
      up_intersection <- sum(joined$up_A & joined$up_B, na.rm = TRUE)
      significant_up_union <- sum(
        joined$significant_up_A | joined$significant_up_B,
        na.rm = TRUE
      )
      significant_up_intersection <- sum(
        joined$significant_up_A & joined$significant_up_B,
        na.rm = TRUE
      )

      data.frame(
        module = module_name,
        comparison_A = comp_a,
        comparison_B = comp_b,
        n_genes_compared = n_genes,
        Spearman_rho_logFC = rho,
        same_direction_fraction = ifelse(
          n_genes > 0L,
          mean(
            sign(joined$logFC_A) == sign(joined$logFC_B),
            na.rm = TRUE
          ),
          NA_real_
        ),
        n_up_in_both = up_intersection,
        n_up_in_either = up_union,
        Jaccard_up = ifelse(
          up_union > 0L,
          up_intersection / up_union,
          NA_real_
        ),
        n_significant_up_in_both = significant_up_intersection,
        n_significant_up_in_either = significant_up_union,
        Jaccard_significant_up = ifelse(
          significant_up_union > 0L,
          significant_up_intersection / significant_up_union,
          NA_real_
        ),
        stringsAsFactors = FALSE
      )
    }
  )
}

pairwise_module_concordance <- purrr::map_dfr(
  module_order,
  function(module_name) {
    compute_pairwise_module_concordance(
      module_name,
      module_concordance_long
    )
  }
)

data.table::fwrite(
  pairwise_module_concordance,
  file.path(
    subdirs$concordance,
    "PAIRWISE_MODULE_LOGFC_AND_UPREGULATED_GENE_CONCORDANCE.csv"
  )
)

make_module_heatmap_df <- function(
    module_name,
    genes,
    comparisons
) {
  genes <- unique(clean_gene_symbols(genes))
  genes <- genes[!is.na(genes) & genes != ""]
  comparisons <- comparisons[comparisons %in% comparison_keys]

  gene_index <- tibble::tibble(
    Gene = genes,
    gene_order = seq_along(genes)
  )

  tidyr::expand_grid(
    comparison = comparisons,
    gene_index
  ) |>
    dplyr::left_join(
      heatmap_de_all |>
        dplyr::select(
          comparison,
          Gene,
          logFC,
          logCPM,
          PValue,
          FDR
        ),
      by = c("comparison", "Gene")
    ) |>
    dplyr::mutate(
      module = module_name,
      comparison_display = factor(
        unname(comparison_display_labels[comparison]),
        levels = unname(comparison_display_labels[comparisons])
      ),
      logFC_colour = pmax(
        pmin(logFC, logfc_colour_limit),
        -logfc_colour_limit
      )
    )
}

build_one_row_module_heatmap <- function(
    df_list,
    comparisons,
    module_order_used = module_order,
    selected_modules = NULL,
    show_legend = TRUE,
    legend_title = "log2FC",
    force_block_height = NULL
) {
  comparisons <- comparisons[comparisons %in% comparison_keys]

  module_names <- module_order_used[
    module_order_used %in% names(df_list)
  ]

  if (!is.null(selected_modules)) {
    module_names <- module_names[module_names %in% selected_modules]
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
      "No genes were available for this heatmap"
    ))
  }

  n_comparisons <- length(comparisons)
  comparison_labels <- unname(
    comparison_display_labels[comparisons]
  )
  comparison_labels_dense <- stringr::str_replace_all(
    comparison_labels,
    "\n",
    " "
  )
  max_label_length <- max(nchar(comparison_labels_dense))

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
  x_gap <- if (n_comparisons >= 6L) 7.75 else 6.75
  if (max_label_length > 24L) x_gap <- x_gap + 0.65
  if (max_label_length > 35L) x_gap <- x_gap + 0.80

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
      title_label = vapply(
        module,
        module_display_title,
        FUN.VALUE = character(1),
        wrap_width = 18L
      ),
      title_lines = stringr::str_count(title_label, "\n") + 1L
    )

  x_label_data <- purrr::map_dfr(
    module_names,
    function(module_name) {
      layout_one <- layout_table |>
        dplyr::filter(module == module_name)

      tibble::tibble(
        module = module_name,
        comparison_label = comparison_labels_dense,
        x_plot = layout_one$x_offset +
          seq_along(comparisons) - 1,
        y_plot = layout_one$y_bottom - 0.52
      )
    }
  )

  x_min <- min(layout_table$x_offset) - 3.05
  x_max <- max(layout_table$x_offset) + n_comparisons - 0.05

  bottom_extra <- if (max_label_length > 24L) 6.6 else 5.8
  if (max_label_length > 35L) bottom_extra <- 7.4

  y_min <- min(x_label_data$y_plot) - bottom_extra
  max_title_lines <- max(module_title_data$title_lines, na.rm = TRUE)
  y_max <- max(layout_table$y_top) + 1.55 +
    0.70 * max(max_title_lines - 1L, 0L)

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
      hjust = module_x_label_hjust,
      vjust = module_x_label_vjust,
      angle = module_x_label_angle,
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

make_full_module_gene_qc <- function(genes, module_name) {
  genes <- unique(clean_gene_symbols(genes))
  genes <- genes[!is.na(genes) & genes != ""]

  signature_table <- tibble::tibble(
    module = module_name,
    Gene = genes,
    signature_order = seq_along(genes)
  )

  gene_stats <- heatmap_de_all |>
    dplyr::filter(
      comparison %in% comparison_keys,
      Gene %in% genes
    ) |>
    dplyr::group_by(Gene) |>
    dplyr::summarise(
      comparisons_present = sum(is.finite(logFC)),
      max_abs_logFC = ifelse(
        any(is.finite(logFC)),
        max(abs(logFC), na.rm = TRUE),
        NA_real_
      ),
      mean_abs_logFC = ifelse(
        any(is.finite(logFC)),
        mean(abs(logFC), na.rm = TRUE),
        NA_real_
      ),
      .groups = "drop"
    )

  signature_table |>
    dplyr::left_join(gene_stats, by = "Gene") |>
    dplyr::mutate(
      comparisons_present = dplyr::coalesce(
        comparisons_present,
        0L
      ),
      present_in_any_of_six_DE_tables = comparisons_present > 0L
    )
}

full_module_gene_qc <- purrr::imap_dfr(
  module_gene_sets,
  make_full_module_gene_qc
)

full_module_gene_lists <- list()

for (module_name in module_order) {
  module_qc <- full_module_gene_qc |>
    dplyr::filter(module == module_name)

  if (
    module_name != "Gastric metaplasia" &&
    !full_hallmark_include_genes_absent_from_all_six_de_tables
  ) {
    module_qc <- module_qc |>
      dplyr::filter(present_in_any_of_six_DE_tables)
  }

  if (module_name == "Gastric metaplasia") {
    module_qc <- module_qc |>
      dplyr::arrange(signature_order)
  } else {
    module_qc <- module_qc |>
      dplyr::arrange(
        dplyr::desc(max_abs_logFC),
        signature_order,
        Gene
      )
  }

  full_module_gene_lists[[module_name]] <- module_qc$Gene
}

data.table::fwrite(
  full_module_gene_qc,
  file.path(
    subdirs$qc,
    "FULL_SIGNATURE_GENE_COVERAGE_AND_ORDER_QC_ALL_SIX.csv"
  )
)

full_module_heatmap_df_list <- list()

for (module_name in module_order) {
  genes <- full_module_gene_lists[[module_name]]

  if (length(genes) == 0L) {
    warning(
      "No genes were available for the full-signature heatmap: ",
      module_name,
      call. = FALSE
    )
    next
  }

  full_module_heatmap_df_list[[module_name]] <-
    make_module_heatmap_df(
      module_name = module_name,
      genes = genes,
      comparisons = comparison_keys
    )
}

full_module_heatmap_long <- if (
  length(full_module_heatmap_df_list) > 0L
) {
  dplyr::bind_rows(
    full_module_heatmap_df_list,
    .id = "module_list_name"
  )
} else {
  tibble::tibble()
}

data.table::fwrite(
  full_module_heatmap_long,
  file.path(
    subdirs$tables,
    "FULL_SIGNATURE_LOG2FC_HEATMAP_LONG_TABLE_ALL_SIX.csv"
  )
)

full_module_plots <- list()

for (module_name in names(full_module_heatmap_df_list)) {
  n_genes <- dplyr::n_distinct(
    full_module_heatmap_df_list[[module_name]]$Gene
  )

  p_one <- build_one_row_module_heatmap(
    df_list = full_module_heatmap_df_list[module_name],
    comparisons = comparison_keys,
    module_order_used = module_name,
    selected_modules = module_name,
    show_legend = TRUE,
    force_block_height = n_genes
  )

  full_module_plots[[module_name]] <- p_one

  height_mm <- max(
    110,
    min(360, 50 + 1.55 * n_genes)
  )

  save_plot_set(
    p_one,
    subdirs$full_modules,
    paste0(
      "FULL_SIGNATURE_",
      safe_filename(module_name),
      "_LOG2FC_HEATMAP_ALL_SIX_COMPARISONS"
    ),
    six_comparison_heatmap_width_mm,
    height_mm
  )
}

if ("Gastric metaplasia" %in% names(full_module_plots)) {
  save_plot_set(
    full_module_plots[["Gastric metaplasia"]],
    subdirs$full_modules,
    "FINAL_24_GENE_GASTRIC_METAPLASIA_LOG2FC_HEATMAP_ALL_SIX_COMPARISONS",
    six_comparison_heatmap_width_mm,
    145
  )
}

select_modules_for_group <- function(group_name) {
  if (!group_name %in% names(comparisons_by_group)) {
    stop("Unknown tissue group: ", group_name)
  }

  reference_comparison <- reference_comparisons[[group_name]]
  reference_de <- heatmap_de_list[[reference_comparison]]

  candidate_qc <- purrr::imap_dfr(
    module_gene_sets,
    function(genes, module_name) {
      cleaned_genes <- unique(clean_gene_symbols(genes))
      cleaned_genes <- cleaned_genes[
        !is.na(cleaned_genes) & cleaned_genes != ""
      ]

      tibble::tibble(
        group = group_name,
        module = module_name,
        Gene = cleaned_genes,
        signature_order = seq_along(cleaned_genes)
      ) |>
        dplyr::left_join(
          reference_de |>
            dplyr::select(
              Gene,
              reference_logFC = logFC,
              reference_logCPM = logCPM,
              reference_PValue = PValue,
              reference_FDR = FDR
            ),
          by = "Gene"
        ) |>
        dplyr::mutate(
          reference_comparison = reference_comparison,
          selection_rule = dplyr::if_else(
            module == "Gastric metaplasia",
            paste0(
              "All final 24 gastric-metaplasia genes; ",
              "no logCPM or logFC threshold"
            ),
            paste0(
              "Reference comparison ",
              reference_comparison,
              ": logCPM >= ",
              module_reference_logcpm_cutoff,
              " and logFC > ",
              module_reference_logfc_cutoff,
              "; top ",
              module_max_genes,
              " by decreasing logFC"
            )
          ),
          passes_reference_filter = dplyr::case_when(
            module == "Gastric metaplasia" ~ TRUE,
            TRUE ~ is.finite(reference_logCPM) &
              reference_logCPM >= module_reference_logcpm_cutoff &
              is.finite(reference_logFC) &
              reference_logFC > module_reference_logfc_cutoff
          )
        )
    }
  )

  selected <- purrr::map_dfr(
    module_order,
    function(module_name) {
      dat <- candidate_qc |>
        dplyr::filter(module == module_name)

      if (module_name == "Gastric metaplasia") {
        dat |>
          dplyr::arrange(signature_order) |>
          dplyr::mutate(selection_rank = dplyr::row_number())
      } else {
        dat |>
          dplyr::filter(passes_reference_filter) |>
          dplyr::arrange(
            dplyr::desc(reference_logFC),
            Gene
          ) |>
          dplyr::slice_head(n = module_max_genes) |>
          dplyr::mutate(selection_rank = dplyr::row_number())
      }
    }
  )

  missing_modules <- setdiff(
    setdiff(module_order, "Gastric metaplasia"),
    unique(selected$module)
  )

  if (length(missing_modules) > 0L) {
    warning(
      group_name,
      ": no genes passed the strict reference rule for ",
      paste(missing_modules, collapse = ", "),
      ". Full-signature heatmaps are still written; the strict threshold ",
      "was not silently relaxed.",
      call. = FALSE
    )
  }

  selected_gene_groups <- split(
    selected$Gene,
    selected$module
  )

  for (module_name in module_order) {
    if (!module_name %in% names(selected_gene_groups)) {
      selected_gene_groups[[module_name]] <- character(0)
    }
  }

  selected_gene_groups <- selected_gene_groups[module_order]

  list(
    group = group_name,
    comparisons = comparisons_by_group[[group_name]],
    reference_comparison = reference_comparison,
    candidate_qc = candidate_qc,
    selected = selected,
    gene_groups = selected_gene_groups
  )
}

selected_module_results_by_group <- lapply(
  names(comparisons_by_group),
  select_modules_for_group
)
names(selected_module_results_by_group) <- names(comparisons_by_group)

selected_module_candidate_qc_all <- dplyr::bind_rows(
  lapply(
    selected_module_results_by_group,
    `[[`,
    "candidate_qc"
  )
)

selected_module_genes_all <- dplyr::bind_rows(
  lapply(
    selected_module_results_by_group,
    `[[`,
    "selected"
  )
)

data.table::fwrite(
  selected_module_candidate_qc_all,
  file.path(
    subdirs$qc,
    "GROUP_SPECIFIC_MODULE_CANDIDATE_AND_REFERENCE_SELECTION_QC.csv"
  )
)

data.table::fwrite(
  selected_module_genes_all,
  file.path(
    subdirs$tables,
    "GROUP_SPECIFIC_GENES_DISPLAYED_IN_SELECTED_HEATMAPS.csv"
  )
)

selected_group_df_lists <- list()
selected_group_panels <- list()
selected_group_individual_plots <- list()
selected_group_long_tables <- list()

for (group_name in names(selected_module_results_by_group)) {
  selection_result <- selected_module_results_by_group[[group_name]]
  group_comparisons <- selection_result$comparisons
  df_list <- list()

  for (module_name in module_order) {
    genes <- selection_result$gene_groups[[module_name]]
    if (length(genes) == 0L) next

    df_list[[module_name]] <- make_module_heatmap_df(
      module_name = module_name,
      genes = genes,
      comparisons = group_comparisons
    )
  }

  selected_group_df_lists[[group_name]] <- df_list

  selected_group_long_tables[[group_name]] <- if (
    length(df_list) > 0L
  ) {
    dplyr::bind_rows(df_list, .id = "module_list_name") |>
      dplyr::mutate(tissue_group = group_name, .before = 1)
  } else {
    tibble::tibble()
  }

  p_base <- build_one_row_module_heatmap(
    df_list = df_list,
    comparisons = group_comparisons,
    module_order_used = module_order,
    show_legend = TRUE
  )

  p_group <- patchwork::wrap_plots(p_base) +
    patchwork::plot_annotation(
      title = paste0(
        group_display_labels[[group_name]],
        " | group-specific ES-reference gene selection"
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          family = font_family,
          face = "bold",
          size = 8,
          hjust = 0
        )
      )
    )

  selected_group_panels[[group_name]] <- p_group

  save_plot_set(
    p_group,
    subdirs$selected_modules,
    paste0(
      "SELECTED_IFNA_INFLAMMATORY_TNFA_GM_ONE_ROW_",
      group_name
    ),
    183,
    125
  )

  if (length(df_list) > 0L) {
    common_block_height <- max(
      vapply(
        df_list,
        function(x) dplyr::n_distinct(x$Gene),
        integer(1)
      )
    )

    group_individual <- list()

    for (module_name in names(df_list)) {
      p_one <- build_one_row_module_heatmap(
        df_list = df_list[module_name],
        comparisons = group_comparisons,
        module_order_used = module_name,
        selected_modules = module_name,
        show_legend = TRUE,
        force_block_height = common_block_height
      )

      group_individual[[module_name]] <- p_one

      individual_height <- if (
        module_name == "Gastric metaplasia"
      ) {
        125
      } else {
        105
      }

      save_plot_set(
        p_one,
        subdirs$selected_modules,
        paste0(
          "SELECTED_",
          group_name,
          "_SINGLE_MODULE_",
          safe_filename(module_name),
          "_LOG2FC_HEATMAP"
        ),
        82,
        individual_height
      )
    }

    selected_group_individual_plots[[group_name]] <- group_individual
  }
}

selected_group_long_all <- dplyr::bind_rows(
  selected_group_long_tables
)

data.table::fwrite(
  selected_group_long_all,
  file.path(
    subdirs$tables,
    "GROUP_SPECIFIC_SELECTED_MODULE_LOG2FC_HEATMAP_LONG_TABLE.csv"
  )
)

p_selected_modules_two_rows <-
  selected_group_panels[["SSL_HP"]] /
  selected_group_panels[["H_ADJ_SSL_HP"]] +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = paste0(
      "Selected gene-level log2FC heatmaps: ",
      "SSL_HP and H_ADJ_SSL_HP"
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
  p_selected_modules_two_rows,
  subdirs$combined,
  "FINAL_SELECTED_GENE_MODULE_HEATMAPS_SSLHP_AND_HADJ_TWO_ROWS",
  183,
  250
)

selected_union_gene_lists <- list()
selected_union_qc <- list()

for (module_name in module_order) {
  selected_one <- selected_module_genes_all |>
    dplyr::filter(module == module_name) |>
    dplyr::arrange(
      factor(group, levels = names(comparisons_by_group)),
      selection_rank,
      Gene
    ) |>
    dplyr::distinct(Gene, .keep_all = TRUE)

  if (module_name == "Gastric metaplasia") {
    selected_one <- selected_one |>
      dplyr::arrange(signature_order)
  }

  selected_union_gene_lists[[module_name]] <- selected_one$Gene
  selected_union_qc[[module_name]] <- selected_one
}

selected_union_qc_all <- dplyr::bind_rows(
  selected_union_qc
)

data.table::fwrite(
  selected_union_qc_all,
  file.path(
    subdirs$tables,
    "UNION_OF_GROUP_SPECIFIC_SELECTED_GENES_FOR_ALL_SIX_HEATMAP.csv"
  )
)

selected_union_df_list <- list()

for (module_name in module_order) {
  genes <- selected_union_gene_lists[[module_name]]
  if (length(genes) == 0L) next

  selected_union_df_list[[module_name]] <- make_module_heatmap_df(
    module_name = module_name,
    genes = genes,
    comparisons = comparison_keys
  )
}

p_selected_union_all6 <- build_one_row_module_heatmap(
  df_list = selected_union_df_list,
  comparisons = comparison_keys,
  module_order_used = module_order,
  show_legend = TRUE
)

save_plot_set(
  p_selected_union_all6,
  subdirs$selected_modules,
  "SELECTED_UNION_IFNA_INFLAMMATORY_TNFA_GM_ALL_SIX_COMPARISONS",
  220,
  145
)

p_final_gsea_summary <- p_hallmark_overlays_final /
  p_metaplasia_curves_2x3 +
  patchwork::plot_layout(heights = c(1, 1)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = paste0(
      "Hallmark and ",
      metaplasia_display_title_wrapped,
      " across SSL_HP and H_ADJ_SSL_HP"
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
  p_final_gsea_summary,
  subdirs$combined,
  "FINAL_GSEA_SUMMARY_HALLMARK_OVERLAYS_AND_GM_ALL_SIX",
  183,
  315
)

p_final_heatmap_summary <- p_hallmark_all6_selected_three /
  p_metaplasia_nes_heatmap /
  p_selected_modules_two_rows +
  patchwork::plot_layout(heights = c(0.55, 0.42, 1.65)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = paste0(
      "Pathway-level and gene-level summaries across all six ",
      "bacterial comparisons"
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
  p_final_heatmap_summary,
  subdirs$combined,
  "FINAL_HEATMAP_SUMMARY_SELECTED_HALLMARK_GM_AND_GENE_MODULES_ALL_SIX",
  183,
  395
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
  "for the four-module heatmap panels..."
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

four_module_candidate_qc <- purrr::map_dfr(
  names(comparisons_by_group),
  function(group_name) {
    reference_comparison <- reference_comparisons[[group_name]]
    reference_de <- heatmap_de_list[[reference_comparison]] |>
      dplyr::select(
        Gene,
        reference_logFC = logFC,
        reference_logCPM = logCPM,
        reference_PValue = PValue,
        reference_FDR = FDR
      )

    purrr::imap_dfr(
      four_module_gene_sets,
      function(genes, module_name) {
        genes <- clean_gene_symbols(genes)
        genes <- unique(genes[!is.na(genes) & genes != ""])

        tibble::tibble(
          group = group_name,
          reference_comparison = reference_comparison,
          module = module_name,
          Gene = genes,
          source_order = seq_along(genes)
        ) |>
          dplyr::left_join(
            reference_de,
            by = "Gene"
          ) |>
          dplyr::mutate(
            present_in_reference_DE = is.finite(reference_logFC),
            passes_selection_filter =
              present_in_reference_DE &
              is.finite(reference_logCPM) &
              reference_logCPM >= module_reference_logcpm_cutoff &
              is.finite(reference_logFC) &
              reference_logFC > module_reference_logfc_cutoff,
            selection_rule = paste0(
              "Reference ", reference_comparison,
              ": logCPM >= ",
              module_reference_logcpm_cutoff,
              " and logFC > ",
              module_reference_logfc_cutoff,
              "; top ",
              module_max_genes,
              " by decreasing logFC"
            )
          )
      }
    )
  }
)

four_module_selected_genes <- purrr::map_dfr(
  names(comparisons_by_group),
  function(group_name) {
    purrr::map_dfr(
      four_module_order,
      function(module_name) {
        four_module_candidate_qc |>
          dplyr::filter(
            group == group_name,
            module == module_name,
            passes_selection_filter
          ) |>
          dplyr::arrange(
            dplyr::desc(reference_logFC),
            reference_FDR,
            Gene
          ) |>
          dplyr::slice_head(n = module_max_genes) |>
          dplyr::mutate(
            selection_rank = dplyr::row_number()
          )
      }
    )
  }
)

expected_four_module_blocks <- tidyr::expand_grid(
  group = names(comparisons_by_group),
  module = four_module_order
)

observed_four_module_blocks <- four_module_selected_genes |>
  dplyr::distinct(group, module)

missing_four_module_blocks <- expected_four_module_blocks |>
  dplyr::anti_join(
    observed_four_module_blocks,
    by = c("group", "module")
  )

if (nrow(missing_four_module_blocks) > 0L) {
  warning(
    "No genes passed the strict four-module reference rule for: ",
    paste0(
      missing_four_module_blocks$group,
      " / ",
      missing_four_module_blocks$module,
      collapse = "; "
    ),
    ". Grey placeholders are retained and the threshold is not relaxed.",
    call. = FALSE
  )
}

data.table::fwrite(
  four_module_candidate_qc,
  file.path(
    subdirs$qc,
    paste0(
      "FOUR_MODULE_PANEL_GROUP_SPECIFIC_REFERENCE_",
      "CANDIDATE_QC.csv"
    )
  )
)

data.table::fwrite(
  four_module_selected_genes,
  file.path(
    subdirs$tables,
    paste0(
      "FOUR_MODULE_PANEL_GROUP_SPECIFIC_SELECTED_",
      "GENES_LOGCPM_GE_2_LOGFC_GT_1_TOP12.csv"
    )
  )
)

make_four_module_gene_lists_for_group <- function(group_name) {
  selected_one <- four_module_selected_genes |>
    dplyr::filter(group == group_name)

  gene_lists <- split(
    selected_one$Gene,
    selected_one$module
  )

  for (module_name in four_module_order) {
    if (!module_name %in% names(gene_lists)) {
      gene_lists[[module_name]] <- character(0)
    }
  }

  gene_lists[four_module_order]
}

make_four_module_heatmap_df_with_placeholder <- function(
    module_name,
    genes,
    comparisons
) {
  genes <- clean_gene_symbols(genes)
  genes <- unique(genes[!is.na(genes) & genes != ""])

  if (length(genes) == 0L) {
    genes <- "NO GENES PASSED"
  }

  make_module_heatmap_df(
    module_name = module_name,
    genes = genes,
    comparisons = comparisons
  ) |>
    dplyr::mutate(
      placeholder = Gene == "NO GENES PASSED"
    )
}

four_module_gene_lists_by_group <- lapply(
  names(comparisons_by_group),
  make_four_module_gene_lists_for_group
)
names(four_module_gene_lists_by_group) <-
  names(comparisons_by_group)

four_module_df_lists_by_group <- list()
four_module_panels_by_group <- list()
four_module_display_long_by_group <- list()

for (group_name in names(comparisons_by_group)) {
  group_comparisons <- comparisons_by_group[[group_name]]
  gene_lists <- four_module_gene_lists_by_group[[group_name]]

  df_list <- lapply(
    four_module_order,
    function(module_name) {
      make_four_module_heatmap_df_with_placeholder(
        module_name = module_name,
        genes = gene_lists[[module_name]],
        comparisons = group_comparisons
      )
    }
  )
  names(df_list) <- four_module_order

  four_module_df_lists_by_group[[group_name]] <- df_list

  four_module_display_long_by_group[[group_name]] <-
    dplyr::bind_rows(df_list, .id = "module_list_name") |>
      dplyr::mutate(tissue_group = group_name, .before = 1)

  p_group_four_module <- build_one_row_module_heatmap(
    df_list = df_list,
    comparisons = group_comparisons,
    module_order_used = four_module_order,
    show_legend = TRUE,
    legend_title = "log2FC"
  )

  four_module_panels_by_group[[group_name]] <-
    p_group_four_module

  save_plot_set(
    p_group_four_module,
    subdirs$four_module,
    paste0(
      "FINAL_ALL4_IFNA_CHEMOKINES_ANTIMICROBIAL_TNFA_",
      "LOG2FC_",
      group_name,
      "_SERR_FUSO_ES_NATURE_EXACT"
    ),
    183,
    82
  )

  common_group_block_height <- max(
    vapply(
      df_list,
      function(x) dplyr::n_distinct(x$Gene),
      integer(1)
    ),
    1L
  )

  for (module_name in four_module_order) {
    p_one_four_module <- build_one_row_module_heatmap(
      df_list = df_list[module_name],
      comparisons = group_comparisons,
      module_order_used = module_name,
      selected_modules = module_name,
      show_legend = TRUE,
      legend_title = "log2FC",
      force_block_height = common_group_block_height
    )

    save_plot_set(
      p_one_four_module,
      subdirs$four_module,
      paste0(
        "FINAL_SINGLE_MODULE_",
        group_name,
        "_",
        safe_filename(module_name),
        "_SERR_FUSO_ES_NATURE_EXACT"
      ),
      82,
      105
    )
  }
}

four_module_display_long_group_specific <- dplyr::bind_rows(
  four_module_display_long_by_group
)

data.table::fwrite(
  four_module_display_long_group_specific,
  file.path(
    subdirs$tables,
    paste0(
      "FOUR_MODULE_PANEL_GROUP_SPECIFIC_DISPLAYED_",
      "LOG2FC_VALUES.csv"
    )
  )
)

p_four_module_two_tissue_rows <-
  four_module_panels_by_group[["SSL_HP"]] /
  four_module_panels_by_group[["H_ADJ_SSL_HP"]] +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = paste0(
      "Interferon-alpha, chemokine, antimicrobial and ",
      "TNFA/NFKB modules"
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
  p_four_module_two_tissue_rows,
  subdirs$four_module,
  paste0(
    "FINAL_ALL4_IFNA_CHEMOKINES_ANTIMICROBIAL_TNFA_",
    "SSLHP_AND_HADJ_TWO_ROWS_NATURE_EXACT"
  ),
  183,
  165
)

save_plot_set(
  p_four_module_two_tissue_rows,
  subdirs$combined,
  paste0(
    "FINAL_ALL4_IFNA_CHEMOKINES_ANTIMICROBIAL_TNFA_",
    "SSLHP_AND_HADJ_TWO_ROWS_NATURE_EXACT"
  ),
  183,
  165
)

four_module_union_selected_genes <- four_module_selected_genes |>
  dplyr::mutate(
    group_sort = match(group, names(comparisons_by_group)),
    module_sort = match(module, four_module_order)
  ) |>
  dplyr::arrange(
    module_sort,
    group_sort,
    selection_rank,
    Gene
  ) |>
  dplyr::distinct(module, Gene, .keep_all = TRUE)

four_module_union_gene_lists <- split(
  four_module_union_selected_genes$Gene,
  four_module_union_selected_genes$module
)

for (module_name in four_module_order) {
  if (!module_name %in% names(four_module_union_gene_lists)) {
    four_module_union_gene_lists[[module_name]] <- character(0)
  }
}

four_module_union_gene_lists <-
  four_module_union_gene_lists[four_module_order]

four_module_all6_df_list <- lapply(
  four_module_order,
  function(module_name) {
    make_four_module_heatmap_df_with_placeholder(
      module_name = module_name,
      genes = four_module_union_gene_lists[[module_name]],
      comparisons = comparison_keys
    )
  }
)
names(four_module_all6_df_list) <- four_module_order

four_module_all6_display_long <- dplyr::bind_rows(
  four_module_all6_df_list,
  .id = "module_list_name"
)

data.table::fwrite(
  four_module_union_selected_genes,
  file.path(
    subdirs$tables,
    paste0(
      "FOUR_MODULE_PANEL_UNION_OF_GROUP_SPECIFIC_",
      "SELECTED_GENES_FOR_ALL_SIX.csv"
    )
  )
)

data.table::fwrite(
  four_module_all6_display_long,
  file.path(
    subdirs$tables,
    paste0(
      "FOUR_MODULE_PANEL_ALL_SIX_DISPLAYED_",
      "LOG2FC_VALUES.csv"
    )
  )
)

p_four_module_all6 <- build_one_row_module_heatmap(
  df_list = four_module_all6_df_list,
  comparisons = comparison_keys,
  module_order_used = four_module_order,
  show_legend = TRUE,
  legend_title = "log2FC"
)

four_module_all6_stem <- paste0(
  "FINAL_ALL4_IFNA_CHEMOKINES_ANTIMICROBIAL_TNFA_",
  "LOG2FC_ALL_SIX_SERR_FUSO_ES_NATURE_EXACT"
)

save_plot_set(
  p_four_module_all6,
  subdirs$four_module,
  four_module_all6_stem,
  183,
  82
)

save_plot_set(
  p_four_module_all6,
  subdirs$four_module,
  paste0(four_module_all6_stem, "_WIDE_220MM"),
  220,
  105
)

save_plot_set(
  p_four_module_all6,
  subdirs$combined,
  four_module_all6_stem,
  183,
  82
)

writeLines(
  c(
    "ORGANOID-MATCHED FOUR-MODULE SIX-COMPARISON HEATMAP PANELS",
    "==========================================================",
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
    "Selection:",
    paste0(
      "  SSL_HP reference: ",
      reference_comparisons[["SSL_HP"]]
    ),
    paste0(
      "  H_ADJ_SSL_HP reference: ",
      reference_comparisons[["H_ADJ_SSL_HP"]]
    ),
    paste0(
      "  Rule: reference logCPM >= ",
      module_reference_logcpm_cutoff,
      " and logFC > ",
      module_reference_logfc_cutoff,
      "; top ",
      module_max_genes,
      " per module by decreasing logFC."
    ),
    "  The same genes are displayed across Serr, Fuso and ES within each tissue compartment.",
    "  The all-six panel uses the union of both compartment-specific selections.",
    "  Empty modules are shown as grey placeholders; thresholds are not relaxed.",
    "",
    "Appearance:",
    paste0("  Font: ", font_family),
    paste0("  Labels: ", shared_label_pt, " pt"),
    paste0(
      "  Standard NES labels: ",
      nes_x_label_angle,
      " degrees; dense gene-heatmap labels: ",
      module_x_label_angle,
      " degrees"
    ),
    paste0(
      "  Colour scale: ", negative_colour,
      " / white / ", positive_colour
    ),
    paste0(
      "  Saturation: -", logfc_colour_limit,
      " to +", logfc_colour_limit, " log2FC"
    ),
    "  White tile borders and right-side vertical log2FC legend.",
    "  Exact-size main exports: 183 x 82 mm."
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
  dplyr::arrange(
    match(comparison, comparison_keys),
    filter_mode
  )

data.table::fwrite(
  rank_method_summary,
  file.path(
    subdirs$qc,
    "FINAL_RANK_METHOD_SUMMARY_ALL_SIX.csv"
  )
)

writeLines(
  c(
    "FINAL SIX-COMPARISON GSEA AND HEATMAP RUN",
    "=========================================",
    "",
    paste0("Run date: ", Sys.time()),
    paste0("Input directory: ", input_dir),
    paste0("Output directory: ", out_dir),
    "",
    "DE tables analysed:",
    paste0(
      "  ",
      de_file_table$group,
      " | ",
      de_file_table$organism,
      " | ",
      de_file_table$file
    ),
    "",
    "Coverage and display annotations:",
    "  SSL_HP and H_ADJ_SSL_HP are both included in every applicable analysis.",
    "  The exact plot annotations are SSL_HP and H_ADJ_SSL_HP throughout.",
    "  Left-to-right plot, heatmap and legend order within each compartment is Serr, Fuso, ES.",
    "  Each compartment includes Serratia, Fusobacterium and Escherichia-Shigella positive-versus-negative DE tables.",
    "  ES remains the pre-specified reference comparison for strict gene selection; display reordering does not change that method.",
    "",
    "Contrast orientation:",
    "  Every table is assumed to represent bacterial-positive minus bacterial-negative.",
    "  Positive logFC/NES therefore means higher or enriched in bacterial-positive samples.",
    "  Verify this in the original DE model/design.",
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
    "  The complete Hallmark collection is tested separately for each of the six comparisons.",
    "  BH FDR is calculated within the Hallmark collection separately for each comparison.",
    "  The final 2 x 3 overlay panel contains SSL_HP in row 1 and H_ADJ_SSL_HP in row 2.",
    "",
    "Final gastric-metaplasia GSEA:",
    paste0("  Signature source: ", metaplasia_signature_file),
    paste0("  Required and used genes: ", length(metaplasia_genes), "."),
    "  Gene-universe filter: none; ALL_VALID_GENES is used.",
    "  One custom pathway is tested, so padj equals pval; report this as nominal P, not Hallmark FDR.",
    "  Six individual curves, two three-comparison group rows, two group overlays and one 2 x 3 panel are exported.",
    "",
    "Hallmark NES heatmaps:",
    "  All-Hallmark, top-Hallmark, key-Hallmark and selected IFNA/inflammatory/TNFA heatmaps contain all six comparisons.",
    "  Stars on Hallmark heatmaps use Hallmark-family BH FDR.",
    "  Stars on the custom gastric-metaplasia NES heatmap use nominal P.",
    "",
    "Full-signature gene-level log2FC heatmaps:",
    "  Separate all-six heatmaps are written for IFNA, inflammatory response, TNFA/NFKB and gastric metaplasia.",
    "  All 24 gastric-metaplasia genes are displayed in supplied signature order, including grey tiles where absent.",
    "  Hallmark genes are ordered by maximum absolute log2FC across the six comparisons.",
    "  These full-signature outputs ensure a pathway is not lost merely because no gene passes the stricter publication threshold.",
    "",
    "Group-specific publication-selected gene heatmaps:",
    paste0(
      "  SSL_HP reference comparison: ",
      reference_comparisons[["SSL_HP"]],
      "."
    ),
    paste0(
      "  H_ADJ_SSL_HP reference comparison: ",
      reference_comparisons[["H_ADJ_SSL_HP"]],
      "."
    ),
    paste0(
      "  Non-GM rule: reference logCPM >= ",
      module_reference_logcpm_cutoff,
      " and logFC > ",
      module_reference_logfc_cutoff,
      "; top ",
      module_max_genes,
      " by decreasing logFC."
    ),
    "  The selected genes are displayed unchanged across Serr, Fuso and ES within that tissue compartment.",
    "  The strict threshold is never silently relaxed; full-signature heatmaps remain available when a module is empty.",
    "  Gastric metaplasia always uses all 24 genes without an expression or effect-size threshold.",
    "",
    "Gene concordance tables:",
    "  Hallmark direction summaries use logCPM >= 2, matching the Hallmark GSEA universe.",
    "  Gastric-metaplasia direction summaries use all valid genes without an expression filter.",
    "  Per-gene counts identify upregulation across all six, within SSL_HP and within H_ADJ_SSL_HP.",
    "  Pairwise tables report Spearman logFC correlation, same-direction fraction and Jaccard overlap of upregulated genes.",
    "",
    "Additional organoid-matched four-module panels:",
    "  Modules are IFNA, chemokines, antimicrobial peptides and TNFA/NFKB in one horizontal row.",
    "  Canonical chemokines come from GOMF_CHEMOKINE_ACTIVITY; antimicrobial genes come from REACTOME_ANTIMICROBIAL_PEPTIDES.",
    paste0(
      "  Within each tissue compartment, the ES reference uses logCPM >= ",
      module_reference_logcpm_cutoff,
      " and logFC > ",
      module_reference_logfc_cutoff,
      "; top ",
      module_max_genes,
      " by decreasing logFC."
    ),
    "  Group-specific panels show Serr, Fuso and ES; the all-six panel uses the union of both group-specific selections.",
    "",
    "Nature-exact appearance:",
    paste0("  Font: ", font_family, "."),
    paste0("  Gene/comparison labels: ", shared_label_pt, " pt."),
    paste0(
      "  Standard NES comparison labels: ",
      nes_x_label_angle,
      " degrees; dense gene-level comparison labels: ",
      module_x_label_angle,
      " degrees."
    ),
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
    "Important interpretation:",
    "  Separate GSEA analyses do not formally test whether one bacterial association is stronger than another.",
    "  Differences in NES between taxa or tissue compartments are descriptive unless tested in a joint interaction model.",
    "  These are association analyses and do not establish microbial causation."
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
message("  02_Tables/SIX_COMPARISONS_FINAL_GASTRIC_METAPLASIA_FGSEA_ALL_VALID_GENES.csv")
message("  03_Hallmark_GSEA_overlays/")
message("  04_Gastric_metaplasia_GSEA_curves/")
message("  05_Hallmark_NES_heatmaps/")
message("  06_Full_signature_gene_log2FC_heatmaps/")
message("  07_Selected_gene_log2FC_heatmaps/")
message("  08_Combined_figures/")
message("  10_Gene_concordance/")
message("  11_Four_module_heatmap_panels/")
