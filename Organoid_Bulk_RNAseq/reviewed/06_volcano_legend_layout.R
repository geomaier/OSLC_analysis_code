# Remake the Chapter 5 infection and Chapter 6 BRAF volcano plots.
# Reads the four existing volcano CSVs from 07_Tables.
# Moves both legends below the plots, retaining the saved points and labels.
# No differential-expression analysis, FDR adjustment or GSEA is run.

local({
  results_dir <- Sys.getenv("OSLC_ORGANOID_RESULTS_DIR")
  tables_dir <- Sys.getenv("OSLC_ORGANOID_TABLES_DIR")
  if (!nzchar(results_dir) && !nzchar(tables_dir)) {
    stop("Set OSLC_ORGANOID_RESULTS_DIR to the completed organoid run, or ",
         "OSLC_ORGANOID_TABLES_DIR to its 07_Tables directory.", call. = FALSE)
  }

  required_packages <- c("ggplot2", "ggrepel")
  missing_packages <- required_packages[!vapply(
    required_packages, requireNamespace, logical(1), quietly = TRUE
  )]
  if (length(missing_packages)) {
    stop(
      "Install the missing plotting packages first:\ninstall.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))\nThen Source this file again.", call. = FALSE
    )
  }

  input_names <- c(
    infection_data = "Infection_volcano_plot_data_logCPM_gt_2.csv",
    infection_labels = "Infection_volcano_labelled_genes_top10_by_PValue.csv",
    braf_data = "BRAF_SSL_HP_only_volcano_plot_data_logCPM_gt_2.csv",
    braf_labels = "BRAF_SSL_HP_only_volcano_labelled_genes_top10_by_PValue.csv"
  )

  if (!nzchar(tables_dir)) {
    original_tables <- file.path(results_dir, "07_Tables")
    if (all(file.exists(file.path(original_tables, input_names)))) {
      tables_dir <- original_tables
    } else if (interactive()) {
      message(
        "Select any one of the four saved volcano CSVs in your original ",
        "07_Tables folder. Both plots will use files from that same folder."
      )
      selected_file <- file.choose()
      if (!basename(selected_file) %in% input_names) {
        stop(
          "Select one of these saved plotting files:\n",
          paste(input_names, collapse = "\n"), call. = FALSE
        )
      }
      tables_dir <- dirname(selected_file)
    } else {
      stop(
        "Set OSLC_ORGANOID_TABLES_DIR to the existing 07_Tables folder, ",
        "then run again.", call. = FALSE
      )
    }
  }

  tables_dir <- normalizePath(tables_dir, winslash = "/", mustWork = TRUE)
  input_paths <- setNames(file.path(tables_dir, input_names), names(input_names))
  missing_files <- input_paths[!file.exists(input_paths)]
  if (length(missing_files)) {
    stop(
      "These saved plotting files are missing from the selected folder:\n",
      paste(basename(missing_files), collapse = "\n"),
      "\nUse the 07_Tables folder from the complete original output.",
      call. = FALSE
    )
  }

  read_saved <- function(path, required_columns) {
    x <- utils::read.csv(
      path, stringsAsFactors = FALSE, check.names = FALSE,
      na.strings = c("", "NA")
    )
    missing_columns <- setdiff(required_columns, names(x))
    if (length(missing_columns)) {
      stop(
        basename(path), " is missing columns: ",
        paste(missing_columns, collapse = ", "), call. = FALSE
      )
    }
    if (!nrow(x)) stop("Empty saved table: ", basename(path), call. = FALSE)
    for (column in c("logFC", "neg_log10_FDR")) {
      if (!is.numeric(x[[column]]) || any(!is.finite(x[[column]]))) {
        stop(
          basename(path), ": ", column,
          " must contain finite numeric coordinates.", call. = FALSE
        )
      }
    }
    x
  }

  regulation_levels <- c("Upregulated", "Downregulated", "Not significant")

  read_volcano <- function(data_file, label_file) {
    dat <- read_saved(
      data_file, c("Gene", "logFC", "neg_log10_FDR", "Regulation")
    )
    labelled <- read_saved(
      label_file, c("Gene", "logFC", "neg_log10_FDR")
    )
    if (anyNA(dat$Regulation) ||
        any(!dat$Regulation %in% regulation_levels)) {
      stop("Unrecognised significance categories in ", basename(data_file),
           call. = FALSE)
    }
    if (anyDuplicated(dat$Gene) || anyDuplicated(labelled$Gene) ||
        anyNA(labelled$Gene) || any(!nzchar(labelled$Gene))) {
      stop("Use the original exported volcano and labelled-gene CSVs.",
           call. = FALSE)
    }
    indices <- match(labelled$Gene, dat$Gene)
    if (anyNA(indices) ||
        !isTRUE(all.equal(labelled$logFC, dat$logFC[indices],
                         tolerance = 1e-10)) ||
        !isTRUE(all.equal(labelled$neg_log10_FDR,
                         dat$neg_log10_FDR[indices], tolerance = 1e-10))) {
      stop("The saved label table does not match ", basename(data_file),
           call. = FALSE)
    }
    list(data = dat, labelled = labelled)
  }

  infection <- read_volcano(
    input_paths[["infection_data"]], input_paths[["infection_labels"]]
  )
  braf <- read_volcano(
    input_paths[["braf_data"]], input_paths[["braf_labels"]]
  )
  if (!"CXCL5" %in% braf$labelled$Gene) {
    stop(
      "CXCL5 is missing from the saved BRAF label table. ",
      "Select the original 07_Tables folder used for Figure 6.7B.",
      call. = FALSE
    )
  }

  # Retain the original shared axes for the two volcano plots.
  common_x_limits <- c(-1, 1) * max(
    2.2, ceiling(max(abs(infection$data$logFC), abs(braf$data$logFC)) * 1.05)
  )
  common_y_limits <- c(0, max(
    5, ceiling(max(infection$data$neg_log10_FDR,
                   braf$data$neg_log10_FDR) * 1.05)
  ))

  volcano_cols <- c(
    "Upregulated" = "#D62728",
    "Downregulated" = "#1F77B4",
    "Not significant" = "#8E8E8E"
  )

  original_theme <- ggplot2::theme_classic(base_family = "Arial", base_size = 7) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Arial", colour = "black"),
      plot.title = ggplot2::element_text(
        family = "Arial", face = "bold", size = 7, hjust = 0
      ),
      axis.title = ggplot2::element_text(family = "Arial", size = 7, colour = "black"),
      axis.text = ggplot2::element_text(family = "Arial", size = 7, colour = "black"),
      legend.title = ggplot2::element_text(family = "Arial", size = 6, colour = "black"),
      legend.text = ggplot2::element_text(family = "Arial", size = 6, colour = "black"),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.35),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = grid::unit(1.2, "mm"),
      panel.grid = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(4, 4, 4, 4, "mm")
    )

  build_volcano <- function(saved, title_text) {
    dat <- saved$data
    dat$Regulation <- factor(dat$Regulation, levels = regulation_levels)

    ggplot2::ggplot(dat, ggplot2::aes(x = logFC, y = neg_log10_FDR)) +
      ggplot2::geom_hline(
        yintercept = -log10(0.05), linetype = "dashed",
        colour = "black", linewidth = 0.35
      ) +
      ggplot2::geom_vline(
        xintercept = c(-1, 1), linetype = "dashed",
        colour = "black", linewidth = 0.35
      ) +
      ggplot2::geom_point(
        ggplot2::aes(colour = Regulation), size = 0.9, stroke = 0, alpha = 0.85
      ) +
      ggplot2::scale_colour_manual(values = volcano_cols, drop = FALSE) +
      ggplot2::coord_cartesian(
        xlim = common_x_limits, ylim = common_y_limits, clip = "off"
      ) +
      ggplot2::labs(
        title = title_text,
        subtitle = NULL,
        x = expression(log[2]~fold~change),
        y = expression(-log[10]~adjusted~italic(P)~value),
        colour = NULL
      ) +
      original_theme +
      ggplot2::theme(
        plot.subtitle = ggplot2::element_blank(),
        # Layout correction: the legend sits below, outside the plotting area.
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.key.height = grid::unit(3.5, "mm"),
        legend.key.width = grid::unit(4.5, "mm"),
        legend.background = ggplot2::element_rect(
          fill = ggplot2::alpha("white", 0.9), colour = NA
        )
      ) +
      ggrepel::geom_text_repel(
        data = saved$labelled, ggplot2::aes(label = Gene),
        size = 2.1, family = "Arial", colour = "black",
        box.padding = 0.28, point.padding = 0.12,
        segment.colour = "grey55", segment.linewidth = 0.28,
        min.segment.length = 0, max.overlaps = Inf, seed = 1
      )
  }

  plots <- list(
    Chapter_5_Infection_volcano_legend_below = build_volcano(
      infection, "Infected vs non-infected organoids"
    ),
    Figure_6_7B_BRAF_volcano_legend_below = build_volcano(
      braf, "BRAF V600E vs WT in SSL_HP organoids"
    )
  )

  # Create a fresh output subfolder on every run; preserve all existing figures.
  output_parent <- file.path(dirname(tables_dir), "Volcano_Legend_Fixes_Chapters_5_and_6")
  dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)
  output_dir <- NULL
  for (i in seq_len(9999L)) {
    candidate <- file.path(output_parent, paste0("r", i))
    if (!file.exists(candidate) && dir.create(candidate, showWarnings = FALSE)) {
      output_dir <- candidate
      break
    }
  }
  if (is.null(output_dir)) stop("Could not create a new output folder.", call. = FALSE)

  for (stem in names(plots)) {
    p <- plots[[stem]]
    prefix <- file.path(output_dir, stem)
    # Original figure dimensions and raster resolution.
    ggplot2::ggsave(
      paste0(prefix, ".png"), plot = p,
      width = 90, height = 85, units = "mm", dpi = 600, bg = "white"
    )
    ggplot2::ggsave(
      paste0(prefix, ".tif"), plot = p,
      width = 90, height = 85, units = "mm", dpi = 600,
      device = grDevices::tiff, compression = "lzw", bg = "white"
    )
    if (capabilities("cairo")) {
      ggplot2::ggsave(
        paste0(prefix, ".pdf"), plot = p,
        width = 90, height = 85, units = "mm",
        device = grDevices::cairo_pdf, bg = "white"
      )
    }
    print(p)
  }

  if (!capabilities("cairo")) {
    message("PNG and TIF saved. PDF skipped because Cairo is unavailable.")
  }
  writeLines(c(
    "Chapter 5 and Chapter 6: existing volcano plots with legends below.",
    paste0("Source tables: ", tables_dir),
    paste0("Chapter 5 genes plotted: ", nrow(infection$data)),
    paste0("Chapter 5 labels: ", paste(infection$labelled$Gene, collapse = ", ")),
    paste0("Chapter 6 genes plotted: ", nrow(braf$data)),
    paste0("Chapter 6 labels: ", paste(braf$labelled$Gene, collapse = ", ")),
    paste0("X limits: ", paste(common_x_limits, collapse = ", ")),
    paste0("Y limits: ", paste(common_y_limits, collapse = ", ")),
    "Original colours, thresholds, selected labels, fonts and dimensions retained.",
    "Dimensions: 90 x 85 mm; PNG/TIF: 600 dpi.",
    "No differential-expression analysis, FDR recalculation or GSEA was run."
  ), file.path(output_dir, "Plot_sources.txt"))

  message(
    "Both volcano plots saved in:\n",
    normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  )
})
