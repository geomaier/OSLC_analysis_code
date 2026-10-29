suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(ggrepel); library(gridExtra)
})
font_family <- "Arial"
font_pt <- 7
legend_pt <- 6
label_pt <- 7
export_dpi <- 600
alpha <- 0.05
alpha_level <- 0.05
volcano_width_mm <- 160
volcano_height_mm <- 110
volcano_colours <- c("Up" = "#A52A2A", "Down" = "green4", "NS" = "grey75")
pt_to_gg <- function(x) x / 2.845276
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

safe_write_csv <- function(x, path) {
  write.csv(
    as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE),
    path,
    row.names = FALSE,
    na = ""
  )
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

wrap_taxon_label <- function(x, width = 34L) {
  vapply(
    as.character(x),
    function(z) paste(strwrap(z, width = width), collapse = "\n"),
    character(1)
  )
}
save_ggplot <- function(
    plot,
    stem,
    width_mm,
    height_mm
) {
  pdf_device <- if (capabilities("cairo")) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }

  ggsave(
    file.path(
      out_dir,
      paste0(stem, ".pdf")
    ),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = pdf_device,
    bg = "white"
  )

  ggsave(
    file.path(
      out_dir,
      paste0(stem, ".png")
    ),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = export_dpi,
    bg = "white"
  )

  ggsave(
    file.path(
      out_dir,
      paste0(stem, ".jpg")
    ),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = export_dpi,
    quality = 100,
    bg = "white"
  )
}
replot_main_volcano <- function(result_df, group_1, group_2, filename_tag, out_dir) {
  volcano_df <- result_df %>%
    filter(is.finite(FDR), FDR >= 0, FDR <= 1, is.finite(log2FoldChange)) %>%
    mutate(
      FDR_plot = pmax(FDR, .Machine$double.xmin),

      Label = ifelse(
        is.na(Genus) | Genus == "",
        Taxon,
        Genus
      ),

      Plot_status = case_when(
        Significant_FDR_0_05 &
          is.finite(log2FoldChange) &
          log2FoldChange > 0 ~ "Up",

        Significant_FDR_0_05 &
          is.finite(log2FoldChange) &
          log2FoldChange < 0 ~ "Down",

        TRUE ~ "NS"
      ),

      Plot_status = factor(
        Plot_status,
        levels = c("Up", "Down", "NS")
      )
    )

  label_df <- volcano_df %>%
    filter(
      Plot_status %in% c("Up", "Down"),
      !is.na(FDR),
      FDR < alpha,
      !is.na(Label),
      Label != "",
      is.finite(log2FoldChange),
      is.finite(FDR_plot)
    ) %>%
    arrange(
      FDR,
      P_value,
      desc(abs(log2FoldChange))
    )

  write.csv(
    label_df %>%
      select(
        any_of(
          c(
            "Label",
            "Genus",
            "Taxon",
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "Wald_statistic",
            "P_value",
            "FDR",
            "Plot_status"
          )
        )
      ),
    file.path(
      out_dir,
      paste0(
        "Volcano_DESeq2_",
        filename_tag,
        "_labelled_genera.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )

  p_volcano <- ggplot(
    volcano_df,
    aes(
      x = log2FoldChange,
      y = -log10(FDR_plot)
    )
  ) +
    geom_point(
      aes(colour = Plot_status),
      size = 1.45,
      alpha = 0.85
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = Label),
      family = font_family,
      size = pt_to_gg(label_pt),
      colour = "black",
      max.overlaps = Inf,
      force = 5,
      force_pull = 0.08,
      box.padding = 0.30,
      point.padding = 0.20,
      segment.size = 0.20,
      min.segment.length = 0,
      seed = 1
    ) +
    scale_colour_manual(
      values = volcano_colours,
      breaks = c("Up", "Down", "NS"),
      limits = c("Up", "Down", "NS"),
      drop = FALSE
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0.10, 0.14))
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.14))
    ) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed",
      linewidth = 0.35,
      colour = "grey50"
    ) +
    geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed",
      linewidth = 0.35,
      colour = "grey50"
    ) +
    labs(
      title = paste0(
        group_1,
        " vs ",
        group_2
      ),
      x = "log2 fold change",
      y = "-log10(FDR)",
      colour = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(
      base_family = font_family,
      base_size = font_pt
    ) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",

      legend.text = element_text(
        family = font_family,
        size = legend_pt,
        colour = "black"
      ),

      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.key.size = grid::unit(3.2, "mm"),
      legend.spacing.x = grid::unit(2, "mm"),

      plot.title = element_text(
        family = font_family,
        face = "bold",
        size = font_pt,
        colour = "black",
        hjust = 0.5
      ),

      axis.title.x = element_text(
        family = font_family,
        size = font_pt,
        colour = "black",
        margin = margin(t = 4)
      ),

      axis.title.y = element_text(
        family = font_family,
        size = font_pt,
        colour = "black",
        margin = margin(r = 4)
      ),

      axis.text.x = element_text(
        family = font_family,
        size = font_pt,
        colour = "black"
      ),

      axis.text.y = element_text(
        family = font_family,
        size = font_pt,
        colour = "black"
      ),

      axis.line = element_line(
        linewidth = 0.35,
        colour = "black"
      ),

      axis.ticks = element_line(
        linewidth = 0.35,
        colour = "black"
      ),

      axis.ticks.length = grid::unit(1.2, "mm"),

      plot.background = element_rect(
        fill = "white",
        colour = NA
      ),

      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),

      plot.margin = margin(
        5,
        10,
        5,
        6,
        "mm"
      )
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(
          size = 2.2,
          alpha = 1
        )
      )
    )

  save_ggplot(
    p_volcano,
    paste0(
      "Volcano_DESeq2_",
      filename_tag
    ),
    width_mm = volcano_width_mm,
    height_mm = volcano_height_mm
  )

  invisible(p_volcano)
}

make_aligned_volcano <- function(
    result_df,
    rank_label,
    output_directory,
    filename_stem
) {
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
locate_saved_run <- function(work_dir, run_name, required) {
  candidates <- list.dirs(work_dir, recursive = TRUE, full.names = TRUE)
  candidates <- candidates[basename(candidates) == run_name]
  candidates <- candidates[vapply(candidates, function(p) all(file.exists(file.path(p, required))), logical(1))]
  if (length(candidates) > 1L) {
    fingerprints <- vapply(candidates, function(p)
      paste(unname(tools::md5sum(file.path(p, required))), collapse = ":"), character(1))
    if (length(unique(fingerprints)) == 1L) candidates <- sort(candidates)[1L]
  }
  if (length(candidates) != 1L) stop("Expected one saved run named ", run_name,
    " containing ", paste(required, collapse = ", "), "; found ", length(candidates),
    ". Set its exact path in the configuration.")
  normalizePath(candidates, mustWork = TRUE)
}

audit_saved_deseq <- function(path, results, output_file) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("DESeq2 is needed to read the saved model.")
  dds <- readRDS(path)
  convergence <- S4Vectors::mcols(dds)$betaConv
  if (is.null(convergence)) stop("Saved DESeq2 model has no betaConv field: ", path)
  audit <- data.frame(Taxon = rownames(dds), Coefficient_converged = as.logical(convergence))
  audit$Reported_significant <- results$Significant_FDR_0_05[match(audit$Taxon, results$Taxon)]
  write.csv(audit, output_file, row.names = FALSE, na = "")
  selected <- !is.na(audit$Reported_significant) & audit$Reported_significant
  if (anyNA(audit$Coefficient_converged[selected]) || any(!audit$Coefficient_converged[selected]))
    stop("A reported DESeq2 hit needs convergence review; inspect ", output_file)
  invisible(audit)
}

run_reporting_corrections <- function(work_dir, deseq_run = NULL, matched_run = NULL) {
  if (is.null(deseq_run)) deseq_run <- locate_saved_run(work_dir,
    "corrected_20261006_164720_527477", c("DESeq2_all_three_comparisons_full_results.csv",
      "DESeq2_SSL_HP_H_COL_dds.rds", "DESeq2_H_ADJ_SSL_HP_H_COL_dds.rds", "DESeq2_SSL_HP_H_ADJ_SSL_HP_dds.rds"))
  if (is.null(matched_run)) matched_run <- locate_saved_run(work_dir,
    "run_20261006_203539_746132", c("05_DESEQ2_GENUS/DESeq2_EXACT_MATCHED_genus_ALL_results.csv",
      "05_DESEQ2_GENUS/DESeq2_EXACT_MATCHED_genus_dds.rds", "07_ANCOMBC2_GENUS/ANCOMBC2_tidy_all_model_based_results.csv"))
  out_dir <- file.path(work_dir, "Reporting_corrected", paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
  if (dir.exists(out_dir)) stop("Output directory exists: ", out_dir)
  dir.create(out_dir, recursive = TRUE)
  input_files <- c(file.path(deseq_run, "DESeq2_all_three_comparisons_full_results.csv"),
    file.path(matched_run, "05_DESEQ2_GENUS/DESeq2_EXACT_MATCHED_genus_ALL_results.csv"),
    file.path(matched_run, "07_ANCOMBC2_GENUS/ANCOMBC2_tidy_all_model_based_results.csv"))
  stopifnot(all(file.exists(input_files)))
  write.csv(data.frame(Path = input_files, MD5 = unname(tools::md5sum(input_files))),
    file.path(out_dir, "input_files.csv"), row.names = FALSE)
  main <- read.csv(input_files[1L], check.names = FALSE)
  main_dir <- file.path(out_dir, "Main_DESeq2"); dir.create(main_dir)
  specifications <- list(c("SSL_HP", "H_COL"), c("H_ADJ_SSL_HP", "H_COL"), c("SSL_HP", "H_ADJ_SSL_HP"))
  plots <- list()
  reporting_changes <- list()
  for (g in specifications) {
    selected <- main$Comparison == paste(g, collapse = " - ")
    one <- main[selected, , drop = FALSE]
    stopifnot(nrow(one) > 0L)
    tag <- paste(g, collapse = "_")
    audit_saved_deseq(file.path(deseq_run, paste0("DESeq2_", tag, "_dds.rds")), one,
      file.path(main_dir, paste0(tag, "_convergence.csv")))
    n_missing <- sum(!is.finite(one$FDR))
    reporting_changes[[tag]] <- data.frame(Analysis = tag, Unavailable_FDR = n_missing,
      Action = if (n_missing) "Replace volcano panel" else "Keep existing volcano panel")
    if (n_missing == 0L) next
    export_env <- new.env(parent = environment())
    export_env$out_dir <- main_dir
    export_env$save_ggplot <- save_ggplot
    environment(export_env$save_ggplot) <- export_env
    export_env$replot_main_volcano <- replot_main_volcano
    environment(export_env$replot_main_volcano) <- export_env
    plots[[tag]] <- export_env$replot_main_volcano(one, g[1], g[2], tag, main_dir)
    write.csv(one[!is.finite(one$FDR), , drop = FALSE],
      file.path(main_dir, paste0(tag, "_unavailable_FDR_not_plotted.csv")), row.names = FALSE, na = "")
  }
  matched <- read.csv(input_files[2L], check.names = FALSE)
  ancom <- read.csv(input_files[3L], check.names = FALSE)
  stopifnot(all(c("Test_available", "Reporting_P", "Reporting_BH_FDR") %in% names(ancom)))
  ancom$Reporting_Holm_P <- ifelse(ancom$Test_available, ancom$q_value_Holm, NA_real_)
  ancom$Reporting_log2_difference <- ifelse(ancom$Test_available, ancom$log2_abundance_difference, NA_real_)
  matched_dir <- file.path(out_dir, "Matched"); dir.create(matched_dir)
  audit_saved_deseq(file.path(matched_run, "05_DESEQ2_GENUS/DESeq2_EXACT_MATCHED_genus_dds.rds"),
    matched, file.path(matched_dir, "DESeq2_convergence.csv"))
  if (any(!is.finite(matched$FDR))) {
    make_aligned_volcano(matched, "Genus", matched_dir, "DESeq2_EXACT_MATCHED_genus_volcano_WIDE_FULL_LABELS")
  } else {
    writeLines("All matched adjusted P values are available. Keep the existing matched DESeq2 volcano panel.",
      file.path(matched_dir, "KEEP_EXISTING_VOLCANO.txt"))
  }
  reporting_changes[["matched"]] <- data.frame(Analysis = "Matched genus", Unavailable_FDR = sum(!is.finite(matched$FDR)),
    Action = "Replace cross-method table; replace volcano only if unavailable FDRs exist")
  write.csv(matched[!is.finite(matched$FDR), , drop = FALSE],
    file.path(matched_dir, "unavailable_FDR_not_plotted.csv"), row.names = FALSE, na = "")
  export_env <- new.env(parent = environment())
  export_env$subdirs <- list(combined = matched_dir)
  export_env$make_cross_method_table <- make_cross_method_table
  environment(export_env$make_cross_method_table) <- export_env
  cross <- export_env$make_cross_method_table(list(result = matched), list(tidy = ancom), "Genus")
  stopifnot(all(is.na(cross$ANCOMBC2_P[!is.na(cross$ANCOMBC2_test_available) & !cross$ANCOMBC2_test_available])))
  write.csv(bind_rows(reporting_changes), file.path(out_dir, "report_changes.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  writeLines(c("Figures and comparison tables regenerated from saved results.",
    "No DESeq2, ANCOM-BC2, CLR or alpha-diversity model was refitted.",
    "Missing adjusted P values remain unavailable and are omitted from volcano plots.",
    "Original P values, adjusted P values and reporting thresholds are unchanged."),
    file.path(out_dir, "README.txt"))
  cat("\nCorrected reporting files saved to:\n", out_dir, "\n")
  invisible(out_dir)
}
