#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

if (capabilities("cairo")) {
  options(bitmapType = "cairo")
}

config <- new.env(parent = baseenv())
if (file.exists("config/config.R")) sys.source("config/config.R", envir = config)
work_dir <- if (exists("PATH_WORK", config, inherits = FALSE)) config$PATH_WORK else "outputs"
if (!dir.exists(work_dir)) stop("Input directory does not exist: ", work_dir)

if (!exists("ANCOMBC2_REVIEWED_RUN", config, inherits = FALSE)) {
  stop("Set ANCOMBC2_REVIEWED_RUN in config/config.R to the reviewed saved run directory.")
}
run_dir <- config$ANCOMBC2_REVIEWED_RUN

figure_dir <- file.path(run_dir, "FINAL_FIGURES")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

output_stem <- file.path(
  figure_dir,
  "Supplementary_Figure_S4_5_3_ANCOMBC2_PANEL_A_ONLY_FINAL"
)

font_family <- "Arial"
base_font_size <- 10
figure_width_mm <- 178
figure_height_mm <- 100
export_dpi <- 600
show_panel_letter <- FALSE

plot_colours <- c(
  "Higher in SSL_HP" = "#FF0000",
  "Higher in H_COL" = "#1616E8"
)

if (!dir.exists(run_dir)) {
  stop("Run directory not found: ", run_dir, call. = FALSE)
}

input_candidates <- c(
  file.path(
    figure_dir,
    "Supplementary_Figure_S4_5_3_panel_A_source_data.csv"
  ),
  file.path(
    run_dir,
    "SSL_HP_H_COL",
    "ANCOMBC2_SSL_HP_H_COL_model_based_BH_FDR_supported.csv"
  )
)

input_files_found <- input_candidates[file.exists(input_candidates)]

if (length(input_files_found) == 0L) {
  stop(
    "Panel A input file not found. Checked:\n",
    paste(input_candidates, collapse = "\n"),
    call. = FALSE
  )
}

input_file <- input_files_found[1]
message("Reading: ", input_file)

plot_data <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "taxon",
  "Comparison",
  "LFC_log2",
  "CI95_low_log2",
  "CI95_high_log2",
  "BH_FDR",
  "Holm_adjusted_P_sensitivity",
  "model_based_BH_and_pseudocount_stable"
)

missing_columns <- setdiff(required_columns, names(plot_data))

if (length(missing_columns) > 0L) {
  stop(
    "Input table is missing: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

to_logical <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(!is.na(x) & x != 0)
  trimws(tolower(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

format_fdr <- function(x) {
  ifelse(
    x < 0.001,
    formatC(x, format = "e", digits = 2),
    formatC(x, format = "f", digits = 5)
  )
}

plot_data <- plot_data %>%
  mutate(
    model_based_BH_and_pseudocount_stable =
      to_logical(model_based_BH_and_pseudocount_stable),
    LFC_log2 = as.numeric(LFC_log2),
    CI95_low_log2 = as.numeric(CI95_low_log2),
    CI95_high_log2 = as.numeric(CI95_high_log2),
    BH_FDR = as.numeric(BH_FDR),
    Holm_adjusted_P_sensitivity =
      as.numeric(Holm_adjusted_P_sensitivity)
  ) %>%
  filter(
    Comparison == "SSL_HP - H_COL",
    model_based_BH_and_pseudocount_stable,
    !is.na(BH_FDR),
    BH_FDR < 0.05
  )

if (nrow(plot_data) == 0L) {
  stop(
    paste0(
      "No SSL_HP-versus-H_COL genera passed BH FDR < 0.05 together ",
      "with the pseudocount-sensitivity criterion."
    ),
    call. = FALSE
  )
}

if (
  any(!is.finite(plot_data$LFC_log2)) ||
  any(!is.finite(plot_data$CI95_low_log2)) ||
  any(!is.finite(plot_data$CI95_high_log2)) ||
  any(!is.finite(plot_data$BH_FDR))
) {
  stop("The plotted result table contains non-finite values.", call. = FALSE)
}

if (
  any(plot_data$CI95_low_log2 > plot_data$LFC_log2) ||
  any(plot_data$CI95_high_log2 < plot_data$LFC_log2)
) {
  stop("At least one confidence interval is inconsistent.", call. = FALSE)
}

plot_data <- plot_data %>%
  mutate(
    Result = ifelse(
      LFC_log2 > 0,
      "Higher in SSL_HP",
      "Higher in H_COL"
    ),
    Holm_pass =
      !is.na(Holm_adjusted_P_sensitivity) &
      is.finite(Holm_adjusted_P_sensitivity) &
      Holm_adjusted_P_sensitivity < 0.05,
    FDR_label = paste0(
      "BH FDR = ",
      format_fdr(BH_FDR),
      ifelse(Holm_pass, " \u2020", "")
    )
  ) %>%
  arrange(LFC_log2) %>%
  mutate(

    taxon = factor(taxon, levels = unique(taxon))
  )

write.csv(
  plot_data,
  paste0(output_stem, "_source_data.csv"),
  row.names = FALSE,
  na = ""
)

minimum_ci <- min(plot_data$CI95_low_log2)
maximum_ci <- max(plot_data$CI95_high_log2)
effect_span <- maximum_ci - minimum_ci

if (!is.finite(effect_span) || effect_span <= 0) {
  effect_span <- 1
}

left_limit <- minimum_ci - 0.06 * effect_span
separator_x <- maximum_ci + 0.08 * effect_span
fdr_text_x <- maximum_ci + 0.15 * effect_span
right_limit <- maximum_ci + 0.72 * effect_span

x_breaks <- pretty(c(left_limit, maximum_ci), n = 6)
x_breaks <- x_breaks[x_breaks >= left_limit & x_breaks <= maximum_ci]

plot_data <- plot_data %>%
  mutate(FDR_text_x = fdr_text_x)

panel_tag <- if (show_panel_letter) "A" else NULL

final_plot <- ggplot(
  plot_data,
  aes(
    x = LFC_log2,
    y = taxon,
    colour = Result
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey55"
  ) +
  geom_vline(
    xintercept = separator_x,
    linewidth = 0.35,
    colour = "grey80"
  ) +
  geom_segment(
    aes(
      x = CI95_low_log2,
      xend = CI95_high_log2,
      yend = taxon
    ),
    linewidth = 0.75,
    colour = "black",
    lineend = "butt"
  ) +
  geom_point(size = 3.0) +
  geom_text(
    aes(
      x = FDR_text_x,
      label = FDR_label
    ),
    hjust = 0,
    colour = "black",
    family = font_family,
    size = 3.45,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = plot_colours,
    breaks = names(plot_colours),
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(left_limit, right_limit),
    breaks = x_breaks,
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    tag = panel_tag,
    x = paste0(
      "ANCOM-BC2 bias-corrected log2 abundance difference\n",
      "(positive values indicate higher abundance in SSL_HP)"
    ),
    y = NULL,
    colour = NULL,
    caption = paste0(
      "\u2020 Holm-adjusted P < 0.05.\n",
      "Horizontal intervals are unadjusted 95% Wald confidence intervals."
    )
  ) +
  theme_classic(
    base_family = font_family,
    base_size = base_font_size
  ) +
  theme(
    axis.text = element_text(
      family = font_family,
      size = base_font_size,
      colour = "black"
    ),
    axis.text.y = element_text(
      family = font_family,
      size = base_font_size + 0.5,
      face = "italic",
      margin = margin(r = 2, unit = "mm")
    ),
    axis.title.x = element_text(
      family = font_family,
      size = base_font_size,
      margin = margin(t = 4, unit = "mm")
    ),
    axis.line = element_line(
      linewidth = 0.45,
      colour = "black"
    ),
    axis.ticks = element_line(
      linewidth = 0.45,
      colour = "black"
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.text = element_text(
      family = font_family,
      size = base_font_size
    ),
    legend.key.width = grid::unit(5, "mm"),
    legend.spacing.x = grid::unit(2, "mm"),
    plot.tag = element_text(
      family = font_family,
      face = "bold",
      size = base_font_size + 4,
      hjust = 0
    ),
    plot.caption = element_text(
      family = font_family,
      size = base_font_size - 1,
      hjust = 0,
      margin = margin(t = 2.5, unit = "mm")
    ),
    plot.margin = margin(
      t = 5,
      r = 7,
      b = 5,
      l = 5,
      unit = "mm"
    )
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 3.0)
    )
  )

pdf_device <- if (capabilities("cairo")) {
  grDevices::cairo_pdf
} else {
  grDevices::pdf
}

ggsave(
  filename = paste0(output_stem, ".pdf"),
  plot = final_plot,
  device = pdf_device,
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = paste0(output_stem, ".png"),
  plot = final_plot,
  device = "png",
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  dpi = export_dpi,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = paste0(output_stem, ".tiff"),
  plot = final_plot,
  device = "tiff",
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  dpi = export_dpi,
  compression = "lzw",
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = paste0(output_stem, ".jpeg"),
  plot = final_plot,
  device = "jpeg",
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  dpi = export_dpi,
  quality = 100,
  bg = "white",
  limitsize = FALSE
)

output_files <- paste0(
  output_stem,
  c(".pdf", ".png", ".tiff", ".jpeg")
)

missing_outputs <- output_files[!file.exists(output_files)]

if (length(missing_outputs) > 0L) {
  stop(
    "The following files were not created: ",
    paste(missing_outputs, collapse = ", "),
    call. = FALSE
  )
}

capture.output(
  sessionInfo(),
  file = paste0(output_stem, "_sessionInfo.txt")
)

message("")
message("FINAL PANEL-A-ONLY FIGURE CREATED SUCCESSFULLY")
message("================================================")
message("Genera plotted: ", nrow(plot_data))
message("Input: ", input_file)
message("")
message("Created files:")
for (output_file in output_files) {
  message(output_file)
}
message("")
message("Panel B is not included, and patchwork is not used.")
