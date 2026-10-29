# Supplementary Figure SF5.6C: cross-sectional gastric-signature enrichment.
# Replots three archived aggregate GSEA rows without refitting any model.
# The first comparison is labelled Early (P1/P2); its archived key remains P1.
# Every run creates a new output directory.

required_packages <- c("ggplot2", "dplyr", "stringr", "scales", "readr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install these missing packages, then rerun: ",
    paste(missing_packages, collapse = ", "),
    "\nFor example: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "), "))"
  )
}
suppressPackageStartupMessages(library(dplyr))

OUTPUT_PARENT <- Sys.getenv("OSLC_FIGURE_OUTPUT_PARENT", getwd())
if (!dir.exists(OUTPUT_PARENT)) {
  stop("OSLC_FIGURE_OUTPUT_PARENT must name an existing output directory.",
       call. = FALSE)
}
output_base <- file.path(
  OUTPUT_PARENT,
  paste0("SF_5_6C_label_correction_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
OUTPUT_DIR <- output_base
output_counter <- 1L
while (file.exists(OUTPUT_DIR)) {
  OUTPUT_DIR <- paste0(output_base, "_", output_counter)
  output_counter <- output_counter + 1L
}

# Exact archived results; numerical values and significance stars are unchanged.
archived_csv_text <- "file,file_name,line_id,condition,numerator_passage,reference_passage,comparison,contrast_definition_assumed,positive_NES_interpretation,pathway,n_ranked_genes_used,intended_signature_size,represented_signature_size,represented_signature_genes,missing_signature_genes,ES,NES,pval,fgsea_padj_within_file,log2err,leading_edge,status,q_cross_sectional_BH,stars,star_basis\nCGG028_First_Passage_SSL_HP_DE_results.txt,CGG028_First_Passage_SSL_HP_DE_results.txt,CGG028,SSL_HP_vs_H_ADJ_SSL_HP,P1,H_ADJ_SSL_HP,P1,SSL_HP - H_ADJ_SSL_HP at P1,enriched in SSL_HP relative to H_ADJ_SSL_HP,Serrated Lesion Gastric Metaplasia Gene Signature,14685,24,20,AQP5;ANXA10;MUC5AC;VSIG1;TCN1;MUC6;PSAPL1;SFTA2;MUCL3;KLK11;DUOXA1;SULT1C2;AADAC;FOXQ1;NR0B2;C13ORF46;FOLR1;LYPD6B;TFR2;SPTSSB,CLDN18;GATA4;IRX3;IGFALS,0.7415008012978949,2.202753477154512,1.278866261666356e-05,1.278866261666356e-05,0.1880573095114043,AQP5;MUC6;ANXA10;MUC5AC;SULT1C2;SPTSSB;NR0B2;FOLR1;C13ORF46;KLK11;SFTA2;MUCL3;FOXQ1;AADAC;VSIG1,OK,3.836598784999068e-05,****,\"Benjamini-Hochberg q-value across the three cross-sectional comparisons P1, P4 and P8\"\nCGG028_Fourth_Passage_SSL_HP_DE_results.txt,CGG028_Fourth_Passage_SSL_HP_DE_results.txt,CGG028,SSL_HP_vs_H_ADJ_SSL_HP,P4,H_ADJ_SSL_HP,P4,SSL_HP - H_ADJ_SSL_HP at P4,enriched in SSL_HP relative to H_ADJ_SSL_HP,Serrated Lesion Gastric Metaplasia Gene Signature,14685,24,20,AQP5;ANXA10;MUC5AC;VSIG1;TCN1;MUC6;PSAPL1;SFTA2;MUCL3;KLK11;DUOXA1;SULT1C2;AADAC;FOXQ1;NR0B2;C13ORF46;FOLR1;LYPD6B;TFR2;SPTSSB,CLDN18;GATA4;IRX3;IGFALS,0.5981174013568318,1.8659238562448208,0.0024001560838551,0.0024001560838551,0.1368317913149987,AQP5;SPTSSB;FOXQ1;ANXA10;FOLR1;AADAC;MUC6;KLK11,OK,0.0036002341257826,**,\"Benjamini-Hochberg q-value across the three cross-sectional comparisons P1, P4 and P8\"\nCGG028_Eighth_Passage_SSL_HP_DE_results.txt,CGG028_Eighth_Passage_SSL_HP_DE_results.txt,CGG028,SSL_HP_vs_H_ADJ_SSL_HP,P8,H_ADJ_SSL_HP,P8,SSL_HP - H_ADJ_SSL_HP at P8,enriched in SSL_HP relative to H_ADJ_SSL_HP,Serrated Lesion Gastric Metaplasia Gene Signature,14685,24,20,AQP5;ANXA10;MUC5AC;VSIG1;TCN1;MUC6;PSAPL1;SFTA2;MUCL3;KLK11;DUOXA1;SULT1C2;AADAC;FOXQ1;NR0B2;C13ORF46;FOLR1;LYPD6B;TFR2;SPTSSB,CLDN18;GATA4;IRX3;IGFALS,0.5642812517553165,1.7450961904720186,0.0080127390367987,0.0080127390367987,0.1206742970889914,ANXA10;AQP5;SPTSSB;TCN1;FOLR1;C13ORF46;SULT1C2;MUC6;AADAC;NR0B2;LYPD6B;VSIG1,OK,0.0080127390367987,**,\"Benjamini-Hochberg q-value across the three cross-sectional comparisons P1, P4 and P8\"\n"
cross_sectional_gsea_results <- utils::read.csv(
  text = archived_csv_text,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
CROSS_SECTIONAL_PASSAGE_ORDER <- c("P1", "P4", "P8")
METAPLASIA_SIGNATURE_TITLE <-
  "Serrated Lesion Gastric Metaplasia Gene Signature"
gastric_metaplasia_24 <- c(
  "AQP5", "ANXA10", "MUC5AC", "VSIG1", "TCN1", "MUC6",
  "PSAPL1", "SFTA2", "MUCL3", "KLK11", "DUOXA1", "SULT1C2",
  "AADAC", "FOXQ1", "NR0B2", "C13ORF46", "FOLR1", "LYPD6B",
  "TFR2", "SPTSSB", "CLDN18", "GATA4", "IRX3", "IGFALS"
)
stopifnot(
  identical(cross_sectional_gsea_results$comparison, CROSS_SECTIONAL_PASSAGE_ORDER),
  all(cross_sectional_gsea_results$represented_signature_size == 20L),
  all(cross_sectional_gsea_results$n_ranked_genes_used == 14685L),
  identical(cross_sectional_gsea_results$stars, c("****", "**", "**"))
)
cross_sectional_represented_size <- 20L
cross_sectional_missing_genes <- sort(unique(unlist(strsplit(
  cross_sectional_gsea_results$missing_signature_genes, ";", fixed = TRUE
))))
primary_pdf_device <- if (capabilities("cairo")) {
  grDevices::cairo_pdf
} else {
  grDevices::pdf
}

if (!dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)) {
  stop("Could not create output folder: ", OUTPUT_DIR)
}

# Plot labels and saved results
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

# Wrap captions to fit the full-width figure.
cross_sectional_subtitle_text <- stringr::str_wrap(
  paste0(
    "Preranked GSEA: SSL_HP vs H_ADJ_SSL_HP organoids at early passage, P4 and P8 (",
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
  "Stars: BH q-values across early passage, P4 and P8 ",
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
    labels = c(P1 = "Early (P1/P2)", P4 = "P4", P8 = "P8"),
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

# Full-width export with space for the caption.
cross_sectional_width_mm <- 183
cross_sectional_height_mm <- 92
cross_sectional_stem <- file.path(
  OUTPUT_DIR,
  "SF_5_6C_early_passage_corrected"
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
    "SF_5_6C_early_passage_corrected_data.csv"
  ),
  na = ""
)


# Also record the corrected display label alongside the archived identifiers.
display_data <- cross_sectional_gsea_results
display_data$display_passage <- c("Early (P1/P2)", "P4", "P8")
readr::write_csv(
  display_data,
  file.path(OUTPUT_DIR, "SF_5_6C_saved_results_with_display_labels.csv")
)
writeLines(
  capture.output(utils::sessionInfo()),
  file.path(OUTPUT_DIR, "sessionInfo.txt")
)
print(cross_sectional_nes_heatmap)
message("Corrected SF_5.6C saved to:\n", normalizePath(OUTPUT_DIR, winslash = "/"))
