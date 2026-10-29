#!/usr/bin/env Rscript
# Build pathway figures from existing GSEA results and differential-expression tables.
# Run from Tissue_Bulk_RNAseq; PATH_GSEA_RESULTS selects the saved GSEA run.
# Gene-selection thresholds, figure dimensions and palettes are specified below.


source("config/config.R")

input_dir <- PATH_DE_TABLES

paper_ready_dir <- PATH_GSEA_RESULTS

fgsea_by_comparison_dir <- file.path(
  paper_ready_dir,
  "tables",
  "fgsea_by_comparison"
)

ranked_dir <- file.path(
  paper_ready_dir,
  "ranked_gene_tables"
)

de_dir <- input_dir

out_dir <- OUTPUT_PATHWAY_FIGURES

ensure_new_output(out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gsea_filter_mode <- "LOGCPM_GE_2"

font_family <- "Arial"

font_base_size <- 6
axis_title_size <- 6
axis_text_size <- 5.0
legend_title_size <- 5.0
legend_text_size <- 5.0
strip_text_size <- 6
plot_title_size <- 7
panel_label_size <- 8

thin_line <- 0.25
regular_line <- 0.40

panelB_reference_comparison <- "SSL_HP_v_H_COL"
panelB_logcpm_threshold <- 2
panelB_logfc_threshold <- 1
max_panelB_genes_per_group <- 12
heatmap_colour_saturation_logfc <- 2

panelA_label_mode <- "stars"
panelA_star_colour <- "white"

panelA_comparisons <- c(
  "SSL_HP_v_H_COL",
  "SSL_HP_v_H_ADJ_SSL_HP",
  "SSL_HP_v_CA"
)

panelB_comparisons <- panelA_comparisons

panelC_comparisons <- c(
  "ES_SSL_HP",
  "ES_H_ADJ_SSL_HP"
)

panelD_comparisons <- panelC_comparisons

es_label_ssl_hp <- "Escherichia-Shigella + vs - SSL_HP"
es_label_h_adj_ssl_hp <- "Escherichia-Shigella + vs - H_ADJ_SSL_HP"

comparison_labels <- c(
  SSL_HP_v_H_COL = "SSL_HP vs H_COL",
  SSL_HP_v_H_ADJ_SSL_HP = "SSL_HP vs H_ADJ_SSL_HP",
  SSL_HP_v_CA = "SSL_HP vs CA",
  H_ADJ_SSL_HP_v_H_COL = "H_ADJ_SSL_HP vs H_COL",
  CA_v_H_COL = "CA vs H_COL",
  ES_SSL_HP = es_label_ssl_hp,
  ES_H_ADJ_SSL_HP = es_label_h_adj_ssl_hp,
  ES_SSL_HP_and_H_ADJ_SSL_HP = "Escherichia-Shigella + vs - all",
  Alt_Infected = "Alternative infected",
  Initial_BRAF_SSL_HP = "Initial BRAF SSL_HP",
  V600E_v_WT_in_SSL_HP = "V600E vs WT in SSL_HP"
)

heatmap_axis_labels <- c(
  SSL_HP_v_H_COL = "SSL_HP vs H_COL",
  SSL_HP_v_H_ADJ_SSL_HP = "SSL_HP vs H_ADJ_SSL_HP",
  SSL_HP_v_CA = "SSL_HP vs CA",
  ES_SSL_HP = "ES+ vs ES- SSL_HP",
  ES_H_ADJ_SSL_HP = "ES+ vs ES- H_ADJ_SSL_HP"
)

panelC_short_labels <- c(
  ES_SSL_HP = "SSL_HP",
  ES_H_ADJ_SSL_HP = "H_ADJ"
)

sample_type_cols <- c(
  SSL_HP = "#FF0000",
  H_ADJ_SSL_HP = "#E9C61D",
  STOMACH = "#F155D0",
  H_COL = "#1313EC",
  CA = "#825CA6"
)

heatmap_positive_colour <- "#B2182B"
heatmap_negative_colour <- "#2166AC"
heatmap_mid_colour <- "white"

heatmap_tile_border_colour <- "white"
heatmap_tile_border_linewidth <- thin_line
heatmap_tile_fraction <- 0.96

shared_heatmap_label_pt <- 5.0
shared_heatmap_x_label_angle <- 55

ggplot_text_pt_per_mm <- 2.845276

geom_text_size_from_pt <- function(pt) {
  pt / ggplot_text_pt_per_mm
}

manual_heatmap_x_label_angle <- shared_heatmap_x_label_angle
manual_heatmap_x_label_size <- geom_text_size_from_pt(shared_heatmap_label_pt)
manual_heatmap_gene_label_size <- geom_text_size_from_pt(shared_heatmap_label_pt)
manual_heatmap_module_title_size <- geom_text_size_from_pt(plot_title_size)

panelA_x_label_angle <- shared_heatmap_x_label_angle
panelA_axis_text_pt <- shared_heatmap_label_pt

panelC_cols <- c(
  ES_SSL_HP = sample_type_cols[["SSL_HP"]],
  ES_H_ADJ_SSL_HP = sample_type_cols[["H_ADJ_SSL_HP"]]
)

panelA_requested_pathways <- c(

  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_COMPLEMENT",
  "REACTOME_NEUTROPHIL_DEGRANULATION",
  "REACTOME_INTERLEUKIN_10_SIGNALING",
  "REACTOME_MYD88_INDEPENDENT_TLR4_CASCADE",
  "REACTOME_TOLL_LIKE_RECEPTOR_CASCADES",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_MEDICUS_PATHOGEN_SHIGELLA_IPGD_TO_ARNO_ARF_ACTB_G_SIGNALING_PATHWAY",
  "KEGG_PATHOGENIC_ESCHERICHIA_COLI_INFECTION"
)

de_file_map <- c(
  CA_v_H_COL = "CGG014_Maier_CA_v_H_COL_DE_results.txt",
  H_ADJ_SSL_HP_v_H_COL = "CGG014_Maier_H_ADJ_SSL_HP_v_H_COL_DE_results.txt",
  SSL_HP_v_CA = "CGG014_Maier_SSL_HP_v_CA_DE_results.txt",
  SSL_HP_v_H_ADJ_SSL_HP = "CGG014_Maier_SSL_HP_v_H_ADJ_SSL_HP_DE_results.txt",
  SSL_HP_v_H_COL = "CGG014_Maier_SSL_HP_v_H_COL_DE_results.txt",
  V600E_v_WT_in_SSL_HP = "CGG014_Maier_V600E_v_WT_in_SSL_HP_DE_results.txt",
  Alt_Infected = "CGG028_Alt_Infected_DE_results.txt",
  Initial_BRAF_SSL_HP = "CGG028_Initial_BRAF_SSL_HP_DE_results.txt",
  ES_H_ADJ_SSL_HP = "CGG014_Maier_ES_H_ADJ_SSL_HP_DE_results.txt",
  ES_SSL_HP_and_H_ADJ_SSL_HP = "CGG014_Maier_ES_SSL_HP_and_H_ADJ_SSL_HP_DE_results.txt",
  ES_SSL_HP = "CGG014_Maier_ES_SSL_HP_DE_results.txt"
)

cran_pkgs <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "purrr",
  "ggplot2",
  "patchwork",
  "scales",
  "msigdbr",
  "tibble"
)

missing_cran <- cran_pkgs[
  !vapply(cran_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_cran) > 0) {
  stop("Install required packages before running: ", paste(missing_cran, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(msigdbr)
  library(tibble)
  library(grid)
})

vertical_heatmap_guide <- function(barheight_mm = 18, barwidth_mm = 2.0) {
  guide_colourbar(
    title.position = "right",
    title.hjust = 0.5,
    title.vjust = 0.5,
    title.theme = element_text(
      family = font_family,
      size = legend_title_size,
      angle = 90,
      hjust = 0.5,
      vjust = 0.5,
      color = "black",
      margin = margin(l = 2, r = 0, unit = "mm")
    ),
    label.position = "right",
    label.hjust = 0.5,
    label.vjust = 0.5,
    label.theme = element_text(
      family = font_family,
      size = legend_text_size,
      color = "black"
    ),
    barwidth = grid::unit(barwidth_mm, "mm"),
    barheight = grid::unit(barheight_mm, "mm"),
    ticks = TRUE
  )
}

if (!dir.exists(paper_ready_dir)) {
  stop("paper_ready_dir does not exist: ", paper_ready_dir)
}
if (!dir.exists(fgsea_by_comparison_dir)) {
  stop("Missing fgsea_by_comparison directory: ", fgsea_by_comparison_dir)
}
if (!dir.exists(ranked_dir)) {
  stop("Missing ranked_gene_tables directory: ", ranked_dir)
}
if (!dir.exists(de_dir)) {
  stop("Missing DE table directory: ", de_dir)
}

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    "NA",
    ifelse(
      p < 1e-4,
      formatC(p, format = "e", digits = 1),
      as.character(signif(p, 2))
    )
  )
}

sig_stars <- function(padj) {
  dplyr::case_when(
    is.na(padj) ~ "",
    padj < 0.001 ~ "***",
    padj < 0.01 ~ "**",
    padj < 0.05 ~ "*",
    TRUE ~ ""
  )
}

safe_filename <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

label_comp <- function(x) {
  out <- comparison_labels[x]
  out[is.na(out)] <- x[is.na(out)]
  unname(out)
}

label_comp_heatmap <- function(x) {
  out <- heatmap_axis_labels[x]
  out[is.na(out)] <- label_comp(x[is.na(out)])
  unname(out)
}

rename_ci <- function(df, target, candidates) {
  if (target %in% names(df)) return(df)

  idx <- match(tolower(candidates), tolower(names(df)))
  idx <- idx[!is.na(idx)]

  if (length(idx) > 0) {
    names(df)[idx[1]] <- target
  }

  df
}

clean_pathway_name <- function(z) {
  toupper(z) %>%
    stringr::str_remove("^HALLMARK_") %>%
    stringr::str_remove("^KEGG_LEGACY_") %>%
    stringr::str_remove("^KEGG_MEDICUS_") %>%
    stringr::str_remove("^KEGG_") %>%
    stringr::str_remove("^REACTOME_") %>%
    stringr::str_remove("^GOMF_")
}

resolve_pathway_name <- function(target, available_pathways) {
  if (is.na(target) || target == "") return(NA_character_)

  available_pathways <- unique(as.character(available_pathways))
  available_pathways <- available_pathways[
    !is.na(available_pathways) & available_pathways != ""
  ]

  if (length(available_pathways) == 0) return(NA_character_)

  target_u <- toupper(as.character(target))
  available_u <- toupper(available_pathways)

  exact_hit <- available_pathways[available_u == target_u]
  if (length(exact_hit) > 0) return(exact_hit[1])

  target_clean <- clean_pathway_name(target_u)
  available_clean <- clean_pathway_name(available_u)

  clean_hit <- available_pathways[available_clean == target_clean]
  if (length(clean_hit) > 0) return(clean_hit[1])

  NA_character_
}

resolve_pathway_list <- function(targets, available_pathways) {
  tibble::tibble(
    requested_pathway = targets,
    resolved_pathway = vapply(
      targets,
      resolve_pathway_name,
      character(1),
      available_pathways = available_pathways
    )
  )
}

find_pathway_by_pattern <- function(named_pathways, patterns) {
  pathway_names <- names(named_pathways)
  pathway_names <- pathway_names[
    !is.na(pathway_names) & pathway_names != ""
  ]

  if (length(pathway_names) == 0) return(NA_character_)

  pathway_upper <- toupper(pathway_names)

  for (pat in patterns) {
    hits <- pathway_names[stringr::str_detect(pathway_upper, pat)]
    if (length(hits) > 0) return(hits[1])
  }

  NA_character_
}

get_required_pathway_info <- function(
    requested_pathway,
    named_pathways,
    loose_patterns = NULL,
    source_label = "MSigDB"
) {
  resolved <- resolve_pathway_name(requested_pathway, names(named_pathways))

  if (is.na(resolved) && !is.null(loose_patterns)) {
    resolved <- find_pathway_by_pattern(named_pathways, loose_patterns)
  }

  if (is.na(resolved) || !resolved %in% names(named_pathways)) {
    stop(
      "Required gene set was not found: ",
      requested_pathway,
      "\nAvailable examples: ",
      paste(utils::head(names(named_pathways), 20), collapse = "; ")
    )
  }

  genes <- sort(unique(toupper(trimws(named_pathways[[resolved]]))))
  genes <- genes[!is.na(genes) & genes != ""]

  if (length(genes) == 0) {
    stop("Required gene set contains no genes: ", requested_pathway)
  }

  list(
    requested_pathway = requested_pathway,
    resolved_pathway = resolved,
    source = source_label,
    genes = genes
  )
}

short_pathway_label <- function(x) {
  x_chr <- as.character(x)
  xu <- toupper(x_chr)

  fallback <- x_chr %>%
    stringr::str_remove("^HALLMARK_") %>%
    stringr::str_remove("^KEGG_LEGACY_") %>%
    stringr::str_remove("^KEGG_MEDICUS_") %>%
    stringr::str_remove("^KEGG_") %>%
    stringr::str_remove("^REACTOME_") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_lower() %>%
    stringr::str_to_sentence()

  dplyr::case_when(
    xu == "KEGG_MEDICUS_PATHOGEN_SHIGELLA_IPGD_TO_ARNO_ARF_ACTB_G_SIGNALING_PATHWAY" ~ "Shigella ARF/ACTB signaling",
    xu == "KEGG_PATHOGENIC_ESCHERICHIA_COLI_INFECTION" ~ "Pathogenic E. coli infection",
    xu == "HALLMARK_INTERFERON_ALPHA_RESPONSE" ~ "Interferon alpha response",
    xu == "HALLMARK_GLYCOLYSIS" ~ "Glycolysis",
    xu == "HALLMARK_OXIDATIVE_PHOSPHORYLATION" ~ "Oxidative phosphorylation",
    xu == "HALLMARK_COMPLEMENT" ~ "Complement",
    xu == "HALLMARK_UNFOLDED_PROTEIN_RESPONSE" ~ "Unfolded protein response",
    xu == "HALLMARK_WNT_BETA_CATENIN_SIGNALING" ~ "Wnt/β-catenin signaling",
    xu == "HALLMARK_INFLAMMATORY_RESPONSE" ~ "Inflammatory response",
    xu == "HALLMARK_TNFA_SIGNALING_VIA_NFKB" ~ "TNFA/NFKB signaling",
    xu == "REACTOME_NEUTROPHIL_DEGRANULATION" ~ "Neutrophil degranulation",
    xu == "REACTOME_INTERLEUKIN_10_SIGNALING" ~ "Interleukin-10 signaling",
    xu == "REACTOME_MYD88_INDEPENDENT_TLR4_CASCADE" ~ "MYD88-independent TLR4 cascade",
    xu == "REACTOME_TOLL_LIKE_RECEPTOR_CASCADES" ~ "Toll-like receptor cascades",
    xu == "REACTOME_ANTIMICROBIAL_PEPTIDES" ~ "Antimicrobial peptides",
    xu == "GOMF_CHEMOKINE_ACTIVITY" ~ "Chemokine activity",
    stringr::str_detect(xu, "SHIGELLA.*ARNO|SHIGELLA.*ACTB|SHIGELLA_IPGD") ~ "Shigella ARF/ACTB signaling",
    stringr::str_detect(xu, "PATHOGENIC_ESCHERICHIA_COLI_INFECTION") ~ "Pathogenic E. coli infection",
    TRUE ~ fallback
  )
}

nature_theme <- function(base_size = font_base_size) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      text = element_text(family = font_family, color = "black"),
      plot.title = element_text(
        family = font_family,
        face = "bold",
        hjust = 0,
        size = plot_title_size
      ),
      plot.subtitle = element_blank(),
      axis.title = element_text(
        family = font_family,
        size = axis_title_size
      ),
      axis.text = element_text(
        family = font_family,
        size = axis_text_size,
        color = "black"
      ),
      legend.title = element_text(
        family = font_family,
        size = legend_title_size
      ),
      legend.text = element_text(
        family = font_family,
        size = legend_text_size
      ),
      strip.background = element_rect(
        fill = "grey92",
        color = "grey70",
        linewidth = thin_line
      ),
      strip.text = element_text(
        family = font_family,
        face = "bold",
        size = strip_text_size
      ),
      panel.border = element_rect(
        fill = NA,
        color = "black",
        linewidth = thin_line
      ),
      axis.line = element_line(linewidth = thin_line),
      axis.ticks = element_line(linewidth = thin_line)
    )
}

ggplot2::theme_set(nature_theme(font_base_size))

save_plot_set <- function(plot, stem, width_mm, height_mm, dpi = 600, jpg_quality = 100) {
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

  ggsave(
    file.path(out_dir, paste0(stem, ".pdf")),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = pdf_device,
    bg = "white"
  )

  ggsave(
    file.path(out_dir, paste0(stem, ".png")),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = dpi,
    bg = "white"
  )

  ggsave(
    file.path(out_dir, paste0(stem, ".tif")),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = dpi,
    device = grDevices::tiff,
    compression = "lzw",
    bg = "white"
  )

  ggsave(
    file.path(out_dir, paste0(stem, ".jpg")),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = dpi,
    device = grDevices::jpeg,
    quality = jpg_quality,
    bg = "white"
  )

  ggsave(
    file.path(out_dir, paste0(stem, ".jpeg")),
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = dpi,
    device = grDevices::jpeg,
    quality = jpg_quality,
    bg = "white"
  )
}

make_empty_plot <- function(title_text, message_text = "No genes passed the selection rule") {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = message_text,
      size = 2.0,
      family = font_family
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_void(base_family = font_family) +
    theme(
      plot.title = element_text(
        family = font_family,
        face = "bold",
        hjust = 0,
        size = plot_title_size
      ),
      plot.subtitle = element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
}

parse_fgsea_filename <- function(path, mode = gsea_filter_mode) {
  name <- basename(path)
  stem <- stringr::str_remove(name, "_FGSEA\\.csv$")

  marker <- paste0("_", mode, "_")
  idx <- regexpr(marker, stem, fixed = TRUE)

  if (idx < 1) {
    return(tibble(
      file = path,
      file_name = name,
      comparison = NA_character_,
      filter_mode = NA_character_,
      gene_set_family = NA_character_
    ))
  }

  comparison <- substr(stem, 1, idx - 1)
  gene_set_family <- substr(stem, idx + nchar(marker), nchar(stem))

  tibble(
    file = path,
    file_name = name,
    comparison = comparison,
    filter_mode = mode,
    gene_set_family = gene_set_family
  )
}

standardize_fgsea_cols <- function(dat, comparison, family, mode, source_file) {
  dat <- rename_ci(dat, "pathway", c("pathway", "Pathway", "gs_name", "geneset", "gene_set", "term"))
  dat <- rename_ci(dat, "NES", c("NES", "nes", "normalized_enrichment_score"))
  dat <- rename_ci(dat, "ES", c("ES", "es", "enrichment_score"))
  dat <- rename_ci(dat, "pval", c("pval", "p.value", "P.Value", "PValue", "pvalue", "p"))
  dat <- rename_ci(dat, "padj", c("padj", "FDR", "fdr", "qvalue", "q.value", "FDR.q.val"))
  dat <- rename_ci(dat, "size", c("size", "setSize", "geneset_size"))
  dat <- rename_ci(dat, "log2err", c("log2err", "log2_error"))
  dat <- rename_ci(dat, "leadingEdge", c("leadingEdge", "leading_edge", "leading.edge"))

  required <- c("pathway", "NES", "pval", "padj")
  missing <- setdiff(required, names(dat))

  if (length(missing) > 0) {
    stop(
      "FGSEA file is missing required columns: ",
      paste(missing, collapse = ", "),
      "\nFile: ",
      source_file
    )
  }

  if (!"ES" %in% names(dat)) dat$ES <- NA_real_
  if (!"size" %in% names(dat)) dat$size <- NA_real_
  if (!"log2err" %in% names(dat)) dat$log2err <- NA_real_
  if (!"leadingEdge" %in% names(dat)) dat$leadingEdge <- NA_character_

  dat %>%
    dplyr::mutate(
      pathway = as.character(pathway),
      NES = suppressWarnings(as.numeric(NES)),
      ES = suppressWarnings(as.numeric(ES)),
      pval = suppressWarnings(as.numeric(pval)),
      padj = suppressWarnings(as.numeric(padj)),
      size = suppressWarnings(as.numeric(size)),
      comparison = comparison,
      filter_mode = mode,
      gene_set_family = family,
      source_file = basename(source_file)
    ) %>%
    dplyr::filter(!is.na(pathway), pathway != "")
}

list_fgsea_files <- function(
    mode = gsea_filter_mode,
    families = c("HALLMARK", "KEGG", "REACTOME")
) {
  files <- list.files(
    fgsea_by_comparison_dir,
    pattern = paste0("_", mode, "_.*_FGSEA\\.csv$"),
    full.names = TRUE,
    recursive = FALSE
  )

  if (length(files) == 0) {
    stop(
      "No ",
      mode,
      " FGSEA files found in: ",
      fgsea_by_comparison_dir
    )
  }

  meta <- purrr::map_dfr(files, parse_fgsea_filename, mode = mode) %>%
    dplyr::filter(!is.na(comparison), !is.na(gene_set_family))

  if (!is.null(families)) {
    meta <- meta %>%
      dplyr::filter(toupper(gene_set_family) %in% toupper(families))
  }

  meta
}

read_fgsea_files <- function(
    mode = gsea_filter_mode,
    families = c("HALLMARK", "KEGG", "REACTOME")
) {
  meta <- list_fgsea_files(mode = mode, families = families)

  purrr::pmap_dfr(
    list(
      meta$file,
      meta$comparison,
      meta$gene_set_family,
      meta$filter_mode
    ),
    function(file, comparison, gene_set_family, filter_mode) {
      dat <- readr::read_csv(file, show_col_types = FALSE)

      standardize_fgsea_cols(
        dat = dat,
        comparison = comparison,
        family = gene_set_family,
        mode = filter_mode,
        source_file = file
      )
    }
  )
}

message("Reading LOGCPM_GE_2 FGSEA files from fgsea_by_comparison...")

combined_fgsea <- read_fgsea_files(
  mode = gsea_filter_mode,
  families = c("HALLMARK", "KEGG", "REACTOME")
)

readr::write_csv(
  combined_fgsea,
  file.path(out_dir, "COMBINED_LOGCPM_GE_2_FGSEA_FROM_PER_COMPARISON_FILES.csv")
)

fgsea_hallmark <- combined_fgsea %>%
  dplyr::filter(toupper(gene_set_family) == "HALLMARK")

fgsea_kegg <- combined_fgsea %>%
  dplyr::filter(toupper(gene_set_family) == "KEGG")

fgsea_reactome <- combined_fgsea %>%
  dplyr::filter(toupper(gene_set_family) == "REACTOME")

if (nrow(fgsea_hallmark) == 0) stop("No LOGCPM_GE_2 HALLMARK FGSEA rows were found.")
if (nrow(fgsea_kegg) == 0) stop("No LOGCPM_GE_2 KEGG FGSEA rows were found.")
if (nrow(fgsea_reactome) == 0) stop("No LOGCPM_GE_2 REACTOME FGSEA rows were found.")

message("Loading MSigDB Hallmark, Reactome and GO gene sets via msigdbr...")

msig_h <- tryCatch(
  msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
  error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H")
)

msig_h_gene_col <- if ("gene_symbol" %in% names(msig_h)) {
  "gene_symbol"
} else if ("human_gene_symbol" %in% names(msig_h)) {
  "human_gene_symbol"
} else {
  stop("No gene-symbol column found in Hallmark msigdbr output.")
}

msig_h_name_col <- if ("gs_name" %in% names(msig_h)) {
  "gs_name"
} else if ("geneset" %in% names(msig_h)) {
  "geneset"
} else if ("gene_set_name" %in% names(msig_h)) {
  "gene_set_name"
} else {
  stop("No gene-set name column found in Hallmark msigdbr output.")
}

msig_h_clean <- msig_h %>%
  dplyr::mutate(
    gene_symbol_clean = toupper(trimws(.data[[msig_h_gene_col]])),
    pathway_name_clean = as.character(.data[[msig_h_name_col]])
  ) %>%
  dplyr::filter(!is.na(gene_symbol_clean), gene_symbol_clean != "")

hallmark_pathways <- split(
  msig_h_clean$gene_symbol_clean,
  msig_h_clean$pathway_name_clean
)

hallmark_pathways <- lapply(hallmark_pathways, function(x) sort(unique(x)))

msig_c2 <- tryCatch(
  msigdbr::msigdbr(species = "Homo sapiens", collection = "C2"),
  error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "C2")
)

msig_c2_gene_col <- if ("gene_symbol" %in% names(msig_c2)) {
  "gene_symbol"
} else if ("human_gene_symbol" %in% names(msig_c2)) {
  "human_gene_symbol"
} else {
  stop("No gene-symbol column found in C2 msigdbr output.")
}

msig_c2_name_col <- if ("gs_name" %in% names(msig_c2)) {
  "gs_name"
} else if ("geneset" %in% names(msig_c2)) {
  "geneset"
} else if ("gene_set_name" %in% names(msig_c2)) {
  "gene_set_name"
} else {
  stop("No gene-set name column found in C2 msigdbr output.")
}

msig_c2_subcollection_col <- if ("gs_subcollection" %in% names(msig_c2)) {
  "gs_subcollection"
} else if ("gs_subcat" %in% names(msig_c2)) {
  "gs_subcat"
} else if ("subcategory" %in% names(msig_c2)) {
  "subcategory"
} else {
  NA_character_
}

msig_c2_clean <- msig_c2 %>%
  dplyr::mutate(
    gene_symbol_clean = toupper(trimws(.data[[msig_c2_gene_col]])),
    pathway_name_clean = as.character(.data[[msig_c2_name_col]]),
    subcollection_clean = if (!is.na(msig_c2_subcollection_col)) {
      as.character(.data[[msig_c2_subcollection_col]])
    } else {
      NA_character_
    }
  ) %>%
  dplyr::filter(!is.na(gene_symbol_clean), gene_symbol_clean != "")

msig_reactome <- msig_c2_clean %>%
  dplyr::filter(
    stringr::str_detect(toupper(pathway_name_clean), "^REACTOME") |
      stringr::str_detect(toupper(subcollection_clean), "REACTOME")
  )

if (nrow(msig_reactome) == 0) {
  stop("No Reactome pathways found in msigdbr C2 output.")
}

reactome_pathways <- split(
  msig_reactome$gene_symbol_clean,
  msig_reactome$pathway_name_clean
)

reactome_pathways <- lapply(reactome_pathways, function(x) sort(unique(x)))

msig_c5 <- tryCatch(
  msigdbr::msigdbr(species = "Homo sapiens", collection = "C5"),
  error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "C5")
)

msig_c5_gene_col <- if ("gene_symbol" %in% names(msig_c5)) {
  "gene_symbol"
} else if ("human_gene_symbol" %in% names(msig_c5)) {
  "human_gene_symbol"
} else {
  stop("No gene-symbol column found in C5 msigdbr output.")
}

msig_c5_name_col <- if ("gs_name" %in% names(msig_c5)) {
  "gs_name"
} else if ("geneset" %in% names(msig_c5)) {
  "geneset"
} else if ("gene_set_name" %in% names(msig_c5)) {
  "gene_set_name"
} else {
  stop("No gene-set name column found in C5 msigdbr output.")
}

msig_c5_clean <- msig_c5 %>%
  dplyr::mutate(
    gene_symbol_clean = toupper(trimws(.data[[msig_c5_gene_col]])),
    pathway_name_clean = as.character(.data[[msig_c5_name_col]])
  ) %>%
  dplyr::filter(!is.na(gene_symbol_clean), gene_symbol_clean != "")

go_pathways <- split(
  msig_c5_clean$gene_symbol_clean,
  msig_c5_clean$pathway_name_clean
)

go_pathways <- lapply(go_pathways, function(x) sort(unique(x)))

chemokine_activity_info_raw <- get_required_pathway_info(
  requested_pathway = "GOMF_CHEMOKINE_ACTIVITY",
  named_pathways = go_pathways,
  loose_patterns = c("^GOMF_CHEMOKINE_ACTIVITY$", "CHEMOKINE_ACTIVITY"),
  source_label = "MSigDB_C5_GOMF_CHEMOKINE_ACTIVITY"
)

canonical_chemokine_regex <- paste0(
  "^(",
  "CCL[0-9]+[A-Z0-9]*|",
  "CXCL[0-9]+[A-Z0-9]*|",
  "CX3CL[0-9]+[A-Z0-9]*|",
  "XCL[0-9]+[A-Z0-9]*|",
  "PF4|PF4V1|PPBP",
  ")$"
)

chemokine_activity_raw_tbl <- tibble::tibble(
  Gene = chemokine_activity_info_raw$genes,
  source_gene_set = chemokine_activity_info_raw$resolved_pathway,
  retained_for_heatmap = stringr::str_detect(Gene, canonical_chemokine_regex)
)

readr::write_csv(
  chemokine_activity_raw_tbl,
  file.path(out_dir, "GOMF_CHEMOKINE_ACTIVITY_raw_and_retained_genes.csv")
)

chemokine_genes <- chemokine_activity_raw_tbl %>%
  dplyr::filter(retained_for_heatmap) %>%
  dplyr::arrange(Gene) %>%
  dplyr::pull(Gene)

if (length(chemokine_genes) == 0) {
  stop("GOMF_CHEMOKINE_ACTIVITY was found, but no canonical chemokine ligands remained after filtering.")
}

chemokine_info <- list(
  requested_pathway = "GOMF_CHEMOKINE_ACTIVITY",
  resolved_pathway = chemokine_activity_info_raw$resolved_pathway,
  source = "MSigDB_C5_GOMF_CHEMOKINE_ACTIVITY_canonical_ligands_only",
  genes = chemokine_genes
)

readr::write_csv(
  tibble::tibble(
    Gene = chemokine_info$genes,
    source_gene_set = chemokine_info$resolved_pathway,
    source_note = "Canonical chemokine ligand subset retained from MSigDB GOMF_CHEMOKINE_ACTIVITY: CCL, CXCL, CX3CL, XCL, PF4, PF4V1, PPBP."
  ),
  file.path(out_dir, "CHEMOKINE_GENES_USED_GOMF_CHEMOKINE_ACTIVITY_CANONICAL_ONLY.csv")
)

antimicrobial_info <- get_required_pathway_info(
  requested_pathway = "REACTOME_ANTIMICROBIAL_PEPTIDES",
  named_pathways = reactome_pathways,
  loose_patterns = c("^REACTOME_ANTIMICROBIAL_PEPTIDES$", "ANTIMICROBIAL_PEPTIDES"),
  source_label = "MSigDB_C2_REACTOME_ANTIMICROBIAL_PEPTIDES"
)

readr::write_csv(
  tibble::tibble(
    Gene = antimicrobial_info$genes,
    source_gene_set = antimicrobial_info$resolved_pathway
  ),
  file.path(out_dir, "REACTOME_ANTIMICROBIAL_PEPTIDES_GENES_USED.csv")
)

read_text_table_auto <- function(file) {
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

  if (!is.null(dat) && ncol(dat) > 1) return(dat)

  dat <- tryCatch(
    readr::read_csv(
      file,
      show_col_types = FALSE,
      progress = FALSE,
      name_repair = "unique"
    ),
    error = function(e) NULL
  )

  if (!is.null(dat) && ncol(dat) > 1) return(dat)

  dat <- tryCatch(
    utils::read.table(
      file,
      header = TRUE,
      sep = "",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    error = function(e) NULL
  )

  if (!is.null(dat) && ncol(dat) > 1) {
    return(tibble::as_tibble(dat, .name_repair = "unique"))
  }

  stop("Could not read DE table: ", file)
}

find_de_file <- function(comparison) {
  if (!comparison %in% names(de_file_map)) {
    return(NA_character_)
  }

  f <- file.path(de_dir, de_file_map[[comparison]])
  if (file.exists(f)) return(f)

  hits <- list.files(
    de_dir,
    pattern = paste0(".*", safe_filename(comparison), ".*DE.*results.*\\.txt$"),
    full.names = TRUE,
    recursive = FALSE
  )

  if (length(hits) > 0) return(hits[1])
  NA_character_
}

standardize_de_cols <- function(dat, comparison, source_file) {
  dat <- tibble::as_tibble(dat, .name_repair = "unique")

  dat <- rename_ci(dat, "Gene", c(
    "Gene", "gene", "SYMBOL", "Symbol", "symbol",
    "gene_symbol", "GeneSymbol", "geneSymbol",
    "genes", "Genes", "feature", "Feature",
    "ID", "id", "Name", "name"
  ))

  if (!"Gene" %in% names(dat)) {
    first_col <- names(dat)[1]
    first_vals <- as.character(dat[[first_col]])

    if (
      length(first_vals) > 0 &&
      any(!is.na(first_vals)) &&
      any(stringr::str_detect(first_vals[!is.na(first_vals)], "[A-Za-z]"))
    ) {
      names(dat)[1] <- "Gene"
    }
  }

  dat <- rename_ci(dat, "logFC", c(
    "logFC", "LogFC", "logfc", "log2FC", "Log2FC",
    "log2FoldChange", "log2_fold_change", "LFC"
  ))

  dat <- rename_ci(dat, "logCPM", c(
    "logCPM", "LogCPM", "logcpm", "log2CPM", "Log2CPM",
    "AveExpr", "avgExpr", "average_expression"
  ))

  dat <- rename_ci(dat, "PValue", c(
    "PValue", "P.Value", "Pvalue", "pvalue", "p.value", "p"
  ))

  dat <- rename_ci(dat, "FDR", c(
    "FDR", "fdr", "padj", "adj.P.Val", "adj_p", "qvalue", "q.value"
  ))

  if (!"Gene" %in% names(dat)) {
    stop("DE table is missing a gene-symbol column for ", comparison, "\nFile: ", source_file)
  }

  if (!"logFC" %in% names(dat)) {
    stop("DE table is missing logFC for ", comparison, "\nFile: ", source_file)
  }

  dat$Gene <- toupper(trimws(as.character(dat$Gene)))
  dat$logFC <- suppressWarnings(as.numeric(dat$logFC))

  if ("logCPM" %in% names(dat)) dat$logCPM <- suppressWarnings(as.numeric(dat$logCPM))
  if ("PValue" %in% names(dat)) dat$PValue <- suppressWarnings(as.numeric(dat$PValue))
  if ("FDR" %in% names(dat)) dat$FDR <- suppressWarnings(as.numeric(dat$FDR))

  dat$comparison <- comparison
  dat$source_file <- basename(source_file)

  dat %>%
    dplyr::filter(
      !is.na(Gene),
      Gene != "",
      is.finite(logFC)
    ) %>%
    dplyr::group_by(Gene) %>%
    dplyr::arrange(dplyr::desc(abs(logFC)), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
}

read_de_table <- function(comparison) {
  f <- find_de_file(comparison)

  if (is.na(f)) {
    stop(
      "Could not find DE table for comparison: ",
      comparison,
      "\nExpected in: ",
      de_dir
    )
  }

  dat <- read_text_table_auto(f)

  standardize_de_cols(
    dat = dat,
    comparison = comparison,
    source_file = f
  )
}

required_de <- unique(c(
  panelB_reference_comparison,
  panelB_comparisons,
  panelD_comparisons
))

missing_de <- required_de[
  is.na(vapply(required_de, find_de_file, character(1)))
]

if (length(missing_de) > 0) {
  stop("Missing individual DE result tables for: ", paste(missing_de, collapse = ", "))
}

de_table_qc <- purrr::map_dfr(required_de, function(comp) {
  dat <- read_de_table(comp)

  tibble(
    comparison = comp,
    source_file = unique(dat$source_file)[1],
    n_genes = nrow(dat),
    has_logCPM = "logCPM" %in% names(dat),
    min_logCPM = if ("logCPM" %in% names(dat)) min(dat$logCPM, na.rm = TRUE) else NA_real_,
    median_logCPM = if ("logCPM" %in% names(dat)) median(dat$logCPM, na.rm = TRUE) else NA_real_,
    max_logFC = max(dat$logFC, na.rm = TRUE),
    min_logFC = min(dat$logFC, na.rm = TRUE)
  )
})

readr::write_csv(
  de_table_qc,
  file.path(out_dir, "DE_tables_used_for_heatmaps_QC.csv")
)

find_ranked_file <- function(comparison, mode = gsea_filter_mode) {
  f <- file.path(
    ranked_dir,
    paste0(safe_filename(comparison), "_", mode, "_ranked_gene_table.csv")
  )

  if (file.exists(f)) return(f)

  hits <- list.files(
    ranked_dir,
    pattern = paste0("^", safe_filename(comparison), ".*", mode, ".*ranked.*gene.*table.*\\.csv$"),
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(hits) > 0) return(hits[1])
  NA_character_
}

standardize_ranked_cols <- function(dat, comparison, mode = gsea_filter_mode) {
  dat <- rename_ci(dat, "Gene", c("Gene", "gene", "SYMBOL", "Symbol", "symbol", "gene_symbol"))
  dat <- rename_ci(dat, "logFC", c("logFC", "LogFC", "logfc", "log2FoldChange", "log2FC"))
  dat <- rename_ci(dat, "logCPM", c("logCPM", "LogCPM", "logcpm", "log2CPM", "Log2CPM", "AveExpr"))
  dat <- rename_ci(dat, "LR", c("LR", "LRT", "likelihood_ratio", "LikelihoodRatio"))
  dat <- rename_ci(dat, "PValue", c("PValue", "P.Value", "Pvalue", "pvalue", "p.value", "p"))
  dat <- rename_ci(dat, "rank_metric", c("rank_metric", "rank", "stat", "STAT"))

  if (!"Gene" %in% names(dat)) stop("Ranked table is missing a Gene/SYMBOL column for ", comparison)
  if (!"logFC" %in% names(dat)) stop("Ranked table is missing logFC for ", comparison)

  dat$Gene <- toupper(trimws(as.character(dat$Gene)))
  dat$logFC <- suppressWarnings(as.numeric(dat$logFC))

  if ("logCPM" %in% names(dat)) dat$logCPM <- suppressWarnings(as.numeric(dat$logCPM))
  if ("LR" %in% names(dat)) dat$LR <- suppressWarnings(as.numeric(dat$LR))
  if ("PValue" %in% names(dat)) dat$PValue <- suppressWarnings(as.numeric(dat$PValue))

  if (!"rank_metric" %in% names(dat)) {
    if ("LR" %in% names(dat)) {
      dat$rank_metric <- sign(dat$logFC) * sqrt(pmax(dat$LR, 0))
    } else if ("PValue" %in% names(dat)) {
      min_nonzero_p <- suppressWarnings(min(dat$PValue[dat$PValue > 0], na.rm = TRUE))
      if (!is.finite(min_nonzero_p)) min_nonzero_p <- 1e-300

      dat$rank_metric <- sign(dat$logFC) *
        -log10(ifelse(dat$PValue <= 0 | is.na(dat$PValue), min_nonzero_p * 0.1, dat$PValue))
    } else {
      stop("Ranked table has no rank_metric, LR, or PValue column for ", comparison)
    }
  }

  dat$rank_metric <- suppressWarnings(as.numeric(dat$rank_metric))
  dat$comparison <- comparison
  dat$filter_mode <- mode

  dat %>%
    dplyr::filter(
      !is.na(Gene),
      Gene != "",
      is.finite(logFC),
      is.finite(rank_metric)
    ) %>%
    dplyr::group_by(Gene) %>%
    dplyr::arrange(
      dplyr::desc(abs(rank_metric)),
      dplyr::desc(abs(logFC)),
      .by_group = TRUE
    ) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
}

read_ranked_table <- function(comparison, mode = gsea_filter_mode) {
  f <- find_ranked_file(comparison, mode)

  if (is.na(f)) {
    stop("Could not find ranked gene table for ", comparison, " / ", mode, " in: ", ranked_dir)
  }

  dat <- readr::read_csv(f, show_col_types = FALSE)

  out <- standardize_ranked_cols(
    dat,
    comparison = comparison,
    mode = mode
  )

  out$source_file <- basename(f)
  out
}

read_ranks_one <- function(comparison, mode = gsea_filter_mode) {
  dat <- read_ranked_table(comparison, mode)

  ranks <- dat$rank_metric
  names(ranks) <- dat$Gene

  ranks <- ranks[
    !is.na(names(ranks)) &
      names(ranks) != "" &
      is.finite(ranks)
  ]

  sort(ranks, decreasing = TRUE)
}

required_ranked <- panelC_comparisons

missing_ranked <- required_ranked[
  is.na(vapply(required_ranked, find_ranked_file, character(1), mode = gsea_filter_mode))
]

if (length(missing_ranked) > 0) {
  stop("Missing LOGCPM_GE_2 ranked gene tables for Panel C: ", paste(missing_ranked, collapse = ", "))
}

ranked_table_qc <- purrr::map_dfr(required_ranked, function(comp) {
  dat <- read_ranked_table(comp, gsea_filter_mode)

  tibble(
    comparison = comp,
    filter_mode = gsea_filter_mode,
    source_file = unique(dat$source_file)[1],
    n_genes = nrow(dat),
    min_logCPM = if ("logCPM" %in% names(dat)) min(dat$logCPM, na.rm = TRUE) else NA_real_,
    median_logCPM = if ("logCPM" %in% names(dat)) median(dat$logCPM, na.rm = TRUE) else NA_real_,
    max_logFC = max(dat$logFC, na.rm = TRUE),
    min_logFC = min(dat$logFC, na.rm = TRUE)
  )
})

readr::write_csv(
  ranked_table_qc,
  file.path(out_dir, "LOGCPM_GE_2_ranked_tables_used_for_panelC_QC.csv")
)

available_panelA_pathways <- combined_fgsea %>%
  dplyr::filter(
    filter_mode == gsea_filter_mode,
    comparison %in% panelA_comparisons
  ) %>%
  dplyr::pull(pathway) %>%
  unique()

panelA_pathway_map <- resolve_pathway_list(
  targets = panelA_requested_pathways,
  available_pathways = available_panelA_pathways
)

readr::write_csv(
  panelA_pathway_map,
  file.path(out_dir, "Panel_A_requested_pathway_resolution_LOGCPM_GE_2.csv")
)

panelA_missing_pathways <- panelA_pathway_map %>%
  dplyr::filter(is.na(resolved_pathway))

if (nrow(panelA_missing_pathways) > 0) {
  stop(
    "The following required Panel A pathways were not found: ",
    paste(panelA_missing_pathways$requested_pathway, collapse = "; ")
  )
}

panelA_resolved_tbl <- panelA_pathway_map %>%
  dplyr::filter(!is.na(resolved_pathway)) %>%
  dplyr::distinct(resolved_pathway, .keep_all = TRUE)

panelA_heatmap_pathways <- panelA_resolved_tbl$resolved_pathway

panelA_results <- combined_fgsea %>%
  dplyr::filter(
    comparison %in% panelA_comparisons,
    filter_mode == gsea_filter_mode,
    pathway %in% panelA_heatmap_pathways
  ) %>%
  dplyr::group_by(comparison, pathway) %>%
  dplyr::arrange(padj, pval, .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

panelA_df <- tidyr::expand_grid(
  comparison = panelA_comparisons,
  pathway = panelA_heatmap_pathways
) %>%
  dplyr::left_join(
    panelA_resolved_tbl %>%
      dplyr::select(pathway = resolved_pathway, requested_pathway),
    by = "pathway"
  ) %>%
  dplyr::left_join(
    panelA_results,
    by = c("comparison", "pathway")
  ) %>%
  dplyr::mutate(
    comparison_label = factor(
      label_comp(comparison),
      levels = label_comp(panelA_comparisons)
    ),
    pathway_label = factor(
      short_pathway_label(pathway),
      levels = rev(short_pathway_label(panelA_heatmap_pathways))
    ),
    star_label = sig_stars(padj),
    nes_star_label = ifelse(
      is.na(NES),
      "",
      paste0(sprintf("%.2f", NES), sig_stars(padj))
    ),
    label = dplyr::case_when(
      panelA_label_mode == "stars" ~ star_label,
      panelA_label_mode == "nes_and_stars" ~ nes_star_label,
      TRUE ~ star_label
    ),
    label_colour = panelA_star_colour
  )

missing_panelA_cells <- panelA_df %>%
  dplyr::filter(is.na(NES))

if (nrow(missing_panelA_cells) > 0) {
  stop("Panel A has missing NES values for required comparison/pathway combinations. Check FGSEA files.")
}

readr::write_csv(
  panelA_results,
  file.path(out_dir, "Panel_A_selected_gene_sets_LOGCPM_GE_2_table.csv")
)

panelA <- ggplot(
  panelA_df,
  aes(x = comparison_label, y = pathway_label, fill = NES)
) +
  geom_tile(
    color = heatmap_tile_border_colour,
    linewidth = heatmap_tile_border_linewidth,
    width = heatmap_tile_fraction,
    height = heatmap_tile_fraction
  ) +
  geom_text(
    aes(label = label, color = label_colour),
    size = geom_text_size_from_pt(5.0),
    family = font_family,
    fontface = "bold"
  ) +
  scale_color_identity() +
  scale_x_discrete(
    drop = FALSE,
    expand = expansion(add = 0.05)
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = function(x) stringr::str_wrap(x, width = 23),
    expand = expansion(add = 0.08)
  ) +
  scale_fill_gradient2(
    low = heatmap_negative_colour,
    mid = heatmap_mid_colour,
    high = heatmap_positive_colour,
    midpoint = 0,
    limits = c(-3, 3),
    oob = scales::squish,
    na.value = "grey95",
    name = "NES",
    guide = vertical_heatmap_guide(barheight_mm = 18, barwidth_mm = 2.0)
  ) +
  coord_fixed(ratio = 1, clip = "off") +
  labs(
    title = NULL,
    x = NULL,
    y = NULL
  ) +
  nature_theme(base_size = font_base_size) +
  theme(
    axis.text.x = element_text(
      family = font_family,
      size = panelA_axis_text_pt,
      angle = panelA_x_label_angle,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      family = font_family,
      size = panelA_axis_text_pt,
      lineheight = 0.9,
      color = "black"
    ),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(
      family = font_family,
      size = legend_title_size,
      angle = 90,
      hjust = 0.5,
      vjust = 0.5,
      color = "black"
    ),
    legend.text = element_text(
      family = font_family,
      size = legend_text_size,
      color = "black"
    ),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  )

save_plot_set(
  panelA,
  "Panel_A_Selected_Gene_Sets_LOGCPM_GE_2_IMAGE_MATCH_CENTERED_LEGEND",
  88,
  115
)

interferon_alpha_info <- get_required_pathway_info(
  requested_pathway = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  named_pathways = hallmark_pathways,
  source_label = "MSigDB_HALLMARK"
)

tnfa_nfkb_info <- get_required_pathway_info(
  requested_pathway = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  named_pathways = hallmark_pathways,
  source_label = "MSigDB_HALLMARK"
)

inflammatory_info <- get_required_pathway_info(
  requested_pathway = "HALLMARK_INFLAMMATORY_RESPONSE",
  named_pathways = hallmark_pathways,
  source_label = "MSigDB_HALLMARK"
)

gene_module_info <- list(
  "Interferon alpha response" = interferon_alpha_info,
  "Chemokines" = chemokine_info,
  "Antimicrobial peptides" = antimicrobial_info,
  "TNFA/NFKB signaling" = tnfa_nfkb_info
)

gene_module_source_tbl <- purrr::imap_dfr(gene_module_info, function(info, group_name) {
  tibble(
    gene_group = group_name,
    requested_pathway = info$requested_pathway,
    resolved_pathway = info$resolved_pathway,
    gene_source = info$source,
    n_candidate_genes = length(info$genes)
  )
})

readr::write_csv(
  gene_module_source_tbl,
  file.path(out_dir, "Panel_BD_gene_module_sources.csv")
)

reference_de_tbl <- read_de_table(panelB_reference_comparison)

if (!"logCPM" %in% names(reference_de_tbl)) {
  stop(
    "Reference DE table is missing logCPM. Cannot enforce logCPM >= ",
    panelB_logcpm_threshold,
    ". Reference comparison: ",
    panelB_reference_comparison
  )
}

reference_for_join <- reference_de_tbl %>%
  dplyr::select(
    Gene,
    SSL_HP_v_H_COL_logFC = logFC,
    SSL_HP_v_H_COL_logCPM = logCPM,
    dplyr::any_of(c("PValue", "FDR"))
  )

if ("PValue" %in% names(reference_for_join)) {
  names(reference_for_join)[names(reference_for_join) == "PValue"] <- "SSL_HP_v_H_COL_PValue"
}

if ("FDR" %in% names(reference_for_join)) {
  names(reference_for_join)[names(reference_for_join) == "FDR"] <- "SSL_HP_v_H_COL_FDR"
}

candidate_qc <- purrr::imap_dfr(gene_module_info, function(info, grp) {
  genes <- unique(toupper(trimws(info$genes)))

  tibble(
    gene_group = grp,
    requested_pathway = info$requested_pathway,
    resolved_pathway = info$resolved_pathway,
    gene_source = info$source,
    Gene = genes
  ) %>%
    dplyr::left_join(reference_for_join, by = "Gene") %>%
    dplyr::mutate(
      passes_panelB_filter =
        !is.na(SSL_HP_v_H_COL_logCPM) &
        is.finite(SSL_HP_v_H_COL_logCPM) &
        SSL_HP_v_H_COL_logCPM >= panelB_logcpm_threshold &
        !is.na(SSL_HP_v_H_COL_logFC) &
        is.finite(SSL_HP_v_H_COL_logFC) &
        SSL_HP_v_H_COL_logFC > panelB_logfc_threshold
    )
})

readr::write_csv(
  candidate_qc,
  file.path(out_dir, "Panel_BD_candidate_gene_filter_QC_DE_logCPM_ge_2_logFC_gt_1.csv")
)

panelB_selected_all_tbl <- candidate_qc %>%
  dplyr::filter(passes_panelB_filter) %>%
  dplyr::arrange(
    gene_group,
    dplyr::desc(SSL_HP_v_H_COL_logFC),
    Gene
  ) %>%
  dplyr::group_by(gene_group) %>%
  dplyr::mutate(selection_rank = dplyr::row_number()) %>%
  dplyr::ungroup()

readr::write_csv(
  panelB_selected_all_tbl,
  file.path(out_dir, "Panel_BD_ALL_selected_genes_FROM_DE_logCPM_ge_2_logFC_gt_1.csv")
)

panelB_display_tbl <- panelB_selected_all_tbl %>%
  dplyr::group_by(gene_group) %>%
  dplyr::filter(selection_rank <= max_panelB_genes_per_group) %>%
  dplyr::ungroup()

readr::write_csv(
  panelB_display_tbl,
  file.path(out_dir, "Panel_BD_DISPLAYED_genes_top_by_logFC.csv")
)

gene_group_order <- names(gene_module_info)

panelB_gene_groups <- split(
  panelB_display_tbl$Gene,
  panelB_display_tbl$gene_group
)

for (grp in gene_group_order) {
  if (!grp %in% names(panelB_gene_groups)) {
    panelB_gene_groups[[grp]] <- character(0)
  }
}

panelB_gene_groups <- panelB_gene_groups[gene_group_order]

make_heatmap_df <- function(comparisons, genes, group_name) {
  de_all <- purrr::map_dfr(comparisons, function(comp) {
    read_de_table(comp) %>%
      dplyr::select(
        comparison,
        Gene,
        logFC,
        dplyr::any_of(c("logCPM", "PValue", "FDR")),
        source_file
      )
  })

  gene_index <- tibble(
    Gene = genes,
    gene_order = seq_along(genes),
    gene_id = paste(group_name, genes, seq_along(genes), sep = "__")
  )

  tidyr::expand_grid(
    comparison = comparisons,
    gene_index
  ) %>%
    dplyr::left_join(
      de_all,
      by = c("comparison", "Gene")
    ) %>%
    dplyr::mutate(
      comparison_label = factor(
        label_comp_heatmap(comparison),
        levels = label_comp_heatmap(comparisons)
      ),
      gene_id = factor(
        gene_id,
        levels = rev(gene_index$gene_id)
      ),
      logFC_colour = pmax(
        pmin(logFC, heatmap_colour_saturation_logfc),
        -heatmap_colour_saturation_logfc
      )
    )
}

build_one_row_heatmap_panel <- function(
    df_list,
    comparisons,
    legend_title = "log2FC",
    selected_modules = NULL,
    force_block_height = NULL,
    show_legend = TRUE
) {
  module_names <- gene_group_order[gene_group_order %in% names(df_list)]

  if (!is.null(selected_modules)) {
    module_names <- module_names[module_names %in% selected_modules]
  }

  module_names <- module_names[
    vapply(
      module_names,
      function(grp) {
        !is.null(df_list[[grp]]) && nrow(df_list[[grp]]) > 0
      },
      logical(1)
    )
  ]

  if (length(module_names) == 0) {
    return(make_empty_plot("Gene-level heatmaps", "No genes available"))
  }

  ncomp <- length(comparisons)
  heatmap_labels <- label_comp_heatmap(comparisons)
  max_label_len <- max(nchar(heatmap_labels))

  module_n_tbl <- purrr::map_dfr(module_names, function(grp) {
    tibble(
      gene_group = grp,
      n_genes = dplyr::n_distinct(df_list[[grp]]$Gene)
    )
  })

  if (!is.null(force_block_height) && is.finite(force_block_height)) {
    block_height <- force_block_height
  } else {
    block_height <- max(module_n_tbl$n_genes, na.rm = TRUE)
  }

  block_height <- max(block_height, 1)

  x_start <- 3.55

  x_gap <- if (ncomp >= 3) 6.75 else 6.25
  if (max_label_len > 18) x_gap <- x_gap + 0.55
  if (max_label_len > 24) x_gap <- x_gap + 0.65

  layout_tbl <- tibble::tibble(
    gene_group = module_names,
    layout_index = seq_along(module_names)
  ) %>%
    dplyr::left_join(module_n_tbl, by = "gene_group") %>%
    dplyr::mutate(
      x_offset = x_start + (layout_index - 1) * (ncomp + x_gap),
      y_top_boundary = block_height,
      y_bottom_boundary = 0
    )

  tile_half_width <- heatmap_tile_fraction / 2

  tile_data <- purrr::map_dfr(module_names, function(grp) {
    df <- df_list[[grp]]
    lay <- layout_tbl %>% dplyr::filter(gene_group == grp)

    genes_in_order <- df %>%
      dplyr::distinct(Gene, gene_order) %>%
      dplyr::arrange(gene_order) %>%
      dplyr::mutate(gene_order_plot = dplyr::row_number())

    comp_tbl <- tibble(
      comparison_label = factor(
        label_comp_heatmap(comparisons),
        levels = label_comp_heatmap(comparisons)
      ),
      x_local = seq_along(comparisons)
    )

    row_height <- block_height / max(lay$n_genes, 1)

    df %>%
      dplyr::left_join(genes_in_order, by = c("Gene", "gene_order")) %>%
      dplyr::left_join(comp_tbl, by = "comparison_label") %>%
      dplyr::mutate(
        x_plot = lay$x_offset + x_local - 1,
        y_plot = lay$y_top_boundary - (gene_order_plot - 0.5) * row_height,
        tile_height = row_height * heatmap_tile_fraction,
        xmin = x_plot - tile_half_width,
        xmax = x_plot + tile_half_width,
        ymin = y_plot - tile_height / 2,
        ymax = y_plot + tile_height / 2,
        gene_group = grp
      )
  })

  gene_label_data <- tile_data %>%
    dplyr::distinct(gene_group, Gene, y_plot) %>%
    dplyr::left_join(layout_tbl, by = "gene_group") %>%
    dplyr::mutate(
      x_plot = x_offset - 0.62
    )

  module_title_data <- layout_tbl %>%
    dplyr::mutate(
      x_plot = x_offset + (ncomp - 1) / 2,
      y_plot = y_top_boundary + 0.82,
      title_label = stringr::str_wrap(gene_group, width = 18)
    )

  x_label_data <- purrr::map_dfr(module_names, function(grp) {
    lay <- layout_tbl %>% dplyr::filter(gene_group == grp)

    tibble(
      gene_group = grp,
      comparison_label = label_comp_heatmap(comparisons),
      x_plot = lay$x_offset + seq_along(comparisons) - 1,
      y_plot = lay$y_bottom_boundary - 0.52
    )
  })

  x_min <- min(layout_tbl$x_offset) - 3.05
  x_max <- max(layout_tbl$x_offset) + ncomp - 0.05

  bottom_extra <- if (max_label_len > 18) 5.8 else 4.8
  if (max_label_len > 24) bottom_extra <- 6.6

  y_min <- min(x_label_data$y_plot) - bottom_extra
  y_max <- max(layout_tbl$y_top_boundary) + 1.55

  legend_position <- if (show_legend) "right" else "none"

  ggplot() +
    geom_rect(
      data = tile_data,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = logFC_colour
      ),
      color = heatmap_tile_border_colour,
      linewidth = heatmap_tile_border_linewidth
    ) +
    geom_text(
      data = gene_label_data,
      aes(x = x_plot, y = y_plot, label = Gene),
      hjust = 1,
      family = font_family,
      size = manual_heatmap_gene_label_size,
      color = "black"
    ) +
    geom_text(
      data = module_title_data,
      aes(x = x_plot, y = y_plot, label = title_label),
      hjust = 0.5,
      vjust = 0,
      family = font_family,
      fontface = "bold",
      size = manual_heatmap_module_title_size,
      lineheight = 0.90,
      color = "black"
    ) +
    geom_text(
      data = x_label_data,
      aes(x = x_plot, y = y_plot, label = comparison_label),
      hjust = 1,
      vjust = 1,
      angle = manual_heatmap_x_label_angle,
      family = font_family,
      size = manual_heatmap_x_label_size,
      color = "black"
    ) +
    scale_fill_gradient2(
      low = heatmap_negative_colour,
      mid = heatmap_mid_colour,
      high = heatmap_positive_colour,
      midpoint = 0,
      limits = c(-heatmap_colour_saturation_logfc, heatmap_colour_saturation_logfc),
      breaks = c(-heatmap_colour_saturation_logfc, 0, heatmap_colour_saturation_logfc),
      labels = c(
        paste0("≤ -", heatmap_colour_saturation_logfc),
        "0",
        paste0("≥ +", heatmap_colour_saturation_logfc)
      ),
      oob = scales::squish,
      na.value = "grey90",
      name = legend_title,
      guide = vertical_heatmap_guide(barheight_mm = 18, barwidth_mm = 2.0)
    ) +
    coord_fixed(
      ratio = 1,
      xlim = c(x_min, x_max),
      ylim = c(y_min, y_max),
      clip = "off"
    ) +
    theme_void(base_family = font_family) +
    theme(
      text = element_text(family = font_family, color = "black"),
      legend.position = legend_position,
      legend.direction = "vertical",
      legend.title = element_text(
        family = font_family,
        size = legend_title_size,
        angle = 90,
        hjust = 0.5,
        vjust = 0.5,
        color = "black"
      ),
      legend.text = element_text(
        family = font_family,
        size = legend_text_size,
        color = "black"
      ),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(2, 3, 2, 3, unit = "mm")
    )
}

save_module_heatmap_sets <- function(
    df_list,
    comparisons,
    stem_prefix,
    legend_title = "log2FC",
    individual_width_mm = 60,
    individual_height_mm = 105
) {
  module_names <- gene_group_order[gene_group_order %in% names(df_list)]

  module_names <- module_names[
    vapply(
      module_names,
      function(grp) {
        !is.null(df_list[[grp]]) && nrow(df_list[[grp]]) > 0
      },
      logical(1)
    )
  ]

  if (length(module_names) == 0) {
    return(list())
  }

  module_n_tbl <- purrr::map_dfr(module_names, function(grp) {
    tibble(
      gene_group = grp,
      n_genes = dplyr::n_distinct(df_list[[grp]]$Gene)
    )
  })

  common_block_height <- max(module_n_tbl$n_genes, na.rm = TRUE)

  plots <- purrr::map(module_names, function(grp) {
    p <- build_one_row_heatmap_panel(
      df_list = df_list[grp],
      comparisons = comparisons,
      legend_title = legend_title,
      selected_modules = grp,
      force_block_height = common_block_height,
      show_legend = TRUE
    )

    save_plot_set(
      p,
      paste0(stem_prefix, "_", safe_filename(grp)),
      individual_width_mm,
      individual_height_mm
    )

    p
  })

  names(plots) <- module_names
  plots
}

panelB_df_list <- list()

for (grp in names(panelB_gene_groups)) {
  genes <- panelB_gene_groups[[grp]]
  if (length(genes) == 0) next

  panelB_df_list[[grp]] <- make_heatmap_df(
    comparisons = panelB_comparisons,
    genes = genes,
    group_name = grp
  )
}

panelB_long <- if (length(panelB_df_list) > 0) {
  bind_rows(panelB_df_list, .id = "gene_group_name")
} else {
  tibble()
}

readr::write_csv(
  panelB_long,
  file.path(out_dir, "Panel_B_gene_heatmaps_FROM_DE_long_table.csv")
)

panelB_module_plots <- save_module_heatmap_sets(
  df_list = panelB_df_list,
  comparisons = panelB_comparisons,
  stem_prefix = "Panel_B_SINGLE_MODULE_HEATMAP_IMAGE_MATCH_CENTERED_LEGEND",
  legend_title = "log2FC",
  individual_width_mm = 60,
  individual_height_mm = 105
)

panelB <- build_one_row_heatmap_panel(
  df_list = panelB_df_list,
  comparisons = panelB_comparisons,
  legend_title = "log2FC",
  show_legend = TRUE
)

save_plot_set(
  panelB,
  "Panel_B_ALL4_Gene_Heatmaps_ONE_ROW_IMAGE_MATCH_CENTERED_LEGEND",
  183,
  82
)

calc_enrichment_curve <- function(ranks, pathway_genes) {
  ranks <- sort(ranks, decreasing = TRUE)

  pathway_genes <- unique(toupper(trimws(pathway_genes)))
  pathway_genes <- intersect(pathway_genes, names(ranks))

  N <- length(ranks)
  Nh <- length(pathway_genes)

  if (Nh < 1 || N <= Nh) {
    return(tibble())
  }

  hits <- names(ranks) %in% pathway_genes
  weights <- abs(ranks)
  NR <- sum(weights[hits])

  if (!is.finite(NR) || NR == 0) {
    weights <- rep(1, length(ranks))
    NR <- sum(weights[hits])
  }

  running_es <- cumsum(ifelse(hits, weights / NR, 0)) -
    cumsum(ifelse(!hits, 1 / (N - Nh), 0))

  tibble(
    position = seq_along(ranks),
    position_pct = 100 * position / N,
    running_es = running_es,
    hit = hits,
    gene = names(ranks)
  )
}

make_overlay_plot <- function(
    pathway_name,
    title_text,
    pathway_genes,
    fgsea_results
) {
  if (is.null(pathway_genes) || length(pathway_genes) == 0) {
    stop("No pathway genes supplied for: ", pathway_name)
  }

  resolved_fgsea_pathway <- resolve_pathway_name(
    pathway_name,
    fgsea_results$pathway
  )

  if (is.na(resolved_fgsea_pathway)) {
    resolved_fgsea_pathway <- pathway_name
  }

  curves <- purrr::map_dfr(panelC_comparisons, function(comp) {
    ranks <- read_ranks_one(comp, mode = gsea_filter_mode)

    curve <- calc_enrichment_curve(
      ranks = ranks,
      pathway_genes = pathway_genes
    )

    fg <- fgsea_results %>%
      dplyr::filter(
        comparison == comp,
        pathway == resolved_fgsea_pathway
      ) %>%
      dplyr::slice(1)

    nes <- if (nrow(fg) > 0) fg$NES[1] else NA_real_
    padj <- if (nrow(fg) > 0) fg$padj[1] else NA_real_

    curve %>%
      dplyr::mutate(
        comparison = comp,
        legend_label = paste0(
          panelC_short_labels[comp],
          "\nNES=",
          sprintf("%.2f", nes),
          ", FDR=",
          fmt_p(padj)
        )
      )
  })

  if (nrow(curves) == 0) {
    return(make_empty_plot(title_text, "No enrichment curve available"))
  }

  hit_df <- curves %>%
    dplyr::filter(hit)

  legend_labs <- curves %>%
    dplyr::distinct(comparison, legend_label)

  lab_vec <- stats::setNames(
    legend_labs$legend_label,
    legend_labs$comparison
  )

  ggplot(curves, aes(x = position_pct, y = running_es, color = comparison)) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = thin_line,
      color = "grey45"
    ) +
    geom_line(linewidth = 0.78) +
    geom_rug(
      data = hit_df,
      aes(x = position_pct, color = comparison),
      sides = "b",
      inherit.aes = FALSE,
      alpha = 1,
      linewidth = 0.18
    ) +
    scale_color_manual(
      values = panelC_cols,
      breaks = panelC_comparisons,
      labels = lab_vec[panelC_comparisons],
      guide = guide_legend(
        ncol = 1,
        byrow = TRUE,
        keyheight = grid::unit(3.5, "mm"),
        keywidth = grid::unit(7, "mm")
      )
    ) +
    scale_x_continuous(
      breaks = c(0, 25, 50, 75, 100),
      limits = c(0, 100),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.08))
    ) +
    labs(
      title = title_text,
      x = "Ranked gene list position (%)",
      y = "Running enrichment score",
      color = "ES+ vs ES−"
    ) +
    nature_theme(base_size = font_base_size) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(family = font_family, size = 5.0),
      legend.text = element_text(family = font_family, size = 4.6, lineheight = 0.85),
      legend.spacing.y = grid::unit(0.4, "mm"),
      legend.box.spacing = grid::unit(1, "mm"),
      plot.title = element_text(size = plot_title_size, face = "bold"),
      plot.margin = margin(2, 4, 8, 2, unit = "mm")
    )
}

panelC_tnfa <- make_overlay_plot(
  pathway_name = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  title_text = "TNFA/NFKB signaling",
  pathway_genes = tnfa_nfkb_info$genes,
  fgsea_results = fgsea_hallmark
)

panelC_ifna <- make_overlay_plot(
  pathway_name = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  title_text = "Interferon alpha response",
  pathway_genes = interferon_alpha_info$genes,
  fgsea_results = fgsea_hallmark
)

panelC_inflammatory <- make_overlay_plot(
  pathway_name = "HALLMARK_INFLAMMATORY_RESPONSE",
  title_text = "Inflammatory response",
  pathway_genes = inflammatory_info$genes,
  fgsea_results = fgsea_hallmark
)

save_plot_set(
  panelC_tnfa,
  "Panel_C_TNFA_NFKB_GSEA_overlay_LOGCPM_GE_2",
  95,
  82
)

save_plot_set(
  panelC_ifna,
  "Panel_C_IFNA_GSEA_overlay_LOGCPM_GE_2",
  95,
  82
)

save_plot_set(
  panelC_inflammatory,
  "Panel_C_INFLAMMATORY_RESPONSE_GSEA_overlay_LOGCPM_GE_2",
  95,
  82
)

panelC_vertical <- patchwork::wrap_plots(
  panelC_tnfa,
  panelC_ifna,
  panelC_inflammatory,
  ncol = 1,
  guides = "keep"
)

panelC_one_row <- patchwork::wrap_plots(
  panelC_tnfa,
  panelC_ifna,
  panelC_inflammatory,
  ncol = 3,
  guides = "keep"
)

save_plot_set(
  panelC_vertical,
  "Panel_C_GSEA_Overlays_TNFA_IFNA_INFLAMMATORY_VERTICAL_LOGCPM_GE_2",
  95,
  235
)

save_plot_set(
  panelC_one_row,
  "Panel_C_GSEA_Overlays_TNFA_IFNA_INFLAMMATORY_ONE_ROW_LOGCPM_GE_2",
  183,
  82
)

panelD_df_list <- list()

for (grp in names(panelB_gene_groups)) {
  genes <- panelB_gene_groups[[grp]]
  if (length(genes) == 0) next

  panelD_df_list[[grp]] <- make_heatmap_df(
    comparisons = panelD_comparisons,
    genes = genes,
    group_name = grp
  )
}

panelD_long <- if (length(panelD_df_list) > 0) {
  bind_rows(panelD_df_list, .id = "gene_group_name")
} else {
  tibble()
}

readr::write_csv(
  panelD_long,
  file.path(out_dir, "Panel_D_ES_gene_heatmaps_FROM_DE_long_table.csv")
)

panelD_module_plots <- save_module_heatmap_sets(
  df_list = panelD_df_list,
  comparisons = panelD_comparisons,
  stem_prefix = "Panel_D_SINGLE_MODULE_HEATMAP_IMAGE_MATCH_CENTERED_LEGEND",
  legend_title = "log2FC",
  individual_width_mm = 58,
  individual_height_mm = 105
)

panelD <- build_one_row_heatmap_panel(
  df_list = panelD_df_list,
  comparisons = panelD_comparisons,
  legend_title = "log2FC",
  show_legend = TRUE
)

save_plot_set(
  panelD,
  "Panel_D_ALL4_ES_Gene_Heatmaps_ONE_ROW_IMAGE_MATCH_CENTERED_LEGEND",
  183,
  82
)

top_row <- patchwork::wrap_plots(
  panelA,
  panelB,
  ncol = 2,
  widths = c(0.62, 1.55),
  guides = "keep"
)

final_figure <- patchwork::wrap_plots(
  top_row,
  panelC_one_row,
  panelD,
  ncol = 1,
  heights = c(1.05, 0.90, 0.90),
  guides = "keep"
) +
  patchwork::plot_annotation(
    tag_levels = "a",
    theme = theme(
      text = element_text(family = font_family, color = "black"),
      plot.tag = element_text(
        family = font_family,
        face = "bold",
        size = panel_label_size
      )
    )
  ) &
  theme(
    text = element_text(family = font_family)
  )

save_plot_set(
  final_figure,
  "M26F_Figure3_FINAL_IMAGE_MATCH_ONE_ROW_HEATMAPS_PANELC_INFLAMMATORY_CENTERED_LEGENDS_183mm",
  183,
  235
)

save_plot_set(
  final_figure,
  "M26F_Figure3_FINAL_IMAGE_MATCH_ONE_ROW_HEATMAPS_PANELC_INFLAMMATORY_CENTERED_LEGENDS_220mm",
  220,
  245
)

writeLines(
  c(
    paste0("Output directory: ", out_dir),
    paste0("Paper-ready analysis directory used: ", paper_ready_dir),
    paste0("FGSEA per-comparison directory used: ", fgsea_by_comparison_dir),
    paste0("Individual DE table directory used for Panels B/D: ", de_dir),
    paste0("Ranked gene table directory used for Panel C only: ", ranked_dir),
    paste0("Filter mode used for FGSEA and GSEA overlays: ", gsea_filter_mode),
    "",
    "No fallback gene lists are used.",
    "If a required Hallmark, Reactome, KEGG or GO/MSigDB gene set is unavailable, this script stops.",
    "",
    "Panel A:",
    "  Uses LOGCPM_GE_2 HALLMARK, KEGG and REACTOME FGSEA files from tables/fgsea_by_comparison.",
    "  Includes all 13 requested pathways.",
    "  Uses the same font, label size, tilted x-label angle, heatmap colours and right-side legend style as Panels B and D.",
    "  The NES legend title is centered beside the vertical colour bar.",
    "",
    "Panels B and D:",
    "  Heatmaps use individual DE result .txt files from the Re-run GSEA folder.",
    paste0("  Reference DE comparison for gene selection: ", panelB_reference_comparison),
    paste0("  Reference DE file: ", de_file_map[[panelB_reference_comparison]]),
    paste0("  Gene selection rule: logCPM >= ", panelB_logcpm_threshold, " and logFC > ", panelB_logfc_threshold, " in ", panelB_reference_comparison, "."),
    paste0("  Main figure displays top ", max_panelB_genes_per_group, " genes per module by SSL_HP_v_H_COL logFC."),
    "  Panel D uses the exact same displayed genes as Panel B.",
    "  Panel B modules are arranged in one horizontal line.",
    "  Panel D modules are arranged in one horizontal line.",
    "  Each Panel B module heatmap is saved individually.",
    "  Each Panel D module heatmap is saved individually.",
    "  The all-four-module Panel B and Panel D one-row heatmaps are saved separately.",
    "  Manual x labels match Panel A font family, label size and tilted orientation.",
    "  The log2FC legend title is centered beside the vertical colour bar.",
    "",
    "Panel C:",
    "  Uses the same enrichment overlay plotting approach as before.",
    "  Now includes TNFA/NFKB signaling, interferon alpha response and Hallmark inflammatory response.",
    "  Panel C is saved as individual plots, a vertical stack and a one-row layout.",
    "",
    "Shared heatmap styling across Panels A, B and D:",
    paste0("  font family = ", font_family),
    paste0("  label size = ", shared_heatmap_label_pt, " pt"),
    paste0("  x-label angle = ", shared_heatmap_x_label_angle, " degrees"),
    paste0("  low = ", heatmap_negative_colour),
    paste0("  mid = ", heatmap_mid_colour),
    paste0("  high = ", heatmap_positive_colour),
    paste0("  tile border colour = ", heatmap_tile_border_colour),
    paste0("  tile border linewidth = ", heatmap_tile_border_linewidth),
    paste0("  tile fill fraction = ", heatmap_tile_fraction),
    "",
    "Centered legend-title fix:",
    "  All heatmap colour bars use vertical_heatmap_guide().",
    "  This sets title.position = 'right', title.hjust = 0.5, title.vjust = 0.5, and a title.theme with angle = 90, hjust = 0.5, vjust = 0.5.",
    "",
    "Exports:",
    "  Each panel and final figure are exported as PDF, PNG, TIFF, JPG and JPEG.",
    "  Panel B single-module heatmaps are exported individually.",
    "  Panel D single-module heatmaps are exported individually.",
    "  Final figure exported at 183 mm width and 220 mm width."
  ),
  con = file.path(out_dir, "RUN_NOTES.txt")
)

message("DONE. Output directory: ", out_dir)
