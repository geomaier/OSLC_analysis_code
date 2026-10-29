# Organoid RNA-seq: infection, BRAF status and serial passage.
# Run the complete script; final score tables and figures are written near the end.

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

# Configuration and reproducible random seed ----
set.seed(1)
options(stringsAsFactors = FALSE)
options(contrasts = c("contr.treatment", "contr.poly"))

base_dir <- required_path("OSLC_ORGANOID_CPM_DIR")

metaplasia_display_name <- "Serrated Lesion Gastric Metaplasia Gene Signature"
metaplasia_file_stem <- "Serrated_Lesion_Gastric_Metaplasia"
metaplasia_pathway_id <- metaplasia_display_name

files <- list(
  infection_cpm = file.path(base_dir, "CGG028_Alt_Infected_cpm_CORRECTED.txt"),
  infection_de = file.path(base_dir, "CGG028_Alt_Infected_DE_results.txt"),
  timecourse_cpm = file.path(base_dir, "CGG028_cpm_table.txt"),
  braf_de = file.path(base_dir, "CGG028_First_Passage_SSL_HP_V600E_WT.txt"),
  infection_meta = file.path(base_dir, "Metadata_Infection_Organoids.csv"),
  timecourse_meta = file.path(base_dir, "Metadata_Timecourse_Organoids.csv")
)

metaplasia_signature_candidates <- c(
  file.path(base_dir, "FINAL_SSL_METAPLASIA_LIST.txt"),
  file.path(base_dir, "FINAL_SSL_METAPLASIA_LIST.csv"),
  Sys.getenv("OSLC_GASTRIC_SIGNATURE_FILE")
)

out_dir <- review_output_dir(file.path(base_dir, "Ch5_Organoid_FINAL"))
overwrite_existing_outputs <- FALSE
install_missing_packages <- FALSE

input_scale_override <- c(
  infection = "AUTO",
  timecourse = "AUTO"
)

pca_n_variable_genes <- 2000L
cpm_pseudocount <- 1
signature_mean_logcpm_cutoff <- 2
minimum_signature_genes <- 5L

# This switch affects the intermediate audit; final submitted scores use paired t-tests.
infection_score_primary_test <- "PAIRED_T"
if (!infection_score_primary_test %in% c("PAIRED_T", "EXACT_WILCOXON")) {
  stop("infection_score_primary_test must be PAIRED_T or EXACT_WILCOXON.")
}

gsea_logcpm_cutoff <- 2
fgsea_min_size <- 15L
fgsea_max_size <- 500L
fgsea_eps <- 0
fgsea_nperm_simple <- 10000L
fgsea_nproc <- 0L
hallmark_bar_top_n <- 12L
hallmark_heatmap_top_n <- 20L

module_heatmap_logcpm_cutoff <- 2
module_heatmap_logfc_cutoff <- 1
module_heatmap_max_genes_per_module <- 20L
module_heatmap_colour_limit <- 2

font_family <- "Arial"
font_pt <- 7
axis_font_pt <- 7
legend_font_pt <- 6
plot_title_pt <- 7
panel_tag_pt <- 8
point_size <- 2.2
point_stroke <- 0.35
export_dpi <- 600

paper_cols <- c(
  SSL_HP = "#FF0000",
  H_ADJ_SSL_HP = "#E9C61D",
  H_COL = "#1313EC",
  CA = "#825CA6",
  STOMACH = "#F155D0"
)

infection_cols <- c(
  "Non-Infected" = "#BDBDBD",
  "Infected" = "#00FF0E"
)

braf_cols <- c(
  "WT" = "#4C78A8",
  "V600E" = "#E45756"
)

comparison_curve_cols <- c(
  "Infected vs Non-Infected" = "#00FF0E",
  "BRAF V600E vs WT in SSL_HP" = "#E45756"
)

# Required packages ----
cran_packages <- c(
  "data.table", "dplyr", "tidyr", "stringr", "purrr", "ggplot2",
  "readr", "tibble", "patchwork", "scales", "msigdbr", "ggrepel"
)

bioc_packages <- c(
  "fgsea", "AnnotationDbi", "org.Hs.eg.db", "singscore", "limma"
)

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

for (pkg in cran_packages) install_if_missing(pkg, bioc = FALSE)
for (pkg in bioc_packages) install_if_missing(pkg, bioc = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(patchwork)
  library(scales)
  library(msigdbr)
  library(ggrepel)
  library(fgsea)
  library(singscore)
  library(grid)
})

if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

if (!dir.exists(base_dir)) {
  stop("Base directory does not exist: ", base_dir)
}

required_files <- unlist(files[c(
  "infection_cpm", "infection_de", "timecourse_cpm", "braf_de",
  "infection_meta", "timecourse_meta"
)])

missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required) > 0) {
  stop(
    "The following required files are missing:\n",
    paste(missing_required, collapse = "\n")
  )
}

if (overwrite_existing_outputs && dir.exists(out_dir)) {
  base_norm <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
  out_norm <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)

  if (identical(base_norm, out_norm)) {
    stop("Refusing to delete out_dir because it is identical to base_dir.")
  }

  unlink(out_dir, recursive = TRUE, force = TRUE)
  if (dir.exists(out_dir)) {
    stop(
      "Could not remove the previous output directory. Close any open output ",
      "files and rerun:\n", out_dir
    )
  }
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

subdirs <- list(
  qc = file.path(out_dir, "01_QC"),
  pca = file.path(out_dir, "02_PCA"),
  scores = file.path(out_dir, "03_Signature_scores"),
  gsea = file.path(out_dir, "04_GSEA"),
  heatmaps = file.path(out_dir, "05_Heatmaps"),
  combined = file.path(out_dir, "06_Combined_figures"),
  tables = file.path(out_dir, "07_Tables")
)

for (d in subdirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

input_role_qc <- data.frame(
  Input_key = names(files),
  File = unname(unlist(files)),
  Biological_role = c(
    "CPM expression matrix: infected vs non-infected organoids",
    "edgeR DE table: infected vs non-infected organoids",
    "CPM expression matrix: full organoid passage series",
    "edgeR DE table: first-passage BRAF V600E vs WT in SSL_HP",
    "Metadata: infection experiment",
    "Metadata: passage series and BRAF annotations"
  ),
  Used_for = c(
    "Infection PCA, sample scores and expression heatmaps",
    "Infection GSEA, enrichment curves and DE-module heatmaps",
    "Time-course PCA/scores and source matrix for first-passage BRAF PCA/scores",
    "BRAF GSEA, enrichment curves and DE-module heatmaps",
    "Infection sample matching and grouping",
    "Time-course matching, passage grouping and BRAF sample grouping"
  ),
  stringsAsFactors = FALSE
)

data.table::fwrite(
  input_role_qc,
  file.path(subdirs$qc, "INPUT_FILE_ROLES.csv")
)

# Input normalisation and plotting helpers ----
stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
warnf <- function(fmt, ...) warning(sprintf(fmt, ...), call. = FALSE)

safe_filename <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

normalise_id <- function(x) {
  x <- clean_character(x)
  x <- basename(x)
  x <- stringr::str_replace(x, "(?i)\\.(bam|sam|txt|csv)$", "")
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_")
  x <- stringr::str_replace_all(x, "_+", "_")
  x <- stringr::str_replace_all(x, "^_|_$", "")
  toupper(x)
}

column_key <- function(x) {
  toupper(gsub("[^A-Za-z0-9]", "", trimws(as.character(x))))
}

detect_col <- function(df, candidates, required = FALSE, label = "column") {
  nms <- names(df)
  if (length(nms) == 0) {
    if (required) stop("No columns are available while detecting ", label, ".")
    return(NA_character_)
  }

  exact <- candidates[candidates %in% nms]
  if (length(exact) > 0) return(exact[1])

  nms_key <- column_key(nms)
  candidate_key <- column_key(candidates)
  idx <- match(candidate_key, nms_key)
  idx <- idx[!is.na(idx)]

  if (length(idx) > 0) return(nms[idx[1]])

  if (required) {
    stop(
      "Could not detect ", label, ". Available columns:\n",
      paste(nms, collapse = "\n")
    )
  }

  NA_character_
}

normalise_type <- function(x) {
  raw <- clean_character(x)
  key <- toupper(gsub("[^A-Za-z0-9]", "", raw))

  out <- rep(NA_character_, length(key))

  out[key %in% c("SSL", "HP", "SSLHP", "SESSILESERRATEDLESION", "HYPERPLASTICPOLYP")] <- "SSL_HP"
  out[grepl("^HADJ", key) & grepl("SSL|HP|SERRATED", key)] <- "H_ADJ_SSL_HP"
  out[key %in% c(
    "HADJ", "HADJSSL", "HADJHP", "HADJSSLHP", "HADJCOL",
    "HADJACENT", "HADJACENTCOLON", "HEALTHYADJACENT",
    "HEALTHYADJACENTCOLON", "ADJACENT"
  )] <- "H_ADJ_SSL_HP"
  out[key %in% c("HCOL", "HEALTHYCOLON", "NORMALCOLON")] <- "H_COL"
  out[key %in% c("CA", "CONVENTIONALADENOMA")] <- "CA"
  out[key %in% c("STOMACH", "ANTRUM")] <- "STOMACH"

  out
}

normalise_braf <- function(x) {
  raw <- clean_character(x)
  key <- toupper(gsub("[^A-Za-z0-9]", "", raw))

  out <- rep(NA_character_, length(key))
  out[grepl("V600E|BRAFMUT|MUTANT", key)] <- "V600E"
  out[grepl("WILDTYPE", key) | key %in% c("WT", "BRAFWT")] <- "WT"
  out
}

normalise_infection <- function(x) {
  raw <- clean_character(x)
  key <- toupper(gsub("[^A-Za-z0-9]", "", raw))

  out <- rep(NA_character_, length(key))

  negative <- grepl(
    "NONINF|NONINFECT|UNINF|UNINFECT|CONTROL|CTRL|MOCK|UNTREAT",
    key
  ) | key %in% c("NI", "NEG", "NEGATIVE")

  positive <- (
    grepl("INFECT|IHE3034|ECOLI|ALTINFECTED", key) |
      grepl("INF$", key) |
      key %in% c("INF", "POS", "POSITIVE")
  ) & !negative

  out[negative] <- "Non-Infected"
  out[positive] <- "Infected"
  out
}

harmonise_infection_sample_id <- function(sample_name, infection_status) {
  original_id <- normalise_id(sample_name)
  corrected_id <- original_id
  status <- clean_character(infection_status)

  infected_idx <- !is.na(status) & status == "Infected" & !is.na(corrected_id)
  noninfected_idx <- !is.na(status) & status == "Non-Infected" & !is.na(corrected_id)

  if (any(infected_idx)) {
    corrected_id[infected_idx] <- stringr::str_replace_all(
      corrected_id[infected_idx],
      "(^|_)(NONINFECTED|NONINFECT|NONINF|UNINFECTED|UNINFECT|UNINF)(?=_|$)",
      "\\1INF"
    )
  }

  if (any(noninfected_idx)) {
    current <- corrected_id[noninfected_idx]

    has_negative_token <- stringr::str_detect(
      current,
      "(^|_)(NONINFECTED|NONINFECT|NONINF|UNINFECTED|UNINFECT|UNINF)(_|$)"
    )

    replace_idx <- !has_negative_token & stringr::str_detect(
      current,
      "(^|_)(INFECTED|INFECT|INF)(_|$)"
    )

    if (any(replace_idx)) {
      current[replace_idx] <- stringr::str_replace_all(
        current[replace_idx],
        "(^|_)(INFECTED|INFECT|INF)(?=_|$)",
        "\\1NONINF"
      )
    }

    corrected_id[noninfected_idx] <- current
  }

  corrected_id
}

parse_passage <- function(x) {
  raw <- clean_character(x)
  suppressWarnings(as.numeric(stringr::str_extract(raw, "[0-9]+(?:\\.[0-9]+)?")))
}

derive_line_id <- function(sample_id) {
  x <- normalise_id(sample_id)

  x <- stringr::str_replace_all(
    x,
    "(^|_)(NONINFECTED|NONINFECT|NONINF|UNINFECTED|UNINFECT|UNINF|INFECTED|INFECT|INF|CONTROL|CTRL|MOCK|UNTREATED|UNTREAT|IHE3034|ALT)($|_)",
    "_"
  )
  x <- stringr::str_replace_all(x, "(^|_)(PASSAGE|P)_?[0-9]+(?:\\.[0-9]+)?($|_)", "_")
  x <- stringr::str_replace_all(x, "(^|_)(REP|R)_?[0-9]+($|_)", "_")
  x <- stringr::str_replace_all(x, "_+", "_")
  x <- stringr::str_replace_all(x, "^_|_$", "")
  x[x == ""] <- normalise_id(sample_id)[x == ""]
  x
}

first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

format_p <- function(p) {
  if (length(p) == 0 || is.na(p)) return("p = NA")
  if (p < 1e-4) return(formatC(p, format = "e", digits = 1))
  paste0("p = ", formatC(p, format = "f", digits = 3))
}

p_to_plot_label <- function(p) {
  if (length(p) == 0 || is.na(p)) return("p = NA")
  if (p < 1e-4) return("****")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  format_p(p)
}

clean_pathway_label <- function(x) {
  x %>%
    stringr::str_remove("^HALLMARK_") %>%
    stringr::str_remove("^REACTOME_") %>%
    stringr::str_remove("^GOMF_") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_sentence()
}

format_signature_plot_title <- function(x, width = 34) {
  x <- as.character(x)
  out <- x
  target <- !is.na(out) & out == metaplasia_display_name
  out[target] <- stringr::str_wrap(out[target], width = width)
  out
}

nature_theme <- function(base_size = font_pt) {
  theme_classic(base_family = font_family, base_size = base_size) +
    theme(
      text = element_text(family = font_family, colour = "black"),
      plot.title = element_text(
        family = font_family,
        face = "bold",
        size = plot_title_pt,
        hjust = 0
      ),
      axis.title = element_text(
        family = font_family,
        size = axis_font_pt,
        colour = "black"
      ),
      axis.text = element_text(
        family = font_family,
        size = axis_font_pt,
        colour = "black"
      ),
      legend.title = element_text(
        family = font_family,
        size = legend_font_pt,
        colour = "black"
      ),
      legend.text = element_text(
        family = font_family,
        size = legend_font_pt,
        colour = "black"
      ),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = unit(1.2, "mm"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(4, 4, 4, 4, "mm")
    )
}

theme_set(nature_theme())

save_plot_set <- function(plot, directory, stem, width_mm, height_mm) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  stem <- safe_filename(stem)
  prefix <- file.path(directory, stem)
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
    dpi = export_dpi, device = grDevices::tiff,
    compression = "lzw", bg = "white"
  )

  ggsave(
    paste0(prefix, ".jpg"), plot,
    width = width_mm, height = height_mm, units = "mm",
    dpi = export_dpi, device = grDevices::jpeg,
    quality = 100, bg = "white"
  )
}

make_empty_plot <- function(title_text, message_text) {
  ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5, label = message_text,
      family = font_family, size = 2.2
    ) +
    xlim(0, 1) + ylim(0, 1) +
    labs(title = format_signature_plot_title(title_text)) +
    theme_void(base_family = font_family) +
    theme(
      plot.title = element_text(
        family = font_family,
        face = "bold",
        size = plot_title_pt,
        hjust = 0
      )
    )
}

extract_passage_from_sample_name <- function(x) {
  id <- normalise_id(x)
  hit <- stringr::str_match(
    id,
    "(?:^|_)(?:PASSAGE|P)_?([0-9]+(?:\\.[0-9]+)?)(?:_|$)"
  )[, 2]
  suppressWarnings(as.numeric(hit))
}

compose_metadata_id <- function(...) {
  parts <- lapply(list(...), clean_character)
  if (length(parts) == 0) return(character())

  n <- max(vapply(parts, length, integer(1)))
  parts <- lapply(parts, rep_len, length.out = n)
  missing_piece <- Reduce(`|`, lapply(parts, is.na))

  out <- do.call(paste, c(parts, sep = "__"))
  out[missing_piece] <- NA_character_
  out
}

unique_nonmissing <- function(x) {
  x <- clean_character(x)
  unique(x[!is.na(x)])
}

consensus_character <- function(x) {
  vals <- unique_nonmissing(x)
  if (length(vals) == 1L) vals[1] else NA_character_
}

consensus_numeric <- function(x) {
  vals <- unique(as.numeric(x[!is.na(x)]))
  if (length(vals) == 1L) vals[1] else NA_real_
}

first_nonmissing_scalar <- function(x) {
  valid <- if (is.character(x)) {
    which(!is.na(x) & trimws(x) != "")
  } else {
    which(!is.na(x))
  }

  if (length(valid) == 0L) return(x[NA_integer_][1])
  x[valid[1]]
}

# Metadata validation and infection-pair assignment ----
read_metadata <- function(
    file,
    dataset_label,
    metadata_role = c("infection", "timecourse")
) {
  metadata_role <- match.arg(metadata_role)

  meta <- data.table::fread(
    file,
    data.table = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "N/A", "na")
  )

  names(meta) <- trimws(names(meta))
  meta$.Metadata_input_row <- seq_len(nrow(meta))

  sample_col <- detect_col(
    meta,
    c("Name", "Sample", "Sample_ID", "SampleID", "Display", "sample_name"),
    required = TRUE,
    label = paste0(dataset_label, " sample-name column")
  )

  type_col <- detect_col(
    meta,
    c("Type", "Sample Type", "Sample_Type", "Group", "Tissue", "Origin"),
    required = FALSE
  )

  braf_col <- detect_col(
    meta,
    c("BRAF Status", "BRAF_Status", "BRAF", "BRAF.status", "BRAFStatus"),
    required = FALSE
  )

  passage_col <- detect_col(
    meta,
    c("Passage", "Passage Number", "Passage_Number", "P"),
    required = FALSE
  )

  infection_col <- detect_col(
    meta,
    c(
      "Infection Status", "Infection_Status", "Infection", "Treatment",
      "Condition", "Infected", "Status"
    ),
    required = FALSE
  )

  line_col <- detect_col(
    meta,
    c(
      "Organoid Line", "Organoid_Line", "Line", "Line_ID", "Organoid",
      "Patient", "Patient_ID", "PatientID", "Donor", "Donor_ID"
    ),
    required = FALSE
  )

  sample_original <- clean_character(meta[[sample_col]])
  sample_id_before_condition_harmonisation <- normalise_id(sample_original)

  type_raw <- if (!is.na(type_col)) meta[[type_col]] else rep(NA_character_, nrow(meta))
  braf_raw <- if (!is.na(braf_col)) meta[[braf_col]] else rep(NA_character_, nrow(meta))
  passage_raw <- if (!is.na(passage_col)) meta[[passage_col]] else rep(NA_character_, nrow(meta))
  infection_raw <- if (!is.na(infection_col)) meta[[infection_col]] else rep(NA_character_, nrow(meta))
  line_raw <- if (!is.na(line_col)) meta[[line_col]] else rep(NA_character_, nrow(meta))

  infection_from_column <- normalise_infection(infection_raw)
  infection_from_type <- normalise_infection(type_raw)
  infection_from_original_name <- normalise_infection(sample_original)

  sample_id <- if (metadata_role == "infection") {
    harmonise_infection_sample_id(
      sample_name = sample_original,
      infection_status = infection_from_column
    )
  } else {
    sample_id_before_condition_harmonisation
  }

  infection_from_name <- normalise_infection(sample_id)

  if (metadata_role == "infection") {
    correction_idx <- which(
      !is.na(sample_id_before_condition_harmonisation) &
        !is.na(sample_id) &
        sample_id_before_condition_harmonisation != sample_id
    )

    if (length(correction_idx) > 0L) {
      correction_audit <- data.frame(
        Metadata_input_row = meta$.Metadata_input_row[correction_idx],
        Sample_original = sample_original[correction_idx],
        Sample_ID_before = sample_id_before_condition_harmonisation[correction_idx],
        Infection_Status_raw = clean_character(infection_raw)[correction_idx],
        Infection_Status_standardised = infection_from_column[correction_idx],
        Sample_ID_after = sample_id[correction_idx],
        Reason = "Sample-name infection token corrected from explicit Infection_Status",
        stringsAsFactors = FALSE
      )

      correction_file <- file.path(
        subdirs$qc,
        paste0(safe_filename(dataset_label), "_metadata_sample_ID_corrections.csv")
      )
      data.table::fwrite(correction_audit, correction_file)

      message(
        "Corrected ", length(correction_idx),
        " infection metadata sample ID(s) using Infection_Status. Audit: ",
        correction_file
      )
    }
  }

  passage_from_column <- parse_passage(passage_raw)
  passage_from_name <- extract_passage_from_sample_name(sample_id)

  line_from_column <- normalise_id(line_raw)
  line_from_name <- derive_line_id(sample_id)

  if (metadata_role == "infection") {
    infection_status <- dplyr::coalesce(
      infection_from_column,
      infection_from_name,
      infection_from_type
    )
    passage_value <- dplyr::coalesce(passage_from_column, passage_from_name)
    line_id <- dplyr::coalesce(line_from_name, line_from_column)
  } else {
    infection_status <- dplyr::coalesce(
      infection_from_column,
      infection_from_type,
      infection_from_name
    )
    passage_value <- dplyr::coalesce(passage_from_column, passage_from_name)
    line_id <- dplyr::coalesce(line_from_column, line_from_name)
  }

  out <- meta %>%
    mutate(
      Dataset = dataset_label,
      Metadata_role = metadata_role,
      Sample_original = sample_original,
      Sample_ID_before_condition_harmonisation = sample_id_before_condition_harmonisation,
      Sample_ID = sample_id,
      Sample_ID_was_corrected = !is.na(sample_id_before_condition_harmonisation) &
        !is.na(sample_id) &
        sample_id_before_condition_harmonisation != sample_id,
      Type_raw = clean_character(type_raw),
      Type = normalise_type(type_raw),
      BRAF_raw = clean_character(braf_raw),
      BRAF_Status = normalise_braf(braf_raw),
      Passage_raw = clean_character(passage_raw),
      Passage = passage_value,
      Infection_Status = infection_status,
      Line_ID = line_id,
      .Infection_from_original_name = infection_from_original_name,
      .Infection_from_name = infection_from_name,
      .Passage_from_name = passage_from_name,
      .Line_from_name = line_from_name
    ) %>%
    filter(!is.na(Sample_ID), Sample_ID != "")

  original_sample_order <- unique(out$Sample_ID)

  if (anyDuplicated(out$Sample_ID)) {
    duplicate_rows <- out %>%
      group_by(Sample_ID) %>%
      filter(n() > 1L) %>%
      ungroup()

    duplicate_qc_file <- file.path(
      subdirs$qc,
      paste0(safe_filename(dataset_label), "_metadata_duplicate_rows_detected.csv")
    )
    data.table::fwrite(duplicate_rows, duplicate_qc_file)

    duplicate_groups <- split(duplicate_rows, duplicate_rows$Sample_ID)
    resolution_log <- list()
    unresolved_log <- list()

    resolved_duplicate_rows <- lapply(names(duplicate_groups), function(id) {
      one_id <- duplicate_groups[[id]]

      completeness_fields <- intersect(
        c(type_col, braf_col, passage_col, infection_col, line_col),
        names(one_id)
      )

      if (length(completeness_fields) > 0L) {
        completeness <- Reduce(
          `+`,
          lapply(
            one_id[completeness_fields],
            function(z) as.integer(!is.na(clean_character(z)))
          )
        )
        representative_index <- which.max(completeness)
      } else {
        representative_index <- 1L
      }

      representative <- one_id[representative_index, , drop = FALSE]

      field_values <- function(field) {
        vals <- unique_nonmissing(one_id[[field]])
        if (length(vals) == 0L) NA_character_ else paste(vals, collapse = " | ")
      }

      if (metadata_role == "infection") {
        name_condition <- unique_nonmissing(one_id$.Infection_from_name)
        row_condition <- unique_nonmissing(one_id$Infection_Status)

        resolved_condition <- if (length(name_condition) == 1L) {
          name_condition[1]
        } else if (length(row_condition) == 1L) {
          row_condition[1]
        } else {
          NA_character_
        }

        name_passage <- unique(as.numeric(one_id$.Passage_from_name[!is.na(one_id$.Passage_from_name)]))
        row_passage <- unique(as.numeric(one_id$Passage[!is.na(one_id$Passage)]))

        resolved_passage <- if (length(name_passage) == 1L) {
          name_passage[1]
        } else if (length(row_passage) == 1L) {
          row_passage[1]
        } else {
          NA_real_
        }

        name_line <- unique_nonmissing(one_id$.Line_from_name)
        row_line <- unique_nonmissing(one_id$Line_ID)

        resolved_line <- if (length(name_line) == 1L) {
          name_line[1]
        } else if (length(row_line) == 1L) {
          row_line[1]
        } else {
          NA_character_
        }

        type_values <- unique_nonmissing(one_id$Type)
        braf_values <- unique_nonmissing(one_id$BRAF_Status)

        representative$Infection_Status <- resolved_condition
        representative$Passage <- resolved_passage
        representative$Line_ID <- resolved_line
        representative$Type <- if (length(type_values) == 1L) type_values[1] else NA_character_
        representative$BRAF_Status <- if (length(braf_values) == 1L) braf_values[1] else NA_character_

        resolution_log[[length(resolution_log) + 1L]] <<- data.frame(
          Sample_ID = id,
          Field = c("Infection_Status", "Passage", "Line_ID", "Type", "BRAF_Status"),
          Values_detected = c(
            field_values("Infection_Status"),
            field_values("Passage"),
            field_values("Line_ID"),
            field_values("Type"),
            field_values("BRAF_Status")
          ),
          Resolution = c(
            ifelse(length(name_condition) == 1L, "Used condition encoded in Name", "Used unique standardised value"),
            ifelse(length(name_passage) == 1L, "Used passage encoded in Name", "Used unique standardised value or NA"),
            ifelse(length(name_line) == 1L, "Used line derived from Name", "Used unique standardised value or NA"),
            ifelse(length(type_values) <= 1L, "Used unique value", "Conflicting optional values set to NA"),
            ifelse(length(braf_values) <= 1L, "Used unique value", "Conflicting optional values set to NA")
          ),
          stringsAsFactors = FALSE
        )

        if (is.na(resolved_condition)) {
          unresolved_log[[length(unresolved_log) + 1L]] <<- data.frame(
            Sample_ID = id,
            Field = "Infection_Status",
            Values_detected = field_values("Infection_Status"),
            Reason = "Condition could not be resolved from Name or metadata",
            stringsAsFactors = FALSE
          )
        }

      } else {
        strict_fields <- c("Type", "BRAF_Status", "Passage", "Line_ID")
        strict_conflicts <- strict_fields[
          vapply(
            strict_fields,
            function(field) length(unique_nonmissing(one_id[[field]])) > 1L,
            logical(1)
          )
        ]

        if (length(strict_conflicts) > 0L) {
          for (field in strict_conflicts) {
            unresolved_log[[length(unresolved_log) + 1L]] <<- data.frame(
              Sample_ID = id,
              Field = field,
              Values_detected = field_values(field),
              Reason = "Conflicting time-course/BRAF metadata",
              stringsAsFactors = FALSE
            )
          }
        }

        representative$Type <- consensus_character(one_id$Type)
        representative$BRAF_Status <- consensus_character(one_id$BRAF_Status)
        representative$Passage <- consensus_numeric(one_id$Passage)
        representative$Line_ID <- consensus_character(one_id$Line_ID)
        representative$Infection_Status <- consensus_character(one_id$Infection_Status)

        resolution_log[[length(resolution_log) + 1L]] <<- data.frame(
          Sample_ID = id,
          Field = c("Type", "BRAF_Status", "Passage", "Line_ID"),
          Values_detected = vapply(
            c("Type", "BRAF_Status", "Passage", "Line_ID"),
            field_values,
            character(1)
          ),
          Resolution = ifelse(
            c("Type", "BRAF_Status", "Passage", "Line_ID") %in% strict_conflicts,
            "Unresolved conflict",
            "Used unique value"
          ),
          stringsAsFactors = FALSE
        )
      }

      representative
    })

    unresolved_tbl <- dplyr::bind_rows(unresolved_log)
    resolution_tbl <- dplyr::bind_rows(resolution_log)

    if (nrow(resolution_tbl) > 0L) {
      data.table::fwrite(
        resolution_tbl,
        file.path(
          subdirs$qc,
          paste0(safe_filename(dataset_label), "_metadata_duplicate_resolution.csv")
        )
      )
    }

    if (nrow(unresolved_tbl) > 0L) {
      conflict_file <- file.path(
        subdirs$qc,
        paste0(safe_filename(dataset_label), "_metadata_duplicate_unresolved_conflicts.csv")
      )
      data.table::fwrite(unresolved_tbl, conflict_file)

      stop(
        "Duplicate metadata rows could not be resolved in ", file, ".\n",
        "Duplicate rows: ", duplicate_qc_file, "\n",
        "Unresolved fields: ", conflict_file
      )
    }

    nonduplicate_rows <- out %>%
      group_by(Sample_ID) %>%
      filter(n() == 1L) %>%
      ungroup()

    out <- bind_rows(nonduplicate_rows, bind_rows(resolved_duplicate_rows)) %>%
      arrange(match(Sample_ID, original_sample_order))

    message(
      "Resolved duplicate annotation rows in ", basename(file),
      " using the ", metadata_role, "-specific metadata rules. Audit: ",
      file.path(
        subdirs$qc,
        paste0(safe_filename(dataset_label), "_metadata_duplicate_resolution.csv")
      )
    )
  }

  out <- out %>%
    mutate(
      Trajectory_ID = if (metadata_role == "infection") {
        compose_metadata_id(Line_ID)
      } else {
        compose_metadata_id(Line_ID, Type)
      },
      Pair_ID = if (metadata_role == "infection") {
        compose_metadata_id(Line_ID, Passage)
      } else {
        compose_metadata_id(Line_ID, Type, Passage)
      }
    ) %>%
    dplyr::select(
      -dplyr::any_of(c(
        ".Infection_from_original_name", ".Infection_from_name",
        ".Passage_from_name", ".Line_from_name"
      ))
    )

  if (anyDuplicated(out$Sample_ID)) {
    stop("Internal error: duplicate Sample_ID values remained after metadata resolution for ", dataset_label, ".")
  }

  if (metadata_role == "infection") {
    unresolved_condition <- out %>%
      filter(is.na(Infection_Status))

    if (nrow(unresolved_condition) > 0L) {
      data.table::fwrite(
        unresolved_condition,
        file.path(subdirs$qc, "Infection_metadata_unresolved_condition_rows.csv")
      )
      warning(
        nrow(unresolved_condition),
        " infection metadata row(s) have no resolvable infection condition and will be excluded."
      )
    }
  }

  data.table::fwrite(
    out,
    file.path(subdirs$qc, paste0(safe_filename(dataset_label), "_metadata_standardised.csv"))
  )

  out
}

infection_meta <- read_metadata(
  files$infection_meta,
  "Infection",
  metadata_role = "infection"
)

timecourse_meta <- read_metadata(
  files$timecourse_meta,
  "Timecourse",
  metadata_role = "timecourse"
)

# Expression import, scale checks and gene-symbol mapping ----
diagnose_expression_scale <- function(mat) {
  if (!is.matrix(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    stop("Expression matrix is empty during scale diagnosis.")
  }

  finite_values <- mat[is.finite(mat)]
  if (length(finite_values) == 0) stop("Expression matrix has no finite values.")

  column_sums <- colSums(mat, na.rm = TRUE)
  median_column_sum <- median(column_sums)
  column_sum_cv <- if (mean(column_sums) == 0) NA_real_ else sd(column_sums) / mean(column_sums)
  negative_fraction <- mean(finite_values < 0)
  integer_fraction <- mean(abs(finite_values - round(finite_values)) < 1e-8)
  q99 <- unname(quantile(finite_values, 0.99, na.rm = TRUE))
  max_value <- max(finite_values, na.rm = TRUE)
  median_value <- median(finite_values, na.rm = TRUE)
  zero_fraction <- mean(finite_values == 0)

  diagnosis <- "AMBIGUOUS"

  if (negative_fraction > 0) {
    diagnosis <- "LOG"
  } else if (
    median_column_sum >= 7e5 && median_column_sum <= 1.3e6 &&
      (is.na(column_sum_cv) || column_sum_cv <= 0.25)
  ) {
    diagnosis <- "CPM"
  } else if (max_value <= 50 && q99 <= 25) {
    diagnosis <- "LOG"
  } else if (integer_fraction >= 0.98 && max_value > 50) {
    diagnosis <- "RAW_COUNTS_OR_AMBIGUOUS"
  } else if (
    median_column_sum >= 3e5 && median_column_sum <= 2e6 && max_value > 50
  ) {
    diagnosis <- "CPM"
  }

  list(
    diagnosis = diagnosis,
    summary = data.frame(
      n_genes = nrow(mat),
      n_samples = ncol(mat),
      min_value = min(finite_values),
      median_value = median_value,
      q99_value = q99,
      max_value = max_value,
      zero_fraction = zero_fraction,
      negative_fraction = negative_fraction,
      integer_fraction = integer_fraction,
      median_column_sum = median_column_sum,
      column_sum_cv = column_sum_cv,
      automatic_scale_diagnosis = diagnosis,
      stringsAsFactors = FALSE
    ),
    column_sums = data.frame(
      Sample_column = colnames(mat),
      Column_sum = column_sums,
      stringsAsFactors = FALSE
    )
  )
}

map_gene_symbols <- function(gene_ids, supplied_symbols = NULL) {
  ids <- clean_character(gene_ids)
  ids_no_version <- toupper(stringr::str_remove(ids, "\\.[0-9]+$"))

  symbols <- if (!is.null(supplied_symbols)) {
    clean_character(supplied_symbols)
  } else {
    rep(NA_character_, length(ids))
  }
  symbols[toupper(symbols) %in% c("NA", "N/A", ".", "NULL")] <- NA_character_

  ensembl_mask <- !is.na(ids_no_version) & grepl("^ENSG[0-9]+$", ids_no_version)
  needs_mapping <- is.na(symbols) & ensembl_mask

  if (any(needs_mapping)) {
    mapped <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = unique(ids_no_version[needs_mapping]),
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )

    symbols[needs_mapping] <- unname(mapped[ids_no_version[needs_mapping]])
  }

  non_ensembl_fallback <- is.na(symbols) & !is.na(ids_no_version) & !ensembl_mask
  symbols[non_ensembl_fallback] <- ids_no_version[non_ensembl_fallback]

  symbols <- toupper(trimws(symbols))
  symbols[symbols == "" | symbols == "NA"] <- NA_character_
  symbols
}

collapse_rows_by_mean <- function(mat, groups) {
  keep <- !is.na(groups) & groups != ""
  mat <- mat[keep, , drop = FALSE]
  groups <- groups[keep]

  if (nrow(mat) == 0) stop("No mapped gene symbols remained after filtering.")

  summed <- rowsum(mat, group = groups, reorder = FALSE, na.rm = TRUE)
  counts <- table(groups)
  divisor <- as.numeric(counts[rownames(summed)])
  sweep(summed, 1, divisor, "/")
}

expression_annotation_keys <- column_key(c(
  "SYMBOL", "GENESYMBOL", "GENE_SYMBOL", "HGNC_SYMBOL",
  "GENENAME", "GENE_NAME", "GENEDESCRIPTION", "GENE_DESCRIPTION",
  "GENETYPE", "GENE_TYPE", "GENEBIOTYPE", "GENE_BIOTYPE", "BIOTYPE",
  "MAP", "MAPLOCATION", "MAP_LOCATION",
  "PATH", "PATHWAY", "PATH_ID", "PATHID",
  "DESCRIPTION", "ANNOTATION",
  "CHR", "CHROM", "CHROMOSOME",
  "START", "END", "STRAND", "LENGTH",
  "ENTREZID", "ENTREZ_ID", "ENTREZGENE", "ENTREZ_GENE",
  "HGNCID", "HGNC_ID", "REFSEQ", "UNIPROT"
))

is_expression_annotation_column <- function(x) {
  column_key(x) %in% expression_annotation_keys
}

numeric_fraction_nonmissing <- function(x) {
  original <- clean_character(x)
  keep <- !is.na(original)

  if (!any(keep)) return(0)

  converted <- suppressWarnings(as.numeric(original[keep]))
  mean(!is.na(converted))
}

read_expression_input <- function(file, dataset_label) {

  imported <- data.table::fread(
    file,
    header = FALSE,
    fill = TRUE,
    data.table = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "N/A", "na")
  )

  if (nrow(imported) < 2L || ncol(imported) < 3L) {
    stop(
      "Expression table has fewer than two rows or fewer than three columns: ",
      basename(file)
    )
  }

  nonempty_columns <- vapply(
    imported,
    function(x) any(!is.na(x) & trimws(as.character(x)) != ""),
    logical(1)
  )
  imported <- imported[, nonempty_columns, drop = FALSE]

  first_row <- clean_character(unlist(imported[1, , drop = TRUE], use.names = FALSE))
  body <- imported[-1, , drop = FALSE]

  valid_body_row <- !is.na(clean_character(body[[1]]))
  body <- body[valid_body_row, , drop = FALSE]

  if (nrow(body) == 0L) {
    stop("No gene rows remained after reading ", basename(file), ".")
  }

  n_body_columns <- ncol(body)
  header_positions <- which(!is.na(first_row))
  header_tokens <- first_row[header_positions]

  gene_header_keys <- column_key(c(
    "genes", "gene", "Gene", "gene_id", "Geneid", "ENSEMBL", "ID"
  ))
  first_token_is_gene_header <- length(header_tokens) > 0L &&
    column_key(header_tokens[1]) %in% gene_header_keys

  header_layout <- NA_character_
  assigned_names <- NULL

  if (
    length(header_tokens) == (n_body_columns - 1L) &&
      !first_token_is_gene_header
  ) {
    header_layout <- "SHIFTED_HEADERS_MISSING_GENE_HEADING"
    assigned_names <- c("Gene_ID", header_tokens)
  } else if (
    length(header_tokens) == n_body_columns &&
      first_token_is_gene_header
  ) {
    header_layout <- "CONVENTIONAL_HEADER"
    assigned_names <- first_row
  } else if (
    length(header_tokens) == (n_body_columns - 1L) &&
      first_token_is_gene_header &&
      is.na(first_row[n_body_columns])
  ) {
    stop(
      "The header in ", basename(file),
      " contains a gene heading but one fewer remaining column name than the data body. ",
      "Inspect row 1 before analysis."
    )
  } else {
    header_debug_file <- file.path(
      subdirs$qc,
      paste0(safe_filename(dataset_label), "_unresolved_expression_header.csv")
    )

    data.table::fwrite(
      data.frame(
        Column_position = seq_len(n_body_columns),
        Header_token = first_row[seq_len(n_body_columns)],
        First_body_value = vapply(
          body,
          function(x) as.character(x[1]),
          character(1)
        ),
        stringsAsFactors = FALSE
      ),
      header_debug_file
    )

    stop(
      "Could not resolve the expression-table header structure for ",
      basename(file), ". Found ", length(header_tokens),
      " non-empty header token(s) for ", n_body_columns,
      " body columns. Diagnostic: ", header_debug_file
    )
  }

  if (length(assigned_names) != n_body_columns || any(is.na(assigned_names))) {
    stop("Internal header-reconstruction error for ", basename(file), ".")
  }

  assigned_names <- trimws(as.character(assigned_names))

  if (any(assigned_names == "") || anyDuplicated(assigned_names)) {
    stop(
      "Blank or duplicated reconstructed column names were detected in ",
      basename(file), ":\n",
      paste(
        unique(assigned_names[duplicated(assigned_names) | assigned_names == ""]),
        collapse = "\n"
      )
    )
  }

  names(body) <- assigned_names

  gene_values <- clean_character(body[[1]])
  gene_numeric <- suppressWarnings(as.numeric(gene_values))
  gene_identifier_fraction <- mean(
    !is.na(gene_values) &
      (stringr::str_detect(
        toupper(gene_values),
        "^ENSG[0-9]+(?:\\.[0-9]+)?$"
      ) | is.na(gene_numeric))
  )

  if (!is.finite(gene_identifier_fraction) || gene_identifier_fraction < 0.90) {
    stop(
      "After header reconstruction, the first column of ", basename(file),
      " does not appear to contain gene identifiers."
    )
  }

  non_gene_columns <- names(body)[-1]
  numeric_fraction <- vapply(
    body[, non_gene_columns, drop = FALSE],
    numeric_fraction_nonmissing,
    numeric(1)
  )
  recognised_annotation <- is_expression_annotation_column(non_gene_columns)

  preliminary_role <- dplyr::case_when(
    recognised_annotation ~ "Recognised_gene_annotation",
    numeric_fraction == 1 ~ "Numeric_expression_candidate",
    TRUE ~ "Text_or_other_annotation"
  )

  header_audit <- data.frame(
    Dataset = dataset_label,
    File = basename(file),
    Header_layout_detected = header_layout,
    Body_columns = n_body_columns,
    Header_tokens_detected = length(header_tokens),
    Gene_column_assigned = assigned_names[1],
    Remaining_columns_assigned = length(non_gene_columns),
    Fully_numeric_remaining_columns = sum(numeric_fraction == 1),
    Recognised_annotation_columns = sum(recognised_annotation),
    Other_text_columns = sum(numeric_fraction < 1 & !recognised_annotation),
    Gene_identifier_fraction = gene_identifier_fraction,
    stringsAsFactors = FALSE
  )

  data.table::fwrite(
    header_audit,
    file.path(
      subdirs$qc,
      paste0(safe_filename(dataset_label), "_expression_header_reconstruction_QC.csv")
    )
  )

  original_token_positions <- if (
    header_layout == "SHIFTED_HEADERS_MISSING_GENE_HEADING"
  ) {
    c(NA_integer_, header_positions)
  } else {
    header_positions
  }

  original_tokens_assigned <- if (
    header_layout == "SHIFTED_HEADERS_MISSING_GENE_HEADING"
  ) {
    c(NA_character_, header_tokens)
  } else {
    header_tokens
  }

  mapping_role <- c("Gene_identifier", preliminary_role)
  mapping_numeric_fraction <- c(NA_real_, unname(numeric_fraction))

  data.table::fwrite(
    data.frame(
      Expression_body_column = seq_len(n_body_columns),
      Reconstructed_name = assigned_names,
      Preliminary_column_role = mapping_role,
      Numeric_fraction_among_nonmissing_values = mapping_numeric_fraction,
      Original_row1_token_position = original_token_positions,
      Original_row1_token = original_tokens_assigned,
      stringsAsFactors = FALSE
    ),
    file.path(
      subdirs$qc,
      paste0(safe_filename(dataset_label), "_expression_header_mapping.csv")
    )
  )

  message(
    "  Header layout: ", header_layout,
    "; restored Gene_ID + ", length(non_gene_columns), " remaining column(s): ",
    sum(numeric_fraction == 1), " fully numeric and ",
    sum(numeric_fraction < 1), " annotation/text column(s)."
  )

  body
}

read_expression_table <- function(file, metadata, dataset_label, scale_override = "AUTO") {
  message("\nReading ", dataset_label, " expression table: ", file)

  raw <- read_expression_input(
    file = file,
    dataset_label = dataset_label
  )

  gene_col <- detect_col(
    raw,
    c("Gene_ID", "genes", "gene", "Gene", "gene_id", "Geneid", "ENSEMBL", "ID"),
    required = FALSE,
    label = "gene identifier column"
  )
  if (is.na(gene_col)) gene_col <- names(raw)[1]

  non_gene_columns <- setdiff(names(raw), gene_col)
  non_gene_ids <- normalise_id(non_gene_columns)
  non_gene_numeric_fraction <- vapply(
    raw[, non_gene_columns, drop = FALSE],
    numeric_fraction_nonmissing,
    numeric(1)
  )

  present_in_metadata <- !is.na(non_gene_ids) &
    non_gene_ids %in% metadata$Sample_ID

  recognised_annotation <- is_expression_annotation_column(non_gene_columns)

  annotation_mask <- !present_in_metadata &
    (recognised_annotation | non_gene_numeric_fraction < 1)

  numeric_unmatched_mask <- !present_in_metadata &
    !annotation_mask &
    non_gene_numeric_fraction == 1

  matched_sample_columns <- non_gene_columns[present_in_metadata]
  matched_sample_ids <- non_gene_ids[present_in_metadata]

  if (any(is.na(matched_sample_ids)) || any(matched_sample_ids == "")) {
    stop("Blank matched sample header(s) remained in ", basename(file), ".")
  }

  if (anyDuplicated(matched_sample_ids)) {
    dup <- unique(matched_sample_ids[duplicated(matched_sample_ids)])
    stop(
      "Duplicate normalised matched sample IDs in ", basename(file), ":\n",
      paste(dup, collapse = "\n")
    )
  }

  column_role <- dplyr::case_when(
    present_in_metadata ~ "Matched_metadata_sample",
    annotation_mask & recognised_annotation ~ "Recognised_gene_annotation",
    annotation_mask ~ "Other_text_annotation",
    numeric_unmatched_mask ~ "Numeric_column_not_in_metadata",
    TRUE ~ "Unclassified"
  )

  expression_header_qc <- data.frame(
    Expression_column_number = match(non_gene_columns, names(raw)),
    Expression_header_original = non_gene_columns,
    Expression_header_normalised = non_gene_ids,
    Numeric_fraction_among_nonmissing_values = unname(non_gene_numeric_fraction),
    Recognised_annotation_name = recognised_annotation,
    Present_in_metadata = present_in_metadata,
    Column_role = column_role,
    stringsAsFactors = FALSE
  )

  metadata_match_qc <- metadata %>%
    dplyr::transmute(
      Sample_original,
      Sample_ID,
      Present_in_expression = Sample_ID %in% matched_sample_ids
    )

  expression_header_qc_file <- file.path(
    subdirs$qc,
    paste0(safe_filename(dataset_label), "_expression_sample_header_match_QC.csv")
  )
  metadata_match_qc_file <- file.path(
    subdirs$qc,
    paste0(safe_filename(dataset_label), "_metadata_sample_match_QC.csv")
  )

  data.table::fwrite(expression_header_qc, expression_header_qc_file)
  data.table::fwrite(metadata_match_qc, metadata_match_qc_file)

  data.table::fwrite(
    expression_header_qc %>%
      dplyr::filter(stringr::str_detect(Column_role, "annotation")),
    file.path(
      subdirs$qc,
      paste0(safe_filename(dataset_label), "_annotation_columns_detected.csv")
    )
  )

  matched_meta <- metadata %>%
    dplyr::filter(Sample_ID %in% matched_sample_ids)

  missing_meta_samples <- setdiff(metadata$Sample_ID, matched_sample_ids)
  numeric_columns_not_in_metadata <- non_gene_ids[numeric_unmatched_mask]

  data.table::fwrite(
    data.frame(Sample_ID = missing_meta_samples),
    file.path(
      subdirs$qc,
      paste0(safe_filename(dataset_label), "_metadata_samples_not_in_expression.csv")
    )
  )

  data.table::fwrite(
    data.frame(
      Expression_column_original = non_gene_columns[numeric_unmatched_mask],
      Expression_column_normalised = numeric_columns_not_in_metadata,
      stringsAsFactors = FALSE
    ),
    file.path(
      subdirs$qc,
      paste0(safe_filename(dataset_label), "_numeric_expression_columns_not_in_metadata.csv")
    )
  )

  if (nrow(matched_meta) < 2) {
    stop(
      "Fewer than two metadata samples matched ", basename(file), ".\n",
      "Expression-header audit: ", expression_header_qc_file, "\n",
      "Metadata matching audit: ", metadata_match_qc_file
    )
  }

  matched_column_idx <- match(matched_meta$Sample_ID, matched_sample_ids)
  sample_columns <- matched_sample_columns[matched_column_idx]

  order_in_table <- order(matched_column_idx)
  matched_meta <- matched_meta[order_in_table, , drop = FALSE]
  sample_columns <- sample_columns[order_in_table]

  expr_df <- raw[, sample_columns, drop = FALSE]
  numeric_expr <- lapply(expr_df, function(x) suppressWarnings(as.numeric(x)))
  numeric_expr <- as.data.frame(numeric_expr, check.names = FALSE)

  nonnumeric_problem <- vapply(
    seq_along(expr_df),
    function(j) {
      original <- clean_character(expr_df[[j]])
      converted <- numeric_expr[[j]]
      any(!is.na(original) & is.na(converted))
    },
    logical(1)
  )

  if (any(nonnumeric_problem)) {
    stop(
      "Non-numeric values were detected in metadata-matched sample columns of ",
      basename(file), ": ",
      paste(sample_columns[nonnumeric_problem], collapse = ", ")
    )
  }

  expr_uncollapsed <- as.matrix(numeric_expr)
  storage.mode(expr_uncollapsed) <- "double"
  colnames(expr_uncollapsed) <- matched_meta$Sample_ID

  if (any(!is.finite(expr_uncollapsed))) {
    stop("Non-finite values were detected in ", basename(file), ".")
  }

  scale_audit <- diagnose_expression_scale(expr_uncollapsed)
  scale_override <- toupper(scale_override)

  if (!scale_override %in% c("AUTO", "CPM", "LOG")) {
    stop("Invalid scale override for ", dataset_label, ": ", scale_override)
  }

  scale_used <- if (scale_override == "AUTO") {
    scale_audit$diagnosis
  } else {
    scale_override
  }

  scale_audit$summary$dataset <- dataset_label
  scale_audit$summary$file <- basename(file)
  scale_audit$summary$scale_override <- scale_override
  scale_audit$summary$scale_used <- scale_used

  data.table::fwrite(
    scale_audit$summary,
    file.path(
      subdirs$qc,
      paste0(safe_filename(dataset_label), "_expression_scale_QC.csv")
    )
  )
  data.table::fwrite(
    scale_audit$column_sums,
    file.path(
      subdirs$qc,
      paste0(safe_filename(dataset_label), "_expression_column_sums_QC.csv")
    )
  )

  if (scale_used %in% c("AMBIGUOUS", "RAW_COUNTS_OR_AMBIGUOUS")) {
    stop(
      "The expression scale for ", basename(file), " was diagnosed as '",
      scale_used, "'. The PCA was not run. Inspect the QC file and set the ",
      "corresponding input_scale_override to CPM or LOG only when justified."
    )
  }

  symbol_col <- detect_col(
    raw,
    c("SYMBOL", "Symbol", "symbol", "gene_symbol", "GeneSymbol"),
    required = FALSE,
    label = "gene symbol column"
  )

  gene_ids <- raw[[gene_col]]
  supplied_symbols <- if (!is.na(symbol_col)) raw[[symbol_col]] else NULL
  gene_symbols <- map_gene_symbols(gene_ids, supplied_symbols)

  mapping_qc <- data.frame(
    dataset = dataset_label,
    input_rows = length(gene_symbols),
    mapped_symbol_rows = sum(!is.na(gene_symbols)),
    mapped_symbol_fraction = mean(!is.na(gene_symbols)),
    unique_mapped_symbols = length(unique(gene_symbols[!is.na(gene_symbols)])),
    gene_id_column = gene_col,
    supplied_symbol_column = ifelse(is.na(symbol_col), NA_character_, symbol_col),
    annotation_columns_detected = paste(
      non_gene_columns[annotation_mask],
      collapse = ";"
    ),
    matched_sample_columns = length(sample_columns),
    stringsAsFactors = FALSE
  )
  data.table::fwrite(
    mapping_qc,
    file.path(
      subdirs$qc,
      paste0(safe_filename(dataset_label), "_gene_mapping_QC.csv")
    )
  )

  expr_collapsed <- collapse_rows_by_mean(expr_uncollapsed, gene_symbols)

  log_expr <- if (scale_used == "CPM") {
    log2(expr_collapsed + cpm_pseudocount)
  } else {
    expr_collapsed
  }

  if (any(!is.finite(log_expr))) {
    stop("Non-finite transformed values were generated for ", dataset_label, ".")
  }

  matched_meta <- matched_meta[
    match(colnames(log_expr), matched_meta$Sample_ID),
    ,
    drop = FALSE
  ]
  stopifnot(identical(colnames(log_expr), matched_meta$Sample_ID))

  data.table::fwrite(
    matched_meta,
    file.path(
      subdirs$qc,
      paste0(safe_filename(dataset_label), "_metadata_aligned_to_expression.csv")
    )
  )

  message(
    "  Metadata-matched sample columns used: ", length(sample_columns),
    "; annotation/text columns excluded: ", sum(annotation_mask), "."
  )

  list(
    dataset = dataset_label,
    source_file = file,
    scale_used = scale_used,
    expression_for_ranking = expr_collapsed,
    log_expr = log_expr,
    meta = matched_meta,
    scale_qc = scale_audit$summary,
    mapping_qc = mapping_qc
  )
}

infection_obj <- read_expression_table(
  files$infection_cpm,
  infection_meta,
  "Infection",
  input_scale_override[["infection"]]
)

timecourse_obj <- read_expression_table(
  files$timecourse_cpm,
  timecourse_meta,
  "Timecourse",
  input_scale_override[["timecourse"]]
)

# Initial-passage and strict-P1 cohorts ----
subset_expression_object <- function(
    expression_object,
    selected_metadata,
    dataset_label,
    selection_note
) {
  sample_ids <- unique(selected_metadata$Sample_ID)

  if (length(sample_ids) < 2L) {
    stop("Fewer than two samples were selected for ", dataset_label, ".")
  }

  missing_from_expression <- setdiff(
    sample_ids,
    colnames(expression_object$log_expr)
  )

  if (length(missing_from_expression) > 0L) {
    stop(
      "Selected ", dataset_label,
      " sample(s) were absent from the source CPM object:\n",
      paste(missing_from_expression, collapse = "\n")
    )
  }

  selected_metadata <- selected_metadata[
    match(sample_ids, selected_metadata$Sample_ID),
    ,
    drop = FALSE
  ]

  log_expr_subset <- expression_object$log_expr[
    ,
    sample_ids,
    drop = FALSE
  ]

  ranking_expr_subset <- expression_object$expression_for_ranking[
    ,
    sample_ids,
    drop = FALSE
  ]

  stopifnot(identical(colnames(log_expr_subset), selected_metadata$Sample_ID))
  stopifnot(identical(colnames(ranking_expr_subset), selected_metadata$Sample_ID))

  derived_scale_qc <- expression_object$scale_qc
  derived_scale_qc$dataset <- dataset_label
  derived_scale_qc$file <- basename(expression_object$source_file)
  if ("n_samples" %in% names(derived_scale_qc)) {
    derived_scale_qc$n_samples <- length(sample_ids)
  }
  derived_scale_qc$derived_subset <- TRUE
  derived_scale_qc$selection_note <- selection_note

  derived_mapping_qc <- expression_object$mapping_qc
  derived_mapping_qc$dataset <- dataset_label
  if ("matched_sample_columns" %in% names(derived_mapping_qc)) {
    derived_mapping_qc$matched_sample_columns <- length(sample_ids)
  }

  list(
    dataset = dataset_label,
    source_file = expression_object$source_file,
    scale_used = expression_object$scale_used,
    expression_for_ranking = ranking_expr_subset,
    log_expr = log_expr_subset,
    meta = selected_metadata,
    scale_qc = derived_scale_qc,
    mapping_qc = derived_mapping_qc,
    derived_from_dataset = expression_object$dataset,
    selection_note = selection_note
  )
}

braf_selection_audit <- timecourse_obj$meta %>%
  dplyr::filter(
    Type %in% c("SSL_HP", "H_ADJ_SSL_HP"),
    is.finite(Passage)
  ) %>%
  dplyr::arrange(Trajectory_ID, Passage, Sample_ID)

if (nrow(braf_selection_audit) < 2L) {
  stop(
    "Fewer than two SSL_HP/H_ADJ_SSL_HP samples with valid passage metadata ",
    "were available in CGG028_cpm_table.txt."
  )
}

if (any(is.na(braf_selection_audit$Trajectory_ID) |
        braf_selection_audit$Trajectory_ID == "")) {
  stop(
    "At least one candidate BRAF sample has no resolvable Trajectory_ID. ",
    "Inspect Timecourse_metadata_standardised.csv."
  )
}

braf_status_conflicts <- braf_selection_audit %>%
  dplyr::group_by(Trajectory_ID) %>%
  dplyr::summarise(
    N_distinct_nonmissing_BRAF = dplyr::n_distinct(
      BRAF_Status[!is.na(BRAF_Status)]
    ),
    BRAF_values = paste(
      sort(unique(BRAF_Status[!is.na(BRAF_Status)])),
      collapse = " | "
    ),
    .groups = "drop"
  ) %>%
  dplyr::filter(N_distinct_nonmissing_BRAF > 1L)

if (nrow(braf_status_conflicts) > 0L) {
  conflict_file <- file.path(
    subdirs$qc,
    "BRAF_status_conflicts_within_timecourse_trajectory.csv"
  )
  data.table::fwrite(braf_status_conflicts, conflict_file)
  stop(
    "Conflicting BRAF annotations were found within at least one organoid ",
    "trajectory. Inspect: ", conflict_file
  )
}

braf_selection_audit <- braf_selection_audit %>%
  dplyr::group_by(Trajectory_ID) %>%
  dplyr::arrange(Passage, Sample_ID, .by_group = TRUE) %>%
  dplyr::mutate(
    BRAF_Status_before_trajectory_fill = BRAF_Status,
    Trajectory_BRAF_Status = consensus_character(BRAF_Status),
    BRAF_Status_filled_from_trajectory =
      is.na(BRAF_Status) & !is.na(Trajectory_BRAF_Status),
    BRAF_Status = dplyr::coalesce(BRAF_Status, Trajectory_BRAF_Status),
    Earliest_available_passage = min(Passage, na.rm = TRUE),
    Has_P1 = any(Passage == 1),
    Order_within_trajectory = dplyr::row_number(),
    Selected_for_BRAF_initial_passage = Order_within_trajectory == 1L,
    Selected_for_BRAF_strict_P1 = Passage == 1,
    Initial_selection_rule = paste0(
      "One earliest available sample per trajectory; P1 preferred and P2 used only when P1 is unavailable. ",
      "Selected passage = P", Earliest_available_passage
    ),
    Strict_P1_selection_rule = ifelse(
      Passage == 1,
      "Selected for strict-P1 sensitivity",
      "Not selected for strict-P1 sensitivity"
    )
  ) %>%
  dplyr::ungroup()

initial_selection_counts <- braf_selection_audit %>%
  dplyr::group_by(Trajectory_ID) %>%
  dplyr::summarise(
    N_selected_initial = sum(Selected_for_BRAF_initial_passage),
    .groups = "drop"
  )

if (any(initial_selection_counts$N_selected_initial != 1L)) {
  bad_file <- file.path(
    subdirs$qc,
    "BRAF_initial_passage_invalid_selection_counts.csv"
  )
  data.table::fwrite(initial_selection_counts, bad_file)
  stop(
    "The initial-passage rule did not select exactly one sample per trajectory. ",
    "Inspect: ", bad_file
  )
}

p1_duplicates <- braf_selection_audit %>%
  dplyr::filter(Passage == 1) %>%
  dplyr::count(Trajectory_ID, name = "N_P1_samples") %>%
  dplyr::filter(N_P1_samples > 1L)

if (nrow(p1_duplicates) > 0L) {
  duplicate_file <- file.path(
    subdirs$qc,
    "BRAF_P1_duplicate_samples_within_trajectory.csv"
  )
  data.table::fwrite(p1_duplicates, duplicate_file)
  stop(
    "More than one P1 sample was found within at least one organoid trajectory. ",
    "Resolve the duplicate P1 samples before analysis: ", duplicate_file
  )
}

braf_initial_meta <- braf_selection_audit %>%
  dplyr::filter(Selected_for_BRAF_initial_passage) %>%
  dplyr::arrange(Type, Line_ID, Passage, Sample_ID)

if (any(!(braf_initial_meta$Passage %in% c(1, 2)))) {
  unexpected_file <- file.path(
    subdirs$qc,
    "BRAF_initial_passage_above_P2.csv"
  )
  data.table::fwrite(
    braf_initial_meta %>% dplyr::filter(!(Passage %in% c(1, 2))),
    unexpected_file
  )
  stop(
    "At least one trajectory has an initial available passage later than P2. ",
    "The approved analysis permits only P1 or P2. Inspect: ", unexpected_file
  )
}

braf_p1_meta <- braf_selection_audit %>%
  dplyr::filter(Selected_for_BRAF_strict_P1) %>%
  dplyr::arrange(Type, Line_ID, Passage, Sample_ID)

if (nrow(braf_initial_meta) < 2L) {
  stop("Fewer than two initial-passage samples remained for the BRAF analysis.")
}
if (nrow(braf_p1_meta) < 2L) {
  stop("Fewer than two strict-P1 samples remained for the BRAF sensitivity analysis.")
}
if (any(braf_p1_meta$Passage != 1)) {
  stop("Internal error: a non-P1 sample entered the strict-P1 BRAF object.")
}

braf_trajectories_without_p1 <- braf_selection_audit %>%
  dplyr::group_by(Trajectory_ID, Line_ID, Type, BRAF_Status) %>%
  dplyr::summarise(
    Available_passages = paste(sort(unique(Passage)), collapse = ";"),
    Earliest_available_passage = min(Passage),
    .groups = "drop"
  ) %>%
  dplyr::filter(Earliest_available_passage != 1)

data.table::fwrite(
  braf_selection_audit,
  file.path(subdirs$qc, "BRAF_initial_and_strict_P1_selection_audit.csv")
)
data.table::fwrite(
  braf_initial_meta,
  file.path(subdirs$qc, "BRAF_initial_passage_P1_or_P2_samples_used.csv")
)
data.table::fwrite(
  braf_p1_meta,
  file.path(subdirs$qc, "BRAF_strict_P1_samples_used.csv")
)
data.table::fwrite(
  braf_trajectories_without_p1,
  file.path(subdirs$qc, "BRAF_strict_P1_trajectories_excluded_no_P1.csv")
)
data.table::fwrite(
  braf_initial_meta %>% dplyr::filter(Passage == 2),
  file.path(subdirs$qc, "BRAF_P2_samples_included_in_initial_analysis.csv")
)

summarise_braf_sample_set <- function(meta, analysis_set) {
  meta %>%
    dplyr::mutate(
      Analysis_set = analysis_set,
      Analysis_group = dplyr::case_when(
        Type == "SSL_HP" & BRAF_Status == "WT" ~ "SSL_HP WT",
        Type == "SSL_HP" & BRAF_Status == "V600E" ~ "SSL_HP V600E",
        Type == "SSL_HP" ~ "SSL_HP BRAF missing",
        Type == "H_ADJ_SSL_HP" ~ "H_ADJ_SSL_HP",
        TRUE ~ "Other"
      )
    ) %>%
    dplyr::count(Analysis_set, Analysis_group, Passage, name = "N") %>%
    dplyr::arrange(Analysis_set, Analysis_group, Passage)
}

braf_sample_set_summary <- dplyr::bind_rows(
  summarise_braf_sample_set(
    braf_initial_meta,
    "Initial passage: P1, or P2 only when P1 unavailable"
  ),
  summarise_braf_sample_set(
    braf_p1_meta,
    "Strict P1 only"
  )
)

data.table::fwrite(
  braf_sample_set_summary,
  file.path(subdirs$qc, "BRAF_initial_vs_strict_P1_sample_counts.csv")
)

braf_obj <- subset_expression_object(
  expression_object = timecourse_obj,
  selected_metadata = braf_initial_meta,
  dataset_label = "BRAF_initial_passage_P1_or_P2",
  selection_note = paste(
    "Derived from CGG028_cpm_table.txt.",
    "One sample per organoid trajectory was retained: P1 where available and P2 only when P1 was unavailable.",
    "The separate CGG028_First_Passage_SSL_HP_V600E_WT.txt file is used only as the BRAF DE table."
  )
)

braf_p1_obj <- subset_expression_object(
  expression_object = timecourse_obj,
  selected_metadata = braf_p1_meta,
  dataset_label = "BRAF_strict_P1",
  selection_note = paste(
    "Derived from CGG028_cpm_table.txt.",
    "Only strict Passage 1 samples were retained; trajectories without P1 were excluded.",
    "This object is a passage-sensitivity analysis for the initial-passage P1/P2 analysis."
  )
)

message(
  "\nCreated initial-passage BRAF expression object from ",
  basename(files$timecourse_cpm), ": ",
  ncol(braf_obj$log_expr), " sample(s) across ",
  dplyr::n_distinct(braf_obj$meta$Trajectory_ID), " trajectory/trajectories; ",
  sum(braf_obj$meta$Passage == 2), " P2 sample(s) included because P1 was unavailable."
)

message(
  "Created strict-P1 BRAF sensitivity object: ",
  ncol(braf_p1_obj$log_expr), " sample(s) across ",
  dplyr::n_distinct(braf_p1_obj$meta$Trajectory_ID), " trajectory/trajectories."
)

# Principal component analysis ----
compute_pca <- function(log_expr, meta, dataset_label, n_variable = pca_n_variable_genes) {
  if (ncol(log_expr) != nrow(meta)) {
    stop("PCA expression/metadata dimension mismatch for ", dataset_label, ".")
  }

  log_expr <- log_expr[rowSums(!is.finite(log_expr)) == 0, , drop = FALSE]
  gene_variance <- apply(log_expr, 1, var)
  keep <- is.finite(gene_variance) & gene_variance > 0
  log_expr <- log_expr[keep, , drop = FALSE]
  gene_variance <- gene_variance[keep]

  if (length(gene_variance) < 2) {
    stop("Too few variable genes remained for PCA: ", dataset_label)
  }

  selected_idx <- head(order(gene_variance, decreasing = TRUE), n_variable)
  pca_input <- log_expr[selected_idx, , drop = FALSE]

  pca <- prcomp(t(pca_input), center = TRUE, scale. = FALSE)
  variance_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)

  scores <- as.data.frame(pca$x, check.names = FALSE) %>%
    tibble::rownames_to_column("Sample_ID") %>%
    left_join(meta, by = "Sample_ID")

  if (anyNA(scores$Dataset)) {
    stop("PCA scores could not be fully joined to metadata for ", dataset_label, ".")
  }

  stem <- safe_filename(dataset_label)

  data.table::fwrite(
    scores,
    file.path(subdirs$tables, paste0(stem, "_PCA_scores.csv"))
  )

  data.table::fwrite(
    data.frame(
      PC = paste0("PC", seq_along(variance_explained)),
      Variance_explained_percent = round(variance_explained, 4)
    ),
    file.path(subdirs$tables, paste0(stem, "_PCA_variance_explained.csv"))
  )

  data.table::fwrite(
    data.frame(
      Rank = seq_along(selected_idx),
      Gene = rownames(pca_input),
      Mean_expression = rowMeans(pca_input),
      Variance_expression = unname(gene_variance[selected_idx])
    ),
    file.path(subdirs$tables, paste0(stem, "_PCA_top_variable_genes.csv"))
  )

  list(
    pca = pca,
    scores = scores,
    variance_explained = variance_explained,
    input = pca_input
  )
}

pca_infection_meta <- infection_obj$meta %>%
  filter(Infection_Status %in% c("Non-Infected", "Infected"))

if (nrow(pca_infection_meta) < 3 || n_distinct(pca_infection_meta$Infection_Status) < 2) {
  stop(
    "Infection metadata must resolve at least three samples across both ",
    "Non-Infected and Infected groups."
  )
}

pca_infection_expr <- infection_obj$log_expr[, pca_infection_meta$Sample_ID, drop = FALSE]
pca_infection <- compute_pca(pca_infection_expr, pca_infection_meta, "Infection_PCA")

pca_infection$scores <- pca_infection$scores %>%
  dplyr::mutate(
    Passage_label = dplyr::if_else(
      is.finite(Passage),
      paste0(
        "P",
        format(Passage, trim = TRUE, scientific = FALSE)
      ),
      NA_character_
    )
  )

if (all(is.na(pca_infection$scores$Passage_label))) {
  warning(
    "No passage values were available for the infection PCA. ",
    "The PCA will be exported without passage labels."
  )
}

data.table::fwrite(
  pca_infection$scores %>%
    dplyr::select(
      Sample_ID,
      Infection_Status,
      Passage,
      Passage_label,
      PC1,
      PC2,
      dplyr::everything()
    ),
  file.path(
    subdirs$tables,
    "Infection_PCA_scores_with_passage_labels.csv"
  )
)

p_infection_pca <- ggplot2::ggplot(
  pca_infection$scores,
  ggplot2::aes(x = PC1, y = PC2, fill = Infection_Status)
) +
  ggplot2::geom_point(
    shape = 21,
    size = point_size,
    stroke = point_stroke,
    colour = "black",
    alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = Passage_label),
    family = font_family,
    size = 5.5 / 2.845276,
    colour = "black",
    box.padding = 0.32,
    point.padding = 0.22,
    min.segment.length = 0,
    segment.colour = "grey55",
    segment.linewidth = 0.25,
    max.overlaps = Inf,
    seed = 1,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  ggplot2::scale_fill_manual(
    values = infection_cols,
    breaks = c("Non-Infected", "Infected"),
    drop = FALSE
  ) +
  ggplot2::labs(
    x = sprintf("PC1 (%.1f%%)", pca_infection$variance_explained[1]),
    y = sprintf("PC2 (%.1f%%)", pca_infection$variance_explained[2]),
    fill = NULL
  ) +
  nature_theme() +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme(
    legend.position = "right",
    plot.margin = ggplot2::margin(5, 8, 5, 5, unit = "mm")
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      override.aes = list(
        shape = 21,
        size = point_size,
        colour = "black",
        alpha = 1
      )
    )
  )

save_plot_set(
  p_infection_pca,
  subdirs$pca,
  "PCA_1_Infected_vs_Non_Infected",
  90,
  70
)

save_plot_set(
  p_infection_pca,
  subdirs$pca,
  "PCA_1_Infected_vs_Non_Infected_with_Passage_Labels",
  90,
  70
)

pca_time_meta <- timecourse_obj$meta %>%
  filter(
    Type %in% c("SSL_HP", "H_ADJ_SSL_HP"),
    is.finite(Passage)
  ) %>%
  arrange(Type, Line_ID, Passage, Sample_ID)

if (nrow(pca_time_meta) < 3 || n_distinct(pca_time_meta$Passage) < 2) {
  stop("Time-course PCA requires at least two distinct passages.")
}

pca_time_expr <- timecourse_obj$log_expr[, pca_time_meta$Sample_ID, drop = FALSE]
pca_time <- compute_pca(pca_time_expr, pca_time_meta, "Timecourse_PCA")

passage_levels <- sort(unique(pca_time$scores$Passage))
passage_labels <- paste0("P", passage_levels)
pca_time$scores$Passage_factor <- factor(
  paste0("P", pca_time$scores$Passage),
  levels = passage_labels
)

filled_shapes <- c(21, 22, 24, 23, 25, 8, 4, 3, 0, 1, 2, 5, 6)
if (length(passage_levels) > length(filled_shapes)) {
  stop("Too many distinct passages for the configured PCA shape palette.")
}
passage_shape_values <- setNames(filled_shapes[seq_along(passage_levels)], passage_labels)

trajectory_lines <- pca_time$scores %>%
  group_by(Trajectory_ID) %>%
  filter(n_distinct(Passage) >= 2) %>%
  arrange(Passage, .by_group = TRUE) %>%
  ungroup()

p_time_pca <- ggplot(pca_time$scores, aes(x = PC1, y = PC2)) +
  geom_path(
    data = trajectory_lines,
    aes(group = Trajectory_ID),
    colour = "grey70",
    linewidth = 0.35,
    alpha = 0.8
  ) +
  geom_point(
    aes(fill = Type, shape = Passage_factor),
    size = point_size,
    stroke = point_stroke,
    colour = "black",
    alpha = 0.95
  ) +
  scale_fill_manual(
    values = paper_cols[c("SSL_HP", "H_ADJ_SSL_HP")],
    breaks = c("SSL_HP", "H_ADJ_SSL_HP"),
    drop = FALSE
  ) +
  scale_shape_manual(values = passage_shape_values, breaks = passage_labels, drop = FALSE) +
  labs(
    x = sprintf("PC1 (%.1f%%)", pca_time$variance_explained[1]),
    y = sprintf("PC2 (%.1f%%)", pca_time$variance_explained[2]),
    fill = "Sample type",
    shape = "Passage"
  ) +
  nature_theme() +
  theme(legend.position = "right") +
  guides(
    fill = guide_legend(order = 1, override.aes = list(shape = 21, size = point_size)),
    shape = guide_legend(order = 2, override.aes = list(fill = "white", size = point_size))
  )

save_plot_set(p_time_pca, subdirs$pca, "PCA_2_Timecourse_SSL_HP_and_H_ADJ", 100, 70)

pca_braf_meta <- braf_obj$meta %>%
  dplyr::filter(
    (Type == "SSL_HP" & BRAF_Status %in% c("WT", "V600E")) |
      Type == "H_ADJ_SSL_HP"
  ) %>%
  dplyr::mutate(
    PCA_Group = dplyr::case_when(
      Type == "SSL_HP" & BRAF_Status == "WT" ~ "SSL_HP WT",
      Type == "SSL_HP" & BRAF_Status == "V600E" ~ "SSL_HP V600E",
      Type == "H_ADJ_SSL_HP" ~ "H_ADJ_SSL_HP",
      TRUE ~ NA_character_
    ),
    Passage_factor = factor(
      paste0("P", Passage),
      levels = paste0("P", sort(unique(Passage)))
    )
  ) %>%
  dplyr::filter(!is.na(PCA_Group))

ssl_braf_groups_present <- pca_braf_meta %>%
  dplyr::filter(Type == "SSL_HP") %>%
  dplyr::distinct(BRAF_Status) %>%
  dplyr::pull(BRAF_Status)

if (!all(c("WT", "V600E") %in% ssl_braf_groups_present)) {
  stop("BRAF PCA requires both WT and V600E SSL/HP samples.")
}

if (!any(pca_braf_meta$Type == "H_ADJ_SSL_HP")) {
  warning(
    "No H_ADJ samples were selected from the time-course CPM table for the ",
    "BRAF PCA; the PCA will contain SSL/HP samples only."
  )
}

pca_braf_expr <- braf_obj$log_expr[, pca_braf_meta$Sample_ID, drop = FALSE]
pca_braf <- compute_pca(
  pca_braf_expr,
  pca_braf_meta,
  "BRAF_PCA_initial_passage_P1_or_P2"
)

braf_pca_levels <- c("SSL_HP WT", "SSL_HP V600E", "H_ADJ_SSL_HP")
braf_pca_levels <- braf_pca_levels[
  braf_pca_levels %in% unique(pca_braf$scores$PCA_Group)
]
pca_braf$scores$PCA_Group <- factor(
  pca_braf$scores$PCA_Group,
  levels = braf_pca_levels
)

braf_pca_cols <- c(
  "SSL_HP WT" = braf_cols[["WT"]],
  "SSL_HP V600E" = braf_cols[["V600E"]],
  "H_ADJ_SSL_HP" = paper_cols[["H_ADJ_SSL_HP"]]
)

braf_passage_levels <- levels(droplevels(pca_braf$scores$Passage_factor))
braf_passage_shapes <- c("P1" = 21, "P2" = 24)
if (!all(braf_passage_levels %in% names(braf_passage_shapes))) {
  stop(
    "The initial-passage BRAF PCA contains a passage other than P1/P2: ",
    paste(setdiff(braf_passage_levels, names(braf_passage_shapes)), collapse = ", ")
  )
}

p_braf_pca <- ggplot2::ggplot(
  pca_braf$scores,
  ggplot2::aes(x = PC1, y = PC2, fill = PCA_Group, shape = Passage_factor)
) +
  ggplot2::geom_point(
    size = point_size,
    stroke = point_stroke,
    colour = "black",
    alpha = 0.95
  ) +
  ggplot2::scale_fill_manual(
    values = braf_pca_cols[braf_pca_levels],
    breaks = braf_pca_levels,
    drop = FALSE
  ) +
  ggplot2::scale_shape_manual(
    values = braf_passage_shapes[braf_passage_levels],
    breaks = braf_passage_levels,
    drop = FALSE
  ) +
  ggplot2::labs(
    x = sprintf("PC1 (%.1f%%)", pca_braf$variance_explained[1]),
    y = sprintf("PC2 (%.1f%%)", pca_braf$variance_explained[2]),
    fill = NULL,
    shape = "Initial passage"
  ) +
  nature_theme() +
  ggplot2::theme(legend.position = "right") +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      order = 1,
      override.aes = list(shape = 21, size = point_size, colour = "black", alpha = 1)
    ),
    shape = ggplot2::guide_legend(
      order = 2,
      override.aes = list(fill = "white", size = point_size, colour = "black", alpha = 1)
    )
  )

save_plot_set(
  p_braf_pca,
  subdirs$pca,
  "PCA_3_BRAF_WT_vs_V600E_Initial_P1_or_P2_with_H_ADJ",
  105,
  72
)

pca_combined <- (p_infection_pca + theme(legend.position = "bottom")) |
  (p_time_pca + theme(legend.position = "bottom")) |
  (p_braf_pca + theme(legend.position = "bottom"))

pca_combined <- pca_combined +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(
        family = font_family,
        face = "bold",
        size = panel_tag_pt
      )
    )
  )

save_plot_set(pca_combined, subdirs$combined, "Combined_three_PCA_plots", 183, 78)

# Pathway definitions and the fixed gastric signature ----
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

standardise_msig <- function(msig_df) {
  gene_col <- detect_col(
    msig_df,
    c("gene_symbol", "human_gene_symbol", "db_gene_symbol", "symbol"),
    required = TRUE,
    label = "MSigDB gene-symbol column"
  )
  pathway_col <- detect_col(
    msig_df,
    c("gs_name", "geneset", "gene_set_name", "gs_exact_source"),
    required = TRUE,
    label = "MSigDB pathway-name column"
  )
  subcollection_col <- detect_col(
    msig_df,
    c("gs_subcollection", "gs_subcat", "subcategory"),
    required = FALSE
  )

  data.frame(
    pathway = as.character(msig_df[[pathway_col]]),
    gene = toupper(trimws(as.character(msig_df[[gene_col]]))),
    subcollection = if (!is.na(subcollection_col)) as.character(msig_df[[subcollection_col]]) else NA_character_,
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(pathway), pathway != "", !is.na(gene), gene != "") %>%
    distinct()
}

pathway_list_from_df <- function(df) {
  split(df$gene, df$pathway) %>%
    lapply(function(x) sort(unique(toupper(trimws(x)))))
}

message("\nLoading MSigDB collections...")
msig_h <- standardise_msig(get_msig_collection("H"))
msig_c2 <- standardise_msig(get_msig_collection("C2"))
msig_c5 <- standardise_msig(get_msig_collection("C5"))

hallmark_pathways <- pathway_list_from_df(msig_h)

reactome_df <- msig_c2 %>%
  filter(
    grepl("^REACTOME_", toupper(pathway)) |
      grepl("REACTOME", toupper(coalesce(subcollection, "")))
  )
reactome_pathways <- pathway_list_from_df(reactome_df)

go_pathways <- pathway_list_from_df(msig_c5)

get_required_pathway <- function(pathways, requested, loose_pattern = NULL) {
  available <- names(pathways)
  exact <- available[toupper(available) == toupper(requested)]

  if (length(exact) == 0 && !is.null(loose_pattern)) {
    exact <- available[grepl(loose_pattern, toupper(available))]
  }

  if (length(exact) == 0) {
    stop("Required pathway was not found in MSigDB: ", requested)
  }

  genes <- sort(unique(toupper(pathways[[exact[1]]])))
  list(name = exact[1], genes = genes)
}

ifna_info <- get_required_pathway(
  hallmark_pathways,
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"
)

tnfa_info <- get_required_pathway(
  hallmark_pathways,
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

inflammatory_info <- get_required_pathway(
  hallmark_pathways,
  "HALLMARK_INFLAMMATORY_RESPONSE"
)

chemokine_info_raw <- get_required_pathway(
  go_pathways,
  "GOMF_CHEMOKINE_ACTIVITY",
  "CHEMOKINE_ACTIVITY"
)

canonical_chemokine_regex <- paste0(
  "^(CCL[0-9]+[A-Z0-9]*|CXCL[0-9]+[A-Z0-9]*|CX3CL[0-9]+[A-Z0-9]*|",
  "XCL[0-9]+[A-Z0-9]*|PF4|PF4V1|PPBP)$"
)
chemokine_genes <- chemokine_info_raw$genes[
  grepl(canonical_chemokine_regex, chemokine_info_raw$genes)
]
if (length(chemokine_genes) == 0) {
  stop("No canonical chemokine ligands remained after filtering GOMF_CHEMOKINE_ACTIVITY.")
}

antimicrobial_info <- get_required_pathway(
  reactome_pathways,
  "REACTOME_ANTIMICROBIAL_PEPTIDES",
  "ANTIMICROBIAL_PEPTIDES"
)

metaplasia_signature_file <- first_existing_file(metaplasia_signature_candidates)
if (is.na(metaplasia_signature_file)) {
  stop(
    "The final SSL metaplasia signature file was not found. Checked:\n",
    paste(metaplasia_signature_candidates, collapse = "\n")
  )
}

read_gene_list <- function(file) {
  if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    dat <- data.table::fread(file, data.table = FALSE, check.names = FALSE)
    gene_col <- detect_col(
      dat,
      c("gene_symbol", "Gene", "gene", "SYMBOL", "Symbol", "symbol"),
      required = FALSE
    )
    if (is.na(gene_col)) gene_col <- names(dat)[1]
    genes <- dat[[gene_col]]
  } else {
    lines <- readLines(file, warn = FALSE)
    if (length(lines) == 1 && grepl("[,;\\t]", lines[1])) {
      genes <- unlist(strsplit(lines, "[,;\\t]"))
    } else {
      genes <- lines
    }
  }

  genes <- unique(toupper(trimws(as.character(genes))))
  genes[genes != "" & !is.na(genes)]
}

metaplasia_genes <- read_gene_list(metaplasia_signature_file)
if (length(metaplasia_genes) < minimum_signature_genes) {
  stop("The metaplasia signature contains too few genes: ", metaplasia_signature_file)
}

if (length(metaplasia_genes) != 24) {
  warning(
    "The final SSL metaplasia file contains ", length(metaplasia_genes),
    " unique genes rather than 24. The script will use the file exactly as supplied."
  )
}

writeLines(
  c(
    paste0("Metaplasia signature file used: ", metaplasia_signature_file),
    paste0("Unique genes: ", length(metaplasia_genes))
  ),
  file.path(subdirs$qc, "metaplasia_signature_source.txt")
)

signature_list <- list(
  "Serrated Lesion Gastric Metaplasia Gene Signature" = metaplasia_genes,
  "Interferon alpha response" = ifna_info$genes,
  "TNFA/NFKB signalling" = tnfa_info$genes,
  "Inflammatory response" = inflammatory_info$genes
)

signature_file_stems <- c(
  "Serrated Lesion Gastric Metaplasia Gene Signature" = metaplasia_file_stem,
  "Interferon alpha response" = "Interferon_alpha_response",
  "TNFA/NFKB signalling" = "TNFA_NFKB_signalling",
  "Inflammatory response" = "Inflammatory_response"
)

# CPM-based singscores; gastric genes are not expression-filtered ----
score_signatures <- function(expr_obj, dataset_label) {
  score_input <- expr_obj$expression_for_ranking
  log_expr <- expr_obj$log_expr
  meta <- expr_obj$meta

  stopifnot(identical(colnames(score_input), meta$Sample_ID))

  message("Scoring signatures for ", dataset_label, " using singscore...")
  rank_data <- singscore::rankGenes(score_input, tiesMethod = "min")

  score_results <- list()
  coverage_results <- list()

  for (signature_name in names(signature_list)) {
    requested_genes <- unique(toupper(signature_list[[signature_name]]))
    present <- intersect(requested_genes, rownames(score_input))

    if (signature_name != "Serrated Lesion Gastric Metaplasia Gene Signature") {
      present <- present[rowMeans(log_expr[present, , drop = FALSE]) >= signature_mean_logcpm_cutoff]
    }

    coverage_results[[signature_name]] <- data.frame(
      Dataset = dataset_label,
      Signature = signature_name,
      Requested_genes = length(requested_genes),
      Present_before_filter = length(intersect(requested_genes, rownames(score_input))),
      Genes_used = length(present),
      Mean_logCPM_filter = ifelse(
        signature_name == "Serrated Lesion Gastric Metaplasia Gene Signature",
        "Not applied",
        paste0(">= ", signature_mean_logcpm_cutoff)
      ),
      stringsAsFactors = FALSE
    )

    if (length(present) < minimum_signature_genes) {
      warning(
        dataset_label, " / ", signature_name, " has only ", length(present),
        " usable genes; scores are returned as NA."
      )

      score_results[[signature_name]] <- meta %>%
        transmute(
          Sample_ID,
          Signature = signature_name,
          Score = NA_real_
        )
      next
    }

    score_obj <- singscore::simpleScore(
      rankData = rank_data,
      upSet = present,
      centerScore = TRUE,
      knownDirection = TRUE
    )

    if (length(score_obj$TotalScore) != ncol(score_input)) {
      stop(
        "singscore returned ", length(score_obj$TotalScore),
        " scores for ", ncol(score_input), " samples in ", dataset_label, "."
      )
    }

    score_df <- data.frame(
      Sample_ID = colnames(score_input),
      Signature = signature_name,
      Score = as.numeric(score_obj$TotalScore),
      stringsAsFactors = FALSE
    )

    score_results[[signature_name]] <- score_df
  }

  scores <- bind_rows(score_results) %>%
    left_join(meta, by = "Sample_ID")

  coverage <- bind_rows(coverage_results)

  data.table::fwrite(
    scores,
    file.path(subdirs$tables, paste0(safe_filename(dataset_label), "_signature_scores_long.csv"))
  )
  data.table::fwrite(
    coverage,
    file.path(subdirs$qc, paste0(safe_filename(dataset_label), "_signature_gene_coverage.csv"))
  )

  list(scores = scores, coverage = coverage)
}

infection_scores <- score_signatures(infection_obj, "Infection")
timecourse_scores <- score_signatures(timecourse_obj, "Timecourse")
braf_scores <- score_signatures(braf_obj, "BRAF_initial_passage_P1_or_P2")

# Intermediate score plots; final statistical exports are written below ----
run_two_group_test <- function(df, group_col, levels, pair_col = NULL, prefer_paired = FALSE) {
  dat <- df %>%
    filter(
      !is.na(.data[[group_col]]),
      .data[[group_col]] %in% levels,
      is.finite(Score)
    ) %>%
    mutate(.Group = factor(.data[[group_col]], levels = levels))

  paired_used <- FALSE
  n_pairs <- 0L
  p_value <- NA_real_
  line_df <- data.frame()

  if (prefer_paired && !is.null(pair_col) && pair_col %in% names(dat)) {
    complete_pairs <- dat %>%
      filter(!is.na(.data[[pair_col]]), .data[[pair_col]] != "") %>%
      group_by(.data[[pair_col]]) %>%
      filter(n_distinct(.Group) == 2, n() == 2) %>%
      ungroup()

    n_pairs <- n_distinct(complete_pairs[[pair_col]])

    if (n_pairs >= 3) {
      wide <- complete_pairs %>%
        dplyr::select(dplyr::all_of(pair_col), .Group, Score) %>%
        tidyr::pivot_wider(names_from = .Group, values_from = Score)

      p_value <- tryCatch(
        wilcox.test(wide[[levels[2]]], wide[[levels[1]]], paired = TRUE, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
      paired_used <- TRUE
      line_df <- complete_pairs
    }
  }

  if (!paired_used) {
    group_counts <- table(dat$.Group)
    if (all(group_counts >= 2)) {
      p_value <- tryCatch(
        wilcox.test(Score ~ .Group, data = dat, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
    }
  }

  list(
    data = dat,
    p_value = p_value,
    paired = paired_used,
    n_pairs = n_pairs,
    line_data = line_df
  )
}

make_two_group_score_plot <- function(
    df,
    signature_name,
    group_col,
    group_levels,
    colours,
    x_labels = group_levels,
    pair_col = NULL,
    prefer_paired = FALSE,
    y_label = "Singscore"
) {
  dat <- df %>% filter(Signature == signature_name)
  test <- run_two_group_test(
    dat,
    group_col = group_col,
    levels = group_levels,
    pair_col = pair_col,
    prefer_paired = prefer_paired
  )

  plot_df <- test$data
  if (nrow(plot_df) == 0) {
    return(make_empty_plot(signature_name, "No valid samples"))
  }

  plot_df$.Group <- factor(plot_df$.Group, levels = group_levels)
  ymax <- max(plot_df$Score, na.rm = TRUE)
  ymin <- min(plot_df$Score, na.rm = TRUE)
  span <- max(ymax - ymin, 0.05)
  label_y <- ymax + 0.15 * span

  p <- ggplot(plot_df, aes(x = .Group, y = Score, fill = .Group))

  if (test$paired && nrow(test$line_data) > 0) {
    line_data <- test$line_data %>%
      mutate(.Group = factor(.data[[group_col]], levels = group_levels))

    p <- p +
      geom_line(
        data = line_data,
        aes(x = .Group, y = Score, group = .data[[pair_col]]),
        inherit.aes = FALSE,
        colour = "grey72",
        linewidth = 0.35,
        alpha = 0.8
      )
  }

  p +
    geom_boxplot(
      width = 0.55,
      outlier.shape = NA,
      linewidth = 0.35,
      alpha = 0.45,
      colour = "black"
    ) +
    geom_point(
      shape = 21,
      size = 2.0,
      stroke = 0.30,
      colour = "black",
      position = position_jitter(width = 0.07, height = 0),
      alpha = 0.95
    ) +
    annotate(
      "segment",
      x = 1, xend = 2,
      y = label_y, yend = label_y,
      linewidth = 0.3,
      colour = "black"
    ) +
    annotate(
      "text",
      x = 1.5,
      y = label_y + 0.04 * span,
      label = p_to_plot_label(test$p_value),
      family = font_family,
      size = 6 / 2.845276
    ) +
    scale_fill_manual(values = colours, breaks = group_levels, drop = FALSE) +
    scale_x_discrete(labels = setNames(x_labels, group_levels), drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    labs(
      title = format_signature_plot_title(signature_name),
      x = NULL,
      y = y_label,
      fill = NULL
    ) +
    nature_theme() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 6)
    )
}

make_timecourse_plot <- function(df, signature_name, y_label = "Singscore") {
  dat <- df %>%
    filter(
      Signature == signature_name,
      Type %in% c("SSL_HP", "H_ADJ_SSL_HP"),
      is.finite(Passage),
      is.finite(Score)
    ) %>%
    arrange(Type, Line_ID, Passage)

  if (nrow(dat) == 0) {
    return(make_empty_plot(signature_name, "No valid time-course samples"))
  }

  mean_df <- dat %>%
    group_by(Type, Passage) %>%
    summarise(
      Mean = mean(Score),
      SEM = ifelse(n() > 1, sd(Score) / sqrt(n()), NA_real_),
      .groups = "drop"
    )

  ggplot(dat, aes(x = Passage, y = Score)) +
    geom_line(
      aes(group = Trajectory_ID, colour = Type),
      linewidth = 0.35,
      alpha = 0.28
    ) +
    geom_point(
      aes(fill = Type),
      shape = 21,
      size = 1.8,
      stroke = 0.28,
      colour = "black",
      alpha = 0.85
    ) +
    geom_line(
      data = mean_df,
      aes(x = Passage, y = Mean, colour = Type, group = Type),
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    geom_errorbar(
      data = mean_df,
      aes(x = Passage, ymin = Mean - SEM, ymax = Mean + SEM, colour = Type),
      width = 0.10,
      linewidth = 0.35,
      inherit.aes = FALSE,
      na.rm = TRUE
    ) +
    geom_point(
      data = mean_df,
      aes(x = Passage, y = Mean, fill = Type),
      shape = 21,
      size = 2.3,
      stroke = 0.35,
      colour = "black",
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = paper_cols[c("SSL_HP", "H_ADJ_SSL_HP")], drop = FALSE) +
    scale_colour_manual(values = paper_cols[c("SSL_HP", "H_ADJ_SSL_HP")], drop = FALSE) +
    scale_x_continuous(breaks = sort(unique(dat$Passage))) +
    labs(
      title = format_signature_plot_title(signature_name),
      x = "Passage",
      y = y_label,
      fill = "Sample type",
      colour = "Sample type"
    ) +
    nature_theme() +
    theme(legend.position = "right") +
    guides(colour = "none")
}

infection_score_plots <- list()
infection_score_stats <- list()

for (sig in names(signature_list)) {
  p <- make_two_group_score_plot(
    infection_scores$scores,
    signature_name = sig,
    group_col = "Infection_Status",
    group_levels = c("Non-Infected", "Infected"),
    colours = infection_cols,
    x_labels = c("Non-Infected", "Infected"),
    pair_col = "Pair_ID",
    prefer_paired = TRUE
  )

  infection_score_plots[[sig]] <- p
  save_plot_set(
    p,
    subdirs$scores,
    paste0("Infection_", signature_file_stems[[sig]], "_singscore"),
    62,
    70
  )

  test <- run_two_group_test(
    infection_scores$scores %>% filter(Signature == sig),
    "Infection_Status",
    c("Non-Infected", "Infected"),
    "Pair_ID",
    TRUE
  )

  infection_score_stats[[sig]] <- data.frame(
    Dataset = "Infection",
    Signature = sig,
    Test = ifelse(test$paired, "Paired Wilcoxon signed-rank", "Wilcoxon rank-sum"),
    Complete_pairs = test$n_pairs,
    P_value = test$p_value,
    stringsAsFactors = FALSE
  )
}

data.table::fwrite(
  bind_rows(infection_score_stats),
  file.path(subdirs$tables, "Infection_signature_score_statistics.csv")
)

infection_score_combined <- wrap_plots(infection_score_plots, ncol = 4) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(family = font_family, face = "bold", size = panel_tag_pt))
  )
save_plot_set(infection_score_combined, subdirs$combined, "Infection_signature_scores_combined", 183, 72)

braf_score_input <- braf_scores$scores %>%
  filter(
    Type == "SSL_HP",
    BRAF_Status %in% c("WT", "V600E")
  )

braf_score_plots <- list()
braf_score_stats <- list()

for (sig in names(signature_list)) {
  p <- make_two_group_score_plot(
    braf_score_input,
    signature_name = sig,
    group_col = "BRAF_Status",
    group_levels = c("WT", "V600E"),
    colours = braf_cols,
    x_labels = c("WT", "V600E"),
    pair_col = NULL,
    prefer_paired = FALSE
  )

  braf_score_plots[[sig]] <- p
  save_plot_set(
    p,
    subdirs$scores,
    paste0("BRAF_SSL_HP_only_", signature_file_stems[[sig]], "_singscore"),
    62,
    70
  )

  test <- run_two_group_test(
    braf_score_input %>% filter(Signature == sig),
    "BRAF_Status",
    c("WT", "V600E"),
    NULL,
    FALSE
  )

  braf_score_stats[[sig]] <- data.frame(
    Dataset = "BRAF_SSL_HP_only",
    Signature = sig,
    Test = "Wilcoxon rank-sum",
    Complete_pairs = NA_integer_,
    P_value = test$p_value,
    stringsAsFactors = FALSE
  )
}

data.table::fwrite(
  bind_rows(braf_score_stats),
  file.path(subdirs$tables, "BRAF_SSL_HP_only_signature_score_statistics.csv")
)

braf_score_combined <- wrap_plots(braf_score_plots, ncol = 4) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(family = font_family, face = "bold", size = panel_tag_pt))
  )
save_plot_set(braf_score_combined, subdirs$combined, "BRAF_SSL_HP_only_signature_scores_combined", 183, 72)

timecourse_score_plots <- list()
timecourse_spearman_stats <- list()

for (sig in names(signature_list)) {
  p <- make_timecourse_plot(timecourse_scores$scores, sig)
  timecourse_score_plots[[sig]] <- p

  save_plot_set(
    p,
    subdirs$scores,
    paste0("Timecourse_", signature_file_stems[[sig]], "_singscore"),
    90,
    70
  )

  stat_one <- timecourse_scores$scores %>%
    filter(
      Signature == sig,
      Type %in% c("SSL_HP", "H_ADJ_SSL_HP"),
      is.finite(Passage),
      is.finite(Score)
    ) %>%
    group_by(Type) %>%
    group_modify(function(.x, .y) {
      tst <- if (nrow(.x) >= 3 && n_distinct(.x$Passage) >= 2) {
        suppressWarnings(cor.test(.x$Passage, .x$Score, method = "spearman", exact = FALSE))
      } else {
        NULL
      }

      data.frame(
        N = nrow(.x),
        Spearman_rho = if (is.null(tst)) NA_real_ else unname(tst$estimate),
        P_value = if (is.null(tst)) NA_real_ else tst$p.value
      )
    }) %>%
    ungroup() %>%
    mutate(Signature = sig, .before = 1)

  timecourse_spearman_stats[[sig]] <- stat_one
}

data.table::fwrite(
  bind_rows(timecourse_spearman_stats),
  file.path(subdirs$tables, "Timecourse_signature_score_Spearman_statistics.csv")
)

timecourse_score_combined <- wrap_plots(timecourse_score_plots, ncol = 2) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(family = font_family, face = "bold", size = panel_tag_pt))
  )
save_plot_set(timecourse_score_combined, subdirs$combined, "Timecourse_signature_scores_combined", 183, 130)

extract_gene_expression <- function(expr_obj, gene) {
  gene <- toupper(gene)
  if (!gene %in% rownames(expr_obj$log_expr)) {
    warning(gene, " is absent from ", expr_obj$dataset, " expression matrix.")
    return(NULL)
  }

  data.frame(
    Sample_ID = colnames(expr_obj$log_expr),
    Signature = gene,
    Score = as.numeric(expr_obj$log_expr[gene, ]),
    stringsAsFactors = FALSE
  ) %>%
    left_join(expr_obj$meta, by = "Sample_ID")
}

muc5ac_infection <- extract_gene_expression(infection_obj, "MUC5AC")
muc5ac_timecourse <- extract_gene_expression(timecourse_obj, "MUC5AC")
muc5ac_braf <- extract_gene_expression(braf_obj, "MUC5AC")

if (!is.null(muc5ac_infection)) {
  p_muc_inf <- make_two_group_score_plot(
    muc5ac_infection,
    "MUC5AC",
    "Infection_Status",
    c("Non-Infected", "Infected"),
    infection_cols,
    c("Non-Infected", "Infected"),
    "Pair_ID",
    TRUE,
    y_label = "MUC5AC log2 expression"
  )
  save_plot_set(p_muc_inf, subdirs$scores, "Infection_MUC5AC_expression", 62, 70)
}

if (!is.null(muc5ac_braf)) {
  p_muc_braf <- make_two_group_score_plot(
    muc5ac_braf %>% filter(Type == "SSL_HP"),
    "MUC5AC",
    "BRAF_Status",
    c("WT", "V600E"),
    braf_cols,
    c("WT", "V600E"),
    NULL,
    FALSE,
    y_label = "MUC5AC log2 expression"
  )
  save_plot_set(p_muc_braf, subdirs$scores, "BRAF_SSL_HP_only_MUC5AC_expression", 62, 70)
}

if (!is.null(muc5ac_timecourse)) {
  p_muc_time <- make_timecourse_plot(
    muc5ac_timecourse,
    "MUC5AC",
    y_label = "MUC5AC log2 expression"
  )
  save_plot_set(p_muc_time, subdirs$scores, "Timecourse_MUC5AC_expression", 90, 70)
}

# Differential-expression tables and ranked GSEA ----
standardise_de_table <- function(file, comparison_label) {

  de <- data.table::fread(
    file,
    data.table = FALSE,
    check.names = FALSE,
    fill = TRUE,
    na.strings = c("", "NA", "N/A", "na")
  )

  names(de) <- trimws(names(de))

  nonempty_columns <- vapply(
    de,
    function(x) any(!is.na(x) & trimws(as.character(x)) != ""),
    logical(1)
  )
  de <- de[, nonempty_columns, drop = FALSE]

  id_col <- detect_col(
    de,
    c("genes", "gene_id", "ENSEMBL", "Ensembl_ID", "Geneid", "ID"),
    required = TRUE,
    label = paste0(comparison_label, " Ensembl gene-ID column")
  )
  symbol_col <- detect_col(
    de,
    c("SYMBOL", "Symbol", "symbol", "gene_symbol", "GeneSymbol", "HGNC.symbol"),
    required = TRUE,
    label = paste0(comparison_label, " gene-symbol column")
  )
  logfc_col <- detect_col(
    de,
    c("logFC", "LogFC", "logfc", "log2FC", "log2FoldChange", "LFC"),
    required = TRUE,
    label = paste0(comparison_label, " logFC column")
  )
  logcpm_col <- detect_col(
    de,
    c("logCPM", "LogCPM", "logcpm", "log2CPM", "AveExpr", "average_expression"),
    required = TRUE,
    label = paste0(comparison_label, " logCPM column")
  )
  lr_col <- detect_col(
    de,
    c("LR", "LRT", "likelihood_ratio", "LikelihoodRatio"),
    required = FALSE,
    label = paste0(comparison_label, " LR column")
  )
  t_col <- detect_col(
    de,
    c("t", "T", "stat", "STAT", "Statistic", "moderated_t"),
    required = FALSE,
    label = paste0(comparison_label, " signed test-statistic column")
  )
  p_col <- detect_col(
    de,
    c("PValue", "P.Value", "Pvalue", "pvalue", "p.value", "p", "P"),
    required = TRUE,
    label = paste0(comparison_label, " P-value column")
  )
  fdr_col <- detect_col(
    de,
    c("FDR", "fdr", "padj", "adj.P.Val", "adj_p", "qvalue", "q.value"),
    required = FALSE,
    label = paste0(comparison_label, " FDR column")
  )

  genename_col <- detect_col(
    de,
    c("GENENAME", "GeneName", "gene_name", "GENE_NAME"),
    required = FALSE
  )
  entrez_col <- detect_col(
    de,
    c("ENTREZID", "ENTREZ_ID", "EntrezID", "entrez_id", "ENTREZGENE"),
    required = FALSE
  )
  genetype_col <- detect_col(
    de,
    c("GENETYPE", "GENE_TYPE", "GeneType", "gene_type", "GENEBIOTYPE", "BIOTYPE"),
    required = FALSE
  )
  map_col <- detect_col(
    de,
    c("MAP", "MAPLOCATION", "MAP_LOCATION", "Map"),
    required = FALSE
  )
  path_col <- detect_col(
    de,
    c("PATH", "Path", "path", "PATHWAY", "PATH_ID"),
    required = FALSE
  )

  column_mapping <- data.frame(
    Standard_role = c(
      "Ensembl_ID", "Gene_symbol", "logFC", "logCPM", "LR", "signed_statistic",
      "PValue", "FDR", "GENENAME", "ENTREZID", "GENETYPE", "MAP", "PATH"
    ),
    Input_column = c(
      id_col, symbol_col, logfc_col, logcpm_col, lr_col, t_col,
      p_col, fdr_col, genename_col, entrez_col, genetype_col, map_col, path_col
    ),
    Used_for = c(
      "Fallback symbol mapping and audit",
      "Primary GSEA/scoring/heatmap identifier",
      "Direction and effect size",
      "Expression filter",
      "Primary edgeR ranking strength",
      "Fallback signed ranking for limma-style input",
      "Nominal significance and ranking fallback",
      "Adjusted significance",
      rep("Annotation only", 5)
    ),
    stringsAsFactors = FALSE
  )

  data.table::fwrite(
    column_mapping,
    file.path(
      subdirs$qc,
      paste0(safe_filename(comparison_label), "_DE_column_mapping.csv")
    )
  )

  convert_numeric_de <- function(x, input_column, role) {
    original <- clean_character(x)
    converted <- suppressWarnings(as.numeric(original))
    invalid <- !is.na(original) & is.na(converted)

    if (any(invalid)) {
      examples <- unique(original[invalid])
      examples <- utils::head(examples, 8)
      stop(
        "Non-numeric values were found in the ", role,
        " field ('", input_column, "') of ", basename(file), ": ",
        paste(examples, collapse = "; ")
      )
    }

    converted
  }

  gene_ids_original <- clean_character(de[[id_col]])
  supplied_symbols <- clean_character(de[[symbol_col]])
  supplied_symbols[
    toupper(supplied_symbols) %in% c("NA", "N/A", ".", "NULL")
  ] <- NA_character_

  mapped_symbols <- map_gene_symbols(
    gene_ids = gene_ids_original,
    supplied_symbols = supplied_symbols
  )

  gene_symbol_source <- ifelse(
    !is.na(supplied_symbols),
    "SYMBOL column",
    ifelse(!is.na(mapped_symbols), "Mapped from Ensembl genes column", NA_character_)
  )

  out <- data.frame(
    Ensembl_ID = toupper(stringr::str_remove(gene_ids_original, "\\.[0-9]+$")),
    Gene = mapped_symbols,
    logFC = convert_numeric_de(de[[logfc_col]], logfc_col, "logFC"),
    logCPM = convert_numeric_de(de[[logcpm_col]], logcpm_col, "logCPM"),
    LR = if (!is.na(lr_col)) {
      convert_numeric_de(de[[lr_col]], lr_col, "LR")
    } else {
      NA_real_
    },
    t = if (!is.na(t_col)) {
      convert_numeric_de(de[[t_col]], t_col, "signed test statistic")
    } else {
      NA_real_
    },
    PValue = convert_numeric_de(de[[p_col]], p_col, "PValue"),
    FDR = if (!is.na(fdr_col)) {
      convert_numeric_de(de[[fdr_col]], fdr_col, "FDR")
    } else {
      NA_real_
    },
    GENENAME = if (!is.na(genename_col)) clean_character(de[[genename_col]]) else NA_character_,
    ENTREZID = if (!is.na(entrez_col)) clean_character(de[[entrez_col]]) else NA_character_,
    GENETYPE = if (!is.na(genetype_col)) clean_character(de[[genetype_col]]) else NA_character_,
    MAP = if (!is.na(map_col)) clean_character(de[[map_col]]) else NA_character_,
    PATH = if (!is.na(path_col)) clean_character(de[[path_col]]) else NA_character_,
    Gene_symbol_source = gene_symbol_source,
    Comparison = comparison_label,
    Source_file = basename(file),
    stringsAsFactors = FALSE
  )

  de_qc <- data.frame(
    Comparison = comparison_label,
    Source_file = basename(file),
    Input_rows = nrow(out),
    Rows_with_nonempty_Ensembl_ID = sum(!is.na(out$Ensembl_ID)),
    Rows_with_supplied_SYMBOL = sum(!is.na(supplied_symbols)),
    Rows_with_final_gene_symbol = sum(!is.na(out$Gene)),
    Rows_mapped_from_Ensembl = sum(out$Gene_symbol_source == "Mapped from Ensembl genes column", na.rm = TRUE),
    Rows_missing_final_gene_symbol = sum(is.na(out$Gene)),
    Finite_logFC = sum(is.finite(out$logFC)),
    Finite_logCPM = sum(is.finite(out$logCPM)),
    Finite_LR = sum(is.finite(out$LR)),
    Finite_PValue = sum(is.finite(out$PValue)),
    Finite_FDR = sum(is.finite(out$FDR)),
    stringsAsFactors = FALSE
  )

  data.table::fwrite(
    de_qc,
    file.path(
      subdirs$qc,
      paste0(safe_filename(comparison_label), "_DE_structure_QC.csv")
    )
  )

  materially_negative_lr <- is.finite(out$LR) & out$LR < -1e-8
  if (any(materially_negative_lr)) {
    negative_rows <- out[materially_negative_lr, , drop = FALSE]
    data.table::fwrite(
      negative_rows,
      file.path(
        subdirs$qc,
        paste0(safe_filename(comparison_label), "_DE_materially_negative_LR_rows.csv")
      )
    )

    stop(
      comparison_label, " contains ", sum(materially_negative_lr),
      " materially negative LR value(s). Minimum LR = ",
      min(out$LR[materially_negative_lr]),
      ". The affected rows were written to the QC folder. Do not use abs(LR)."
    )
  }

  out$LR[is.finite(out$LR) & out$LR < 0] <- 0

  out <- out %>%
    dplyr::filter(
      !is.na(Gene), Gene != "",
      is.finite(logFC),
      is.finite(logCPM),
      is.finite(PValue)
    )

  if (nrow(out) == 0) {
    stop("No valid DE rows remained after parsing ", basename(file), ".")
  }

  if (all(is.na(out$FDR))) {
    out$FDR <- p.adjust(out$PValue, method = "BH")
  }

  positive_p <- out$PValue[is.finite(out$PValue) & out$PValue > 0]
  minimum_positive_p <- if (length(positive_p) > 0) min(positive_p) else 1e-300
  p_safe <- ifelse(
    !is.finite(out$PValue) | out$PValue <= 0,
    minimum_positive_p * 0.1,
    out$PValue
  )

  out$.collapse_rank <- dplyr::case_when(
    is.finite(out$LR) ~ sign(out$logFC) * sqrt(out$LR),
    is.finite(out$t) ~ out$t,
    TRUE ~ sign(out$logFC) * -log10(p_safe)
  )

  duplicate_symbol_rows <- out %>%
    dplyr::group_by(Gene) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup()

  if (nrow(duplicate_symbol_rows) > 0) {
    data.table::fwrite(
      duplicate_symbol_rows,
      file.path(
        subdirs$qc,
        paste0(safe_filename(comparison_label), "_DE_duplicate_symbol_rows_before_collapse.csv")
      )
    )
  }

  out <- out %>%
    dplyr::group_by(Gene) %>%
    dplyr::arrange(
      dplyr::desc(abs(.collapse_rank)),
      PValue,
      dplyr::desc(abs(logFC)),
      .by_group = TRUE
    ) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.collapse_rank)

  out
}

infection_de <- standardise_de_table(
  files$infection_de,
  "Infected vs Non-Infected"
)

data.table::fwrite(
  infection_de,
  file.path(subdirs$tables, "Infection_DE_standardised.csv")
)

braf_de <- standardise_de_table(
  files$braf_de,
  "BRAF V600E vs WT in SSL_HP"
)

braf_de_source <- paste0(
  "Supplied edgeR DE table: ",
  files$braf_de,
  ". Contrast interpreted as V600E - WT in SSL_HP."
)

data.table::fwrite(
  braf_de,
  file.path(subdirs$tables, "BRAF_SSL_HP_only_DE_standardised.csv")
)

writeLines(
  braf_de_source,
  file.path(subdirs$qc, "BRAF_DE_source.txt")
)

make_ranked_vector <- function(
    de,
    comparison_label,
    apply_logcpm_filter = TRUE
) {
  dat <- de

  if (apply_logcpm_filter && "logCPM" %in% names(dat) && any(is.finite(dat$logCPM))) {
    dat <- dat %>%
      dplyr::filter(is.finite(logCPM), logCPM >= gsea_logcpm_cutoff)
  }

  rank_method <- NA_character_

  if (sum(is.finite(dat$LR)) >= 100) {
    materially_negative <- dat$LR < -1e-8 & is.finite(dat$LR)
    if (any(materially_negative)) {
      stop(
        comparison_label, " contains ", sum(materially_negative),
        " materially negative LR values. Minimum LR = ",
        min(dat$LR[materially_negative]),
        ". LR must not be made positive with abs(). Inspect the DE table."
      )
    }

    dat$LR[is.finite(dat$LR) & dat$LR < 0] <- 0
    dat$rank_metric <- sign(dat$logFC) * sqrt(dat$LR)
    rank_method <- "sign(logFC) * sqrt(LR)"
  } else if (sum(is.finite(dat$t)) >= 100) {
    dat$rank_metric <- dat$t
    rank_method <- "signed limma t statistic"
  } else {
    positive_p <- dat$PValue[is.finite(dat$PValue) & dat$PValue > 0]
    min_positive_p <- if (length(positive_p) > 0) min(positive_p) else 1e-300
    safe_p <- ifelse(
      dat$PValue <= 0 | !is.finite(dat$PValue),
      min_positive_p * 0.1,
      dat$PValue
    )
    dat$rank_metric <- sign(dat$logFC) * -log10(safe_p)
    rank_method <- "sign(logFC) * -log10(PValue)"
  }

  dat <- dat %>%
    dplyr::filter(is.finite(rank_metric), !is.na(Gene), Gene != "") %>%
    dplyr::group_by(Gene) %>%
    dplyr::arrange(dplyr::desc(abs(rank_metric)), PValue, .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(rank_metric), Gene)

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
    method = rank_method,
    filter_mode = ifelse(apply_logcpm_filter, "LOGCPM_GE_2", "ALL_VALID_GENES")
  )
}

run_fgsea_one <- function(ranks, pathways, comparison_label, family_label, min_size = fgsea_min_size) {
  pathways <- lapply(pathways, function(x) unique(toupper(x)))
  overlap <- vapply(pathways, function(x) sum(x %in% names(ranks)), integer(1))
  pathways <- pathways[overlap >= min_size & overlap <= fgsea_max_size]

  if (length(pathways) == 0) {
    warning("No pathways passed overlap filters for ", comparison_label, " / ", family_label)
    return(data.frame())
  }

  set.seed(1)
  fg <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = min_size,
    maxSize = fgsea_max_size,
    eps = fgsea_eps,
    nproc = fgsea_nproc,
    nPermSimple = fgsea_nperm_simple,
    scoreType = "std"
  )

  fg <- as.data.frame(fg)
  if (!"padj" %in% names(fg)) fg$padj <- p.adjust(fg$pval, method = "BH")
  if (!"log2err" %in% names(fg)) fg$log2err <- NA_real_
  if (!"leadingEdge" %in% names(fg)) fg$leadingEdge <- vector("list", nrow(fg))

  fg %>%
    mutate(
      Comparison = comparison_label,
      Gene_set_family = family_label,
      leadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
    ) %>%
    arrange(padj, desc(abs(NES)))
}

infection_rank <- make_ranked_vector(
  infection_de,
  "Infected vs Non-Infected",
  apply_logcpm_filter = TRUE
)
braf_rank <- make_ranked_vector(
  braf_de,
  "BRAF V600E vs WT in SSL_HP",
  apply_logcpm_filter = TRUE
)

infection_metaplasia_rank <- make_ranked_vector(
  infection_de,
  "Infected vs Non-Infected",
  apply_logcpm_filter = FALSE
)
braf_metaplasia_rank <- make_ranked_vector(
  braf_de,
  "BRAF V600E vs WT in SSL_HP",
  apply_logcpm_filter = FALSE
)

data.table::fwrite(
  infection_rank$table,
  file.path(subdirs$tables, "Infection_LOGCPM_GE_2_ranked_gene_table.csv")
)
data.table::fwrite(
  braf_rank$table,
  file.path(subdirs$tables, "BRAF_SSL_HP_only_LOGCPM_GE_2_ranked_gene_table.csv")
)
data.table::fwrite(
  infection_metaplasia_rank$table,
  file.path(subdirs$tables, "Infection_ALL_VALID_GENES_ranked_gene_table_for_metaplasia.csv")
)
data.table::fwrite(
  braf_metaplasia_rank$table,
  file.path(subdirs$tables, "BRAF_SSL_HP_only_ALL_VALID_GENES_ranked_gene_table_for_metaplasia.csv")
)

infection_hallmark_fgsea <- run_fgsea_one(
  infection_rank$ranks,
  hallmark_pathways,
  "Infected vs Non-Infected",
  "HALLMARK"
)

braf_hallmark_fgsea <- run_fgsea_one(
  braf_rank$ranks,
  hallmark_pathways,
  "BRAF V600E vs WT in SSL_HP",
  "HALLMARK"
)

custom_metaplasia_pathway <- stats::setNames(list(metaplasia_genes), metaplasia_pathway_id)

infection_metaplasia_fgsea <- run_fgsea_one(
  infection_metaplasia_rank$ranks,
  custom_metaplasia_pathway,
  "Infected vs Non-Infected",
  "CUSTOM_METAPLASIA_NO_LOGCPM_FILTER",
  min_size = 3L
) %>%
  dplyr::mutate(
    Filter_mode = "ALL_VALID_GENES",
    N_ranked_genes = length(infection_metaplasia_rank$ranks),
    P_value_interpretation = paste0(
      "Nominal custom-set p value. One custom pathway was tested; ",
      "padj equals pval and is not a Hallmark-family FDR."
    )
  )

braf_metaplasia_fgsea <- run_fgsea_one(
  braf_metaplasia_rank$ranks,
  custom_metaplasia_pathway,
  "BRAF V600E vs WT in SSL_HP",
  "CUSTOM_METAPLASIA_NO_LOGCPM_FILTER",
  min_size = 3L
) %>%
  dplyr::mutate(
    Filter_mode = "ALL_VALID_GENES",
    N_ranked_genes = length(braf_metaplasia_rank$ranks),
    P_value_interpretation = paste0(
      "Nominal custom-set p value. One custom pathway was tested; ",
      "padj equals pval and is not a Hallmark-family FDR."
    )
  )

combined_fgsea <- dplyr::bind_rows(
  infection_hallmark_fgsea %>%
    dplyr::mutate(Filter_mode = "LOGCPM_GE_2"),
  braf_hallmark_fgsea %>%
    dplyr::mutate(Filter_mode = "LOGCPM_GE_2"),
  infection_metaplasia_fgsea,
  braf_metaplasia_fgsea
)

data.table::fwrite(
  combined_fgsea,
  file.path(subdirs$tables, "Organoid_combined_FGSEA_results.csv")
)

data.table::fwrite(
  dplyr::bind_rows(infection_hallmark_fgsea, braf_hallmark_fgsea),
  file.path(subdirs$tables, "Organoid_Hallmark_FGSEA_results.csv")
)

data.table::fwrite(
  dplyr::bind_rows(infection_metaplasia_fgsea, braf_metaplasia_fgsea),
  file.path(subdirs$tables, "Organoid_custom_metaplasia_FGSEA_ALL_VALID_GENES.csv")
)

writeLines(
  c(
    paste0("Infection Hallmark ranking: ", infection_rank$method),
    paste0("BRAF Hallmark ranking: ", braf_rank$method),
    paste0("Hallmark GSEA expression filter: logCPM >= ", gsea_logcpm_cutoff),
    paste0(metaplasia_display_name, " GSEA expression filter: none; ALL_VALID_GENES ranking is used."),
    paste0("BRAF DE source: ", braf_de_source),
    "Positive infection NES/logFC means higher in infected organoids.",
    "Positive BRAF NES/logFC means higher in V600E SSL/HP organoids than WT SSL/HP organoids.",
    "H_ADJ samples were excluded from BRAF DE, GSEA and individual-score comparisons."
  ),
  file.path(subdirs$qc, "GSEA_orientation_and_ranking_notes.txt")
)

# Pathway and expression figures ----
make_hallmark_barplot <- function(fgsea_df, comparison_label) {
  dat <- fgsea_df %>%
    filter(
      Comparison == comparison_label,
      Gene_set_family == "HALLMARK",
      is.finite(NES),
      is.finite(padj)
    ) %>%
    arrange(padj, desc(abs(NES))) %>%
    slice_head(n = hallmark_bar_top_n) %>%
    mutate(
      Pathway_label = clean_pathway_label(pathway),
      Pathway_label = factor(Pathway_label, levels = rev(Pathway_label)),
      Significance = case_when(
        padj < 0.001 ~ "***",
        padj < 0.01 ~ "**",
        padj < 0.05 ~ "*",
        TRUE ~ ""
      ),
      Text_hjust = ifelse(NES >= 0, -0.15, 1.15)
    )

  if (nrow(dat) == 0) {
    return(make_empty_plot(comparison_label, "No Hallmark results"))
  }

  ggplot(dat, aes(x = Pathway_label, y = NES, fill = NES)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
    geom_col(colour = "black", linewidth = 0.25, width = 0.72) +
    geom_text(
      aes(label = Significance, hjust = Text_hjust),
      family = font_family,
      size = 5.8 / 2.845276,
      show.legend = FALSE
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "NES"
    ) +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.18))) +
    labs(
      title = comparison_label,
      x = NULL,
      y = "Normalised enrichment score"
    ) +
    nature_theme() +
    theme(
      legend.position = "right",
      axis.text.y = element_text(size = 5.8),
      plot.margin = margin(4, 8, 4, 4, "mm")
    )
}

p_bar_infection <- make_hallmark_barplot(combined_fgsea, "Infected vs Non-Infected")
p_bar_braf <- make_hallmark_barplot(combined_fgsea, "BRAF V600E vs WT in SSL_HP")

save_plot_set(p_bar_infection, subdirs$gsea, "Infection_Hallmark_GSEA_top_barplot", 100, 95)
save_plot_set(p_bar_braf, subdirs$gsea, "BRAF_SSL_HP_only_Hallmark_GSEA_top_barplot", 100, 95)

barplots_combined <- p_bar_infection | p_bar_braf
barplots_combined <- barplots_combined +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(family = font_family, face = "bold", size = panel_tag_pt))
  )
save_plot_set(barplots_combined, subdirs$combined, "Hallmark_GSEA_barplots_combined", 183, 100)

make_hallmark_heatmap <- function(fgsea_df, top_only = FALSE) {
  dat <- fgsea_df %>%
    filter(Gene_set_family == "HALLMARK", is.finite(NES))

  if (top_only) {
    selected <- dat %>%
      group_by(Comparison) %>%
      arrange(padj, desc(abs(NES)), .by_group = TRUE) %>%
      slice_head(n = hallmark_heatmap_top_n) %>%
      ungroup() %>%
      distinct(pathway) %>%
      pull(pathway)
    dat <- dat %>% filter(pathway %in% selected)
  }

  dat <- dat %>%
    mutate(
      Comparison = factor(
        Comparison,
        levels = c("Infected vs Non-Infected", "BRAF V600E vs WT in SSL_HP")
      ),
      Pathway_label = clean_pathway_label(pathway),
      Significance = case_when(
        padj < 0.001 ~ "***",
        padj < 0.01 ~ "**",
        padj < 0.05 ~ "*",
        pval < 0.05 ~ ".",
        TRUE ~ ""
      )
    )

  pathway_order <- dat %>%
    group_by(Pathway_label) %>%
    summarise(Max_abs_NES = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
    arrange(Max_abs_NES) %>%
    pull(Pathway_label)

  dat$Pathway_label <- factor(dat$Pathway_label, levels = pathway_order)

  ggplot(dat, aes(x = Comparison, y = Pathway_label, fill = NES)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label = Significance), family = font_family, size = 5.5 / 2.845276) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "NES"
    ) +
    labs(x = NULL, y = NULL) +
    nature_theme() +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 6),
      axis.text.y = element_text(size = ifelse(top_only, 6, 5)),
      legend.position = "right"
    )
}

p_hallmark_heatmap_top <- make_hallmark_heatmap(combined_fgsea, top_only = TRUE)
p_hallmark_heatmap_all <- make_hallmark_heatmap(combined_fgsea, top_only = FALSE)

save_plot_set(p_hallmark_heatmap_top, subdirs$heatmaps, "Hallmark_NES_heatmap_top_pathways", 120, 110)
save_plot_set(p_hallmark_heatmap_all, subdirs$heatmaps, "Hallmark_NES_heatmap_all_50_pathways", 130, 190)

make_single_hallmark_heatmap <- function(
    fgsea_df,
    comparison_label,
    title_text,
    top_n = hallmark_heatmap_top_n
) {
  dat <- fgsea_df %>%
    dplyr::filter(
      Comparison == comparison_label,
      Gene_set_family == "HALLMARK",
      is.finite(NES)
    ) %>%
    dplyr::arrange(padj, dplyr::desc(abs(NES))) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      Pathway_label = clean_pathway_label(pathway),
      Pathway_label = factor(Pathway_label, levels = rev(Pathway_label)),
      Comparison_display = factor(title_text, levels = title_text),
      Significance = dplyr::case_when(
        padj < 0.001 ~ "***",
        padj < 0.01 ~ "**",
        padj < 0.05 ~ "*",
        pval < 0.05 ~ ".",
        TRUE ~ ""
      )
    )

  if (nrow(dat) == 0) {
    return(make_empty_plot(title_text, "No Hallmark GSEA results"))
  }

  ggplot2::ggplot(
    dat,
    ggplot2::aes(x = Comparison_display, y = Pathway_label, fill = NES)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::geom_text(
      ggplot2::aes(label = Significance),
      family = font_family,
      size = 5.5 / 2.845276
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "NES"
    ) +
    ggplot2::labs(title = format_signature_plot_title(title_text), x = NULL, y = NULL) +
    nature_theme() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = font_family,
        face = "bold",
        hjust = 0,
        size = plot_title_pt
      ),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 5.7),
      axis.text.y = ggplot2::element_text(size = 5.5),
      legend.position = "bottom"
    )
}

p_hallmark_heatmap_infection_single <- make_single_hallmark_heatmap(
  combined_fgsea,
  "Infected vs Non-Infected",
  "Hallmark GSEA: Infected vs Non-Infected"
)

p_hallmark_heatmap_braf_single <- make_single_hallmark_heatmap(
  combined_fgsea,
  "BRAF V600E vs WT in SSL_HP",
  "Hallmark GSEA: BRAF V600E vs WT in SSL/HP"
)

save_plot_set(
  p_hallmark_heatmap_infection_single,
  subdirs$heatmaps,
  "Infection_Hallmark_NES_heatmap_top_pathways",
  95,
  120
)

save_plot_set(
  p_hallmark_heatmap_braf_single,
  subdirs$heatmaps,
  "BRAF_SSL_HP_only_Hallmark_NES_heatmap_top_pathways",
  95,
  120
)

calc_enrichment_curve <- function(ranks, pathway_genes) {
  ranks <- sort(ranks, decreasing = TRUE)
  pathway_genes <- intersect(unique(toupper(pathway_genes)), names(ranks))

  N <- length(ranks)
  Nh <- length(pathway_genes)
  if (Nh < 1 || N <= Nh) return(data.frame())

  hits <- names(ranks) %in% pathway_genes
  weights <- abs(ranks)
  NR <- sum(weights[hits])

  if (!is.finite(NR) || NR == 0) {
    weights <- rep(1, length(ranks))
    NR <- sum(weights[hits])
  }

  running_es <- cumsum(ifelse(hits, weights / NR, 0)) -
    cumsum(ifelse(!hits, 1 / (N - Nh), 0))

  data.frame(
    Position = seq_along(ranks),
    Position_percent = 100 * seq_along(ranks) / N,
    Running_ES = running_es,
    Hit = hits,
    Gene = names(ranks),
    stringsAsFactors = FALSE
  )
}

find_fgsea_result <- function(fgsea_df, comparison, pathway_name) {

  hit <- dplyr::filter(
    fgsea_df,
    .data$Comparison == comparison,
    toupper(.data$pathway) == toupper(pathway_name)
  )

  if (nrow(hit) == 0L) return(NULL)
  hit[1L, , drop = FALSE]
}

make_gsea_curve_plot <- function(
    ranks,
    pathway_genes,
    pathway_name,
    title_text,
    comparison_label,
    fgsea_df
) {
  curve <- calc_enrichment_curve(ranks, pathway_genes)
  if (nrow(curve) == 0) {
    return(make_empty_plot(title_text, "No overlapping genes"))
  }

  result <- find_fgsea_result(fgsea_df, comparison_label, pathway_name)
  nes <- if (is.null(result)) NA_real_ else result$NES[1]
  padj <- if (is.null(result)) NA_real_ else result$padj[1]
  pval <- if (is.null(result)) NA_real_ else result$pval[1]
  family <- if (is.null(result) || !"Gene_set_family" %in% names(result)) {
    NA_character_
  } else {
    as.character(result$Gene_set_family[1])
  }

  is_custom_metaplasia <- !is.na(family) && grepl("CUSTOM_METAPLASIA", family)
  displayed_p <- if (is_custom_metaplasia) pval else padj
  displayed_p_label <- if (is_custom_metaplasia) "Nominal p" else "FDR"

  annotation <- paste0(
    "NES = ", ifelse(is.na(nes), "NA", sprintf("%.2f", nes)),
    "\n", displayed_p_label, " = ",
    ifelse(
      is.na(displayed_p),
      "NA",
      ifelse(
        displayed_p < 1e-4,
        formatC(displayed_p, format = "e", digits = 1),
        sprintf("%.3f", displayed_p)
      )
    )
  )

  ggplot(curve, aes(x = Position_percent, y = Running_ES)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey45") +
    geom_line(
      colour = comparison_curve_cols[[comparison_label]],
      linewidth = 0.75
    ) +
    geom_rug(
      data = curve %>% filter(Hit),
      aes(x = Position_percent),
      sides = "b",
      inherit.aes = FALSE,
      colour = comparison_curve_cols[[comparison_label]],
      linewidth = 0.18
    ) +
    annotate(
      "text",
      x = 98,
      y = Inf,
      label = annotation,
      hjust = 1,
      vjust = 1.25,
      family = font_family,
      size = 5.5 / 2.845276
    ) +
    scale_x_continuous(
      breaks = c(0, 25, 50, 75, 100),
      limits = c(0, 100),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.14))) +
    labs(
      title = format_signature_plot_title(title_text),
      x = "Ranked gene list position (%)",
      y = "Running enrichment score"
    ) +
    nature_theme() +
    theme(legend.position = "none")
}

selected_gsea_pathways <- list(
  "Interferon alpha response" = list(
    pathway = ifna_info$name,
    genes = ifna_info$genes
  ),
  "TNFA/NFKB signalling" = list(
    pathway = tnfa_info$name,
    genes = tnfa_info$genes
  ),
  "Inflammatory response" = list(
    pathway = inflammatory_info$name,
    genes = inflammatory_info$genes
  ),
  "Serrated Lesion Gastric Metaplasia Gene Signature" = list(
    pathway = metaplasia_pathway_id,
    genes = metaplasia_genes
  )
)

curve_plot_sets <- list()

for (comparison_label in c("Infected vs Non-Infected", "BRAF V600E vs WT in SSL_HP")) {
  hallmark_ranks <- if (comparison_label == "Infected vs Non-Infected") {
    infection_rank$ranks
  } else {
    braf_rank$ranks
  }
  metaplasia_ranks <- if (comparison_label == "Infected vs Non-Infected") {
    infection_metaplasia_rank$ranks
  } else {
    braf_metaplasia_rank$ranks
  }
  prefix <- if (comparison_label == "Infected vs Non-Infected") "Infection" else "BRAF_SSL_HP_only"

  plots <- list()

  for (display_name in names(selected_gsea_pathways)) {
    info <- selected_gsea_pathways[[display_name]]
    ranks_for_plot <- if (display_name == "Serrated Lesion Gastric Metaplasia Gene Signature") {
      metaplasia_ranks
    } else {
      hallmark_ranks
    }

    p <- make_gsea_curve_plot(
      ranks = ranks_for_plot,
      pathway_genes = info$genes,
      pathway_name = info$pathway,
      title_text = display_name,
      comparison_label = comparison_label,
      fgsea_df = combined_fgsea
    )

    plots[[display_name]] <- p

    save_plot_set(
      p,
      subdirs$gsea,
      paste0(prefix, "_GSEA_curve_", safe_filename(display_name)),
      90,
      70
    )
  }

  combined <- wrap_plots(plots, ncol = 2) +
    plot_annotation(
      tag_levels = "a",
      title = comparison_label,
      theme = theme(
        plot.title = element_text(family = font_family, face = "bold", size = 8, hjust = 0),
        plot.tag = element_text(family = font_family, face = "bold", size = panel_tag_pt)
      )
    )

  curve_plot_sets[[comparison_label]] <- combined
  save_plot_set(
    combined,
    subdirs$combined,
    paste0(prefix, "_selected_GSEA_curves_combined"),
    183,
    145
  )
}

gene_modules <- list(
  "Interferon alpha response" = ifna_info$genes,
  "Chemokines" = chemokine_genes,
  "Antimicrobial peptides" = antimicrobial_info$genes,
  "TNFA/NFKB signalling" = tnfa_info$genes,
  "Serrated Lesion Gastric Metaplasia Gene Signature" = metaplasia_genes
)

heatmap_comparisons <- list(
  Infection = list(
    comparison = "Infected vs Non-Infected",
    display = "Infected vs Non-Infected",
    selection = "INFECTION_UP"
  ),
  BRAF_SSL_HP_only = list(
    comparison = "BRAF V600E vs WT in SSL_HP",
    display = "BRAF V600E vs WT in SSL/HP",
    selection = "BRAF_ABSOLUTE"
  )
)

make_individual_module_heatmap_data <- function(
    de,
    module_name,
    genes,
    comparison_key,
    comparison_display,
    selection_mode
) {
  requested <- unique(toupper(genes))
  dat <- tibble::tibble(
    Gene = requested,
    Signature_order = seq_along(requested)
  ) %>%
    dplyr::left_join(
      de %>%
        dplyr::select(Gene, logFC, logCPM, PValue, FDR),
      by = "Gene"
    ) %>%
    dplyr::mutate(
      Comparison_key = comparison_key,
      Comparison_display = comparison_display,
      Module = module_name,
      Present_in_DE_table = is.finite(logFC)
    )

  if (module_name == "Serrated Lesion Gastric Metaplasia Gene Signature") {
    selected <- dat %>%
      dplyr::mutate(
        Selected = TRUE,
        Selection_rule = paste0("All genes in ", metaplasia_display_name, "; no logCPM or logFC threshold")
      ) %>%
      dplyr::arrange(Signature_order)
  } else if (selection_mode == "INFECTION_UP") {
    selected <- dat %>%
      dplyr::mutate(
        Selected = is.finite(logCPM) &
          logCPM >= module_heatmap_logcpm_cutoff &
          is.finite(logFC) &
          logFC > module_heatmap_logfc_cutoff,
        Selection_rule = paste0(
          "logCPM >= ", module_heatmap_logcpm_cutoff,
          " and logFC > ", module_heatmap_logfc_cutoff
        )
      ) %>%
      dplyr::filter(Selected) %>%
      dplyr::arrange(dplyr::desc(logFC), Gene) %>%
      dplyr::slice_head(n = module_heatmap_max_genes_per_module)
  } else {
    selected <- dat %>%
      dplyr::mutate(
        Selected = is.finite(logCPM) &
          logCPM >= module_heatmap_logcpm_cutoff &
          is.finite(logFC) &
          abs(logFC) > module_heatmap_logfc_cutoff,
        Selection_rule = paste0(
          "logCPM >= ", module_heatmap_logcpm_cutoff,
          " and abs(logFC) > ", module_heatmap_logfc_cutoff
        )
      ) %>%
      dplyr::filter(Selected) %>%
      dplyr::arrange(dplyr::desc(abs(logFC)), Gene) %>%
      dplyr::slice_head(n = module_heatmap_max_genes_per_module)
  }

  list(candidate = dat, selected = selected)
}

make_individual_module_heatmap_plot <- function(
    selected,
    module_name,
    comparison_display
) {
  title_text <- paste0(format_signature_plot_title(module_name, width = 31), "\n", comparison_display)

  if (nrow(selected) == 0) {
    return(make_empty_plot(title_text, "No genes passed the comparison-specific selection rule"))
  }

  if (module_name == "Serrated Lesion Gastric Metaplasia Gene Signature") {
    gene_order <- rev(selected$Gene)
  } else {
    gene_order <- rev(selected$Gene)
  }

  plot_df <- selected %>%
    dplyr::mutate(
      Gene = factor(Gene, levels = gene_order),
      Comparison_display = factor(Comparison_display, levels = comparison_display),
      logFC_plot = pmax(
        pmin(logFC, module_heatmap_colour_limit),
        -module_heatmap_colour_limit
      ),
      Sig = dplyr::case_when(
        is.finite(FDR) & FDR < 0.001 ~ "***",
        is.finite(FDR) & FDR < 0.01 ~ "**",
        is.finite(FDR) & FDR < 0.05 ~ "*",
        TRUE ~ ""
      )
    )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = Comparison_display, y = Gene, fill = logFC_plot)
  ) +
    ggplot2::geom_tile(
      colour = "white",
      linewidth = 0.25,
      width = 0.94,
      height = 0.82
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = Sig),
      family = font_family,
      size = 5.0 / 2.845276
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-module_heatmap_colour_limit, module_heatmap_colour_limit),
      oob = scales::squish,
      na.value = "grey90",
      name = "log2FC"
    ) +
    ggplot2::labs(title = format_signature_plot_title(title_text), x = NULL, y = NULL) +
    nature_theme() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = font_family,
        face = "bold",
        hjust = 0,
        size = plot_title_pt
      ),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 0,
        hjust = 0.5,
        size = 5.6
      ),
      axis.text.y = ggplot2::element_text(size = 5.5),
      legend.position = "bottom"
    )
}

individual_module_heatmaps <- list()
module_candidate_tables <- list()
module_selected_tables <- list()

for (comparison_key in names(heatmap_comparisons)) {
  settings <- heatmap_comparisons[[comparison_key]]
  de_one <- if (comparison_key == "Infection") infection_de else braf_de
  plots_one_comparison <- list()

  for (module_name in names(gene_modules)) {
    data_one <- make_individual_module_heatmap_data(
      de = de_one,
      module_name = module_name,
      genes = gene_modules[[module_name]],
      comparison_key = comparison_key,
      comparison_display = settings$display,
      selection_mode = settings$selection
    )

    key <- paste(comparison_key, module_name, sep = "__")
    module_candidate_tables[[key]] <- data_one$candidate
    module_selected_tables[[key]] <- data_one$selected

    p <- make_individual_module_heatmap_plot(
      selected = data_one$selected,
      module_name = module_name,
      comparison_display = settings$display
    )

    individual_module_heatmaps[[key]] <- p
    plots_one_comparison[[module_name]] <- p

    n_genes <- nrow(data_one$selected)
    height_mm <- max(68, min(170, 44 + 3.2 * max(n_genes, 1)))

    save_plot_set(
      p,
      subdirs$heatmaps,
      paste0(
        comparison_key,
        "_DE_logFC_heatmap_",
        safe_filename(module_name)
      ),
      82,
      height_mm
    )
  }

  comparison_panel <- patchwork::wrap_plots(
    plots_one_comparison,
    ncol = 2,
    guides = "collect"
  ) &
    ggplot2::theme(legend.position = "bottom")

  comparison_panel <- comparison_panel +
    patchwork::plot_annotation(
      tag_levels = "a",
      title = settings$display,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          family = font_family,
          face = "bold",
          hjust = 0,
          size = 8
        ),
        plot.tag = ggplot2::element_text(
          family = font_family,
          face = "bold",
          size = panel_tag_pt
        )
      )
    )

  n_max <- max(
    vapply(
      module_selected_tables[grepl(paste0("^", comparison_key, "__"), names(module_selected_tables))],
      nrow,
      integer(1)
    ),
    1L
  )
  panel_height_mm <- max(170, min(300, 130 + 4.0 * n_max))

  save_plot_set(
    comparison_panel,
    subdirs$combined,
    paste0(comparison_key, "_DE_module_logFC_heatmaps_combined"),
    183,
    panel_height_mm
  )
}

data.table::fwrite(
  dplyr::bind_rows(module_candidate_tables, .id = "Heatmap_key"),
  file.path(subdirs$tables, "DE_module_heatmap_candidate_gene_table_by_comparison.csv")
)

data.table::fwrite(
  dplyr::bind_rows(module_selected_tables, .id = "Heatmap_key"),
  file.path(subdirs$tables, "DE_module_heatmap_selected_genes_by_comparison.csv")
)

heatmap_selection_rules <- tibble::tribble(
  ~Comparison, ~Module, ~Selection_rule,
  "Infected vs Non-Infected", "Hallmark/chemokine/antimicrobial modules", "logCPM >= 2 and logFC > 1",
  "BRAF V600E vs WT in SSL/HP", "Hallmark/chemokine/antimicrobial modules", "logCPM >= 2 and abs(logFC) > 1",
  "Both comparisons", metaplasia_display_name, "No logCPM or logFC threshold; all requested signature genes displayed"
)

data.table::fwrite(
  heatmap_selection_rules,
  file.path(subdirs$qc, "DE_module_heatmap_selection_rules.csv")
)

make_signature_expression_heatmap <- function(
    expr_obj,
    metadata_subset,
    genes,
    facet_col,
    title_text,
    stem
) {
  present <- intersect(unique(toupper(genes)), rownames(expr_obj$log_expr))
  if (length(present) < 2 || nrow(metadata_subset) < 2) {
    p <- make_empty_plot(title_text, "Insufficient genes or samples")
    save_plot_set(p, subdirs$heatmaps, stem, 120, 80)
    return(p)
  }

  metadata_subset <- metadata_subset %>%
    filter(Sample_ID %in% colnames(expr_obj$log_expr))

  if (facet_col == "Type_Passage") {
    metadata_subset <- metadata_subset %>%
      mutate(Type_Passage = paste0(Type, " | P", Passage)) %>%
      arrange(Type, Passage, Line_ID, Sample_ID)
  } else {
    if (!facet_col %in% names(metadata_subset)) {
      stop("Heatmap facet column is absent: ", facet_col)
    }
    metadata_subset <- metadata_subset %>%
      arrange(.data[[facet_col]], Line_ID, Passage, Sample_ID)
  }

  facet_values <- metadata_subset[[facet_col]]
  facet_levels <- if (is.factor(facet_values)) {
    levels(droplevels(facet_values))
  } else {
    unique(as.character(facet_values))
  }
  metadata_subset$.Facet <- factor(as.character(facet_values), levels = facet_levels)

  mat <- expr_obj$log_expr[present, metadata_subset$Sample_ID, drop = FALSE]
  row_sd <- apply(mat, 1, sd)
  keep <- is.finite(row_sd) & row_sd > 0
  mat <- mat[keep, , drop = FALSE]

  if (nrow(mat) < 2) {
    p <- make_empty_plot(title_text, "Fewer than two variable signature genes")
    save_plot_set(p, subdirs$heatmaps, stem, 120, 80)
    return(p)
  }

  zmat <- t(scale(t(mat)))
  zmat <- pmax(pmin(zmat, 2.5), -2.5)

  gene_order <- rownames(zmat)
  if (nrow(zmat) >= 3) {
    hc <- hclust(dist(zmat), method = "complete")
    gene_order <- rownames(zmat)[hc$order]
  }

  long <- as.data.frame(zmat, check.names = FALSE) %>%
    tibble::rownames_to_column("Gene") %>%
    pivot_longer(-Gene, names_to = "Sample_ID", values_to = "Z") %>%
    left_join(metadata_subset, by = "Sample_ID") %>%
    mutate(
      Gene = factor(Gene, levels = rev(gene_order)),
      Sample_ID = factor(Sample_ID, levels = metadata_subset$Sample_ID)
    )

  p <- ggplot(long, aes(x = Sample_ID, y = Gene, fill = Z)) +
    geom_tile(colour = "white", linewidth = 0.12) +
    facet_grid(
      cols = vars(.Facet),
      scales = "free_x",
      space = "free_x"
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = scales::squish,
      name = "Row z-score"
    ) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_bw(base_family = font_family, base_size = font_pt) +
    theme(
      text = element_text(family = font_family, colour = "black"),
      plot.title = element_text(face = "bold", size = plot_title_pt, hjust = 0),
      panel.grid = element_blank(),
      panel.spacing.x = unit(0.8, "mm"),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.25),
      strip.text = element_text(size = 5.5, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.5),
      axis.text.y = element_text(size = 5.2),
      axis.ticks = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(4, 4, 4, 4, "mm")
    )

  height <- max(85, min(170, 45 + 3.2 * nrow(zmat)))
  save_plot_set(p, subdirs$heatmaps, stem, 183, height)

  data.table::fwrite(
    long,
    file.path(subdirs$tables, paste0(safe_filename(stem), "_long_table.csv"))
  )

  p
}

infection_heatmap_meta <- infection_obj$meta %>%
  filter(Infection_Status %in% c("Non-Infected", "Infected")) %>%
  mutate(Infection_Status = factor(Infection_Status, levels = c("Non-Infected", "Infected")))

p_expr_heatmap_infection <- make_signature_expression_heatmap(
  infection_obj,
  infection_heatmap_meta,
  metaplasia_genes,
  "Infection_Status",
  paste0(metaplasia_display_name, ": infected vs non-infected organoids (sample-level row z-scores)"),
  paste0("Infection_", metaplasia_file_stem, "_sample_level_row_zscore_heatmap")
)

braf_heatmap_meta <- braf_obj$meta %>%
  filter(Type == "SSL_HP", BRAF_Status %in% c("WT", "V600E")) %>%
  mutate(BRAF_Status = factor(BRAF_Status, levels = c("WT", "V600E")))

p_expr_heatmap_braf <- make_signature_expression_heatmap(
  braf_obj,
  braf_heatmap_meta,
  metaplasia_genes,
  "BRAF_Status",
  paste0(metaplasia_display_name, ": BRAF WT vs V600E SSL/HP organoids at initial passage (sample-level row z-scores)"),
  paste0("BRAF_SSL_HP_only_", metaplasia_file_stem, "_sample_level_row_zscore_heatmap")
)

timecourse_heatmap_meta <- timecourse_obj$meta %>%
  filter(Type %in% c("SSL_HP", "H_ADJ_SSL_HP"), is.finite(Passage)) %>%
  mutate(Type_Passage = paste0(Type, " | P", Passage))

p_expr_heatmap_timecourse <- make_signature_expression_heatmap(
  timecourse_obj,
  timecourse_heatmap_meta,
  metaplasia_genes,
  "Type_Passage",
  paste0(metaplasia_display_name, ": organoid passage series (sample-level row z-scores)"),
  paste0("Timecourse_", metaplasia_file_stem, "_sample_level_row_zscore_heatmap")
)

summary_top <- (p_infection_pca + theme(legend.position = "bottom")) |
  (p_braf_pca + theme(legend.position = "bottom"))

summary_middle <- p_bar_infection | p_bar_braf
summary_bottom <- p_hallmark_heatmap_top

chapter5_summary <- summary_top / summary_middle / summary_bottom +
  plot_layout(heights = c(0.8, 1.2, 1.0)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(family = font_family, face = "bold", size = panel_tag_pt))
  )

save_plot_set(chapter5_summary, subdirs$combined, "Chapter5_organoid_analysis_summary", 183, 255)

all_scale_qc <- bind_rows(
  infection_obj$scale_qc,
  timecourse_obj$scale_qc,
  braf_obj$scale_qc,
  braf_p1_obj$scale_qc
)

data.table::fwrite(
  all_scale_qc,
  file.path(subdirs$qc, "ALL_expression_scale_QC_summary.csv")
)

writeLines(
  c(
    paste0("Run date: ", Sys.time()),
    paste0("Base directory: ", base_dir),
    paste0("Output directory: ", out_dir),
    "",
    "PCA:",
    paste0("  Variable genes: top ", pca_n_variable_genes, " by variance"),
    "  Centred: TRUE",
    "  Unit-variance scaled: FALSE",
    "  CPM transformation: log2(CPM + 1)",
    "  CPM import: missing Gene_ID/shifted sample headers are reconstructed when detected.",
    "  Infection metadata: Sample_ID infection token reconciled with Infection_Status before matching.",
    "  Infection PCA: samples are coloured by infection status and labelled by passage using non-overlapping labels.",
    "",
    "BRAF restrictions:",
    paste0("  PCA and sample-level scores use: ", basename(files$timecourse_cpm), "."),
    "  Primary sample-level set: one initial sample per trajectory; P1 where present, otherwise P2.",
    "  Strict-P1 PCA, scores, MUC5AC and CPM heatmaps are exported as a sensitivity analysis.",
    "  No trajectory contributes both P1 and P2 to either sample-level analysis.",
    "  Missing passage-level BRAF labels are filled only from a non-conflicting trajectory-level consensus.",
    "  PCA: SSL_HP WT, SSL_HP V600E and H_ADJ_SSL_HP included.",
    "  Scores, DE and GSEA: SSL/HP only; H_ADJ excluded.",
    paste0("  DE/GSEA source: ", braf_de_source),
    "  The DE summary table does not contain its sample manifest; its exact upstream P1/P2 membership must be checked from the edgeR design/sample list.",
    "",
    "GSEA:",
    paste0("  Infection ranking: ", infection_rank$method),
    paste0("  BRAF ranking: ", braf_rank$method),
    paste0("  Hallmark GSEA logCPM filter: >= ", gsea_logcpm_cutoff),
    paste0("  ", metaplasia_display_name, " GSEA logCPM filter: none (ALL_VALID_GENES)."),
    paste0("  minSize: ", fgsea_min_size),
    paste0("  maxSize: ", fgsea_max_size),
    paste0("  eps: ", fgsea_eps),
    "  Hallmark FDR is calculated within the Hallmark family for each comparison.",
    paste0("  ", metaplasia_display_name, " is analysed separately as the custom set."),
    "",
    "Signature scoring:",
    "  Method: singscore, centred single-direction score.",
    paste0("  Hallmark genes require mean logCPM >= ", signature_mean_logcpm_cutoff, "."),
    paste0("  ", metaplasia_display_name, " is never mean-logCPM filtered."),
    paste0("  Metaplasia signature source: ", metaplasia_signature_file),
    paste0("  Metaplasia signature genes: ", length(metaplasia_genes)),
    "",
    "Outputs are exported as PDF, PNG, TIFF and JPG at 600 dpi for raster formats."
  ),
  file.path(out_dir, "RUN_NOTES.txt")
)

package_versions <- data.frame(
  Package = c(
    "R", "data.table", "dplyr", "ggplot2", "fgsea",
    "msigdbr", "singscore", "org.Hs.eg.db"
  ),
  Version = c(
    paste0(R.version$major, ".", R.version$minor),
    as.character(packageVersion("data.table")),
    as.character(packageVersion("dplyr")),
    as.character(packageVersion("ggplot2")),
    as.character(packageVersion("fgsea")),
    as.character(packageVersion("msigdbr")),
    as.character(packageVersion("singscore")),
    as.character(packageVersion("org.Hs.eg.db"))
  ),
  stringsAsFactors = FALSE
)

data.table::fwrite(package_versions, file.path(subdirs$qc, "PACKAGE_VERSIONS.csv"))
capture.output(sessionInfo(), file = file.path(subdirs$qc, "SESSION_INFO.txt"))

message("\nDONE.")
message("Output directory: ", out_dir)
message("Inspect first:")
message("  01_QC/ALL_expression_scale_QC_summary.csv")
message("  02_PCA/")
message("  03_Signature_scores/")
message("  04_GSEA/")
message("  05_Heatmaps/")
message("  06_Combined_figures/")

message("\nRunning post-run statistical audit and figure corrections...")

type_display_labels <- c(
  SSL_HP = "SSL/HP",
  H_ADJ_SSL_HP = "H-ADJ SSL/HP"
)

audited_passage_shapes <- c(21, 22, 24, 23, 25, 8, 4, 3, 7, 0, 1, 2, 5)
if (length(passage_levels) > length(audited_passage_shapes)) {
  stop("Too many distinct passages for the audited PCA shape palette.")
}
audited_passage_shape_values <- stats::setNames(
  audited_passage_shapes[seq_along(passage_levels)],
  passage_labels
)

p_time_pca <- ggplot2::ggplot(pca_time$scores, ggplot2::aes(x = PC1, y = PC2)) +
  ggplot2::geom_path(
    data = trajectory_lines,
    ggplot2::aes(group = Trajectory_ID),
    colour = "grey74",
    linewidth = 0.32,
    alpha = 0.75
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      colour = Type,
      fill = Type,
      shape = Passage_factor
    ),
    size = 2.35,
    stroke = 0.48,
    alpha = 0.96
  ) +
  ggplot2::scale_colour_manual(
    values = paper_cols[c("SSL_HP", "H_ADJ_SSL_HP")],
    breaks = c("SSL_HP", "H_ADJ_SSL_HP"),
    labels = type_display_labels[c("SSL_HP", "H_ADJ_SSL_HP")],
    drop = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = paper_cols[c("SSL_HP", "H_ADJ_SSL_HP")],
    breaks = c("SSL_HP", "H_ADJ_SSL_HP"),
    labels = type_display_labels[c("SSL_HP", "H_ADJ_SSL_HP")],
    drop = FALSE,
    guide = "none"
  ) +
  ggplot2::scale_shape_manual(
    values = audited_passage_shape_values,
    breaks = passage_labels,
    drop = FALSE
  ) +
  ggplot2::labs(
    x = sprintf("PC1 (%.1f%%)", pca_time$variance_explained[1]),
    y = sprintf("PC2 (%.1f%%)", pca_time$variance_explained[2]),
    colour = "Sample type",
    shape = "Passage"
  ) +
  nature_theme() +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.justification = "left",
    legend.spacing.y = grid::unit(1.0, "mm"),
    legend.margin = ggplot2::margin(t = 2, r = 1, b = 0, l = 1, unit = "mm")
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(
      order = 1,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(shape = 16, size = 2.4, alpha = 1)
    ),
    shape = ggplot2::guide_legend(
      order = 2,
      nrow = 2,
      byrow = TRUE,
      override.aes = list(colour = "black", fill = "white", size = 2.25, alpha = 1)
    )
  )

save_plot_set(
  p_time_pca,
  subdirs$pca,
  "PCA_2_Timecourse_SSL_HP_and_H_ADJ",
  125,
  92
)

pca_combined <-
  (p_infection_pca + ggplot2::theme(legend.position = "bottom")) |
  p_time_pca |
  (p_braf_pca + ggplot2::theme(legend.position = "bottom"))

pca_combined <- pca_combined +
  patchwork::plot_annotation(
    tag_levels = "a",
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(
        family = font_family,
        face = "bold",
        size = panel_tag_pt
      )
    )
  )

save_plot_set(
  pca_combined,
  subdirs$combined,
  "Combined_three_PCA_plots",
  183,
  100
)

# Sensitivity tests and longitudinal mixed models ----
safe_wilcox_p <- function(...) {
  suppressWarnings(
    tryCatch(
      stats::wilcox.test(...)$p.value,
      error = function(e) NA_real_
    )
  )
}

run_two_group_test_audited <- function(
    df,
    group_col,
    levels,
    pair_col = NULL,
    prefer_paired = FALSE
) {
  dat <- df %>%
    dplyr::filter(
      !is.na(.data[[group_col]]),
      .data[[group_col]] %in% levels,
      is.finite(Score)
    ) %>%
    dplyr::mutate(.Group = factor(.data[[group_col]], levels = levels))

  paired_used <- FALSE
  n_pairs <- 0L
  line_df <- data.frame()
  p_exact <- NA_real_
  p_asym_no_correction <- NA_real_
  p_asym_continuity_corrected <- NA_real_
  p_paired_t <- NA_real_
  paired_t_ci_lower <- NA_real_
  paired_t_ci_upper <- NA_real_
  exact_available <- FALSE
  mean_difference <- NA_real_
  median_difference <- NA_real_
  n_higher <- NA_integer_
  n_lower <- NA_integer_
  n_equal <- NA_integer_
  n_nonzero <- NA_integer_
  group_1_n <- sum(dat$.Group == levels[1], na.rm = TRUE)
  group_2_n <- sum(dat$.Group == levels[2], na.rm = TRUE)

  if (prefer_paired && !is.null(pair_col) && pair_col %in% names(dat)) {
    complete_pairs <- dat %>%
      dplyr::filter(!is.na(.data[[pair_col]]), .data[[pair_col]] != "") %>%
      dplyr::group_by(.data[[pair_col]]) %>%
      dplyr::filter(dplyr::n_distinct(.Group) == 2, dplyr::n() == 2) %>%
      dplyr::ungroup()

    n_pairs <- dplyr::n_distinct(complete_pairs[[pair_col]])

    if (n_pairs >= 2) {
      wide <- complete_pairs %>%
        dplyr::select(dplyr::all_of(pair_col), .Group, Score) %>%
        tidyr::pivot_wider(names_from = .Group, values_from = Score)

      x_level_2 <- wide[[levels[2]]]
      x_level_1 <- wide[[levels[1]]]
      differences <- x_level_2 - x_level_1
      finite_differences <- differences[is.finite(differences)]
      nonzero_differences <- finite_differences[finite_differences != 0]

      mean_difference <- mean(finite_differences)
      median_difference <- stats::median(finite_differences)
      n_higher <- sum(finite_differences > 0)
      n_lower <- sum(finite_differences < 0)
      n_equal <- sum(finite_differences == 0)
      n_nonzero <- length(nonzero_differences)

      exact_available <-
        n_nonzero > 0 &&
        n_equal == 0 &&
        !anyDuplicated(abs(nonzero_differences)) &&
        n_nonzero <= 50

      if (exact_available) {
        p_exact <- safe_wilcox_p(
          x_level_2,
          x_level_1,
          paired = TRUE,
          exact = TRUE,
          correct = FALSE,
          alternative = "two.sided"
        )
      }

      p_asym_no_correction <- safe_wilcox_p(
        x_level_2,
        x_level_1,
        paired = TRUE,
        exact = FALSE,
        correct = FALSE,
        alternative = "two.sided"
      )

      p_asym_continuity_corrected <- safe_wilcox_p(
        x_level_2,
        x_level_1,
        paired = TRUE,
        exact = FALSE,
        correct = TRUE,
        alternative = "two.sided"
      )

      paired_t_result <- tryCatch(
        stats::t.test(
          x_level_2,
          x_level_1,
          paired = TRUE,
          alternative = "two.sided",
          conf.level = 0.95
        ),
        error = function(e) NULL
      )

      if (!is.null(paired_t_result)) {
        p_paired_t <- paired_t_result$p.value
        paired_t_ci_lower <- unname(paired_t_result$conf.int[1])
        paired_t_ci_upper <- unname(paired_t_result$conf.int[2])
      }

      paired_used <- TRUE
      line_df <- complete_pairs
    }
  }

  if (!paired_used) {
    x_level_1 <- dat$Score[dat$.Group == levels[1]]
    x_level_2 <- dat$Score[dat$.Group == levels[2]]

    if (length(x_level_1) >= 2 && length(x_level_2) >= 2) {
      mean_difference <- mean(x_level_2) - mean(x_level_1)
      median_difference <- stats::median(x_level_2) - stats::median(x_level_1)

      exact_available <-
        !anyDuplicated(c(x_level_1, x_level_2)) &&
        length(x_level_1) * length(x_level_2) <= 200

      if (exact_available) {
        p_exact <- safe_wilcox_p(
          x_level_2,
          x_level_1,
          paired = FALSE,
          exact = TRUE,
          correct = FALSE,
          alternative = "two.sided"
        )
      }

      p_asym_no_correction <- safe_wilcox_p(
        x_level_2,
        x_level_1,
        paired = FALSE,
        exact = FALSE,
        correct = FALSE,
        alternative = "two.sided"
      )

      p_asym_continuity_corrected <- safe_wilcox_p(
        x_level_2,
        x_level_1,
        paired = FALSE,
        exact = FALSE,
        correct = TRUE,
        alternative = "two.sided"
      )
    }
  }

  if (paired_used) {
    if (infection_score_primary_test == "PAIRED_T" && is.finite(p_paired_t)) {
      p_primary <- p_paired_t
      primary_method <- "Two-sided paired t-test"
    } else if (is.finite(p_exact)) {
      p_primary <- p_exact
      primary_method <- "Exact two-sided paired Wilcoxon signed-rank"
    } else {
      p_primary <- p_asym_no_correction
      primary_method <- "Asymptotic paired Wilcoxon signed-rank without continuity correction"
    }
  } else {
    p_primary <- if (is.finite(p_exact)) p_exact else p_asym_no_correction
    primary_method <- if (is.finite(p_exact)) {
      "Exact two-sided Wilcoxon rank-sum"
    } else {
      "Asymptotic Wilcoxon rank-sum without continuity correction"
    }
  }

  list(
    data = dat,
    p_value = p_primary,
    p_value_exact = p_exact,
    p_value_asymptotic_no_correction = p_asym_no_correction,
    p_value_asymptotic_continuity_corrected = p_asym_continuity_corrected,
    p_value_paired_t = p_paired_t,
    paired_t_CI_lower = paired_t_ci_lower,
    paired_t_CI_upper = paired_t_ci_upper,
    exact_available = exact_available,
    primary_method = primary_method,
    paired = paired_used,
    n_pairs = n_pairs,
    line_data = line_df,
    group_1_n = group_1_n,
    group_2_n = group_2_n,
    mean_difference_level_2_minus_level_1 = mean_difference,
    median_difference_level_2_minus_level_1 = median_difference,
    pairs_higher_in_level_2 = n_higher,
    pairs_lower_in_level_2 = n_lower,
    pairs_equal = n_equal,
    nonzero_pairs = n_nonzero
  )
}

make_two_group_score_plot_audited <- function(
    df,
    signature_name,
    group_col,
    group_levels,
    colours,
    x_labels = group_levels,
    pair_col = NULL,
    prefer_paired = FALSE,
    y_label = "Singscore"
) {
  dat <- df %>% dplyr::filter(Signature == signature_name)
  test <- run_two_group_test_audited(
    dat,
    group_col = group_col,
    levels = group_levels,
    pair_col = pair_col,
    prefer_paired = prefer_paired
  )

  plot_df <- test$data
  if (nrow(plot_df) == 0) {
    return(make_empty_plot(signature_name, "No valid samples"))
  }

  plot_df$.Group <- factor(plot_df$.Group, levels = group_levels)
  ymax <- max(plot_df$Score, na.rm = TRUE)
  ymin <- min(plot_df$Score, na.rm = TRUE)
  span <- max(ymax - ymin, 0.05)
  label_y <- ymax + 0.17 * span

  p_label <- format_p(test$p_value)
  if (test$paired && is.finite(test$nonzero_pairs)) {
    if (test$pairs_higher_in_level_2 > test$pairs_lower_in_level_2) {
      direction_label <- paste0(
        test$pairs_higher_in_level_2,
        "/",
        test$nonzero_pairs,
        " pairs higher in ",
        x_labels[2]
      )
    } else if (test$pairs_lower_in_level_2 > test$pairs_higher_in_level_2) {
      direction_label <- paste0(
        test$pairs_lower_in_level_2,
        "/",
        test$nonzero_pairs,
        " pairs lower in ",
        x_labels[2]
      )
    } else {
      direction_label <- paste0(test$nonzero_pairs, " non-zero pairs")
    }
    method_short <- if (grepl("paired t-test", test$primary_method, ignore.case = TRUE)) {
      "paired t-test"
    } else if (grepl("Exact", test$primary_method, fixed = TRUE)) {
      "exact paired Wilcoxon"
    } else {
      "paired Wilcoxon"
    }
    p_label <- paste0(direction_label, "
", method_short, " ", format_p(test$p_value))
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .Group, y = Score, fill = .Group))

  if (test$paired && nrow(test$line_data) > 0) {
    line_data <- test$line_data %>%
      dplyr::mutate(.Group = factor(.data[[group_col]], levels = group_levels))

    p <- p +
      ggplot2::geom_line(
        data = line_data,
        ggplot2::aes(x = .Group, y = Score, group = .data[[pair_col]]),
        inherit.aes = FALSE,
        colour = "grey72",
        linewidth = 0.35,
        alpha = 0.8
      )
  }

  p +
    ggplot2::geom_boxplot(
      width = 0.55,
      outlier.shape = NA,
      linewidth = 0.35,
      alpha = 0.45,
      colour = "black"
    ) +
    ggplot2::geom_point(
      shape = 21,
      size = 2.0,
      stroke = 0.30,
      colour = "black",
      position = ggplot2::position_jitter(width = 0.055, height = 0),
      alpha = 0.95
    ) +
    ggplot2::annotate(
      "segment",
      x = 1,
      xend = 2,
      y = label_y,
      yend = label_y,
      linewidth = 0.3,
      colour = "black"
    ) +
    ggplot2::annotate(
      "text",
      x = 1.5,
      y = label_y + 0.035 * span,
      label = p_label,
      family = font_family,
      size = 5.4 / 2.845276,
      lineheight = 0.95
    ) +
    ggplot2::scale_fill_manual(values = colours, breaks = group_levels, drop = FALSE) +
    ggplot2::scale_x_discrete(labels = stats::setNames(x_labels, group_levels), drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.30))) +
    ggplot2::labs(
      title = format_signature_plot_title(signature_name),
      x = NULL,
      y = y_label,
      fill = NULL
    ) +
    nature_theme() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(size = 6)
    )
}

infection_score_plots_audited <- list()
infection_score_stats_audited <- list()

for (sig in names(signature_list)) {
  p <- make_two_group_score_plot_audited(
    infection_scores$scores,
    signature_name = sig,
    group_col = "Infection_Status",
    group_levels = c("Non-Infected", "Infected"),
    colours = infection_cols,
    x_labels = c("Non-Infected", "Infected"),
    pair_col = "Pair_ID",
    prefer_paired = TRUE
  )

  infection_score_plots_audited[[sig]] <- p
  save_plot_set(
    p,
    subdirs$scores,
    paste0("Infection_", signature_file_stems[[sig]], "_singscore"),
    68,
    76
  )

  test <- run_two_group_test_audited(
    infection_scores$scores %>% dplyr::filter(Signature == sig),
    "Infection_Status",
    c("Non-Infected", "Infected"),
    "Pair_ID",
    TRUE
  )

  infection_score_stats_audited[[sig]] <- data.frame(
    Dataset = "Infection",
    Signature = sig,
    Level_1 = "Non-Infected",
    Level_2 = "Infected",
    Test_primary = test$primary_method,
    Complete_pairs = test$n_pairs,
    P_value = test$p_value,
    P_value_exact = test$p_value_exact,
    P_value_asymptotic_no_correction = test$p_value_asymptotic_no_correction,
    P_value_asymptotic_continuity_corrected = test$p_value_asymptotic_continuity_corrected,
    P_value_paired_t = test$p_value_paired_t,
    Paired_t_95CI_lower = test$paired_t_CI_lower,
    Paired_t_95CI_upper = test$paired_t_CI_upper,
    Mean_difference_level_2_minus_level_1 = test$mean_difference_level_2_minus_level_1,
    Median_difference_level_2_minus_level_1 = test$median_difference_level_2_minus_level_1,
    Pairs_higher_in_level_2 = test$pairs_higher_in_level_2,
    Pairs_lower_in_level_2 = test$pairs_lower_in_level_2,
    Pairs_equal = test$pairs_equal,
    stringsAsFactors = FALSE
  )
}

infection_score_stats_audited_df <- dplyr::bind_rows(infection_score_stats_audited) %>%
  dplyr::mutate(
    P_value_BH_across_four_prespecified_scores = stats::p.adjust(P_value, method = "BH")
  )
data.table::fwrite(
  infection_score_stats_audited_df,
  file.path(subdirs$tables, "Infection_signature_score_statistics.csv")
)

infection_score_combined <- patchwork::wrap_plots(infection_score_plots_audited, ncol = 4) +
  patchwork::plot_annotation(
    tag_levels = "a",
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(
        family = font_family,
        face = "bold",
        size = panel_tag_pt
      )
    )
  )

save_plot_set(
  infection_score_combined,
  subdirs$combined,
  "Infection_signature_scores_combined",
  183,
  78
)

braf_score_plots_audited <- list()
braf_score_stats_audited <- list()

for (sig in names(signature_list)) {
  p <- make_two_group_score_plot_audited(
    braf_score_input,
    signature_name = sig,
    group_col = "BRAF_Status",
    group_levels = c("WT", "V600E"),
    colours = braf_cols,
    x_labels = c("WT", "V600E"),
    pair_col = NULL,
    prefer_paired = FALSE
  )

  braf_score_plots_audited[[sig]] <- p
  save_plot_set(
    p,
    subdirs$scores,
    paste0("BRAF_SSL_HP_only_", signature_file_stems[[sig]], "_singscore"),
    64,
    72
  )

  test <- run_two_group_test_audited(
    braf_score_input %>% dplyr::filter(Signature == sig),
    "BRAF_Status",
    c("WT", "V600E"),
    NULL,
    FALSE
  )

  braf_score_stats_audited[[sig]] <- data.frame(
    Dataset = "BRAF_SSL_HP_only",
    Signature = sig,
    Level_1 = "WT",
    Level_2 = "V600E",
    Test_primary = test$primary_method,
    N_WT = test$group_1_n,
    N_V600E = test$group_2_n,
    P_value = test$p_value,
    P_value_exact = test$p_value_exact,
    P_value_asymptotic_no_correction = test$p_value_asymptotic_no_correction,
    P_value_asymptotic_continuity_corrected = test$p_value_asymptotic_continuity_corrected,
    Mean_difference_V600E_minus_WT = test$mean_difference_level_2_minus_level_1,
    Median_difference_V600E_minus_WT = test$median_difference_level_2_minus_level_1,
    stringsAsFactors = FALSE
  )
}

braf_score_stats_audited_df <- dplyr::bind_rows(braf_score_stats_audited)
data.table::fwrite(
  braf_score_stats_audited_df,
  file.path(subdirs$tables, "BRAF_SSL_HP_only_signature_score_statistics.csv")
)

braf_score_combined <- patchwork::wrap_plots(braf_score_plots_audited, ncol = 4) +
  patchwork::plot_annotation(
    tag_levels = "a",
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(
        family = font_family,
        face = "bold",
        size = panel_tag_pt
      )
    )
  )

save_plot_set(
  braf_score_combined,
  subdirs$combined,
  "BRAF_SSL_HP_only_signature_scores_combined",
  183,
  74
)

derive_patient_id_audited <- function(x) {
  out <- stringr::str_extract(toupper(as.character(x)), "^[A-Z]+_[A-Z0-9]+")
  out[is.na(out) | out == ""] <- toupper(as.character(x))[is.na(out) | out == ""]
  out
}

braf_sample_audit <- braf_score_input %>%
  dplyr::mutate(Patient_ID = derive_patient_id_audited(Line_ID)) %>%
  dplyr::distinct(
    Sample_ID,
    Patient_ID,
    Line_ID,
    Trajectory_ID,
    Passage,
    Type,
    BRAF_Status
  ) %>%
  dplyr::arrange(BRAF_Status, Patient_ID, Sample_ID)

data.table::fwrite(
  braf_sample_audit,
  file.path(subdirs$qc, "BRAF_SSL_HP_only_sample_independence_and_passage_audit.csv")
)

braf_patient_collapsed <- braf_score_input %>%
  dplyr::mutate(Patient_ID = derive_patient_id_audited(Line_ID)) %>%
  dplyr::group_by(Signature, Patient_ID, BRAF_Status) %>%
  dplyr::summarise(
    Score = mean(Score, na.rm = TRUE),
    N_organoid_lines = dplyr::n_distinct(Line_ID),
    Passages = paste(sort(unique(Passage)), collapse = ";"),
    .groups = "drop"
  )

patient_collapsed_stats <- purrr::map_dfr(names(signature_list), function(sig) {
  test <- run_two_group_test_audited(
    braf_patient_collapsed %>% dplyr::filter(Signature == sig),
    "BRAF_Status",
    c("WT", "V600E"),
    NULL,
    FALSE
  )

  data.frame(
    Analysis = "Patient-collapsed sensitivity",
    Signature = sig,
    N_WT_patients = test$group_1_n,
    N_V600E_patients = test$group_2_n,
    P_value = test$p_value,
    P_value_exact = test$p_value_exact,
    Mean_difference_V600E_minus_WT = test$mean_difference_level_2_minus_level_1,
    Median_difference_V600E_minus_WT = test$median_difference_level_2_minus_level_1,
    stringsAsFactors = FALSE
  )
})

data.table::fwrite(
  braf_patient_collapsed,
  file.path(subdirs$tables, "BRAF_SSL_HP_only_patient_collapsed_signature_scores.csv")
)
data.table::fwrite(
  patient_collapsed_stats,
  file.path(subdirs$tables, "BRAF_SSL_HP_only_patient_collapsed_signature_statistics.csv")
)

if (!requireNamespace("nlme", quietly = TRUE)) {
  if (!install_missing_packages) {
    stop("Package 'nlme' is required for the time-course mixed-effects audit.")
  }
  install.packages("nlme", repos = "https://cloud.r-project.org")
}

fit_timecourse_slope_by_type <- function(dat, signature_name, type_name) {
  x <- dat %>%
    dplyr::filter(
      Signature == signature_name,
      Type == type_name,
      is.finite(Passage),
      is.finite(Score),
      !is.na(Trajectory_ID),
      Trajectory_ID != ""
    ) %>%
    dplyr::mutate(
      Passage_centered = Passage - min(Passage, na.rm = TRUE),
      Trajectory_ID = factor(Trajectory_ID)
    )

  n_trajectories <- dplyr::n_distinct(x$Trajectory_ID)

  if (nrow(x) < 6 || n_trajectories < 3 || dplyr::n_distinct(x$Passage) < 2) {
    return(data.frame(
      Signature = signature_name,
      Type = type_name,
      N_samples = nrow(x),
      N_trajectories = n_trajectories,
      Term = "Passage_centered",
      Estimate = NA_real_,
      SE = NA_real_,
      DF = NA_real_,
      t_value = NA_real_,
      P_value = NA_real_,
      Model_status = "Insufficient repeated-measures data",
      stringsAsFactors = FALSE
    ))
  }

  fit <- tryCatch(
    nlme::lme(
      fixed = Score ~ Passage_centered,
      random = ~1 | Trajectory_ID,
      data = x,
      method = "REML",
      na.action = stats::na.omit,
      control = nlme::lmeControl(returnObject = TRUE, opt = "optim")
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(data.frame(
      Signature = signature_name,
      Type = type_name,
      N_samples = nrow(x),
      N_trajectories = n_trajectories,
      Term = "Passage_centered",
      Estimate = NA_real_,
      SE = NA_real_,
      DF = NA_real_,
      t_value = NA_real_,
      P_value = NA_real_,
      Model_status = paste0("Model failed: ", conditionMessage(fit)),
      stringsAsFactors = FALSE
    ))
  }

  tt <- summary(fit)$tTable
  term <- "Passage_centered"

  data.frame(
    Signature = signature_name,
    Type = type_name,
    N_samples = nrow(x),
    N_trajectories = n_trajectories,
    Term = term,
    Estimate = unname(tt[term, "Value"]),
    SE = unname(tt[term, "Std.Error"]),
    DF = unname(tt[term, "DF"]),
    t_value = unname(tt[term, "t-value"]),
    P_value = unname(tt[term, "p-value"]),
    Model_status = "OK: random-intercept mixed model",
    stringsAsFactors = FALSE
  )
}

fit_timecourse_interaction <- function(dat, signature_name) {
  x <- dat %>%
    dplyr::filter(
      Signature == signature_name,
      Type %in% c("H_ADJ_SSL_HP", "SSL_HP"),
      is.finite(Passage),
      is.finite(Score),
      !is.na(Trajectory_ID),
      Trajectory_ID != ""
    ) %>%
    dplyr::mutate(
      Passage_centered = Passage - min(Passage, na.rm = TRUE),
      Type = factor(Type, levels = c("H_ADJ_SSL_HP", "SSL_HP")),
      Trajectory_ID = factor(Trajectory_ID)
    )

  n_trajectories <- dplyr::n_distinct(x$Trajectory_ID)

  fit <- tryCatch(
    nlme::lme(
      fixed = Score ~ Passage_centered * Type,
      random = ~1 | Trajectory_ID,
      data = x,
      method = "REML",
      na.action = stats::na.omit,
      control = nlme::lmeControl(returnObject = TRUE, opt = "optim")
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(data.frame(
      Signature = signature_name,
      N_samples = nrow(x),
      N_trajectories = n_trajectories,
      Term = NA_character_,
      Estimate = NA_real_,
      SE = NA_real_,
      DF = NA_real_,
      t_value = NA_real_,
      P_value = NA_real_,
      Model_status = paste0("Model failed: ", conditionMessage(fit)),
      stringsAsFactors = FALSE
    ))
  }

  tt <- as.data.frame(summary(fit)$tTable)
  tt$Term <- rownames(tt)
  rownames(tt) <- NULL

  data.frame(
    Signature = signature_name,
    N_samples = nrow(x),
    N_trajectories = n_trajectories,
    Term = tt$Term,
    Estimate = tt$Value,
    SE = tt$Std.Error,
    DF = tt$DF,
    t_value = tt$`t-value`,
    P_value = tt$`p-value`,
    Model_status = "OK: random-intercept mixed model",
    stringsAsFactors = FALSE
  )
}

timecourse_mixed_by_type <- purrr::map_dfr(names(signature_list), function(sig) {
  purrr::map_dfr(c("H_ADJ_SSL_HP", "SSL_HP"), function(type_name) {
    fit_timecourse_slope_by_type(timecourse_scores$scores, sig, type_name)
  })
})

timecourse_mixed_interaction <- purrr::map_dfr(names(signature_list), function(sig) {
  fit_timecourse_interaction(timecourse_scores$scores, sig)
})

data.table::fwrite(
  timecourse_mixed_by_type,
  file.path(subdirs$tables, "Timecourse_signature_score_mixed_models_by_type.csv")
)
data.table::fwrite(
  timecourse_mixed_interaction,
  file.path(subdirs$tables, "Timecourse_signature_score_mixed_model_type_interactions.csv")
)

timecourse_spearman_descriptive <- data.table::fread(
  file.path(subdirs$tables, "Timecourse_signature_score_Spearman_statistics.csv"),
  data.table = FALSE
)
timecourse_spearman_descriptive$Inference_note <- paste0(
  "Descriptive only: samples are repeated within organoid trajectories; ",
  "use the mixed-effects model tables for inferential p values."
)
data.table::fwrite(
  timecourse_spearman_descriptive,
  file.path(subdirs$tables, "Timecourse_signature_score_Spearman_statistics.csv")
)

infection_rank_all <- infection_metaplasia_rank
braf_rank_all <- braf_metaplasia_rank

custom_metaplasia_fgsea_audited <- dplyr::bind_rows(
  infection_metaplasia_fgsea,
  braf_metaplasia_fgsea
) %>%
  dplyr::mutate(
    Filter_mode = "ALL_VALID_GENES",
    P_value_interpretation = paste0(
      "Nominal custom-set p value. Only one custom pathway was tested, ",
      "so padj equals pval and is not comparable with Hallmark-family FDR."
    )
  )

hallmark_fgsea_audited <- dplyr::bind_rows(
  infection_hallmark_fgsea,
  braf_hallmark_fgsea
) %>%
  dplyr::mutate(
    Filter_mode = "LOGCPM_GE_2",
    N_ranked_genes = ifelse(
      Comparison == "Infected vs Non-Infected",
      length(infection_rank$ranks),
      length(braf_rank$ranks)
    ),
    P_value_interpretation = "BH FDR within the Hallmark family for this comparison."
  )

combined_fgsea_audited <- dplyr::bind_rows(
  hallmark_fgsea_audited,
  custom_metaplasia_fgsea_audited
)

data.table::fwrite(
  combined_fgsea_audited,
  file.path(subdirs$tables, "Organoid_combined_FGSEA_results.csv")
)

data.table::fwrite(
  hallmark_fgsea_audited,
  file.path(subdirs$tables, "Organoid_Hallmark_FGSEA_results.csv")
)

data.table::fwrite(
  custom_metaplasia_fgsea_audited,
  file.path(subdirs$tables, "Organoid_custom_metaplasia_FGSEA_ALL_VALID_GENES.csv")
)

make_metaplasia_gsea_gene_audit <- function(
    de,
    comparison_label,
    rank_all,
    custom_results
) {
  result_all <- custom_results %>%
    dplyr::filter(Comparison == comparison_label) %>%
    dplyr::slice(1)

  leading_all <- if (nrow(result_all) == 1 && !is.na(result_all$leadingEdge[1])) {
    strsplit(result_all$leadingEdge[1], ";", fixed = TRUE)[[1]]
  } else {
    character(0)
  }

  data.frame(
    Gene = metaplasia_genes,
    Signature_order = seq_along(metaplasia_genes),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(
      de %>%
        dplyr::select(Gene, logFC, logCPM, LR, PValue, FDR),
      by = "Gene"
    ) %>%
    dplyr::mutate(
      Comparison = comparison_label,
      Present_in_DE_table = !is.na(logFC),
      LogCPM_filter_applied = FALSE,
      Rank_metric_ALL_VALID_GENES = unname(rank_all$ranks[Gene]),
      Leading_edge_ALL_VALID_GENES = Gene %in% leading_all
    ) %>%
    dplyr::arrange(Signature_order)
}

infection_metaplasia_gene_audit <- make_metaplasia_gsea_gene_audit(
  infection_de,
  "Infected vs Non-Infected",
  infection_rank_all,
  custom_metaplasia_fgsea_audited
)

braf_metaplasia_gene_audit <- make_metaplasia_gsea_gene_audit(
  braf_de,
  "BRAF V600E vs WT in SSL_HP",
  braf_rank_all,
  custom_metaplasia_fgsea_audited
)

data.table::fwrite(
  infection_metaplasia_gene_audit,
  file.path(subdirs$qc, "Infection_custom_metaplasia_GSEA_24_gene_audit.csv")
)
data.table::fwrite(
  braf_metaplasia_gene_audit,
  file.path(subdirs$qc, "BRAF_custom_metaplasia_GSEA_24_gene_audit.csv")
)

custom_metaplasia_summary <- custom_metaplasia_fgsea_audited %>%
  dplyr::transmute(
    Comparison,
    Filter_mode,
    Requested_signature_genes = length(metaplasia_genes),
    Genes_used_by_fgsea = size,
    N_ranked_genes,
    ES,
    NES,
    Nominal_P_value = pval,
    Reported_padj_equals_nominal_P = padj,
    Leading_edge_genes = leadingEdge,
    Interpretation = P_value_interpretation
  )

data.table::fwrite(
  custom_metaplasia_summary,
  file.path(subdirs$tables, "Custom_metaplasia_GSEA_ALL_VALID_GENES_summary.csv")
)

make_custom_gsea_curve_audited <- function(
    rank_obj,
    fgsea_result,
    title_text,
    comparison_label,
    p_label = "Nominal p"
) {
  curve <- calc_enrichment_curve(rank_obj$ranks, metaplasia_genes)
  if (nrow(curve) == 0) {
    return(make_empty_plot(title_text, "No overlapping genes"))
  }

  result <- fgsea_result %>% dplyr::slice(1)
  nes <- if (nrow(result) == 0) NA_real_ else result$NES[1]
  p_value <- if (nrow(result) == 0) NA_real_ else result$pval[1]
  overlap_size <- if (nrow(result) == 0) NA_integer_ else result$size[1]

  annotation <- paste0(
    "NES = ", ifelse(is.na(nes), "NA", sprintf("%.2f", nes)),
    "\n", p_label, " = ",
    ifelse(
      is.na(p_value),
      "NA",
      ifelse(
        p_value < 1e-4,
        formatC(p_value, format = "e", digits = 1),
        sprintf("%.3f", p_value)
      )
    ),
    "\nGenes = ", ifelse(is.na(overlap_size), "NA", overlap_size)
  )

  ggplot2::ggplot(curve, ggplot2::aes(x = Position_percent, y = Running_ES)) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.25,
      colour = "grey45"
    ) +
    ggplot2::geom_line(
      colour = comparison_curve_cols[[comparison_label]],
      linewidth = 0.75
    ) +
    ggplot2::geom_rug(
      data = curve %>% dplyr::filter(Hit),
      ggplot2::aes(x = Position_percent),
      sides = "b",
      inherit.aes = FALSE,
      colour = comparison_curve_cols[[comparison_label]],
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
      size = 5.3 / 2.845276
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
      title = format_signature_plot_title(title_text),
      x = "Ranked gene list position (%)",
      y = "Running enrichment score"
    ) +
    nature_theme() +
    ggplot2::theme(legend.position = "none")
}

for (comparison_label in c("Infected vs Non-Infected", "BRAF V600E vs WT in SSL_HP")) {
  prefix <- if (comparison_label == "Infected vs Non-Infected") {
    "Infection"
  } else {
    "BRAF_SSL_HP_only"
  }

  hallmark_rank_obj <- if (comparison_label == "Infected vs Non-Infected") {
    infection_rank
  } else {
    braf_rank
  }

  custom_rank_all <- if (comparison_label == "Infected vs Non-Infected") {
    infection_rank_all
  } else {
    braf_rank_all
  }

  custom_result_all <- custom_metaplasia_fgsea_audited %>%
    dplyr::filter(Comparison == comparison_label)

  hallmark_results_one <- hallmark_fgsea_audited %>%
    dplyr::filter(Comparison == comparison_label)

  p_ifna <- make_gsea_curve_plot(
    ranks = hallmark_rank_obj$ranks,
    pathway_genes = ifna_info$genes,
    pathway_name = ifna_info$name,
    title_text = "Interferon alpha response",
    comparison_label = comparison_label,
    fgsea_df = hallmark_results_one
  )

  p_tnfa <- make_gsea_curve_plot(
    ranks = hallmark_rank_obj$ranks,
    pathway_genes = tnfa_info$genes,
    pathway_name = tnfa_info$name,
    title_text = "TNFA/NFKB signalling",
    comparison_label = comparison_label,
    fgsea_df = hallmark_results_one
  )

  p_inflammatory <- make_gsea_curve_plot(
    ranks = hallmark_rank_obj$ranks,
    pathway_genes = inflammatory_info$genes,
    pathway_name = inflammatory_info$name,
    title_text = "Inflammatory response",
    comparison_label = comparison_label,
    fgsea_df = hallmark_results_one
  )

  p_metaplasia_all <- make_custom_gsea_curve_audited(
    rank_obj = custom_rank_all,
    fgsea_result = custom_result_all,
    title_text = "Serrated Lesion Gastric Metaplasia Gene Signature",
    comparison_label = comparison_label
  )

  save_plot_set(
    p_metaplasia_all,
    subdirs$gsea,
    paste0(prefix, "_GSEA_curve_", metaplasia_file_stem),
    90,
    70
  )

  audited_curve_set <- patchwork::wrap_plots(
    list(
      "Interferon alpha response" = p_ifna,
      "TNFA/NFKB signalling" = p_tnfa,
      "Inflammatory response" = p_inflammatory,
      "Serrated Lesion Gastric Metaplasia Gene Signature" = p_metaplasia_all
    ),
    ncol = 2
  ) +
    patchwork::plot_annotation(
      tag_levels = "a",
      title = comparison_label,
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
          size = panel_tag_pt
        )
      )
    )

  save_plot_set(
    audited_curve_set,
    subdirs$combined,
    paste0(prefix, "_selected_GSEA_curves_combined"),
    183,
    145
  )
}

infection_n_pairs <- unique(infection_score_stats_audited_df$Complete_pairs)
braf_samples_unique <- braf_sample_audit %>%
  dplyr::count(BRAF_Status, name = "N")

writeLines(
  c(
    "POST-RUN RESULTS AUDIT",
    "======================",
    "",
    "Infection sample scores:",
    paste0("  Complete infected/non-infected pairs: ", paste(infection_n_pairs, collapse = ", "), "."),
    paste0("  Primary plotted test: ", infection_score_primary_test, "."),
    "  The two-sided paired t-test is used for continuous passage-matched scores when PAIRED_T is selected.",
    "  Exact and asymptotic paired Wilcoxon p values are exported as sensitivity analyses.",
    "  With only three non-zero pairs, the exact two-sided signed-rank p value cannot be below 0.25.",
    "  All infection samples are from one organoid/donor series, so inference is within-line and exploratory.",
    "",
    "Signature filtering:",
    paste0("  Hallmark singscores: genes require mean logCPM >= ", signature_mean_logcpm_cutoff, "."),
    paste0("  ", metaplasia_display_name, " singscore: all present signature genes are retained in the primary score."),
    paste0("  Hallmark GSEA: DE-table logCPM >= ", gsea_logcpm_cutoff, "."),
    paste0("  ", metaplasia_display_name, " GSEA: ALL_VALID_GENES only; no logCPM threshold is ever applied."),
    "  The custom-set padj equals its nominal p because only one custom pathway is tested; do not call this Hallmark-family FDR.",
    "",
    "BRAF sample-level analysis:",
    paste0(
      "  Unique SSL/HP samples: ",
      paste(paste0(braf_samples_unique$BRAF_Status, " = ", braf_samples_unique$N), collapse = "; "),
      "."
    ),
    "  Initial available passage is the primary sample-level analysis: P1 where present, otherwise P2; strict P1 is exported as a sensitivity analysis.",
    "  A patient-collapsed score sensitivity is exported.",
    "  The supplied edgeR BRAF DE/GSEA table cannot be validly regenerated from CPM alone; raw counts/design would be required for a count-model sensitivity analysis.",
    "",
    "Time course:",
    "  Pooled Spearman results are labelled descriptive because repeated samples are not independent.",
    "  Random-intercept mixed-effects models by trajectory are exported for inference.",
    "",
    "PCA:",
    "  The infection PCA is coloured by infection status and labelled with passage using non-overlapping labels.",
    "  The time-course plot has been overwritten with a bottom legend and colour-preserving passage shapes."
  ),
  file.path(subdirs$qc, "RESULTS_AUDIT_NOTES.txt")
)

cat(
  paste0(
    "\n\nPOST-RUN AUDIT (final outputs)\n",
    "  Time-course PCA legend corrected.\n",
    "  Paired t-test plus exact-Wilcoxon sensitivity tables written.\n",
    "  Time-course mixed-effects tables written.\n",
    "  Unfiltered 24-gene metaplasia GSEA and gene audits written.\n",
    "  Initial-passage BRAF PCA/scores plus strict-P1 sensitivity analyses written.\n",
    "  Individual titled infection and BRAF heatmaps written.\n"
  ),
  file = file.path(out_dir, "RUN_NOTES.txt"),
  append = TRUE
)

audit_package_versions <- data.frame(
  Package = c("nlme"),
  Version = c(as.character(utils::packageVersion("nlme"))),
  stringsAsFactors = FALSE
)

package_versions_final <- dplyr::bind_rows(package_versions, audit_package_versions) %>%
  dplyr::distinct(Package, .keep_all = TRUE)

data.table::fwrite(
  package_versions_final,
  file.path(subdirs$qc, "PACKAGE_VERSIONS.csv")
)

message("Post-run audit complete. Final audited output directory: ", out_dir)

message("\nRunning BRAF initial-passage versus strict-P1 sensitivity analysis...")

pca_braf_p1_meta <- braf_p1_obj$meta %>%
  dplyr::filter(
    (Type == "SSL_HP" & BRAF_Status %in% c("WT", "V600E")) |
      Type == "H_ADJ_SSL_HP"
  ) %>%
  dplyr::mutate(
    PCA_Group = dplyr::case_when(
      Type == "SSL_HP" & BRAF_Status == "WT" ~ "SSL_HP WT",
      Type == "SSL_HP" & BRAF_Status == "V600E" ~ "SSL_HP V600E",
      Type == "H_ADJ_SSL_HP" ~ "H_ADJ_SSL_HP",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(PCA_Group))

p1_groups_present <- pca_braf_p1_meta %>%
  dplyr::filter(Type == "SSL_HP") %>%
  dplyr::distinct(BRAF_Status) %>%
  dplyr::pull(BRAF_Status)

if (!all(c("WT", "V600E") %in% p1_groups_present)) {
  stop("Strict-P1 BRAF PCA requires both WT and V600E SSL/HP samples.")
}

pca_braf_p1_expr <- braf_p1_obj$log_expr[
  ,
  pca_braf_p1_meta$Sample_ID,
  drop = FALSE
]
pca_braf_p1 <- compute_pca(
  pca_braf_p1_expr,
  pca_braf_p1_meta,
  "BRAF_PCA_strict_P1"
)

p1_pca_levels <- c("SSL_HP WT", "SSL_HP V600E", "H_ADJ_SSL_HP")
p1_pca_levels <- p1_pca_levels[
  p1_pca_levels %in% unique(pca_braf_p1$scores$PCA_Group)
]
pca_braf_p1$scores$PCA_Group <- factor(
  pca_braf_p1$scores$PCA_Group,
  levels = p1_pca_levels
)

p_braf_p1_pca <- ggplot2::ggplot(
  pca_braf_p1$scores,
  ggplot2::aes(x = PC1, y = PC2, fill = PCA_Group)
) +
  ggplot2::geom_point(
    shape = 21,
    size = point_size,
    stroke = point_stroke,
    colour = "black",
    alpha = 0.95
  ) +
  ggplot2::scale_fill_manual(
    values = braf_pca_cols[p1_pca_levels],
    breaks = p1_pca_levels,
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "Strict P1 only",
    x = sprintf("PC1 (%.1f%%)", pca_braf_p1$variance_explained[1]),
    y = sprintf("PC2 (%.1f%%)", pca_braf_p1$variance_explained[2]),
    fill = NULL
  ) +
  nature_theme() +
  ggplot2::theme(
    legend.position = "right",
    plot.title = ggplot2::element_text(
      family = font_family,
      face = "bold",
      size = plot_title_pt,
      hjust = 0
    )
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      override.aes = list(shape = 21, size = point_size, colour = "black", alpha = 1)
    )
  )

save_plot_set(
  p_braf_p1_pca,
  subdirs$pca,
  "PCA_3B_BRAF_WT_vs_V600E_Strict_P1_with_H_ADJ",
  100,
  72
)

p_braf_initial_for_comparison <- p_braf_pca +
  ggplot2::labs(title = "Initial passage: P1, or P2 if P1 unavailable") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      family = font_family,
      face = "bold",
      size = plot_title_pt,
      hjust = 0
    )
  )

braf_pca_passage_comparison <- patchwork::wrap_plots(
  p_braf_initial_for_comparison,
  p_braf_p1_pca,
  ncol = 2,
  guides = "keep"
) +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = "BRAF PCA passage sensitivity",
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
        size = panel_tag_pt
      )
    )
  )

save_plot_set(
  braf_pca_passage_comparison,
  subdirs$combined,
  "BRAF_PCA_Initial_P1_or_P2_vs_Strict_P1_comparison",
  183,
  92
)

# Strict-P1 score sensitivity ----
braf_p1_scores <- score_signatures(braf_p1_obj, "BRAF_strict_P1")

braf_p1_score_input <- braf_p1_scores$scores %>%
  dplyr::filter(
    Type == "SSL_HP",
    BRAF_Status %in% c("WT", "V600E")
  )

braf_p1_score_plots <- list()
braf_p1_score_stats <- list()

for (sig in names(signature_list)) {
  p <- make_two_group_score_plot_audited(
    braf_p1_score_input,
    signature_name = sig,
    group_col = "BRAF_Status",
    group_levels = c("WT", "V600E"),
    colours = braf_cols,
    x_labels = c("WT", "V600E"),
    pair_col = NULL,
    prefer_paired = FALSE
  ) +
    ggplot2::labs(subtitle = "Strict P1 only") +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(
        family = font_family,
        size = 5.5,
        hjust = 0
      )
    )

  braf_p1_score_plots[[sig]] <- p
  save_plot_set(
    p,
    subdirs$scores,
    paste0("BRAF_Strict_P1_SSL_HP_only_", signature_file_stems[[sig]], "_singscore"),
    64,
    74
  )

  test <- run_two_group_test_audited(
    braf_p1_score_input %>% dplyr::filter(Signature == sig),
    "BRAF_Status",
    c("WT", "V600E"),
    NULL,
    FALSE
  )

  braf_p1_score_stats[[sig]] <- data.frame(
    Analysis_code = "Strict_P1",
    Analysis_set = "Strict P1 only",
    Signature = sig,
    Level_1 = "WT",
    Level_2 = "V600E",
    Test_primary = test$primary_method,
    N_WT = test$group_1_n,
    N_V600E = test$group_2_n,
    P_value = test$p_value,
    P_value_exact = test$p_value_exact,
    P_value_asymptotic_no_correction = test$p_value_asymptotic_no_correction,
    P_value_asymptotic_continuity_corrected = test$p_value_asymptotic_continuity_corrected,
    Mean_difference_V600E_minus_WT = test$mean_difference_level_2_minus_level_1,
    Median_difference_V600E_minus_WT = test$median_difference_level_2_minus_level_1,
    stringsAsFactors = FALSE
  )
}

braf_p1_score_stats_df <- dplyr::bind_rows(braf_p1_score_stats) %>%
  dplyr::mutate(
    P_value_BH_across_four_prespecified_scores = stats::p.adjust(
      P_value,
      method = "BH"
    )
  )

data.table::fwrite(
  braf_p1_score_stats_df,
  file.path(
    subdirs$tables,
    "BRAF_Strict_P1_SSL_HP_only_signature_score_statistics.csv"
  )
)

braf_p1_score_combined <- patchwork::wrap_plots(
  braf_p1_score_plots,
  ncol = 4
) +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = "BRAF WT versus V600E: strict P1 only",
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
        size = panel_tag_pt
      )
    )
  )

save_plot_set(
  braf_p1_score_combined,
  subdirs$combined,
  "BRAF_Strict_P1_SSL_HP_only_signature_scores_combined",
  183,
  80
)

initial_stats_for_comparison <- braf_score_stats_audited_df %>%
  dplyr::mutate(
    Analysis_code = "Initial_P1_or_P2",
    Analysis_set = "Initial passage: P1, or P2 if P1 unavailable",
    P_value_BH_across_four_prespecified_scores = stats::p.adjust(
      P_value,
      method = "BH"
    )
  ) %>%
  dplyr::select(
    Analysis_set,
    dplyr::everything()
  )

braf_passage_score_stats_long <- dplyr::bind_rows(
  initial_stats_for_comparison,
  braf_p1_score_stats_df
)

data.table::fwrite(
  braf_passage_score_stats_long,
  file.path(
    subdirs$tables,
    "BRAF_Initial_P1_or_P2_vs_Strict_P1_signature_statistics_long.csv"
  )
)

braf_passage_score_stats_wide <- braf_passage_score_stats_long %>%
  dplyr::select(
    Analysis_code,
    Signature,
    N_WT,
    N_V600E,
    P_value,
    P_value_exact,
    P_value_BH_across_four_prespecified_scores,
    Mean_difference_V600E_minus_WT,
    Median_difference_V600E_minus_WT
  ) %>%
  tidyr::pivot_wider(
    names_from = Analysis_code,
    values_from = c(
      N_WT,
      N_V600E,
      P_value,
      P_value_exact,
      P_value_BH_across_four_prespecified_scores,
      Mean_difference_V600E_minus_WT,
      Median_difference_V600E_minus_WT
    ),
    names_glue = "{.value}__{Analysis_code}"
  )

data.table::fwrite(
  braf_passage_score_stats_wide,
  file.path(
    subdirs$tables,
    "BRAF_Initial_P1_or_P2_vs_Strict_P1_signature_statistics_wide.csv"
  )
)

score_comparison_initial_plots <- lapply(names(signature_list), function(sig) {
  braf_score_plots_audited[[sig]] +
    ggplot2::labs(subtitle = "Initial passage: P1/P2") +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(
        family = font_family,
        size = 5.5,
        hjust = 0
      )
    )
})

score_comparison_p1_plots <- lapply(names(signature_list), function(sig) {
  braf_p1_score_plots[[sig]]
})

score_comparison_plots <- c(
  score_comparison_initial_plots,
  score_comparison_p1_plots
)

braf_score_passage_comparison <- patchwork::wrap_plots(
  score_comparison_plots,
  ncol = 4
) +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = "BRAF score sensitivity: initial passage versus strict P1",
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
        size = panel_tag_pt
      )
    )
  )

save_plot_set(
  braf_score_passage_comparison,
  subdirs$combined,
  "BRAF_Initial_P1_or_P2_vs_Strict_P1_signature_score_comparison",
  183,
  150
)

muc5ac_braf_initial_sensitivity <- extract_gene_expression(braf_obj, "MUC5AC")
muc5ac_braf_p1_sensitivity <- extract_gene_expression(braf_p1_obj, "MUC5AC")

if (!is.null(muc5ac_braf_initial_sensitivity) &&
    !is.null(muc5ac_braf_p1_sensitivity)) {
  p_muc_initial_sensitivity <- make_two_group_score_plot_audited(
    muc5ac_braf_initial_sensitivity %>% dplyr::filter(Type == "SSL_HP"),
    "MUC5AC",
    "BRAF_Status",
    c("WT", "V600E"),
    braf_cols,
    c("WT", "V600E"),
    NULL,
    FALSE,
    y_label = "MUC5AC log2 expression"
  ) +
    ggplot2::labs(subtitle = "Initial passage: P1/P2")

  p_muc_p1_sensitivity <- make_two_group_score_plot_audited(
    muc5ac_braf_p1_sensitivity %>% dplyr::filter(Type == "SSL_HP"),
    "MUC5AC",
    "BRAF_Status",
    c("WT", "V600E"),
    braf_cols,
    c("WT", "V600E"),
    NULL,
    FALSE,
    y_label = "MUC5AC log2 expression"
  ) +
    ggplot2::labs(subtitle = "Strict P1 only")

  save_plot_set(
    p_muc_p1_sensitivity,
    subdirs$scores,
    "BRAF_Strict_P1_SSL_HP_only_MUC5AC_expression",
    64,
    74
  )

  save_plot_set(
    p_muc_initial_sensitivity | p_muc_p1_sensitivity,
    subdirs$combined,
    "BRAF_Initial_P1_or_P2_vs_Strict_P1_MUC5AC_comparison",
    130,
    78
  )
}

braf_initial_heatmap_meta <- braf_obj$meta %>%
  dplyr::filter(
    Type == "SSL_HP",
    BRAF_Status %in% c("WT", "V600E")
  ) %>%
  dplyr::mutate(
    BRAF_Status = factor(BRAF_Status, levels = c("WT", "V600E"))
  )

braf_p1_heatmap_meta <- braf_p1_obj$meta %>%
  dplyr::filter(
    Type == "SSL_HP",
    BRAF_Status %in% c("WT", "V600E")
  ) %>%
  dplyr::mutate(
    BRAF_Status = factor(BRAF_Status, levels = c("WT", "V600E"))
  )

p_braf_initial_metaplasia_heatmap <- make_signature_expression_heatmap(
  braf_obj,
  braf_initial_heatmap_meta,
  metaplasia_genes,
  "BRAF_Status",
  paste0(metaplasia_display_name, ": initial passage P1/P2 (sample-level row z-scores)"),
  paste0("BRAF_Initial_P1_or_P2_SSL_HP_only_", metaplasia_file_stem, "_sample_level_row_zscore_heatmap")
)

p_braf_p1_metaplasia_heatmap <- make_signature_expression_heatmap(
  braf_p1_obj,
  braf_p1_heatmap_meta,
  metaplasia_genes,
  "BRAF_Status",
  paste0(metaplasia_display_name, ": strict P1 only (sample-level row z-scores)"),
  paste0("BRAF_Strict_P1_SSL_HP_only_", metaplasia_file_stem, "_sample_level_row_zscore_heatmap")
)

braf_metaplasia_heatmap_passage_comparison <-
  p_braf_initial_metaplasia_heatmap /
  p_braf_p1_metaplasia_heatmap +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = paste0("BRAF ", metaplasia_display_name, ": sample-level row z-score passage sensitivity"),
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
        size = panel_tag_pt
      )
    )
  ) &
  ggplot2::theme(legend.position = "bottom")

save_plot_set(
  braf_metaplasia_heatmap_passage_comparison,
  subdirs$combined,
  paste0("BRAF_Initial_P1_or_P2_vs_Strict_P1_", metaplasia_file_stem, "_sample_level_row_zscore_heatmaps"),
  183,
  250
)

braf_gsea_passage_scope <- data.frame(
  Analysis = c(
    "BRAF sample-level PCA/scores/CPM heatmaps",
    "BRAF edgeR DE and GSEA"
  ),
  Initial_P1_or_P2_status = c(
    "Generated directly from the sample-level time-course CPM matrix",
    "Existing supplied first-passage DE table used; upstream P1/P2 membership is not encoded in the summary table"
  ),
  Strict_P1_status = c(
    "Generated directly from the sample-level time-course CPM matrix",
    "No separate count-model edgeR result; exploratory CPM-based limma-trend GSEA is generated in Section 20"
  ),
  Reason = c(
    "Both sample subsets can be selected directly from the time-course CPM matrix.",
    paste(
      "The supplied BRAF input is an already fitted edgeR DE summary table.",
      "Separate strict-P1 DE/GSEA requires the original raw-count matrix and edgeR design/sample manifest; CPM must not be used to refit edgeR."
    )
  ),
  stringsAsFactors = FALSE
)

data.table::fwrite(
  braf_gsea_passage_scope,
  file.path(subdirs$qc, "BRAF_initial_vs_strict_P1_DE_GSEA_scope.csv")
)

cat(
  paste0(
    "\n\nBRAF PASSAGE SENSITIVITY\n",
    "  Primary sample-level analysis: one initial sample per trajectory; P1 where present, otherwise P2.\n",
    "  Sensitivity analysis: strict P1 only.\n",
    "  No trajectory contributes both P1 and P2 to either analysis.\n",
    "  Separate strict-P1 DE/GSEA was not generated because the available BRAF input is a fitted DE summary table rather than raw counts plus design.\n"
  ),
  file = file.path(out_dir, "RUN_NOTES.txt"),
  append = TRUE
)

message(
  "BRAF passage sensitivity complete. Compare:\n",
  "  06_Combined_figures/BRAF_PCA_Initial_P1_or_P2_vs_Strict_P1_comparison.*\n",
  "  07_Tables/BRAF_Initial_P1_or_P2_vs_Strict_P1_signature_statistics_wide.csv\n",
  "  06_Combined_figures/BRAF_Initial_P1_or_P2_vs_Strict_P1_signature_score_comparison.*"
)

message("\nRunning exploratory CPM-based BRAF GSEA for initial P1/P2 and strict P1...")

install_if_missing("limma", bioc = TRUE)

# Exploratory limma-trend sensitivity using CPM inputs ----
run_braf_limma_trend <- function(
    expression_object,
    analysis_code,
    analysis_label,
    unit_of_analysis = c("ORGANOID_LINE", "PATIENT_COLLAPSED")
) {
  unit_of_analysis <- match.arg(unit_of_analysis)

  meta <- expression_object$meta %>%
    dplyr::filter(
      Type == "SSL_HP",
      BRAF_Status %in% c("WT", "V600E")
    ) %>%
    dplyr::mutate(
      Patient_ID = derive_patient_id_audited(Line_ID)
    ) %>%
    dplyr::arrange(BRAF_Status, Patient_ID, Line_ID, Passage, Sample_ID)

  if (nrow(meta) < 4L) {
    stop("Too few SSL/HP samples for exploratory BRAF limma analysis: ", analysis_label)
  }

  missing_samples <- setdiff(meta$Sample_ID, colnames(expression_object$log_expr))
  if (length(missing_samples) > 0L) {
    stop(
      "These samples were selected for ", analysis_label,
      " but were absent from the expression object:\n",
      paste(missing_samples, collapse = "\n")
    )
  }

  expr <- expression_object$log_expr[
    ,
    meta$Sample_ID,
    drop = FALSE
  ]

  if (unit_of_analysis == "PATIENT_COLLAPSED") {
    patient_conflicts <- meta %>%
      dplyr::group_by(Patient_ID) %>%
      dplyr::summarise(
        N_BRAF_groups = dplyr::n_distinct(BRAF_Status),
        BRAF_values = paste(sort(unique(BRAF_Status)), collapse = ";"),
        .groups = "drop"
      ) %>%
      dplyr::filter(N_BRAF_groups > 1L)

    if (nrow(patient_conflicts) > 0L) {
      data.table::fwrite(
        patient_conflicts,
        file.path(
          subdirs$qc,
          paste0(safe_filename(analysis_code), "_patient_BRAF_conflicts.csv")
        )
      )
      stop("Conflicting BRAF groups were found within a patient for ", analysis_label, ".")
    }

    patient_order <- unique(meta$Patient_ID)

    collapsed_expr <- vapply(
      patient_order,
      function(pid) {
        sample_ids <- meta$Sample_ID[meta$Patient_ID == pid]
        rowMeans(expr[, sample_ids, drop = FALSE], na.rm = TRUE)
      },
      FUN.VALUE = numeric(nrow(expr))
    )

    if (is.null(dim(collapsed_expr))) {
      collapsed_expr <- matrix(
        collapsed_expr,
        nrow = nrow(expr),
        dimnames = list(rownames(expr), patient_order)
      )
    } else {
      rownames(collapsed_expr) <- rownames(expr)
      colnames(collapsed_expr) <- patient_order
    }

    model_meta <- meta %>%
      dplyr::group_by(Patient_ID) %>%
      dplyr::summarise(
        BRAF_Status = consensus_character(BRAF_Status),
        N_organoid_lines = dplyr::n_distinct(Line_ID),
        N_samples_collapsed = dplyr::n(),
        Sample_IDs = paste(Sample_ID, collapse = ";"),
        Passages = paste(sort(unique(Passage)), collapse = ";"),
        .groups = "drop"
      ) %>%
      dplyr::mutate(Model_sample_ID = Patient_ID)

    model_meta <- model_meta[
      match(patient_order, model_meta$Patient_ID),
      ,
      drop = FALSE
    ]

    expr_model <- collapsed_expr[, model_meta$Model_sample_ID, drop = FALSE]
  } else {
    model_meta <- meta %>%
      dplyr::transmute(
        Model_sample_ID = Sample_ID,
        Patient_ID,
        Line_ID,
        Passage,
        BRAF_Status,
        N_organoid_lines = 1L,
        N_samples_collapsed = 1L,
        Sample_IDs = Sample_ID,
        Passages = as.character(Passage)
      )

    expr_model <- expr[, model_meta$Model_sample_ID, drop = FALSE]
  }

  group_counts <- table(model_meta$BRAF_Status)
  n_wt <- unname(group_counts["WT"])
  n_v600e <- unname(group_counts["V600E"])
  if (length(n_wt) == 0L || is.na(n_wt)) n_wt <- 0L
  if (length(n_v600e) == 0L || is.na(n_v600e)) n_v600e <- 0L

  if (n_wt < 2L || n_v600e < 2L) {
    stop(
      analysis_label, " requires at least two WT and two V600E model units. ",
      "Observed WT = ", n_wt, "; V600E = ", n_v600e, "."
    )
  }

  group <- factor(model_meta$BRAF_Status, levels = c("WT", "V600E"))
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  if (qr(design)$rank < ncol(design)) {
    stop("The exploratory limma design is rank deficient for ", analysis_label, ".")
  }

  fit <- limma::lmFit(expr_model, design)
  contrast_matrix <- limma::makeContrasts(
    V600E_vs_WT = V600E - WT,
    levels = design
  )
  fit <- limma::contrasts.fit(fit, contrast_matrix)
  fit <- limma::eBayes(fit, trend = TRUE, robust = FALSE)

  tt <- limma::topTable(
    fit,
    coef = "V600E_vs_WT",
    number = Inf,
    sort.by = "none",
    adjust.method = "BH"
  )

  de <- data.frame(
    Gene = rownames(tt),
    logFC = as.numeric(tt$logFC),
    AveExpr = as.numeric(tt$AveExpr),
    t = as.numeric(tt$t),
    PValue = as.numeric(tt$P.Value),
    FDR = as.numeric(tt$adj.P.Val),
    B = as.numeric(tt$B),
    Analysis_code = analysis_code,
    Analysis_set = analysis_label,
    Unit_of_analysis = unit_of_analysis,
    Contrast = "V600E - WT",
    N_WT = n_wt,
    N_V600E = n_v600e,
    Residual_df = unname(fit$df.residual[1]),
    Model = "limma-trend on log2(CPM + 1)",
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(
      !is.na(Gene),
      Gene != "",
      is.finite(t),
      is.finite(logFC),
      is.finite(PValue)
    ) %>%
    dplyr::group_by(Gene) %>%
    dplyr::arrange(
      dplyr::desc(abs(t)),
      PValue,
      dplyr::desc(abs(logFC)),
      .by_group = TRUE
    ) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  sample_manifest <- model_meta %>%
    dplyr::mutate(
      Analysis_code = analysis_code,
      Analysis_set = analysis_label,
      Unit_of_analysis = unit_of_analysis,
      .before = 1
    )

  list(
    de = de,
    meta = model_meta,
    sample_manifest = sample_manifest,
    expression = expr_model,
    analysis_code = analysis_code,
    analysis_label = analysis_label,
    unit_of_analysis = unit_of_analysis,
    n_wt = n_wt,
    n_v600e = n_v600e
  )
}

make_limma_t_rank <- function(de, apply_mean_expression_filter = TRUE) {
  dat <- de

  if (apply_mean_expression_filter) {
    dat <- dat %>%
      dplyr::filter(
        is.finite(AveExpr),
        AveExpr >= gsea_logcpm_cutoff
      )
  }

  dat <- dat %>%
    dplyr::filter(
      !is.na(Gene),
      Gene != "",
      is.finite(t)
    ) %>%
    dplyr::arrange(dplyr::desc(t), Gene)

  ranks <- dat$t
  names(ranks) <- dat$Gene

  if (anyDuplicated(ranks)) {
    ranks <- ranks - seq_along(ranks) * 1e-12
    dat$t_rank_used <- ranks
  } else {
    dat$t_rank_used <- ranks
  }

  ranks <- sort(ranks, decreasing = TRUE)

  list(
    ranks = ranks,
    table = dat,
    method = "limma moderated t statistic for V600E - WT",
    filter_mode = ifelse(
      apply_mean_expression_filter,
      paste0("MEAN_LOG2_CPM_GE_", gsea_logcpm_cutoff),
      "ALL_VALID_GENES"
    )
  )
}

run_cpm_limma_fgsea <- function(limma_result) {
  hall_rank <- make_limma_t_rank(
    limma_result$de,
    apply_mean_expression_filter = TRUE
  )

  metaplasia_rank <- make_limma_t_rank(
    limma_result$de,
    apply_mean_expression_filter = FALSE
  )

  comparison_label <- paste0(
    "BRAF CPM limma: ",
    limma_result$analysis_label,
    ifelse(
      limma_result$unit_of_analysis == "PATIENT_COLLAPSED",
      " (patient-collapsed)",
      ""
    )
  )

  hallmark_fgsea <- run_fgsea_one(
    hall_rank$ranks,
    hallmark_pathways,
    comparison_label,
    "HALLMARK"
  ) %>%
    dplyr::mutate(
      Analysis_code = limma_result$analysis_code,
      Analysis_set = limma_result$analysis_label,
      Unit_of_analysis = limma_result$unit_of_analysis,
      Filter_mode = hall_rank$filter_mode,
      Ranking_method = hall_rank$method,
      N_WT = limma_result$n_wt,
      N_V600E = limma_result$n_v600e,
      Result_scope = "Exploratory CPM-based passage-sensitivity GSEA"
    )

  metaplasia_fgsea <- run_fgsea_one(
    metaplasia_rank$ranks,
    custom_metaplasia_pathway,
    comparison_label,
    "CUSTOM_METAPLASIA_NO_LOGCPM_FILTER",
    min_size = 3L
  ) %>%
    dplyr::mutate(
      Analysis_code = limma_result$analysis_code,
      Analysis_set = limma_result$analysis_label,
      Unit_of_analysis = limma_result$unit_of_analysis,
      Filter_mode = "ALL_VALID_GENES",
      Ranking_method = metaplasia_rank$method,
      N_WT = limma_result$n_wt,
      N_V600E = limma_result$n_v600e,
      N_ranked_genes = length(metaplasia_rank$ranks),
      Signature_genes_requested = length(unique(metaplasia_genes)),
      Signature_genes_in_rank_universe = sum(
        unique(metaplasia_genes) %in% names(metaplasia_rank$ranks)
      ),
      P_value_interpretation = paste0(
        "Nominal custom-set p value. One custom pathway was tested; ",
        "padj equals pval and is not a Hallmark-family FDR."
      ),
      Result_scope = "Exploratory CPM-based passage-sensitivity GSEA"
    )

  list(
    hallmark = hallmark_fgsea,
    metaplasia = metaplasia_fgsea,
    hall_rank = hall_rank,
    metaplasia_rank = metaplasia_rank,
    comparison_label = comparison_label,
    limma_result = limma_result
  )
}

braf_cpm_limma_specs <- list(
  Initial_P1_or_P2_line = list(
    object = braf_obj,
    analysis_code = "Initial_P1_or_P2_line",
    analysis_label = "Initial P1/P2",
    unit = "ORGANOID_LINE"
  ),
  Strict_P1_line = list(
    object = braf_p1_obj,
    analysis_code = "Strict_P1_line",
    analysis_label = "Strict P1",
    unit = "ORGANOID_LINE"
  ),
  Initial_P1_or_P2_patient = list(
    object = braf_obj,
    analysis_code = "Initial_P1_or_P2_patient_collapsed",
    analysis_label = "Initial P1/P2",
    unit = "PATIENT_COLLAPSED"
  ),
  Strict_P1_patient = list(
    object = braf_p1_obj,
    analysis_code = "Strict_P1_patient_collapsed",
    analysis_label = "Strict P1",
    unit = "PATIENT_COLLAPSED"
  )
)

braf_cpm_limma_results <- list()
braf_cpm_limma_fgsea_results <- list()

for (spec_name in names(braf_cpm_limma_specs)) {
  spec <- braf_cpm_limma_specs[[spec_name]]

  limma_one <- run_braf_limma_trend(
    expression_object = spec$object,
    analysis_code = spec$analysis_code,
    analysis_label = spec$analysis_label,
    unit_of_analysis = spec$unit
  )

  fgsea_one <- run_cpm_limma_fgsea(limma_one)

  braf_cpm_limma_results[[spec_name]] <- limma_one
  braf_cpm_limma_fgsea_results[[spec_name]] <- fgsea_one

  prefix <- paste0("BRAF_CPM_limma_", safe_filename(spec$analysis_code))

  data.table::fwrite(
    limma_one$sample_manifest,
    file.path(subdirs$qc, paste0(prefix, "_sample_manifest.csv"))
  )
  data.table::fwrite(
    limma_one$de,
    file.path(subdirs$tables, paste0(prefix, "_DE_results.csv"))
  )
  data.table::fwrite(
    fgsea_one$hall_rank$table,
    file.path(subdirs$tables, paste0(prefix, "_Hallmark_ranked_genes.csv"))
  )
  data.table::fwrite(
    fgsea_one$metaplasia_rank$table,
    file.path(subdirs$tables, paste0(prefix, "_Metaplasia_ALL_VALID_GENES_ranked_genes.csv"))
  )
  data.table::fwrite(
    fgsea_one$hallmark,
    file.path(subdirs$tables, paste0(prefix, "_Hallmark_FGSEA.csv"))
  )
  data.table::fwrite(
    fgsea_one$metaplasia,
    file.path(subdirs$tables, paste0(prefix, "_Custom_metaplasia_FGSEA_NO_FILTER.csv"))
  )
}

braf_cpm_limma_hallmark_all <- dplyr::bind_rows(
  lapply(braf_cpm_limma_fgsea_results, `[[`, "hallmark")
)

braf_cpm_limma_metaplasia_all <- dplyr::bind_rows(
  lapply(braf_cpm_limma_fgsea_results, `[[`, "metaplasia")
)

braf_cpm_limma_de_all <- dplyr::bind_rows(
  lapply(braf_cpm_limma_results, `[[`, "de")
)

braf_cpm_limma_manifests_all <- dplyr::bind_rows(
  lapply(braf_cpm_limma_results, `[[`, "sample_manifest")
)

data.table::fwrite(
  braf_cpm_limma_hallmark_all,
  file.path(
    subdirs$tables,
    "BRAF_CPM_limma_Initial_vs_Strict_P1_Hallmark_FGSEA_ALL_SENSITIVITIES.csv"
  )
)

data.table::fwrite(
  braf_cpm_limma_metaplasia_all,
  file.path(
    subdirs$tables,
    "BRAF_CPM_limma_Initial_vs_Strict_P1_custom_metaplasia_FGSEA_NO_FILTER_ALL_SENSITIVITIES.csv"
  )
)

data.table::fwrite(
  braf_cpm_limma_de_all,
  file.path(
    subdirs$tables,
    "BRAF_CPM_limma_Initial_vs_Strict_P1_DE_results_ALL_SENSITIVITIES.csv"
  )
)

data.table::fwrite(
  braf_cpm_limma_manifests_all,
  file.path(
    subdirs$qc,
    "BRAF_CPM_limma_Initial_vs_Strict_P1_sample_manifest_ALL_SENSITIVITIES.csv"
  )
)

initial_line_code <- "Initial_P1_or_P2_line"
strict_line_code <- "Strict_P1_line"

initial_line_fgsea <- braf_cpm_limma_fgsea_results[["Initial_P1_or_P2_line"]]
strict_line_fgsea <- braf_cpm_limma_fgsea_results[["Strict_P1_line"]]

p_braf_cpm_gsea_initial <- make_hallmark_barplot(
  braf_cpm_limma_hallmark_all,
  initial_line_fgsea$comparison_label
) +
  ggplot2::labs(title = "Initial P1/P2: exploratory limma GSEA")

p_braf_cpm_gsea_strict <- make_hallmark_barplot(
  braf_cpm_limma_hallmark_all,
  strict_line_fgsea$comparison_label
) +
  ggplot2::labs(title = "Strict P1: exploratory limma GSEA")

save_plot_set(
  p_braf_cpm_gsea_initial,
  subdirs$gsea,
  "BRAF_CPM_limma_Initial_P1_or_P2_Hallmark_GSEA_top_barplot",
  100,
  95
)

save_plot_set(
  p_braf_cpm_gsea_strict,
  subdirs$gsea,
  "BRAF_CPM_limma_Strict_P1_Hallmark_GSEA_top_barplot",
  100,
  95
)

p_braf_cpm_gsea_bars_combined <-
  (p_braf_cpm_gsea_initial | p_braf_cpm_gsea_strict) +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = "BRAF Hallmark GSEA passage sensitivity",
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
        size = panel_tag_pt
      )
    )
  )

save_plot_set(
  p_braf_cpm_gsea_bars_combined,
  subdirs$combined,
  "BRAF_CPM_limma_Initial_P1_or_P2_vs_Strict_P1_Hallmark_GSEA_barplots",
  183,
  102
)

make_cpm_limma_hallmark_comparison_heatmap <- function(fgsea_df) {
  dat <- fgsea_df %>%
    dplyr::filter(
      Unit_of_analysis == "ORGANOID_LINE",
      Analysis_code %in% c(initial_line_code, strict_line_code),
      Gene_set_family == "HALLMARK",
      is.finite(NES)
    )

  selected_pathways <- dat %>%
    dplyr::group_by(Analysis_code) %>%
    dplyr::arrange(padj, dplyr::desc(abs(NES)), .by_group = TRUE) %>%
    dplyr::slice_head(n = hallmark_heatmap_top_n) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(pathway) %>%
    dplyr::pull(pathway)

  dat <- dat %>%
    dplyr::filter(pathway %in% selected_pathways) %>%
    dplyr::mutate(
      Analysis_display = dplyr::case_when(
        Analysis_code == initial_line_code ~ "Initial P1/P2",
        Analysis_code == strict_line_code ~ "Strict P1",
        TRUE ~ Analysis_set
      ),
      Analysis_display = factor(
        Analysis_display,
        levels = c("Initial P1/P2", "Strict P1")
      ),
      Pathway_label = clean_pathway_label(pathway),
      Significance = dplyr::case_when(
        padj < 0.001 ~ "***",
        padj < 0.01 ~ "**",
        padj < 0.05 ~ "*",
        pval < 0.05 ~ ".",
        TRUE ~ ""
      )
    )

  pathway_order <- dat %>%
    dplyr::group_by(Pathway_label) %>%
    dplyr::summarise(
      Max_abs_NES = max(abs(NES), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Max_abs_NES) %>%
    dplyr::pull(Pathway_label)

  dat$Pathway_label <- factor(dat$Pathway_label, levels = pathway_order)

  ggplot2::ggplot(
    dat,
    ggplot2::aes(x = Analysis_display, y = Pathway_label, fill = NES)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::geom_text(
      ggplot2::aes(label = Significance),
      family = font_family,
      size = 5.4 / 2.845276
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "NES"
    ) +
    ggplot2::labs(
      title = "Exploratory BRAF Hallmark GSEA: passage sensitivity",
      subtitle = "limma-trend on log2(CPM + 1); V600E - WT",
      x = NULL,
      y = NULL
    ) +
    nature_theme() +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 6.2),
      axis.text.y = ggplot2::element_text(size = 5.6),
      plot.subtitle = ggplot2::element_text(
        family = font_family,
        size = 5.5,
        hjust = 0
      ),
      legend.position = "right"
    )
}

p_braf_cpm_gsea_heatmap <- make_cpm_limma_hallmark_comparison_heatmap(
  braf_cpm_limma_hallmark_all
)

save_plot_set(
  p_braf_cpm_gsea_heatmap,
  subdirs$heatmaps,
  "BRAF_CPM_limma_Initial_P1_or_P2_vs_Strict_P1_Hallmark_NES_heatmap",
  120,
  130
)

comparison_curve_cols <- c(
  comparison_curve_cols,
  stats::setNames("#4C78A8", initial_line_fgsea$comparison_label),
  stats::setNames("#E45756", strict_line_fgsea$comparison_label)
)

braf_cpm_line_curve_sets <- list(
  Initial_P1_or_P2 = initial_line_fgsea,
  Strict_P1 = strict_line_fgsea
)

for (analysis_name in names(braf_cpm_line_curve_sets)) {
  result_one <- braf_cpm_line_curve_sets[[analysis_name]]
  fgsea_for_curves <- dplyr::bind_rows(
    result_one$hallmark,
    result_one$metaplasia
  )

  curve_plots <- list()

  for (display_name in names(selected_gsea_pathways)) {
    pathway_info <- selected_gsea_pathways[[display_name]]
    ranks_for_plot <- if (display_name == "Serrated Lesion Gastric Metaplasia Gene Signature") {
      result_one$metaplasia_rank$ranks
    } else {
      result_one$hall_rank$ranks
    }

    p_curve <- make_gsea_curve_plot(
      ranks = ranks_for_plot,
      pathway_genes = pathway_info$genes,
      pathway_name = pathway_info$pathway,
      title_text = display_name,
      comparison_label = result_one$comparison_label,
      fgsea_df = fgsea_for_curves
    )

    curve_plots[[display_name]] <- p_curve

    save_plot_set(
      p_curve,
      subdirs$gsea,
      paste0(
        "BRAF_CPM_limma_",
        analysis_name,
        "_GSEA_curve_",
        safe_filename(display_name)
      ),
      90,
      70
    )
  }

  p_curves_combined <- patchwork::wrap_plots(curve_plots, ncol = 2) +
    patchwork::plot_annotation(
      tag_levels = "a",
      title = paste0(
        "BRAF ",
        ifelse(analysis_name == "Strict_P1", "strict P1", "initial P1/P2"),
        ": exploratory limma GSEA"
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
          size = panel_tag_pt
        )
      )
    )

  save_plot_set(
    p_curves_combined,
    subdirs$combined,
    paste0("BRAF_CPM_limma_", analysis_name, "_selected_GSEA_curves_combined"),
    183,
    145
  )
}

braf_cpm_limma_gsea_summary <- dplyr::bind_rows(
  braf_cpm_limma_hallmark_all %>%
    dplyr::filter(
      Unit_of_analysis == "ORGANOID_LINE",
      Analysis_code %in% c(initial_line_code, strict_line_code)
    ) %>%
    dplyr::transmute(
      Analysis_code,
      Analysis_set,
      Unit_of_analysis,
      Gene_set_family,
      pathway,
      NES,
      pval,
      padj,
      size,
      leadingEdge,
      Filter_mode,
      N_WT,
      N_V600E
    ),
  braf_cpm_limma_metaplasia_all %>%
    dplyr::filter(
      Unit_of_analysis == "ORGANOID_LINE",
      Analysis_code %in% c(initial_line_code, strict_line_code)
    ) %>%
    dplyr::transmute(
      Analysis_code,
      Analysis_set,
      Unit_of_analysis,
      Gene_set_family,
      pathway,
      NES,
      pval,
      padj,
      size,
      leadingEdge,
      Filter_mode,
      N_WT,
      N_V600E
    )
)

data.table::fwrite(
  braf_cpm_limma_gsea_summary,
  file.path(
    subdirs$tables,
    "BRAF_CPM_limma_Initial_P1_or_P2_vs_Strict_P1_GSEA_direct_comparison.csv"
  )
)

braf_cpm_limma_scope <- data.frame(
  Analysis = c(
    "Primary supplied BRAF edgeR GSEA",
    "Exploratory initial P1/P2 CPM GSEA",
    "Exploratory strict-P1 CPM GSEA",
    "Patient-collapsed CPM sensitivity"
  ),
  Input = c(
    "Supplied count-model DE summary table",
    "Initial P1/P2 log2(CPM + 1) matrix",
    "Strict P1 log2(CPM + 1) matrix",
    "Patient-mean log2(CPM + 1) matrix"
  ),
  Gene_level_model = c(
    "Existing edgeR likelihood-ratio model",
    "limma-trend moderated linear model",
    "limma-trend moderated linear model",
    "limma-trend moderated linear model"
  ),
  Ranking = c(
    "sign(logFC) * sqrt(LR)",
    "moderated t for V600E - WT",
    "moderated t for V600E - WT",
    "moderated t for V600E - WT"
  ),
  Interpretation = c(
    "Primary pathway result, subject to the original model/design",
    "Method-matched passage-sensitivity analysis",
    "Method-matched passage-sensitivity analysis",
    "Sensitivity to repeated organoid lines from the same patient"
  ),
  stringsAsFactors = FALSE
)

data.table::fwrite(
  braf_cpm_limma_scope,
  file.path(subdirs$qc, "BRAF_CPM_limma_GSEA_scope_and_methods.csv")
)

writeLines(
  c(
    "STRICT-P1 BRAF GSEA SENSITIVITY",
    "=================================",
    "",
    "The supplied edgeR BRAF DE table remains the primary GSEA input.",
    "Strict-P1 edgeR cannot be reconstructed from CPM because raw counts, library sizes and the original design are unavailable.",
    "",
    "For passage sensitivity, both initial P1/P2 and strict P1 were therefore reanalysed with the same limma-trend model on log2(CPM + 1).",
    "Genes were ranked by the moderated t statistic for V600E - WT.",
    paste0("Hallmark GSEA used genes with mean log2(CPM + 1) >= ", gsea_logcpm_cutoff, "."),
    paste0(metaplasia_display_name, " GSEA used ALL_VALID_GENES and no expression threshold."),
    "H_ADJ samples were excluded from all BRAF limma/GSEA models.",
    "",
    "The line-level strict-P1 analysis contains two WT and four V600E organoid lines.",
    "Because GM_7699 contributes two V600E lines, patient-collapsed results are exported as a dependence sensitivity.",
    "",
    "These fgsea p values describe gene-set enrichment in the fitted rank and do not create additional biological replicates.",
    "With the small number of independent lines/patients, use concordance of NES direction and leading-edge genes across sample selections as the main sensitivity criterion, not whether one analysis crosses p = 0.05."
  ),
  file.path(subdirs$qc, "BRAF_CPM_limma_GSEA_INTERPRETATION.txt")
)

cat(
  paste0(
    "\n\nSTRICT-P1 EXPLORATORY GSEA ADDED\n",
    "  Both initial P1/P2 and strict P1 were analysed with limma-trend on log2(CPM + 1).\n",
    "  Hallmark ranking uses mean log2(CPM + 1) >= ", gsea_logcpm_cutoff, ".\n",
    paste0("  ", metaplasia_display_name, " ranking uses ALL_VALID_GENES with no expression threshold.\n"),
    "  Patient-collapsed sensitivity results were also exported.\n",
    "  These outputs are exploratory and do not replace the supplied edgeR GSEA.\n"
  ),
  file = file.path(out_dir, "RUN_NOTES.txt"),
  append = TRUE
)

message(
  "Exploratory strict-P1 BRAF GSEA complete. Main files:\n",
  "  07_Tables/BRAF_CPM_limma_Initial_P1_or_P2_vs_Strict_P1_GSEA_direct_comparison.csv\n",
  "  05_Heatmaps/BRAF_CPM_limma_Initial_P1_or_P2_vs_Strict_P1_Hallmark_NES_heatmap.*\n",
  "  06_Combined_figures/BRAF_CPM_limma_Strict_P1_selected_GSEA_curves_combined.*\n",
  "  01_QC/BRAF_CPM_limma_GSEA_INTERPRETATION.txt"
)

# Volcano plots and saved plotting tables ----
volcano_dir <- file.path(out_dir, "03_Volcano_plots")
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

volcano_logcpm_cutoff <- 2
volcano_logfc_cutoff <- 1
volcano_fdr_cutoff <- 0.05
volcano_n_labels <- 10

volcano_cols <- c(
  "Upregulated" = "#D62728",
  "Downregulated" = "#1F77B4",
  "Not significant" = "#8E8E8E"
)

prepare_volcano_table <- function(
    de,
    comparison_label,
    direction_positive = "Upregulated",
    direction_negative = "Downregulated"
) {
  input_n <- nrow(de)

  dat <- de %>%
    dplyr::mutate(
      Display_pass = is.finite(logCPM) & logCPM > volcano_logcpm_cutoff,
      Significant = Display_pass & is.finite(FDR) & FDR < volcano_fdr_cutoff &
        is.finite(logFC) & abs(logFC) >= volcano_logfc_cutoff,
      Regulation = dplyr::case_when(
        Significant & logFC >= volcano_logfc_cutoff ~ direction_positive,
        Significant & logFC <= -volcano_logfc_cutoff ~ direction_negative,
        TRUE ~ "Not significant"
      ),
      neg_log10_FDR = -log10(pmax(FDR, .Machine$double.xmin)),
      neg_log10_PValue = -log10(pmax(PValue, .Machine$double.xmin)),
      Comparison = comparison_label
    ) %>%
    dplyr::filter(Display_pass, is.finite(logFC), is.finite(neg_log10_FDR)) %>%
    dplyr::arrange(PValue, dplyr::desc(abs(logFC)))

  labelled <- dat %>%
    dplyr::filter(Regulation != "Not significant", !is.na(Gene), Gene != "") %>%
    dplyr::arrange(PValue, dplyr::desc(abs(logFC))) %>%
    dplyr::slice_head(n = volcano_n_labels)

  if (nrow(labelled) < volcano_n_labels) {
    labelled <- dat %>%
      dplyr::filter(!is.na(Gene), Gene != "") %>%
      dplyr::arrange(PValue, dplyr::desc(abs(logFC))) %>%
      dplyr::slice_head(n = volcano_n_labels)
  }

  labelled <- labelled %>%
    dplyr::group_by(Gene) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  dat$Label <- ifelse(dat$Gene %in% labelled$Gene, dat$Gene, NA_character_)

  summary_table <- data.frame(
    Comparison = comparison_label,
    Genes_in_input = input_n,
    Genes_after_logCPM_display_filter = nrow(dat),
    Significant_total = sum(dat$Significant, na.rm = TRUE),
    Upregulated_significant = sum(dat$Regulation == direction_positive, na.rm = TRUE),
    Downregulated_significant = sum(dat$Regulation == direction_negative, na.rm = TRUE),
    logCPM_display_cutoff = volcano_logcpm_cutoff,
    abs_logFC_cutoff = volcano_logfc_cutoff,
    FDR_cutoff = volcano_fdr_cutoff,
    Labelled_genes = paste(labelled$Gene, collapse = "; "),
    stringsAsFactors = FALSE
  )

  list(data = dat, labelled = labelled, summary = summary_table)
}

build_volcano_plot <- function(
    volcano_object,
    title_text,
    subtitle_text = NULL,
    x_limits = NULL,
    y_limits = NULL
) {
  dat <- volcano_object$data
  labelled <- volcano_object$labelled

  if (nrow(dat) == 0) {
    return(make_empty_plot(title_text, "No genes remained after volcano display filtering."))
  }

  if (is.null(x_limits)) {
    max_abs_fc <- max(abs(dat$logFC), na.rm = TRUE)
    x_limits <- c(-1, 1) * max(2.2, ceiling(max_abs_fc * 1.05))
  }

  if (is.null(y_limits)) {
    max_y <- max(dat$neg_log10_FDR, na.rm = TRUE)
    y_limits <- c(0, max(5, ceiling(max_y * 1.05)))
  }

  dat$Regulation <- factor(
    dat$Regulation,
    levels = c("Upregulated", "Downregulated", "Not significant")
  )

  p <- ggplot(dat, aes(x = logFC, y = neg_log10_FDR)) +
    geom_hline(
      yintercept = -log10(volcano_fdr_cutoff),
      linetype = "dashed", colour = "black", linewidth = 0.35
    ) +
    geom_vline(
      xintercept = c(-volcano_logfc_cutoff, volcano_logfc_cutoff),
      linetype = "dashed", colour = "black", linewidth = 0.35
    ) +
    geom_point(
      aes(colour = Regulation),
      size = 0.9, stroke = 0, alpha = 0.85
    ) +
    scale_colour_manual(values = volcano_cols, drop = FALSE) +
    coord_cartesian(xlim = x_limits, ylim = y_limits, clip = "off") +
    labs(
      title = format_signature_plot_title(title_text),
      subtitle = subtitle_text,
      x = expression(log[2]~fold~change),
      y = expression(-log[10]~adjusted~italic(P)~value),
      colour = NULL
    ) +
    nature_theme() +
    theme(
      plot.subtitle = element_blank(),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.direction = "vertical",
      legend.key.height = unit(3.5, "mm"),
      legend.key.width = unit(4.5, "mm"),
      legend.background = element_rect(fill = alpha("white", 0.9), colour = NA)
    )

  if (nrow(labelled) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = labelled,
      aes(label = Gene),
      size = 2.1,
      family = font_family,
      colour = "black",
      box.padding = 0.28,
      point.padding = 0.12,
      segment.colour = "grey55",
      segment.linewidth = 0.28,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 1
    )
  }

  p
}

infection_volcano <- prepare_volcano_table(
  infection_de,
  comparison_label = "Infected vs Non-Infected",
  direction_positive = "Upregulated",
  direction_negative = "Downregulated"
)

braf_volcano <- prepare_volcano_table(
  braf_de,
  comparison_label = "BRAF V600E vs WT in SSL/HP",
  direction_positive = "Upregulated",
  direction_negative = "Downregulated"
)

data.table::fwrite(
  infection_volcano$data,
  file.path(subdirs$tables, "Infection_volcano_plot_data_logCPM_gt_2.csv")
)
data.table::fwrite(
  infection_volcano$labelled,
  file.path(subdirs$tables, "Infection_volcano_labelled_genes_top10_by_PValue.csv")
)
data.table::fwrite(
  infection_volcano$summary,
  file.path(subdirs$tables, "Infection_volcano_summary.csv")
)

data.table::fwrite(
  braf_volcano$data,
  file.path(subdirs$tables, "BRAF_SSL_HP_only_volcano_plot_data_logCPM_gt_2.csv")
)
data.table::fwrite(
  braf_volcano$labelled,
  file.path(subdirs$tables, "BRAF_SSL_HP_only_volcano_labelled_genes_top10_by_PValue.csv")
)
data.table::fwrite(
  braf_volcano$summary,
  file.path(subdirs$tables, "BRAF_SSL_HP_only_volcano_summary.csv")
)

volcano_xmax <- max(
  abs(infection_volcano$data$logFC),
  abs(braf_volcano$data$logFC),
  na.rm = TRUE
)
volcano_ymax <- max(
  infection_volcano$data$neg_log10_FDR,
  braf_volcano$data$neg_log10_FDR,
  na.rm = TRUE
)
common_x_limits <- c(-1, 1) * max(2.2, ceiling(volcano_xmax * 1.05))
common_y_limits <- c(0, max(5, ceiling(volcano_ymax * 1.05)))

infection_volcano_plot <- build_volcano_plot(
  infection_volcano,
  title_text = "Infected vs non-infected organoids",
  subtitle_text = NULL,
  x_limits = common_x_limits,
  y_limits = common_y_limits
)

braf_volcano_plot <- build_volcano_plot(
  braf_volcano,
  title_text = "BRAF V600E vs WT in SSL_HP organoids",
  subtitle_text = NULL,
  x_limits = common_x_limits,
  y_limits = common_y_limits
)

save_plot_set(
  infection_volcano_plot,
  volcano_dir,
  "Volcano_Infected_vs_Non_Infected",
  width_mm = 90,
  height_mm = 85
)

save_plot_set(
  braf_volcano_plot,
  volcano_dir,
  "Volcano_BRAF_V600E_vs_WT_in_SSL_HP",
  width_mm = 90,
  height_mm = 85
)

volcano_combined <- patchwork::wrap_plots(
  infection_volcano_plot,
  braf_volcano_plot,
  ncol = 2,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "Volcano plots of organoid differential-expression analyses"
  ) &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    plot.title = element_text(
      family = font_family,
      face = "bold",
      size = plot_title_pt,
      hjust = 0
    )
  )

save_plot_set(
  volcano_combined,
  subdirs$combined,
  "Volcano_plots_combined",
  width_mm = 183,
  height_mm = 90
)

writeLines(
  c(
    "VOLCANO PLOT NOTES",
    "==================",
    "",
    "Input DE tables:",
    paste0("- Infection: ", basename(files$infection_de)),
    paste0("- BRAF: ", basename(files$braf_de)),
    "",
    paste0("Display filter applied before plotting: logCPM > ", volcano_logcpm_cutoff, "."),
    paste0("Display significance category: FDR < ", volcano_fdr_cutoff,
           " and |logFC| >= ", volcano_logfc_cutoff, "."),
    paste0("Top labelled genes: top ", volcano_n_labels,
           " genes by nominal PValue after the display filter; prioritised among significant genes."),
    "The DE tables themselves retain all rows; the logCPM filter is for volcano display only.",
    "Y-axis uses -log10(FDR).",
    ""
  ),
  file.path(subdirs$qc, "Volcano_plot_methods.txt")
)

cat(
  paste0(
    "\n\nVOLCANO PLOTS ADDED\n",
    "  Infection and BRAF edgeR DE tables were used directly.\n",
    "  Display filter: logCPM > ", volcano_logcpm_cutoff, ".\n",
    "  Significance: FDR < ", volcano_fdr_cutoff,
    " and |log2FC| >= ", volcano_logfc_cutoff, ".\n",
    "  Labels: top ", volcano_n_labels,
    " genes by nominal P value after filtering, prioritised among significant genes.\n"
  ),
  file = file.path(out_dir, "RUN_NOTES.txt"),
  append = TRUE
)

message(
  "Publication-style volcano plots complete. Main files:\n",
  "  03_Volcano_plots/Volcano_Infected_vs_Non_Infected.*\n",
  "  03_Volcano_plots/Volcano_BRAF_V600E_vs_WT_in_SSL_HP.*\n",
  "  06_Combined_figures/Volcano_plots_combined.*\n",
  "  07_Tables/Infection_volcano_summary.csv\n",
  "  07_Tables/BRAF_SSL_HP_only_volcano_summary.csv\n",
  "  01_QC/Volcano_plot_methods.txt"
)

message("\nRunning final Nature-exact heatmaps and final singscore tests...")

# Final heatmap layouts; earlier renderings are retained in LEGACY_STYLE ----
final_nature_heatmap_dir <- file.path(subdirs$heatmaps, "FINAL_Nature_exact")
final_nature_combined_dir <- file.path(subdirs$combined, "FINAL_Nature_exact")
final_legacy_heatmap_dir <- file.path(subdirs$heatmaps, "LEGACY_STYLE")
final_legacy_combined_dir <- file.path(subdirs$combined, "LEGACY_STYLE")
final_legacy_score_dir <- file.path(subdirs$scores, "LEGACY_STYLE")
final_legacy_table_dir <- file.path(subdirs$tables, "LEGACY_STYLE")

for (d in c(
  final_nature_heatmap_dir,
  final_nature_combined_dir,
  final_legacy_heatmap_dir,
  final_legacy_combined_dir,
  final_legacy_score_dir,
  final_legacy_table_dir
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

copy_existing_plot_set <- function(
    source_directory,
    source_stem,
    destination_directory,
    destination_stem = paste0(source_stem, "_LEGACY_STYLE")
) {
  dir.create(destination_directory, recursive = TRUE, showWarnings = FALSE)

  copied <- character(0)
  for (extension in c("pdf", "png", "tif", "jpg", "jpeg")) {
    source_file <- file.path(
      source_directory,
      paste0(source_stem, ".", extension)
    )
    destination_file <- file.path(
      destination_directory,
      paste0(destination_stem, ".", extension)
    )

    if (file.exists(source_file)) {
      ok <- file.copy(
        source_file,
        destination_file,
        overwrite = TRUE
      )
      if (isTRUE(ok)) copied <- c(copied, destination_file)
    }
  }

  invisible(copied)
}

copy_existing_table <- function(
    source_directory,
    source_name,
    destination_directory,
    destination_name = source_name
) {
  source_file <- file.path(source_directory, source_name)
  destination_file <- file.path(destination_directory, destination_name)

  if (file.exists(source_file)) {
    file.copy(source_file, destination_file, overwrite = TRUE)
  }

  invisible(destination_file)
}

for (stem in c(
  "Hallmark_NES_heatmap_top_pathways",
  "Hallmark_NES_heatmap_all_50_pathways",
  "Infection_Hallmark_NES_heatmap_top_pathways",
  "BRAF_SSL_HP_only_Hallmark_NES_heatmap_top_pathways"
)) {
  copy_existing_plot_set(
    subdirs$heatmaps,
    stem,
    final_legacy_heatmap_dir
  )
}

for (stem in c(
  "Infection_DE_module_logFC_heatmaps_combined",
  "BRAF_SSL_HP_only_DE_module_logFC_heatmaps_combined",
  "Infection_signature_scores_combined",
  "BRAF_SSL_HP_only_signature_scores_combined",
  "BRAF_Strict_P1_SSL_HP_only_signature_scores_combined",
  "BRAF_Initial_P1_or_P2_vs_Strict_P1_signature_score_comparison"
)) {
  copy_existing_plot_set(
    subdirs$combined,
    stem,
    final_legacy_combined_dir
  )
}

for (signature_name in names(signature_list)) {
  signature_stem <- signature_file_stems[[signature_name]]

  for (stem in c(
    paste0("Infection_", signature_stem, "_singscore"),
    paste0("BRAF_SSL_HP_only_", signature_stem, "_singscore"),
    paste0("BRAF_Strict_P1_SSL_HP_only_", signature_stem, "_singscore")
  )) {
    copy_existing_plot_set(
      subdirs$scores,
      stem,
      final_legacy_score_dir
    )
  }
}

for (module_stem in c(
  "Interferon_alpha_response",
  "Chemokines",
  "Antimicrobial_peptides",
  "TNFA_NFKB_signalling",
  "TNFA_NFKB_signaling"
)) {
  for (comparison_prefix in c(
    "Infection_DE_logFC_heatmap_",
    "BRAF_SSL_HP_only_DE_logFC_heatmap_"
  )) {
    copy_existing_plot_set(
      subdirs$heatmaps,
      paste0(comparison_prefix, module_stem),
      final_legacy_heatmap_dir
    )
  }
}

for (table_name in c(
  "Infection_signature_score_statistics.csv",
  "BRAF_SSL_HP_only_signature_score_statistics.csv",
  "BRAF_Strict_P1_SSL_HP_only_signature_score_statistics.csv"
)) {
  copy_existing_table(
    subdirs$tables,
    table_name,
    final_legacy_table_dir,
    paste0(
      stringr::str_remove(table_name, "\\.csv$"),
      "_LEGACY_STYLE.csv"
    )
  )
}

final_heatmap_font_base_size <- 6
final_heatmap_axis_title_size <- 6
final_heatmap_axis_text_size <- 5.0
final_heatmap_legend_title_size <- 5.0
final_heatmap_legend_text_size <- 5.0
final_heatmap_plot_title_size <- 7
final_heatmap_thin_line <- 0.25

final_heatmap_positive_colour <- "#B2182B"
final_heatmap_negative_colour <- "#2166AC"
final_heatmap_mid_colour <- "white"
final_heatmap_tile_border_colour <- "white"
final_heatmap_tile_border_linewidth <- final_heatmap_thin_line
final_heatmap_tile_fraction <- 0.96
final_heatmap_shared_label_pt <- 5.0
final_heatmap_x_label_angle <- 55
final_heatmap_logfc_limit <- 2
final_heatmap_nes_limit <- 3
final_heatmap_max_genes_per_module <- 12L

final_ggplot_text_pt_per_mm <- 2.845276

final_geom_text_size_from_pt <- function(pt) {
  pt / final_ggplot_text_pt_per_mm
}

final_vertical_heatmap_guide <- function(
    barheight_mm = 18,
    barwidth_mm = 2.0
) {
  ggplot2::guide_colourbar(
    title.position = "right",
    title.hjust = 0.5,
    title.vjust = 0.5,
    title.theme = ggplot2::element_text(
      family = font_family,
      size = final_heatmap_legend_title_size,
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
      size = final_heatmap_legend_text_size,
      colour = "black"
    ),
    barwidth = grid::unit(barwidth_mm, "mm"),
    barheight = grid::unit(barheight_mm, "mm"),
    ticks = TRUE
  )
}

final_nature_heatmap_theme <- function(
    base_size = final_heatmap_font_base_size
) {
  ggplot2::theme_classic(
    base_size = base_size,
    base_family = font_family
  ) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = font_family,
        colour = "black"
      ),
      plot.title = ggplot2::element_text(
        family = font_family,
        face = "bold",
        hjust = 0,
        size = final_heatmap_plot_title_size
      ),
      axis.title = ggplot2::element_text(
        family = font_family,
        size = final_heatmap_axis_title_size
      ),
      axis.text = ggplot2::element_text(
        family = font_family,
        size = final_heatmap_axis_text_size,
        colour = "black"
      ),
      legend.title = ggplot2::element_text(
        family = font_family,
        size = final_heatmap_legend_title_size
      ),
      legend.text = ggplot2::element_text(
        family = font_family,
        size = final_heatmap_legend_text_size
      ),
      panel.border = ggplot2::element_rect(
        fill = NA,
        colour = "black",
        linewidth = final_heatmap_thin_line
      ),
      axis.line = ggplot2::element_line(
        linewidth = final_heatmap_thin_line
      ),
      axis.ticks = ggplot2::element_line(
        linewidth = final_heatmap_thin_line
      )
    )
}

final_significance_stars <- function(padj) {
  dplyr::case_when(
    is.na(padj) ~ "",
    padj < 0.001 ~ "***",
    padj < 0.01 ~ "**",
    padj < 0.05 ~ "*",
    TRUE ~ ""
  )
}

final_short_hallmark_label <- function(x) {
  x_upper <- toupper(as.character(x))

  fallback <- x %>%
    as.character() %>%
    stringr::str_remove("^HALLMARK_") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_lower() %>%
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
    x_upper == "HALLMARK_GLYCOLYSIS" ~
      "Glycolysis",
    x_upper == "HALLMARK_COMPLEMENT" ~
      "Complement",
    TRUE ~ fallback
  )
}

final_hallmark_fgsea <- if (
  exists("hallmark_fgsea_audited", inherits = FALSE)
) {
  hallmark_fgsea_audited
} else {
  combined_fgsea %>%
    dplyr::filter(Gene_set_family == "HALLMARK")
}

final_hallmark_comparisons <- c(
  "Infected vs Non-Infected",
  "BRAF V600E vs WT in SSL_HP"
)

final_hallmark_comparison_labels <- c(
  "Infected vs Non-Infected" = "Infected vs Non-Infected",
  "BRAF V600E vs WT in SSL_HP" = "BRAF V600E vs WT in SSL/HP"
)

make_final_nature_hallmark_heatmap <- function(
    fgsea_df,
    comparisons,
    pathways = NULL,
    top_n_per_comparison = NULL,
    title_text = NULL,
    legend_title = "NES",
    nes_limit = final_heatmap_nes_limit
) {
  dat <- fgsea_df %>%
    dplyr::filter(
      Comparison %in% comparisons,
      Gene_set_family == "HALLMARK",
      is.finite(NES)
    ) %>%
    dplyr::group_by(Comparison, pathway) %>%
    dplyr::arrange(
      padj,
      pval,
      dplyr::desc(abs(NES)),
      .by_group = TRUE
    ) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  if (nrow(dat) == 0L) {
    return(make_empty_plot(
      ifelse(is.null(title_text), "Hallmark GSEA", title_text),
      "No Hallmark FGSEA results were available"
    ))
  }

  if (!is.null(pathways)) {
    selected_pathways <- unique(pathways[pathways %in% dat$pathway])
  } else if (!is.null(top_n_per_comparison)) {
    selected_pathways <- dat %>%
      dplyr::group_by(Comparison) %>%
      dplyr::arrange(
        padj,
        dplyr::desc(abs(NES)),
        .by_group = TRUE
      ) %>%
      dplyr::slice_head(n = top_n_per_comparison) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(pathway) %>%
      dplyr::pull(pathway)
  } else {
    selected_pathways <- unique(dat$pathway)
  }

  if (length(selected_pathways) == 0L) {
    return(make_empty_plot(
      ifelse(is.null(title_text), "Hallmark GSEA", title_text),
      "No requested Hallmark pathways were available"
    ))
  }

  pathway_order <- dat %>%
    dplyr::filter(pathway %in% selected_pathways) %>%
    dplyr::group_by(pathway) %>%
    dplyr::summarise(
      Max_abs_NES = max(abs(NES), na.rm = TRUE),
      Best_FDR = min(padj, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      dplyr::desc(Max_abs_NES),
      Best_FDR,
      pathway
    ) %>%
    dplyr::pull(pathway)

  comparison_display_levels <- unname(
    final_hallmark_comparison_labels[comparisons]
  )
  comparison_display_levels[
    is.na(comparison_display_levels)
  ] <- comparisons[is.na(comparison_display_levels)]

  plot_df <- tidyr::expand_grid(
    Comparison = comparisons,
    pathway = pathway_order
  ) %>%
    dplyr::left_join(
      dat %>%
        dplyr::select(
          Comparison,
          pathway,
          NES,
          pval,
          padj
        ),
      by = c("Comparison", "pathway")
    ) %>%
    dplyr::mutate(
      Comparison_display_text = unname(
        final_hallmark_comparison_labels[Comparison]
      ),
      Comparison_display_text = dplyr::if_else(
        is.na(Comparison_display_text),
        Comparison,
        Comparison_display_text
      ),
      Comparison_display = factor(
        Comparison_display_text,
        levels = comparison_display_levels
      ),
      Pathway_label_text = final_short_hallmark_label(pathway),
      Pathway_label = factor(
        Pathway_label_text,
        levels = rev(final_short_hallmark_label(pathway_order))
      ),
      Significance = final_significance_stars(padj)
    )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = Comparison_display,
      y = Pathway_label,
      fill = NES
    )
  ) +
    ggplot2::geom_tile(
      colour = final_heatmap_tile_border_colour,
      linewidth = final_heatmap_tile_border_linewidth,
      width = final_heatmap_tile_fraction,
      height = final_heatmap_tile_fraction
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = Significance),
      family = font_family,
      fontface = "bold",
      colour = "white",
      size = final_geom_text_size_from_pt(
        final_heatmap_shared_label_pt
      )
    ) +
    ggplot2::scale_x_discrete(
      drop = FALSE,
      expand = ggplot2::expansion(add = 0.05)
    ) +
    ggplot2::scale_y_discrete(
      drop = FALSE,
      labels = function(x) stringr::str_wrap(x, width = 23),
      expand = ggplot2::expansion(add = 0.08)
    ) +
    ggplot2::scale_fill_gradient2(
      low = final_heatmap_negative_colour,
      mid = final_heatmap_mid_colour,
      high = final_heatmap_positive_colour,
      midpoint = 0,
      limits = c(-nes_limit, nes_limit),
      oob = scales::squish,
      na.value = "grey95",
      name = legend_title,
      guide = final_vertical_heatmap_guide(
        barheight_mm = 18,
        barwidth_mm = 2.0
      )
    ) +
    ggplot2::coord_fixed(ratio = 1, clip = "off") +
    ggplot2::labs(
      title = format_signature_plot_title(title_text),
      x = NULL,
      y = NULL
    ) +
    final_nature_heatmap_theme(
      base_size = final_heatmap_font_base_size
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        family = font_family,
        size = final_heatmap_shared_label_pt,
        angle = final_heatmap_x_label_angle,
        hjust = 1,
        vjust = 1,
        colour = "black"
      ),
      axis.text.y = ggplot2::element_text(
        family = font_family,
        size = final_heatmap_shared_label_pt,
        lineheight = 0.9,
        colour = "black"
      ),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      legend.direction = "vertical",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(2, 2, 2, 2, unit = "mm")
    )
}

final_hallmark_key_pathways <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING"
)

p_final_hallmark_top <- make_final_nature_hallmark_heatmap(
  final_hallmark_fgsea,
  comparisons = final_hallmark_comparisons,
  top_n_per_comparison = hallmark_heatmap_top_n,
  title_text = NULL
)

p_final_hallmark_all <- make_final_nature_hallmark_heatmap(
  final_hallmark_fgsea,
  comparisons = final_hallmark_comparisons,
  pathways = sort(unique(final_hallmark_fgsea$pathway)),
  title_text = NULL
)

p_final_hallmark_key <- make_final_nature_hallmark_heatmap(
  final_hallmark_fgsea,
  comparisons = final_hallmark_comparisons,
  pathways = final_hallmark_key_pathways,
  title_text = NULL
)

p_final_hallmark_infection <- make_final_nature_hallmark_heatmap(
  final_hallmark_fgsea,
  comparisons = "Infected vs Non-Infected",
  top_n_per_comparison = hallmark_heatmap_top_n,
  title_text = "Hallmark GSEA: Infected vs Non-Infected"
)

p_final_hallmark_braf <- make_final_nature_hallmark_heatmap(
  final_hallmark_fgsea,
  comparisons = "BRAF V600E vs WT in SSL_HP",
  top_n_per_comparison = hallmark_heatmap_top_n,
  title_text = "Hallmark GSEA: BRAF V600E vs WT in SSL/HP"
)

save_plot_set(
  p_final_hallmark_top,
  final_nature_heatmap_dir,
  "FINAL_Hallmark_NES_heatmap_top_pathways_NATURE_EXACT",
  125,
  125
)
save_plot_set(
  p_final_hallmark_all,
  final_nature_heatmap_dir,
  "FINAL_Hallmark_NES_heatmap_all_pathways_NATURE_EXACT",
  135,
  205
)
save_plot_set(
  p_final_hallmark_key,
  final_nature_heatmap_dir,
  "FINAL_Hallmark_NES_heatmap_key_pathways_NATURE_EXACT",
  110,
  92
)
save_plot_set(
  p_final_hallmark_infection,
  final_nature_heatmap_dir,
  "FINAL_Infection_Hallmark_NES_heatmap_NATURE_EXACT",
  88,
  125
)
save_plot_set(
  p_final_hallmark_braf,
  final_nature_heatmap_dir,
  "FINAL_BRAF_Hallmark_NES_heatmap_NATURE_EXACT",
  88,
  125
)

save_plot_set(
  p_final_hallmark_top,
  subdirs$heatmaps,
  "Hallmark_NES_heatmap_top_pathways",
  125,
  125
)
save_plot_set(
  p_final_hallmark_all,
  subdirs$heatmaps,
  "Hallmark_NES_heatmap_all_50_pathways",
  135,
  205
)
save_plot_set(
  p_final_hallmark_infection,
  subdirs$heatmaps,
  "Infection_Hallmark_NES_heatmap_top_pathways",
  88,
  125
)
save_plot_set(
  p_final_hallmark_braf,
  subdirs$heatmaps,
  "BRAF_SSL_HP_only_Hallmark_NES_heatmap_top_pathways",
  88,
  125
)

final_gene_module_info <- list(
  "Interferon alpha response" = ifna_info$genes,
  "Chemokines" = chemokine_genes,
  "Antimicrobial peptides" = antimicrobial_info$genes,
  "TNFA/NFKB signaling" = tnfa_info$genes
)

final_gene_group_order <- names(final_gene_module_info)

final_module_comparison_labels <- c(
  Infection = "Infected vs Non-Infected",
  BRAF_initial = "BRAF V600E vs WT in SSL/HP"
)

final_module_de_tables <- list(
  Infection = infection_de %>%
    dplyr::mutate(comparison = "Infection"),
  BRAF_initial = braf_de %>%
    dplyr::mutate(comparison = "BRAF_initial")
)

final_reference_de <- final_module_de_tables[["Infection"]]

final_reference_for_join <- final_reference_de %>%
  dplyr::select(
    Gene,
    Reference_logFC = logFC,
    Reference_logCPM = logCPM,
    dplyr::any_of(c("PValue", "FDR"))
  )

if ("PValue" %in% names(final_reference_for_join)) {
  names(final_reference_for_join)[
    names(final_reference_for_join) == "PValue"
  ] <- "Reference_PValue"
}
if ("FDR" %in% names(final_reference_for_join)) {
  names(final_reference_for_join)[
    names(final_reference_for_join) == "FDR"
  ] <- "Reference_FDR"
}

final_module_candidate_qc <- purrr::imap_dfr(
  final_gene_module_info,
  function(genes, group_name) {
    tibble::tibble(
      gene_group = group_name,
      Gene = unique(toupper(trimws(genes)))
    ) %>%
      dplyr::filter(!is.na(Gene), Gene != "") %>%
      dplyr::left_join(
        final_reference_for_join,
        by = "Gene"
      ) %>%
      dplyr::mutate(
        passes_reference_filter =
          !is.na(Reference_logCPM) &
          is.finite(Reference_logCPM) &
          Reference_logCPM >= module_heatmap_logcpm_cutoff &
          !is.na(Reference_logFC) &
          is.finite(Reference_logFC) &
          Reference_logFC > module_heatmap_logfc_cutoff
      )
  }
)

final_reference_selected_all <- final_module_candidate_qc %>%
  dplyr::filter(passes_reference_filter) %>%
  dplyr::arrange(
    factor(gene_group, levels = final_gene_group_order),
    dplyr::desc(Reference_logFC),
    Gene
  ) %>%
  dplyr::group_by(gene_group) %>%
  dplyr::mutate(selection_rank = dplyr::row_number()) %>%
  dplyr::ungroup()

final_reference_display_genes <- final_reference_selected_all %>%
  dplyr::group_by(gene_group) %>%
  dplyr::filter(
    selection_rank <= final_heatmap_max_genes_per_module
  ) %>%
  dplyr::ungroup()

final_missing_modules <- setdiff(
  final_gene_group_order,
  unique(final_reference_display_genes$gene_group)
)

if (length(final_missing_modules) > 0L) {
  warning(
    "No infection-reference genes passed logCPM >= ",
    module_heatmap_logcpm_cutoff,
    " and logFC > ",
    module_heatmap_logfc_cutoff,
    " for: ",
    paste(final_missing_modules, collapse = ", "),
    ". The affected module cannot be populated without relaxing the ",
    "pre-specified Nature/tissue selection rule."
  )
}

final_gene_groups <- split(
  final_reference_display_genes$Gene,
  final_reference_display_genes$gene_group
)

for (group_name in final_gene_group_order) {
  if (!group_name %in% names(final_gene_groups)) {
    final_gene_groups[[group_name]] <- character(0)
  }
}

final_gene_groups <- final_gene_groups[final_gene_group_order]

final_braf_for_join <- final_module_de_tables[["BRAF_initial"]] %>%
  dplyr::select(
    Gene,
    BRAF_logFC = logFC,
    BRAF_logCPM = logCPM,
    dplyr::any_of(c("PValue", "FDR"))
  )

if ("PValue" %in% names(final_braf_for_join)) {
  names(final_braf_for_join)[
    names(final_braf_for_join) == "PValue"
  ] <- "BRAF_PValue"
}
if ("FDR" %in% names(final_braf_for_join)) {
  names(final_braf_for_join)[
    names(final_braf_for_join) == "FDR"
  ] <- "BRAF_FDR"
}

final_braf_module_candidate_qc <- purrr::imap_dfr(
  final_gene_module_info,
  function(genes, group_name) {
    tibble::tibble(
      gene_group = group_name,
      Gene = unique(toupper(trimws(genes)))
    ) %>%
      dplyr::filter(!is.na(Gene), Gene != "") %>%
      dplyr::left_join(
        final_braf_for_join,
        by = "Gene"
      ) %>%
      dplyr::mutate(
        passes_braf_filter =
          !is.na(BRAF_logCPM) &
          is.finite(BRAF_logCPM) &
          BRAF_logCPM >= module_heatmap_logcpm_cutoff &
          !is.na(BRAF_logFC) &
          is.finite(BRAF_logFC) &
          abs(BRAF_logFC) > module_heatmap_logfc_cutoff
      )
  }
)

final_braf_selected_all <- final_braf_module_candidate_qc %>%
  dplyr::filter(passes_braf_filter) %>%
  dplyr::arrange(
    factor(gene_group, levels = final_gene_group_order),
    dplyr::desc(abs(BRAF_logFC)),
    Gene
  ) %>%
  dplyr::group_by(gene_group) %>%
  dplyr::mutate(selection_rank = dplyr::row_number()) %>%
  dplyr::ungroup()

final_braf_display_genes <- final_braf_selected_all %>%
  dplyr::group_by(gene_group) %>%
  dplyr::filter(
    selection_rank <= final_heatmap_max_genes_per_module
  ) %>%
  dplyr::ungroup()

final_braf_gene_groups <- split(
  final_braf_display_genes$Gene,
  final_braf_display_genes$gene_group
)

for (group_name in final_gene_group_order) {
  if (!group_name %in% names(final_braf_gene_groups)) {
    final_braf_gene_groups[[group_name]] <- character(0)
  }
}

final_braf_gene_groups <- final_braf_gene_groups[
  final_gene_group_order
]

make_final_heatmap_df <- function(
    comparisons,
    genes,
    group_name
) {
  de_all <- purrr::map_dfr(
    comparisons,
    function(comparison_key) {
      if (!comparison_key %in% names(final_module_de_tables)) {
        stop("Unknown module-heatmap comparison: ", comparison_key)
      }

      final_module_de_tables[[comparison_key]] %>%
        dplyr::select(
          comparison,
          Gene,
          logFC,
          dplyr::any_of(c("logCPM", "PValue", "FDR"))
        )
    }
  )

  gene_index <- tibble::tibble(
    Gene = genes,
    gene_order = seq_along(genes),
    gene_id = paste(
      group_name,
      genes,
      seq_along(genes),
      sep = "__"
    )
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
        unname(final_module_comparison_labels[comparison]),
        levels = unname(
          final_module_comparison_labels[comparisons]
        )
      ),
      gene_id = factor(
        gene_id,
        levels = rev(gene_index$gene_id)
      ),
      logFC_colour = pmax(
        pmin(logFC, final_heatmap_logfc_limit),
        -final_heatmap_logfc_limit
      )
    )
}

build_final_one_row_heatmap_panel <- function(
    df_list,
    comparisons,
    legend_title = "log2FC",
    selected_modules = NULL,
    force_block_height = NULL,
    show_legend = TRUE
) {
  module_names <- final_gene_group_order[
    final_gene_group_order %in% names(df_list)
  ]

  if (!is.null(selected_modules)) {
    module_names <- module_names[
      module_names %in% selected_modules
    ]
  }

  module_names <- module_names[
    vapply(
      module_names,
      function(group_name) {
        !is.null(df_list[[group_name]]) &&
          nrow(df_list[[group_name]]) > 0
      },
      logical(1)
    )
  ]

  if (length(module_names) == 0L) {
    return(make_empty_plot(
      "Gene-level heatmaps",
      "No genes passed the selection rule"
    ))
  }

  ncomp <- length(comparisons)
  heatmap_labels <- unname(
    final_module_comparison_labels[comparisons]
  )
  max_label_len <- max(nchar(heatmap_labels))

  module_n_tbl <- purrr::map_dfr(
    module_names,
    function(group_name) {
      tibble::tibble(
        gene_group = group_name,
        n_genes = dplyr::n_distinct(
          df_list[[group_name]]$Gene
        )
      )
    }
  )

  if (!is.null(force_block_height) &&
      is.finite(force_block_height)) {
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
    dplyr::left_join(
      module_n_tbl,
      by = "gene_group"
    ) %>%
    dplyr::mutate(
      x_offset = x_start +
        (layout_index - 1) * (ncomp + x_gap),
      y_top_boundary = block_height,
      y_bottom_boundary = 0
    )

  tile_half_width <- final_heatmap_tile_fraction / 2

  tile_data <- purrr::map_dfr(
    module_names,
    function(group_name) {
      df <- df_list[[group_name]]
      lay <- layout_tbl %>%
        dplyr::filter(gene_group == group_name)

      genes_in_order <- df %>%
        dplyr::distinct(Gene, gene_order) %>%
        dplyr::arrange(gene_order) %>%
        dplyr::mutate(
          gene_order_plot = dplyr::row_number()
        )

      comp_tbl <- tibble::tibble(
        comparison_label = factor(
          heatmap_labels,
          levels = heatmap_labels
        ),
        x_local = seq_along(comparisons)
      )

      row_height <- block_height / max(lay$n_genes, 1)

      df %>%
        dplyr::left_join(
          genes_in_order,
          by = c("Gene", "gene_order")
        ) %>%
        dplyr::left_join(
          comp_tbl,
          by = "comparison_label"
        ) %>%
        dplyr::mutate(
          x_plot = lay$x_offset + x_local - 1,
          y_plot = lay$y_top_boundary -
            (gene_order_plot - 0.5) * row_height,
          tile_height =
            row_height * final_heatmap_tile_fraction,
          xmin = x_plot - tile_half_width,
          xmax = x_plot + tile_half_width,
          ymin = y_plot - tile_height / 2,
          ymax = y_plot + tile_height / 2,
          gene_group = group_name
        )
    }
  )

  gene_label_data <- tile_data %>%
    dplyr::distinct(
      gene_group,
      Gene,
      y_plot
    ) %>%
    dplyr::left_join(
      layout_tbl,
      by = "gene_group"
    ) %>%
    dplyr::mutate(
      x_plot = x_offset - 0.62
    )

  module_title_data <- layout_tbl %>%
    dplyr::mutate(
      x_plot = x_offset + (ncomp - 1) / 2,
      y_plot = y_top_boundary + 0.82,
      title_label = stringr::str_wrap(
        gene_group,
        width = 18
      )
    )

  x_label_data <- purrr::map_dfr(
    module_names,
    function(group_name) {
      lay <- layout_tbl %>%
        dplyr::filter(gene_group == group_name)

      tibble::tibble(
        gene_group = group_name,
        comparison_label = heatmap_labels,
        x_plot = lay$x_offset +
          seq_along(comparisons) - 1,
        y_plot = lay$y_bottom_boundary - 0.52
      )
    }
  )

  x_min <- min(layout_tbl$x_offset) - 3.05
  x_max <- max(layout_tbl$x_offset) + ncomp - 0.05

  bottom_extra <- if (max_label_len > 18) 5.8 else 4.8
  if (max_label_len > 24) bottom_extra <- 6.6

  y_min <- min(x_label_data$y_plot) - bottom_extra
  y_max <- max(layout_tbl$y_top_boundary) + 1.55

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
      colour = final_heatmap_tile_border_colour,
      linewidth = final_heatmap_tile_border_linewidth
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
      size = final_geom_text_size_from_pt(
        final_heatmap_shared_label_pt
      ),
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
      size = final_geom_text_size_from_pt(
        final_heatmap_plot_title_size
      ),
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
      angle = final_heatmap_x_label_angle,
      family = font_family,
      size = final_geom_text_size_from_pt(
        final_heatmap_shared_label_pt
      ),
      colour = "black"
    ) +
    ggplot2::scale_fill_gradient2(
      low = final_heatmap_negative_colour,
      mid = final_heatmap_mid_colour,
      high = final_heatmap_positive_colour,
      midpoint = 0,
      limits = c(
        -final_heatmap_logfc_limit,
        final_heatmap_logfc_limit
      ),
      breaks = c(
        -final_heatmap_logfc_limit,
        0,
        final_heatmap_logfc_limit
      ),
      labels = c(
        paste0("≤ -", final_heatmap_logfc_limit),
        "0",
        paste0("≥ +", final_heatmap_logfc_limit)
      ),
      oob = scales::squish,
      na.value = "grey90",
      name = legend_title,
      guide = final_vertical_heatmap_guide(
        barheight_mm = 18,
        barwidth_mm = 2.0
      )
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
        size = final_heatmap_legend_title_size,
        angle = 90,
        hjust = 0.5,
        vjust = 0.5,
        colour = "black"
      ),
      legend.text = ggplot2::element_text(
        family = font_family,
        size = final_heatmap_legend_text_size,
        colour = "black"
      ),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(
        2,
        3,
        2,
        3,
        unit = "mm"
      )
    )
}

save_final_module_heatmap_sets <- function(
    df_list,
    comparisons,
    stem_prefix,
    directory = final_nature_heatmap_dir,
    legend_title = "log2FC",
    individual_width_mm = 60,
    individual_height_mm = 105
) {
  module_names <- final_gene_group_order[
    final_gene_group_order %in% names(df_list)
  ]

  module_names <- module_names[
    vapply(
      module_names,
      function(group_name) {
        !is.null(df_list[[group_name]]) &&
          nrow(df_list[[group_name]]) > 0
      },
      logical(1)
    )
  ]

  if (length(module_names) == 0L) {
    return(list())
  }

  module_n_tbl <- purrr::map_dfr(
    module_names,
    function(group_name) {
      tibble::tibble(
        gene_group = group_name,
        n_genes = dplyr::n_distinct(
          df_list[[group_name]]$Gene
        )
      )
    }
  )

  common_block_height <- max(
    module_n_tbl$n_genes,
    na.rm = TRUE
  )

  plots <- purrr::map(
    module_names,
    function(group_name) {
      p <- build_final_one_row_heatmap_panel(
        df_list = df_list[group_name],
        comparisons = comparisons,
        legend_title = legend_title,
        selected_modules = group_name,
        force_block_height = common_block_height,
        show_legend = TRUE
      )

      save_plot_set(
        p,
        directory,
        paste0(
          stem_prefix,
          "_",
          safe_filename(group_name)
        ),
        individual_width_mm,
        individual_height_mm
      )

      p
    }
  )

  names(plots) <- module_names
  plots
}

final_module_df_list_both <- list()
final_module_df_list_infection <- list()
final_module_df_list_braf_reference <- list()
final_module_df_list_braf_specific <- list()

for (group_name in names(final_gene_groups)) {
  genes <- final_gene_groups[[group_name]]
  if (length(genes) == 0L) next

  final_module_df_list_both[[group_name]] <-
    make_final_heatmap_df(
      comparisons = c("Infection", "BRAF_initial"),
      genes = genes,
      group_name = group_name
    )

  final_module_df_list_infection[[group_name]] <-
    make_final_heatmap_df(
      comparisons = "Infection",
      genes = genes,
      group_name = group_name
    )

  final_module_df_list_braf_reference[[group_name]] <-
    make_final_heatmap_df(
      comparisons = "BRAF_initial",
      genes = genes,
      group_name = group_name
    )
}

for (group_name in names(final_braf_gene_groups)) {
  genes <- final_braf_gene_groups[[group_name]]
  if (length(genes) == 0L) next

  final_module_df_list_braf_specific[[group_name]] <-
    make_final_heatmap_df(
      comparisons = "BRAF_initial",
      genes = genes,
      group_name = group_name
    )
}

p_final_modules_both <- build_final_one_row_heatmap_panel(
  final_module_df_list_both,
  comparisons = c("Infection", "BRAF_initial"),
  legend_title = "log2FC",
  show_legend = TRUE
)

p_final_modules_infection <- build_final_one_row_heatmap_panel(
  final_module_df_list_infection,
  comparisons = "Infection",
  legend_title = "log2FC",
  show_legend = TRUE
)

p_final_modules_braf_reference <- build_final_one_row_heatmap_panel(
  final_module_df_list_braf_reference,
  comparisons = "BRAF_initial",
  legend_title = "log2FC",
  show_legend = TRUE
)

p_final_modules_braf_specific <- build_final_one_row_heatmap_panel(
  final_module_df_list_braf_specific,
  comparisons = "BRAF_initial",
  legend_title = "log2FC",
  show_legend = TRUE
)

save_plot_set(
  p_final_modules_both,
  final_nature_heatmap_dir,
  "FINAL_ALL4_Module_log2FC_Heatmaps_Infection_and_BRAF_NATURE_EXACT",
  183,
  82
)
save_plot_set(
  p_final_modules_infection,
  final_nature_heatmap_dir,
  "FINAL_ALL4_Module_log2FC_Heatmaps_Infection_NATURE_EXACT",
  183,
  82
)
save_plot_set(
  p_final_modules_braf_reference,
  final_nature_heatmap_dir,
  "FINAL_ALL4_Module_log2FC_Heatmaps_BRAF_Initial_INFECTION_REFERENCE_GENES_NATURE_EXACT",
  183,
  82
)
save_plot_set(
  p_final_modules_braf_specific,
  final_nature_heatmap_dir,
  "FINAL_ALL4_Module_log2FC_Heatmaps_BRAF_Initial_BRAF_SELECTED_NATURE_EXACT",
  183,
  82
)

save_plot_set(
  p_final_modules_both,
  subdirs$heatmaps,
  "ALL4_Module_log2FC_Heatmaps_Infection_and_BRAF_FINAL",
  183,
  82
)
save_plot_set(
  p_final_modules_infection,
  subdirs$heatmaps,
  "ALL4_Module_log2FC_Heatmaps_Infection_FINAL",
  183,
  82
)
save_plot_set(
  p_final_modules_braf_specific,
  subdirs$heatmaps,
  "ALL4_Module_log2FC_Heatmaps_BRAF_Initial_FINAL",
  183,
  82
)
save_plot_set(
  p_final_modules_braf_reference,
  subdirs$heatmaps,
  "ALL4_Module_log2FC_Heatmaps_BRAF_Initial_INFECTION_REFERENCE_FINAL",
  183,
  82
)

save_plot_set(
  p_final_modules_infection,
  subdirs$combined,
  "Infection_DE_module_logFC_heatmaps_combined",
  183,
  82
)
save_plot_set(
  p_final_modules_braf_specific,
  subdirs$combined,
  "BRAF_SSL_HP_only_DE_module_logFC_heatmaps_combined",
  183,
  82
)

final_individual_module_plots_both <- save_final_module_heatmap_sets(
  final_module_df_list_both,
  comparisons = c("Infection", "BRAF_initial"),
  stem_prefix = "FINAL_BOTH_COMPARISONS_SINGLE_MODULE_NATURE_EXACT",
  directory = final_nature_heatmap_dir,
  individual_width_mm = 70,
  individual_height_mm = 105
)

final_individual_module_plots_infection <- save_final_module_heatmap_sets(
  final_module_df_list_infection,
  comparisons = "Infection",
  stem_prefix = "FINAL_INFECTION_SINGLE_MODULE_NATURE_EXACT",
  directory = final_nature_heatmap_dir,
  individual_width_mm = 60,
  individual_height_mm = 105
)

final_individual_module_plots_braf_reference <- save_final_module_heatmap_sets(
  final_module_df_list_braf_reference,
  comparisons = "BRAF_initial",
  stem_prefix = "FINAL_BRAF_INITIAL_INFECTION_REFERENCE_SINGLE_MODULE_NATURE_EXACT",
  directory = final_nature_heatmap_dir,
  individual_width_mm = 60,
  individual_height_mm = 105
)

final_individual_module_plots_braf <- save_final_module_heatmap_sets(
  final_module_df_list_braf_specific,
  comparisons = "BRAF_initial",
  stem_prefix = "FINAL_BRAF_INITIAL_BRAF_SELECTED_SINGLE_MODULE_NATURE_EXACT",
  directory = final_nature_heatmap_dir,
  individual_width_mm = 60,
  individual_height_mm = 105
)

final_standard_module_stem <- c(
  "Interferon alpha response" = "Interferon_alpha_response",
  "Chemokines" = "Chemokines",
  "Antimicrobial peptides" = "Antimicrobial_peptides",
  "TNFA/NFKB signaling" = "TNFA_NFKB_signalling"
)

for (group_name in names(final_standard_module_stem)) {
  standard_stem <- final_standard_module_stem[[group_name]]

  if (group_name %in% names(final_individual_module_plots_infection)) {
    save_plot_set(
      final_individual_module_plots_infection[[group_name]],
      subdirs$heatmaps,
      paste0(
        "Infection_DE_logFC_heatmap_",
        standard_stem
      ),
      60,
      105
    )
  }

  if (group_name %in% names(final_individual_module_plots_braf)) {
    save_plot_set(
      final_individual_module_plots_braf[[group_name]],
      subdirs$heatmaps,
      paste0(
        "BRAF_SSL_HP_only_DE_logFC_heatmap_",
        standard_stem
      ),
      60,
      105
    )
  }
}

if ("TNFA/NFKB signaling" %in%
    names(final_individual_module_plots_infection)) {
  save_plot_set(
    final_individual_module_plots_infection[[
      "TNFA/NFKB signaling"
    ]],
    subdirs$heatmaps,
    "Infection_DE_logFC_heatmap_TNFA_NFKB_signaling",
    60,
    105
  )
}
if ("TNFA/NFKB signaling" %in%
    names(final_individual_module_plots_braf)) {
  save_plot_set(
    final_individual_module_plots_braf[[
      "TNFA/NFKB signaling"
    ]],
    subdirs$heatmaps,
    "BRAF_SSL_HP_only_DE_logFC_heatmap_TNFA_NFKB_signaling",
    60,
    105
  )
}

final_module_long_both <- if (
  length(final_module_df_list_both) > 0L
) {
  dplyr::bind_rows(
    final_module_df_list_both,
    .id = "Module_list_name"
  )
} else {
  tibble::tibble()
}

data.table::fwrite(
  final_module_candidate_qc,
  file.path(
    subdirs$tables,
    "FINAL_NATURE_EXACT_module_candidate_gene_QC.csv"
  )
)

data.table::fwrite(
  final_reference_selected_all,
  file.path(
    subdirs$tables,
    "FINAL_NATURE_EXACT_ALL_genes_passing_infection_reference_rule.csv"
  )
)

data.table::fwrite(
  final_reference_display_genes,
  file.path(
    subdirs$tables,
    "FINAL_NATURE_EXACT_DISPLAYED_genes_top12_from_infection_reference.csv"
  )
)

data.table::fwrite(
  final_braf_module_candidate_qc,
  file.path(
    subdirs$tables,
    "FINAL_NATURE_EXACT_BRAF_module_candidate_gene_QC.csv"
  )
)

data.table::fwrite(
  final_braf_display_genes,
  file.path(
    subdirs$tables,
    "FINAL_NATURE_EXACT_BRAF_DISPLAYED_genes_top12_by_abs_logFC.csv"
  )
)

data.table::fwrite(
  final_module_long_both,
  file.path(
    subdirs$tables,
    "FINAL_NATURE_EXACT_ALL4_heatmap_long_Infection_and_BRAF.csv"
  )
)

final_heatmap_summary <- p_final_hallmark_key /
  p_final_modules_both +
  patchwork::plot_layout(
    heights = c(0.82, 1.18)
  ) +
  patchwork::plot_annotation(
    tag_levels = "a",
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(
        family = font_family,
        face = "bold",
        size = panel_tag_pt
      )
    )
  )

save_plot_set(
  final_heatmap_summary,
  final_nature_combined_dir,
  "FINAL_Organoid_Hallmark_and_ALL4_gene_modules_NATURE_EXACT",
  183,
  170
)

# Final score statistics: paired infection tests and unpaired BRAF tests ----
safe_final_t_test <- function(...) {
  tryCatch(
    stats::t.test(...),
    error = function(e) NULL
  )
}

safe_final_wilcox <- function(...) {
  suppressWarnings(
    tryCatch(
      stats::wilcox.test(...)$p.value,
      error = function(e) NA_real_
    )
  )
}

prepare_final_infection_test <- function(
    score_data,
    signature_name
) {
  dat <- score_data %>%
    dplyr::filter(
      Signature == signature_name,
      Infection_Status %in% c("Non-Infected", "Infected"),
      is.finite(Score),
      !is.na(Pair_ID),
      Pair_ID != ""
    )

  event_condition_scores <- dat %>%
    dplyr::group_by(
      Pair_ID,
      Infection_Status
    ) %>%
    dplyr::summarise(
      Score = mean(Score, na.rm = TRUE),
      N_technical_replicates = dplyr::n(),
      Line_ID = first_nonmissing_scalar(Line_ID),
      Trajectory_ID = first_nonmissing_scalar(Trajectory_ID),
      Passage = first_nonmissing_scalar(Passage),
      .groups = "drop"
    )

  complete_event_ids <- event_condition_scores %>%
    dplyr::group_by(Pair_ID) %>%
    dplyr::summarise(
      N_conditions = dplyr::n_distinct(Infection_Status),
      N_rows = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      N_conditions == 2L,
      N_rows == 2L
    ) %>%
    dplyr::pull(Pair_ID)

  paired_long <- event_condition_scores %>%
    dplyr::filter(Pair_ID %in% complete_event_ids) %>%
    dplyr::mutate(
      Infection_Status = factor(
        Infection_Status,
        levels = c("Non-Infected", "Infected")
      )
    )

  paired_wide <- paired_long %>%
    dplyr::select(
      Pair_ID,
      Infection_Status,
      Score,
      Line_ID,
      Trajectory_ID,
      Passage
    ) %>%
    tidyr::pivot_wider(
      id_cols = c(
        Pair_ID,
        Line_ID,
        Trajectory_ID,
        Passage
      ),
      names_from = Infection_Status,
      values_from = Score
    ) %>%
    dplyr::filter(
      is.finite(.data[["Non-Infected"]]),
      is.finite(.data[["Infected"]])
    ) %>%
    dplyr::mutate(
      Difference_infected_minus_noninfected =
        .data[["Infected"]] -
        .data[["Non-Infected"]]
    )

  n_events <- nrow(paired_wide)
  noninfected <- paired_wide[["Non-Infected"]]
  infected <- paired_wide[["Infected"]]

  paired_t <- if (n_events >= 2L) {
    safe_final_t_test(
      infected,
      noninfected,
      paired = TRUE,
      alternative = "two.sided",
      conf.level = 0.95
    )
  } else {
    NULL
  }

  exact_wilcoxon <- if (n_events >= 2L) {
    safe_final_wilcox(
      infected,
      noninfected,
      paired = TRUE,
      exact = TRUE,
      correct = FALSE,
      alternative = "two.sided"
    )
  } else {
    NA_real_
  }

  asymptotic_wilcoxon <- if (n_events >= 2L) {
    safe_final_wilcox(
      infected,
      noninfected,
      paired = TRUE,
      exact = FALSE,
      correct = FALSE,
      alternative = "two.sided"
    )
  } else {
    NA_real_
  }

  primary_p <- if (is.null(paired_t)) {
    NA_real_
  } else {
    paired_t$p.value
  }

  statistics <- data.frame(
    Dataset = "Infection",
    Signature = signature_name,
    Experimental_unit =
      "Passage-specific independent infection event (Pair_ID)",
    Primary_test = "Two-sided paired t-test",
    N_independent_infection_events = n_events,
    N_noninfected = n_events,
    N_infected = n_events,
    Mean_noninfected = ifelse(
      n_events > 0L,
      mean(noninfected),
      NA_real_
    ),
    Mean_infected = ifelse(
      n_events > 0L,
      mean(infected),
      NA_real_
    ),
    Mean_difference_infected_minus_noninfected = ifelse(
      n_events > 0L,
      mean(
        paired_wide$Difference_infected_minus_noninfected
      ),
      NA_real_
    ),
    Median_difference_infected_minus_noninfected = ifelse(
      n_events > 0L,
      stats::median(
        paired_wide$Difference_infected_minus_noninfected
      ),
      NA_real_
    ),
    P_value = primary_p,
    Paired_t_95CI_lower = if (is.null(paired_t)) {
      NA_real_
    } else {
      unname(paired_t$conf.int[1])
    },
    Paired_t_95CI_upper = if (is.null(paired_t)) {
      NA_real_
    } else {
      unname(paired_t$conf.int[2])
    },
    Exact_paired_Wilcoxon_P_value = exact_wilcoxon,
    Asymptotic_paired_Wilcoxon_P_value =
      asymptotic_wilcoxon,
    Points_connected_in_plot = FALSE,
    stringsAsFactors = FALSE
  )

  list(
    plot_data = paired_long,
    paired_wide = paired_wide,
    statistics = statistics
  )
}

prepare_final_braf_test <- function(
    score_data,
    signature_name,
    analysis_label
) {
  dat <- score_data %>%
    dplyr::filter(
      Signature == signature_name,
      BRAF_Status %in% c("WT", "V600E"),
      is.finite(Score)
    ) %>%
    dplyr::mutate(
      BRAF_Status = factor(
        BRAF_Status,
        levels = c("WT", "V600E")
      )
    )

  wt <- dat$Score[dat$BRAF_Status == "WT"]
  v600e <- dat$Score[dat$BRAF_Status == "V600E"]
  n_wt <- length(wt)
  n_v600e <- length(v600e)

  welch_t <- if (n_wt >= 2L && n_v600e >= 2L) {
    safe_final_t_test(
      v600e,
      wt,
      paired = FALSE,
      var.equal = FALSE,
      alternative = "two.sided",
      conf.level = 0.95
    )
  } else {
    NULL
  }

  exact_wilcoxon <- if (
    n_wt >= 2L &&
    n_v600e >= 2L &&
    !anyDuplicated(c(wt, v600e)) &&
    n_wt * n_v600e <= 200L
  ) {
    safe_final_wilcox(
      v600e,
      wt,
      paired = FALSE,
      exact = TRUE,
      correct = FALSE,
      alternative = "two.sided"
    )
  } else {
    NA_real_
  }

  asymptotic_wilcoxon <- if (
    n_wt >= 2L &&
    n_v600e >= 2L
  ) {
    safe_final_wilcox(
      v600e,
      wt,
      paired = FALSE,
      exact = FALSE,
      correct = FALSE,
      alternative = "two.sided"
    )
  } else {
    NA_real_
  }

  primary_p <- if (is.null(welch_t)) {
    NA_real_
  } else {
    welch_t$p.value
  }

  statistics <- data.frame(
    Dataset = analysis_label,
    Signature = signature_name,
    Experimental_unit = "Organoid line",
    Primary_test =
      "Two-sided unpaired Welch t-test (V600E versus WT)",
    N_WT = n_wt,
    N_V600E = n_v600e,
    Mean_WT = ifelse(n_wt > 0L, mean(wt), NA_real_),
    Mean_V600E = ifelse(
      n_v600e > 0L,
      mean(v600e),
      NA_real_
    ),
    Mean_difference_V600E_minus_WT = ifelse(
      n_wt > 0L && n_v600e > 0L,
      mean(v600e) - mean(wt),
      NA_real_
    ),
    Median_difference_V600E_minus_WT = ifelse(
      n_wt > 0L && n_v600e > 0L,
      stats::median(v600e) - stats::median(wt),
      NA_real_
    ),
    P_value = primary_p,
    Welch_t_95CI_lower = if (is.null(welch_t)) {
      NA_real_
    } else {
      unname(welch_t$conf.int[1])
    },
    Welch_t_95CI_upper = if (is.null(welch_t)) {
      NA_real_
    } else {
      unname(welch_t$conf.int[2])
    },
    Welch_degrees_of_freedom = if (is.null(welch_t)) {
      NA_real_
    } else {
      unname(welch_t$parameter)
    },
    Exact_Wilcoxon_rank_sum_P_value = exact_wilcoxon,
    Asymptotic_Wilcoxon_rank_sum_P_value =
      asymptotic_wilcoxon,
    Points_connected_in_plot = FALSE,
    stringsAsFactors = FALSE
  )

  list(
    plot_data = dat,
    statistics = statistics
  )
}

make_final_two_group_plot <- function(
    test_object,
    signature_name,
    group_col,
    group_levels,
    colours,
    x_labels,
    test_label,
    n_label,
    y_label = "Singscore"
) {
  dat <- test_object$plot_data
  stats_one <- test_object$statistics[1, , drop = FALSE]

  if (nrow(dat) == 0L) {
    return(make_empty_plot(signature_name, "No valid samples"))
  }

  dat[[group_col]] <- factor(
    dat[[group_col]],
    levels = group_levels
  )

  y_max <- max(dat$Score, na.rm = TRUE)
  y_min <- min(dat$Score, na.rm = TRUE)
  y_span <- max(y_max - y_min, 0.05)

  bracket_y <- y_max + 0.09 * y_span
  text_y <- y_max + 0.31 * y_span

  annotation_text <- if (
    is.finite(stats_one$P_value)
  ) {
    paste0(
      test_label,
      " ",
      format_p(stats_one$P_value),
      "\n",
      n_label
    )
  } else {
    paste0(
      test_label,
      ": p = NA\n",
      n_label
    )
  }

  ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = .data[[group_col]],
      y = Score,
      fill = .data[[group_col]]
    )
  ) +

    ggplot2::geom_boxplot(
      width = 0.55,
      outlier.shape = NA,
      linewidth = 0.35,
      alpha = 0.45,
      colour = "black"
    ) +
    ggplot2::geom_point(
      shape = 21,
      size = 2.0,
      stroke = 0.30,
      colour = "black",
      position = ggplot2::position_jitter(
        width = 0.055,
        height = 0,
        seed = 1
      ),
      alpha = 0.95
    ) +
    ggplot2::annotate(
      "segment",
      x = 1,
      xend = 2,
      y = bracket_y,
      yend = bracket_y,
      linewidth = 0.30,
      colour = "black"
    ) +
    ggplot2::annotate(
      "text",
      x = 1.5,
      y = text_y,
      label = annotation_text,
      family = font_family,
      size = 5.2 / 2.845276,
      lineheight = 0.94,
      vjust = 0.5
    ) +
    ggplot2::scale_fill_manual(
      values = colours,
      breaks = group_levels,
      drop = FALSE
    ) +
    ggplot2::scale_x_discrete(
      labels = stats::setNames(
        x_labels,
        group_levels
      ),
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(0.05, 0.45)
      )
    ) +
    ggplot2::labs(
      title = format_signature_plot_title(signature_name),
      x = NULL,
      y = y_label,
      fill = NULL
    ) +
    nature_theme() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(size = 6)
    )
}

final_infection_tests <- purrr::map(
  names(signature_list),
  function(signature_name) {
    prepare_final_infection_test(
      infection_scores$scores,
      signature_name
    )
  }
)
names(final_infection_tests) <- names(signature_list)

final_infection_plots <- purrr::imap(
  final_infection_tests,
  function(test_object, signature_name) {
    n_events <- test_object$statistics$
      N_independent_infection_events[1]

    p <- make_final_two_group_plot(
      test_object = test_object,
      signature_name = signature_name,
      group_col = "Infection_Status",
      group_levels = c("Non-Infected", "Infected"),
      colours = infection_cols,
      x_labels = c("Non-Infected", "Infected"),
      test_label = "paired t-test",
      n_label = paste0(
        "n = ",
        n_events,
        " independent infection events"
      )
    )

    save_plot_set(
      p,
      subdirs$scores,
      paste0(
        "FINAL_Infection_",
        signature_file_stems[[signature_name]],
        "_singscore_NO_CONNECTING_LINES"
      ),
      72,
      82
    )

    save_plot_set(
      p,
      subdirs$scores,
      paste0(
        "Infection_",
        signature_file_stems[[signature_name]],
        "_singscore"
      ),
      72,
      82
    )

    p
  }
)

final_infection_statistics <- dplyr::bind_rows(
  lapply(final_infection_tests, `[[`, "statistics")
) %>%
  dplyr::mutate(
    P_value_BH_across_four_prespecified_scores =
      stats::p.adjust(P_value, method = "BH")
  )

final_infection_pair_differences <- dplyr::bind_rows(
  purrr::imap(
    final_infection_tests,
    function(test_object, signature_name) {
      test_object$paired_wide %>%
        dplyr::mutate(
          Signature = signature_name,
          .before = 1
        )
    }
  )
)

data.table::fwrite(
  final_infection_statistics,
  file.path(
    subdirs$tables,
    "Infection_signature_score_statistics.csv"
  )
)

data.table::fwrite(
  final_infection_statistics,
  file.path(
    subdirs$tables,
    "FINAL_Infection_paired_t_singscore_statistics.csv"
  )
)

data.table::fwrite(
  final_infection_pair_differences,
  file.path(
    subdirs$tables,
    "FINAL_Infection_event_level_singscore_differences.csv"
  )
)

final_infection_combined <- patchwork::wrap_plots(
  final_infection_plots,
  ncol = 4
) +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = "Infection signature scores",
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
        size = panel_tag_pt
      )
    )
  )

save_plot_set(
  final_infection_combined,
  subdirs$combined,
  "Infection_signature_scores_combined",
  183,
  84
)

save_plot_set(
  final_infection_combined,
  final_nature_combined_dir,
  "FINAL_Infection_signature_scores_NO_CONNECTING_LINES_PAIRED_T",
  183,
  84
)

final_braf_tests <- purrr::map(
  names(signature_list),
  function(signature_name) {
    prepare_final_braf_test(
      braf_score_input,
      signature_name,
      "BRAF_SSL_HP_only_initial_passage"
    )
  }
)
names(final_braf_tests) <- names(signature_list)

final_braf_plots <- purrr::imap(
  final_braf_tests,
  function(test_object, signature_name) {
    n_wt <- test_object$statistics$N_WT[1]
    n_v600e <- test_object$statistics$N_V600E[1]

    p <- make_final_two_group_plot(
      test_object = test_object,
      signature_name = signature_name,
      group_col = "BRAF_Status",
      group_levels = c("WT", "V600E"),
      colours = braf_cols,
      x_labels = c("WT", "V600E"),
      test_label = "Welch t-test",
      n_label = paste0(
        "WT n = ",
        n_wt,
        "; V600E n = ",
        n_v600e
      )
    )

    save_plot_set(
      p,
      subdirs$scores,
      paste0(
        "FINAL_BRAF_SSL_HP_only_",
        signature_file_stems[[signature_name]],
        "_singscore_WELCH_T"
      ),
      68,
      78
    )

    save_plot_set(
      p,
      subdirs$scores,
      paste0(
        "BRAF_SSL_HP_only_",
        signature_file_stems[[signature_name]],
        "_singscore"
      ),
      68,
      78
    )

    p
  }
)

final_braf_statistics <- dplyr::bind_rows(
  lapply(final_braf_tests, `[[`, "statistics")
) %>%
  dplyr::mutate(
    P_value_BH_across_four_prespecified_scores =
      stats::p.adjust(P_value, method = "BH")
  )

data.table::fwrite(
  final_braf_statistics,
  file.path(
    subdirs$tables,
    "BRAF_SSL_HP_only_signature_score_statistics.csv"
  )
)

data.table::fwrite(
  final_braf_statistics,
  file.path(
    subdirs$tables,
    "FINAL_BRAF_SSL_HP_only_Welch_t_singscore_statistics.csv"
  )
)

final_braf_combined <- patchwork::wrap_plots(
  final_braf_plots,
  ncol = 4
) +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = "BRAF WT versus V600E signature scores",
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
        size = panel_tag_pt
      )
    )
  )

save_plot_set(
  final_braf_combined,
  subdirs$combined,
  "BRAF_SSL_HP_only_signature_scores_combined",
  183,
  80
)

save_plot_set(
  final_braf_combined,
  final_nature_combined_dir,
  "FINAL_BRAF_signature_scores_UNPAIRED_WELCH_T",
  183,
  80
)

final_braf_p1_tests <- purrr::map(
  names(signature_list),
  function(signature_name) {
    prepare_final_braf_test(
      braf_p1_score_input,
      signature_name,
      "BRAF_SSL_HP_only_strict_P1"
    )
  }
)
names(final_braf_p1_tests) <- names(signature_list)

final_braf_p1_plots <- purrr::imap(
  final_braf_p1_tests,
  function(test_object, signature_name) {
    n_wt <- test_object$statistics$N_WT[1]
    n_v600e <- test_object$statistics$N_V600E[1]

    p <- make_final_two_group_plot(
      test_object = test_object,
      signature_name = signature_name,
      group_col = "BRAF_Status",
      group_levels = c("WT", "V600E"),
      colours = braf_cols,
      x_labels = c("WT", "V600E"),
      test_label = "Welch t-test",
      n_label = paste0(
        "WT n = ",
        n_wt,
        "; V600E n = ",
        n_v600e,
        " (strict P1)"
      )
    ) +
      ggplot2::labs(subtitle = "Strict P1 only") +
      ggplot2::theme(
        plot.subtitle = ggplot2::element_text(
          family = font_family,
          size = 5.5,
          hjust = 0
        )
      )

    save_plot_set(
      p,
      subdirs$scores,
      paste0(
        "BRAF_Strict_P1_SSL_HP_only_",
        signature_file_stems[[signature_name]],
        "_singscore"
      ),
      68,
      80
    )

    p
  }
)

final_braf_p1_statistics <- dplyr::bind_rows(
  lapply(final_braf_p1_tests, `[[`, "statistics")
) %>%
  dplyr::mutate(
    P_value_BH_across_four_prespecified_scores =
      stats::p.adjust(P_value, method = "BH")
  )

data.table::fwrite(
  final_braf_p1_statistics,
  file.path(
    subdirs$tables,
    "BRAF_Strict_P1_SSL_HP_only_signature_score_statistics.csv"
  )
)

data.table::fwrite(
  final_braf_p1_statistics,
  file.path(
    subdirs$tables,
    "FINAL_BRAF_Strict_P1_Welch_t_singscore_statistics.csv"
  )
)

final_braf_p1_combined <- patchwork::wrap_plots(
  final_braf_p1_plots,
  ncol = 4
) +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = "BRAF WT versus V600E signature scores: strict P1",
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
        size = panel_tag_pt
      )
    )
  )

save_plot_set(
  final_braf_p1_combined,
  subdirs$combined,
  "BRAF_Strict_P1_SSL_HP_only_signature_scores_combined",
  183,
  82
)

final_braf_passage_comparison <- patchwork::wrap_plots(
  c(final_braf_plots, final_braf_p1_plots),
  ncol = 4
) +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = "BRAF score sensitivity: initial passage versus strict P1",
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
        size = panel_tag_pt
      )
    )
  )

save_plot_set(
  final_braf_passage_comparison,
  subdirs$combined,
  "BRAF_Initial_P1_or_P2_vs_Strict_P1_signature_score_comparison",
  183,
  154
)

final_braf_patient_scores <- braf_score_input %>%
  dplyr::mutate(
    Patient_ID = derive_patient_id_audited(Line_ID)
  ) %>%
  dplyr::group_by(
    Signature,
    Patient_ID,
    BRAF_Status
  ) %>%
  dplyr::summarise(
    Score = mean(Score, na.rm = TRUE),
    N_organoid_lines = dplyr::n_distinct(Line_ID),
    Lines = paste(
      sort(unique(Line_ID)),
      collapse = ";"
    ),
    Passages = paste(
      sort(unique(Passage)),
      collapse = ";"
    ),
    .groups = "drop"
  )

final_braf_patient_tests <- purrr::map(
  names(signature_list),
  function(signature_name) {
    prepare_final_braf_test(
      final_braf_patient_scores,
      signature_name,
      "BRAF_SSL_HP_patient_collapsed"
    )
  }
)
names(final_braf_patient_tests) <- names(signature_list)

final_braf_patient_statistics <- dplyr::bind_rows(
  lapply(final_braf_patient_tests, `[[`, "statistics")
) %>%
  dplyr::mutate(
    Experimental_unit = "Patient-collapsed mean score",
    P_value_BH_across_four_prespecified_scores =
      stats::p.adjust(P_value, method = "BH")
  )

data.table::fwrite(
  final_braf_patient_scores,
  file.path(
    subdirs$tables,
    "FINAL_BRAF_patient_collapsed_singscores.csv"
  )
)

data.table::fwrite(
  final_braf_patient_statistics,
  file.path(
    subdirs$tables,
    "FINAL_BRAF_patient_collapsed_Welch_t_statistics.csv"
  )
)

if (exists("muc5ac_infection", inherits = FALSE) &&
    !is.null(muc5ac_infection)) {
  muc5ac_infection_final_input <- muc5ac_infection %>%
    dplyr::mutate(Signature = "MUC5AC")

  muc5ac_infection_final_test <- prepare_final_infection_test(
    muc5ac_infection_final_input,
    "MUC5AC"
  )

  p_muc5ac_infection_final <- make_final_two_group_plot(
    test_object = muc5ac_infection_final_test,
    signature_name = "MUC5AC",
    group_col = "Infection_Status",
    group_levels = c("Non-Infected", "Infected"),
    colours = infection_cols,
    x_labels = c("Non-Infected", "Infected"),
    test_label = "paired t-test",
    n_label = paste0(
      "n = ",
      muc5ac_infection_final_test$statistics$
        N_independent_infection_events[1],
      " independent infection events"
    ),
    y_label = "MUC5AC log2 expression"
  )

  save_plot_set(
    p_muc5ac_infection_final,
    subdirs$scores,
    "Infection_MUC5AC_expression",
    72,
    82
  )
}

if (exists("muc5ac_braf", inherits = FALSE) &&
    !is.null(muc5ac_braf)) {
  muc5ac_braf_final_input <- muc5ac_braf %>%
    dplyr::filter(Type == "SSL_HP") %>%
    dplyr::mutate(Signature = "MUC5AC")

  muc5ac_braf_final_test <- prepare_final_braf_test(
    muc5ac_braf_final_input,
    "MUC5AC",
    "BRAF_SSL_HP_only_initial_passage"
  )

  p_muc5ac_braf_final <- make_final_two_group_plot(
    test_object = muc5ac_braf_final_test,
    signature_name = "MUC5AC",
    group_col = "BRAF_Status",
    group_levels = c("WT", "V600E"),
    colours = braf_cols,
    x_labels = c("WT", "V600E"),
    test_label = "Welch t-test",
    n_label = paste0(
      "WT n = ",
      muc5ac_braf_final_test$statistics$N_WT[1],
      "; V600E n = ",
      muc5ac_braf_final_test$statistics$N_V600E[1]
    ),
    y_label = "MUC5AC log2 expression"
  )

  save_plot_set(
    p_muc5ac_braf_final,
    subdirs$scores,
    "BRAF_SSL_HP_only_MUC5AC_expression",
    68,
    78
  )
}

final_singscore_method_table <- data.frame(
  Experiment = c(
    "Infection",
    "BRAF initial passage",
    "BRAF strict P1",
    "BRAF patient-collapsed sensitivity",
    "Passage time course"
  ),
  Groups = c(
    "Non-Infected versus Infected",
    "WT versus V600E",
    "WT versus V600E",
    "WT versus V600E",
    "Repeated passage measurements"
  ),
  Statistical_design = c(
    "Matched within each separately conducted passage-specific infection event",
    "Unpaired organoid lines",
    "Unpaired organoid lines",
    "Unpaired patients after averaging multiple lines within patient",
    "Repeated measures within trajectory"
  ),
  Primary_test = c(
    "Two-sided paired t-test",
    "Two-sided unpaired Welch t-test",
    "Two-sided unpaired Welch t-test",
    "Two-sided unpaired Welch t-test sensitivity",
    "Random-intercept mixed-effects model from the main pipeline"
  ),
  Plot_connecting_lines = c(
    "No",
    "No",
    "No",
    "No",
    "Time-course trajectories retain biologically meaningful passage lines"
  ),
  Sensitivity_test = c(
    "Exact and asymptotic paired Wilcoxon",
    "Exact/asymptotic Wilcoxon rank-sum",
    "Exact/asymptotic Wilcoxon rank-sum",
    "Exact/asymptotic Wilcoxon rank-sum",
    "Descriptive Spearman table"
  ),
  stringsAsFactors = FALSE
)

data.table::fwrite(
  final_singscore_method_table,
  file.path(
    subdirs$qc,
    "FINAL_SINGSCORE_STATISTICAL_METHODS.csv"
  )
)

writeLines(
  c(
    "FINAL NATURE-EXACT FIGURE AND STATISTICS NOTES",
    "===============================================",
    "",
    "Singscore plots:",
    "  Final two-group plots contain no point-connecting lines.",
    "  Pairing is nevertheless retained in the infection paired t-test.",
    "  Infection uses passage-specific Pair_ID as the independent infection event.",
    "  Technical replicates within one Pair_ID/condition are averaged before testing.",
    "  BRAF WT versus V600E is not paired; final BRAF plots use a two-sided Welch t-test.",
    "  Exact/asymptotic Wilcoxon results are retained as sensitivity columns.",
    "  BRAF initial and strict-P1 analyses both use the same unpaired Welch framework.",
    "  A patient-collapsed Welch sensitivity is also exported.",
    "",
    "Four-module log2FC heatmap:",
    "  Exact visual construction copied from the final Nature-style tissue script.",
    "  Modules: interferon alpha response, chemokines, antimicrobial peptides and TNFA/NFKB signaling.",
    paste0(
      "  Reference gene selection: infection logCPM >= ",
      module_heatmap_logcpm_cutoff,
      " and logFC > ",
      module_heatmap_logfc_cutoff,
      "."
    ),
    paste0(
      "  Up to ",
      final_heatmap_max_genes_per_module,
      " genes per module, ranked by decreasing infection logFC."
    ),
    "  The same selected genes are displayed unchanged for infection and initial BRAF.",
    "  Four modules are arranged in one horizontal row.",
    paste0("  Font: ", font_family, "."),
    paste0("  Gene/comparison labels: ", final_heatmap_shared_label_pt, " pt."),
    paste0("  Comparison-label angle: ", final_heatmap_x_label_angle, " degrees."),
    paste0(
      "  Heatmap colours: ",
      final_heatmap_negative_colour,
      " / white / ",
      final_heatmap_positive_colour,
      "."
    ),
    paste0(
      "  log2FC colour saturation: -",
      final_heatmap_logfc_limit,
      " to +",
      final_heatmap_logfc_limit,
      "."
    ),
    "  White tile borders and a right-side centred vertical legend are used.",
    "",
    "Hallmark NES heatmaps:",
    "  Rewritten with the same Arial, 5-point labels, 55-degree x labels,",
    "  white tile borders, blue-white-red scale and centred vertical legend.",
    "",
    "The immediately preceding renderings were copied to LEGACY_STYLE folders",
    "before the standard output stems were overwritten with these final files."
  ),
  file.path(
    subdirs$qc,
    "FINAL_NATURE_EXACT_AND_SINGSCORE_NOTES.txt"
  )
)

cat(
  paste0(
    "\n\nFINAL NATURE-EXACT OVERRIDE\n",
    "  Final two-group singscore plots have no connecting lines.\n",
    "  Infection primary test: two-sided paired t-test across independent Pair_ID events.\n",
    "  BRAF primary test: two-sided unpaired Welch t-test for WT versus V600E.\n",
    "  Four-module log2FC heatmaps use the exact Nature-style one-row construction.\n",
    "  Standard Hallmark NES heatmap stems now use the matching Nature style.\n",
    "  Previous renderings were copied to LEGACY_STYLE folders.\n"
  ),
  file = file.path(out_dir, "RUN_NOTES.txt"),
  append = TRUE
)

message(
  "Final override complete. Inspect:\n",
  "  05_Heatmaps/FINAL_Nature_exact/FINAL_ALL4_Module_log2FC_Heatmaps_Infection_and_BRAF_NATURE_EXACT.*\n",
  "  06_Combined_figures/Infection_signature_scores_combined.*\n",
  "  06_Combined_figures/BRAF_SSL_HP_only_signature_scores_combined.*\n",
  "  07_Tables/Infection_signature_score_statistics.csv\n",
  "  07_Tables/BRAF_SSL_HP_only_signature_score_statistics.csv\n",
  "  01_QC/FINAL_NATURE_EXACT_AND_SINGSCORE_NOTES.txt"
)

message("\nWriting pooled comparison-level log2FC heatmaps for ", metaplasia_display_name, "...")

# Gastric-signature log2 fold-change heatmaps ----
prepare_pooled_metaplasia_logfc <- function(
    de,
    comparison_key,
    comparison_label
) {
  signature_order <- unique(toupper(metaplasia_genes))

  tibble::tibble(
    Gene = signature_order,
    Signature_order = seq_along(signature_order)
  ) %>%
    dplyr::left_join(
      de %>%
        dplyr::select(
          Gene,
          logFC,
          logCPM,
          PValue,
          FDR
        ),
      by = "Gene"
    ) %>%
    dplyr::mutate(
      Comparison_key = comparison_key,
      Comparison_label = comparison_label,
      Present_in_DE_table = is.finite(logFC),
      logFC_colour = pmax(
        pmin(logFC, final_heatmap_logfc_limit),
        -final_heatmap_logfc_limit
      ),
      Significance = final_significance_stars(FDR)
    ) %>%
    dplyr::arrange(Signature_order)
}

pooled_metaplasia_logfc_values <- dplyr::bind_rows(
  prepare_pooled_metaplasia_logfc(
    infection_de,
    comparison_key = "Infection",
    comparison_label = "Infected vs Non-Infected"
  ),
  prepare_pooled_metaplasia_logfc(
    braf_de,
    comparison_key = "BRAF_initial",
    comparison_label = "BRAF V600E vs WT in SSL/HP"
  )
)

make_pooled_metaplasia_logfc_heatmap <- function(
    pooled_data,
    comparison_keys,
    subtitle_text
) {
  signature_order <- unique(toupper(metaplasia_genes))

  plot_df <- pooled_data %>%
    dplyr::filter(Comparison_key %in% comparison_keys) %>%
    dplyr::mutate(
      Comparison_label = factor(
        Comparison_label,
        levels = unique(Comparison_label[match(comparison_keys, Comparison_key)])
      ),
      Gene = factor(
        Gene,
        levels = rev(signature_order)
      ),
      Star_colour = dplyr::if_else(
        is.finite(logFC_colour) & abs(logFC_colour) >= 0.85,
        "white",
        "black"
      )
    )

  if (nrow(plot_df) == 0L) {
    return(make_empty_plot(
      metaplasia_display_name,
      "No pooled differential-expression values were available"
    ))
  }

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = Comparison_label,
      y = Gene,
      fill = logFC_colour
    )
  ) +
    ggplot2::geom_tile(
      colour = final_heatmap_tile_border_colour,
      linewidth = final_heatmap_tile_border_linewidth,
      width = final_heatmap_tile_fraction,
      height = final_heatmap_tile_fraction
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = Significance,
        colour = Star_colour
      ),
      family = font_family,
      fontface = "bold",
      size = final_geom_text_size_from_pt(
        final_heatmap_shared_label_pt
      ),
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_x_discrete(
      drop = FALSE,
      expand = ggplot2::expansion(add = 0.05)
    ) +
    ggplot2::scale_y_discrete(
      drop = FALSE,
      expand = ggplot2::expansion(add = 0.08)
    ) +
    ggplot2::scale_fill_gradient2(
      low = final_heatmap_negative_colour,
      mid = final_heatmap_mid_colour,
      high = final_heatmap_positive_colour,
      midpoint = 0,
      limits = c(
        -final_heatmap_logfc_limit,
        final_heatmap_logfc_limit
      ),
      breaks = c(
        -final_heatmap_logfc_limit,
        0,
        final_heatmap_logfc_limit
      ),
      labels = c(
        paste0("≤ -", final_heatmap_logfc_limit),
        "0",
        paste0("≥ +", final_heatmap_logfc_limit)
      ),
      oob = scales::squish,
      na.value = "grey90",
      name = "log2FC",
      guide = final_vertical_heatmap_guide(
        barheight_mm = 18,
        barwidth_mm = 2.0
      )
    ) +
    ggplot2::coord_fixed(ratio = 1, clip = "off") +
    ggplot2::labs(
      title = format_signature_plot_title(
        metaplasia_display_name,
        width = 46
      ),
      subtitle = subtitle_text,
      x = NULL,
      y = NULL
    ) +
    final_nature_heatmap_theme(
      base_size = final_heatmap_font_base_size
    ) +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(
        family = font_family,
        size = 5.5,
        hjust = 0,
        colour = "black"
      ),
      axis.text.x = ggplot2::element_text(
        family = font_family,
        size = final_heatmap_shared_label_pt,
        angle = final_heatmap_x_label_angle,
        hjust = 1,
        vjust = 1,
        colour = "black"
      ),
      axis.text.y = ggplot2::element_text(
        family = font_family,
        size = final_heatmap_shared_label_pt,
        colour = "black"
      ),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      legend.direction = "vertical",
      plot.margin = ggplot2::margin(3, 3, 3, 3, unit = "mm")
    )
}

p_metaplasia_logfc_infection <- make_pooled_metaplasia_logfc_heatmap(
  pooled_metaplasia_logfc_values,
  comparison_keys = "Infection",
  subtitle_text = "Pooled differential-expression log2FC: Infected vs Non-Infected"
)

p_metaplasia_logfc_braf <- make_pooled_metaplasia_logfc_heatmap(
  pooled_metaplasia_logfc_values,
  comparison_keys = "BRAF_initial",
  subtitle_text = "Pooled differential-expression log2FC: BRAF V600E vs WT in SSL/HP"
)

p_metaplasia_logfc_combined <- make_pooled_metaplasia_logfc_heatmap(
  pooled_metaplasia_logfc_values,
  comparison_keys = c("Infection", "BRAF_initial"),
  subtitle_text = "Pooled differential-expression log2FC"
)

save_plot_set(
  p_metaplasia_logfc_infection,
  final_nature_heatmap_dir,
  paste0(
    "FINAL_Infection_",
    metaplasia_file_stem,
    "_pooled_log2FC_heatmap"
  ),
  88,
  150
)

save_plot_set(
  p_metaplasia_logfc_braf,
  final_nature_heatmap_dir,
  paste0(
    "FINAL_BRAF_",
    metaplasia_file_stem,
    "_pooled_log2FC_heatmap"
  ),
  88,
  150
)

save_plot_set(
  p_metaplasia_logfc_combined,
  final_nature_heatmap_dir,
  paste0(
    "FINAL_Infection_and_BRAF_",
    metaplasia_file_stem,
    "_pooled_log2FC_heatmap"
  ),
  112,
  150
)

save_plot_set(
  p_metaplasia_logfc_infection,
  subdirs$heatmaps,
  paste0(
    "Infection_DE_logFC_heatmap_",
    metaplasia_file_stem
  ),
  88,
  150
)

save_plot_set(
  p_metaplasia_logfc_braf,
  subdirs$heatmaps,
  paste0(
    "BRAF_SSL_HP_only_DE_logFC_heatmap_",
    metaplasia_file_stem
  ),
  88,
  150
)

save_plot_set(
  p_metaplasia_logfc_combined,
  subdirs$heatmaps,
  paste0(
    "Infection_and_BRAF_DE_logFC_heatmap_",
    metaplasia_file_stem
  ),
  112,
  150
)

data.table::fwrite(
  pooled_metaplasia_logfc_values,
  file.path(
    subdirs$tables,
    paste0(
      metaplasia_file_stem,
      "_pooled_DE_log2FC_heatmap_values.csv"
    )
  )
)

writeLines(
  c(
    "POOLED METAPLASIA LOG2FC HEATMAPS",
    "=================================",
    "",
    paste0("Signature display name: ", metaplasia_display_name),
    "The infection heatmap uses the pooled edgeR log2FC for Infected versus Non-Infected.",
    "The BRAF heatmap uses the pooled edgeR log2FC for V600E versus WT in SSL/HP.",
    "All requested signature genes are displayed in signature-file order.",
    "No expression, log2FC or significance threshold is applied to gene inclusion.",
    "Missing genes are grey; stars indicate FDR < 0.05, 0.01 or 0.001.",
    "The original sample-level row-z-score heatmaps are retained under explicit sample_level_row_zscore filenames."
  ),
  file.path(
    subdirs$qc,
    paste0(
      metaplasia_file_stem,
      "_pooled_log2FC_heatmap_methods.txt"
    )
  )
)

cat(
  paste0(
    "\n\nPOOLED METAPLASIA LOG2FC HEATMAPS ADDED\n",
    "  Plot name: ", metaplasia_display_name, ".\n",
    "  Infection and BRAF comparison-level log2FC heatmaps were written.\n",
    "  Sample-level row-z-score heatmaps were retained separately.\n"
  ),
  file = file.path(out_dir, "RUN_NOTES.txt"),
  append = TRUE
)

message(
  "Pooled metaplasia log2FC heatmaps complete. Main files:\n",
  "  05_Heatmaps/Infection_DE_logFC_heatmap_", metaplasia_file_stem, ".*\n",
  "  05_Heatmaps/BRAF_SSL_HP_only_DE_logFC_heatmap_", metaplasia_file_stem, ".*\n",
  "  05_Heatmaps/Infection_and_BRAF_DE_logFC_heatmap_", metaplasia_file_stem, ".*\n",
  "  07_Tables/", metaplasia_file_stem, "_pooled_DE_log2FC_heatmap_values.csv"
)
