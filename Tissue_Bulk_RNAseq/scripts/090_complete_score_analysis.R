#!/usr/bin/env Rscript
# Tissue scores, patient-aware comparisons and correlation figures.
# Run the complete file in a fresh R session; inputs are configured below.
# Default: CPM -> TPM/scores -> comparison with saved scores -> tests and figures.
# OSLC_SCORE_INPUT_MODE=saved_scores explicitly reuses the verified private scores.
# Spearman correlations use patient-block bootstrap inference; mixed-model
# association tests are exported separately. Original outputs are preserved.
# Scope: Figures 3.5.2, 3.9B, 3.10, SF3.10 and SF4.6.3, plus sensitivities.

BASE_DIR <- Sys.getenv("OSLC_BASE_DIR", Sys.getenv("OSLC_BULK_DIR"))
OUTPUT_PARENT <- Sys.getenv("OSLC_RERUN_OUTPUT_PARENT", file.path(getwd(), "outputs"))
REFERENCE_DIR <- Sys.getenv("OSLC_SCORE_REFERENCE_DIR", file.path(getwd(), "private_inputs", "score"))
read_reference_text <- function(name) {
  path <- file.path(REFERENCE_DIR, name)
  if (!file.exists(path)) stop("Missing private reference input: ", path,
    ". Set OSLC_SCORE_REFERENCE_DIR to the supplied private score input folder.")
  paste(readLines(path, warn=FALSE), collapse="\n")
}

# Input paths and run options.
SCORE_INPUT_MODE <- Sys.getenv('OSLC_SCORE_INPUT_MODE', 'from_cpm')
CPM_FILE <- Sys.getenv('OSLC_CPM_FILE',file.path(BASE_DIR,'Re-run GSEA','CGG014_Maier_main_cpm_table.txt'))
GENE_LENGTH_FILE <- Sys.getenv('OSLC_GENE_LENGTH_FILE',file.path(BASE_DIR,'GM_publication_ready','GRCh38.113.gene.lengths'))
METADATA_FILE <- Sys.getenv('OSLC_METADATA_FILE',file.path(BASE_DIR,'Metadata_March26.csv'))
BOOTSTRAP_REPLICATES <- as.integer(Sys.getenv('OSLC_BOOTSTRAP_REPLICATES','9999'))
if(is.na(BOOTSTRAP_REPLICATES)||BOOTSTRAP_REPLICATES<199L) stop('At least 199 bootstrap replicates are required; use 9999 for the thesis.')

INSTALL_MISSING <- identical(Sys.getenv("OSLC_INSTALL_MISSING"), "1")
MAKE_PLOTS <- !identical(Sys.getenv('OSLC_RERUN_MAKE_PLOTS'), '0')
FIGURE_FORMATS <- strsplit(Sys.getenv('OSLC_RERUN_FORMATS', 'pdf,png,jpg,jpeg,tiff'), ',')[[1]]
FIGURE_DPI <- as.integer(Sys.getenv('OSLC_RERUN_DPI', '600'))
SET_SEED <- 1L
THIS_SCRIPT <- local({
  args <- commandArgs(trailingOnly = FALSE)
  command_file <- sub('^--file=', '', args[grepl('^--file=', args)])
  source_files <- unlist(lapply(sys.frames(), function(f) if (!is.null(f$ofile)) as.character(f$ofile) else NULL))
  candidates <- unique(c(command_file, rev(source_files)))
  candidates <- candidates[file.exists(candidates)]
  if(length(candidates)) normalizePath(candidates[1],winslash='/',mustWork=TRUE) else NA_character_
})

params <- list(
  font_family = 'Arial', font_pt = 7, x_font_pt = 6, label_pt = 6,
  figure_width_mm = 183, figure_dpi = FIGURE_DPI,
  score_plot_groups = c('SSL_HP','H_ADJ_SSL_HP','H_COL','CA'),
  group_order = c('SSL_HP','H_ADJ_SSL_HP','H_COL','CA'),
  group_labels = c(SSL_HP='SSL_HP',H_ADJ_SSL_HP='H_ADJ_SSL_HP',H_COL='H_COL',CA='CA'),
  group_colours = c(SSL_HP='#FF0000',H_ADJ_SSL_HP='#E9C61D',H_COL='#1313EC',CA='#825CA6'),
  signature_display_names = c(GASTRIC_METAPLASIA_24_GENE='Gastric metaplasia',
    INTERFERON_ALPHA_RESPONSE='Interferon alpha response',
    INFLAMMATORY_RESPONSE='Inflammatory response',TNFA_SIGNALING_VIA_NFKB='TNFA/NFKB signalling'),
  method_axis_labels = c(Zscore='Mean gene-wise z-score',Singscore='Centred singscore',SSGSEA='ssGSEA score'))
SIGNATURES <- names(params$signature_display_names)
METHODS <- c('Singscore','Zscore','SSGSEA')
SUBSETS <- list(SSL_HP_ONLY='SSL_HP',ALL_COLORECTAL=params$score_plot_groups,SSL_HP_CA=c('SSL_HP','CA'))
pt_to_mm <- function(x) x / 2.845276


# Private reference inputs from the exact saved analysis; exported into each local run.
score_csv <- read_reference_text("score_csv.csv")

original_singscore_csv <- read_reference_text("original_singscore_csv.csv")

score_coverage_csv <- read_reference_text("score_coverage_csv.csv")

metadata_csv <- read_reference_text("metadata_csv.csv")

master_csv <- read_reference_text("master_csv.csv")

ifna_filter_csv <- read_reference_text("ifna_filter_csv.csv")

inflam_filter_csv <- read_reference_text("inflam_filter_csv.csv")

tnf_filter_csv <- read_reference_text("tnf_filter_csv.csv")

gastric_genes_txt <- paste(c(
  "AQP5",
  "ANXA10",
  "MUC5AC",
  "VSIG1",
  "TCN1",
  "MUC6",
  "PSAPL1",
  "CLDN18",
  "SFTA2",
  "MUCL3",
  "GATA4",
  "KLK11",
  "DUOXA1",
  "SULT1C2",
  "IRX3",
  "AADAC",
  "FOXQ1",
  "IGFALS",
  "NR0B2",
  "C13ORF46",
  "FOLR1",
  "LYPD6B",
  "TFR2",
  "SPTSSB"
), collapse = "\n")

reference_alpha_csv <- read_reference_text("reference_alpha_csv.csv")

reference_gastric_csv <- read_reference_text("reference_gastric_csv.csv")


# Plot functions retained from the original v11 script.
paper_theme <- function(base_size = params$font_pt) {
  theme_classic(base_family = params$font_family, base_size = base_size) +
    theme(
      text = element_text(family = params$font_family, colour = "black"),
      axis.title = element_text(size = params$font_pt, colour = "black"),
      axis.text = element_text(size = params$font_pt, colour = "black"),
      axis.text.x = element_text(size = params$x_font_pt, colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = unit(1.2, "mm"),
      plot.title = element_text(
        size = params$font_pt, face = "bold", hjust = 0, colour = "black"
      ),
      legend.title = element_text(size = params$font_pt, colour = "black"),
      legend.text = element_text(size = params$font_pt, colour = "black"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(4, 4, 4, 4, "mm")
    )
}

score_method_y_label <- function(method) {
  label <- params$method_axis_labels[[method]]
  if (is.null(label) || is.na(label)) method else label
}

signature_pretty_label <- function(signature) {
  label <- params$signature_display_names[[signature]]
  if (is.null(label) || is.na(label)) signature else label
}

pair_key <- function(a, b) {
  paste(sort(c(as.character(a), as.character(b))), collapse = "__")
}

p_label <- function(p) {
  if (!is.finite(p)) return("p = NA")
  if (p < 0.0001) return("****")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  paste0("p = ", formatC(p, format = "fg", digits = 3))
}

make_generic_score_distribution_plot <- function(df, stats_obj, signature, method) {
  plot_df <- df %>%
    dplyr::filter(
      Group_final %in% params$score_plot_groups,
      is.finite(score)
    ) %>%
    dplyr::mutate(
      Group_plot = factor(Group_final, levels = params$score_plot_groups)
    )

  p <- ggplot(plot_df, aes(x = Group_plot, y = score)) +
    geom_violin(aes(fill = Group_plot), trim = FALSE, alpha = 0.25, colour = NA) +
    geom_boxplot(
      aes(fill = Group_plot), width = 0.16, outlier.shape = NA,
      alpha = 0.85, colour = "black", linewidth = 0.35
    ) +
    geom_jitter(
      aes(colour = Group_plot), position = position_jitter(width = 0.10, height = 0, seed = 1), size = 1.2, alpha = 0.85
    ) +
    scale_fill_manual(values = params$group_colours, drop = FALSE) +
    scale_colour_manual(values = params$group_colours, drop = FALSE) +
    scale_x_discrete(
      breaks = params$score_plot_groups,
      labels = params$group_labels[params$score_plot_groups],
      drop = FALSE
    ) +
    labs(
      title = signature_pretty_label(signature),
      x = NULL,
      y = score_method_y_label(method)
    ) +
    paper_theme() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        family = params$font_family, size = params$x_font_pt,
        angle = 25, hjust = 1, vjust = 1, colour = "black"
      )
    )

  pairwise <- stats_obj$pairwise
  if (is.data.frame(pairwise) && nrow(pairwise) > 0) {
    y_max <- max(plot_df$score, na.rm = TRUE)
    y_min <- min(plot_df$score, na.rm = TRUE)
    y_span <- y_max - y_min
    if (!is.finite(y_span) || y_span == 0) y_span <- 0.1

    group_positions <- stats::setNames(
      seq_along(params$score_plot_groups), params$score_plot_groups
    )
    # Use an absolute minimum row spacing as well as a proportional spacing.
    # This prevents overlap when the score range is narrow (for example, the
    # inflammatory-response singscore panel). Exact p-value labels receive
    # extra clearance before the next bracket; asterisks remain close to their
    # own bracket.
    base_row_gap <- max(0.16 * y_span, 0.012)
    p_row_extra <- max(0.08 * y_span, 0.008)
    tick_depth <- max(0.028 * y_span, 0.002)
    p_label_gap <- max(0.055 * y_span, 0.005)
    star_label_gap <- max(0.014 * y_span, 0.002)

    pairwise <- pairwise %>%
      dplyr::mutate(
        x1 = unname(group_positions[group1]),
        x2 = unname(group_positions[group2]),
        x_left = pmin(x1, x2),
        x_right = pmax(x1, x2),
        span = x_right - x_left
      ) %>%
      dplyr::arrange(span, x_left, x_right) %>%
      dplyr::mutate(
        is_exact_p = grepl("^p[[:space:]]*=", label),
        row_gap = base_row_gap + dplyr::lag(
          dplyr::if_else(is_exact_p, p_row_extra, 0),
          default = 0
        ),
        y = y_max + cumsum(row_gap),
        y_tick = y - tick_depth,
        y_text = y + dplyr::if_else(is_exact_p, p_label_gap, star_label_gap)
      )

    p <- p +
      geom_segment(
        data = pairwise,
        aes(x = x_left, xend = x_right, y = y, yend = y),
        inherit.aes = FALSE, linewidth = 0.35, colour = "black"
      ) +
      geom_segment(
        data = pairwise,
        aes(x = x_left, xend = x_left, y = y, yend = y_tick),
        inherit.aes = FALSE, linewidth = 0.35, colour = "black"
      ) +
      geom_segment(
        data = pairwise,
        aes(x = x_right, xend = x_right, y = y, yend = y_tick),
        inherit.aes = FALSE, linewidth = 0.35, colour = "black"
      ) +
      geom_text(
        data = pairwise,
        aes(x = (x_left + x_right) / 2, y = y_text, label = label),
        inherit.aes = FALSE, family = params$font_family,
        size = pt_to_mm(params$label_pt), colour = "black",
        vjust = 0
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  }

  p
}

make_score_distribution_plot <- function(score_df, pairwise_for_plot) {
  plot_df <- score_df %>%
    dplyr::filter(
      Group_final %in% params$score_plot_groups,
      is.finite(centered_singscore)
    ) %>%
    dplyr::mutate(
      Group_plot = factor(Group_final, levels = params$score_plot_groups),
      Group_num = as.numeric(Group_plot)
    )

  group_positions <- setNames(
    seq_along(params$score_plot_groups), params$score_plot_groups
  )

  p <- ggplot(
    plot_df,
    aes(x = Group_num, y = centered_singscore, group = Group_plot)
  ) +
    geom_violin(
      aes(fill = Group_plot),
      trim = FALSE,
      alpha = 0.25,
      colour = NA
    ) +
    geom_boxplot(
      aes(fill = Group_plot),
      width = 0.16,
      outlier.shape = NA,
      alpha = 0.85,
      colour = "black",
      linewidth = 0.35
    ) +
    geom_jitter(
      aes(colour = Group_plot),
      position = position_jitter(width = 0.10, height = 0, seed = 1),
      size = 1.2,
      alpha = 0.85
    ) +
    scale_fill_manual(values = params$group_colours, drop = FALSE) +
    scale_colour_manual(values = params$group_colours, drop = FALSE) +
    scale_x_continuous(
      breaks = seq_along(params$score_plot_groups),
      labels = params$group_labels[params$score_plot_groups],
      limits = c(0.5, length(params$score_plot_groups) + 0.5)
    ) +
    labs(x = NULL, y = "Centred gastric metaplasia singscore") +
    paper_theme() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        family = params$font_family,
        size = params$x_font_pt,
        angle = 25,
        hjust = 1,
        vjust = 1,
        colour = "black"
      )
    )

  if (nrow(pairwise_for_plot) > 0) {
    y_max <- max(plot_df$centered_singscore, na.rm = TRUE)
    y_min <- min(plot_df$centered_singscore, na.rm = TRUE)
    y_span <- y_max - y_min
    if (!is.finite(y_span) || y_span == 0) y_span <- 0.1

    # Apply the same collision-safe bracket geometry used for all generic
    # score-distribution panels. The absolute minimum spacing is important for
    # signatures whose scores occupy a narrow numerical range.
    base_row_gap <- max(0.16 * y_span, 0.012)
    p_row_extra <- max(0.08 * y_span, 0.008)
    tick_depth <- max(0.028 * y_span, 0.002)
    p_label_gap <- max(0.055 * y_span, 0.005)
    star_label_gap <- max(0.014 * y_span, 0.002)

    brackets <- pairwise_for_plot %>%
      dplyr::mutate(
        x1 = unname(group_positions[group1]),
        x2 = unname(group_positions[group2]),
        x_left = pmin(x1, x2),
        x_right = pmax(x1, x2),
        is_exact_p = grepl("^p[[:space:]]*=", label),
        row_gap = base_row_gap + dplyr::lag(
          dplyr::if_else(is_exact_p, p_row_extra, 0),
          default = 0
        ),
        y = y_max + cumsum(row_gap),
        y_tick = y - tick_depth,
        y_text = y + dplyr::if_else(is_exact_p, p_label_gap, star_label_gap)
      )

    p <- p +
      geom_segment(
        data = brackets,
        aes(x = x_left, xend = x_right, y = y, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.35,
        colour = "black"
      ) +
      geom_segment(
        data = brackets,
        aes(x = x_left, xend = x_left, y = y, yend = y_tick),
        inherit.aes = FALSE,
        linewidth = 0.35,
        colour = "black"
      ) +
      geom_segment(
        data = brackets,
        aes(x = x_right, xend = x_right, y = y, yend = y_tick),
        inherit.aes = FALSE,
        linewidth = 0.35,
        colour = "black"
      ) +
      geom_text(
        data = brackets,
        aes(x = (x_left + x_right) / 2, y = y_text, label = label),
        inherit.aes = FALSE,
        family = params$font_family,
        size = pt_to_mm(params$label_pt),
        colour = "black",
        vjust = 0
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  }

  p
}

format_p_text <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  paste0("= ", formatC(p, format = "fg", digits = 3))
}


.rerun <- new.env(parent = emptyenv())
read_embedded <- function(text) read.csv(text = text, check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = c('', 'NA', 'N/A'))
write_table <- function(x, relative) {
  p <- file.path(.rerun$out, relative)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, p, row.names = FALSE, na = '')
  invisible(p)
}
write_note <- function(x, relative) {
  p <- file.path(.rerun$out, relative)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  writeLines(x, p, useBytes = TRUE)
}
norm_sample <- function(x) {
  x <- gsub('[^A-Z0-9]+', '_', toupper(trimws(as.character(x))))
  gsub('^_+|_+$', '', x)
}
pool_group <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c('SSL','HP')] <- 'SSL_HP'
  x[x %in% c('H_ADJ_SSL','H_ADJ_HP')] <- 'H_ADJ_SSL_HP'
  x[x %in% c('H_COL_D','H_COL_ND')] <- 'H_COL'
  x
}
bh_family <- function(p, n = length(p), method = 'BH') {
  ans <- rep(NA_real_, length(p)); ok <- is.finite(p)
  if (any(ok)) ans[ok] <- p.adjust(p[ok], method = method, n = n)
  ans
}
standardise <- function(x) {
  if (any(!is.finite(x)) || !is.finite(sd(x)) || sd(x) <= 0)
    stop('Non-finite or constant variable in model input.')
  as.numeric((x - mean(x)) / sd(x))
}
score_column <- function(method, signature) paste(method, signature, sep = '__')

validate_inputs <- function() {
  wide <- read_embedded(score_csv)
  original <- read_embedded(original_singscore_csv)
  meta <- read_embedded(metadata_csv)
  master <- read_embedded(master_csv)
  coverage <- read_embedded(score_coverage_csv)
  expected_n <- c(SSL_HP=40L,H_ADJ_SSL_HP=13L,H_COL=14L,CA=12L)
  wide$Sample_column <- trimws(wide$Sample_column)
  wide$Patient_ID <- trimws(as.character(wide$Patient_ID))
  original$Sample_column <- trimws(original$Sample_column)
  required_scores <- as.vector(outer(METHODS, SIGNATURES, paste, sep='__'))
  stopifnot(nrow(wide)==79L, length(unique(wide$Patient_ID))==42L,
    !anyDuplicated(norm_sample(wide$Sample_column)), !anyNA(wide$Patient_ID),
    all(nzchar(wide$Patient_ID)), !any(grepl('UNMATCHED',wide$Patient_ID)),
    all(required_scores %in% names(wide)),
    all(as.integer(table(factor(wide$Group_final,levels=names(expected_n)))) == expected_n))
  if (any(!is.finite(as.matrix(wide[,required_scores])))) stop('Missing/non-finite saved scores.')
  mi <- match(norm_sample(wide$Sample_column), norm_sample(meta$Display))
  if (anyNA(mi) || anyDuplicated(norm_sample(meta$Display))) stop('Metadata sample matching is not unique.')
  group_meta <- pool_group(meta$Group[mi])
  override <- norm_sample(wide$Sample_column)=='B11_CA_DESC_6947'
  group_meta[override] <- 'CA'
  alignment <- data.frame(Sample=wide$Sample_column, Patient_score=wide$Patient_ID,
    Patient_metadata=trimws(as.character(meta$Patient[mi])), Group_score=wide$Group_final,
    Group_metadata=group_meta, CA6947_override=override)
  write_table(alignment,'qc/metadata_alignment.csv')
  if (any(alignment$Patient_score!=alignment$Patient_metadata) || any(alignment$Group_score!=alignment$Group_metadata))
    stop('Score metadata disagree with the supplied authoritative metadata.')
  j <- match(wide$Sample_column,original$Sample_column)
  stopifnot(!anyNA(j))
  for (sig in SIGNATURES) {
    if (max(abs(wide[[score_column('Singscore',sig)]] - original[[sig]][j])) > 1e-12)
      stop('Saved singscore tables disagree: ', sig)
  }
  audit_inputs <- list(INTERFERON_ALPHA_RESPONSE=ifna_filter_csv,
    INFLAMMATORY_RESPONSE=inflam_filter_csv,TNFA_SIGNALING_VIA_NFKB=tnf_filter_csv)
  expected_genes <- c(GASTRIC_METAPLASIA_24_GENE=24L,INTERFERON_ALPHA_RESPONSE=88L,
    INFLAMMATORY_RESPONSE=111L,TNFA_SIGNALING_VIA_NFKB=144L)
  for (sig in names(audit_inputs)) {
    a <- read_embedded(audit_inputs[[sig]])
    stopifnot(all(a$retained == (a$present_in_matrix & is.finite(a$mean_logCPM) & a$mean_logCPM>=2)),
      sum(a$retained)==expected_genes[[sig]])
    write_table(a,paste0('inputs/',sig,'_original_expression_filter_audit.csv'))
  }
  stopifnot(nrow(coverage)==12L,
    all(coverage$n_signature_used==unname(expected_genes[coverage$signature])))
  genes <- trimws(strsplit(gastric_genes_txt,'\n')[[1]])
  genes <- genes[nzchar(genes)]
  stopifnot(length(unique(genes))==24L)
  write_note(genes,'inputs/FINAL_24_GASTRIC_GENES.txt')
  write_table(coverage,'inputs/score_gene_coverage.csv')
  # These five historical SSL/HP aliases have the same patient, site, pooled
  # group and all four saved scores. No fuzzy matching or patient inference.
  aliases <- c(B2_HP_SIG_6830='B2_SSL_SIG_6830',B2_HP_SIG_7154='B2_SSL_SIG_7154',
    B2_HP_SIG_7784='B2_SSL_SIG_7784',B3_SSL_DESC_2102='B3_HP_DESC_2102',
    B3_SSL_TRVS_7241='B3_HP_TRVS_7241')
  master$Sample_Name <- trimws(master$Sample_Name)
  master$Patient <- trimws(as.character(master$Patient))
  mk <- norm_sample(master$Sample_Name)
  hit <- mk %in% names(aliases); mk[hit] <- unname(aliases[mk[hit]])
  match_wide <- match(mk,norm_sample(wide$Sample_column))
  if (anyNA(match_wide) || anyDuplicated(mk)) stop('Microbiome score samples do not match uniquely.')
  master_map <- c(Gastric_Metaplasia_Score='GASTRIC_METAPLASIA_24_GENE',
    Interferon_Alpha_Score='INTERFERON_ALPHA_RESPONSE',
    Inflammatory_Response_Score='INFLAMMATORY_RESPONSE',TNFA_NFKB_Score='TNFA_SIGNALING_VIA_NFKB')
  score_differences <- sapply(names(master_map),function(nm)
    master[[nm]] - wide[[score_column('Singscore',master_map[[nm]])]][match_wide])
  stopifnot(nrow(master)==63L,!any(master$Patient=='7699'),
    all(master$Patient==wide$Patient_ID[match_wide]),
    all(master$Group==wide$Group_final[match_wide]), max(abs(score_differences))<1e-8)
  write_table(data.frame(Sample_master=master$Sample_Name,Sample_score=wide$Sample_column[match_wide],
    Patient=master$Patient,Historical_SSL_HP_alias=hit,score_differences,check.names=FALSE),
    'qc/microbiome_and_clinical_score_agreement.csv')
  write_table(wide,'inputs/FROZEN_ALL_METHODS_SCORES.csv')
  write_table(master,'inputs/FROZEN_MASTER_METADATA_CH4.csv')
  write_table(meta,'inputs/FROZEN_RNA_METADATA.csv')
  write_table(as.data.frame(wide %>% group_by(Group_final) %>%
    summarise(Biopsies=n(),Patients=n_distinct(Patient_ID),.groups='drop')),
    'qc/colorectal_cohort_counts.csv')
  list(wide=wide,master=master,master_map=master_map)
}

fit_checked <- function(formula, data, id, REML=TRUE) {
  data$Patient_ID <- factor(data$Patient_ID)
  if (anyNA(data$Patient_ID) || any(!nzchar(as.character(data$Patient_ID)))) stop('Missing patient in ',id)
  if (!anyDuplicated(data$Patient_ID)) stop('No repeated patients: random-intercept model not identifiable in ',id)
  fit_once <- function(optimizer) {
    warnings <- messages <- character()
    m <- withCallingHandlers(lmerTest::lmer(formula,data=data,REML=REML,
      na.action=na.fail,control=lme4::lmerControl(optimizer=optimizer,
        optCtrl=if(optimizer=='bobyqa') list(maxfun=200000L) else list(maxeval=200000L),
        check.rankX='stop.deficient')),
      warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')},
      message=function(w){messages <<- c(messages,conditionMessage(w));invokeRestart('muffleMessage')})
    conv <- m@optinfo$conv$lme4$messages
    # Boundary messages are reported as singularity, not optimiser failure.
    conv <- conv[!grepl('boundary.*singular',conv)]
    opt <- m@optinfo$conv$opt
    list(model=m,warnings=warnings,messages=messages,conv=conv,
      converged=!length(conv) && (!length(opt)||all(opt==0)),optimizer=optimizer)
  }
  attempt <- tryCatch(fit_once('bobyqa'),error=function(e) list(error=conditionMessage(e),converged=FALSE))
  first_issue <- if (!is.null(attempt$error)) attempt$error else paste(attempt$conv,collapse='; ')
  if (!isTRUE(attempt$converged)) {
    retry <- tryCatch(fit_once('nloptwrap'),error=function(e) list(error=conditionMessage(e),converged=FALSE))
    if (isTRUE(retry$converged)) attempt <- retry
    else if (is.null(attempt$model) && !is.null(retry$model)) attempt <- retry
  }
  if (is.null(attempt$model)) stop('Model failed: ',id,'; ',first_issue)
  m <- attempt$model
  singular <- lme4::isSingular(m,tol=1e-4)
  v <- as.data.frame(lme4::VarCorr(m)); pv <- v$vcov[v$grp=='Patient_ID'][1]
  issues <- unique(c(attempt$warnings,attempt$conv))
  status <- if (!attempt$converged) 'UNAVAILABLE_NONCONVERGENCE' else
    if(length(issues)) 'REVIEW_WARNINGS' else if(singular) 'REVIEW_SINGULAR' else 'OK'
  qc <- data.frame(Model_ID=id,Formula=paste(deparse(formula),collapse=' '),REML=REML,
    N_biopsies=nrow(data),N_patients=nlevels(data$Patient_ID),
    N_patients_repeated=sum(table(data$Patient_ID)>1),Patient_variance=pv,
    Residual_variance=sigma(m)^2,ICC=pv/(pv+sigma(m)^2),
    Singular=singular,Converged=attempt$converged,Status=status,
    Optimizer=attempt$optimizer,First_attempt_issue=first_issue,
    Warnings=paste(attempt$warnings,collapse='; '),Messages=paste(attempt$messages,collapse='; '),
    Convergence_messages=paste(attempt$conv,collapse='; '),stringsAsFactors=FALSE)
  .rerun$models[[id]] <- list(model=m,data=data,qc=qc)
  write_note(capture.output(summary(m)),paste0('models/',id,'_summary.txt'))
  input_cols <- unique(c(intersect(c('Sample_column','Sample_Name','Patient_ID','Group_final','Group'),names(data)),
    all.vars(formula)))
  write_table(data[,input_cols,drop=FALSE],paste0('models/',id,'_inputs.csv'))
  list(model=m,data=data,qc=qc,available=attempt$converged)
}

focal_result <- function(object, term='x') {
  cf <- coef(summary(object$model,ddf='Satterthwaite'))
  if (!term %in% rownames(cf)) stop('Focal term missing: ',term)
  r <- cf[term,]; estimate <- unname(r['Estimate']); se <- unname(r['Std. Error']); df <- unname(r['df'])
  ok <- object$available && all(is.finite(c(estimate,se,df))) && se>0 && df>0
  ci <- if(ok) estimate+c(-1,1)*qt(.975,df)*se else c(NA_real_,NA_real_)
  data.frame(Estimate=estimate,SE=se,DF=df,CI_low=ci[1],CI_high=ci[2],
    P_model=if(ok) unname(r['Pr(>|t|)']) else NA_real_,
    N_biopsies=nrow(object$data),N_patients=length(unique(object$data$Patient_ID)),
    Singular=object$qc$Singular,Status=object$qc$Status,stringsAsFactors=FALSE)
}

run_groups <- function(wide) {
  result <- list(); all_pairs <- all_overall <- all_means <- legacy <- list()
  k<-0L
  for(method in METHODS) for(sig in SIGNATURES) {
    k<-k+1L; id<-sprintf('GROUP_%02d',k)
    d<-wide;d$Group_final<-factor(d$Group_final,levels=params$score_plot_groups)
    raw<-d[[score_column(method,sig)]]; d$y<-standardise(raw)
    obj<-fit_checked(y~Group_final+(1|Patient_ID),d,id,REML=TRUE)
    emm<-emmeans::emmeans(obj$model,~Group_final,lmer.df='satterthwaite')
    contrasts<-as.data.frame(summary(emmeans::contrast(emm,'pairwise',adjust='none'),
      infer=c(TRUE,TRUE),adjust='none'))
    parts<-strsplit(as.character(contrasts$contrast),' - ',fixed=TRUE)
    pairs<-data.frame(method=method,signature=sig,Model_ID=id,
      group1=vapply(parts,`[`,character(1),1L),group2=vapply(parts,`[`,character(1),2L),
      estimate=contrasts$estimate*sd(raw),SE=contrasts$SE*sd(raw),df=contrasts$df,
      CI_low=contrasts$lower.CL*sd(raw),CI_high=contrasts$upper.CL*sd(raw),
      statistic=contrasts$t.ratio,p_value=if(obj$available) contrasts$p.value else rep(NA_real_,6),
      Status=obj$qc$Status,Singular=obj$qc$Singular,stringsAsFactors=FALSE)
    pairs$p_adj<-bh_family(pairs$p_value,6L)
    pairs$label<-vapply(pairs$p_adj,p_label,character(1))
    pairs$key<-mapply(pair_key,pairs$group1,pairs$group2)
    # Independent-biopsy tests are comparison records, never plot fallback.
    oldp<-mapply(function(a,b) suppressWarnings(wilcox.test(raw[d$Group_final==a],raw[d$Group_final==b],
      exact=FALSE)$p.value),pairs$group1,pairs$group2)
    pairs$OLD_unpaired_P_BH_DO_NOT_REPORT<-bh_family(oldp,6L)
    pairs$Changed_significance_vs_old_unpaired<-ifelse(is.finite(pairs$p_adj),
      (pairs$p_adj<.05)!=(pairs$OLD_unpaired_P_BH_DO_NOT_REPORT<.05),NA)
    an<-as.data.frame(anova(obj$model,type=3,ddf='Satterthwaite'))
    aa<-an['Group_final',,drop=FALSE]
    overall<-data.frame(method=method,signature=sig,Model_ID=id,F=aa[['F value']],
      Numerator_DF=aa$NumDF,Denominator_DF=aa$DenDF,
      P=if(obj$available) aa[['Pr(>F)']] else NA_real_,Status=obj$qc$Status)
    means<-as.data.frame(confint(emm,adjust='none'))
    means$emmean<-means$emmean*sd(raw)+mean(raw);means$SE<-means$SE*sd(raw)
    means$lower.CL<-means$lower.CL*sd(raw)+mean(raw);means$upper.CL<-means$upper.CL*sd(raw)+mean(raw)
    means$method<-method;means$signature<-sig
    key<-paste(method,sig,sep='__')
    result[[key]]<-list(object=obj,pairwise=pairs,overall=overall,means=means)
    all_pairs[[key]]<-pairs;all_overall[[key]]<-overall;all_means[[key]]<-means
  }
  pairs<-bind_rows(all_pairs);overall<-bind_rows(all_overall)
  overall<-overall %>% group_by(method) %>% mutate(P_BH_four_scores_sensitivity=bh_family(P,4L)) %>% ungroup()
  write_table(pairs,'tables/GROUP_COMPARISONS_ALL_METHODS.csv')
  write_table(filter(pairs,method=='Singscore'),'tables/PRIMARY_SINGSCORE_GROUP_COMPARISONS.csv')
  write_table(overall,'tables/GROUP_OVERALL_TESTS.csv')
  write_table(bind_rows(all_means),'tables/GROUP_ESTIMATED_MARGINAL_MEANS.csv')
  result
}

run_model_associations <- function(wide) {
  rows<-old_rows<-mean_rows<-list(); k<-0L
  for(subset in names(SUBSETS)) for(xm in METHODS) for(ym in METHODS) for(sig in SIGNATURES[-1]) {
    k<-k+1L;id<-sprintf('ASSOC_%03d',k)
    d<-wide[wide$Group_final %in% SUBSETS[[subset]],,drop=FALSE]
    xc<-score_column(xm,SIGNATURES[1]);yc<-score_column(ym,sig)
    d$x_raw<-d[[xc]];d$y_raw<-d[[yc]]
    d$x<-standardise(d$x_raw);d$y<-standardise(d$y_raw)
    d$Group_final<-droplevels(factor(d$Group_final,levels=params$group_order))
    f<-if(nlevels(d$Group_final)>1L) y~x+Group_final+(1|Patient_ID) else y~x+(1|Patient_ID)
    obj<-fit_checked(f,d,id,REML=TRUE)
    meta<-data.frame(Model_ID=id,subset=subset,gastric_method=xm,inflammatory_method=ym,
      signature=sig,x_column=xc,y_column=yc,Group_adjusted=nlevels(d$Group_final)>1L)
    r<-cbind(meta,focal_result(obj))
    r$Pearson_r_descriptive<-cor(d$x_raw,d$y_raw,method='pearson')
    r$Spearman_rho_descriptive<-cor(d$x_raw,d$y_raw,method='spearman')
    r$Raw_scale_beta<-r$Estimate*sd(d$y_raw)/sd(d$x_raw)
    rows[[id]]<-r
    for(cm in c('pearson','spearman')) {
      ct<-suppressWarnings(cor.test(d$x_raw,d$y_raw,method=cm,exact=FALSE))
      old_rows[[length(old_rows)+1L]]<-cbind(meta,data.frame(correlation=cm,
        Coefficient=unname(ct$estimate),OLD_INDEPENDENT_BIOPSY_P_DO_NOT_REPORT=ct$p.value))
    }
    # Patient means answer a between-patient question within one tissue. Do not
    # pool patient-tissue means across groups and call them independent patients.
    if(subset=='ALL_COLORECTAL') for(gr in levels(d$Group_final)) {
      pm<-d[d$Group_final==gr,,drop=FALSE] %>% group_by(Patient_ID) %>%
        summarise(x=mean(x_raw),y=mean(y_raw),Biopsies=n(),.groups='drop')
      for(cm in c('pearson','spearman')) {
        ct<-suppressWarnings(cor.test(pm$x,pm$y,method=cm,exact=FALSE))
        mean_rows[[length(mean_rows)+1L]]<-data.frame(Group=gr,gastric_method=xm,
          inflammatory_method=ym,signature=sig,correlation=cm,N_patients=nrow(pm),
          Coefficient=unname(ct$estimate),P=ct$p.value)
      }
    }
  }
  tab<-bind_rows(rows) %>% group_by(subset,gastric_method,inflammatory_method) %>%
    mutate(P_BH_three=bh_family(P_model,3L)) %>%
    ungroup() %>% group_by(subset) %>% mutate(P_BH_all27_sensitivity=bh_family(P_model,27L)) %>% ungroup()
  old<-bind_rows(old_rows) %>% group_by(subset,gastric_method,inflammatory_method,correlation) %>%
    mutate(OLD_BH_THREE_DO_NOT_REPORT=bh_family(OLD_INDEPENDENT_BIOPSY_P_DO_NOT_REPORT,3L)) %>% ungroup()
  pm<-bind_rows(mean_rows) %>% group_by(Group,gastric_method,inflammatory_method,correlation) %>%
    mutate(P_BH_three=bh_family(P,3L)) %>% ungroup()
  write_table(tab,'sensitivity/SEPARATE_MIXED_MODEL_ASSOCIATIONS_ALL_METHODS.csv')
  write_table(filter(tab,gastric_method=='Singscore',inflammatory_method=='Singscore'),
    'sensitivity/SEPARATE_MIXED_MODEL_SINGSCORE_ASSOCIATIONS.csv')
  write_table(old,'comparison_only/OLD_NAIVE_BIOPSY_CORRELATIONS_DO_NOT_REPORT.csv')
  write_table(pm,'sensitivity/PATIENT_MEAN_CORRELATIONS_WITHIN_EACH_TISSUE.csv')
  write_note('ALL_WITH_STOMACH was not regenerated: the verified saved scoring dataset contains no stomach samples. It would duplicate ALL_COLORECTAL under a misleading label.',
    'qc/OLD_ALL_WITH_STOMACH_SUBSET.txt')
  tab
}

run_alpha <- function(master,score_map) {
  sets<-list(ALL_THREE_GROUPS=c('SSL_HP','H_ADJ_SSL_HP','H_COL'),
    SERRATED_COHORT=c('SSL_HP','H_ADJ_SSL_HP'),SSL_HP_ONLY='SSL_HP')
  rows<-list();k<-0L
  for(subset in names(sets)) for(metric in c('Shannon','Observed','Simpson')) for(score in names(score_map)) {
    k<-k+1L;id<-sprintf('ALPHA_%02d',k)
    d<-master[master$Group %in% sets[[subset]],,drop=FALSE]
    d$Patient_ID<-d$Patient;d$Group_final<-droplevels(factor(d$Group,levels=params$group_order))
    d$BRAF_01<-ifelse(trimws(d$BRAF)=='V600E',1,ifelse(trimws(d$BRAF)=='WT',0,NA_real_))
    keep<-is.finite(d[[score]]) & is.finite(d[[metric]]) &
      is.finite(d$Reads_After_Selective_Rarefaction_Pre5pct)
    if(subset=='SSL_HP_ONLY') keep<-keep & is.finite(d$BRAF_01)
    write_table(data.frame(Sample=d$Sample_Name,Patient=d$Patient,Included=keep),paste0('qc/',id,'_inclusion.csv'))
    d<-d[keep,,drop=FALSE];d$y<-standardise(d[[score]]);d$x<-standardise(d[[metric]])
    d$Reads_z<-standardise(d$Reads_After_Selective_Rarefaction_Pre5pct)
    full<-if(subset=='SSL_HP_ONLY') y~x+BRAF_01+Reads_z+(1|Patient_ID) else y~x+Group_final+Reads_z+(1|Patient_ID)
    null<-if(subset=='SSL_HP_ONLY') y~BRAF_01+Reads_z+(1|Patient_ID) else y~Group_final+Reads_z+(1|Patient_ID)
    fit<-fit_checked(full,d,id,REML=FALSE)
    reduced<-fit_checked(null,d,paste0(id,'_NULL'),REML=FALSE)
    lr<-as.data.frame(anova(reduced$model,fit$model,refit=FALSE))
    plrt<-if(fit$available && reduced$available) lr[['Pr(>Chisq)']][2] else NA_real_
    r<-cbind(data.frame(Model_ID=id,subset=subset,signature=unname(score_map[[score]]),
      score_column=score,metric=metric,P_LRT=plrt),focal_result(fit))
    r$Spearman_rho_descriptive<-cor(d[[metric]],d[[score]],method='spearman')
    rows[[id]]<-r
  }
  tab<-bind_rows(rows) %>% group_by(subset) %>% mutate(P_Holm12=bh_family(P_LRT,12L,'holm')) %>% ungroup()
  write_table(tab,'tables/ALPHA_DIVERSITY_ALL_36_PATIENT_MODELS.csv')
  write_table(filter(tab,subset=='ALL_THREE_GROUPS',metric=='Shannon'),'tables/SF_4_6_3_SHANNON_RESULTS.csv')
  ref<-read_embedded(reference_alpha_csv)
  old<-data.frame(score_column=ref$Outcome,metric=ref$Predictor,
    OLD_beta=ref$Standardized_beta,OLD_P_LRT=ref$P_LRT,OLD_P_Holm12=ref$P_Holm_analysis_family,
    OLD_N_biopsies=ref$N_samples,OLD_N_patients=ref$N_patients)
  comparison<-merge(tab[tab$subset=='ALL_THREE_GROUPS',],old,by=c('score_column','metric'),all=TRUE)
  comparison$Beta_difference<-comparison$Estimate-comparison$OLD_beta
  comparison$P_difference<-comparison$P_LRT-comparison$OLD_P_LRT
  comparison$Counts_agree<-comparison$N_biopsies==comparison$OLD_N_biopsies &
    comparison$N_patients==comparison$OLD_N_patients
  comparison$Numerically_reproduced<-abs(comparison$Beta_difference)<1e-6 &
    abs(comparison$P_difference)<1e-6 & comparison$Counts_agree
  write_table(comparison,'qc/ALPHA_MODELS_AGREEMENT_WITH_SAVED_RESULTS.csv')
  if(nrow(comparison)!=12L || anyNA(comparison$Numerically_reproduced) || !all(comparison$Numerically_reproduced))
    warning('Alpha models differ from the saved reference beyond 1e-6: review qc/ALPHA_MODELS_AGREEMENT_WITH_SAVED_RESULTS.csv.')
  tab
}

save_plot_set <- function(plot,relative,width_mm,height_mm) {
  if(!MAKE_PLOTS) return(invisible(NULL))
  stem<-file.path(.rerun$out,relative);dir.create(dirname(stem),recursive=TRUE,showWarnings=FALSE)
  for(ext in FIGURE_FORMATS) {
    args<-list(filename=paste0(stem,'.',ext),plot=plot,width=width_mm,height=height_mm,
      units='mm',bg='white',limitsize=FALSE)
    if(ext=='pdf') args$device<-if(capabilities('cairo'))grDevices::cairo_pdf else grDevices::pdf
    else {
      args$dpi<-FIGURE_DPI
      if(ext=='png') args$device<-ragg::agg_png
      if(ext %in% c('jpg','jpeg')) {args$device<-ragg::agg_jpeg;args$quality<-100}
      if(ext=='tiff') {args$device<-ragg::agg_tiff;args$compression<-'lzw'}
    }
    do.call(ggplot2::ggsave,args)
  }
}

make_association_plot <- function(wide,r,tag=NULL) {
  d<-wide[wide$Group_final %in% SUBSETS[[r$subset]],,drop=FALSE]
  d$Group_plot<-factor(d$Group_final,levels=params$group_order)
  present<-params$group_order[params$group_order %in% as.character(d$Group_plot)]
  same_singscore<-r$gastric_method=='Singscore' && r$inflammatory_method=='Singscore'
  pretty<-signature_pretty_label(r$signature)
  short<-c(Singscore='singscore',Zscore='z-score',SSGSEA='ssGSEA')
  title<-if(same_singscore) pretty else paste0(pretty,': ',short[[r$gastric_method]],' vs ',short[[r$inflammatory_method]])
  compact_p<-function(p) if(!is.finite(p)) '=NA' else if(p<.001) '<.001' else
    paste0('=',sub('^0','',formatC(p,format='fg',digits=2)))
  annotation<-paste0('n=',r$N_biopsies,'\nN=',r$N_patients,
    '\n\u03c1=',sub('^(-?)0','\\1',sprintf('%.3f',r$Spearman_rho)),
    '\nP',compact_p(r$Spearman_P),'\nq',compact_p(r$Spearman_q))
  scatter<-ggplot(d,aes(x=.data[[r$x_column]],y=.data[[r$y_column]],colour=Group_plot))+
    geom_point(size=1.5,alpha=.88)+
    scale_colour_manual(values=params$group_colours[present],breaks=present,
      labels=params$group_labels[present],drop=TRUE)+
    labs(title=title,x=paste('Gastric metaplasia',short[[r$gastric_method]]),
      y=paste(pretty,short[[r$inflammatory_method]]),colour=NULL,tag=tag)+paper_theme()+
    theme(legend.position='bottom',legend.direction='horizontal',legend.margin=margin(0,0,0,0),
      legend.text=element_text(size=6),legend.key.size=grid::unit(2.8,'mm'),
      legend.spacing.x=grid::unit(.8,'mm'),plot.margin=margin(4,1,4,4,'mm'),
      plot.tag=element_text(size=8))+
    guides(colour=guide_legend(nrow=1,byrow=TRUE,override.aes=list(size=2,alpha=1)))
  # Preserve the original grey descriptive OLS lines. They do not supply
  # any model coefficient, P value or interval reported in the side panel.
  scatter<-scatter+geom_smooth(method='lm',formula=y~x,se=FALSE,colour='grey40',linewidth=.45)
  side<-ggplot()+annotate('text',x=0,y=1,label=annotation,hjust=0,vjust=1,
    family=params$font_family,size=pt_to_mm(params$label_pt),lineheight=1.05)+
    xlim(0,1)+ylim(0,1)+theme_void(base_family=params$font_family)+
    theme(plot.margin=margin(6,1,4,0,'mm'))
  (scatter|side)+plot_layout(widths=c(1,.31))
}

make_alpha_plot <- function(master,r) {
  saved<-.rerun$models[[r$Model_ID]];d<-saved$data;m<-saved$model
  x<-r$metric;y<-r$score_column
  # The same black line/grey ribbon format now uses population-level mixed-
  # model predictions averaged equally across tissue groups, at mean depth.
  xx<-seq(min(d[[x]]),max(d[[x]]),length.out=100L)
  pred<-as.data.frame(confint(emmeans::emmeans(m,~x,
    at=list(x=(xx-mean(d[[x]]))/sd(d[[x]]),Reads_z=0),
    weights='equal',lmer.df='satterthwaite'),adjust='none'))
  line<-data.frame(x=pred$x*sd(d[[x]])+mean(d[[x]]),
    y=pred$emmean*sd(d[[y]])+mean(d[[y]]),
    lo=pred$lower.CL*sd(d[[y]])+mean(d[[y]]),hi=pred$upper.CL*sd(d[[y]])+mean(d[[y]]))
  ylabel<-c(GASTRIC_METAPLASIA_24_GENE='Gastric metaplasia singscore',
    INTERFERON_ALPHA_RESPONSE='Interferon-alpha singscore',INFLAMMATORY_RESPONSE='Inflammatory-response singscore',
    TNFA_SIGNALING_VIA_NFKB='TNFA/NFKB singscore')[[r$signature]]
  caption<-paste0('Patient-aware beta = ',sprintf('%.3f',r$Estimate),
    '; 95% CI [',sprintf('%.3f',r$CI_low),', ',sprintf('%.3f',r$CI_high),']',
    '\nP (LRT) ',format_p_text(r$P_LRT),'; Holm P (12 tests) ',format_p_text(r$P_Holm12),
    '\nn = ',r$N_biopsies,' biopsies; ',r$N_patients,' patients',
    if(r$Status!='OK') '\nFit flagged for review; see diagnostics.' else '')
  ggplot(d,aes(x=.data[[x]],y=.data[[y]],colour=Group_final))+
    geom_ribbon(data=line,aes(x=x,ymin=lo,ymax=hi),inherit.aes=FALSE,fill='grey80',alpha=.35)+
    geom_line(data=line,aes(x=x,y=y),inherit.aes=FALSE,colour='black',linewidth=.65)+
    geom_point(size=2.1,alpha=.80)+
    scale_colour_manual(values=params$group_colours,breaks=params$group_order[1:3],drop=TRUE)+
    labs(title=paste0(ylabel,' vs Shannon diversity'),subtitle='SSL_HP, H_ADJ_SSL_HP and H_COL',
      x='Shannon diversity',y=ylabel,colour='Tissue group',caption=caption)+
    theme_classic(base_family=params$font_family,base_size=10)+
    theme(text=element_text(family=params$font_family,colour='black'),
      plot.title=element_text(face='bold',size=11,hjust=0),plot.subtitle=element_text(size=9,colour='black'),
      axis.title=element_text(size=10,colour='black'),axis.text=element_text(size=9,colour='black'),
      legend.position='bottom',legend.title=element_text(size=9),legend.text=element_text(size=8),
      plot.caption.position='plot',plot.caption=element_text(size=3.2*72.27/25.4,hjust=1,lineheight=1.15,
        margin=margin(t=3,unit='mm')),plot.margin=margin(5,6,5,5,unit='mm'))
}

export_figures <- function(wide,master,groups,associations,alpha) {
  plots<-list()
  for(method in METHODS) for(sig in SIGNATURES) {
    key<-paste(method,sig,sep='__');d<-wide;d$score<-d[[score_column(method,sig)]]
    p<-make_generic_score_distribution_plot(d,groups[[key]],sig,method)
    plots[[key]]<-p
    if(method=='Singscore') save_plot_set(p,paste0('figures/individual/',sig,'__Singscore_distribution'),82,88)
  }
  inf_order<-c('INFLAMMATORY_RESPONSE','INTERFERON_ALPHA_RESPONSE','TNFA_SIGNALING_VIA_NFKB')
  panel<-wrap_plots(plots[paste('Singscore',inf_order,sep='__')],ncol=3)
  save_plot_set(panel,'figures/thesis/Figure_3_5_2_inflammatory_singscores_PATIENT_AWARE',183,96)
  d<-wide;d$centered_singscore<-d[[score_column('Singscore',SIGNATURES[1])]]
  gp<-make_score_distribution_plot(d,groups[[paste('Singscore',SIGNATURES[1],sep='__')]]$pairwise)
  save_plot_set(gp,'figures/thesis/Figure_3_9B_gastric_singscore_PATIENT_AWARE',105,95)
  save_plot_set(gp,'figures/individual/Panel_D_metaplasia_singscore_distribution',90,88)
  # Method-specific grids keep the original 4-column score comparison layout.
  for(method in METHODS[-1]) {
    pp<-wrap_plots(plots[paste(method,SIGNATURES,sep='__')],ncol=4)+
      plot_annotation(title=score_method_y_label(method))
    save_plot_set(pp,paste0('figures/sensitivity/COMBINED_ALL_SIGNATURES__',method),220,96)
  }
  for(subset in names(SUBSETS)) {
    pp<-list()
    for(sig in SIGNATURES[-1]) {
      r<-associations[associations$subset==subset & associations$signature==sig &
        associations$gastric_method=='Singscore' & associations$inflammatory_method=='Singscore',,drop=FALSE]
      p<-make_association_plot(wide,r,tag=NULL);pp[[sig]]<-p
      save_plot_set(p,paste0('figures/correlations/',subset,'_METAPLASIA_vs_',sig),108,78)
    }
    combined<-wrap_plots(pp,ncol=3)
    stem<-switch(subset,SSL_HP_ONLY='Figure_3_10_SSL_HP_gastric_inflammation_PATIENT_AWARE',
      ALL_COLORECTAL='SF_3_10_all_colorectal_gastric_inflammation_PATIENT_AWARE',
      SSL_HP_CA='SUPPLEMENT_SSL_HP_CA_gastric_inflammation_PATIENT_AWARE')
    save_plot_set(combined,paste0('figures/thesis/',stem),183,67)
    # All 3x3 scoring-method combinations are recalculated and exported as
    # grids, without regenerating unrelated figures from the original pipeline.
    for(sig in SIGNATURES[-1]) {
      pp<-list()
      for(xm in c('Zscore','Singscore','SSGSEA')) for(ym in c('Zscore','Singscore','SSGSEA')) {
        r<-associations[associations$subset==subset & associations$signature==sig &
          associations$gastric_method==xm & associations$inflammatory_method==ym,,drop=FALSE]
        pp[[paste(xm,ym)]]<-make_association_plot(wide,r,tag=NULL)
      }
      grid<-wrap_plots(pp,ncol=3)+plot_annotation(title=paste(subset,signature_pretty_label(sig)))
      save_plot_set(grid,paste0('figures/sensitivity/',subset,'__3x3__',sig),255,215)
    }
  }
  ap<-list()
  for(sig in SIGNATURES) {
    r<-alpha[alpha$subset=='ALL_THREE_GROUPS' & alpha$metric=='Shannon' & alpha$signature==sig,,drop=FALSE]
    p<-make_alpha_plot(master,r);ap[[sig]]<-p
    save_plot_set(p,paste0('figures/thesis/SF_4_6_3_',sig,'_PATIENT_AWARE'),6.4*25.4,6.2*25.4)
  }
  if(MAKE_PLOTS) {
    path<-file.path(.rerun$out,'figures/thesis/SF_4_6_3_four_pages_PATIENT_AWARE.pdf')
    grDevices::cairo_pdf(path,width=6.4,height=6.2,family=params$font_family)
    tryCatch(for(p in ap) print(p),finally=grDevices::dev.off())
  }
}

export_diagnostics <- function() {
  qc<-bind_rows(lapply(.rerun$models,`[[`,'qc'))
  write_table(qc,'qc/ALL_MODEL_FIT_CHECKS.csv')
  write_table(filter(qc,Status!='OK'),'qc/MODELS_REQUIRING_REVIEW.csv')
  saveRDS(.rerun$models,file.path(.rerun$out,'models/ALL_FITTED_MODELS_AND_INPUTS.rds'),compress='xz')
  if(!MAKE_PLOTS) return(invisible(qc))
  path<-file.path(.rerun$out,'qc/PRIMARY_MODEL_DIAGNOSTICS.pdf')
  grDevices::cairo_pdf(path,width=9,height=4.5,family=params$font_family)
  on.exit(grDevices::dev.off(),add=TRUE)
  for(id in .rerun$diagnostic_ids) {
    z<-.rerun$models[[id]];m<-z$model
    par(mfrow=c(1,2),mar=c(4,4,4,1))
    plot(fitted(m),residuals(m),pch=16,cex=.7,xlab='Fitted value',ylab='Residual',main=paste(id,z$qc$Status))
    abline(h=0,lty=2,col='grey50');lines(lowess(fitted(m),residuals(m)),col='red',lwd=1)
    qqnorm(residuals(m),pch=16,cex=.7,main=paste(id,'Residual Q-Q plot'));qqline(residuals(m),col='red')
  }
  invisible(qc)
}

run_reanalysis <- function() {
  set.seed(SET_SEED)
  dir.create(OUTPUT_PARENT,recursive=TRUE,showWarnings=FALSE)
  run_id<-paste0('OSLC_COMPLETE_SCORE_ANALYSIS_',format(Sys.time(),'%Y%m%d_%H%M%S'))
  out<-file.path(OUTPUT_PARENT,run_id);i<-0L
  while(dir.exists(out)) {i<-i+1L;out<-file.path(OUTPUT_PARENT,paste0(run_id,'_',i))}
  dir.create(out,recursive=TRUE);.rerun$out<-normalizePath(out,winslash='/',mustWork=TRUE)
  .rerun$models<-list();.rerun$diagnostic_ids<-character()
  for(d in c('inputs','qc','tables','models','figures','sensitivity','comparison_only')) dir.create(file.path(out,d))
  if(!is.na(THIS_SCRIPT)) file.copy(THIS_SCRIPT,file.path(out,'OSLC_complete_score_analysis.R'),overwrite=FALSE)
  err<-NULL
  tryCatch({
    required<-c('ggplot2','dplyr','tidyr','patchwork','lme4','lmerTest','emmeans','zip','ragg')
    missing<-required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
    if(length(missing) && INSTALL_MISSING) install.packages(missing,repos='https://cloud.r-project.org')
    missing<-required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
    if(length(missing)) stop('Missing packages: ',paste(missing,collapse=', '),'. Install these and run this file again.')
    suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(patchwork)})
    if(packageVersion('ggplot2')<'3.4.0') stop('ggplot2 >= 3.4.0 is required.')
    if(.Platform$OS.type=='windows') grDevices::windowsFonts(Arial=grDevices::windowsFont('Arial'))
    if(requireNamespace('systemfonts',quietly=TRUE) && !'Arial' %in% systemfonts::system_fonts()$family) {
      params$font_family <<- if('Liberation Sans' %in% systemfonts::system_fonts()$family) 'Liberation Sans' else 'sans'
    }
    if(!capabilities('cairo')) stop('Cairo PDF support is required for matching figure typography.')
    write_note(paste('Font used:',params$font_family),'qc/FONT_USED.txt')
    message('1/6 Validating inputs and using score entry point: ',SCORE_INPUT_MODE)
    dat<-validate_inputs()
    message('2/6 Fitting the same patient-aware group models for all four signatures...')
    gr<-run_groups(dat$wide)
    old_gastric<-read_embedded(reference_gastric_csv)
    new_gastric<-gr[[paste('Singscore',SIGNATURES[1],sep='__')]]$pairwise
    gastric_comparison<-merge(new_gastric,old_gastric,by=c('group1','group2'),suffixes=c('_new','_old'))
    gastric_comparison$Same_significance<-(gastric_comparison$p_adj_new<.05)==(gastric_comparison$p_adj_old<.05)
    write_table(gastric_comparison,'qc/GASTRIC_COMPARISON_WITH_OLD_PATIENT_ADJUSTED_TABLE.csv')
    write_note('The old emmeans call left its degrees-of-freedom method at the package default. This rerun explicitly uses Satterthwaite for every signature. Gastric estimates should agree, but SEs, degrees of freedom and P values can differ from older default-based output. Use the new table consistently.',
      'qc/GASTRIC_EMMEANS_DF_NOTE.txt')
    message('3/6 Calculating clustered Spearman and Pearson tests, plus separate model results...')
    aa<-run_associations(dat$wide)
    message('4/6 Reproducing the existing patient-aware alpha-diversity models...')
    al<-run_alpha(dat$master,dat$master_map)
    .rerun$diagnostic_ids<-c(vapply(gr[paste('Singscore',SIGNATURES,sep='__')],function(x)x$object$qc$Model_ID,character(1)),
      aa$Model_ID[aa$gastric_method=='Singscore' & aa$inflammatory_method=='Singscore'],
      al$Model_ID[al$subset=='ALL_THREE_GROUPS'])
    message('5/6 Saving figures in the original formats and recording model diagnostics...')
    export_figures(dat$wide,dat$master,gr,aa,al)
    qc<-export_diagnostics()
    write_records(aa,al)
    status<-if(any(qc$Status!='OK')) 'COMPLETED_WITH_FITS_TO_REVIEW' else 'COMPLETED_PENDING_DIAGNOSTIC_REVIEW'
    write_note(c(status,paste('Models:',nrow(qc)),paste('Flagged:',sum(qc$Status!='OK')),
      'The run completed. Review the saved model diagnostics before using the results.',
      paste('Figures generated:',MAKE_PLOTS)), 'RUN_STATUS.txt')
    write_table(data.frame(Package=required,Version=vapply(required,function(p)as.character(packageVersion(p)),character(1))),
      'qc/PACKAGE_VERSIONS.csv')
  },error=function(e) {
    err <<- conditionMessage(e)
    write_note(c('FAILED',err),'RUN_STATUS.txt')
    if(length(.rerun$models)) {
      write_table(do.call(rbind,lapply(.rerun$models,`[[`,'qc')),'qc/ALL_MODEL_FIT_CHECKS.csv')
      saveRDS(.rerun$models,file.path(out,'models/PARTIAL_FITTED_MODELS.rds'))
    }
  })
  write_note(capture.output(sessionInfo()),'qc/SESSION_INFO.txt')
  files<-list.files(out,recursive=TRUE,full.names=TRUE)
  write_table(data.frame(File=substring(files,nchar(out)+2),MD5=unname(tools::md5sum(files))), 'qc/OUTPUT_CHECKSUMS.csv')
  message('6/6 Creating the review ZIP...')
  zip_path<-paste0(.rerun$out,'.zip')
  if(requireNamespace('zip',quietly=TRUE)) {
    zip::zipr(zipfile=zip_path,files=basename(.rerun$out),root=dirname(.rerun$out))
    message('Review ZIP: ',zip_path)
  } else message('ZIP package unavailable. Zip this output folder manually: ',.rerun$out)
  message('Output folder: ',.rerun$out)
  if(!is.null(err)) stop('Rerun stopped: ',err,'\nThe partial output ZIP records the error.',call.=FALSE)
  invisible(.rerun$out)
}



# Original score construction, isolated from figure/test functions.
make_original_scoring_environment <- function() {
 e<-new.env(parent=globalenv())
 evalq({
 params<-list(cpm_file=CPM_FILE,gene_length_file=GENE_LENGTH_FILE,metadata_file=METADATA_FILE,
 out_dir='',hallmark_mean_logcpm_cutoff=2,tpm_min_gene_length_coverage=.80,
 require_all_final_signature_genes_measured=TRUE,cpm_sum_expected=1e6,cpm_sum_tolerance_fraction=.05,
 ssgsea_norm=TRUE,ssgsea_alpha=.25,ssgsea_max_genes=Inf,groups_to_exclude='DELETE',
 group_order=c('SSL_HP','H_ADJ_SSL_HP','H_COL','CA','STOMACH'),
 patient_id_candidates=c('Patient_ID','PatientID','patient_id','Patient','patient'),
 manual_group_overrides=data.frame(sample_key='B11_CA_DESC_6947',forced_group='CA'))
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

safe_mean <- function(x) {
  x <- to_num(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_max <- function(x) {
  x <- to_num(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

safe_sd <- function(x) {
  x <- to_num(x)
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  stats::sd(x)
}

normalise_name <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

normalise_sample_id <- function(x) {
  x <- basename(as.character(x))
  x <- toupper(trimws(x))
  x <- gsub("\\.(BAM|TXT|TSV|CSV|GZ|FASTQ|FQ)$", "", x)
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

strip_batch_prefix <- function(x) {
  gsub("^B[0-9]+_", "", normalise_sample_id(x))
}

collapse_group <- function(x) {
  x <- normalise_name(x)
  dplyr::case_when(
    x %in% c("SSL", "HP", "SSL_HP", "SSLHP") ~ "SSL_HP",
    x %in% c("H_ADJ_SSL", "H_ADJ_HP", "H_ADJ_SSL_HP", "HADJSSLHP") ~ "H_ADJ_SSL_HP",
    x %in% c("H_COL", "HCOL", "HEALTHY_COLON") ~ "H_COL",
    x %in% c("CA", "CONVENTIONAL_ADENOMA") ~ "CA",
    x %in% c("ANTRUM", "STOMACH", "GASTRIC_ANTRUM") ~ "STOMACH",
    TRUE ~ NA_character_
  )
}

parse_group_from_sample <- function(x) {
  nx <- normalise_sample_id(x)
  dplyr::case_when(
    grepl("H_COL", nx) ~ "H_COL",
    grepl("H_ADJ", nx) ~ "H_ADJ_SSL_HP",
    grepl("(^|_)CA(_|$)", nx) ~ "CA",
    grepl("(^|_)SSL(_|$)", nx) | grepl("(^|_)HP(_|$)", nx) ~ "SSL_HP",
    grepl("ANTRUM", nx) | grepl("STOMACH", nx) ~ "STOMACH",
    TRUE ~ NA_character_
  )
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

write_tsv_safe <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(df, path, sep = "\t")
}

write_lines_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(as.character(x), con = path)
}

pt_to_mm <- function(x) x / 2.845276
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
    # Preserve unmapped Ensembl genes as stable background identifiers. They
    # cannot match HGNC signatures, but they should remain in the all-gene
    # within-sample ranking background used by singscore.
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

collapse_expression_by_symbol <- function(expr_mat, gene_symbols) {
  gene_symbols <- toupper(trimws(as.character(gene_symbols)))
  keep <- !is.na(gene_symbols) & nzchar(gene_symbols)

  df <- data.frame(
    gene_symbol = gene_symbols[keep],
    as.data.frame(expr_mat[keep, , drop = FALSE], check.names = FALSE),
    check.names = FALSE
  )

  df <- df %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), safe_mean), .groups = "drop")

  mat <- as.matrix(df[, -1, drop = FALSE])
  rownames(mat) <- df$gene_symbol
  storage.mode(mat) <- "double"
  mat
}

# =============================================================================
# 5. HEADER-SAFE CPM READER
# =============================================================================

is_blank_field <- function(x) {
  is.na(x) | !nzchar(trimws(as.character(x)))
}

count_delimiter <- function(line, delimiter) {
  hits <- gregexpr(delimiter, line, fixed = TRUE)[[1]]
  if (length(hits) == 1 && hits[1] == -1) return(0L)
  length(hits)
}

split_fields <- function(line, delimiter) {
  line <- sub(paste0("^", intToUtf8(65279)), "", line)
  line <- sub("\r$", "", line)
  strsplit(line, delimiter, fixed = TRUE)[[1]]
}

trim_trailing_blank_fields <- function(x) {
  while (length(x) > 0 && is_blank_field(x[length(x)])) {
    x <- x[-length(x)]
  }
  x
}

header_field_is_gene_label <- function(x) {
  normalise_name(x) %in% c(
    "GENE", "GENE_ID", "GENEID", "ID", "ENSEMBL", "ENSEMBL_ID",
    "ENSEMBL_GENE_ID", "GENE_SYMBOL", "SYMBOL", "GENE_NAME"
  )
}

looks_like_ensembl_gene_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("\\..*$", "", x)
  grepl("^ENSG[0-9]+$", x)
}

read_cpm_header_safe <- function(path, gene_col_name = "gene_id") {
  probe <- readLines(path, n = 200, warn = FALSE)
  nonempty <- which(nzchar(trimws(probe)))
  noncomment <- nonempty[!grepl("^\\s*#", probe[nonempty])]

  if (length(noncomment) < 2) {
    stop("The CPM file contains fewer than two non-empty data lines: ", path)
  }

  header_idx <- noncomment[1]
  data_idx <- noncomment[2]
  header_line <- probe[header_idx]
  data_line <- probe[data_idx]

  tab_count <- count_delimiter(header_line, "\t") + count_delimiter(data_line, "\t")
  comma_count <- count_delimiter(header_line, ",") + count_delimiter(data_line, ",")
  sep <- if (tab_count >= comma_count) "\t" else ","

  header_raw <- split_fields(header_line, sep)
  header_trimmed <- trim_trailing_blank_fields(header_raw)
  data_fields <- split_fields(data_line, sep)

  first_header_blank <- length(header_raw) > 0 && is_blank_field(header_raw[1])
  first_header_gene <- length(header_raw) > 0 && header_field_is_gene_label(header_raw[1])
  first_data_ensembl <- length(data_fields) > 0 && looks_like_ensembl_gene_id(data_fields[1])

  mode <- "HEADER_OK"

  if (length(header_raw) == length(data_fields) - 1) {
    mode <- "MISSING_GENE_HEADER_PREPENDED"
    df <- fread(
      path, sep = sep, header = FALSE, skip = header_idx,
      col.names = c(gene_col_name, header_raw),
      data.table = FALSE, fill = TRUE, check.names = FALSE
    )
  } else if (
    length(header_trimmed) == length(data_fields) - 1 &&
      !first_header_blank && !first_header_gene && first_data_ensembl
  ) {
    mode <- "MISSING_GENE_HEADER_PREPENDED_TRAILING_BLANK_DROPPED"
    df <- fread(
      path, sep = sep, header = FALSE, skip = header_idx,
      col.names = c(gene_col_name, header_trimmed),
      data.table = FALSE, fill = TRUE, check.names = FALSE
    )
  } else {
    df <- fread(
      path, sep = sep, header = TRUE, skip = header_idx - 1,
      data.table = FALSE, fill = TRUE, check.names = FALSE
    )

    if (is_blank_field(names(df)[1]) || header_field_is_gene_label(names(df)[1])) {
      names(df)[1] <- gene_col_name
      mode <- "GENE_HEADER_PRESENT_OR_BLANK_RENAMED"
    }
  }

  if (names(df)[1] != gene_col_name) {
    first_values <- as.character(df[[1]][seq_len(min(25, nrow(df)))])
    if (mean(looks_like_ensembl_gene_id(first_values), na.rm = TRUE) > 0.5) {
      names(df)[1] <- gene_col_name
      mode <- paste0(mode, "_FIRST_COLUMN_RENAMED_BY_ENSEMBL_CONTENT")
    }
  }

  if (ncol(df) < 2) {
    stop("The corrected CPM file has fewer than two columns.")
  }

  qc <- data.frame(
    file = path,
    delimiter = ifelse(sep == "\t", "tab", "comma"),
    reader_mode = mode,
    header_field_count = length(header_raw),
    first_data_field_count = length(data_fields),
    first_corrected_column = names(df)[1],
    first_corrected_sample = names(df)[2],
    stringsAsFactors = FALSE
  )

  list(df = df, qc = qc)
}

# =============================================================================
# 6. METADATA ALIGNMENT
# =============================================================================

load_and_align_metadata <- function(path, sample_columns) {
  meta <- read_table_flexible(path)

  sample_col <- pick_first_existing(
    c("Display", "Sample", "sample", "SampleID", "sample_id", "expr_name_norm"),
    meta
  )
  group_col <- pick_first_existing(
    c("Group", "GROUP", "Group2", "group"),
    meta
  )
  patient_col <- pick_first_existing(params$patient_id_candidates, meta)

  if (is.na(sample_col)) stop("No sample-ID column was detected in the metadata.")
  if (is.na(group_col)) stop("No group column was detected in the metadata.")

  meta_core <- data.frame(
    Sample_raw = as.character(meta[[sample_col]]),
    Sample_norm = normalise_sample_id(meta[[sample_col]]),
    Sample_strip_batch = strip_batch_prefix(meta[[sample_col]]),
    Group_raw = as.character(meta[[group_col]]),
    Group_final = collapse_group(meta[[group_col]]),
    Patient_ID = if (!is.na(patient_col)) as.character(meta[[patient_col]]) else NA_character_,
    stringsAsFactors = FALSE
  )

  meta_core$exclude <- normalise_name(meta_core$Group_raw) %in%
    normalise_name(params$groups_to_exclude)
  meta_keep <- meta_core[!meta_core$exclude, , drop = FALSE]
  meta_excluded <- meta_core[meta_core$exclude, , drop = FALSE]

  if (anyDuplicated(meta_keep$Sample_norm)) {
    dup <- unique(meta_keep$Sample_norm[duplicated(meta_keep$Sample_norm)])
    stop("Duplicated metadata sample IDs after normalisation: ", paste(dup, collapse = ", "))
  }

  sample_meta <- data.frame(
    Sample_column = sample_columns,
    Sample_norm = normalise_sample_id(sample_columns),
    Sample_strip_batch = strip_batch_prefix(sample_columns),
    stringsAsFactors = FALSE
  )

  exact <- meta_keep %>%
    dplyr::select(
      Sample_norm,
      Meta_sample_exact = Sample_raw,
      Group_raw_exact = Group_raw,
      Group_final_exact = Group_final,
      Patient_ID_exact = Patient_ID
    )

  strip_counts <- meta_keep %>%
    dplyr::count(Sample_strip_batch, name = "n_strip")

  stripped <- meta_keep %>%
    dplyr::left_join(strip_counts, by = "Sample_strip_batch") %>%
    dplyr::filter(n_strip == 1) %>%
    dplyr::select(
      Sample_strip_batch,
      Meta_sample_strip = Sample_raw,
      Group_raw_strip = Group_raw,
      Group_final_strip = Group_final,
      Patient_ID_strip = Patient_ID
    )

  sample_meta <- sample_meta %>%
    dplyr::left_join(exact, by = "Sample_norm") %>%
    dplyr::left_join(stripped, by = "Sample_strip_batch") %>%
    dplyr::mutate(
      Sample_raw = dplyr::coalesce(Meta_sample_exact, Meta_sample_strip, Sample_column),
      Group_raw = dplyr::coalesce(Group_raw_exact, Group_raw_strip),
      Group_final = dplyr::coalesce(
        Group_final_exact,
        Group_final_strip,
        parse_group_from_sample(Sample_column)
      ),
      Patient_ID = dplyr::coalesce(Patient_ID_exact, Patient_ID_strip),
      metadata_match_mode = dplyr::case_when(
        !is.na(Meta_sample_exact) ~ "exact_normalised_match",
        is.na(Meta_sample_exact) & !is.na(Meta_sample_strip) ~ "batch_prefix_stripped_match",
        TRUE ~ "group_parsed_from_sample_name"
      )
    )

  # Exclude samples explicitly marked DELETE in the full metadata.
  excluded_exact <- unique(meta_excluded$Sample_norm)
  excluded_strip <- unique(meta_excluded$Sample_strip_batch)
  sample_meta$excluded_by_metadata <-
    sample_meta$Sample_norm %in% excluded_exact |
    sample_meta$Sample_strip_batch %in% excluded_strip
  sample_meta <- sample_meta[!sample_meta$excluded_by_metadata, , drop = FALSE]

  # Apply only explicitly declared manual overrides.
  sample_meta$manual_override <- FALSE
  sample_meta$manual_override_to <- NA_character_

  if (!is.null(params$manual_group_overrides) && nrow(params$manual_group_overrides) > 0) {
    overrides <- params$manual_group_overrides
    overrides$sample_key_norm <- normalise_sample_id(overrides$sample_key)
    overrides$sample_key_strip <- strip_batch_prefix(overrides$sample_key)
    overrides$forced_group <- collapse_group(overrides$forced_group)

    for (i in seq_len(nrow(sample_meta))) {
      keys <- unique(c(
        sample_meta$Sample_norm[i],
        sample_meta$Sample_strip_batch[i],
        normalise_sample_id(sample_meta$Sample_raw[i]),
        strip_batch_prefix(sample_meta$Sample_raw[i])
      ))

      hits <- which(
        overrides$sample_key_norm %in% keys |
        overrides$sample_key_strip %in% keys
      )

      if (length(hits) > 0) {
        forced <- unique(overrides$forced_group[hits])
        forced <- forced[!is.na(forced)]
        if (length(forced) > 1) {
          stop("Conflicting manual group overrides for ", sample_meta$Sample_column[i])
        }
        if (length(forced) == 1) {
          sample_meta$Group_final[i] <- forced
          sample_meta$manual_override[i] <- TRUE
          sample_meta$manual_override_to[i] <- forced
        }
      }
    }
  }

  # Preserve genuine patient identifiers where available. Assign a unique ID to
  # samples lacking patient information so that optional mixed models do not
  # incorrectly link unrelated observations.
  sample_meta$Patient_ID <- trimws(as.character(sample_meta$Patient_ID))
  missing_patient <- is.na(sample_meta$Patient_ID) | !nzchar(sample_meta$Patient_ID)
  sample_meta$Patient_ID[missing_patient] <- paste0(
    "UNMATCHED_", sample_meta$Sample_norm[missing_patient]
  )

  sample_meta$Group_plot <- factor(sample_meta$Group_final, levels = params$group_order)

  if (!any(sample_meta$Group_final == "H_COL", na.rm = TRUE)) {
    stop("No H_COL samples were identified after metadata alignment.")
  }

  attr(sample_meta, "patient_column_detected") <- patient_col
  sample_meta
}

# =============================================================================
# 7. DIFFERENTIAL-EXPRESSION TABLE HELPERS
# =============================================================================

load_cpm_and_metadata <- function() {
  message("[5/9] Loading the CPM matrix and aligning metadata...")

  obj <- read_cpm_header_safe(params$cpm_file, gene_col_name = "gene_id")
  df <- obj$df
  write_csv_safe(obj$qc, file.path(params$out_dir, "qc", "CPM_header_detection.csv"))

  annotation_names <- normalise_name(c(
    "gene_id", "Length", "gene_length", "gene_length_bp", "exon_length",
    "effective_length", "gene_biotype", "biotype", "description",
    "chr", "chromosome", "start", "end", "strand"
  ))
  current_names <- normalise_name(names(df))
  annotation_hit <- current_names %in% annotation_names

  sample_cols <- names(df)[!annotation_hit]
  sample_cols <- base::setdiff(sample_cols, "gene_id")
  expr_df <- df[, sample_cols, drop = FALSE]
  expr_num <- as.data.frame(lapply(expr_df, to_num), check.names = FALSE)

  nonnumeric <- names(expr_num)[
    vapply(expr_num, function(x) all(is.na(x)), logical(1))
  ]
  if (length(nonnumeric) > 0) {
    warning("Dropping nonnumeric CPM columns: ", paste(nonnumeric, collapse = ", "))
    sample_cols <- base::setdiff(sample_cols, nonnumeric)
    expr_num <- expr_num[, sample_cols, drop = FALSE]
  }

  if (length(sample_cols) == 0) stop("No numeric CPM sample columns remained.")

  cpm_raw <- as.matrix(expr_num)
  storage.mode(cpm_raw) <- "double"
  rownames(cpm_raw) <- make.unique(as.character(df$gene_id))

  if (any(!is.finite(cpm_raw))) stop("The CPM matrix contains non-finite values.")
  if (any(cpm_raw < 0)) stop("The CPM matrix contains negative values.")
  if (anyDuplicated(colnames(cpm_raw))) {
    stop("The CPM matrix contains duplicated sample names.")
  }

  cpm_totals <- colSums(cpm_raw)
  cpm_qc <- data.frame(
    Sample_column = colnames(cpm_raw),
    total_CPM = as.numeric(cpm_totals),
    expected_if_complete_matrix = params$cpm_sum_expected,
    relative_deviation =
      (cpm_totals - params$cpm_sum_expected) / params$cpm_sum_expected,
    within_tolerance = abs(cpm_totals - params$cpm_sum_expected) <=
      params$cpm_sum_expected * params$cpm_sum_tolerance_fraction,
    stringsAsFactors = FALSE
  )
  write_csv_safe(cpm_qc, file.path(params$out_dir, "qc", "CPM_value_sanity_QC.csv"))

  gene_ann <- prepare_gene_symbols(df$gene_id)
  write_csv_safe(
    gene_ann,
    file.path(params$out_dir, "qc", "CPM_gene_ID_to_symbol_mapping.csv")
  )

  cpm_symbol <- collapse_expression_by_symbol(cpm_raw, gene_ann$gene_symbol)
  logcpm_symbol <- log2(cpm_symbol + 1)

  sample_meta <- load_and_align_metadata(
    params$metadata_file,
    colnames(cpm_symbol)
  )

  cpm_symbol <- cpm_symbol[, sample_meta$Sample_column, drop = FALSE]
  logcpm_symbol <- logcpm_symbol[, sample_meta$Sample_column, drop = FALSE]
  cpm_raw <- cpm_raw[, sample_meta$Sample_column, drop = FALSE]

  stopifnot(identical(colnames(cpm_symbol), sample_meta$Sample_column))

  write_csv_safe(
    sample_meta,
    file.path(params$out_dir, "qc", "metadata_aligned_to_CPM.csv")
  )
  write_csv_safe(
    sample_meta[sample_meta$manual_override %in% TRUE, , drop = FALSE],
    file.path(params$out_dir, "qc", "manual_group_overrides_applied.csv")
  )
  write_tsv_safe(
    data.frame(
      gene_symbol = rownames(cpm_symbol),
      as.data.frame(cpm_symbol, check.names = FALSE),
      check.names = FALSE
    ),
    file.path(params$out_dir, "tables", "CPM_gene_symbol_matrix.tsv")
  )
  write_tsv_safe(
    data.frame(
      gene_symbol = rownames(logcpm_symbol),
      as.data.frame(logcpm_symbol, check.names = FALSE),
      check.names = FALSE
    ),
    file.path(params$out_dir, "tables", "log2_CPM_plus_1_gene_symbol_matrix.tsv")
  )

  list(
    cpm_raw = cpm_raw,
    cpm_symbol = cpm_symbol,
    logcpm_symbol = logcpm_symbol,
    gene_ann = gene_ann,
    sample_meta = sample_meta,
    cpm_qc = cpm_qc
  )
}

read_gene_lengths <- function(path, gene_ann) {
  df <- read_table_flexible(path)

  gene_col <- pick_first_existing(
    c(
      "gene_id", "Geneid", "GeneID", "ENSEMBL", "ensembl_gene_id",
      "gene", "Gene", "gene_symbol", "SYMBOL", "symbol", "V1"
    ),
    df
  )
  length_col <- pick_first_existing(
    c(
      "gene_length_bp", "length_bp", "exon_length_bp", "effective_length_bp",
      "gene_length", "length", "Length", "exon_length", "effective_length", "V2"
    ),
    df
  )

  # Re-read a two-column, no-header file if named columns were not recognised.
  if (is.na(gene_col) || is.na(length_col)) {
    no_header <- tryCatch(
      fread(
        path, header = FALSE, data.table = FALSE,
        fill = TRUE, check.names = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(no_header) && ncol(no_header) >= 2) {
      candidate_lengths <- to_num(no_header[[2]])
      if (sum(is.finite(candidate_lengths) & candidate_lengths > 0) > 0) {
        df <- no_header
        gene_col <- names(df)[1]
        length_col <- names(df)[2]
      }
    }
  }

  if (is.na(gene_col) || is.na(length_col)) {
    stop("Gene and length columns could not be identified in: ", path)
  }

  keys <- toupper(trimws(as.character(df[[gene_col]])))
  keys <- sub("\\..*$", "", keys)
  lengths <- to_num(df[[length_col]])

  valid <- !is.na(keys) & nzchar(keys) & is.finite(lengths) & lengths > 0
  keys <- keys[valid]
  lengths <- lengths[valid]

  unit_inferred <- "bp"
  if (median(lengths, na.rm = TRUE) < 1000) {
    lengths <- lengths * 1000
    unit_inferred <- "kb_converted_to_bp"
  }

  lookup <- data.frame(
    key = keys,
    gene_length_bp = lengths,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(key) %>%
    dplyr::summarise(gene_length_bp = mean(gene_length_bp), .groups = "drop")

  length_map <- setNames(lookup$gene_length_bp, lookup$key)
  id_key <- toupper(gene_ann$gene_id_novers)
  symbol_key <- toupper(gene_ann$gene_symbol)

  matched <- unname(length_map[id_key])
  missing <- !is.finite(matched) | matched <= 0
  matched[missing] <- unname(length_map[symbol_key[missing]])

  qc <- data.frame(
    gene_id = gene_ann$gene_id,
    gene_id_novers = gene_ann$gene_id_novers,
    gene_symbol = gene_ann$gene_symbol,
    gene_length_bp = matched,
    usable_for_TPM = is.finite(matched) & matched > 0,
    source_gene_column = gene_col,
    source_length_column = length_col,
    length_unit_inferred = unit_inferred,
    stringsAsFactors = FALSE
  )

  qc
}

convert_cpm_to_tpm <- function(cpm_raw, gene_ann) {
  message("[6/9] Converting CPM to TPM for singscore...")

  length_qc <- read_gene_lengths(params$gene_length_file, gene_ann)
  write_csv_safe(
    length_qc,
    file.path(params$out_dir, "qc", "CPM_to_TPM_gene_length_mapping.csv")
  )

  valid <- length_qc$usable_for_TPM
  coverage <- mean(valid)

  coverage_table <- data.frame(
    n_genes_total = length(valid),
    n_genes_with_usable_length = sum(valid),
    usable_length_fraction = coverage,
    required_fraction = params$tpm_min_gene_length_coverage,
    stringsAsFactors = FALSE
  )
  write_csv_safe(
    coverage_table,
    file.path(params$out_dir, "qc", "CPM_to_TPM_gene_length_coverage.csv")
  )

  if (coverage < params$tpm_min_gene_length_coverage) {
    stop(
      "Gene-length coverage for CPM-to-TPM conversion was ",
      round(coverage, 3),
      ", below the required ", params$tpm_min_gene_length_coverage, "."
    )
  }

  cpm_valid <- cpm_raw[valid, , drop = FALSE]
  lengths_kb <- length_qc$gene_length_bp[valid] / 1000

  rpk <- sweep(cpm_valid, 1, lengths_kb, "/")
  denominators <- colSums(rpk)
  if (any(!is.finite(denominators) | denominators <= 0)) {
    stop("At least one sample had a non-positive CPM/length denominator.")
  }

  tpm_raw <- sweep(rpk, 2, denominators, "/") * 1e6
  rownames(tpm_raw) <- rownames(cpm_raw)[valid]
  storage.mode(tpm_raw) <- "double"

  gene_ann_valid <- gene_ann[valid, , drop = FALSE]
  tpm_symbol <- collapse_expression_by_symbol(
    tpm_raw,
    gene_ann_valid$gene_symbol
  )

  write_tsv_safe(
    data.frame(
      gene_symbol = rownames(tpm_symbol),
      as.data.frame(tpm_symbol, check.names = FALSE),
      check.names = FALSE
    ),
    file.path(params$out_dir, "tables", "TPM_gene_symbol_matrix_for_singscore.tsv")
  )

  raw_tpm_qc <- data.frame(
    Sample_column = colnames(tpm_raw),
    TPM_sum_before_symbol_collapse = colSums(tpm_raw),
    stringsAsFactors = FALSE
  )
  write_csv_safe(
    raw_tpm_qc,
    file.path(params$out_dir, "qc", "TPM_column_sums_before_symbol_collapse.csv")
  )

  list(
    tpm_raw = tpm_raw,
    tpm_symbol = tpm_symbol,
    length_qc = length_qc,
    coverage = coverage
  )
}

# =============================================================================
# 13. HALLMARK SIGNATURES AND SINGSCORE
# =============================================================================

filter_hallmark_by_mean_logcpm <- function(gene_set, signature_name, logcpm_mat) {
  stats <- data.frame(
    Gene = rownames(logcpm_mat),
    mean_logCPM = rowMeans(logcpm_mat, na.rm = TRUE),
    median_logCPM = apply(logcpm_mat, 1, median, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  audit <- data.frame(Gene = gene_set, stringsAsFactors = FALSE) %>%
    dplyr::left_join(stats, by = "Gene") %>%
    dplyr::mutate(
      present_in_matrix = Gene %in% rownames(logcpm_mat),
      retained = present_in_matrix & is.finite(mean_logCPM) &
        mean_logCPM >= params$hallmark_mean_logcpm_cutoff
    ) %>%
    dplyr::arrange(dplyr::desc(retained), dplyr::desc(mean_logCPM), Gene)

  retained <- audit %>% dplyr::filter(retained) %>% dplyr::pull(Gene) %>% unique()
  if (length(retained) == 0) {
    stop(signature_name, " has no genes after mean logCPM filtering.")
  }

  write_csv_safe(
    audit,
    file.path(
      params$out_dir, "tables",
      paste0(signature_name, "_mean_logCPM_filter_audit.csv")
    )
  )
  write_lines_safe(
    retained,
    file.path(
      params$out_dir, "tables",
      paste0(signature_name, "_genes_used.txt")
    )
  )

  retained
}

score_one_signature <- function(rank_data, tpm_mat, genes, signature_name, sample_meta) {
  genes <- standardise_gene_vector(genes)
  genes_present <- base::intersect(genes, rownames(tpm_mat))

  if (length(genes_present) == 0) {
    stop("No genes from ", signature_name, " were present in the TPM matrix.")
  }

  centered <- singscore::simpleScore(
    rankData = rank_data,
    upSet = genes_present,
    knownDirection = TRUE,
    centerScore = TRUE
  )
  uncentered <- singscore::simpleScore(
    rankData = rank_data,
    upSet = genes_present,
    knownDirection = TRUE,
    centerScore = FALSE
  )

  dispersion <- if ("TotalDispersion" %in% names(centered)) {
    as.numeric(centered$TotalDispersion)
  } else if ("Dispersion" %in% names(centered)) {
    as.numeric(centered$Dispersion)
  } else {
    rep(NA_real_, length(centered$TotalScore))
  }

  out <- data.frame(
    Sample_column = colnames(tpm_mat),
    signature = signature_name,
    centered_singscore = as.numeric(centered$TotalScore),
    uncentered_singscore = as.numeric(uncentered$TotalScore),
    rank_dispersion = dispersion,
    n_signature_intended = length(genes),
    n_signature_present = length(genes_present),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(sample_meta, by = "Sample_column")

  hcol_median <- median(
    out$centered_singscore[out$Group_final == "H_COL"],
    na.rm = TRUE
  )
  out$delta_HCOL_median <- out$centered_singscore - hcol_median
  out$Group_plot <- factor(out$Group_final, levels = params$group_order)

  list(scores = out, genes_present = genes_present)
}

calculate_all_singscores <- function(final_genes, expression_obj, tpm_obj) {
  message("[7/9] Calculating sample-level metaplasia and inflammatory singscores...")

  tpm_mat <- tpm_obj$tpm_symbol[, expression_obj$sample_meta$Sample_column, drop = FALSE]
  rank_data <- singscore::rankGenes(tpm_mat, tiesMethod = "min")

  final_present <- base::intersect(final_genes, rownames(tpm_mat))
  final_missing <- base::setdiff(final_genes, rownames(tpm_mat))

  write_lines_safe(
    final_present,
    file.path(params$out_dir, "tables", "FINAL_SIGNATURE_genes_present_in_TPM.txt")
  )
  write_lines_safe(
    final_missing,
    file.path(params$out_dir, "qc", "FINAL_SIGNATURE_genes_missing_from_TPM.txt")
  )

  if (
    params$require_all_final_signature_genes_measured &&
      length(final_missing) > 0
  ) {
    stop(
      "Not all 24 final-signature genes were present in the TPM matrix. ",
      "Inspect FINAL_SIGNATURE_genes_missing_from_TPM.txt."
    )
  }

  hallmark_all <- load_hallmark_gene_sets()
  hallmark_filtered <- list(
    INTERFERON_ALPHA_RESPONSE = filter_hallmark_by_mean_logcpm(
      hallmark_all$INTERFERON_ALPHA_RESPONSE,
      "HALLMARK_INTERFERON_ALPHA_RESPONSE",
      expression_obj$logcpm_symbol
    ),
    INFLAMMATORY_RESPONSE = filter_hallmark_by_mean_logcpm(
      hallmark_all$INFLAMMATORY_RESPONSE,
      "HALLMARK_INFLAMMATORY_RESPONSE",
      expression_obj$logcpm_symbol
    ),
    TNFA_SIGNALING_VIA_NFKB = filter_hallmark_by_mean_logcpm(
      hallmark_all$TNFA_SIGNALING_VIA_NFKB,
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      expression_obj$logcpm_symbol
    )
  )

  signature_list <- c(
    list(GASTRIC_METAPLASIA_24_GENE = final_genes),
    hallmark_filtered
  )

  scored <- lapply(names(signature_list), function(sig) {
    score_one_signature(
      rank_data = rank_data,
      tpm_mat = tpm_mat,
      genes = signature_list[[sig]],
      signature_name = sig,
      sample_meta = expression_obj$sample_meta
    )
  })
  names(scored) <- names(signature_list)

  score_long <- dplyr::bind_rows(lapply(scored, `[[`, "scores"))

  score_wide <- score_long %>%
    dplyr::select(
      Sample_column, Sample_norm, Patient_ID, Group_final,
      signature, centered_singscore
    ) %>%
    tidyr::pivot_wider(
      names_from = signature,
      values_from = centered_singscore
    )

  coverage <- data.frame(
    signature = names(scored),
    n_intended = vapply(signature_list, length, integer(1)),
    n_present_in_TPM = vapply(
      scored,
      function(x) length(x$genes_present),
      integer(1)
    ),
    final_signature_mean_logCPM_filter_applied = c(
      FALSE,
      rep(TRUE, length(scored) - 1)
    ),
    stringsAsFactors = FALSE
  )

  write_csv_safe(
    score_long,
    file.path(params$out_dir, "tables", "ALL_SIGNATURES_singscore_long.csv")
  )
  write_csv_safe(
    score_wide,
    file.path(params$out_dir, "tables", "ALL_SIGNATURES_singscore_wide.csv")
  )
  write_csv_safe(
    coverage,
    file.path(params$out_dir, "qc", "ALL_SIGNATURES_singscore_gene_coverage.csv")
  )

  list(
    long = score_long,
    wide = score_wide,
    scored = scored,
    signature_list = signature_list,
    coverage = coverage
  )
}


# =============================================================================
# 13B. COMPLETE SCORE SUITE: Z-SCORE, SINGSCORE AND ssGSEA
# =============================================================================

compute_hcol_gene_z_matrix <- function(logcpm_mat, sample_meta) {
  hcol_samples <- sample_meta$Sample_column[sample_meta$Group_final == "H_COL"]
  hcol_samples <- base::intersect(hcol_samples, colnames(logcpm_mat))
  if (length(hcol_samples) < 2) {
    stop("At least two H_COL samples are required for gene-wise z-score calculation.")
  }

  ref <- logcpm_mat[, hcol_samples, drop = FALSE]
  ref_mean <- rowMeans(ref, na.rm = TRUE)
  ref_sd <- apply(ref, 1, stats::sd, na.rm = TRUE)
  usable <- is.finite(ref_mean) & is.finite(ref_sd) & ref_sd > 0

  z_mat <- matrix(
    NA_real_, nrow = nrow(logcpm_mat), ncol = ncol(logcpm_mat),
    dimnames = dimnames(logcpm_mat)
  )
  if (any(usable)) {
    tmp <- sweep(logcpm_mat[usable, , drop = FALSE], 1, ref_mean[usable], "-")
    tmp <- sweep(tmp, 1, ref_sd[usable], "/")
    z_mat[usable, ] <- tmp
  }

  gene_stats <- data.frame(
    gene_symbol = rownames(logcpm_mat),
    H_COL_reference_mean_log2CPM_plus_1 = ref_mean,
    H_COL_reference_sd_log2CPM_plus_1 = ref_sd,
    usable_for_z_score = usable,
    stringsAsFactors = FALSE
  )
  write_csv_safe(
    gene_stats,
    file.path(params$out_dir, "qc", "GENEWISE_ZSCORE_HCOL_REFERENCE_STATS.csv")
  )
  write_tsv_safe(
    data.frame(
      gene_symbol = rownames(z_mat),
      as.data.frame(z_mat, check.names = FALSE),
      check.names = FALSE
    ),
    file.path(params$out_dir, "tables", "GENEWISE_ZSCORE_HCOL_MATRIX.tsv")
  )

  list(z_mat = z_mat, gene_stats = gene_stats, hcol_samples = hcol_samples)
}

score_one_z_signature <- function(z_mat, genes, signature_name, sample_meta) {
  genes <- standardise_gene_vector(genes)
  genes_present <- base::intersect(genes, rownames(z_mat))
  if (length(genes_present) > 0) {
    genes_present <- genes_present[
      rowSums(is.finite(z_mat[genes_present, , drop = FALSE])) > 0
    ]
  }
  if (length(genes_present) == 0) {
    stop("No genes from ", signature_name, " were usable in the z-score matrix.")
  }

  sub <- z_mat[genes_present, , drop = FALSE]
  n_used <- colSums(is.finite(sub))
  if (any(n_used == 0)) {
    stop("At least one sample had no usable genes for z-score signature ", signature_name, ".")
  }

  score <- colSums(sub, na.rm = TRUE) / n_used
  stouffer <- colSums(sub, na.rm = TRUE) / sqrt(n_used)

  out <- data.frame(
    Sample_column = colnames(z_mat),
    signature = signature_name,
    method = "Zscore",
    score = as.numeric(score),
    score_stouffer = as.numeric(stouffer),
    n_signature_intended = length(genes),
    n_signature_present = length(genes_present),
    n_signature_used = as.numeric(n_used),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(sample_meta, by = "Sample_column")

  hcol_median <- median(out$score[out$Group_final == "H_COL"], na.rm = TRUE)
  out$delta_HCOL_median <- out$score - hcol_median
  out$Group_plot <- factor(out$Group_final, levels = params$group_order)

  list(scores = out, genes_present = genes_present)
}

calculate_all_zscores <- function(signature_list, expression_obj) {
  z_obj <- compute_hcol_gene_z_matrix(
    expression_obj$logcpm_symbol,
    expression_obj$sample_meta
  )

  scored <- lapply(names(signature_list), function(sig) {
    score_one_z_signature(
      z_mat = z_obj$z_mat,
      genes = signature_list[[sig]],
      signature_name = sig,
      sample_meta = expression_obj$sample_meta
    )
  })
  names(scored) <- names(signature_list)

  long <- dplyr::bind_rows(lapply(scored, `[[`, "scores"))
  wide <- long %>%
    dplyr::select(
      Sample_column, Sample_norm, Patient_ID, Group_final,
      signature, score
    ) %>%
    tidyr::pivot_wider(names_from = signature, values_from = score)

  write_csv_safe(long, file.path(params$out_dir, "tables", "ALL_SIGNATURES_Zscore_long.csv"))
  write_csv_safe(wide, file.path(params$out_dir, "tables", "ALL_SIGNATURES_Zscore_wide.csv"))

  list(long = long, wide = wide, scored = scored, z_obj = z_obj)
}

run_ssgsea_gsva_safe <- function(expr_mat, gene_sets) {
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "double"

  gene_sets <- lapply(gene_sets, standardise_gene_vector)
  gene_sets <- lapply(gene_sets, function(g) base::intersect(g, rownames(expr_mat)))
  empty <- names(gene_sets)[vapply(gene_sets, length, integer(1)) == 0]
  if (length(empty) > 0) {
    stop("The following signatures had no genes in the ssGSEA matrix: ", paste(empty, collapse = ", "))
  }

  max_size_effective <- if (is.infinite(params$ssgsea_max_genes)) {
    nrow(expr_mat)
  } else {
    min(params$ssgsea_max_genes, nrow(expr_mat))
  }

  out <- NULL
  errors <- character(0)

  if (exists("ssgseaParam", where = asNamespace("GSVA"), mode = "function")) {
    attempts <- list(
      function() GSVA::ssgseaParam(
        exprData = expr_mat,
        geneSets = gene_sets,
        minSize = 1,
        maxSize = max_size_effective,
        alpha = params$ssgsea_alpha,
        normalize = params$ssgsea_norm
      ),
      function() GSVA::ssgseaParam(
        exprData = expr_mat,
        geneSets = gene_sets,
        minSize = 1,
        maxSize = max_size_effective,
        normalize = params$ssgsea_norm
      ),
      function() GSVA::ssgseaParam(
        exprData = expr_mat,
        geneSets = gene_sets,
        minSize = 1,
        maxSize = max_size_effective
      )
    )

    for (make_param in attempts) {
      candidate <- tryCatch({
        param <- make_param()
        GSVA::gsva(param, verbose = FALSE)
      }, error = function(e) {
        errors <<- c(errors, paste0("new_api: ", conditionMessage(e)))
        NULL
      })
      if (!is.null(candidate)) {
        out <- candidate
        break
      }
    }
  }

  if (is.null(out)) {
    out <- tryCatch(
      GSVA::gsva(
        expr_mat,
        gene_sets,
        method = "ssgsea",
        kcdf = "Gaussian",
        min.sz = 1,
        max.sz = max_size_effective,
        abs.ranking = FALSE,
        ssgsea.norm = params$ssgsea_norm,
        verbose = FALSE
      ),
      error = function(e) {
        errors <<- c(errors, paste0("old_api: ", conditionMessage(e)))
        NULL
      }
    )
  }

  if (is.null(out)) {
    stop("ssGSEA failed with both GSVA APIs: ", paste(unique(errors), collapse = " | "))
  }

  out <- as.matrix(out)
  if (!all(colnames(expr_mat) %in% colnames(out)) &&
      all(colnames(expr_mat) %in% rownames(out))) {
    out <- t(out)
  }
  if (!all(colnames(expr_mat) %in% colnames(out))) {
    stop("ssGSEA output sample names did not match the input expression matrix.")
  }
  out[, colnames(expr_mat), drop = FALSE]
}

calculate_all_ssgsea <- function(signature_list, expression_obj) {
  matrix_scores <- run_ssgsea_gsva_safe(
    expression_obj$logcpm_symbol,
    signature_list
  )
  matrix_scores <- matrix_scores[names(signature_list), , drop = FALSE]

  scored <- lapply(names(signature_list), function(sig) {
    genes <- standardise_gene_vector(signature_list[[sig]])
    genes_present <- base::intersect(genes, rownames(expression_obj$logcpm_symbol))
    out <- data.frame(
      Sample_column = colnames(expression_obj$logcpm_symbol),
      signature = sig,
      method = "SSGSEA",
      score = as.numeric(matrix_scores[sig, ]),
      n_signature_intended = length(genes),
      n_signature_present = length(genes_present),
      n_signature_used = length(genes_present),
      ssgsea_norm = params$ssgsea_norm,
      ssgsea_alpha = params$ssgsea_alpha,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::left_join(expression_obj$sample_meta, by = "Sample_column")

    hcol_median <- median(out$score[out$Group_final == "H_COL"], na.rm = TRUE)
    out$delta_HCOL_median <- out$score - hcol_median
    out$Group_plot <- factor(out$Group_final, levels = params$group_order)
    list(scores = out, genes_present = genes_present)
  })
  names(scored) <- names(signature_list)

  long <- dplyr::bind_rows(lapply(scored, `[[`, "scores"))
  wide <- long %>%
    dplyr::select(
      Sample_column, Sample_norm, Patient_ID, Group_final,
      signature, score
    ) %>%
    tidyr::pivot_wider(names_from = signature, values_from = score)

  write_csv_safe(
    data.frame(
      signature = rownames(matrix_scores),
      as.data.frame(matrix_scores, check.names = FALSE),
      check.names = FALSE
    ),
    file.path(params$out_dir, "tables", "ALL_SIGNATURES_SSGSEA_score_matrix.csv")
  )
  write_csv_safe(long, file.path(params$out_dir, "tables", "ALL_SIGNATURES_SSGSEA_long.csv"))
  write_csv_safe(wide, file.path(params$out_dir, "tables", "ALL_SIGNATURES_SSGSEA_wide.csv"))

  list(long = long, wide = wide, scored = scored, score_matrix = matrix_scores)
}

calculate_complete_score_suite <- function(final_genes, expression_obj, tpm_obj) {
  message("[7/9] Calculating z-score, singscore and ssGSEA for all four signatures...")

  singscore_obj <- calculate_all_singscores(final_genes, expression_obj, tpm_obj)
  signature_list <- singscore_obj$signature_list
  zscore_obj <- calculate_all_zscores(signature_list, expression_obj)
  ssgsea_obj <- calculate_all_ssgsea(signature_list, expression_obj)

  singscore_long <- singscore_obj$long %>%
    dplyr::transmute(
      Sample_column, Sample_norm, Patient_ID, Group_final, Group_plot,
      signature,
      method = "Singscore",
      score = centered_singscore,
      auxiliary_score = uncentered_singscore,
      rank_dispersion,
      n_signature_intended,
      n_signature_present,
      n_signature_used = n_signature_present,
      delta_HCOL_median
    )

  zscore_long <- zscore_obj$long %>%
    dplyr::transmute(
      Sample_column, Sample_norm, Patient_ID, Group_final, Group_plot,
      signature, method, score,
      auxiliary_score = score_stouffer,
      rank_dispersion = NA_real_,
      n_signature_intended,
      n_signature_present,
      n_signature_used,
      delta_HCOL_median
    )

  ssgsea_long <- ssgsea_obj$long %>%
    dplyr::transmute(
      Sample_column, Sample_norm, Patient_ID, Group_final, Group_plot,
      signature, method, score,
      auxiliary_score = NA_real_,
      rank_dispersion = NA_real_,
      n_signature_intended,
      n_signature_present,
      n_signature_used,
      delta_HCOL_median
    )

  all_long <- dplyr::bind_rows(zscore_long, singscore_long, ssgsea_long)
  all_wide <- all_long %>%
    dplyr::mutate(score_key = paste(method, signature, sep = "__")) %>%
    dplyr::select(
      Sample_column, Sample_norm, Patient_ID, Group_final,
      score_key, score
    ) %>%
    tidyr::pivot_wider(names_from = score_key, values_from = score)

  coverage <- dplyr::bind_rows(
    data.frame(
      method = "Zscore",
      signature = names(zscore_obj$scored),
      n_intended = vapply(signature_list, length, integer(1)),
      n_present = vapply(zscore_obj$scored, function(x) length(x$genes_present), integer(1)),
      stringsAsFactors = FALSE
    ),
    data.frame(
      method = "Singscore",
      signature = names(singscore_obj$scored),
      n_intended = vapply(signature_list, length, integer(1)),
      n_present = vapply(singscore_obj$scored, function(x) length(x$genes_present), integer(1)),
      stringsAsFactors = FALSE
    ),
    data.frame(
      method = "SSGSEA",
      signature = names(ssgsea_obj$scored),
      n_intended = vapply(signature_list, length, integer(1)),
      n_present = vapply(ssgsea_obj$scored, function(x) length(x$genes_present), integer(1)),
      stringsAsFactors = FALSE
    )
  ) %>%
    dplyr::mutate(
      mean_logCPM_filter_applied = signature != "GASTRIC_METAPLASIA_24_GENE",
      mean_logCPM_cutoff = ifelse(
        mean_logCPM_filter_applied,
        params$hallmark_mean_logcpm_cutoff,
        NA_real_
      )
    )

  write_csv_safe(
    all_long,
    file.path(params$out_dir, "tables", "ALL_SIGNATURES_ALL_METHODS_long.csv")
  )
  write_csv_safe(
    all_wide,
    file.path(params$out_dir, "tables", "ALL_SIGNATURES_ALL_METHODS_wide.csv")
  )
  write_csv_safe(
    coverage,
    file.path(params$out_dir, "qc", "ALL_SIGNATURES_ALL_METHODS_gene_coverage.csv")
  )

  # One convenient table per signature containing all three score methods.
  for (signature_name in names(signature_list)) {
    cols <- c(
      "Sample_column", "Sample_norm", "Patient_ID", "Group_final",
      paste0("Zscore__", signature_name),
      paste0("Singscore__", signature_name),
      paste0("SSGSEA__", signature_name)
    )
    missing_cols <- base::setdiff(cols, names(all_wide))
    if (length(missing_cols) > 0) {
      stop(
        "Missing score columns while exporting combined table for ",
        signature_name, ": ", paste(missing_cols, collapse = ", ")
      )
    }
    write_csv_safe(
      all_wide[, cols, drop = FALSE],
      file.path(
        params$out_dir, "tables",
        paste0(signature_name, "__combined_Zscore_Singscore_SSGSEA.csv")
      )
    )
  }

  # Top-level aliases preserve compatibility with the established primary
  # singscore figures and correlations in the remainder of the pipeline.
  list(
    long = singscore_obj$long,
    wide = singscore_obj$wide,
    scored = singscore_obj$scored,
    signature_list = signature_list,
    coverage = coverage,
    singscore = singscore_obj,
    zscore = zscore_obj,
    ssgsea = ssgsea_obj,
    all_long = all_long,
    all_wide = all_wide
  )
}

 },envir=e)
 e
}

# Clustered correlation computation.
# Patient-block studentized bootstrap for a biopsy-weighted Spearman correlation.
# Null: population marginal Spearman rho = 0. No covariate adjustment.
# Every sampled patient contributes all their biopsies; ranks are recalculated.
# A delete-one-patient jackknife estimates the SE of Fisher's atanh(rho),
# both in the original sample and each bootstrap sample. Repeated bootstrap
# copies of a patient are distinct resampling units, not merged into one unit.
make_correlation_engine <- function(d, xcols, ycols, correlation='spearman') {
  z<-as.matrix(d[,unique(c(xcols,ycols)),drop=FALSE]);storage.mode(z)<-'double'
  if(any(!is.finite(z)) || anyNA(d$Patient_ID)) stop('Incomplete correlation inputs.')
  ids<-sort(unique(as.character(d$Patient_ID)))
  cl<-match(as.character(d$Patient_ID),ids);G<-length(ids);n<-nrow(z)
  pairs<-expand.grid(x_column=xcols,y_column=ycols,stringsAsFactors=FALSE)
  xi<-match(pairs$x_column,colnames(z));yi<-match(pairs$y_column,colnames(z))
  # Midrank kernel: rows are observed values; columns are reference values.
  # Integer row weights exactly reproduce the expanded resampled data,
  # including average ranks for ties introduced by resampling.
  if(correlation=='spearman') {
    K<-lapply(seq_len(ncol(z)),function(j)
      outer(z[,j],z[,j],function(a,b) (a>b)+0.5*(a==b)))
  }
  point<-function(mult) {
    w<-mult[cl];tot<-sum(w)
    if(tot<3) return(rep(NA_real_,nrow(pairs)))
    W<-w/tot
    zz<-if(correlation=='spearman') vapply(K,function(k) as.vector(k%*%W),numeric(n)) else z
    zz<-sweep(zz,2,colSums(zz*W),'-')
    cv<-crossprod(zz,zz*W);vv<-diag(cv)
    r<-cv[cbind(xi,yi)]/sqrt(vv[xi]*vv[yi])
    r[!is.finite(r)|abs(r)>1+1e-10]<-NA_real_
    pmin(1,pmax(-1,r))
  }
  fisher<-function(r) {r[abs(r)>=1-1e-12]<-NA_real_;atanh(r)}
  summary<-function(mult) {
    rho<-point(mult);zz<-fisher(rho);active<-which(mult>0);total_units<-sum(mult)
    dels<-vapply(active,function(g) {m<-mult;m[g]<-m[g]-1L;fisher(point(m))},numeric(length(zz)))
    if(is.null(dim(dels))) dels<-matrix(dels,nrow=length(zz))
    mu<-as.vector(dels%*%mult[active]/total_units)
    v<-(total_units-1)/total_units*rowSums(sweep((dels-mu)^2,2,mult[active],'*'))
    list(rho=rho,z=zz,se=sqrt(v))
  }
  list(point=point,summary=summary,G=G,n=n,pairs=pairs,ids=ids,cluster=cl)
}

cluster_bootstrap_correlations <- function(d,xcols,ycols,B=9999L,seed=1L,correlation='spearman') {
  engine<-make_correlation_engine(d,xcols,ycols,correlation)
  G<-engine$G;if(G<10L) stop('Fewer than 10 independent patients in a requested correlation subset.')
  observed<-engine$summary(rep(1L,G));P<-length(observed$rho)
  if(any(!is.finite(observed$se)|observed$se<=0)) stop('Correlation or patient-jackknife SE is not estimable.')
  set.seed(seed)
  multiplicities<-rmultinom(B,G,rep(1/G,G))
  roots<-matrix(NA_real_,B,P);rho_boot<-matrix(NA_real_,B,P)
  for(b in seq_len(B)) {
    z<-engine$summary(multiplicities[,b])
    roots[b,]<-(z$z-observed$z)/z$se;rho_boot[b,]<-z$rho
    if(b%%1000L==0L) message('  ',correlation,': ',b,'/',B,' patient-block bootstrap samples')
  }
  valid<-colSums(is.finite(roots));failed<-B-valid
  # Do not quietly replace invalid replicates or switch inferential methods.
  if(any(valid<0.99*B)) stop('More than 1% of bootstrap replicates are non-estimable; inspect inputs.')
  t0<-observed$z/observed$se
  p<-vapply(seq_len(P),function(j) {
    v<-roots[is.finite(roots[,j]),j]
    (1+sum(abs(v)>=abs(t0[j])))/(1+length(v))
  },numeric(1))
  crit<-vapply(seq_len(P),function(j) quantile(abs(roots[,j]),.95,na.rm=TRUE,type=1),numeric(1))
  tab<-cbind(engine$pairs,data.frame(Correlation=correlation,Coefficient=observed$rho,
    CI_low=tanh(observed$z-crit*observed$se),CI_high=tanh(observed$z+crit*observed$se),
    Fisher_SE_patient_jackknife=observed$se,Test_statistic=t0,P=p,
    Bootstrap_replicates=B,Valid_replicates=valid,Invalid_replicates=failed,
    Seed=seed,N_biopsies=engine$n,N_patients=G,
    Estimand='Marginal biopsy-weighted correlation',
    Test='Two-sided patient-block studentized bootstrap; Fisher transform',
    Tissue_group_adjusted=FALSE,stringsAsFactors=FALSE))
  list(table=tab,roots=roots,bootstrap_coefficients=rho_boot,multiplicities=multiplicities,
    observed=observed,patients=engine$ids,pairs=engine$pairs)
}

# CPM-to-score entry point using the original scoring functions above.
recompute_original_scores <- function(reference) {
  needed<-c(CPM_FILE,GENE_LENGTH_FILE,METADATA_FILE)
  if(any(!file.exists(needed))) stop('CPM-to-score mode requires these original files: ',
    paste(needed[!file.exists(needed)],collapse='; '),
    '. This script never silently substitutes saved scores. To explicitly rerun from the private verified score tables, set SCORE_INPUT_MODE to saved_scores.')
  cran<-c('data.table');bioc<-c('AnnotationDbi','org.Hs.eg.db','singscore','GSVA')
  missing<-cran[!vapply(cran,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)&&INSTALL_MISSING) install.packages(missing,repos='https://cloud.r-project.org')
  missing<-bioc[!vapply(bioc,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)&&INSTALL_MISSING) {
    if(!requireNamespace('BiocManager',quietly=TRUE)) install.packages('BiocManager',repos='https://cloud.r-project.org')
    BiocManager::install(missing,ask=FALSE,update=FALSE)
  }
  absent<-c(cran,bioc)[!vapply(c(cran,bioc),requireNamespace,logical(1),quietly=TRUE)]
  if(length(absent)) stop('Missing scoring packages: ',paste(absent,collapse=', '))
  e<-make_original_scoring_environment()
  e$params$cpm_file<-CPM_FILE;e$params$gene_length_file<-GENE_LENGTH_FILE
  e$params$metadata_file<-METADATA_FILE;e$params$out_dir<-file.path(.rerun$out,'scoring')
  for(p in c('','qc','tables')) dir.create(file.path(e$params$out_dir,p),recursive=TRUE,showWarnings=FALSE)
  e$fread<-data.table::fread;e$org.Hs.eg.db<-org.Hs.eg.db::org.Hs.eg.db
  # Freeze the original full Hallmark memberships; do not fetch a newer MSigDB.
  audits<-list(INTERFERON_ALPHA_RESPONSE=read_embedded(ifna_filter_csv),
    INFLAMMATORY_RESPONSE=read_embedded(inflam_filter_csv),TNFA_SIGNALING_VIA_NFKB=read_embedded(tnf_filter_csv))
  e$load_hallmark_gene_sets<-local({a<-audits;function() lapply(a,function(z) unique(as.character(z$Gene)))})
  # Keep alpha and normalisation explicit; never drop an argument on an API error.
  e$run_ssgsea_gsva_safe<-function(expr_mat,gene_sets) {
    if(!exists('ssgseaParam',envir=asNamespace('GSVA'),inherits=FALSE))
      stop('This scoring entry point requires GSVA with ssgseaParam. Do not silently substitute a different API.')
    par<-GSVA::ssgseaParam(exprData=as.matrix(expr_mat),geneSets=gene_sets,
      minSize=1,maxSize=nrow(expr_mat),alpha=.25,normalize=TRUE)
    out<-as.matrix(GSVA::gsva(par,verbose=FALSE))
    if(!all(colnames(expr_mat)%in%colnames(out))) stop('ssGSEA sample mismatch.')
    out[,colnames(expr_mat),drop=FALSE]
  }
  ex<-e$load_cpm_and_metadata()
  if(!setequal(ex$sample_meta$Sample_column,reference$Sample_column)||nrow(ex$sample_meta)!=79L)
    stop('CPM scoring cohort does not match the original 79 biopsies. No cohort change was accepted.')
  md<-ex$sample_meta[match(reference$Sample_column,ex$sample_meta$Sample_column),]
  if(any(as.character(md$Patient_ID)!=reference$Patient_ID)||any(md$Group_final!=reference$Group_final))
    stop('Rebuilt patient or tissue metadata differ from the verified reference.')
  tpm<-e$convert_cpm_to_tpm(ex$cpm_raw,ex$gene_ann)
  genes<-trimws(strsplit(gastric_genes_txt,'\n')[[1]]);genes<-genes[nzchar(genes)]
  scored<-e$calculate_complete_score_suite(genes,ex,tpm)
  out<-as.data.frame(scored$all_wide)
  out<-out[match(reference$Sample_column,out$Sample_column),names(reference),drop=FALSE]
  out$Patient_ID<-as.character(out$Patient_ID)
  keys<-grep('__',names(reference),value=TRUE)
  differences<-vapply(keys,function(k) max(abs(out[[k]]-reference[[k]])),numeric(1))
  write_table(data.frame(Score=keys,Maximum_absolute_difference=differences,Tolerance=1e-8),
    'qc/RECOMPUTED_VS_ORIGINAL_SCORES.csv')
  for(s in names(audits)) {
    expected<-as.character(audits[[s]]$Gene[audits[[s]]$retained])
    if(!setequal(scored$signature_list[[s]],expected)) stop('The retained gene list changed: ',s)
  }
  if(any(!is.finite(differences)|differences>1e-8))
    stop('Recomputed scores differ from the original beyond tolerance. Inspect scoring/mapping/package versions; no changed figures were accepted.')
  write_table(data.frame(Package=c(cran,bioc),Version=vapply(c(cran,bioc),function(p)as.character(packageVersion(p)),character(1))),
    'qc/SCORING_PACKAGE_VERSIONS.csv')
  write_table(data.frame(Path=needed,MD5=unname(tools::md5sum(needed))),'qc/RAW_INPUT_CHECKSUMS.csv')
  write_table(out,'inputs/RECOMPUTED_ALL_METHODS_SCORES.csv')
  out
}

validate_saved_inputs<-validate_inputs
validate_inputs<-function() {
  dat<-validate_saved_inputs()
  if(!SCORE_INPUT_MODE%in%c('from_cpm','saved_scores')) stop('SCORE_INPUT_MODE must be from_cpm or saved_scores.')
  if(SCORE_INPUT_MODE=='from_cpm') dat$wide<-recompute_original_scores(dat$wide)
  write_note(c(paste('Score input mode:',SCORE_INPUT_MODE),
    if(SCORE_INPUT_MODE=='from_cpm') 'Scores recalculated from CPM and matched against the original reference.' else
      'The exact verified saved scores are used. CPM-to-TPM and scoring are not rerun in this mode.'), 'SCORE_INPUT_MODE.txt')
  dat
}

run_associations<-function(wide) {
  # Retain the old tissue-adjusted mixed slopes only as a separately labelled
  # model analysis. Their tests are never used as Spearman P values.
  model_rows<-run_model_associations(wide)
  rows<-list()
  for(s in seq_along(SUBSETS)) {
    subset<-names(SUBSETS)[s]
    d<-wide[wide$Group_final%in%SUBSETS[[subset]],,drop=FALSE]
    # Stable order makes the seeded resampling independent of the input row order.
    d<-d[order(as.character(d$Patient_ID),d$Sample_column),,drop=FALSE]
    xc<-vapply(METHODS,function(m)score_column(m,SIGNATURES[1]),character(1))
    yc<-unlist(lapply(METHODS,function(m)vapply(SIGNATURES[-1],function(g)score_column(m,g),character(1))))
    for(cm in c('spearman','pearson')) {
      message('Patient-block ',cm,' analysis: ',subset)
      obj<-cluster_bootstrap_correlations(d,xc,yc,B=BOOTSTRAP_REPLICATES,
        seed=SET_SEED+1000L*s+if(cm=='spearman')0L else 500L,correlation=cm)
      tab<-obj$table;tab$subset<-subset
      tab$gastric_method<-sub('__.*','',tab$x_column)
      tab$inflammatory_method<-sub('__.*','',tab$y_column)
      tab$signature<-sub('^[^_]+__','',tab$y_column)
      rows[[paste(subset,cm)]]<-tab
      saveRDS(obj,file.path(.rerun$out,'models',paste0('BOOTSTRAP_',subset,'_',cm,'.rds')),compress='xz')
    }
  }
  all<-bind_rows(rows) %>% group_by(subset,gastric_method,inflammatory_method,Correlation) %>%
    mutate(P_BH_three=bh_family(P,3L)) %>% ungroup() %>% group_by(subset,Correlation) %>%
    mutate(P_BH_all27_sensitivity=bh_family(P,27L)) %>% ungroup()
  write_table(all,'tables/PATIENT_CLUSTERED_CORRELATIONS_ALL_METHODS.csv')
  primary<-filter(all,gastric_method=='Singscore',inflammatory_method=='Singscore',Correlation=='spearman')
  write_table(primary,'tables/PRIMARY_SINGSCORE_SPEARMAN.csv')
  write_table(filter(all,Correlation=='pearson'),'sensitivity/PATIENT_CLUSTERED_PEARSON.csv')
  sp<-filter(all,Correlation=='spearman')
  # Preserve names used by the plotting loop, while carrying the matching
  # Spearman inference explicitly in new, unambiguous fields.
  keycols<-c('subset','gastric_method','inflammatory_method','signature','x_column','y_column')
  model_rows<-model_rows %>% rename(Model_P_BH_three=P_BH_three)
  sp<-sp %>% rename(Spearman_rho=Coefficient,Spearman_P=P,Spearman_q=P_BH_three,
    Spearman_CI_low=CI_low,Spearman_CI_high=CI_high)
  joined<-left_join(model_rows,sp[,c(keycols,'Spearman_rho','Spearman_P','Spearman_q',
    'Spearman_CI_low','Spearman_CI_high','Valid_replicates','Invalid_replicates')],by=keycols)
  stopifnot(nrow(joined)==81L,!anyNA(joined$Spearman_P))
  write_table(joined,'tables/CORRELATION_PLOT_DATA_AND_SEPARATE_MODEL_RESULTS.csv')
  joined
}

write_records<-function(associations,alpha) {
  write_table(data.frame(Figure=c('3.5.2','3.9B','3.10','SF_3.10','SF_4.6.3'),
    New_file_stem=c('Figure_3_5_2_inflammatory_singscores_PATIENT_AWARE',
      'Figure_3_9B_gastric_singscore_PATIENT_AWARE','Figure_3_10_SSL_HP_gastric_inflammation_PATIENT_AWARE',
      'SF_3_10_all_colorectal_gastric_inflammation_PATIENT_AWARE','SF_4_6_3_four_pages_PATIENT_AWARE'),
    Test=c('Patient-intercept model; BH6','Patient-intercept model; BH6',
      'Biopsy Spearman; patient-block bootstrap; BH3','Marginal biopsy Spearman; patient-block bootstrap; BH3',
      'Existing patient/group/depth-adjusted mixed models; LRT; Holm12'),
    Panel_letters_added=FALSE),'FIGURE_REPLACEMENT_MAP.csv')
  write_note(c('WHAT THIS COMPLETE SCRIPT DOES',
    paste('Score entry point:',SCORE_INPUT_MODE),
    'Run this entire file in a fresh R session. No earlier R script or workspace object is required.',
    'Default from_cpm mode reads the original CPM, gene lengths and metadata from the three paths at the top.',
    'The 24-gene signature and original Hallmark memberships are fixed. This is not a re-run of signature discovery, DE or fgsea.',
    'Saved-scores mode explicitly uses private verified score tables; it is the mode validated in the supplied execution.',
    'From-CPM mode checks sample/patient/group identity, retained genes and all 12 scores against the originals. It stops on disagreement.',
    'All four signatures use the same group model and same correlation method. The expression filter applies only to the Hallmarks.',
    'Singscore uses TPM. Z-scores and ssGSEA retain their original log2(CPM+1) inputs.',
    'Figures have no A/B/C or other automatic panel tags. Individual and combined panels are exported.',
    'Correlation panels show Spearman rho, its P, BH q, n biopsies and N patients. There is no beta or rho asterisk.',
    'The coefficient gives each biopsy equal weight. Whole patients, including all their biopsies, are resampled for uncertainty.',
    'P tests zero marginal Spearman correlation. The bootstrap is approximate, not an exact permutation test.',
    'The primary SSL_HP analysis concerns biopsies within that tissue group. Pooled analyses are marginal and not adjusted for tissue group.',
    'Tissue-adjusted mixed slopes remain separate sensitivity/model results; their P values never label Spearman rho.',
    'Alpha-diversity models retain their existing covariates and inferential question; their figure labels still identify model slopes.',
    'Clinical-score tests, DE/fgsea, organoid analyses and qPCR are outside this score-figure rerun.',
    'Do not reuse an older legend describing beta as the main correlation statistic.',
    'A ZIP is created automatically with inputs, tables, model objects, bootstrap draws and figures.'),'READ_ME_FIRST.txt')
  write_note(c('METHODS: SCORE ASSOCIATIONS',
    'I assessed associations between gastric metaplasia and inflammatory scores using Spearman correlations across biopsies. To account for repeated biopsies, I resampled patients with replacement, retaining all biopsies from each sampled patient and recalculating ranks in every resample. Two-sided P values for zero Spearman correlation and 95% confidence intervals used a studentized patient-block bootstrap on the Fisher-transformed correlation, with a delete-one-patient jackknife standard error. I used 9,999 resamples with fixed seeds. BH correction was applied across the three inflammatory signatures within each sample subset and scoring-method combination. Pearson correlations were calculated separately as a sensitivity analysis.',
    'Each biopsy contributed equally to the correlation estimate. Pooled tissue analyses were marginal associations and did not adjust for tissue group. Tissue-adjusted patient-intercept models were retained as separate analyses.',
    'Gastric and inflammatory group comparisons used the same patient-intercept mixed models, fitted by REML, Satterthwaite tests and BH correction across six comparisons per signature and scoring method. The only signature-specific expression filter was the original Hallmark mean log2(CPM+1)>=2 filter.',
    'Alpha-diversity associations retained ML patient-intercept models, tissue group/BRAF and depth covariates, likelihood-ratio tests and Holm correction across 12 tests per subset.',
    'These methods support inference conditional on the already selected signature and gene sets. This rerun does not create an independent validation cohort.'),'METHODS_RECORD.txt')
  write_note(c('CORRELATION LEGEND',
    'Each point represents one biopsy. Spearman rho describes the association between the two scores. P is the two-sided patient-block bootstrap P value for this correlation, accounting for repeated biopsies from the same patient; q is the BH-adjusted P value across the three inflammatory signatures. n and N give the numbers of biopsies and patients. Grey lines show descriptive linear fits.',
    'For pooled tissue plots add: correlations describe the pooled association and are not adjusted for tissue group.',
    'No stars or asterisks are appended to rho. Automatic panel letters are disabled.',
    'Alpha-diversity panels remain model-based and must retain their own slope/LRT legend, not the correlation legend.'),'LEGEND_CHANGES.txt')
  write_note(c('BOOTSTRAP IMPLEMENTATION',
    'Independent resampling unit: patient. Each draw samples G patients with replacement from G observed patients.',
    'Every selected patient contributes their complete paired biopsy rows; x and y are never resampled separately.',
    'Point estimate: ordinary biopsy-weighted Spearman correlation (average ranks for ties).',
    'Fisher transform z=atanh(rho). SE from delete-one-patient jackknife; in bootstrap samples duplicate patient copies count as distinct drawn blocks.',
    'Studentized bootstrap root: (z_star-z_observed)/SE_star; observed null statistic z_observed/SE_observed.',
    'Two-sided P = (1 + count(abs(root)>=abs(observed statistic)))/(1 + number of valid replicates).',
    'Symmetric 95% bootstrap-t interval uses the 95th percentile of abs(root), then inverse Fisher transformation.',
    'This is a bootstrap approximation for marginal correlation, not a test for pure within-patient change, not a mixed-model test, and not an exact permutation test.',
    'Failure if more than 1% of roots are non-estimable; no silent method substitution. All root arrays and resampling multiplicities are saved.',
    'Methods references: Davison and Hinkley (1997), Bootstrap Methods and Their Application, Chapter 5.',
    'https://stat.ethz.ch/R-manual/R-devel/library/boot/html/boot.ci.html',
    'Huang (2018), Using Cluster Bootstrapping to Analyze Nested Data With a Few Clusters. DOI:10.1177/0013164416678980.',
    'https://pmc.ncbi.nlm.nih.gov/articles/PMC5965657/',
    'The scoring and gene-selection definitions are inherited from the original OSLC scripts, not from these bootstrap references.'),'BOOTSTRAP_METHOD_AND_REFERENCES.txt')
}

if (!identical(Sys.getenv("OSLC_RERUN_DEFINE_ONLY"), "1")) run_reanalysis()
