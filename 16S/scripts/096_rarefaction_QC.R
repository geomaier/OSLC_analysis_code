suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
})

input_file <- Sys.getenv("OSLC_RAREFACTION_INPUT")
output_parent <- Sys.getenv("OSLC_RAREFACTION_OUTPUT", file.path(getwd(), "outputs", "rarefaction_QC"))
if (!nzchar(input_file)) stop("Set OSLC_RAREFACTION_INPUT to the decontaminated, prevalence-unfiltered biological phyloseq RDS.")
excluded_patients <- "7699"

# Use decontaminated, prevalence-unfiltered integer ASV counts.
# This script produces descriptive QC; it does not alter the input RDS.
if (!file.exists(input_file)) stop("Input file not found: ", input_file)
physeq <- readRDS(input_file)
if (!inherits(physeq, "phyloseq")) stop("Input is not a phyloseq object.")

otu <- as(otu_table(physeq), "matrix")
if (taxa_are_rows(physeq)) otu <- t(otu)
if (!is.numeric(otu) || any(!is.finite(otu)) ||
    any(otu < 0) || any(otu != round(otu))) {
  stop("The OTU table must contain finite, non-negative integer counts.")
}
if (nrow(otu) == 0L || ncol(otu) == 0L) stop("The count table is empty.")
if (is.null(rownames(otu)) || anyDuplicated(rownames(otu))) {
  stop("Sample identifiers must be present and unique.")
}

sd <- as.data.frame(sample_data(physeq), check.names = FALSE)
if (!"Sample_Type" %in% names(sd)) stop("Metadata must contain Sample_Type.")
if (anyDuplicated(rownames(sd)) || !all(rownames(otu) %in% rownames(sd))) {
  stop("Sample metadata do not match the count table.")
}
sd <- sd[rownames(otu), , drop = FALSE]

group_map <- c(
  SSL = "SSL_HP", HP = "SSL_HP", SSL_HP = "SSL_HP",
  H_ADJ_SSL = "H_ADJ_SSL_HP", H_ADJ_HP = "H_ADJ_SSL_HP",
  H_ADJ_SSL_HP = "H_ADJ_SSL_HP",
  H_COL = "H_COL", H_COL_ND = "H_COL", H_COL_D = "H_COL"
)
sample_ids <- rownames(otu)
sample_type <- trimws(as.character(sd$Sample_Type))
group <- unname(group_map[sample_type])

# In the supplied sample identifiers, the final underscore-separated
# numeric component is the patient ID, e.g. B10-SSL_TRVS_7699.
patient <- rep(NA_character_, length(sample_ids))
has_patient_suffix <- grepl("_[0-9]+$", sample_ids)
patient[has_patient_suffix] <- sub(".*_", "", sample_ids[has_patient_suffix])
if (any(is.na(patient) & !is.na(group))) {
  stop("Cannot extract a patient ID from: ",
       paste(sample_ids[is.na(patient) & !is.na(group)], collapse = ", "))
}

reads <- rowSums(otu)
reason <- vapply(seq_along(sample_ids), function(i) {
  reasons <- character()
  if (is.na(group[i])) reasons <- c(reasons, "Outside target tissue groups")
  if (!is.na(patient[i]) && patient[i] %in% excluded_patients) {
    reasons <- c(reasons, "Excluded patient")
  }
  if (reads[i] < 2) reasons <- c(reasons, "Fewer than two retained reads")
  paste(reasons, collapse = "; ")
}, character(1))

audit <- data.frame(
  Sample = sample_ids, Patient = patient, Sample_Type = sample_type,
  Group = group, Reads = reads, Included = reason == "",
  Exclusion_Reason = reason, stringsAsFactors = FALSE
)
if ("Batch" %in% names(sd)) audit$Batch <- as.character(sd$Batch)

dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)
out_dir <- tempfile(
  pattern = paste0("reviewed_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
  tmpdir = output_parent
)
if (!dir.create(out_dir)) stop("Could not create output directory: ", out_dir)
write.csv(audit, file.path(out_dir, "sample_inclusion_audit.csv"), row.names = FALSE)
if (!any(audit$Included)) stop("No eligible samples. See ", out_dir)

meta <- audit[audit$Included, , drop = FALSE]
otu <- otu[meta$Sample, , drop = FALSE]
otu <- otu[, colSums(otu) > 0, drop = FALSE]
manifest_columns <- c("Sample", "Patient", "Sample_Type", "Group", "Reads")
if ("Batch" %in% names(meta)) manifest_columns <- c(manifest_columns, "Batch")
write.csv(meta[, manifest_columns, drop = FALSE],
          file.path(out_dir, "included_sample_manifest.csv"), row.names = FALSE)

build_one_curve <- function(counts, sample_name) {
  counts <- counts[counts > 0]
  n <- sum(counts)
  depth_90 <- max(1, floor(0.90 * n))
  depth_99 <- max(1, floor(0.99 * n))
  grid_size <- min(200L, n)
  depths <- sort(unique(c(
    round(seq(1, n, length.out = grid_size)),
    round(exp(seq(0, log(n), length.out = grid_size))),
    depth_90, depth_99, n
  )))
  depths <- depths[depths >= 1 & depths <= n]

  warning_messages <- character()
  richness <- withCallingHandlers(
    as.numeric(vegan::rarefy(counts, sample = depths)),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  observed <- length(counts)
  if (length(richness) != length(depths) ||
      any(!is.finite(richness)) || any(diff(richness) < -1e-7) ||
      abs(tail(richness, 1) - observed) > 1e-7) {
    stop("Unexpected rarefaction output for sample: ", sample_name)
  }

  rich_90 <- richness[match(depth_90, depths)]
  rich_99 <- richness[match(depth_99, depths)]
  rich_final <- richness[match(n, depths)]
  p <- counts / n
  singletons <- sum(counts == 1)

  list(
    curve = data.frame(
      Sample = sample_name, Depth = depths,
      Expected_ASV_Richness = richness, stringsAsFactors = FALSE
    ),
    summary = data.frame(
      Sample = sample_name, Reads = n, Observed_ASVs = observed,
      Depth_90pct = depth_90, Expected_Richness_90pct = rich_90,
      Expected_Richness_Final = rich_final,
      Gain_Last_10pct = max(0, rich_final - rich_90),
      Depth_99pct = depth_99,
      Average_Slope_Last_1pct_per_1000_Reads =
        1000 * max(0, rich_final - rich_99) / (n - depth_99),
      Singleton_ASVs = singletons,
      Goods_Coverage = 1 - singletons / n,
      Shannon = -sum(p * log(p)),
      Simpson = 1 - sum(p^2),
      stringsAsFactors = FALSE
    ),
    warnings = data.frame(
      Sample = rep(sample_name, length(warning_messages)),
      Warning = warning_messages, stringsAsFactors = FALSE
    )
  )
}

results <- lapply(seq_len(nrow(otu)), function(i) {
  build_one_curve(otu[i, ], rownames(otu)[i])
})
curve_df <- do.call(rbind, lapply(results, function(x) x$curve))
metrics <- do.call(rbind, lapply(results, function(x) x$summary))
warning_df <- do.call(rbind, lapply(results, function(x) x$warnings))
stopifnot(identical(as.character(metrics$Sample), as.character(meta$Sample)))

summary_df <- cbind(
  meta[, setdiff(manifest_columns, "Reads"), drop = FALSE],
  metrics[, setdiff(names(metrics), "Sample"), drop = FALSE]
)
curve_df$Group <- meta$Group[match(curve_df$Sample, meta$Sample)]

write.csv(summary_df, file.path(out_dir, "per_sample_rarefaction_summary.csv"),
          row.names = FALSE)
write.csv(curve_df, file.path(out_dir, "per_sample_rarefaction_curve_data.csv"),
          row.names = FALSE)
write.csv(warning_df, file.path(out_dir, "rarefaction_warnings.csv"),
          row.names = FALSE)
if (nrow(warning_df) > 0) {
  message("vegan warnings recorded in rarefaction_warnings.csv:\n",
          paste(unique(warning_df$Warning), collapse = "\n"))
}

group_levels <- c("H_COL", "H_ADJ_SSL_HP", "SSL_HP")
cohort_summary <- do.call(rbind, lapply(group_levels, function(g) {
  keep <- summary_df$Group == g
  data.frame(Group = g, Samples = sum(keep),
             Patients = length(unique(summary_df$Patient[keep])))
}))
write.csv(cohort_summary, file.path(out_dir, "cohort_summary.csv"),
          row.names = FALSE)

# Descriptive associations across biopsies; no independence-based P values.
descriptive_rho <- function(x, y) {
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(NA_real_)
  }
  cor(x, y, method = "spearman")
}
depth_associations <- data.frame(
  Metric = c("Observed_ASVs", "Shannon"),
  Spearman_rho = c(
    descriptive_rho(summary_df$Reads, summary_df$Observed_ASVs),
    descriptive_rho(summary_df$Reads, summary_df$Shannon)
  ),
  Biopsies = nrow(summary_df)
)
write.csv(depth_associations,
          file.path(out_dir, "descriptive_depth_associations.csv"),
          row.names = FALSE)

palette <- c(H_COL = "#0072B2", H_ADJ_SSL_HP = "#009E73", SSL_HP = "#D55E00")
curve_df$Group <- factor(curve_df$Group, levels = group_levels)
summary_df$Group <- factor(summary_df$Group, levels = group_levels)

p_all <- ggplot(curve_df, aes(Depth, Expected_ASV_Richness,
                             group = Sample, colour = Group)) +
  geom_line(linewidth = 0.45, alpha = 0.7) +
  scale_colour_manual(values = palette, drop = FALSE) +
  labs(x = "Reads subsampled", y = "Expected ASV richness",
       colour = "Tissue group") +
  theme_classic(base_size = 12)
ggsave(file.path(out_dir, "per_sample_rarefaction_curves_all_together.pdf"),
       p_all, width = 9, height = 6, useDingbats = FALSE)
ggsave(file.path(out_dir, "per_sample_rarefaction_curves_all_together.png"),
       p_all, width = 9, height = 6, dpi = 300)

write_faceted_pdf <- function() {
  grDevices::pdf(
    file.path(out_dir, "per_sample_rarefaction_curves_faceted.pdf"),
    width = 12, height = 9, onefile = TRUE, useDingbats = FALSE
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  ids <- summary_df$Sample[order(summary_df$Group, summary_df$Reads)]
  for (start in seq(1L, length(ids), by = 16L)) {
    page_ids <- ids[start:min(start + 15L, length(ids))]
    page_df <- curve_df[curve_df$Sample %in% page_ids, , drop = FALSE]
    page_df$Sample <- factor(page_df$Sample, levels = page_ids)
    p <- ggplot(page_df, aes(Depth, Expected_ASV_Richness, colour = Group)) +
      geom_line(linewidth = 0.45) +
      facet_wrap(~ Sample, scales = "free", ncol = 4) +
      scale_colour_manual(values = palette, drop = FALSE) +
      scale_x_continuous(n.breaks = 3) +
      scale_y_continuous(n.breaks = 3) +
      labs(x = "Reads subsampled", y = "Expected ASV richness",
           colour = "Tissue group") +
      theme_classic(base_size = 10) +
      theme(strip.text = element_text(size = 7), legend.position = "bottom")
    print(p)
  }
}
write_faceted_pdf()

p_depth <- ggplot(summary_df, aes(Reads, Observed_ASVs, colour = Group)) +
  geom_point(size = 2.3, alpha = 0.85) +
  scale_colour_manual(values = palette, drop = FALSE) +
  labs(x = "Retained reads", y = "Observed ASVs", colour = "Tissue group") +
  theme_classic(base_size = 12)
ggsave(file.path(out_dir, "observed_asvs_vs_reads.pdf"),
       p_depth, width = 8, height = 6, useDingbats = FALSE)
ggsave(file.path(out_dir, "observed_asvs_vs_reads.png"),
       p_depth, width = 8, height = 6, dpi = 300)

run_notes <- c(
  paste("Input:", normalizePath(input_file)),
  paste("Input MD5:", unname(tools::md5sum(input_file))),
  paste("Output:", normalizePath(out_dir)),
  paste("Excluded patient IDs:", paste(excluded_patients, collapse = ", ")),
  paste("Included biopsies:", nrow(summary_df)),
  paste("Unique patients across included biopsies:",
        length(unique(summary_df$Patient))),
  paste("Samples with no singleton ASVs:", sum(summary_df$Singleton_ASVs == 0)),
  "Only the target tissue groups were retained; all exclusions are in the audit.",
  "No prevalence or abundance filter was applied to retained samples by this script.",
  "Input counts must already have the intended upstream processing and decontamination.",
  "rarefy calculates expected richness; no randomly rarefied analysis dataset was created.",
  "The 90% and 99% depths were rounded down to whole reads and included in the grid.",
  "Tail gains and slopes are descriptive; no asymptote or adequacy classification was used.",
  "Good's coverage is calculated from retained counts and is sensitive to singleton removal.",
  "High coverage does not estimate the fraction of all species that has been detected.",
  "Spearman coefficients are descriptive biopsy-level summaries, without P values.",
  "The overall plot uses common axes; the faceted plots use separate axes for each sample.",
  "A plateau in the processed count table does not rule out depth-related richness differences."
)
writeLines(run_notes, file.path(out_dir, "analysis_notes.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

print(cohort_summary, row.names = FALSE)
print(depth_associations, row.names = FALSE)
message("Included ", nrow(summary_df), " biopsies from ",
        length(unique(summary_df$Patient)), " patients.")
message("Samples with no singleton ASVs: ", sum(summary_df$Singleton_ASVs == 0))
message("Saved outputs to: ", out_dir)
