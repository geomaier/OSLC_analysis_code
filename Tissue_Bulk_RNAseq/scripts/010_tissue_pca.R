# PCA of the supplied VSN-transformed tissue expression matrix.
# Run from Tissue_Bulk_RNAseq after configuring config/config.R.
# The matrix is centred without another log transformation or variance scaling.

source("config/config.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

expr_file <- PATH_PCA_MATRIX
meta_file <- PATH_METADATA
out_dir <- OUTPUT_PCA

ensure_new_output(out_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

font_family <- "Arial"
font_pt <- 7

groups <- c(
  "SSL_HP",
  "H_ADJ_SSL_HP",
  "H_COL",
  "CA"
)

paper_cols <- c(
  "SSL_HP"        = "#FF0000",
  "H_ADJ_SSL_HP" = "#E9C61D",
  "H_COL"         = "#1313EC",
  "CA"            = "#825CA6DE",
  "STOMACH"       = "#F155D0"
)

raw <- fread(expr_file, check.names = FALSE)

gene_ids <- trimws(as.character(raw[[1]]))

if (anyNA(gene_ids) || any(!nzchar(gene_ids)) || anyDuplicated(gene_ids)) {
  stop("Gene identifiers must be non-missing and unique.")
}

sample_columns <- names(raw)[-1]
sample_ids <- gsub("-", "_", trimws(sample_columns))

if (anyDuplicated(sample_ids)) {
  stop("Duplicate sample identifiers were found in the expression matrix.")
}

expr <- suppressWarnings(
  as.matrix(
    raw[, lapply(.SD, as.numeric), .SDcols = sample_columns]
  )
)

if (anyNA(expr)) {
  stop("The expression matrix contains missing or non-numeric values.")
}

rownames(expr) <- gene_ids
colnames(expr) <- sample_ids

rm(raw)

meta <- fread(meta_file)
setnames(meta, trimws(names(meta)))

if (!all(c("Display", "Group") %in% names(meta))) {
  stop("Metadata must contain columns named 'Display' and 'Group'.")
}

meta[, `:=`(
  Sample_ID = gsub("-", "_", trimws(as.character(Display))),
  Group = gsub(
    "[[:space:]-]+",
    "_",
    toupper(trimws(as.character(Group)))
  )
)]

meta <- meta[
  !is.na(Sample_ID) &
    nzchar(Sample_ID) &
    !is.na(Group) &
    Group != "DELETE"
]

meta[Group == "ANTRUM", Group := "STOMACH"]

if (anyDuplicated(meta$Sample_ID)) {
  stop("Duplicate sample identifiers were found in the metadata.")
}

meta[, Plot_Group := fcase(
  Group %chin% c("SSL", "HP", "SSL_HP"),
  "SSL_HP",

  Group %chin% c("H_ADJ_SSL", "H_ADJ_HP", "H_ADJ_SSL_HP"),
  "H_ADJ_SSL_HP",

  Group == "H_COL",
  "H_COL",

  Group == "CA",
  "CA",

  Group == "STOMACH",
  "STOMACH",

  default = Group
)]

meta <- meta[Plot_Group %chin% groups]

missing_samples <- setdiff(meta$Sample_ID, colnames(expr))

if (length(missing_samples) > 0) {
  fwrite(
    data.table(Sample_ID = missing_samples),
    file.path(
      out_dir,
      "Panel_D_metadata_samples_not_in_expression_matrix.csv"
    )
  )

  warning(
    length(missing_samples),
    " selected metadata sample(s) were absent from the expression matrix."
  )
}

meta <- meta[Sample_ID %chin% colnames(expr)]

if (nrow(meta) == 0) {
  stop("No metadata samples matched the expression matrix.")
}

missing_groups <- setdiff(groups, unique(meta$Plot_Group))

if (length(missing_groups) > 0) {
  stop(
    "No matched samples were found for: ",
    paste(missing_groups, collapse = ", ")
  )
}

meta[, Plot_Group := factor(Plot_Group, levels = groups)]

expr <- expr[, meta$Sample_ID, drop = FALSE]

stopifnot(
  identical(colnames(expr), meta$Sample_ID)
)

expr <- expr[
  rowSums(!is.finite(expr)) == 0,
  ,
  drop = FALSE
]

gene_variance <- apply(expr, 1, var)

keep <- is.finite(gene_variance) & gene_variance > 0

expr <- expr[keep, , drop = FALSE]
gene_variance <- gene_variance[keep]

if (length(gene_variance) < 2) {
  stop("Too few variable genes remained for PCA.")
}

selected <- head(
  order(gene_variance, decreasing = TRUE),
  2000
)

pca_input <- expr[selected, , drop = FALSE]

pca <- prcomp(
  t(pca_input),
  center = TRUE,
  scale. = FALSE
)

variance_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)

scores <- as.data.table(
  pca$x,
  keep.rownames = "Sample_ID"
)

scores[, Plot_Group := meta$Plot_Group[
  match(Sample_ID, meta$Sample_ID)
]]

stopifnot(!anyNA(scores$Plot_Group))

fwrite(
  scores,
  file.path(out_dir, "Panel_D_PCA_scores.csv")
)

fwrite(
  data.table(
    PC = paste0("PC", seq_along(variance_explained)),
    Variance_explained_percent = round(variance_explained, 3)
  ),
  file.path(out_dir, "Panel_D_PCA_variance_explained.csv")
)

fwrite(
  data.table(
    Rank = seq_len(nrow(pca_input)),
    Gene = rownames(pca_input),
    Mean_VSN_expression = rowMeans(pca_input),
    Variance_VSN_expression = unname(gene_variance[selected])
  ),
  file.path(out_dir, "Panel_D_PCA_top_2000_variable_genes_VSN.csv")
)

fwrite(
  data.table(
    Plot_Group = groups,
    n = as.integer(table(meta$Plot_Group))
  ),
  file.path(out_dir, "Panel_D_PCA_group_sample_numbers.csv")
)

p <- ggplot(
  scores,
  aes(x = PC1, y = PC2, fill = Plot_Group)
) +
  geom_point(
    shape = 21,
    size = 2.2,
    stroke = 0.35,
    colour = "black",
    alpha = 0.95
  ) +
  scale_fill_manual(
    values = paper_cols[groups],
    breaks = groups,
    limits = groups,
    labels = groups,
    drop = FALSE
  ) +
  labs(
    x = sprintf("PC1 (%.1f%%)", variance_explained[1]),
    y = sprintf("PC2 (%.1f%%)", variance_explained[2])
  ) +
  theme_classic(
    base_family = font_family,
    base_size = font_pt
  ) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(
      family = font_family,
      size = font_pt,
      colour = "black"
    ),
    legend.key.size = unit(3.2, "mm"),
    legend.spacing.y = unit(0.5, "mm"),
    axis.title = element_text(
      family = font_family,
      size = font_pt,
      colour = "black"
    ),
    axis.title.x = element_text(
      margin = margin(t = 4)
    ),
    axis.title.y = element_text(
      margin = margin(r = 4)
    ),
    axis.text = element_text(
      family = font_family,
      size = font_pt,
      colour = "black"
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.35
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.35
    ),
    axis.ticks.length = unit(1.2, "mm"),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    plot.margin = margin(4, 4, 4, 4, "mm")
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        size = 2.2,
        colour = "black",
        alpha = 1
      )
    )
  )

print(p)

output_prefix <- file.path(
  out_dir,
  "Panel_D_bulk_RNAseq_PCA_final"
)

pdf_device <- if (capabilities("cairo")) {
  grDevices::cairo_pdf
} else {
  grDevices::pdf
}

ggsave(
  paste0(output_prefix, ".pdf"),
  p,
  width = 160,
  height = 60,
  units = "mm",
  device = pdf_device
)

for (extension in c("png", "jpg")) {
  ggsave(
    paste0(output_prefix, ".", extension),
    p,
    width = 160,
    height = 60,
    units = "mm",
    dpi = 600,
    bg = "white"
  )
}

writeLines(
  c(
    "Input: direct VSN-transformed expression matrix",
    "Additional log transformation: none",
    "Non-finite and zero-variance genes: excluded",
    paste0("Genes used for PCA: ", nrow(pca_input)),
    "PCA centring: TRUE",
    "PCA unit-variance scaling: FALSE",
    sprintf("PC1 variance explained: %.3f%%", variance_explained[1]),
    sprintf("PC2 variance explained: %.3f%%", variance_explained[2])
  ),
  file.path(out_dir, "Panel_D_PCA_analysis_settings.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(out_dir, "Panel_D_PCA_R_sessionInfo.txt")
)

cat("\nPCA completed successfully.\n")
cat("Samples per group:\n")
print(table(meta$Plot_Group))
cat(
  sprintf(
    "\nPC1: %.1f%% | PC2: %.1f%%\n",
    variance_explained[1],
    variance_explained[2]
  )
)
