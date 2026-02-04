#!/usr/bin/env Rscript
# scripts/120_alpha_diversity.R
# Alpha diversity: Observed, Chao1, Shannon, Simpson
# Group pooling + pairwise LM vs H_COL baseline (adjusting for Reads)

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(data.table)
})

source("config/config.R")

## ============================================================
## 0) INPUTS / OUTPUTS
## ============================================================
phyloseq_in <- file.path("outputs", "phyloseq_after_reseq_noNC_5_percent.rds")
stopifnot(file.exists(phyloseq_in))

out_dir <- file.path("outputs", "alpha_4metrics_PAIRWISELM_vs_HCOL")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ============================================================
## 1) Load phyloseq + clean metadata
## ============================================================
physeq <- readRDS(phyloseq_in)

meta <- as.data.frame(sample_data(physeq))
stopifnot(ALPHA_GROUP_COL %in% colnames(meta))

meta[[ALPHA_GROUP_COL]] <- trimws(as.character(meta[[ALPHA_GROUP_COL]]))
meta[[ALPHA_GROUP_COL]][meta[[ALPHA_GROUP_COL]] == ""] <- NA
sample_data(physeq) <- sample_data(meta)

## ============================================================
## 2) Alpha diversity
## ============================================================
alpha_df <- estimate_richness(
  physeq,
  measures = c("Observed", "Chao1", "Shannon", "Simpson")
)

alpha_df$SampleID <- rownames(alpha_df)

# Optional: fix dots -> hyphens if needed
if (exists("ALPHA_DOT_TO_HYPHEN") && isTRUE(ALPHA_DOT_TO_HYPHEN)) {
  alpha_df$SampleID <- gsub("\\.", "-", alpha_df$SampleID)
}

# Read depth (match by SampleID)
alpha_df$Reads <- sample_sums(physeq)[alpha_df$SampleID]

# Attach group column
meta_df <- as.data.frame(sample_data(physeq))
stopifnot(all(alpha_df$SampleID %in% rownames(meta_df)))

alpha_df$New_Sample_Type <- meta_df[[ALPHA_GROUP_COL]][match(alpha_df$SampleID, rownames(meta_df))]
stopifnot(any(!is.na(alpha_df$New_Sample_Type)))

## ============================================================
## 3) GROUPING (H_COL = H_COL_ND only, drop H_COL_D)
## ============================================================
alpha_df$Group <- as.character(alpha_df$New_Sample_Type)

alpha_df$Group[alpha_df$Group %in% c("SSL", "HP")] <- "SSL_HP"
alpha_df$Group[alpha_df$Group %in% c("H_ADJ_SSL", "H_ADJ_HP")] <- "H_ADJ_SSL_HP"
alpha_df$Group[alpha_df$Group == "H_COL_ND"] <- "H_COL"
alpha_df$Group[alpha_df$Group == "H_COL_D"]  <- NA

alpha_df$Group <- factor(alpha_df$Group)

cat("\nGroup counts after pooling (ND-only baseline):\n")
print(table(alpha_df$Group, useNA = "ifany"))

## ============================================================
## 4) Palette + group order
## ============================================================
pal <- c(
  "H_COL"        = "#00A65A",
  "H_ADJ_SSL_HP" = "#F39C12",
  "H_ADJ_CA"     = "#2980B9",
  "CA"           = "#8E44AD",
  "SSL_HP"       = "#E74C3C"
)

groups_to_plot <- c("H_COL", "H_ADJ_SSL_HP", "SSL_HP", "H_ADJ_CA", "CA")
groups_to_plot <- groups_to_plot[groups_to_plot %in% levels(alpha_df$Group)]
stopifnot("H_COL" %in% groups_to_plot)

dfp <- alpha_df[alpha_df$Group %in% groups_to_plot, , drop = FALSE]
dfp$Group <- factor(dfp$Group, levels = groups_to_plot)

cat("\nN per plotted group:\n")
print(table(dfp$Group))

## ============================================================
## 5) Metrics + helpers
## ============================================================
metrics <- c("Observed", "Chao1", "Shannon", "Simpson")
metric_labels <- c(
  Observed = "Observed richness",
  Chao1    = "Chao1 richness",
  Shannon  = "Shannon index",
  Simpson  = "Simpson index"
)

p_to_star <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ""))))
}

pairwise_lm_p <- function(df, metric, g1, g2) {
  sub <- df[df$Group %in% c(g1, g2), , drop = FALSE]
  sub$Group <- droplevels(factor(sub$Group))
  if (any(table(sub$Group) < 2)) return(NA_real_)
  sub$Group <- relevel(sub$Group, ref = g2)
  anova(lm(as.formula(paste(metric, "~ Group + Reads")), data = sub))["Group", "Pr(>F)"]
}

make_metric_plot <- function(df, metric, baseline = "H_COL") {
  
  ann <- do.call(rbind, lapply(setdiff(levels(df$Group), baseline), function(g) {
    p <- pairwise_lm_p(df, metric, g, baseline)
    s <- p_to_star(p)
    if (s == "") return(NULL)
    data.frame(g = g, p = p, star = s)
  }))
  
  y <- df[[metric]]
  yrng <- diff(range(y, na.rm = TRUE)); if (yrng == 0) yrng <- 1
  
  if (!is.null(ann)) {
    ann <- ann[order(ann$p), ]
    ann$x1 <- match(baseline, levels(df$Group))
    ann$x2 <- match(ann$g, levels(df$Group))
    ann$y  <- max(y, na.rm = TRUE) + seq_len(nrow(ann)) * 0.08 * yrng
    ann$y_text <- ann$y + 0.02 * yrng
  }
  
  p <- ggplot(df, aes(Group, .data[[metric]], fill = Group)) +
    geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85, linewidth = 0.3) +
    geom_jitter(aes(color = Group), width = 0.12, size = 1.4, alpha = 0.65) +
    scale_fill_manual(values = pal) +
    scale_color_manual(values = pal) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none") +
    labs(title = metric, y = metric_labels[[metric]])
  
  if (!is.null(ann)) {
    p <- p +
      geom_segment(data = ann, aes(x = x1, xend = x2, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_text(data = ann, aes(x = (x1 + x2)/2, y = y_text, label = star),
                inherit.aes = FALSE, size = 4)
  }
  
  p
}

## ============================================================
## 6) Save plots
## ============================================================
for (m in metrics) {
  ggsave(
    file.path(out_dir, paste0("Alpha_", m, "_PAIRWISELM_vs_HCOL.png")),
    make_metric_plot(dfp, m),
    width = 3, height = 5.2, dpi = 600
  )
}

# Save the plotted table too (paper traceability)
fwrite(dfp, file.path(out_dir, "alpha_table_plotted.csv"))

cat("\nDONE. Figures saved to:\n", out_dir, "\n")
