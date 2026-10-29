# MUC5AC qPCR: patient-intercept model for 89 verified tissue biopsies.
# Sources: FINAL_FINAL_APRIL26.xlsx (Sheet1) and MUC5ACPCR(2).csv, 10 Oct 2026.
# Read the verified private snapshot via OSLC_PCR_SNAPSHOT_FILE; original
# expression values were checked against workbook log2 fold changes (3 d.p.).
# Primary model: log2FC ~ Group + (1 | Patient), fitted by REML.
# Satterthwaite tests; Holm correction across all ten tissue comparisons.
# Ct=40 denotes nondetection. This analysis does not resolve qPCR censoring.
# Bars show model-estimated means +/- one model-based SE.

# SHA256 FINAL_FINAL_APRIL26.xlsx: 060128d83773387edc8c2229764e492b66ee993eb02aca7d8bb9e4e1a94f79f0
# SHA256 MUC5ACPCR(2).csv: f1268187dda8d7f5c1dc423f1b0fa862dfb08a892b6a0acbbead66393a39de4e

# Inputs and output settings

OUTPUT_PARENT <- Sys.getenv("OSLC_PCR_OUTPUT_PARENT", file.path(getwd(), "outputs"))
SNAPSHOT_FILE <- Sys.getenv("OSLC_PCR_SNAPSHOT_FILE", file.path(getwd(), "private_inputs", "pcr", "snapshot_csv.csv"))
read_snapshot_text <- function(name) {
  if (!file.exists(SNAPSHOT_FILE)) stop("Missing private qPCR snapshot: ", SNAPSHOT_FILE)
  paste(readLines(SNAPSHOT_FILE, warn=FALSE), collapse="\n")
}

# Optional verification against the unchanged wide CSV. Leave NULL to use the
# verified private snapshot. A different/reordered file stops the script.
INPUT_CSV <- NULL

# Confirmed source-label correction (10 Oct 2026): 1544 and 1554 are one patient.
# Preserve the original source label and use 1554 for the analysis.
PATIENT_ID_CORRECTIONS <- c("1544" = "1554")

# Ct=40 substitutions are retained as nondetections, not measured cycles.
PCR_QC_NOTE <- "Confirmed on 10 Oct 2026: undetected reactions were assigned Ct=40 for the expression calculations. Ct=40 therefore denotes nondetection, not measured amplification."
PATIENT_ID_NOTE <- "Resolved on 10 Oct 2026: George confirmed 1544 is a typo for 1554; both refer to the same patient."
COMPARE_MODE <- "vs_SSL_HP"  # or "all"; Holm always includes all 10 contrasts.
FONT_FAMILY <- if (.Platform$OS.type == "windows") "Arial" else "sans"
SEED <- 1L

required <- c("dplyr", "ggplot2", "lme4", "lmerTest", "emmeans", "zip")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_pkgs)) stop("Install missing packages, then rerun: install.packages(c(",
    paste(sprintf('"%s"',missing_pkgs),collapse=", "), "))")
suppressPackageStartupMessages({library(dplyr);library(ggplot2);library(grid)})
if (.Platform$OS.type == "windows" && FONT_FAMILY == "Arial") windowsFonts(Arial=windowsFont("Arial"))
set.seed(SEED)
if (!dir.exists(OUTPUT_PARENT)) dir.create(OUTPUT_PARENT,recursive=TRUE)
run_base <- file.path(OUTPUT_PARENT,paste0("MUC5AC_PATIENT_AWARE_RERUN_",format(Sys.time(),"%Y%m%d_%H%M%S")))
outdir <- run_base
suffix <- 1L
while (dir.exists(outdir)) {outdir <- paste0(run_base,"_",suffix); suffix <- suffix+1L}
dir.create(outdir,recursive=TRUE)
write_table <- function(x,name) write.csv(x,file.path(outdir,name),row.names=FALSE,na="")
plot_order <- c("SSL_HP","H_ADJ_SSL_HP","H_COL","CA","STOMACH")
paper_cols <- c(SSL_HP="#FF0000",H_ADJ_SSL_HP="#E9C61D",H_COL="#1313EC",CA="#825CA6",STOMACH="#F155D0")
group_labels <- setNames(plot_order,plot_order)
font_family <- FONT_FAMILY
font_pt <- 7; x_font_pt <- 6; label_pt <- 6
pt_to_gg <- function(x) x / 2.845276
p_to_label <- function(p) {
  if (!is.finite(p)) return("unavailable")
  if (p<0.0001) return("****")
  if (p<0.001) return("***")
  if (p<0.01) return("**")
  if (p<0.05) return("*")
  paste0("p=",sprintf("%.3f",p))
}
normalise_header <- function(x) {
  keys <- toupper(gsub("[^A-Za-z0-9]","",trimws(x)))
  lookup <- c(SSLHP="SSL_HP",HADJSSLHP="H_ADJ_SSL_HP",HCOL="H_COL",CA="CA",STOMACH="STOMACH",ANTRUM="STOMACH")
  found <- unname(lookup[keys]);ifelse(is.na(found),trimws(x),found)
}

# Load and verify the source snapshot
snapshot_csv <- read_snapshot_text("snapshot_csv.csv")
df <- read.csv(text=snapshot_csv,stringsAsFactors=FALSE,check.names=FALSE,
               colClasses=c(Patient_original="character"))
stopifnot(nrow(df)==89L,!anyDuplicated(df$Workbook_row),
          identical(as.integer(table(factor(df$Group,levels=plot_order))),c(40L,13L,14L,12L,10L)))
if (!is.null(INPUT_CSV)) {
  raw <- read.csv(INPUT_CSV,check.names=FALSE,stringsAsFactors=FALSE)
  names(raw) <- normalise_header(names(raw))
  if (anyDuplicated(names(raw)) || !all(plot_order %in% names(raw))) stop("CSV group headers differ from verified snapshot.")
  for (g in plot_order) {
    original <- raw[[g]]
    nonblank <- !is.na(original) & trimws(as.character(original))!=""
    x <- suppressWarnings(as.numeric(as.character(original[nonblank])))
    expected <- df[df$Group==g,];expected <- expected[order(expected$CSV_row),]
    if (anyNA(x) || length(x)!=nrow(expected) || any(abs(x-expected$MUC5AC_expression)>1e-8))
      stop("CSV values/order differ for ",g,". Rebuild the sample mapping; do not assign IDs by row alone.")
  }
}
ct <- as.matrix(df[paste0("MUC5AC_Ct",1:3)])
ref <- as.matrix(df[paste0("GAPDH_Ct",1:3)])
df$Mean_deltaCt <- rowMeans(ct-ref)
df$Recomputed_log2FC <- 15.913-df$Mean_deltaCt
stopifnot(all(abs(df$Recomputed_log2FC-df$Workbook_log2FC)<1e-8),
          all(abs(df$MUC5AC_expression-df$Workbook_log2FC)<0.000501))
df$MUC5AC_Ct40_wells <- rowSums(abs(ct-40)<1e-8)
df$GAPDH_Ct40_wells <- rowSums(abs(ref-40)<1e-8)
df$Current_Ct_SD <- apply(ct,1,sd)
df$Recorded_SD_differs <- !is.na(df$Recorded_Ct_SD) & abs(df$Recorded_Ct_SD-df$Current_Ct_SD)>0.001
df$Replacement_label <- grepl("replace",df$Sample,ignore.case=TRUE)
df$Patient <- df$Patient_original
if (length(PATIENT_ID_CORRECTIONS)) {
  if (is.null(names(PATIENT_ID_CORRECTIONS)) || anyDuplicated(names(PATIENT_ID_CORRECTIONS)) ||
      anyNA(PATIENT_ID_CORRECTIONS) || any(!grepl("^[0-9]+$",PATIENT_ID_CORRECTIONS)) ||
      !all(names(PATIENT_ID_CORRECTIONS) %in% df$Patient_original)) stop("Invalid patient correction mapping.")
  changed <- df$Patient_original %in% names(PATIENT_ID_CORRECTIONS)
  df$Patient[changed] <- unname(PATIENT_ID_CORRECTIONS[df$Patient_original[changed]])
}
df$Patient <- factor(df$Patient)
df$Group <- factor(df$Group,levels=plot_order)
contrasts(df$Group) <- contr.treatment(length(plot_order))
df$Group_num <- as.numeric(df$Group)
write_table(df,"MUC5AC_patient_linked_input.csv")
write_table(df[df$MUC5AC_Ct40_wells>0 | df$GAPDH_Ct40_wells>0 | df$Recorded_SD_differs |
                 df$Replacement_label | df$Patient_original=="1544",],"MUC5AC_input_QC_flags.csv")
counts <- df %>% group_by(Group) %>% summarise(Biopsies=n(),Patients=n_distinct(Patient),
    Any_MUC5AC_Ct40=sum(MUC5AC_Ct40_wells>0),All_MUC5AC_Ct40=sum(MUC5AC_Ct40_wells==3),.groups="drop")
write_table(counts,"MUC5AC_sample_patient_counts.csv")
write_table(as.data.frame.matrix(table(df$Patient,df$Group)) %>% tibble::rownames_to_column("Patient"),
            "MUC5AC_biopsies_per_patient_and_group.csv")
input_notes <- c("INPUT RECORD: patient-aware analysis of the existing calculation.",
  PATIENT_ID_NOTE,PCR_QC_NOTE,
  "Ct=40 values are retained exactly as in the original calculation, not reclassified as measured amplification.",
  "No SD threshold, failed-well exclusion, or replacement rule has been invented.",
  "Tissue groups come from workbook column B, not the sample-name text.",
  "The embedded CSV values (3 decimal places) are the response; workbook precision is saved for auditing.",
  "The H_COL calibrator is the original biopsy-weighted mean deltaCt, rounded to 15.913.",
  "Subtracting this common calibrator shifts all values equally; group contrasts are unchanged by that shift.",
  "The SD discrepancy flag is an arithmetic audit tolerance, not a sample-exclusion threshold.",
  paste("Patient corrections:",if(length(PATIENT_ID_CORRECTIONS)) paste(names(PATIENT_ID_CORRECTIONS),PATIENT_ID_CORRECTIONS,sep=" -> ",collapse="; ") else "none"),
  paste("Biopsies:",nrow(df),"; distinct patient IDs:",nlevels(df$Patient)))
writeLines(input_notes,file.path(outdir,"INPUT_REVIEW.txt"))

# ----------------------- PATIENT-AWARE MODEL -----------------------
fit_warnings <- character()
model <- withCallingHandlers(lmerTest::lmer(MUC5AC_expression ~ Group + (1|Patient),data=df,
        REML=TRUE,na.action=na.fail,
        control=lme4::lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=200000))),
        warning=function(w) {fit_warnings <<- c(fit_warnings,conditionMessage(w));invokeRestart("muffleWarning")})
saveRDS(model,file.path(outdir,"MUC5AC_patient_random_intercept_model.rds"))
convergence <- model@optinfo$conv$lme4$messages
optimizer_code <- model@optinfo$conv$opt
singular <- lme4::isSingular(model,tol=1e-4)
var_components <- as.data.frame(lme4::VarCorr(model))
write_table(var_components,"MUC5AC_variance_components.csv")
capture.output(summary(model),file=file.path(outdir,"MUC5AC_model_summary.txt"))
writeLines(c(paste("Singular fit:",singular),paste("Optimizer code:",paste(optimizer_code,collapse="; ")),
             paste("Convergence messages:",paste(convergence,collapse="; ")),
             paste("Warnings:",paste(unique(fit_warnings),collapse="; "))),file.path(outdir,"MUC5AC_model_status.txt"))
hard_convergence <- convergence[!grepl("boundary (singular) fit",convergence,fixed=TRUE)]
if (length(hard_convergence) || any(optimizer_code!=0)) stop("Model convergence issue. Review saved status before inference.")
if (singular) warning("Patient variance is on the boundary. Review the saved model status; no model was substituted.")

overall <- as.data.frame(anova(model,ddf="Satterthwaite"))
overall$Term <- rownames(overall);rownames(overall) <- NULL
write_table(overall,"MUC5AC_overall_patient_adjusted_F_test.csv")
emm <- emmeans::emmeans(model,~Group,lmer.df="satterthwaite")
means <- as.data.frame(summary(emm,infer=c(TRUE,FALSE),adjust="none"))
write_table(means,"MUC5AC_model_estimated_group_means.csv")
pair_list <- combn(plot_order,2,simplify=FALSE)
contrast_list <- lapply(pair_list,function(pair) {
    z <- rep(0,length(plot_order));z[match(pair[1],plot_order)]<-1;z[match(pair[2],plot_order)]<--1;z})
names(contrast_list) <- vapply(pair_list,paste,character(1),collapse=" - ")
pair_emm <- emmeans::contrast(emm,method=contrast_list,adjust="none")
pair_est <- as.data.frame(summary(pair_emm,infer=c(TRUE,TRUE),adjust="none"))
pairwise_stats <- data.frame(group1=vapply(pair_list,`[`,character(1),1),
                            group2=vapply(pair_list,`[`,character(1),2),pair_est)
names(pairwise_stats)[names(pairwise_stats)=="p.value"] <- "p_value"
pairwise_stats$p_adj <- p.adjust(pairwise_stats$p_value,method="holm",n=10L)
pairwise_stats$p_label <- vapply(pairwise_stats$p_adj,p_to_label,character(1))
pairwise_stats$p_adjust_method <- "Holm across all 10 tissue comparisons"
pairwise_stats$effect_direction <- "group1 minus group2, in log2FC units"
write_table(pairwise_stats,"MUC5AC_pairwise_patient_adjusted_Holm.csv")

# Reproduce the former tests solely to identify changes; these are not the
# patient-adjusted results and are not used for the new figure annotations.
old <- data.frame(group1=pairwise_stats$group1,group2=pairwise_stats$group2,
    old_unpaired_p=vapply(pair_list,function(pair) wilcox.test(
      df$MUC5AC_expression[df$Group==pair[1]],df$MUC5AC_expression[df$Group==pair[2]],exact=FALSE)$p.value,numeric(1)))
old$old_unpaired_Holm <- p.adjust(old$old_unpaired_p,"holm",n=10L)
old$new_patient_model_Holm <- pairwise_stats$p_adj
old$Different_5percent_decision <- (old$old_unpaired_Holm<0.05)!=(old$new_patient_model_Holm<0.05)
write_table(old,"MUC5AC_old_vs_patient_adjusted_tests.csv")

# Diagnostic plots are for review, not an automatic claim of validity.
diag_df <- data.frame(Sample=df$Sample,Patient=df$Patient,Group=df$Group,
    Fitted=fitted(model),Residual=residuals(model),Standardized=residuals(model)/sigma(model))
write_table(diag_df,"MUC5AC_model_residuals.csv")
if (capabilities("cairo")) {
  cairo_pdf(file.path(outdir,"MUC5AC_model_diagnostics.pdf"),width=8,height=7,family="sans")
} else pdf(file.path(outdir,"MUC5AC_model_diagnostics.pdf"),width=8,height=7)
par(mfrow=c(2,2),mar=c(7,4,3,1),cex.main=0.9)
plot(diag_df$Fitted,diag_df$Residual,xlab="Fitted log2FC",ylab="Residual",main="Residuals versus fitted");abline(h=0,lty=2)
qqnorm(diag_df$Standardized,main="Residual Q-Q plot");qqline(diag_df$Standardized)
boxplot(Residual~Group,data=diag_df,las=2,main="Residual spread by tissue",ylab="Residual",xlab="",cex.axis=0.65)
qqnorm(lme4::ranef(model)$Patient[,1],main="Patient random-intercept Q-Q plot");qqline(lme4::ranef(model)$Patient[,1])
dev.off()

# The original plot geometry is retained. Means and error bars are now derived
# from the patient-aware model. Points remain individual biopsy observations.
summary_df <- data.frame(Group=means$Group,mean=means$emmean,sem=means$SE,
                         x=match(as.character(means$Group),plot_order))
if (COMPARE_MODE == "vs_SSL_HP") {
  
  ref_group <- "SSL_HP"
  
  comparisons_to_show <- pairwise_stats %>%
    filter(group1 == ref_group | group2 == ref_group) %>%
    mutate(
      other_group = ifelse(group1 == ref_group, group2, group1),
      group1 = ref_group,
      group2 = other_group
    ) %>%
    select(-other_group)
  
} else if (COMPARE_MODE == "all") {
  
  comparisons_to_show <- pairwise_stats
  
} else {
  
  stop("COMPARE_MODE must be either 'vs_SSL_HP' or 'all'.")
}

comparisons_to_show <- comparisons_to_show %>%
  mutate(
    x1 = match(group1, plot_order),
    x2 = match(group2, plot_order),
    x_left = pmin(x1, x2),
    x_right = pmax(x1, x2),
    span = x_right - x_left
  ) %>%
  filter(!is.na(x1), !is.na(x2)) %>%
  arrange(span, x_left, x_right)

# Bracket positions
y_min <- min(df$MUC5AC_expression, summary_df$mean-summary_df$sem, na.rm = TRUE)
y_max <- max(df$MUC5AC_expression, summary_df$mean+summary_df$sem, na.rm = TRUE)
y_range <- y_max - y_min

if (!is.finite(y_range) || y_range == 0) {
  y_range <- 1
}

comparisons_to_show <- comparisons_to_show %>%
  mutate(
    bracket_number = row_number(),
    y_bracket = y_max + (0.12 + 0.13 * (bracket_number - 1)) * y_range,
    y_tick = y_bracket - 0.04 * y_range,
    y_text = y_bracket + 0.035 * y_range
  )

write.csv(
  comparisons_to_show,
  file.path(outdir, "MUC5AC_pairwise_comparisons_displayed_on_plot.csv"),
  row.names = FALSE
)

p <- ggplot(
  df,
  aes(x = Group_num, y = MUC5AC_expression, colour = Group)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_point(
    position = position_jitter(width=0.13, height=0, seed=SEED),
    size = 1.6,
    alpha = 0.9
  ) +
  geom_errorbar(
    data = summary_df,
    aes(
      x = x,
      ymin = mean - sem,
      ymax = mean + sem
    ),
    inherit.aes = FALSE,
    width = 0.15,
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_segment(
    data = summary_df,
    aes(
      x = x - 0.25,
      xend = x + 0.25,
      y = mean,
      yend = mean
    ),
    inherit.aes = FALSE,
    linewidth = 0.45,
    colour = "black"
  )

# Add brackets
if (nrow(comparisons_to_show) > 0) {
  p <- p +
    geom_segment(
      data = comparisons_to_show,
      aes(
        x = x_left,
        xend = x_right,
        y = y_bracket,
        yend = y_bracket
      ),
      inherit.aes = FALSE,
      linewidth = 0.35,
      colour = "black"
    ) +
    geom_segment(
      data = comparisons_to_show,
      aes(
        x = x_left,
        xend = x_left,
        y = y_bracket,
        yend = y_tick
      ),
      inherit.aes = FALSE,
      linewidth = 0.35,
      colour = "black"
    ) +
    geom_segment(
      data = comparisons_to_show,
      aes(
        x = x_right,
        xend = x_right,
        y = y_bracket,
        yend = y_tick
      ),
      inherit.aes = FALSE,
      linewidth = 0.35,
      colour = "black"
    ) +
    geom_text(
      data = comparisons_to_show,
      aes(
        x = (x_left + x_right) / 2,
        y = y_text,
        label = p_label
      ),
      inherit.aes = FALSE,
      family = font_family,
      size = pt_to_gg(label_pt),
      colour = "black"
    )
}

top_y <- if (nrow(comparisons_to_show) > 0) {
  max(comparisons_to_show$y_text, y_max, na.rm = TRUE) + 0.08 * y_range
} else {
  y_max + 0.15 * y_range
}

p <- p +
  scale_colour_manual(
    values = paper_cols,
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = seq_along(plot_order),
    labels = group_labels[plot_order],
    limits = c(0.5, length(plot_order) + 0.5)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.08, 0.04))
  ) +
  coord_cartesian(
    ylim = c(
      y_min - 0.08 * y_range,
      top_y
    ),
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = expression(italic(MUC5AC)~"expression (log"[2]~"fold change)")
  ) +
  theme_classic(
    base_family = font_family,
    base_size = font_pt
  ) +
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(
      family = font_family,
      size = x_font_pt,
      colour = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      family = font_family,
      size = font_pt,
      colour = "black"
    ),
    axis.title.y = element_text(
      family = font_family,
      size = font_pt,
      colour = "black",
      margin = margin(r = 4)
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
    
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(5, 5, 5, 5, "mm")
  )

if (interactive()) print(p)

# ============================================================
# SAVE
# ============================================================

out_pdf <- file.path(outdir, "MUC5AC_qPCR_dotplot_patient_adjusted.pdf")
out_png <- file.path(outdir, "MUC5AC_qPCR_dotplot_patient_adjusted.png")
out_jpg <- file.path(outdir, "MUC5AC_qPCR_dotplot_patient_adjusted.jpg")

ggsave(
  filename = out_pdf,
  plot = p,
  width = 82,
  height = 92,
  units = "mm",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = out_png,
  plot = p,
  width = 82,
  height = 92,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = out_jpg,
  plot = p,
  width = 82,
  height = 92,
  units = "mm",
  dpi = 600,
  bg = "white"
)

cat("\nSaved:\n")
cat(out_pdf, "\n")
cat(out_png, "\n")
cat(out_jpg, "\n")


writeLines(c(
 "MUC5AC qPCR rerun: read INPUT_REVIEW.txt before replacing thesis results.",
 "Points are individual biopsies. Black bars are model-estimated means +/- one SE.",
 "These bars replace the original descriptive biopsy means +/- naive SEM.",
 "Primary model: log2 fold change ~ tissue group + (1 | patient), fitted by REML.",
 "Overall group effect: Satterthwaite F test. Pairwise tests: estimated marginal means with Satterthwaite degrees of freedom.",
 "Figure brackets show Holm-adjusted primary P values across all ten tissue comparisons, even when only four brackets are displayed.",
 "**** P<0.0001; *** P<0.001; ** P<0.01; * P<0.05. Other values are printed numerically.",
 "Group mean and contrast 95% confidence intervals in the CSV tables are pointwise, not simultaneous.",
 "Check model diagnostics and experimental QC before treating this output as final.",
 "The old-versus-new table also reflects a change from rank tests to a model of mean log2FC; it is not a pure isolated test of patient adjustment.",
 "Patient 1544 is corrected to 1554; undetected reactions were assigned Ct=40, as confirmed by George.",
 "These clarifications establish the recorded input handling; the model does not explicitly model nondetection.",
 "Methods wording: Undetected reactions were assigned a Ct value of 40 before calculating relative expression."
),file.path(outdir,"READ_ME_FIRST.txt"))
capture.output(sessionInfo(),file=file.path(outdir,"sessionInfo.txt"))
write_table(data.frame(Package=required,Version=vapply(required,function(p)as.character(packageVersion(p)),character(1))),"package_versions.csv")
source_candidate <- tryCatch(sys.frame(1)$ofile,error=function(e)NULL)
if (is.null(source_candidate)) {
  arg <- grep("^--file=",commandArgs(),value=TRUE)
  if (length(arg)) source_candidate <- sub("^--file=","",arg[1])
}
if (!is.null(source_candidate) && file.exists(source_candidate)) file.copy(source_candidate,file.path(outdir,"MUC5AC_patient_aware_rerun.R"),overwrite=FALSE)
zipfile <- paste0(outdir,".zip")
zip::zipr(zipfile,files=list.files(outdir,full.names=TRUE),root=outdir)
cat("\nCompleted patient-aware rerun. Confirmed input handling is recorded in INPUT_REVIEW.txt.\n")
print(counts)
print(pairwise_stats[,c("group1","group2","estimate","p_adj")])
cat("\nResults folder: ",normalizePath(outdir,winslash="/"),"\nZIP to upload: ",normalizePath(zipfile,winslash="/"),"\n",sep="")
