# Organoid analyses

These scripts reproduce the organoid RNA-seq, gastric-signature, qPCR and tissue–organoid concordance analyses. Input data are external to this repository. Each run creates a new output subfolder and preserves earlier results.

## Scripts and inputs

| Script | Purpose | Input setting |
|---|---|---|
| `reviewed/01_Organoid_RNAseq.R` | Infection and BRAF analyses, passage scores, sensitivity analyses and final figures | `OSLC_ORGANOID_CPM_DIR` |
| `reviewed/02_Passage_GSEA.R` | Gastric-signature GSEA across passages and SSL_HP–H_ADJ_SSL_HP contrasts | `OSLC_ORGANOID_PASSAGE_DIR` |
| `reviewed/03_MUC5AC_qPCR.R` | Technical-replicate QC, MUC5AC trajectories and mixed models | `OSLC_ORGANOID_PCR_DIR` |
| `reviewed/04_Concordance.R` | Concordance of existing tissue and organoid DE results | `OSLC_CONCORDANCE_DIR` |
| `reviewed/05_SF5_6C_early_passage_labels.R` | Submitted SF5.6C labels from three embedded aggregate GSEA rows | Optional `OSLC_FIGURE_OUTPUT_PARENT`; defaults to the working directory |
| `reviewed/06_volcano_legend_layout.R` | Submitted infection and BRAF volcano layouts from saved plotting tables | `OSLC_ORGANOID_RESULTS_DIR` or `OSLC_ORGANOID_TABLES_DIR` |

For script 01, the input directory must contain these files:

- `CGG028_Alt_Infected_cpm_CORRECTED.txt`
- `CGG028_Alt_Infected_DE_results.txt`
- `CGG028_cpm_table.txt`
- `CGG028_First_Passage_SSL_HP_V600E_WT.txt`
- `Metadata_Infection_Organoids.csv`
- `Metadata_Timecourse_Organoids.csv`

Put `FINAL_SSL_METAPLASIA_LIST.txt` in the script 01 and script 04 input directories. The fixed 24-gene list is supplied in `provenance/`. Script 01 also accepts its path through `OSLC_GASTRIC_SIGNATURE_FILE`.

Script 02 expects passage DE files matching `*_Passage_*_P<number>_v_P<number>_DE_results.txt` directly in its input directory, plus `CGG028_First_Passage_SSL_HP_DE_results.txt`, `CGG028_Fourth_Passage_SSL_HP_DE_results.txt` and `CGG028_Eighth_Passage_SSL_HP_DE_results.txt` for the cross-sectional comparisons. It contains the same fixed 24-gene signature.

Script 03 expects `Organoids_MUC5AC_ANALYSIS.csv` in its original ten-column export format, including sample name, passage, organoid type and `deltaCT1`, `deltaCT2`, `deltaCT3`. It checks column names, retains sample rows beginning with `GM_`, and records discarded export artefacts. Its pseudonymous sample IDs identify QC checks against the study's original data.

Script 04 reads the 12 DE-result files listed in its `de_file_table` comparison manifest. All 12 are required by default. It does not generate these upstream DE tables or depend on objects left in memory by another script. See [upstream provenance](../upstream_provenance/README.md) for the supplied DE models and their limitations.

## Run instructions

Install the packages listed in [SOFTWARE.md](../SOFTWARE.md), using the archived session information in `provenance/` when recreating the original environment. Script 01 also needs `nlme` for its later time-course models. Automatic package installation is disabled by default. The MSigDB-dependent scripts may require network access when obtaining gene sets through `msigdbr`. Arial and a Cairo-capable R graphics device are needed to match the original typography and PDF exports.

Run each analysis in a fresh R session. Scripts 01, 03 and 04 clear their workspace at startup. For example, from the repository root:

```r
Sys.setenv(OSLC_ORGANOID_CPM_DIR = "/path/to/organoid_inputs")
source("Organoid_Bulk_RNAseq/reviewed/01_Organoid_RNAseq.R")
```

Scripts 01–04 can be run independently once their inputs are available. Script 05 is a standalone replot of the saved aggregate results. Run script 06 after script 01, or point it at the saved `07_Tables` directory from the completed original analysis.

**Run script 01 to completion.** Its later sections write the submitted score tables and figures. Earlier Wilcoxon exports and layouts are retained in `LEGACY_STYLE` folders for traceability; they are not the final primary score results. The introductory `infection_score_primary_test` switch controls an intermediate audit; the final submitted infection analysis explicitly uses paired t-tests.

| Script | Output parent; each run uses an unused child directory |
|---|---|
| 01 | `<CPM input>/Ch5_Organoid_FINAL/rN` |
| 02 | `<passage input>/Final_24_gene_gastric_metaplasia_GSEA_THESIS_FINAL/rN` |
| 03 | `<qPCR input>/MUC5AC_Organoids_Final_Analysis_Aug26/rN` |
| 04 | `<concordance input>/M26F_DE_PATHWAY_GENE_CONCORDANCE_FINAL_ALL_DE_TABLES/rN` |

Here `rN` is the next unused directory (`r1`, `r2`, and so on). Missing required inputs cause a reported error. Keep generated sample-level tables outside the public repository.

## Statistical scope

Final infection scores use paired tests across three passage-specific infection events from one donor. These passages do not represent three independent donors. The primary BRAF score comparisons use unpaired Welch tests across seven lesion lines, with patient-collapsed and strict-P1 sensitivities retained. The longitudinal cohort contains 13 lines (seven lesion and six adjacent) from six patients; the infection line is separate.

Organoid singscores use CPM ranks; tissue scores use TPM ranks. The gastric signature is not subject to the Hallmark expression filter. In scripts 01 and 04 its single-set GSEA P value is nominal; Hallmark FDR is calculated within the tested Hallmark family. Script 02 additionally reports BH correction across its passage comparisons and separately across the three cross-sectional comparisons. Concordance correction across all 264 comparison–module tests is retained.

The exploratory limma-trend analyses in script 01 use the available CPM inputs and remain labelled as sensitivities. They do not replace the supplied count-based DE analyses. Scientific thresholds, sample exclusions, contrasts, plot settings and final test choices are preserved in this release.

## Submitted figure layouts

- Script 05 retains the archived `P1` comparison key but displays **Early (P1/P2)**, followed by P4 and P8. It replots the embedded results without refitting any model.
- Script 06 reads the four saved volcano CSVs from `07_Tables` and places the legends below the plots. It does not refit differential expression or GSEA.

Use these two plotting scripts for the submitted layouts. The main scripts retain their earlier plot exports for traceability.

## Verification records

The files in `audit/` and the session records in `provenance/` describe earlier checks against saved results. They do not establish that the complete pipeline has been rerun from this release. The repository-level `FILE_MANIFEST.json` identifies the current files.

The 26 October readability review added concise section headings, shortened plotting-edit commentary, and corrected a reporting-only qPCR flag: wells with no excluded replicate now have explicit `FALSE`/`TRUE` exclusion/use indicators in the long QC table. This correction does not alter sample means, QC inclusion, model inputs or statistical results. See the repository verification record for the scope of current checks.
