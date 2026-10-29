# Tissue bulk RNA sequencing

These scripts reproduce the tissue transcriptomic analyses and figures. They are individual entry points with explicit inputs; the script numbers are an index, not an automatic execution order. The private count, expression and clinical input files are not included in this repository.

## Configure and run

For scripts 010–070, run R from `Tissue_Bulk_RNAseq`, copy `config/config.example.R` to `config/config.R`, and set `OSLC_BULK_DIR` to the existing `March_2026_Final_Analysis` input directory. Review the paths in the copied config before running a script:

```r
Sys.setenv(OSLC_BULK_DIR = "/absolute/path/to/March_2026_Final_Analysis")
source("scripts/010_tissue_pca.R")
```

The config keeps the original input layout and uses a timestamped output suffix. Existing non-empty output directories are rejected. Install the recorded dependencies before analysis; see the repository's `SOFTWARE.md` and the historical package-version files in `provenance/`.

Some scripts consume saved outputs. For example, `030_pathway_figures.R` reads the directory selected by `PATH_GSEA_RESULTS`. To plot a new run of script 020, point that config value to the directory created by that run. Likewise, `PATH_GASTRIC_SIGNATURE` selects the fixed signature text file used by downstream scripts. These paths do not automatically switch to the newest output directory.

## Entry points

| Script | Purpose | Main inputs |
| --- | --- | --- |
| `010_tissue_pca.R` | PCA and sample-group figure | VSN-transformed expression matrix and metadata |
| `020_pathway_gsea.R` | General pathway GSEA | Final tissue differential-expression tables |
| `030_pathway_figures.R` | Pathway and gene-module figures | Saved GSEA results, ranked tables and differential-expression tables |
| `040_gastric_signature_and_scores.R` | Signature discovery, verification, heatmap and gastric GSEA | GTEx/HPA tables and tissue differential-expression tables |
| `045_signature_threshold_sensitivity.R` | Signature-selection sensitivity | Discovery inputs and the fixed final signature |
| `050_external_signature_validation.R` | External validation in GSE76987 and GSE198692 | Configured external-data directory and public GEO inputs |
| `060_bacterial_detection_gsea.R` | GSEA by bacterial detection | Six genus/tissue differential-expression tables and the signature |
| `065_braf_overall_gsea.R` | Overall BRAF-associated GSEA | BRAF V600E versus WT differential-expression table and the signature |
| `070_braf_es_gsea.R` | Five BRAF/detection-stratified GSEA comparisons | Original differential-expression tables and recorded gene sets |
| `080_tissue_MUC5AC_qPCR.R` | Tissue MUC5AC qPCR models and figure | Verified private qPCR snapshot |
| `090_complete_score_analysis.R` | Final score models, correlations and figures | Private reference scores, expression inputs and matched metadata |
| `110_clinical_host_scores.R` | Clinical-score models using corrected ages | Corrected clinical master table and separate age-validation CSV |
| `115_braf_status_figure.R` | Aggregate BRAF-status figure | `BRAF_STATUS_TISSUE.csv` |

Script 040 retains its historical filename but performs only signature discovery and gastric GSEA. Its unused score-analysis helpers have been removed. Script 090 is the final entry point for score calculations, tests and score figures. Script 045 is a sensitivity analysis and is not needed to regenerate the final score figures.

## Final score analysis

Set `OSLC_SCORE_REFERENCE_DIR` to the score directory in the private input archive, and set `OSLC_BASE_DIR` or `OSLC_BULK_DIR` to the original bulk-analysis directory. The default `from_cpm` mode reads the original CPM matrix, GRCh38.113 gene lengths and metadata, then checks all recalculated scores against the saved reference. `OSLC_CPM_FILE`, `OSLC_GENE_LENGTH_FILE` and `OSLC_METADATA_FILE` override individual paths.

To reuse the verified scores explicitly, set `OSLC_SCORE_INPUT_MODE=saved_scores`. This changes the input route; the group models, patient bootstrap and alpha-score models are the same. Use `OSLC_RERUN_OUTPUT_PARENT` for the destination and retain `OSLC_BOOTSTRAP_REPLICATES=9999` for the thesis analysis. Package installation is opt-in with `OSLC_INSTALL_MISSING=1`. Figures retain their submitted layouts and have no automatic panel letters.

Singscore uses TPM. Z-score and ssGSEA analyses use log2(CPM+1). The gastric signature has no expression filter; the three Hallmark signatures retain the original mean log2(CPM+1) >= 2 filter. Patient adjustment changes inference, not the assigned score values. Spearman panels use patient-block bootstrap inference; tissue-adjusted mixed-model association tests are exported separately. Pooled correlations describe marginal associations across the included tissue groups.

## Tissue qPCR

Set `OSLC_PCR_SNAPSHOT_FILE` to the private `pcr/snapshot_csv.csv` and `OSLC_PCR_OUTPUT_PARENT` to an output destination, then run script 080. It retains all 89 biopsies, the documented 1544-to-1554 patient-label correction, Ct=40 nondetections, existing log2 fold changes and Holm correction across all ten tissue comparisons. The fitted patient variance was effectively zero in the saved source analysis; the script does not automatically substitute another model.

## Corrected clinical ages

Script 110 reads a corrected 63-row master table and a separate 63-row age-validation CSV. It retains the 37-lesion, 22-patient clinical analysis and its original sample labels:

```sh
Rscript scripts/110_clinical_host_scores.R /absolute/path/to/MASTER_METADATA_CH4_AGE_CORRECTED.csv /absolute/path/to/output
```

Place `clinical_age_validation.csv` beside the master CSV or set `OSLC_CLINICAL_AGE_KEY` to its location. RStudio users can set `OSLC_CLINICAL_INPUT`, `OSLC_CLINICAL_OUTPUT` and, if needed, `OSLC_CLINICAL_AGE_KEY`, then source the script. `PATH_MASTER_CH4` is not used by this entry point. Original ML models and BH correction across 16 primary tests are retained; batch checks are separate. This exploratory Chapter 6 analysis is separate from the Chapter 4 alpha-diversity age sensitivity.

## BRAF count figure and self-test

For submitted Figure 6.2, set `OSLC_BRAF_STATUS_FILE` to the aggregate `BRAF_STATUS_TISSUE.csv` and run script 115. `OSLC_BRAF_FIGURE_OUTPUT` optionally selects its output parent. Each run creates a new subdirectory and retains the submitted pooling and ANTRUM label.

The existing script 050 self-test runs without private configuration or input data:

```sh
Rscript scripts/050_external_signature_validation.R --self-test
```

It checks identifier harmonisation and sample-selection guards. It does not rerun the external-cohort differential-expression or enrichment analyses.

## Interpretation and reproducibility

General, bacterial-detection and overall BRAF GSEA use signed square-root LR ranks for the recorded inputs. The three original tissue gastric-signature contrasts use signed negative-log10 nominal gene P values; the five stratified gastric-signature comparisons use signed negative-log10 gene FDR ranks. Preserve these distinct definitions when describing the analyses. Each script records the ranking statistic actually used, including any fallback. GSEA does not add patient adjustment to supplied differential-expression tables.

Signature discovery and SSL_HP-versus-H_COL enrichment share the discovery contrast; that enrichment is descriptive rather than independent validation. The external datasets provide separate validation. Historical provenance and package versions describe the saved source runs; see `verification/STATUS.md` at repository level for the scope and limits of the current code checks. Parsing and self-tests do not establish that the private-data analyses have been reproduced.
