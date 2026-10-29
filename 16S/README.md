# 16S microbiome analysis

Run these scripts from the `16S` directory. Copy `config/config_example.R` to `config/config.R` and set local paths. Script 80 applies the configured identifier aliases, records the applied mappings in `outputs/metadata_name_corrections_applied.csv` and retains all tissue metadata. The example records the original metadata-name correction and the DADA2, control-ASV removal, selective 95th-percentile subsampling and 5% prevalence settings. Exclusions remain as implemented in the source analyses.

Scripts 05 and 10 represent different recorded filtering stages. Scripts 20–110 reconstruct subsequent processing from their original input layout. Existing final phyloseq objects are sufficient for the downstream analyses; preprocessing is not needed merely to refresh thesis reporting. Earlier processing stages can write fixed intermediate filenames, so use a copy of the original data tree if reconstructing them.

## Primary alpha and corrected age sensitivity

Scripts 120 and 180 retain the primary age-free model `Metric ~ Group + Reads + (1 | Patient)`. Both scripts use the same pre-prevalence-filter read-depth covariate for both prevalence versions. Both scripts read `PATH_WORK` from `config/config.R`. Script 180's covariate-sensitivity output follows the ages in the supplied phyloseq object and is the historical reference when that object still contains the earlier ages.

**Use 181_alpha_age_corrected.R for corrected ST4.2.3. The saved source run completed on 20 October and its table matches the submitted thesis.** Set:

```r
Sys.setenv(
  ALPHA_AGE_INPUT_DIR = "/absolute/path/to/private/alpha",
  ALPHA_AGE_OUTPUT_DIR = "/absolute/path/to/existing/alpha/output/parent"
)
source("scripts/181_alpha_age_corrected.R")
```

The input folder must contain:

- `Alpha_input_PRIMARY_no_prevalence_filter_COMMON_READS.csv`: exact saved 64-row primary input.
- `alpha_corrected_age_key.csv`: 64 corrected ages with retained original analysis IDs and the table-label mapping.
- `FINAL_alpha_covariate_sensitivity_PRIMARY_unfiltered_COMMON_READS.csv`: nine original sensitivity results.

The script checks the old-age results, fits the three corrected age models, and reapplies BH across Age, Gender and Location within each metric. The six unchanged Gender/Location raw P values are retained. Their adjusted values can change when the age P value changes. The script writes a new timestamped directory and does not refit or replace the primary age-free analysis. Its saved corrected table and session information are in `../verification/historical/`. This confirms the historical result, not a new execution of the packaged input-loading revision.

## Final beta workflow

`144_beta_final.R` runs the final three-comparison workflow with a shared figure function. Set:

```r
Sys.setenv(
  BETA_WORK_DIR = "/absolute/path/to/16S/work",
  BETA_MATCHED_RUN = "/absolute/path/to/the/reviewed/exact/matched/run",
  BETA_OUTPUT_DIR = "/absolute/path/to/new/results/parent"
)
source("scripts/144_beta_final.R")
```

The work directory contains `phyloseq_final_March26_5pct_without7699.rds`. The matched run must contain `00_PAIRING_AND_QC/07_FINAL_MATCHED_SAMPLE_MANIFEST.csv` and the saved Aitchison and Bray–Curtis objects/statistics from script 175. The script preserves the recorded run name as a default; supply `BETA_MATCHED_RUN` explicitly if the directory changes. It does not choose a new pairing silently.

There are three tests per distance metric: lesion/control and adjacent/control comparisons use patient-tissue mean profiles; lesion/adjacent uses 24 exactly matched biopsies, 12 site-specific pairs from 10 patients. Synchronous patient swaps preserve multiple pairs from the same patient. Holm correction is across the three tissue tests within each distance metric. There is no additional overall three-group inferential test. This final matched beta output supersedes the interim beta inference emitted by script 175.

The archived 18 October beta run supports the statistical workflow. The current v7 plotting revision has been inspected but its figure output has not been independently rendered in this packaging environment.

## Other entry points

- 096: descriptive rarefaction QC. Set `OSLC_RAREFACTION_INPUT` to `phyloseq_March26_prevalence_unfiltered_biological.rds`; optionally set `OSLC_RAREFACTION_OUTPUT`. It reads the object without altering it.
- 130: relative abundance; 150: targeted genus CLR/prevalence.
- 160 and 161: DESeq2 and ANCOM-BC2; 162: figure from an explicitly selected ANCOM-BC2 run. Set `ANCOMBC2_REVIEWED_RUN` in config. Its `FINAL_FIGURES` output location is inside that run; use a copy if preserving previous figures.
- 170: final 22 October species workflow, configured with `OSLC_16S_WORK_DIR`; 175: exact-site paired analyses and saved distances.
- 176: final differential-abundance reporting from existing fit objects. Optional `DESEQ2_REVIEWED_RUN` and `MATCHED_REVIEWED_RUN` config values select explicit saved runs. No new differential-abundance models are fitted.

Fresh runs of 160 and 175 use the final reporting rules already implemented in 176: missing adjusted P values are omitted from volcano plots, and unavailable ANCOM-BC2 tests remain `NA` in cross-method tables. Complete fitted-result tables retain their original values. For existing saved runs, set `DESEQ2_REVIEWED_RUN` and `MATCHED_REVIEWED_RUN` explicitly before running 176.

The version-specific ANCOMBC 2.4.0 missing-value alignment guard is preserved. The final species source replaces the earlier ASV-derived DESeq2 size factors with poscounts size factors estimated from the tested species counts, matching the supplied 22 October analysis. Species labels remain provisional. See `../verification/STATUS.md` for the limits of the current checks.

## Final species source and diagnostics

Set `OSLC_16S_WORK_DIR` to the directory containing `phyloseq_final_March26_5pct_without7699.rds`, then run `scripts/170_species_analysis.R`. The script creates a new timestamped results folder. Primary comparisons retain all available patient-tissue profiles; exact-site and 10% prevalence sensitivities remain separate.

The saved 22 October run contains `REVIEW_REQUIRED` statuses. In the two primary control comparisons, the DESeq2 convergence field is unavailable for 13 and 7 all-zero species, respectively; none of these rows is called significant. This is not evidence of failed convergence for a fitted non-zero taxon. ANCOM-BC2 additionally reports fewer than 50 taxa for bias estimation, numerical warnings and some unavailable contrasts. Keep these limitations visible in interpretation. Species candidate tables record rows passing numerical criteria; they do not override a `REVIEW_REQUIRED` fit status. Do not describe the entire species run as warning-free or independently validated.

Preprocessing scripts write phyloseq intermediates under `16S/outputs`. For downstream use, set the module's work-directory setting (or `OSLC_16S_WORK_DIR` for species) to the directory holding the selected final objects; copying an example configuration does not move those objects.
