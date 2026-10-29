# Local configuration

Copy `config_example.R` to `config.R` and configure the original input paths. The local file is ignored by Git. `PATH_WORK` selects the downstream phyloseq objects; script 80 applies the explicit `METADATA_PATCHES` aliases and records the applied mappings. Set `ANCOMBC2_REVIEWED_RUN`, `DESEQ2_REVIEWED_RUN` and `MATCHED_REVIEWED_RUN` when using saved fits for reporting.

Script 144 uses `BETA_WORK_DIR`/`BETA_MATCHED_RUN` and optional `BETA_OUTPUT_DIR`; script 170 uses `OSLC_16S_WORK_DIR`; script 181 uses `ALPHA_AGE_INPUT_DIR`/`ALPHA_AGE_OUTPUT_DIR`. See the module README for the required inputs.
