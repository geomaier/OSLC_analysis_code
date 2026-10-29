# 26 October 2026 readability and reporting review

- Removed unused legacy scoring helpers from tissue signature-discovery script 040; score analysis remains in script 090. Retained discovery calculations are unchanged.
- Replaced revision/chat comments with concise purpose and section descriptions, and clarified per-module configuration and run instructions.
- Connected the documented 16S metadata-name mapping and alpha-analysis work-directory configuration to their entry points.
- Aligned older DESeq2 volcano and matched ANCOM-BC2 reporting paths with the existing missing-test reporting corrections. Model formulas, fitted results, filtering thresholds and multiplicity rules are preserved.
- Corrected qPCR well-level QC flags when no replicate was excluded; the sample means and model calculations are unchanged.
- Corrected the upstream shell unstranded-condition test and made launchers refuse existing run directories and propagate failures through their logging pipelines.
- Redacted a personal source-directory prefix from historical audit metadata, retained contributor credits, and separated historical verification evidence from checks on this reviewed snapshot.

See `verification/review_changes_20261026.json` for changed-file hashes and `verification/STATUS.md` for the completed checks and remaining execution limits.

# 26 October 2026 source preparation (before readability review)

- Selected the final 22 October species analysis, including species-level DESeq2 normalisation and all available primary profiles. Preserved its diagnostics and modelling rules.
- Added submitted Figure 6.2, SF5.6C early-passage labels and Chapter 5/6 volcano-layout scripts. Replaced their personal machine paths with explicit configuration and protected new BRAF figure outputs.
- Made the existing external-validation self-test usable without private configuration; removed automatic dependency installation from that entry point and repaired misleading dependency messages elsewhere.
- Added the completed alpha-age evidence from 20 October, checked its nine-row table against the submitted thesis, and updated the HPA version record from the submitted Methods.
- Parsed all 55 R files and passed the existing external-validation logical tests. Updated script mapping, documentation and exact verification status. No complete statistical rerun, GitHub publication, licence selection or DOI creation is claimed.

Source changes are recorded in `verification/source_update_20261026.json`. The following earlier entry is retained as history; its then-pending alpha-age result was subsequently completed.

# 19 October 2026 release candidate

- Replaced the 10 October bundle's beta entry point with the current 18 October three-comparison script, v7 figures. Only its author-specific work-directory default was changed.
- Replaced the earlier clinical entry point with the corrected-age source from the completed 11 October analysis. Paths are configurable; the 63 individual validation ages are now an external private input.
- Added the pending corrected alpha-age sensitivity as script 181, with its three original input tables moved into private CSV files. The model and comparison code is unchanged.
- Added the saved-fit reporting helper/entry point 176 and descriptive rarefaction QC 096.
- Removed an interactive R history and the superseded beta diagnostic entry point 143 from the public tree.
- Replaced unconditional package-installation blocks in pathway figures and passage GSEA with missing-dependency errors; no statistical functions were changed.
- Updated configuration examples, documented inputs and pending work, and added integrity/syntax-check entry points.
- Preserved the older aggregate verification records as historical evidence. No private clinical data, raw data, or new result tables were added to the public tree.

At the 19 October preparation stage, the package had not been executed in R and the alpha-age correction was pending. Later entries above supersede that historical status. See `verification/release_sources.json` and `verification/STATUS.md` for the current evidence.
