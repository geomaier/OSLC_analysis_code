# Verification status — 26 October 2026

## Checks on the reviewed snapshot

- All **55 R files** passed the actual R parser in **R 4.6.0 via webR 0.6.0**. Among the 54 R files shared with the upload, 46 have identical parsed expressions. Eight contain documented input, configuration, reporting or unused-code changes. One unused configuration fragment was removed and one regression-check entry point was added.
- The existing external-validation base-R self-tests passed without study data or private configuration. They check identifier harmonisation and sample selection.
- `tools/check_reporting_logic.R` passed eight 16S checks and nine qPCR-well cases by evaluating the actual reviewed source expressions. These cover missing/invalid FDR values, unavailable ANCOM-BC2 tests, cross-method direction rules, identifier aliases and replicate flags.
- `tools/check_shell_guards.py` passed **19 checks**, including syntax for all seven shell files, preservation of existing directories, Picard strand labels and logging/queue pipeline handling. It used temporary files and submitted no cluster jobs.
- An independent R syntax-tree review confirmed that all 25 functions reachable from tissue script 040's original entry point remain. All 24 retained helper bodies are identical; the entry point only drops three input checks for data it does not read. The 53 removed function definitions were unreachable. No dynamic function lookup or removed-package dependency was found in the retained path.
- The final species source, including its wrapper and scientific settings, has identical parsed expressions to the uploaded version.
- Current file hashes, the script index, local documentation links and targeted private-input/credential scans are recorded with this snapshot. The original contributor credits and necessary pseudonymous analysis IDs remain.

The runtime checks establish syntax and the stated logical behavior. They do not validate package installation or numerical execution in each historical R/Bioconductor environment.

## Current evidence files

- `r_parse_and_self_test_20261026.json`: current file hashes, R parsing, expression comparisons and test results.
- `shell_review_checks_20261026.json`: the 19 shell checks.
- `tissue_040_cleanup_20261026.json`: removed-function evidence and the independent check of the retained execution path.
- `review_changes_20261026.json`: file-level changes relative to the uploaded archive.
- `current_packaging_checks.json`: current packaging, source and link checks.

## Historical result evidence

Records in `historical/`, the original source-update ledgers and the organoid `audit/` folder describe earlier runs or preparations. Their old hashes identify those earlier sources; they are not hashes of the current files. The pre-review 26 October parsing and packaging records are retained in `historical/`.

The saved corrected alpha-age run completed on 20 October. The supplied comparison records its nine-row agreement with submitted ST4.2.3. This review did not repeat that thesis-table comparison. Saved primary alpha models include singular random-effect fits for Shannon and Simpson; those diagnostics remain.

The saved 22 October species run has small taxon sets, numerical warnings and unavailable ANCOM-BC2 tests. Its `REVIEW_REQUIRED` statuses remain. Numerically passing a threshold does not clear those fit-level limitations. The two control DESeq2 comparisons have unavailable convergence values for all-zero species, none called significant in the recorded result evidence.

Earlier score/qPCR agreement files, the 11 October clinical run and the 18 October beta run remain historical evidence. No saved numerical result has been relabelled as a new execution.

## Work still needed for a reproducible manuscript release

The archive does not include the original sequencing, expression, clinical or qPCR inputs. Targeted execution of the configured loaders and final exports in the recorded environments is still required. The supplied upstream DE models retain their documented patient/passage limitations, and the species warnings still require appropriate interpretation.

An approved data-access statement, a permitted data-backed example where required, confirmed licence/citation metadata, and an archived release identifier remain outstanding. This review did not push to GitHub or create a release tag or DOI. See [the release checklist](../RELEASE_CHECKLIST.md).
