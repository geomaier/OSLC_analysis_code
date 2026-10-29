# Historical organoid checks

These JSON files record earlier reviews of the supplied organoid scripts and saved results. The source and reviewed hashes in `code_changes.json` identify those earlier versions; use the repository-level `FILE_MANIFEST.json` for the current release.

- `code_changes.json`: earlier input/output edits and the source filenames to which they applied. One personal Windows path has been replaced with a descriptive placeholder and explicitly marked as redacted.
- `saved_result_checks.json`: checks performed on previously supplied saved-result tables.
- `syntax_and_token_checks.json`: earlier syntax-tree and scientific-token comparisons, including the stated parser limitation for `switch`.

These records are retained for traceability. They are not evidence of a fresh end-to-end execution of the current release.
