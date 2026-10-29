# Upstream handover provenance

These supplied handover scripts document the origin of the DE tables used downstream. The R analysis code and its model choices are preserved. They require the original ARNAQ configuration, count matrices, reference files and initial R session objects; the handover scripts also contain stage-selection flags. They are not a standalone installation or an automatic end-to-end pipeline.

The tissue DE models did not account for repeated patients. The infection DE model omitted the passage term; the primary organoid BRAF model omitted patient. Downstream patient-aware score tests do not alter these DE tables or their dependent GSEA results. This limitation remains in the release.

The 26 October readability review made three shell-only corrections: the Picard unstranded-condition test now has valid shell spacing; the two launchers require unused working/output directories instead of deleting existing ones; and launcher failures propagate through their logging pipelines. A conversational comment was replaced with a description of the version-command handling. These changes do not alter alignment, counting parameters or statistical models.

The shell scripts target the original Linux/Slurm module environment. Use paths without spaces and supply its `config.txt`. Historical count collation uses `join` and requires consistently sorted, matching gene keys. Reconstruct and validate that environment before using these handover scripts to regenerate count matrices. Shell syntax checks alone do not validate read processing or count collation.

Credits supplied by Ben Wright: CGG028_analysis.R — Ben Wright; edgeR_analysis.R — Ben Wright, Eshita Sharma and Helen Lockstone. Tissue analysis.R was supplied by Ben Wright. ARNAQ: https://github.com/Queex/arnaq. Read-processing scripts retain their original Avi Taylor / Ben Wright attribution. Preserve these credits when sharing the code. No new licence is asserted for third-party code.
