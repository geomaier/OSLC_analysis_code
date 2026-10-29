# Analysis entry points

Thesis chapter/figure labels refer to the submitted thesis. Manuscript figure numbers have not been assigned here. Run each entry point in the directory specified by its module README.

| Thesis analysis | Entry point(s) | Inputs / output role |
|---|---|---|
| Chapter 3 tissue PCA | `Tissue_Bulk_RNAseq/scripts/010_tissue_pca.R` | Normalised expression matrix and metadata |
| Chapter 3 pathway enrichment and plots | Tissue scripts `020`, `030` | Supplied DE tables and saved GSEA outputs |
| Gastric signature and threshold sensitivity | Tissue scripts `040`, `045` | GTEx/HPA references and tissue DE tables; 040 is discovery-only |
| GSE76987/GSE198692 validation | Tissue script `050` | Public cohort inputs and validated metadata; `--self-test` checks logic only |
| Bacterial detection and BRAF enrichment | Tissue scripts `060`, `065`, `070` | Bacterial/overall BRAF/stratified BRAF-ES DE tables |
| Tissue MUC5AC qPCR | Tissue script `080` | External private Ct snapshot; patient-aware inference |
| Tissue scores, correlations and alpha-score associations | Tissue script `090` | Expression/metadata or verified saved scores; private reference snapshots |
| Corrected-age clinical models | Tissue script `110` | Corrected clinical master table and validation key |
| Figure 6.2 BRAF tissue counts | Tissue script `115` | Aggregate `BRAF_STATUS_TISSUE.csv`; includes ANTRUM label |
| 16S preprocessing | `16S/scripts/05` and `10`–`110` | Original read/metadata/reference layout; stages 05 and 10 are distinct recorded filters |
| 16S rarefaction QC | 16S script `096` | Prevalence-unfiltered biological phyloseq object |
| Primary alpha diversity / supplementary tables | 16S scripts `120`, `180` | No-prevalence-filter primary data; 5% filter sensitivity; common read-depth covariate |
| Corrected alpha-age sensitivity, ST4.2.3 | 16S script `181` | Three private CSV inputs; saved corrected run completed 20 October |
| Final beta diversity | 16S script `144` | Final phyloseq object and explicit exact-site matched run; replaces interim beta inference from 175 |
| Relative abundance and targeted genus analyses | 16S scripts `130`, `150` | Final phyloseq object |
| Genus DESeq2 / ANCOM-BC2 | 16S scripts `160`, `161`, `162`, `176` | Final objects or explicitly selected saved fits; 162/176 are reporting entry points |
| Species differential abundance | 16S script `170` | Final 22 October source; patient-tissue profiles and species-level size factors; diagnostic restrictions retained |
| Exact-site matched analyses | 16S script `175` | Final phyloseq and preserved matching definition; final beta reporting uses 144 |
| Organoid infection/BRAF analyses | `Organoid_Bulk_RNAseq/reviewed/01_Organoid_RNAseq.R` | Supplied CPM/DE/metadata inputs; use final score exports |
| Organoid passage / MUC5AC / concordance | Organoid scripts `02`, `03`, `04` | Corresponding original data directories; see module instructions |
| SF5.6C passage-label correction | Organoid script `05` | Replots three archived aggregate rows; Early (P1/P2), P4 and P8 labels |
| Chapter 5/6 volcano legend layout | Organoid script `06` | Four saved plot/label CSVs in the completed organoid run's `07_Tables` folder |

The contributed upstream scripts require the original ARNAQ setup and session inputs. They are retained as provenance with their original author credits, rather than advertised as a new automatic end-to-end pipeline.
