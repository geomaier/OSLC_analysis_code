# Copy to config/config.R and set the local data paths.

PATH_WORK <- "/path/to/work_directory"

PATH_RAW_FASTQ_NOV2025 <- "/path/to/pooled_host_filtered_fastq"
PATH_FILT_NOV2025 <- file.path(PATH_WORK, "filterandtrimmed_final")

PATH_RAW_FASTQ <- file.path(PATH_WORK, "filterandtrimmed_final_March26")
PATH_FILT <- file.path(PATH_WORK, "filterandtrimmed_final_March26_filtered")

# FilterAndTrim parameters
FILTER_truncLen    <- 0
FILTER_minQ        <- 3
FILTER_minLen      <- 1000
FILTER_maxLen      <- 1600
FILTER_maxN        <- 0
FILTER_maxEE       <- 2
FILTER_rm_phix     <- FALSE
FILTER_multithread <- FALSE

# Error model learning
PACBIO_BAND_SIZE     <- 32
LEARNERR_multithread <- TRUE
ERROR_DEREP_FILES <- c(
  "B3-HP_TRVS_7241_filt_derep.rds",
  "B4-H_ADJ_DESC_1554_filt_derep.rds",
  "B5-H_COL_TRVS_7852_filt_derep.rds",
  "B2-SSL_SIG_7154_filt_derep.rds"
)

DADA_multithread <- TRUE

CHIMERA_method <- "consensus"
CHIMERA_multithread <- TRUE

SILVA_TRAIN_FA <- file.path(PATH_WORK, "silva_nr99_v138.1_train_set.fa.gz")
SILVA_SPEC_FA <- file.path(PATH_WORK, "silva_species_assignment_v138.1.fa.gz")
TAX_multithread <- TRUE


METADATA_CSV <- file.path(PATH_WORK, "Final_Metadata_Table_16S_March26_cleaned.csv")
METADATA_SAMPLE_COL <- "Sample_clean"

NC_PATTERN1 <- "NC1"
NC_PATTERN2 <- "NC2"
MIN_READS_KEEP <- 1   # Samples must have more than this number of reads.

RAREFY_seed <- 123

PREV_FRAC <- 0.05

# Existing count exports used by script 143.
BETA_ASV_EXPORT_DIR <- file.path(PATH_WORK, "Beta_diversity_BrayCurtis_ASV_FINAL")
BETA_GENUS_EXPORT_DIR <- file.path(PATH_WORK, "Beta_diversity_BrayCurtis_Genus_FINAL")
