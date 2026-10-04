# config/config_example.R
# Copy this file to config/config.R and edit paths locally.
# Do NOT commit config/config.R (it contains private paths).

# Paths (examples)
PATH_RAW_FASTQ <- "/path/to/raw_fastq_on_cluster"
PATH_WORK      <- "/path/to/workdir_on_cluster"
PATH_FILT      <- file.path(PATH_WORK, "filterandtrimmed_final")

# FilterAndTrim parameters (examples)
FILTER_truncLen    <- 0
FILTER_minQ        <- 3
FILTER_minLen      <- 1000
FILTER_maxLen      <- 1600
FILTER_maxN        <- 0
FILTER_maxEE       <- 2
FILTER_rm_phix     <- FALSE
FILTER_multithread <- FALSE

# Error model learning (examples)
PACBIO_BAND_SIZE     <- 32
LEARNERR_multithread <- TRUE
ERROR_DEREP_FILES <- c("SAMPLE_A_derep.rds", "SAMPLE_B_derep.rds")

DADA_multithread <- TRUE

CHIMERA_method <- "consensus"
CHIMERA_multithread <- TRUE

SILVA_TRAIN_FA  <- "/path/to/silva_nr99_v138.1_train_set.fa.gz"
SILVA_SPEC_FA   <- "/path/to/silva_species_assignment_v138.1.fa.gz"
TAX_multithread <- TRUE


METADATA_CSV        <- "/path/to/metadata.csv"
METADATA_SAMPLE_COL <- "Old_Sample_Name"
METADATA_PATCHES <- data.frame(from = character(0), to = character(0))

NC_PATTERN1 <- "NC1"
NC_PATTERN2 <- "NC2"
MIN_READS_KEEP <- 1   # this keeps samples with >1 reads (i.e., drops 0 or 1)

RAREFY_seed <- 123

PREV_FRAC <- 0.05
RAREFY_seed <- 123

ALPHA_GROUP_COL <- "New_Sample_Type"
ALPHA_DOT_TO_HYPHEN <- FALSE
