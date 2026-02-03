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
