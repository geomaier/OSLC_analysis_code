#!/bin/bash

##############################
# RNAseq trim/map/assign pipeline collation launcher
# Adapted by Ben Wright
##############################

set -e

JOB_ID=$1

if [[ -z ${JOB_ID} ]]; then
    DEPENDENCY=""
else
    DEPENDENCY="--dependency=afterok:$JOB_ID"
fi

set -u

# Expects a file called config in the same directory with the following lines:
# WORKING_DIR
# OUTPUT_DIR
# FASTQ_DIR
# SAMPLE_INFO
# PAIRED
# NUM_CORES
# GENOME_SUPER_DIR
# STRANDEDNESS
# FC_OPTS
# MEM_PER_CPU
# NODES

# Read config

. config.txt

# Submit collation tasks

sbatch ${DEPENDENCY} ./collate_featureCounts.sh .fC 7 ${OUTPUT_DIR}/count_table.txt
sbatch ${DEPENDENCY} ./collate_featureCounts.sh .fC.summary 2 ${OUTPUT_DIR}/assignment_summary_table.txt
sbatch ${DEPENDENCY} ./collate_picard.sh .metrics ${OUTPUT_DIR}/rna_metrics.txt
sbatch ${DEPENDENCY} ./collate_picard.sh .dup ${OUTPUT_DIR}/duplication_metrics.txt

