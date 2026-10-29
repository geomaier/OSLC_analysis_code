#!/bin/bash

##############################
# RNAseq trim/map/assign pipeline batch job submitter for CCG
# Original by Avi Taylor
# Adapted by Ben Wright
##############################

set -e
set -u
set -o pipefail

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

LOG="launch_trim_map_assign.sh.out"

{ # This allows us to redirect all output from this script

echo "----------------------------------"
echo "Launching trim/map/assign pipeline"
echo "----------------------------------"

# Read config

. config.txt

NUM_SAMPLES=`cat ${SAMPLE_INFO} | tail -n+2 | wc -l`
echo "Found this many samples: ${NUM_SAMPLES}"

# Require unused run directories to preserve existing data and results.
mkdir -- "$WORKING_DIR"
mkdir -- "$OUTPUT_DIR"

# Set up subdirectories

CUTADAPT_DIR="${WORKING_DIR}/cutadapt_output/"
STAR_DIR="${WORKING_DIR}/STAR_output/"
FC_DIR="${OUTPUT_DIR}/featureCounts_output/"
PICARD_DIR="${WORKING_DIR}/picard_output/"
mkdir $CUTADAPT_DIR
mkdir $STAR_DIR
mkdir $FC_DIR
mkdir $PICARD_DIR

# Submit jobs

echo "Submitting jobs..."
echo ""

JOB_ID=$(sbatch --parsable --array 1-${NUM_SAMPLES} -n ${NUM_CORES} -p long \
--export CUTADAPT_DIR=${CUTADAPT_DIR},STAR_DIR=${STAR_DIR},FC_DIR=${FC_DIR},PICARD_DIR=${PICARD_DIR} \
--mem-per-cpu=${MEM_PER_CPU} --nodes=${NODES} trim_map_assign.sh)
echo "Submitted job ID: ${JOB_ID}"
echo ""

# Check for completion

FINISHED=-1
NEW_FINISHED=0
touch $OUTPUT_DIR/dummy.fC

while [ $FINISHED -lt $NUM_SAMPLES ]; do
	NEW_FINISHED=`ls $OUTPUT_DIR/*.fC | wc -l`
	set +e
	let "NEW_FINISHED=$NEW_FINISHED - 1"
	set -e
	if [ $NEW_FINISHED -ne $FINISHED ]; then
		FINISHED=$NEW_FINISHED
		echo `date`"  Trim/Map/Assign jobs complete: $FINISHED"
	else
		RUNNING=$(squeue | awk -v job_id="$JOB_ID" 'index($0, job_id) { n++ } END { print n + 0 }')
		if [ $RUNNING -lt 1 ]; then
			echo `date`"  ERROR: No running jobs left!"
			exit 1
		fi 
		sleep 2
		echo -ne '+\r'
		sleep 2
		echo -ne 'x\r'
	fi
done

echo "-----------------------------"
echo "All jobs complete"
echo "-----------------------------"

rm ${OUTPUT_DIR}/dummy.fC # remove dummy file

# Compile outputs into one file

echo "Collating outputs into tables"
echo "-----------------------------"

./collate_featureCounts.sh .fC 7 ${OUTPUT_DIR}/count_table.txt
./collate_featureCounts.sh .fC.summary 2 ${OUTPUT_DIR}/assignment_summary_table.txt
./collate_picard.sh .metrics ${OUTPUT_DIR}/rna_metrics.txt
./collate_picard.sh .dup ${OUTPUT_DIR}/duplication_metrics.txt

# Remove unwanted files

rm ${OUTPUT_DIR}/*.dup
rm ${OUTPUT_DIR}/*.fC
rm ${OUTPUT_DIR}/*.fC.summary
rm ${OUTPUT_DIR}/*.metrics

echo "ALL FINISHED"
echo "-----------------------------"

} 2>&1 | tee >(sed 's/.*\r//' > $LOG) # This redirects to both standard output and a file
