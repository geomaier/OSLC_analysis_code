#!/bin/bash

##############################
# RNAseq trim/map/assign check strandedness pipeline batch job submitter
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
# NUM_SAMPLES_FOR_STRANDEDNESS_CHECK
# NUM_CORES
# MEM_PER_CPU
# GENOME_SUPER_DIR
# NODES

LOG="launch_check_strandedness.sh.out"

{ # This allows us to redirect all output from this script

echo "--------------------------------"
echo "Launching check for strandedness"
echo "--------------------------------"

# Read config

. config.txt

# Require unused run directories to preserve existing data and results.
mkdir -- "$WORKING_DIR"
mkdir -- "$OUTPUT_DIR"

# Set up subdirectories

CUTADAPT_DIR="${WORKING_DIR}/cutadapt_output/"
STAR_DIR="${WORKING_DIR}/STAR_output/"
FC_DIR="${OUTPUT_DIR}/featureCounts_output/"
mkdir $CUTADAPT_DIR
mkdir $STAR_DIR
mkdir $FC_DIR

# Output file

SUMMARY="${OUTPUT_DIR}/Strandedness_summaries.txt"
DUMMY="${OUTPUT_DIR}/featureCounts_dummy_summaries.txt"

# Submit jobs

echo "Submitting jobs..."
echo ""

JOB_ID=$(sbatch --parsable --array 1-${NUM_SAMPLES_FOR_STRANDEDNESS_CHECK} -n ${NUM_CORES} \
--export CUTADAPT_DIR=${CUTADAPT_DIR},STAR_DIR=${STAR_DIR},FC_DIR=${FC_DIR} \
--mem-per-cpu=${MEM_PER_CPU} --nodes=${NODES} check_strandedness.sh)
echo "Submitted job ID: ${JOB_ID}"
echo ""
# Check for completion

FINISHED=-1
NEW_FINISHED=0
touch ${DUMMY}

while [ $FINISHED -lt $NUM_SAMPLES_FOR_STRANDEDNESS_CHECK ]; do
	NEW_FINISHED=`ls $OUTPUT_DIR/featureCounts_*_summaries.txt | wc -l`
    set +e
	let "NEW_FINISHED=$NEW_FINISHED - 1"
    set -e
	if [ $NEW_FINISHED -ne $FINISHED ]; then
		FINISHED=$NEW_FINISHED
		echo `date`"  Check strandedness jobs complete: $FINISHED"
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

# Compile outputs into one file

echo ""
echo "Collating output summaries..."

rm ${DUMMY}

cat $OUTPUT_DIR/featureCounts_*_summaries.txt | sed -n '1p' > ${SUMMARY}
cat $OUTPUT_DIR/featureCounts_*_summaries.txt | sed -n 'n;p' >> ${SUMMARY}

# Remove individual summaries

rm $OUTPUT_DIR/featureCounts_*_summaries.txt
rm -r $OUTPUT_DIR/featureCounts_output

echo "-----------------------"
echo "All jobs complete"
echo "-----------------------"

} 2>&1 | tee >(sed 's/.*\r//' > $LOG) # This redirects to both standard output and a file
