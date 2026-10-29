#!/bin/bash

##################################################
# RNA seq trim/map/assign pipeline for CCG
# Original by Avi Taylor
# Updated by Ben Wright
##################################################

# Set to abort early on error, and to abort on missing params

set -e
set -u
set -v

. config.txt # load in parameters for this run

# Expected params:
# GENOME_SUPER_DIR - Directory containing the different organism ref genomes
# FASTQ_DIR - Directory where all the fastqs can be found
# SAMPLE_INFO - File with a table of 'fastq filename', 'readable filename' pairs
# WORKING_DIR - Directory to put working data in
# OUTPUT_DIR - Directory to put outputs in
# NUM_CORES - number of cores to allow
# PAIRED - 1 is paired-end, otherwise 0
# STRANDEDNESS - 1 or 2 for featureCounts, as per the check_strandedness script
# FC_OPTS - Some combination of -B and -C, or blank, for featureCounts
# GENOME - named genome in GENOME_SUPER_DIR

# Task summary

echo "----------------------------------------"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Slurm array ID: $SLURM_ARRAY_TASK_ID of $SLURM_ARRAY_TASK_COUNT"
echo "Hostname: "`hostname`
echo "OS: "`uname -s`
echo "User: "`whoami`
echo "Started at: "`date`
echo "----------------------------------------"

# Set up variables

SAMPLE_ID=`cut -f1 ${SAMPLE_INFO} | tail -n+2 | tail -n+${SLURM_ARRAY_TASK_ID} | head -n1`
SAMPLE_NAME=`cut -f2 ${SAMPLE_INFO} | tail -n+2 | tail -n+${SLURM_ARRAY_TASK_ID} | head -n1`

echo "Sample ID: $SAMPLE_ID"
echo "Sample name: $SAMPLE_NAME"
echo "Paired-end: $PAIRED"
echo "Strandedness: $STRANDEDNESS"
echo "Genome: $GENOME"
echo "----------------------------------------"

# Set up working and output dirs and files

SUMMARY_FILE="${OUTPUT_DIR}/featureCounts_${SAMPLE_NAME}_summaries.txt"

if [ $PAIRED -eq 1 ]; then
	INPUT1="${FASTQ_DIR}/${SAMPLE_ID}_R1.fastq.gz"
	INPUT2="${FASTQ_DIR}/${SAMPLE_ID}_R2.fastq.gz"
else
	INPUT1="${FASTQ_DIR}/${SAMPLE_ID}.fastq.gz"
fi
CUT_OUTPUT1="${CUTADAPT_DIR}/${SAMPLE_ID}_${SAMPLE_NAME}_R1_trimmed.fq.gz"
CUT_OUTPUT2="${CUTADAPT_DIR}/${SAMPLE_ID}_${SAMPLE_NAME}_R2_trimmed.fq.gz"
CUT_OUTPUT1_U=`echo "${CUT_OUTPUT1}" | sed 's/...$//'`
CUT_OUTPUT2_U=`echo "${CUT_OUTPUT2}" | sed 's/...$//'`
STAR_PREFIX="${STAR_DIR}/${SAMPLE_NAME}"
STAR_FILE="${STAR_PREFIX}.Aligned.out.bam"

############
# Cutadapt #
############

echo "--- Running Cutadapt --- "`date`

module load cutadapt/4.9-GCCcore-12.3.0

cutadapt --version

# Need to check for correct adapter sequences here? When coming from different sources?

if [ $PAIRED -eq 1 ]; then
	ADAPTERS="-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
	INPUTS="${INPUT1} ${INPUT2}"
	OUTPUTS="${CUT_OUTPUT1} -p ${CUT_OUTPUT2}"
else
	ADAPTERS="-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
	INPUTS="${INPUT1}"
	OUTPUTS="${CUT_OUTPUT1}"
fi

cutadapt --cores=${NUM_CORES} --nextseq-trim=20 ${ADAPTERS} -m 20 -o ${OUTPUTS} ${INPUTS}  

module purge

########
# STAR #
########

echo "--- Running STAR --- "`date`

gunzip ${CUT_OUTPUT1}
if [ $PAIRED -eq 1 ]; then
	gunzip ${CUT_OUTPUT2}
	INPUTS="${CUT_OUTPUT1_U} ${CUT_OUTPUT2_U}"
else
	INPUTS="${CUT_OUTPUT1_U}"
fi

module load STAR/2.7.11b-GCC-12.3.0

STAR --version

STAR --genomeDir ${GENOME_SUPER_DIR}/${GENOME}/index/ \
--runThreadN ${NUM_CORES} \
--readFilesIn ${INPUTS} \
--outFileNamePrefix ${STAR_PREFIX}. \
--outSAMtype BAM Unsorted \
--outSAMattributes Standard \
--outSAMunmapped Within \
--outReadsUnmapped Fastx

rm ${CUT_OUTPUT1_U}
if [ $PAIRED -eq 1 ]; then
	rm ${CUT_OUTPUT2_U}
fi

module purge

##########
# Picard #
##########

THIS_PIC="${PICARD_DIR}/${SAMPLE_NAME}"
PIC_METRICS="${THIS_PIC}/${SAMPLE_NAME}.metrics"
PIC_DUP="${THIS_PIC}/${SAMPLE_NAME}.dup"
PIC_MARKED="${THIS_PIC}/${SAMPLE_NAME}.marked.bam"
if [ $STRANDEDNESS -eq 1 ]; then
  PIC_STRAND="FIRST_READ_TRANSCRIPTION_STRAND"
else
  if [ "$STRANDEDNESS" -eq 0 ]; then
    PIC_STRAND="NONE"
  else
  	PIC_STRAND="SECOND_READ_TRANSCRIPTION_STRAND"
  fi
fi
SORT_BAM="${PICARD_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}.sorted.bam"

mkdir ${THIS_PIC}

module load SAMtools/1.18-GCC-12.3.0

echo "--- Running samtools to sort bam --- "`date`

samtools sort ${STAR_FILE} -o ${SORT_BAM}

module load picard/3.0.0-Java-17
PICARD=${EBROOTPICARD}/picard.jar

echo "--- Running Picard for metrics --- "`date`

set +e # Allow a nonzero exit status from the version command.
java -jar ${PICARD} CollectRnaSeqMetrics --version
set -e

echo "STAR_FILE: ${STAR_FILE}"
echo "PIC_METRICS: ${PIC_METRICS}"
echo "PIC_STRAND: ${PIC_STRAND}"

java -jar ${PICARD} CollectRnaSeqMetrics -I ${STAR_FILE} -O ${PIC_METRICS} \
  -STRAND ${PIC_STRAND} \
  --REF_FLAT ${GENOME_SUPER_DIR}/${GENOME}/refFlat/refFlat.txt

echo "--- Running Picard for duplicates --- "`date`

java -jar ${PICARD} MarkDuplicates -I ${SORT_BAM} -O ${PIC_MARKED} -M ${PIC_DUP}

rm ${PIC_MARKED}
rm ${SORT_BAM}

module purge

#################
# featureCounts #
#################

echo "--- Running featureCounts --- "`date`

module load Subread/2.0.6-GCC-12.3.0

THIS_FC="${FC_DIR}/${SAMPLE_NAME}/"
FC_COUNT="${THIS_FC}/${SAMPLE_NAME}.fC"
FC_SUMMARY="${FC_COUNT}.summary"
mkdir ${THIS_FC}    

featureCounts -v

if [ $PAIRED -eq 1 ]; then
	featureCounts -T ${NUM_CORES} -p -s ${STRANDEDNESS} ${FC_OPTS} \
		--countReadPairs \
		-a ${GENOME_SUPER_DIR}/${GENOME}/gtf/*.gtf \
		--largestOverlap \
		-o ${FC_COUNT} \
		${STAR_FILE}
else
	featureCounts -T ${NUM_CORES} -s ${STRANDEDNESS} \
		-a ${GENOME_SUPER_DIR}/${GENOME}/gtf/*.gtf \
		--largestOverlap \
		-o ${FC_COUNT} \
		${STAR_FILE}
fi

module purge

################################
# Move files to final location #
################################

echo "--- Moving output files --- "`date`

mv ${PIC_DUP} ${OUTPUT_DIR}/
mv ${PIC_METRICS} ${OUTPUT_DIR}/
mv ${FC_COUNT} ${OUTPUT_DIR}/
mv ${FC_SUMMARY} ${OUTPUT_DIR}/ # This is done last as it is the canary that tells the launch script this run is completed

#######
# END #
#######

echo "----------------------------------------"
echo "Finished at: "`date`
echo "----------------------------------------"
