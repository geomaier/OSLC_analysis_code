#!/bin/bash

##################################################
# RNA seq pipeline for CGG
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
# NUM_SAMPLES_FOR_STRANDEDNESS_CHECK - number of samples for perform strandedness check for, rather than entire project
# PAIRED - 1 is paired-end, otherwise 0
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
echo "Genome: $GENOME"
echo "----------------------------------------"

# Set up working and output dirs and files

SUMMARY_FILE="${OUTPUT_DIR}/featureCounts_${SAMPLE_NAME}_summaries.txt"

if [ $PAIRED -eq 1 ]; then
	INPUT1="${FASTQ_DIR}/${SAMPLE_ID}_1.fq.gz"
	INPUT2="${FASTQ_DIR}/${SAMPLE_ID}_2.fq.gz"
else
	INPUT1="${FASTQ_DIR}/${SAMPLE_ID}.fastq.gz"
fi
CUT_OUTPUT1="${FASTQ_DIR}/${SAMPLE_ID}_R1.fq.gz" # Changed as trimming already done
CUT_OUTPUT2="${FASTQ_DIR}/${SAMPLE_ID}_R2.fq.gz"
CUT_OUTPUT1_U=`echo "${CUT_OUTPUT1}" | sed 's/...$//'`
CUT_OUTPUT2_U=`echo "${CUT_OUTPUT2}" | sed 's/...$//'`
STAR_PREFIX="${STAR_DIR}/${SAMPLE_NAME}"
STAR_FILE="${STAR_PREFIX}.Aligned.out.bam"

############
# Cutadapt #
############

echo "--- NOT Running Cutadapt --- "`date`

#module load cutadapt/4.9-GCCcore-12.3.0

# Need to check for correct adapter sequences here? When coming from different sources?

#if [ $PAIRED -eq 1 ]; then
#	ADAPTERS="-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
#	INPUTS="${INPUT1} ${INPUT2}"
#	OUTPUTS="${CUT_OUTPUT1} -p ${CUT_OUTPUT2}"
#else
#	ADAPTERS="-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
#	INPUTS="${INPUT1}"
#	OUTPUTS="${CUT_OUTPUT1}"
#fi

#cutadapt --cores=${NUM_CORES} --nextseq-trim=20 ${ADAPTERS} -m 20 -o ${OUTPUTS} ${INPUTS}  

#module purge

########
# STAR #
########

echo "--- Running STAR --- "`date`

gunzip -k -f ${CUT_OUTPUT1}
if [ $PAIRED -eq 1 ]; then
	gunzip -k -f ${CUT_OUTPUT2}
	INPUTS="${CUT_OUTPUT1_U} ${CUT_OUTPUT2_U}"
else
	INPUTS="${CUT_OUTPUT1_U}"
fi

module load STAR/2.7.11b-GCC-12.3.0

STAR --genomeDir ${GENOME_SUPER_DIR}/${GENOME}/index/ \
--runThreadN ${NUM_CORES} \
--readFilesIn ${INPUTS} \
--outFileNamePrefix ${STAR_PREFIX}. \
--outSAMtype BAM Unsorted \
--outSAMattributes Standard

rm ${CUT_OUTPUT1_U}
if [ $PAIRED -eq 1 ]; then
	rm ${CUT_OUTPUT2_U}
fi

module purge

#################
# featureCounts #
#################

echo "--- Running featureCounts --- "`date`

module load Subread/2.0.6-GCC-12.3.0

TMP1=${SUMMARY_FILE}.tmp
TMP2=${SUMMARY_FILE}.tmp2

echo -n "Sample" > $TMP1
echo -n "${SAMPLE_NAME}" > $TMP2

for S in 1 2
do
	if [ $PAIRED -eq 1 ]; then
		for BC in "" "-B" "-C" "-B -C";	do
			FIXED=`echo "${BC}" | sed 's/ //g' | sed 's/-/_/g'`
			SUFFIX="S${S}${FIXED}"
			THIS_RUN="${FC_DIR}/${SAMPLE_NAME}_${SUFFIX}"
			mkdir ${THIS_RUN}

			featureCounts -T ${NUM_CORES} -p -s ${S} ${BC} \
			--countReadPairs \
			-a ${GENOME_SUPER_DIR}/${GENOME}/gtf/*.gtf \
			--largestOverlap \
			-o ${THIS_RUN}/${SAMPLE_NAME}.fC \
			${STAR_FILE}

			VAL=`grep Assigned ${THIS_RUN}/${SAMPLE_NAME}.fC.summary | cut -f2`
			echo -n " ${SUFFIX}" >> ${TMP1}
			echo -n " ${VAL}" >> ${TMP2}
		done
	else
		SUFFIX="S${S}"
		THIS_RUN="${FC_DIR}/${SAMPLE_ID}_${SAMPLE_NAME}_${SUFFIX}"
		mkdir ${THIS_RUN}

		featureCounts -T ${NUM_CORES} -s ${S} \
		-a ${GENOME_SUPER_DIR}/${GENOME}/gtf/*.gtf \
		--largestOverlap \
		-o ${THIS_RUN}/${SAMPLE_NAME}.fC \
		${STAR_FILE}

		VAL=`grep Assigned ${THIS_RUN}/${SAMPLE_NAME}.fC.summary | cut -f2`
		echo -n " ${SUFFIX}" >> ${TMP1}
		echo -n " ${VAL}" >> ${TMP2}
	fi
done

echo "" >> ${TMP1}
echo "" >> ${TMP2}
cat ${TMP2} >> ${TMP1}
cat ${TMP1} | sed 's/ /\t/g' > ${TMP2}
rm ${TMP1}


# This makes the appearance of the combined summary an atomic operation
mv ${TMP2} ${SUMMARY_FILE} 

rm ${STAR_FILE}

module purge

#######
# END #
#######

echo "----------------------------------------"
echo "Finished at: "`date`
echo "----------------------------------------"

