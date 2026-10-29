#!/bin/bash

##########################
# Collates featureCounts output from multiple samples into a single table
##########################

set -u

# Collate 

. config.txt

# Expects a file called config in the same directory with the following lines:
# OUTPUT_DIR

SUFFIX=$1 # file suffix to collate on (either .fC or .fC.summary)
COL=$2 # The column to extract
OUT=$3 # The output file

echo "Suffix: ${SUFFIX}"
echo "Column: ${COL}"
echo "Output: ${OUT}"

FILES=`ls ${OUTPUT_DIR}/*${SUFFIX}`
echo "Files: ${FILES}"

FIRST=`echo ${FILES} | sed 's/ .*//'`
cat ${FIRST} | grep "^[^#]" | cut -f 1 > tmp

for F in $FILES; do
	cat ${F} | grep "^[^#]" | cut -f 1,${COL} | sed 's~\t.*/~\t~' | sed 's/.Aligned.*//' > tmp2
	join tmp tmp2 > tmp3
	mv tmp3 tmp
done

mv tmp ${OUT}
rm tmp2
