#!/bin/bash

##########################
# Collates Picard output from multiple samples into a single table
##########################

set -u

# Collate 

. config.txt

# Expects a file called config in the same directory with the following lines:
# OUTPUT_DIR

SUFFIX=$1 # file suffix to collate on (either .dup or .metrics)
OUT=$2 # The output file

FILES=`ls ${OUTPUT_DIR}/*${SUFFIX}`
FIRST=`echo ${FILES} | sed 's/ .*//'`
echo -n -e "SAMPLE\t" > ${OUT}
cat ${FIRST} | grep -A2 "METRICS CLASS" | head -n 2 | tail -n 1 >> ${OUT}

for F in $FILES; do
  echo -n -e `basename ${F} ${SUFFIX}`"\t" >> ${OUT}
	cat ${F} | grep -A2 "METRICS CLASS" | tail -n 1 >> ${OUT}
done
