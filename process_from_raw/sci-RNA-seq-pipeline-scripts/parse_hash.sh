#$ -S /bin/bash
#$ -l mfree=2G
#$ -l h_rt=8:0:0

INPUT_DIR=$1
R1_FILE_LIST=$2
SCRIPTS_DIR=$3
RT_OLIGO_LIST=$4
INDEXING_KEY=$5
OUTPUT_DIR=$6
LOGS_DIR=$7

BATCH_ID=`basename "$R1_FILE_LIST"`

cat $R1_FILE_LIST | while read R1_FILE; do
    PCR_COMBO=`echo "$R1_FILE" | cut -d '.' -f 1`

    zcat $INPUT_DIR/$R1_FILE \
    | awk -f $SCRIPTS_DIR/parseHash.awk -v PCR_COMBO="$PCR_COMBO" \
        $RT_OLIGO_LIST $INDEXING_KEY - \
    | sed -e 's/|/,/g' \
    | awk 'BEGIN {FS=","; OFS="\t";} {print $2,$3"_"$4"_"$5,$6,$7,$8}' \
    | sort -S 16G -k1,1 -k2,2 -k4,4 -k3,3 \
    | gzip > $OUTPUT_DIR/$PCR_COMBO.hash.gz

    echo "Processed $PCR_COMBO"
done 2>$LOGS_DIR/$BATCH_ID

