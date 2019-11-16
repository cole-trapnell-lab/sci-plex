#$ -S /bin/bash
#$ -l mfree=2G
#$ -l h_rt=8:0:0

INPUT_DIR=$1
R1_FILE_LIST=$2
SCRIPTS_DIR=$3
RT_OLIGO_LIST=$4
LIG_OLIGO_LIST=$5
INDEXING_KEY=$6
OUTPUT_DIR=$7
LOGS_DIR=$8

BATCH_ID=`basename "$R1_FILE_LIST"`

cat $R1_FILE_LIST | while read R1_FILE; do
    R2_FILE=`echo "$R1_FILE" | sed 's/_R1_/_R2_/'`
    PCR_COMBO=`echo "$R1_FILE" | cut -d '_' -f 1`

    paste \
        <(gunzip <$INPUT_DIR/$R1_FILE) \
        <(gunzip <$INPUT_DIR/$R2_FILE) \
    | awk -f $SCRIPTS_DIR/put-read1-info-in-read2_3Level.awk -v PCR_COMBO="$PCR_COMBO" \
        $RT_OLIGO_LIST $LIG_OLIGO_LIST $INDEXING_KEY - \
    | gzip >$OUTPUT_DIR/$PCR_COMBO.fastq.gz

    echo "Processed $PCR_COMBO"
done 2>$LOGS_DIR/$BATCH_ID

