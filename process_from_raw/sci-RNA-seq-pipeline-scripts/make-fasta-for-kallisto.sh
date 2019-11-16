#$ -S /bin/bash
#$ -l mfree=2G
#$ -l h_rt=8:0:0

INPUT_DIR=$1
FILE_LIST=$2
OUTPUT_DIR=$3

BATCH_ID=`basename "$FILE_LIST"`
mkdir $OUTPUT_DIR/$BATCH_ID

cat $FILE_LIST | while read FILE; do
    gunzip <$INPUT_DIR/$FILE | awk -v OUTPUT_STEM="$OUTPUT_DIR/$BATCH_ID" '
        BEGIN { FS="|"; }
        {
            x = NR % 4;
            if (x == 1) {
                sample = $2;
                if (sample == "?") next;
                sub(/^@/, ">")
            }
            if (x == 1 || x == 2) {
                print | "gzip >> " OUTPUT_STEM "/" sample ".fa.gz"
            }
        }'

    echo "$FILE"
done

