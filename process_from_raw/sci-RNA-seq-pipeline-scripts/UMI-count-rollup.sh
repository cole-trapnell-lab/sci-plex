#$ -S /bin/bash
#$ -l mfree=3G
#$ -l h_rt=8:0:0

INPUT_DIR=$1
FILE_LIST=$2
OUTPUT_DIR=$3

cat $FILE_LIST | while read FILE; do
    awk '$3 == "exonic" || $3 == "intronic" {
        split($1, arr, "|");
        printf "%s|%s_%s_%s\t%s\n",
            arr[2], arr[3], arr[4], arr[5], $2;
    }' $INPUT_DIR/$FILE \
    | sort -k1,1 -k2,2 -S 2G \
    | /net/trapnell/vol1/jspacker/datamash/datamash -g 1,2 count 2 \
    >$OUTPUT_DIR/$FILE

    echo "Processed $FILE"
done

