#$ -S /bin/bash
#$ -l mfree=2G
#$ -l h_rt=8:0:0

INPUT_DIR=$1
FILE_LIST=$2
OUTPUT_DIR=$3

cat $FILE_LIST | while read FILE; do
    FILE_STEM=`echo "$FILE" | cut -d '.' -f 1`

    gunzip <$INPUT_DIR/$FILE | grep "^@" | awk '
        BEGIN { FS="|"; }
        $2 != "?" {
            printf "%s\t%s_%s_%s\t%s\n",
                $2, $3, $4, $5, $6;
    }' | sort -k1,1 -k2,2 -S 1G --compress-program=/bin/gzip \
    | /net/trapnell/vol1/jspacker/datamash/datamash -g 1,2 countunique 3 count 3 \
    >$OUTPUT_DIR/$FILE_STEM

    echo "Processed $FILE"
done

