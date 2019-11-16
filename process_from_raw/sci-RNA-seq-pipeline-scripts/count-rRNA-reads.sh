#$ -S /bin/bash
#$ -l mfree=2G
#$ -l h_rt=8:0:0

INPUT_DIR=$1
FILE_LIST=$2
RRNA_BED=$3
OUTPUT_DIR=$4

module load bedtools/2.26.0

cat $FILE_LIST | while read FILE; do
    PCR_WELL=`basename "$FILE" .Aligned.out.bam`

    bedtools intersect -a $INPUT_DIR/$FILE -b $RRNA_BED -c -nonamecheck -bed \
    | awk '{
        split($4, arr, "|");
        if (!seen[arr[1]]) {
            if ($NF > 0)
                rrna_count[arr[2]]++;
            total_count[arr[2]]++;
            seen[arr[1]] = 1;
        }
    } END {
        for (sample in total_count)
            printf "%s\t%d\t%d\n",
                sample, rrna_count[sample], total_count[sample];
    }' \
    >$OUTPUT_DIR/$PCR_WELL

    echo "Processed $FILE"
done

