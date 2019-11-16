#$ -S /bin/bash
#$ -l mfree=3G
#$ -l h_rt=8:0:0

FILTERED_READS_DIR=$1
RMDUP_SPLIT_BED_DIR=$2
BATCH=$3
OUTPUT_DIR=$4

module load samtools/1.4

cat $BATCH | while read PCR_WELL; do
    awk '{
        split($4, arr, "|");
        if (!seen[arr[1]]) {
            seen[arr[1]] = 1;
            count[arr[2]]++;
        }
    } END {
        for (sample in count)
            print sample "\t" count[sample];
    }' $RMDUP_SPLIT_BED_DIR/$PCR_WELL.bed \
    | sort -k1,1 \
    >$OUTPUT_DIR/$PCR_WELL.UMI.count

    samtools view $FILTERED_READS_DIR/$PCR_WELL.bam \
    | cut -d '|' -f 2 \
    | /net/trapnell/vol1/jspacker/datamash/datamash -g 1 count 1 \
    | sort -k1,1 -S 2G \
    | /net/trapnell/vol1/jspacker/datamash/datamash -g 1 sum 2 \
    >$OUTPUT_DIR/$PCR_WELL.read.count

    echo "Processed $PCR_WELL"
done

