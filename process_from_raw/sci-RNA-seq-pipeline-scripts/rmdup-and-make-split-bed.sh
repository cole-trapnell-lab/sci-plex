#$ -S /bin/bash
#$ -l mfree=4G
#$ -l h_rt=8:0:0

INPUT_DIR=$1
FILE_LIST=$2
SCRIPTS_DIR=$3
OUTPUT_DIR=$4

module load samtools/1.4
module load bedtools/2.26.0

# this makes "sort" case sensitive
export LC_ALL=C
START_TIME=$SECONDS
cat $FILE_LIST | while read FILE; do
    PCR_WELL=`basename "$FILE" .bam`

    samtools view -h $INPUT_DIR/$FILE \
    | awk -f $SCRIPTS_DIR/rmdup.awk \
    | samtools view -bh \
    | bedtools bamtobed -i - -split \
    | sort -k1,1 -k2,2n -k3,3n -S 3G \
    >$OUTPUT_DIR/$PCR_WELL.bed

    ELAPSED_TIME=$(($SECONDS - $START_TIME))
    echo "Processed $FILE in $ELAPSED_TIME seconds"

done

