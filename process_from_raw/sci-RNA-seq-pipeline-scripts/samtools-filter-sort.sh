#$ -S /bin/bash
#$ -l mfree=1G
#$ -l h_rt=8:0:0
#$ -pe serial 10

INPUT_DIR=$1
FILE_LIST=$2
OUTPUT_DIR=$3

module load samtools/1.4

cat $FILE_LIST | while read FILE; do
    SAMPLE=`basename "$FILE" .Aligned.out.bam`

    samtools view -bh -q 30 -F 4 $INPUT_DIR/$FILE \
    | samtools sort -@ 10 - \
    >$OUTPUT_DIR/$SAMPLE.bam

    echo "Processed $FILE"
done

