#$ -S /bin/bash
#$ -l mfree=4G
#$ -l h_rt=8:0:0

INPUT_DIR=$1
FILE_LIST=$2
EXON_BED=$3
GENE_BED=$4
SCRIPTS_DIR=$5
OUTPUT_DIR=$6

module load bedtools/2.26.0

cat $FILE_LIST | while read FILE; do
    PCR_WELL=`basename "$FILE" .bed`

    bedtools map \
        -a $INPUT_DIR/$FILE \
        -b $EXON_BED \
        -nonamecheck -s -f 0.95 -c 7 -o distinct -delim '|' \
    | bedtools map \
        -a - -b $GENE_BED \
        -nonamecheck -s -f 0.95 -c 4 -o distinct -delim '|' \
    | sort -k4,4 -k2,2n -k3,3n -S 3G \
    | /net/trapnell/vol1/jspacker/datamash/datamash \
        -g 4 first 1 first 2 last 3 first 5 first 6 collapse 7 collapse 8 \
    >$OUTPUT_DIR/$PCR_WELL

    echo "Processed $FILE"
done

