#$ -S /bin/bash
#$ -l mfree=32G
#$ -l h_rt=24:0:0

module load bcl2fastq/2.20

INPUT_PATH=$1
OUTPUT_PATH=$2
SAMPLE_SHEET=$3

bcl2fastq \
-R $INPUT_PATH \
-o $OUTPUT_PATH \
--sample-sheet $SAMPLE_SHEET \
--barcode-mismatches=1 \
--ignore-missing-positions \
--ignore-missing-controls \
--ignore-missing-filter \
--ignore-missing-bcls \
--no-lane-splitting \
--minimum-trimmed-read-length 15 \
--mask-short-adapter-reads 15

