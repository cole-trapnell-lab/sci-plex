#$ -S /bin/bash
#$ -l mfree=2G
#$ -l h_rt=8:0:0

INPUT_DIR=$1
FILE_LIST=$2
OUTPUT_DIR=$3

module load perl/5.24.0
module load python/2.7.3
module load cutadapt/1.8.3
module load trim_galore/0.4.1

cat $FILE_LIST | while read FILE; do
    trim_galore $INPUT_DIR/$FILE -a AAAAAAAA --three_prime_clip_R1 1 --gzip -o $OUTPUT_DIR
    echo "Processed $FILE"
done

