#$ -S /bin/bash
#$ -l mfree=24G
#$ -l h_rt=72:0:0

fasta=$1
index=$2
out=$3
fragsize_mean=$4
fragsize_sd=$5

START_TIME=$(date +%s)

/net/trapnell/vol1/jspacker/kallisto/kallisto quant \
-i $index -o $out --bias --single -l $fragsize_mean -s $fragsize_sd $fasta

END_TIME=$(date +%s)
echo "Finished in $(($END_TIME-$START_TIME)) seconds"

