#$ -S /bin/bash
#$ -l mfree=8G
#$ -l h_rt=24:0:0

module load boost/latest

read_fastq=$1
index=$2
output_file=$3

/net/trapnell/vol1/jspacker/hisat2/hisat2 \
-t --dta-cufflinks \
-x $index          \
-U $read_fastq     \
-S $output_file
