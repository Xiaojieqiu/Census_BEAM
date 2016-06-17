#$ -S /bin/bash
#$ -l mfree=8G
#$ -l h_rt=72:0:0

id_file=$1
working_dir=/net/trapnell/vol1/jspacker/sandberg/
data_dir=$working_dir/data/
index=/net/trapnell/vol1/jspacker/mouse-genome/grcm38_plus_cast_snp_gencode_tran/grcm38_plus_cast_snp_gencode_tran

module load boost/latest
module load samtools/latest

cat $id_file | while read id
do

START_TIME=$(date +%s)

/net/trapnell/vol1/jspacker/hisat2/hisat2 \
-t --dta-cufflinks      \
-x $index               \
-U $data_dir/$id.fastq  \
-S $data_dir/$id.sam

CHECKPOINT_1=$(date +%s)
echo "Mapped reads with HiSAT2 in $(($CHECKPOINT_1 - $START_TIME)) seconds"

samtools view -b $data_dir/$id.sam >$data_dir/$id.bam
samtools sort \
-o $data_dir/$id.sorted.bam  \
-O bam -T $data_dir/$id.tmp  \
$data_dir/$id.bam

rm $data_dir/$id.sam
rm $data_dir/$id.bam

END_TIME=$(date +%s)
echo "Converted SAM to sorted BAM in $(($END_TIME-$CHECKPOINT_1)) seconds"
echo "Processed $id in $(($END_TIME-$START_TIME)) seconds"

done

