#$ -S /bin/bash
#$ -l mfree=8G
#$ -l h_rt=72:0:0

id_file=$1

cat $id_file | while read id
do

START_TIME=$(date +%s)

/net/trapnell/vol1/jspacker/kallisto/kallisto quant        \
-o /net/trapnell/vol1/jspacker/sandberg/kallisto/$id       \
-i /net/trapnell/vol1/jspacker/kallisto/CASTxC57.index     \
-l 200 -s 80 --bias --single                               \
/net/trapnell/vol1/jspacker/sandberg/data/$id.fastq

END_TIME=$(date +%s)
echo "Processed $id in $(($END_TIME-$START_TIME)) seconds"

done

