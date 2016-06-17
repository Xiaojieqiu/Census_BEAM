#$ -S /bin/bash
#$ -l mfree=8G
#$ -l h_rt=3:0:0:0

WORKING_DIR=/net/trapnell/vol1/jspacker/kallisto/

$WORKING_DIR/prepare-ASE-transcriptome.py                    \
$WORKING_DIR/../mouse-genome/gencode.vM9.transcripts.fa.bgz  \
$WORKING_DIR/CAST.transcriptome.snps.gz                      \
$WORKING_DIR/CASTxC57.transcriptome.fa.gz

