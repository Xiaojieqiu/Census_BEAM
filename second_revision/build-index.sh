#$ -S /bin/bash
#$ -l mfree=32G
#$ -l h_rt=8:0:0

WORKING_DIR=/net/trapnell/vol1/jspacker/kallisto/

$WORKING_DIR/kallisto index     \
-i $WORKING_DIR/CASTxC57.index  \
$WORKING_DIR/CASTxC57.transcriptome.fa.gz

