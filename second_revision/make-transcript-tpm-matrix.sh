#$ -S /bin/bash
#$ -l mfree=16G
#$ -l h_rt=24:0:0

WORKING_DIR=/net/trapnell/vol1/jspacker/sandberg/

gunzip <$WORKING_DIR/kallisto-results/transcript.tpm.gz  \
| awk -v GROUPING="transcript" -f $WORKING_DIR/make.results.matrix.awk  \
| gzip >$WORKING_DIR/kallisto-results/transcript.tpm.matrix.gz

