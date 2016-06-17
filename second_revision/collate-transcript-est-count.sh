#$ -S /bin/bash
#$ -l mfree=8G
#$ -l h_rt=24:0:0

WORKING_DIR=/net/trapnell/vol1/jspacker/sandberg/

ls $WORKING_DIR/kallisto/ | while read CELL; do
    awk -v CELL=$CELL 'NR > 1 { print $1 "\t" CELL "\t" $4; }' \
    $WORKING_DIR/kallisto/$CELL/abundance.tsv
done | gzip \
>$WORKING_DIR/kallisto-results/transcript.est.count.gz

gunzip <$WORKING_DIR/kallisto-results/transcript.est.count.gz \
| sort -k1,1 -k2,2 --compress-program=/bin/gzip \
| gzip >$WORKING_DIR/kallisto-results/.transcript.est.count.gz \

mv $WORKING_DIR/kallisto-results/.transcript.est.count.gz \
   $WORKING_DIR/kallisto-results/transcript.est.count.gz

