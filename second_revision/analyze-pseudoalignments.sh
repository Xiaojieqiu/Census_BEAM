#$ -S /bin/bash
#$ -l mfree=8G
#$ -l h_rt=24:0:0

id_file=$1
WORKING_DIR=/net/trapnell/vol1/jspacker/sandberg/kallisto-pseudoalignments/

cat $id_file | while read CELL
do

START_TIME=$(date +%s)

gunzip <$WORKING_DIR/$CELL.pseudoalignment.gz \
| grep -v "^@" | tr '_' '\t' | cut -f 1,3 \
| /net/trapnell/vol1/jspacker/datamash/datamash -g 1 countunique 2 \
| awk '$2 == 1 { print $1; }' \
>$WORKING_DIR/$CELL.allele.informative.reads

gunzip <$WORKING_DIR/$CELL.pseudoalignment.gz \
| grep -v "^@" \
| tr '_' '\t' \
| awk '{
    if (ARGIND == 1)
        informative[$1]=1;
    else {
        n[$2]++;
        if ($1 in informative)
            x[$2]++;
    }
} END {
    for (transcript in n) {
        printf "%s\t%d\t%d\t%.4g\n",
                transcript, x[transcript], n[transcript],
                x[transcript]/n[transcript];
    }
}' $WORKING_DIR/$CELL.allele.informative.reads - \
>$WORKING_DIR/$CELL.transcript.allele.informativeness

#gunzip <$WORKING_DIR/$CELL.pseudoalignment.gz \
#| grep -v "^@" \
#| awk -v CELL=$CELL '{
#    if (ARGIND == 1)
#        r[$1] = 1;
#    else if (ARGIND == 2)
#        t[$1] = $2;
#    else {
#        if (!($1 in seen1)) {
#            n_tot++;
#            seen1[$1] = 1;
#        }
#        if ($1 in r && $2 in t && !($1 in seen2)) {
#            n[t[$2]]++;
#            seen2[$1] = 1;
#        }
#    }
#} END {
#    for (x in n)
#        printf "%s\t%s\t%d\n", CELL, x, n[x];
#}' \
#$WORKING_DIR/$CELL.allele.informative.reads \
#/net/trapnell/vol1/jspacker/mouse-genome/sex.chr.allelic.transcripts \
#- >$WORKING_DIR/$CELL.sex.chr.stats


END_TIME=$(date +%s)
echo "Processed $id in $(($END_TIME-$START_TIME)) seconds"

done

