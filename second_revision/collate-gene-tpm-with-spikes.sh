#$ -S /bin/bash
#$ -l mfree=16G
#$ -l h_rt=24:0:0

WORKING_DIR=/net/trapnell/vol1/jspacker/sandberg/

gunzip <$WORKING_DIR/kallisto-results/with-spikes.transcript.tpm.gz \
| awk '{
    if (ARGIND == 1) {
        gene[$4 "_M"] = $5 "_M";
        gene[$4 "_P"] = $5 "_P";
    } else {
        if ($1 in gene)
            print gene[$1] "\t" $2 "\t" $3;
        else
            print;
    }
}' $WORKING_DIR/../mouse-genome/gencode.vM9.transcripts.bed - \
| gzip >$WORKING_DIR/kallisto-results/with-spikes.gene.tpm.gz

gunzip <$WORKING_DIR/kallisto-results/with-spikes.gene.tpm.gz \
| sort -k1,1 -k2,2 \
| /net/trapnell/vol1/jspacker/datamash/datamash -g 1,2 sum 3 \
| awk '{ printf "%s\t%s\t%.6g\n", $1, $2, $3; }' \
| gzip >$WORKING_DIR/kallisto-results/.with-spikes.gene.tpm.gz

mv $WORKING_DIR/kallisto-results/.with-spikes.gene.tpm.gz \
$WORKING_DIR/kallisto-results/with-spikes.gene.tpm.gz

