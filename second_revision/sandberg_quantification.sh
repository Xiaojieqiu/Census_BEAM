#
# Parse GENCODE vM9 transcript model
#

cd /net/trapnell/vol1/jspacker/mouse-genome

wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/gencode.vM9.annotation.gtf.gz

gunzip <gencode.vM9.annotation.gtf.gz \
| grep -v "^#" \
| awk 'BEGIN { FS="\t"; }
       $3 == "transcript" { print $1 "\t" $4 "\t" $5 "\t" $9; }' \
| sed 's/^chr//' \
| sed 's/gene_id "\(ENSMUSG[0-9]\+\)[.][0-9]\+".*transcript_id "\(ENSMUST[0-9]\+\)[.][0-9]\+".*gene_name "\([A-Za-z0-9.-]\+\).*/\2\t\1\t\3/' \
| sort -k1,1 -k2,2n \
>gencode.vM9.transcripts.bed

awk '{ if (($1 == "X" && $2 >= 169969759) || ($1 == "Y" && $2 >= 90745845))
           pseudo=1;
       else
           pseudo=0;
       printf "%s_M\t%s\t%d\tM\t%s\t%s\n", $4, $1, pseudo, $5, $6;
       printf "%s_P\t%s\t%d\tP\t%s\t%s\n", $4, $1, pseudo, $5, $6; }' \
gencode.vM9.transcripts.bed \
| sort -k1,1 \
| awk 'BEGIN { print "transcript\tchr\tis.pseudo.autosomal\tallele\tgene\tsymbol" } { print; }' \
>transcript.attr

#
# Extract info on CAST strain-specific SNPs
# Run Ensembl Variant Effect Predictor to get SNP cDNA positions relative to GENCODE vM9 transcripts
# Build allele-specific transcriptome index
#

cd /net/trapnell/vol1/jspacker/kallisto

wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz.tbi

bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz \
| awk '$6 == "1/1" { id = $1; if (id == ".") id = "snp-" NR;
                     printf "%s\t%d\t%d\t%s/%s\t+\t%s\n", $2, $3, $3, $4, $5, id; }' \
| gzip >CAST.snps.for.vep.gz

gunzip <CAST.snps.for.vep.gz \
| awk '{ printf "%s\t%d\t%d\t%s\t%s\t%s\n", $1, $2-1, $3, $4, $5, $6; }' \
| bedtools intersect -a - -b ../mouse-genome/gencode.vM9.exons.bed -u \
| awk '{ printf "%s\t%s\t%s\t%s\t%s\t%s\n", $1, $3, $3, $4, $5, $6 }' \
| gzip >tmp && mv tmp CAST.snps.for.vep.gz

gunzip CAST.snps.for.vep.gz

perl variant_effect_predictor.pl \
--species mus_musculus \
-i /net/trapnell/vol1/jspacker/kallisto/CAST.snps.for.vep   \
-o /net/trapnell/vol1/jspacker/kallisto/CAST.snp.vep.output \
--no_stats --cache --tab

grep -v "^#" CAST.snp.vep.output \
| cut -f 1,5,7,8,16 \
| grep -P -i "synonymous|missense|stop|start|UTR|exon|frame" \
| cut -f 1,2,4,5 \
| sort -k1,1 -k2,2 >CAST.snp.to.transcript

awk '{ split($4, arr, "/"); printf "%s\t%s\t%s\n", $6, arr[1], arr[2]; }' CAST.snps.for.vep \
| sort -k1,1 | join - CAST.snp.to.transcript | tr ' ' '\t' \
| awk '{ print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t.\t" $6; }' \
| sort -k4,4 -k5,5n | gzip >CAST.transcriptome.snps.gz

qsub -P sage prepare-ASE-transcriptome.sh
qsub -P sage build-index.sh

#
# Run Kallisto and aggregate results into a matrices.
#

cd /net/trapnell/vol1/jspacker/sandberg

ls batches | grep "^batch-" | while read BATCH; do
    qsub -P sage kallisto.sh /net/trapnell/vol1/jspacker/sandberg/batches/$BATCH
done

grep lateblast cell.info | cut -f 1 | shuf | split -l 2 -d - batches/lateblast-

ls batches | grep "^lateblast-" | while read BATCH; do
    qsub -P sage kallisto-with-spikes.sh /net/trapnell/vol1/jspacker/sandberg/batches/$BATCH
done

qsub -P sage collate-transcript-est-count.sh
qsub -P sage collate-transcript-tpm.sh
qsub -P sage collate-gene-est-count.sh
qsub -P sage collate-gene-tpm.sh

qsub -P sage make-transcript-est-count-matrix.sh
qsub -P sage make-transcript-tpm-matrix.sh
qsub -P sage make-gene-est-count-matrix.sh
qsub -P sage make-gene-tpm-matrix.sh

#
# Make CDS from read count and TPM matrices, save as RData files
#

cd /net/trapnell/vol1/jspacker/sandberg/monocle

R

library(monocle)

cell.attr <- new("AnnotatedDataFrame", data=read.table("cell.attr", header=T, row.names=1))
transcript.attr <- new("AnnotatedDataFrame", data=read.table("transcript.attr", header=T, row.names=1))
gene.attr <- new("AnnotatedDataFrame", data=read.table("gene.attr", header=T, row.names=1))

with.spikes.cell.attr <- new("AnnotatedDataFrame", data=read.table("with-spikes.cell.attr", header=T, row.names=1))
with.spikes.transcript.attr <- new("AnnotatedDataFrame", data=read.table("with-spikes.transcript.attr", header=T, row.names=1))
with.spikes.gene.attr <- new("AnnotatedDataFrame", data=read.table("with-spikes.gene.attr", header=T, row.names=1))

sandberg.transcript.read.counts <- read.table(
    gzfile("../kallisto-results/transcript.est.count.matrix.gz"), header=T, row.names=1)
sandberg.transcript.read.counts <- as(as.matrix(sandberg.transcript.read.counts), "sparseMatrix")
sandberg.transcript.read.counts <- newCellDataSet(sandberg.transcript.read.counts,
       phenoData = cell.attr,
       featureData = transcript.attr,
       expressionFamily = negbinomial.size())
sandberg.transcript.read.counts <- sandberg.transcript.read.counts[,row.names(subset(pData(sandberg.transcript.read.counts), !is.na(embryo)))]
save(sandberg.transcript.read.counts, file="sandberg.transcript.read.counts.RData")

sandberg.transcript.tpm <- read.table(
    gzfile("../kallisto-results/transcript.tpm.matrix.gz"), header=T, row.names=1)
sandberg.transcript.tpm <- as(as.matrix(sandberg.transcript.tpm), "sparseMatrix")
sandberg.transcript.tpm <- newCellDataSet(sandberg.transcript.tpm,
       phenoData = cell.attr,
       featureData = transcript.attr,
       lowerDetectionLimit=1,
       expressionFamily = negbinomial.size())
sandberg.transcript.tpm <- sandberg.transcript.tpm[,row.names(subset(pData(sandberg.transcript.tpm), !is.na(embryo)))]
sandberg.transcript.tpm <- detectGenes(sandberg.transcript.tpm, min_expr=1)
save(sandberg.transcript.tpm, file="sandberg.transcript.tpm.RData")

sandberg.gene.read.counts <- read.table(
    gzfile("../kallisto-results/gene.est.count.matrix.gz"), header=T, row.names=1)
sandberg.gene.read.counts <- as(as.matrix(sandberg.gene.read.counts), "sparseMatrix")
sandberg.gene.read.counts <- newCellDataSet(sandberg.gene.read.counts,
       phenoData = cell.attr,
       featureData = gene.attr,
       expressionFamily = negbinomial.size())
sandberg.gene.read.counts <- sandberg.gene.read.counts[,row.names(subset(pData(sandberg.gene.read.counts), !is.na(embryo)))]
save(sandberg.gene.read.counts, file="sandberg.gene.read.counts.RData")

sandberg.gene.tpm <- read.table(
    gzfile("../kallisto-results/gene.tpm.matrix.gz"), header=T, row.names=1)
sandberg.gene.tpm <- as(as.matrix(sandberg.gene.tpm), "sparseMatrix")
sandberg.gene.tpm <- newCellDataSet(sandberg.gene.tpm,
       phenoData = cell.attr,
       featureData = gene.attr,
       lowerDetectionLimit=1,
       expressionFamily = negbinomial.size())
sandberg.gene.tpm <- sandberg.gene.tpm[,row.names(subset(pData(sandberg.gene.tpm), !is.na(embryo)))]
sandberg.gene.tpm <- detectGenes(sandberg.gene.tpm, min_expr=1)
save(sandberg.gene.tpm, file="sandberg.gene.tpm.RData")

with.spikes.sandberg.transcript.read.counts <- read.table(
    gzfile("../kallisto-results/with-spikes.transcript.est.count.matrix.gz"), header=T, row.names=1)
with.spikes.sandberg.transcript.read.counts <- as(as.matrix(with.spikes.sandberg.transcript.read.counts), "sparseMatrix")
with.spikes.sandberg.transcript.read.counts <- newCellDataSet(with.spikes.sandberg.transcript.read.counts,
       phenoData = with.spikes.cell.attr,
       featureData = with.spikes.transcript.attr,
       expressionFamily = negbinomial.size())
save(with.spikes.sandberg.transcript.read.counts, file="with-spikes.sandberg.transcript.read.counts.RData")

with.spikes.sandberg.transcript.tpm <- read.table(
    gzfile("../kallisto-results/with-spikes.transcript.tpm.matrix.gz"), header=T, row.names=1)
with.spikes.sandberg.transcript.tpm <- as(as.matrix(with.spikes.sandberg.transcript.tpm), "sparseMatrix")
with.spikes.sandberg.transcript.tpm <- newCellDataSet(with.spikes.sandberg.transcript.tpm,
       phenoData = with.spikes.cell.attr,
       featureData = with.spikes.transcript.attr,
       lowerDetectionLimit=1,
       expressionFamily = negbinomial.size())
with.spikes.sandberg.transcript.tpm <- detectGenes(with.spikes.sandberg.transcript.tpm, min_expr=1)
save(with.spikes.sandberg.transcript.tpm, file="with-spikes.sandberg.transcript.tpm.RData")

with.spikes.sandberg.gene.read.counts <- read.table(
    gzfile("../kallisto-results/with-spikes.gene.est.count.matrix.gz"), header=T, row.names=1)
with.spikes.sandberg.gene.read.counts <- as(as.matrix(with.spikes.sandberg.gene.read.counts), "sparseMatrix")
with.spikes.sandberg.gene.read.counts <- newCellDataSet(with.spikes.sandberg.gene.read.counts,
       phenoData = with.spikes.cell.attr,
       featureData = with.spikes.gene.attr,
       expressionFamily = negbinomial.size())
save(with.spikes.sandberg.gene.read.counts, file="with-spikes.sandberg.gene.read.counts.RData")

with.spikes.sandberg.gene.tpm <- read.table(
    gzfile("../kallisto-results/with-spikes.gene.tpm.matrix.gz"), header=T, row.names=1)
with.spikes.sandberg.gene.tpm <- as(as.matrix(with.spikes.sandberg.gene.tpm), "sparseMatrix")
with.spikes.sandberg.gene.tpm <- newCellDataSet(with.spikes.sandberg.gene.tpm,
       phenoData = with.spikes.cell.attr,
       featureData = with.spikes.gene.attr,
       lowerDetectionLimit=1,
       expressionFamily = negbinomial.size())
with.spikes.sandberg.gene.tpm <- detectGenes(with.spikes.sandberg.gene.tpm, min_expr=1)
save(with.spikes.sandberg.gene.tpm, file="with-spikes.sandberg.gene.tpm.RData")

system.time({ with.spikes.sandberg.transcript.rpc <- relative2abs(with.spikes.sandberg.transcript.tpm, verbose=T) })
with.spikes.sandberg.transcript.rpc <- newCellDataSet(
    as(as.matrix(with.spikes.sandberg.transcript.rpc), "sparseMatrix"),
    phenoData = with.spikes.cell.attr,
    featureData = with.spikes.transcript.attr,
    lowerDetectionLimit = 1,
    expressionFamily = negbinomial.size())
save(with.spikes.sandberg.transcript.rpc, file="with-spikes.sandberg.transcript.rpc.RData")

quit("no")

#
# run rel2abs, convert resulting matrices into CDS, and store as RData files
#

qsub -P sage rel2abs.sh sandberg.transcript.tpm.RData sandberg.transcript.rpc.matrix.nbt-3rd-submission.gz
qsub -P sage rel2abs.sh sandberg.gene.tpm.RData sandberg.gene.rpc.matrix.nbt-3rd-submission.gz
qsub -P sage rel2abs.sh with-spikes.sandberg.gene.tpm.RData with.spikes.sandberg.gene.rpc.matrix.nbt-3rd-submission.gz

R

sandberg.transcript.rpc <- as.matrix(read.table(gzfile("sandberg.transcript.rpc.matrix.nbt-3rd-submission.gz"), header=T, row.names=1))

sandberg.gene.rpc <- as.matrix(read.table(gzfile("sandberg.gene.rpc.matrix.nbt-3rd-submission.gz"), header=T, row.names=1))

with.spikes.sandberg.gene.rpc <- as.matrix(read.table(gzfile("with.spikes.sandberg.gene.rpc.matrix.nbt-3rd-submission.gz"), header=T, row.names=1))

sandberg.transcript.rpc <- newCellDataSet(
    as(as.matrix(sandberg.transcript.rpc), "sparseMatrix"),
    phenoData = cell.attr[!is.na(cell.attr$embryo),],
    featureData = transcript.attr,
    lowerDetectionLimit = 1,
    expressionFamily = negbinomial.size())

sandberg.gene.rpc <- newCellDataSet(
    as(as.matrix(sandberg.gene.rpc), "sparseMatrix"),
    phenoData = cell.attr[!is.na(cell.attr$embryo),],
    featureData = gene.attr,
    lowerDetectionLimit = 1,
    expressionFamily = negbinomial.size())

with.spikes.sandberg.gene.rpc <- newCellDataSet(
    as(as.matrix(with.spikes.sandberg.gene.rpc), "sparseMatrix"),
    phenoData = with.spikes.cell.attr,
    featureData = with.spikes.gene.attr,
    lowerDetectionLimit = 1,
    expressionFamily = negbinomial.size())

save(sandberg.transcript.rpc, file="sandberg.transcript.rpc.nbt-3rd-submission.RData")
save(sandberg.gene.rpc, file="sandberg.gene.rpc.nbt-3rd-submission.RData")
save(with.spikes.sandberg.gene.rpc, file="with.spikes.sandberg.gene.rpc.nbt-3rd-submission.RData")

quit("no")

