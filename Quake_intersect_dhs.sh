# #!/bin/bash

#intersect peak file from three Lung DHS replicates:
bedtools intersect -a wgEncodeUwDnaseLungC57bl6MAdult8wksPkRep1.narrowPeak -b wgEncodeUwDnaseLungC57bl6MAdult8wksPkRep2.narrowPeak -sorted > tmp
bedtools intersect -a tmp -b wgEncodeUwDnaseLungC57bl6MAdult8wksPkRep3.narrowPeak -sorted >  Lung_EncodeIntersect_PKs.bed
rm tmp

# On trapnell cluster:
bedtools getfasta -fi mm9.fa -bed Lung_EncodeIntersect_PKs.bed -fo Lung_EncodeIntersect_PKs.fasta
fimo --text --parse-genomic-coord motif_databases/JASPAR_CORE_2014_vertebrates.meme Lung_EncodeIntersect_PKs.fasta > Lung_EncodeIntersect_PKs_fimo_output.txt

# bedtools getfasta -fi /net/trapnell/vol1/genomes/human/hg19/fasta/hg19.fa -bed Lung_EncodeIntersect_PKs.bed -fo HSMMtube_EncodeUniform_PKs.fasta
# fimo --text --parse-genomic-coord motif_databases/JASPAR_CORE_2014_vertebrates.meme HSMMtube_EncodeUniform_PKs.fasta > HSMMtube_EncodeUniform_PKs_fimo_output.txt

python dhs_fimo_to_bed.py motif_databases/JASPAR_CORE_2014_vertebrates.meme Lung_EncodeIntersect_PKs_fimo_output.txt > Lung_EncodeIntersect_PKs_fimo_output.bed
# python dhs_fimo_to_bed.py motif_databases/JASPAR_CORE_2014_vertebrates.meme HSMMemb_EncodeUniform_PKs_fimo_output.txt > HSMMemb_EncodeUniform_PKs_fimo_output.bed

cat *_fimo_output.bed > fimo_hits.bed
Rscript filter_fimo.q fimo_hits.bed fimo_hits_filtered.bed

awk '{print $4}' fimo_hits.bed | sort | uniq > tf_motif_list.txt

#convert genecode gtf file to bed file with bedops:
gtf2bed < gencode.vM6.annotation.gtf > sorted-gencode.vM6.annotation.gtf.bed
gtf2bed < gencode.vM1.annotation.gtf > sorted-gencode.vM1.annotation.gtf.bed

bedtools intersect -a sorted-gencode.vM1.annotation.gtf.bed -b fimo_hits_filtered.bed > fimo_hits_filtered.gencode.bed

sortBed -i fimo_hits_filtered.bed > fimo_hits_filtered.bed.sorted

#make the correct gmt file: 
closestBed -d -a fimo_hits_filtered.bed.sorted -b sorted-gencode.vM1.annotation.gtf.bed | awk '$NF >-5000 && $NF < 5000 { print $10    , $4; }'| sort | uniq > Lung_JASPAR_5kb_hits_olap.txt

#change all names into gene short name:
Rscript ensemble2gene_short_name.R Lung_JASPAR_5kb_hits_olap.txt Lung_JASPAR_5kb_hits_olap_gene_short_name.txt
./make_gmt.py Lung_JASPAR_5kb_hits_olap_gene_short_name.txt tf_motif_list.txt Lung_JASPAR_5kb_hits_olap.gmt
