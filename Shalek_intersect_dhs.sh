#!/bin/bash

# On trapnell cluster:
bedtools getfasta -fi /net/trapnell/vol1/home/xqiu/absolute_Rtranscript/DHS_analysis/mm9.fa -bed mouse_ATAC_seq_union.bed.gz -fo mouse_encode_dnase_lymphatic_cell_union.fasta
fimo --text --parse-genomic-coord /net/trapnell/vol1/home/xqiu/absolute_Rtranscript/DHS_analysis/motif_databases/JASPAR_CORE_2014_vertebrates.meme mouse_encode_dnase_lymphatic_cell_union.fasta > mouse_encode_dnase_lymphatic_cell_union_PKs_fimo_output.txt

# bedtools getfasta -fi /net/trapnell/vol1/genomes/human/hg19/fasta/hg19.fa -bed Lung_EncodeIntersect_PKs.bed -fo HSMMtube_EncodeUniform_PKs.fasta
# fimo --text --parse-genomic-coord motif_databases/JASPAR_CORE_2014_vertebrates.meme HSMMtube_EncodeUniform_PKs.fasta > HSMMtube_EncodeUniform_PKs_fimo_output.txt

python dhs_fimo_to_bed.py ../motif_databases/JASPAR_CORE_2014_vertebrates.meme mouse_encode_dnase_lymphatic_cell_union_PKs_fimo_output.txt > mouse_encode_dnase_lymphatic_cell_union_PKs_fimo_output.bed
# python dhs_fimo_to_bed.py motif_databases/JASPAR_CORE_2014_vertebrates.meme HSMMemb_EncodeUniform_PKs_fimo_output.txt > HSMMemb_EncodeUniform_PKs_fimo_output.bed

cat *_fimo_output.bed > fimo_hits.bed
Rscript filter_fimo.q fimo_hits.bed fimo_hits_filtered.bed

awk '{print $4}' fimo_hits.bed | sort | uniq > tf_motif_list.txt

# gtf2bed < gencode.vM6.annotation.gtf > sorted-gencode.vM6.annotation.gtf.bed
gtf2bed < gencode.vM1.annotation.gtf > sorted-gencode.vM1.annotation.gtf.bed

#bedtools intersect -a sorted-gencode.vM1.annotation.gtf.bed -b fimo_hits_filtered.bed > fimo_hits_filtered.gencode.bed

sortBed -i fimo_hits_filtered.bed > fimo_hits_filtered.bed.sorted

#use the approach I used for the quake data:
closestBed -d -a fimo_hits_filtered.bed.sorted -b sorted-gencode.vM1.annotation.gtf.bed | awk '$NF >-5000 && $NF < 5000 { print $10    , $4; }'| sort | uniq > DC_JASPAR_5kb_hits_olap.txt

#change all names into gene short name: 
Rscript ensemble2gene_short_name.R DC_JASPAR_5kb_hits_olap.txt DC_JASPAR_5kb_hits_olap_gene_short_name.txt
./make_gmt.py DC_JASPAR_5kb_hits_olap_gene_short_name.txt tf_motif_list.txt DC_JASPAR_5kb_hits_olap.gmt
