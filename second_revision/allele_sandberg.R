library(monocle) 
# library(devtools)
# load_all('~/Projects/monocle-dev')
# load_all('~/Dropbox (Personal)/Projects/monocle-dev/')
library(stringr)
library(dplyr)
library(pheatmap)
library(matrixStats)
library(reshape2)
library(zoo)

#########################

#source allele related scripts: 
source('./allele_helper.R')

script.dir <- dirname(sys.frame(1)$ofile)
source(paste(script.dir, '/allele_helper.R', sep = ''))

# Set global ggplot2 properties for making print-scaled PDF panels
theme_set(theme_bw(base_size=6))

time_colors = c("#D7191C", "#FDAE61", "#ABD9E9", "#2C7BB6")

#/net/trapnell/vol1/jspacker/sandberg/monocle/sandberg.transcript.rpc.nbt-3rd-submission.RData
#/net/trapnell/vol1/jspacker/sandberg/monocle/with-spikes.sandberg.transcript.rpc.nbt-3rd-submission.RData
#/net/trapnell/vol1/jspacker/sandberg/monocle/with-spikes.sandberg.transcript.read.counts.RData

# load("./data/sandberg.transcript.rpc.nbt-3rd-submission.RData")
load('sandberg.gene.rpc.nbt-3rd-submission.RData')
sandberg <- sandberg.gene.rpc

# row.names(sandberg) <- str_replace(row.names(sandberg), 'ENSMUST', 'ENSMUSG')

fData(sandberg)$allele_id <- row.names(fData(sandberg))
fData(sandberg)$allele_id <- str_split_fixed(fData(sandberg)$allele_id, "_", 2)[,1]
fData(sandberg)$gene <- fData(sandberg)$allele_id
fData(sandberg)$gene_id <- fData(sandberg)$gene

sandberg <- estimateSizeFactors(sandberg)
sandberg <- detectGenes(sandberg, min_expr=1)
sandberg <- estimateDispersions(sandberg)
 
#select only protein coding and lincRNA for downstream analysis: 
gene_metadata_vM9 <- read.table('./gencode.vM9.annotation.gtf.processed')
row.names(gene_metadata_vM9) <- gene_metadata_vM9$V1
valid_genes <- row.names(subset(gene_metadata_vM9[fData(sandberg)$gene, ], as.character(V2) %in% c('protein_coding', 'lincRNA')))
valid_genes_ids <- row.names(subset(fData(sandberg), gene %in% valid_genes))
sandberg <- sandberg[valid_genes_ids, ]

fData(sandberg)$gene_short_name <- fData(sandberg)$symbol
pData(sandberg)$stage <- factor(pData(sandberg)$stage, levels = c("zygote",
                                  "early-2-cell", 
                                "mid-2-cell",
                                "late-2-cell",
                                "4-cell",
                                "8-cell",
                                "16-cell",
                                "early-blast",
                                "mid-blast",
                                "late-blast"))

qplot(stage, fill=sex, geom="bar", data=pData(sandberg))

target_genes <-row.names(subset(fData(sandberg), symbol %in% c("Gnai3")))

target_genes <- row.names(fData(sandberg)[1:20,])

iso_switch_test_res <- allelic_bias_test(sandberg[target_genes,], fullModelFormulaStr="~stage", cores=detectCores() / 2)
closeAllConnections()
time_de_res <- differentialGeneTest(sandberg[row.names(subset(fData(sandberg), num_cells_expressed >= 15)),], fullModelFormulaStr="~stage", cores=detectCores() / 2, verbose=T)
sig_genes <- subset(time_de_res, qval < 0.01)$gene_id

script.dir <- dirname(sys.frame(1)$ofile)
source(paste(script.dir, '/allele_fig6ab.R', sep = ''))

spike_ids <-  c("RNA_SPIKE_1_EC2_Gene_Bank_AE000151",                             
                "RNA_SPIKE_2_EC12_Gene_Bank_D86239",                               
                "RNA_SPIKE_3_EC3_Gene_Bank_AE000151",                                
                "RNA_SPIKE_4_EC15_ACCESSION_AE000185_U00096",                        
                "RNA_SPIKE_5_EC17_Accession_Number_V00307_J01654",                   
                "RNA_SPIKE_6_EC13_Gene_Bank_M34825",                                 
                "RNA_SPIKE_7_EC18_Accession_Number_M15263",                           
                "RNA_SPIKE_8_EC5_Accession_Number_AE000184_U00096")

load('with.spikes.sandberg.gene.rpc.nbt-3rd-submission.RData')
sandberg_rpc <- with.spikes.sandberg.gene.rpc
# row.names(sandberg_rpc) <- str_replace(row.names(sandberg_rpc), 'ENSMUST', 'ENSMUSG')

# fData(sandberg_rpc)$gene_id <- str_split_fixed(row.names(fData(sandberg_rpc)), "_", 2)[,1]
fData(sandberg_rpc)$allele_id <- row.names(fData(sandberg_rpc))
fData(sandberg_rpc)$allele_id <- str_split_fixed(fData(sandberg_rpc)$allele_id, "_", 2)[,1]
fData(sandberg_rpc)$gene <- fData(sandberg_rpc)$allele_id
fData(sandberg_rpc)$gene_id <- fData(sandberg_rpc)$gene

sandberg_rpc <- detectGenes(sandberg_rpc)
sandberg_rpc <- estimateSizeFactors(sandberg_rpc)
sandberg_rpc <- estimateDispersions(sandberg_rpc)

#select only protein coding and lincRNA for downstream analysis: 
valid_genes <- row.names(subset(gene_metadata_vM9[fData(sandberg_rpc)$gene, ], as.character(V2) %in% c('protein_coding', 'lincRNA')))
valid_genes_ids <- row.names(subset(fData(sandberg_rpc), gene %in% valid_genes))
sandberg_rpc <- sandberg_rpc[c(valid_genes_ids, spike_ids), ]

spike_rpcs <- sandberg_rpc[spike_ids,]
spike_rpcs <- estimateSizeFactors(spike_rpcs)



load('with-spikes.sandberg.gene.read.counts.RData')
 
sandberg_count <- with.spikes.sandberg.gene.read.counts
# row.names(sandberg_count) <- str_replace(row.names(sandberg_count), 'ENSMUST', 'ENSMUSG')
fData(sandberg_count)$allele_id <- row.names(fData(sandberg_count))
fData(sandberg_count)$allele_id <- str_split_fixed(fData(sandberg_count)$allele_id, "_", 2)[,1]
fData(sandberg_count)$gene <- fData(sandberg_count)$allele_id
fData(sandberg_count)$gene_id <- fData(sandberg_count)$gene

sandberg_count <- detectGenes(sandberg_count)
sandberg_count <- estimateSizeFactors(sandberg_count)
sandberg_count <- estimateDispersions(sandberg_count)

valid_genes <- row.names(subset(gene_metadata_vM9[fData(sandberg_count)$gene, ], as.character(V2) %in% c('protein_coding', 'lincRNA')))
valid_genes_ids <- row.names(subset(fData(sandberg_count), gene %in% valid_genes))
sandberg_count <- sandberg_count[c(valid_genes_ids, spike_ids), ]

spike_counts <- sandberg_count[spike_ids,]
spike_counts <- estimateSizeFactors(spike_counts)

pData(sandberg_count)$Spike_factor <- pData(spike_counts)$Size_Factor

#qplot(rowMeans(exprs(spike_rpcs))

pData(sandberg_rpc)$Spike_factor <- pData(spike_rpcs)$Size_Factor

autosomal_informative_genes <- as.character(subset(fData(sandberg_rpc), is.na(chr) == FALSE & chr %in% c("X", "Y") == FALSE & is.allele.informative == 1)$gene)


gencode_gene_metadata <- read.delim("gencode_biotypes.txt")

#locus_cols <- str_split_fixed(gencode_gene_metadata$locus, ":|-", 3) 
gencode_gene_metadata$gene_id <- str_split_fixed(gencode_gene_metadata$gene_id, "\\.", 2)[,1]


allowed_gene_types <- c("protein_coding", "lincRNA")
ribosomal_gene_ids <- subset(fData(sandberg_rpc), grepl("^Rps|^Rpl", symbol))$gene

autosomal_gene_records <- row.names(subset(fData(sandberg_rpc), gene %in% autosomal_informative_genes 
                                           & #gene %in% subset(gencode_gene_metadata, biotype %in% allowed_gene_types)$gene_id &
                                             gene %in% ribosomal_gene_ids == FALSE))
#autosomal_gene_records <- setdiff(autosomal_gene_records, ribosomal_gene_ids)

disp_calibration_sim<- sim_escape(sandberg_count[autosomal_gene_records,], 
                                        formulaStr="~embryo + Spike_factor", 
                                        allele_thresh=0.1,
                                        dispersion=0,
                                        nsims=100,
                                        cores=detectCores() / 2)
qplot(mean_expression, link_phi, data=subset(disp_calibration_sim, link_phi < 1e4), log="xy") + geom_smooth()

source(paste(script.dir, '/allele_fig6e.R', sep = ''))
source(paste(script.dir, '/allele_fig6cd.R', sep = ''))
source(paste(script.dir, '/allele_fig6f.R', sep = ''))

################################

# Xist_maternal <- "ENSMUSG00000086503_M"
# Xist_paternal <- "ENSMUSG00000086503_P"

# Jpx_maternal <- "ENSMUSG00000097571_M"
# Jpx_paternal <- "ENSMUSG00000097571_P"

# paternal_ncrnas <- Deng_sexed_cells[c(Xist_paternal, Xist_maternal, Jpx_paternal, Jpx_maternal), pData(Deng_sexed_cells)$Deng_Sex == "Female"]
# fData(paternal_ncrnas)$gene_short_name = fData(paternal_ncrnas)$symbol

# qplot(as.numeric(exprs(Deng_sexed_cells[c(Xist_paternal), pData(Deng_sexed_cells)$Deng_Sex == "Female"])),
#       as.numeric(exprs(Deng_sexed_cells[c(Jpx_paternal), pData(Deng_sexed_cells)$Deng_Sex == "Female"])))

# monocle:::plot_coexpression_matrix(paternal_ncrnas, rowgenes=c("Xist"), colgenes=c("Jpx"))


all_informative_genes <- as.character(subset(fData(Deng_sexed_cells), is.allele.informative == 1)$gene)

all_informative_gene_records <- row.names(subset(fData(Deng_sexed_cells), gene %in% all_informative_genes 
  & gene %in% subset(gencode_gene_metadata, biotype %in% allowed_gene_types)$gene_id))

all_informative_sim_res<- sim_escape(Deng_sexed_cells[all_informative_gene_records,], 
                                        formulaStr="~stage + Deng_Sex + num_genes_expressed", 
                                        nsims=100,
                                        cores=detectCores() / 2)

nrow(subset(all_informative_sim_res, (maternal_expr_qval < 0.05 | paternal_expr_qval < 0.05) &  monoallele_qval < 0.05))




escape_sim_res <- rbind(sandberg_escape_sim_res_4_cell, sandberg_escape_sim_res_16_cell, sandberg_escape_sim_res_early_blast)

escape_sim_res$stage <- factor(escape_sim_res$stage, levels=c("4-cell", "16-cell", "early-blast"))
escape_sim_res$gene_id <- row.names(escape_sim_res)


allelic_bias_genes <- subset(escape_sim_res, allelic_bias_qval < 0.05)$gene_id

subset(fData(Deng_sexed_cells), gene %in% allelic_bias_genes)

#nrow(all_informative_sim_res, )



#subset(fData(sandberg_escape_cds), gene_id %in% subset(sandberg_escape_sim_res, paternal_expr_qval < 0.05)$gene_id)


# escape_test_res <- allelic_bias_test(sandberg_escape_cds, fullModelFormulaStr="~embryo", cores=8)
# escape_test_res <- merge(distinct(fData(sandberg)[,c("gene_id", "symbol", "chr")]), escape_test_res, by.x="gene_id", by.y="row.names")
# row.names(escape_test_res) <- escape_test_res$gene_id
# escape_test_res$qval <- 1.0
# escape_test_res$qval[escape_test_res$status == "OK"] <- p.adjust(subset(escape_test_res, status == "OK")$pval, method="fdr")

# x_escape_candidate_gene_ids <- setdiff(x_escape_candidate_gene_ids, subset(escape_test_res, qval < 0.05)$gene_id)

# sandberg_escape_cds <- Deng_sexed_cells[row.names(subset(fData(Deng_sexed_cells), gene_id %in% x_escape_candidate_gene_ids)),]

save.image('./RData/analysis_allele_switch.RData')

