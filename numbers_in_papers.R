################################################################################################
###the following code is used to calculate all the number used in the paper: 
#the following script is missed from the first submission. It was later added 
#in the second submission and originally included as part of the `remaining_review_figures.R` 
################################################################################################

load('./RData/gen_lung_figures_mc.RData')
load('./RData/gen_shalek_figures.RData')
load('./RData/deg_benchmark_analysis.RData')
load('./RData/cmpr_three_packages.RData')
load('./RData/gen_supplementary_figure.RData')

#1. consensus between five tools: 
length(readcount_overlap)
#[1] 1324
length(readcount_union)
#[1] 14386
length(abs_overlap)
#[1] 1957
length(abs_union)
#[1] 5946

#2. false discovery rate for the DESEq: 
fdr_sensitivity_df
subset(fdr_sensitivity_df, Type == 'DESeq1')
#read_count transcript_counts
#fpr: 0.12194463 0.02481542
#sensitivity: 0.3019018 0.6183206

subset(fdr_sensitivity_df, Type == 'Monocle' & data_type %in% c('Spikein transcripts', 'Read counts'))[, c('data_type', 'fpr')]
#read_counts transcript_counts: 
#fpr: 0.018550392 0.008000911   
#3. number of enriched TFs: 
nrow(valid_hyper_df)
# 56

#number of branch genes
length(row.names(subset(weihgted_relative_abs_AT12_cds_subset_all_gene, qval <= 0.05)))
#[1] 1720
length(row.names(subset(ko_branching_genes, qval < 0.05)))
#[1] 1062
length(row.names(subset(golgi_branching_genes, qval < 0.05)))

#number of genes related to interferon signaling: 
table(Shalek_abs_subset_ko_LPS_heatmap_annotations$annotation_row$Cluster)

#  1   2   3   4   5   6
# 92 306 155 206 127 117

#debug the permutation test: 
#the error comes from the change of the permutation algorithm. When the original approach used. we get good results 

#reproduce the previous result: 
#we got clear explannation on the changes of genes but also confirm our algorithm is correct: 
std_HSMM_myo_pseudotime_res_ori <- differentialGeneTest((std_HSMM[, ]), relative_expr = F, 
                                                        fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", #log10(Total_mRNAs) + spike_total_mRNAs
                                                        reducedModelFormulaStr = "~1", cores = detectCores())


HSMM_myo_size_norm_res <- differentialGeneTest(HSMM_myo[, ], 
                                               fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", #log10(Total_mRNAs) + spike_total_mRNAs
                                               reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

#update the supplementary data: 
lung_go_enrichment <- read.delim('./supplementary_data/lung_hyper_df.xls', header = T)
ko_reactome_enrichment <- read.delim('./supplementary_data/ko_hyper_df.xls', header = T)
golgiplug_reactome_enrichment <- read.delim('./supplementary_data/golgiplug_hyper_df.xls', header = T)

lung_go_enrichment$cluster_id <- as.character(lung_go_enrichment$cluster_id)
lung_go_enrichment$cluster_id <- revalue(lung_go_enrichment$cluster_id, c("1" = "cluster_1", "2" = "cluster_6", "3" = "cluster_5", "4" = "cluster_4", "5" = "cluster_2", "6" = "cluster_3"))

ko_reactome_enrichment$cluster_id <- as.character(ko_reactome_enrichment$cluster_id)
ko_reactome_enrichment$cluster_id <- revalue(ko_reactome_enrichment$cluster_id, c("1" = "cluster_5", "2" = "cluster_6", "3" = "cluster_2", "4" = "cluster_4", "5" = "cluster_1", "6" = "cluster_3"))

golgiplug_reactome_enrichment$cluster_id <- as.character(golgiplug_reactome_enrichment$cluster_id)
golgiplug_reactome_enrichment$cluster_id <- revalue(golgiplug_reactome_enrichment$cluster_id, c("1" = "cluster_5", "2" = "cluster_4", "3" = "cluster_6", "4" = "cluster_2", "5" = "cluster_3", "6" = "cluster_1"))

write.table(lung_go_enrichment, sep = '\t', quote = F, row.names = F, file = 'lung_go_enrichment')
write.table(ko_reactome_enrichment, sep = '\t', quote = F, row.names = F, file = 'ko_reactome_enrichment')
write.table(golgiplug_reactome_enrichment, sep = '\t', quote = F, row.names = F, file = 'golgiplug_reactome_enrichment')

#save it to file: 
permutation_res <- data.frame(spikein_free = mc_mode_size_norm_permutate_ratio_by_geometric_mean, spikein = mode_size_norm_permutate_ratio_by_geometric_mean, read_count = readcount_permutate_pval[names(std_permutate_pval)], fpkm = std_permutate_pval)
write.table(permutation_res, file = './permutation.txt', sep = '\t', row.names = T, quote = F)

#number of cells in the tree downsampling figure: 
load('./RData/analysis_cell_downsampling.RData')
unlist(lapply(c(1, 6, 8, 10, 14, 18, 19, 23, 27, 30, 33, 36), function(x) ncol(cds_downsampled_cells_ordered[[x]])))
# Samples Samples Samples Samples Samples Samples Samples Samples Samples Samples
#      37      78     117     150     196     230     273     304     325     342
# Samples Samples
#     371     387

#figure 5: 
fData(HSMM_myo_isoform[row.names(subset(fData(HSMM_myo_isoform), 
                                   gene_short_name %in% c("ACTA1", "ACTA2", "ACTB", "ACTG1", "ACTG2", "ACTC1") & num_cells_expressed >= 15)),])
#                              gene_id gene_short_name        biotype                     tss_id num_cells_expressed
# ENST00000366684.3 ENSG00000143632.10           ACTA1 protein_coding ENSG00000143632.10_TSS9173                  18
# ENST00000224784.6  ENSG00000107796.8           ACTA2 protein_coding ENSG00000107796.8_TSS12411                 107
# ENST00000290378.4  ENSG00000159251.6           ACTC1 protein_coding ENSG00000159251.6_TSS32091                 164
# ENST00000557860.1  ENSG00000159251.6           ACTC1 protein_coding ENSG00000159251.6_TSS32092                 168
# ENST00000560563.1  ENSG00000159251.6           ACTC1 protein_coding ENSG00000159251.6_TSS32093                  63
# ENST00000331925.2  ENSG00000184009.5           ACTG1 protein_coding ENSG00000184009.5_TSS44058                 174
# ENST00000576209.1  ENSG00000184009.5           ACTG1 protein_coding ENSG00000184009.5_TSS44060                 173
# ENST00000464611.1  ENSG00000075624.9            ACTB protein_coding ENSG00000075624.9_TSS85735                 165
# ENST00000331789.5  ENSG00000075624.9            ACTB protein_coding ENSG00000075624.9_TSS85736                 174

table(cutree(iso_switch_heatmap$tree_row, 6))
#  1  2  3  4  5  6 
# 44 58 29 15 21 56 

iso_names <- names(cutree(iso_switch_heatmap$tree_row, 6))
length(unique(fData(HSMM_myo_isoform)[iso_names[cutree(iso_switch_heatmap$tree_row, 6) == 1], 'gene_short_name']))
length(unique(fData(HSMM_myo_isoform)[iso_names[cutree(iso_switch_heatmap$tree_row, 6) == 2], 'gene_short_name']))
length(unique(fData(HSMM_myo_isoform)[iso_names[cutree(iso_switch_heatmap$tree_row, 6) == 3], 'gene_short_name']))
length(unique(fData(HSMM_myo_isoform)[iso_names[cutree(iso_switch_heatmap$tree_row, 6) == 4], 'gene_short_name']))
length(unique(fData(HSMM_myo_isoform)[iso_names[cutree(iso_switch_heatmap$tree_row, 6) == 5], 'gene_short_name']))
length(unique(fData(HSMM_myo_isoform)[iso_names[cutree(iso_switch_heatmap$tree_row, 6) == 6], 'gene_short_name']))





dim(iso_switch_test_res$qval < 0.01)
dim(count_iso_switch_test_res$qval < 0.1)

#figure 6: 
1-	nrow(subset(genes_interval_test_df, value > up)) / nrow(genes_interval_test_df)
#315/ 6629

nrow(subset(genes_interval_test_df, value > up)) / nrow(genes_interval_test_df)
# 4502/ 21988

#numbers in supplementary texts: 

#number of cells selected for analysis: 
absolute_cds

Shalek_abs_subset_ko_LPS #check the code (analysis_shalek) for the qval threshold for making the tree 

HSMM_myo
HSMM_myo #in figure 5 comparing to the supplementary file is a little different: 174 cells

