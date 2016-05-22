################################################################################################
###the following code is used to calculate all the number used in the paper: 
#the following script is missed from the first submission. It was later added 
#in the second submission and originally included as part of the `remaining_review_figures.R` 
################################################################################################

load('./RData/gen_supplementary_figure.RData')
load('./RData/gen_lung_figures.RData')
load('./RData/gen_shalek_figures.RData')
load('./RData/deg_benchmark_analysis.RData')
load('./RData/cmpr_three_packages.RData')
load('./RData/cmpr_three_packages.RData')

#1. consensus between five tools: 
length(readcount_overlap)
#[1] 946
length(readcount_union)
#[1] 13364
length(abs_overlap)
#[1] 1649
length(abs_union)
#[1] 7209

#2. false discovery rate for the DESEq: 
fdr_sensitivity_df
subset(fdr_sensitivity_df, Type == 'DESeq1')
#read_count transcript_counts
#fpr: 0.1490958 0.0490389
#sensitivity: 0.2721133 0.3959332

subset(fdr_sensitivity_df, Type == 'Monocle' & data_type %in% c('Spikein transcripts', 'Read counts'))[, c('data_type', 'fpr')]
#read_counts transcript_counts: 
#fpr: 0.04736089   0.02127532 
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

