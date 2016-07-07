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
load('./RData/umi_normalization.RData')

#1. consensus between five tools: 
length(readcount_overlap)
#2035
length(readcount_union)
#6235
#2035 / 6235
#0.3263833
length(abs_overlap)
#2437
length(abs_union)
#4220
# > 2437 / 4220
# [1] 0.5774882
length(census_overlap)
#[1] 2119
length(census_union)
#[1] 3875
# > 2119 / 3875
# [1] 0.5468387

#2. false discovery rate for the DESEq: 
fdr_sensitivity_cal(default_deseq2_p, mode_size_norm_permutate_ratio_by_geometric_mean, 'fdr')
# [1] 0.5821187
fdr_sensitivity_cal(abs_default_deseq2_p, mode_size_norm_permutate_ratio_by_geometric_mean, 'fdr')
# [1] 0.2618625
fdr_sensitivity_cal(default_deseq2_p, mode_size_norm_permutate_ratio_by_geometric_mean, 'sensitivity')
# [1] 0.902727
fdr_sensitivity_cal(abs_default_deseq2_p, mode_size_norm_permutate_ratio_by_geometric_mean, 'sensitivity')
# [1] 0.935901
fdr_sensitivity_cal(mc_default_deseq2_p, mode_size_norm_permutate_ratio_by_geometric_mean, 'fdr')
# [1] 0.2579564
fdr_sensitivity_cal(mc_default_deseq2_p, mode_size_norm_permutate_ratio_by_geometric_mean, 'sensitivity')
# [1] 0.8742599

fdr_sensitivity_cal(monocle_p_readcount, mode_size_norm_permutate_ratio_by_geometric_mean, 'fdr')
# [1] 0.494898
fdr_sensitivity_cal(new_abs_size_norm_monocle_p_ratio_by_geometric_mean, mode_size_norm_permutate_ratio_by_geometric_mean, 'fdr')
# [1] 0.09314034
fdr_sensitivity_cal(mc_mode_size_norm_permutate_ratio_by_geometric_mean, mode_size_norm_permutate_ratio_by_geometric_mean, 'fdr')
# [1] 0.04165144

#3. number of enriched TFs: 
nrow(valid_hyper_df)
# [1] 82
length(branch_motif_Tfs_id)
# [1] 17

#number of branch genes
# length(row.names(subset(weihgted_relative_abs_AT12_cds_subset_all_gene, qval <= 0.05)))
#[1] 1720
length(row.names(subset(weihgted_mc_AT12_cds_subset_all_gene, qval <= 0.05)))
#[1] 1685
length(row.names(subset(ko_branching_genes, qval < 0.05)))
#[1] 1811
# length(row.names(subset(golgi_branching_genes, qval < 0.05)))

#number of genes related to interferon signaling: 
table(Shalek_abs_subset_ko_LPS_heatmap_annotations$annotation_row$Cluster)
#   1   2   3   4   5   6
# 575 219 185 532 153 144

#debug the permutation test: 
#the error comes from the change of the permutation algorithm. When the original approach used. we get good results 

#reproduce the previous result: 
#we got clear explannation on the changes of genes but also confirm our algorithm is correct: 
# std_HSMM_myo_pseudotime_res_ori <- differentialGeneTest((std_HSMM[, ]), relative_expr = F, 
#                                                         fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", #log10(Total_mRNAs) + spike_total_mRNAs
#                                                         reducedModelFormulaStr = "~1", cores = detectCores())


# HSMM_myo_size_norm_res <- differentialGeneTest(HSMM_myo[, ], 
#                                                fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", #log10(Total_mRNAs) + spike_total_mRNAs
#                                                reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

#update the supplementary data: 
lung_go_enrichment <- read.delim('./supplementary_data/lung_hyper_df.xls', header = T)
ko_reactome_enrichment <- read.delim('./supplementary_data/ko_hyper_df.xls', header = T)
# golgiplug_reactome_enrichment <- read.delim('./supplementary_data/golgiplug_hyper_df.xls', header = T)

lung_go_enrichment$cluster_id <- as.character(lung_go_enrichment$cluster_id)
lung_go_enrichment$cluster_id <- revalue(lung_go_enrichment$cluster_id, c("1" = "cluster_3", "2" = "cluster_4", "3" = "cluster_6", "4" = "cluster_5", "5" = "cluster_1", "6" = "cluster_2"))

ko_reactome_enrichment$cluster_id <- as.character(ko_reactome_enrichment$cluster_id)
ko_reactome_enrichment$cluster_id <- revalue(ko_reactome_enrichment$cluster_id, c("1" = "cluster_4", "2" = "cluster_6", "3" = "cluster_5", "4" = "cluster_1", "5" = "cluster_6", "6" = "cluster_3"))

# golgiplug_reactome_enrichment$cluster_id <- as.character(golgiplug_reactome_enrichment$cluster_id)
# golgiplug_reactome_enrichment$cluster_id <- revalue(golgiplug_reactome_enrichment$cluster_id, c("1" = "cluster_5", "2" = "cluster_4", "3" = "cluster_6", "4" = "cluster_2", "5" = "cluster_3", "6" = "cluster_1"))

write.table(lung_go_enrichment, sep = '\t', quote = F, row.names = F, file = 'lung_go_enrichment.txt')
write.table(ko_reactome_enrichment, sep = '\t', quote = F, row.names = F, file = 'ko_reactome_enrichment.txt')
# write.table(golgiplug_reactome_enrichment, sep = '\t', quote = F, row.names = F, file = 'golgiplug_reactome_enrichment')

#save it to file: 
# permutation_res <- data.frame(spikein_free = mc_mode_size_norm_permutate_ratio_by_geometric_mean, spikein = mode_size_norm_permutate_ratio_by_geometric_mean, read_count = readcount_permutate_pval[names(std_permutate_pval)], fpkm = std_permutate_pval)
# write.table(permutation_res, file = './permutation.txt', sep = '\t', row.names = T, quote = F)

permutation_res <- data.frame(E14_18 = mode_size_norm_permutate_ratio_by_geometric_mean, E16_18 = mode_size_norm_permutate_ratio_by_geometric_mean)
write.table(permutation_res, file = './permutation.txt', sep = '\t', row.names = T, quote = F)


#number of cells in the tree downsampling figure: 
load('./RData/analysis_cell_downsampling.RData')
unlist(lapply(c(1, 6, 8, 10, 14, 18, 19, 23, 27, 30, 33, 36), function(x) ncol(cds_downsampled_cells_ordered[[x]])))
# Samples Samples Samples Samples Samples Samples Samples Samples Samples Samples
#      37      78     117     150     196     230     273     304     325     342
# Samples Samples
#     371     387

#figure 5: 
unique(subset(iso_switch_test_res, qval < 0.01)$gene_short_name

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

#PSI
cluster_summaries

dim(iso_switch_test_res$qval < 0.01)
dim(count_iso_switch_test_res$qval < 0.1)

#figure 6: 
load('./RData/analysis_allele_switch.RData')
fig6e_a
fig6e_b
fig6e_c

fig6f_a
fig6f_b
fig6f_c

#numbers in supplementary texts: 

#number of cells selected for analysis: 
absolute_cds
#183

Shalek_abs_subset_ko_LPS #check the code (analysis_shalek) for the qval threshold for making the tree 
#461

HSMM_myo
#168

HSMM_myo #in figure 5 comparing to the supplementary file is a little different: 174 cells
UMI_cds_76_55
#131


