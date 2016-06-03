# # load('./RData/prepare_lung_data.RData')
library(devtools)
load_all('~/Projects/monocle-dev')
 
library(xacHelper)
library(MAST)

# library(scde)
#########################run the following scripts in remote sever#############################

#generate the pvals from the statistical test (permutation based or from software: monocle/DESeq/SCDE)

#double check the cells selected are the same as we used in the paper:  identical(colnames(new_abs_cds_14_18), cells_used)
new_abs_cds_14_18 <- absolute_cds[, intersect(valid_cells, row.names(subset(pData(absolute_cds), Time %in% c('E18.5', 'E14.5'))))]
new_abs_cds_14_18 <- estimateSizeFactors(new_abs_cds_14_18) #calculate the size factor for performing the relative absolute expression tests

table(pData(new_abs_cds_14_18)$Time) #check the number of used cells at each Time Point

new_mc_cds_14_18 <- mc_adj_cds[, colnames(new_abs_cds_14_18)]
new_mc_cds_14_18 <- estimateSizeFactors(new_mc_cds_14_18)

new_std_cds_14_18 <- standard_cds[, colnames(new_abs_cds_14_18)]
# new_std_cds_14_18 <- estimateSizeFactors(new_std_cds_14_18) #FPKM values are already on relative scale

#prepare the readcount data for DESeq / SCDE: 
read_countdata <- read_countdata[row.names(new_abs_cds_14_18), colnames(new_abs_cds_14_18)] 

#create a cds for readcount data to perform the default DEG tests for DESeq and SCDE : 
count_cds <- newCellDataSet(read_countdata[row.names(new_abs_cds_14_18), colnames(new_abs_cds_14_18)],
                            phenoData = new("AnnotatedDataFrame", data = pData(new_abs_cds_14_18)),
                            featureData = new("AnnotatedDataFrame", data = fData(new_abs_cds_14_18)),
                            expressionFamily = negbinomial(),
                            lowerDetectionLimit = 1)

count_cds <- estimateSizeFactors(count_cds)

# # calculate the pval with the readcount with scde: (calculate the scde associate DEG test result LOCALLY) 
# std_scde_res_list <- scde_DEG(dir = NULL, count_cds = count_cds, DEG_attribute = 'Time', contrast = c('E14.5', 'E18.5'), n.cores = detectCores())

# #calculate the pval with the normalized transcripts with scde: 
# abs_scde_res_list <- scde_DEG(dir = NULL, count_cds = new_abs_cds_14_18, DEG_attribute = 'Time', contrast = c('E14.5', 'E18.5'), n.cores = detectCores(), normalize = T)
# abs_scde_res_list_no_normalize <- scde_DEG(dir = NULL, count_cds = new_abs_cds_14_18, DEG_attribute = 'Time', contrast = c('E14.5', 'E18.5'), n.cores = detectCores())

#load all other necessary packages: 
load_all_libraries()

# #perform  the stastical tests on the data: 
new_std_diff_test_res <- differentialGeneTest(new_std_cds_14_18[1:transcript_num, ], 
                                              fullModelFormulaStr = "~Time", 
                                              reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)  

new_size_norm_abs_diff_test_res <- differentialGeneTest(new_abs_cds_14_18[1:transcript_num, ], 
                                                        fullModelFormulaStr = "~Time", 
                                                        reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

new_size_norm_mc_diff_test_res <- differentialGeneTest(new_mc_cds_14_18[1:transcript_num, ], 
                                                       fullModelFormulaStr = "~Time", 
                                                       reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

new_cds_diff_test_res <- differentialGeneTest(count_cds[1:transcript_num, ], 
                                              fullModelFormulaStr = "~Time", 
                                              reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

new_cds_diff_test_res_no_relative <- differentialGeneTest(count_cds[1:transcript_num, ], 
                                                          fullModelFormulaStr = "~Time", 
                                                          reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)

#pseudotime test on lineage 2/3:    
std_pseudotime_test_lineage2_res <- differentialGeneTest(std_AT12_cds_subset_all_gene[1:transcript_num, pData(std_AT12_cds_subset_all_gene)$State %in% c(1, 2)], 
                                              fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                              reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)  
std_pseudotime_test_lineage3_res <- differentialGeneTest(std_AT12_cds_subset_all_gene[1:transcript_num, pData(std_AT12_cds_subset_all_gene)$State %in% c(1, 3)], 
                                                      fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                                      reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)
abs_pseudotime_test_lineage2_res <- differentialGeneTest(abs_AT12_cds_subset_all_gene[1:transcript_num, pData(abs_AT12_cds_subset_all_gene)$State %in% c(1, 2)], 
                                                fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                                reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)  
abs_pseudotime_test_lineage3_res <- differentialGeneTest(abs_AT12_cds_subset_all_gene[1:transcript_num, pData(abs_AT12_cds_subset_all_gene)$State %in% c(1, 3)], 
                                                      fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                                      reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)
mc_pseudotime_test_lineage2_res <- differentialGeneTest(mc_AT12_cds_subset_all_gene[1:transcript_num, pData(mc_AT12_cds_subset_all_gene)$State %in% c(1, 2)], 
                                                fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                                reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)    
mc_pseudotime_test_lineage3_res <- differentialGeneTest(mc_AT12_cds_subset_all_gene[1:transcript_num, pData(mc_AT12_cds_subset_all_gene)$State %in% c(1, 3)], 
                                                      fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                                      reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)
Time_ori <- pData(new_abs_cds_14_18)$Time #Time as input for newCountDataSet

# #calculate the pval with the readcount with DESeq: 

# read_count_d <- newCountDataSet(round(read_countdata)[, ], Time_ori) 
# # std_dtable_pool_max_nbinomTest <- DESeq1_test(read_count_d, disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomTest', scale = T) 
# # row.names(std_dtable_pool_max_nbinomTest$dtalbe) <- std_dtable_pool_max_nbinomTest$dtalbe$id

# #calculate the pval with the normalized transcripts with DESeq: 
# abs_count_d <- newCountDataSet(round(t(t(exprs(new_abs_cds_14_18)) / sizeFactors(new_abs_cds_14_18))), (Time_ori)) #normalized the data by size factor
# # abs_dtable_pool_max_nbinomTest <- DESeq1_test(abs_count_d, disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomTest') 
# # row.names(abs_dtable_pool_max_nbinomTest$dtalbe) <- abs_dtable_pool_max_nbinomTest$dtalbe$id

# #DESeq glm: (GLM tests are more relevant to our software)
# std_dtable_pool_max_nbinomGLMTest <- DESeq1_test(read_count_d[1:transcript_num, ], disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomGLMTest', scale = T) 
# row.names(std_dtable_pool_max_nbinomGLMTest$dtalbe) <-  row.names(read_count_d[1:transcript_num, ])
# abs_dtable_pool_max_nbinomGLMTest <- DESeq1_test(abs_count_d[1:transcript_num, ], disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomGLMTest') 
# row.names(abs_dtable_pool_max_nbinomGLMTest$dtalbe) <- row.names(new_abs_cds_14_18[1:transcript_num, ])

# #permutation results for two group test on the FPKM value (gold standard for monocle on FPKM values with tobit model): 
# Time_order <- order(Time_ori) #order the E14.5 cell at the begining (ensure E14.5 cells number is 43 while E18.5 is 74)
# std_split_cds <- split(t(exprs(new_std_cds_14_18[1:transcript_num, Time_order])), col(t(exprs(new_std_cds_14_18[1:transcript_num, Time_order])), as.factor = T))
# std_fc <- esApply(new_std_cds_14_18[1:transcript_num, ], 1, mean_fc, grp0 = 'E14.5', grp1 = 'E18.5')
# std_split_fc <- split(t(std_fc), col(t(std_fc), as.factor = T))

# permuation_pval <- function (x, fc, alpha = table(pData(new_abs_cds_14_18)$Time)['E14.5'], beta = table(pData(new_abs_cds_14_18)$Time)['E18.5'], permutate_num = 10000,
#     return_fc = F)
# {
#     fc_vec <- rep(0, permutate_num)
#     if (is.na(fc))
#         return(1)
#     for (i in 1:permutate_num) {
#         x <- sample(x, length(x))
#         x_0 <- (mean(x[1:alpha]))
#         x_1 <- (mean(x[(alpha + 1):(alpha + beta)]))
#         mean_fc <- log2(x_1/x_0)
#         fc_vec[i] <- mean_fc
#     }
#     if (fc > 0)
#         pval <- sum(fc <= fc_vec)/length(fc_vec)
#     else pval <- sum(fc >= fc_vec)/length(fc_vec)
#     if (return_fc)
#         return(list(pval = pval, fc_vec = fc_vec))
#     else return(pval)
# }
# std_permutate_pval <- mcmapply(permuation_pval, std_split_cds, std_split_fc, mc.cores = detectCores()) #multiple cores 

# closeAllConnections()
# # permutation results for two group test on the readcounts value (gold standard for DESEeq / SCDE on readcount data): 
# readcount_sf_mat <- matrix(rep(sizeFactors(count_cds), nrow(count_cds)), nrow = nrow(count_cds), byrow = T)
# readcount_split_cds <- split(t(exprs(count_cds[1:transcript_num, Time_order]) / readcount_sf_mat[1:transcript_num, Time_order]), col(t(exprs(count_cds[1:transcript_num, Time_order])), as.factor = T))
# #readcount_split_cds <- split(t(exprs(count_cds[1:transcript_num, Time_order])), col(t(exprs(count_cds[1:transcript_num, Time_order])), as.factor = T))

# readcount_fc <- apply(exprs(count_cds[1:transcript_num, Time_order]) / readcount_sf_mat[1:transcript_num, Time_order], 1, mean_fc, grp0 = 'E14.5', grp1 = 'E18.5', grp = pData(count_cds)$Time) #valid_gene_id
# # readcount_fc <- esApply(count_cds[1:transcript_num, ], 1, mean_fc, grp0 = 'E14.5', grp1 = 'E18.5') #valid_gene_id
# readcount_split_fc <- split(t(readcount_fc), col(t(readcount_fc), as.factor = T))
# ## 

# readcount_permutate_pval <- mcmapply(permuation_pval, readcount_split_cds, readcount_split_fc, mc.cores = detectCores()) #multiple core
# closeAllConnections()

# # two-group permutation tests (the same as above)
# # readcount_permutate_pval <- permu_two_group_gen(count_cds[, Time_order], alpha = 43, beta = 60, grp0 = 'E14.5', grp1 = 'E18.5', group = pData(count_cds)$Time)

# # permutation results for two group test on the normalized absolute value (gold standard for DESEeq / SCDE on readcount data): 
# # mode_size_norm_permutate_ratio_by_geometric_mean <- cal_perm_pval_size_norm(new_abs_cds_14_18[1:transcript_num, Time_order])
# # mc_mode_size_norm_permutate_ratio_by_geometric_mean <- cal_perm_pval_size_norm(new_mc_cds_14_18[1:transcript_num, Time_order])

# #new implementation: 
# abs_sf_mat <- matrix(rep(sizeFactors(new_abs_cds_14_18), nrow(new_abs_cds_14_18)), nrow = nrow(new_abs_cds_14_18), byrow = T)
# abs_split_cds <- split(t(exprs(new_abs_cds_14_18[1:transcript_num, Time_order]) / abs_sf_mat[1:transcript_num, Time_order]), col(t(exprs(new_abs_cds_14_18[1:transcript_num, Time_order])), as.factor = T))
# abs_fc <- apply(exprs(new_abs_cds_14_18[1:transcript_num, Time_order]) / abs_sf_mat[1:transcript_num, Time_order], 1, mean_fc, grp0 = 'E14.5', grp1 = 'E18.5', grp = pData(new_abs_cds_14_18)$Time) #valid_gene_id
# abs_split_fc <- split(t(abs_fc), col(t(abs_fc), as.factor = T))
# mode_size_norm_permutate_ratio_by_geometric_mean <- mcmapply(permuation_pval, abs_split_cds, abs_split_fc, mc.cores = detectCores()) #multiple core
# closeAllConnections()

# mc_sf_mat <- matrix(rep(sizeFactors(new_mc_cds_14_18), nrow(new_mc_cds_14_18)), nrow = nrow(new_mc_cds_14_18), byrow = T)
# mc_split_cds <- split(t(exprs(new_mc_cds_14_18[1:transcript_num, Time_order]) / mc_sf_mat[1:transcript_num, Time_order]), col(t(exprs(new_mc_cds_14_18[1:transcript_num, Time_order])), as.factor = T))
# mc_fc <- esApply(new_mc_cds_14_18[1:transcript_num, ], 1, mean_fc, grp0 = 'E14.5', grp1 = 'E18.5') #valid_gene_id
# mc_fc <- apply(exprs(new_mc_cds_14_18[1:transcript_num, Time_order]) / mc_sf_mat[1:transcript_num, Time_order], 1, mean_fc, grp0 = 'E14.5', grp1 = 'E18.5', grp = pData(new_mc_cds_14_18)$Time) #valid_gene_id
# mc_split_fc <- split(t(mc_fc), col(t(mc_fc), as.factor = T))
# mc_mode_size_norm_permutate_ratio_by_geometric_mean <- mcmapply(permuation_pval, mc_split_cds, mc_split_fc, mc.cores = detectCores()) #multiple core
# closeAllConnections()

#add the DEG tests using edgeR / DESeq2: 
# edgeR_res <- edgeR_test(counts = exprs(count_cds)[1:transcript_num, ], glm = F)
# abs_edgeR_res <- edgeR_test(exprs(new_abs_cds_14_18[1:transcript_num, ]), group = Time_ori, glm = F)

# edgeR_res_glm <- edgeR_test(counts = exprs(count_cds)[1:transcript_num, ], glm = T)
# abs_edgeR_res_glm <- edgeR_test(exprs(new_abs_cds_14_18[1:transcript_num, ]), group = Time_ori, glm = T)
  
# deseq2_res <- DESeq2_deg(dir = NULL, count_cds[1:transcript_num, ], Time = Time_ori, pd = pData(count_cds))
# abs_deseq2_res <- DESeq2_deg(dir = NULL, new_abs_cds_14_18[1:transcript_num, ], Time = Time_ori, pd = pData(new_abs_cds_14_18))

#prepare to generate the data for create the precision/recall/F1 score data.frame: 
#new monocle_p: 
monocle_p <- new_std_diff_test_res[, 'pval'] 
names(monocle_p) <- row.names(new_std_diff_test_res)
#use cds: 
monocle_p_readcount <- new_cds_diff_test_res[, 'pval'] 
names(monocle_p_readcount) <- row.names(new_cds_diff_test_res)

#abs
new_abs_size_norm_monocle_p_ratio_by_geometric_mean <- new_size_norm_abs_diff_test_res[, 'pval']
names(new_abs_size_norm_monocle_p_ratio_by_geometric_mean) <- row.names(new_size_norm_abs_diff_test_res)
#mc 
new_mc_size_norm_monocle_p_ratio_by_geometric_mean <- new_size_norm_mc_diff_test_res[, 'pval']
names(new_mc_size_norm_monocle_p_ratio_by_geometric_mean) <- row.names(new_size_norm_mc_diff_test_res)
#deseq
default_deseq_p <- std_dtable_pool_max_nbinomGLMTest$dtalbe[, 'pval'] #std_dtable_pool_max_nbinomTest
names(default_deseq_p) <- row.names(std_dtable_pool_max_nbinomGLMTest$dtalbe) #std_dtable_pool_max_nbinomGLMTest
abs_default_deseq_p <- abs_dtable_pool_max_nbinomGLMTest$dtalbe[, 'pval'] #abs_dtable_pool_max_nbinomTest
names(abs_default_deseq_p) <- row.names(abs_dtable_pool_max_nbinomGLMTest$dtalbe)
# scde
# save.image('tmp_benchmark_analysis.RData')

scde_p <- std_scde_res_list$pval 
abs_scde_p <- abs_scde_res_list$pval #_no_normalize

# generate df3.1 for making the comparision: 
#select only high expressed genes: 
high_gene_list <- esApply(absolute_cds[, colnames(AT12_cds_subset_all_gene)], 1, function(x) sum(x > 1) > 15)
#overlap genes: 
permutation_sets <- list(FPKM = names(std_permutate_pval[which(std_permutate_pval < .01)]), 
                         transcript = names(mode_size_norm_permutate_ratio_by_geometric_mean[mode_size_norm_permutate_ratio_by_geometric_mean < 0.01]),
                         estimate_transcript = names(mc_mode_size_norm_permutate_ratio_by_geometric_mean[mc_mode_size_norm_permutate_ratio_by_geometric_mean < 0.01]),
                         readcount = names(readcount_permutate_pval[readcount_permutate_pval < 0.01]))
overlap_genes <- Reduce(intersect,  permutation_sets) 

#prepare the data for the DESeq2 / edgeR: 
default_edgeR_p_glm = edgeR_res_glm$lrt$table$PValue
names(default_edgeR_p_glm) <- row.names(edgeR_res_glm$lrt$table)
abs_default_edgeR_p_glm = abs_edgeR_res_glm$lrt$table$PValue
names(abs_default_edgeR_p_glm) <- row.names(abs_edgeR_res_glm$lrt$table)

#non_glm: 
default_edgeR_p = edgeR_res$et$table$PValue
names(default_edgeR_p) <- row.names(edgeR_res$et$table)
abs_default_edgeR_p = abs_edgeR_res$et$table$PValue
names(abs_default_edgeR_p) <- row.names(abs_edgeR_res$et$table)


default_deseq2_p = results(deseq2_res)$pvalue
names(default_deseq2_p) <- row.names(results(deseq2_res))
abs_default_deseq2_p = results(abs_deseq2_res)$pvalue 
names(abs_default_deseq2_p) <- row.names(results(abs_deseq2_res))
 
#generate the dataframe for making the benchmarking plots: 
############NOTE: mode_size_norm_permutate_ratio_by_geometric_mean may be changed into the same as DESEQ size normalization############ 
df3 <- plot_pre_rec_f1(test_p_list = list(monocle_p = monocle_p, 
                                        monocle_p_readcount = monocle_p_readcount,
                                        mode_size_norm_permutate_ratio_by_geometric_mean = new_abs_size_norm_monocle_p_ratio_by_geometric_mean,
                                        mc_mode_size_norm_permutate_ratio_by_geometric_mean = new_mc_size_norm_monocle_p_ratio_by_geometric_mean, 
                                        default_edgeR_p = default_edgeR_p, 
                                        abs_default_edgeR_p = abs_default_edgeR_p,         
                                        default_deseq2_p = default_deseq2_p, 
                                        abs_default_deseq2_p = abs_default_deseq2_p, 
                                        default_deseq_p = default_deseq_p, 
                                        abs_default_deseq_p = abs_default_deseq_p, 
                                        scde_p = scde_p, 
                                        abs_scde_p = abs_scde_p),
                     permutate_pval = list(monocle_p = std_permutate_pval, #readcount_permutate_pval, #std_permutate_pval, 
                                           monocle_p_readcount = readcount_permutate_pval, 
                                           mode_size_norm_permutate_ratio_by_geometric_mean = mode_size_norm_permutate_ratio_by_geometric_mean,
                                           mc_mode_size_norm_permutate_ratio_by_geometric_mean = mc_mode_size_norm_permutate_ratio_by_geometric_mean,
                                           default_edgeR_p = readcount_permutate_pval, #readcount_permutate_pval use readcount instead of std_permutate_pval
                                           abs_default_edgeR_p = mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           default_deseq2_p = readcount_permutate_pval, #readcount_permutate_pval use readcount instead of std_permutate_pval
                                           abs_default_deseq2_p = mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           default_deseq_p = readcount_permutate_pval, #readcount_permutate_pval use readcount instead of std_permutate_pval
                                           abs_default_deseq_p = mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           # abs_default_deseq_p_new_norm = mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           scde_p = readcount_permutate_pval, 
                                           abs_scde_p = mode_size_norm_permutate_ratio_by_geometric_mean),
                     row.names(absolute_cds), #gene_list, overlap_genes, high_gene_list
                     return_df = T, #na.rm = T, 
                     p_thrsld = 0.01, #0.05
                     rownames = c('monocle (FPKM)', 'monocle (readcount)', 'monocle (New size normalization)', 'monocle (New size normalization, Estimate transcript)', 
                        'edgeR (edgeR size normalization)', 'edgeR (New Size normalization)', 'DESeq2 (DESeq2 size normalization)', 'DESeq2 (New Size normalization)',
                        'DESeq (DESeq size normalization)', "DESeq (New Size normalization)", 'SCDE (Read Counts)', 'SCDE (New size normalization)'))
df3$data_type = c("Spikein transcripts", "Read counts", "Spikein transcripts", "Read counts", 

"Spikein transcripts", "Read counts", "Spikein transcripts", "Read counts", 

"MC transcripts", "Spikein transcripts", "Read counts", "FPKM")

df3$class = '3relative'


# only show new size normalization: 
df3.1 <- df3
df3.1[, 'Type'] <- c('SCDE', 'SCDE', 'DESeq1', 'DESeq1', 'DESeq2', 'DESeq2', 'edgeR', 'edgeR', 'Monocle', 'Monocle', 'Monocle', 'Monocle') # geom_bar(stat = 'identity', position = 'dodge') 

pdf('./supplementary_figures/fig2a_si.pdf', width = 3, height = 2)
qplot(factor(Type), value, stat = "identity", geom = 'bar', position = 'dodge', fill = data_type, data = melt(df3.1)) + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + monocle_theme_opts() + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) + nm_theme()
dev.off()

pdf('./supplementary_figures/fig2a_si_helper.pdf', width = 3, height = 2)
qplot(factor(Type), value, stat = "identity", geom = 'bar', position = 'dodge', fill = data_type, data = melt(df3.1)) + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + monocle_theme_opts() + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) #+ nm_theme()
dev.off()

# distribution of pval vs the mean expression: 
pval_mean <- data.frame(std = std_permutate_pval, abs = mode_size_norm_permutate_ratio_by_geometric_mean, 
                        mc = mc_mode_size_norm_permutate_ratio_by_geometric_mean, read = readcount_permutate_pval, 
                        fpkm_mean = apply(new_std_cds_14_18[1:transcript_num, ], 1, mean), 
                        spike_mean = apply(new_abs_cds_14_18[1:transcript_num, ], 1, mean), 
                        mc_mean = apply(new_mc_cds_14_18[1:transcript_num, ], 1, mean), 
                        read_mean = apply(count_cds[1:transcript_num, ], 1, mean) )


# # fig 3h: 
# # show only the spike-in / mc algorithm test: 
mc_spikein_df <- plot_pre_rec_f1(test_p_list = list(mode_size_norm_permutate_ratio_by_geometric_mean = new_abs_size_norm_monocle_p_ratio_by_geometric_mean,
                                    mc_mode_size_norm_permutate_ratio_by_geometric_mean = new_mc_size_norm_monocle_p_ratio_by_geometric_mean),
                 permutate_pval = list(mode_size_norm_permutate_ratio_by_geometric_mean = mode_size_norm_permutate_ratio_by_geometric_mean,
                                       mc_mode_size_norm_permutate_ratio_by_geometric_mean = mc_mode_size_norm_permutate_ratio_by_geometric_mean),
                 row.names(absolute_cds), #gene_list, overlap_genes, high_gene_list
                 return_df = T, #na.rm = T, 
                 p_thrsld = 0.01, #0.05
                 rownames = c('monocle (New size normalization)', 'monocle (New size normalization, Estimate transcript)'))
mc_spikein_df$data_type = c("Spikein transcripts", "estimated transcripts")

mc_spikein_df[, 'Type'] <- c('Monocle', 'Monocle') # geom_bar(stat = 'identity', position = 'dodge') 
colnames(mc_spikein_df)[1:3] <- c('Precision', 'Recall', 'F1 score')

pdf('./main_figures/fig3h.pdf', width = 1.7, height = 1.9)
ggplot(aes(factor(Type), value,  fill = data_type), data = melt(mc_spikein_df)) + geom_bar(position = position_dodge(), stat = 'identity') + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) + nm_theme() + xlab('') + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

pdf('./tmp/fig3h_helper.pdf', width = 3, height = 2)
ggplot(aes(factor(Type), value,  fill = data_type), data = melt(mc_spikein_df)) + geom_bar(position = position_dodge(), stat = 'identity') + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank())
dev.off()


# save.image('./RData/deg_benchmark_analysis.RData')
