# library(monocle)
library(devtools)
load_all('~/Projects/monocle-dev')
library(xacHelper)
library(MAST)
library(ROCR)

# source('roc_curves.R')
#load the data: 
load('./RData/analysis_HSMM_data.RData')
HSMM_bulk <- read.table("./data/HSMM_data/bulk_cuffdiff/gene_exp.diff", header = T, row.names = NULL)

HSMM_bulk_T0_T72 <-  subset(HSMM_bulk, sample_1 == "T0" & sample_2 == "T24")
order_stats_HSMM_bulk_T0_T72 <-  HSMM_bulk_T0_T72[order(abs(HSMM_bulk_T0_T72$test_stat), decreasing=T), ]

#add read counts: 
HSMM_readcounts_matrix <- read.delim("./data/HSMM_data/muscle/HSMM/HSMM_cuffnorm_out/genes.count_table")
row.names(HSMM_readcounts_matrix) <- HSMM_readcounts_matrix$tracking_id
HSMM_readcounts_matrix <- HSMM_readcounts_matrix[,-1]

pd <- new("AnnotatedDataFrame", data = pData(HSMM_myo))
fd <- new("AnnotatedDataFrame", data = fData(HSMM_myo))
HSMM_readcounts_matrix_cds <-  newCellDataSet(as.matrix(HSMM_fpkm_matrix[row.names(HSMM_myo), colnames(HSMM_myo)]), 
                                   phenoData = pd, 
                                   featureData = fd, 
                                   expressionFamily=negbinomial(), 
                                   lowerDetectionLimit=1)

##integrate into the benchmark analysis: 

##########################based on the deg_benchmark_analysis.R#############################

library(scde)
#########################run the following scripts in remote sever#############################

HSMM_valid_cells <- colnames(HSMM_myo)
#generate the pvals from the statistical test (permutation based or from software: monocle/DESeq/SCDE)

#double check the cells selected are the same as we used in the paper:  identical(colnames(new_mc_cds_T0_T72), cells_used)
new_mc_cds_T0_T72 <- HSMM_myo[, intersect(HSMM_valid_cells, row.names(subset(pData(HSMM_myo), Time %in% c('T72_CT_', 'T0_CT_'))))]
new_mc_cds_T0_T72 <- estimateSizeFactors(new_mc_cds_T0_T72) #calculate the size factor for performing the relative absolute expression tests
new_mc_cds_T0_T72 <- estimateDispersions(new_mc_cds_T0_T72) #calculate the size factor for performing the relative absolute expression tests

new_std_cds_T0_T72 <- std_HSMM[, colnames(new_mc_cds_T0_T72)]
# new_std_cds_T0_T72 <- estimateSizeFactors(new_std_cds_T0_T72) #FPKM values are already on relative scale

#create a cds for readcount data to perform the default DEG tests for DESeq and SCDE : 
count_cds_HSMM_bulk <- HSMM_readcounts_matrix_cds[row.names(new_mc_cds_T0_T72), colnames(new_mc_cds_T0_T72)] 

count_cds_HSMM_bulk <- estimateSizeFactors(count_cds_HSMM_bulk)
count_cds_HSMM_bulk <- estimateDispersions(count_cds_HSMM_bulk)

# calculate the pval with the readcount with scde: (calculate the scde associate DEG test result LOCALLY) 
std_scde_res_list <- scde_DEG(dir = NULL, count_cds = count_cds_HSMM_bulk, DEG_attribute = 'Time', contrast = c('T0_CT_', 'T72_CT_'), n.cores = 1)

#calculate the pval with the normalized transcripts with scde: 
mc_scde_res_list <- scde_DEG(dir = NULL, count_cds = new_mc_cds_T0_T72, DEG_attribute = 'Time', contrast = c('T0_CT_', 'T72_CT_'), n.cores = 1, normalize = T)
mc_scde_res_list_no_normalize <- scde_DEG(dir = NULL, count_cds = new_mc_cds_T0_T72, DEG_attribute = 'Time', contrast = c('T0_CT_', 'T72_CT_'), n.cores = 1)

#load all other necessary packages: 
load_all_libraries()

#perform  the stastical tests on the data: 
new_std_diff_test_res <- differentialGeneTest(new_std_cds_T0_T72[, ], 
                                              fullModelFormulaStr = "~Time", 
                                              reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)  

new_size_norm_mc_diff_test_res <- differentialGeneTest(new_mc_cds_T0_T72[, ], 
                                                       fullModelFormulaStr = "~Time", 
                                                       reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

new_cds_diff_test_res <- differentialGeneTest(count_cds_HSMM_bulk[, ], 
                                              fullModelFormulaStr = "~Time", 
                                              reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

new_cds_diff_test_res_no_relative <- differentialGeneTest(count_cds_HSMM_bulk[, ], 
                                                          fullModelFormulaStr = "~Time", 
                                                          reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)

Time_ori <- pData(new_mc_cds_T0_T72)$Time #Time as input for newCountDataSet
Time_ori <- as.character(Time_ori) #Time as input for newCountDataSet


#benchmark with edgeR / DESeq2:
edgeR_res <- edgeR_test(counts = exprs(count_cds_HSMM_bulk), glm = T)
mc_edgeR_res <- edgeR_test(exprs(new_mc_cds_T0_T72), group = Time_ori, glm = T)

edgeR_res_glm <- edgeR_test(counts = exprs(count_cds_HSMM_bulk))
mc_edgeR_res_glm <- edgeR_test(exprs(new_mc_cds_T0_T72), group = Time_ori)

#calculate the pval with the readcount with DESeq: 

read_count_d <- newCountDataSet(round(exprs(count_cds_HSMM_bulk)), Time_ori) 
std_dtable_pool_max_nbinomTest <- DESeq1_test(read_count_d, disp_method = 'pooled', contrast = c("T0_CT_", "T72_CT_"), sharing = 'maximum', test_type = 'nbinomTest', scale = T) 
row.names(std_dtable_pool_max_nbinomTest$dtalbe) <- std_dtable_pool_max_nbinomTest$dtalbe$id

# # calculate the pval with the normalized transcripts with DESeq: 
mc_count_d <- newCountDataSet(round(t(t(exprs(new_mc_cds_T0_T72)) / sizeFactors(new_mc_cds_T0_T72))), (Time_ori)) #normalized the data by size factor
mc_dtable_pool_max_nbinomTest <- DESeq1_test(mc_count_d, disp_method = 'pooled', contrast = c("T0_CT_", "T72_CT_"), sharing = 'maximum', test_type = 'nbinomTest') 
row.names(mc_dtable_pool_max_nbinomTest$dtalbe) <- mc_dtable_pool_max_nbinomTest$dtalbe$id

# DESeq glm: (GLM tests are more relevant to our software)
std_dtable_pool_max_nbinomGLMTest <- DESeq1_test(read_count_d[, ], disp_method = 'pooled', contrast = c("T0_CT_", "T72_CT_"), sharing = 'maximum', test_type = 'nbinomGLMTest', scale = T) 
row.names(std_dtable_pool_max_nbinomGLMTest$dtalbe) <-  row.names(read_count_d[, ])
mc_dtable_pool_max_nbinomGLMTest <- DESeq1_test(mc_count_d[, ], disp_method = 'pooled', contrast = c("T0_CT_", "T72_CT_"), sharing = 'maximum', test_type = 'nbinomGLMTest') 
row.names(mc_dtable_pool_max_nbinomGLMTest$dtalbe) <- row.names(new_mc_cds_T0_T72[, ])

#add the DEG tests using edgeR / DESeq2: 
edgeR_res <- edgeR_test(counts = exprs(count_cds_HSMM_bulk), glm = T)
mc_edgeR_res <- edgeR_test(exprs(new_mc_cds_T0_T72[, ]), group = Time_ori, glm = F)

# edgeR_res_glm <- edgeR_test()
mc_edgeR_res_glm <- edgeR_test(exprs(new_mc_cds_T0_T72[, ]), group = Time_ori, glm = T)

deseq2_res <- DESeq2_deg(dir = NULL, count_cds_HSMM_bulk[, ], Time = Time_ori, pd = pData(count_cds_HSMM_bulk))
mc_deseq2_res <- DESeq2_deg(dir = NULL, new_mc_cds_T0_T72[, ], Time = Time_ori, pd = pData(new_mc_cds_T0_T72))

#prepare to generate the data for create the precision/recall/F1 score data.frame: 
# #new monocle_p: 
monocle_p_HSMM_bulk <- new_std_diff_test_res[, 'pval'] 
names(monocle_p_HSMM_bulk) <- row.names(new_std_diff_test_res)
#use cds: 
monocle_p_readcount_HSMM_bulk <- new_cds_diff_test_res[, 'pval'] 
names(monocle_p_readcount_HSMM_bulk) <- row.names(new_cds_diff_test_res)

#mc
new_mc_size_norm_monocle_p_ratio_by_geometric_mean_HSMM_bulk <- new_size_norm_mc_diff_test_res[, 'pval']
names(new_mc_size_norm_monocle_p_ratio_by_geometric_mean_HSMM_bulk) <- row.names(new_size_norm_mc_diff_test_res)
#mc 
new_mc_size_norm_monocle_p_ratio_by_geometric_mean_HSMM_bulk <- new_size_norm_mc_diff_test_res[, 'pval']
names(new_mc_size_norm_monocle_p_ratio_by_geometric_mean_HSMM_bulk) <- row.names(new_size_norm_mc_diff_test_res)
#deseq
default_deseq_p_HSMM_bulk <- std_dtable_pool_max_nbinomGLMTest$dtalbe[, 'pval'] #std_dtable_pool_max_nbinomTest
names(default_deseq_p_HSMM_bulk) <- row.names(std_dtable_pool_max_nbinomGLMTest$dtalbe) #std_dtable_pool_max_nbinomGLMTest
mc_default_deseq_p_HSMM_bulk <- mc_dtable_pool_max_nbinomGLMTest$dtalbe[, 'pval'] #abs_dtable_pool_max_nbinomTest
names(mc_default_deseq_p_HSMM_bulk) <- row.names(mc_dtable_pool_max_nbinomGLMTest$dtalbe)
# scde
#save.image('tmp_benchmark_analysis.RData')

scde_p_HSMM_bulk <- std_scde_res_list$pval 
mc_scde_p_HSMM_bulk <- mc_scde_res_list$pval #_no_normalize

# generate df3.1 for making the comparision: 
#select only high expressed genes: 
# high_gene_list <- esApply(HSMM_myo[, colnames(AT12_cds_subset_all_gene)], 1, function(x) sum(x > 1) > 15)

#prepare the data for the DESeq2 / edgeR: 
default_edgeR_p_HSMM_bulk = edgeR_res$lrt$table$PValue
names(default_edgeR_p_HSMM_bulk) <- row.names(edgeR_res$lrt$table)
mc_default_edgeR_p = mc_edgeR_res$et$table$PValue
names(mc_default_edgeR_p) <- row.names(mc_edgeR_res$et$table)

default_deseq2_p_HSMM_bulk = results(deseq2_res)$pvalue
names(default_deseq2_p_HSMM_bulk) <- row.names(results(deseq2_res))
mc_default_deseq2_p_HSMM_bulk = results(mc_deseq2_res)$pvalue 
names(mc_default_deseq2_p_HSMM_bulk) <- row.names(results(mc_deseq2_res))

#mast function: 
mast_mc_pval_no_norm_HSMM_bulk <- MAST_deg(new_mc_cds_T0_T72)
mast_std_pval_no_norm_HSMM_bulk <- MAST_deg(new_std_cds_T0_T72)
mast_count_pval_no_norm_HSMM_bulk <- MAST_deg(count_cds_HSMM_bulk)

mast_mc_pval_norm <- MAST_deg(new_mc_cds_T0_T72, normalization = T)
# mast_std_pval_norm <- MAST_deg(new_std_cds_T0_T72, normalization = T)
mast_count_pval_norm <- MAST_deg(count_cds_HSMM_bulk, normalization = T)

#generate the dataframe for making the benchmarking plots: 
############NOTE: mode_size_norm_permutate_ratio_by_geometric_mean may be changed into the same as DESEQ size normalization############ 
df3_HSMM_bulk <- plot_pre_rec_f1(test_p_list = list(monocle_p = monocle_p_HSMM_bulk, 
                                        monocle_p_readcount = monocle_p_readcount_HSMM_bulk,
                                        mc_mode_size_norm_permutate_ratio_by_geometric_mean = new_mc_size_norm_monocle_p_ratio_by_geometric_mean_HSMM_bulk, 
                                        default_edgeR_p = default_edgeR_p_HSMM_bulk, 
                                        mc_default_edgeR_p = mc_default_edgeR_p,         
                                        default_deseq2_p = default_deseq2_p_HSMM_bulk, 
                                        mc_default_deseq2_p = mc_default_deseq2_p_HSMM_bulk, 
                                        default_deseq_p = default_deseq_p_HSMM_bulk, 
                                        mc_default_deseq_p = mc_default_deseq_p_HSMM_bulk, 
                                        scde_p = scde_p_HSMM_bulk, 
                                        mc_scde_p = mc_scde_p_HSMM_bulk, 
                                        # mast_abs_pval_no_norm = mast_abs_pval_no_norm_HSMM_bulk, 
                                        mast_mc_pval_no_norm = mast_mc_pval_no_norm_HSMM_bulk, 
                                        mast_std_pval_no_norm = mast_std_pval_no_norm_HSMM_bulk, 
                                        mast_count_pval_no_norm = mast_count_pval_no_norm_HSMM_bulk),
                     permutate_pval = list(monocle_p = top_1k_HSMM_bulk_T0_T72_pval, #readcount_permutate_pval, #std_permutate_pval, 
                                           monocle_p_readcount = top_1k_HSMM_bulk_T0_T72_pval, 
                                           mc_mode_size_norm_permutate_ratio_by_geometric_mean = top_1k_HSMM_bulk_T0_T72_pval,
                                           default_edgeR_p = top_1k_HSMM_bulk_T0_T72_pval, #readcount_permutate_pval use readcount instead of std_permutate_pval
                                           mc_default_edgeR_p = top_1k_HSMM_bulk_T0_T72_pval, 
                                           default_deseq2_p = top_1k_HSMM_bulk_T0_T72_pval, #readcount_permutate_pval use readcount instead of std_permutate_pval
                                           mc_default_deseq2_p = top_1k_HSMM_bulk_T0_T72_pval, 
                                           default_deseq_p = top_1k_HSMM_bulk_T0_T72_pval, #readcount_permutate_pval use readcount instead of std_permutate_pval
                                           mc_default_deseq_p = top_1k_HSMM_bulk_T0_T72_pval, 
                                           # abs_default_deseq_p_new_norm = mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           scde_p = top_1k_HSMM_bulk_T0_T72_pval, 
                                           mc_scde_p = top_1k_HSMM_bulk_T0_T72_pval, 
                                           # mast_abs_pval_no_norm = top_1k_HSMM_bulk_T0_T72_pval, 
                                           mast_mc_pval_no_norm = top_1k_HSMM_bulk_T0_T72_pval, 
                                           mast_std_pval_no_norm = top_1k_HSMM_bulk_T0_T72_pval, 
                                           mast_count_pval_no_norm = top_1k_HSMM_bulk_T0_T72_pval),
                     names(top_1k_HSMM_bulk_T0_T72_pval), #gene_list, overlap_genes, high_gene_list
                     return_df = T, #na.rm = T, 
                     p_thrsld = 0.01, #0.05
                     rownames = c('monocle (FPKM)', 'monocle (readcount)', 'monocle (New size normalization, Estimate transcript)', 
                        'edgeR (edgeR size normalization)', 'edgeR (New Size normalization)', 'DESeq2 (DESeq2 size normalization)', 'DESeq2 (New Size normalization)',
                        'DESeq (DESeq size normalization)', "DESeq (New Size normalization)", 'SCDE (Read Counts)', 'SCDE (New size normalization)', 
                        'MAST (mc no normalization)', 'MAST (FPKM no normalization)', 'MAST (readcount no normalization)'))
df3_HSMM_bulk$data_type = c("Read counts", "FPKM", "MC transcripts", "MC transcripts", "Read counts", "MC transcripts", "Read counts", 

"MC transcripts", "Read counts", "MC transcripts", "Read counts", 

"MC transcripts", "Read counts", "FPKM")

df3_HSMM_bulk$class = '3relative'


# only show new size normalization: 
df3.1_HSMM_bulk <- df3_HSMM_bulk
df3.1_HSMM_bulk[, 'Type'] <- c('MAST', 'MAST', 'MAST', 'SCDE', 'SCDE', 'DESeq1', 'DESeq1', 'DESeq2', 'DESeq2', 'edgeR', 'edgeR', 'Monocle', 'Monocle', 'Monocle') # geom_bar(stat = 'identity', position = 'dodge') 

pdf('./supplementary_figures/fig2a_si_HSMM_bulk.pdf', width = 3, height = 2)
qplot(factor(Type), value, stat = "identity", geom = 'bar', position = 'dodge', fill = data_type, data = melt(df3.1_HSMM_bulk)) + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + monocle_theme_opts() + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) + nm_theme()
dev.off()

pdf('./supplementary_figures/fig2a_si_HSMM_bulk_helper.pdf', width = 6, height = 2)
qplot(factor(Type), value, stat = "identity", geom = 'bar', position = 'dodge', fill = data_type, data = melt(df3.1_HSMM_bulk)) + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + monocle_theme_opts() + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) 
dev.off()

HSMM_bulk_pval_df <- data.frame(monocle_p = monocle_p_HSMM_bulk, 
                                        monocle_p_readcount = monocle_p_readcount_HSMM_bulk,
                                        mc_mode_size_norm_permutate_ratio_by_geometric_mean = new_mc_size_norm_monocle_p_ratio_by_geometric_mean_HSMM_bulk, 
                                        default_edgeR_p = default_edgeR_p_HSMM_bulk[names(monocle_p_HSMM_bulk)], 
                                        mc_default_edgeR_p = mc_default_edgeR_p[names(monocle_p_HSMM_bulk)],         
                                        default_deseq2_p = default_deseq2_p_HSMM_bulk[names(monocle_p_HSMM_bulk)], 
                                        mc_default_deseq2_p = mc_default_deseq2_p_HSMM_bulk[names(monocle_p_HSMM_bulk)], 
                                        default_deseq_p = default_deseq_p_HSMM_bulk[names(monocle_p_HSMM_bulk)], 
                                        mc_default_deseq_p = mc_default_deseq_p_HSMM_bulk[names(monocle_p_HSMM_bulk)], 
                                        scde_p = scde_p_HSMM_bulk[names(monocle_p_HSMM_bulk)], 
                                        mc_scde_p = mc_scde_p_HSMM_bulk[names(monocle_p_HSMM_bulk)], 
                                        # mast_abs_pval_no_norm = mast_abs_pval_no_norm_HSMM_bulk, 
                                        mast_mc_pval_no_norm = mast_mc_pval_no_norm_HSMM_bulk, 
                                        mast_std_pval_no_norm = mast_std_pval_no_norm_HSMM_bulk, 
                                        mast_count_pval_no_norm = mast_count_pval_no_norm_HSMM_bulk)

top_1k_genes <- order_stats_HSMM_bulk_T0_T72[order_stats_HSMM_bulk_T0_T72$significant == 'yes', ][1:1000, 'gene_id']

reverse_order_stats_HSMM_bulk_T0_T72 <-  HSMM_bulk_T0_T72[order(abs(HSMM_bulk_T0_T72$q_value), decreasing=T), ]
bottom_1k_genes <- reverse_order_stats_HSMM_bulk_T0_T72[reverse_order_stats_HSMM_bulk_T0_T72$significant == 'no', ][1:1000, 'gene_id']

#top/bottom 1k genes: 
select_genes <- c(as.character(top_1k_genes), as.character(bottom_1k_genes))
true_data <- rep(0, length(select_genes))
true_data[1:length(top_1k_genes)] <- 1
true_data[(length(top_1k_genes) + 1):(length(select_genes))] <- 0

true_df <- data.frame(Type = true_data, pval = top_1k_HSMM_bulk_T0_T72_pval[select_genes])

generate_roc_df <-function(p_value, classification, type = 'fpr') {
	p_value[is.na(p_value)] <- 1
	pred_p_value <- prediction(p_value, classification)
	perf_tpr_fpr <- performance(pred_p_value, "tpr", "fpr")
	
    fpr = perf_tpr_fpr@x.values

    tpr = perf_tpr_fpr@y.values
    
	perf_auc <- performance(pred_p_value, "auc")
	auc <- perf_auc@y.values

    data.frame(tpr = tpr, fpr = fpr, auc = auc)
}

roc_df_list <- apply(HSMM_bulk_pval_df[select_genes, ], 2, function(x) generate_roc_df(x, true_data != 1))

roc_df_list <- lapply(roc_df_list, function(x) {colnames(x) <- c('tpr', 'fpr', 'auc'); x} )
roc_df <- do.call(rbind, roc_df_list)
roc_df$method <-  str_split_fixed(row.names(roc_df), '\\.', 2)[, 1]

auc <- unique(roc_df[, c('method', 'auc')])
row.names(auc) <- auc$method
hmcols <- blue2green2red(nrow(auc))

test <- c(paste(sort(colnames(HSMM_bulk_pval_df)), ':', " AUC = ", auc[order(colnames(HSMM_bulk_pval_df)), 2]))

roc_df$software <- revalue(roc_df$method, c("monocle_p" = 'Monocle', "monocle_p_readcount" = 'Monocle', "mc_mode_size_norm_permutate_ratio_by_geometric_mean" = 'Monocle', 
                                            "default_edgeR_p" = 'edgeR', "mc_default_edgeR_p" = 'edgeR',
                                            "default_deseq2_p" = 'DESeq2', "mc_default_deseq2_p" = 'DESeq2',
                                            "default_deseq_p" = 'DESeq', "mc_default_deseq_p" = 'DESeq',
                                            "scde_p" = 'SCDE', "mc_scde_p" = 'SCDE',
                                            "mast_mc_pval_no_norm" = 'MAST', 
                                            "mast_std_pval_no_norm" = "MAST", 
                                            "mast_count_pval_no_norm" = "MAST"))

roc_df$Type <- revalue(roc_df$method, c("monocle_p" = 'FPKM', "monocle_p_readcount" = 'read counts', "mc_mode_size_norm_permutate_ratio_by_geometric_mean" = 'estimated transcript counts', 
                                            "default_edgeR_p" = 'read counts', "mc_default_edgeR_p" = 'estimated transcript counts',
                                            "default_deseq2_p" = 'read counts', "mc_default_deseq2_p" = 'estimated transcript counts',
                                            "default_deseq_p" = 'FPKM', "mc_default_deseq_p" = 'estimated transcript counts',
                                            "scde_p" = 'read counts', "mc_scde_p" = 'estimated transcript counts',
                                            "mast_mc_pval_no_norm" = 'estimated transcript counts', 
                                            "mast_std_pval_no_norm" = "FPKM", 
                                            "mast_count_pval_no_norm" = "estimated transcript counts"))

pdf('./supplementary_figures/roc.pdf', height = 3, width = 4)
qplot(fpr, tpr, data=roc_df, geom="line", linetype = Type) + 
    xlab("False positive rate") +
    ylab("True positive rate") +
      ylim(c(0, 1.0)) + facet_wrap(~software) + 
  xlim(c(0, 1.0)) + monocle_theme_opts() +
   scale_color_manual(values = hmcols, name = "Test", label = test)  + nm_theme()
dev.off()

pdf('./supplementary_figures/roc_helper.pdf', height = 3, width = 4)
qplot(fpr, tpr, data=roc_df, geom="line", linetype = Type) + 
    xlab("False positive rate") +
    ylab("True positive rate") +
      ylim(c(0, 1.0)) + facet_wrap(~software) + 
  xlim(c(0, 1.0)) + monocle_theme_opts() +
   scale_color_manual(values = hmcols, name = "Test", label = test) 
dev.off()


pdf('./supplementary_figures/roc_auc_bar.pdf', height = 3, width = 4)
qplot(method, auc, stat = 'identity', data=roc_df, fill=method, geom="bar") + 
    xlab("") +
    ylim(c(0, 1.0)) + 
    monocle_theme_opts() +  theme(axis.text.x=element_text(angle=30, hjust=1)) + 
    scale_fill_manual(values = hmcols, name = "Software", label = test)  + nm_theme()
dev.off()

pdf('./supplementary_figures/roc_helper.pdf', height = 3, width = 14)
qplot(fpr, tpr, data=roc_df, color=method, geom="line") + 
    xlab("False positive rate") +
    ylab("True positive rate") +
      ylim(c(0, 1.0)) + 
  xlim(c(0, 1.0)) + monocle_theme_opts() +
   scale_color_manual(values = hmcols, name = "Test", label = test)  # + nm_theme()
dev.off()

save.image('./RData/deg_benchmark_analysis_HSMM_bulk.RData')

