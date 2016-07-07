library(ROCR)
library(monocle)
# library(devtools)
# load_all('~/Projects/monocle-dev')
library(plyr)
library(xacHelper)
library(grid)


# source(paste(script.dir, '/roc_curves.R', sep = ''))
#load the data: 
# load(paste('./RData/cmpr_three_packages', conditions[1], conditions[2], '.RData', sep = ''))

lung_pval_df <- data.frame(monocle_p = monocle_p, 
                                        monocle_p_readcount = monocle_p_readcount,
                                        mode_size_norm_permutate_ratio_by_geometric_mean = new_abs_size_norm_monocle_p_ratio_by_geometric_mean,
                                        mc_mode_size_norm_permutate_ratio_by_geometric_mean = new_mc_size_norm_monocle_p_ratio_by_geometric_mean, 
                                        default_edgeR_p = default_edgeR_p_glm, 
                                        abs_default_edgeR_p = abs_default_edgeR_p_glm,         
                                        default_deseq2_p = default_deseq2_p, 
                                        abs_default_deseq2_p = abs_default_deseq2_p, 
                                        default_deseq_p = default_deseq_p, 
                                        abs_default_deseq_p = abs_default_deseq_p, 
                                        scde_p = scde_p[names(monocle_p)], 
                                        abs_scde_p = abs_scde_p[names(monocle_p)] ,

	                                    mast_abs_pval_no_norm = mast_abs_pval_norm[names(monocle_p)], #no_norm
	                                    mast_mc_pval_no_norm = mast_mc_pval_norm[names(monocle_p)], 
	                                    mast_std_pval_no_norm = mast_std_pval_norm[names(monocle_p)], 
	                                    mast_count_pval_no_norm = mast_count_pval_norm[names(monocle_p)], 

                                        # tpm_count_monocle_p = tpm_count_monocle_p, 
                                        # tpm_count_edgeR_p_glm = tpm_count_edgeR_p_glm, 
                                        # tpm_count_deseq2_p = tpm_count_deseq2_p, 
                                        # tpm_count_deseq_p = tpm_count_deseq_p, 
                                        # tpm_count_scde_p = tpm_count_scde_p, 

                                        # mc_count_monocle_p = new_mc_size_norm_monocle_p_ratio_by_geometric_mean, 
                                        mc_count_edgeR_p_glm = mc_edgeR_p_glm, 
                                        mc_count_deseq2_p = mc_default_deseq2_p, 
                                        mc_count_deseq_p = mc_default_deseq_p, 
                                        mc_count_scde_p = mc_scde_p[names(monocle_p)] 
                                        )

row.names(lung_pval_df) <- names(monocle_p)
permutation_lung_pval_df <- data.frame(monocle_p = std_permutate_pval, #readcount_permutate_pval, #std_permutate_pval, 
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
                                           abs_scde_p = mode_size_norm_permutate_ratio_by_geometric_mean,

                                           mast_abs_pval_no_norm = mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           mast_mc_pval_no_norm = mc_mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           mast_std_pval_no_norm = std_permutate_pval, 
                                           mast_count_pval_no_norm = readcount_permutate_pval,

                                           # tpm_count_monocle_p = TPM_permutate_pval, #readcount_permutate_pval, #std_permutate_pval, 
                                           # tpm_count_edgeR_p_glm = TPM_permutate_pval, 
                                           # tpm_count_deseq2_p = TPM_permutate_pval,
                                           # tpm_count_deseq_p = TPM_permutate_pval,
                                           # tpm_count_scde_p = TPM_permutate_pval, 

                                           # mc_count_monocle_p = mc_mode_size_norm_permutate_ratio_by_geometric_mean, #readcount_permutate_pval, #std_permutate_pval, 
                                           mc_count_edgeR_p_glm = mc_mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           mc_count_deseq2_p = mc_mode_size_norm_permutate_ratio_by_geometric_mean,
                                           mc_count_deseq_p = mc_mode_size_norm_permutate_ratio_by_geometric_mean,
                                           mc_count_scde_p = mc_mode_size_norm_permutate_ratio_by_geometric_mean)
row.names(permutation_lung_pval_df) <- names(monocle_p)

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

####
plot_roc_df <-function(p_value, classification, type = 'fpr') {
  p_value[is.na(p_value)] <- 1
  pred_p_value <- prediction(p_value, classification)
  perf_tpr_fpr <- performance(pred_p_value, "tpr", "fpr")
  
  pdf('roc_plot_ori.pdf')
  plot(perf_tpr_fpr)
  dev.off()
}

perm_pvals <- permutation_lung_pval_df[select_genes, 1]
# select_genes <- names[perm_pvals[!is.na(perm_pvals)]]
software_pvals <- lung_pval_df[select_genes, 1] 

select_genes <- row.names(new_std_cds_14_18[1:transcript_num])[esApply(new_std_cds_14_18[1:transcript_num], 1, function(x) sum(x > 0.5) > 10)]
p_thrsld <- 0.1

perm_pvals[is.na(perm_pvals)] <- 1
software_pvals[is.na(software_pvals)] <- 1
plot_roc_df(software_pvals, perm_pvals > p_thrsld)
####

# select_genes <- row.names(standard_cds[1:transcript_num])
# > names(lung_pval_df)
#  [1] "monocle_p"
#  [2] "monocle_p_readcount"
#  [3] "mode_size_norm_permutate_ratio_by_geometric_mean"
#  [4] "mc_mode_size_norm_permutate_ratio_by_geometric_mean"
#  [5] "default_edgeR_p"
#  [6] "abs_default_edgeR_p"
#  [7] "default_deseq2_p"
#  [8] "abs_default_deseq2_p"
#  [9] "default_deseq_p"
# [10] "abs_default_deseq_p"
# [11] "scde_p"
# [12] "abs_scde_p"
# [13] "mast_abs_pval_no_norm"
# [14] "mast_mc_pval_no_norm"
# [15] "mast_count_pval_no_norm"
# [16] "mc_count_edgeR_p_glm"
# [17] "mc_count_deseq2_p"
# [18] "mc_count_deseq_p"
# [19] "mc_count_scde_p"
# >

lung_roc_df_list <- lapply(colnames(lung_pval_df), function(x) {
	print(x)
	
  # if(x %in% c('mode_size_norm_permutate_ratio_by_geometric_mean', 'abs_default_edgeR_p', 'abs_default_deseq2_p', 'abs_scde_p', 'mast_abs_pval_no_norm'))
  #   select_genes <- row.names(new_abs_cds_14_18[1:transcript_num])[esApply(new_abs_cds_14_18[1:transcript_num], 1, function(x) sum(x > 0.5) > 10)]
  # if(x %in% c('monocle_p', 'mast_std_pval'))
  #   select_genes <- row.names(new_std_cds_14_18[1:transcript_num])[esApply(new_std_cds_14_18[1:transcript_num], 1, function(x) sum(x > 0.5) > 10)]
  # if(x %in% c('monocle_p_readcount', 'default_edgeR_p', 'default_deseq2_p', 'default_deseq_p', 'scde_p', 'mast_count_pval_no_norm'))
  #   select_genes <- row.names(count_cds[1:transcript_num])[esApply(count_cds[1:transcript_num], 1, function(x) sum(x > 0.5) > 10)]
  # if(x %in% c('mc_mode_size_norm_permutate_ratio_by_geometric_mean', 'mast_mc_pval_no_norm', 'mc_count_edgeR_p_glm', 'mc_count_deseq2_p', 'mc_count_deseq_p', 'mc_count_scde_p'))
  #   select_genes <- row.names(mc_adj_cds[1:transcript_num])[esApply(mc_adj_cds[1:transcript_num], 1, function(x) sum(x > 0.5) > 10)]

  perm_pvals <- permutation_lung_pval_df[select_genes, 3] #use the spike-in as the gold standard
	software_pvals <- lung_pval_df[select_genes, x] 

	perm_pvals[is.na(perm_pvals)] <- 1
	software_pvals[is.na(software_pvals)] <- 1
	res <- generate_roc_df(software_pvals, perm_pvals > p_thrsld)
	colnames(res) <- c('tpr', 'fpr', 'auc')
	cbind(res, method = x)
})

lung_roc_df_list <- lapply(lung_roc_df_list, function(x) {colnames(x) <- c('tpr', 'fpr', 'auc', 'method'); x} )
lung_roc_df <- do.call(rbind, lung_roc_df_list)
lung_roc_df[1:5, ]
# str_split_fixed(row.names(lung_roc_df), '\\.', 2)[, 1]



lung_auc <- unique(lung_roc_df[, c('method', 'auc')])
row.names(lung_auc) <- lung_auc$method

save(file = './lung_roc_df', lung_roc_df)

lung_roc_df$software <- revalue(lung_roc_df$method, c("monocle_p" = 'Monocle', "monocle_p_readcount" = 'Monocle', "mode_size_norm_permutate_ratio_by_geometric_mean" = "Monocle", "mc_mode_size_norm_permutate_ratio_by_geometric_mean" = 'Monocle', 
                                            "default_edgeR_p" = 'edgeR', "mc_count_edgeR_p_glm" = 'edgeR', "abs_default_edgeR_p" = "edgeR", 
                                            "default_deseq2_p" = 'DESeq2', "mc_count_deseq2_p" = 'DESeq2', "abs_default_deseq2_p" = "DESeq2", 
                                            "default_deseq_p" = 'DESeq', "mc_count_deseq_p" = 'DESeq', "abs_default_deseq_p" = "DESeq",
                                            "scde_p" = 'SCDE', "mc_count_scde_p" = 'SCDE', "abs_scde_p" = "SCDE",
                                            "mast_std_pval_no_norm" = "MAST", 
                                            "mast_abs_pval_no_norm" = 'MAST', 
                                            "mast_mc_pval_no_norm" = "MAST", 
                                            "mast_count_pval_no_norm" = "MAST"))

lung_roc_df$Type <- revalue(lung_roc_df$method, c("monocle_p" = 'FPKM', "monocle_p_readcount" = 'Read counts', "mode_size_norm_permutate_ratio_by_geometric_mean" = "Transcript counts", "mc_mode_size_norm_permutate_ratio_by_geometric_mean" = 'Estimated transcript counts', 
                                            "default_edgeR_p" = 'Read counts', "mc_count_edgeR_p_glm" = 'Estimated transcript counts', "abs_default_edgeR_p" = "Transcript counts", 
                                            "default_deseq2_p" = 'Read counts', "mc_count_deseq2_p" = 'Estimated transcript counts', "abs_default_deseq2_p" = "Transcript counts", 
                                            "default_deseq_p" = 'Read counts', "mc_count_deseq_p" = 'Estimated transcript counts', "abs_default_deseq_p" = "Transcript counts",
                                            "scde_p" = 'Read counts', "mc_count_scde_p" = 'Estimated transcript counts', "abs_scde_p" = "Transcript counts",
                                            "mast_abs_pval_no_norm" = 'Estimated transcript counts', 
                                            "mast_std_pval_no_norm" = "FPKM", 
                                            "mast_mc_pval_no_norm" = "Estimated transcript counts", 
                                            "mast_count_pval_no_norm" = "Read counts"))

cols <- c("FPKM" = "#F2756D", "Read counts" = "#6F94CC", "Transcript counts" = "#000202", "Estimated transcript counts" = "#7BAE41") ##A680B9
pdf(paste('./supplementary_figures/lung_roc', condotions[1], condotions[2], '.pdf', sep = ''), height = 2, width = 3)
qplot(fpr, tpr, data= subset(lung_roc_df, software %in% c('Monocle', 'edgeR', 'DESeq2', 'SCDE')), geom="step", size = 0.5, color = Type) + #linetype = Type, 
    xlab("False positive rate") +
    ylab("True positive rate") +
      ylim(c(0, 1.0)) + facet_wrap(~software) + scale_size(range = c(0.1, 0.5)) + 
   scale_color_manual(values = cols, name = "Type") + 
  xlim(c(0, 1.0)) + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pdf('./supplementary_figures/lung_roc_helper.pdf', height = 13, width = 14)
qplot(fpr, tpr, data= lung_roc_df, geom="step", color = Type) + #linetype = Type, 
    xlab("False positive rate") +
    ylab("True positive rate") +
      ylim(c(0, 1.0)) + facet_wrap(~software) + 
   scale_color_manual(values = cols, name = "Type") #+ nm_theme()
  xlim(c(0, 1.0)) #+ nm_theme()
dev.off()

#ensure the color is the same: 
#"tpr"      "fpr"      "auc"      "method"   "software" "Type"
tmp <- data.frame(tpr = NA, fpr = NA, auc = NA, method = NA, 
                  software = c('DESeq', 'DESeq2', 'edgeR', 'MAST', 'SCDE'), 
                  Type = c('FPKM', 'FPKM', 'FPKM', 'FPKM', 'FPKM'))
lung_roc_df <- rbind(lung_roc_df, tmp) 

#roc_df <-  melt(df_res)
#                              software
# Type                          DESeq DESeq2 edgeR MAST Monocle SCDE
#   estimated transcript counts  1139   1239  1234 2484    1249 1207
#   FPKM                         1229      0     0 1242     996    0
#   read counts                     0   1289  1424    0    1313 1275

pdf('./supplementary_figures/lung_roc_dfroc_auc_bar.pdf', height = 3, width = 3)
ggplot(aes(software, auc), data = lung_roc_df) + geom_bar(position = 'dodge', stat = 'identity', aes(fill=Type)) + 
    xlab("") +
    # ylim(0.5, 1.0) + 
    # monocle_theme_opts() +  theme(axis.text.x=element_text(angle=30, hjust=1)) + 
    # scale_fill_manual(values = cols, name = "Software", label = test)  + 
    nm_theme()
dev.off()

pdf('./supplementary_figures/lung_roc_dfroc_auc_bar_helper.pdf', height = 6, width = 9)
ggplot(aes(software, auc), data = lung_roc_df) + geom_bar(position = 'dodge', stat = 'identity', aes(fill=Type)) + 
    xlab("") #+
    # ylim(0.5, 1.0) + 
    # monocle_theme_opts() +  theme(axis.text.x=element_text(angle=30, hjust=1)) + 
    # scale_fill_manual(values = cols, name = "Software", label = test)  + 
    # nm_theme()
dev.off()


#show the values of auc 
unique(lung_roc_df[, c('method','auc')])

# #save the result: 
save.image(paste('./RData/deg_benchmark_lung_roc_auc_', conditions[1], conditions[2], '.RData', sep = ''))

# subset(lung_roc_df, software == 'Monocle' & fpr < 0.1 & Type == 'FPKM')

# subset(lung_roc_df, software == 'Monocle' & fpr < 0.1)

# differentialGeneTest(new_std_cds_14_18[names(monocle_p)[std_permutate_pval > 0.5 & monocle_p < 0.001], ], 
#                                               fullModelFormulaStr = "~Time", 
#                                               reducedModelFormulaStr = "~1", cores = 1, relative = F, verbose = T) 

# esApply(new_std_cds_14_18[names(readcount_permutate_pval)[readcount_permutate_pval > 0.5 & default_deseq_p < 0.001], ], 2, function(x) sum(x > 0))

# names(readcount_permutate_pval)[readcount_permutate_pval > 0.5 & default_deseq_p < 0.001]

# new_std_diff_test_res[names(fits)]

# new_std_diff_test_res['ENSMUSG00000082196.2', ]

# dim(differentialGeneTest(new_std_cds_14_18[c('ENSMUSG00000082196.2', 'ENSMUSG00000082196.2'), ], 
#                                               fullModelFormulaStr = "~Time", 
#                                               reducedModelFormulaStr = "~1", cores = 2, relative = F))  

# differentialGeneTest(new_std_cds_14_18[c('ENSMUSG00000082196.2', 'ERCC-00077'), ], 
#                                               fullModelFormulaStr = "~Time", 
#                                               reducedModelFormulaStr = "~1", cores = 1, relative = F)

# qplot(readcount_permutate_pval, monocle_p_readcount) + ggsave('scatter_pval.pdf', heigth = 10, width = 10)

# which((readcount_permutate_pval > 0.5 & default_deseq_p < 0.001) == T)
# plot_genes_jitter(count_cds[which((readcount_permutate_pval > 0.5 & default_deseq_p < 0.001) == T)[1:10], ], color_by = 'Time', grouping = 'Time', nrow = 5, ncol = 2) + 
#   ggsave('example_genes.pdf')

# qplot(readcount_permutate_pval, default_deseq2_p) + ggsave('deseq2_scatter_pval.pdf') #, heigth = 10, width = 10

# qplot(std_permutate_pval, monocle_p) + ggsave('std_monocle_scatter_pval.pdf') #, heigth = 10, width = 10

# std_dtable_pool_max_nbinomGLMTest$dtalbe[which((readcount_permutate_pval > 0.5 & default_deseq_p < 0.001) == T)[1:10], ]




# default_deseq_p <- std_dtable_pool_max_nbinomGLMTest$dtalbe[, 'pval'] #std_dtable_pool_max_nbinomTest
# names(default_deseq_p) <- row.names(std_dtable_pool_max_nbinomGLMTest$dtalbe) #std_dtable_pool_max_nbinomGLMTest
# default_deseq_p[which(std_dtable_pool_max_nbinomGLMTest$dtalbe$converged)]









