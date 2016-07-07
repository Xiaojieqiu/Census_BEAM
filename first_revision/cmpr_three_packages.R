#function to compare the performance of the three packages
# load(paste('./RData/deg_benchmark_analysis', conditions[1], conditions[2], '.RData', sep = ''))
# #1. MAST 
# library(devtools)
# # install_github('RGLab/MAST')
# # library(monocle)
# library(devtools)
# load_all('~/Projects/monocle-dev')
# # library(MAST)
# library(xacHelper)
# library(grid)

# MAST_deg <- function(cds, grp = 'Time', test.type = 'hurdle', pseudo_cnt = 1, normalization = F) {
#   data <- as.matrix(exprs(cds))
  
#   if(normalization)
#     data <- log2(t(t(data) / sizeFactors(cds)) + pseudo_cnt)
#   else
#     data <- log2(data + pseudo_cnt)
  
#   data <- melt(data)
#   colnames(data) <- c('Gene', 'Cell', 'exprs')
#   data[, grp] <- c(pData(cds)[data[, 2], grp])
#   data$ncells <- 1 
  
#   mast_cds <- FluidigmAssay(data, idvars="Cell", 
#                             primerid='Gene', measurement='exprs', geneid="Gene", 
#                             ncells = 'ncells', phenovars=grp)
  
#   zlm.output <- zlm.SingleCellAssay(as.formula(paste("~ ", grp)), mast_cds, method='glm', ebayes=TRUE) #
#   # show(zlm.output)
#   zlm.lr <- lrTest(zlm.output, grp)
#   # dimnames(zlm.lr)
  
#   pval_df <- zlm.lr[,,'Pr(>Chisq)']
#   pval <- pval_df[, test.type]
  
#   return(pval)
# }

#run the following separately because MAST requires a different version of ggplot2 we are using
# mast_abs_pval_no_norm <- MAST_deg(new_abs_cds_14_18)
# mast_mc_pval_no_norm <- MAST_deg(new_mc_cds_14_18)
# mast_std_pval_no_norm <- MAST_deg(new_std_cds_14_18)
# mast_count_pval_no_norm <- MAST_deg(count_cds)
# mast_TPM_pval_no_norm <- MAST_deg(lung_TPM_count_cds_14_18)

# mast_abs_pval_norm <- MAST_deg(new_abs_cds_14_18, normalization = T)
# mast_mc_pval_norm <- MAST_deg(new_mc_cds_14_18, pseudo_cnt = 0.1, normalization = T)
# mast_std_pval_norm <- MAST_deg(new_std_cds_14_18, pseudo_cnt = 0.1, normalization = T)
# mast_count_pval_norm <- MAST_deg(count_cds, pseudo_cnt = 0.1, normalization = T)
# mast_TPM_pval_norm <- MAST_deg(lung_TPM_count_cds_14_18, normalization = T)

# # save(mast_abs_pval_no_norm, mast_mc_pval_no_norm, mast_std_pval_no_norm, mast_count_pval_no_norm, mast_TPM_pval_no_norm, mast_abs_pval_norm, mast_mc_pval_norm, mast_std_pval_norm, mast_count_pval_norm, mast_TPM_pval_norm, file = './RData/cmpr_three_package_mast_res')


# MAST_readcount_split_cds <- split(log2(t(round(exprs(count_cds[1:transcript_num, Time_order]) / readcount_sf_mat[1:transcript_num, Time_order])) + 1), col(t(exprs(count_cds[1:transcript_num, Time_order])), as.factor = T))
# MAST_readcount_fc <- apply(log2(round(exprs(count_cds)[1:transcript_num, Time_order] / readcount_sf_mat[1:transcript_num, Time_order]) + 1), 1, mean_fc, grp0 = 'E14.5', grp1 = 'E18.5', grp = pData(count_cds[1:transcript_num, Time_order])$Time) #valid_gene_id
# MAST_readcount_split_fc <- split(t(MAST_readcount_fc), col(t(MAST_readcount_fc), as.factor = T))
# MAST_readcount_permutate_pval <- mcmapply(permuation_pval, MAST_readcount_split_cds, MAST_readcount_fc, mc.cores = detectCores()) #multiple core
# closeAllConnections()

# MAST_abs_split_cds <- split(log2(round(t(exprs(new_abs_cds_14_18[1:transcript_num, Time_order]) / abs_sf_mat[1:transcript_num, Time_order])) + 0.1), col(t(exprs(new_abs_cds_14_18[1:transcript_num, Time_order])), as.factor = T))
# MAST_abs_fc <- apply(log2(round(exprs(new_abs_cds_14_18)[1:transcript_num, Time_order] / abs_sf_mat[1:transcript_num, Time_order]) + 0.1), 1, mean_fc, grp0 = 'E14.5', grp1 = 'E18.5', grp = pData(new_abs_cds_14_18[1:transcript_num, Time_order])$Time) #valid_gene_id
# MAST_abs_split_fc <- split(t(MAST_abs_fc), col(t(MAST_abs_fc), as.factor = T))
# MAST_mode_size_norm_permutate_ratio_by_geometric_mean <- mcmapply(permuation_pval, MAST_abs_split_cds, MAST_abs_split_fc, mc.cores = detectCores()) #multiple core
# closeAllConnections()

# MAST_mc_split_cds <- split(log2(round(t(exprs(new_mc_cds_14_18[1:transcript_num, Time_order]) / mc_sf_mat[1:transcript_num, Time_order])) + 0.1), col(t(exprs(new_mc_cds_14_18[1:transcript_num, Time_order])), as.factor = T))
# MAST_mc_fc <- apply(log2(round(exprs(new_mc_cds_14_18)[1:transcript_num, Time_order] / mc_sf_mat[1:transcript_num, Time_order]) + 0.1), 1, mean_fc, grp0 = 'E14.5', grp1 = 'E18.5', grp = pData(new_mc_cds_14_18[1:transcript_num, Time_order])$Time) #valid_gene_id
# MAST_mc_split_fc <- split(t(MAST_mc_fc), col(t(MAST_mc_fc), as.factor = T))
# MAST_mc_mode_size_norm_permutate_ratio_by_geometric_mean <- mcmapply(permuation_pval, MAST_mc_split_cds, MAST_mc_split_fc, mc.cores = detectCores()) #multiple core
# closeAllConnections()

# MAST_std_split_cds <- split(log2(t(exprs(new_std_cds_14_18[1:transcript_num, Time_order])) + 1), col(t(exprs(new_std_cds_14_18[1:transcript_num, Time_order])), as.factor = T))
# MAST_std_fc <- apply(log2(exprs(new_std_cds_14_18[1:transcript_num, ]) + 1), 1, mean_fc, grp0 = 'E14.5', grp1 = 'E18.5', grp = pData(count_cds[1:transcript_num, Time_order])$Time)
# MAST_std_split_fc <- split(t(MAST_std_fc), col(t(MAST_std_fc), as.factor = T))
# MAST_std_permutate_pval <- mcmapply(permuation_pval, MAST_std_split_cds, MAST_std_split_fc, mc.cores = detectCores()) #multiple cores 

# save(file = 'MAST_permutation_res', MAST_readcount_permutate_pval, MAST_mode_size_norm_permutate_ratio_by_geometric_mean, 
#   MAST_mc_mode_size_norm_permutate_ratio_by_geometric_mean, MAST_std_permutate_pval)

load('./RData/cmpr_three_package_mast_res')

#select genes for benchmarking the performance: (this should matach with the roc_auc plot) 
select_genes <- row.names(new_std_cds_14_18[1:transcript_num])[esApply(new_std_cds_14_18[1:transcript_num], 1, function(x) sum(x > 0.5) > 10)]

#integrate into the benchmark analysis:  
df3 <- plot_pre_rec_f1(test_p_list = list(monocle_p = monocle_p, 
                                        monocle_p_readcount = monocle_p_readcount,
                                        mode_size_norm_permutate_ratio_by_geometric_mean = new_abs_size_norm_monocle_p_ratio_by_geometric_mean,
                                        mc_mode_size_norm_permutate_ratio_by_geometric_mean = new_mc_size_norm_monocle_p_ratio_by_geometric_mean, 
                                        default_edgeR_p = default_edgeR_p_glm, 
                                        abs_default_edgeR_p = abs_default_edgeR_p_glm,           
                                        default_deseq2_p = default_deseq2_p, 
                                        abs_default_deseq2_p = abs_default_deseq2_p, 
                                        default_deseq_p = default_deseq_p, 
                                        abs_default_deseq_p = abs_default_deseq_p, 
                                        scde_p = scde_p,  
                                        abs_scde_p = abs_scde_p, 
                                        # mast_abs_pval_norm = mast_abs_pval_norm, 
                                        # mast_mc_pval_norm = mast_mc_pval_norm,  
                                        # #mast_std_pval_no_norm = mast_std_pval_no_norm, 
                                        # mast_count_pval_norm = mast_count_pval_norm
                                        mast_abs_pval_no_norm = mast_abs_pval_norm, #no_norm
                                        mast_mc_pval_no_norm = mast_mc_pval_norm, 
                                        mast_std_pval_no_norm = mast_std_pval_no_norm, 
                                        mast_count_pval_no_norm = mast_count_pval_norm, 

                                        # mc_count_monocle_p = new_mc_size_norm_monocle_p_ratio_by_geometric_mean, 
                                        mc_count_edgeR_p_glm = mc_edgeR_p_glm, 
                                        mc_count_deseq2_p = mc_default_deseq2_p, 
                                        mc_count_deseq_p = mc_default_deseq_p, 
                                        mc_count_scde_p = mc_scde_p                                        
                                        ),
                     permutate_pval = list(monocle_p = mode_size_norm_permutate_ratio_by_geometric_mean,#std_permutate_pval, #readcount_permutate_pval, #std_permutate_pval, 
                                           monocle_p_readcount = mode_size_norm_permutate_ratio_by_geometric_mean,#readcount_permutate_pval, 
                                           mode_size_norm_permutate_ratio_by_geometric_mean = mode_size_norm_permutate_ratio_by_geometric_mean,
                                           mc_mode_size_norm_permutate_ratio_by_geometric_mean = mode_size_norm_permutate_ratio_by_geometric_mean,#mc_mode_size_norm_permutate_ratio_by_geometric_mean,
                                           default_edgeR_p = mode_size_norm_permutate_ratio_by_geometric_mean,#readcount_permutate_pval, #readcount_permutate_pval use readcount instead of std_permutate_pval
                                           abs_default_edgeR_p = mode_size_norm_permutate_ratio_by_geometric_mean,#mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           default_deseq2_p = mode_size_norm_permutate_ratio_by_geometric_mean,#readcount_permutate_pval, #readcount_permutate_pval use readcount instead of std_permutate_pval
                                           abs_default_deseq2_p = mode_size_norm_permutate_ratio_by_geometric_mean,#mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           default_deseq_p = mode_size_norm_permutate_ratio_by_geometric_mean,#readcount_permutate_pval, #readcount_permutate_pval use readcount instead of std_permutate_pval
                                           abs_default_deseq_p = mode_size_norm_permutate_ratio_by_geometric_mean,#mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           # abs_default_deseq_p_new_norm = mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           scde_p = mode_size_norm_permutate_ratio_by_geometric_mean,#readcount_permutate_pval, 
                                           abs_scde_p = mode_size_norm_permutate_ratio_by_geometric_mean,#mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           
                                           mast_abs_pval_no_norm = mode_size_norm_permutate_ratio_by_geometric_mean,#MAST_mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           mast_mc_pval_no_norm = mode_size_norm_permutate_ratio_by_geometric_mean,#MAST_mc_mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           mast_std_pval_no_norm = mode_size_norm_permutate_ratio_by_geometric_mean,#MAST_std_permutate_pval, 
                                           mast_count_pval_no_norm = mode_size_norm_permutate_ratio_by_geometric_mean,#MAST_readcount_permutate_pval,
                                           # mast_abs_pval_norm = mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           # mast_mc_pval_norm = mc_mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           # #mast_std_pval_no_norm = std_permutate_pval, 
                                           # mast_count_pval_norm = readcount_permutate_pval

                                           # mc_count_monocle_p = mc_mode_size_norm_permutate_ratio_by_geometric_mean, #readcount_permutate_pval, #std_permutate_pval, 
                                           mc_count_edgeR_p_glm = mode_size_norm_permutate_ratio_by_geometric_mean,#mc_mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           mc_count_deseq2_p = mode_size_norm_permutate_ratio_by_geometric_mean,#mc_mode_size_norm_permutate_ratio_by_geometric_mean,
                                           mc_count_deseq_p = mode_size_norm_permutate_ratio_by_geometric_mean,#mc_mode_size_norm_permutate_ratio_by_geometric_mean,
                                           mc_count_scde_p = mode_size_norm_permutate_ratio_by_geometric_mean#mc_mode_size_norm_permutate_ratio_by_geometric_mean
                                           ),
                     row.names(absolute_cds), #gene_list, overlap_genes, high_gene_list
                     return_df = T, #na.rm = T, 
                     p_thrsld = 0.05, #0.05
                     rownames = c('monocle (FPKM)', 'monocle (readcount)', 'monocle (New size normalization)', 'monocle (New size normalization, Estimate transcript)', 
                        'edgeR (edgeR size normalization)', 'edgeR (New Size normalization)', 'DESeq2 (DESeq2 size normalization)', 'DESeq2 (New Size normalization)',
                        'DESeq (DESeq size normalization)', "DESeq (New Size normalization)", 'SCDE (Read Counts)', 'SCDE (New size normalization)', 
                        'MAST (absolute no normalization)', 'MAST (mc no normalization)', 'MAST (FPKM)', 'MAST (readcount no normalization)', 
                        'edgeR (New size normalization, Estimate transcript)', 'DESeq2 (New size normalization, Estimate transcript)', 'DESeq (New size normalization, Estimate transcript)', 'SCDE (New size normalization, Estimate transcript)'

                        )) # 'MAST (FPKM no normalization)',
df3$data_type = c("MC transcripts", "MC transcripts", "MC transcripts", "MC transcripts", 
"Read counts", "FPKM", "MC transcripts", "Spikein transcripts", "Spikein transcripts", "Read counts", "Spikein transcripts", "Read counts", 

"Spikein transcripts", "Read counts", "Spikein transcripts", "Read counts", 

"MC transcripts", "Spikein transcripts", "Read counts", "FPKM")

df3$class = '3relative'


# only show new size normalization: 
df3.1 <- df3
df3.1[, 'Type'] <- c('SCDE', 'DESeq1', 'DESeq2', 'edgeR', 'MAST', 'MAST', 'MAST', 'MAST', 'SCDE', 'SCDE', 'DESeq1', 'DESeq1', 'DESeq2', 'DESeq2', 'edgeR', 'edgeR', 'Monocle', 'Monocle', 'Monocle', 'Monocle') # geom_bar(stat = 'identity', position = 'dodge') 

tmp <- data.frame(Type = c('SCDE', 'DESeq1', 'DESeq2', 'edgeR'), 
                  data_type = c('FPKM', 'FPKM', 'FPKM', 'FPKM'),
                  class = '3relative', 
                  pre = NA, rec = NA, f1 = NA)
df_res <- rbind(df3.1, tmp) 
df_res <-  melt(df_res)

pdf('./supplementary_figures/fig2a_si_test.pdf', width = 4.5, height = 2)
qplot(factor(Type), value, stat = "identity", geom = 'bar', position = 'dodge', fill = data_type, data = subset(df_res, Type %in% c("DESeq2", "edgeR", "Monocle", "SCDE"))) + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + monocle_theme_opts() + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) + nm_theme()
dev.off()

pdf('./tmp/fig2a_si_test_helper.pdf', width = 3, height = 2)
qplot(factor(Type), value, stat = "identity", geom = 'bar', position = 'dodge', fill = data_type, data = df_res) + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + monocle_theme_opts() + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank())
dev.off()

#save the permutation result as a supplementary file: 
permutate_df <- data.frame(FPKM = std_permutate_pval, read_count = readcount_permutate_pval, spike_free = mc_mode_size_norm_permutate_ratio_by_geometric_mean, spike = mode_size_norm_permutate_ratio_by_geometric_mean)
row.names(permutate_df) <- names(std_permutate_pval)
write.table(permutate_df, sep = "\t", quote = F, row.names = T, file = './supplementary_data/supplementary_file_1.txt')

# save.image('./RData/cmpr_three_packages.RData')
save.image(paste('./RData/cmpr_three_packages', conditions[1], conditions[2], '.RData', sep = ''))


#2. SUBRA (matlab code) this is also a very different tool and we are not going to test that too


#3. BASiCs (very different tools and we are not going to test this) 
# install_github('catavallejos/BASiCS')



