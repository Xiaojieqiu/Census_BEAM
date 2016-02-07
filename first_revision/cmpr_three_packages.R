#function to compare the performance of the three packages

#1. MAST 
# library(devtools)
# install_github('RGLab/MAST')

library(MAST)

mast_abs_pval_no_norm <- MAST_deg(new_abs_cds_14_18)
mast_mc_pval_no_norm <- MAST_deg(new_mc_cds_14_18)
mast_std_pval_no_norm <- MAST_deg(new_std_cds_14_18)
mast_count_pval_no_norm <- MAST_deg(count_cds)

mast_abs_pval_norm <- MAST_deg(new_abs_cds_14_18, normalization = T)
mast_mc_pval_norm <- MAST_deg(new_mc_cds_14_18, normalization = T)
mast_std_pval_norm <- MAST_deg(new_std_cds_14_18, normalization = T)
mast_count_pval_norm <- MAST_deg(count_cds, normalization = T)

#integrate into the benchmark analysis: 
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
                                        abs_scde_p = abs_scde_p, 
                                        mast_abs_pval_no_norm = mast_abs_pval_no_norm, 
                                        mast_mc_pval_no_norm = mast_mc_pval_no_norm, 
                                        mast_std_pval_no_norm = mast_std_pval_no_norm, 
                                        mast_count_pval_no_norm = mast_count_pval_no_norm),
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
                                           abs_scde_p = mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           mast_abs_pval_no_norm = mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           mast_mc_pval_no_norm = mc_mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           mast_std_pval_no_norm = std_permutate_pval, 
                                           mast_count_pval_no_norm = readcount_permutate_pval),
                     row.names(absolute_cds), #gene_list, overlap_genes, high_gene_list
                     return_df = T, #na.rm = T, 
                     p_thrsld = 0.01, #0.05
                     rownames = c('monocle (FPKM)', 'monocle (readcount)', 'monocle (New size normalization)', 'monocle (New size normalization, Estimate transcript)', 
                        'edgeR (edgeR size normalization)', 'edgeR (New Size normalization)', 'DESeq2 (DESeq2 size normalization)', 'DESeq2 (New Size normalization)',
                        'DESeq (DESeq size normalization)', "DESeq (New Size normalization)", 'SCDE (Read Counts)', 'SCDE (New size normalization)', 
                        'MAST (absolute no normalization)', 'MAST (mc no normalization)', 'MAST (FPKM no normalization)', 'MAST (readcount no normalization)'))
df3$data_type = c("Read counts", "FPKM", "MC transcripts", "Spikein transcripts", "Spikein transcripts", "Read counts", "Spikein transcripts", "Read counts", 

"Spikein transcripts", "Read counts", "Spikein transcripts", "Read counts", 

"MC transcripts", "Spikein transcripts", "Read counts", "FPKM")

df3$class = '3relative'


# only show new size normalization: 
df3.1 <- df3
df3.1[, 'Type'] <- c('MAST', 'MAST', 'MAST', 'MAST', 'SCDE', 'SCDE', 'DESeq1', 'DESeq1', 'DESeq2', 'DESeq2', 'edgeR', 'edgeR', 'Monocle', 'Monocle', 'Monocle', 'Monocle') # geom_bar(stat = 'identity', position = 'dodge') 

pdf('./supplementary_figures/fig2a_si.pdf', width = 3, height = 2)
qplot(factor(Type), value, stat = "identity", geom = 'bar', position = 'dodge', fill = data_type, data = melt(df3.1)) + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + monocle_theme_opts() + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) + nm_theme()
dev.off()

save.image('./RData/cmpr_three_packages.RData')
#use LRT test
# mast_abs_cds_14_18.split<-split(mast_abs_cds_14_18,"Cell")

# two.sample <- LRT(mast_abs_cds_14_18, Time, referent='E14.5')
# car::some(two.sample)

#2. SUBRA (matlab code)


#3. BASiCs (probably not necessary) 
# install_github('catavallejos/BASiCS')



