library(monocle)
# library(devtools)
# load_all('~/Projects/monocle-dev')
library(xacHelper)

load_all_libraries()

load('./RData/deg_benchmark_analysis.RData')
load('./RData/analysis_lung_data_mc.RData') #load('./RData/analysis_lung_data.RData')
load('./RData/analysis_HSMM_data.RData')
load('./RData/analysis_distribution_fitting.RData')
load('./RData/analysis_shalek_data.RData')
load('./RData/cmpr_three_packages.RData')

##################################################fig1_si############################################################

three_cells_cds <- data.frame(transcript = round(as.vector(exprs(absolute_cds)[1:38919, c("SRR1033854_0", "SRR1033855_0", "SRR1033856_0")])), 
                              read_counts = as.vector(exprs(read_countdata_cds)[1:38919, c("SRR1033854_0", "SRR1033855_0", "SRR1033856_0")]), 
                              cell = rep(c("SRR1033854_0", "SRR1033855_0", "SRR1033856_0"), each = 38919)
)

pdf('./supplementary_figures/fig1a_si.pdf', width = 2, height = 3)
qplot(value, data = melt(three_cells_cds), log = 'x') + geom_histogram(aes(fill = variable)) + facet_wrap(~cell+variable, nrow = 3, scale = 'free') + xlab('') + ylab('Gene number') +  
  theme(strip.background = element_blank(), strip.text.x = element_blank())+ nm_theme()
dev.off()


abs_gd_fit_res <- cal_gd_statistics(abs_gd_fit_df[, c('nb_pvalue', 'zinb_pvalue')], percentage = F, type = 'absolute', gene_list = valid_genes)
readcount_gd_fit_res <- cal_gd_statistics(read_gd_fit_df[, c('nb_pvalue', 'zinb_pvalue')], percentage = F,  type = 'readcount', gene_list = valid_genes)
gd_fit_res <- rbind(abs_gd_fit_res, readcount_gd_fit_res)
gd_fit_res <- cbind(gd_fit_res, data_type = row.names(gd_fit_res))
row.names(gd_fit_res) <- NULL
gd_fit_res <- as.data.frame(gd_fit_res)

gd_fit_res_num <- subset(gd_fit_res, data_type == 'gd_fit_num')
gd_fit_res_success_num <- subset(gd_fit_res, data_type == 'success_fit_num')
# 

#generate the result of goodness of fit for each gene: 
colnames(gd_fit_res_num)[1:2] <- c('NB', 'ZINB')
test <- melt(gd_fit_res_num[, 1:3], id.vars = 'type')
p1 <- qplot(as.factor(variable), as.numeric(as.character(value)), geom = 'bar', stat = 'identity', data = test, fill = type) + facet_wrap('type') + nm_theme() + 
  theme(legend.position = 'none') + xlab('Fit types') + ylab('number of genes') + theme(strip.background = element_blank(),
                                                                                        strip.text.x = element_blank()) + theme(axis.text.x = element_text(angle = 30, hjust = .9))
pdf('./supplementary_figures/fig1b_si.pdf', height = 1.5, width = 1)
p1 + xlab('')
dev.off()

colnames(gd_fit_res_success_num)[1:2] <- c('NB', 'ZINB')
test <- melt(gd_fit_res_success_num[, 1:3], id.vars = 'type')

p2 <- qplot(as.factor(variable), as.numeric(as.character(value)), geom = 'bar', stat = 'identity', data = test, fill = type) + facet_wrap('type') + nm_theme() + 
  theme(legend.position = 'none') + xlab('Fit types') + ylab('number of genes') + theme(strip.background = element_blank(),
                                                                                        strip.text.x = element_blank()) + theme(axis.text.x = element_text(angle = 30, hjust = .9))

pdf('./supplementary_figures/fig1c_si.pdf', height = 1.5, width = 1)
p2 + xlab('')
dev.off()

# ##################################################fig2_si############################################################
#add the benchmark plot with precision/recall/F1 score (for FPR/sensitivity plot, see the very bottom)
#fig2a_si has already generated in the benchmark_analysis.R 

#add the barplot for the overlapping genes: 
q_thrsld <- 0.1
element_all_list <- list(
  intersect(select_genes, names(default_edgeR_p_glm[p.adjust(default_edgeR_p_glm, method = "BH") < q_thrsld])), 
  intersect(select_genes, names(default_deseq2_p[p.adjust(default_deseq2_p, method = "BH") < q_thrsld])), 
  intersect(select_genes, names(mode_size_norm_permutate_ratio_by_geometric_mean[p.adjust(mode_size_norm_permutate_ratio_by_geometric_mean, method = "BH") < q_thrsld]))#, 
  # names(default_deseq_p[default_deseq_p < p_thrsld]), 
  # intersect(select_genes, names(monocle_p_readcount[p.adjust(monocle_p_readcount, method = "BH") < q_thrsld])),  
  # intersect(select_genes, names(scde_p[p.adjust(scde_p, method = "BH") < q_thrsld])) 
  # names(mast_count_pval_no_norm[mast_count_pval_no_norm < p_thrsld])
  )

census_element_all_list <- list(
  intersect(select_genes, names(mc_edgeR_p_glm[p.adjust(mc_edgeR_p_glm, method = "BH") < q_thrsld])), 
  intersect(select_genes, names(mc_default_deseq2_p[p.adjust(mc_default_deseq2_p, method = "BH") < q_thrsld])), 
  intersect(select_genes, names(new_mc_size_norm_monocle_p_ratio_by_geometric_mean[p.adjust(new_mc_size_norm_monocle_p_ratio_by_geometric_mean, method = "BH") < q_thrsld]))#, 
  # intersect(select_genes, names(mode_size_norm_permutate_ratio_by_geometric_mean[p.adjust(mode_size_norm_permutate_ratio_by_geometric_mean, method = "BH") < q_thrsld])), 

  # intersect(select_genes, names(mc_scde_p[mc_scde_p < q_thrsld]))
 
  # names(abs_default_deseq_p[abs_default_deseq_p < p_thrsld]), 
  # names(mode_size_norm_permutate_ratio_by_geometric_mean[mode_size_norm_permutate_ratio_by_geometric_mean < p_thrsld])
  )

abs_element_all_list <- list(
  intersect(select_genes, names(abs_default_edgeR_p_glm[p.adjust(abs_default_edgeR_p_glm, method = "BH") < q_thrsld])), 
  intersect(select_genes, names(abs_default_deseq2_p[p.adjust(abs_default_deseq2_p, method = "BH") < q_thrsld])), 
  intersect(select_genes, names(new_abs_size_norm_monocle_p_ratio_by_geometric_mean[p.adjust(new_abs_size_norm_monocle_p_ratio_by_geometric_mean, method = "BH") < q_thrsld]))#, 
  # names(abs_default_deseq_p[abs_default_deseq_p < p_thrsld]), 
  # intersect(select_genes, names(mode_size_norm_permutate_ratio_by_geometric_mean[p.adjust(mode_size_norm_permutate_ratio_by_geometric_mean, method = "BH") < q_thrsld])), 
  
  # intersect(select_genes, names(abs_scde_p[p.adjust(abs_scde_p, method = "BH") < q_thrsld])
  )
  
readcount_overlap <- Reduce(intersect, element_all_list)
readcount_union <- Reduce(union, element_all_list)
census_overlap <- Reduce(intersect, census_element_all_list)
census_union <- Reduce(union, census_element_all_list)
abs_overlap <- Reduce(intersect, abs_element_all_list)
abs_union <- Reduce(union, abs_element_all_list)

overlap_union_df <- data.frame(Measurement = factor(rep(c('Read counts', 'Census counts', 'Spikein Transcripts'), each = 2), levels = c('Read counts', 'Census counts', 'Spikein Transcripts')), 
			Number = c(length(readcount_overlap), length(setdiff(readcount_union, readcount_overlap)), length(census_overlap), length(setdiff(census_union, census_overlap)), length(abs_overlap), length(setdiff(abs_union, abs_overlap))),
			Type = rep(c('Overlapping', 'Union'), 3))

#levels(overlap_union_df$Measurement) <- c('Read counts', 'Census counts', 'Spikein Transcripts')
cols <- c("Read counts" = "#00BFC4", "Read.counts" = "#00BFC4", "read.counts" = "#00BFC4", "MC transcripts" = "#7CAE00", 'Census counts' = "#7CAE00", 'Census.counts' = "#7CAE00", "Spikein transcripts" = "#C77CFF", "Spikein.transcripts" = "#C77CFF", "transcript.counts" = "#C77CFF", "FPKM" = "#F8766D", 'read_counts' = "#00BFC4", transcript_counts = '#C77CFF')

cols <- c("Read_overlap" = "#7093CB", "Read_union" = "#A1BCE2", "spike_overlap" = "#010101", "spike_union" = "#7B7C7C", 'census_overlap' = "#7BAE41", 'census_union' = "#A0C174")
overlap_union_df$col <- c('Read_overlap', 'Read_union', 'census_overlap', 'census_union', 'spike_overlap', 'spike_union')

pdf('./supplementary_figures/fig2b_three_measurements_fraction.pdf', width = 1, height = 2) 
qplot(Measurement, Number, geom = 'bar', stat = 'identity', position = 'fill', fill = col, data = overlap_union_df)  + 
	theme(axis.text.x = element_text(angle = 30, hjust = 1)) + nm_theme() + xlab('') + scale_fill_manual(values = cols)
dev.off()

pdf('./supplementary_figures/fig2b_three_measurements_number.pdf', width = 1, height = 2) 
qplot(Measurement, Number, geom = 'bar', stat = 'identity', position = 'stack', fill = col, data = overlap_union_df)  + 
	theme(axis.text.x = element_text(angle = 30, hjust = 1)) + nm_theme() + xlab('') + scale_fill_manual(values = cols) + ylab('DE genes')
dev.off()

pdf('./supplementary_figures/fig2b_three_measurements.pdf', width = 4, height = 4) 
qplot(Measurement, Number, geom = 'bar', stat = 'identity', position = 'fill', fill = Type, data = overlap_union_df)  + ylab('DE genes') #+ scale_fill_manual(values = cols) 
dev.off()

overlap_df <- data.frame("read.counts" = length(readcount_overlap), "transcript.counts" = length(abs_overlap)) 
union_df <- data.frame("read.counts" = length(readcount_union), "transcript.counts" = length(abs_union))  

pdf('./supplementary_figures/fig2b.1.pdf', width = 0.75, height = 1.5)
qplot(variable, value, geom = 'bar', stat = 'identity', fill = variable, data = melt(overlap_df)) + xlab('') + ylab('number') + nm_theme() + 
scale_fill_manual(values = cols) + theme(axis.text.x = element_text(angle = 30, hjust = .9)) #+ theme(axis.ticks = element_blank(), axis.text.x = element_blank()) #
dev.off()

pdf('./supplementary_figures/fig2b.2.pdf', width = 0.75, height = 1.5)
qplot(variable, value, geom = 'bar', stat = 'identity', fill = variable, data = melt(union_df)) + xlab('') + 
  ylab('number') + nm_theme() +  scale_fill_manual(values = cols) + theme(axis.text.x = element_text(angle = 30, hjust = .9)) #+ theme(axis.ticks = element_blank(), axis.text.x = element_blank()) #
dev.off()

##UMI data: 
#fig2c_si has already generated in the analysis_UMI_data.R

#################################################fig4_si############################################################
#read the read count data for the genes: 
df <- data.frame(Time = pData(read_countdata_cds)$Time, sum_readcounts = esApply(read_countdata_cds[fData(read_countdata_cds)$biotype == 'spike', ], 2, sum))
# qplot(sum_readcounts, fill = Time, log = 'x') + facet_wrap(~Time)
pdf('./supplementary_figures/fig_4a_si.pdf', height = 3, width = 2)
qplot(sum_readcounts, fill = Time, log = 'x', data = df) + facet_wrap(~Time, ncol = 1, scales = 'free_y') + nm_theme() + xlab('Total reads') + ylab('Cells')
dev.off()

#number of ERCC spike-in detected in each cell
loss_ercc_spikein <- esApply(ercc_controls, 2, function(x) as.numeric(x < 10e-4))

ercc_controls_detected_df <- data.frame(detected = esApply(ercc_controls, 2, function(x) sum(x > 0)), Time = pData(absolute_cds[, colnames(loss_ercc_spikein)])$Time)

pdf('./supplementary_figures/fig_4b_si.pdf', height = 3, width = 2)
qplot(detected, fill = Time, data = ercc_controls_detected_df) + facet_wrap(~Time, ncol = 1) + xlab('Spike-ins detected') + ylab('Cells') + nm_theme()
dev.off()

dir = "./data/Aviv_data/cuffnorm_output_files"
Shalek_sample_table <- read.delim(paste(dir, "/samples.table", sep = ''))
Shalek_norm_count <- read.delim(paste(dir, "/genes.count_table", sep = ''))
row.names(Shalek_norm_count) <- Shalek_norm_count$tracking_id
Shalek_norm_count <- Shalek_norm_count[, -1]

Shalek_read_countdata <- round(t(t(Shalek_norm_count) * Shalek_sample_table$internal_scale)) #convert back to the raw counts 
Shalek_read_countdata <- Shalek_read_countdata[row.names(Shalek_abs), colnames(Shalek_abs)]
colnames(Shalek_read_countdata) <- colnames(Shalek_abs)

Shalek_read_countdata_cds <- newCellDataSet(as.matrix(Shalek_read_countdata),
                                            phenoData = new("AnnotatedDataFrame", data = pData(Shalek_abs)),
                                            featureData = new("AnnotatedDataFrame", data = fData(Shalek_abs)),
                                            expressionFamily = negbinomial(),
                                            lowerDetectionLimit = 1)

pData(Shalek_read_countdata_cds)$Total_mRNAs <- esApply(Shalek_read_countdata_cds, 2, sum)
pData(Shalek_read_countdata_cds)$endogenous_RNA <- esApply(Shalek_read_countdata_cds, 2, sum)

# Shalek_golgi_update <- Shalek_abs[,pData(Shalek_abs)$experiment_name %in% c("LPS_GolgiPlug", "LPS", "Unstimulated_Replicate")] #Shalek_golgi_update is commented in the current analysis, added back here 
Shalek_gene_df <- data.frame(experiment_name = pData(Shalek_read_countdata_cds[, c(colnames(Shalek_abs_subset_ko_LPS))])$experiment_name, #, colnames(Shalek_golgi_update)
                             sum_readcounts = esApply(Shalek_read_countdata_cds[, c(colnames(Shalek_abs_subset_ko_LPS))], 2, sum)) #, colnames(Shalek_golgi_update)


Shalek_gene_df <- data.frame(experiment_name = pData(Shalek_read_countdata_cds[, c(colnames(Shalek_abs_subset_ko_LPS))])$experiment_name, #, colnames(Shalek_golgi_update)
                             sum_readcounts = esApply(Shalek_read_countdata_cds[, c(colnames(Shalek_abs_subset_ko_LPS))], 2, sum)) #, colnames(Shalek_golgi_update)

pdf('./supplementary_figures/fig_4c_si.pdf', height = 1.5, width = 2)
qplot(sum_readcounts, fill = experiment_name, log = 'x', data = Shalek_gene_df) + facet_wrap(~~experiment_name, ncol = 1, scales = 'free_y') + nm_theme() + xlab('Total reads') + ylab('Cells')
dev.off()

################################################fig5_si##############################################################
####generate the SI figures for HSMM data: 

pdf('./supplementary_figures/fig5a_si.pdf', width = 1.5, height = 1.2)
plot_spanning_tree(std_HSMM, color_by="Time", show_backbone=T, backbone_color = 'black',
                   markers=NULL, show_cell_names = F, cell_size = 1, cell_link_size = 0.2) + nm_theme() #+ scale_size(range = c(0.5, .5)) 
dev.off()


pdf('./supplementary_figures/fig5b_si.pdf', width = 1.5, height = 1.2)
plot_spanning_tree(HSMM_myo, color_by="Time", show_backbone=T, backbone_color = 'black',
                   markers=NULL, show_cell_names = F, cell_size = 1, cell_link_size = 0.2) + nm_theme() #+ scale_size(range = c(0.5, .5)) 
dev.off()

plot_tree_pairwise_cor2 <- function (std_tree_cds, absolute_tree_cds) 
{
  maturation_df <- data.frame(cell = rep(colnames(std_tree_cds), 
                                         2), maturation_level = 100 * c(pData(std_tree_cds)$Pseudotime/max(pData(std_tree_cds)$Pseudotime), 
                                                                        pData(absolute_tree_cds)$Pseudotime/max(pData(absolute_tree_cds)$Pseudotime)), 
                              Type = rep(c("FPKM", "Transcript counts (vst)"), each = ncol(std_tree_cds)), rownames = colnames(absolute_tree_cds))
  cor.coeff <- cor(pData(absolute_tree_cds)$Pseudotime, pData(std_tree_cds)$Pseudotime, 
                   method = "spearman")
  message(cor.coeff)
  p <- ggplot(aes(x = maturation_level, y = Type, group = cell), 
              data = maturation_df) + geom_point(size = 1) + geom_line(color = "blue", alpha = .3) + 
    xlab("Pseudotime") + ylab("Type of tree construction") + monocle_theme_opts()
  return(p)
}

std_HSMM <- orderCells(std_HSMM, reverse = F) #reverse the tree
pdf('./supplementary_figures/fig5c_si.pdf', width = 4, height = 2)
plot_tree_pairwise_cor2(std_HSMM, HSMM_myo) + nm_theme()
dev.off()

#muscle_std_glm_perm_pval, muscle_size_normalized_mc_glm_perm_pval
HSMM_muscle_abs_pseudotime_deg <- intersect(row.names(HSMM_myo_size_norm_res[HSMM_myo_size_norm_res$pval <0.05, ]), names(muscle_size_normalized_mc_glm_perm_pval[!is.na(muscle_size_normalized_mc_glm_perm_pval)]))
HSMM_muscle_fpkm_pseudotime_deg <- intersect(row.names(std_HSMM_myo_pseudotime_res_ori[std_HSMM_myo_pseudotime_res_ori$pval <0.05, ]), names(muscle_std_glm_perm_pval[!is.na(muscle_std_glm_perm_pval)]))
element_all <- c(HSMM_muscle_abs_pseudotime_deg, 
                 HSMM_muscle_fpkm_pseudotime_deg)
sets_all <- c(rep(paste('Transcript counts (Size + VST)', sep = ''), length(HSMM_muscle_abs_pseudotime_deg)), 
              rep(paste('FPKM', sep = ''), length(HSMM_muscle_fpkm_pseudotime_deg)))
table(sets_all) #number of genes

pdf('./supplementary_figures/fig5d_si.pdf')
venneuler_venn(element_all, sets_all)
dev.off()
table(sets_all) #number of genes

# #precision: 
# length(intersect(row.names(subset(HSMM_myo_size_norm_res, pval <= 0.05)), names(muscle_size_normalized_mc_glm_perm_pval[which(muscle_size_normalized_mc_glm_perm_pval < 0.05)]))) / 
#       length(row.names(subset(HSMM_myo_size_norm_res, pval <= 0.05)))
# length(intersect(row.names(subset(HSMM_myo_size_norm_res, pval <= 0.05)), names(muscle_std_glm_perm_pval[which(muscle_std_glm_perm_pval < 0.05)]))) / 
#       length(row.names(subset(HSMM_myo_size_norm_res, pval <= 0.05)))
# #recall: 
# length(intersect(row.names(subset(HSMM_myo_size_norm_res, pval <= 0.05)), names(muscle_size_normalized_mc_glm_perm_pval[which(muscle_size_normalized_mc_glm_perm_pval < 0.05)]))) / 
#       length(names(muscle_size_normalized_mc_glm_perm_pval[which(muscle_size_normalized_mc_glm_perm_pval < 0.05)]))
# length(intersect(row.names(subset(HSMM_myo_size_norm_res, pval <= 0.05)), names(muscle_std_glm_perm_pval[which(muscle_std_glm_perm_pval < 0.05)]))) / 
#       length(names(muscle_std_glm_perm_pval[which(muscle_std_glm_perm_pval < 0.05)]))

# #F1 score: 

names(muscle_size_normalized_mc_glm_perm_pval < 0.05)

cols <- c("FPKM" = "#F2756D", "Read counts" = "#6F94CC", "Transcript counts" = "#000202", "Estimated transcript counts" = "#7AAE41") ##A680B9
muscle_df$data_type = c("Estimated transcript counts", "Estimated transcript counts", "FPKM", "FPKM")

muscle_df$class = '3relative'
muscle_df.1 <- muscle_df
muscle_df.1 <- muscle_df[1:2, ]

muscle_df[, 'Type'] <- c('Monocle', 'DESeq', 'DESeq', 'Monocle')
colnames(muscle_df)[1:3] <- c('Precision', 'Recall', 'F1')

pdf('./supplementary_figures/fig5e_si_ori.pdf', width = 1.8, height = 1.5)
qplot(factor(Type), value, stat = "identity", geom = 'bar', position = 'dodge', fill = data_type, data = melt(muscle_df), log = 'y') + #facet_wrap(~variable) + 
  ggtitle(title) + scale_fill_manual(values = cols, name = "Type")  + xlab('') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
  ggtitle('') + monocle_theme_opts() + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) + ylim(0, 1) + 
  theme(strip.background = element_blank(), strip.text.x = element_blank()) + nm_theme() + theme(strip.background = element_blank(), strip.text.x = element_blank())
dev.off()

pdf('./supplementary_figures/fig5e_si_ori_helper.pdf', width = 1.8, height = 1.5)
qplot(factor(Type), value, stat = "identity", geom = 'bar', position = 'dodge', fill = data_type, data = melt(muscle_df), log = 'y') + #facet_wrap(~variable) + 
  ggtitle(title) + scale_fill_manual(values = cols, name = "Type")  + xlab('') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
  ggtitle('') + monocle_theme_opts() + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) + ylim(0, 1) + 
  theme(strip.background = element_blank(), strip.text.x = element_blank()) + theme(strip.background = element_blank(), strip.text.x = element_blank())
dev.off()

####4. enrichment for GO/motif analysis for the muscle permuation test ####
abs_HSMM_ids <- row.names(HSMM_myo_size_norm_res[HSMM_myo_size_norm_res$pval <0.05, ]) 
std_HSMM_ids <- row.names(std_HSMM_myo_pseudotime_res_ori[std_HSMM_myo_pseudotime_res_ori$pval <0.05, ])

abs_HSMM_name <- fData(HSMM_myo[abs_HSMM_ids, ])$gene_short_name
std_HSMM_name <- fData(HSMM_myo[std_HSMM_ids, ])$gene_short_name

write.table(std_HSMM_name, file = '../std_muscle_pseudotime_gene.txt', sep = '\t', row.names = F, quote = F)
write.table(abs_HSMM_name, file = '../abs_muscle_pseudotime_gene.txt', sep = '\t', row.names = F, quote = F)

abs_gsaRes <- runGSAhyper(unique(abs_HSMM_name), gsc=human_go_gsc, universe=unique(as.vector(fData(HSMM_myo)$gene_short_name)))
std_gsaRes <- runGSAhyper(unique(std_HSMM_name), gsc=human_go_gsc, universe=unique(as.vector(fData(HSMM_myo)$gene_short_name)))

abs_gsaRes_reactome <- runGSAhyper(unique(abs_HSMM_name), gsc=human_reactome_gsc, universe=unique(as.vector(fData(HSMM_myo)$gene_short_name)))
std_gsaRes_reactome <- runGSAhyper(unique(std_HSMM_name), gsc=human_reactome_gsc, universe=unique(as.vector(fData(HSMM_myo)$gene_short_name)))

#plot the terms and significance: 
gsa_results <- list('Transcript_counts' = abs_gsaRes, 'FPKM' = std_gsaRes)
plot_gsa_hyper_heatmap(HSMM_myo, gsa_results, significance=1e-3)

gsa_results_reactome <- list('Transcript_counts' = abs_gsaRes_reactome, 'FPKM' = std_gsaRes_reactome)
plot_gsa_hyper_heatmap(HSMM_myo, gsa_results_reactome, significance=1e-3)

benchmark_pseudotime_test <- function(abs_gsaRes , std_gsaRes ) {
    abs_hyper_df_all <- data.frame(gene_set = names(abs_gsaRes$pvalues), pval = abs_gsaRes$pvalues, qval = abs_gsaRes$p.adj)
    abs_hyper_df_all$qval <- p.adjust(abs_hyper_df_all$pval, method = 'fdr')
    colnames(abs_hyper_df_all)[1] <- "cluster_id"
    
    std_hyper_df_all <- data.frame(gene_set = names(std_gsaRes$pvalues), pval = std_gsaRes$pvalues, qval = std_gsaRes$p.adj)
    std_hyper_df_all$qval <- p.adjust(std_hyper_df_all$pval, method = 'fdr')
    colnames(std_hyper_df_all)[1] <- "cluster_id"
    
    abs_hyper_df <- subset(abs_hyper_df_all, abs_hyper_df_all[, 'qval'] <= 0.01)
    std_hyper_df <- subset(std_hyper_df_all, std_hyper_df_all[, 'qval'] <= 0.01)
    
    gsa_results <- list('1' = abs_gsaRes, '2' = std_gsaRes)
    
    #performance comparision based on Cole's idea: 
    element_all_list <- c(
        abs_hyper_df$cluster_id, 
        std_hyper_df$cluster_id
    )
    
    sets_all <- c(
        rep(paste('transcript counts', sep = ''), nrow(abs_hyper_df)),
        rep(paste('FPKM values', sep = ''), nrow(std_hyper_df))
    )
    
    #show the number for the overlapping:
    print(table(sets_all))
    print(intersect(abs_hyper_df$cluster_id, std_hyper_df$cluster_id))

    pdf(file = './supplementary_figures/overlap_enriched_muscle_term.pdf')
    venneuler_venn(element_all_list, sets_all)
    dev.off()
    
    #p-val plot show the gene names relevant for muscle differentiation 
    abs_hyper_df_all$cluster_id
    muscle_term_ids <- c(grep(pattern = 'MUSCLE', abs_hyper_df_all$cluster_id, ignore.case = T), grep(pattern = 'Myogenesis', abs_hyper_df_all$cluster_id, ignore.case = T))
    # muscle_term_ids <- Reduce(intersect, list(which(abs_hyper_df_all$qval < 0.01), which(std_hyper_df_all$qval < 0.01), muscle_term_ids))
    Types <- rep('Term without muscle', nrow(abs_hyper_df_all))
    Types[muscle_term_ids] <- 'Term with muscle'
    muscle_gs_df <- data.frame(cluster_id = abs_hyper_df_all$cluster_id, abs = abs_hyper_df_all$pval, std = std_hyper_df_all$pval, Type = Types) #muscle related/ non-muscle related 
    muscle_gs_df$ratio <- muscle_gs_df$std /  muscle_gs_df$abs
    valid_muscle_gs_df <- muscle_gs_df[c(muscle_term_ids, which(muscle_gs_df$abs < 0.01 & muscle_gs_df$std < 0.01 )), ]
    pdf(file = './supplementary_figures/muscle_pseudotime_benchmark_qval.pdf', width = 2.5, height = 1.5)
    ggplot(data = valid_muscle_gs_df, aes(Type, log(ratio))) + 
        geom_boxplot(aes(fill = Type), alpha = 0.3, outlier.size = 0.25, lwd = 0.25, fatten = 0.5) + 
        geom_jitter(size = 0.5) + nm_theme() + geom_vline(xintercept = 0) + xlab('log(P(FPKM) / P(transcript counts))') + ylab('') + scale_size(range = c(0.25, 0.5))
    dev.off()
    
    pdf(file = './supplementary_figures/muscle_pseudotime_benchmark_qval2.pdf')
    qplot(log(ratio), data = muscle_gs_df[c(muscle_term_ids, which(muscle_gs_df$abs < 0.01 & muscle_gs_df$std < 0.01 )), ], 
          fill = Type, geom = 'density', log = 'x', alpha = 0.3) + nm_theme() + geom_vline(xintercept = 0) + xlab('log(P(FPKM) / P(transcript counts))') + ylab('')
    #qplot(abs, std, data = muscle_gs_df[c(muscle_term_ids, which(muscle_gs_df$abs < 0.01 | muscle_gs_df$std < 0.01 )), ], color = Type, log = 'xy') + nm_theme() + xlab('transcript counts') + ylab('FPKM')
    dev.off()

    return(valid_muscle_gs_df)
}

valid_muscle_gs_df <- benchmark_pseudotime_test(abs_gsaRes_reactome, std_gsaRes_reactome)
pdf(file = './supplementary_figures/muscle_pseudotime_benchmark_qval.pdf', width = 2.5, height = 1.5)
ggplot(data = valid_muscle_gs_df, aes(Type, log(ratio))) + 
    geom_boxplot(aes(fill = Type), alpha = 0.3, outlier.size = 0.25, lwd = 0.25, fatten = 0.5) + 
    geom_jitter(size = 0.5) + nm_theme() + geom_vline(xintercept = 0) + xlab('log(P(FPKM) / P(transcript counts))') + ylab('') + scale_size(range = c(0.25, 0.5))
dev.off()

##############################################################################################################
#make the tree plot with quake annotation:     
pData(AT12_cds_subset_all_gene)$Cell_type <- as.character(pData(AT12_cds_subset_all_gene)$Cell_type)
pData(AT12_cds_subset_all_gene)$Cell_type[c(14, 67:183)] <- 'no_avail'

lung_custom_color_scale_plus_states <- c('no_avail' = 'gray', 'AT1' = '#40A43A', 'AT2' = '#CB1B1E', 'BP' = '#3660A5', 'bulk' = 'black')

pdf('./supplementary_figures/fig6_si.pdf', height = 2, width = 2.5)
plot_spanning_tree(AT12_cds_subset_all_gene, color_by="Cell_type", show_backbone=T, backbone_color = 'black', 
                   markers=NULL, show_cell_names = F, cell_link_size = 0.2) + scale_size(range = c(0.1, 2.5)) + 
  scale_color_manual(values=lung_custom_color_scale_plus_states) + nm_theme() + scale_y_reverse()
dev.off()


pdf('./tmp/fig6_helper.pdf', height = 2, width = 2.5)
plot_spanning_tree(AT12_cds_subset_all_gene, color_by="Cell_type", show_backbone=T, backbone_color = 'black', 
                   markers=NULL, show_cell_names = F, cell_link_size = 0.2) + scale_size(range = c(0.1, 2.5)) + 
  scale_color_manual(values=lung_custom_color_scale_plus_states) 
dev.off()

#make kinetic plots: 
# Nkx2-1, Hopx, Sox9, Foxa2,  and Gata6
mc_AT12_cds_subset_all_gene@dim_reduce_type <- 'ICA'
important_tf_ids <- row.names(subset(fData(mc_AT12_cds_subset_all_gene), gene_short_name %in% c('Hopx', 'Sox9', 'Foxa2', 'Gata6', 'Nkx2-1')))

new_cds <- buildLineageBranchCellDataSet(mc_AT12_cds_subset_all_gene[1:10, ], lineage_labels = c('AT1', 'AT2'))

colour_cell <- rep(0, length(new_cds$Lineage))
names(colour_cell) <- as.character(new_cds$Time)
colour_cell[names(colour_cell) == 'E14.5'] <- "#7CAF42"
colour_cell[names(colour_cell) == 'E16.5'] <- "#00BCC3"
colour_cell[names(colour_cell) == 'E18.5'] <- "#A680B9"
colour_cell[names(colour_cell) == 'Adult'] <- "#F3756C"

colour <- rep(0, length(new_cds$Lineage))
names(colour) <- as.character(new_cds$Lineage)
colour[names(colour) == 'AT1'] <- AT1_Lineage
colour[names(colour) == 'AT2'] <- AT2_Lineage

pdf('./supplementary_figures/fig6b_si_kinetic_plots.pdf', height = 1.5, width = 5)
plot_genes_branched_pseudotime2(mc_AT12_cds_subset_all_gene[important_tf_ids, ], cell_color_by = "Time", lineage_labels = c('AT1', 'AT2'), 
                                trajectory_color_by = "Lineage", trend_formula = '~sm.ns(Pseudotime, df = 3)*Lineage', normalize = F, stretch = T,
                                cell_size = 1, ncol = 3, reducedModelFormulaStr = "~ sm.ns(Pseudotime, df=3)", add_pval = T) + xlab('Pseudotime') + 
  ylab('Transcript counts') + nm_theme()
dev.off()

##############################################################################################################

#FDR, sensitivity: 
plot_fdr_sensitivity <- function (test_p_list = list(monocle_p = monocle_p, deseq_p = deseq_p,
                                                     scde_p = scde_p), permutate_pval, gene_list, na.rm = T, p_thrsld = 0.05,
                                  title = "Comparison of the two-group DEG tests", rownames = NULL,
                                  fill = NULL, return_df = F)
{
  if (is.list(permutate_pval)) {
    if (length(permutate_pval) != length(test_p_list))
      stop("Error: if using permuation pval for each bar, make sure the permuation list matches with the DEG pval list")
    fdr_sensitivity_df <- NULL
    for (ind in 1:length(permutate_pval)) {
      gene_list_true_data_list <- gene_list_true_data(p_thrsld = p_thrsld,
                                                      permutate_pval = permutate_pval[[ind]][gene_list],
                                                      na.rm = na.rm)
      gene_list_new <- gene_list_true_data_list$gene_list
      true_data <- gene_list_true_data_list$true_data
      test_p_vec <- test_p_list[[ind]][gene_list_new]
      TF_PN <- TF_PN_vec(true_data, test_p_vec)
      mc_abs_two_grp <- fdr_sensitivity(test_p_vec, permutate_pval[[ind]][gene_list],
                                        TF_PN)
      fdr_sensitivity_df <- rbind(mc_abs_two_grp, fdr_sensitivity_df)
    }
  }
  else {
    gene_list_true_data_list <- gene_list_true_data(p_thrsld = p_thrsld,
                                                    permutate_pval = permutate_pval[gene_list])
    gene_list <- gene_list_true_data_list$gene_list
    true_data <- gene_list_true_data_list$true_data
    test_p_list <- lapply(test_p_list, function(x) x[gene_list])
    pre_rec_f1_list <- lapply(test_p_list, function(x, true_data,
                                                    permutate_pval) {
      TF_PN <- TF_PN_vec(true_data, x)
      mc_abs_two_grp <- fdr_sensitivity(x, permutate_pval, TF_PN)
    }, true_data = true_data, permutate_pval = permutate_pval[gene_list])
    fdr_sensitivity_df <- do.call(rbind.data.frame, pre_rec_f1_list)
  }
  print(fdr_sensitivity_df)
  if (!is.null(rownames)) {
    row.names(fdr_sensitivity_df) <- rownames
    names(test_p_list) <- rownames
  }
  print(row.names(fdr_sensitivity_df))
  if (is.list(permutate_pval))
    fdr_sensitivity_df$Type <- rev(names(test_p_list))
  else fdr_sensitivity_df$Type <- row.names(fdr_sensitivity_df)
  print(fdr_sensitivity_df)
  if (return_df)
    return(fdr_sensitivity_df)
  if (is.null(fill))
    p1 <- qplot(factor(Type), value, stat = "identity", geom = "bar",
                position = "dodge", fill = Type, data = melt(fdr_sensitivity_df)) +
    facet_wrap(~variable) + ggtitle(title) + scale_fill_discrete("Type") +
    xlab("Type") + ylab("")
  else p1 <- qplot(factor(Type), value, stat = "identity",
                   geom = "bar", position = "dodge", fill = fill, data = melt(fdr_sensitivity_df)) +
    facet_wrap(~variable) + ggtitle(title) + scale_fill_discrete("Type") +
    xlab("Type") + ylab("")
  return(p1)
}

fdr_sensitivity_cal <- function (est_pval, true_pval, type = c("precision", "recall", "fpr", 'sensitivity', 'fdr'),
                                 q_thrsld = 0.1)
{
  qval <- p.adjust(est_pval, method = "BH")
  names(qval) <- names(est_pval)
  true_qval <- p.adjust(true_pval, method = "BH")
  fp <- setdiff(names(qval[qval <= q_thrsld]), names(true_pval[true_qval <=
                                                                 q_thrsld]))
  tp <- intersect(names(qval[qval <= q_thrsld]), names(true_qval[true_qval <=
                                                                   q_thrsld]))
  fn <- setdiff(names(qval[qval > q_thrsld]), names(true_qval[true_qval >
                                                                q_thrsld]))
  tn <- intersect(names(qval[qval >= q_thrsld]), names(true_qval[true_qval >=
                                                                   q_thrsld]))
  if (type == "precision")
    length(tp)/length(union(fp, tp))
  else if (type == "recall")
    length(tp)/length(union(tp, fn))
  else if (type == 'fpr')
    length(fp) / length(union(fp, tn))
  else if (type == 'sensitivity')
    length(tp) / length(union(tp, fn))
  else if (type == 'fdr')
    length(fp) / length(union(tp, fp))
}

fdr_sensitivity <- function (est_pval, true_pval, TF_PN_vec, q_thrsld = 0.1, beta = 1) {
  fdr <- fdr_sensitivity_cal(est_pval, true_pval, type = c("fdr"),
                             q_thrsld = 0.1)
  sensitivity <- fdr_sensitivity_cal(est_pval, true_pval, type = c("sensitivity"),
                                     q_thrsld = 0.1)
  data.frame(fdr = fdr, sensitivity = sensitivity)
}

fdr_sensitivity_cal(monocle_p_readcount, readcount_permutate_pval, 'fdr')
fdr_sensitivity_cal(new_abs_size_norm_monocle_p_ratio_by_geometric_mean, mode_size_norm_permutate_ratio_by_geometric_mean, 'fdr')
fdr_sensitivity_cal(default_deseq2_p, readcount_permutate_pval, 'fdr')
fdr_sensitivity_cal(abs_default_deseq2_p, mode_size_norm_permutate_ratio_by_geometric_mean, 'fdr')

fdr_sensitivity_df <- plot_fdr_sensitivity(test_p_list = list(monocle_p = monocle_p, 
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
                                           permutate_pval = list(monocle_p = mode_size_norm_permutate_ratio_by_geometric_mean, #std_permutate_pval, #readcount_permutate_pval, #std_permutate_pval, 
                                                                 monocle_p_readcount = mode_size_norm_permutate_ratio_by_geometric_mean, #readcount_permutate_pval, 
                                                                 mode_size_norm_permutate_ratio_by_geometric_mean = mode_size_norm_permutate_ratio_by_geometric_mean, #mode_size_norm_permutate_ratio_by_geometric_mean,
                                                                 mc_mode_size_norm_permutate_ratio_by_geometric_mean = mode_size_norm_permutate_ratio_by_geometric_mean, #mc_mode_size_norm_permutate_ratio_by_geometric_mean,
                                                                 default_edgeR_p = mode_size_norm_permutate_ratio_by_geometric_mean, #readcount_permutate_pval, #readcount_permutate_pval use readcount instead of std_permutate_pval
                                                                 abs_default_edgeR_p = mode_size_norm_permutate_ratio_by_geometric_mean, #mode_size_norm_permutate_ratio_by_geometric_mean, 
                                                                 default_deseq2_p = mode_size_norm_permutate_ratio_by_geometric_mean, #readcount_permutate_pval, #readcount_permutate_pval use readcount instead of std_permutate_pval
                                                                 abs_default_deseq2_p = mode_size_norm_permutate_ratio_by_geometric_mean, #mode_size_norm_permutate_ratio_by_geometric_mean, 
                                                                 default_deseq_p = mode_size_norm_permutate_ratio_by_geometric_mean, #readcount_permutate_pval, #readcount_permutate_pval use readcount instead of std_permutate_pval
                                                                 abs_default_deseq_p = mode_size_norm_permutate_ratio_by_geometric_mean, #mode_size_norm_permutate_ratio_by_geometric_mean, 
                                                                 # abs_default_deseq_p_new_norm = mode_size_norm_permutate_ratio_by_geometric_mean, 
                                                                 scde_p = mode_size_norm_permutate_ratio_by_geometric_mean, #readcount_permutate_pval, 
                                                                 abs_scde_p = mode_size_norm_permutate_ratio_by_geometric_mean, #mode_size_norm_permutate_ratio_by_geometric_mean, 
                                                                 mast_abs_pval_no_norm = mode_size_norm_permutate_ratio_by_geometric_mean, #mode_size_norm_permutate_ratio_by_geometric_mean, 
                                                                 mast_mc_pval_no_norm = mode_size_norm_permutate_ratio_by_geometric_mean, #mc_mode_size_norm_permutate_ratio_by_geometric_mean, 
                                                                 mast_std_pval_no_norm = mode_size_norm_permutate_ratio_by_geometric_mean, #std_permutate_pval, 
                                                                 mast_count_pval_no_norm = mode_size_norm_permutate_ratio_by_geometric_mean), #readcount_permutate_pval),
                                           row.names(absolute_cds), #gene_list, overlap_genes, high_gene_list
                                           return_df = T, #na.rm = T, 
                                           p_thrsld = 0.05, #0.05
                                           rownames = c('monocle (FPKM)', 'monocle (readcount)', 'monocle (New size normalization)', 'monocle (New size normalization, Estimate transcript)', 
                                                        'edgeR (edgeR size normalization)', 'edgeR (New Size normalization)', 'DESeq2 (DESeq2 size normalization)', 'DESeq2 (New Size normalization)',
                                                        'DESeq (DESeq size normalization)', "DESeq (New Size normalization)", 'SCDE (Read Counts)', 'SCDE (New size normalization)', 
                                                        'MAST (absolute no normalization)', 'MAST (mc no normalization)', 'MAST (FPKM no normalization)', 'MAST (readcount no normalization)'))
fdr_sensitivity_df$data_type = c("Read counts", "FPKM", "MC transcripts", "Spikein transcripts", "Spikein transcripts", "Read counts", "Spikein transcripts", "Read counts", 
                                 
                                 "Spikein transcripts", "Read counts", "Spikein transcripts", "Read counts", 
                                 
                                 "MC transcripts", "Spikein transcripts", "Read counts", "FPKM")

fdr_sensitivity_df$class = '3relative'

fdr_sensitivity_df[, 'Type'] <- c('MAST', 'MAST', 'MAST', 'MAST', 'SCDE', 'SCDE', 'DESeq1', 'DESeq1', 'DESeq2', 'DESeq2', 'edgeR', 'edgeR', 'Monocle', 'Monocle', 'Monocle', 'Monocle') # geom_bar(stat = 'identity', position = 'dodge') 


tmp <- data.frame(Type = c('SCDE', 'SCDE', 'DESeq1', 'DESeq1', 'DESeq2', 'DESeq2', 'edgeR', 'edgeR'), 
                  data_type = c('MC transcripts', 'FPKM', 'MC transcripts', 'FPKM', 'MC transcripts', 'FPKM', 'MC transcripts', 'FPKM'),
                  class = '3relative', 
                  fpr = NA, sensitivity = NA)
df_res <- rbind(fdr_sensitivity_df, tmp) 
df_res <-  melt(df_res)

pdf('./tmp/cmpr_fpr_sensitivity.pdf', width = 3, height = 2)
ggplot(aes(factor(Type), value), data = df_res) + geom_bar(stat = "identity", position=position_dodge(), aes(fill = data_type)) + #facet_wrap(~variable) + 
  ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
  ggtitle('') + monocle_theme_opts() + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) + nm_theme()
dev.off()

pdf('./tmp/cmpr_fpr_sensitivity_helper.pdf', width = 3, height = 2)
ggplot(aes(factor(Type), value), data = df_res) + geom_bar(stat = "identity", position=position_dodge(), aes(fill = data_type)) + #facet_wrap(~variable) + 
  ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
  ggtitle('') + monocle_theme_opts() + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) + nm_theme()
dev.off()


#generate the figure for calculating capture efficiency: 
Time <- pData(absolute_cds)$Time
ercc_exprs <- exprs(ercc_controls)

#fit the data: 
detected_ercc_spikein <- apply(ercc_exprs, 2, function(x) as.numeric(x > 10e-4))

#total number of detected spikein for each cell
num_time_spike_in_detect_df <- data.frame(Time = rep(unique(pData(absolute_cds)$Time), each = length(unique(spike_df$rounded_numMolecules))), 
                                          rounded_numMolecules = rep(sort(unique(spike_df$rounded_numMolecules)), 4),
                                          detected = 0, tau = 0)

#calculate number of observed cases for a particular spikein in a particular time point,
#also calculate tau, the probability for a particular spikein to be detected 
for(t in unique(Time)){
  for(round_numMolecules in sort(unique(spike_df$rounded_numMolecules))){
    print(t)
    print(round_numMolecules)
    ind <- num_time_spike_in_detect_df$Time == t & num_time_spike_in_detect_df$rounded_numMolecules == round_numMolecules
    num_time_spike_in_detect_df[ind, 3] <- 
      sum(detected_ercc_spikein[spike_df$rounded_numMolecules == round_numMolecules, Time == t])
    
    num_spikein <- sum(spike_df$rounded_numMolecules == round_numMolecules)
    num_time_spike_in_detect_df[ind, 4] <- num_time_spike_in_detect_df[ind, 3] / (table(Time)[t] * num_spikein)
  }
  
}

#write the objective function: 
optim_p_from_tau <- function(x, tau, rounded_numMolecules) {
  res <- 0
  for(i in 1:length(tau)){
    tmp <- (tau[i] - (1 - (1 - x)^rounded_numMolecules[i]))^2
    res <- res + tmp
  }
  res
}

optim_p <- function(num_time_spike_in_detect_df, Time, rounded_numMolecules) {
  tau <- num_time_spike_in_detect_df[num_time_spike_in_detect_df$Time == Time, 4]
  message('initial guess for tau, ', tau[2])
  optim_res <- optim(par = c(x = tau[2]), optim_p_from_tau,
                     gr = NULL, tau = tau[2:9],
                     rounded_numMolecules = rounded_numMolecules[2:9],
                     method = c("Brent"), 
                     lower = c(0), 
                     upper = c(1),
                     hessian = FALSE)
  optim_p_val <- optim_res$par
  predict_tau <- rep(0, length(tau))
  optim_predict_tau <- rep(0, length(tau))
  empirical_p <- rep(0, length(tau))
  
  for(i in 1:length(tau)) {
    optim_predict_tau[i] <- (1 - (1 - optim_p_val)^rounded_numMolecules[i])
  }
  
  for(i in 1:length(tau)) {
    predict_tau[i] <- (1 - (1 - tau[2])^rounded_numMolecules[i])
  }
  
  for(i in 1:length(tau)) {
    empirical_p[i] <- (1 - (1 - tau)^(1/rounded_numMolecules[i]))
  }
  
  list(optim_p_val = optim_p_val, optim_predict_tau = optim_predict_tau, predict_tau = predict_tau, empirical_p = empirical_p)
}

#show the results: 
optim_res_14 <- optim_p(num_time_spike_in_detect_df, 'E14.5', unique(num_time_spike_in_detect_df$rounded_numMolecules))
optim_res_16 <- optim_p(num_time_spike_in_detect_df, 'E16.5', unique(num_time_spike_in_detect_df$rounded_numMolecules))
optim_res_18 <- optim_p(num_time_spike_in_detect_df, 'E18.5', unique(num_time_spike_in_detect_df$rounded_numMolecules))
optim_res_adult <- optim_p(num_time_spike_in_detect_df, 'Adult', unique(num_time_spike_in_detect_df$rounded_numMolecules))

#plot the results: 
#detected number for each trnscripts 
num_time_spike_in_detect_df$optim_predicted_tau <- c(optim_res_18$optim_predict_tau, optim_res_14$optim_predict_tau, optim_res_adult$optim_predict_tau, optim_res_16$optim_predict_tau)
num_time_spike_in_detect_df$predicted_tau <- c(optim_res_18$predict_tau, optim_res_14$predict_tau, optim_res_adult$predict_tau, optim_res_16$predict_tau)
num_time_spike_in_detect_df$optim_p_val <- rep(c(optim_res_18$optim_p_val, optim_res_14$optim_p_val, optim_res_adult$optim_p_val, optim_res_16$optim_p_val), each = length(unique(spike_df$rounded_numMolecules)))
num_time_spike_in_detect_df$empirical_p <- c(optim_res_18$empirical_p, optim_res_14$empirical_p, optim_res_adult$empirical_p, optim_res_16$empirical_p)

ggplot(aes(rounded_numMolecules, tau), data = num_time_spike_in_detect_df) + geom_point(aes(color = Time), size = 4)+ geom_line(aes(rounded_numMolecules, optim_predicted_tau, color = Time)) + scale_x_log10() + facet_wrap(~Time + optim_p_val)
ggplot(aes(rounded_numMolecules, tau), data = num_time_spike_in_detect_df) + geom_point(aes(color = Time), size = 4)+ geom_line(aes(rounded_numMolecules, predicted_tau, color = Time)) + scale_x_log10() + facet_wrap(~Time)

pdf('./supplementary_figures//sequencing_efficiency.pdf', width = 2, height = 1.4)
ggplot(aes(rounded_numMolecules, tau), data = num_time_spike_in_detect_df) + geom_point(aes(color = Time), alpha = 0.5, size = 1) + 
    geom_line(aes(rounded_numMolecules, optim_predicted_tau, color = Time), size = 0.2) + scale_x_log10() + nm_theme() + xlab('Spike-in molecules') + ylab('Probability of observation (P)”')
dev.off()

pdf('./supplementary_figures//sequencing_efficiency_helper.pdf', width = 2, height = 1.4)
ggplot(aes(rounded_numMolecules, tau), data = num_time_spike_in_detect_df) + geom_point(aes(color = Time), size = 1) + 
  geom_line(aes(rounded_numMolecules, optim_predicted_tau, color = Time), size = 0.2) + scale_x_log10() 
dev.off()

#num_time_spike_in_detect_df$Time <- factor(num_time_spike_in_detect_df$Time, levels = c('E14.5', 'E16.5', 'E18.5', 'Adult'))
pdf('./supplementary_figures/capture_rate.pdf', width = 2, height = 1.4)
ggplot(aes(Time, unique(optim_p_val)), data = unique(num_time_spike_in_detect_df[, c('Time', 'optim_p_val')])) + geom_bar(stat = 'identity', aes(fill = unique(Time))) + ylab('Estimated capture rate') + nm_theme()
dev.off()

# #generate the perfect recovery figure: 
# test <- mapply(function(cell_dmode, model) {
#   predict(model, newdata = data.frame(log_fpkm = log10(cell_dmode)), type = 'response')
# }, as.list(estimate_t(exprs(standard_cds)[1:transcript_num, ])), molModels_select) #molModels_select

# df <- pData(absolute_cds)
# df$mode_transcript <- 10^test
# df$estimate_mode <- estimate_t(exprs(standard_cds))

# abs_mat_1 <- relative2abs(standard_cds, verbose = T, return_all = T, reads_per_cell = colSums(read_countdata), expected_total_mRNAs = exp(dmode(log(pData(absolute_cds)$endo ))), use_fixed_intercept = T, expected_mRNA_mode = df$mode_transcript,
#                           expected_capture_rate = c(rep(0.2326, sum(Time == 'E18.5')), rep(0.2987, sum(Time == 'E14.5')), rep(0.3093, sum(Time == 'Adult')), rep(0.09174, sum(Time == 'E16.5'))), 
#                           cores = detectCores(), kb_intercept = 2.268025, kb_slope = -3.655828, kb_slope_rng = c(-3.665828, -3.645828))
# script.dir <- dirname(sys.frame(1)$ofile)
# source(paste(script.dir, '/second_revision/transcript_count_cmpr.R', sep = ''))

# file.rename(c('./main_figures/fig3f_new.pdf', './main_figures/fig3g2_new.pdf'), c('./main_figures/fig3f_perferct_recovery.pdf', './main_figures/fig3g_perferct_recovery.pdf'))
save.image('./RData/gen_supplementary_figure.RData') 
