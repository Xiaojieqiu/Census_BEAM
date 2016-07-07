#scripts for responding the third rebuttal submission: 
# setwd('~/Dropbox (Personal)/Projects/BEAM/')

library(monocle)
# library(devtools)
# load_all('~/Projects/monocle-dev') 
library(xacHelper)

load_all_libraries()
library(sp)
library(igraph)
library(grid)

#load the data: 
# load('./RData/analysis_lung_data.RData')

pData(absolute_cds)$spike_in_RNA <-  esApply(absolute_cds, 2, function(x) sum(x[transcript_num:nrow(absolute_cds)]))
pData(mc_adj_cds)$spike_in_RNA <-  esApply(mc_adj_cds, 2, function(x) sum(x[transcript_num:nrow(mc_adj_cds)]))

#use spike-in transcript counts
#fig 3f_ercc_genes: 
pdf('./main_figures/fig3f_ercc_genes.pdf', width = 2, height = 1.7)
qplot(pData(absolute_cds)$spike_in_RNA[pData(absolute_cds)$endogenous_RNA > 1e3], 
    pData(mc_adj_cds)$spike_in_RNA[pData(absolute_cds)$endogenous_RNA > 1e3], 
    log="xy", color=pData(absolute_cds)$Time[pData(absolute_cds)$endogenous_RNA > 1e3], size = I(1)) + 
   geom_smooth(method="lm", color="black", size = .1) + geom_abline(color="red") +  
  xlab("spike-in transcript counts") +
  ylab("spike-in transcript counts (algorithm)") + #scale_size(range = c(0.25, 0.25)) + 
  scale_color_discrete(name = "Time points") + nm_theme() + xlim(1e4, 2e4) + ylim(1e4, 2e4)
dev.off()

#fig 3g_ercc_genes:Time <- pData(abs_AT12_cds_subset_all_gene)$Time
E14.5_cell <- rowMeans(exprs(abs_AT12_cds_subset_all_gene[transcript_num:nrow(mc_adj_cds), which(Time == 'E14.5')]))
E16.5_cell <- rowMeans(exprs(abs_AT12_cds_subset_all_gene[transcript_num:nrow(mc_adj_cds), which(Time == 'E16.5')]))
E18.5_cell <- rowMeans(exprs(abs_AT12_cds_subset_all_gene[transcript_num:nrow(mc_adj_cds), which(Time == 'E18.5')]))
Adult_cell <- rowMeans(exprs(abs_AT12_cds_subset_all_gene[transcript_num:nrow(mc_adj_cds), which(Time == 'Adult')]))

mc_E14.5_cell <- rowMeans(exprs(mc_adj_cds[transcript_num:nrow(mc_adj_cds), colnames(abs_AT12_cds_subset_all_gene)[which(Time == 'E14.5')]]))
mc_E16.5_cell <- rowMeans(exprs(mc_adj_cds[transcript_num:nrow(mc_adj_cds), colnames(abs_AT12_cds_subset_all_gene)[which(Time == 'E16.5')]]))
mc_E18.5_cell <- rowMeans(exprs(mc_adj_cds[transcript_num:nrow(mc_adj_cds), colnames(abs_AT12_cds_subset_all_gene)[which(Time == 'E18.5')]]))
mc_Adult_cell <- rowMeans(exprs(mc_adj_cds[transcript_num:nrow(mc_adj_cds), colnames(abs_AT12_cds_subset_all_gene)[which(Time == 'Adult')]]))

mc_abs_exprs_df <- data.frame(spikein = as.vector(c(E14.5_cell, E16.5_cell, E18.5_cell, Adult_cell)), 
                              mc_algorithm = as.vector(c(mc_E14.5_cell, mc_E16.5_cell, mc_E18.5_cell, mc_Adult_cell)),
                              cell = rep(c("E14.5_cell", "E16.5_cell", "E18.5_cell", "Adult_cell"), each = length(transcript_num:nrow(mc_adj_cds))))

levels(mc_abs_exprs_df$cell) <- c("Adult", "E14.5", "E16.5", "E18.5")
pdf('./main_figures/fig3g2_ercc_genes.pdf', width = 3, height = 2) #main_figures/
qplot(spikein + 1, mc_algorithm + 1, log = 'xy', 
      color = cell, data = mc_abs_exprs_df, size  = 1.5) + facet_wrap(~cell, scales = 'free', ncol = 2) +  geom_smooth(method = 'rlm', aes(group = 199), size = .1) + 
    scale_size(range = c(1.5, 1)) +  geom_abline(size = .1) + xlab('Transcript counts (Spike-in)') + #scale_size(range = c(0.25, 0.25)) + 
    ylab('Transcript counts (Recovery algorithm)') + nm_theme()
dev.off()

pdf('./main_figures/fig3g_ercc_genes.pdf', width = 3, height = 2)
ggplot(mc_abs_exprs_df) + aes(x=spikein + 1, y= mc_algorithm + 1) + scale_x_log10() + scale_y_log10() + facet_wrap(~cell, scales = 'free', ncol = 2) + 
    xlab('Transcript counts (Spike-in)') + scale_size(range = c(1.5, 1)) + 
    ylab('Transcript counts (Recovery algorithm)')  + 
    #stat_density2d(geom="tile", aes(fill=..density..^1, alpha=1), contour=FALSE) + 
    geom_point(size=0.5, aes(color = cell)) + geom_abline(size = .1) + nm_theme()
# stat_density2d(geom="tile", aes(fill=..density..^1, alpha=ifelse(..density..^1<0.4,0,1)), contour=FALSE) 
dev.off()

#figure 8c: 
#use UMI counts
load('./RData/umi_normalization.RData')

#save.image('./RData/umi_normalization.RData')
pData(UMI_cds)$Total_mRNAs <- esApply(UMI_cds, 2, sum)
pData(UMI_cds)$endogenous_RNA <- esApply(UMI_cds[-ERCC_ids, ], 2, sum)

pdf('./supplementary_figures/UMI_algorithm_mean.pdf', width = 2, height = 2)
qplot(apply(UMI_norm_recovery_all_correct_mc$norm_cds[-ERCC_ids, ], 1, mean), esApply(UMI_cds[-ERCC_ids,  colnames(UMI_TPM_cds)], 1, mean), log = 'xy') + 
    geom_smooth(method = 'rlm') + xlab('Mean endogenous mRNA count \n (algorithm)') + ylab('Mean endogenous mRNA count (regression)') + 
    nm_theme()
dev.off()


#use fpkm for performing all the similar analysis? 
#change all data for downstream analysis to fpkm values: 
#run the fpkm values in the shalek data, etc. 























