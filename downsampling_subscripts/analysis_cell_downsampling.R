#######################################
# Downsampling the number of cells
# ########################################
library(monocle)
# library(devtools)
# load_all('~/Projects/monocle-dev')
library(xacHelper)
# source("monocle_helper_functions.R")
library(plyr)
library(stringr)
library(dplyr) 
library(grid)
library(gridExtra)
load("RData/analysis_shalek_data.RData")

#functions used in the scripts
order_shalek_cells_by_original_states <- function(cds_subset, original_cds, root_state, cells_state_2, cells_state_3) {
  cds_subset = estimateSizeFactors(cds_subset)
  cds_subset = estimateDispersions(cds_subset)
  cds_subset = detectGenes(cds_subset)

  closeAllConnections()
  expressed_genes <- row.names(subset(fData(original_cds), num_cells_expressed > 50))
  genes_in_range <- selectGenesInExpressionRange(original_cds[expressed_genes,], 2, Inf, 0.1, stat_fun=function(x) { median(round(x)) })

  differential_genes = differentialGeneTest(cds_subset[genes_in_range, ], fullModelFormulaStr = '~experiment_name + stim_time', cores =detectCores())
  closeAllConnections()
  
  # Filter out high mRNAs
  cds_subset <- cds_subset[, pData(cds_subset)$Total_mRNAs < 75000]
  order_genes = subset(differential_genes, qval < 0.05) %>% top_n(34, -qval)
  order_genes = order_genes$gene_id
  cds_subset = setOrderingFilter(cds_subset, order_genes)

  cds_subset = reduceDimension(cds_subset, use_irlba = F, use_vst=T, method="ICA", scaling=F, pseudo_expr=0) 
  cds_subset = orderCells(cds_subset, num_paths = 2, reverse = T)
  
  #determine the mapping between original state 1/2/3 and new state 1/2/3:  
  overlap_state_1 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 1, ]), root_state))
  overlap_state_2 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 2, ]), root_state))
  overlap_state_3 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 3, ]), root_state))
  
  #find the state corresponding to the original root state
  overlap_vec <- c(overlap_state_1, overlap_state_2, overlap_state_3)
  max_ind <- which(overlap_vec == max(overlap_vec))

  cds_subset = orderCells(cds_subset, max_ind, num_paths = 2, reverse = T)
  overlap_state_2 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 2, ]), cells_state_2))
  overlap_state_3 <- length(intersect(row.names(pData(cds_subset)[pData(cds_subset)$State == 3, ]), cells_state_2))

  #find the new state corresponding to the original state 2: 
  overlap_state_2_vec <- c(overlap_state_2, overlap_state_3)
  state_2_max_ind <- which(overlap_state_2_vec == max(overlap_state_2_vec)) #root state
  
  if(state_2_max_ind != 1) {
    State_vec <- pData(cds_subset)$State
    pData(cds_subset)$State[State_vec == 2] <- 3
    pData(cds_subset)$State[State_vec == 3] <- 2
  }

  cds_subset
}

plot_monocle_spanning_tree_vectorized = function(cds) {
    # monocle::plot_spanning_tree(cds, color_by="interaction(experiment_name, time)", cell_size=2) + 
    #     scale_color_manual(values=shalek_custom_color_scale_plus_states) + nm_theme()
    monocle::plot_spanning_tree(cds, color_by="State", cell_size=2) + nm_theme()
        
}

# Get the states assignment used in the manuscript
root_state <- row.names(subset(pData(Shalek_abs_subset_ko_LPS), State == 1))
cells_state_2 <- row.names(subset(pData(Shalek_abs_subset_ko_LPS), State == 2))
cells_state_3 <- row.names(subset(pData(Shalek_abs_subset_ko_LPS), State == 3))

# Now get the original Shalek KO transcript counts for subsetting
Shalek_abs_subset_ko_LPS <- Shalek_abs[, pData(Shalek_abs)$experiment_name %in% c('Ifnar1_KO_LPS', 'Stat1_KO_LPS',  "LPS", "Unstimulated_Replicate")]
pData(Shalek_abs_subset_ko_LPS)[, 'stim_time'] <- as.character(pData(Shalek_abs_subset_ko_LPS)$time)
pData(Shalek_abs_subset_ko_LPS)$stim_time[pData(Shalek_abs_subset_ko_LPS)$stim_time == ''] <- 0
pData(Shalek_abs_subset_ko_LPS)$stim_time <- as.integer(revalue(pData(Shalek_abs_subset_ko_LPS)$stim_time, c("1h" = 1, "2h" = 2, "4h" = 4, "6h" = 6)))
Shalek_abs_subset_ko_LPS <- detectGenes(Shalek_abs_subset_ko_LPS, min_expr = 0.1)

load('./RData/abs_cds_downsampled_cells_branch_genes_1')
load('./RData/abs_cds_downsampled_cells_branch_genes_2')
load('./RData/abs_cds_downsampled_cells_branch_genes_3')
load('./RData/abs_cds_downsampled_cells_branch_genes_4')
load('./RData/abs_cds_downsampled_cells_branch_genes_5')
load('./RData/abs_cds_downsampled_cells_branch_genes_6')
load('./RData/abs_cds_downsampled_cells_branch_genes_7')
load('./RData/abs_cds_downsampled_cells_branch_genes_8')

cds_downsampled_cells_ordered <- c(lapply(abs_cds_downsampled_cells_branch_genes_1, function(x) x$order_cds), 
                                lapply(abs_cds_downsampled_cells_branch_genes_2, function(x) x$order_cds), 
                                lapply(abs_cds_downsampled_cells_branch_genes_3, function(x) x$order_cds), 
                                lapply(abs_cds_downsampled_cells_branch_genes_4, function(x) x$order_cds),
                                lapply(abs_cds_downsampled_cells_branch_genes_5, function(x) x$order_cds),
                                lapply(abs_cds_downsampled_cells_branch_genes_6, function(x) x$order_cds),
                                lapply(abs_cds_downsampled_cells_branch_genes_7, function(x) x$order_cds),
                                lapply(abs_cds_downsampled_cells_branch_genes_8, function(x) x$order_cds))

cds_downsampled_cells_branch_genes <- c(lapply(abs_cds_downsampled_cells_branch_genes_1, function(x) x$res), 
                                lapply(abs_cds_downsampled_cells_branch_genes_2, function(x) x$res), 
                                lapply(abs_cds_downsampled_cells_branch_genes_3, function(x) x$res), 
                                lapply(abs_cds_downsampled_cells_branch_genes_4, function(x) x$res),
                                lapply(abs_cds_downsampled_cells_branch_genes_5, function(x) x$res),
                                lapply(abs_cds_downsampled_cells_branch_genes_6, function(x) x$res),
                                lapply(abs_cds_downsampled_cells_branch_genes_7, function(x) x$res),
                                lapply(abs_cds_downsampled_cells_branch_genes_8, function(x) x$res))

# Generate downsampled sets of cells
set.seed(5)
MIN_PROPORTION = 0.1
MAX_PROPORTION = 1
STEP = 0.1
REPS_PER = 3
EXTRA_PROPORTIONS = c(0.85, 0.95) 

downsampled_proportions = sort(rep(c(seq(MIN_PROPORTION, MAX_PROPORTION, by=STEP), EXTRA_PROPORTIONS) , REPS_PER))
names(downsampled_proportions) = downsampled_proportions  # will tie CDS objects to proportion for later

## Plot the reordered trajectories
cds_downsampled_cells_ordered_trajectories = lapply(cds_downsampled_cells_ordered, plot_monocle_spanning_tree_vectorized)

#save the tree structure across different depth: 
for(i in 1:length(cds_downsampled_cells_ordered)) {
    p <- monocle::plot_spanning_tree(cds_downsampled_cells_ordered[[i]], color_by="interaction(experiment_name, time)", cell_size=0.5, cell_link_size = 0.01) + 
         scale_color_manual(values=shalek_custom_color_scale_plus_states) + nm_theme() + 
         scale_size(range = c(0.01, 0.5),  limits = c(0.01, 0.5))
    ggsave(p, filename = paste("tmp/", i, downsampled_proportions[i], '.pdf', sep = '_'), height = 1.1, width = 1.1)
}

#get the number of cells: 
lapply(cds_downsampled_cells_ordered[c(3, 6, 9, 12, 15, 17, 21, 24, 27, 30, 33, 36)], ncol)
lapply(cds_downsampled_cells_ordered[c(1, 5, 8, 12, 15, 16, 20, 22, 27, 30, 33, 36)], ncol)

# # Get branching genes for all these subsets
# #parallel on multiple clusters (think about MPI, etc): 
# cds_downsampled_cells_branch_genes = lapply(cds_downsampled_cells_ordered, function(cds) { 
#   branchTest(cds[, ], fullModelFormulaStr = '~sm.ns(Pseudotime, df = 3)*Lineage', cores = detectCores(), relative_expr = T, weighted = T) })

# save(cds_downsampled_cells_ordered, cds_downsampled_cells_branch_genes, Shalek_abs_subset_ko_LPS, file="analysis_cell_downsampling.gz")

#BEAM gene overlap, FPR, p-value spearman correlation (with number_gene_removed)
cds_downsampled_cells_significant_genes= lapply(cds_downsampled_cells_branch_genes, function(branching_genes) row.names(subset(branching_genes, qval <= 0.01 )))
original_branch_genes = cds_downsampled_cells_significant_genes[[36]] 

statistics_per_depth <- lapply(cds_downsampled_cells_branch_genes, function(branching_genes) {
  true_positive <- intersect(row.names(branching_genes[branching_genes$qval < 0.01, ]), 
                            original_branch_genes)
  data.frame(true_positive = length(true_positive), 
    precision = length(true_positive) / nrow(branching_genes[branching_genes$qval < 0.01, ]))

}
)

recall <- lapply(cds_downsampled_cells_branch_genes, function(branching_genes) {
  true_positive <- intersect(row.names(branching_genes[branching_genes$qval < 0.01, ]), 
                            original_branch_genes)
 length(true_positive) / length(original_branch_genes)
}
)

original_branch_genes_res = cds_downsampled_cells_branch_genes[[36]] 

spearman_correlation = lapply(cds_downsampled_cells_branch_genes, function(branching_genes) {
    cor(original_branch_genes_res[row.names(branching_genes), 'pval'], branching_genes[, 'pval'], method = 'spearman') 
  }
  )

statistics_per_depth <- do.call(rbind.data.frame, statistics_per_depth)
statistics_per_depth$proportion_original_cells <- downsampled_proportions
statistics_per_depth$recall = unlist(recall)
statistics_per_depth$spearman_correlation <- unlist(spearman_correlation)

#number of cells assigned correct states: 
load('./RData/Shalek_abs_subset_ko_LPS')
state_correct_fraction <- lapply(cds_downsampled_cells_ordered, function(cds) {
  state_correct_num <- sum(as.numeric(pData(cds)$State) - as.numeric(pData(Shalek_abs_subset_ko_LPS)[colnames(cds), 'State']) == 0) #cds_downsampled_cells_ordered[[36]]
  state_correct_num / ncol(cds)
}
)

statistics_per_depth$state_correct_fraction <- unlist(state_correct_fraction)
pdf("./supplementary_figures/fig10c_state_corr_fraction.pdf", height=1.6, width=2.2)
ggplot(statistics_per_depth, aes(proportion_original_cells, state_correct_fraction)) + ylim(0, 1) +
    geom_point(color = 'red', size = 0.75) + xlab('Proportion of original cells') + ylab('Fraction of correctly \n assigned states') + 
    nm_theme()
dev.off()

pseudotime_correlation <- lapply(cds_downsampled_cells_ordered, function(cds) {
  state_correlation <- cor(as.numeric(pData(cds)$Pseudotime), as.numeric(pData(Shalek_abs_subset_ko_LPS)[colnames(cds), 'Pseudotime'])) ##cds_downsampled_cells_ordered[[36]]
}
)

statistics_per_depth$pseudotime_correlation <- unlist(pseudotime_correlation)
pdf("./supplementary_figures/fig10c_pseudotime_cor.pdf", height=1.6, width=2.2)
ggplot(statistics_per_depth, aes(proportion_original_cells, pseudotime_correlation)) + ylim(0, 1) +
    geom_point(color = 'red', size = 0.75) + xlab('Proportion of original cells') + ylab('Pseudotime correlation') + 
    nm_theme()
dev.off()

#BEAM gene overlap, FPR, p-value spearman correlation without tree re-ordering: 
#######just sample the cells: 
# Generate downsampled sets of cells

################################################Parallel the above analysis in three individual runs#####################################################
load('./RData/PILOT_SAMPLING_branch_genes_A.RData')
load('./RData/PILOT_SAMPLING_branch_genes_B.RData')
load('./RData/PILOT_SAMPLING_branch_genes_C.RData')
################################################Parallel the above analysis in three individual runs#####################################################
#calculate the result: 
PILOT_SAMPLING_branch_genes <- Reduce(append, list(PILOT_SAMPLING_branch_genes_A, PILOT_SAMPLING_branch_genes_B, PILOT_SAMPLING_branch_genes_C))

#calculate the pilot experiment: 
statistics_per_depth_fix_tree <- lapply(PILOT_SAMPLING_branch_genes, function(branching_genes) {
  true_positive <- intersect(row.names(branching_genes[branching_genes$qval < 0.01, ]), 
                            row.names(ko_branching_genes[ko_branching_genes$qval < 0.05, ]))
  precision = length(true_positive) / nrow(branching_genes[branching_genes$qval < 0.01, ])
  recall = length(true_positive) / nrow(ko_branching_genes[ko_branching_genes$qval < 0.01, ])
  data.frame(true_positive = length(true_positive), precision = precision, recall = recall)
}
)

spearman_correlation = lapply(PILOT_SAMPLING_branch_genes, function(branching_genes) {
    cor(ko_branching_genes[row.names(branching_genes), 'pval'], branching_genes[, 'pval'], method = 'spearman') 
  }
  )


PILOT_SAMPLING_A = rep(c(0.99, 0.98, 0.96, 0.94, 0.9), each = 3)
PILOT_SAMPLING_B = rep(c(0.85, 0.8, 0.75, 0.7, 0.6), each = 3)
PILOT_SAMPLING_C = rep(c(0.5, 0.4, 0.3, 0.2, 0.1), each = 3)

statistics_per_depth_fix_tree <- do.call(rbind.data.frame, statistics_per_depth_fix_tree)
statistics_per_depth_fix_tree$proportion_original_cells = c(PILOT_SAMPLING_A, PILOT_SAMPLING_B, PILOT_SAMPLING_C)
statistics_per_depth_fix_tree$spearman_correlation <- unlist(spearman_correlation)

statistics_per_depth_all <- rbind(statistics_per_depth_fix_tree[, c('proportion_original_cells', 'precision', 'recall', 'spearman_correlation')], 
      statistics_per_depth[, c('proportion_original_cells', 'precision', 'recall', 'spearman_correlation')])

statistics_per_depth_all$Type <- c(rep('Cell ordering fixed', nrow(statistics_per_depth_fix_tree)), 
                                  rep('Cell reordering', nrow(statistics_per_depth)))

subset_statistics_per_depth_all <- subset(statistics_per_depth_all, proportion_original_cells != 1)

pdf("./supplementary_figures/fig10c_si_recall_test.pdf", height=1.5, width=1.5)
ggplot(subset_statistics_per_depth_all, aes(proportion_original_cells, recall)) +
    geom_point(aes(color = Type), size = 0.75) + xlab('Proportion of original cells') + ylab('recall') + xlim(0, 1) + ylim(0, 1) + 
    nm_theme() 
dev.off()

pdf("./supplementary_figures/fig10c_si_recall_helper_test.pdf", height=1.5, width=1.5)
ggplot(subset_statistics_per_depth_all, aes(proportion_original_cells, recall)) +
    geom_point(aes(color = Type), size = 0.75) + xlab('Proportion of original cells') + ylab('recall') + xlim(0, 1) + ylim(0, 1) 
dev.off()

pdf("./supplementary_figures/fig10c_si_precision_test.pdf", height=1.5, width=1.5)
ggplot(subset_statistics_per_depth_all, aes(proportion_original_cells, precision)) +
    geom_point(aes(color = Type), size = 0.75) + xlab('Proportion of original cells') + ylab('Precision') + xlim(0, 1) + ylim(0, 1) +
    nm_theme()
dev.off()

pdf("./supplementary_figures/fig10c_si_spearman_test.pdf", height=1.5, width=1.5)
ggplot(subset_statistics_per_depth_all, aes(proportion_original_cells, spearman_correlation)) +
  geom_point(aes(color = Type), size = 0.75) + xlab('Proportion of original cells') +  ylab('Spearman rank correlation \n (P value)') + ylim(0, 1) +
  nm_theme()
dev.off()

abs_statistics_per_depth_all <- statistics_per_depth_all
save(file = 'abs_statistics_per_depth_all', abs_statistics_per_depth_all)

save.image('./RData/analysis_cell_downsampling.RData')
