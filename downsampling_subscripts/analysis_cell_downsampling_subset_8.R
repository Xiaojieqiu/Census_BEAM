########################################
# Downsampling the number of cells
# ########################################
# library(monocle)
library(devtools)
load_all('~/Projects/monocle-dev')
library(xacHelper)
# source("monocle_helper_functions.R")
library(plyr)
library(stringr)
library(dplyr) 
library(grid)
library(gridExtra)
load("RData/analysis_shalek_data.RData")

#functions used in the scripts
order_shalek_cells_by_original_states <- function(cds_subset, original_cds, root_state, cells_state_2, cells_state_3, genes_in_range_vec = genes_in_range) {
  cds_subset = estimateSizeFactors(cds_subset)
  cds_subset = estimateDispersions(cds_subset)
  cds_subset = detectGenes(cds_subset)

  closeAllConnections()
  # expressed_genes <- row.names(subset(fData(original_cds), num_cells_expressed > 50))
  # genes_in_range <- selectGenesInExpressionRange(original_cds[expressed_genes,], 2, Inf, 0.1, stat_fun=function(x) { median(round(x)) })

  differential_genes = differentialGeneTest(cds_subset[genes_in_range_vec, ], fullModelFormulaStr = '~experiment_name + stim_time', cores =detectCores())
  closeAllConnections()
  
  # Filter out high mRNAs
  # cds_subset <- cds_subset[, pData(cds_subset)$Total_mRNAs < 170000]  
  top_gene_num <- sum(fData(cds_subset)$use_for_ordering)
  print(paste('top_gen_num is', top_gene_num))
  order_genes = subset(differential_genes, qval < 0.05) %>% top_n(top_gene_num, -qval)
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

# # Now get the original Shalek KO transcript counts for subsetting
# Shalek_abs_subset_ko_LPS <- Shalek_abs[, pData(Shalek_abs)$experiment_name %in% c('Ifnar1_KO_LPS', 'Stat1_KO_LPS',  "LPS", "Unstimulated_Replicate")]
# pData(Shalek_abs_subset_ko_LPS)[, 'stim_time'] <- as.character(pData(Shalek_abs_subset_ko_LPS)$time)
# pData(Shalek_abs_subset_ko_LPS)$stim_time[pData(Shalek_abs_subset_ko_LPS)$stim_time == ''] <- 0
# pData(Shalek_abs_subset_ko_LPS)$stim_time <- as.integer(revalue(pData(Shalek_abs_subset_ko_LPS)$stim_time, c("1h" = 1, "2h" = 2, "4h" = 4, "6h" = 6)))
# Shalek_abs_subset_ko_LPS <- detectGenes(Shalek_abs_subset_ko_LPS, min_expr = 0.1)

print(paste('number of cells for the ko data is: ', ncol(Shalek_abs_subset_ko_LPS)))
set.seed(2016)

# Generate downsampled sets of cells
set.seed(5)
MIN_PROPORTION = 0.1
MAX_PROPORTION = 1
STEP = 0.1
REPS_PER = 3
EXTRA_PROPORTIONS = c(0.85, 0.95) 

downsampled_proportions = sort(rep(c(seq(MIN_PROPORTION, MAX_PROPORTION, by=STEP), EXTRA_PROPORTIONS) , REPS_PER))
names(downsampled_proportions) = downsampled_proportions  # will tie CDS objects to proportion for later

cds_downsampled_cells = lapply(downsampled_proportions, function(x) { Shalek_abs_subset_ko_LPS[, sample(ncol(Shalek_abs_subset_ko_LPS), round(ncol(Shalek_abs_subset_ko_LPS) * x))] })

# # Reorder those subsets
# cds_downsampled_cells_ordered = lapply(cds_downsampled_cells, order_shalek_cells_by_original_states, Shalek_abs_subset_ko_LPS, root_state, cells_state_2, cells_state_3)
# save(cds_downsampled_cells_ordered, Shalek_abs_subset_ko_LPS, file="analysis_cell_downsampling.gz")

#################################################Parallel the above analysis#####################################################
abs_cds_downsampled_cells_branch_genes = lapply(cds_downsampled_cells[33:36], function(cds) { 
  cds <- order_shalek_cells_by_original_states(cds, Shalek_abs_subset_ko_LPS, root_state, cells_state_2, cells_state_3)
  res <- branchTest(cds[, ], fullModelFormulaStr = '~sm.ns(Pseudotime, df = 3)*Lineage', cores = detectCores(), relative_expr = T, weighted = T)
  return(list(order_cds = cds, res = res))
   })
abs_cds_downsampled_cells_branch_genes_8 <- abs_cds_downsampled_cells_branch_genes
save(abs_cds_downsampled_cells_branch_genes_8, file = 'abs_cds_downsampled_cells_branch_genes_8')
#################################################Parallel the above analysis#####################################################

