######################################################################
# Downsampling the number of cells to benchmark BEAM performance 
# ####################################################################
library(monocle)
# library(devtools)
# load_all('~/Projects/monocle-dev')
library(xacHelper)
# source("monocle_helper_functions.R")
library(dplyr) 
library(grid)
library(gridExtra)

# load('./RData/analysis_cell_downsampling.RData')
load("RData/analysis_shalek_data.RData")

#functions used in the scripts
branchTest_downsampling <- function(cds_subset, cell_id_list, relative_expr = T, cores = detectCores()) {
  branching_genes = branchTest(cds_subset, fullModelFormulaStr =  '~sm.ns(Pseudotime, df = 3)*Lineage', 
    cores = cores, relative_expr = relative_expr, weighted = T,  cell_id_list = cell_id_list)

  branching_genes
}

set.seed(5)
MIN_PROPORTION = 0.1
MAX_PROPORTION = 1
STEP = 0.1 #0.1
REPS_PER = 3
EXTRA_PROPORTIONS = c(0.85, 0.95) 

downsampled_proportions = sort(rep(c(seq(MIN_PROPORTION, MAX_PROPORTION, by=STEP), EXTRA_PROPORTIONS) , REPS_PER))
names(downsampled_proportions) = downsampled_proportions  # will tie CDS objects to proportion for later

# cds_downsampled_cells_branch_genes = lapply(downsampled_proportions, function(x) { 
#   cell_id_list <- sample(ncol(Shalek_abs_subset_ko_LPS), round(ncol(Shalek_abs_subset_ko_LPS) * x))
#   res <- branchTest_downsampling(Shalek_abs_subset_ko_LPS[1:10, ], cell_id_list, cores = 1)

#   res
#      })

#do the same thing for the FPKM values: 
Shalek_std_subset_KO_LPS <- Shalek_std[, colnames(Shalek_abs_subset_ko_LPS)]

#assign pseudotime associated data calculated with FPKM values to all datasets: USE THE ORIGINAL CELL ORDERING
pData(Shalek_std_subset_KO_LPS) <- pData(Shalek_abs_subset_ko_LPS[, colnames(Shalek_std_subset_KO_LPS)])

#pass cell ordering information as well: 
Shalek_std_subset_KO_LPS@auxOrderingData <- Shalek_abs_subset_ko_LPS@auxOrderingData
Shalek_std_subset_KO_LPS@dim_reduce_type <- 'ICA'

#calculate the branch genes for knockout experiment, using fpkm values
std_ko_branching_genes <- branchTest(Shalek_std_subset_KO_LPS[, ], fullModelFormulaStr = full_model_string, cores = 1, relative_expr = F, weighted = T) #detectCores()
closeAllConnections()

closeAllConnections()

#distributed into eight clusters: 
#  [1] "0.1"  "0.1"  "0.1"  "0.2"  "0.2"  "0.2"  | "0.3"  "0.3"  "0.3"  "0.4"
# [11] "0.4"  "0.4" |  "0.5"  "0.5"  "0.5"  "| 0.6"  "0.6"  "0.6"  "0.7"  "0.7"
# [21] "0.7"  "0.8"  "0.8"  "0.8"  "0.85" "0.85" "0.85" "0.9"  "0.9"  "0.9"
# [31] "0.95" "0.95" "0.95" "1"    "1"    "1"

# ######################################################################################################
# std_cds_downsampled_cells_branch_genes = lapply(downsampled_proportions[c(1, 15, 31)], function(x) { 
#   cell_id_list <- sample(ncol(Shalek_abs_subset_ko_LPS), round(ncol(Shalek_std_subset_KO_LPS) * x))
#   res <- branchTest_downsampling(Shalek_std_subset_KO_LPS[, ], cell_id_list, relative_expr = F, cores = 1)

#   return(list(cell_id_list = cell_id_list, res = res))
#    })
# std_cds_downsampled_cells_branch_genes_test <- std_cds_downsampled_cells_branch_genes
# save(std_cds_downsampled_cells_branch_genes_test, file = 'std_cds_downsampled_cells_branch_genes_test')
######################################################################################################

std_cds_downsampled_cells_branch_genes = lapply(downsampled_proportions[1:6], function(x) { 
  cell_id_list <- sample(ncol(Shalek_abs_subset_ko_LPS), round(ncol(Shalek_std_subset_KO_LPS) * x))
  res <- branchTest_downsampling(Shalek_std_subset_KO_LPS[, ], cell_id_list, relative_expr = F, cores = 1)

  return(list(cell_id_list = cell_id_list, res = res))
   })
std_cds_downsampled_cells_branch_genes_1 <- std_cds_downsampled_cells_branch_genes
save(std_cds_downsampled_cells_branch_genes_1, file = 'std_cds_downsampled_cells_branch_genes_1')

std_cds_downsampled_cells_branch_genes = lapply(downsampled_proportions[7:12], function(x) { 
  cell_id_list <- sample(ncol(Shalek_abs_subset_ko_LPS), round(ncol(Shalek_std_subset_KO_LPS) * x))
  res <- branchTest_downsampling(Shalek_std_subset_KO_LPS[, ], cell_id_list, relative_expr = F, cores = 1)

  return(list(cell_id_list = cell_id_list, res = res))
   })
std_cds_downsampled_cells_branch_genes_2 <- std_cds_downsampled_cells_branch_genes
save(std_cds_downsampled_cells_branch_genes_2, file = 'std_cds_downsampled_cells_branch_genes_2')

std_cds_downsampled_cells_branch_genes = lapply(downsampled_proportions[13:19], function(x) { 
  cell_id_list <- sample(ncol(Shalek_abs_subset_ko_LPS), round(ncol(Shalek_std_subset_KO_LPS) * x))
  res <- branchTest_downsampling(Shalek_std_subset_KO_LPS[, ], cell_id_list, relative_expr = F, cores = 1)

  return(list(cell_id_list = cell_id_list, res = res))
   })
std_cds_downsampled_cells_branch_genes_3 <- std_cds_downsampled_cells_branch_genes
save(std_cds_downsampled_cells_branch_genes_3, file = 'std_cds_downsampled_cells_branch_genes_3')

std_cds_downsampled_cells_branch_genes = lapply(downsampled_proportions[20:21], function(x) { 
  cell_id_list <- sample(ncol(Shalek_abs_subset_ko_LPS), round(ncol(Shalek_std_subset_KO_LPS) * x))
  res <- branchTest_downsampling(Shalek_std_subset_KO_LPS[, ], cell_id_list, relative_expr = F, cores = 1)

  return(list(cell_id_list = cell_id_list, res = res))
   })
std_cds_downsampled_cells_branch_genes_4 <- std_cds_downsampled_cells_branch_genes
save(std_cds_downsampled_cells_branch_genes_4, file = 'std_cds_downsampled_cells_branch_genes_4')

std_cds_downsampled_cells_branch_genes = lapply(downsampled_proportions[22:24], function(x) { 
  cell_id_list <- sample(ncol(Shalek_abs_subset_ko_LPS), round(ncol(Shalek_std_subset_KO_LPS) * x))
  res <- branchTest_downsampling(Shalek_std_subset_KO_LPS[, ], cell_id_list, relative_expr = F, cores = 1)

  return(list(cell_id_list = cell_id_list, res = res))
   })
std_cds_downsampled_cells_branch_genes_5 <- std_cds_downsampled_cells_branch_genes
save(std_cds_downsampled_cells_branch_genes_5, file = 'std_cds_downsampled_cells_branch_genes_5')

std_cds_downsampled_cells_branch_genes = lapply(downsampled_proportions[25:28], function(x) { 
  cell_id_list <- sample(ncol(Shalek_abs_subset_ko_LPS), round(ncol(Shalek_std_subset_KO_LPS) * x))
  res <- branchTest_downsampling(Shalek_std_subset_KO_LPS[, ], cell_id_list, relative_expr = F, cores = 1)

  return(list(cell_id_list = cell_id_list, res = res))
   })
std_cds_downsampled_cells_branch_genes_6 <- std_cds_downsampled_cells_branch_genes
save(std_cds_downsampled_cells_branch_genes_6, file = 'std_cds_downsampled_cells_branch_genes_6')

std_cds_downsampled_cells_branch_genes = lapply(downsampled_proportions[29:32], function(x) { 
  cell_id_list <- sample(ncol(Shalek_abs_subset_ko_LPS), round(ncol(Shalek_std_subset_KO_LPS) * x))
  res <- branchTest_downsampling(Shalek_std_subset_KO_LPS[, ], cell_id_list, relative_expr = F, cores = 1)

  return(list(cell_id_list = cell_id_list, res = res))
   })
std_cds_downsampled_cells_branch_genes_7 <- std_cds_downsampled_cells_branch_genes
save(std_cds_downsampled_cells_branch_genes_7, file = 'std_cds_downsampled_cells_branch_genes_7')

std_cds_downsampled_cells_branch_genes = lapply(downsampled_proportions[33:36], function(x) { 
  cell_id_list <- sample(ncol(Shalek_abs_subset_ko_LPS), round(ncol(Shalek_std_subset_KO_LPS) * x))
  res <- branchTest_downsampling(Shalek_std_subset_KO_LPS[, ], cell_id_list, relative_expr = F, cores = 1)

  return(list(cell_id_list = cell_id_list, res = res))
   })
std_cds_downsampled_cells_branch_genes_8 <- std_cds_downsampled_cells_branch_genes
save(std_cds_downsampled_cells_branch_genes_8, file = 'std_cds_downsampled_cells_branch_genes_8')

load('std_cds_downsampled_cells_branch_genes_1')
load('std_cds_downsampled_cells_branch_genes_2')
load('std_cds_downsampled_cells_branch_genes_3')
load('std_cds_downsampled_cells_branch_genes_4')
load('std_cds_downsampled_cells_branch_genes_5')
load('std_cds_downsampled_cells_branch_genes_6')
load('std_cds_downsampled_cells_branch_genes_7')
load('std_cds_downsampled_cells_branch_genes_8')

std_cds_downsampled_cells_branch_genes <- c(std_cds_downsampled_cells_branch_genes_1, std_cds_downsampled_cells_branch_genes_2, 
                                std_cds_downsampled_cells_branch_genes_3, std_cds_downsampled_cells_branch_genes_4, 
                                std_cds_downsampled_cells_branch_genes_5, std_cds_downsampled_cells_branch_genes_6, 
                                std_cds_downsampled_cells_branch_genes_7, std_cds_downsampled_cells_branch_genes_8)

std_cds_downsampled_cells_ordered <- c(lapply(std_cds_downsampled_cells_branch_genes_1, function(x) x$order_cds), 
                                lapply(std_cds_downsampled_cells_branch_genes_2, function(x) x$order_cds), 
                                lapply(std_cds_downsampled_cells_branch_genes_3, function(x) x$order_cds), 
                                lapply(std_cds_downsampled_cells_branch_genes_4, function(x) x$order_cds),
                                lapply(std_cds_downsampled_cells_branch_genes_5, function(x) x$order_cds),
                                lapply(std_cds_downsampled_cells_branch_genes_6, function(x) x$order_cds),
                                lapply(std_cds_downsampled_cells_branch_genes_7, function(x) x$order_cds),
                                lapply(std_cds_downsampled_cells_branch_genes_8, function(x) x$order_cds))

std_cds_downsampled_cells_branch_genes <- c(lapply(std_cds_downsampled_cells_branch_genes_1, function(x) x$res), 
                                lapply(std_cds_downsampled_cells_branch_genes_2, function(x) x$res), 
                                lapply(std_cds_downsampled_cells_branch_genes_3, function(x) x$res), 
                                lapply(std_cds_downsampled_cells_branch_genes_4, function(x) x$res),
                                lapply(std_cds_downsampled_cells_branch_genes_5, function(x) x$res),
                                lapply(std_cds_downsampled_cells_branch_genes_6, function(x) x$res),
                                lapply(std_cds_downsampled_cells_branch_genes_7, function(x) x$res),
                                lapply(std_cds_downsampled_cells_branch_genes_8, function(x) x$res))

# Perform the BEAM test on those genes: 
cds_downsampled_cells_branch_genes <- std_cds_downsampled_cells_branch_genes
PILOT_SAMPLING_branch_genes <- lapply(cds_downsampled_cells_branch_genes, function(x) x$res) 
# std_PILOT_SAMPLING_branch_genes <- lapply(std_cds_downsampled_cells_branch_genes, function(x) x$res) #std_cds_downsampled_cells_branch_genes
save(cds_downsampled_cells_branch_genes, file="cds_downsampled_cells_branch_genes.gz")

#calculate the pilot experiment: 
cal_cell_downsampling_stats <- function(PILOT_SAMPLING_branch_genes, ori_branching_genes = ko_branching_genes, proportion_original_cells = as.numeric(names(downsampled_proportions))) {
  stat <- lapply(PILOT_SAMPLING_branch_genes, function(branching_genes) {
    true_positive <- intersect(row.names(branching_genes[branching_genes$qval < 0.01, ]), 
                              row.names(ori_branching_genes[ori_branching_genes$qval < 0.05, ]))
    precision = length(true_positive) / nrow(branching_genes[branching_genes$qval < 0.01, ])
    recall = length(true_positive) / nrow(ori_branching_genes[ori_branching_genes$qval < 0.01, ])
    data.frame(true_positive = length(true_positive), precision = precision, recall = recall)
  }
  )

  spearman_correlation <- lapply(PILOT_SAMPLING_branch_genes, function(branching_genes) {
      cor(ori_branching_genes[row.names(branching_genes), 'pval'], branching_genes[, 'pval'], method = 'spearman') 
    }
    )
  
  stat <- do.call(rbind.data.frame, stat)
  stat$proportion_original_cells = proportion_original_cells
  stat$spearman_correlation <- unlist(spearman_correlation)

  return(stat)
}

# abs_statistics_per_depth <- cal_cell_downsampling_stats(PILOT_SAMPLING_branch_genes)
abs_statistics_per_depth <- abs_statistics_per_depth_all

std_statistics_per_depth <- cal_cell_downsampling_stats(std_cds_downsampled_cells_branch_genes, std_ko_branching_genes)
std_statistics_per_depth$Type <- rep('FPKM', nrow(std_statistics_per_depth))

statistics_per_depth_all <- rbind(abs_statistics_per_depth[, c('proportion_original_cells', 'precision', 'recall', 'spearman_correlation', 'Type')], 
      std_statistics_per_depth[, c('proportion_original_cells', 'precision', 'recall', 'spearman_correlation', 'Type')])

# statistics_per_depth_all$Type <- c(rep('Transcript counts', nrow(abs_statistics_per_depth)), 
#                                   rep('FPKM', nrow(std_statistics_per_depth)))

subset_statistics_per_depth_all <- subset(statistics_per_depth_all, proportion_original_cells != 1)

pdf("./supplementary_figures/beam_fpkm_abs_recall.pdf", height=1.5, width=1.5)
ggplot(subset_statistics_per_depth_all, aes(proportion_original_cells, recall)) +
    geom_point(aes(color = Type), size = 0.75) + xlab('Proportion of original cells') + ylab('recall') + xlim(0, 1) + ylim(0, 1) + nm_theme() 
dev.off()

pdf("./supplementary_figures/beam_fpkm_abs_recall_helper.pdf", height=1.5, width=1.5)
ggplot(subset_statistics_per_depth_all, aes(proportion_original_cells, recall)) +
    geom_point(aes(color = Type), size = 0.75) + xlab('Proportion of original cells') + ylab('recall') + xlim(0, 1) + ylim(0, 1) 
dev.off()

pdf("./supplementary_figures/beam_fpkm_abs_precision.pdf", height=1.5, width=1.5)
ggplot(subset_statistics_per_depth_all, aes(proportion_original_cells, precision)) +
    geom_point(aes(color = Type), size = 0.75) + xlab('Proportion of original cells') + ylab('Precision') + xlim(0, 1) + ylim(0, 1) +
    nm_theme()
dev.off()

pdf("./supplementary_figures/beam_fpkm_abs_spearman.pdf", height=1.5, width=1.5)
ggplot(subset_statistics_per_depth_all, aes(proportion_original_cells, spearman_correlation)) +
  geom_point(aes(color = Type), size = 0.75) + xlab('Proportion of original cells') +  ylab('Spearman rank correlation \n (P value)') + ylim(0, 1) +
  nm_theme()
dev.off()

save.image('./RData/analysis_BEAM_cell_downsampling.RData')
