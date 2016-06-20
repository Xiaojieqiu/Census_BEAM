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

#################################################Parallel the above analysis#####################################################

std_cds_downsampled_cells_branch_genes = lapply(downsampled_proportions[29:32], function(x) { 
  cell_id_list <- sample(ncol(Shalek_abs_subset_ko_LPS), round(ncol(Shalek_std_subset_KO_LPS) * x))
  res <- branchTest_downsampling(Shalek_std_subset_KO_LPS[, ], cell_id_list, relative_expr = F, cores = 1)

  return(list(cell_id_list = cell_id_list, res = res))
   })
std_cds_downsampled_cells_branch_genes_7 <- std_cds_downsampled_cells_branch_genes
save(std_cds_downsampled_cells_branch_genes_7, file = 'std_cds_downsampled_cells_branch_genes_7')

