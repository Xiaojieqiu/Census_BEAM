################################################################################################################################################
 # The following script is based on Adnrew's work #
 ################################################################################################################################################
#color scheme for shalek tree analysis: 
shalek_custom_color_scale = c("Unstimulated."="#eff3ff", "On_Chip_Unstim."="#eff3ff", "LPS.1h"="#bdd7e7", "LPS.2h"="#6baed6", "LPS.4h"="#3182bd", "LPS.6h"="#08519c", "PIC.1h"="#bae4b3", "PIC.2h"="#74c476", "PIC.4h"="#31a354", "PIC.6h"="#006d2c", "PAM.1h"="#fdbe85", "PAM.2h"="#fd8d3c", "PAM.4h"="#e6550d", "PAM.6h"="#a63603", "LPS_GolgiPlug.4h_0h"="#fcae91", "LPS_GolgiPlug.4h_1h"="#fb6a4a", "LPS_GolgiPlug.4h_2h"="#cb181d", "Ifnar1_KO_LPS.2h"="#df65b0", "Ifnar1_KO_LPS.4h"="#dd1c77", "Stat1_KO_LPS.2h"="#969696", "Stat1_KO_LPS.4h"="#252525")
shalek_custom_color_scale_plus_states= c(shalek_custom_color_scale, c('1'='#40A43A', '2'='#CB1B1E', '3'='#3660A5', 'Unstimulated_Replicate.' = 'gray'))

library(stringr)
library(plyr)
library(xacHelper)
library(igraph)
library(pheatmap)
library(RColorBrewer)

#use monocle2: 
library(devtools)
load_all('~/Projects/monocle-dev')

#   #load all the go/reactome/kegg datasets for the analysis: 
root_directory <- "./data/Quake_data"

human_go_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Human/symbol/Human_GO_AllPathways_with_GO_iea_June_20_2014_symbol.gmt", sep=""), encoding="latin1")
names(human_go_gsc$gsc) <- str_split_fixed(names(human_go_gsc$gsc), "%", 2)[,1]

human_reactome_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Human/symbol/Pathways/Human_Reactome_June_20_2014_symbol.gmt", sep=""), encoding="latin1")
names(human_reactome_gsc$gsc) <- str_split_fixed(names(human_reactome_gsc$gsc), "%", 2)[,1]

mouse_go_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/GO/MOUSE_GO_bp_with_GO_iea_symbol.gmt", sep=""), encoding="latin1")
names(mouse_go_gsc$gsc) <- str_split_fixed(names(mouse_go_gsc$gsc), "%", 2)[,1]

mouse_reactome_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/Pathways/Mouse_Reactome_June_20_2014_symbol.gmt", sep=""), encoding="latin1")
names(mouse_reactome_gsc$gsc) <- str_split_fixed(names(mouse_reactome_gsc$gsc), "%", 2)[,1]

mouse_kegg_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/Pathways/Mouse_Human_KEGG_June_20_2014_symbol.gmt", sep=""), encoding="latin1")
names(mouse_kegg_gsc$gsc) <- str_split_fixed(names(mouse_kegg_gsc$gsc), "%", 2)[,1]

mouse_go_gsc_cc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/GO/MOUSE_GO_cc_with_GO_iea_symbol.gmt", sep=""), encoding="latin1")
names(mouse_go_gsc_cc$gsc) <- str_split_fixed(names(mouse_go_gsc_cc$gsc), "%", 2)[,1]
mouse_go_gsc_mf <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/GO/MOUSE_GO_mf_with_GO_iea_symbol.gmt", sep=""), encoding="latin1")
names(mouse_go_gsc_mf$gsc) <- str_split_fixed(names(mouse_go_gsc_mf$gsc), "%", 2)[,1]

mouse_reactome_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/Pathways/Mouse_Reactome_June_20_2014_symbol.gmt", sep=""), encoding="latin1")
names(mouse_reactome_gsc$gsc) <- str_split_fixed(names(mouse_reactome_gsc$gsc), "%", 2)[,1]

mouse_kegg_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/Pathways/Mouse_Human_KEGG_June_20_2014_symbol.gmt", sep=""), encoding="latin1")
names(mouse_kegg_gsc$gsc) <- str_split_fixed(names(mouse_kegg_gsc$gsc), "%", 2)[,1]

#set the directory: 
prog_cell_state = "#979797"
AT1_cell_state = "#F05662" 
AT2_cell_state = "#7990C8" 
AT1_Lineage = "#BD1C7C"
AT2_Lineage = "#337DB9" 

#  #load the data: 
Shalek_valid_genes <- read.table('./data/Aviv_data/valid_genes_for_analyis.txt', header = T)
Shalek_exprs_mat <- read.table("./data/Aviv_data/cuffnorm_output_files/genes.fpkm_table", row.names = 1, header = T)
Shalek_fd <- read.table("./data/Aviv_data/cuffnorm_output_files/genes.attr_table", row.names = 1, header = T)
Shalek_pd <- read.table("./data/Aviv_data/sample_metadata_table.txt", sep = '\t', row.names = 1, header = T)
rownames(Shalek_pd) <- paste(rownames(Shalek_pd), "_0", sep = "")
Shalek_exprs_mat <- Shalek_exprs_mat[row.names(Shalek_fd), row.names(Shalek_pd)]
Shalek_std <- newCellDataSet(as.matrix(Shalek_exprs_mat), 
                         phenoData = new("AnnotatedDataFrame", data = Shalek_pd), 
                         featureData = new("AnnotatedDataFrame", data = Shalek_fd), 
                         expressionFamily=tobit(), 
                         lowerDetectionLimit=1)
Shalek_std <- Shalek_std[, which(pData(Shalek_std)$used_in_study == T)]
Shalek_std <- Shalek_std[row.names(Shalek_std) %in% Shalek_valid_genes$gene_id, ] #27386 * 1787 cells
#check the consistency with the current Shalek_abs data: 

Shalek_isoform_fpkm_matrix <- read.table("./data/Aviv_data/cuffnorm_output_files/isoforms.fpkm_table", row.names = 1, header = T)
Shalek_isoform_fpkm_matrix <- Shalek_isoform_fpkm_matrix[, colnames(Shalek_std)]

#load read counts data for recovery algorithm: 
dir = "./data/Aviv_data/cuffnorm_output_files"
Shalek_sample_table <- read.delim(paste(dir, "/samples.table", sep = ''))
Shalek_norm_count <- read.delim(paste(dir, "/genes.count_table", sep = ''))
row.names(Shalek_norm_count) <- Shalek_norm_count$tracking_id
Shalek_norm_count <- Shalek_norm_count[, -1]

Shalek_read_countdata <- round(t(t(Shalek_norm_count) * Shalek_sample_table$internal_scale)) #convert back to the raw counts 
Shalek_read_countdata <- Shalek_read_countdata[row.names(Shalek_std), colnames(Shalek_std)]
colnames(Shalek_read_countdata) <- colnames(Shalek_std)

Shalek_read_countdata_cds <- newCellDataSet(as.matrix(Shalek_read_countdata),
                                            phenoData = new("AnnotatedDataFrame", data = pData(Shalek_std)),
                                            featureData = new("AnnotatedDataFrame", data = fData(Shalek_std)),
                                            expressionFamily = negbinomial(),
                                            lowerDetectionLimit = 1)

pData(Shalek_read_countdata_cds)$Total_mRNAs <- esApply(Shalek_read_countdata_cds, 2, sum)
pData(Shalek_read_countdata_cds)$endogenous_RNA <- esApply(Shalek_read_countdata_cds, 2, sum)

# Convert expression measurements from FPKM to absolute transcript counts, using the isoforms object to estimate the t parameter
Shalek_abs= relative2abs(Shalek_std, estimate_t(Shalek_std), modelFormulaStr = "~1", cores=detectCores(), reads_per_cell = pData(Shalek_read_countdata_cds[, colnames(Shalek_std)])$Total_mRNAs, verbose = T) #
# Shalek_abs_37500= relative2abs(Shalek_std, estimate_t(Shalek_std), modelFormulaStr = "~1", cores=detectCores(), expected_total_mRNAs = 37500, reads_per_cell = pData(Shalek_read_countdata_cds[, colnames(Shalek_std)])$Total_mRNAs, verbose = T) #
# Shalek_abs_50000= relative2abs(Shalek_std, estimate_t(Shalek_std), modelFormulaStr = "~1", cores=detectCores(), expected_total_mRNAs = 50000, reads_per_cell = pData(Shalek_read_countdata_cds[, colnames(Shalek_std)])$Total_mRNAs, verbose = T) #
# Shalek_abs_50000_subset= relative2abs(Shalek_std[, pData(Shalek_abs)$experiment_name %in% c('Ifnar1_KO_LPS', 'Stat1_KO_LPS',  "LPS", "Unstimulated_Replicate")], 
#                     estimate_t(Shalek_std), modelFormulaStr = "~1", cores=detectCores(), expected_total_mRNAs = 50000, reads_per_cell = pData(Shalek_read_countdata_cds[, colnames(Shalek_std)])$Total_mRNAs, verbose = T) #
# Shalek_abs_37500_subset= relative2abs(Shalek_std[, pData(Shalek_abs)$experiment_name %in% c('Ifnar1_KO_LPS', 'Stat1_KO_LPS',  "LPS", "Unstimulated_Replicate")], 
#                     estimate_t(Shalek_std), modelFormulaStr = "~1", cores=detectCores(), expected_total_mRNAs = 37500, reads_per_cell = pData(Shalek_read_countdata_cds[, colnames(Shalek_std)])$Total_mRNAs, verbose = T) #

pd <- new("AnnotatedDataFrame", data = pData(Shalek_std))
fd <- new("AnnotatedDataFrame", data = fData(Shalek_std))
Shalek_abs <-  newCellDataSet(as.matrix(Shalek_abs), 
                             phenoData = pd, 
                             featureData = fd, 
                             expressionFamily=negbinomial(), 
                             lowerDetectionLimit=1)
pData(Shalek_abs)$Total_mRNAs <- colSums(exprs(Shalek_abs))

#Calculate size factors and dispersions
Shalek_abs = estimateSizeFactors(Shalek_abs)
Shalek_abs = estimateDispersions(Shalek_abs)

# Filter to only expressed genes (feel free to change this as needed, I have tried several methods)
Shalek_abs = detectGenes(Shalek_abs, min_expr = 1)

###################################################################################################################################
### performing the DEG tests to obtain the genes used for ordering the cells #####

# LPS: 
# Shalek_LPS <- Shalek_abs[, pData(Shalek_abs)$experiment_name %in% c('LPS', 'Unstimulated_Replicate')]
# pData(Shalek_LPS)[, 'stim_time'] <- as.character(pData(Shalek_LPS)$time)

# pData(Shalek_LPS)$stim_time[pData(Shalek_LPS)$stim_time == ''] <- 0
# Shalek_LPS <- detectGenes(Shalek_LPS, min_expr = 0.1)
# expressed_genes <- row.names(subset(fData(Shalek_LPS), num_cells_expressed > 50))
# genes_in_range <- selectGenesInExpressionRange(Shalek_LPS[expressed_genes,], 2, Inf, 0.1, stat_fun=function(x) { median(round(x)) })
# Shalek_LPS_subset_DEG_res <- differentialGeneTest(Shalek_LPS[genes_in_range, ], fullModelFormulaStr = '~stim_time', cores = detectCores())
# closeAllConnections()

#ko: include all LPS cells
Shalek_abs_subset_ko_LPS <- Shalek_abs[, pData(Shalek_abs)$experiment_name %in% c('Ifnar1_KO_LPS', 'Stat1_KO_LPS',  "LPS", "Unstimulated_Replicate")]
pData(Shalek_abs_subset_ko_LPS)[, 'stim_time'] <- as.character(pData(Shalek_abs_subset_ko_LPS)$time)

pData(Shalek_abs_subset_ko_LPS)$stim_time[pData(Shalek_abs_subset_ko_LPS)$stim_time == ''] <- 0
pData(Shalek_abs_subset_ko_LPS)$stim_time <- as.integer(revalue(pData(Shalek_abs_subset_ko_LPS)$stim_time, c("1h" = 1, "2h" = 2, "4h" = 4, "6h" = 6)))
Shalek_abs_subset_ko_LPS <- detectGenes(Shalek_abs_subset_ko_LPS, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(Shalek_abs_subset_ko_LPS), num_cells_expressed > 50))
genes_in_range <- selectGenesInExpressionRange(Shalek_abs_subset_ko_LPS[expressed_genes,], 2, Inf, 0.1, stat_fun=function(x) { median(round(x)) })
Shalek_abs_subset_ko_LPS_subset_DEG_res <- differentialGeneTest(Shalek_abs_subset_ko_LPS[genes_in_range, ], fullModelFormulaStr = '~experiment_name + stim_time', cores = detectCores())
closeAllConnections()

#make spanning trees (select subset of the CDS for downstream analysis): 
pData(Shalek_abs_subset_ko_LPS)$Total_mRNAs <- colSums(exprs(Shalek_abs_subset_ko_LPS))
Shalek_abs_subset_ko_LPS <- Shalek_abs_subset_ko_LPS[, pData(Shalek_abs_subset_ko_LPS)$Total_mRNAs < 170000]

order_genes <- c(row.names(subset(Shalek_abs_subset_ko_LPS_subset_DEG_res, qval < 1e-30)))

Shalek_abs_subset_ko_LPS <- setOrderingFilter(Shalek_abs_subset_ko_LPS, order_genes)
Shalek_abs_subset_ko_LPS <- reduceDimension(Shalek_abs_subset_ko_LPS, pseudo_expr = 0, method = "ICA", scaling = F) #
# Shalek_abs_subset_ko_LPS <- reduceDimension(Shalek_abs_subset_ko_LPS, use_vst = T, use_irlba=F, pseudo_expr = 0, residualModelFormulaStr = "~num_genes_expressed", scaling = F, method = "ICA")
Shalek_abs_subset_ko_LPS <- orderCells(Shalek_abs_subset_ko_LPS, num_path = 2)

pdf('tmp/Shalek_LPS_tree.pdf')
monocle::plot_spanning_tree(Shalek_abs_subset_ko_LPS, color_by="interaction(experiment_name, time)", cell_size=4) + 
  scale_color_manual(values=shalek_custom_color_scale_plus_states)
dev.off()

save.image('./RData/analysis_shalek_data_tmp.RData')
# # state_1_cell <- 'Unstimulated_Replicate_S47_0'
# # # State_2_cell <- 'LPS_4h_S59_0'
# # State_3_cell <- 'Stat1_KO_LPS_4h_S22_0'

# # #ko data: 
# # root_state <- pData(Shalek_abs_subset_ko_LPS[, state_1_cell])$State
# # Shalek_abs_subset_ko_LPS <- orderCells(Shalek_abs_subset_ko_LPS, root_state = root_state, num_path = 2)
# # if (pData(Shalek_abs_subset_ko_LPS)[State_3_cell, 'State'] != 3) {
# #   State <- pData(Shalek_abs_subset_ko_LPS)$State 
# #   pData(Shalek_abs_subset_ko_LPS)$State[State == 3] <- 2
# #   pData(Shalek_abs_subset_ko_LPS)$State[State == 2] <- 3
# # }

# # Figure 5C -- Heatmap
# # Detect branching genes and calulate ABCs and ILRs
# full_model_string = '~sm.ns(Pseudotime, df = 3)*Lineage'

# ko_branching_genes = branchTest(Shalek_abs_subset_ko_LPS, fullModelFormulaStr = full_model_string, cores = detectCores(), relative_expr = T) #, weighted = T
# closeAllConnections()

# # Figure 5B annotations -- Enrichment analysis on clusters

# # make venn diagram for the genes for figure 5/6:
# # figure 5: 
# # pseudotime test for the WT cells

# ##two group tests: 
# #comparign with time 4h:
# ko_Ifnar1_wt4 <- differentialGeneTest(Shalek_abs_subset_ko_LPS[, c(pData(Shalek_abs_subset_ko_LPS)$experiment_name %in% c('LPS') & 
#                                          pData(Shalek_abs_subset_ko_LPS)$time %in% '4h') | c(pData(Shalek_abs_subset_ko_LPS)$experiment_name %in% c('Ifnar1_KO_LPS') & 
#                                          pData(Shalek_abs_subset_ko_LPS)$time %in% '4h')
#                                          ], fullModelFormulaStr="~experiment_name", reducedModelFormulaStr="~1", cores=detectCores())
# closeAllConnections()
# ko_stat1_wt4 <- differentialGeneTest(Shalek_abs_subset_ko_LPS[, c(pData(Shalek_abs_subset_ko_LPS)$experiment_name %in% c('LPS') & 
#                                          pData(Shalek_abs_subset_ko_LPS)$time %in% '4h') | c(pData(Shalek_abs_subset_ko_LPS)$experiment_name %in% c('Stat1_KO_LPS') & 
#                                          pData(Shalek_abs_subset_ko_LPS)$time %in% '4h')
#                                          ], fullModelFormulaStr="~experiment_name", reducedModelFormulaStr="~1", cores=detectCores())
# closeAllConnections()

# #####################golgi: with all LPS cells: ######################
# Shalek_golgi_update <- Shalek_abs[,pData(Shalek_abs)$experiment_name %in% c("LPS_GolgiPlug", "LPS", "Unstimulated_Replicate")]

# #add both the golgi time and stim time: 
# split_cols <- str_split_fixed(pData(Shalek_golgi_update)$time, '_', 2)
# pData(Shalek_golgi_update)[, 'stim_time'] <- split_cols[, 1]
# pData(Shalek_golgi_update)$stim_time[pData(Shalek_golgi_update)$stim_time == ''] <- 0
# pData(Shalek_golgi_update)$stim_time <- as.numeric(revalue(pData(Shalek_golgi_update)$stim_time, c("1h" = 1, "2h" = 2, "4h" = 4, "6h" = 6)))

# #the predictor cannot be Inf
# pData(Shalek_golgi_update)[, 'golgi_time'] <- split_cols[, 2]
# pData(Shalek_golgi_update)$golgi_time[pData(Shalek_golgi_update)$golgi_time == ''] <- 'NEVER' 

# Shalek_golgi_update <- detectGenes(Shalek_golgi_update, min_expr = 0.1)
# expressed_genes <- row.names(subset(fData(Shalek_golgi_update), num_cells_expressed > 50))
# genes_in_range <- selectGenesInExpressionRange(Shalek_golgi_update[expressed_genes,], 2, Inf, 0.1, stat_fun=function(x) { median(round(x)) })

# Shalek_golgi_update_subset_DEG_res <- differentialGeneTest(Shalek_golgi_update[genes_in_range, ], fullModelFormulaStr = '~stim_time + golgi_time', cores = detectCores())
# closeAllConnections()

# #make spanning trees for golgi-plug: 
# pData(Shalek_golgi_update)$Total_mRNAs <- colSums(exprs(Shalek_golgi_update))
# Shalek_golgi_update <- Shalek_golgi_update[, pData(Shalek_golgi_update)$Total_mRNAs < 75000]

# #select genes for ordering cells: 
# golgi_order_genes <- c(row.names(subset(Shalek_golgi_update_subset_DEG_res, qval < 1e-40)))

# Shalek_golgi_update <- setOrderingFilter(Shalek_golgi_update, golgi_order_genes)
# Shalek_golgi_update <- reduceDimension(Shalek_golgi_update, pseudo_expr = 0, residualModelFormulaStr = "~num_genes_expressed", reduction_method = "ICA")
# Shalek_golgi_update <- orderCells(Shalek_golgi_update, num_path = 2)

# #gologi data: 
# state_1_cell <- 'Unstimulated_Replicate_S75_0'
# State_2_cell <- 'LPS_6h_S48_0'
# # State_3_cell <- 'Stat1_KO_LPS_4h_S22_0'
# root_state <- pData(Shalek_golgi_update[, state_1_cell])$State
# Shalek_golgi_update <- orderCells(Shalek_golgi_update, root_state = root_state, num_path = 2)
# if (pData(Shalek_golgi_update)[State_2_cell, 'State'] != 2) {
#   State <- pData(Shalek_golgi_update)$State 
#   pData(Shalek_golgi_update)$State[State == 3] <- 2
#   pData(Shalek_golgi_update)$State[State == 2] <- 3
# }

# # Figure 6C -- Heatmap of trajectory from 6A

# # Perform branch test and calculate ABCs
# full_model_string = '~sm.ns(Pseudotime, df = 3)*Lineage'

# golgi_branching_genes = branchTest(Shalek_golgi_update, fullModelFormulaStr = full_model_string, cores= detectCores(), relative_expr = T) #, weighted = T
# closeAllConnections()

# #figure 6: 
# #pseudotime test for the WT cells
# golgi_wt_0to4_pseudo <- differentialGeneTest(Shalek_golgi_update[, pData(Shalek_golgi_update)$experiment_name %in% c('LPS', 'Unstimulated_Replicate') & pData(Shalek_golgi_update)$time %in% c('', '1h', '2h', '4h')], fullModelFormulaStr="~sm.ns(Pseudotime, df = 3)", reducedModelFormulaStr="~1", cores=detectCores() / 1)
# closeAllConnections()
# golgi_wt_0to4_pseudo_gene_ids = row.names(subset(golgi_wt_0to4_pseudo, qval < 1e-2))
# golgi_wt_0to6_pseudo <- differentialGeneTest(Shalek_golgi_update[, pData(Shalek_golgi_update)$experiment_name %in% c('LPS', 'Unstimulated_Replicate')], fullModelFormulaStr="~sm.ns(Pseudotime, df = 3)", reducedModelFormulaStr="~1", cores=detectCores() )
# closeAllConnections()
# golgi_wt_0to6_pseudo_gene_ids = row.names(subset(golgi_wt_0to6_pseudo, qval < 1e-2))

# ##two group tests: 
# #test all Golgi plug cells at once: 
# all_golgi_plug0_wt4 <- differentialGeneTest(Shalek_golgi_update[, c(pData(Shalek_golgi_update)$experiment_name %in% c('LPS') & 
#                                          pData(Shalek_golgi_update)$time %in% '4h') | c(pData(Shalek_golgi_update)$experiment_name %in% c('LPS_GolgiPlug'))
#                                          ], fullModelFormulaStr="~time", reducedModelFormulaStr="~1", cores=detectCores())
# closeAllConnections()
# all_golgi_plug0_wt4_gene_ids = row.names(subset(all_golgi_plug0_wt4, qval < 1e-2))

# all_golgi_plug0_wt0 <- differentialGeneTest(Shalek_golgi_update[, c(pData(Shalek_golgi_update)$experiment_name %in% c('Unstimulated_Replicate')) | c(pData(Shalek_golgi_update)$experiment_name %in% c('LPS_GolgiPlug'))
#                                          ], fullModelFormulaStr="~time", reducedModelFormulaStr="~1", cores=detectCores())
# closeAllConnections()
# all_golgi_plug0_wt0_gene_ids = row.names(subset(all_golgi_plug0_wt0, qval < 1e-2))

# #different time comparing to WT 4h
# golgi_plug0_wt4 <- differentialGeneTest(Shalek_golgi_update[, c(pData(Shalek_golgi_update)$experiment_name %in% c('LPS') & 
#                                          pData(Shalek_golgi_update)$time %in% '4h') | c(pData(Shalek_golgi_update)$experiment_name %in% c('LPS_GolgiPlug') & 
#                                          pData(Shalek_golgi_update)$time %in% '4h_0h')
#                                          ], fullModelFormulaStr="~time", reducedModelFormulaStr="~1", cores=detectCores())
# closeAllConnections()
# golgi_plug0_wt4_gene_ids = row.names(subset(golgi_plug0_wt4, qval < 1e-2))

# # #different time comparing to WT 0h
# golgi_plug0_wt0 <- differentialGeneTest(Shalek_golgi_update[, c(pData(Shalek_golgi_update)$experiment_name %in% c('Unstimulated_Replicate')) | c(pData(Shalek_golgi_update)$experiment_name %in% c('LPS_GolgiPlug') & 
#                                          pData(Shalek_golgi_update)$time %in% '4h_0h')
#                                          ], fullModelFormulaStr="~time", reducedModelFormulaStr="~1", cores=detectCores())
# closeAllConnections()
# golgi_plug0_wt0_gene_ids = row.names(subset(golgi_plug0_wt0, qval < 1e-2))

save.image('./RData/analysis_shalek_data.RData')