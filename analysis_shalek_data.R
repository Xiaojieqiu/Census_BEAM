# ################################################################################################################################################
#  # The following script is based on Adnrew's work #
#  ################################################################################################################################################

#  #install all packages: 
#  packages = c("ggplot2", "VGAM", "igraph", "pRlyr", "combinat", "fastICA", "irlba", "matrixStats", "reshape2", "R.utils", "snow", 
#             "stringr", "modeest", "Hmisc", "boot", "doMC", "data.table", "fitdistrplus", "ggdendro", "gplots", "princurve", "sp",
#             "lmtest", "MASS", "mixsmsn", "pheatmap", "plyr", "pscl", "RColorBrewer", "VennDiagram", "zoo", "raster", "colorRamps", "grid")
#  install.packages(packages, repo = 'http://cran.fhcrc.org/')

#  bio_packages = c("Biobase", "BiocGenerics",  "limma", "edgeR", "DESeq", "DESeq2", "piano")
#  source("http://bioconductor.org/biocLite.R")
#  biocLite(bio_packages)

 # go to https://github.com/settings/tokens and generate personal tokens for install the private monocle / devtree package: 
 # install_github("cole-trapnell-lab/monocle-dev", auth_token = "2b5f9747e17c8512f1ecd2bf76f5df4730be21e2")
 # install_github("cole-trapnell-lab/branch-diff", auth_token = "2b5f9747e17c8512f1ecd2bf76f5df4730be21e2")

 # install.packages('./xacHelper_0.0.0.9000.tar.gz', dependencies = TRUE)
 # install.packages('./monocle_1.99.0.tar.gz', dependencies = TRUE)
 library(monocle)
 library(xacHelper)
 library(igraph)

 load_all_libraries()
 
 #   #load all the go/reactome/kegg datasets for the analysis: 
 root_directory <- "./Quake_data"

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
#  # load('xiaojie_test_data.gz') #make the dataset
 source('./monocle_helper_functions.R')
 Shalek_valid_genes <- read.table('./Aviv_data/valid_genes_for_analyis.txt', header = T)
 Shalek_exprs_mat <- read.table("./Aviv_data/cuffnorm_output_files/genes.fpkm_table", row.names = 1, header = T)
 Shalek_fd <- read.table("./Aviv_data/cuffnorm_output_files/genes.attr_table", row.names = 1, header = T)
 Shalek_pd <- read.table("./Aviv_data/sample_metadata_table.txt", sep = '\t', row.names = 1, header = T)
 rownames(Shalek_pd) <- paste(rownames(Shalek_pd), "_0", sep = "")
 Shalek_exprs_mat <- Shalek_exprs_mat[row.names(Shalek_fd), row.names(Shalek_pd)]
 Shalek_std <- newCellDataSet(Shalek_exprs_mat, 
                             phenoData = new("AnnotatedDataFrame", data = Shalek_pd), 
                             featureData = new("AnnotatedDataFrame", data = Shalek_fd), 
                             expressionFamily=tobit(), 
                             lowerDetectionLimit=1)
 Shalek_std <- Shalek_std[, which(pData(Shalek_std)$used_in_study == T)]
 Shalek_std <- Shalek_std[row.names(Shalek_std) %in% Shalek_valid_genes$gene_id, ] #27386 * 1787 cells
 #check the consistency with the current Shalek_abs data: 

 Shalek_isoform_fpkm_matrix <- read.table("./Aviv_data/cuffnorm_output_files/isoforms.fpkm_table", row.names = 1, header = T)
 # colnames(Shalek_isoform_fpkm_matrix) <- str_replace(colnames(Shalek_isoform_fpkm_matrix), "_0$", "")
 #colnames(isoform_fpkm_matrix) <- str_replace(colnames(isoform_fpkm_matrix), "GolgiPlugh", "GolgiPlug")
 # row.names(Shalek_isoform_fpkm_matrix) <- Shalek_isoform_fpkm_matrix$tracking_id
 Shalek_isoform_fpkm_matrix <- Shalek_isoform_fpkm_matrix[, colnames(Shalek_std)]

 # save(Shalek_exprs_mat, Shalek_pd, Shalek_fd, Shalek_isoform_fpkm_matrix, Shalek_valid_genes, file = 'AvivDC_cell.RData') #data for making the help package 
 # save(Shalek_valid_genes, file = 'Shalek_valid_genes')
 
 # Convert expression measurements from FPKM to absolute transcript counts, using the isoforms object to estimate the t parameter
 Shalek_abs= relative2abs(Shalek_std, estimate_t(Shalek_isoform_fpkm_matrix), modelFormulaStr = "~1", cores=detectCores())

 pd <- new("AnnotatedDataFrame", data = pData(Shalek_std))
 fd <- new("AnnotatedDataFrame", data = fData(Shalek_std))
 Shalek_abs <-  newCellDataSet(Shalek_abs, 
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
 # expressed_genes = row.names(subset(fData(Shalek_abs), num_cells_expressed >= 50))
 # Shalek_abs_cp <- Shalek_abs
 # Shalek_abs = Shalek_abs[expressed_genes,]
####################################################################################################################################

 Shalek_abs <- estimateSizeFactors(Shalek_abs)
 Shalek_abs <- estimateDispersions(Shalek_abs)

 # fd <- read.delim("/Users/xqiu/Dropbox (Cole Trapnell's Lab)/Shared Data/Regev DC/fData.txt")
 # exprs_mat <- read.delim("/Users/xqiu/Dropbox (Cole Trapnell's Lab)/Shared Data/Regev DC/gene_exprs_counts.txt")
 # pd <-  read.delim("/Users/xqiu/Dropbox (Cole Trapnell's Lab)/Shared Data/Regev DC/pData.txt")

 # fd <- read.delim("./Aviv_data/cuffnorm_output_files/fData.txt")
 # exprs_mat <- read.delim("./Aviv_data/cuffnorm_output_files/gene_exprs_counts.txt")
 # pd <-  read.delim("./Aviv_data/cuffnorm_output_files/pData.txt")


 # Shalek_abs <- newCellDataSet(exprs_mat, 
 #                               phenoData = new("AnnotatedDataFrame", data = pd), 
 #                               featureData = new("AnnotatedDataFrame", data = fd), 
 #                               expressionFamily=negbinomial(), 
 #                               lowerDetectionLimit=1)
 # Shalek_abs <- estimateSizeFactors(Shalek_abs)
 # Shalek_abs <- estimateDispersions(Shalek_abs)

 ###################################################################################################################################
 ### performing the DEG tests to obtain the genes used for ordering the cells #####

 # LPS: 
 Shalek_LPS <- Shalek_abs[, pData(Shalek_abs)$experiment_name %in% c('LPS', 'Unstimulated_Replicate')]
 pData(Shalek_LPS)[, 'stim_time'] <- as.character(pData(Shalek_LPS)$time)

 pData(Shalek_LPS)$stim_time[pData(Shalek_LPS)$stim_time == ''] <- 0
 Shalek_LPS <- detectGenes(Shalek_LPS, min_expr = 0.1)
 expressed_genes <- row.names(subset(fData(Shalek_LPS), num_cells_expressed > 50))
 genes_in_range <- selectGenesInExpressionRange(Shalek_LPS[expressed_genes,], 2, Inf, 0.1, stat_fun=function(x) { median(round(x)) })
 Shalek_LPS_subset_DEG_res <- differentialGeneTest(Shalek_LPS[genes_in_range, ], fullModelFormulaStr = '~stim_time', cores = detectCores() / 2)

 #ko: include all LPS cells
 Shalek_abs_subset_ko_LPS <- Shalek_abs[, pData(Shalek_abs)$experiment_name %in% c('Ifnar1_KO_LPS', 'Stat1_KO_LPS',  "LPS", "Unstimulated_Replicate")]
 pData(Shalek_abs_subset_ko_LPS)[, 'stim_time'] <- as.character(pData(Shalek_abs_subset_ko_LPS)$time)

 pData(Shalek_abs_subset_ko_LPS)$stim_time[pData(Shalek_abs_subset_ko_LPS)$stim_time == ''] <- 0
 pData(Shalek_abs_subset_ko_LPS)$stim_time <- as.integer(revalue(pData(Shalek_abs_subset_ko_LPS)$stim_time, c("1h" = 1, "2h" = 2, "4h" = 4, "6h" = 6)))
 Shalek_abs_subset_ko_LPS <- detectGenes(Shalek_abs_subset_ko_LPS, min_expr = 0.1)
 expressed_genes <- row.names(subset(fData(Shalek_abs_subset_ko_LPS), num_cells_expressed > 50))
 genes_in_range <- selectGenesInExpressionRange(Shalek_abs_subset_ko_LPS[expressed_genes,], 2, Inf, 0.1, stat_fun=function(x) { median(round(x)) })
 Shalek_abs_subset_ko_LPS_subset_DEG_res <- differentialGeneTest(Shalek_abs_subset_ko_LPS[genes_in_range, ], fullModelFormulaStr = '~experiment_name + stim_time', cores = detectCores() / 2)

 #make spanning trees (select subset of the CDS for downstream analysis): 
 pData(Shalek_abs_subset_ko_LPS)$Total_mRNAs <- colSums(exprs(Shalek_abs_subset_ko_LPS))
 Shalek_abs_subset_ko_LPS <- Shalek_abs_subset_ko_LPS[, pData(Shalek_abs_subset_ko_LPS)$Total_mRNAs < 75000]

 order_genes <- c(row.names(subset(Shalek_abs_subset_ko_LPS_subset_DEG_res, qval < 1e-40)))

 Shalek_abs_subset_ko_LPS <- setOrderingFilter(Shalek_abs_subset_ko_LPS, order_genes)
 Shalek_abs_subset_ko_LPS <- reduceDimension(Shalek_abs_subset_ko_LPS, use_vst = T, use_irlba=F, pseudo_expr = 0, covariates = as.vector(pData(Shalek_abs_subset_ko_LPS)$num_genes_expressed) )
 #save.image('~/Projects/BEAM/Parallel_the_reproduce/tmp_analysis_shalek_data.RData')
 Shalek_abs_subset_ko_LPS <- orderCells(Shalek_abs_subset_ko_LPS, num_path = 2)

 # Figure 5C -- Heatmap
 ## Detect branching genes and calulate ABCs and ILRs
 full_model_string = '~sm.ns(Pseudotime, df = 3)*Lineage'

 ko_branching_genes = branchTest(Shalek_abs_subset_ko_LPS, fullModelFormulaStr = full_model_string, cores = detectCores(), relative_expr = T, weighted = T)
 # Shalek_abs_subset_ko_LPS_abcs = calABCs(Shalek_abs_subset_ko_LPS, fullModelFormulaStr=full_model_string, cores=1) 

 ## Generate a plot
 regev_cat <- read.table(file = './Aviv_data/cuffnorm_output_files/study_gene_categories.txt', header = T, sep = '\t') #pass this to the plot_genes_branched_heatmap function

 # Figure 5B annotations -- Enrichment analysis on clusters
 
#save_hyper_df(Shalek_abs_subset_ko_LPS_tree_hyper_geometric_results_reactome, 'ko_hyper_df.xls') 

 #make venn diagram for the genes for figure 5/6:
 #figure 5: 
 #pseudotime test for the WT cells

 # ##two group tests: 
 # #comparign with time 4h:
 ko_Ifnar1_wt4 <- differentialGeneTest(Shalek_abs_subset_ko_LPS[, c(pData(Shalek_abs_subset_ko_LPS)$experiment_name %in% c('LPS') & 
                                             pData(Shalek_abs_subset_ko_LPS)$time %in% '4h') | c(pData(Shalek_abs_subset_ko_LPS)$experiment_name %in% c('Ifnar1_KO_LPS') & 
                                             pData(Shalek_abs_subset_ko_LPS)$time %in% '4h')
                                             ], fullModelFormulaStr="~experiment_name", reducedModelFormulaStr="~1", cores=detectCores() / 2)
 ko_stat1_wt4 <- differentialGeneTest(Shalek_abs_subset_ko_LPS[, c(pData(Shalek_abs_subset_ko_LPS)$experiment_name %in% c('LPS') & 
                                             pData(Shalek_abs_subset_ko_LPS)$time %in% '4h') | c(pData(Shalek_abs_subset_ko_LPS)$experiment_name %in% c('Stat1_KO_LPS') & 
                                             pData(Shalek_abs_subset_ko_LPS)$time %in% '4h')
                                             ], fullModelFormulaStr="~experiment_name", reducedModelFormulaStr="~1", cores=detectCores() / 2)

 #####################golgi: with all LPS cells: ######################
 Shalek_golgi_update <- Shalek_abs[,pData(Shalek_abs)$experiment_name %in% c("LPS_GolgiPlug", "LPS", "Unstimulated_Replicate")]

 #add both the golgi time and stim time: 
 split_cols <- str_split_fixed(pData(Shalek_golgi_update)$time, '_', 2)
 pData(Shalek_golgi_update)[, 'stim_time'] <- split_cols[, 1]
 pData(Shalek_golgi_update)$stim_time[pData(Shalek_golgi_update)$stim_time == ''] <- 0
 pData(Shalek_golgi_update)$stim_time <- as.numeric(revalue(pData(Shalek_golgi_update)$stim_time, c("1h" = 1, "2h" = 2, "4h" = 4, "6h" = 6)))

 #the predictor cannot be Inf
 pData(Shalek_golgi_update)[, 'golgi_time'] <- split_cols[, 2]
 pData(Shalek_golgi_update)$golgi_time[pData(Shalek_golgi_update)$golgi_time == ''] <- 'NEVER' 

 Shalek_golgi_update <- detectGenes(Shalek_golgi_update, min_expr = 0.1)
 expressed_genes <- row.names(subset(fData(Shalek_golgi_update), num_cells_expressed > 50))
 genes_in_range <- selectGenesInExpressionRange(Shalek_golgi_update[expressed_genes,], 2, Inf, 0.1, stat_fun=function(x) { median(round(x)) })

 Shalek_golgi_update_subset_DEG_res <- differentialGeneTest(Shalek_golgi_update[genes_in_range, ], fullModelFormulaStr = '~stim_time + golgi_time', cores = detectCores() / 2)

 #make spanning trees for golgi-plug: 
 pData(Shalek_golgi_update)$Total_mRNAs <- colSums(exprs(Shalek_golgi_update))
 Shalek_golgi_update <- Shalek_golgi_update[, pData(Shalek_golgi_update)$Total_mRNAs < 75000]

 #select genes for ordering cells: 
 golgi_order_genes <- c(row.names(subset(Shalek_golgi_update_subset_DEG_res, qval < 1e-40)))

 Shalek_golgi_update <- setOrderingFilter(Shalek_golgi_update, golgi_order_genes)
 Shalek_golgi_update <- reduceDimension(Shalek_golgi_update, use_vst = T, use_irlba=F, pseudo_expr = 0, covariates = as.vector(pData(Shalek_golgi_update)$num_genes_expressed) )
 Shalek_golgi_update <- orderCells(Shalek_golgi_update, num_path = 2)

  # Figure 6C -- Heatmap of trajectory from 6A

 ## Perform branch test and calculate ABCs
 full_model_string = '~sm.ns(Pseudotime, df = 3)*Lineage'

 golgi_branching_genes = branchTest(Shalek_golgi_update, fullModelFormulaStr = full_model_string, cores= detectCores(), relative_expr = T, weighted = T)
 # ABCs_golgi = calABCs(Shalek_golgi_update, fullModelFormulaStr=full_model_string, cores=1)  # change to calABCs once fixed

 #figure 6: 
 #pseudotime test for the WT cells
 golgi_wt_0to4_pseudo <- differentialGeneTest(Shalek_golgi_update[, pData(Shalek_golgi_update)$experiment_name %in% c('LPS', 'Unstimulated_Replicate') & pData(Shalek_golgi_update)$time %in% c('', '1h', '2h', '4h')], fullModelFormulaStr="~sm.ns(Pseudotime, df = 3)", reducedModelFormulaStr="~1", cores=detectCores() / 1)
 golgi_wt_0to4_pseudo_gene_ids = row.names(subset(golgi_wt_0to4_pseudo, qval < 1e-2))
 golgi_wt_0to6_pseudo <- differentialGeneTest(Shalek_golgi_update[, pData(Shalek_golgi_update)$experiment_name %in% c('LPS', 'Unstimulated_Replicate')], fullModelFormulaStr="~sm.ns(Pseudotime, df = 3)", reducedModelFormulaStr="~1", cores=detectCores() / 2)
 golgi_wt_0to6_pseudo_gene_ids = row.names(subset(golgi_wt_0to6_pseudo, qval < 1e-2))

 ##two group tests: 
 #test all Golgi plug cells at once: 
 all_golgi_plug0_wt4 <- differentialGeneTest(Shalek_golgi_update[, c(pData(Shalek_golgi_update)$experiment_name %in% c('LPS') & 
                                             pData(Shalek_golgi_update)$time %in% '4h') | c(pData(Shalek_golgi_update)$experiment_name %in% c('LPS_GolgiPlug'))
                                             ], fullModelFormulaStr="~time", reducedModelFormulaStr="~1", cores=detectCores() - 1)
 all_golgi_plug0_wt4_gene_ids = row.names(subset(all_golgi_plug0_wt4, qval < 1e-2))

 all_golgi_plug0_wt0 <- differentialGeneTest(Shalek_golgi_update[, c(pData(Shalek_golgi_update)$experiment_name %in% c('Unstimulated_Replicate')) | c(pData(Shalek_golgi_update)$experiment_name %in% c('LPS_GolgiPlug'))
                                             ], fullModelFormulaStr="~time", reducedModelFormulaStr="~1", cores=detectCores() - 1)
 all_golgi_plug0_wt0_gene_ids = row.names(subset(all_golgi_plug0_wt0, qval < 1e-2))

 #different time comparing to WT 4h
 golgi_plug0_wt4 <- differentialGeneTest(Shalek_golgi_update[, c(pData(Shalek_golgi_update)$experiment_name %in% c('LPS') & 
                                             pData(Shalek_golgi_update)$time %in% '4h') | c(pData(Shalek_golgi_update)$experiment_name %in% c('LPS_GolgiPlug') & 
                                             pData(Shalek_golgi_update)$time %in% '4h_0h')
                                             ], fullModelFormulaStr="~time", reducedModelFormulaStr="~1", cores=detectCores())
 golgi_plug0_wt4_gene_ids = row.names(subset(golgi_plug0_wt4, qval < 1e-2))

 # #different time comparing to WT 0h
 golgi_plug0_wt0 <- differentialGeneTest(Shalek_golgi_update[, c(pData(Shalek_golgi_update)$experiment_name %in% c('Unstimulated_Replicate')) | c(pData(Shalek_golgi_update)$experiment_name %in% c('LPS_GolgiPlug') & 
                                             pData(Shalek_golgi_update)$time %in% '4h_0h')
                                             ], fullModelFormulaStr="~time", reducedModelFormulaStr="~1", cores=detectCores())
 golgi_plug0_wt0_gene_ids = row.names(subset(golgi_plug0_wt0, qval < 1e-2))

save.image('shalek_data_analysis.RData')
