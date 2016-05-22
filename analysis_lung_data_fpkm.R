# library(monocle)
# library(xacHelper)

# load_all_libraries()

# load('./RData/prepare_lung_data.RData')

# # specify states used for branchTest and calABCs
# progenitor_state <- 1
# lineage_states <- c(2, 3)

# ########################################################prepare the HSMM dataset########################################################
# #   #load all the go/reactome/kegg datasets for the analysis: 
# root_directory <- "./data/Quake_data"

# human_go_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Human/symbol/Human_GO_AllPathways_with_GO_iea_June_20_2014_symbol.gmt", sep=""), encoding="latin1")
# names(human_go_gsc$gsc) <- str_split_fixed(names(human_go_gsc$gsc), "%", 2)[,1]

# human_reactome_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Human/symbol/Pathways/Human_Reactome_June_20_2014_symbol.gmt", sep=""), encoding="latin1")
# names(human_reactome_gsc$gsc) <- str_split_fixed(names(human_reactome_gsc$gsc), "%", 2)[,1]

# mouse_go_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/GO/MOUSE_GO_bp_with_GO_iea_symbol.gmt", sep=""), encoding="latin1")
# names(mouse_go_gsc$gsc) <- str_split_fixed(names(mouse_go_gsc$gsc), "%", 2)[,1]

# mouse_reactome_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/Pathways/Mouse_Reactome_June_20_2014_symbol.gmt", sep=""), encoding="latin1")
# names(mouse_reactome_gsc$gsc) <- str_split_fixed(names(mouse_reactome_gsc$gsc), "%", 2)[,1]

# mouse_kegg_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/Pathways/Mouse_Human_KEGG_June_20_2014_symbol.gmt", sep=""), encoding="latin1")
# names(mouse_kegg_gsc$gsc) <- str_split_fixed(names(mouse_kegg_gsc$gsc), "%", 2)[,1]

# mouse_go_gsc_cc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/GO/MOUSE_GO_cc_with_GO_iea_symbol.gmt", sep=""), encoding="latin1")
# names(mouse_go_gsc_cc$gsc) <- str_split_fixed(names(mouse_go_gsc_cc$gsc), "%", 2)[,1]
# mouse_go_gsc_mf <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/GO/MOUSE_GO_mf_with_GO_iea_symbol.gmt", sep=""), encoding="latin1")
# names(mouse_go_gsc_mf$gsc) <- str_split_fixed(names(mouse_go_gsc_mf$gsc), "%", 2)[,1]

# mouse_reactome_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/Pathways/Mouse_Reactome_June_20_2014_symbol.gmt", sep=""), encoding="latin1")
# names(mouse_reactome_gsc$gsc) <- str_split_fixed(names(mouse_reactome_gsc$gsc), "%", 2)[,1]

# mouse_kegg_gsc <- loadGSCSafe(paste(root_directory,"/GMT/EM_pathways/Mouse/by_symbol/Pathways/Mouse_Human_KEGG_June_20_2014_symbol.gmt", sep=""), encoding="latin1")
# names(mouse_kegg_gsc$gsc) <- str_split_fixed(names(mouse_kegg_gsc$gsc), "%", 2)[,1]

# pseudotime_abs_AT12_cds_subset <- differentialGeneTest(abs_AT12_cds_subset_all_gene[1:transcript_num, ], cores = detectCores(), relative_expr = T)
# branch_pseudotime_abs_AT12_cds_subset <- branchTest(abs_AT12_cds_subset_all_gene[1:transcript_num, ], fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", reducedModelFormulaStr = "~1", cores = detectCores(), relative_expr = T, weighted = F)

# abs_AT12_cds_subset_all_gene@dim_reduce_type <- 'ICA'
# weihgted_relative_abs_AT12_cds_subset_all_gene <- branchTest(abs_AT12_cds_subset_all_gene[1:transcript_num, ], cores = detectCores(), relative_expr = T, weighted = T)
# weihgted_relative_abs_AT12_cds_subset_quake_gene <- branchTest(abs_AT12_cds_subset_all_gene[add_quake_gene_all_marker_ids, ], cores = detectCores(), relative_expr = T, weighted = T)

# closeAllConnections()

# abs_AT12_cds_subset_all_gene_ILRs_list <- calILRs(cds = abs_AT12_cds_subset_all_gene[1:transcript_num, ], lineage_states = c(2, 3), stretch = T, cores = detectCores(), 
#   trend_formula = "~sm.ns(Pseudotime, df = 3) * Lineage", ILRs_limit = 3, 
#   relative_expr = T, weighted = FALSE, label_by_short_name = F, 
#   useVST = T, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "abs_AT12_cds_subset_all_gene_ILRs", return_all = T)

# weighted_abs_AT12_cds_subset_all_gene_ILRs_list <- calILRs(cds = abs_AT12_cds_subset_all_gene[1:transcript_num, ], lineage_states = c(2, 3), stretch = T, cores = detectCores(), 
#   trend_formula = "~sm.ns(Pseudotime, df = 3) * Lineage", ILRs_limit = 3, 
#   relative_expr = T, weighted = T, label_by_short_name = F, 
#   useVST = T, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "weighted_abs_AT12_cds_subset_all_gene_ILRs_list", return_all = T)

# #with progenitor duplication
# duplicate_progenitors_weighted_abs_AT12_cds_subset_all_gene_ILRs_list <- calILRs(cds = abs_AT12_cds_subset_all_gene[1:transcript_num, ], lineage_states = c(2, 3), stretch = T, cores = detectCores(), 
#     trend_formula = "~sm.ns(Pseudotime, df = 3) * Lineage", ILRs_limit = 3, progenitor_method = 'duplicate', 
#     relative_expr = T, weighted = T, label_by_short_name = F, 
#     useVST = T, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "duplicate_progenitors_weighted_abs_AT12_cds_subset_all_gene_ILRs_list", return_all = T)

# abs_bifurcation_time <- detectBifurcationPoint(weighted_abs_AT12_cds_subset_all_gene_ILRs_list$norm_str_logfc_df[1:transcript_num, 27:100], div_threshold = 0.3)
# names(abs_bifurcation_time) <- fData(absolute_cds[names(abs_bifurcation_time), ])$gene_short_name

# ####################################################################################################################################################################################################################

# ########################################FPKM value based results: ########################################
# pseudotime_std_AT12_cds_subset <- differentialGeneTest(standard_cds[1:transcript_num, ], cores = detectCores(), relative_expr = T)
# branch_pseudotime_std_AT12_cds_subset <- branchTest(standard_cds[1:transcript_num, ], fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", reducedModelFormulaStr = "~1", cores = detectCores(), relative_expr = T, weighted = F)

# standard_cds@dim_reduce_type <- 'ICA'
# weihgted_std_AT12_cds_subset_all_gene <- branchTest(standard_cds[1:transcript_num, ], fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", reducedModelFormulaStr = "~1", cores = detectCores(), relative_expr = T, weighted = T)
# weihgted_std_AT12_cds_subset_quake_gene <- branchTest(std_AT12_cds_subset_all_gene[add_quake_gene_all_marker_ids, ], fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", reducedModelFormulaStr = "~1", cores = detectCores(), relative_expr = T, weighted = T)

closeAllConnections()

standard_cds
std_AT12_cds_subset_all_gene_ILRs_list <- calILRs(cds = std_AT12_cds_subset_all_gene[1:transcript_num, ], lineage_states = c(2, 3), stretch = T, cores = detectCores(), 
  trend_formula = "~sm.ns(Pseudotime, df = 3) * Lineage", ILRs_limit = 3, 
  relative_expr = T, weighted = FALSE, label_by_short_name = F, 
  useVST = T, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "std_AT12_cds_subset_all_gene_ILRs", return_all = T)

weighted_std_AT12_cds_subset_all_gene_ILRs_list <- calILRs(cds = std_AT12_cds_subset_all_gene[1:transcript_num, ], lineage_states = c(2, 3), stretch = T, cores = detectCores(), 
  trend_formula = "~sm.ns(Pseudotime, df = 3) * Lineage", ILRs_limit = 3, 
  relative_expr = T, weighted = T, label_by_short_name = F, 
  useVST = T, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "weighted_std_AT12_cds_subset_all_gene_ILRs_list", return_all = T)

#with progenitor duplication
duplicate_progenitors_weighted_std_AT12_cds_subset_all_gene_ILRs_list <- calILRs(cds = std_AT12_cds_subset_all_gene[1:transcript_num, ], lineage_states = c(2, 3), stretch = T, cores = detectCores(), 
    trend_formula = "~sm.ns(Pseudotime, df = 3) * Lineage", ILRs_limit = 3, progenitor_method = 'duplicate', 
    relative_expr = T, weighted = T, label_by_short_name = F, 
    useVST = T, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "duplicate_progenitors_weighted_std_AT12_cds_subset_all_gene_ILRs_list", return_all = T)

std_bifurcation_time <- detectBifurcationPoint(weighted_std_AT12_cds_subset_all_gene_ILRs_list$norm_str_logfc_df[1:transcript_num, 27:100], div_threshold = 0.3)
names(abs_bifurcation_time) <- fData(absolute_cds[names(std_bifurcation_time), ])$gene_short_name

########################################mc transcript count value based results: ########################################
pseudotime_mc_AT12_cds_subset <- differentialGeneTest(mc_AT12_cds_subset_all_gene[1:transcript_num, ], cores = detectCores(), relative_expr = T)
branch_pseudotime_mc_AT12_cds_subset <- branchTest(mc_AT12_cds_subset_all_gene[1:transcript_num, ], fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", reducedModelFormulaStr = "~1", cores = detectCores(), relative_expr = T, weighted = F)

mc_adj_cds@dim_reduce_type <- 'ICA'
weihgted_mc_AT12_cds_subset_all_gene <- branchTest(mc_AT12_cds_subset_all_gene[1:transcript_num, ], cores = detectCores(), relative_expr = T, weighted = T)
weihgted_mc_AT12_cds_subset_quake_gene <- branchTest(mc_AT12_cds_subset_all_gene[add_quake_gene_all_marker_ids, ], cores = detectCores(), relative_expr = T, weighted = T)

closeAllConnections()

mc_AT12_cds_subset_all_gene_ILRs_list <- calILRs(cds = mc_AT12_cds_subset_all_gene[1:transcript_num, ], lineage_states = c(2, 3), stretch = T, cores = detectCores(), 
  trend_formula = "~sm.ns(Pseudotime, df = 3) * Lineage", ILRs_limit = 3, 
  relative_expr = T, weighted = FALSE, label_by_short_name = F, 
  useVST = T, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "mc_AT12_cds_subset_all_gene_ILRs", return_all = T)

weighted_mc_AT12_cds_subset_all_gene_ILRs_list <- calILRs(cds = mc_AT12_cds_subset_all_gene[1:transcript_num, ], lineage_states = c(2, 3), stretch = T, cores = detectCores(), 
  trend_formula = "~sm.ns(Pseudotime, df = 3) * Lineage", ILRs_limit = 3, 
  relative_expr = T, weighted = T, label_by_short_name = F, 
  useVST = T, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "weighted_mc_AT12_cds_subset_all_gene_ILRs_list", return_all = T)

#with progenitor duplication
duplicate_progenitors_weighted_mc_AT12_cds_subset_all_gene_ILRs_list <- calILRs(cds = mc_AT12_cds_subset_all_gene[1:transcript_num, ], lineage_states = c(2, 3), stretch = T, cores = detectCores(), 
    trend_formula = "~sm.ns(Pseudotime, df = 3) * Lineage", ILRs_limit = 3, progenitor_method = 'duplicate', 
    relative_expr = T, weighted = T, label_by_short_name = F, 
    useVST = T, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "duplicate_progenitors_weighted_mc_AT12_cds_subset_all_gene_ILRs_list", return_all = T)

mc_bifurcation_time <- detectBifurcationPoint(weighted_mc_AT12_cds_subset_all_gene_ILRs_list$norm_str_logfc_df[1:transcript_num, 27:100], div_threshold = 0.3)
names(mc_bifurcation_time) <- fData(absolute_cds[names(mc_bifurcation_time), ])$gene_short_name

save.image('./RData/analysis_lung_data_fpkm.RData')






