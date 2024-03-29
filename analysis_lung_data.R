  library(monocle)
  library(xacHelper)

  load_all_libraries()

  load('prepare_lung_data.RData')
  
  # specify states used for branchTest and calABCs
  progenitor_state <- 1
  lineage_states <- c(2, 3)

  # #   ########################################################prepare the HSMM dataset########################################################
  # #   
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

  pseudotime_abs_AT12_cds_subset <- differentialGeneTest(abs_AT12_cds_subset_all_gene[1:transcript_num, ], cores = detectCores(), relative_expr = T)
  branch_pseudotime_abs_AT12_cds_subset <- branchTest(abs_AT12_cds_subset_all_gene[1:transcript_num, ], fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", reducedModelFormulaStr = "~1", cores = detectCores(), relative_expr = T, weighted = F)
  
  # abs_AT12_cds_subset_all_gene_res <- branchTest(abs_AT12_cds_subset_all_gene[1:transcript_num, ], cores = detectCores(), relative_expr = F, weighted = F)
  # mc_AT12_cds_subset_all_gene_res <- branchTest(mc_AT12_cds_subset_all_gene[1:transcript_num, ], cores = detectCores(), relative_expr = F, weighted = F)

  #perform relative expression analysis: 

  # relative_abs_AT12_cds_subset_all_gene <- branchTest(abs_AT12_cds_subset_all_gene[1:transcript_num, ], cores = detectCores(), relative_expr = T, weighted = F)
  # relative_abs_AT12_cds_subset_quake_gene <- branchTest(abs_AT12_cds_subset_all_gene[add_quake_gene_all_marker_ids, ], cores = 1, relative_expr = T, weighted = F)

  # relative_mc_AT12_cds_subset_all_gene <- branchTest(mc_AT12_cds_subset_all_gene[1:transcript_num, ], cores = detectCores(), relative_expr = T, weighted = F)
  # relative_mc_AT12_cds_subset_quake_gene <- branchTest(mc_AT12_cds_subset_all_gene[add_quake_gene_all_marker_ids, ], cores = 1, relative_expr = T, weighted = F)

  weihgted_relative_abs_AT12_cds_subset_all_gene <- branchTest(abs_AT12_cds_subset_all_gene[1:transcript_num, ], cores = detectCores(), relative_expr = T, weighted = T)
  weihgted_relative_abs_AT12_cds_subset_quake_gene <- branchTest(abs_AT12_cds_subset_all_gene[add_quake_gene_all_marker_ids, ], cores = 1, relative_expr = T, weighted = T)

  # weihgted_relative_mc_AT12_cds_subset_all_gene <- branchTest(mc_AT12_cds_subset_all_gene[1:transcript_num, ], cores = detectCores(), relative_expr = T, weighted = T)
  # weihgted_relative_mc_AT12_cds_subset_quake_gene <- branchTest(mc_AT12_cds_subset_all_gene[add_quake_gene_all_marker_ids, ], cores = 1, relative_expr = T, weighted = T)
  
  # save(relative_abs_AT12_cds_subset_all_gene, relative_abs_AT12_cds_subset_quake_gene, relative_mc_AT12_cds_subset_all_gene, relative_mc_AT12_cds_subset_quake_gene, 
  #   weihgted_relative_abs_AT12_cds_subset_all_gene, weihgted_relative_abs_AT12_cds_subset_quake_gene, weihgted_relative_mc_AT12_cds_subset_all_gene, weihgted_relative_mc_AT12_cds_subset_quake_gene,
  #     file = 'branchTest_res_update')

  ## calILRs: 
  # std_AT12_cds_subset_all_gene_ILRs <- monocle::calILRs(cds = std_AT12_cds_subset_all_gene[add_quake_gene_all_marker_ids, ], lineage_states = c(2, 3), stretch = T, cores = 1, 
  #   trend_formula = "~sm.ns(Pseudotime, df = 3)", ILRs_limit = 3, 
  #   relative_expr = F, weighted = FALSE, label_by_short_name = F, 
  #   useVST = FALSE, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "std_AT12_cds_subset_all_gene_ILRs")
  
  # abs_AT12_cds_subset_all_gene_ILRs <- calILRs(cds = abs_AT12_cds_subset_all_gene[, ], lineage_states = c(2, 3), stretch = T, cores = 1, 
  #   trend_formula = "~sm.ns(Pseudotime, df = 3) * Lineage", ILRs_limit = 3, 
  #   relative_expr = T, weighted = FALSE, label_by_short_name = F, 
  #   useVST = T, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "abs_AT12_cds_subset_all_gene_ILRs")

  abs_AT12_cds_subset_all_gene_ILRs_list <- calILRs(cds = abs_AT12_cds_subset_all_gene[1:transcript_num, ], lineage_states = c(2, 3), stretch = T, cores = 1, 
    trend_formula = "~sm.ns(Pseudotime, df = 3) * Lineage", ILRs_limit = 3, 
    relative_expr = T, weighted = FALSE, label_by_short_name = F, 
    useVST = T, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "abs_AT12_cds_subset_all_gene_ILRs", return_all = T)

  weighted_abs_AT12_cds_subset_all_gene_ILRs_list <- calILRs(cds = abs_AT12_cds_subset_all_gene[1:transcript_num, ], lineage_states = c(2, 3), stretch = T, cores = 1, 
    trend_formula = "~sm.ns(Pseudotime, df = 3) * Lineage", ILRs_limit = 3, 
    relative_expr = T, weighted = T, label_by_short_name = F, 
    useVST = T, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "weighted_abs_AT12_cds_subset_all_gene_ILRs_list", return_all = T)

  # abs_AT12_cds_subset_all_gene_ILRs <- monocle::calILRs(cds = abs_AT12_cds_subset_all_gene[add_quake_gene_all_marker_ids, ], lineage_states = c(2, 3), stretch = T, cores =1, 
  #   trend_formula = "~sm.ns(Pseudotime, df = 3)", ILRs_limit = 3, 
  #   relative_expr = F, weighted = FALSE, label_by_short_name = F, 
  #   useVST = FALSE, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "abs_AT12_cds_subset_all_gene_ILRs")

  # mc_AT12_cds_subset_all_gene_ILRs <- calILRs(cds = mc_AT12_cds_subset_all_gene[1:transcript_num, ], lineage_states = c(2, 3), stretch = T, cores = 1, 
  #   trend_formula = "~sm.ns(Pseudotime, df = 3) * Lineage", ILRs_limit = 3, 
  #   relative_expr = T, weighted = FALSE, label_by_short_name = F, 
  #   useVST = FALSE, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "mc_AT12_cds_subset_all_gene_ILRs")

  # calculate the priority:  
  # std_bifurcation_time <- detectBifurcationPoint(std_AT12_cds_subset_all_gene_ILRs[1:transcript_num, 27:100], div_threshold = 0.3)
  # names(std_bifurcation_time) <- fData(absolute_cds[names(std_bifurcation_time), ])$gene_short_name

  abs_bifurcation_time <- detectBifurcationPoint(weighted_abs_AT12_cds_subset_all_gene_ILRs_list$norm_str_logfc_df[1:transcript_num, 27:100], div_threshold = 0.3)
  names(abs_bifurcation_time) <- fData(absolute_cds[names(abs_bifurcation_time), ])$gene_short_name

  # mc_bifurcation_time <- detectBifurcationPoint(mc_AT12_cds_subset_all_gene_ILRs[1:transcript_num, 27:100], div_threshold = 0.3)
  # names(mc_bifurcation_time) <- fData(absolute_cds[names(abs_bifurcation_time), ])$gene_short_name

  # calculate the ABCs: 
  # std_AT12_cds_subset_all_gene_ABCs <- calABCs(std_AT12_cds_subset_all_gene[add_quake_gene_all_marker_ids, ], trend_formula = "~sm.ns(Pseudotime, df = 3)*Lineage",
  #                                              reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)",
  #                                              branchTest = FALSE,
  #                                              lineage_states = c(2, 3),
  #                                              relative_expr = FALSE,
  #                                              stretch = TRUE,
  #                                              pseudocount=0,
  #                                              cores = 1,
  #                                              weighted = TRUE,
  #                                              min_expr = 0.5,
  #                                              integer_expression = FALSE,
  #                                              num = 5000)
  # abs_AT12_cds_subset_all_gene_ABCs <- calABCs(abs_AT12_cds_subset_all_gene[1:transcript_num, ], trend_formula = "~sm.ns(Pseudotime, df = 3)*Lineage",
  #                                              reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)",
  #                                              branchTest = FALSE,
  #                                              lineage_states = c(2, 3),
  #                                              relative_expr = TRUE,
  #                                              stretch = TRUE,
  #                                              pseudocount=0,
  #                                              cores = 1,
  #                                              weighted = TRUE,
  #                                              min_expr = 0.5,
  #                                              integer_expression = F,
  #                                              num = 5000)
  # mc_AT12_cds_subset_all_gene_ABCs <- calABCs(mc_AT12_cds_subset_all_gene[1:transcript_num, ], trend_formula = "~sm.ns(Pseudotime, df = 3)*Lineage",
  #                                             reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)",
  #                                             branchTest = FALSE,
  #                                             lineage_states = c(2, 3),
  #                                             relative_expr = TRUE,
  #                                             stretch = TRUE,
  #                                             pseudocount=0,
  #                                             cores = 1,
  #                                             weighted = TRUE,
  #                                             min_expr = 0.5,
  #                                             integer_expression = FALSE,
  #                                             num = 5000)
  
  save.image('analysis_lung_data.RData')

