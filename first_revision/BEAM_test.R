#two group test on states 2/3:    
std_grp_test_lineage23_res <- differentialGeneTest(std_AT12_cds_subset_all_gene[1:transcript_num, pData(std_AT12_cds_subset_all_gene)$State %in% c(2, 3)], 
                                              fullModelFormulaStr = "~State", 
                                              reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)  
abs_grp_test_lineage23_res <- differentialGeneTest(abs_AT12_cds_subset_all_gene[1:transcript_num, pData(abs_AT12_cds_subset_all_gene)$State %in% c(2, 3)], 
                                                fullModelFormulaStr = "~State", 
                                                reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)  
mc_grp_test_lineage23_res <- differentialGeneTest(mc_AT12_cds_subset_all_gene[1:transcript_num, pData(mc_AT12_cds_subset_all_gene)$State %in% c(2, 3)], 
                                                fullModelFormulaStr = "State", 
                                                reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)    

#make venn diagram: 
abs_pseudotime_test_lineage23_res <- c(abs_pseudotime_test_lineage2_res, abs_pseudotime_test_lineage3_res)
beam_element_all <- c(
            row.names(abs_pseudotime_test_lineage23_res[std_grp_test_lineage23_res$qval < .05, ]), #0to6
            row.names(abs_grp_test_lineage23_res[abs_grp_test_lineage23_res$qval < .05, ]), 
            row.names(quake_branching_genes[quake_branching_genes$qval < .05, ]) 
            )
beam_sets_all <- c(
         rep(paste('two group test', sep = ''), length(row.names(abs_pseudotime_test_lineage23_res[abs_pseudotime_test_lineage23_res$qval < .05, ]))),
         rep(paste('union of pseudotime', sep = ''), length(row.names(abs_grp_test_lineage23_res[abs_grp_test_lineage23_res$qval < .05, ]))),
         rep(paste('BEAM test', sep = ''), length(row.names(quake_branching_genes[quake_branching_genes$qval < .05, ])))
         )

pdf(file = paste(fig_root_dir, 'lung_beam_test_benchmark.pdf', sep = ''))
venneuler_venn(beam_element_all, beam_sets_all)
dev.off()

##########################shalek ko data: ##########################
#pseudotime test on lineage 2/3:    
# std_pseudotime_test_lineage2_res <- differentialGeneTest(Shalek_abs_subset_ko_LPS[, pData(Shalek_abs_subset_ko_LPS)$State %in% c(1, 2)], 
#                                               fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
#                                               reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)  
# std_pseudotime_test_lineage3_res <- differentialGeneTest(Shalek_abs_subset_ko_LPS[1:transcript_num, pData(Shalek_abs_subset_ko_LPS)$State %in% c(1, 3)], 
#                                                       fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
#                                                       reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)
ko_abs_pseudotime_test_lineage2_res <- differentialGeneTest(Shalek_abs_subset_ko_LPS[, pData(Shalek_abs_subset_ko_LPS)$State %in% c(1, 2)], 
                                                fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                                reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)  
ko_abs_pseudotime_test_lineage3_res <- differentialGeneTest(Shalek_abs_subset_ko_LPS[, pData(Shalek_abs_subset_ko_LPS)$State %in% c(1, 3)], 
                                                      fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                                      reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)
   

# ko_std_grp_test_lineage2_res <- differentialGeneTest(Shalek_abs_subset_ko_LPS[, pData(Shalek_abs_subset_ko_LPS)$State %in% c(2, 3)], 
#                                               fullModelFormulaStr = "~State", 
#                                               reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)  
ko_abs_grp_test_lineage23_res <- differentialGeneTest(Shalek_abs_subset_ko_LPS[, pData(Shalek_abs_subset_ko_LPS)$State %in% c(2, 3)], 
                                                fullModelFormulaStr = "State", 
                                                reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)    

#make venn diagram: 
ko_abs_pseudotime_test_lineage23_res <- c(ko_abs_pseudotime_test_lineage2_res, ko_abs_pseudotime_test_lineage3_res)
ko_beam_element_all <- c(
            row.names(ko_abs_grp_test_lineage23_res[ko_abs_grp_test_lineage23_res$qval < .05, ]), #0to6
            row.names(ko_abs_pseudotime_test_lineage23_res[ko_abs_pseudotime_test_lineage23_res$qval < .05, ]), 
            row.names(ko_branching_genes[ko_branching_genes$qval < .05, ]) 
            )
ko_beam_sets_all <- c(
         rep(paste('two group test', sep = ''), length(row.names(ko_abs_grp_test_lineage23_res[ko_abs_grp_test_lineage23_res$qval < .05, ]))),
         rep(paste('union of pseudotime', sep = ''), length(row.names(ko_abs_pseudotime_test_lineage23_res[ko_abs_pseudotime_test_lineage23_res$qval < .05, ]))),
         rep(paste('BEAM test', sep = ''), length(row.names(ko_branching_genes[ko_branching_genes$qval < .05, ])))
         )

pdf(file = paste(fig_root_dir, 'ko_beam_test_benchmark.pdf', sep = ''))
venneuler_venn(ko_beam_element_all, ko_beam_sets_all)
dev.off()

##########################shalek golgi data: ##########################
golgi_abs_pseudotime_test_lineage2_res <- differentialGeneTest(Shalek_golgi_update[, pData(Shalek_abs_subset_ko_LPS)$State %in% c(1, 2)], 
                                                fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                                reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)  
golgi_abs_pseudotime_test_lineage3_res <- differentialGeneTest(Shalek_golgi_update[, pData(Shalek_abs_subset_ko_LPS)$State %in% c(1, 3)], 
                                                      fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                                      reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)
   
#two group test: 
# golgi_std_grp_test_lineage23_res <- differentialGeneTest(Shalek_golgi_update[, pData(Shalek_golgi_update)$State %in% c(2, 3)], 
#                                               fullModelFormulaStr = "~State", 
#                                               reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)  
golgi_abs_grp_test_lineage23_res <- differentialGeneTest(Shalek_golgi_update[, pData(Shalek_golgi_update)$State %in% c(2, 3)], 
                                                fullModelFormulaStr = "State", 
                                                reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)    

#make venn diagram: 
golgi_abs_pseudotime_test_lineage23_res <- c(golgi_abs_pseudotime_test_lineage2_res, golgi_abs_pseudotime_test_lineage3_res)
golgi_beam_element_all <- c(
            row.names(golgi_abs_grp_test_lineage23_res[golgi_abs_grp_test_lineage23_res$qval < .05, ]), #0to6
            row.names(golgi_abs_pseudotime_test_lineage23_res[golgi_abs_pseudotime_test_lineage23_res$qval < .05, ]), 
            row.names(golgi_branching_genes[golgi_branching_genes$qval < .05, ]) 
            )
golgi_beam_sets_all <- c(
         rep(paste('two group test', sep = ''), length(row.names(golgi_abs_grp_test_lineage23_res[golgi_abs_grp_test_lineage23_res$qval < .05, ]))),
         rep(paste('union of pseudotime', sep = ''), length(row.names(golgi_abs_pseudotime_test_lineage23_res[golgi_abs_pseudotime_test_lineage23_res$qval < .05, ]))),
         rep(paste('BEAM test', sep = ''), length(row.names(golgi_branching_genes[golgi_branching_genes$qval < .05, ])))
         )

pdf(file = paste(fig_root_dir, 'golgi_beam_test_benchmark.pdf', sep = ''))
venneuler_venn(ko_beam_element_all, ko_beam_sets_all)
dev.off()

