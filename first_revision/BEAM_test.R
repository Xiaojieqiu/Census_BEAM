load('./RData/analysis_shalek_data.RData')
load('./RData/analysis_lung_data.RData')
library(monocle)
library(xacHelper)

fig_root_dir <- './supplementary_figures/'

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

abs_pseudotime_test_lineage2_res <- differentialGeneTest(abs_AT12_cds_subset_all_gene[1:transcript_num, pData(abs_AT12_cds_subset_all_gene)$State %in% c(1, 2)], 
                                                fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                                reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)  
abs_pseudotime_test_lineage3_res <- differentialGeneTest(abs_AT12_cds_subset_all_gene[1:transcript_num, pData(abs_AT12_cds_subset_all_gene)$State %in% c(1, 3)], 
                                                      fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                                      reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

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
                                                fullModelFormulaStr = "~State", 
                                                reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)    

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
                                                fullModelFormulaStr = "~State", 
                                                reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)    

#make venn diagram: 
abs_pseudotime_test_lineage23_res_gene_names <- intersect(row.names(abs_pseudotime_test_lineage2_res[abs_pseudotime_test_lineage2_res$qval < 0.05, ]),
									 row.names(abs_pseudotime_test_lineage3_res[abs_pseudotime_test_lineage3_res$qval < 0.05, ]))
beam_element_all <- c(
   # abs_pseudotime_test_lineage23_res_gene_names, 
            row.names(abs_grp_test_lineage23_res[abs_grp_test_lineage23_res$qval < .05, ]), 
            quake_branch_genes #  %in% valid_expressed_genes
            )
beam_sets_all <- c(
         # rep(paste('intersection of pseudotime test', sep = ''), length(abs_pseudotime_test_lineage23_res_gene_names)),
         rep(paste('two group test', sep = ''), sum(abs_grp_test_lineage23_res$qval < .05)),
         rep(paste('BEAM test', sep = ''), length(quake_branch_genes))
         )

pdf(file = paste(fig_root_dir, 'lung_beam_test_benchmark.pdf', sep = ''))
venneuler_venn(beam_element_all, beam_sets_all)
dev.off()

intersect(row.names(abs_grp_test_lineage23_res[abs_grp_test_lineage23_res$qval < .05, ]), 
            quake_branch_genes )
table(andrew_sets_all)
#make venn diagram: 
ko_abs_pseudotime_test_lineage23_res <- intersect(row.names(ko_abs_pseudotime_test_lineage2_res[ko_abs_pseudotime_test_lineage2_res$qval < 0.05, ]),
									 row.names(ko_abs_pseudotime_test_lineage3_res[ko_abs_pseudotime_test_lineage3_res$qval < 0.05, ]))
ko_beam_element_all <- c(
            row.names(ko_abs_grp_test_lineage23_res[ko_abs_grp_test_lineage23_res$qval < .05, ]), #0to6
            # ko_abs_pseudotime_test_lineage23_res, 
            row.names(ko_branching_genes[ko_branching_genes$qval < .05, ]) 
            )
ko_beam_sets_all <- c(
         rep(paste('two group test', sep = ''), length(row.names(ko_abs_grp_test_lineage23_res[ko_abs_grp_test_lineage23_res$qval < .05, ]))),
         # rep(paste('intersection of pseudotime test', sep = ''), length(ko_abs_pseudotime_test_lineage23_res)),
         rep(paste('BEAM test', sep = ''), length(row.names(ko_branching_genes[ko_branching_genes$qval < .05, ])))
         )

pdf(file = paste(fig_root_dir, 'ko_beam_test_benchmark.pdf', sep = ''))
venneuler_venn(ko_beam_element_all, ko_beam_sets_all)
dev.off()

intersect(row.names(ko_abs_grp_test_lineage23_res[ko_abs_grp_test_lineage23_res$qval < .05, ]), #0to6
            # ko_abs_pseudotime_test_lineage23_res, 
            row.names(ko_branching_genes[ko_branching_genes$qval < .05, ]) )
table(andrew_sets_all)

#make venn diagram: 
golgi_abs_pseudotime_test_lineage23_res <- intersect(row.names(golgi_abs_pseudotime_test_lineage2_res[golgi_abs_pseudotime_test_lineage2_res$qval < 0.05, ]),
									 row.names(golgi_abs_pseudotime_test_lineage3_res[golgi_abs_pseudotime_test_lineage3_res$qval < 0.05, ]))

golgi_beam_element_all <- c(
            row.names(golgi_abs_grp_test_lineage23_res[golgi_abs_grp_test_lineage23_res$qval < .05, ]), #0to6
            # golgi_abs_pseudotime_test_lineage23_res, 
            row.names(golgi_branching_genes[golgi_branching_genes$qval < .05, ]) 
            )
golgi_beam_sets_all <- c(
         rep(paste('two group test', sep = ''), length(row.names(golgi_abs_grp_test_lineage23_res[golgi_abs_grp_test_lineage23_res$qval < .05, ]))),
         # rep(paste('intersection of pseudotime test', sep = ''), length(golgi_abs_pseudotime_test_lineage23_res)),
         rep(paste('BEAM test', sep = ''), length(row.names(golgi_branching_genes[golgi_branching_genes$qval < .05, ])))
         )

intersect(row.names(golgi_abs_grp_test_lineage23_res[golgi_abs_grp_test_lineage23_res$qval < .05, ]), #0to6
            # golgi_abs_pseudotime_test_lineage23_res, 
            row.names(golgi_branching_genes[golgi_branching_genes$qval < .05, ]) )
table(andrew_sets_all)

pdf(file = paste(fig_root_dir, 'golgi_beam_test_benchmark.pdf', sep = ''))
venneuler_venn(golgi_beam_element_all, golgi_beam_sets_all)
dev.off()

#union: 
#make venn diagram: 
abs_pseudotime_test_lineage23_res_gene_names <- union(row.names(abs_pseudotime_test_lineage2_res[abs_pseudotime_test_lineage2_res$qval < 0.05, ]),
									 row.names(abs_pseudotime_test_lineage3_res[abs_pseudotime_test_lineage3_res$qval < 0.05, ]))
beam_element_all <- c(#abs_pseudotime_test_lineage23_res_gene_names, 
            row.names(abs_grp_test_lineage23_res[abs_grp_test_lineage23_res$qval < .05, ]), 
            quake_branch_genes #  %in% valid_expressed_genes
            )
beam_sets_all <- c(
         # rep(paste('union of pseudotime test', sep = ''), length(abs_pseudotime_test_lineage23_res_gene_names)),
         rep(paste('two group test', sep = ''), sum(abs_grp_test_lineage23_res$qval < .05)),
         rep(paste('BEAM test', sep = ''), length(quake_branch_genes))
         )

pdf(file = paste(fig_root_dir, 'lung_beam_test_benchmark_union.pdf', sep = ''))
venneuler_venn(beam_element_all, beam_sets_all)
dev.off()

#make venn diagram: 
ko_abs_pseudotime_test_lineage23_res <- union(row.names(ko_abs_pseudotime_test_lineage2_res[ko_abs_pseudotime_test_lineage2_res$qval < 0.05, ]),
									 row.names(ko_abs_pseudotime_test_lineage3_res[ko_abs_pseudotime_test_lineage3_res$qval < 0.05, ]))
ko_beam_element_all <- c(
            row.names(ko_abs_grp_test_lineage23_res[ko_abs_grp_test_lineage23_res$qval < .05, ]), #0to6
            # ko_abs_pseudotime_test_lineage23_res, 
            row.names(ko_branching_genes[ko_branching_genes$qval < .05, ]) 
            )
ko_beam_sets_all <- c(
         rep(paste('two group test', sep = ''), length(row.names(ko_abs_grp_test_lineage23_res[ko_abs_grp_test_lineage23_res$qval < .05, ]))),
         # rep(paste('union of pseudotime test', sep = ''), length(ko_abs_pseudotime_test_lineage23_res)),
         rep(paste('BEAM test', sep = ''), length(row.names(ko_branching_genes[ko_branching_genes$qval < .05, ])))
         )

pdf(file = paste(fig_root_dir, 'ko_beam_test_benchmark_union.pdf', sep = ''))
venneuler_venn(ko_beam_element_all, ko_beam_sets_all)
dev.off()

#make venn diagram: 
golgi_abs_pseudotime_test_lineage23_res <- union(row.names(golgi_abs_pseudotime_test_lineage2_res[golgi_abs_pseudotime_test_lineage2_res$qval < 0.05, ]),
									 row.names(golgi_abs_pseudotime_test_lineage3_res[golgi_abs_pseudotime_test_lineage3_res$qval < 0.05, ]))

golgi_beam_element_all <- c(
            row.names(golgi_abs_grp_test_lineage23_res[golgi_abs_grp_test_lineage23_res$qval < .05, ]), #0to6
            # golgi_abs_pseudotime_test_lineage23_res, 
            row.names(golgi_branching_genes[golgi_branching_genes$qval < .05, ]) 
            )
golgi_beam_sets_all <- c(
         rep(paste('two group test', sep = ''), length(row.names(golgi_abs_grp_test_lineage23_res[golgi_abs_grp_test_lineage23_res$qval < .05, ]))),
         # rep(paste('union of pseudotime test', sep = ''), length(golgi_abs_pseudotime_test_lineage23_res)),
         rep(paste('BEAM test', sep = ''), length(row.names(golgi_branching_genes[golgi_branching_genes$qval < .05, ])))
         )

pdf(file = paste(fig_root_dir, 'golgi_beam_test_benchmark_union.pdf', sep = ''))
venneuler_venn(golgi_beam_element_all, golgi_beam_sets_all)
dev.off()

save.image('./RData/BEAM_test.RData')
