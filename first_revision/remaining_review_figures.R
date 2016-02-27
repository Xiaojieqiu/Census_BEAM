#this script will generate all the figures we have not done yet for the review: 

#####1. same heatmap#####
ph_res <- Shalek_abs_subset_ko_LPS_heatmap_annotations$ph_res  #Shalek_golgi_update_heatmap_annotations

heatmap_gene_list <- ph_res$tree_row$labels #get the label
lineage_states=c(2, 3)
cds_subset <- Shalek_abs_subset_ko_LPS[heatmap_gene_list, ]
new_cds <- buildLineageBranchCellDataSet(cds_subset, 
                                         lineage_states=c(2, 3), 
                                         branch_point=NULL,
                                         stretch = T)

LineageA_exprs <- exprs(new_cds[, pData(new_cds)$Lineage == as.numeric(unique(as.character(pData(new_cds)$Lineage))[1])])[, sort(pData(new_cds[, pData(new_cds)$Lineage == as.numeric(unique(as.character(pData(new_cds)$Lineage))[1])])$Pseudotime, index.return = T)$ix]
LineageB_exprs <- exprs(new_cds[, pData(new_cds)$Lineage == as.numeric(unique(as.character(pData(new_cds)$Lineage))[2])])[, sort(pData(new_cds[, pData(new_cds)$Lineage == as.numeric(unique(as.character(pData(new_cds)$Lineage))[2])])$Pseudotime, index.return = T)$ix]

col_gap_ind <- sum(pData(new_cds)$Lineage == as.numeric(unique(as.character(pData(new_cds)$Lineage))[1])) + 1 

newdataA <- data.frame(Pseudotime = sort(pData(new_cds[, pData(new_cds)$Lineage == as.numeric(unique(as.character(pData(new_cds)$Lineage))[1])])$Pseudotime))
newdataB <- data.frame(Pseudotime = sort(pData(new_cds[, pData(new_cds)$Lineage == as.numeric(unique(as.character(pData(new_cds)$Lineage))[2])])$Pseudotime))

#change to half of the number of the progenitor cell 
#common_ancestor_cells <- row.names(pData(new_cds)[duplicated(pData(new_cds)$original_cell_id),])
common_ancestor_cells <- row.names(pData(cds_subset)[pData(cds_subset)$State == setdiff(pData(cds_subset)$State, lineage_states),])

#account for the changes in the buildLineageBranchCellDataSet
  LineageP_num <- length(common_ancestor_cells) / 2
  LineageA_num <- sum(pData(cds_subset)$State == as.numeric(unique(as.character(pData(new_cds)$Lineage))[1]))
  LineageB_num <- sum(pData(cds_subset)$State == as.numeric(unique(as.character(pData(new_cds)$Lineage))[2])) + 1

  LineageA_exprs <- log10(LineageA_exprs + 1)
  LineageB_exprs <- log10(LineageB_exprs + 1)

heatmap_matrix <- cbind(LineageA_exprs[, (col_gap_ind - 1):1], LineageB_exprs)

row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1

scaling <- T

if(scaling) {
  heatmap_matrix <- Matrix::t(scale(Matrix::t(heatmap_matrix)))
  heatmap_matrix[heatmap_matrix > 3] <- 3
  heatmap_matrix[heatmap_matrix < -3] <- -3     
}
heatmap_matrix[heatmap_matrix > 3] <- 3
heatmap_matrix[heatmap_matrix < -3] <- -3     
heatmap_matrix_ori <- heatmap_matrix
heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[, 1]) & is.finite(heatmap_matrix[, col_gap_ind]), ] #remove the NA fitting failure genes for each lineage 

exp_rng <- range(heatmap_matrix) #bks is based on the expression range
bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by=0.1)
hmcols <- blue2green2red(length(bks) - 1)

ph_res <- pheatmap(heatmap_matrix[, ], #ph$tree_row$order
                   useRaster = T,
                   cluster_cols=FALSE, 
                   cluster_rows=F, 
                   show_rownames=F, 
                   show_colnames=F, 
                   #scale="row",
                   #clustering_distance_rows=row_dist, #row_dist
                   #clustering_method = "ward", #ward.D2
                   #cutree_rows=6,
                   # cutree_cols = 2,
#                    annotation_row=annotation_row,
#                    annotation_col=annotation_col,
#                    annotation_colors=annotation_colors,
                   gaps_col = col_gap_ind,
                   #gaps_row = row_gap_ind,

                   #treeheight_row = 1.5, 
                   breaks=bks,
                   fontsize = 6,
                   color=hmcols, 
                   height = 2.45, 
                   width = 3.4,
                   filename="./supplementary_figures/raw_expression_pheatmap_ko.pdf"
)

pdf('../BEAM/supplementary_figures/raw_expression_pheatmap_clustering.pdf')
test <- plot_genes_branched_heatmap(cds_subset, norm_method = 'log', use_fitting_curves = F)
dev.off()

#####2. bifurcation time for lung regulators ##### (see the branch time point analysis)
#this doesn't really work since the lung has less data at transition time point 

####3. mode estimate of log10(FPKM) distribution ####
cell_dmode_df <- data.frame(dmode = estimate_t(exprs(isoform_count_cds)[1:119469, ]), Time = pData(abs_AT12_cds_subset_all_gene)$Time)
pdf('./supplementary_figures/mode_distribution.pdf', width = 2.2, height = 1.4)
qplot(dmode, fill = Time, log = 'x', data = cell_dmode_df) + facet_wrap(~Time)  + xlab('Most frequent log10(FPKM)') + ylab('Cells') + nm_theme() #+ geom_vline(x = 1, linetype = 'longdash', color = I('blue'), size = .1)
dev.off()

#mode of fpkm values VS sequence depth: 
t_kb <- lapply(molModels_select, function(x) {
  coefs <- coef(x)
  10^(-coefs[1] / coefs[2])
}
)

cell_dmode_df$t_kb <- unlist(t_kb)
pdf('./supplementary_figures//true_mode_distribution.pdf', width = 2.2, height = 1.4)
qplot(t_kb, fill = Time, data = cell_dmode_df, log = 'x') + facet_wrap(~Time)  + xlab('t estimate calculated from regression') + ylab('Cells') + nm_theme() #+ geom_vline(x = 1, linetype = 'longdash', color = I('blue'), size = .1)
dev.off()

cell_dmode_df$total_counts <- apply(read_countdata, 2, sum)
pdf('./supplementary_figures/total_counts_vs_true_mode.pdf', width = 2.2, height = 1.4)
qplot(total_counts, t_kb, color = Time, data = cell_dmode_df, log = 'xy') + nm_theme()
dev.off()

cell_dmode_df$total_counts <- apply(read_countdata, 2, sum)
pdf('./supplementary_figures/total_counts_vs_estimated_mode.pdf', width = 2.2, height = 1.4)
qplot(total_counts, dmode, color = Time, data = cell_dmode_df, log = 'xy') + nm_theme()
dev.off()

cell_dmode_df$dmode_transcript_count <- estimate_t(absolute_cds[1:transcript_num, ])

####4. enrichment for GO/motif analysis for the muscle permuation test ####
abs_HSMM_ids <- row.names(HSMM_myo_size_norm_res[HSMM_myo_size_norm_res$pval <0.05, ]) 
std_HSMM_ids <- row.names(std_HSMM_myo_pseudotime_res_ori[std_HSMM_myo_pseudotime_res_ori$pval <0.05, ])

abs_HSMM_name <- fData(HSMM_myo[abs_HSMM_ids, ])$gene_short_name
std_HSMM_name <- fData(HSMM_myo[std_HSMM_ids, ])$gene_short_name

write.table(std_HSMM_name, file = '../std_muscle_pseudotime_gene.txt', sep = '\t', row.names = F, quote = F)
write.table(abs_HSMM_name, file = '../abs_muscle_pseudotime_gene.txt', sep = '\t', row.names = F, quote = F)

# abs_gsaRes <- runGSAhyper(unique(abs_HSMM_name), gsc=human_go_gsc, universe=unique(as.vector(fData(HSMM_myo)$gene_short_name)))
# std_gsaRes <- runGSAhyper(unique(std_HSMM_name), gsc=human_go_gsc, universe=unique(as.vector(fData(HSMM_myo)$gene_short_name)))

abs_gsaRes_reactome <- runGSAhyper(unique(abs_HSMM_name), gsc=human_reactome_gsc, universe=unique(as.vector(fData(HSMM_myo)$gene_short_name)))
std_gsaRes_reactome <- runGSAhyper(unique(std_HSMM_name), gsc=human_reactome_gsc, universe=unique(as.vector(fData(HSMM_myo)$gene_short_name)))

#plot the terms and significance: 
gsa_results <- list('Transcript_counts' = abs_gsaRes, 'FPKM' = std_gsaRes)
plot_gsa_hyper_heatmap(HSMM_myo, gsa_results, significance=1e-3)

gsa_results_reactome <- list('Transcript_counts' = abs_gsaRes_reactome, 'FPKM' = std_gsaRes_reactome)
plot_gsa_hyper_heatmap(HSMM_myo, gsa_results_reactome, significance=1e-3)

benchmark_pseudotime_test <- function(abs_gsaRes , std_gsaRes ) {
  abs_hyper_df_all <- data.frame(gene_set = names(abs_gsaRes$pvalues), pval = abs_gsaRes$pvalues, qval = abs_gsaRes$p.adj)
  abs_hyper_df_all$qval <- p.adjust(abs_hyper_df_all$pval, method = 'fdr')
  colnames(abs_hyper_df_all)[1] <- "cluster_id"
  
  std_hyper_df_all <- data.frame(gene_set = names(std_gsaRes$pvalues), pval = std_gsaRes$pvalues, qval = std_gsaRes$p.adj)
  std_hyper_df_all$qval <- p.adjust(std_hyper_df_all$pval, method = 'fdr')
  colnames(std_hyper_df_all)[1] <- "cluster_id"
  
  abs_hyper_df <- subset(abs_hyper_df_all, abs_hyper_df_all[, 'qval'] <= 0.01)
  std_hyper_df <- subset(std_hyper_df_all, std_hyper_df_all[, 'qval'] <= 0.01)
  
  gsa_results <- list('1' = abs_gsaRes, '2' = std_gsaRes)
  
  #performance comparision based on Cole's idea: 
  element_all_list <- c(
    abs_hyper_df$cluster_id, 
    std_hyper_df$cluster_id
  )
  
  sets_all <- c(
    rep(paste('transcript counts', sep = ''), nrow(abs_hyper_df)),
    rep(paste('FPKM values', sep = ''), nrow(std_hyper_df))
  )
  
  pdf(file = './nbt_2nd_sub_reviewers/overlap_enriched_muscle_term.pdf')
  venneuler_venn(element_all_list, sets_all)
  dev.off()
  
  #p-val plot show the gene names relevant for muscle differentiation 
  abs_hyper_df_all$cluster_id
  muscle_term_ids <- c(grep(pattern = 'MUSCLE', abs_hyper_df_all$cluster_id, ignore.case = T), grep(pattern = 'Myogenesis', abs_hyper_df_all$cluster_id, ignore.case = T))
  # muscle_term_ids <- Reduce(intersect, list(which(abs_hyper_df_all$qval < 0.01), which(std_hyper_df_all$qval < 0.01), muscle_term_ids))
  Types <- rep('Term without muscle', nrow(abs_hyper_df_all))
  Types[muscle_term_ids] <- 'Term with muscle'
  muscle_gs_df <- data.frame(cluster_id = abs_hyper_df_all$cluster_id, abs = abs_hyper_df_all$pval, std = std_hyper_df_all$pval, Type = Types) #muscle related/ non-muscle related 
  muscle_gs_df$ratio <- muscle_gs_df$std /  muscle_gs_df$abs
  pdf(file = './nbt_2nd_sub_reviewers/muscle_pseudotime_benchmark_qval.pdf')
    ggplot(data = muscle_gs_df[c(muscle_term_ids, which(muscle_gs_df$abs < 0.01 & muscle_gs_df$std < 0.01 )), ], aes(Type, log(ratio))) + 
        geom_boxplot(aes(fill = Type), alpha = 0.3, size = 0.5,  outlier.size = 0.5, lwd = 0.5, fatten = 0.5) + 
        geom_jitter() + nm_theme() + geom_vline(xintercept = 0) + xlab('log(P(FPKM) / P(transcript counts))') + ylab('') + scale_size(range = c(0.5, 0.5))
  dev.off()
  
  qplot(log(ratio), data = muscle_gs_df[c(muscle_term_ids, which(muscle_gs_df$abs < 0.01 & muscle_gs_df$std < 0.01 )), ], 
        fill = Type, geom = 'density', log = 'x', alpha = 0.3) + nm_theme() + geom_vline(xintercept = 0) + xlab('log(P(FPKM) / P(transcript counts))') + ylab('')
  #qplot(abs, std, data = muscle_gs_df[c(muscle_term_ids, which(muscle_gs_df$abs < 0.01 | muscle_gs_df$std < 0.01 )), ], color = Type, log = 'xy') + nm_theme() + xlab('transcript counts') + ylab('FPKM')
  dev.off()
}

benchmark_pseudotime_test(abs_gsaRes_reactome, std_gsaRes_reactome)
####5. overlap the genes from golgi-plug test VS those from the paper ####
#perform the dimension reduction with the new list of genes: 
golgi_order_genes_from_paper_df <- read.table('/Users/xqiu/Dropbox (Personal)/Projects/BEAM/study_gene_categories.txt', row.names = 1, header = T)
golgi_order_genes_from_paper <- row.names(subset(golgi_order_genes_from_paper_df, cluster %in% c('Id', 'IIIc', 'IIId'))) #core antiviral, peak inf, sustained inf
golgi_order_genes_from_paper <- capitalize(tolower(golgi_order_genes_from_paper))
golgi_order_genes_from_paper_id <- row.names(subset(fData(Shalek_golgi_update), gene_short_name %in% golgi_order_genes_from_paper))

Shalek_golgi_update_from_paper <- setOrderingFilter(Shalek_golgi_update, golgi_order_genes_from_paper_id)
Shalek_golgi_update_from_paper <- reduceDimension(Shalek_golgi_update_from_paper, use_vst = T, use_irlba=F, pseudo_expr = 0, residualModelFormulaStr = "~num_genes_expressed", scaling = F, method = "ICA")
Shalek_golgi_update_from_paper <- orderCells(Shalek_golgi_update_from_paper, num_path = 2)

shalek_custom_color_scale_plus_states= c(shalek_custom_color_scale, c('1'='#40A43A', '2'='#CB1B1E', '3'='#3660A5', 'Unstimulated_Replicate.' = 'gray'))

pdf(file = paste(fig_root_dir, 'fig6b_from_paper.pdf', sep = ''), height = 3, width = 3)
monocle::plot_spanning_tree(Shalek_golgi_update_from_paper, color_by="interaction(experiment_name, time)", cell_size=1) + 
  scale_color_manual(values=shalek_custom_color_scale_plus_states) + illustrator_theme() 
dev.off()

#overlap the genes: 
andrew_element_all <- c(
  golgi_order_genes_from_paper_id, #0to6
  row.names(subset(fData(Shalek_golgi_update), use_for_ordering == T))
)
andrew_sets_all <- c(
  rep(paste('Paper', sep = ''), length(golgi_order_genes_from_paper_id)),
  rep(paste('Order genes', sep = ''), length(row.names(subset(fData(Shalek_golgi_update), use_for_ordering == T))))
)

pdf(file = paste('./tmp/', 'overlap_shalek_golgi_plug_gene_order_gene.pdf', sep = ''))
venneuler_venn(andrew_element_all, andrew_sets_all)
dev.off()

####6. overlap the degs with strectching and not stretch####
save(no_stretch_weihgted_relative_abs_AT12_cds_subset_all_gene, no_stretch_ko_branching_genes, no_stretch_golgi_branching_genes, file = './tmp/no_stretch_branch_genes')
#lung
no_stretch_weihgted_relative_abs_AT12_cds_subset_all_gene <- branchTest(abs_AT12_cds_subset_all_gene[1:transcript_num, ], cores = detectCores(), relative_expr = T, weighted = T, stretch = F)
andrew_element_all <- c(
  row.names(weihgted_relative_abs_AT12_cds_subset_all_gene[weihgted_relative_abs_AT12_cds_subset_all_gene$qval < .05, ]), #0to6
  row.names(no_stretch_weihgted_relative_abs_AT12_cds_subset_all_gene[no_stretch_weihgted_relative_abs_AT12_cds_subset_all_gene$qval < .05, ])
)
andrew_sets_all <- c(
  rep(paste('Lung (stretch)', sep = ''), length(row.names(weihgted_relative_abs_AT12_cds_subset_all_gene[weihgted_relative_abs_AT12_cds_subset_all_gene$qval < .05, ]))),
  rep(paste('Lung (no stretch)', sep = ''), length(row.names(no_stretch_weihgted_relative_abs_AT12_cds_subset_all_gene[no_stretch_weihgted_relative_abs_AT12_cds_subset_all_gene$qval < .01, ])))
)

pdf(file = paste('./tmp/', 'lung_stretch_no_stretch_cmpr.pdf', sep = ''))
venneuler_venn(andrew_element_all, andrew_sets_all)
dev.off()

intersect(row.names(weihgted_relative_abs_AT12_cds_subset_all_gene[weihgted_relative_abs_AT12_cds_subset_all_gene$qval < .05, ]), #0to6
          row.names(no_stretch_weihgted_relative_abs_AT12_cds_subset_all_gene[no_stretch_weihgted_relative_abs_AT12_cds_subset_all_gene$qval < .05, ]))
table(andrew_sets_all)
#ko
no_stretch_ko_branching_genes = branchTest(Shalek_abs_subset_ko_LPS, fullModelFormulaStr = full_model_string, cores = detectCores(), relative_expr = T, weighted = T, stretch = F)
andrew_element_all <- c(
  row.names(ko_branching_genes[ko_branching_genes$qval < .05, ]), #0to6
  row.names(no_stretch_ko_branching_genes[no_stretch_ko_branching_genes$qval < .05, ])
)
andrew_sets_all <- c(
  rep(paste('ko (stretch)', sep = ''), length(row.names(ko_branching_genes[ko_branching_genes$qval < .05, ]))),
  rep(paste('ko (no stretch)', sep = ''), length(row.names(no_stretch_ko_branching_genes[no_stretch_ko_branching_genes$qval < .05, ])))
)

pdf(file = paste('./tmp/', 'ko_stretch_no_stretch_cmpr.pdf', sep = ''))
venneuler_venn(andrew_element_all, andrew_sets_all)
dev.off()

intersect(row.names(ko_branching_genes[ko_branching_genes$qval < .05, ]), #0to6
          row.names(no_stretch_ko_branching_genes[no_stretch_ko_branching_genes$qval < .05, ]))
table(andrew_sets_all)
#golgi-plug
no_stretch_golgi_branching_genes = branchTest(Shalek_golgi_update, fullModelFormulaStr = full_model_string, cores= detectCores(), relative_expr = T, weighted = T, stretch = F)
andrew_element_all <- c(
  row.names(golgi_branching_genes[golgi_branching_genes$qval < .05, ]), #0to6
  row.names(no_stretch_golgi_branching_genes[no_stretch_golgi_branching_genes$qval < .05, ])
)
andrew_sets_all <- c(
  rep(paste('golgi-plug (stretch)', sep = ''), length(row.names(golgi_branching_genes[golgi_branching_genes$qval < .05, ]))),
  rep(paste('golgi-plug (no stretch)', sep = ''), length(row.names(no_stretch_golgi_branching_genes[no_stretch_golgi_branching_genes$qval < .05, ])))
)

pdf(file = paste('./tmp/', 'golig_plug_stretch_no_stretch_cmpr.pdf', sep = ''))
venneuler_venn(andrew_element_all, andrew_sets_all)
dev.off()

intersect(row.names(golgi_branching_genes[golgi_branching_genes$qval < .05, ]), #0to6
          row.names(no_stretch_golgi_branching_genes[no_stretch_golgi_branching_genes$qval < .05, ]))
table(andrew_sets_all)
save(no_stretch_weihgted_relative_abs_AT12_cds_subset_all_gene, no_stretch_ko_branching_genes, no_stretch_golgi_branching_genes, file = '/tmp/branch_gene.RData')
####7. number number of detected genes affect the normalization approach ####

#function to calculate all the number used in the paper: 
#1. consensus between five tools: 
length(readcount_overlap)
#[1] 946
length(readcount_union)
#[1] 13364
length(abs_overlap)
#[1] 1649
length(abs_union)
#[1] 7209

#2. false discovery rate for the DESEq: 
fdr_sensitivity_df

#3. number of enriched TFs: 
dim(valid_hyper_df)

#number of branch genes
length(row.names(subset(weihgted_relative_abs_AT12_cds_subset_all_gene, qval <= 0.05)))
#[1] 1720
length(row.names(subset(ko_branching_genes, qval < 0.05)))
#[1] 1062

#number of genes related to interferon signaling: 
table(Shalek_abs_subset_ko_LPS_heatmap_annotations$annotation_row$Cluster)

#debug the permutation test: 
#the error comes from the change of the permutation algorithm. When the original approach used. we get good results 

#reproduce the previous result: 
#we got clear explannation on the changes of genes but also confirm our algorithm is correct: 
std_HSMM_myo_pseudotime_res_ori <- differentialGeneTest((std_HSMM[, ]), relative_expr = F, 
                                                        fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", #log10(Total_mRNAs) + spike_total_mRNAs
                                                        reducedModelFormulaStr = "~1", cores = detectCores())


HSMM_myo_size_norm_res <- differentialGeneTest(HSMM_myo[, ], 
                                               fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", #log10(Total_mRNAs) + spike_total_mRNAs
                                               reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

#update the supplementary data: 
lung_go_enrichment <- read.delim('./supplementary_data/lung_hyper_df.xls', header = T)
ko_reactome_enrichment <- read.delim('./supplementary_data/ko_hyper_df.xls', header = T)
golgiplug_reactome_enrichment <- read.delim('./supplementary_data/golgiplug_hyper_df.xls', header = T)

lung_go_enrichment$cluster_id <- as.character(lung_go_enrichment$cluster_id)
lung_go_enrichment$cluster_id <- revalue(lung_go_enrichment$cluster_id, c("1" = "cluster_1", "2" = "cluster_6", "3" = "cluster_5", "4" = "cluster_4", "5" = "cluster_2", "6" = "cluster_3"))

ko_reactome_enrichment$cluster_id <- as.character(ko_reactome_enrichment$cluster_id)
ko_reactome_enrichment$cluster_id <- revalue(ko_reactome_enrichment$cluster_id, c("1" = "cluster_5", "2" = "cluster_6", "3" = "cluster_2", "4" = "cluster_4", "5" = "cluster_1", "6" = "cluster_3"))

golgiplug_reactome_enrichment$cluster_id <- as.character(golgiplug_reactome_enrichment$cluster_id)
golgiplug_reactome_enrichment$cluster_id <- revalue(golgiplug_reactome_enrichment$cluster_id, c("1" = "cluster_5", "2" = "cluster_4", "3" = "cluster_6", "4" = "cluster_2", "5" = "cluster_3", "6" = "cluster_1"))

write.table(lung_go_enrichment, sep = '\t', quote = F, row.names = F, file = 'lung_go_enrichment')
write.table(ko_reactome_enrichment, sep = '\t', quote = F, row.names = F, file = 'ko_reactome_enrichment')
write.table(golgiplug_reactome_enrichment, sep = '\t', quote = F, row.names = F, file = 'golgiplug_reactome_enrichment')

#save it to file: 
permutation_res <- data.frame(spikein_free = mc_mode_size_norm_permutate_ratio_by_geometric_mean, spikein = mode_size_norm_permutate_ratio_by_geometric_mean, read_count = readcount_permutate_pval[names(std_permutate_pval)], fpkm = std_permutate_pval)
write.table(permutation_res, file = './permutation.txt', sep = '\t', row.names = T, quote = F)

