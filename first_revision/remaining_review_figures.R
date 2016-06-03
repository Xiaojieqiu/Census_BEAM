#this script will generate all the figures we have not included in other scripts yet for the review: 

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
