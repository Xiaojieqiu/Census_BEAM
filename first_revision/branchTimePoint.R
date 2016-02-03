#code for performing branching time point calculation

#Quake data: 
# 1. add branching time point for the casual genes
bifurcation_time  = abs(abs_bifurcation_time[fData(abs_AT12_cds_subset_all_gene[branch_motif_Tfs_id, ])$gene_short_name])

plot_genes_branched_pseudotime2(abs_AT12_cds_subset_all_gene[branch_motif_Tfs_id, ], cell_color_by = "State", bifurcation_time  = abs(bif_time), 
                trajectory_color_by = "Lineage", fullModelFormulaStr = '~sm.ns(Pseudotime, df = 3)*Lineage', normalize = F, stretch = T,
                lineage_labels = c('AT1', 'AT2'), cell_size = 1, ncol = 4, reducedModelFormulaStr = "~sm.ns(Pseudotime, df=3)", add_pval = T) + 
                ylab('Transcript counts') + nm_theme() + xlab('Pseudotime')

# # 2. GSEA on all the branching genes?
# get_gene_set_enrichment <- function(cds, branching_genes, branch_genes_abs_bifurcation_time, gene_set, cores=1) {
# 	# Normalize the gene names in GO enrichment set so cases match
# 	# for (name in names(gene_set$gsc)) {
# 	#     gene_set$gsc[[name]] = toupper(gene_set$gsc[[name]])
# 	# }

# 	# Clean up dataframe so no NA values
# 	branch_genes_abs_bifurcation_time = branch_genes_abs_bifurcation_time[!is.na(branch_genes_abs_bifurcation_time)]

# 	# Get the p-values from branch test for GSEA
# 	overlapping_branching_genes = branching_genes[as.character(branching_genes[, 'gene_short_name']) %in% names(branch_genes_abs_bifurcation_time), ]
# 	pval = overlapping_branching_genes$pval
# 	names(pval) = overlapping_branching_genes$gene_short_name
# 	pval = pval[!is.na(pval)]
# 	branch_genes_abs_bifurcation_time = branch_genes_abs_bifurcation_time[names(pval)]

# 	# Now run the GSEA using the GO enrichment terms
# 	gsaRes_go <- runGSA(pval, branch_genes_abs_bifurcation_time, gsc = gene_set, ncpus = cores)

# 	return(gsaRes_go)
# }

branch_genes_abs_bifurcation_time <- abs_bifurcation_time[as.character(fData(abs_AT12_cds_subset_all_gene[quake_branch_genes[quake_branch_genes %in% valid_expressed_genes], ])$gene_short_name)]
qval_df <- weihgted_relative_abs_AT12_cds_subset_all_gene[quake_branch_genes[quake_branch_genes %in% valid_expressed_genes], ]

# gsaRes_go <- get_gene_set_enrichment(abs_AT12_cds_subset_all_gene, 
# 	qval_df, 
# 	branch_genes_abs_bifurcation_time, mouse_go_gsc, cores = 50)

#use branching time point as the geneLevelStats: 
branch_genes_abs_bifurcation_time <- branch_genes_abs_bifurcation_time[!is.na(branch_genes_abs_bifurcation_time)]
# gsaRes_go <- runGSA(branch_genes_abs_bifurcation_time, gsc = mouse_go_gsc, ncpus = 50) #floor(detectCores() / 2) * 2

# pdf(paste(submission_directory, "submission_fig4b_go_enrichment.pdf", sep = ''), height=100, width=15)
# plot_gsa_hyper_heatmap(abs_AT12_cds_subset_all_gene, branchGenes_gsa_results_branch_genes, significance = 1e-2)
# dev.off()

# #make the heatmap plot for the enrichment: 
# plot_gene_set_enrichment <- function (gsa_res, gene_num = 20, selection = F,
#     custom_p_adjust = F, add_terms = F)
# {
# 	# Extract gene set enrichment data for more genes than needed (ensures will get equal AT1 and AT2 numbers)
#     go_data <- make_enrichment_df(gsa_res, 1000000, custom_p_adjust, add_terms)
#     go_data <- unique(go_data[, ])

#     # Get only the significant values and assign to lineages based on adjusted pvalues
#     AT1_values = go_data[go_data$"p adj (dist.dir.up)" < 0.05, ]
#     AT1_values$lineage = 'Lineage 2'

#     AT2_values = go_data[go_data$"p adj (dist.dir.dn)" < 0.05, ]
#     AT2_values$lineage = 'Lineage 3'

#     # Get the top significant values by stat
#     top_AT1 = tail(AT1_values, gene_num)
#     top_AT2 = head(AT2_values, gene_num)

#     top_go_data = rbind(top_AT2, top_AT1)

#     print(AT1_values$stat)
#     print(AT2_values$stat)

#     gsea_plot <- qplot(x = 1:nrow(top_go_data), y = abs(as.numeric(stat)), data = top_go_data, geom = "bar", stat = "identity", fill = lineage) +
#         coord_flip() + 
#         scale_x_discrete(limits = 1:nrow(top_go_data), labels = str_split_fixed(top_go_data$Name, "%", 2)[, 1]) +
#         xlab("") + 
#         ylab("Normalized Enrichment Score") + 
#         scale_fill_discrete(name = "Type", labels = c("Lineage 2", "Lineage 3"))

#     return(gsea_plot)
# }

#Shalek data: 
#
weighted_abs_Shalek_KO_subset_all_gene_ILRs_list <- calILRs(cds = Shalek_abs_subset_ko_LPS[, ], lineage_states = c(2, 3), stretch = T, cores = 1, 
    trend_formula = "~sm.ns(Pseudotime, df = 3) * Lineage", ILRs_limit = 3, 
    relative_expr = T, weighted = T, label_by_short_name = F, 
    useVST = T, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "weighted_abs_Shalek_KO_subset_all_gene_ILRs_list", return_all = T)

#find the MST branching point: 
test <- buildLineageBranchCellDataSet(Shalek_abs_subset_ko_LPS, stretch=T)
range(pData(test)[488:688, 'Pseudotime'])
subset_MST_branch_abs_ko_bifurcation_time <- detectBifurcationPoint(weighted_abs_Shalek_KO_subset_all_gene_ILRs_list$norm_str_logfc_df[, ])
all_abs_ko_bifurcation_time <- detectBifurcationPoint(weighted_abs_Shalek_KO_subset_all_gene_ILRs_list$norm_str_logfc_df[, 44:100])


# weighted_abs_Shalek_Golgi_subset_all_gene_ILRs_list <- calILRs(cds = Shalek_golgi_update[, ], lineage_states = c(2, 3), stretch = T, cores = 1, 
# 	trend_formula = "~sm.ns(Pseudotime, df = 3) * Lineage", ILRs_limit = 3, 
# 	relative_expr = T, weighted = T, label_by_short_name = F, 
# 	useVST = T, round_exprs = T, pseudocount = 0, output_type = "all", file = "weighted_abs_Shalek_Golgi_subset_all_gene_ILRs_list", return_all = T)

# #find the MST branching point: 
# test <- buildLineageBranchCellDataSet(Shalek_golgi_update, stretch=T)
# range(pData(test)[496:641, 'Pseudotime'])
# subset_MST_branch_abs_Golgi_bifurcation_time <- detectBifurcationPoint(weighted_abs_Shalek_Golgi_subset_all_gene_ILRs_list$norm_str_logfc_df[, 33:100])
# all_abs_Golgi_bifurcation_time <- detectBifurcationPoint(weighted_abs_Shalek_Golgi_subset_all_gene_ILRs_list$norm_str_logfc_df[, ])

# #find the MST branching point: 
# test <- buildLineageBranchCellDataSet(Shalek_golgi_update, stretch=T)
# range(pData(test)[496:641, 'Pseudotime'])
# subset_MST_branch_abs_Golgi_bifurcation_time <- detectBifurcationPoint(weighted_abs_Shalek_Golgi_subset_all_gene_ILRs_list_df4$norm_str_logfc_df[, 33:100])
# all_abs_Golgi_bifurcation_time <- detectBifurcationPoint(weighted_abs_Shalek_Golgi_subset_all_gene_ILRs_list_df4$norm_str_logfc_df[, ])

#collect list of grn from Aviv science paper: 
grn_fig4a <- c("STAT2", "STAT1", "HHEX", "RBL1", "Timeless", "Ifnb1", "Cxcl9", "Ifit3", "Tsc22d1", "Tnfsf8", "Isg20", "Nmi", "Iigp1", "Irf7", "Lhx2", "Bbx", "Fus", "Tcf4", "Pml", "Usp12", "Irf7", "Ifit1", "Isg15", "Rbl1", "Usp25", "Daxx", "Cd40", "Atm", "Lrf1", "Ligp2", "Mertk", "Cxcl11", "Trim12", "Trim21", "NfkbIz")
grn_fig4b <- c("NFKBIZ", "ATF4", "ATF4", "IIf2b", "FUS", "RBL1", "RBL1", "ligp1", "CBX4", "DNMT3A", "DNMT3A", "Ifnb1", "NFKBIZ", "FOS", "FOS", "II12b", "JUN", "FUS", "FUS", "IFNB1", "FUS", "NFkbIz", "NFKBIZ", "FOS", "CBX4", "STAT1", "STAT1", "IFnb1", "CEBP2", "HHEX", "HHEX", "IL12b")
grn_fig4c <- c("Tlr3", "STAT1", "Tlr4", "STAT2", "Tlr2", "IRf8", "ETV8", "JUN", "STAT4", "IRF9", "ATF3", "RBL1", "PLAGL2", "NFKB1", "NFKB2", "RUX1", "CEBPB", "ATF4", "HAT1", "RELA", "PNRC2", "CXCL10", "IFNb1", "IL15", "HMGN3", "FUS", "SAP30", "IL12b", "CXCL2", "CXCL1", "IL1B", "CCL3", "NFE2L2", "IRF4", "FOS")

grn_all_genes <- c(grn_fig4a, grn_fig4b, grn_fig4c) 

grn_ids <- row.names(subset(fData(Shalek_abs_subset_ko_LPS), toupper(gene_short_name) %in% toupper(grn_fig4a)))
fData(Shalek_golgi_update)[grn_ids, 'gene_short_name']
grn_branch_time <- subset_MST_branch_abs_ko_bifurcation_time[grn_ids]
names(grn_branch_time) <- fData(Shalek_abs_subset_ko_LPS)[grn_ids, 'gene_short_name']
sort(abs(grn_branch_time), index.return = T)

# grn_branch_time <- all_abs_ko_bifurcation_time[grn_ids]
# names(grn_branch_time) <- fData(Shalek_abs_subset_ko_LPS)[grn_ids, 'gene_short_name']
# sort(abs(grn_branch_time), index.return = T)

# grn_branch_time <- subset_MST_branch_abs_Golgi_bifurcation_time[grn_ids]
# names(grn_branch_time) <- fData(Shalek_abs_subset_ko_LPS)[grn_ids, 'gene_short_name']
# sort(abs(grn_branch_time), index.return = T)

# grn_branch_time <- all_abs_Golgi_bifurcation_time[grn_ids]
# names(grn_branch_time) <- fData(Shalek_abs_subset_ko_LPS)[grn_ids, 'gene_short_name']
# sort(abs(grn_branch_time), index.return = T)

#make the new branching time plot with the network structure: 

#two group comparison: 
top_group <- c("STAT2", "STAT1", "HHEX", "RBL1", "Timeless", "FUS")
bottom_group <- c("Ifnb1", "Cxcl9", "Ifit3", "Tsc22d1", "Tnfsf8", "Isg20", "Nmi", "Iigp1", "Irf7", "Lhx2", "Bbx", "Fus", "Tcf4", "Pml", "Usp12", "Irf7", "Ifit1", "Isg15", "Usp25", "Daxx", "Cd40", "Atm", "Lrf1", "Lgp2", "Mertk", "Cxcl11", "Trim12", "Trim21")

top_group_ids <- row.names(subset(fData(Shalek_abs_subset_ko_LPS), toupper(gene_short_name) %in% toupper(top_group)))
fData(Shalek_golgi_update)[grn_ids, 'gene_short_name']
bottom_group_ids <- row.names(subset(fData(Shalek_abs_subset_ko_LPS), toupper(gene_short_name) %in% toupper(bottom_group)))
fData(Shalek_abs_subset_ko_LPS)[bottom_group_ids, 'gene_short_name']

#select only the significant genes: 
valid_ko_branching_genes_list <- intersect(row.names(subset(ko_branching_genes, qval < 0.05)), ko_valid_expressed_genes)
valid_top_group_ids <- intersect(top_group_ids, valid_ko_branching_genes_list)
top_group_branch_time <- subset_MST_branch_abs_ko_bifurcation_time[valid_top_group_ids]
names(top_group_branch_time) <- fData(Shalek_golgi_update)[valid_top_group_ids, 'gene_short_name']
sort(abs(top_group_branch_time), index.return = T)

valid_bottom_group_ids <- intersect(bottom_group_ids, valid_ko_branching_genes_list)
bottom_group_branch_time <- subset_MST_branch_abs_ko_bifurcation_time[valid_bottom_group_ids]
names(bottom_group_branch_time) <- fData(Shalek_abs_subset_ko_LPS)[valid_bottom_group_ids, 'gene_short_name']
sort(abs(bottom_group_branch_time), index.return = T)

df <- data.frame(bifurcation_time = c(top_group_branch_time, bottom_group_branch_time), 
		type = c(rep("Antiviral regulators", length(top_group_branch_time)), rep("Targets", length(bottom_group_branch_time))))

pdf('tmp/DC_regulation_hierarchy.pdf', width = 2.25, height = 2.25)
qplot(type, abs(bifurcation_time), color = type, geom = c('jitter', 'boxplot'), data = df, alpha = I(0.7), log = 'y') + 
    xlab('') + ylab('bifurcation time point') + #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`)) 
    nm_theme()
dev.off()

#figure 4c: 
grn_fig4c <- c("Tlr3", "STAT1", "Tlr4", "STAT2", "Tlr2", "IRf8", "ETV8", "JUN", "STAT4", "IRF9", "ATF3", "RBL1", "PLAGL2", "NFKB1", "NFKB2", "RUX1", "CEBPB", "ATF4", "HAT1", "RELA", "PNRC2", "CXCL10", "IFNb1", "IL15", "HMGN3", "FUS", "SAP30", "IL12b", "CXCL2", "CXCL1", "IL1B", "CCL3", "NFE2L2", "IRF4", "FOS")
fig4c_layer1 <- c("Tlr2", "Tlr3", "Tlr4")
fig4c_layer2.1 <- c("STAT1", "STAT2", "IRF8", "ETV6", "JUN", "STAT4", "IRF9", "ATF3", "RBL1")
fig4c_layer2.2 <- c("PLAGL2", "NFKB1", "RUNX1", "CEBPB", "ATF4", "HAT1", "RELA", "PNRC2")
fig4c_layer3.1 <- c("HMGN3", "FUS", "SAP30")
fig4c_layer3.2 <- c("NFE2L2", "IRF4", "FOS")
fig4c_layer4.1 <- c("CXCL10", "Ifnb1", "IL15")
fig4c_layer4.2 <- c("IL12b", "CXCL2", "CXCL1", "IL1B", "CCL3")

# #find ids for genes at each layers
# fig4c_layer1_ids <- row.names(subset(fData(Shalek_abs_subset_ko_LPS), toupper(gene_short_name) %in% toupper(fig4c_layer1)))
# fData(Shalek_golgi_update)[fig4c_layer1_ids, 'gene_short_name']
# fig4c_layer2.1_ids <- row.names(subset(fData(Shalek_abs_subset_ko_LPS), toupper(gene_short_name) %in% toupper(fig4c_layer2.1)))
# fData(Shalek_golgi_update)[fig4c_layer2.1_ids, 'gene_short_name']
# fig4c_layer2.2_ids <- row.names(subset(fData(Shalek_abs_subset_ko_LPS), toupper(gene_short_name) %in% toupper(fig4c_layer2.2)))
# fData(Shalek_golgi_update)[fig4c_layer2.2_ids, 'gene_short_name']
# fig4c_layer3.1_ids <- row.names(subset(fData(Shalek_abs_subset_ko_LPS), toupper(gene_short_name) %in% toupper(fig4c_layer3.1)))
# fData(Shalek_golgi_update)[fig4c_layer3.1_ids, 'gene_short_name']
# fig4c_layer3.2_ids <- row.names(subset(fData(Shalek_abs_subset_ko_LPS), toupper(gene_short_name) %in% toupper(fig4c_layer3.2)))
# fData(Shalek_golgi_update)[fig4c_layer3.2_ids, 'gene_short_name']
# fig4c_layer4.1_ids <- row.names(subset(fData(Shalek_abs_subset_ko_LPS), toupper(gene_short_name) %in% toupper(fig4c_layer4.1)))
# fData(Shalek_golgi_update)[fig4c_layer4.1_ids, 'gene_short_name']
# fig4c_layer4.2_ids <- row.names(subset(fData(Shalek_abs_subset_ko_LPS), toupper(gene_short_name) %in% toupper(fig4c_layer4.2)))
# fData(Shalek_golgi_update)[fig4c_layer4.2_ids, 'gene_short_name']

# valid_ko_branching_genes_list2 <- row.names(subset(ko_branching_genes, qval < 0.1))

# valid_fig4c_layer1_ids<- intersect(fig4c_layer1_ids, valid_ko_branching_genes_list2)
# valid_fig4c_layer2.1_ids <- intersect(fig4c_layer2.1_ids, valid_ko_branching_genes_list2)
# valid_fig4c_layer2.2_ids <- intersect(fig4c_layer2.2_ids, valid_ko_branching_genes_list2)
# valid_fig4c_layer3.1_ids <- intersect(fig4c_layer3.1_ids, valid_ko_branching_genes_list2)
# valid_fig4c_layer3.2_ids <- intersect(fig4c_layer3.2_ids, valid_ko_branching_genes_list2)
# valid_fig4c_layer4.1_ids <- intersect(fig4c_layer4.1_ids, valid_ko_branching_genes_list2)
# valid_fig4c_layer4.2_ids <- intersect(fig4c_layer4.2_ids, valid_ko_branching_genes_list2)

# fig4c_layer1_ids_branching_time <- subset_MST_branch_abs_ko_bifurcation_time[valid_fig4c_layer1_ids]
# names(fig4c_layer1_ids_branching_time) <- fData(Shalek_golgi_update)[valid_fig4c_layer1_ids, 'gene_short_name']
# sort(abs(fig4c_layer1_ids_branching_time), index.return = T)

# fig4c_layer2.1_ids_branching_time <- subset_MST_branch_abs_ko_bifurcation_time[valid_fig4c_layer2.1_ids]
# names(fig4c_layer2.1_ids_branching_time) <- fData(Shalek_golgi_update)[valid_fig4c_layer2.1_ids, 'gene_short_name']
# sort(abs(fig4c_layer2.1_ids_branching_time), index.return = T)

# fig4c_layer2.2_ids_branching_time <- subset_MST_branch_abs_ko_bifurcation_time[valid_fig4c_layer2.2_ids]
# names(fig4c_layer2.2_ids_branching_time) <- fData(Shalek_golgi_update)[valid_fig4c_layer2.2_ids, 'gene_short_name']
# sort(abs(fig4c_layer2.2_ids_branching_time ), index.return = T)

# fig4c_layer3.1_ids_branching_time <- subset_MST_branch_abs_ko_bifurcation_time[valid_fig4c_layer3.1_ids]
# names(fig4c_layer3.1_ids_branching_time ) <- fData(Shalek_golgi_update)[valid_fig4c_layer3.1_ids, 'gene_short_name']
# sort(abs(fig4c_layer3.1_ids_branching_time ), index.return = T)

# fig4c_layer3.2_ids_branching_time <- subset_MST_branch_abs_ko_bifurcation_time[valid_fig4c_layer3.2_ids]
# names(fig4c_layer3.2_ids_branching_time ) <- fData(Shalek_golgi_update)[valid_fig4c_layer3.2_ids, 'gene_short_name']
# sort(abs(fig4c_layer3.2_ids_branching_time ), index.return = T)

# fig4c_layer4.1_ids_branching_time <- subset_MST_branch_abs_ko_bifurcation_time[valid_fig4c_layer4.1_ids]
# names(fig4c_layer4.1_ids_branching_time ) <- fData(Shalek_golgi_update)[valid_fig4c_layer4.1_ids, 'gene_short_name']
# sort(abs(fig4c_layer4.1_ids_branching_time ), index.return = T)

# fig4c_layer4.2_ids_branching_time <- subset_MST_branch_abs_ko_bifurcation_time[valid_fig4c_layer4.2_ids]
# names(fig4c_layer4.2_ids_branching_time ) <- fData(Shalek_golgi_update)[valid_fig4c_layer4.2_ids, 'gene_short_name']
# sort(abs(fig4c_layer4.2_ids_branching_time ), index.return = T)

# #plot the bifurcation time for each genes at each layers: 
# df <- data.frame(bifurcation_time = c(fig4c_layer1_ids_branching_time, fig4c_layer2.1_ids_branching_time, fig4c_layer2.2_ids_branching_time, fig4c_layer3.1_ids_branching_time, fig4c_layer3.2_ids_branching_time, fig4c_layer4.1_ids_branching_time, fig4c_layer4.2_ids_branching_time), 
# 		type = c(rep("Layer1", length(fig4c_layer1_ids_branching_time)), rep("Layer2.1", length(fig4c_layer2.1_ids_branching_time)), rep("Layer2.2", length(fig4c_layer2.2_ids_branching_time)), 
# 			rep("Layer3.1", length(fig4c_layer3.1_ids_branching_time)), rep("Layer3.2", length(fig4c_layer3.2_ids_branching_time)), rep("Layer4.1", length(fig4c_layer4.1_ids_branching_time)), rep("Layer4.2", length(fig4c_layer4.2_ids_branching_time))))
# df$gene_short_name <- row.names(df)
# pdf('tmp/valid_fig4c_DC_regulation_hierarchy.pdf', width = 2.25, height = 2.25)
# qplot(type, abs(bifurcation_time), color = type, geom = c('jitter', 'boxplot'), data = df, alpha = I(0.7), log = 'y') + 
#     xlab('') + ylab('bifurcation time point') + #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`)) 
#     nm_theme() + coord_flip() + geom_text(aes(label = df$gene_short_name), size = 2, angle = 45, color = 'black')
# dev.off()

#enriched TFs and the corresponding targets from Lung data: 
#branching time between TFs and their targets: 
branch_motif_Tfs_id_time <- all_abs_bifurcation_time[branch_motif_Tfs_id]

#find the targets of each enriched branching TF: 
enrich_branching_TF_gsc_names <- names(TF_5k_enrichment_gsc$gsc)[which(names(TF_5k_enrichment_gsc$gsc) %in% toupper(fData(abs_AT12_cds_subset_all_gene[branch_motif_Tfs_id, ])$gene_short_name))]
valid_quake_branch_genes_names <- as.character(fData(absolute_cds[valid_quake_branch_genes, ])$gene_short_name)

#intersect the enriched TFs' targets with ALL significant branch genes 
branching_TF_TF_5k_enrichment_gsc <- lapply(enrich_branching_TF_gsc_names, function(x) intersect(TF_5k_enrichment_gsc$gsc[[x]], valid_quake_branch_genes_names))

#targets of genes for each enriched branching TF at the CORRESPONDING enriched clusters 
significant_branching_TF_TF_5k_enrichment_gsc <- lapply(enrich_branching_TF_gsc_names, function(x) {
	if(toupper(x) %in% c('TCF7L2', 'MAFF')){
		res <- intersect(TF_5k_enrichment_gsc$gsc[[x]], names(clusters[clusters == 6]))
	}
	else{	
		res <- intersect(TF_5k_enrichment_gsc$gsc[[x]], names(clusters[clusters %in% c(4, 6)]))
	}
	data.frame(regulators = x, targets = res, branching_time = all_abs_bifurcation_time_gene_names[res])
})

#obtain the branching time for the genes: 
all_abs_bifurcation_time_gene_names <- all_abs_bifurcation_time
names(all_abs_bifurcation_time_gene_names) <- as.character(fData(absolute_cds[names(all_abs_bifurcation_time), ])$gene_short_name)

branching_TF_TF_5k_enrichment_gsc_branch_time <- all_abs_bifurcation_time_gene_names[unique(unlist(branching_TF_TF_5k_enrichment_gsc))]

significant_branching_TF_TF_5k_enrichment_gsc_branch_time <- all_abs_bifurcation_time_gene_names[unique(unlist(significant_branching_TF_TF_5k_enrichment_gsc))]

# df <- data.frame(bifurcation_time = c(branch_motif_Tfs_id_time, significant_branching_TF_TF_5k_enrichment_gsc_branch_time), 
# 		type = c(rep("enriched branching TFs", length(branch_motif_Tfs_id_time)), rep("targets", length(branching_TF_TF_5k_enrichment_gsc_branch_time))))

# df <- do.call(rbind.data.frame, significant_branching_TF_TF_5k_enrichment_gsc)

# pdf('tmp/lung_regulation_hierarchy.pdf', width = 3, height = 3)
# qplot(regulators, abs(branching_time), color = regulators, geom = c('jitter', 'boxplot'), data = subset(df, abs(branching_time) > 10), alpha = I(0.7), log = 'y') + 
#     ylab('bifurcation time point') + #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`)) 
#     nm_theme()
# dev.off()

# pdf('tmp/lung_regulation_hierarchy.pdf', width = 3, height = 3)
# qplot(regulators, abs(branching_time), color = regulators, geom = c('jitter', 'boxplot'), data = df, alpha = I(0.7), log = 'y') + 
#     ylab('bifurcation time point') + #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`)) 
#     nm_theme()
# dev.off()

# df$regulators_branching_time <- all_abs_bifurcation_time_gene_names[df$regulators]
# pdf('tmp/lung_regulator_pseudotime.pdf', width = 3, height = 3)
# qplot(regulators, abs(regulators_branching_time), fill = regulators, geom = 'bar', stat = 'identity', data = unique(df[, c(1, 4)]), alpha = I(0.7), log = 'y') + 
#     ylab('bifurcation time point') + #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`)) 
#     nm_theme()
# dev.off()

#Shalek hierarchy: 
#add the branch time points for three TFs for the barplot and their targets: 
#add the genes from different layers: 
annotated_TFs <- c('Stat1', 'Stat2', 'IRF1', 'IRF2', 'PRDM1') #Stat1 Stat2 IRF1 IRF2 PRDM1
annotated_TF_ids <- row.names(subset(fData(Shalek_abs_subset_ko_LPS), toupper(gene_short_name) %in% toupper(annotated_TFs)))

enrich_branching_TF_gsc_names <- names(TF_5k_enrichment_gsc$gsc)[names(TF_5k_enrichment_gsc$gsc) %in% toupper(annotated_TFs)]
enrich_branching_TF_gsc_names <- c(enrich_branching_TF_gsc_names, 'STAT2::STAT1')

ko_branch_gene_ids <- row.names(ko_branching_genes)[ko_branching_genes$qval < 0.05]

valid_ko_branch_gene_ids <- ko_branch_gene_ids[ko_branch_gene_ids %in% ko_valid_expressed_genes]
valid_ko_branch_gene_names <- as.character(fData(Shalek_abs_subset_ko_LPS[valid_ko_branch_gene_ids, ])$gene_short_name)
branching_TF_TF_5k_enrichment_gsc <- lapply(enrich_branching_TF_gsc_names, function(x) intersect(toupper(TF_5k_enrichment_gsc$gsc[[x]]), toupper(valid_ko_branch_gene_names)))

#targets of genes for each enriched branching TF at the CORRESPONDING enriched clusters 
#calculate the 
significant_branching_TF_TF_5k_enrichment_gsc <- lapply(enrich_branching_TF_gsc_names, function(x) {
	res <- intersect(toupper(TF_5k_enrichment_gsc$gsc[[x]]), 
		toupper(names(Shalek_abs_subset_ko_LPS_tree_heatmap_clusters[Shalek_abs_subset_ko_LPS_tree_heatmap_clusters %in% c(5)])))
	data.frame(regulators = x, targets = res, branching_time = all_abs_bifurcation_time_gene_names[res])
})

df <- do.call(rbind.data.frame, significant_branching_TF_TF_5k_enrichment_gsc)
branching_tf_df <- data.frame(regulators = enrich_branching_TF_gsc_names, 
								all_abs_bifurcation_time_gene_names[enrich_branching_TF_gsc_names]))

pdf('tmp/ko_regulation_hierarchy.pdf', width = 3, height = 3)
qplot(regulators, abs(branching_time), color = regulators, geom = c('jitter', 'boxplot'), 
	data = subset(df, abs(branching_time) > 10), alpha = I(0.5), log = 'y') + 
    ylab('bifurcation time point') + 
    geom_point(aes(regulators, abs(branching_time), color = 'black'), data = branching_tf_df) + 
    #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`)) 
    nm_theme()
dev.off()












