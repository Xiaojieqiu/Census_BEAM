load('./RData/analysis_shalek_data.RData')
load('./RData/gen_lung_figures.RData')
library(monocle)
library(xacHelper)
library(grid)

#code for performing branching time point calculation

#Quake data: 
# 1. add branching time point for the casual genes
bifurcation_time  = abs(all_abs_bifurcation_time[branch_motif_Tfs_id])

names(bifurcation_time) <- fData(abs_AT12_cds_subset_all_gene[branch_motif_Tfs_id, ])$gene_short_name
# pdf('./main_figures/fig4d_branch_time.pdf', height = 4, width = 6)
# plot_genes_branched_pseudotime2(abs_AT12_cds_subset_all_gene[branch_motif_Tfs_id, ], color_by = "State", trajectory_color_by = 'Lineage', #panel_order = names(bifurcation_time),
#   fullModelFormulaStr = '~sm.ns(Pseudotime, df = 3)*Lineage', reducedModelFormulaStr = '~sm.ns(Pseudotime, df = 3)', normalize = T, stretch = T,
#   lineage_labels = c('AT1', 'AT2'), cell_size = 1, ncol = 2, bifurcation_time  = abs(bifurcation_time)) + nm_theme()+ ylab('Transcript counts') + xlab('Pseudotime')
# dev.off()

#Shalek data: 
#
weighted_abs_Shalek_KO_subset_all_gene_ILRs_list <- calILRs(cds = Shalek_abs_subset_ko_LPS[, ], lineage_states = c(2, 3), stretch = T, cores = detectCores(), 
    trend_formula = "~sm.ns(Pseudotime, df = 3) * Lineage", ILRs_limit = 3, 
    relative_expr = T, weighted = T, label_by_short_name = F, 
    useVST = T, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "weighted_abs_Shalek_KO_subset_all_gene_ILRs_list", return_all = T)

duplicate_progenitors_weighted_abs_Shalek_KO_subset_all_gene_ILRs_list <- calILRs(cds = Shalek_abs_subset_ko_LPS[, ], lineage_states = c(2, 3), stretch = T, cores = detectCores(), 
    trend_formula = "~sm.ns(Pseudotime, df = 3) * Lineage", ILRs_limit = 3, progenitor_method = 'duplicate', 
    relative_expr = T, weighted = T, label_by_short_name = F, 
    useVST = T, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "duplicate_progenitors_weighted_abs_Shalek_KO_subset_all_gene_ILRs_list", return_all = T)

#find the MST branching point: 
new_cds <- buildLineageBranchCellDataSet(Shalek_abs_subset_ko_LPS, stretch=T)

MST_branch_range <- range(pData(Shalek_abs_subset_ko_LPS[, pData(Shalek_abs_subset_ko_LPS)$State == 1])[, 'Pseudotime'])
test <- (pData(Shalek_abs_subset_ko_LPS)$Pseudotime - MST_branch_range[2])
first_branching_cell <- colnames(Shalek_abs_subset_ko_LPS)[which(abs(test) == min(abs(test)))]
MST_branch_time <- round(pData(new_cds)[first_branching_cell, 'Pseudotime'])
subset_MST_branch_abs_ko_bifurcation_time <- detectBifurcationPoint(duplicate_progenitors_weighted_abs_Shalek_KO_subset_all_gene_ILRs_list$norm_str_logfc_df[, ]) #weighted_abs_Shalek_KO_subset_all_gene_ILRs_list
all_abs_ko_bifurcation_time <- detectBifurcationPoint(duplicate_progenitors_weighted_abs_Shalek_KO_subset_all_gene_ILRs_list$norm_str_logfc_df[, MST_branch_time:100]) #weighted_abs_Shalek_KO_subset_all_gene_ILRs_list

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
#function to create a data frame storing all the branching time for the regulators and their targets:  
gen_branchTime_df <- function(cds = Shalek_golgi_update, branch_time = all_abs_ko_bifurcation_time, gene_clusters  = Shalek_abs_subset_ko_LPS_tree_heatmap_clusters, clusters_id = 5, 
	enrich_branching_TF_gsc_names) {
	all_abs_bifurcation_time_gene_names <- branch_time
	names(all_abs_bifurcation_time_gene_names) <- toupper(fData(cds[names(branch_time), ])$gene_short_name)
	significant_branching_TF_TF_5k_enrichment_gsc <- lapply(enrich_branching_TF_gsc_names, function(x) {
		res <- intersect(toupper(TF_5k_enrichment_gsc$gsc[[x]]), 
			toupper(names(gene_clusters[gene_clusters %in% clusters_id])))
		message('res is', res)
		regulator_time <- NA
		if(length(res) <= 1)
			return(data.frame(regulators = x, targets = NA, branching_time = NA, regulators_time = regulator_time))

		res <- intersect(res, names(all_abs_bifurcation_time_gene_names))

		regulator_time <- all_abs_bifurcation_time_gene_names[toupper(x)]
		if(length(regulator_time) != 1)
			regulator_time <- NA
		message('regulator_time is', res)
		
		data.frame(regulators = x, targets = res, branching_time = all_abs_bifurcation_time_gene_names[res], regulators_time = regulator_time)
	})

	df <- do.call(rbind.data.frame, significant_branching_TF_TF_5k_enrichment_gsc)
	# branching_tf_df <- data.frame(regulators = enrich_branching_TF_gsc_names, 
	# 								all_abs_bifurcation_time_gene_names[c(enrich_branching_TF_gsc_names[1:4], 'STAT1')])
	
	return(df)
}

test <- gen_branchTime_df(enrich_branching_TF_gsc_names = enrich_branching_TF_gsc_names)
df[df$regulators == 'STAT2::STAT1', 'regulators_time'] <- unique(df[df$regulators == 'STAT1', 'regulators_time']) 
pdf('tmp/ko_regulation_hierarchy.pdf', width = 1.5, height = 1.5)
qplot(regulators, abs(branching_time), color = regulators, geom = c('jitter', 'boxplot'), 
	data = df, alpha = I(0.5), log = 'y') +  
    ylab('bifurcation time point') + 
    geom_point(aes(regulators, abs(regulators_time)), color = 'black', data = df) + 
    coord_flip() + scale_y_continuous(breaks = round(seq(abs(min(df$regulators_time, na.rm = T)), max(df$branching_time, na.rm = T), by = 12),1)) + #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`)) 
    nm_theme()
dev.off()

#all the TFs enriched from the barplot: 
all_enriched_genes <- hyper_df_cluster2_order$gene_set
all_enriched_TF_branching_time_df1 <- gen_branchTime_df(enrich_branching_TF_gsc_names = all_enriched_genes, clusters_id = 1)
all_enriched_TF_branching_time_df2 <- gen_branchTime_df(enrich_branching_TF_gsc_names = all_enriched_genes, clusters_id = 2)
all_enriched_TF_branching_time_df3 <- gen_branchTime_df(enrich_branching_TF_gsc_names = all_enriched_genes, clusters_id = 3)
all_enriched_TF_branching_time_df4 <- gen_branchTime_df(enrich_branching_TF_gsc_names = all_enriched_genes, clusters_id = 4)
all_enriched_TF_branching_time_df5 <- gen_branchTime_df(enrich_branching_TF_gsc_names = all_enriched_genes, clusters_id = 5)
all_enriched_TF_branching_time_df6 <- gen_branchTime_df(enrich_branching_TF_gsc_names = all_enriched_genes, clusters_id = 6)

df[df$regulators == 'STAT2::STAT1', 'regulators_time'] <- unique(df[df$regulators == 'STAT1', 'regulators_time']) 
pdf('tmp/ko_regulation_hierarchy_cluster_1.pdf', width = 5, height = 5)
qplot(regulators, abs(branching_time), color = regulators, geom = c('jitter', 'boxplot'), 
	data = all_enriched_TF_branching_time_df1, alpha = I(0.5), log = 'y') +  
    ylab('bifurcation time point') + 
    geom_point(aes(regulators, abs(regulators_time)), color = 'black', data = all_enriched_TF_branching_time_df1) + 
    coord_flip() + scale_y_continuous(breaks = round(seq(abs(min(df$regulators_time, na.rm = T)), max(df$branching_time, na.rm = T), by = 12),1)) + #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`)) 
    nm_theme()
dev.off()

pdf('tmp/ko_regulation_hierarchy_cluster_2.pdf', width = 5, height = 5)
qplot(regulators, abs(branching_time), color = regulators, geom = c('jitter', 'boxplot'), 
	data = all_enriched_TF_branching_time_df2, alpha = I(0.5), log = 'y') +  
    ylab('bifurcation time point') + 
    geom_point(aes(regulators, abs(regulators_time)), color = 'black', data = all_enriched_TF_branching_time_df2) + 
    coord_flip() + scale_y_continuous(breaks = round(seq(abs(min(df$regulators_time, na.rm = T)), max(df$branching_time, na.rm = T), by = 12),1)) + #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`)) 
    nm_theme()
dev.off()

pdf('tmp/ko_regulation_hierarchy_cluster_3.pdf', width = 5, height = 5)
qplot(regulators, abs(branching_time), color = regulators, geom = c('jitter', 'boxplot'), 
	data = all_enriched_TF_branching_time_df3, alpha = I(0.5), log = 'y') +  
    ylab('bifurcation time point') + 
    geom_point(aes(regulators, abs(regulators_time)), color = 'black', data = all_enriched_TF_branching_time_df3) + 
    coord_flip() + scale_y_continuous(breaks = round(seq(abs(min(df$regulators_time, na.rm = T)), max(df$branching_time, na.rm = T), by = 12),1)) + #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`)) 
    nm_theme()
dev.off()

pdf('tmp/ko_regulation_hierarchy_cluster_4.pdf', width = 5, height = 5)
qplot(regulators, abs(branching_time), color = regulators, geom = c('jitter', 'boxplot'), 
	data = all_enriched_TF_branching_time_df4, alpha = I(0.5), log = 'y') +  
    ylab('bifurcation time point') + 
    geom_point(aes(regulators, abs(regulators_time)), color = 'black', data = all_enriched_TF_branching_time_df4) + 
    coord_flip() + scale_y_continuous(breaks = round(seq(abs(min(df$regulators_time, na.rm = T)), max(df$branching_time, na.rm = T), by = 12),1)) + #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`)) 
    nm_theme()
dev.off()

pdf('tmp/ko_regulation_hierarchy_cluster_5.pdf', width = 5, height = 5)
qplot(regulators, abs(branching_time), color = regulators, geom = c('jitter', 'boxplot'), 
	data = all_enriched_TF_branching_time_df5, alpha = I(0.5), log = 'y') +  
    ylab('bifurcation time point') + 
    geom_point(aes(regulators, abs(regulators_time)), color = 'black', data = all_enriched_TF_branching_time_df5) + 
    coord_flip() + scale_y_continuous(breaks = round(seq(abs(min(df$regulators_time, na.rm = T)), max(df$branching_time, na.rm = T), by = 12),1)) + #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`)) 
    nm_theme()
dev.off()

pdf('tmp/ko_regulation_hierarchy_cluster_6.pdf', width = 5, height = 8)
qplot(regulators, abs(branching_time), color = regulators, geom = c('jitter', 'boxplot'), 
	data = all_enriched_TF_branching_time_df6, alpha = I(0.5), log = 'y') +  
    ylab('bifurcation time point') + 
    geom_point(aes(regulators, abs(regulators_time)), color = 'black', data = all_enriched_TF_branching_time_df6) + 
    coord_flip() + scale_y_continuous(breaks = round(seq(abs(min(df$regulators_time, na.rm = T)), max(df$branching_time, na.rm = T), by = 12),1)) + #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`)) 
    nm_theme()
dev.off()

save.image('./RData/branchTimePoint.RData')





