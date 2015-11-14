load('shalek_data_analysis.RData')
library(monocle)
library(xacHelper)
library(grid)
library(colorRamps)
library(R.utils)
library(piano)
library(venneuler)
fig_root_dir = './main/'
shalek_custom_color_scale_plus_states= c(shalek_custom_color_scale, c('1'='#40A43A', '2'='#CB1B1E', '3'='#3660A5', 'Unstimulated_Replicate.' = 'gray'))

#########################################################################################################
#figure 5: 

#########################################################################################################
#panel b: 
pdf(file =paste(fig_root_dir, 'figure_5A_new.pdf', sep = ''), height = 3, width = 3)
monocle::plot_spanning_tree(Shalek_abs_subset_ko_LPS, color_by="interaction(experiment_name, time)", cell_size=1) + 
scale_color_manual(values=shalek_custom_color_scale_plus_states) + illustrator_theme() 
dev.off()

#for making the legends: 
pdf(file =paste(fig_root_dir, 'figure_5A_new_helpA.pdf', sep = ''), height = 3, width = 3)
monocle::plot_spanning_tree(Shalek_abs_subset_ko_LPS, color_by="interaction(experiment_name, time)", cell_size=1) + 
scale_color_manual(values=shalek_custom_color_scale_plus_states) + illustrator_theme() 
dev.off()

pdf(file =paste(fig_root_dir, 'figure_5A_new_helpB.pdf', sep = ''), height = 3, width = 3)
monocle::plot_spanning_tree(Shalek_abs_subset_ko_LPS, color_by="State", cell_size=1) + 
scale_color_manual(values=shalek_custom_color_scale_plus_states) + illustrator_theme() 
dev.off()

#########################################################################################################
#panel c: 
fData(Shalek_abs_subset_ko_LPS)$num_cell_expressed <- esApply(Shalek_abs_subset_ko_LPS[, ], 1, function(x) sum(round(x) > 0))
ko_valid_expressed_genes <- row.names(subset(fData(Shalek_abs_subset_ko_LPS), num_cell_expressed > 5))

Shalek_abs_subset_ko_LPS_heatmap_annotations = plot_genes_branched_heatmap(Shalek_abs_subset_ko_LPS[intersect(row.names(subset(ko_branching_genes, qval < 0.05)), ko_valid_expressed_genes) ,], num_clusters=6, norm_method = "vstExprs", file_name=paste(fig_root_dir, 'figure_5C_new.pdf', sep = ''), cores=1, ABC_df=NULL, branchTest_df=ko_branching_genes, hmcols=NULL, lineage_labels = c('Normal cells', 'Knockout cells'))

# #test the NA: (branchTest and the kinetic plots)
# ko_branching_genes[row.names(subset(Shalek_abs_subset_ko_LPS_heatmap_annotations$annotation_row,  Cluster == 6)), ]

#new function: 
# cds_exprs <- buildLineageBranchCellDataSet(Shalek_abs_subset_ko_LPS[1:10, ])
# colour_cell <- rep(0, length(cds_exprs$experiment_name))
# names(colour_cell) <- as.character(cds_exprs$experiment_name)

# color_scheme <- brewer.pal(4, "Set2")

# colour_cell[names(colour_cell) == 'Ifnar1_KO_LPS'] <- color_scheme[1]
# colour_cell[names(colour_cell) ==  'Stat1_KO_LPS'] <- color_scheme[2]
# colour_cell[names(colour_cell) == 'LPS'] <- color_scheme[3]
# colour_cell[names(colour_cell) ==  'On_Chip_Unstim'] <- color_scheme[4]

# AT1_Lineage = "#BD1C7C"
# AT2_Lineage = "#337DB9" 

# colour <- rep(0, length(pData(cds_exprs)$Lineage))
# names(colour) <- as.character(pData(cds_exprs)$Lineage)
# names(colour)[pData(cds_exprs)$Lineage == '1'] <- 'Progenitor'
# names(colour)[pData(cds_exprs)$Lineage == '2'] <- 'Normal'
# names(colour)[pData(cds_exprs)$Lineage == '3'] <- 'Knockout'

# plot_genes_branched_pseudotime2(Shalek_abs_subset_ko_LPS[row.names(subset(Shalek_abs_subset_ko_LPS_heatmap_annotations$annotation_row,  Cluster == 6)), ], cell_color_by = "experiment_name", lineage_labels = c('Normal', 'Knockout'), 
# trajectory_color_by = "Lineage", fullModelFormulaStr = '~sm.ns(Pseudotime, df = 3)*Lineage', normalize = F, stretch = T,
# cell_size = 1, ncol = 2, reducedModelFormulaStr = "~ sm.ns(Pseudotime, df=3)", add_pval = T) + 
# ylab('Transcript counts') 

# # Get hyper geometric GSA test results for different enrichment sets (GO, KEGG, reactome, etc.)
# Shalek_abs_subset_ko_LPS_tree_heatmap_clusters <- as.numeric(Shalek_abs_subset_ko_LPS_heatmap_annotations$annotation_row$Cluster)
# names(Shalek_abs_subset_ko_LPS_tree_heatmap_clusters) = fData(Shalek_abs_subset_ko_LPS[row.names(Shalek_abs_subset_ko_LPS_heatmap_annotations$annotation_row), ])$gene_short_name
names(Shalek_abs_subset_ko_LPS_tree_heatmap_clusters) <- capitalize(tolower(names(Shalek_abs_subset_ko_LPS_tree_heatmap_clusters))) # normalize reactome terms names to all uppercase

#Shalek_abs_subset_ko_LPS_tree_hyper_geometric_results_go <- collect_gsa_hyper_results(Shalek_abs_subset_ko_LPS, gsc = mouse_go_gsc, Shalek_abs_subset_ko_LPS_tree_heatmap_clusters)
Shalek_abs_subset_ko_LPS_tree_hyper_geometric_results_reactome <- collect_gsa_hyper_results(Shalek_abs_subset_ko_LPS, gsc = mouse_reactome_gsc, Shalek_abs_subset_ko_LPS_tree_heatmap_clusters)
# Shalek_abs_subset_ko_LPS_tree_hyper_geometric_results_kegg <- collect_gsa_hyper_results(Shalek_abs_subset_ko_LPS, gsc = mouse_kegg_gsc, Shalek_abs_subset_ko_LPS_tree_heatmap_clusters)
# Shalek_abs_subset_ko_LPS_tree_hyper_geometric_results_go_cc <- collect_gsa_hyper_results(Shalek_abs_subset_ko_LPS, gsc = mouse_go_gsc_cc, Shalek_abs_subset_ko_LPS_tree_heatmap_clusters)
# Shalek_abs_subset_ko_LPS_tree_hyper_geometric_results_go_mf <- collect_gsa_hyper_results(Shalek_abs_subset_ko_LPS, gsc = mouse_go_gsc_mf, Shalek_abs_subset_ko_LPS_tree_heatmap_clusters)

## Plot each of the enrichment heatmaps (Reactome terms are the best)
# plot_gsa_hyper_heatmap(Shalek_abs_subset_ko_LPS, Shalek_abs_subset_ko_LPS_tree_hyper_geometric_results_go, significance=1e-25) + illustrator_theme() + 
# ggsave(filename = paste(fig_root_dir, 'figure_5D_go.pdf', sep = ''), height=15, width = 7)
pdf(file =paste(fig_root_dir, 'figure_5D_reactome.pdf', sep = ''), height = 15, width = 7)
plot_gsa_hyper_heatmap(Shalek_abs_subset_ko_LPS, Shalek_abs_subset_ko_LPS_tree_hyper_geometric_results_reactome, significance=1e-3) + illustrator_theme() 
dev.off()

save_hyper_df(Shalek_abs_subset_ko_LPS_tree_hyper_geometric_results_reactome, './supplementary_data/ko_hyper_df.xls') 
# plot_gsa_hyper_heatmap(Shalek_abs_subset_ko_LPS, Shalek_abs_subset_ko_LPS_tree_hyper_geometric_results_kegg, significance=1e-2) + illustrator_theme() + 
# ggsave(filename = paste(fig_root_dir, 'figure_5D_kegg.pdf', sep = ''), height=12, width = 7)
# plot_gsa_hyper_heatmap(Shalek_abs_subset_ko_LPS, Shalek_abs_subset_ko_LPS_tree_hyper_geometric_results_go_cc, significance=1e-10) + illustrator_theme() + 
# ggsave(filename =paste(fig_root_dir, 'figure_5D_go_cc.pdf', sep = ''), height=12, width = 7)
# plot_gsa_hyper_heatmap(Shalek_abs_subset_ko_LPS, Shalek_abs_subset_ko_LPS_tree_hyper_geometric_results_go_mf, significance=1e-10) + illustrator_theme() + 
# ggsave(filename =paste(fig_root_dir, 'figure_5D_go_mf.pdf', sep = ''), height=12, width = 7)

#overlap between branching genes: 
#make the venn diagram for all the above analysis (Fig 6):
andrew_element_all <- c(
            row.names(ko_Ifnar1_wt4[ko_Ifnar1_wt4$qval < .01, ]), #0to6
            row.names(ko_stat1_wt4[ko_stat1_wt4$qval < .01, ]), 
            row.names(subset(ko_branching_genes, qval < 0.01))
            )
andrew_sets_all <- c(
         rep(paste('ko_Ifnar1_wt4', sep = ''), length(row.names(ko_Ifnar1_wt4[ko_Ifnar1_wt4$qval < .01, ]))),
         rep(paste('ko_stat1_wt4', sep = ''), length(row.names(ko_stat1_wt4[ko_stat1_wt4$qval < .01, ]))),
         rep(paste('ko_branching_genes', sep = ''), length(row.names(subset(ko_branching_genes, qval < 0.01))))
         )

list(A = row.names(ko_Ifnar1_wt4[ko_Ifnar1_wt4$qval < .01, ]), #0to6
     B = row.names(ko_stat1_wt4[ko_stat1_wt4$qval < .01, ]), 
     C = row.names(subset(ko_branching_genes, qval < 0.01)))
# save(branch_pseudotime_element_all, branch_pseudotime_sets_all, file = 'branchTest_cmpr_subset')

pdf(file = paste(fig_root_dir, 'fig_SI_branchTest_cmpr.pdf', sep = ''), height = 2, width = 3)
#pdf(file = paste(fig_root_dir, 'fig6_SI_ko_overlapping.pdf', sep = ''))
venneuler_venn(andrew_element_all, andrew_sets_all)
dev.off()

#########################################################################################################
##panel d: 
##motif enrichment: 
#comparing with lung data: 
#run the fimo analysis on the DHS site
TF_5k_enrichment_gsc <- loadGSCSafe("./data/DC_JASPAR_5kb_hits_olap.gmt", encoding="latin1") 

for(i in 1:length(TF_5k_enrichment_gsc$gsc)) {
    TF_5k_enrichment_gsc$gsc[[i]] <- toupper(TF_5k_enrichment_gsc$gsc[[i]])
}

names(Shalek_abs_subset_ko_LPS_tree_heatmap_clusters) <- toupper(names(Shalek_abs_subset_ko_LPS_tree_heatmap_clusters))
subset(ko_branching_genes, toupper(gene_short_name) %in% names(Shalek_abs_subset_ko_LPS_tree_heatmap_clusters[Shalek_abs_subset_ko_LPS_tree_heatmap_clusters == 6]))

TF_enrichment_results_5k <- collect_gsa_hyper_results(Shalek_abs_subset_ko_LPS[, ], TF_5k_enrichment_gsc, Shalek_abs_subset_ko_LPS_tree_heatmap_clusters)

motif_enrich_plot <- plot_gsa_hyper_heatmap(Shalek_abs_subset_ko_LPS, TF_enrichment_results_5k, significance = 1e-1)

load('hyper_df')

max_cluster_id <- as.vector(which(table(Shalek_abs_subset_ko_LPS_tree_heatmap_clusters)  == max(table(Shalek_abs_subset_ko_LPS_tree_heatmap_clusters))))
hyper_df_cluster2 <- subset(hyper_df, cluster_id == max_cluster_id) #
hyper_df_cluster2_order <- hyper_df_cluster2[order(hyper_df_cluster2$pval), ]
hyper_df_cluster2_order$label <- hyper_df_cluster2_order$gene_set
STAT1_id <- which(hyper_df_cluster2_order$label == 'STAT1')
hyper_df_cluster2_order$label[c(5:STAT1_id, STAT1_id:nrow(hyper_df_cluster2_order))] <- ""

#add the barplot with the horizontal line: 
pdf(file = paste(fig_root_dir, 'figure_5_enrichment.pdf', sep = ''), height = 2, width = 3)
qplot(1:nrow(hyper_df_cluster2_order), - log10(qval), data=hyper_df_cluster2_order, geom = c('bar'), stat = 'identity', fill = 'red') +  
illustrator_theme() + geom_hline(yintercept = 1, linetype = 2, size  = 0.2, color = 'blue') + xlab('')
dev.off()

#########################################################################################################
#figure 6: 

#########################################################################################################
#panel b: 
pdf(file = paste(fig_root_dir, 'figure_6A_new.pdf', sep = ''), height = 3, width = 3)
monocle::plot_spanning_tree(Shalek_golgi_update, color_by="interaction(experiment_name, time)", cell_size=1) + 
scale_color_manual(values=shalek_custom_color_scale_plus_states) + illustrator_theme() 
dev.off()

#help figures for adding the legend: 
pdf(file = paste(fig_root_dir, 'figure_6A_new_helpA.pdf', sep = ''), height = 12, width = 12)
monocle::plot_spanning_tree(Shalek_golgi_update, color_by="interaction(experiment_name, time)", cell_size=1) + 
scale_color_manual(values=shalek_custom_color_scale_plus_states) 
dev.off()

pdf(file = paste(fig_root_dir, 'figure_6A_new_helpB.pdf', sep = ''), height = 12, width = 12)
monocle::plot_spanning_tree(Shalek_golgi_update, color_by="State", cell_size=1) + 
scale_color_manual(values=shalek_custom_color_scale_plus_states)
dev.off()

#########################################################################################################
#panel d
## Plot Heatmap
fData(Shalek_golgi_update)$num_cell_expressed <- esApply(Shalek_golgi_update[, ], 1, function(x) sum(round(x) > 0))
golgi_valid_expressed_genes <- row.names(subset(fData(Shalek_golgi_update), num_cell_expressed > 5))

Shalek_golgi_update_heatmap_annotations = plot_genes_branched_heatmap(Shalek_golgi_update[intersect(row.names(subset(golgi_branching_genes, qval < 0.05)), golgi_valid_expressed_genes),], num_clusters=6, norm_method = "vstExprs", lineage_labels = c('Normal', 'Golgi Plug'), 
file_name=paste(fig_root_dir, 'figure_6C.pdf', sep = ''), cores=1, ABC_df=NULL, branchTest_df=golgi_branching_genes, hmcols=NULL)

# Figure 6B annotations -- Enrichment analysis on clusters
Shalek_golgi_update_heatmap_clusters = as.numeric(Shalek_golgi_update_heatmap_annotations$annotation_row$Cluster)
names(Shalek_golgi_update_heatmap_clusters) = fData(Shalek_golgi_update[row.names(Shalek_golgi_update_heatmap_annotations$annotation_row), ])$gene_short_name
names(Shalek_golgi_update_heatmap_clusters) <- capitalize(tolower(names(Shalek_golgi_update_heatmap_clusters))) # normalize reactome terms names to all uppercase

pdf(file = paste(fig_root_dir, 'figure_6D.pdf', sep = ''), height = 20, width = 9)
Shalek_golgi_hyper_geometric_results = collect_gsa_hyper_results(Shalek_golgi_update, mouse_go_gsc, Shalek_golgi_update_heatmap_clusters)
plot_gsa_hyper_heatmap(Shalek_golgi_update, Shalek_golgi_hyper_geometric_results, significance=1e-2) 
dev.off()

## Get hyper geometric GSA test results for different enrichment sets (GO, KEGG, reactome, etc.)
Shalek_golgi_hyper_geometric_results_reactome <- collect_gsa_hyper_results(Shalek_golgi_update, gsc = mouse_reactome_gsc, Shalek_golgi_update_heatmap_clusters)
# Shalek_golgi_hyper_geometric_results_kegg <- collect_gsa_hyper_results(Shalek_golgi_update, gsc = mouse_kegg_gsc, Shalek_golgi_update_heatmap_clusters)
# Shalek_golgi_hyper_geometric_results_go_cc <- collect_gsa_hyper_results(Shalek_golgi_update, gsc = mouse_go_gsc_cc, Shalek_golgi_update_heatmap_clusters)
# Shalek_golgi_hyper_geometric_results_go_mf <- collect_gsa_hyper_results(Shalek_golgi_update, gsc = mouse_go_gsc_mf, Shalek_golgi_update_heatmap_clusters)

## Generate the heatmaps
pdf(file = paste(fig_root_dir, 'figure_6D_reactome.pdf', sep = ''), height = 30, width = 7)
plot_gsa_hyper_heatmap(Shalek_golgi_update_heatmap_annotations, Shalek_golgi_hyper_geometric_results_reactome, significance=1e-4) + illustrator_theme() 
dev.off()
# plot_gsa_hyper_heatmap(Shalek_golgi_update_heatmap_annotations, Shalek_golgi_hyper_geometric_results_kegg, significance=1e-2) + illustrator_theme() + 
# ggsave(paste(fig_root_dir, 'figure_6D_kegg.pdf', sep = ''), height=12, width = 7)
# plot_gsa_hyper_heatmap(Shalek_golgi_update_heatmap_annotations, Shalek_golgi_hyper_geometric_results_go_cc, significance=1e-10) + illustrator_theme() + 
# ggsave(paste(fig_root_dir, 'figure_6D_go_cc.pdf', sep = '') , height=12, width = 7)
# plot_gsa_hyper_heatmap(Shalek_golgi_update_heatmap_annotations, Shalek_golgi_hyper_geometric_results_go_mf, significance=1e-10) + illustrator_theme() + 
# ggsave(paste(fig_root_dir, 'figure_6D_go_mf.pdf', sep = '') , height=12, width = 7)

save_hyper_df(Shalek_golgi_hyper_geometric_results_reactome, './supplementary_data/golgiplug_hyper_df.xls') 
# save(branch_pseudotime_element_all, branch_pseudotime_sets_all, file = 'branchTest_cmpr_subset')

# pdf(file = paste(fig_root_dir, 'fig_SI_branchTest_cmpr.pdf', sep = ''), height = 2, width = 3)
#pdf(file = paste(fig_root_dir, 'fig6_SI_blocking.pdf', sep = ''))
#venneuler_venn(andrew_element_all, andrew_sets_all)
#dev.off()

#########################################################################################################
#panel d
#make the golgi-plug enrichment plot: 
##motif enrichment: 
#comparing with lung data: 
#run the fimo analysis on the DHS site
TF_5k_enrichment_gsc <- loadGSCSafe("./data/DC_JASPAR_5kb_hits_olap.gmt", encoding="latin1") 

for(i in 1:length(TF_5k_enrichment_gsc$gsc)) {
    TF_5k_enrichment_gsc$gsc[[i]] <- toupper(TF_5k_enrichment_gsc$gsc[[i]])
}

names(Shalek_golgi_update_heatmap_clusters) <- toupper(names(Shalek_golgi_update_heatmap_clusters))

golgi_TF_enrichment_results_5k <- collect_gsa_hyper_results(Shalek_golgi_update[, ], TF_5k_enrichment_gsc, Shalek_golgi_update_heatmap_clusters)

pdf(file = paste(fig_root_dir, 'figure_6D_motif_enrichment.pdf', sep = ''), height = 30, width = 7)
golgi_motif_enrich_plot <- plot_gsa_hyper_heatmap(Shalek_golgi_update, golgi_TF_enrichment_results_5k, significance = 1e-1)
dev.off()

#add kinetic plots for the branch genes in the TFs with motifs: 
load('hyper_df')

hyper_df[, 'first'] <- str_split_fixed(hyper_df[, 1], "::", 2)[,1]
hyper_df[, 'second'] <- str_split_fixed(hyper_df[, 1], "::", 2)[,2]

golgi_valid_gene_id_20_cell <- row.names(Shalek_golgi_update[which(rowSums(exprs(Shalek_golgi_update) >= 1) > 40), ])
golgi_valid_hyper_df <- subset(hyper_df, (toupper(first) %in% toupper(fData(Shalek_golgi_update[golgi_valid_gene_id_20_cell, ])$gene_short_name)) | 
(toupper(second) %in% toupper(fData(Shalek_golgi_update[golgi_valid_gene_id_20_cell, ])$gene_short_name))
) 

golgi_valid_hyper_df$sig <- T
pdf(file = paste(fig_root_dir, 'fig6_SI_enrichment.pdf', sep = ''), height = 1.2, width = 8)
qplot(cluster_id, gene_set, fill=sig, geom="tile", data=golgi_valid_hyper_df) + scale_fill_manual(values='black') + nm_theme() +
coord_flip() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab('') + ylab('')
dev.off()

# golgi_new_cds <- buildLineageBranchCellDataSet(Shalek_golgi_update[1:10, ], lineage_labels = c('Normal cells', 'GolgiPlug'))

# colour_cell <- rep(0, length(golgi_new_cds$Lineage))
# names(colour_cell) <- as.character(golgi_new_cds$State)
# colour_cell[names(colour_cell) == '1'] <- prog_cell_state
# colour_cell[names(colour_cell) == '2'] <- AT1_cell_state
# colour_cell[names(colour_cell) == '3'] <- AT2_cell_state

# colour <- rep(0, length(golgi_new_cds$Lineage))
# names(colour) <- as.character(pData(golgi_new_cds)$Lineage)
# colour[names(colour) == 'Normal cells'] <- AT1_Lineage
# colour[names(colour) ==  'GolgiPlug'] <- AT2_Lineage

# golgi_gene_grn_list <- infer_branch_gene_grn(cds = Shalek_golgi_update, branchTest_res = golgi_branching_genes, p_thrsld = 0.05,
# TF_enrichment_gsc = TF_5k_enrichment_gsc, motif_tfs_all = c(motif_TFs, motif_TFs_add),
# file = 'golgi_gene_regulatory_net_up5k')

# branch_motif_Tfs <- golgi_gene_grn_list$branch_tfs[toupper(golgi_gene_grn_list$branch_tfs) %in% toupper(golgi_valid_hyper_df$first) | 
#            toupper(golgi_gene_grn_list$branch_tfs) %in% toupper(golgi_valid_hyper_df$second)]
# branch_motif_Tfs_id <- row.names(subset(fData(Shalek_golgi_update), toupper(gene_short_name) %in% branch_motif_Tfs)) 

# pdf(file = paste(fig_root_dir, 'fig6_SI_branching_genes.pdf', sep = ''), height = 4, width = 15)
# plot_genes_branched_pseudotime2(Shalek_golgi_update[branch_motif_Tfs_id, ], cell_color_by = "State", 
#            trajectory_color_by = "Lineage", fullModelFormulaStr = '~sm.ns(Pseudotime, df = 3)*Lineage', normalize = F, stretch = T,
#            lineage_labels = c('Normal cells', 'GolgiPlug'), cell_size = 1, ncol = 8, reducedModelFormulaStr = "~sm.ns(Pseudotime, df=3)", add_pval = T) + 
#            ylab('Transcript counts') + nm_theme()
# dev.off()

########################################################################################################
#panel C: 
#make the venn diagram for all the above analysis (Fig 6):
andrew_element_all <- c(
            row.names(golgi_wt_0to6_pseudo[golgi_wt_0to6_pseudo$qval < .05, ]), #0to6
            row.names(all_golgi_plug0_wt4[all_golgi_plug0_wt4$qval < .05, ]), 
            row.names(all_golgi_plug0_wt0[all_golgi_plug0_wt0$qval < .05, ])
            )
andrew_sets_all <- c(
         rep(paste('golgi_wt_0to6_pseudo', sep = ''), length(row.names(golgi_wt_0to6_pseudo[golgi_wt_0to6_pseudo$qval < .05, ]))),
         rep(paste('all_golgi_plug0_wt4', sep = ''), length(row.names(all_golgi_plug0_wt4[all_golgi_plug0_wt4$qval < .05, ]))),
         rep(paste('all_golgi_plug_wt0', sep = ''), length(row.names(all_golgi_plug0_wt0[all_golgi_plug0_wt0$qval < .05, ])))
         )

pdf(file = paste(fig_root_dir, 'fig6_SI_blocking.pdf', sep = ''))
venneuler_venn(andrew_element_all, andrew_sets_all)
dev.off()

#this one is used in the paper: 
andrew_element_all <- c(
            row.names(golgi_wt_0to4_pseudo[golgi_wt_0to4_pseudo$qval < .05, ]), #0to6
            row.names(golgi_plug0_wt4[golgi_plug0_wt4$qval < .05, ]), 
            row.names(golgi_plug0_wt0[golgi_plug0_wt0$qval < .05, ])
            )
andrew_sets_all <- c(
         rep(paste('golgi_wt_0to4_pseudo', sep = ''), length(row.names(golgi_wt_0to4_pseudo[golgi_wt_0to4_pseudo$qval < .05, ]))),
         rep(paste('all_golgi_plug0_wt4', sep = ''), length(row.names(golgi_plug0_wt4[golgi_plug0_wt4$qval < .05, ]))),
         rep(paste('all_golgi_plug_wt0', sep = ''), length(row.names(golgi_plug0_wt0[golgi_plug0_wt0$qval < .05, ])))
         )

pdf(file = paste(fig_root_dir, 'fig6_SI_blocking_4h.pdf', sep = ''))
venneuler_venn(andrew_element_all, andrew_sets_all)
dev.off()

# test the branch gene number to two group test DEGs:  
andrew_element_all <- c(
            row.names(golgi_branching_genes[golgi_branching_genes$qval < .01, ]), #0to6
            row.names(all_golgi_plug0_wt4[all_golgi_plug0_wt4$qval < .01, ]), 
            row.names(all_golgi_plug0_wt0[all_golgi_plug0_wt0$qval < .01, ])
            )
andrew_sets_all <- c(
         rep(paste('golgi_branching_genes', sep = ''), length(row.names(golgi_branching_genes[golgi_branching_genes$qval < .01, ]))),
         rep(paste('all_golgi_plug0_wt4', sep = ''), length(row.names(all_golgi_plug0_wt4[all_golgi_plug0_wt4$qval < .01, ]))),
         rep(paste('all_golgi_plug_wt0', sep = ''), length(row.names(all_golgi_plug0_wt0[all_golgi_plug0_wt0$qval < .01, ])))
         )

pdf(file = paste(fig_root_dir, 'fig6_golgi_branch_overlapping.pdf', sep = ''))
venneuler_venn(andrew_element_all, andrew_sets_all)
dev.off()

# #TFs for the two group tests, jitter plots, violin plot and the corresponding kinetic plots: 
# DEG_tf_list <- intersect(subset(golgi_plug0_wt0, qval < 1e-1)[, 'gene_short_name'],  c(motif_TFs, motif_TFs_add))
# write.table(subset(golgi_plug0_wt0, qval < 1e-1)[, 'gene_short_name'], file = 'golgi_plug0_wt0_deg', quote = F, row.names= F, col.names=F)
# write.table(subset(golgi_plug0_wt4, qval < 1e-1)[, 'gene_short_name'], file = 'golgi_plug0_wt4_deg', quote = F, row.names= F, col.names=F)

# DEG_tf_ids <- row.names(subset(fData(Shalek_golgi_update), toupper(gene_short_name) %in% c(DEG_tf_list, 'TNF', 'IL6', 'IL-1')))

# plot_genes_violin(Shalek_golgi_update[DEG_tf_ids, pData(Shalek_golgi_update)$time %in% c('4h_0h', '4h', '1h', '2h', '')], pvalue = golgi_plug0_wt0[DEG_tf_ids, 'pval'], grouping = "time", min_expr = 0.1,
# cell_size = 0.75, nrow = NULL, ncol = 1, panel_order = NULL, color_by = 'time',
# plot_trend = FALSE, label_by_short_name = TRUE,
# outlier_rm_type = NULL)
# ggsave('time_DEG.pdf', height = 10, width = 5)

# plot_genes_violin(Shalek_golgi_update[DEG_tf_ids, ], pvalue = golgi_plug0_wt0[DEG_tf_ids, 'pval'], grouping = "time", min_expr = 0.1,
# cell_size = 0.75, nrow = NULL, ncol = 1, panel_order = NULL, color_by = 'time',
# plot_trend = FALSE, label_by_short_name = TRUE,
# outlier_rm_type = NULL)
# ggsave('time_DEG_all_time.pdf', height = 10, width = 8)

# plot_genes_jitter(Shalek_golgi_update[DEG_tf_ids, pData(Shalek_golgi_update)$time %in% c('4h_0h', '4h', '1h', '2h', '')], #pvalue = golgi_plug0_wt0[DEG_tf_ids, 'pval'], 
# grouping = "time", min_expr = 0.1,
# cell_size = 0.75, nrow = NULL, ncol = 1, panel_order = NULL,
# color_by = NULL, plot_trend = FALSE, label_by_short_name = TRUE,
# )
# # DEG_tf_ids <- row.names(subset(fData(Shalek_golgi_update), toupper(gene_short_name) %in% DEG_tf_list)) 
# plot_genes_branched_pseudotime2(Shalek_golgi_update[DEG_tf_ids, ], cell_color_by = "State", 
#            trajectory_color_by = "Lineage", fullModelFormulaStr = '~sm.ns(Pseudotime, df = 3)*Lineage', normalize = F, stretch = T,
#            lineage_labels = c('Normal cells', 'GolgiPlug'), cell_size = 1, ncol = 8, reducedModelFormulaStr = "~sm.ns(Pseudotime, df=3)", add_pval = T) + 
#            ylab('Transcript counts') #+ nm_theme()
# ggsave(paste(fig_root_dir, 'fig6_SI_DEG_genes.pdf', sep = ''), height = 4, width = 15)

#figure 6: for golgi-plug (Note that the GSC are the same): 
# TF_enrichment_results_5k_golgi <- collect_gsa_hyper_results(Shalek_golgi_update[, ], TF_5k_enrichment_gsc, Shalek_golgi_update_heatmap_clusters)

# motif_enrich_plot_golgi <- plot_gsa_hyper_heatmap(Shalek_golgi_update, TF_enrichment_results_5k_golgi, significance = 1e-1)

# motif_enrich_plot_golgi + illustrator_theme() + ggsave(paste(fig_root_dir, "/figure_6_enrichment_heatmap.pdf", sep = ''), height = 2, width = 3)

save.image('shalek_data_analysis_figs_Rdata')
