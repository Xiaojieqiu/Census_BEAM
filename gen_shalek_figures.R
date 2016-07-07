# library(devtools)
# load_all('~/Projects/monocle-dev')
library(monocle)
library(xacHelper)
library(grid)
library(igraph)
library(RColorBrewer)
library(colorRamps) 
library(R.utils)
library(piano)
library(venneuler)
library(pheatmap)
library(plyr)
library(stringr)

load('./RData/analysis_shalek_data.RData')
fig_root_dir = './main_figures/'

#########################################################################################################
#figure 5: 

#########################################################################################################
#panel b: 
pdf(file =paste(fig_root_dir, 'fig5b.pdf', sep = ''), height = 3, width = 3)
plot_spanning_tree(Shalek_abs_subset_ko_LPS, color_by="interaction(experiment_name, time)", cell_size=2) + 
scale_color_manual(values=shalek_custom_color_scale_plus_states) + nm_theme()
dev.off()

#for making the legends: 
pdf(file =paste('./tmp/', 'fig5b_helpA.pdf', sep = ''), height = 3, width = 3)
monocle::plot_spanning_tree(Shalek_abs_subset_ko_LPS, color_by="interaction(experiment_name, time)", cell_size=2) + 
scale_color_manual(values=shalek_custom_color_scale_plus_states) + nm_theme()
dev.off()

pdf(file =paste('./tmp/', 'fig5b_helpB.pdf', sep = ''), height = 3, width = 3)
monocle::plot_spanning_tree(Shalek_abs_subset_ko_LPS, color_by="State", cell_size=2) + 
scale_color_manual(values=shalek_custom_color_scale_plus_states) + nm_theme()
dev.off()

#########################################################################################################
#panel c: 
fData(Shalek_abs_subset_ko_LPS)$num_cell_expressed <- esApply(Shalek_abs_subset_ko_LPS[, ], 1, function(x) sum(round(x) > 0))
ko_valid_expressed_genes <- row.names(subset(fData(Shalek_abs_subset_ko_LPS), num_cell_expressed > 5))

rm(plot_genes_branched_heatmap)
pdf(paste(fig_root_dir, 'fig5c_ori.pdf', sep = ''))
Shalek_abs_subset_ko_LPS_heatmap_annotations = plot_genes_branched_heatmap(Shalek_abs_subset_ko_LPS[intersect(row.names(subset(ko_branching_genes, qval < 0.05)), ko_valid_expressed_genes) ,], num_clusters=6, norm_method = "vstExprs", 
    cores=detectCores(), ABC_df=NULL, branchTest_df=ko_branching_genes, hmcols=NULL, lineage_labels = c('Normal cells', 'Knockout cells'), return_all = T)
dev.off()

# Get hyper geometric GSA test results for different enrichment sets (GO, KEGG, reactome, etc.)
Shalek_abs_subset_ko_LPS_tree_heatmap_clusters <- as.numeric(Shalek_abs_subset_ko_LPS_heatmap_annotations$annotation_row$Cluster)
names(Shalek_abs_subset_ko_LPS_tree_heatmap_clusters) = fData(Shalek_abs_subset_ko_LPS[row.names(Shalek_abs_subset_ko_LPS_heatmap_annotations$annotation_row), ])$gene_short_name
names(Shalek_abs_subset_ko_LPS_tree_heatmap_clusters) <- capitalize(tolower(names(Shalek_abs_subset_ko_LPS_tree_heatmap_clusters))) # normalize reactome terms names to all uppercase

Shalek_abs_subset_ko_LPS_tree_hyper_geometric_results_reactome <- collect_gsa_hyper_results(Shalek_abs_subset_ko_LPS, gsc = mouse_reactome_gsc, Shalek_abs_subset_ko_LPS_tree_heatmap_clusters)

plot_gsa_hyper_heatmap <- function (cds, gsa_results, significance = 0.05, sign_type = "qval")
{
    hyper_df <- ldply(gsa_results, function(gsa_res) {
        data.frame(gene_set = names(gsa_res$pvalues), pval = gsa_res$pvalues,
            qval = gsa_res$p.adj)
    })
    hyper_df$qval <- p.adjust(hyper_df$pval, method = 'fdr')
    colnames(hyper_df)[1] <- "cluster_id"
    hyper_df <- subset(hyper_df, hyper_df[, sign_type] <= significance)
    hyper_df <- merge(hyper_df, ddply(hyper_df, .(gene_set),
        function(x) {
            nrow(x)
        }), by = "gene_set")
    save(hyper_df, file = "hyper_df")
    hyper_df$gene_set <- factor(hyper_df$gene_set, levels = unique(arrange(hyper_df,
        V1, cluster_id)$gene_set))
    qplot(cluster_id, gene_set, fill = -log10(qval), geom = "tile",
        data = hyper_df) + scale_fill_gradientn(colours = rainbow(7))
}

pdf(file =paste(fig_root_dir, 'fig5c_reactome.pdf', sep = ''), height = 15, width = 7)
plot_gsa_hyper_heatmap(Shalek_abs_subset_ko_LPS, Shalek_abs_subset_ko_LPS_tree_hyper_geometric_results_reactome, significance=1e-1) + nm_theme()
dev.off()

save_hyper_df(Shalek_abs_subset_ko_LPS_tree_hyper_geometric_results_reactome, './supplementary_data/ko_hyper_df.xls') 

Shalek_abs_subset_ko_LPS_tree_hyper_geometric_results_go <- collect_gsa_hyper_results(Shalek_abs_subset_ko_LPS, gsc = mouse_go_gsc, Shalek_abs_subset_ko_LPS_tree_heatmap_clusters)

pdf(file =paste(fig_root_dir, 'fig5c_go.pdf', sep = ''), height = 70, width = 12)
plot_gsa_hyper_heatmap(Shalek_abs_subset_ko_LPS, Shalek_abs_subset_ko_LPS_tree_hyper_geometric_results_go, significance=5e-2) + nm_theme()
dev.off()

save_hyper_df(Shalek_abs_subset_ko_LPS_tree_hyper_geometric_results_go, './supplementary_data/go_ko_hyper_df.xls') 

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

pdf(file = paste('./tmp/', 'fig_SI_branchTest_cmpr.pdf', sep = ''), height = 2, width = 3)
venneuler_venn(andrew_element_all, andrew_sets_all)
dev.off()

#########################################################################################################
##panel d: 
##motif enrichment: 
#comparing with lung data: 
#run the fimo analysis on the DHS site
shalek_TF_5k_enrichment_gsc <- loadGSCSafe("./data/DC_JASPAR_5kb_hits_olap.gmt", encoding="latin1") 

for(i in 1:length(shalek_TF_5k_enrichment_gsc$gsc)) {
    shalek_TF_5k_enrichment_gsc$gsc[[i]] <- toupper(shalek_TF_5k_enrichment_gsc$gsc[[i]])
}

names(Shalek_abs_subset_ko_LPS_tree_heatmap_clusters) <- toupper(names(Shalek_abs_subset_ko_LPS_tree_heatmap_clusters))
subset(ko_branching_genes, toupper(gene_short_name) %in% names(Shalek_abs_subset_ko_LPS_tree_heatmap_clusters[Shalek_abs_subset_ko_LPS_tree_heatmap_clusters == 6]))

TF_enrichment_results_5k <- collect_gsa_hyper_results(Shalek_abs_subset_ko_LPS[, ], shalek_TF_5k_enrichment_gsc, Shalek_abs_subset_ko_LPS_tree_heatmap_clusters)

motif_enrich_plot <- plot_gsa_hyper_heatmap(Shalek_abs_subset_ko_LPS, TF_enrichment_results_5k, significance = 1e-1)

load('hyper_df')

max_cluster_id <- as.vector(which(table(Shalek_abs_subset_ko_LPS_tree_heatmap_clusters)  == max(table(Shalek_abs_subset_ko_LPS_tree_heatmap_clusters))))
hyper_df_cluster2 <- subset(hyper_df, cluster_id == max_cluster_id) #
hyper_df_cluster2 <- hyper_df_cluster2[hyper_df_cluster2$qval !=0 , ] #remove outlier

hyper_df_cluster2_order <- hyper_df_cluster2[order(hyper_df_cluster2$pval), ]
hyper_df_cluster2_order$label <- hyper_df_cluster2_order$gene_set
STAT1_id <- which(hyper_df_cluster2_order$label == 'STAT1')
hyper_df_cluster2_order$label[c(5:STAT1_id, STAT1_id:nrow(hyper_df_cluster2_order))] <- ""


#add the barplot with the horizontal line: 
pdf(file = paste(fig_root_dir, 'fig5d_ori.pdf', sep = ''), height = 2, width = 3)
qplot(1:nrow(hyper_df_cluster2_order), - log10(qval), data=hyper_df_cluster2_order, geom = c('bar'), stat = 'identity', fill = 'red') +  
nm_theme()+ geom_hline(yintercept = 1, linetype = 2, size  = 0.2, color = 'blue') + xlab('')
dev.off()

# #########################################################################################################
# #figure 6: 

# #########################################################################################################
# #panel b: 
# pdf(file = paste(fig_root_dir, 'fig6b.pdf', sep = ''), height = 3, width = 3)
# monocle::plot_spanning_tree(Shalek_golgi_update, color_by="interaction(experiment_name, time)", cell_size=2) + 
# scale_color_manual(values=shalek_custom_color_scale_plus_states) + nm_theme()
# dev.off()

# #help figures for adding the legend: 
# pdf(file = paste('./supplementary_figures/', 'fig6b_heplerA.pdf', sep = ''), height = 12, width = 12)
# monocle::plot_spanning_tree(Shalek_golgi_update, color_by="interaction(experiment_name, time)", cell_size=2) + 
# scale_color_manual(values=shalek_custom_color_scale_plus_states) 
# dev.off()

# pdf(file = paste('./supplementary_figures/', 'fig6b_heplerB.pdf', sep = ''), height = 12, width = 12)
# monocle::plot_spanning_tree(Shalek_golgi_update, color_by="State", cell_size=2) + 
# scale_color_manual(values=shalek_custom_color_scale_plus_states)
# dev.off()

# #########################################################################################################
# #panel d
# ## Plot Heatmap
# fData(Shalek_golgi_update)$num_cell_expressed <- esApply(Shalek_golgi_update[, ], 1, function(x) sum(round(x) > 0))
# golgi_valid_expressed_genes <- row.names(subset(fData(Shalek_golgi_update), num_cell_expressed > 5))

# pdf(paste(fig_root_dir, 'fig6c.pdf', sep = ''))
# Shalek_golgi_update_heatmap_annotations = monocle::plot_genes_branched_heatmap(Shalek_golgi_update[intersect(row.names(subset(golgi_branching_genes, qval < 0.05)), golgi_valid_expressed_genes),], num_clusters=6, norm_method = "vstExprs", lineage_labels = c('Normal', 'Golgi Plug'), 
# cores=detectCores(), ABC_df=NULL, branchTest_df=golgi_branching_genes, hmcols=NULL, return_all = T)
# dev.off()

# # Figure 6B annotations -- Enrichment analysis on clusters
# Shalek_golgi_update_heatmap_clusters = as.numeric(Shalek_golgi_update_heatmap_annotations$annotation_row$Cluster)
# names(Shalek_golgi_update_heatmap_clusters) = fData(Shalek_golgi_update[row.names(Shalek_golgi_update_heatmap_annotations$annotation_row), ])$gene_short_name
# names(Shalek_golgi_update_heatmap_clusters) <- capitalize(tolower(names(Shalek_golgi_update_heatmap_clusters))) # normalize reactome terms names to all uppercase

# ## Get hyper geometric GSA test results for different enrichment sets (GO, KEGG, reactome, etc.)
# Shalek_golgi_hyper_geometric_results_reactome <- collect_gsa_hyper_results(Shalek_golgi_update, gsc = mouse_reactome_gsc, Shalek_golgi_update_heatmap_clusters)

# ## Generate the heatmaps
# pdf(file = paste(fig_root_dir, 'fig6c_reactome.pdf', sep = ''), height = 30, width = 7)
# plot_gsa_hyper_heatmap(Shalek_golgi_update_heatmap_annotations, Shalek_golgi_hyper_geometric_results_reactome, significance=1e-1) + nm_theme()
# dev.off()

# save_hyper_df(Shalek_golgi_hyper_geometric_results_reactome, './supplementary_data/golgiplug_hyper_df.xls') 

# #########################################################################################################
# #panel d
# #make the golgi-plug enrichment plot: 
# ##motif enrichment: 
# #comparing with lung data: 
# #run the fimo analysis on the DHS site
# shalek_TF_5k_enrichment_gsc <- loadGSCSafe("./data/DC_JASPAR_5kb_hits_olap.gmt", encoding="latin1") 

# for(i in 1:length(shalek_TF_5k_enrichment_gsc$gsc)) {
#     shalek_TF_5k_enrichment_gsc$gsc[[i]] <- toupper(shalek_TF_5k_enrichment_gsc$gsc[[i]])
# }
# names(Shalek_golgi_update_heatmap_clusters) <- toupper(names(Shalek_golgi_update_heatmap_clusters))

# golgi_TF_enrichment_results_5k <- collect_gsa_hyper_results(Shalek_golgi_update[, ], shalek_TF_5k_enrichment_gsc, Shalek_golgi_update_heatmap_clusters)

# pdf(file = paste('./tmp/', 'fig6d_motif_enrichment.pdf', sep = ''), height = 30, width = 7)
# golgi_motif_enrich_plot <- plot_gsa_hyper_heatmap(Shalek_golgi_update, golgi_TF_enrichment_results_5k, significance = 1e-1)
# dev.off()

# #add kinetic plots for the branch genes in the TFs with motifs: 
# load('hyper_df')

# hyper_df[, 'first'] <- str_split_fixed(hyper_df[, 1], "::", 2)[,1]
# hyper_df[, 'second'] <- str_split_fixed(hyper_df[, 1], "::", 2)[,2]

# golgi_valid_gene_id_20_cell <- row.names(Shalek_golgi_update[which(rowSums(exprs(Shalek_golgi_update) >= 1) > 40), ])
# golgi_valid_hyper_df <- subset(hyper_df, (toupper(first) %in% toupper(fData(Shalek_golgi_update[golgi_valid_gene_id_20_cell, ])$gene_short_name)) | 
# (toupper(second) %in% toupper(fData(Shalek_golgi_update[golgi_valid_gene_id_20_cell, ])$gene_short_name))
# ) 

# golgi_valid_hyper_df$sig <- T
# pdf(file = paste('./tmp/', 'fig6_SI_enrichment.pdf', sep = ''), height = 1.2, width = 8)
# qplot(cluster_id, gene_set, fill=sig, geom="tile", data=golgi_valid_hyper_df) + scale_fill_manual(values='black') + nm_theme() +
# coord_flip() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab('') + ylab('')
# dev.off()

# ########################################################################################################
# #panel C: 
# #make the venn diagram for all the above analysis (Fig 6):
# andrew_element_all <- c(
#             row.names(golgi_wt_0to6_pseudo[golgi_wt_0to6_pseudo$qval < .05, ]), #0to6
#             row.names(all_golgi_plug0_wt4[all_golgi_plug0_wt4$qval < .05, ]), 
#             row.names(all_golgi_plug0_wt0[all_golgi_plug0_wt0$qval < .05, ])
#             )
# andrew_sets_all <- c(
#          rep(paste('golgi_wt_0to6_pseudo', sep = ''), length(row.names(golgi_wt_0to6_pseudo[golgi_wt_0to6_pseudo$qval < .05, ]))),
#          rep(paste('all_golgi_plug0_wt4', sep = ''), length(row.names(all_golgi_plug0_wt4[all_golgi_plug0_wt4$qval < .05, ]))),
#          rep(paste('all_golgi_plug_wt0', sep = ''), length(row.names(all_golgi_plug0_wt0[all_golgi_plug0_wt0$qval < .05, ])))
#          )

# pdf(file = paste(fig_root_dir, 'fig6d_6h.pdf', sep = ''))
# venneuler_venn(andrew_element_all, andrew_sets_all)
# dev.off()

# #this one is used in the paper: 
# andrew_element_all <- c(
#             row.names(golgi_wt_0to4_pseudo[golgi_wt_0to4_pseudo$qval < .05, ]), #0to6
#             row.names(golgi_plug0_wt4[golgi_plug0_wt4$qval < .05, ]), 
#             row.names(golgi_plug0_wt0[golgi_plug0_wt0$qval < .05, ])
#             )
# andrew_sets_all <- c(
#          rep(paste('golgi_wt_0to4_pseudo', sep = ''), length(row.names(golgi_wt_0to4_pseudo[golgi_wt_0to4_pseudo$qval < .05, ]))),
#          rep(paste('all_golgi_plug0_wt4', sep = ''), length(row.names(golgi_plug0_wt4[golgi_plug0_wt4$qval < .05, ]))),
#          rep(paste('all_golgi_plug_wt0', sep = ''), length(row.names(golgi_plug0_wt0[golgi_plug0_wt0$qval < .05, ])))
#          )

# pdf(file = paste(fig_root_dir, 'fig6d.pdf', sep = ''))
# venneuler_venn(andrew_element_all, andrew_sets_all)
# dev.off()

# # test the branch gene number to two group test DEGs:  
# andrew_element_all <- c(
#             row.names(golgi_branching_genes[golgi_branching_genes$qval < .01, ]), #0to6
#             row.names(all_golgi_plug0_wt4[all_golgi_plug0_wt4$qval < .01, ]), 
#             row.names(all_golgi_plug0_wt0[all_golgi_plug0_wt0$qval < .01, ])
#             )
# andrew_sets_all <- c(
#          rep(paste('golgi_branching_genes', sep = ''), length(row.names(golgi_branching_genes[golgi_branching_genes$qval < .01, ]))),
#          rep(paste('all_golgi_plug0_wt4', sep = ''), length(row.names(all_golgi_plug0_wt4[all_golgi_plug0_wt4$qval < .01, ]))),
#          rep(paste('all_golgi_plug_wt0', sep = ''), length(row.names(all_golgi_plug0_wt0[all_golgi_plug0_wt0$qval < .01, ])))
#          )

# pdf(file = paste('./tmp/', 'fig6_golgi_branch_overlapping.pdf', sep = ''))
# venneuler_venn(andrew_element_all, andrew_sets_all)
# dev.off()

save.image('./RData/gen_shalek_figures.RData')
