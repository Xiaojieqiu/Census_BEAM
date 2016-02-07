# load('./RData/analysis_shalek_data.RData')

library(monocle)
library(xacHelper) #plot_clustering_heatmap
library(pheatmap)
library(colorRamps)

#select unstimulated cells: 
shalek_unstimulated_cells <- row.names(subset(pData(Shalek_golgi_update), experiment_name == 'Unstimulated_Replicate'))
shalek_unstimulated_cells_cds <- Shalek_golgi_update[, shalek_unstimulated_cells]
dim(shalek_unstimulated_cells_cds) #only 49 cells left

#use the cluster cells function: 

#use the list of genes from the paper: 
fig6a_mac <- c('Ptplad2', '1810011H11Rik', 'Tlr4', 'Fgd4', 'Sqrdl', 
	'Csf3r', 'Plod1', 'Tom1', 'Pld3', 'Tpp1', 'Ctsd', 'Lamp2', 'Pla2g4a', 
	'Fcgr1', 'Mr1', 'Mertk', 'Cd14', 'Tbxas1', 'Fcgr1', 'Sepp1', 'Cd164', 'Tcn2',
	'Dok3', 'Ctsl', 'Tspan14')

fig6a_dc <- c('Adam19', 'Ccr7', 'Flt3', 'Gpr132', 'H2-Eb2', 'Hmgn3', 'Kit', 'Klri1', 
	'Kmo', 'P2ry10', 'Pvrl1', 'Rab30', 'Slamf7', 'Traf1', 'Zbtb46', 'Dpp4', 'Runx3')

fig6c_lps_response <- c('Ccl4', 'Ccl3', 'Cxcl2', 'Ccrl2', 'Tank', 'Cflar')

fig6g <- c('Ifit2', 'Oasl2', 'Ifi204', 'Ddx58', 'Clec4e', 'Ifit3', 'iNos', 'Il19', 'Il12a', 'Ccl7', 'Ifnb1', 'Il33')

marker_genes_paper <- c(fig6a_mac, fig6a_dc) #, fig6c_lps_response, fig6g

#plot_clustering_heatmap from xacHelper package
plot_clustering_heatmap <- function (cds, num_clusters = 3, height = 3, width = 3)
{
    exprs_df <- exprs(cds)
    exprs_df <- exprs_df[apply(exprs_df, 1, function(x) !any(is.na(x))),
        ]
    exprs_df <- t(scale(t(exprs_df)))
    exprs_df[exprs_df > 3] <- 3
    exprs_df[exprs_df < -3] <- -3
    row_dist <- as.dist((1 - cor(t(exprs_df)))/2)
    row_dist[is.na(row_dist)] <- 1
    bk <- seq(-3.1, 3.1, by = 0.1)
    hmcols <- blue2green2red(length(bk) - 1)
    row.names(exprs_df) <- fData(cds[row.names(exprs_df), ])$gene_short_name
    ph_res <- pheatmap(exprs_df, cluster_cols = T, cluster_rows = TRUE,
        show_rownames = T, show_colnames = T, clustering_distance_rows = row_dist,
        clustering_method = "ward.D2", color = hmcols, fontsize = 6, treeheight_row = 3, treeheight_col = 3, 
        cutree_cols = num_clusters, height = height, width = width)

    grid::grid.rect(gp = grid::gpar("fill"))
    grid::grid.draw(ph_res$gtable)

    return(ph_res)
}

pdf('./nbt_2nd_sub_reviewers/subpopulation_heatmap.pdf')
clustering_results <- plot_clustering_heatmap(shalek_unstimulated_cells_cds[fData(Shalek_golgi_update)$gene_short_name %in% marker_genes_paper, ])
dev.off()
as.numeric(clustering_results$annotation_col$Cluster)

Cluster_df <- data.frame(Cluster = factor(cutree(clustering_results$tree_col, 3)))
pData(Shalek_golgi_update)$Clusters <- 0
pData(Shalek_golgi_update)[row.names(Cluster_df), 'Clusters'] <- as.vector(Cluster_df$Cluster)

#color the cells by the clusters from the heatmap: 
pdf('./tmp/cluster_unstimulted_cells.pdf', height = 2, width = 3)
plot_spanning_tree(Shalek_golgi_update, color = 'Clusters') + nm_theme()
dev.off()

#compare with the results from the immunity paper: 
# shalek_unstimulated_cells_cds <- reduceDimension(shalek_unstimulated_cells_cds, method = 'DDRTree', ncenter = 2)
# pData(shalek_unstimulated_cells_cds)$Clusters <- 0
# pData(shalek_unstimulated_cells_cds)[row.names(Cluster_df), 'Clusters'] <- as.vector(Cluster_df$Cluster)
# pdf('./tmp/unstimulated_DDRTree.pdf', height = 2, width = 3)
# plot_spanning_tree(shalek_unstimulated_cells_cds, color = 'Clusters') #
# dev.off()

save.image('./RData/subpopulation.RData')
