shalek_custom_color_scale = c("Unstimulated."="#eff3ff", "On_Chip_Unstim."="#eff3ff", "LPS.1h"="#bdd7e7", "LPS.2h"="#6baed6", "LPS.4h"="#3182bd", "LPS.6h"="#08519c", "PIC.1h"="#bae4b3", "PIC.2h"="#74c476", "PIC.4h"="#31a354", "PIC.6h"="#006d2c", "PAM.1h"="#fdbe85", "PAM.2h"="#fd8d3c", "PAM.4h"="#e6550d", "PAM.6h"="#a63603", "LPS_GolgiPlug.4h_0h"="#fcae91", "LPS_GolgiPlug.4h_1h"="#fb6a4a", "LPS_GolgiPlug.4h_2h"="#cb181d", "Ifnar1_KO_LPS.2h"="#df65b0", "Ifnar1_KO_LPS.4h"="#dd1c77", "Stat1_KO_LPS.2h"="#969696", "Stat1_KO_LPS.4h"="#252525")

get_monocle_hsmm <- function(fpkm_matrix_file, sample_sheet_file, gene_annotation_file) {
	# Load the data files
	fpkm_matrix = read.table(fpkm_matrix_file, sep='\t', header=TRUE, row.names=1)
	sample_sheet = read.table(sample_sheet_file, sep='\t', header=TRUE, row.names=1)
	gene_ann = read.table(gene_annotation_file, sep='\t', header=TRUE, row.names=1)

	# Normalize the gene short names to all uppercase
	gene_ann$gene_short_name = toupper(gene_ann$gene_short_name)

	# Cuffnorm sometimes adds _0 to each sample name, strip this away
	colnames(fpkm_matrix) = gsub('_0$', '', colnames(fpkm_matrix))

	# Make sure the order of each dataframe is the same
	fpkm_matrix = fpkm_matrix[, row.names(sample_sheet)]
	gene_ann = gene_ann[row.names(fpkm_matrix), ]

	# Initialize data structures for monocle and return HSMM
	pd <- new("AnnotatedDataFrame", data = sample_sheet)
	fd <- new("AnnotatedDataFrame", data = gene_ann)
	HSMM <- new("CellDataSet", exprs = as.matrix(fpkm_matrix), phenoData = pd, featureData = fd, expressionFamily=tobit(), lowerDetectionLimit=0.1)

	# Initialize expression detection 
	HSMM <- detectGenes(HSMM, min_expr = 0.1)

	return(HSMM)
}

subset_HSMM_experiment <- function(HSMM, stimulants, experiment_names) {
	HSMM_subset =  HSMM[, pData(HSMM)$stimulation %in% stimulants & pData(HSMM)$experiment_name %in% experiment_names]
	return (HSMM_subset)
}

reduce_dimension_and_order_cells <- function(HSMM, ordering_genes, paths=1, root_cell=NULL, use_vst=FALSE, reverse=T) {
	HSMM <- setOrderingFilter(HSMM, ordering_genes)

	HSMM <- reduceDimension(HSMM, use_irlba = F, use_vst=use_vst, pseudo_expr=0.1, fun="exp")
	HSMM <- orderCells(HSMM, num_paths = paths, reverse = reverse, root_cell=root_cell)

	return(HSMM)
}

get_gene_cluster_ensembl_ids <- function(HSMM, gene_categories, clusters_of_interest) {
	cluster_subset = subset(gene_categories, cluster %in% clusters_of_interest)
	cluster_gene_ids = row.names(subset(fData(HSMM), gene_short_name %in% cluster_subset$gene_short_name))
	return(cluster_gene_ids)
}

get_gene_set_enrichment <- function(branch_test_results, ABCs, gene_set, ABCs_threshold=5, cores=1) {
	# Normalize the gene names in GO enrichment set so cases match
	for (name in names(gene_set$gsc)) {
	    gene_set$gsc[[name]] = toupper(gene_set$gsc[[name]])
	}

	# Clean up dataframe so no NA values
	ABCs = subset(ABCs, !is.na(ABCs))
	ABCs_top_values = subset(ABCs, abs(ABCs) > ABCs_threshold)

	# We only need the ABCs values here as a list, not the full dataframe
	ABCs_GSEA = ABCs[, 'ABCs']
	names(ABCs_GSEA) = ABCs[, 'gene_short_name']

	# Get the p-values from branch test for GSEA
	overlapping_branching_genes = branching_genes[ABCs[, 'gene_id'], ]
	pval = overlapping_branching_genes$pval
	names(pval) = overlapping_branching_genes$gene_short_name
	pval = pval[!is.na(pval)]
	ABCs_GSEA = ABCs_GSEA[names(pval)]

	# Now run the GSEA using the GO enrichment terms
	gsaRes_go <- runGSA(pval, ABCs_GSEA, gsc = gene_set, ncpus = cores)

	return(gsaRes_go)
}

plot_gene_set_enrichment <- function (gsa_res, gene_num = 20, selection = F,
    custom_p_adjust = F, add_terms = F)
{
	# Extract gene set enrichment data for more genes than needed (ensures will get equal AT1 and AT2 numbers)
    go_data <- make_enrichment_df(gsa_res, 1000000, custom_p_adjust, add_terms)
    go_data <- unique(go_data[, ])

    # Get only the significant values and assign to lineages based on adjusted pvalues
    AT1_values = go_data[go_data$"p adj (dist.dir.up)" < 0.05, ]
    AT1_values$lineage = 'Lineage 2'

    AT2_values = go_data[go_data$"p adj (dist.dir.dn)" < 0.05, ]
    AT2_values$lineage = 'Lineage 3'

    # Get the top significant values by stat
    top_AT1 = tail(AT1_values, gene_num)
    top_AT2 = head(AT2_values, gene_num)

    top_go_data = rbind(top_AT2, top_AT1)

    print(AT1_values$stat)
    print(AT2_values$stat)

    gsea_plot <- qplot(x = 1:nrow(top_go_data), y = abs(as.numeric(stat)), data = top_go_data, geom = "bar", stat = "identity", fill = lineage) +
        coord_flip() + 
        scale_x_discrete(limits = 1:nrow(top_go_data), labels = str_split_fixed(top_go_data$Name, "%", 2)[, 1]) +
        xlab("") + 
        ylab("Normalized Enrichment Score") + 
        scale_fill_discrete(name = "Type", labels = c("Lineage 2", "Lineage 3"))

    return(gsea_plot)
}

# This is simply an alternative version of Xiaojie's heatmap function that allows direct saving to a file with defined plot widths, etc.
plot_heatmap <- function (cds, ILRs_df, ABC_df, branching_pval_df, filename, height=10, width=10, rownames_type = c("gene_short_name",
    "ensemble_id"), ABC_type = c("positive", "negative", "all"),
    dist_method = "euclidean", hclust_method = "ward", ILRs_limit = 3,
    cluster_num = 4, ...)
{
    ILRs_df <- ILRs_df[!is.na(ILRs_df[, 1]), ]
    ILRs_df <- ILRs_df[row.names(ILRs_df) != "-", ]
    ILRs_df[which(ILRs_df <= -ILRs_limit)] <- -ILRs_limit
    ILRs_df[which(ILRs_df >= ILRs_limit)] <- ILRs_limit
    if (rownames_type == "rownames_type")
        branch_gene_ABCs <- subset(ABC_df, gene_short_name %in%
            row.names(ILRs_df))
    else if (rownames_type == "ensemble_id")
        branch_gene_ABCs <- ABC_df[row.names(ILRs_df), ]
    if (ABC_type == "positive")
        ILRs_df <- ILRs_df[unique(as.character(subset(branch_gene_ABCs,
            ABCs > 0)[, "gene_short_name"])), ]
    else if (ABC_type == "negative")
        ILRs_df <- ILRs_df[unique(as.character(subset(branch_gene_ABCs,
            ABCs < 0)[, "gene_short_name"])), ]
    else if (ABC_type == "all")
        ILRs_df <- ILRs_df
        
    test <- pheatmap(ILRs_df, cluster_cols = FALSE, clustering_distance_rows = dist_method,
        clustering_method = hclust_method, ...)

    annotation <- data.frame(class = as.factor(cutree(test$tree_row,
        cluster_num)), row.names = names(cutree(test$tree_row,
        cluster_num)))

    gene_names <- row.names(ILRs_df[test$tree_row$order, ])
    if (rownames_type == "gene_short_name") {
        
        ensemble_names <- row.names(subset(fData(cds), gene_short_name %in%
            gene_names))
        ensemble_names <- ensemble_names[!duplicated(fData(cds[ensemble_names,
            ])$gene_short_name)]
        log_qval <- log10(branching_pval_df[ensemble_names, "qval"])

    }
    else {
        log_qval <- log10(branching_pval_df[gene_names, "qval"])
    }
    annotation$log_qval <- -log_qval
    
    pheatmap(t(ILRs_df[test$tree_row$order, ]), cluster_cols = T,
        cluster_rows = F, show_rownames = F,
        border_color = NA, clustering_distance_cols = dist_method,
        clustering_method = hclust_method, annotation = annotation,
        annotation_legend = T, filename=filename, height=height, width=width, ...)

    return(annotation)
}

illustrator_theme <- function() {
 theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
   theme(panel.border = element_blank(), axis.line = element_line()) +
   theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
   theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
   theme(panel.background = element_rect(fill='white')) + 
   #theme(text = element_text(size=6)) + 
    theme(axis.text.y=element_text(size=6)) + 
    theme(axis.text.x=element_text(size=6)) +
    theme(axis.title.y=element_text(size=6)) + 
    theme(axis.title.x=element_text(size=6)) +
    theme(panel.border = element_blank(), axis.line = element_line(size = .1), axis.ticks = element_line(size = .1)) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(legend.position = "none") +
    theme(strip.text.x = element_text(colour="black", size=6)) + 
    theme(strip.text.y = element_text(colour="black", size=6)) + 
    theme(legend.title = element_text(colour="black", size = 6)) + 
    theme(legend.text = element_text(colour="black", size = 6)) + 
    theme(plot.margin=unit(c(0,0,0,0), "lines")) 
}

illustrator_theme2 <- function() {
    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    #theme(text = element_text(size=6)) +
    theme(axis.text.y=element_text(size=6)) +
    theme(axis.text.x=element_text(size=6)) +
    theme(axis.title.y=element_text(size=6)) +
    theme(axis.title.x=element_text(size=6)) +
    theme(panel.border = element_blank(), axis.line = element_line(size = .1), axis.ticks = element_line(size = .1)) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
    # theme(legend.position = "none") +
    theme(strip.text.x = element_text(colour="black", size=6)) +
    theme(strip.text.y = element_text(colour="black", size=6)) +
    theme(legend.title = element_text(colour="black", size = 6)) +
    theme(legend.text = element_text(colour="black", size = 6)) +
    theme(plot.margin=unit(c(0,0,0,0), "lines"))
}


