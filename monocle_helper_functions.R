# library(grid)
# # library(DevTree)
# library(monocle)

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

plot_genes_branched_heatmap <- function(cds_subset, num_clusters = 6,
ABC_df = abs_AT12_cds_subset_all_gene_ABCs, branchTest_df = relative_abs_AT12_cds_subset_all_gene,
lineage_labels = c("AT1", "AT2"), stretch = T, scaling = T,
norm_method = c("vstExprs", "log"), use_fitting_curves = T,
dist_method = NULL, hclust_method = "ward", heatmap_height = 3,
heatmap_width = 4, ABC_lowest_thrsd = 0, ABC_highest_thrsd = 2,
qval_lowest_thrsd = 1, qval_highest_thrsd = 5, hmcols = NULL,
Cell_type_color = c('#979797', '#F05662', '#7990C8'),
file_name = "genes_branched_heatmap.pdf", pseudo_cnt = 1, cores = 1) {
    new_cds <- buildLineageBranchCellDataSet(cds_subset, stretch = stretch)
    new_cds@dispFitInfo <- cds_subset@dispFitInfo
    
    if(use_fitting_curves) {
        col_gap_ind <- 101
        newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100))
        newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100))
        
        LineageA_exprs <- monocle::responseMatrix(monocle::fitModel(new_cds[, pData(new_cds)$Lineage == 2],
        modelFormulaStr="~sm.ns(Pseudotime, df=3)", cores=cores), newdataA)
        LineageB_exprs <- monocle::responseMatrix(monocle::fitModel(new_cds[, pData(new_cds)$Lineage == 3],
        modelFormulaStr="~sm.ns(Pseudotime, df=3)", cores=cores), newdataB)
        LineageP_num <- 100 - floor(max(pData(new_cds)[pData(new_cds)$State == 1, 'Pseudotime']))
        LineageA_num <- floor(max(pData(new_cds)[pData(new_cds)$State == 1, 'Pseudotime']))
        LineageB_num <- LineageA_num
    }
    else {
        LineageA_exprs <- exprs(new_cds[, pData(new_cds)$Lineage == 2])[, sort(pData(new_cds[, pData(new_cds)$Lineage == 2])$Pseudotime, index.return = T)$ix]
        LineageB_exprs <- exprs(new_cds[, pData(new_cds)$Lineage == 3])[, sort(pData(new_cds[, pData(new_cds)$Lineage == 3])$Pseudotime, index.return = T)$ix]
        
        col_gap_ind <- sum(pData(new_cds)$Lineage == 2) + 1
        
        newdataA <- data.frame(Pseudotime = sort(pData(new_cds[, pData(new_cds)$Lineage == 2])$Pseudotime))
        newdataB <- data.frame(Pseudotime = sort(pData(new_cds[, pData(new_cds)$Lineage == 3])$Pseudotime))
        
        LineageP_num <- sum(pData(new_cds)$State == 1) / 2
        LineageA_num <- sum(pData(new_cds)$State == 2)
        LineageB_num <- sum(pData(new_cds)$State == 3)
    }
    
    if(norm_method == 'vstExprs') {
        LineageA_exprs <- vstExprs(new_cds, expr_matrix=LineageA_exprs)
        LineageB_exprs <- vstExprs(new_cds, expr_matrix=LineageB_exprs)
    }
    else if(norm_method == 'log') {
        LineageA_exprs <- log2(LineageA_exprs + pseudo_cnt)
        LineageB_exprs <- log2(LineageB_exprs + pseudo_cnt)
    }
    
    heatmap_matrix <- cbind(LineageA_exprs[, (col_gap_ind - 1):1], LineageB_exprs)
    heatmap_matrix <- heatmap_matrix[apply(heatmap_matrix, 1, function(x) !any(is.na(x))), ]
    
    if(scaling) {
        heatmap_matrix <- t(scale(t(heatmap_matrix)))
        heatmap_matrix[heatmap_matrix > 3] <- 3
        heatmap_matrix[heatmap_matrix < -3] <- -3
    }
    
    row_dist <- as.dist((1 - cor(t(heatmap_matrix)))/2)
    row_dist[is.na(row_dist)] <- 1
    
    if(is.null(hmcols)) {
        bk <- seq(-3.1,3.1, by=0.1)
        hmcols <- blue2green2red(length(bk) - 1)
    }
    
    # print(hmcols)
    ph <- pheatmap(heatmap_matrix,
    cluster_cols=FALSE,
    cluster_rows=TRUE,
    show_rownames=F,
    show_colnames=F,
    #scale="row",
    clustering_distance_rows=row_dist,
    clustering_method = hclust_method,
    cutree_rows=num_clusters,
    #breaks=bks,
    color=hmcols#,
    # filename="expression_pseudotime_pheatmap.pdf",
    )
    annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
    annotation_row[, "-log 10(qval)"] <- - log10(branchTest_df[row.names(annotation_row), 'qval'])
    annotation_row[which(annotation_row[, "-log 10(qval)"] < qval_lowest_thrsd), "-log 10(qval)"] <- qval_lowest_thrsd
    annotation_row[which(annotation_row[, "-log 10(qval)"] > qval_highest_thrsd), "-log 10(qval)"] <- qval_highest_thrsd
    
    # annotation_row[, "log10(abs(ABCs))"] <- log10(abs(ABC_df[row.names(annotation_row), 'ABCs']))
    # annotation_row[which(annotation_row[, "log10(abs(ABCs))"] < ABC_lowest_thrsd), "log10(abs(ABCs))"] <- ABC_lowest_thrsd
    # annotation_row[which(annotation_row[, "log10(abs(ABCs))"] > ABC_highest_thrsd), "log10(abs(ABCs))"] <- ABC_highest_thrsd
    
    # pData(AT12_cds_subset_all_gene)$cell_type <- "Progenitor"
    # pData(AT12_cds_subset_all_gene)$cell_type[pData(AT12_cds_subset_all_gene)$State == 2] <- lineage_labels[1]
    # pData(AT12_cds_subset_all_gene)$cell_type[pData(AT12_cds_subset_all_gene)$State == 3] <- lineage_labels[2]
    # AT1_num <- (pData(AT12_cds_subset_all_gene)$State == 2)
    
    colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
    annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)), "Cell Type" = c(rep(lineage_labels[1], LineageA_num),
    rep("Progenitor",  2 * LineageP_num),
    rep(lineage_labels[2], LineageB_num)))
    
    colnames(annotation_col) <- "Cell Type"
    
    Cluster_color <- brewer.pal(length(unique(annotation_row$Cluster)),"Set1")
    names(Cluster_color) <- 1:length(unique(annotation_row$Cluster))
    annotation_colors=list("Cell Type"=c(Progenitor=Cell_type_color[1], AT1=Cell_type_color[2], AT2=Cell_type_color[3]),
    'Cluster' = Cluster_color
    )
    names(annotation_colors$`Cell Type`) = c('Progenitor', lineage_labels)
    names(annotation_colors$`Cell Type`) <- c('Progenitor', lineage_labels)

    # pdf(paste(elife_directory, 'AT2_branch_gene_str_norm_div_df_heatmap_cole.pdf', sep = ''))#, height = 4, width = 3)
    # save(heatmap_matrix, hmcols, annotation_row, annotation_col, annotation_colors, row_dist, hclust_method, num_clusters, col_gap_ind, file = 'heatmap_matrix')
    pdf(file_name, height = heatmap_height, width = heatmap_width)
    pheatmap(heatmap_matrix[, ], #ph$tree_row$order
    cluster_cols=FALSE,
    cluster_rows=TRUE,
    show_rownames=F,
    show_colnames=F,
    #scale="row",
    clustering_distance_rows=row_dist, #row_dist
    clustering_method = hclust_method, #ward.D2
    cutree_rows=num_clusters,
    cutree_cols = 2,
    annotation_row=annotation_row,
    annotation_col=annotation_col,
    annotation_colors=annotation_colors,
    gaps_col = col_gap_ind,
    treeheight_row = 1.5,
    #breaks=bks,
    fontsize = 6,
    color=hmcols
    # filename="expression_pseudotime_pheatmap.pdf",
    )
    dev.off()
    
    return(list(LineageA_exprs = LineageA_exprs, LineageB_exprs = LineageB_exprs, heatmap_matrix = heatmap_matrix, ph = ph, annotation_row = annotation_row, annotation_col = annotation_col))
}

plot_genes_branched_pseudotime2 <- function (cds, lineage_states = c(2, 3), lineage_labels = NULL,
method = "fitting", stretch = TRUE, min_expr = NULL, cell_size = 0.75,
nrow = NULL, ncol = 1, panel_order = NULL, cell_color_by = "State", trajectory_color_by = "State",
fullModelFormulaStr = "~ sm.ns(Pseudotime, df=3) * Lineage", reducedModelFormulaStr = NULL,
label_by_short_name = TRUE, weighted = TRUE, add_ABC = FALSE, add_pval = FALSE, normalize = TRUE, bifurcation_time = NULL, ...)
{
    if (add_ABC) {
        ABCs_df <- calABCs(cds, fullModelFormulaStr = fullModelFormulaStr,
        lineage_states = lineage_states, stretch = stretch,
        weighted = weighted, min_expr = min_expr, lineage_labels = lineage_labels,
        ...)
        fData(cds)[, "ABCs"] <- ABCs_df$ABCs
    }
    if (add_pval) {
        pval_df <- branchTest(cds, fullModelFormulaStr = fullModelFormulaStr,
        reducedModelFormulaStr = "~ sm.ns(Pseudotime, df=3)")
        fData(cds)[, "pval"] <- pval_df[row.names(cds), 'pval']
    }
    if("Lineage" %in% all.vars(terms(as.formula(fullModelFormulaStr)))) {
        cds_subset <- buildLineageBranchCellDataSet(cds = cds, lineage_states = lineage_states,
        lineage_labels = lineage_labels, method = method, stretch = stretch,
        weighted = weighted, ...)
    }
    else {
        cds_subset <- cds
        pData(cds_subset)$Lineage <- pData(cds_subset)$State
    }
    
    if (cds_subset@expressionFamily@vfamily %in% c("zanegbinomialff",
    "negbinomial", "poissonff", "quasipoissonff")) {
        integer_expression <- TRUE
    }
    else {
        integer_expression <- FALSE
    }
    if (integer_expression) {
        CM <- exprs(cds_subset)
        if (normalize)
        CM <- t(t(CM)/sizeFactors(cds_subset))
        cds_exprs <- reshape2::melt(round(CM))
    }
    else {
        cds_exprs <- reshape2::melt(exprs(cds_subset))
    }
    if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    colnames(cds_exprs) <- c("f_id", "Cell", "expression")
    cds_pData <- pData(cds_subset)
    if (add_ABC)
    cds_fData <- ABCs_df
    else cds_fData <- fData(cds_subset)
    cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
    cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
    if (integer_expression) {
        cds_exprs$adjusted_expression <- round(cds_exprs$expression)
    }
    else {
        cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
    }
    if (label_by_short_name == TRUE) {
        if (is.null(cds_exprs$gene_short_name) == FALSE) {
            cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
            cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
        }
        else {
            cds_exprs$feature_label <- cds_exprs$f_id
        }
    }
    else {
        cds_exprs$feature_label <- cds_exprs$f_id
    }
    cds_exprs$feature_label <- as.factor(cds_exprs$feature_label)
    # trend_formula <- paste("adjusted_expression", trend_formula,
    #     sep = "")
    cds_exprs$Lineage <- as.factor(cds_exprs$Lineage)
    # merged_df_with_vgam <- plyr::ddply(cds_exprs, .(feature_label),
    #     function(x) {
    #         fit_res <- tryCatch({
    #             expressionFamily <- cds@expressionFamily
    #             if (expressionFamily@vfamily == "negbinomial") {
    #               if (!is.null(cds@dispFitInfo[["blind"]]$disp_func)) {
    #                 disp_guess <- calulate_NB_dispersion_hint(cds@dispFitInfo[["blind"]]$disp_func,
    #                   x$adjusted_expression)
    #                 if (is.null(disp_guess) == FALSE && disp_guess >
    #                   0 && is.na(disp_guess) == FALSE) {
    #                   expressionFamily <- negbinomial(isize = 1/disp_guess)
    #                 }
    #               }
    #             }
    #             vg <- suppressWarnings(vglm(formula = as.formula(trend_formula),
    #               family = expressionFamily, data = x, maxit = 30,
    #               checkwz = FALSE, trace = T))
    #             if (integer_expression) {
    #               res <- predict(vg, type = "response")
    #               res[res < min_expr] <- min_expr
    #             }
    #             else {
    #               res <- 10^(predict(vg, type = "response"))
    #               res[res < log10(min_expr)] <- log10(min_expr)
    #             }
    #             res
    #         }, error = function(e) {
    #             print("Error!")
    #             print(e)
    #             res <- rep(NA, nrow(x))
    #             res
    #         })
    #         df <- cbind(x, expectation = fit_res)
    
    #         df
    #     })
    
    new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime, Lineage = pData(cds_subset)$Lineage)
    full_model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = fullModelFormulaStr,
    relative_expr = T, pseudocount = 0, new_data = new_data)
    colnames(full_model_expectation) <- colnames(cds_subset)
    
    cds_exprs$full_model_expectation <- apply(cds_exprs,1, function(x) full_model_expectation[x[2], x[1]])
    if(!is.null(reducedModelFormulaStr)){
        reduced_model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = reducedModelFormulaStr,
        relative_expr = T, pseudocount = 0, new_data = new_data)
        colnames(reduced_model_expectation) <- colnames(cds_subset)
        cds_exprs$reduced_model_expectation <- apply(cds_exprs,1, function(x) reduced_model_expectation[x[2], x[1]])
    }
    
    if(!is.null(bifurcation_time)){
        cds_exprs$bifurcation_time <- bifurcation_time[as.vector(cds_exprs$gene_short_name)]
    }
    # merged_df_with_vgam <- cds_exprs
    # save(full_model_expectation, new_data, file = 'full_model_expectation')
    if (method == "loess")
    cds_exprs$expression <- cds_exprs$expression + cds@lowerDetectionLimit
    if (label_by_short_name == TRUE) {
        if (is.null(cds_exprs$gene_short_name) == FALSE) {
            cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
            cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
        }
        else {
            cds_exprs$feature_label <- cds_exprs$f_id
        }
    }
    else {
        cds_exprs$feature_label <- cds_exprs$f_id
    }
    cds_exprs$feature_label <- factor(cds_exprs$feature_label)
    if (is.null(panel_order) == FALSE) {
        cds_exprs$feature_label <- factor(cds_exprs$feature_label,
        levels = panel_order)
    }
    
    save(cds_exprs, cds_subset, file = 'cds_exprs')
    #fix the bug to have NA or 0 values in the expression / expectation values:
    cds_exprs$expression[is.na(cds_exprs$expression)] <- min_expr
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_exprs$full_model_expectation[is.na(cds_exprs$full_model_expectation)] <- min_expr
    cds_exprs$full_model_expectation[cds_exprs$full_model_expectation < min_expr] <- min_expr
    
    if(!is.null(reducedModelFormulaStr)){
        cds_exprs$reduced_model_expectation[is.na(cds_exprs$reduced_model_expectation)] <- min_expr
        cds_exprs$reduced_model_expectation[cds_exprs$reduced_model_expectation < min_expr] <- min_expr
    }
    
    cds_exprs$State <- as.factor(cds_exprs$State)
    cds_exprs$Lineage <- as.factor(cds_exprs$Lineage)
    # q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
    # cds_exprs <- cbind(cds_exprs, cds_exprs)
    # cds_exprs$State <- as.factor(cds_exprs$State)
    # cds_exprs$Lineage <- as.factor(cds_exprs$Lineage)
    q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
    
    if (!is.null(bifurcation_time)){
        q <- q + geom_vline(aes(xintercept = bifurcation_time), color = 'black', linetype = "longdash")
    }
    
    if (is.null(cell_color_by) == FALSE) {
        q <- q + geom_point(aes_string(color = cell_color_by), size = I(cell_size)) + scale_color_discrete(name = "Cell")
    }
    if (add_ABC)
    q <- q + scale_y_log10() + facet_wrap(~feature_label +
    ABCs, nrow = nrow, ncol = ncol, scales = "free_y")
    else if (add_pval)
    q <- q + scale_y_log10() + facet_wrap(~feature_label +
    pval, nrow = nrow, ncol = ncol, scales = "free_y")
    else q <- q + scale_y_log10() + facet_wrap(~feature_label,
    nrow = nrow, ncol = ncol, scales = "free_y")
    save(q, cds_exprs, file = 'test_color_scheme')

    if (method == "loess")
    q <- q + stat_smooth(aes(fill = Lineage, color = Lineage),
    method = "loess")
    else if (method == "fitting") {
        save(q, file = 'q')
        q <- q + geom_line(aes_string(x = "Pseudotime", y = "full_model_expectation",
        color = trajectory_color_by), data = cds_exprs) #+ scale_color_manual(name = "Lineage", values = c(colour_cell, colour), labels = c("Ifnar1_KO_LPS", "Stat1_KO_LPS", "LPS", "On_Chip_Unstim", lineage_labels))
    }
    
    if(!is.null(reducedModelFormulaStr)) {
        q <- q + geom_line(aes_string(x = "Pseudotime", y = "reduced_model_expectation"),
        color = 'black', linetype = 2, data =  cds_exprs)
    }
    if (stretch)
    q <- q + ylab("Expression") + xlab("Maturation levels")
    else q <- q + ylab("Expression") + xlab("Pseudotime")
    q <- q + monocle_theme_opts()
    # q + expand_limits(y = min_expr)
}
