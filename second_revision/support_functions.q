so_switch_test <- function(isoform_cds, cores=1){

  platform <- Sys.info()[['sysname']]
  if (platform == "Windows")
    cl <- makeCluster(cores)
  if (platform %in% c("Linux", "Darwin")) 
    cl <- makeCluster(cores)
  
  cleanup <- function(){
    stopCluster(cl)
  }
  on.exit(cleanup)
  
  required_packages <- c("BiocGenerics", "Biobase", "VGAM", "plyr")
  if (is.null(required_packages) == FALSE){
    clusterCall(cl, function(pkgs) {
      for (req in pkgs) {
        library(req, character.only=TRUE)
      }
    }, required_packages)
  }

  gene_ids <- unique(fData(isoform_cds)$gene_id)
  res <- parLapply(cl, gene_ids, function(gene, isoform_cds){
  #print (head(iso_exprs))
  fd_isoforms <- subset(fData(isoform_cds), gene_id == gene)
  #print (head(iso_exprs))
  tryCatch({
    fd_isoforms <- subset(fd_isoforms, num_cells_expressed >= 15)
    if (nrow(fd_isoforms) > 1){
      iso_names <- row.names(fd_isoforms)
      #print (iso_names)
      iso_exprs <- exprs(isoform_cds[iso_names,])
      
      iso_exprs <- round(t(iso_exprs) / pData(isoform_cds)$Size_Factor)

      #iso_exprs <- iso_exprs[,colSums(iso_exprs > 1) >= 15]
      #print (head(iso_exprs))

      iso_names <- colnames(iso_exprs)
      iso_exprs <- cbind(iso_exprs, pData(isoform_cds))

      iso_exprs <- iso_exprs[rowSums(iso_exprs[,iso_names]) > 0,]

      # fit_time <- vglm(as.formula(paste("cbind(",toString(iso_names),") ~ sm.ns(Pseudotime, df=3)", sep="")), dirmultinomial(lphi="logit"), data = iso_exprs, epsilon=0.1)
      # fit_null<- vglm(as.formula(paste("cbind(",toString(iso_names),") ~ 1", sep="")), dirmultinomial(lphi="logit"), data = iso_exprs, epsilon=0.1)
      # link_phi <- predict(fit_time)
      # link_phi <- unique(link_phi[,ncol(link_phi)])[1]

      fit_time <- vglm(as.formula(paste("cbind(",toString(iso_names),") ~ sm.ns(Pseudotime, df=3)", sep="")), dirmultinomial(), data = iso_exprs, epsilon=0.1)
      fit_null<- vglm(as.formula(paste("cbind(",toString(iso_names),") ~ 1", sep="")), dirmultinomial(), data = iso_exprs, epsilon=0.1)


      #print(summary(fit_time))
      
      lrt <- lrtest(fit_time,fit_null) 
      #print(lrt)
      pval=lrt@Body["Pr(>Chisq)"][2,]


      #print(pval)
      test_res <- data.frame(link_phi=NA, status = "OK", pval=pval)
    }else{
      test_res <- data.frame(link_phi=NA, status = "NOTEST", pval=1.0)
    }
  }, error = function(e) { print (e); data.frame(link_phi=NA, status = "LOWDATA", pval=1.0)})
  #print(test_res)
  }, isoform_cds)
  res <- do.call(rbind, res)
  row.names(res) <- gene_ids
  res
}

compare_cell_types_in_pseudotime <-function(cds_subset, 
                                            trend_formula="~ sm.ns(Pseudotime, df=3)*Cell.Type",
                      										  min_expr=NULL, 
                      											cell_size=0.75, 
                      											cell_alpha=1.0, 
                      											nrow=NULL, 
                      											ncol=1, 
                      											panel_order=NULL, 
                      											color_by="Cell.Type", 
                      											shade_by=NULL, 
                      											df=3, 
                      											maxit=300,
                      											relative_expr=TRUE, 
                                            pseudocount=0){
  
  cds_pData <- pData(cds_subset)
  cds_fData <- fData(cds_subset)
    
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial")){
    integer_expression <- TRUE
  }else{
    integer_expression <- FALSE
	  relative_expr <- TRUE
  }
  
  if (integer_expression)
  {
    cds_exprs <- exprs(cds_subset)
    cds_exprs <- cds_exprs + pseudocount
    if (relative_expr){
      if (is.null(sizeFactors(cds_subset)))
      {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs <- t(t(cds_exprs) / sizeFactors(cds_subset))
    }
    cds_exprs <- melt(ceiling(cds_exprs))
  }else{
    cds_exprs <- exprs(cds_subset)
    cds_exprs <- cds_exprs + pseudocount
    cds_exprs <- melt(cds_exprs)
  }

  
  
  colnames(cds_exprs) <- c("gene_id", "Cell", "expression")
 
  cds_exprs <- merge(cds_exprs, cds_fData, by.x="gene_id", by.y="row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x="Cell", by.y="row.names")
 
  if (is.null(min_expr)){
    min_expr <- cds_subset@lowerDetectionLimit
  }
  
  if (integer_expression)
  {
    cds_exprs$expression <- cds_exprs$expression
  }else{
    cds_exprs$expression <- log10(cds_exprs$expression)
  }
  
  if (is.null(cds_exprs$gene_short_name) == FALSE){
    cds_exprs$gene_label <- as.character(cds_exprs$gene_short_name)
    cds_exprs$gene_label[is.na(cds_exprs$gene_label)]  <- cds_exprs$gene_id
  }else{
    cds_exprs$gene_label <- cds_exprs$gene_id
  }
  cds_exprs$gene_label <- factor(cds_exprs$gene_label)
  
  # merged_df_with_vgam <- ddply(cds_exprs, .(gene_label), function(x) { 
  #   fit_res <- tryCatch({
  #     #Extra <- list(leftcensored = with(x, fpkm <= min_fpkm), rightcencored = rep(FALSE, nrow(x)))
  #     #vg <- vgam(formula = fpkm ~ s(pseudo_time), family = cennormal, data = x, extra=Extra, maxit=30, trace = TRUE) 
  #     vg <- vglm(formula = expression ~ sm.ns(Pseudotime, df=df) * Cell.Type, family = negbinomial(), data = x, maxit=maxit, checkwz=F) 
  #     if (integer_expression){
  #       res <- predict(vg, type="response")
  #       res[res < min_expr] <- min_expr
  #     }else{
  #       res <- 10^(predict(vg, type="response"))
  #       res[res < log10(min_expr)] <- log10(min_expr)
  #     }
  #     res
  #   }
  #   ,error = function(e) {
  #   					traceback()
  #   					print(e)
  #   					res <- rep(NA, nrow(x))
  #   					res
  #   			}
  #   )
    
  #   expectation = fit_res

  #   data.frame(Pseudotime=x$Pseudotime, expectation, Cell.Type=x$Cell.Type)
  # })
  
  new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime, Cell.Type=pData(cds_subset)$Cell.Type)
  model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = trend_formula,
                      relative_expr = T, pseudocount = 0, new_data = new_data, weights = pData(cds_subset)$weight)
  colnames(model_expectation) <- colnames(cds_subset)

  cds_exprs$expectation <- apply(cds_exprs,1, function(x) model_expectation[x[2], x[1]])

  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr

  if (is.null(panel_order) == FALSE)
  {
    cds_subset$gene_label <- factor(cds_subset$gene_label, levels=panel_order)
  }
  
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  
  q <- ggplot(aes(Pseudotime, expression), data=cds_exprs) 
  
  if (is.null(color_by) == FALSE){
    q <- q + geom_point(aes_string(color=color_by), size=I(cell_size), alpha=I(cell_alpha))
  }else{
    q <- q + geom_point(size=I(cell_size))
  }
    
  q <- q + geom_line(aes(Pseudotime, expectation, color=Cell.Type))
  
  #q <- q + geom_ribbon(aes(x=pseudo_time, ymin=conf_lo, ymax=conf_hi), alpha=I(0.15), data=merged_df_with_vgam)
  
  q <- q + scale_y_log10() + facet_wrap(~gene_label, nrow=nrow, ncol=ncol, scales="free_y")
  
  q <- q + ylab("Expression") + xlab("Pseudo-time")
  q
  
}

# qp <- compare_cell_types_in_pseudotime(combined_cds[row.names(subset(fData(combined_cds),
# 								 gene_short_name %in% missing_myo_genes
# 								 )),], color_by="Cell.Type", ncol=3, df=4) + scale_color_brewer(palette="Set1")
# pdf("missing_myogenesis_genes.pdf", width=6, height=8)
# print (qp)
# dev.off()
#

# We're gonna want to use this implementation of Mahalanobis distance for various clustering applications
fastPwMahal = function(x1) {
    require(corpcor)
    #SQRT = with(irlba(invCovMat, min(dc)-1, min(dc)-1), u %*% diag(d^0.5) %*% t(v))
    SQRT = with(fast.svd(ginv(cov(x1))), u %*% diag(d^0.5) %*% t(v))
    dist(x1 %*% SQRT)
}


loadGSCSafe <- function (file, type = "auto", addInfo, sep="\t", encoding="latin1") 
{
    if (missing(addInfo)) {
        addUserInfo <- "skip"
        addInfo <- "none"
    }
    else {
        addUserInfo <- "yes"
    }
    tmp <- try(type <- match.arg(type, c("auto", "gmt", "sbml", 
        "sif", "data.frame"), several.ok = FALSE), silent = TRUE)
    if (class(tmp) == "try-error") {
        stop("argument type set to unknown value")
    }
    if (type == "auto") {
        if (class(file) == "character") {
            tmp <- unlist(strsplit(file, "\\."))
            type <- tolower(tmp[length(tmp)])
            if (!type %in% c("gmt", "sif", "sbml", "xml")) 
                stop(paste("can not handle .", type, " file extension, read manually using e.g. read.delim() and load as data.frame", 
                  sep = ""))
        }
        else {
            type <- "data.frame"
        }
    }
    if (type == "gmt") {
        con <- file(file, encoding=encoding)
        tmp <- try(suppressWarnings(open(con)), silent = TRUE)
        if (class(tmp) == "try-error") 
            stop("file could not be read")
        if (addUserInfo == "skip") 
            addInfo <- vector()
        gscList <- list()
        i <- 1
        tmp <- try(suppressWarnings(while (length(l <- scan(con, 
            nlines = 1, what = "character", quiet = T, sep=sep)) > 0) {
            if (addUserInfo == "skip") 
                addInfo <- rbind(addInfo, l[1:2])
            tmp <- l[3:length(l)]
            gscList[[l[1]]] <- unique(tmp[tmp != "" & tmp != 
                " " & !is.na(tmp)])
            i <- i + 1
        }), silent = TRUE)
        if (class(tmp) == "try-error") 
            stop("file could not be read")
        close(con)
        gsc <- gscList[!duplicated(names(gscList))]
        if (addUserInfo == "skip") 
            addInfo <- unique(addInfo)
    }
    else if (type %in% c("sbml", "xml")) {
        require(rsbml)
        tmp <- try(sbml <- rsbml_read(file))
        if (class(tmp) == "try-error") {
            stop("file could not be read by rsbml_read()")
        }
        gsc <- list()
        for (iReaction in 1:length(reactions(model(sbml)))) {
            metIDs <- names(c(reactants(reactions(model(sbml))[[iReaction]]), 
                products(reactions(model(sbml))[[iReaction]])))
            geneIDs <- names(modifiers(reactions(model(sbml))[[iReaction]]))
            if (length(geneIDs) > 0) {
                geneNames <- rep(NA, length(geneIDs))
                for (iGene in 1:length(geneIDs)) {
                  geneNames[iGene] <- name(species(model(sbml))[[geneIDs[iGene]]])
                }
                for (iMet in 1:length(metIDs)) {
                  gsc[[metIDs[iMet]]] <- c(gsc[[metIDs[iMet]]], 
                    geneNames)
                }
            }
        }
        if (length(gsc) == 0) {
            stop("no gene association found")
        }
        else {
            for (iMet in 1:length(gsc)) {
                tmp1 <- name(species(model(sbml))[[names(gsc)[iMet]]])
                tmp2 <- compartment(species(model(sbml))[[names(gsc)[iMet]]])
                names(gsc)[iMet] <- paste(tmp1, " (", tmp2, ")", 
                  sep = "")
            }
        }
    }
    else if (type == "sif") {
        tmp <- try(gsc <- as.data.frame(read.delim(file, header = FALSE, 
            quote = "", as.is = TRUE), stringsAsFactors = FALSE), 
            silent = TRUE)
        if (class(tmp) == "try-error") {
            stop("argument file could not be read and converted into a data.frame")
        }
        if (ncol(gsc) != 3) {
            stop("sif file should contain three columns")
        }
        if (addUserInfo == "skip") 
            addInfo <- gsc[, c(1, 2)]
        gsc <- gsc[, c(3, 1)]
        tmp <- nrow(gsc)
        gsc <- unique(gsc)
        geneSets <- unique(gsc[, 2])
        gscList <- list()
        for (iGeneSet in 1:length(geneSets)) {
            gscList[[iGeneSet]] <- gsc[gsc[, 2] == geneSets[iGeneSet], 
                1]
        }
        names(gscList) <- geneSets
        gsc <- gscList
    }
    else if (type == "data.frame") {
        tmp <- try(gsc <- as.data.frame(file, stringsAsFactors = FALSE), 
            silent = TRUE)
        if (class(tmp) == "try-error") {
            stop("argument file could not be converted into a data.frame")
        }
        for (i in 1:ncol(gsc)) {
            gsc[, i] <- as.character(gsc[, i])
        }
        if (ncol(gsc) != 2) {
            stop("argument file has to contain exactly two columns")
        }
        tmp <- nrow(gsc)
        gsc <- unique(gsc)
        geneSets <- unique(gsc[, 2])
        gscList <- list()
        for (iGeneSet in 1:length(geneSets)) {
            gscList[[iGeneSet]] <- gsc[gsc[, 2] == geneSets[iGeneSet], 
                1]
        }
        names(gscList) <- geneSets
        gsc <- gscList
    }
    if (addUserInfo == "yes") {
        tmp <- try(addInfo <- as.data.frame(addInfo, stringsAsFactors = FALSE), 
            silent = TRUE)
        if (class(tmp) == "try-error") {
            stop("failed to convert additional info in argument 'addInfo' into a data.frame")
        }
    }
    if (class(addInfo) == "data.frame") {
        if (ncol(addInfo) != 2) 
            stop("additional info in argument 'file' or 'addInfo' has to contain 2 columns")
        tmp <- nrow(addInfo)
        addInfo <- unique(addInfo[addInfo[, 1] %in% names(gsc), 
            ])
    }
    else {
    }
    res <- list(gsc, addInfo)
    names(res) <- c("gsc", "addInfo")
    class(res) <- "GSC"
    return(res)
}



plot_gsa_hyper_heatmap <- function(cds, gsa_results, significance=0.05)
{
	hyper_df <- ldply(gsa_results, function(gsa_res)
		{
			data.frame(gene_set = names(gsa_res$pvalues), pval = gsa_res$pvalues, qval = gsa_res$p.adj)
		})
	colnames(hyper_df)[1] <- "cluster_id"
	
	#hyper_df 
	
	hyper_df <- subset(hyper_df, qval <= significance)
	print (head(hyper_df))
	hyper_df <- merge(hyper_df, ddply(hyper_df, .(gene_set), function(x) { nrow(x) }), by="gene_set")
	#print (hyper_df)
	hyper_df$gene_set <- factor(hyper_df$gene_set, levels=unique(arrange(hyper_df, V1, cluster_id)$gene_set))
	
	qplot(cluster_id, gene_set, fill=-log10(qval), geom="tile", data=hyper_df) + scale_fill_gradientn(colours=rainbow(7))
}



collect_gsa_hyper_results <- function(cds, gsc, clusters, gsSizeLim=c(1, Inf))
{
    gene_universe <- unique(as.character(fData(cds)$gene_short_name))
    gsa_results <- list()
    cluster_ids <- unique(clusters)
    for (i in (1:length(cluster_ids))) {
        cluster_genes <- unique(fData(cds[names(clusters[clusters == i]),])$gene_short_name)
        gsaRes <- runGSAhyper(cluster_genes, gsc=gsc, universe=gene_universe, gsSizeLim=gsSizeLim)
        gsa_results[[length(gsa_results) + 1]] <- gsaRes
    }
    names(gsa_results) <- cluster_ids
    gsa_results
}

filter_gsc <- function(gsc, cds, alias_table, allowed_genes, gsc_whitelist)
{
	gene_gsc_names_allowed <- sapply(names(gsc$gsc), function(gsc_name, cds, alias_table, allowed_genes) { 
	
		allowed_gene_symbols <- fData(cds[allowed_genes,])$gene_short_name
	
		gsc_matches_to_official_symbol <- sum(grepl(paste("^",gsc_name, sep=""), allowed_gene_symbols))
		gsc_matches_to_alias <- sum(grepl(paste("^", gsc_name, sep=""), subset(alias_table, symbol %in% allowed_gene_symbols)$alias_symbol))
		allowed_gsc_names <- gsc_matches_to_official_symbol + gsc_matches_to_alias > 0
		allowed_gsc_names
		},
		cds, alias_table, allowed_genes
		)
	gsc$gsc <- gsc$gsc[union(names(gene_gsc_names_allowed)[gene_gsc_names_allowed], gsc_whitelist)]
	gsc
}

clusterGenesRelative<-function(expr_matrix, k, method=function(x){as.dist((1 - cor(t(x)))/2)}, ...){
  n<-method(expr_matrix)
	  	  	
  clusters<-pam(n,k, ...)
  class(clusters)<-"list"
  clusters$exprs<-expr_matrix
  clusters
}

plot_clusters_relative<-function(cds, 
                       	 		 clustering,
                        		 drawSummary=TRUE, 
                        		 sumFun=mean_cl_normal,
                        		 ncol=NULL, 
                        		 nrow=NULL, 
                        		 row_samples=NULL, 
                        		 callout_ids=NULL,
                                 conf_int=0.68,
                                 draw_individual_genes=FALSE){
  m <- as.data.frame(clustering$exprs)

  baseline_expr <- mean(rowMeans(m))

  m$ids <- rownames(clustering$exprs)
  if (is.null(clustering$labels) == FALSE)
  {
    m$cluster = factor(clustering$labels[clustering$clustering], levels = levels(clustering$labels))
  }else{
    m$cluster <- factor(clustering$clustering)
  }
  
  cluster_sizes <- as.data.frame(table(m$cluster))    
  
  cluster_sizes$Freq <- paste("(", cluster_sizes$Freq, ")")   
  facet_labels <- str_join(cluster_sizes$Var1, cluster_sizes$Freq, sep=" ")
 
  facet_wrap_labeller <- function(gg.plot,labels=NULL) {
    #works with R 3.0.1 and ggplot2 0.9.3.1
    require(gridExtra)

    g <- ggplotGrob(gg.plot)
    gg <- g$grobs      
    strips <- grep("strip_t", names(gg))

    for(ii in seq_along(labels))  {
      modgrob <- getGrob(gg[[strips[ii]]], "strip.text", 
                         grep=TRUE, global=TRUE)
      gg[[strips[ii]]]$children[[modgrob$name]] <- editGrob(modgrob,label=labels[ii])
    }

    g$grobs <- gg
    class(g) = c("arrange", "ggplot",class(g)) 
    g
  }
  
  
  m.melt <- melt(m, id.vars = c("ids", "cluster"))

  m.melt <- merge(m.melt, pData(cds), by.x="variable", by.y="row.names")
 
  
  if (is.null(row_samples) == FALSE){
    m.melt <- m.melt[sample(nrow(m.melt), row_samples),]
  }
  
  c <- ggplot(m.melt) + facet_wrap("cluster", ncol=ncol, nrow=nrow, scales="free_y")
  #c <- c + stat_density2d(aes(x = Pseudotime, y = value), geom="polygon", fill="white", color="black", size=I(0.1)) + facet_wrap("cluster", ncol=ncol, nrow=nrow)

  if (draw_individual_genes){
      c <- c + geom_line(aes(x=Pseudotime, y=value, group=ids), alpha=I(0.1), size=I(0.5))
  }
  
  if (drawSummary) {
    c <- c + stat_summary(aes(x = Pseudotime, y = value, group = 1),
                          fun.data = sumFun, color = "red",
                          alpha = 0.2, size = 0.5, geom = "smooth")
    c <- c + stat_summary(aes(x = Pseudotime, y = value, group = 1), 
                          fun.data = "median_hilow", fill = "black", 
                          alpha = 0.2, size = 0.5, geom = "ribbon", conf.int=conf_int)
  }
  
  #cluster_medians <- subset(m.melt, ids %in% clustering$medoids)
  
  #c <- c + geom_line()
  #c <- c + geom_line(aes(x=Pseudotime, y=value), data=cluster_medians, color=I("red"))
  c <- c + scale_color_hue(l = 50, h.start = 200) + theme(axis.text.x = element_text(angle = 0, 
                                                                                     hjust = 0)) + xlab("Pseudo-time") + ylab("Expression")
  c <- c + theme(strip.background = element_rect(colour = 'white', fill = 'white')) + 
    theme(panel.border = element_blank()) +
    theme(legend.position="none") +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank())
  	
  c <- c + geom_hline(yintercept=baseline_expr, color="steelblue")
#   if (draw_cluster_size){
#     cluster_sizes <- as.data.frame(table(m$cluster))
#     colnames(cluster_sizes) <- c("cluster", "Freq")
#     cluster_sizes <- cbind (cluster_sizes, Pseudotime = cluster_label_text_x, value = cluster_label_text_y)
#     c <- c + geom_text(aes(x=Pseudotime, y=value, label=Freq), data=cluster_sizes, size=cluster_label_text_size)
#   }
  
  if (is.null(callout_ids) == FALSE)
  {
    callout_melt <- subset(m.melt, ids %in% callout_ids)
    c <- c + geom_line(aes(x=Pseudotime, y=value), data=callout_melt, color=I("steelblue"))
  }
  #c <- c + monocle_theme_opts()
  #c <- facet_wrap_labeller(c, facet_labels)
  c
}

# selectHighDispersionGenes <- function(cds, 
#                                       detectionLimit=-Inf, 
#                                       relative_expr=TRUE){
#   disp_df<-esApply(cds,1,
#     function(f_expression) { 
#       if (relative_expr && cds@expressionFamily@vfamily == "negbinomial"){
#         f_expression <- f_expression / Size_Factor
#       }
#       f_expression <- f_expression[f_expression > detectionLimit]
#       expr_mean <- mean(f_expression[f_expression > 0])
#       if (is.null(expr_mean) == FALSE) {
#         disp_guess_fit <- cds@dispFitInfo[["blind"]]$disp_func(expr_mean)
        
#         # For NB: Var(Y)=mu*(1+mu/k)
#         f_expression_var <- var(f_expression)
#         f_expression_mean <- mean(f_expression)
        
#         disp_guess_meth_moments <- f_expression_var - f_expression_mean 
#         disp_guess_meth_moments <- disp_guess_meth_moments / (f_expression_mean^2) #fix the calculation of k 
        
#         return (data.frame(mean_exp=expr_mean, disp_fit=disp_guess_fit, disp_empirical=disp_guess_meth_moments))
#       }
#       return (NULL)
#   } )
#   do.call(rbind,disp_df)
# }

# NOTE: our definition of "reference" and "query" are reversed w.r.t what DTW uses.
align_cells <- function(query_cell_matrix, ref_cell_matrix, 
            step_pattern=rabinerJuangStepPattern(3, "d", smoothed=TRUE),
            open.end=FALSE,
            open.begin=FALSE)
{
    # we want the query and reference to have the same number of cells, so downsample the larger
    # one to match the smaller one
    
    
    if (ncol(ref_cell_matrix) < ncol(query_cell_matrix)){
        ref_cell_matrix_sub <- ref_cell_matrix
        sampled_idxs <- c(1, sort(sample(seq(2,ncol(query_cell_matrix)-1), ncol(ref_cell_matrix)-2)), ncol(query_cell_matrix))
        query_cell_matrix_sub <- query_cell_matrix[,sampled_idxs]
    }else if (ncol(ref_cell_matrix) > ncol(query_cell_matrix)){
        query_cell_matrix_sub <- query_cell_matrix
        sampled_idxs <- c(1, sort(sample(seq(2,ncol(ref_cell_matrix)-1), ncol(query_cell_matrix)-2)), ncol(ref_cell_matrix))
        ref_cell_matrix_sub <- ref_cell_matrix[,sampled_idxs]
    }else{
        query_cell_matrix_sub <- query_cell_matrix
        ref_cell_matrix_sub <- ref_cell_matrix
    }

    ref_cell_matrix_sub <- t(scale(t(ref_cell_matrix_sub)))
    query_cell_matrix_sub <- t(scale(t(query_cell_matrix_sub)))
    
  valid_genes <- rowSums(is.na(ref_cell_matrix_sub)) == 0 & rowSums(is.na(query_cell_matrix_sub)) == 0
  ref_cell_matrix_sub <- ref_cell_matrix_sub[valid_genes,]
  query_cell_matrix_sub <- query_cell_matrix_sub[valid_genes,]
    
  print (dim(ref_cell_matrix_sub))
    #dist_matrix <- (1 - cor(ref_cell_matrix_sub, query_cell_matrix_sub, use="pairwise.complete.obs"))^2
    dist_matrix <- (1 - cor(ref_cell_matrix_sub, query_cell_matrix_sub, use="pairwise.complete.obs"))
    
  cell_alignment_dtw <- dtw(dist_matrix,
                    step.pattern=step_pattern, 
                    open.end=open.end,
                    open.begin=open.begin,
                    keep=T,
                              keep.internals=T)
  cell_alignment_dtw
}


customDtwPlotDensity <- function (d, normalize = FALSE, xlab = "Query index", ylab = "Reference index", 
    ...) 
{
    cm <- d$costMatrix
    if (is.null(cm)) 
        stop("dtwPlotDensity requires dtw internals (set keep.internals=TRUE on dtw() call)")
    if (normalize) {
        norm <- attr(d$stepPattern, "norm")
        if (is.na(norm)) 
            stop("No normalization known for step pattern used")
        if (norm == "N") {
            cm <- cm/row(cm)
        }
        else if (norm == "N+M") {
            cm <- cm/(row(cm) + col(cm))
        }
        else if (norm == "M") {
            cm <- cm/col(cm)
        }
    }
    xd <- dim(cm)[1]
    yd <- dim(cm)[2]
    image(cm, col = terrain.colors(100), x = 1:xd, y = 1:yd, 
        xlab = xlab, ylab = ylab, useRaster=TRUE, ...)
    contour(cm, x = 1:xd, y = 1:yd, add = TRUE)
    lines(d$index1, d$index2, col = "cyan", lwd = 2)
}

warp_pseudotime <- function(ref_cds, query_cds, alignment_dtw)
{
    ref_pseudotimes <- pData(ref_cds)$Pseudotime
    names(ref_pseudotimes) <- row.names(pData(ref_cds))
    ref_pseudotimes <- sort(ref_pseudotimes)
    
    query_pseudotimes <- pData(query_cds)$Pseudotime
    names(query_pseudotimes) <- row.names(pData(query_cds))
    query_pseudotimes <- sort(query_pseudotimes)
    
    query_alignment_time <- ref_pseudotimes[warp(alignment_dtw, index.reference=T)]
    
    #print (query_alignment_time)
    pData(ref_cds)$Alignment_Pseudotime <- pData(ref_cds)$Pseudotime
    pData(query_cds)$Alignment_Pseudotime <- NA
  
    #names(query_alignment_time) <- rownames(alignment_dtw$localCostMatrix)
  
    #pData(query_cds)[names(query_alignment_time),]$Alignment_Pseudotime <- query_alignment_time
    
    fc <- approxfun(seq(1, 101, by=1), query_alignment_time, rule=2)
    pData(query_cds)$Alignment_Pseudotime  <- fc(query_cds$Pseudotime)
    return(list(ref_cds=ref_cds, query_cds=query_cds))
    
    # BJ_MYO_aligned <- BJ_MYO_selected
    # pData(BJ_MYO_aligned)[colnames(smoothed_BJ_MYO_exprs),]$Pseudotime <- BJ_MYO_alignment_time
    # pData(BJ_MYO_selected)[colnames(smoothed_BJ_MYO_exprs),]$Pseudotime <- BJ_MYO_alignment_time
}

calABCs_from_alignment <- function(cds, 
                                   trend_formula = "~sm.ns(Pseudotime, df = 3)*Cell.Type",
                                   method = 'fitting', 
                                   ABC_method = 'integral', 
                                   points_num = 1000, 
                                   fc_limit = 3, 
                                   branchTest = FALSE, 
                                   cell_types = c("BJ-MYO", "HSMM"), 
                                   relative_expr = TRUE,
                                   stretch = TRUE,
                                   pseudocount=0,
                                   cores = 1, 
                                   weighted = TRUE, 
                                   min_expr = 0.5, 
                                   integer_expression = FALSE, 
                                   num = 100, 
                                   cell_type_labels = NULL, ...) {


  new_data = list()
  for (cell_type in unique(pData(cds)$Cell.Type)){
    new_data[[length(new_data) + 1]] = data.frame(Pseudotime=seq(1, 100, length.out=num),
                                                  Cell.Type=cell_type)
  }
  new_data = do.call(rbind, new_data)
  
  str_branchAB_expression_curve_matrix <- genSmoothCurves(cds, 
                                                          cores=cores, trend_formula = trend_formula,
                                                          relative_expr = relative_expr, 
                                                          pseudocount = pseudocount, 
                                                          new_data = new_data)
  
  str_branchA_expression_curve_matrix <- str_branchAB_expression_curve_matrix[, 1:num]
  str_branchB_expression_curve_matrix <- str_branchAB_expression_curve_matrix[, (num + 1):(2 * num)]
  
  ABCs_res <- str_branchA_expression_curve_matrix - str_branchB_expression_curve_matrix
  ILR_res <- log2(str_branchA_expression_curve_matrix / (str_branchB_expression_curve_matrix + 0.1))
  
  ABCs_res <- apply(ABCs_res, 1, function(x, num, ABC_method) {
    avg_delta_x <- (x[1:(num - 1)] + x[2:(num)])/2
    step <- (100/(num - 1))
    
    if (ABC_method == "integral") {
      res <- round(sum(avg_delta_x * step), 3)
    }
    else if (ABC_method == "global_normalization") {
      max <- max(max(predictBranchOri), max(x))
      res <- round(sum(avg_delta_x/max * step), 3)
    }
    else if (ABC_method == "local_normalization") {
      pair_wise_max <- apply(data.frame(x = x, y = predictBranchOri),
                             1, max)
      res <- round(sum((((predictBranchOri - x)/pair_wise_max)[1:(num -
                                                                    1)] + ((predictBranchOri - x)/pair_wise_max)[2:(num)])/2 *
                         step), 3)
    }
    else if (ABC_method == "four_values") {
      ori_ABCs <- round(sum((x[1:(num - 1)] + x[2:(num)])/2 *
                              step), 3)
      other_ABCs <- round(sum((predictBranchOri[1:(num -
                                                     1)] + predictBranchOri[2:(num)])/2 * step),
                          3)
      ori_ABCs_H <- round(sum(avg_delta_x[avg_delta_x >
                                            0] * step), 3)
      other_ABCs_H <- round(sum(avg_delta_x[avg_delta_x <
                                              0] * step), 3)
      res <- c(ori_ABCs = ori_ABCs, other_ABCs = other_ABCs,
               ori_ABCs_H = ori_ABCs_H, other_ABCs_H = other_ABCs_H)
    }
    else if (ABC_method == "ILRs") {
      str_logfc_df <- log2((predictBranchOri + 1)/(x +
                                                     1))
      res <- sum(str_logfc_df)
    }
    return(res)}, num = num, ABC_method = ABC_method
  )
  
  ABCs_res <- cbind(ABCs_res, ILR_res[,ncol(ILR_res)])
  colnames(ABCs_res)<- c("ABCs", "Endpoint_ILR")

  ABCs_res <- merge(ABCs_res, fData(cds), by = "row.names")
  row.names(ABCs_res) <- ABCs_res[, 1]
  ABCs_res[, 1] <- NULL
  colnames(ABCs_res)[1] <- "ABCs"
  
  
  return(ABCs_res)
}


#The following code is swipped from colorRamps package which is used to make the pallette
table.ramp <- function(n, mid = 0.5, sill = 0.5, base = 1, height = 1)
{
    x <- seq(0, 1, length.out = n)
    y <- rep(0, length(x))
    sill.min <- max(c(1, round((n - 1) * (mid - sill / 2)) + 1))
    sill.max <- min(c(n, round((n - 1) * (mid + sill / 2)) + 1))
    y[sill.min:sill.max] <- 1
    base.min <- round((n - 1) * (mid - base / 2)) + 1
    base.max <- round((n - 1) * (mid + base / 2)) + 1
    xi <- base.min:sill.min
    yi <- seq(0, 1, length.out = length(xi))
    i <- which(xi > 0 & xi <= n)
    y[xi[i]] <- yi[i]
    xi <- sill.max:base.max
    yi <- seq(1, 0, length.out = length(xi))
    i <- which(xi > 0 & xi <= n)
    y[xi[i]] <- yi[i]
    height * y
}


rgb.tables <- function(n,
red = c(0.75, 0.25, 1),
green = c(0.5, 0.25, 1),
blue = c(0.25, 0.25, 1))
{
    rr <- do.call("table.ramp", as.list(c(n, red)))
    gr <- do.call("table.ramp", as.list(c(n, green)))
    br <- do.call("table.ramp", as.list(c(n, blue)))
    rgb(rr, gr, br)
}

matlab.like <- function(n) rgb.tables(n)

matlab.like2 <- function(n)
rgb.tables(n,
red = c(0.8, 0.2, 1),
green = c(0.5, 0.4, 0.8),
blue = c(0.2, 0.2, 1))

blue2green2red <- matlab.like2


#beam test on the isoform switch: 
beam_iso_switch_test <- function(isoform_cds, cores=1){

  platform <- Sys.info()[['sysname']]
  if (platform == "Windows")
    cl <- makeCluster(cores)
  if (platform %in% c("Linux", "Darwin")) 
    cl <- makeCluster(cores)
  
  cleanup <- function(){
    stopCluster(cl)
  }
  on.exit(cleanup)
  
  required_packages <- c("BiocGenerics", "Biobase", "VGAM", "plyr")
  if (is.null(required_packages) == FALSE){
    clusterCall(cl, function(pkgs) {
      for (req in pkgs) {
        library(req, character.only=TRUE)
      }
    }, required_packages)
  }

  isoform_cds <- buildLineageBranchCellDataSet(isoform_cds)
  gene_ids <- unique(fData(isoform_cds)$gene_id)
  res <- parLapply(cl, gene_ids, function(gene, isoform_cds){
  #print (head(iso_exprs))
  fd_isoforms <- subset(fData(isoform_cds), gene_id == gene)
  #print (head(iso_exprs))
  tryCatch({
    fd_isoforms <- subset(fd_isoforms, num_cells_expressed >= 15)
    if (nrow(fd_isoforms) > 1){
      iso_names <- row.names(fd_isoforms)
      #print (iso_names)
      iso_exprs <- exprs(isoform_cds[iso_names,])
      
      iso_exprs <- round(t(iso_exprs) / pData(isoform_cds)$Size_Factor)

      #iso_exprs <- iso_exprs[,colSums(iso_exprs > 1) >= 15]
      #print (head(iso_exprs))

      iso_names <- colnames(iso_exprs)
      iso_exprs <- cbind(iso_exprs, pData(isoform_cds))

      iso_exprs <- iso_exprs[rowSums(iso_exprs[,iso_names]) > 0,]

      # fit_time <- vglm(as.formula(paste("cbind(",toString(iso_names),") ~ sm.ns(Pseudotime, df=3)", sep="")), dirmultinomial(lphi="logit"), data = iso_exprs, epsilon=0.1)
      # fit_null<- vglm(as.formula(paste("cbind(",toString(iso_names),") ~ 1", sep="")), dirmultinomial(lphi="logit"), data = iso_exprs, epsilon=0.1)
      # link_phi <- predict(fit_time)
      # link_phi <- unique(link_phi[,ncol(link_phi)])[1]

      fit_time <- vglm(as.formula(paste("cbind(",toString(iso_names),") ~ sm.ns(Pseudotime, df=3) * Lineage", sep="")), dirmultinomial(), data = iso_exprs, epsilon=0.1)
      fit_null<- vglm(as.formula(paste("cbind(",toString(iso_names),") ~ sm.ns(Pseudotime, df=3)", sep="")), dirmultinomial(), data = iso_exprs, epsilon=0.1)


      #print(summary(fit_time))
      
      lrt <- lrtest(fit_time,fit_null) 
      #print(lrt)
      pval=lrt@Body["Pr(>Chisq)"][2,]


      #print(pval)
      test_res <- data.frame(link_phi=NA, status = "OK", pval=pval)
    }else{
      test_res <- data.frame(link_phi=NA, status = "NOTEST", pval=1.0)
    }
  }, error = function(e) { print (e); data.frame(link_phi=NA, status = "LOWDATA", pval=1.0)})
  #print(test_res)
  }, isoform_cds)
  res <- do.call(rbind, res)
  row.names(res) <- gene_ids
  res
}




