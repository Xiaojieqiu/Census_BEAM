#beam test on the isoform switch: 
iso_switch_test <- function(isoform_cds, cores=1){

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

plot_isoform_expression_heatmap <- function(isoform_cds, gene_names, num_clusters=6, add_annotation_row=NULL){
  gene_cds <- isoform_cds[row.names(subset(fData(isoform_cds), gene_short_name %in% gene_names & num_cells_expressed >= 15)),]
  isoform_exprs <- vstExprs(gene_cds)
  isoform_exprs <- t(scale(t(isoform_exprs)))
  isoform_exprs[isoform_exprs > 3] <- 3
  isoform_exprs[isoform_exprs < -3] <- -3
  isoform_exprs[is.na(isoform_exprs)] <- 0
  #pheatmap(isoform_exprs, cluster_cols=F, show_rownames=F, show_colnames=F)
  ph_res <- plot_pseudotime_heatmap(gene_cds, 
                                    use_gene_short_name=FALSE, 
                                    num_clusters=num_clusters, 
                                    add_annotation_row=add_annotation_row,
                                    return_heatmap=TRUE)
  return(ph_res)
}





