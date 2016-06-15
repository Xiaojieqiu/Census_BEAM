
bration_sim<- sim_escape(sandberg_rpc[autosomal_gene_records,], 
                                  formulaStr="~embryo + Spike_factor", 
                                  allele_thresh=0.1,
                                  dispersion=0,
                                  nsims=100,
                                  cores=detectCores() / 2)
rpc_dispersion <- median(disp_calibration_sim$link_phi, na.rm=T)

qplot(mean_expression, link_phi, data=disp_calibration_sim, log="xy") + geom_smooth()

#qplot(link_phi, data=disp_calibration_sim)

autosome_sim_res<- sim_monoallelic_expr(sandberg_rpc[autosomal_gene_records,], 
                                        formulaStr="~embryo + Spike_factor", 
                                        allele_thresh=0.1,
                                        dispersion=rpc_dispersion, #maybe from rpc_dispersion
                                        nsims=100,
                                        cores=detectCores() / 2)
autosome_sim_res$Region <- "Autosomal"

autosome_sim_res$expr_interval <- cut_interval(log(autosome_sim_res$mean_expression), 25)
expected_autosome_sim_res <- subset(autosome_sim_res, variable=="expected_monoallelic_calls" & is.na(value) == FALSE)
obs_autosome_sim_res <- subset(autosome_sim_res, variable=="obs_monoallelic_calls" & is.na(value) == FALSE)


exp_med_df <- expected_autosome_sim_res %>% 
  group_by(expr_interval) %>% 
  mutate(interval_expression=median(mean_expression),
         med_expr=median(value), 
         up = quantile(value, probs = 0.95, na.rm=TRUE),
         lo = quantile(value, probs = 0.05, na.rm=TRUE)) %>% 
  select(interval_expression, expr_interval, med_expr, up, lo) %>% distinct()


obs_med_df <- obs_autosome_sim_res %>% 
  group_by(expr_interval) %>% 
  mutate(interval_expression=median(mean_expression),
         med_expr=median(value), 
         up = quantile(value, probs = 0.95, na.rm=TRUE),
         lo = quantile(value, probs = 0.05, na.rm=TRUE)) %>% 
  select(interval_expression, expr_interval, med_expr, up, lo) %>% distinct()

#qplot(x=interval_expression, y=value, geom="boxplot", data=exp_med_df) + geom_line(aes(x=interval_expression, y=exp_med_df, group=1), color="red", data=obs_med_df) 


genes_interval_test_df <- inner_join(obs_autosome_sim_res, exp_med_df)
genes_above_interval <- subset(genes_interval_test_df, value > up)

pdf("fig6e.pdf", width=1.5, height=1)
ggplot(aes(x=interval_expression, y=med_expr, ymax=up, ymin=lo), data=exp_med_df) + 
  scale_x_continuous(trans="log2", limits=c(1,32)) + 
  scale_y_continuous(labels=scales::percent) + 
  geom_linerange() + 
  xlab("Transcript counts") + 
  ylab("Monoallelic genes\n(transcript counts)") + 
  geom_point(aes(x=interval_expression, y=value, group=1), data=genes_above_interval, color="red", position="jitter", size=0.01, alpha=0.05) + 
  geom_line(aes(x=interval_expression, y=med_expr, group=1), color="steelblue", size=I(0.25)) +
  geom_line(aes(x=interval_expression, y=med_expr, group=1), color="red", , size=I(0.25), data=obs_med_df) +
  monocle:::monocle_theme_opts()
dev.off()

###########
allelic_bias_test <- function(isoform_cds, 
              fullModelFormulaStr="~ stage",
              reducedModelFormulaStr="~1",
              cores=1){

  platform <- Sys.info()[['sysname']]
  if (platform == "Windows")
    cl <- makeCluster(cores)
  if (platform %in% c("Linux", "Darwin")) 
    cl <- makeCluster(cores, outfile="./debug.txt")
  
  cleanup <- function(){
    stopCluster(cl)
  }
  on.exit(cleanup)
  
  required_packages <- c("BiocGenerics", "Biobase", "VGAM", "plyr", "dplyr", "stats", "base")
  if (is.null(required_packages) == FALSE){
    clusterCall(cl, function(pkgs) {
      for (req in pkgs) {
        library(req, character.only=TRUE)
      }
    }, required_packages)
  }

  gene_ids <- unique(fData(isoform_cds)$gene_id)
  res <-parLapply(cl, gene_ids, function(gene_id_to_match, isoform_cds){
  #print (head(iso_exprs))
  fd_isoforms <- subset(fData(isoform_cds), gene_id == gene_id_to_match)
  print (head(fd_isoforms))
  # save(gene_ids, gene, fd_isoforms, file = 'test_allelic_bias_test')
  tryCatch({
    #fd_isoforms <- subset(fd_isoforms, num_cells_expressed >= 15)
    if (nrow(fd_isoforms) > 1){
      iso_names <- row.names(fd_isoforms)
      #print (iso_names)
      iso_exprs <- as.matrix(exprs(isoform_cds[iso_names,]))
      
      iso_exprs <- round(t(iso_exprs) / pData(isoform_cds)$Size_Factor)

      #iso_exprs <- iso_exprs[,colSums(iso_exprs > 1) >= 15]
      #print (head(iso_exprs))
      iso_names <- colnames(iso_exprs)
      iso_exprs <- cbind(iso_exprs, pData(isoform_cds))

      iso_exprs <- iso_exprs[rowSums(iso_exprs[,iso_names]) > 0,]

      fit_time <- glm(as.formula(paste("cbind(",toString(iso_names),")", fullModelFormulaStr, sep="")), quasibinomial(), data = iso_exprs, epsilon=0.1)
      fit_null<- glm(as.formula(paste("cbind(",toString(iso_names),")", reducedModelFormulaStr, sep="")), quasibinomial(), data = iso_exprs, epsilon=0.1)

      test_res <- anova(fit_time, fit_null, test="F")
      #print (str(test_res))
      pval=test_res$`Pr(>F)`[2]

      print(pval)
      test_res <- data.frame(link_phi=NA, status = "OK", pval=pval)
    }else{
      test_res <- data.frame(link_phi=NA, status = "NOTEST", pval=1.0)
    }
  }, error = function(e) { print (e); data.frame(link_phi=NA, status = "FAIL", pval=1.0)})
  #print(test_res)
  }, isoform_cds)
  res <- do.call(rbind, res)
  row.names(res) <- gene_ids
  res
}

calculate_allele_ratio_matrix <- function(isoform_cds, gene_names, num_cells_expressed = 10, cores=1){


  platform <- Sys.info()[['sysname']]
  if (platform == "Windows")
    cl <- makeCluster(cores)
  if (platform %in% c("Linux", "Darwin")) 
    cl <- makeCluster(cores, outfile="./debug.txt")
  
  cleanup <- function(){
    stopCluster(cl)
  }
  on.exit(cleanup)
  
  required_packages <- c("BiocGenerics", "Biobase", "VGAM", "plyr", "dplyr")
  if (is.null(required_packages) == FALSE){
    clusterCall(cl, function(pkgs) {
      for (req in pkgs) {
        library(req, character.only=TRUE)
      }
    }, required_packages)
  }
  inv.logit <- function(x) exp(x)/(1+exp(x))
  gene_names <- setdiff(gene_names, NA)
  gene_ids <- row.names(subset(fData(isoform_cds), as.character(symbol) %in% as.character(gene_names) & num_cells_expressed >= num_cells_expressed))
  #print (length(gene_ids))
  gene_cds <- isoform_cds[gene_ids,]
  #print (gene_cds)
  gene_ids <- unique(fData(gene_cds)$gene_id)

  res <- parLapply(cl, gene_ids, function(gene_id_to_match, isoform_cds){
  #print (head(iso_exprs))
  fd_isoforms <- subset(fData(isoform_cds), gene_id == gene_id_to_match)
#   print (head(iso_exprs))
  tryCatch({
    #fd_isoforms <- subset(fd_isoforms, num_cells_expressed >= 15)
    if (nrow(fd_isoforms) > 1){
      iso_names <- row.names(fd_isoforms)
      iso_names <- sort(iso_names)
      #print (iso_names)
      iso_exprs <- as.matrix(exprs(isoform_cds[iso_names,]))
      
      iso_exprs <- round(t(iso_exprs) / pData(isoform_cds)$Size_Factor)

      #iso_exprs <- iso_exprs[,colSums(iso_exprs > 1) >= 15]
      #print (head(iso_exprs))
      iso_names <- colnames(iso_exprs)
      iso_exprs <- cbind(iso_exprs, pData(isoform_cds))

      iso_exprs <- iso_exprs[rowSums(iso_exprs[,iso_names]) > 0,]

      fit_time <- glm(as.formula(paste("cbind(",toString(iso_names),") ~ stage", sep="")), quasibinomial(), data = iso_exprs, epsilon=0.1)
      
      #print(summary(fit_time))
      #print (head(coef(fit_time)))
      newdata=data.frame(stage = levels(pData(isoform_cds)$stage))
      id <- which(!(paste("stage",newdata$stage,sep="") %in% names(coef(fit_time))))
      newdata$stage[id] <- NA
      
      res_mat <- as.matrix(rep(NA, length(levels(pData(isoform_cds)$stage))))
      row.names(res_mat) <- levels(pData(isoform_cds)$stage)
      colnames(res_mat) <- rep(gene_id_to_match, ncol(res_mat))
      res_df <- distinct(data.frame(stage=iso_exprs$stage, percent_maternal = signif(predict(fit_time, type="response"))))
      res_mat[as.character(res_df$stage),] <- res_df$percent_maternal
      #print (head(res_mat))

      #print(head(res_mat))
      res_mat
    }else{
      res_mat <- as.matrix(rep(NA, length(levels(pData(isoform_cds)$stage))))
      colnames(res_mat) <- rep(gene_id_to_match, ncol(res_mat))
      res_mat
    }
  }, error = function(e) { 
    print (e); 
    res_mat <- as.matrix(rep(NA, length(levels(pData(isoform_cds)$stage))))
    colnames(res_mat) <- rep(gene_id_to_match, ncol(res_mat))
    res_mat
    })
  #print(test_res)
  #
  }, isoform_cds)
  print (head(res))

  res <- as.matrix(t(do.call(cbind, res)))
  res[is.na(res)] <- 0.5
  #res <- res[rowSums(is.na(res)) == 0,]
  res <- res[rowSums(is.nan(res)) == 0,]
  res <- res[rowSums(is.finite(res)) == ncol(res),]
  res <- res[rowSds(res) > 0,]

  res
}

plot_allele_ratio_heatmap <- function(isoform_cds, res, num_clusters=3, return_heatmap=TRUE){

  ph <- pheatmap(res, 
                 useRaster = T,
                 cluster_cols=FALSE, 
                 cluster_rows=TRUE, 
                 show_rownames=F, 
                 show_colnames=F, 
                 #clustering_distance_rows="correlation",
                 clustering_method = "ward.D2",
                 cutree_rows=num_clusters,
                 silent=TRUE,
                 filename=NA)

  annotation_row <- data.frame(#Cluster=factor(cutree(ph$tree_row, num_clusters)),
                               Region=factor(rep("Autosomal", length(row.names(res))), levels=c("Autosomal", "X")))


  annotation_row$Region[fData(isoform_cds)[paste(row.names(res), "_M", sep=""),]$chr =="X"] <- "X"
  annotation_colors <- list("Region"=c("Autosomal"="gray", "X"="red"))
  feature_label <- row.names(res)
  row_ann_labels <- row.names(res)
  
  
  row.names(res) <- feature_label
  row.names(annotation_row) <- row_ann_labels
  print (head(annotation_row))
  
  #colnames(res) <- c(1:ncol(res))
  

  
  #hmcols <- colorRampPalette(c("green4", "green", "white","violet","purple"))(100)
  hmcols <- colorRampPalette(rev(brewer.pal(n = 7, name ="PiYG")))(100)
  ph_res <- pheatmap(res, 
                     useRaster = T,
                     cluster_cols=FALSE, 
                     cluster_rows=TRUE, 
                     show_rownames=F, 
                     show_colnames=T, 
                     #clustering_distance_rows="correlation",
                     clustering_method = "ward.D2",
                     cutree_rows=num_clusters,
                     # cutree_cols = 2,
                     annotation_row=annotation_row,
                     annotation_colors=annotation_colors,
                     treeheight_row = 20, 
                     fontsize = 6,
                     silent=TRUE,
                     color=hmcols,
                     filename=NA,
                     lwd=0.5
  )
  
  grid::grid.rect(gp=grid::gpar("fill", col=NA))
  grid::grid.draw(ph_res$gtable)

  if (return_heatmap){
    return(ph_res)
  }
}


simulate_mRNA_counts <- function(isoform_cds, 
                                 fullModelFormulaStr,
                                 allele_thresh=0.1,
                                 dispersion=0.0,
                                 nsims=100,
                                 cores=1){

  platform <- Sys.info()[['sysname']]
  if (platform == "Windows")
    cl <- makeCluster(cores)
  if (platform %in% c("Linux", "Darwin")) 
    cl <- makeCluster(cores, outfile="debug_simulate_mRNA_counts.txt")
  
  cleanup <- function(){
    stopCluster(cl)
  }
  on.exit(cleanup)
  
  required_packages <- c("BiocGenerics", "Biobase", "VGAM", "plyr", "dplyr", "monocle")
  if (is.null(required_packages) == FALSE){
    clusterCall(cl, function(pkgs) {
      for (req in pkgs) {
        library(req, character.only=TRUE)
      }
    }, required_packages)
  }

  gene_ids <- unique(fData(isoform_cds)$gene)
  res <- parLapply(cl, gene_ids, function(gene_id_to_match, isoform_cds, dispersion){
    #print (head(iso_exprs))
    fd_isoforms <- subset(fData(isoform_cds), gene == gene_id_to_match)
    #print (head(iso_exprs))
    tryCatch({
      #fd_isoforms <- subset(fd_isoforms, num_cells_expressed >= 15)
      if (nrow(fd_isoforms) > 1){
        #print ("******")
        iso_names <- sort(row.names(fd_isoforms))
        #print (iso_names)
        #print (iso_names)
        iso_exprs <- as.matrix(exprs(isoform_cds[iso_names,]))
        #print (dim(iso_exprs))
        iso_exprs <- round(t(iso_exprs) / pData(isoform_cds)$Size_Factor)
        #print (dim(iso_exprs))
        #iso_exprs <- iso_exprs[,colSums(iso_exprs > 1) >= 15]
        #print (head(iso_exprs))

        iso_names <- colnames(iso_exprs)
        iso_exprs <- cbind(iso_exprs, pData(isoform_cds))

        #print (dim(iso_exprs))
        iso_exprs <- iso_exprs[rowSums(iso_exprs[,iso_names]) > 0,]

        #print (dim(iso_exprs))
        #print (head(iso_exprs))
        total_per_gene <- rowSums(as.matrix(iso_exprs[,c(1,2)]))
      
        full_formula_str <- paste("cbind(",toString(iso_names),")", fullModelFormulaStr, sep="")
        #print (full_formula_str)
        fit_time <- vglm(as.formula(full_formula_str), binomialff(dispersion=dispersion), data = iso_exprs, epsilon=0.1)
        #print (summary(fit_time))
        #print (total_per_gene)
        sims <- simulate(fit_time, nsim=nsims)
        #print (sims)
        #sims <- sims/total_per_gene
        # print (ncol(fitted(fit_time)))
        # print (dim(sims))
        # print (sims)
        # print (dim(iso_exprs))
        # print (colnames(iso_exprs))
        # print (row.names(iso_exprs))
        sims <- round(sims * total_per_gene)/total_per_gene #sims can be 0 or total_per_gene which gives 0/1

        expected_monoallelic_calls <- colSums(sims == 0) + colSums(sims == 1)

        obs_monoallelic_calls <- sum(iso_exprs[,1]/total_per_gene == 0) + sum(iso_exprs[,1]/total_per_gene == 1)
        #print (qplot(expected_monoallelic_calls))
        #print (expected_monoallelic_calls)
        #print (obs_monoallelic_calls)
        #print (sum(obs_monoallelic_calls > expected_monoallelic_calls))
        empirical_p_val <- signif(1 - (sum(obs_monoallelic_calls >= expected_monoallelic_calls) + 1)/(ncol(sims) + 1), 1/ncol(sims))
        #print(dim(iso_exprs))
        


        # sim_paternal_expr <- 1 - sims
        # obs_paternal_nonzero <- sum(iso_exprs[,2] > 0)
        # expected_paternal_nonzero <- colSums(sim_paternal_expr > 0)
        # paternal_expr_pval <- signif((sum(expected_paternal_nonzero == 0) + 1)/(ncol(sims) + 1), 1/ncol(sims))
        
        # sim_maternal_expr <- sims
        # obs_maternal_nonzero <- sum(iso_exprs[,1] > 0)
        # expected_maternal_nonzero <- colSums(sim_maternal_expr > 0)
        # maternal_expr_pval <- signif((sum(expected_maternal_nonzero == 0) + 1)/(ncol(sims) + 1), 1/ncol(sims))
        #print(sims[,1:5])
        sim_paternal_expr <- 1 - sims
        #print (colMeans(sim_paternal_expr)[1:5])
        obs_paternal_nonzero <- sum(iso_exprs[,2] > 0)
        expected_paternal_below_thresh <- sum(colMeans(sim_paternal_expr) <= allele_thresh)
        #print (expected_paternal_below_thresh)
        paternal_expr_pval <- signif((expected_paternal_below_thresh + 1)/(ncol(sims) + 1), 1/ncol(sims))
        
        sim_maternal_expr <- sims
        #print (colMeans(sim_maternal_expr)[1:5])
        obs_paternal_nonzero <- sum(iso_exprs[,1] > 0)
        expected_maternal_below_thresh <- sum(colMeans(sim_maternal_expr) <= allele_thresh)
        #print(expected_maternal_below_thresh)
        maternal_expr_pval <- signif((expected_maternal_below_thresh + 1)/(ncol(sims) + 1), 1/ncol(sims))


        #obs_allele_bias <- mean(sum(iso_exprs[,1]/total_per_gene)) 
        mean_bias <- colMeans(sims)
        num_samples_biased_toward_paternal <- sum(colMeans(sims) > 0.5)
        #print (num_samples_biased_toward_paternal)

        num_samples_biased_toward_maternal <- sum(colMeans(sims) < 0.5)
        #print (num_samples_biased_toward_maternal)

        paternal_bias_pval <- (num_samples_biased_toward_paternal + 1)/(ncol(sims) + 1)
        
        maternal_bias_pval <- (num_samples_biased_toward_paternal + 1)/(ncol(sims) + 1)
        
        bias_pval <- min(paternal_bias_pval, maternal_bias_pval)

        print(bias_pval)
        allele_bias_pval <- signif(bias_pval, 1/ncol(sims))

        #link_phi <- predict(fit_time)
        #link_phi <- unique(link_phi[,ncol(link_phi)])[1]
        link_phi <- fit_time@misc$dispersion
        test_res <- data.frame(link_phi=link_phi, 
                               status = "OK", 
                               mean_expression=mean(total_per_gene),
                               expected_monoallelic_calls=mean(expected_monoallelic_calls)/length(total_per_gene), 
                               obs_monoallelic_calls=obs_monoallelic_calls/length(total_per_gene), 
                               monoallele_pval=empirical_p_val,
                               maternal_expr_pval=maternal_expr_pval,
                               paternal_expr_pval=paternal_expr_pval,
                               allele_bias_pval=allele_bias_pval)

      }else{
        test_res <- data.frame(link_phi=NA, 
                               status = "NOTEST", 
                               mean_expression=NA,
                               expected_monoallelic_calls=NA, 
                               obs_monoallelic_calls=NA,
                               monoallele_pval=1.0,
                               maternal_expr_pval=1.0,
                               paternal_expr_pval=1.0,
                               allele_bias_pval=1.0)
      }
    }, error = function(e) { print (e); data.frame(link_phi=NA, 
                                                   status = "FAIL", 
                                                   mean_expression=NA,
                                                   expected_monoallelic_calls=NA, 
                                                   obs_monoallelic_calls=NA,
                                                   monoallele_pval=1.0,
                                                   maternal_expr_pval=1.0,
                                                   paternal_expr_pval=1.0, 
                                                   allele_bias_pval=1.0)})
  #print(test_res)
  }, isoform_cds, dispersion)
  res <- do.call(rbind, res)
  row.names(res) <- gene_ids
  res
}

sim_escape <- function(sandberg_records, formulaStr, allele_thresh, dispersion=0, expr_bins=10, nsims=100, cores=1){
  type_grouping <- pData(sandberg_records) %>% add_rownames() #%>% group_by(type)
  sim_res <- type_grouping %>% do(simulate_mRNA_counts(sandberg_records[,.$rowname], 
                                fullModelFormulaStr=formulaStr,
                                allele_thresh=allele_thresh,
                                dispersion=dispersion,
                                nsims=nsims,
                                cores=cores))
  print (dim(sim_res))

  sim_res <- as.data.frame(sim_res)
  sim_res$gene_id <- row.names(sim_res)
  sim_res$monoallele_qval <- 1.0
  sim_res$monoallele_qval[sim_res$status == "OK"] <- p.adjust(subset(sim_res, status == "OK")$monoallele_pval, method="fdr")

  sim_res$maternal_expr_qval <- 1.0
  sim_res$maternal_expr_qval[sim_res$status == "OK"] <- p.adjust(subset(sim_res, status == "OK")$maternal_expr_pval, method="fdr")

  sim_res$paternal_expr_qval <- 1.0
  sim_res$paternal_expr_qval[sim_res$status == "OK"] <- p.adjust(subset(sim_res, status == "OK")$paternal_expr_pval, method="fdr")

  sim_res$allelic_bias_qval <- 1.0
  sim_res$allelic_bias_qval[sim_res$status == "OK"] <- p.adjust(subset(sim_res, status == "OK")$allele_bias_pval, method="fdr")

  return(sim_res)
}

sim_monoallelic_expr <- function(sandberg_records, formulaStr, allele_thresh, dispersion=0, expr_bins=10, nsims=100, cores=1){
  type_grouping <- pData(sandberg_records) %>% add_rownames() #%>% group_by(type)
  sim_res <- type_grouping %>% do(simulate_mRNA_counts(sandberg_records[,.$rowname], 
                                fullModelFormulaStr=formulaStr,
                                allele_thresh=allele_thresh,
                                dispersion=dispersion,
                                nsims=nsims,
                                cores=cores))
  print (dim(sim_res))

  sim_res <- as.data.frame(sim_res)
  sim_res$gene_id <- row.names(sim_res)
  sim_res$monoallele_qval <- 1.0
  sim_res$monoallele_qval[sim_res$status == "OK"] <- p.adjust(subset(sim_res, status == "OK")$monoallele_pval, method="fdr")

  sim_res$maternal_expr_qval <- 1.0
  sim_res$maternal_expr_qval[sim_res$status == "OK"] <- p.adjust(subset(sim_res, status == "OK")$maternal_expr_pval, method="fdr")

  sim_res$paternal_expr_qval <- 1.0
  sim_res$paternal_expr_qval[sim_res$status == "OK"] <- p.adjust(subset(sim_res, status == "OK")$paternal_expr_pval, method="fdr")

  sim_res$allelic_bias_qval <- 1.0
  sim_res$allelic_bias_qval[sim_res$status == "OK"] <- p.adjust(subset(sim_res, status == "OK")$allele_bias_pval, method="fdr")


  sim_res_selected <- sim_res #subset(sim_res, exp(link_phi)/(1+exp(link_phi)) == 0)

  sim_res_expected_vs_obs <- rbind(melt(sim_res_selected[,c("gene_id", "expected_monoallelic_calls")]), 
                                 melt(sim_res_selected[,c("gene_id", "obs_monoallelic_calls")]))
  sim_res_expected_vs_obs <- merge(sim_res_expected_vs_obs, sim_res_selected[,c("gene_id", "mean_expression")])

  sim_res_expected_vs_obs$expr_quantile <- as.numeric(cut_interval(log(sim_res_expected_vs_obs$mean_expression), n = expr_bins))

  return(sim_res_expected_vs_obs)
}

# # obtain the ladder of the spike-in they added in: 
# spike_exprs <- exprs(with.spikes.sandberg.transcript.tpm[spike_ids, ])
# total_spike_exprs <- matrix(rep(colSums(spike_exprs), each = nrow(spike_exprs)), nrow = nrow(spike_exprs))
# spike_exprs_ladder <- spike_exprs / total_spike_exprs


