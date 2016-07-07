# # library(monocle)
# library(devtools)
# # load_all('~/Projects/monocle-dev')
# load_all('~/Dropbox (Personal)/Projects/monocle-dev/')
# library(stringr)
# library(dplyr)
# library(pheatmap)
# library(matrixStats)
# library(reshape2)
# library(zoo)
# 
# #########################
# theme_set(theme_bw())
# 
# rpc_dispersion <- 1.05777
# 
# # Set global ggplot2 properties for making print-scaled PDF panels
# theme_set(theme_bw(base_size=6))
# 
# time_colors = c("#D7191C", "#FDAE61", "#ABD9E9", "#2C7BB6")
# 
# 
# load('sandberg.gene.rpc.nbt-3rd-submission.RData')
# sandberg <- sandberg.gene.rpc
# 
# fData(sandberg)$gene_id <- str_split_fixed(row.names(fData(sandberg)), "_", 2)[,1]
# fData(sandberg)$allele_id <- row.names(fData(sandberg))
# 
# sandberg <- estimateSizeFactors(sandberg)
# sandberg <- detectGenes(sandberg, min_expr=1)
# sandberg <- estimateDispersions(sandberg)
# fData(sandberg)$gene_short_name <- fData(sandberg)$symbol
# pData(sandberg)$stage <- factor(pData(sandberg)$stage, levels = c("zygote",
#                                                                   "early-2-cell", 
#                                                                   "mid-2-cell",
#                                                                   "late-2-cell",
#                                                                   "4-cell",
#                                                                   "8-cell",
#                                                                   "16-cell",
#                                                                   "early-blast",
#                                                                   "mid-blast",
#                                                                   "late-blast"))
# 
# qplot(stage, fill=sex, geom="bar", data=pData(sandberg))
# 
# allelic_bias_test <- function(isoform_cds, 
#                               fullModelFormulaStr="~ stage",
#                               reducedModelFormulaStr="~1",
#                               cores=1){
#   
#   platform <- Sys.info()[['sysname']]
#   if (platform == "Windows")
#     cl <- makeCluster(cores)
#   if (platform %in% c("Linux", "Darwin")) 
#     cl <- makeCluster(cores)
#   
#   cleanup <- function(){
#     stopCluster(cl)
#   }
#   on.exit(cleanup)
#   
#   required_packages <- c("BiocGenerics", "Biobase", "VGAM", "plyr", "dplyr")
#   if (is.null(required_packages) == FALSE){
#     clusterCall(cl, function(pkgs) {
#       for (req in pkgs) {
#         library(req, character.only=TRUE)
#       }
#     }, required_packages)
#   }
#   
#   gene_ids <- unique(fData(isoform_cds)$gene_id)
#   res <-parLapply(cl, gene_ids, function(gene, isoform_cds){
#     #print (head(iso_exprs))
#     fd_isoforms <- subset(fData(isoform_cds), gene_id == gene)
#     #print (head(iso_exprs))
#     tryCatch({
#       #fd_isoforms <- subset(fd_isoforms, num_cells_expressed >= 15)
#       if (nrow(fd_isoforms) > 1){
#         iso_names <- row.names(fd_isoforms)
#         #print (iso_names)
#         iso_exprs <- as.matrix(exprs(isoform_cds[iso_names,]))
#         
#         iso_exprs <- round(t(iso_exprs) / pData(isoform_cds)$Size_Factor)
#         
#         #iso_exprs <- iso_exprs[,colSums(iso_exprs > 1) >= 15]
#         #print (head(iso_exprs))
#         iso_names <- colnames(iso_exprs)
#         iso_exprs <- cbind(iso_exprs, pData(isoform_cds))
#         
#         iso_exprs <- iso_exprs[rowSums(iso_exprs[,iso_names]) > 0,]
#         
#         fit_time <- glm(as.formula(paste("cbind(",toString(iso_names),")", fullModelFormulaStr, sep="")), quasibinomial(), data = iso_exprs, epsilon=0.1)
#         fit_null<- glm(as.formula(paste("cbind(",toString(iso_names),")", reducedModelFormulaStr, sep="")), quasibinomial(), data = iso_exprs, epsilon=0.1)
#         
#         #pred_out <- predict(fit_time, se.fit = TRUE)
#         #print (pred_out)
#         #param_df <- data.frame(fitted_vals = pred_out$fit, se_fit = pred_out$se.fit)
#         #param_df <- distinct(param_df)
#         #print(head(param_df))
#         #print(summary(fit_time))
#         #print(logit(0.5))
#         # pval_lower <- pnorm(0, mean = param_df$fitted_vals, sd = param_df$se_fit)
#         # pval_upper <- pnorm(0, mean = param_df$fitted_vals, sd = param_df$se_fit, lower.tail = FALSE)
#         # pvals <- rep(0, length(pval_lower))
#         # pvals <- p.adjust(pvals, method="fdr")
#         # pvals[param_df$fitted_vals > 0] <- pval_lower[param_df$fitted_vals > 0]
#         # pvals[param_df$fitted_vals <= 0] <- pval_upper[param_df$fitted_vals <= 0]
#         
#         #link_phi <- predict(fit_time)
#         #link_phi <- unique(link_phi[,ncol(link_phi)])[1]
#         
#         test_res <- anova(fit_time, fit_null, test="F")
#         #print (str(test_res))
#         pval=test_res$`Pr(>F)`[2]
#         
#         #print(pval)
#         test_res <- data.frame(link_phi=NA, status = "OK", pval=pval)
#       }else{
#         test_res <- data.frame(link_phi=NA, status = "NOTEST", pval=1.0)
#       }
#     }, error = function(e) { print (e); data.frame(link_phi=NA, status = "FAIL", pval=1.0)})
#     #print(test_res)
#   }, isoform_cds)
#   res <- do.call(rbind, res)
#   row.names(res) <- gene_ids
#   res
# }
# 
# 
# target_genes <-row.names(subset(fData(sandberg), symbol %in% c("Gnai3")))
# 
# target_genes <- row.names(fData(sandberg)[1:20,])
# 
# iso_switch_test_res <- allelic_bias_test(sandberg[target_genes,], fullModelFormulaStr="~stage", cores=detectCores() / 2)
# 
# time_de_res <- differentialGeneTest(sandberg[row.names(subset(fData(sandberg), num_cells_expressed >= 15)),], fullModelFormulaStr="~stage", cores=detectCores() / 2, verbose=T)
# sig_genes <- subset(time_de_res, qval < 0.1)$gene_id
# 
# iso_switch_test_res <- allelic_bias_test(sandberg[row.names(subset(fData(sandberg), gene_id %in% sig_genes)),], fullModelFormulaStr="~stage", cores=detectCores() / 2)
# iso_switch_test_res <- merge(distinct(fData(sandberg)[,c("gene_id", "symbol", "chr")]), iso_switch_test_res, by.x="gene_id", by.y="row.names")
# row.names(iso_switch_test_res) <- iso_switch_test_res$gene_id
# iso_switch_test_res$qval <- 1.0
# iso_switch_test_res$qval[iso_switch_test_res$status == "OK"] <- p.adjust(subset(iso_switch_test_res, status == "OK")$pval, method="fdr")
# 
# allele_switch_genes <- subset(iso_switch_test_res, qval < 0.01 & chr == "X")
# allele_switch_records <- row.names(subset(fData(sandberg), symbol %in% allele_switch_genes$symbol))
# plot_genes_jitter(sandberg[allele_switch_records[1:16],], 
#                   grouping="stage", color_by="allele", plot_trend=TRUE) + 
#   facet_wrap( ~ feature_label, scales="free_y", ncol=1)
# 
# allele_switch_records <- row.names(subset(fData(sandberg), symbol %in% c("Nanog")))
# 
# plot_genes_jitter(sandberg[allele_switch_records,], 
#                   grouping="stage", color_by="allele", plot_trend=TRUE, min_expr=0.1) + 
#   facet_wrap( ~ feature_label, scales="free_y", ncol=1)
# 
# 
# 
# allele_switch_records <- row.names(subset(fData(sandberg), symbol %in% c("Xist", "Tsix", "Jpx", "Firre")))
# #iso_switch_test_res <- allelic_bias_test(sandberg[allele_switch_records,], fullModelFormulaStr="~stage", cores=2)
# plot_genes_jitter(sandberg[allele_switch_records,is.na(pData(sandberg)$sex) == FALSE & pData(sandberg)$sex == "F"], 
#                   grouping="stage", color_by="allele", plot_trend=TRUE, min_expr=0.1) + 
#   facet_wrap( ~ feature_label, scales="free_y", ncol=1)
# 
# 
# allele_switch_records <- row.names(subset(fData(sandberg), symbol %in% c("Gm44279", "Mrpl33", "Cbx7", "Dock5")))
# plot_genes_jitter(sandberg[allele_switch_records,], 
#                   grouping="stage", color_by="allele", plot_trend=TRUE, min_expr=0.1) + 
#   facet_wrap( ~ feature_label, scales="free_y", ncol=1)
# 
# # 
# calculate_allele_ratio_matrix <- function(isoform_cds, gene_names, cores=1){
#   
#   
#   platform <- Sys.info()[['sysname']]
#   if (platform == "Windows")
#     cl <- makeCluster(cores)
#   if (platform %in% c("Linux", "Darwin")) 
#     cl <- makeCluster(cores)
#   
#   cleanup <- function(){
#     stopCluster(cl)
#   }
#   on.exit(cleanup)
#   
#   required_packages <- c("BiocGenerics", "Biobase", "VGAM", "plyr", "dplyr")
#   if (is.null(required_packages) == FALSE){
#     clusterCall(cl, function(pkgs) {
#       for (req in pkgs) {
#         library(req, character.only=TRUE)
#       }
#     }, required_packages)
#   }
#   inv.logit <- function(x) exp(x)/(1+exp(x))
#   gene_names <- setdiff(gene_names, NA)
#   allele_ids <- row.names(subset(fData(isoform_cds), as.character(symbol) %in% as.character(gene_names) & num_cells_expressed >= 15))
#   #print (length(allele_ids))
#   gene_cds <- isoform_cds[allele_ids,]
#   #print (gene_cds)
#   gene_ids <- unique(fData(gene_cds)$gene_id)
#   
#   res <- parLapply(cl, gene_ids, function(gene, isoform_cds){
#     #print (head(iso_exprs))
#     fd_isoforms <- subset(fData(isoform_cds), gene_id == gene)
#     #print (head(iso_exprs))
#     tryCatch({
#       #fd_isoforms <- subset(fd_isoforms, num_cells_expressed >= 15)
#       if (nrow(fd_isoforms) > 1){
#         iso_names <- row.names(fd_isoforms)
#         iso_names <- sort(iso_names)
#         #print (iso_names)
#         iso_exprs <- as.matrix(exprs(isoform_cds[iso_names,]))
#         
#         iso_exprs <- round(t(iso_exprs) / pData(isoform_cds)$Size_Factor)
#         
#         #iso_exprs <- iso_exprs[,colSums(iso_exprs > 1) >= 15]
#         #print (head(iso_exprs))
#         iso_names <- colnames(iso_exprs)
#         iso_exprs <- cbind(iso_exprs, pData(isoform_cds))
#         
#         iso_exprs <- iso_exprs[rowSums(iso_exprs[,iso_names]) > 0,]
#         
#         fit_time <- glm(as.formula(paste("cbind(",toString(iso_names),") ~ stage", sep="")), quasibinomial(), data = iso_exprs, epsilon=0.1)
#         
#         #print(summary(fit_time))
#         #print (head(coef(fit_time)))
#         newdata=data.frame(stage = levels(pData(isoform_cds)$stage))
#         id <- which(!(paste("stage",newdata$stage,sep="") %in% names(coef(fit_time))))
#         newdata$stage[id] <- NA
#         
#         res_mat <- as.matrix(rep(NA, length(levels(pData(isoform_cds)$stage))))
#         row.names(res_mat) <- levels(pData(isoform_cds)$stage)
#         colnames(res_mat) <- rep(gene, ncol(res_mat))
#         res_df <- distinct(data.frame(stage=iso_exprs$stage, percent_maternal = signif(predict(fit_time, type="response"))))
#         res_mat[as.character(res_df$stage),] <- res_df$percent_maternal
#         #print (head(res_mat))
#         
#         #print(head(res_mat))
#         res_mat
#       }else{
#         res_mat <- as.matrix(rep(NA, length(levels(pData(isoform_cds)$stage))))
#         colnames(res_mat) <- rep(gene, ncol(res_mat))
#         res_mat
#       }
#     }, error = function(e) { 
#       print (e); 
#       res_mat <- as.matrix(rep(NA, length(levels(pData(isoform_cds)$stage))))
#       colnames(res_mat) <- rep(gene, ncol(res_mat))
#       res_mat
#     })
#     #print(test_res)
#     #
#   }, isoform_cds)
#   #print (head(res))
#   
#   res <- as.matrix(t(do.call(cbind, res)))
#   res[is.na(res)] <- 0.5
#   #res <- res[rowSums(is.na(res)) == 0,]
#   res <- res[rowSums(is.nan(res)) == 0,]
#   res <- res[rowSums(is.finite(res)) == ncol(res),]
#   res <- res[rowSds(res) > 0,]
#   res
# }
# 
# # X_inact_allele_matrix <- calculate_allele_ratio_matrix(sandberg_female_X, 
# #                                                        subset(fData(sandberg_female_X), num_cells_expressed >= 1)$symbol,
# #                                                        cores=detectCores() / 2)
# 
# 
# # 
# # allele_matrix <- calculate_allele_ratio_matrix(sandberg, 
# #                                                c("Xist", "Tsix", "Jpx", "Firre"),
# #                                                cores=2)
# 
# 
# 
# 
# 
# plot_allele_ratio_heatmap <- function(isoform_cds, res, num_clusters=3, return_heatmap=TRUE){
#   
#   ph <- pheatmap(res, 
#                  useRaster = T,
#                  cluster_cols=FALSE, 
#                  cluster_rows=TRUE, 
#                  show_rownames=F, 
#                  show_colnames=F, 
#                  #clustering_distance_rows="correlation",
#                  clustering_method = "ward.D2",
#                  cutree_rows=num_clusters,
#                  silent=TRUE,
#                  filename=NA)
#   
#   annotation_row <- data.frame(#Cluster=factor(cutree(ph$tree_row, num_clusters)),
#     Region=factor(rep("Autosomal", length(row.names(res))), levels=c("Autosomal", "X")))
#   
#   
#   annotation_row$Region[fData(isoform_cds)[paste(row.names(res), "_M", sep=""),]$chr =="X"] <- "X"
#   annotation_colors <- list("Region"=c("Autosomal"="gray", "X"="red"))
#   feature_label <- row.names(res)
#   row_ann_labels <- row.names(res)
#   
#   
#   row.names(res) <- feature_label
#   row.names(annotation_row) <- row_ann_labels
#   print (head(annotation_row))
#   
#   #colnames(res) <- c(1:ncol(res))
#   
#   
#   
#   #hmcols <- colorRampPalette(c("green4", "green", "white","violet","purple"))(100)
#   hmcols <- colorRampPalette(rev(brewer.pal(n = 7, name ="PiYG")))(100)
#   ph_res <- pheatmap(res, 
#                      useRaster = T,
#                      cluster_cols=FALSE, 
#                      cluster_rows=TRUE, 
#                      show_rownames=F, 
#                      show_colnames=T, 
#                      #clustering_distance_rows="correlation",
#                      clustering_method = "ward.D2",
#                      cutree_rows=num_clusters,
#                      # cutree_cols = 2,
#                      annotation_row=annotation_row,
#                      annotation_colors=annotation_colors,
#                      treeheight_row = 20, 
#                      fontsize = 6,
#                      silent=TRUE,
#                      color=hmcols,
#                      filename=NA,
#                      lwd=0.5
#   )
#   
#   grid::grid.rect(gp=grid::gpar("fill", col=NA))
#   grid::grid.draw(ph_res$gtable)
#   
#   if (return_heatmap){
#     return(ph_res)
#   }
# }
# 
# 
# allele_switch_genes <- subset(iso_switch_test_res, qval < 1e-4)
# allele_matrix <- calculate_allele_ratio_matrix(sandberg, 
#                                                allele_switch_genes$symbol,
#                                                cores=detectCores() / 2)
# 
# num_allele_switch_clusters <- 8
# pdf ("fig6a.pdf", height=4, width=3)
# #iso_switch_heatmap <- 
# plot_allele_ratio_heatmap(sandberg, 
#                           allele_matrix, 
#                           num_clusters=num_allele_switch_clusters)
# dev.off()
# 
# 
# 
# Deng_sexed_cells <- sandberg[,pData(sandberg)$stage %in% c("4-cell", "16-cell", "early-blast")]
# Deng_sexed_cells <- detectGenes(Deng_sexed_cells, min_expr=1)
# 
# Deng_Y_markers <- row.names(subset(fData(sandberg), symbol %in% c("Eif2s3y", "Zfy2", "Uba1y", "Kdm5d")))
# all_Y_genes <- row.names(subset(fData(sandberg), chr == "Y" & num_cells_expressed >= 5))

plot_genes_jitter(female_cds[row.names(subset(fData(Deng_sexed_cells), gene_short_name %in%  c("Eif2s3y", "Zfy2", "Uba1y", "Kdm5d"))),], 
                  grouping="stage", color_by="allele", plot_trend=TRUE, min_expr=0.5) 

Y_marker_exprs <- limma:::removeBatchEffect(vstExprs(Deng_sexed_cells[all_Y_genes,], round_vals=F), batch=pData(Deng_sexed_cells)$stage)

sex_cluster <- hclust(dist(t(Y_marker_exprs)), method="ward.D2")
plot(sex_cluster)
sex_cluster_cut <- cutree(sex_cluster, k=2)
#Y_marker_exprs <- colMeans(round(Y_marker_exprs))
sex_cluster_cut[sex_cluster_cut == 1] <- "Male"
sex_cluster_cut[sex_cluster_cut == 2] <- "Female"

pData(Deng_sexed_cells)$Sex_cluster <- sex_cluster_cut[row.names(pData(Deng_sexed_cells))]
sex_embryo_table <- table(pData(Deng_sexed_cells)$Sex_cluster, pData(Deng_sexed_cells)$embryo)
sex_embryo_table <- as.matrix(sex_embryo_table)
embryo_sexes <- sex_embryo_table["Female",] >  sex_embryo_table["Male",]
male_embryos <- names(embryo_sexes[embryo_sexes == FALSE])

Deng_sexed_cells$Deng_Sex <- factor(rep("Female", ncol(Deng_sexed_cells)), levels=c("Male", "Female"))
Deng_sexed_cells$Deng_Sex[pData(Deng_sexed_cells)$embryo %in% male_embryos] <- "Male"
#pData(sandberg)$Deng_Sex[Y_marker_exprs > 0] <- "Male"

gene_metadata <- read.delim("genes.attr_table")
gene_metadata <- subset(gene_metadata, grepl("ERCC", tracking_id) == FALSE)

locus_cols <- str_split_fixed(gene_metadata$locus, ":|-", 3) 
gene_metadata$gene_id <- str_split_fixed(gene_metadata$gene_id, "\\.", 2)[,1]

locus_df <- data.frame(gene_id=gene_metadata$gene_id, start=as.integer(locus_cols[,2]), stop=as.integer(locus_cols[,3]))


allele_switch_records <- row.names(subset(fData(Deng_sexed_cells), chr == "X" & is.allele.informative == 1))


X_inact_allele_matrix <- calculate_allele_ratio_matrix(Deng_sexed_cells[allele_switch_records, pData(Deng_sexed_cells)$Deng_Sex == "Female"], 
                                                       subset(fData(Deng_sexed_cells), num_cells_expressed >= 1)$symbol,
                                                       cores=detectCores() / 2)

X_inact_allele_matrix_melted <- melt(X_inact_allele_matrix)
colnames(X_inact_allele_matrix_melted) <- c("gene_id", "stage", "expression")
X_inact_allele_matrix_melted <- merge (X_inact_allele_matrix_melted, locus_df, by="gene_id")
X_inact_allele_matrix_melted <- subset(X_inact_allele_matrix_melted, stage %in% c("4-cell", "16-cell", "early-blast"))

X_inact_df <- X_inact_allele_matrix_melted %>%
  group_by(stage) %>%
  arrange(stage, start) %>%
  mutate(moving_avg_expr = rollmean(x = expression, 10, align = "right", fill = NA))

pdf("fig6b.pdf", width=1.35, height=1.5)
qplot(as.numeric(start) / 1e6, moving_avg_expr, data=X_inact_df, color=stage, geom="step", size=I(0.25)) + 
  xlab("X chr position (Mb)") + 
  ylab("Maternal expression") + 
  theme(legend.position = "none") + 
  geom_hline(yintercept=0.5, linetype=1) + 
  geom_vline(xintercept=103.414466, linetype=2) + 
  scale_color_brewer(palette="Set1") + 
  scale_y_continuous(limits=c(0,1), label=scales::percent) + monocle:::monocle_theme_opts()
dev.off()

plot_genes_jitter(Deng_sexed_cells[c("ENSMUSG00000000787_M", "ENSMUSG00000000787_P"),], 
                  grouping="stage", color_by="Deng_Sex", plot_trend=TRUE, min_expr=0.1) + 
  facet_wrap( ~ feature_label, scales="free_y", ncol=1)

############


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
    cl <- makeCluster(cores, outfile="xxx.txt")
  
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
  
  gene_ids <- unique(fData(isoform_cds)$gene_id)
  res <- parLapply(cl, gene_ids, function(gene, isoform_cds, dispersion){
    #print (head(iso_exprs))
    fd_isoforms <- subset(fData(isoform_cds), gene_id == gene)
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
        sims <- round(sims * total_per_gene)/total_per_gene
        
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
  }, isoform_cds, dispersion = dispersion)
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


X_inact_targets <- Deng_sexed_cells[row.names(subset(fData(Deng_sexed_cells), gene_short_name %in%   c("Xist", "Rlim", "Kif4", "Huwe1", "Atrx"))),
                                    pData(Deng_sexed_cells)$stage == "early-blast" &
                                      pData(Deng_sexed_cells)$Deng_Sex == "Female" ]
sim_escape(X_inact_targets, formulaStr="~1", allele_thresh=0.1, dispersion=rpc_dispersion, nsims=1000, cores=detectCores() / 2)




all_informative_sim_res<- sim_escape(Deng_sexed_cells[c("ENSMUSG00000086503_M", "ENSMUSG00000086503_P"),pData(Deng_sexed_cells)$Deng_Sex == "Female"], 
                                     formulaStr="~num_genes_expressed", 
                                     allele_thresh=0.1,
                                     dispersion=0.1,
                                     nsims=100,
                                     cores=1)




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


spike_ids <-  c("RNA_SPIKE_1_EC2_Gene_Bank_AE000151",                             
                "RNA_SPIKE_2_EC12_Gene_Bank_D86239",                               
                "RNA_SPIKE_3_EC3_Gene_Bank_AE000151",                                
                "RNA_SPIKE_4_EC15_ACCESSION_AE000185_U00096",                        
                "RNA_SPIKE_5_EC17_Accession_Number_V00307_J01654",                   
                "RNA_SPIKE_6_EC13_Gene_Bank_M34825",                                 
                "RNA_SPIKE_7_EC18_Accession_Number_M15263",                           
                "RNA_SPIKE_8_EC5_Accession_Number_AE000184_U00096")

load('with.spikes.sandberg.gene.rpc.nbt-3rd-submission.RData')
sandberg_rpc <- with.spikes.sandberg.gene.rpc

fData(sandberg_rpc)$gene_id <- str_split_fixed(row.names(fData(sandberg_rpc)), "_", 2)[,1]
fData(sandberg_rpc)$allele_id <- row.names(fData(sandberg_rpc))

sandberg_rpc <- detectGenes(sandberg_rpc)
sandberg_rpc <- estimateSizeFactors(sandberg_rpc)
sandberg_rpc <- estimateDispersions(sandberg_rpc)
spike_rpcs <- sandberg_rpc[spike_ids,]
spike_rpcs <- estimateSizeFactors(spike_rpcs)

#qplot(rowMeans(exprs(spike_rpcs))

pData(sandberg_rpc)$Spike_factor <- pData(spike_rpcs)$Size_Factor

autosomal_informative_genes <- as.character(subset(fData(sandberg_rpc), is.na(chr) == FALSE & chr %in% c("X", "Y") == FALSE & is.allele.informative == 1)$gene_id)


gencode_gene_metadata <- read.delim("gencode_biotypes.txt")

#locus_cols <- str_split_fixed(gencode_gene_metadata$locus, ":|-", 3) 
gencode_gene_metadata$gene_id <- str_split_fixed(gencode_gene_metadata$gene_id, "\\.", 2)[,1]


allowed_gene_types <- c("protein_coding", "lincRNA")
ribosomal_gene_ids <- subset(fData(sandberg_rpc), grepl("^Rps|^Rpl", symbol))$gene_id

autosomal_gene_records <- row.names(subset(fData(sandberg_rpc), gene_id %in% autosomal_informative_genes 
                                           & gene_id %in% subset(gencode_gene_metadata, biotype %in% allowed_gene_types)$gene_id &
                                             gene_id %in% ribosomal_gene_ids == FALSE))
#autosomal_gene_records <- setdiff(autosomal_gene_records, ribosomal_gene_ids)

disp_calibration_sim<- sim_escape(sandberg_rpc[autosomal_gene_records,], 
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
                                        dispersion=rpc_dispersion,
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

fig6e_a <- nrow(subset(genes_interval_test_df, value > up))
fig6e_b <- nrow(genes_interval_test_df)
fig6e_c <- 1 - nrow(subset(genes_interval_test_df, value > up)) / nrow(genes_interval_test_df)

print(nrow(subset(genes_interval_test_df, value > up)))
print(nrow(genes_interval_test_df))
print(1 - nrow(subset(genes_interval_test_df, value > up)) / nrow(genes_interval_test_df))

###########

load('with-spikes.sandberg.gene.read.counts.RData')

sandberg_count <- with.spikes.sandberg.gene.read.counts

fData(sandberg_count)$gene_id <- str_split_fixed(row.names(fData(sandberg_count)), "_", 2)[,1]
fData(sandberg_count)$allele_id <- row.names(fData(sandberg_count))

sandberg_count <- detectGenes(sandberg_count)
sandberg_count <- estimateSizeFactors(sandberg_count)
sandberg_count <- estimateDispersions(sandberg_count)
spike_counts <- sandberg_count[spike_ids,]
spike_counts <- estimateSizeFactors(spike_counts)

pData(sandberg_count)$Spike_factor <- pData(spike_counts)$Size_Factor

disp_calibration_sim<- sim_escape(sandberg_count[autosomal_gene_records,], 
                                  formulaStr="~embryo + Spike_factor", 
                                  allele_thresh=0.1,
                                  dispersion=0,
                                  nsims=100,
                                  cores=detectCores() / 2)
qplot(mean_expression, link_phi, data=subset(disp_calibration_sim, link_phi < 1e4), log="xy") + geom_smooth()

#qplot(link_phi, data=disp_calibration_sim)


autosome_sim_res_counts<- sim_monoallelic_expr(sandberg_count[autosomal_gene_records,], 
                                               formulaStr="~embryo + Spike_factor", 
                                               allele_thresh=0.1,
                                               dispersion=0,
                                               cores=detectCores() / 2)
autosome_sim_res_counts$Region <- "Autosomal"

rpc_expr_df <- autosome_sim_res %>% select(gene_id, mean_expression) %>% rename(mean_rpc = mean_expression)
autosome_sim_res_counts <- inner_join(autosome_sim_res_counts, rpc_expr_df)

autosome_sim_res_counts$expr_interval <- cut_interval(log(autosome_sim_res_counts$mean_rpc), 25)
expected_autosome_sim_res_counts <- subset(autosome_sim_res_counts, variable=="expected_monoallelic_calls" & is.na(value) == FALSE)
obs_autosome_sim_res_counts <- subset(autosome_sim_res_counts, variable=="obs_monoallelic_calls" & is.na(value) == FALSE)



exp_med_df <- expected_autosome_sim_res_counts %>% 
  group_by(expr_interval) %>% 
  mutate(interval_expression=median(mean_rpc),
         med_expr=median(value), 
         up = quantile(value, probs = 0.97, na.rm=TRUE),
         lo = quantile(value, probs = 0.05, na.rm=TRUE)) %>% 
  select(interval_expression, expr_interval, med_expr, up, lo) %>% distinct()


obs_med_df <- obs_autosome_sim_res_counts %>% 
  group_by(expr_interval) %>% 
  mutate(interval_expression=median(mean_rpc),
         med_expr=median(value), 
         up = quantile(value, probs = 0.95, na.rm=TRUE),
         lo = quantile(value, probs = 0.05, na.rm=TRUE)) %>% 
  select(interval_expression, expr_interval, med_expr, up, lo) %>% distinct()

#qplot(x=interval_expression, y=value, geom="boxplot", data=exp_med_df) + geom_line(aes(x=interval_expression, y=exp_med_df, group=1), color="red", data=obs_med_df) 
genes_interval_test_df <- inner_join(obs_autosome_sim_res_counts, exp_med_df)
genes_above_interval <- subset(genes_interval_test_df, value > up)

pdf("fig6f.pdf",  width=1.5, height=1)
ggplot(aes(x=interval_expression, y=med_expr, ymax=up, ymin=lo), data=exp_med_df) + 
  scale_x_continuous(trans="log2", limits=c(1,32)) + 
  scale_y_continuous(labels=scales::percent) + 
  geom_linerange() + 
  xlab("Transcript counts") + 
  ylab("Monoallelic genes\n(read counts)") + 
  geom_point(aes(x=interval_expression, y=value, group=1), data=genes_above_interval, color="red", position="jitter", size=0.01, alpha=0.05) + 
  geom_line(aes(x=interval_expression, y=med_expr, group=1), color="steelblue") +
  geom_line(aes(x=interval_expression, y=med_expr, group=1), color="red", data=obs_med_df) +
  monocle:::monocle_theme_opts()
dev.off()

fig6f_a <- nrow(subset(genes_interval_test_df, value > up))
fig6f_b <- nrow(genes_interval_test_df)
fig6f_c <- 1 - nrow(subset(genes_interval_test_df, value > up)) / nrow(genes_interval_test_df)

print(nrow(subset(genes_interval_test_df, value > up)))
print(nrow(genes_interval_test_df))
print(1 - nrow(subset(genes_interval_test_df, value > up)) / nrow(genes_interval_test_df))

################################



x_escape_candidate_gene_ids <- unique(as.character(subset(fData(Deng_sexed_cells), chr == "X" & is.allele.informative == 1)$gene_id))

#x_escape_candidate_gene_ids <- intersect(x_escape_candidate_gene_ids, subset(gencode_gene_metadata, biotype %in% c("protein_coding", "lincRNA", "antisense"))$gene_id)
x_escape_candidate_paternal_allele_ids <- row.names(subset(fData(Deng_sexed_cells), gene_id %in% x_escape_candidate_gene_ids & allele == "P"))

Deng_sexed_cells_paternal_alleles <- Deng_sexed_cells[x_escape_candidate_paternal_allele_ids,]
Deng_sexed_cells_paternal_alleles <- detectGenes(Deng_sexed_cells_paternal_alleles, min_exp=0.1)
x_escape_candidate_gene_ids <- unique(as.character(subset(fData(Deng_sexed_cells_paternal_alleles), num_cells_expressed >= 1)$gene_id))


#1.05777

sandberg_escape_cds <- Deng_sexed_cells[row.names(subset(fData(Deng_sexed_cells), gene_id %in% x_escape_candidate_gene_ids)),
                                        pData(Deng_sexed_cells)$stage == "early-blast" &
                                          pData(Deng_sexed_cells)$Deng_Sex == "Female" ]
sandberg_escape_sim_res_early_blast<- sim_escape(sandberg_escape_cds, formulaStr="~1", allele_thresh=0.1, dispersion=1.05777, nsims=1000, cores=detectCores() / 2)
sandberg_escape_sim_res_early_blast$stage= "early-blast"

#subset(fData(sandberg_escape_cds), gene_id %in% subset(sandberg_escape_sim_res, paternal_expr_qval < 0.05)$gene_id)

sandberg_escape_cds <- Deng_sexed_cells[row.names(subset(fData(Deng_sexed_cells), gene_id %in% x_escape_candidate_gene_ids)),
                                        pData(Deng_sexed_cells)$stage == "16-cell" &
                                          pData(Deng_sexed_cells)$Deng_Sex == "Female"]
sandberg_escape_sim_res_16_cell<- sim_escape(sandberg_escape_cds, formulaStr="~1", allele_thresh=0.1, dispersion=1.05777, nsims=1000, cores=detectCores() / 2)
sandberg_escape_sim_res_16_cell$stage= "16-cell"

sandberg_escape_cds <- Deng_sexed_cells[row.names(subset(fData(Deng_sexed_cells), gene_id %in% x_escape_candidate_gene_ids)),
                                        pData(Deng_sexed_cells)$stage == "4-cell" &
                                          pData(Deng_sexed_cells)$Deng_Sex == "Female"]
sandberg_escape_sim_res_4_cell<- sim_escape(sandberg_escape_cds, formulaStr="~1", allele_thresh=0.1, dispersion=1.05777, nsims=1000, cores=detectCores() / 2)
sandberg_escape_sim_res_4_cell$stage= "4-cell"

escape_sim_res <- rbind(sandberg_escape_sim_res_4_cell, sandberg_escape_sim_res_16_cell, sandberg_escape_sim_res_early_blast)

escape_sim_res$stage <- factor(escape_sim_res$stage, levels=c("4-cell", "16-cell", "early-blast"))
escape_sim_res$gene_id <- row.names(escape_sim_res)

allelic_bias_genes <- subset(escape_sim_res, allelic_bias_qval < 0.05)$gene_id
subset(fData(Deng_sexed_cells), gene_id %in% allelic_bias_genes)

escape_sim_res <- escape_sim_res[,c("gene_id", "stage", "maternal_expr_qval", "paternal_expr_qval")]
#Xist results:
subset(escape_sim_res, gene_id == "ENSMUSG00000086503")

#
#escape_sim_res <- subset(escape_sim_res, monoallele_qval )

escape_sim_res <- melt(escape_sim_res)

escape_df <- escape_sim_res %>% group_by(stage) %>%
  mutate(n_detected = sum(value < 0.01)) %>% 
  group_by(stage, variable) %>% 
  mutate(percent_allele = sum(value < 0.01) / n_detected) %>%
  select(stage, variable, percent_allele) %>%
  distinct()


#escape_sim_res <- merge(escape_sim_res, fData(Deng_sexed_cells)[,c("gene_id", "gene_short_name", "num_cells_expressed")], by.x="row.names", by.y="gene_id")
pdf("fig6c.pdf", width=1.35, height=1)
ggplot(aes(stage, percent_allele, fill=variable), data=escape_df) + 
  geom_bar(position="dodge", stat="identity") + scale_fill_manual(values=c("#C51B7D", "#4D9221")) + 
  ylab("Detectable alleles")+
  scale_y_continuous(label=scales::percent) +
  theme(axis.title.x=element_blank()) +
  theme(legend.position="none") + 
  #scale_fill_brewer(palette="Set1") + 
  monocle:::monocle_theme_opts()
dev.off()

early_blast_escape_genes <- subset(sandberg_escape_sim_res_early_blast, paternal_expr_qval < 0.05)$gene_id
subset(fData(Deng_sexed_cells), gene_id %in% early_blast_escape_genes)

known_escape_gene_names <- c("Xist", "Mid1", "Firre", "Kdm6a", "Pbdc1", "Eif2s3x", "Ddx3x", "Shroom4", "BC022960", "Kdm5c", "Car5b", "Bgn", "1810030O07Rik")
known_escape_genes <- unique(subset(fData(Deng_sexed_cells), symbol %in% c(known_escape_gene_names))$gene_id)
missing_escapees <- setdiff(known_escape_genes, early_blast_escape_genes)
missing_escapees <- subset(fData(Deng_sexed_cells), gene_id %in% missing_escapees)

early_blast_maternal <- subset(sandberg_escape_sim_res_early_blast, maternal_expr_qval < 0.01)$gene_id
subset(fData(Deng_sexed_cells), gene_id %in% early_blast_maternal)

female_cds <- Deng_sexed_cells[, pData(Deng_sexed_cells)$Deng_Sex == "Female" ]

plot_genes_jitter(female_cds[row.names(subset(fData(Deng_sexed_cells), gene_short_name %in%  c("Xist", "Jpx", "Eif2s3x",  "Ddx3x", "Kdm6a", "Pbdc1",  "Kdm5c", "Car5b"))),], 
                  grouping="stage", color_by="allele", plot_trend=TRUE, min_expr=0.5) 

plot_genes_jitter(female_cds[row.names(subset(fData(Deng_sexed_cells), gene_short_name %in%  c("Pgk1", "Gprasp1", "Klf8", "Mid1", "Fam3a", "Cetn2", "Dlg3"))),], 
                  grouping="stage", color_by="allele", plot_trend=TRUE, min_expr=0.5) 

early_blast_inactivated_genes <- subset(sandberg_escape_sim_res_early_blast, paternal_expr_qval > 0.05 & maternal_expr_qval < 0.05)$gene_id
subset(fData(Deng_sexed_cells), gene_id %in% early_blast_inactivated_genes)


subset(sandberg_escape_sim_res_early_blast, gene_id %in% subset(fData(Deng_sexed_cells), 
                                                                gene_short_name %in%  c("Xist", "Rlim", "Kif4", "Huwe1", "Atrx"))$gene_id)

fData(female_cds[row.names(subset(fData(Deng_sexed_cells), 
                                  gene_short_name %in%  c("Xist", "Rlim", "Kif4", "Huwe1", "Atrx"))),])


X_inact_targets <- Deng_sexed_cells[row.names(subset(fData(Deng_sexed_cells), gene_short_name %in%  c("Xist", "Rlim", "Kif4", "Huwe1", "Atrx"))),
                                    pData(Deng_sexed_cells)$stage == "early-blast" &
                                      pData(Deng_sexed_cells)$Deng_Sex == "Female" ]
sim_escape(X_inact_targets, formulaStr="~num_genes_expressed", allele_thresh=0.1, dispersion=0, nsims=1000, cores=detectCores() / 2)


pdf ("fig6d.pdf", width=1.5, height=4.25)
plot_genes_jitter(female_cds[row.names(subset(fData(Deng_sexed_cells), 
                                              gene_short_name %in%  c("Xist", "Rlim", "Kif4", "Huwe1", "Atrx"))),],
                  grouping="stage", color_by="allele", plot_trend=FALSE, min_expr=0, cell_size=I(0.15)) + 
  stat_summary(aes_string(x="stage", y="expression", color="allele", group="allele"), fun.data = "mean_cl_boot", size=0.35, geom="line") +
  stat_summary(aes_string(group="allele"), color="black", fun.data = "mean_cl_boot", size=0.75, geom="point") +
  stat_summary(aes_string(color="allele"), fun.data = "mean_cl_boot", size=0.25, geom="point") +
  scale_color_manual(values=c("#C51B7D", "#4D9221")) +
  scale_y_continuous() +  
  theme(legend.position="none") + 
  monocle:::monocle_theme_opts()
dev.off()

# Xist_maternal <- "ENSMUSG00000086503_M"
# Xist_paternal <- "ENSMUSG00000086503_P"

# Jpx_maternal <- "ENSMUSG00000097571_M"
# Jpx_paternal <- "ENSMUSG00000097571_P"

# paternal_ncrnas <- Deng_sexed_cells[c(Xist_paternal, Xist_maternal, Jpx_paternal, Jpx_maternal), pData(Deng_sexed_cells)$Deng_Sex == "Female"]
# fData(paternal_ncrnas)$gene_short_name = fData(paternal_ncrnas)$symbol
# 
# qplot(as.numeric(exprs(Deng_sexed_cells[c(Xist_paternal), pData(Deng_sexed_cells)$Deng_Sex == "Female"])),
#       as.numeric(exprs(Deng_sexed_cells[c(Jpx_paternal), pData(Deng_sexed_cells)$Deng_Sex == "Female"])))
# 
# monocle:::plot_coexpression_matrix(paternal_ncrnas, rowgenes=c("Xist"), colgenes=c("Jpx"))


all_informative_genes <- as.character(subset(fData(Deng_sexed_cells), is.allele.informative == 1)$gene_id)

all_informative_gene_records <- row.names(subset(fData(Deng_sexed_cells), gene_id %in% all_informative_genes 
                                                 & gene_id %in% subset(gencode_gene_metadata, biotype %in% allowed_gene_types)$gene_id))

all_informative_sim_res<- sim_escape(Deng_sexed_cells[all_informative_gene_records,], 
                                     formulaStr="~stage + Deng_Sex + num_genes_expressed", 
                                     nsims=100,
                                     cores=detectCores() / 2)

nrow(subset(all_informative_sim_res, (maternal_expr_qval < 0.05 | paternal_expr_qval < 0.05) &  monoallele_qval < 0.05))




escape_sim_res <- rbind(sandberg_escape_sim_res_4_cell, sandberg_escape_sim_res_16_cell, sandberg_escape_sim_res_early_blast)

escape_sim_res$stage <- factor(escape_sim_res$stage, levels=c("4-cell", "16-cell", "early-blast"))
escape_sim_res$gene_id <- row.names(escape_sim_res)


allelic_bias_genes <- subset(escape_sim_res, allelic_bias_qval < 0.05)$gene_id

subset(fData(Deng_sexed_cells), gene_id %in% allelic_bias_genes)



#nrow(all_informative_sim_res, )



#subset(fData(sandberg_escape_cds), gene_id %in% subset(sandberg_escape_sim_res, paternal_expr_qval < 0.05)$gene_id)


# escape_test_res <- allelic_bias_test(sandberg_escape_cds, fullModelFormulaStr="~embryo", cores=detectCores() / 2)
# escape_test_res <- merge(distinct(fData(sandberg)[,c("gene_id", "symbol", "chr")]), escape_test_res, by.x="gene_id", by.y="row.names")
# row.names(escape_test_res) <- escape_test_res$gene_id
# escape_test_res$qval <- 1.0
# escape_test_res$qval[escape_test_res$status == "OK"] <- p.adjust(subset(escape_test_res, status == "OK")$pval, method="fdr")

# x_escape_candidate_gene_ids <- setdiff(x_escape_candidate_gene_ids, subset(escape_test_res, qval < 0.05)$gene_id)

# sandberg_escape_cds <- Deng_sexed_cells[row.names(subset(fData(Deng_sexed_cells), gene_id %in% x_escape_candidate_gene_ids)),]


