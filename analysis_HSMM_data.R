# # use monocle2: 
# library(devtools)
# load_all('~/Projects/monocle-dev')
# #library(monocle)
# library(xacHelper)
# library(plyr)
# load_all_libraries()

# #export this function in xacHelper package: 
# optim_mc_func_fix_c <- function (m, c, t_estimate = estimate_t(TPM_isoform_count_cds),
#           relative_expr_matrix = relative_expr_matrix, split_relative_expr_matrix = split_relative_exprs,
#           alpha = rep(1, ncol(relative_expr_matrix)), total_RNAs = rep(50000, ncol(relative_expr_matrix)),
#           cores = 1, weight = 0.5, add_kl_divergence  = T, verbose = F,  ...) {
#   data('spike_df') #add the spikein dataset

#   if(is.null(spike_df$log_numMolecule))
#     spike_df$log_numMolecules <- log10(spike_df$numMolecules)
  
#   m_val <- m
#   c_val <- c
  
#   cell_num <- ncol(relative_expr_matrix)
#   names(t_estimate) <- colnames(relative_expr_matrix)
#   split_t <- split(t(t_estimate), col(as.matrix(t(t_estimate)), as.factor = T))
  
#   total_rna_df <- data.frame(Cell = colnames(relative_expr_matrix), t_estimate = t_estimate)
  
#   t_k_b_solution <- tryCatch({
#     k_b_solution <- plyr::ddply(total_rna_df, .(Cell), function(x) {
#       a_matrix <- matrix(c(log10(x[, "t_estimate"]), 1,
#                            m_val, -1), ncol = 2, nrow = 2, byrow = T)
#       colnames(a_matrix) <- c("k", "b")
#       b_matrix <- matrix(c(0, -c_val), nrow = 2, byrow = T)
#       k_b_solution <- t(solve(a_matrix, b_matrix))
#     })
#     k_b_solution},
#     error = function(e) {print(e); c(NA, NA)}
#   )
  
#   if(any(is.na(t_k_b_solution)))
#     return(NA)
  
#   cell_dmode <- tryCatch({
#     if(cores > 1){
#       cell_dmode <- mcmapply(opt_norm_t, split_t, split_relative_expr_matrix, m = m_val, c = c_val, pseudocnt = 0.01, mc.cores = cores)
#       adj_est_std_cds <- mcmapply(opt_norm_t, split_t, split_relative_expr_matrix, m = m_val, c = c_val, pseudocnt = 0.01, return_norm = T, mc.cores = cores)
#     }
#     else {
#       cell_dmode <- mapply(opt_norm_t, split_t, split_relative_expr_matrix, m = m_val, c = c_val, pseudocnt = 0.01)
#       adj_est_std_cds <- mapply(opt_norm_t, split_t, split_relative_expr_matrix, m = m_val, c = c_val, pseudocnt = 0.01, return_norm = T)
#     }
#     cell_dmode},
#     error = function(e) {print(e); NA}
#   )
  
#   if(any(is.na(cell_dmode)))
#     return(NA)
  
#   #adj_est_std_cds <- mapply(opt_norm_t, split_t, split_fpkm, m = m_val, c = c_val, return_norm = T)
#   sum_total_cells_rna <- colSums(adj_est_std_cds)
  
#   #minimization function:
#   #8.
#   dmode_rmse_weight_total <- mean(weight*((cell_dmode - alpha)/alpha)^2 + (1 - weight)*((sum_total_cells_rna -  total_RNAs)/total_RNAs)^2)
#   #add the JS distance measure:
#   split_relative_expr_matrix <- split(t(adj_est_std_cds), 1:ncol(adj_est_std_cds))
#   round_split_relative_expr_matrix <- split(t(round(adj_est_std_cds)), 1:ncol(adj_est_std_cds))
  
#   p_df <- makeprobs(relative_expr_matrix) #relative expression
#   p_list <- split(t(p_df), 1:ncol(p_df))
#   q_df_round <- makeprobs(round(adj_est_std_cds)) #round
#   q_df <- makeprobs(adj_est_std_cds) #no rounding
#   q_list <- split(t(q_df), 1:ncol(q_df))
#   q_list_round <- split(t(q_df_round), 1:ncol(q_df_round))
  
#   dist_divergence <- mcmapply(function(x, y) {
#     JSdistVec(x, y)
#   }, p_list, q_list, mc.cores = cores)
  
#   dist_divergence_round <- mcmapply(function(x, y) {
#     JSdistVec(x, y)
#   }, p_list, q_list_round, mc.cores = cores)
  
#   gm_dist_divergence <- exp(mean(log(dist_divergence)))
  
#   if(add_kl_divergence)
#     res <- 0.25 * log10(dmode_rmse_weight_total + 1) + 0.75 * gm_dist_divergence
#   else
#     res <- log10(dmode_rmse_weight_total + 1)
  
#   #use the algorithm:
#   if(add_kl_divergence)
#     res <- (weight * (mean(((cell_dmode - alpha)/alpha)^2) - 0)) + (1 - weight) * (gm_dist_divergence - 0.0) + dmode_rmse_weight_total
#   else
#     res <- log10(dmode_rmse_weight_total + 1)
  
#   if(verbose){
#     message('current m, c values are ', paste(m, c, sep = ', '))
#     message('dmode_rmse_weight_total is ', mean(((cell_dmode - alpha)/alpha)^2) - 0)
#     message('gm_dist_divergence is ', gm_dist_divergence)
#   }
#   #   return(list(m = m_val, c = c_val, dmode_rmse_weight_total = dmode_rmse_weight_total, gm_dist_divergence = gm_dist_divergence, dist_divergence_round = dist_divergence_round,
#   #               cell_dmode = cell_dmode, t_k_b_solution = t_k_b_solution, sum_total_cells_rna = sum_total_cells_rna, optim_res = res))
#   #
#   if(is.finite(dmode_rmse_weight_total))
#     return(res)
#   else
#     return(10)
# }

# ###################the muscle data####################
# #Cole's code to order the muscle cells##
# HSMM_fpkm_matrix <- read.delim("./data/HSMM_data/muscle/HSMM/HSMM_cuffnorm_out/genes.fpkm_table")
# row.names(HSMM_fpkm_matrix) <- HSMM_fpkm_matrix$tracking_id
# HSMM_fpkm_matrix <- HSMM_fpkm_matrix[,-1]

# HSMM_isoform_fpkm_matrix <- read.delim("./data/HSMM_data/muscle/HSMM/HSMM_cuffnorm_out/genes.fpkm_table")
# row.names(HSMM_isoform_fpkm_matrix) <- HSMM_isoform_fpkm_matrix$tracking_id
# HSMM_isoform_fpkm_matrix <- HSMM_isoform_fpkm_matrix[,-1]

# HSMM_readcount_matrix <- read.delim("./data/HSMM_data/muscle/HSMM/HSMM_cuffnorm_out/genes.count_table")
# row.names(HSMM_readcount_matrix) <- HSMM_readcount_matrix$tracking_id
# HSMM_readcount_matrix <- HSMM_readcount_matrix[,-1]

# sample_sheet <- read.delim("./data/HSMM_data/muscle/HSMM/sample_sheet.txt")
# sample_sheet$cell_id <- paste(sample_sheet$cell_id, "_0", sep="")
# row.names(sample_sheet) <- sample_sheet$cell_id
# sample_sheet <- sample_sheet[colnames(HSMM_fpkm_matrix),]

# cell_is_valid_singleton <- row.names(subset(sample_sheet, Control == FALSE & Unusual.Shape == FALSE & Debris == FALSE & Clump == FALSE & Cells.in.Well == 1)) 
# sample_sheet <- sample_sheet[cell_is_valid_singleton,]

# HSMM_fpkm_matrix <- HSMM_fpkm_matrix[,row.names(sample_sheet)]
# HSMM_isoform_fpkm_matrix <- HSMM_isoform_fpkm_matrix[,row.names(sample_sheet)]
# HSMM_readcount_matrix <- HSMM_readcount_matrix[,row.names(sample_sheet)]

# #use the new recovery algorithm: 
# TPM_HSMM_isoform_fpkm_matrix <- apply(HSMM_isoform_fpkm_matrix[1:5, 1:5], 2, function(x) x / sum(x) * 10^6)

# gene_ann <- read.delim("./data/HSMM_data/muscle/HSMM/gene_annotations.txt")

# mito_genes <- subset(gene_ann, grepl("chrM", gene_ann$locus))$gene_short_name

# gencode_biotypes <- read.delim("./data/HSMM_data/muscle/HSMM/gencode_biotypes.txt")

# gene_ann <- merge(gene_ann, gencode_biotypes, by = "gene_id")
# row.names(gene_ann) <- gene_ann$gene_id
# gene_ann <- gene_ann[,c("gene_short_name", "biotype")]
# gene_ann <- gene_ann[row.names(HSMM_fpkm_matrix),]
# #sample_sheet <- sample_sheet[colnames(fpkm_matrix_adj),]

# pd <- new("AnnotatedDataFrame", data = sample_sheet)
# fd <- new("AnnotatedDataFrame", data = gene_ann)

# HSMM_fpkm_matrix_cds <-  newCellDataSet(as.matrix(HSMM_fpkm_matrix), 
#                                    phenoData = pd, 
#                                    featureData = fd, 
#                                    expressionFamily=negbinomial(), 
#                                    lowerDetectionLimit=1)

# # HSMM_fpkm_matrix_adj_select <- relative2abs(HSMM_fpkm_matrix_cds, t_estimate = estimate_t(TPM_HSMM_isoform_fpkm_matrix, relative_expr_thresh = 0.1), 
# #                                                         alpha_v = 1, total_RNAs = 50000, weight = 0.01, verbose = T, return_all = T, cores = detectCores(), m =  -4.864207, c = 2.77514, c_rng = c(2.77514, 2.77514)) # mean(mean_m_c_select[1, ])

# HSMM_fpkm_matrix_adj <- relative2abs(HSMM_fpkm_matrix_cds, t_estimate = estimate_t(HSMM_isoform_fpkm_matrix), cores=detectCores())


# HSMM <-  newCellDataSet(as.matrix(HSMM_fpkm_matrix_adj), 
#                         phenoData = pd, 
#                         featureData = fd, 
#                         expressionFamily=negbinomial(), 
#                         lowerDetectionLimit=1)
# # HSMM_abs_select <-  newCellDataSet(as.matrix(HSMM_fpkm_matrix_adj_select$norm_cds), 
# #                                    phenoData = pd, 
# #                                    featureData = fd, 
# #                                    expressionFamily=negbinomial(), 
# #                                    lowerDetectionLimit=1)

# # exprs(HSMM) <- as.matrix(HSMM_fpkm_matrix_adj_select$norm_cds)

# ######################################analysis code for the HSMM dataset######################################
# pData(HSMM)$Total_mRNAs <- colSums(exprs(HSMM))

# #HSMM <- HSMM[,row.names(subset(pData(HSMM), Total_mRNAs >= 10000 & Total_mRNAs <= 100000))]
# HSMM <- detectGenes(HSMM, min_expr = 0.1)

# PDGFRA_expr <- exprs(HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("PDGFRA"))),])
# SPHK1_expr <- exprs(HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("SPHK1"))),])
# ANPEP_expr <- exprs(HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("ANPEP"))),])
# MEF2C_expr <- exprs(HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("MEF2C"))),])
# MYF5_expr <- exprs(HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("MYF5"))),])


# pData(HSMM)$PDGFRA <- as.vector(PDGFRA_expr)
# pData(HSMM)$SPHK1 <- as.vector(SPHK1_expr)
# pData(HSMM)$ANPEP <- as.vector(ANPEP_expr)
# pData(HSMM)$MEF2C <- as.vector(MEF2C_expr)
# pData(HSMM)$MYF5 <- as.vector(MYF5_expr)


# SPHK1_thresh <- 1
# PDGFRA_thresh <- 1
# ANPEP_thresh <- Inf #0.1
# MEF2C_thresh <- 5 #0.1
# MYF5_thresh <- 1 #0.1

# #cell_is_hsmm <- (pData(HSMM)$MEF2C > MEF2C_thresh | pData(HSMM)$MYF5 > MYF5_thresh) 
# cell_is_hsmm <- (pData(HSMM)$MEF2C > MEF2C_thresh | pData(HSMM)$MYF5 > MYF5_thresh) 

# names(cell_is_hsmm) <- row.names(pData(HSMM))
# pData(HSMM)$CellType <- "Myoblast"
# pData(HSMM)$CellType[cell_is_hsmm == FALSE] <-  "Fibroblast"

# HSMM_myo <- HSMM[,pData(HSMM)$CellType == "Myoblast" ] #& pData(HSMM)$Total_mRNAs > 10000 & pData(HSMM)$num_genes_expressed > 1000 & pData(HSMM)$Total_mRNAs < 40000  
# # HSMM_myo <- HSMM[, valid_HSMM_cell]

# HSMM_expressed_genes <- row.names(subset(fData(HSMM_myo), 
#                                          num_cells_expressed >= 15 & 
#                                            biotype %in% c("protein_coding", "lincRNA")))

# #HSMM_myo <- setOrderingFilter(HSMM_myo, HSMM_ordering_genes)

# HSMM_myo <- estimateSizeFactors(HSMM_myo)

# HSMM_myo <- estimateDispersions(HSMM_myo)

# #HSMM_myo <- computeVarianceStabilizedValues(HSMM_myo, model_row_names=HSMM_expressed_genes, modelFormulaStr="expression~Media", cores=12)
# HSMM_myo_grp_res <- differentialGeneTest(HSMM_myo[, ], 
#                                              fullModelFormulaStr = "~Media", #log10(Total_mRNAs) + spike_total_mRNAs
#                                              reducedModelFormulaStr = "~1", cores =detectCores(), relative = T)

# std_HSMM <-  newCellDataSet(as.matrix(HSMM_fpkm_matrix[row.names(HSMM_myo), colnames(HSMM_myo)]), 
#                             phenoData = new("AnnotatedDataFrame", data = pData(HSMM_myo)), 
#                             featureData = new("AnnotatedDataFrame", data = fData(HSMM_myo)), 
#                             expressionFamily=tobit(), 
#                             lowerDetectionLimit=1)
# readcounts_HSMM <-  newCellDataSet(as.matrix(HSMM_readcount_matrix[row.names(HSMM_myo), colnames(HSMM_myo)]), 
#                             phenoData = new("AnnotatedDataFrame", data = pData(HSMM_myo)), 
#                             featureData = new("AnnotatedDataFrame", data = fData(HSMM_myo)), 
#                             expressionFamily=tobit(), 
#                             lowerDetectionLimit=1)
# readcounts_HSMM <- detectGenes(readcounts_HSMM)
# readcounts_HSMM <- estimateSizeFactors(readcounts_HSMM)
# readcounts_HSMM <- estimateDispersions(readcounts_HSMM)

# std_HSMM <- detectGenes(std_HSMM, min_expr = 0.1)
# std_HSMM_grp_res <- differentialGeneTest(std_HSMM[, ], 
#                                              fullModelFormulaStr = "~Media", #log10(Total_mRNAs) + spike_total_mRNAs
#                                              reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

# HSMM_ordering_genes <- rownames(subset(HSMM_myo_grp_res, qval < 1e-5))  
# std_HSMM_ordering_genes <- rownames(subset(std_HSMM_grp_res, qval < 1e-5)) 

# # HSMM_ordering_genes <- selectGenesInExpressionRange(HSMM_myo, 5, Inf, 0.1, stat_fun=function(x) { median(round(x)) })

# HSMM_myo <- setOrderingFilter(HSMM_myo, c(HSMM_ordering_genes, std_HSMM_ordering_genes))

# HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, use_irlba=T, use_vst=T, pseudo_expr=0, fun="exp", scaling = F, method = "ICA")
# #HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, use_irlba=T, use_vst=F, pseudo_expr=0.5, fun="exp", scaling = F, method = "ICA")

# HSMM_myo <- orderCells(HSMM_myo, reverse=F, num_paths=1) #, root_cell = 'T48_CT_G10_0'

# #std of HSMM: 

# # HSMM_myo <- computeVarianceStabilizedValues(HSMM_myo, model_row_names=HSMM_expressed_genes, modelFormulaStr="expression~Media", cores=12)
# std_HSMM <- setOrderingFilter(std_HSMM, c(HSMM_ordering_genes, std_HSMM_ordering_genes))
 
# std_HSMM <- reduceDimension(std_HSMM, max_components = 2, use_irlba=T, use_vst=F, pseudo_expr=0.1, scaling = F, method = "ICA")
# #HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, use_irlba=T, use_vst=F, pseudo_expr=0.5, fun="exp", scaling = F, method = "ICA")

# std_HSMM <- orderCells(std_HSMM, reverse=F, num_paths=1) #, root_cell = 'T48_CT_G10_0'

# # #readcounts of HSMM: 

# # # HSMM_myo <- computeVarianceStabilizedValues(HSMM_myo, model_row_names=HSMM_expressed_genes, modelFormulaStr="expression~Media", cores=12)
# # std_HSMM <- setOrderingFilter(std_HSMM, c(HSMM_ordering_genes, std_HSMM_ordering_genes))
 
# # std_HSMM <- reduceDimension(std_HSMM, max_components = 2, use_irlba=T, use_vst=F, pseudo_expr=0.1, scaling = F, method = "ICA")
# # #HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, use_irlba=T, use_vst=F, pseudo_expr=0.5, fun="exp", scaling = F, method = "ICA")

# # std_HSMM <- orderCells(std_HSMM, reverse=F, num_paths=1) #, root_cell = 'T48_CT_G10_0'


# ################################################data for generating the muscle p-val for calculating the result################################################
# std_HSMM_myo_pseudotime_res_ori <- differentialGeneTest((std_HSMM[, ]), relative_expr = F, 
#                                                         fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", #log10(Total_mRNAs) + spike_total_mRNAs
#                                                         reducedModelFormulaStr = "~1", cores = detectCores())

# std_HSMM_myo <- std_HSMM[, colnames(HSMM_myo)]
# std_count_d <- DESeq::newCountDataSet(round(exprs(std_HSMM_myo)), pData(std_HSMM_myo)$Time) 
# #we should use read counts if we need to compare the relative expression performances
# std_HSMM_myo_DEseq_pseudotime_res <- DESeq1_pseudotime_test(std_count_d[, ], scale = F, Pseudotime = pData(std_HSMM_myo)$Pseudotime) 

# HSMM_myo_size_norm_res <- differentialGeneTest(HSMM_myo[, ], 
#                                                fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", #log10(Total_mRNAs) + spike_total_mRNAs
#                                                reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

# #create size normalized cds for permutation tests: 
# HSMM_myo_size_norm <- newCellDataSet(as.matrix(t(t(exprs(HSMM_myo[, ])) / sizeFactors(HSMM_myo))), 
#                                      phenoData = new("AnnotatedDataFrame", data = pData(HSMM_myo)), 
#                                      featureData = new("AnnotatedDataFrame", data = fData(HSMM_myo[, ])), 
#                                      expressionFamily = negbinomial(), 
#                                      lowerDetectionLimit = HSMM_myo@lowerDetectionLimit)
# HSMM_myo_size_norm <- estimateSizeFactors(HSMM_myo_size_norm)
# HSMM_myo_size_norm_count_d <- DESeq::newCountDataSet(round(exprs(HSMM_myo_size_norm)), pData(HSMM_myo_size_norm)$Time) 
# size_norm_HSMM_myo_DEseq_pseudotime_res <- DESeq1_pseudotime_test(HSMM_myo_size_norm_count_d[, ], scale = F, Pseudotime = pData(HSMM_myo_size_norm)$Pseudotime)

# muscle_std_glm_perm <- cal_glm_perm(std_HSMM[, ]) #for DESeq: the read count need to be used if we want to compare the relative expression values
# muscle_std_glm_perm <- do.call(rbind, muscle_std_glm_perm)
# muscle_std_glm_perm_pval <- muscle_std_glm_perm[, 1]

# muscle_size_normalized_mc_glm_perm <- cal_glm_perm(HSMM_myo_size_norm)
# muscle_size_normalized_mc_glm_perm <- do.call(rbind, muscle_size_normalized_mc_glm_perm)
# muscle_size_normalized_mc_glm_perm_pval <- muscle_size_normalized_mc_glm_perm[, 1]
# #save.image()

# #then generate the data for the plot: 
# monocle_p <- std_HSMM_myo_pseudotime_res_ori[, 'pval']
# names(monocle_p) <- row.names(std_HSMM_myo_pseudotime_res_ori)
# deseq_p <- std_HSMM_myo_DEseq_pseudotime_res$dtalbe$pval
# names(deseq_p) <- row.names(std_HSMM_myo_DEseq_pseudotime_res$dtalbe)
# HSMM_myo_size_norm_pval <- HSMM_myo_size_norm_res[, 'pval']
# names(HSMM_myo_size_norm_pval) <- row.names(HSMM_myo_size_norm_res)
# size_norm_HSMM_myo_DEseq_pseudotime_pval <- size_norm_HSMM_myo_DEseq_pseudotime_res$dtalbe$pval 
# names(size_norm_HSMM_myo_DEseq_pseudotime_pval) <- row.names(size_norm_HSMM_myo_DEseq_pseudotime_res$dtalbe)
# #
# size_norm_HSMM_myo_DEseq_pseudotime_pval <- size_norm_HSMM_myo_DEseq_pseudotime_res$dtalbe$pval 
# names(size_norm_HSMM_myo_DEseq_pseudotime_pval) <- row.names(size_norm_HSMM_myo_DEseq_pseudotime_res$dtalbe)

# muscle_std_glm_perm_pval <- muscle_std_glm_perm[, 1]
# muscle_mode_size_norm_permutate_ratio_by_geometric_mean <- muscle_size_normalized_mc_glm_perm[, 1]

# #choose the gene set: 
# element_all <- c(names(monocle_p[monocle_p <0.1]), 
#                names(HSMM_myo_size_norm_pval[HSMM_myo_size_norm_pval <0.1]))

# muscle_df <- plot_pre_rec_f1(test_p_list = list(std_monocle_p = monocle_p, std_deseq_p = deseq_p, size_norm_deseq_p = size_norm_HSMM_myo_DEseq_pseudotime_pval, size_norm_monocle_pval = HSMM_myo_size_norm_pval), 
#                              permutate_pval = list(std_monocle_p = muscle_std_glm_perm_pval, std_deseq_p = muscle_std_glm_perm_pval, size_norm_deseq_p = muscle_size_normalized_mc_glm_perm_pval, size_norm_monocle_pval = muscle_size_normalized_mc_glm_perm_pval), names(monocle_p),
#                              title = '', return_df = T) #+ ylim(0, 1) + ggtitle('') + theme(legend.position = 'NULL') + monocle_theme_opts() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
# muscle_df$data_type = c("Transcript (size normalization)", "Transcript (size normalization)", "FPKM", "FPKM")

# muscle_df$class = '3relative'


####################################################################################################################################
###make the pre/rec/F1 score and the ROC curves for the HSMM_pseudotime dataset: 

HSMM_pseudotime_pval_df <- data.frame(std_monocle_p = monocle_p, 
                          std_deseq_p = deseq_p[names(monocle_p)], 
                          size_norm_deseq_p = size_norm_HSMM_myo_DEseq_pseudotime_pval[names(monocle_p)], 
                          size_norm_monocle_pval = HSMM_myo_size_norm_pval[names(monocle_p)])

row.names(HSMM_pseudotime_pval_df) <- names(monocle_p)
permutation_HSMM_pseudotime_pval_df <- data.frame(std_monocle_p = muscle_std_glm_perm_pval[names(monocle_p)], 
                                      std_deseq_p = muscle_std_glm_perm_pval[names(monocle_p)], 
                                      size_norm_deseq_p = muscle_size_normalized_mc_glm_perm_pval[names(monocle_p)], 
                                      size_norm_monocle_pval = muscle_size_normalized_mc_glm_perm_pval[names(monocle_p)])
row.names(permutation_HSMM_pseudotime_pval_df) <- names(monocle_p)

generate_roc_df <-function(p_value, classification, type = 'fpr') {
  p_value[is.na(p_value)] <- 1
  pred_p_value <- prediction(p_value, classification)
  perf_tpr_fpr <- performance(pred_p_value, "tpr", "fpr")
  
    fpr = perf_tpr_fpr@x.values

    tpr = perf_tpr_fpr@y.values
    
  perf_auc <- performance(pred_p_value, "auc")
  auc <- perf_auc@y.values

    data.frame(tpr = tpr, fpr = fpr, auc = auc)
}

####
plot_roc_df <-function(p_value, classification, type = 'fpr') {
  p_value[is.na(p_value)] <- 1
  pred_p_value <- prediction(p_value, classification)
  perf_tpr_fpr <- performance(pred_p_value, "tpr", "fpr")
  
  pdf('roc_plot_ori.pdf')
  plot(perf_tpr_fpr)
  dev.off()
}

select_genes <- row.names(HSMM_myo)[esApply(HSMM_myo, 1, function(x) sum(x > 1) > 20)]

p_thrsld <- 0.01

perm_pvals <- permutation_HSMM_pseudotime_pval_df[select_genes, 1]
software_pvals <- HSMM_pseudotime_pval_df[select_genes, 1] 

perm_pvals[is.na(perm_pvals)] <- 1
software_pvals[is.na(software_pvals)] <- 1
plot_roc_df(software_pvals, perm_pvals > p_thrsld)
####

HSMM_pseudotime_roc_df_list <- lapply(colnames(HSMM_pseudotime_pval_df), function(x) {
  print(x)
  
  # if(x %in% c('mode_size_norm_permutate_ratio_by_geometric_mean', 'abs_default_edgeR_p', 'abs_default_deseq2_p', 'abs_scde_p', 'mast_abs_pval_no_norm'))
  #   select_genes <- row.names(new_abs_cds_14_18[1:transcript_num])[esApply(new_abs_cds_14_18[1:transcript_num], 1, function(x) sum(x > 0.5) > 10)]
  # if(x %in% c('monocle_p', 'mast_std_pval'))
  #   select_genes <- row.names(new_std_cds_14_18[1:transcript_num])[esApply(new_std_cds_14_18[1:transcript_num], 1, function(x) sum(x > 0.5) > 10)]
  # if(x %in% c('monocle_p_readcount', 'default_edgeR_p', 'default_deseq2_p', 'default_deseq_p', 'scde_p', 'mast_count_pval_no_norm'))
  #   select_genes <- row.names(count_cds[1:transcript_num])[esApply(count_cds[1:transcript_num], 1, function(x) sum(x > 0.5) > 10)]
  # if(x %in% c('mc_mode_size_norm_permutate_ratio_by_geometric_mean', 'mast_mc_pval_no_norm', 'mc_count_edgeR_p_glm', 'mc_count_deseq2_p', 'mc_count_deseq_p', 'mc_count_scde_p'))
  #   select_genes <- row.names(mc_adj_cds[1:transcript_num])[esApply(mc_adj_cds[1:transcript_num], 1, function(x) sum(x > 0.5) > 10)]

  perm_pvals <- permutation_HSMM_pseudotime_pval_df[select_genes, 1] #use the fpkm results as the gold standard
  software_pvals <- HSMM_pseudotime_pval_df[select_genes, x] 

  perm_pvals[is.na(perm_pvals)] <- 1
  software_pvals[is.na(software_pvals)] <- 1
  res <- generate_roc_df(software_pvals, perm_pvals > p_thrsld)
  colnames(res) <- c('tpr', 'fpr', 'auc')
  cbind(res, method = x)
})

HSMM_pseudotime_roc_df_list <- lapply(HSMM_pseudotime_roc_df_list, function(x) {colnames(x) <- c('tpr', 'fpr', 'auc', 'method'); x} )
HSMM_pseudotime_roc_df <- do.call(rbind, HSMM_pseudotime_roc_df_list)
HSMM_pseudotime_roc_df[1:5, ]
# str_split_fixed(row.names(HSMM_pseudotime_roc_df), '\\.', 2)[, 1]

cols <- c("FPKM" = "#F2756D", "Read counts" = "#6F94CC", "Transcript counts" = "#000202", "Estimated transcript counts" = "#7BAE41") ##A680B9

HSMM_pseudotime_auc <- unique(HSMM_pseudotime_roc_df[, c('method', 'auc')])
row.names(HSMM_pseudotime_auc) <- HSMM_pseudotime_auc$method

save(file = './HSMM_pseudotime_roc_df', HSMM_pseudotime_roc_df)

HSMM_pseudotime_roc_df$software <- revalue(HSMM_pseudotime_roc_df$method, c("std_monocle_p" = 'Monocle', "std_deseq_p" = 'DESeq1',
                                            "size_norm_deseq_p" = 'DESeq1', "size_norm_monocle_pval" = 'Monocle'))

HSMM_pseudotime_roc_df$Type <- revalue(HSMM_pseudotime_roc_df$method, c("std_monocle_p" = 'FPKM', "std_deseq_p" = 'FPKM',  
                                            "size_norm_deseq_p" = 'Estimated transcript counts', "size_norm_monocle_pval" = 'Estimated transcript counts'))

pdf('./supplementary_figures/HSMM_pseudotime_roc.pdf', height = 1.5, width = 3)
qplot(fpr, tpr, data= subset(HSMM_pseudotime_roc_df, software %in% c('Monocle', 'DESeq1')), geom="step", size = 0.5, color = Type) + #linetype = Type, 
    xlab("False positive rate") +
    ylab("True positive rate") +
      ylim(c(0, 1.0)) + facet_wrap(~software) + scale_size(range = c(0.1, 0.5)) + 
   scale_color_manual(values = cols, name = "Type") + 
  xlim(c(0, 1.0)) + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pdf('./supplementary_figures/HSMM_pseudotime_roc_helper.pdf', height = 13, width = 14)
qplot(fpr, tpr, data= HSMM_pseudotime_roc_df, geom="step", color = Type) + #linetype = Type, 
    xlab("False positive rate") +
    ylab("True positive rate") +
      ylim(c(0, 1.0)) + facet_wrap(~software) + 
   scale_color_manual(values = cols, name = "Type") #+ nm_theme()
  xlim(c(0, 1.0)) #+ nm_theme()
dev.off()


# HSMM_pseudotime_roc_df <- rbind(HSMM_pseudotime_roc_df, tmp) 

#roc_df <-  melt(df_res)
#                              software
# Type                          DESeq DESeq2 edgeR MAST Monocle SCDE
#   estimated transcript counts  1139   1239  1234 2484    1249 1207
#   FPKM                         1229      0     0 1242     996    0
#   read counts                     0   1289  1424    0    1313 1275

pdf('./supplementary_figures/HSMM_pseudotime_roc_dfroc_auc_bar.pdf', height = 3, width = 3)
ggplot(aes(software, auc), data = HSMM_pseudotime_roc_df) + geom_bar(position = 'dodge', stat = 'identity', aes(fill=Type)) + 
    xlab("") +
    # ylim(0.5, 1.0) + 
    # monocle_theme_opts() +  theme(axis.text.x=element_text(angle=30, hjust=1)) + 
    # scale_fill_manual(values = cols, name = "Software", label = test)  + 
    nm_theme()
dev.off()

pdf('./supplementary_figures/HSMM_pseudotime_roc_dfroc_auc_bar_helper.pdf', height = 6, width = 9)
ggplot(aes(software, auc), data = HSMM_pseudotime_roc_df) + geom_bar(position = 'dodge', stat = 'identity', aes(fill=Type)) + 
    xlab("") #+
    # ylim(0.5, 1.0) + 
    # monocle_theme_opts() +  theme(axis.text.x=element_text(angle=30, hjust=1)) + 
    # scale_fill_manual(values = cols, name = "Software", label = test)  + 
    # nm_theme()
dev.off()


#show the values of auc 
unique(HSMM_pseudotime_roc_df[, c('method','auc')])


# save.image('./RData/analysis_HSMM_data.RData')
