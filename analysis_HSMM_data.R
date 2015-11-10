 ###################the muscle data####################
  #Cole's code to order the muscle cells##
  HSMM_fpkm_matrix <- read.delim("/net/trapnell/vol1/home/xqiu/Projects/BEAM/Xiaojie_reproduce/muscle/HSMM/sc-RNA-Seq/HSMM_cuffnorm_out/genes.fpkm_table")
  row.names(HSMM_fpkm_matrix) <- HSMM_fpkm_matrix$tracking_id
  HSMM_fpkm_matrix <- HSMM_fpkm_matrix[,-1]
  
  HSMM_isoform_fpkm_matrix <- read.delim("/net/trapnell/vol1/home/xqiu/Projects/BEAM/Xiaojie_reproduce/muscle/HSMM/sc-RNA-Seq/HSMM_cuffnorm_out/isoforms.fpkm_table")
  row.names(HSMM_isoform_fpkm_matrix) <- HSMM_isoform_fpkm_matrix$tracking_id
  HSMM_isoform_fpkm_matrix <- HSMM_isoform_fpkm_matrix[,-1]
  
  sample_sheet <- read.delim("/net/trapnell/vol1/home/xqiu/Projects/BEAM/Xiaojie_reproduce/muscle/HSMM/sc-RNA-Seq/sample_sheet.txt")
  sample_sheet$cell_id <- paste(sample_sheet$cell_id, "_0", sep="")
  row.names(sample_sheet) <- sample_sheet$cell_id
  sample_sheet <- sample_sheet[colnames(HSMM_fpkm_matrix),]
  
  cell_is_valid_singleton <- row.names(subset(sample_sheet, Control == FALSE & Unusual.Shape == FALSE & Debris == FALSE & Clump == FALSE & Cells.in.Well == 1)) 
  sample_sheet <- sample_sheet[cell_is_valid_singleton,]
  
  HSMM_fpkm_matrix <- HSMM_fpkm_matrix[,row.names(sample_sheet)]
  HSMM_isoform_fpkm_matrix <- HSMM_isoform_fpkm_matrix[,row.names(sample_sheet)]
    
  #use the new recovery algorithm: 
  TPM_HSMM_isoform_fpkm_matrix <- apply(HSMM_isoform_fpkm_matrix, 2, function(x) x / sum(x) * 10^6)
  
  HSMM_fpkm_matrix_adj_select <- relative2abs_optim_fix_c(HSMM_fpkm_matrix, t_estimate = estimate_t(TPM_HSMM_isoform_fpkm_matrix, relative_expr_thresh = 0.1), 
                                                          alpha_v = 1, total_RNAs = 50000, weight = 0.01, verbose = T, return_all = T, cores = 2, m =  -4.864207, c = 2.77514) # mean(mean_m_c_select[1, ])
  gene_ann <- read.delim("/net/trapnell/vol1/home/xqiu/Projects/BEAM/Xiaojie_reproduce/Muscle/HSMM/sc-RNA-Seq/gene_annotations.txt")
  
  mito_genes <- subset(gene_ann, grepl("chrM", gene_ann$locus))$gene_short_name
  
  gencode_biotypes <- read.delim("/net/trapnell/vol1/home/xqiu/Projects/BEAM/Xiaojie_reproduce/Muscle/HSMM/sc-RNA-Seq/gencode_biotypes.txt")
  
  gene_ann <- merge(gene_ann, gencode_biotypes, by = "gene_id")
  row.names(gene_ann) <- gene_ann$gene_id
  gene_ann <- gene_ann[,c("gene_short_name", "biotype")]
  gene_ann <- gene_ann[row.names(HSMM_fpkm_matrix),]
  #sample_sheet <- sample_sheet[colnames(fpkm_matrix_adj),]
  
  pd <- new("AnnotatedDataFrame", data = sample_sheet)
  fd <- new("AnnotatedDataFrame", data = gene_ann)
  
  HSMM_fpkm_matrix_cds <-  newCellDataSet(HSMM_fpkm_matrix, 
                                     phenoData = pd, 
                                     featureData = fd, 
                                     expressionFamily=negbinomial(), 
                                     lowerDetectionLimit=1)

  HSMM_fpkm_matrix_adj <- relative2abs(HSMM_fpkm_matrix_cds, t_estimate = estimate_t(HSMM_isoform_fpkm_matrix), cores=1)


  HSMM <-  newCellDataSet(HSMM_fpkm_matrix_adj, 
                          phenoData = pd, 
                          featureData = fd, 
                          expressionFamily=negbinomial(), 
                          lowerDetectionLimit=1)
  HSMM_abs_select <-  newCellDataSet(HSMM_fpkm_matrix_adj_select$norm_cds, 
                                     phenoData = pd, 
                                     featureData = fd, 
                                     expressionFamily=negbinomial(), 
                                     lowerDetectionLimit=1)
  qplot(esApply(HSMM, 2, sum), apply(HSMM_fpkm_matrix_adj_select$norm_cds, 2, sum), log = 'xy') + xlab('Previous recovery algorithm') +
    ylab('New recovery algorithm') + ggtitle('Compare the result from previous and current recovery algorithm (HSMM data)') + geom_smooth(method = 'rlm') + geom_abline()
  
  # valid_HSMM_cell <- load('valid_HSMM_cell') #load the cells
  
  if(use_select_algorithm)
    exprs(HSMM) <- HSMM_fpkm_matrix_adj_select$norm_cds
  
    ######################################analysis code for the HSMM dataset######################################
  pData(HSMM)$Total_mRNAs <- colSums(exprs(HSMM))
  
  #HSMM <- HSMM[,row.names(subset(pData(HSMM), Total_mRNAs >= 10000 & Total_mRNAs <= 100000))]
  HSMM <- detectGenes(HSMM, min_expr = 0.1)
  
  # myo_and_fibro_expressed_genes <- row.names(subset(fData(HSMM),
  # num_cells_expressed >= 15 &
  # biotype %in% c("protein_coding", "lincRNA")))
  # vstExprs(HSMM) <- computeVarianceStabilizedValues(HSMM, model_row_names=myo_and_fibro_expressed_genes, modelFormulaStr="expression~Media", cores=12)
  
  PDGFRA_expr <- exprs(HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("PDGFRA"))),])
  SPHK1_expr <- exprs(HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("SPHK1"))),])
  ANPEP_expr <- exprs(HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("ANPEP"))),])
  MEF2C_expr <- exprs(HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("MEF2C"))),])
  MYF5_expr <- exprs(HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("MYF5"))),])
  
  
  pData(HSMM)$PDGFRA <- as.vector(PDGFRA_expr)
  pData(HSMM)$SPHK1 <- as.vector(SPHK1_expr)
  pData(HSMM)$ANPEP <- as.vector(ANPEP_expr)
  pData(HSMM)$MEF2C <- as.vector(MEF2C_expr)
  pData(HSMM)$MYF5 <- as.vector(MYF5_expr)
  
  
  SPHK1_thresh <- 1
  PDGFRA_thresh <- 1
  ANPEP_thresh <- Inf #0.1
  MEF2C_thresh <- 5 #0.1
  MYF5_thresh <- 1 #0.1
  
  #cell_is_hsmm <- (pData(HSMM)$MEF2C > MEF2C_thresh | pData(HSMM)$MYF5 > MYF5_thresh) 
  cell_is_hsmm <- (pData(HSMM)$MEF2C > MEF2C_thresh | pData(HSMM)$MYF5 > MYF5_thresh) 
  
  names(cell_is_hsmm) <- row.names(pData(HSMM))
  pData(HSMM)$CellType <- "Myoblast"
  pData(HSMM)$CellType[cell_is_hsmm == FALSE] <-  "Fibroblast"
  
  HSMM_myo <- HSMM[,pData(HSMM)$CellType == "Myoblast" ] #& pData(HSMM)$Total_mRNAs > 10000 & pData(HSMM)$num_genes_expressed > 1000 & pData(HSMM)$Total_mRNAs < 40000  
  # HSMM_myo <- HSMM[, valid_HSMM_cell]
  
  HSMM_expressed_genes <- row.names(subset(fData(HSMM_myo), 
                                           num_cells_expressed >= 15 & 
                                             biotype %in% c("protein_coding", "lincRNA")))
  
  #HSMM_myo <- setOrderingFilter(HSMM_myo, HSMM_ordering_genes)
  
  HSMM_myo <- estimateSizeFactors(HSMM_myo)
  
  HSMM_myo <- estimateDispersions(HSMM_myo)
  
  #HSMM_myo <- computeVarianceStabilizedValues(HSMM_myo, model_row_names=HSMM_expressed_genes, modelFormulaStr="expression~Media", cores=12)
  HSMM_myo_grp_res <- differentialGeneTest(HSMM_myo[, ], 
                                               fullModelFormulaStr = "~Media", #log10(Total_mRNAs) + spike_total_mRNAs
                                               reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

  std_HSMM <-  newCellDataSet(HSMM_fpkm_matrix[row.names(HSMM_myo), colnames(HSMM_myo)], 
                              phenoData = new("AnnotatedDataFrame", data = pData(HSMM_myo)), 
                              featureData = new("AnnotatedDataFrame", data = fData(HSMM_myo)), 
                              expressionFamily=tobit(), 
                              lowerDetectionLimit=1)
  std_HSMM <- detectGenes(std_HSMM, min_expr = 0.1)
  std_HSMM_grp_res <- differentialGeneTest(std_HSMM[, ], 
                                               fullModelFormulaStr = "~Media", #log10(Total_mRNAs) + spike_total_mRNAs
                                               reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

  HSMM_ordering_genes <- rownames(subset(HSMM_myo_grp_res, qval < 1e-5))  
  std_HSMM_ordering_genes <- rownames(subset(std_HSMM_grp_res, qval < 1e-5)) 

  # HSMM_ordering_genes <- selectGenesInExpressionRange(HSMM_myo, 5, Inf, 0.1, stat_fun=function(x) { median(round(x)) })
  
  HSMM_myo <- setOrderingFilter(HSMM_myo, c(HSMM_ordering_genes, std_HSMM_ordering_genes))
  
  HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, use_irlba=T, use_vst=T, pseudo_expr=0, fun="exp")
  #HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, use_irlba=T, use_vst=F, pseudo_expr=0.5, fun="exp")
  
  HSMM_myo <- orderCells(HSMM_myo, reverse=F, num_paths=1, root_cell = 'T48_CT_G10_0')
  
  #std of HSMM: 

  
  #order the cells with the fpkm data: 
#  std_HSMM <- estimateSizeFactors(std_HSMM)
  
#  std_HSMM <- estimateDispersions(std_HSMM, cores=detectCores())

  # std_HSMM_ordering_genes <- selectGenesInExpressionRange(std_HSMM, 5, Inf, 0.1, stat_fun=function(x) { median(round(x)) })
  # HSMM_myo <- computeVarianceStabilizedValues(HSMM_myo, model_row_names=HSMM_expressed_genes, modelFormulaStr="expression~Media", cores=12)
  std_HSMM <- setOrderingFilter(std_HSMM, c(HSMM_ordering_genes, std_HSMM_ordering_genes))
   
  std_HSMM <- reduceDimension(std_HSMM, max_components = 2, use_irlba=T, use_vst=F, pseudo_expr=0.1)
  #HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, use_irlba=T, use_vst=F, pseudo_expr=0.5, fun="exp")
  
  std_HSMM <- orderCells(std_HSMM, reverse=F, num_paths=1, root_cell = 'T48_CT_G10_0')
  
  ################################################data for generating the muscle p-val for calculating the result################################################
  std_HSMM_myo_pseudotime_res_ori <- differentialGeneTest((std_HSMM[, ]), relative_expr = F, 
                                                          fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", #log10(Total_mRNAs) + spike_total_mRNAs
                                                          reducedModelFormulaStr = "~1", cores = detectCores())
  
  std_HSMM_myo <- std_HSMM[, colnames(HSMM_myo)]
  std_count_d <- DESeq::newCountDataSet(round(exprs(std_HSMM_myo)), pData(std_HSMM_myo)$Time) 
  #we should use read counts if we need to compare the relative expression performances
  std_HSMM_myo_DEseq_pseudotime_res <- DESeq1_pseudotime_test(std_count_d[, ], scale = F, Pseudotime = pData(std_HSMM_myo)$Pseudotime) 
  
  HSMM_myo_size_norm_res <- differentialGeneTest(HSMM_myo[, ], 
                                                 fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", #log10(Total_mRNAs) + spike_total_mRNAs
                                                 reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)
  
  # save.image()
  #create size normalized cds for permutation tests: 
  HSMM_myo_size_norm <- newCellDataSet(t(t(exprs(HSMM_myo[, ])) / sizeFactors(HSMM_myo)), 
                                       phenoData = new("AnnotatedDataFrame", data = pData(HSMM_myo)), 
                                       featureData = new("AnnotatedDataFrame", data = fData(HSMM_myo[, ])), 
                                       expressionFamily = negbinomial(), 
                                       lowerDetectionLimit = HSMM_myo@lowerDetectionLimit)
  HSMM_myo_size_norm <- estimateSizeFactors(HSMM_myo_size_norm)
  HSMM_myo_size_norm_count_d <- DESeq::newCountDataSet(round(exprs(HSMM_myo_size_norm)), pData(HSMM_myo_size_norm)$Time) 
  size_norm_HSMM_myo_DEseq_pseudotime_res <- DESeq1_pseudotime_test(HSMM_myo_size_norm_count_d[, ], scale = F, Pseudotime = pData(HSMM_myo_size_norm)$Pseudotime)
  
  muscle_std_glm_perm <- cal_glm_perm(std_HSMM[, ]) #for DESeq: the read count need to be used if we want to compare the relative expression values
  muscle_std_glm_perm <- do.call(rbind, muscle_std_glm_perm)
  muscle_std_glm_perm_pval <- muscle_std_glm_perm[, 1]
  
  muscle_size_normalized_mc_glm_perm <- cal_glm_perm(HSMM_myo_size_norm)
  muscle_size_normalized_mc_glm_perm <- do.call(rbind, muscle_size_normalized_mc_glm_perm)
  muscle_size_normalized_mc_glm_perm_pval <- muscle_size_normalized_mc_glm_perm[, 1]
  save.image()
  
  #then generate the data for the plot: 
  monocle_p <- std_HSMM_myo_pseudotime_res_ori[, 'pval']
  names(monocle_p) <- row.names(std_HSMM_myo_pseudotime_res_ori)
  deseq_p <- std_HSMM_myo_DEseq_pseudotime_res$dtalbe$pval
  names(deseq_p) <- row.names(std_HSMM_myo_DEseq_pseudotime_res$dtalbe)
  HSMM_myo_size_norm_pval <- HSMM_myo_size_norm_res[, 'pval']
  names(HSMM_myo_size_norm_pval) <- row.names(HSMM_myo_size_norm_res)
  size_norm_HSMM_myo_DEseq_pseudotime_pval <- size_norm_HSMM_myo_DEseq_pseudotime_res$dtalbe$pval 
  names(size_norm_HSMM_myo_DEseq_pseudotime_pval) <- row.names(size_norm_HSMM_myo_DEseq_pseudotime_res$dtalbe)
  #
  size_norm_HSMM_myo_DEseq_pseudotime_pval <- size_norm_HSMM_myo_DEseq_pseudotime_res$dtalbe$pval 
  names(size_norm_HSMM_myo_DEseq_pseudotime_pval) <- row.names(size_norm_HSMM_myo_DEseq_pseudotime_res$dtalbe)
  
  muscle_std_glm_perm_pval <- muscle_std_glm_perm[, 1]
  muscle_mode_size_norm_permutate_ratio_by_geometric_mean <- muscle_size_normalized_mc_glm_perm[, 1]
  
  #choose the gene set: 
  element_all <- c(names(monocle_p[monocle_p <0.1]), 
                 names(HSMM_myo_size_norm_pval[HSMM_myo_size_norm_pval <0.1]))
  
  muscle_df <- plot_pre_rec_f1(test_p_list = list(std_monocle_p = monocle_p, std_deseq_p = deseq_p, size_norm_deseq_p = size_norm_HSMM_myo_DEseq_pseudotime_pval, size_norm_monocle_pval = HSMM_myo_size_norm_pval), 
                               permutate_pval = list(std_monocle_p = muscle_std_glm_perm_pval, std_deseq_p = muscle_std_glm_perm_pval, size_norm_deseq_p = muscle_size_normalized_mc_glm_perm_pval, size_norm_monocle_pval = muscle_size_normalized_mc_glm_perm_pval), names(monocle_p),
                               title = '', return_df = T) #+ ylim(0, 1) + ggtitle('') + theme(legend.position = 'NULL') + monocle_theme_opts() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
    #why the precision is so low? 
  muscle_df$data_type = c("Transcript (size normalization)", "Transcript (size normalization)", "FPKM", "FPKM")
  
  muscle_df$class = '3relative'