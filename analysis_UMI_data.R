  #UMI data: 
  #from the paper GSE54695_analysis_script.R: 
  library(monocle)
  library(xacHelper)

  load_all_libraries()

  input.ERCC.annotation<-read.delim("./data/Quake_data/quake_lung/ERCC_specification.txt", header=T)
  colnames(input.ERCC.annotation)<-c("Resort_ID",
                                     "ERCC_ID",
                                     "subgroup",
                                     "conc_attomoles_ul_Mix1",
                                     "conc_attomoles_ul_Mix2",
                                     "exp_fch_ratio",
                                     "log2_Mix1_Mix2")
  
  # So we can index this data frame by ERCC transcript ID                  
  rownames(input.ERCC.annotation)<-input.ERCC.annotation[,"ERCC_ID"]

  ##rm(list = setdiff(ls(), lsf.str()))
  # load('/Users/xqiu/Dropbox (Personal)/Quake/scRNA-seq_confirm_algorithm/SRP030617/input.ERCC.annotation')
  umi_matrix <- read.delim('./data/UMI_data/GSE54695_data_transcript_counts.txt', row.names="GENENAME")
  
  # pdf('umi_sum_dist.pdf')
  # qplot(apply(umi_matrix, 2, sum), log = 'x')
  # dev.off()

  ERCC_ids <- (grep('^ERCC', row.names(umi_matrix)))
  input.ERCC.annotation[row.names(umi_matrix)[ERCC_ids], ]
  input.ERCC.annotation$numMolecules <- input.ERCC.annotation$conc_attomoles_ul_Mix1*(20*10^(-3)*1/50000*10^(-18)*6.02214179*10^(23))
  input.ERCC.annotation[setdiff(row.names(input.ERCC.annotation), row.names(umi_matrix)[ERCC_ids]), ] #anything below 3.528598705 are removed
  sort(unique(input.ERCC.annotation[setdiff(row.names(input.ERCC.annotation), row.names(umi_matrix)[ERCC_ids]), 'numMolecules'])) 
  #(ERCC-00123, ERCC-00104, ERCC-00142, ERCC-00024, ERCC-00016, ERCC-00081, ERCC-00041)
  #input.ERCC.annotation[c('ERCC-00123', 'ERCC-00104', 'ERCC-00142', 'ERCC-00024', 'ERCC-00016', 'ERCC-00081', 'ERCC-00041'), ]
  # umi_matrix[c('ERCC-00123', 'ERCC-00104', 'ERCC-00142', 'ERCC-00024', 'ERCC-00016', 'ERCC-00081', 'ERCC-00041'), ]
  # rowSums(umi_matrix[c('ERCC-00123', 'ERCC-00104', 'ERCC-00142', 'ERCC-00024', 'ERCC-00016', 'ERCC-00081', 'ERCC-00041'), ]) #very small number of counts

  sample_sheet_tmp <- do.call(rbind.data.frame, strsplit(colnames(umi_matrix), '_'))[, 1:2]
  colnames(sample_sheet_tmp) <- c('Type', "Condition")
  sample_sheet <- cbind(sample_sheet_tmp, group = apply(sample_sheet_tmp[ , 1:2] , 1 , paste , collapse = "_" ))
  row.names(sample_sheet) <- colnames(umi_matrix)

  pd <- new("AnnotatedDataFrame", data = sample_sheet)
  type <- rep('gene', nrow(umi_matrix))
  type[ERCC_ids] <- 'ERCC'
  fd <- new("AnnotatedDataFrame", data = data.frame(gene = row.names(umi_matrix), type = type, row.names = row.names(umi_matrix)))
  umi_matrix <- as.matrix(umi_matrix)
  umi_matrix <- umi_matrix[,row.names(sample_sheet)]

  UMI_cds <- newCellDataSet(umi_matrix, 
                            expressionFamily=negbinomial(), 
                            phenoData = pd, 
                            featureData = fd,
                            lowerDetectionLimit=1)

  pData(UMI_cds)$Total_mRNAs <- esApply(UMI_cds, 2, sum)
  summary(pData(UMI_cds)$Total_mRNAs)

  #show capture efficiency distribution: 
  UMI_cds_total_mRNAs <- UMI_cds[, pData(UMI_cds)$Total_mRNAs > 4800]
  umi_spike_df <- data.frame(spikein = input.ERCC.annotation[row.names(umi_matrix)[ERCC_ids], 'numMolecules'], UMI = exprs(UMI_cds_total_mRNAs)[ERCC_ids, 1])
  
  # pdf('recovery_efficiency.pdf')
  # qplot(UMI, spikein, data = umi_spike_df, log = 'xy') #recovery efficiency is around 0.1
  # dev.off()

  all_umi_spike_df <- data.frame(spikein = input.ERCC.annotation[row.names(umi_matrix)[ERCC_ids], 'numMolecules'], UMI = as.vector(exprs(UMI_cds)[ERCC_ids, ]))
  
  # pdf('recovery_efficiency_loess.pdf')
  # qplot(UMI, spikein, data = all_umi_spike_df, log = 'xy') + geom_smooth(method = 'loess') #recovery efficiency is around 0.1 (mixture 1/2?)
  # dev.off()

  UMI_cds <- UMI_cds[, pData(UMI_cds)$Total_mRNAs > 0]
  pData(UMI_cds)$dmode <- estimate_t(UMI_cds)

  # pdf('mode_cell_umi.pdf')
  # qplot(unlist(pData(UMI_cds)$dmode), log = 'x') + xlab('Mode of UMI for each cell') #confirm the first assumption
  # dev.off()
  # 
  # fig3b <- c('Cd19', 'Cd79b', 'Cd22', 'Cd37', 'Ctsd', 'Apoe', 'C1qa', 'C1qb', 'C1qc', 'Csf1r', 'Slpi', 'Tlr2', 'Mmp13', 'Marco', 'Ifng', 'Gzmb', 'Myc', 'Xcl1', 'Ccl5', 'Gzma', 'Nkg7', 'Spic', 'Cebpb', 'Lyz2', 'Sfpi1', 'Nfkbiz', 'Bst2', 'Siglech', 'Ly6d', 'Irf8', 'Cst3', 'Naaa', 'Ccr7', 'Cxcl9', 'Traf1', 'Relb', 'Itgax', 'Tmem176b', 'Tnf', 'Tnfaip3', 'Nfkbia', 'Il15', 'Cxcl10', 'Ifit1', 'Isg15', 'Irf7')
  # fig4c <- c('Irf1', 'Ly6e', 'Ppt1', 'Stat1', 'Irf2', 'Ccl4', 'Stat2')
  # figs10 <- c('V00821', 'Hmox1', 'Cd7', 'Tyrobp')
  # figs12 <- c("Isg15","Cxcl10","Irf7","Stat2", "Oas3", "Ifit1", "Ifit2", "Ifit3", "Mx1", "Mx2", "Cd69", "Ly6a", "Cd274", "Tnf", "Nfkbia", "Myd88", "Nmi", "Oas2", "Traf1", "Nfkb2", "Etv3", "Tnfaip3", "Ccl5", "Jak1", "Jak2", "Irf1", "Il15", "Ifi205","Stat1", "Ly6c2", "Cebpb", "Irf9", "Il1b", "Tlr2", "Rnd3", "Id2", "Cdh1", "Jarid2", "Relb", "Ccr5", "Fscn1","Ccl22", "Nfkbib", "Lyz2", "Marco", "Csf1", "Stat3", "Junb", "Bst2", "Nfkbiz", "Stat5a", "Il4i1", "Ccr7", "Vcam1","Adap2","Atf4","Irf2","Xbp1","Myc","Ifng","Runx1","Apoe","Cd9","Tnfsf14","Tnfrsf9","Nkg7","Emr4","Hmox1")
  # figs13 <- c(“Plbd1”, “Psap”, “Crip1”, “Snora41”, “Txn1”, “Ppt1”, “Irf8”, “Lgals3”, “Atpif1”, “Ifi30”, “CCl5”, “Fscn1”, “II4i1”, “Nfkbia”, “Marcksl1”, “Gadd45b”, “Relb”, “Mdh2”, “Ifi27l2a”, “Itgb2”, “Cotl1”, “Ffar2”, “Gsn”, “Pglyrp1”, “Ppp1r14a”, “Dtx1”, “Lyz2”)
  # 
  UMI_cds_subset <-  UMI_cds[, pData(UMI_cds)$Total_mRNAs > 10000 & pData(UMI_cds)$group %in% c('SC_2i', 'SC_serum')]
  UMI_cds_subset <- estimateSizeFactors(UMI_cds_subset)
  UMI_cds_subset@expressionFamily <- tobit() #fix the bug of estimateDispersions
  UMI_ordering_genes <- selectGenesInExpressionRange(UMI_cds_subset, 5, Inf, 0.1, stat_fun=function(x) { median(round(x)) })
  UMI_cds_subset <- setOrderingFilter(UMI_cds_subset, UMI_ordering_genes)
  UMI_cds_subset <- reduceDimension(UMI_cds_subset, use_irlba = F, use_vst = F) 
  UMI_cds_subset <- orderCells(UMI_cds_subset, num_paths = 2, reverse = F) #SRR1033962_0 
  UMI_cds_subset@expressionFamily <- negbinomial()

  # pdf('UMI_tree.pdf')
  # plot_spanning_tree(UMI_cds_subset, color_by="group", show_backbone=T, show_cell_names = F)
  # dev.off()

  #permform the permuation test and benchmarking? 
  #use CD11c+CD8+CD86+ VS CD11c+CD8+pDC as control: 
  UMI_cds_76_55 <- UMI_cds_subset[, pData(UMI_cds_subset)$group %in% c('SC_2i', 'SC_serum')]
  #convert UMI counts to TPM values: 
  UMI_cds_76_55_relative <- newCellDataSet(apply(exprs(UMI_cds_76_55), 2, function(x) x / sum(x) * 10^6 ), 
                                           expressionFamily=tobit(), 
                                           phenoData = new("AnnotatedDataFrame", pData(UMI_cds_76_55)), 
                                           featureData = new("AnnotatedDataFrame", fData(UMI_cds_76_55)),
                                           lowerDetectionLimit=1)
  Marker_ori <- as.character(pData(UMI_cds_76_55)$group)
  Marker_order <- order(Marker_ori) #order the E14.5 cell at the begining (ensure E14.5 cells number is 43 while E18.5 is 74)

  UMI_cds_76_55_perm_res <- cal_perm_pval_size_norm(UMI_cds_76_55[, Marker_order], alpha = 76, beta = 55, grp0 = 'SC_2i', grp1 = 'SC_serum', group = Marker_ori)
  UMI_cds_76_55_relative_perm_res <- cal_perm_pval_size_norm(UMI_cds_76_55_relative[, Marker_order], alpha = 76, beta = 55, grp0 = 'SC_2i', grp1 = 'SC_serum', group = Marker_ori, size_norm = F)

  # UMI_split_cds <- split(t(exprs(UMI_cds_76_55[, Time_order])), col(t(exprs(UMI_cds_76_55[, Time_order])), as.factor = T))
  # UMI_fc <- esApply(UMI_cds_76_55[, ], 1, mean_fc, grp0 = 'SC_2i', grp1 = 'SC_serum', grp = Marker_ori)
  # UMI_split_fc <- split(t(UMI_fc), col(t(UMI_fc), as.factor = T))
  # UMI_permutate_pval <- mcmapply(permuation_pval, std_split_cds, std_split_fc, alpha = 76, beta = 55, mc.cores = detectCores()) #multiple cores 
  UMI_cds_76_55 <- estimateSizeFactors(UMI_cds_76_55)
  UMI_cds_76_55_count_d <- newCountDataSet(round(t(t(exprs(UMI_cds_76_55)) / sizeFactors(UMI_cds_76_55))), (Marker_ori)) #normalized the data by size factor
  # abs_dtable_pool_max_nbinomTest <- DESeq1_test(abs_count_d, disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomTest') 
  # row.names(abs_dtable_pool_max_nbinomTest$dtalbe) <- abs_dtable_pool_max_nbinomTest$dtalbe$id

  #DESeq glm: (GLM tests are more relevant to our software)
  UMI_cds_76_55_dtable_pool_max_nbinomGLMTest <- DESeq1_test(UMI_cds_76_55_count_d, disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomGLMTest', scale = T) 
  row.names(UMI_cds_76_55_dtable_pool_max_nbinomGLMTest$dtalbe) <-  row.names(UMI_cds_76_55_count_d)

  # calculate the pval with the readcount with scde: (calculate the scde associate DEG test result LOCALLY) 
  UMI_cds_76_55_scde_res_list <- scde_DEG(dir = NULL, count_cds = UMI_cds_76_55, DEG_attribute = 'group', contrast = c('SC_2i', 'SC_serum'), n.cores = detectCores())

  UMI_cds_76_55_diff_test_res <- differentialGeneTest(UMI_cds_76_55[, ], 
                                                      fullModelFormulaStr = "~group", 
                                                      reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)
  UMI_cds_76_55_relative_diff_test_res <- differentialGeneTest(UMI_cds_76_55_relative[, ], 
                                                               fullModelFormulaStr = "~group", 
                                                               reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)
  #calculate the pval with the normalized transcripts with scde: 
  #   abs_scde_res_list <- scde_DEG(dir = NULL, count_cds = UMI_cds_76_55, DEG_attribute = 'group', contrast = c('SC_2i', 'SC_serum'), n.cores = detectCores(), normalize = T)
  #   abs_scde_res_list_no_normalize <- scde_DEG(dir = NULL, count_cds = UMI_cds_76_55, DEG_attribute = 'Time', contrast = c('SC_2i', 'SC_serum'), n.cores = detectCores())
  #calculate the pval with DESeq 1/2, edgeR: 
  #DESeq2/edgeR: 

  # UMI_cds_76_55_deseq2_res <- DESeq2_deg(dir = NULL, UMI_cds_76_55, Time = Marker_ori, pd = pData(UMI_cds_76_55), design = as.formula("~ group"))
  UMI_cds_76_55_edgeR_res <- edgeR_test(exprs(UMI_cds_76_55), group = Marker_ori, glm = T)

  #benchmark the performance: 
  abs_monocle_p <- UMI_cds_76_55_diff_test_res[, 'pval'] 
  names(abs_monocle_p) <- row.names(UMI_cds_76_55_diff_test_res)

  monocle_p <- UMI_cds_76_55_relative_diff_test_res[, 'pval'] 
  names(monocle_p) <- row.names(UMI_cds_76_55_relative_diff_test_res)

  #deseq
  abs_default_deseq_p <- UMI_cds_76_55_dtable_pool_max_nbinomGLMTest$dtalbe[, 'pval'] #abs_dtable_pool_max_nbinomTest
  names(abs_default_deseq_p) <- row.names(UMI_cds_76_55_dtable_pool_max_nbinomGLMTest$dtalbe)
  #scde
  abs_scde_p <- UMI_cds_76_55_scde_res_list$pval #_no_normalize

  df3 <- plot_pre_rec_f1(test_p_list = list(monocle_p = monocle_p,
                                            abs_monocle_p = abs_monocle_p, 
                                            abs_default_deseq_p = abs_default_deseq_p, 
                                            abs_scde_p = abs_scde_p),
                         permutate_pval = list(monocle_p = UMI_cds_76_55_relative_perm_res,   
                                               abs_monocle_p = UMI_cds_76_55_perm_res,
                                               abs_default_deseq_p = UMI_cds_76_55_perm_res, 
                                               # abs_default_deseq_p_new_norm = mode_size_norm_permutate_ratio_by_geometric_mean, 
                                               abs_scde_p = UMI_cds_76_55_perm_res),
                         row.names(UMI_cds_76_55), #gene_list, overlap_genes, high_gene_list
                         return_df = T, #na.rm = T, 
                         title = 'Comparison of the two-group DEG tests on relative gene expression', 
                         rownames = c('monocle (TPM)', 'monocle (UMI)', 'DESeq (UMI)', 'SCDE (UMI)'))

  df3$class = '3relative'
  colnames(df3)[1:3] <- c('Precision', 'Recall', 'F1 score')
  df3$class <- c('UMI', 'UMI', 'UMI', 'TPM')


  pdf(file = "./data/supplementary_figure/figsc_UMI_GSE54695.pdf", width = 2.5, height = 2)
  ggplot(aes(factor(Type), value,  fill = class), data = melt(df3)) + geom_bar(position = position_dodge(), stat = 'identity') + #facet_wrap(~variable) + 
    ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
    ggtitle('') + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) + nm_theme() + xlab('') + ylim(0, 1)
  dev.off()


