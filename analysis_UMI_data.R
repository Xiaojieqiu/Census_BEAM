# cal_perm_pval_size_norm <- function (cds = new_abs_cds_14_18_2, alpha = sum(pData(cds_size_norm_cds)$Time ==
#     "E14.5"), beta = sum(pData(cds_size_norm_cds)$Time == "E18.5"),
#     grp0 = "E14.5", grp1 = "E18.5", group = pData(cds_size_norm_cds)$Time,
#     size_norm = T)
# {
#     if (size_norm) {
#         cds_norm_mat <- t(t(round(exprs(cds)/sizeFactors(cds)))) #round after the size factor normalization 
#     }
#     else cds_norm_mat = round(exprs(cds))
#     cds_size_norm_cds <- newCellDataSet(cds_norm_mat, phenoData = new("AnnotatedDataFrame",
#         data = pData(cds)), featureData = new("AnnotatedDataFrame",
#         data = fData(cds)), expressionFamily = negbinomial(),
#         lowerDetectionLimit = 0.1)
#     size_norm_permutate_pval <- permu_two_group_gen(cds_size_norm_cds[,
#         ], round = T, alpha = alpha, beta = beta, grp0 = grp0,
#         grp1 = grp1, group = group)
#     return(size_norm_permutate_pval)
# }

# edgeR_test <- function (counts = exprs(count_cds), sf = sizeFactors(count_cds), group = Time_ori, glm = F, normalize = T) {
#     print(paste('size factors are ', sf))
#     y <- DGEList(counts = round(counts), group = group)
#     if (!glm) {
        
#         if(normalize){
#           y$sample$norm.factors <- sf
#         }
#         else
#             y <- calcNormFactors(y)

#         y <- estimateCommonDisp(y)
#         y <- estimateTagwiseDisp(y)
#         et <- exactTest(y)
#         return(list(et = et, topTags = topTags(et)))
#     }
#     else {
#         # y <- calcNormFactors(y)
#         if(normalize){
#           y$sample$norm.factors <- sf
#         }
#         design <- model.matrix(~group)
#         y <- estimateGLMCommonDisp(y, design)
#         y <- estimateGLMTrendedDisp(y, design)
#         y <- estimateGLMTagwiseDisp(y, design)
#         fit <- glmFit(y, design)
#         lrt <- glmLRT(fit, coef = 2)
#         topTags(lrt)
#         return(list(lrt = lrt, topTags = topTags(lrt)))
#     }
# }

# #update DESeq 1: 
# DESeq1_test <- function (count_d, pd = pData(new_abs_cds_14_18[, colnames(count_d)]), condition = Time_ori, disp_method = "blind",
#     sharing = "fit-only", test_type = "nbinomGLMTest", scale = F,
#     contrast = conditions, fullModelFormulaStr = "count~condition",
#     reducedModelFormulaStr = "count~1"){
#     message(class(count_d))
#     if (class(count_d)[1] != "CountDataSet")
#         count_d <- newCountDataSet(round(exprs(cds)), condition)
#     d <- estimateSizeFactors(count_d)
#     if (!scale)
#         pData(d)$sizeFactor <- 1
#     else
#         pData(d)$sizeFactor <- pd[colnames(count_d), 'Size_Factor']

#     d2 <- d
#     sizeFactor <- pData(d)$sizeFactor
#     print(sizeFactor)
#     d <- DESeq::estimateDispersions(d2, method = disp_method,
#         sharingMode = sharing)
#     message("pass here")
#     message("estimate dispersion")
#     tmp <- do.call(cbind.data.frame, lapply(d@fitInfo, function(x) x[c(1,
#         3)]))
#     disp_func <- lapply(ls(d@fitInfo), function(x) fitInfo(d,
#         name = x))
#     message("all disp functions extracted")
#     if (test_type == "nbinomGLMTest") {
#         message("test is ", test_type)
#         print(DESeq::sizeFactors(d))
#         dfit1 <- fitNbinomGLMs(d, as.formula(fullModelFormulaStr))
#         dfit0 <- fitNbinomGLMs(d, as.formula(reducedModelFormulaStr))
#         dpval <- nbinomGLMTest(dfit1, dfit0)
#         dpadj <- p.adjust(dpval, method = "BH")
#         dtable <- transform(dfit1, pval = dpval, padj = dpadj)
#         dtable <- cbind(dtable, tmp)
#         head(dtable)
#         message("pass glm test")
#         return(list(dfit1 = dfit1, dfit0 = dfit0, full_df.residual = attr(dfit1,
#             "df.residual"), reduced_df.residual = attr(dfit0,
#             "df.residual"), full_deviance = dfit1$deviance, reduced_deviance = dfit0$deviance,
#             dtalbe = dtable, disp_func = disp_func))
#     }
#     else if (test_type == "nbinomTest") {
#         dtable <- nbinomTest(d, contrast[1], contrast[2])
#         dtable <- cbind(dtable, tmp)
#         message("pass nb test")
#         return(list(dtalbe = dtable, disp_func = disp_func))
#     }
# }

# #update DESeq 2: 
# DESeq2_deg <- function (dir, cds, Time, test_type = "LRT", design = as.formula("~ Time"),
#     reduced = as.formula("~ 1"), pd = pData(cds)) {
#     pd$sizeFactor <- pd$Size_Factor
#     if (!is.null(dir)) {
#         sample_table <- read.delim(paste(dir, "/samples.table",
#             sep = ""))
#         norm_count <- read.delim(paste(dir, "/genes.count_table",
#             sep = ""))
#         row.names(norm_count) <- norm_count$tracking_id
#         norm_count <- norm_count[, -1]
#         countdata <- round(t(t(norm_count) * sample_table$internal_scale))
#         if (exists("sample_sheet"))
#             countdata <- countdata[, row.names(sample_sheet)]
#     }
#     else countdata <- round(exprs(cds))
#     DEseq2_cnt_cds <- DESeqDataSetFromMatrix(countData = countdata[,
#         ], colData = pd, design = design)
#     DEseq2_cnt_cds <- DESeq(DEseq2_cnt_cds, test = test_type,
#         reduced = reduced)
#     return(DEseq2_cnt_cds)
# }

# #UMI data: 
# #from the paper GSE54695_analysis_script.R: 
# library(devtools)
# load_all('~/Projects/monocle-dev')
# # library(monocle)
# library(xacHelper)
# library(ROCR)

# load_all_libraries()

# input.ERCC.annotation<-read.delim("./data/Quake_data/quake_lung/ERCC_specification.txt", header=T)
# colnames(input.ERCC.annotation)<-c("Resort_ID",
#                                    "ERCC_ID",
#                                    "subgroup",
#                                    "conc_attomoles_ul_Mix1",
#                                    "conc_attomoles_ul_Mix2",
#                                    "exp_fch_ratio",
#                                    "log2_Mix1_Mix2")

# # So we can index this data frame by ERCC transcript ID                  
# rownames(input.ERCC.annotation)<-input.ERCC.annotation[,"ERCC_ID"]

# umi_matrix <- read.delim('./data/UMI_data/GSE54695_data_transcript_counts.txt', row.names="GENENAME")

# ERCC_ids <- (grep('^ERCC', row.names(umi_matrix)))
# input.ERCC.annotation <- input.ERCC.annotation[row.names(umi_matrix)[ERCC_ids], ]
# input.ERCC.annotation$numMolecules <- input.ERCC.annotation$conc_attomoles_ul_Mix1*(1000*10^(-3)*1/2500000*10^(-18)*6.02214179*10^(23))
# input.ERCC.annotation[setdiff(row.names(input.ERCC.annotation), row.names(umi_matrix)[ERCC_ids]), ] #anything below 3.528598705 are removed
# sort(unique(input.ERCC.annotation[setdiff(row.names(input.ERCC.annotation), row.names(umi_matrix)[ERCC_ids]), 'numMolecules'])) 
# #(ERCC-00123, ERCC-00104, ERCC-00142, ERCC-00024, ERCC-00016, ERCC-00081, ERCC-00041)
# #input.ERCC.annotation[c('ERCC-00123', 'ERCC-00104', 'ERCC-00142', 'ERCC-00024', 'ERCC-00016', 'ERCC-00081', 'ERCC-00041'), ]
# # umi_matrix[c('ERCC-00123', 'ERCC-00104', 'ERCC-00142', 'ERCC-00024', 'ERCC-00016', 'ERCC-00081', 'ERCC-00041'), ]
# # rowSums(umi_matrix[c('ERCC-00123', 'ERCC-00104', 'ERCC-00142', 'ERCC-00024', 'ERCC-00016', 'ERCC-00081', 'ERCC-00041'), ]) #very small number of counts

# sample_sheet_tmp <- do.call(rbind.data.frame, strsplit(colnames(umi_matrix), '_'))[, 1:2]
# colnames(sample_sheet_tmp) <- c('Type', "Condition")
# sample_sheet <- cbind(sample_sheet_tmp, group = apply(sample_sheet_tmp[ , 1:2] , 1 , paste , collapse = "_" ))
# row.names(sample_sheet) <- colnames(umi_matrix)

# pd <- new("AnnotatedDataFrame", data = sample_sheet)
# type <- rep('gene', nrow(umi_matrix))
# type[ERCC_ids] <- 'ERCC'
# fd <- new("AnnotatedDataFrame", data = data.frame(gene = row.names(umi_matrix), type = type, row.names = row.names(umi_matrix)))
# umi_matrix <- as.matrix(umi_matrix)
# umi_matrix <- umi_matrix[,row.names(sample_sheet)]

# UMI_cds <- newCellDataSet(umi_matrix, 
#                           expressionFamily=negbinomial(), 
#                           phenoData = pd, 
#                           featureData = fd,
#                           lowerDetectionLimit=1)

# pData(UMI_cds)$Total_mRNAs <- esApply(UMI_cds, 2, sum)
# summary(pData(UMI_cds)$Total_mRNAs)

# #show capture efficiency distribution: 
# UMI_cds_total_mRNAs <- UMI_cds[, pData(UMI_cds)$Total_mRNAs > 4800]
# umi_spike_df <- data.frame(spikein = input.ERCC.annotation[row.names(umi_matrix)[ERCC_ids], 'numMolecules'], UMI = exprs(UMI_cds_total_mRNAs)[ERCC_ids, 1])

# all_umi_spike_df <- data.frame(spikein = input.ERCC.annotation[row.names(umi_matrix)[ERCC_ids], 'numMolecules'], UMI = as.vector(exprs(UMI_cds)[ERCC_ids, ]))

# UMI_cds <- UMI_cds[, pData(UMI_cds)$Total_mRNAs > 0]
# pData(UMI_cds)$dmode <- estimate_t(exprs(UMI_cds))

# UMI_cds_subset <-  UMI_cds[, pData(UMI_cds)$Total_mRNAs > 10000 & pData(UMI_cds)$group %in% c('SC_2i', 'SC_serum')]
# UMI_cds_subset <- estimateSizeFactors(UMI_cds_subset)
# UMI_cds_subset@expressionFamily <- tobit() #fix the bug of estimateDispersions
# UMI_ordering_genes <- selectGenesInExpressionRange(UMI_cds_subset, 5, Inf, 0.1, stat_fun=function(x) { median(round(x)) })
# UMI_cds_subset <- setOrderingFilter(UMI_cds_subset, UMI_ordering_genes)
# UMI_cds_subset <- reduceDimension(UMI_cds_subset, use_irlba = F, use_vst = F, method = "ICA", scaling = F) 
# UMI_cds_subset <- orderCells(UMI_cds_subset, num_paths = 2, reverse = F) #SRR1033962_0 
# UMI_cds_subset@expressionFamily <- negbinomial()

# #permform the permuation test and benchmarking? 
# #use CD11c+CD8+CD86+ VS CD11c+CD8+pDC as control: 
# UMI_cds_76_55 <- UMI_cds_subset[, pData(UMI_cds_subset)$group %in% c('SC_2i', 'SC_serum')]

# #convert UMI counts to TPM values: 
# UMI_cds_76_55_relative <- newCellDataSet(apply(exprs(UMI_cds_76_55), 2, function(x) x / sum(x) * 10^6 ), 
#                                          expressionFamily=tobit(), 
#                                          phenoData = new("AnnotatedDataFrame", pData(UMI_cds_76_55)), 
#                                          featureData = new("AnnotatedDataFrame", fData(UMI_cds_76_55)),
#                                          lowerDetectionLimit=1)
# Marker_ori <- as.character(pData(UMI_cds_76_55)$group)
# Marker_order <- order(Marker_ori) #order the E14.5 cell at the begining (ensure E14.5 cells number is 43 while E18.5 is 74)

# UMI_cds_76_55 <- estimateSizeFactors(UMI_cds_76_55)
# UMI_cds_76_55 <- estimateDispersions(UMI_cds_76_55)

# UMI_sf_mat <- matrix(rep(sizeFactors(UMI_cds_76_55), nrow(UMI_cds_76_55)), nrow = nrow(UMI_cds_76_55), byrow = T, dimnames = dimnames(UMI_cds_76_55))

# UMI_cds_76_55_perm_res <- cal_perm_pval_size_norm(UMI_cds_76_55[, Marker_order], alpha = 76, beta = 55, grp0 = 'SC_2i', grp1 = 'SC_serum', group = Marker_ori)
# UMI_cds_76_55_relative_perm_res <- cal_perm_pval_size_norm(UMI_cds_76_55_relative[, Marker_order], alpha = 76, beta = 55, grp0 = 'SC_2i', grp1 = 'SC_serum', group = Marker_ori, size_norm = F)

# UMI_cds_76_55_count_d <- newCountDataSet(round(exprs(UMI_cds_76_55)), (Marker_ori)) # / sizeFactors(UMI_cds_76_55)) normalized the data by size factor

# # abs_dtable_pool_max_nbinomTest <- DESeq1_test(abs_count_d, disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomTest') 
# # row.names(abs_dtable_pool_max_nbinomTest$dtalbe) <- abs_dtable_pool_max_nbinomTest$dtalbe$id

# #DESeq glm: (GLM tests are more relevant to our software)
# UMI_cds_76_55_dtable_pool_max_nbinomGLMTest <- DESeq1_test(UMI_cds_76_55_count_d,  pd = pData(UMI_cds_76_55), disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomGLMTest', scale = T) 
# row.names(UMI_cds_76_55_dtable_pool_max_nbinomGLMTest$dtalbe) <-  row.names(UMI_cds_76_55_count_d)

# # calculate the pval with the readcount with scde: (calculate the scde associate DEG test result LOCALLY) 
# UMI_cds_76_55_scde_res_list <- scde_DEG(dir = NULL, count_cds = UMI_cds_76_55, DEG_attribute = 'group', contrast = c('SC_2i', 'SC_serum'), n.cores = detectCores(), normalize = T)

# UMI_cds_76_55_diff_test_res <- differentialGeneTest(UMI_cds_76_55[, ], 
#                                                     fullModelFormulaStr = "~group", 
#                                                     reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)
# UMI_cds_76_55_relative_diff_test_res <- differentialGeneTest(UMI_cds_76_55_relative[, ], 
#                                                              fullModelFormulaStr = "~group", 
#                                                              reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)
# #calculate the pval with the normalized transcripts with scde: 
# #   abs_scde_res_list <- scde_DEG(dir = NULL, count_cds = UMI_cds_76_55, DEG_attribute = 'group', contrast = c('SC_2i', 'SC_serum'), n.cores = detectCores(), normalize = T)
# #   abs_scde_res_list_no_normalize <- scde_DEG(dir = NULL, count_cds = UMI_cds_76_55, DEG_attribute = 'Time', contrast = c('SC_2i', 'SC_serum'), n.cores = detectCores())
# #calculate the pval with DESeq 1/2, edgeR: 
# #DESeq2/edgeR: 

# UMI_cds_76_55_deseq2_res <- DESeq2_deg(dir = NULL, UMI_cds_76_55, Time = Marker_ori, pd = pData(UMI_cds_76_55), design = as.formula("~ group"))
# UMI_cds_76_55_edgeR_res <- edgeR_test(exprs(UMI_cds_76_55), sf = sizeFactors(UMI_cds_76_55), group = Marker_ori, glm = T)

# #benchmark the performance: 
# abs_monocle_p <- UMI_cds_76_55_diff_test_res[, 'pval'] 
# names(abs_monocle_p) <- row.names(UMI_cds_76_55_diff_test_res)

# monocle_p <- UMI_cds_76_55_relative_diff_test_res[, 'pval'] 
# names(monocle_p) <- row.names(UMI_cds_76_55_relative_diff_test_res)

# #deseq
# abs_default_deseq_p <- UMI_cds_76_55_dtable_pool_max_nbinomGLMTest$dtalbe[, 'pval'] #abs_dtable_pool_max_nbinomTest
# names(abs_default_deseq_p) <- row.names(UMI_cds_76_55_dtable_pool_max_nbinomGLMTest$dtalbe)

# #deseq2
# abs_default_deseq2_p <- results(UMI_cds_76_55_deseq2_res)$pvalue #abs_dtable_pool_max_nbinomTest
# names(abs_default_deseq2_p) <- row.names(results(UMI_cds_76_55_deseq2_res))

# #edgeR:
# abs_default_edgeR_p <- UMI_cds_76_55_edgeR_res$lrt$table$PValue #abs_dtable_pool_max_nbinomTest
# names(abs_default_edgeR_p) <- row.names(UMI_cds_76_55_edgeR_res$rt$table)

# #scde
# abs_scde_p <- UMI_cds_76_55_scde_res_list$pval #_no_normalize

# df3 <- plot_pre_rec_f1(test_p_list = list(monocle_p = monocle_p,
#                                           abs_monocle_p = abs_monocle_p, 
#                                           abs_default_deseq_p = abs_default_deseq_p, 
#                                           abs_scde_p = abs_scde_p),
#                        permutate_pval = list(monocle_p = UMI_cds_76_55_relative_perm_res,   
#                                              abs_monocle_p = UMI_cds_76_55_perm_res,
#                                              abs_default_deseq_p = UMI_cds_76_55_perm_res, 
#                                              # abs_default_deseq_p_new_norm = mode_size_norm_permutate_ratio_by_geometric_mean, 
#                                              abs_scde_p = UMI_cds_76_55_perm_res),
#                        row.names(UMI_cds_76_55), #gene_list, overlap_genes, high_gene_list
#                        return_df = T, #na.rm = T, 
#                        title = 'Comparison of the two-group DEG tests on relative gene expression', 
#                        rownames = c('monocle (TPM)', 'monocle (UMI)', 'DESeq (UMI)', 'SCDE (UMI)'))

# df3$class = '3relative'
# colnames(df3)[1:3] <- c('Precision', 'Recall', 'F1 score')
# df3$class <- c('UMI', 'UMI', 'UMI', 'TPM')

# cols <- c("Read counts" = "#00BFC4","MC transcripts" = "#7CAE00", "Spikein transcripts" = "#C77CFF", "FPKM" = "#F8766D", 'read_counts' = "#00BFC4", transcript_counts = '#C77CFF', UMI = '#C77CFF', TPM = '#F8766D')
# pdf(file = "./supplementary_figures/fig2c_si.pdf", width = 2.5, height = 2)
# ggplot(aes(factor(Type), value,  fill = class), data = melt(df3)) + geom_bar(position = position_dodge(), stat = 'identity') + #facet_wrap(~variable) + 
#   ggtitle(title) +  scale_fill_manual(values = cols) + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') + 
#   ggtitle('') + theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 30, hjust = .9)) + theme(strip.background = element_blank()) + nm_theme() + xlab('') + ylim(0, 1)
# dev.off()


# #########add the test based on spike-in transcript counts: 

# pData(UMI_cds)$Total_mRNAs <- esApply(UMI_cds, 2, sum)
# pData(UMI_cds)$endogenous_RNA <- esApply(UMI_cds[-ERCC_ids, ], 2, sum)

# #1. test the normalization algorithm by using the TPM dataset converted from UMI data

# #convert UMI data to TPM (a relative abundance measure) for the following analysis: 
# UMI_TPM <- esApply(UMI_cds, 2, function(x) x / sum(x) * 10^6)
# # UMI_TPM <- UMI_TPM[, esApply(UMI_cds, 2, sum) > 5e4] #remove cells have terrible recovery as well as TPM distribution 
# UMI_tpm_align <- melt(UMI_TPM)
# UMI_tpm_align <- UMI_tpm_align[UMI_tpm_align$value > 0.1, ]

# pdf('./tmp/UMI_tpm_align_distr.pdf', width = 20, height = 20) 
# qplot(value, log = 'x', geom = 'density', color = Var2, data = UMI_tpm_align) + facet_wrap(~Var2) + theme_bw() + theme(legend.position = 'none') + xlab('TPM')
# dev.off()

# UMI_TPM_cds <- newCellDataSet(UMI_TPM, 
#                               phenoData = new("AnnotatedDataFrame", data = pData(UMI_cds[, colnames(UMI_TPM)])),
#                               featureData = new("AnnotatedDataFrame", data = fData(UMI_cds)),
#                               expressionFamily=tobit(), 
#                               lowerDetectionLimit=1)
# pData(UMI_TPM_cds)$Total_mRNAs <- esApply(UMI_TPM_cds, 2, sum)
# pData(UMI_TPM_cds)$endogenous_RNA <- esApply(UMI_TPM_cds[-ERCC_ids, ], 2, sum)

# UMI_TPM_ercc_controls <- newCellDataSet(exprs(UMI_TPM_cds[ERCC_ids, ]), 
#                                         phenoData = new("AnnotatedDataFrame", data = pData(UMI_cds[, colnames(UMI_TPM)])),
#                                         featureData = new("AnnotatedDataFrame", data = fData(UMI_cds[ERCC_ids, ])), 
#                                         expressionFamily=tobit(), 
#                                         lowerDetectionLimit=1)

# #We use the UMI counts as the ladder instead of ERCC spike-in transcripts: 
# #1. Only 50 ERCC species detected by their experiment (92 species was added)
# #2. Cell-seq with UMI is lower in recovery efficiency 
# #3. Authors from the paper mentioned low recovery efficiency based on spike-in regression (~ 0.025) and potential ERCC degradation

# #subset dataframe for ERCC annotation: 
# ERCC_ids_names <- row.names(UMI_cds)[ERCC_ids]
# subset_spike_df <- spike_df[ERCC_ids_names, ] #spike_df from the monocle package 
# identical(ERCC_ids_names, row.names(subset_spike_df))

# #there are variations between total UMI counts of the detected ERCC species for each cell, so median is used 
# subset_spike_df$numMolecules <- input.ERCC.annotation[row.names(subset_spike_df), 'numMolecules'] #as.numeric(esApply(UMI_cds[ERCC_ids, ], 1, median)) + 1

# #regression model between TPM and UMI counts 
# UMI_molModels <- esApply(UMI_TPM_ercc_controls, 2, function(cell_exprs, input.ERCC.annotation) {
#   #print (cell_exprs)
#   spike_df <- input.ERCC.annotation 
#   spike_df <- cbind(spike_df, cell_exprs[row.names(spike_df)]) 
#   colnames(spike_df)[length(colnames(spike_df))] <- "FPKM"
#   #spike_df$numMolecules <- spike_df$conc_attomoles_ul_Mix1*(1000*10^(-3)*1/2500000*10^(-18)*6.02214179*10^(23))
#   #spike_df$rounded_numMolecules <- round(spike_df$conc_attomoles_ul_Mix1*(1000*10^(-3)*1/2500000*10^(-18)*6.02214179*10^(23)))
#   # print (mean(spike_df$rounded_numMolecules))
  
#   spike_df <- subset(spike_df, FPKM >= 1e-10)
#   spike_df$log_fpkm <- log10(spike_df$FPKM) 
#   spike_df$log_numMolecules <- log10(spike_df$numMolecules)
  
#   # print (paste("geometric mean of log numMol ", mean(spike_df$log_numMolecules), "geometric mean of log FPKM ", mean(spike_df$log_fpkm)))
  
#   molModel <- tryCatch({
#     #molModel <- vgam (rounded_numMolecules ~ sm.ns(log_fpkm, df=3), data=spike_df, family=negbinomial(zero=NULL))
#     molModel <- rlm(log_numMolecules ~ log_fpkm, data=spike_df)
    
#     #
#     #        qp <- qplot(FPKM, numMolecules, data=spike_df, log="xy") +
#     #      geom_abline(color="green") +
#     # geom_smooth(aes(y=10^log_numMolecules), method="lm") +
#     # geom_smooth(aes(y=10^log_numMolecules), method="rlm", color="red") +
#     # geom_line(aes(y=10^predict(molModel,type="response")))
#     # print(qp)
#     molModel
#   }, 
#   error = function(e) { print(e); NULL })
#   molModel
# }, subset_spike_df)

# # #calculate the corresponding transcript counts of the mode with the model learnt from the UMI data (from TPM to UMI)
# t_estimate <- estimate_t(UMI_TPM_cds)
# UMI_mode_absolute_counts <- mapply(function(cell_exprs, molModel) {
#   tryCatch({
#     norm_df <- data.frame(log_fpkm=log10(as.numeric(cell_exprs)))
#     res <- 10^predict(molModel, type="response", newdata=norm_df)
#   }, 
#   error = function(e) {
#     rep(NA, length(cell_exprs))
#   })
# }, 
# split(as.vector(t_estimate), names(t_estimate)), 
# UMI_molModels) 

# UMI_absolute_counts <- mapply(function(cell_exprs, molModel) {
#   tryCatch({
#     norm_df <- data.frame(log_fpkm=log10(as.numeric(cell_exprs)))
#     res <- 10^predict(molModel, type="response", newdata=norm_df)
#   }, 
#   error = function(e) {
#     rep(NA, length(cell_exprs))
#   })
# }, 
# split(exprs(UMI_TPM_cds), rep(1:ncol(exprs(UMI_TPM_cds)), each = nrow(exprs(UMI_TPM_cds)))), 
# UMI_molModels) 

# dimnames(UMI_absolute_counts) <- dimnames(UMI_TPM_cds)
# pdf('./supplementary_figures/UMI_mode_dist.pdf', width = 2, height = 2)
# qplot(UMI_mode_absolute_counts) + xlab('Transcript count for most frequent log10(TPM)') + nm_theme() + ylab('Cells')
# dev.off()

# #obtain the regression parameters from the UMI_molModels
# UMI_models_kd_df <- do.call(rbind.data.frame, lapply(UMI_molModels, function(x) data.frame(b = coef(x)[1], k = coef(x)[2])))

# pdf('./supplementary_figures/UMI_kb_line.pdf', width = 2, height = 2)
# qplot(k, b, data = UMI_models_kd_df) + geom_smooth(method = 'rlm') + 
#   xlab('Slope from TPM vs. Median UMI counts \n for ERCC spike-in transcripts') + ylab('Intercept from TPM vs. Median UMI counts \n for ERCC spike-in transcripts') +  nm_theme()
# dev.off()

# rlm(b ~ k, data = UMI_models_kd_df) 

# ####################################################################################################################################

# ###perform the benchmark analysis: 
# spike_cds_76_55 <- newCellDataSet(UMI_absolute_counts[row.names(UMI_cds_76_55), colnames(UMI_cds_76_55)], 
#                                          expressionFamily=negbinomial(), 
#                                          phenoData = new("AnnotatedDataFrame", pData(UMI_cds_76_55)), 
#                                          featureData = new("AnnotatedDataFrame", fData(UMI_cds_76_55)),
#                                          lowerDetectionLimit=1)

# spike_cds_76_55 <- estimateSizeFactors(spike_cds_76_55)
# spike_cds_76_55 <- estimateDispersions(spike_cds_76_55)

# spike_sf_mat <- matrix(rep(sizeFactors(spike_cds_76_55), nrow(spike_cds_76_55)), nrow = nrow(spike_cds_76_55), byrow = T, dimnames = dimnames(spike_cds_76_55))

# spike_cds_76_55_perm_res <- cal_perm_pval_size_norm(spike_cds_76_55[, Marker_order], alpha = 76, beta = 55, grp0 = 'SC_2i', grp1 = 'SC_serum', group = Marker_ori)
# spike_cds_76_55_relative_perm_res <- cal_perm_pval_size_norm(spike_cds_76_55_relative[, Marker_order], alpha = 76, beta = 55, grp0 = 'SC_2i', grp1 = 'SC_serum', group = Marker_ori, size_norm = F)

# spike_cds_76_55_count_d <- newCountDataSet(round(exprs(spike_cds_76_55)), (Marker_ori)) # / sizeFactors(spike_cds_76_55)) normalized the data by size factor

# # abs_dtable_pool_max_nbinomTest <- DESeq1_test(abs_count_d, disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomTest') 
# # row.names(abs_dtable_pool_max_nbinomTest$dtalbe) <- abs_dtable_pool_max_nbinomTest$dtalbe$id

# #DESeq glm: (GLM tests are more relevant to our software)
# spike_cds_76_55_dtable_pool_max_nbinomGLMTest <- DESeq1_test(spike_cds_76_55_count_d,  pd = pData(spike_cds_76_55), disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomGLMTest', scale = T) 
# row.names(spike_cds_76_55_dtable_pool_max_nbinomGLMTest$dtalbe) <-  row.names(spike_cds_76_55_count_d)

# # calculate the pval with the readcount with scde: (calculate the scde associate DEG test result LOCALLY) 
# spike_cds_76_55_scde_res_list <- scde_DEG(dir = NULL, count_cds = spike_cds_76_55, DEG_attribute = 'group', contrast = c('SC_2i', 'SC_serum'), n.cores = detectCores(), normalize = T)

# spike_cds_76_55_diff_test_res <- differentialGeneTest(spike_cds_76_55[, ], 
#                                                     fullModelFormulaStr = "~group", 
#                                                     reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)
# # spike_cds_76_55_relative_diff_test_res <- differentialGeneTest(spike_cds_76_55_relative[, ], 
# #                                                              fullModelFormulaStr = "~group", 
# #                                                              reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)
# #calculate the pval with the normalized transcripts with scde: 
# #   abs_scde_res_list <- scde_DEG(dir = NULL, count_cds = spike_cds_76_55, DEG_attribute = 'group', contrast = c('SC_2i', 'SC_serum'), n.cores = detectCores(), normalize = T)
# #   abs_scde_res_list_no_normalize <- scde_DEG(dir = NULL, count_cds = spike_cds_76_55, DEG_attribute = 'Time', contrast = c('SC_2i', 'SC_serum'), n.cores = detectCores())
# #calculate the pval with DESeq 1/2, edgeR: 
# #DESeq2/edgeR: 

# spike_cds_76_55_deseq2_res <- DESeq2_deg(dir = NULL, spike_cds_76_55, Time = Marker_ori, pd = pData(spike_cds_76_55), design = as.formula("~ group"))
# spike_cds_76_55_edgeR_res <- edgeR_test(exprs(spike_cds_76_55), sf = sizeFactors(spike_cds_76_55), group = Marker_ori, glm = T)


# ####################################################################################################################################
# ###make the pre/rec/F1 score and the ROC curves for the UMI dataset: 

# #benchmark the performance: 
# spike_monocle_p <- spike_cds_76_55_diff_test_res[, 'pval'] 
# names(spike_monocle_p) <- row.names(spike_cds_76_55_diff_test_res)

# #deseq
# spike_default_deseq_p <- spike_cds_76_55_dtable_pool_max_nbinomGLMTest$dtalbe[, 'pval'] #abs_dtable_pool_max_nbinomTest
# names(spike_default_deseq_p) <- row.names(spike_cds_76_55_dtable_pool_max_nbinomGLMTest$dtalbe)

# #deseq2
# spike_default_deseq2_p <- results(spike_cds_76_55_deseq2_res)$pvalue #abs_dtable_pool_max_nbinomTest
# names(spike_default_deseq2_p) <- row.names(results(spike_cds_76_55_deseq2_res))

# #edgeR:
# spike_default_edgeR_p <- spike_cds_76_55_edgeR_res$lrt$table$PValue #abs_dtable_pool_max_nbinomTest
# names(spike_default_edgeR_p) <- row.names(spike_cds_76_55_edgeR_res$lrt$table)

# #scde
# spike_scde_p <- spike_cds_76_55_scde_res_list$pval #_no_normalize

df3 <- plot_pre_rec_f1(test_p_list = list(monocle_p = monocle_p,
                                          abs_monocle_p = abs_monocle_p, 
                                          spike_monocle_p = spike_monocle_p, 
                                          abs_default_edgeR_p = spike_default_edgeR_p, 
                                          spike_default_edgeR_p = spike_default_edgeR_p, 
                                          abs_default_deseq_p = abs_default_deseq_p, 
                                          spike_default_deseq_p = spike_default_deseq_p, 
                                          abs_default_deseq2_p = abs_default_deseq2_p, 
                                          spike_default_deseq2_p = spike_default_deseq2_p, 
                                          abs_scde_p = abs_scde_p,
                                          spike_scde_p = spike_scde_p),
                       permutate_pval = list(monocle_p = UMI_cds_76_55_perm_res,
                                          abs_monocle_p = UMI_cds_76_55_perm_res, 
                                          abs_default_edgeR_p = UMI_cds_76_55_perm_res, 
                                          spike_default_edgeR_p = UMI_cds_76_55_perm_res, 
                                          spike_monocle_p = UMI_cds_76_55_perm_res, 
                                          abs_default_deseq_p = UMI_cds_76_55_perm_res, 
                                          spike_default_deseq_p = UMI_cds_76_55_perm_res, 
                                          abs_default_deseq2_p = UMI_cds_76_55_perm_res, 
                                          spike_default_deseq2_p = UMI_cds_76_55_perm_res, 
                                          abs_scde_p = UMI_cds_76_55_perm_res,
                                          spike_scde_p = UMI_cds_76_55_perm_res),
                       row.names(UMI_cds_76_55), #gene_list, overlap_genes, high_gene_list
                       return_df = T, #na.rm = T, 
                       title = 'Comparison of the two-group DEG tests on relative gene expression', 
                       rownames = c('monocle (TPM)', 'monocle (UMI)', 'monocle (spike-in)', 
                        'edgeR (UMI)', 'edgeR (spike-in)',
                        'DESeq (UMI)', 'DESeq (spike-in)', 'DESeq2 (UMI)', 'DESeq2 (spike-in)', 
                        'SCDE (UMI)', 'SCDE (spike-in)'
                        ))

df3$class = '3relative'
colnames(df3)[1:3] <- c('Precision', 'Recall', 'F1 score')
df3$class <- c('spike-in', 'UMI', 'spike-in', 'UMI', 'spike-in', 'UMI', 'spike-in', 'UMI', 'spike-in', 'UMI', 'TPM')
df3[, 'Type'] <- c('SCDE', 'SCDE', 'DESeq2', 'DESeq2', 'DESeq1', 'DESeq1', 'edgeR', 'edgeR', 'monocle', 'monocle', 'monocle') # geom_bar(stat = 'identity', position = 'dodge') 

tmp <- data.frame(Precision = NA, Recall = NA, 'F1 score' = NA,
                  Type = c('SCDE', 'DESeq1', 'DESeq2', 'edgeR'), 
                  class = c('FPKM', 'FPKM', 'FPKM', 'FPKM'))
colnames(tmp)[3] <- 'F1 score'
df_res <- rbind(df3, tmp) 
df_res <-  melt(df_res)

cols <- c("Read counts" = "#00BFC4","UMI" = "#7CAE00", "spike-in" = "#C77CFF", "FPKM" = "#F8766D", 'read_counts' = "#00BFC4", transcript_counts = '#C77CFF', TPM = '#F8766D')
pdf(file = "./supplementary_figures/fig2c_si.pdf", width = 3, height = 2.5)
ggplot(aes(factor(Type), value,  fill = class), data = subset(df_res, Type %in% c('SCDE', 'monocle', 'DESeq2', 'edgeR'))) + geom_bar(position = position_dodge(), stat = 'identity') + #facet_wrap(~variable) + 
  ggtitle(title) +  scale_fill_manual(values = cols) + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') + 
  ggtitle('') + theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 30, hjust = .9)) + theme(strip.background = element_blank()) + nm_theme() + xlab('') + ylim(0, 1)
dev.off()

pdf(file = "./supplementary_figures/fig2c_si_helper.pdf", width = 3, height = 2)
ggplot(aes(factor(Type), value,  fill = class), data = melt(df3)) + geom_bar(position = position_dodge(), stat = 'identity') + #facet_wrap(~variable) + 
  ggtitle(title) +  scale_fill_manual(values = cols) + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') + 
  ggtitle('') + theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 30, hjust = .9)) + theme(strip.background = element_blank())  + xlab('') + ylim(0, 1) #+ nm_theme()
dev.off()
# ###make the roc curves: 


UMI_pval_df <- data.frame(monocle_p = monocle_p,
                                          abs_monocle_p = abs_monocle_p, 
                                          spike_monocle_p = spike_monocle_p, 
                                          abs_default_edgeR_p = abs_default_edgeR_p, 
                                          spike_default_edgeR_p = spike_default_edgeR_p, 
                                          abs_default_deseq_p = abs_default_deseq_p, 
                                          spike_default_deseq_p = spike_default_deseq_p, 
                                          abs_default_deseq2_p = abs_default_deseq2_p, 
                                          spike_default_deseq2_p = spike_default_deseq2_p, 
                                          abs_scde_p = abs_scde_p,
                                          spike_scde_p = spike_scde_p)

row.names(UMI_pval_df) <- names(monocle_p)
permutation_UMI_pval_df <- data.frame(monocle_p = UMI_cds_76_55_perm_res,
                                          abs_monocle_p = UMI_cds_76_55_perm_res, 
                                          abs_default_edgeR_p = UMI_cds_76_55_perm_res, 
                                          spike_default_edgeR_p = UMI_cds_76_55_perm_res, 
                                          spike_monocle_p = UMI_cds_76_55_perm_res, 
                                          abs_default_deseq_p = UMI_cds_76_55_perm_res, 
                                          spike_default_deseq_p = UMI_cds_76_55_perm_res, 
                                          abs_default_deseq2_p = UMI_cds_76_55_perm_res, 
                                          spike_default_deseq2_p = UMI_cds_76_55_perm_res, 
                                          abs_scde_p = UMI_cds_76_55_perm_res,
                                          spike_scde_p = UMI_cds_76_55_perm_res)
row.names(permutation_UMI_pval_df) <- names(monocle_p)

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

select_genes <- row.names(permutation_UMI_pval_df)
p_thrsld <- 0.1

perm_pvals <- permutation_UMI_pval_df[select_genes, 1]
# select_genes <- names[perm_pvals[!is.na(perm_pvals)]]
software_pvals <- UMI_pval_df[select_genes, 1] 

perm_pvals[is.na(perm_pvals)] <- 1
software_pvals[is.na(software_pvals)] <- 1
plot_roc_df(software_pvals, perm_pvals > p_thrsld)
####


select_genes <- names(monocle_p)

UMI_roc_df_list <- lapply(colnames(UMI_pval_df), function(x) {
  print(x)
  
  # if(x %in% c('mode_size_norm_permutate_ratio_by_geometric_mean', 'abs_default_edgeR_p', 'abs_default_deseq2_p', 'abs_scde_p', 'mast_abs_pval_no_norm'))
  #   select_genes <- row.names(new_abs_cds_14_18[1:transcript_num])[esApply(new_abs_cds_14_18[1:transcript_num], 1, function(x) sum(x > 0.5) > 10)]
  # if(x %in% c('monocle_p', 'mast_std_pval'))
  #   select_genes <- row.names(new_std_cds_14_18[1:transcript_num])[esApply(new_std_cds_14_18[1:transcript_num], 1, function(x) sum(x > 0.5) > 10)]
  # if(x %in% c('monocle_p_readcount', 'default_edgeR_p', 'default_deseq2_p', 'default_deseq_p', 'scde_p', 'mast_count_pval_no_norm'))
  #   select_genes <- row.names(count_cds[1:transcript_num])[esApply(count_cds[1:transcript_num], 1, function(x) sum(x > 0.5) > 10)]
  # if(x %in% c('mc_mode_size_norm_permutate_ratio_by_geometric_mean', 'mast_mc_pval_no_norm', 'mc_count_edgeR_p_glm', 'mc_count_deseq2_p', 'mc_count_deseq_p', 'mc_count_scde_p'))
  #   select_genes <- row.names(mc_adj_cds[1:transcript_num])[esApply(mc_adj_cds[1:transcript_num], 1, function(x) sum(x > 0.5) > 10)]

  perm_pvals <- permutation_UMI_pval_df[select_genes, 1] #use the spike-in as the gold standard
  software_pvals <- UMI_pval_df[select_genes, x] 

  perm_pvals[is.na(perm_pvals)] <- 1
  software_pvals[is.na(software_pvals)] <- 1
  res <- generate_roc_df(software_pvals, perm_pvals > p_thrsld)
  colnames(res) <- c('tpr', 'fpr', 'auc')
  cbind(res, method = x)
})

UMI_roc_df_list <- lapply(UMI_roc_df_list, function(x) {colnames(x) <- c('tpr', 'fpr', 'auc', 'method'); x} )
UMI_roc_df <- do.call(rbind, UMI_roc_df_list)
UMI_roc_df[1:5, ]
# str_split_fixed(row.names(UMI_roc_df), '\\.', 2)[, 1]



UMI_auc <- unique(UMI_roc_df[, c('method', 'auc')])
row.names(UMI_auc) <- UMI_auc$method

save(file = './UMI_roc_df', UMI_roc_df)

UMI_roc_df$software <- revalue(UMI_roc_df$method, c("monocle_p" = 'Monocle', "abs_monocle_p" = 'Monocle', "spike_monocle_p" = "Monocle", 
                                            "abs_default_edgeR_p" = 'edgeR', "spike_default_edgeR_p" = 'edgeR',
                                            "abs_default_deseq2_p" = 'DESeq2', "spike_default_deseq2_p" = 'DESeq2',
                                            "abs_default_deseq_p" = 'DESeq', "spike_default_deseq_p" = 'DESeq',
                                            "abs_scde_p" = 'SCDE', "spike_scde_p" = 'SCDE'))

UMI_roc_df$Type <- revalue(UMI_roc_df$method, c("monocle_p" = 'TPM', "abs_monocle_p" = 'UMI', "spike_monocle_p" = "spike-in", 
                                            "abs_default_edgeR_p" = 'UMI', "spike_default_edgeR_p" = 'spike-in',
                                            "abs_default_deseq2_p" = 'UMI', "spike_default_deseq2_p" = 'spike-in',
                                            "abs_default_deseq_p" = 'UMI', "spike_default_deseq_p" = 'spike-in',
                                            "abs_scde_p" = 'UMI', "spike_scde_p" = 'spike-in'))

pdf('./supplementary_figures/UMI_roc.pdf', height = 2.5, width = 3)
qplot(fpr, tpr, data= subset(UMI_roc_df, software %in% c('Monocle', 'edgeR', 'DESeq2', 'SCDE')), geom="step", size = 0.5, color = Type) + #linetype = Type, 
    xlab("False positive rate") +
    ylab("True positive rate") +
      ylim(c(0, 1.0)) + facet_wrap(~software) + scale_size(range = c(0.1, 0.5)) + 
   scale_color_manual(values = cols, name = "Type") + 
  xlim(c(0, 1.0)) + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

pdf('./supplementary_figures/UMI_roc_helper.pdf', height = 13, width = 14)
qplot(fpr, tpr, data= UMI_roc_df, geom="step", color = Type) + #linetype = Type, 
    xlab("False positive rate") +
    ylab("True positive rate") +
      ylim(c(0, 1.0)) + facet_wrap(~software) + 
   scale_color_manual(values = cols, name = "Type") #+ nm_theme()
  xlim(c(0, 1.0)) #+ nm_theme()
dev.off()

#ensure the color is the same: 
#"tpr"      "fpr"      "auc"      "method"   "software" "Type"
tmp <- data.frame(tpr = NA, fpr = NA, auc = NA, method = NA, 
                  software = c('DESeq', 'DESeq2', 'edgeR', 'MAST', 'SCDE'), 
                  Type = c('FPKM', 'FPKM', 'FPKM', 'FPKM', 'FPKM'))
UMI_roc_df <- rbind(UMI_roc_df, tmp) 

#roc_df <-  melt(df_res)
#                              software
# Type                          DESeq DESeq2 edgeR MAST Monocle SCDE
#   estimated transcript counts  1139   1239  1234 2484    1249 1207
#   FPKM                         1229      0     0 1242     996    0
#   read counts                     0   1289  1424    0    1313 1275

pdf('./supplementary_figures/UMI_roc_dfroc_auc_bar.pdf', height = 3, width = 3)
ggplot(aes(software, auc), data = UMI_roc_df) + geom_bar(position = 'dodge', stat = 'identity', aes(fill=Type)) + 
    xlab("") +
    # ylim(0.5, 1.0) + 
    # monocle_theme_opts() +  theme(axis.text.x=element_text(angle=30, hjust=1)) + 
    # scale_fill_manual(values = cols, name = "Software", label = test)  + 
    nm_theme()
dev.off()

pdf('./supplementary_figures/UMI_roc_dfroc_auc_bar_helper.pdf', height = 6, width = 9)
ggplot(aes(software, auc), data = UMI_roc_df) + geom_bar(position = 'dodge', stat = 'identity', aes(fill=Type)) + 
    xlab("") #+
    # ylim(0.5, 1.0) + 
    # monocle_theme_opts() +  theme(axis.text.x=element_text(angle=30, hjust=1)) + 
    # scale_fill_manual(values = cols, name = "Software", label = test)  + 
    # nm_theme()
dev.off()

#show the values of auc 
unique(UMI_roc_df[, c('method','auc')])

# save.image('./RData/analysis_UMI_data.RData')
