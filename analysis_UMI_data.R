#UMI data: 
#from the paper GSE54695_analysis_script.R: 
# library(devtools)
# load_all('~/Projects/monocle-dev')
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

umi_matrix <- read.delim('./data/UMI_data/GSE54695_data_transcript_counts.txt', row.names="GENENAME")

ERCC_ids <- (grep('^ERCC', row.names(umi_matrix)))
input.ERCC.annotation <- input.ERCC.annotation[row.names(umi_matrix)[ERCC_ids], ]
input.ERCC.annotation$numMolecules <- input.ERCC.annotation$conc_attomoles_ul_Mix1*(1000*10^(-3)*1/2500000*10^(-18)*6.02214179*10^(23))
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

all_umi_spike_df <- data.frame(spikein = input.ERCC.annotation[row.names(umi_matrix)[ERCC_ids], 'numMolecules'], UMI = as.vector(exprs(UMI_cds)[ERCC_ids, ]))

UMI_cds <- UMI_cds[, pData(UMI_cds)$Total_mRNAs > 0]
pData(UMI_cds)$dmode <- estimate_t(exprs(UMI_cds))

UMI_cds_subset <-  UMI_cds[, pData(UMI_cds)$Total_mRNAs > 10000 & pData(UMI_cds)$group %in% c('SC_2i', 'SC_serum')]
UMI_cds_subset <- estimateSizeFactors(UMI_cds_subset)
UMI_cds_subset@expressionFamily <- tobit() #fix the bug of estimateDispersions
UMI_ordering_genes <- selectGenesInExpressionRange(UMI_cds_subset, 5, Inf, 0.1, stat_fun=function(x) { median(round(x)) })
UMI_cds_subset <- setOrderingFilter(UMI_cds_subset, UMI_ordering_genes)
UMI_cds_subset <- reduceDimension(UMI_cds_subset, use_irlba = F, use_vst = F, method = "ICA", scaling = F) 
UMI_cds_subset <- orderCells(UMI_cds_subset, num_paths = 2, reverse = F) #SRR1033962_0 
UMI_cds_subset@expressionFamily <- negbinomial()

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

cols <- c("Read counts" = "#00BFC4","MC transcripts" = "#7CAE00", "Spikein transcripts" = "#C77CFF", "FPKM" = "#F8766D", 'read_counts' = "#00BFC4", transcript_counts = '#C77CFF', UMI = '#C77CFF', TPM = '#F8766D')
pdf(file = "./supplementary_figures/fig2c_si.pdf", width = 2.5, height = 2)
ggplot(aes(factor(Type), value,  fill = class), data = melt(df3)) + geom_bar(position = position_dodge(), stat = 'identity') + #facet_wrap(~variable) + 
  ggtitle(title) +  scale_fill_manual(values = cols) + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') + 
  ggtitle('') + theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 30, hjust = .9)) + theme(strip.background = element_blank()) + nm_theme() + xlab('') + ylim(0, 1)
dev.off()

save.image('./RData/analysis_UMI_data.RData')
