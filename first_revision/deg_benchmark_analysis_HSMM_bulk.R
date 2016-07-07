edgeR_test <- function (counts = exprs(count_cds), sf = sizeFactors(count_cds), group = Time_ori, glm = F, normalize = T) {
    print(paste('size factors are ', sf))
    y <- DGEList(counts = round(counts), group = group)
    if (!glm) {
        # y <- calcNormFactors(y)
        if(normalize){
          y$sample$norm.factors <- sf
        }

        y <- estimateCommonDisp(y)
        y <- estimateTagwiseDisp(y)
        et <- exactTest(y)
        return(list(et = et, topTags = topTags(et)))
    }
    else {
        # y <- calcNormFactors(y)
        if(normalize){
          y$sample$norm.factors <- sf
        }
        design <- model.matrix(~group)
        y <- estimateGLMCommonDisp(y, design)
        y <- estimateGLMTrendedDisp(y, design)
        y <- estimateGLMTagwiseDisp(y, design)
        fit <- glmFit(y, design)
        lrt <- glmLRT(fit, coef = 2)
        topTags(lrt)
        return(list(lrt = lrt, topTags = topTags(lrt)))
    }
}

#update DESeq 1: 
DESeq1_test <- function (count_d, pd = pData(new_abs_cds_14_18[, colnames(count_d)]), condition = Time_ori, disp_method = "blind",
    sharing = "fit-only", test_type = "nbinomGLMTest", scale = F,
    contrast = c("E14.5", "E18.5"), fullModelFormulaStr = "count~condition",
    reducedModelFormulaStr = "count~1"){
    message(class(count_d))
    if (class(count_d)[1] != "CountDataSet")
        count_d <- newCountDataSet(round(exprs(cds)), condition)
    d <- estimateSizeFactors(count_d)
    if (!scale)
        pData(d)$sizeFactor <- 1
    else
        pData(d)$sizeFactor <- pd[colnames(count_d), 'Size_Factor']

    d2 <- d
    sizeFactor <- pData(d)$sizeFactor
    print(sizeFactor)
    d <- DESeq::estimateDispersions(d2, method = disp_method,
        sharingMode = sharing)
    message("pass here")
    message("estimate dispersion")
    tmp <- do.call(cbind.data.frame, lapply(d@fitInfo, function(x) x[c(1,
        3)]))
    disp_func <- lapply(ls(d@fitInfo), function(x) fitInfo(d,
        name = x))
    message("all disp functions extracted")
    if (test_type == "nbinomGLMTest") {
        message("test is ", test_type)
        print(DESeq::sizeFactors(d))
        dfit1 <- fitNbinomGLMs(d, as.formula(fullModelFormulaStr))
        dfit0 <- fitNbinomGLMs(d, as.formula(reducedModelFormulaStr))
        dpval <- nbinomGLMTest(dfit1, dfit0)
        dpadj <- p.adjust(dpval, method = "BH")
        dtable <- transform(dfit1, pval = dpval, padj = dpadj)
        dtable <- cbind(dtable, tmp)
        head(dtable)
        message("pass glm test")
        return(list(dfit1 = dfit1, dfit0 = dfit0, full_df.residual = attr(dfit1,
            "df.residual"), reduced_df.residual = attr(dfit0,
            "df.residual"), full_deviance = dfit1$deviance, reduced_deviance = dfit0$deviance,
            dtalbe = dtable, disp_func = disp_func))
    }
    else if (test_type == "nbinomTest") {
        dtable <- nbinomTest(d, contrast[1], contrast[2])
        dtable <- cbind(dtable, tmp)
        message("pass nb test")
        return(list(dtalbe = dtable, disp_func = disp_func))
    }
}

#update DESeq 2: 
DESeq2_deg <- function (dir, cds, Time, test_type = "LRT", design = as.formula("~ Time"),
    reduced = as.formula("~ 1"), pd = pData(cds)) {
    pd$sizeFactor <- pd$Size_Factor
    if (!is.null(dir)) {
        sample_table <- read.delim(paste(dir, "/samples.table",
            sep = ""))
        norm_count <- read.delim(paste(dir, "/genes.count_table",
            sep = ""))
        row.names(norm_count) <- norm_count$tracking_id
        norm_count <- norm_count[, -1]
        countdata <- round(t(t(norm_count) * sample_table$internal_scale))
        if (exists("sample_sheet"))
            countdata <- countdata[, row.names(sample_sheet)]
    }
    else countdata <- round(exprs(cds))
    DEseq2_cnt_cds <- DESeqDataSetFromMatrix(countData = countdata[,
        ], colData = pd, design = design)
    DEseq2_cnt_cds <- DESeq(DEseq2_cnt_cds, test = test_type,
        reduced = reduced)
    return(DEseq2_cnt_cds)
}

MAST_deg <- function(cds, grp = 'Time', test.type = 'hurdle', pseudo_cnt = 1, normalization = F) {
  data <- as.matrix(exprs(cds))
  
  if(normalization)
    data <- log2(t(t(data) / sizeFactors(cds)) + pseudo_cnt)
  else
    data <- log2(data + pseudo_cnt)
  
  data <- melt(data)
  colnames(data) <- c('Gene', 'Cell', 'exprs')
  data[, grp] <- c(pData(cds)[data[, 2], grp])
  data$ncells <- 1 
  
  mast_cds <- FluidigmAssay(data, idvars="Cell", 
                            primerid='Gene', measurement='exprs', geneid="Gene", 
                            ncells = 'ncells', phenovars=grp)
  
  zlm.output <- zlm.SingleCellAssay(as.formula(paste("~ ", grp)), mast_cds, method='glm', ebayes=TRUE) #
  # show(zlm.output)
  zlm.lr <- lrTest(zlm.output, grp)
  # dimnames(zlm.lr)
  
  pval_df <- zlm.lr[,,'Pr(>Chisq)']
  pval <- pval_df[, test.type]
  
  return(pval)
}

# library(monocle)
library(devtools)
load_all('~/Projects/monocle-dev')
library(xacHelper) 
# library(MAST)
library(ROCR)
script.dir <- dirname(sys.frame(1)$ofile)

source(paste(script.dir, '/roc_curves.R', sep = ''))
#load the data: 
load('./RData/analysis_HSMM_data.RData')
HSMM_bulk <- read.table("./data/HSMM_data/bulk_cuffdiff/gene_exp.diff", header = T, row.names = NULL)

HSMM_bulk_T0_T72 <-  subset(HSMM_bulk, sample_1 == "T0" & sample_2 == "T24")
order_stats_HSMM_bulk_T0_T72 <-  HSMM_bulk_T0_T72[order(abs(HSMM_bulk_T0_T72$test_stat), decreasing=T), ]

#add read counts: 
HSMM_readcounts_matrix <- read.delim("./data/HSMM_data/muscle/HSMM/HSMM_cuffnorm_out/genes.count_table")
row.names(HSMM_readcounts_matrix) <- HSMM_readcounts_matrix$tracking_id
HSMM_readcounts_matrix <- HSMM_readcounts_matrix[,-1]

pd <- new("AnnotatedDataFrame", data = pData(HSMM_myo))
fd <- new("AnnotatedDataFrame", data = fData(HSMM_myo))
HSMM_readcounts_matrix_cds <-  newCellDataSet(as.matrix(HSMM_fpkm_matrix[row.names(HSMM_myo), colnames(HSMM_myo)]), 
                                   phenoData = pd, 
                                   featureData = fd, 
                                   expressionFamily=negbinomial(), 
                                   lowerDetectionLimit=1)

##integrate into the benchmark analysis: 

##########################based on the deg_benchmark_analysis.R#############################

library(scde)
#########################run the following scripts in remote sever#############################

HSMM_valid_cells <- colnames(HSMM_myo)
#generate the pvals from the statistical test (permutation based or from software: monocle/DESeq/SCDE)

#double check the cells selected are the same as we used in the paper:  identical(colnames(new_mc_cds_T0_T72), cells_used)
new_mc_cds_T0_T72 <- HSMM_myo[, intersect(HSMM_valid_cells, row.names(subset(pData(HSMM_myo), Time %in% c('T72_CT_', 'T0_CT_'))))]
new_mc_cds_T0_T72 <- estimateSizeFactors(new_mc_cds_T0_T72) #calculate the size factor for performing the relative absolute expression tests
new_mc_cds_T0_T72 <- estimateDispersions(new_mc_cds_T0_T72) #calculate the size factor for performing the relative absolute expression tests

new_std_cds_T0_T72 <- std_HSMM[, colnames(new_mc_cds_T0_T72)]
# new_std_cds_T0_T72 <- estimateSizeFactors(new_std_cds_T0_T72) #FPKM values are already on relative scale

#create a cds for readcount data to perform the default DEG tests for DESeq and SCDE : 
count_cds_HSMM_bulk <- HSMM_readcounts_matrix_cds[row.names(new_mc_cds_T0_T72), colnames(new_mc_cds_T0_T72)] 

count_cds_HSMM_bulk <- estimateSizeFactors(count_cds_HSMM_bulk)
count_cds_HSMM_bulk <- estimateDispersions(count_cds_HSMM_bulk)

# # calculate the pval with the readcount with scde: (calculate the scde associate DEG test result LOCALLY) 
std_scde_res_list <- scde_DEG(dir = NULL, count_cds = count_cds_HSMM_bulk, DEG_attribute = 'Time', contrast = c('T0_CT_', 'T72_CT_'), n.cores = detectCores(), normalize = T)
# std_scde_res_list_no_normalize <- scde_DEG(dir = NULL, count_cds = count_cds_HSMM_bulk, DEG_attribute = 'Time', contrast = c('T0_CT_', 'T72_CT_'), n.cores = detectCores())

# #calculate the pval with the normalized transcripts with scde: 
# mc_scde_res_list <- scde_DEG(dir = NULL, count_cds = new_mc_cds_T0_T72, DEG_attribute = 'Time', contrast = c('T0_CT_', 'T72_CT_'), n.cores = detectCores(), normalize = T)
# mc_scde_res_list_no_normalize <- scde_DEG(dir = NULL, count_cds = new_mc_cds_T0_T72, DEG_attribute = 'Time', contrast = c('T0_CT_', 'T72_CT_'), n.cores = detectCores())

#load all other necessary packages: 
load_all_libraries()

#perform  the stastical tests on the data: 
new_std_diff_test_res <- differentialGeneTest(new_std_cds_T0_T72[, ], 
                                              fullModelFormulaStr = "~Time", 
                                              reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)  

new_size_norm_mc_diff_test_res <- differentialGeneTest(new_mc_cds_T0_T72[, ], 
                                                       fullModelFormulaStr = "~Time", 
                                                       reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

new_cds_diff_test_res <- differentialGeneTest(count_cds_HSMM_bulk[, ], 
                                              fullModelFormulaStr = "~Time", 
                                              reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

new_cds_diff_test_res_no_relative <- differentialGeneTest(count_cds_HSMM_bulk[, ], 
                                                          fullModelFormulaStr = "~Time", 
                                                          reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)

Time_ori <- pData(new_mc_cds_T0_T72)$Time #Time as input for newCountDataSet
Time_ori <- as.character(Time_ori) #Time as input for newCountDataSet

calculate the pval with the readcount with DESeq: 
read_count_d <- newCountDataSet(round(exprs(count_cds_HSMM_bulk)), Time_ori) 
std_dtable_pool_max_nbinomTest <- DESeq1_test(read_count_d, pd = pData(count_cds_HSMM_bulk[, colnames(read_count_d)]), disp_method = 'pooled', contrast = c("T0_CT_", "T72_CT_"), sharing = 'maximum', test_type = 'nbinomTest', scale = T) 
row.names(std_dtable_pool_max_nbinomTest$dtalbe) <- std_dtable_pool_max_nbinomTest$dtalbe$id

# # calculate the pval with the normalized transcripts with DESeq: 
mc_count_d <- newCountDataSet(round(exprs(new_mc_cds_T0_T72)),  Time_ori) #normalized the data by size factor
mc_dtable_pool_max_nbinomTest <- DESeq1_test(mc_count_d, pd = pData(new_mc_cds_T0_T72[, colnames(mc_count_d)]), disp_method = 'pooled', contrast = c("T0_CT_", "T72_CT_"), sharing = 'maximum', test_type = 'nbinomTest', scale = T) 
row.names(mc_dtable_pool_max_nbinomTest$dtalbe) <- mc_dtable_pool_max_nbinomTest$dtalbe$id

# DESeq glm: (GLM tests are more relevant to our software)
std_dtable_pool_max_nbinomGLMTest <- DESeq1_test(read_count_d[, ], pd = pData(count_cds_HSMM_bulk[, colnames(read_count_d)]), disp_method = 'pooled', contrast = c("T0_CT_", "T72_CT_"), sharing = 'maximum', test_type = 'nbinomGLMTest', scale = T) 
row.names(std_dtable_pool_max_nbinomGLMTest$dtalbe) <-  row.names(read_count_d[, ])
mc_dtable_pool_max_nbinomGLMTest <- DESeq1_test(mc_count_d[, ], pd = pData(new_mc_cds_T0_T72[, colnames(mc_count_d)]), disp_method = 'pooled', contrast = c("T0_CT_", "T72_CT_"), sharing = 'maximum', test_type = 'nbinomGLMTest', scale = T) 
row.names(mc_dtable_pool_max_nbinomGLMTest$dtalbe) <- row.names(new_mc_cds_T0_T72[, ])

#add the DEG tests using edgeR / DESeq2: 
edgeR_res <- edgeR_test(counts = exprs(count_cds_HSMM_bulk), group = Time_ori, sf = sizeFactors(count_cds_HSMM_bulk), glm = F)
mc_edgeR_res <- edgeR_test(exprs(new_mc_cds_T0_T72), group = Time_ori, sf = sizeFactors(new_mc_cds_T0_T72), glm = F)

edgeR_res_glm <- edgeR_test(counts = exprs(count_cds_HSMM_bulk), group = Time_ori, sf = sizeFactors(count_cds_HSMM_bulk), glm = T)
mc_edgeR_res_glm <- edgeR_test(exprs(new_mc_cds_T0_T72), group = Time_ori, sf = sizeFactors(new_mc_cds_T0_T72), glm = T)

deseq2_res <- DESeq2_deg(dir = NULL, count_cds_HSMM_bulk[, ], Time = Time_ori, pd = pData(count_cds_HSMM_bulk))
mc_deseq2_res <- DESeq2_deg(dir = NULL, new_mc_cds_T0_T72[, ], Time = Time_ori, pd = pData(new_mc_cds_T0_T72))

prepare to generate the data for create the precision/recall/F1 score data.frame: 
#new monocle_p: 
monocle_p_HSMM_bulk <- new_std_diff_test_res[, 'pval'] 
names(monocle_p_HSMM_bulk) <- row.names(new_std_diff_test_res)
#use cds: 
monocle_p_readcount_HSMM_bulk <- new_cds_diff_test_res[, 'pval'] 
names(monocle_p_readcount_HSMM_bulk) <- row.names(new_cds_diff_test_res)

#mc
new_mc_size_norm_monocle_p_ratio_by_geometric_mean_HSMM_bulk <- new_size_norm_mc_diff_test_res[, 'pval']
names(new_mc_size_norm_monocle_p_ratio_by_geometric_mean_HSMM_bulk) <- row.names(new_size_norm_mc_diff_test_res)
#mc 
new_mc_size_norm_monocle_p_ratio_by_geometric_mean_HSMM_bulk <- new_size_norm_mc_diff_test_res[, 'pval']
names(new_mc_size_norm_monocle_p_ratio_by_geometric_mean_HSMM_bulk) <- row.names(new_size_norm_mc_diff_test_res)
#deseq
default_deseq_p_HSMM_bulk <- std_dtable_pool_max_nbinomGLMTest$dtalbe[, 'pval'] #std_dtable_pool_max_nbinomTest
names(default_deseq_p_HSMM_bulk) <- row.names(std_dtable_pool_max_nbinomGLMTest$dtalbe) #std_dtable_pool_max_nbinomGLMTest
mc_default_deseq_p_HSMM_bulk <- mc_dtable_pool_max_nbinomGLMTest$dtalbe[, 'pval'] #abs_dtable_pool_max_nbinomTest
names(mc_default_deseq_p_HSMM_bulk) <- row.names(mc_dtable_pool_max_nbinomGLMTest$dtalbe)
# scde
#save.image('tmp_benchmark_analysis.RData')

scde_p_HSMM_bulk <- std_scde_res_list$pval 
mc_scde_p_HSMM_bulk <- mc_scde_res_list$pval #_no_normalize

# generate df3.1 for making the comparision: 
#select only high expressed genes: 
# high_gene_list <- esApply(HSMM_myo[, colnames(AT12_cds_subset_all_gene)], 1, function(x) sum(x > 1) > 15)

#prepare the data for the DESeq2 / edgeR: 
default_edgeR_p_HSMM_bulk = edgeR_res_glm$lr$table$PValue
names(default_edgeR_p_HSMM_bulk) <- row.names(edgeR_res_glm$lrt$table)
mc_default_edgeR_p = mc_edgeR_res_glm$lrt$table$PValue
names(mc_default_edgeR_p) <- row.names(mc_edgeR_res_glm$lrt$table)

default_deseq2_p_HSMM_bulk = results(deseq2_res)$pvalue
names(default_deseq2_p_HSMM_bulk) <- row.names(results(deseq2_res))
mc_default_deseq2_p_HSMM_bulk = results(mc_deseq2_res)$pvalue 
names(mc_default_deseq2_p_HSMM_bulk) <- row.names(results(mc_deseq2_res))

#mast function: 
# mast_mc_pval_no_norm_HSMM_bulk <- MAST_deg(new_mc_cds_T0_T72)
# mast_std_pval_no_norm_HSMM_bulk <- MAST_deg(new_std_cds_T0_T72)
# mast_count_pval_no_norm_HSMM_bulk <- MAST_deg(count_cds_HSMM_bulk)

# mast_mc_pval_norm <- MAST_deg(new_mc_cds_T0_T72, normalization = T)
# # mast_std_pval_norm <- MAST_deg(new_std_cds_T0_T72, normalization = T)
# mast_count_pval_norm <- MAST_deg(count_cds_HSMM_bulk, normalization = T)

# save(file = './RData/MAST_res', mast_mc_pval_no_norm_HSMM_bulk, mast_std_pval_no_norm_HSMM_bulk, mast_count_pval_no_norm_HSMM_bulk, mast_mc_pval_norm, mast_count_pval_norm)

load('./RData/MAST_res')
top_1k_HSMM_bulk_T0_T72_pval <- HSMM_bulk_T0_T72$p_value
names(top_1k_HSMM_bulk_T0_T72_pval) <- HSMM_bulk_T0_T72$gene_id

# 1) signficant by bulk

# 2) AND have a abs(log2(fold_change)) > 2 in bulk

# 3) AND are detectable in at least 15 single cells

# or however many cells with expression >= 1  you used for the ROC plots in the lung data

select_genes <- row.names(new_std_cds_T0_T72)[esApply(new_std_cds_T0_T72, 1, function(x) sum(x > 0.5) > 10)]
# select_genes <- as.character(subset(order_stats_HSMM_bulk_T0_T72, value_1 > 10 & value_2 > 10)[, 'gene_id'])

top_1k_genes <- order_stats_HSMM_bulk_T0_T72[which(order_stats_HSMM_bulk_T0_T72[order_stats_HSMM_bulk_T0_T72$significant == 'no', 'gene_id'] %in% select_genes), ][1:1000, 'gene_id']

reverse_order_stats_HSMM_bulk_T0_T72 <-  HSMM_bulk_T0_T72[order(abs(HSMM_bulk_T0_T72$q_value), decreasing=T), ]
bottom_1k_genes <- reverse_order_stats_HSMM_bulk_T0_T72[which(reverse_order_stats_HSMM_bulk_T0_T72[reverse_order_stats_HSMM_bulk_T0_T72$significant == 'no', 'gene_id'] %in% select_genes), ][1:1000, 'gene_id']

# bottom_1k_genes <- setdiff(select_genes, bottom_1k_genes)
# #top/bottom 1k genes: 
select_genes <- c(as.character(top_1k_genes), as.character(bottom_1k_genes))
true_data <- rep(0, length(select_genes))
names(true_data) <- select_genes
true_data[1:length(top_1k_genes)] <- 1
true_data[(length(top_1k_genes) + 1):(length(select_genes))] <- 0

true_df <- data.frame(Type = true_data, pval = top_1k_HSMM_bulk_T0_T72_pval[select_genes])

#new version: 
# bulk_sig_genes <- as.character(subset(order_stats_HSMM_bulk_T0_T72, significant == 'yes' & abs(log2.fold_change.) > 1)[, 'gene_id'])

# true_data <- rep(0, nrow(new_std_cds_T0_T72))
# names(true_data) <- row.names(new_std_cds_T0_T72)
# true_data[bulk_sig_genes] <- 1

# true_df <- data.frame(Type = true_data[select_genes], pval = top_1k_HSMM_bulk_T0_T72_pval[select_genes])

#generate the dataframe for making the benchmarking plots: 
############NOTE: mode_size_norm_permutate_ratio_by_geometric_mean may be changed into the same as DESEQ size normalization############ 
df3_HSMM_bulk <- plot_pre_rec_f1(test_p_list = list(monocle_p = monocle_p_HSMM_bulk, 
                                        monocle_p_readcount = monocle_p_readcount_HSMM_bulk,
                                        mc_mode_size_norm_permutate_ratio_by_geometric_mean = new_mc_size_norm_monocle_p_ratio_by_geometric_mean_HSMM_bulk, 
                                        default_edgeR_p = default_edgeR_p_HSMM_bulk, 
                                        mc_default_edgeR_p = mc_default_edgeR_p,         
                                        default_deseq2_p = default_deseq2_p_HSMM_bulk, 
                                        mc_default_deseq2_p = mc_default_deseq2_p_HSMM_bulk, 
                                        default_deseq_p = default_deseq_p_HSMM_bulk, 
                                        mc_default_deseq_p = mc_default_deseq_p_HSMM_bulk, 
                                        scde_p = scde_p_HSMM_bulk, 
                                        mc_scde_p = mc_scde_p_HSMM_bulk, 
                                        # mast_abs_pval_no_norm = mast_abs_pval_no_norm_HSMM_bulk, 
                                        mast_mc_pval_no_norm = mast_mc_pval_no_norm_HSMM_bulk, 
                                        mast_std_pval_no_norm = mast_std_pval_no_norm_HSMM_bulk, 
                                        mast_count_pval_no_norm = mast_count_pval_no_norm_HSMM_bulk),
                     permutate_pval = list(monocle_p = top_1k_HSMM_bulk_T0_T72_pval, #readcount_permutate_pval, #std_permutate_pval, 
                                           monocle_p_readcount = top_1k_HSMM_bulk_T0_T72_pval, 
                                           mc_mode_size_norm_permutate_ratio_by_geometric_mean = top_1k_HSMM_bulk_T0_T72_pval,
                                           default_edgeR_p = top_1k_HSMM_bulk_T0_T72_pval, #readcount_permutate_pval use readcount instead of std_permutate_pval
                                           mc_default_edgeR_p = top_1k_HSMM_bulk_T0_T72_pval, 
                                           default_deseq2_p = top_1k_HSMM_bulk_T0_T72_pval, #readcount_permutate_pval use readcount instead of std_permutate_pval
                                           mc_default_deseq2_p = top_1k_HSMM_bulk_T0_T72_pval, 
                                           default_deseq_p = top_1k_HSMM_bulk_T0_T72_pval, #readcount_permutate_pval use readcount instead of std_permutate_pval
                                           mc_default_deseq_p = top_1k_HSMM_bulk_T0_T72_pval, 
                                           # abs_default_deseq_p_new_norm = mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           scde_p = top_1k_HSMM_bulk_T0_T72_pval, 
                                           mc_scde_p = top_1k_HSMM_bulk_T0_T72_pval, 
                                           # mast_abs_pval_no_norm = top_1k_HSMM_bulk_T0_T72_pval, 
                                           mast_mc_pval_no_norm = top_1k_HSMM_bulk_T0_T72_pval, 
                                           mast_std_pval_no_norm = top_1k_HSMM_bulk_T0_T72_pval, 
                                           mast_count_pval_no_norm = top_1k_HSMM_bulk_T0_T72_pval),
                     names(true_data), #gene_list, overlap_genes, high_gene_list
                     return_df = T, #na.rm = T, 
                     p_thrsld = 0.05, #0.05
                     rownames = c('monocle (FPKM)', 'monocle (readcount)', 'monocle (New size normalization, Estimate transcript)', 
                        'edgeR (edgeR size normalization)', 'edgeR (New Size normalization)', 'DESeq2 (DESeq2 size normalization)', 'DESeq2 (New Size normalization)',
                        'DESeq (DESeq size normalization)', "DESeq (New Size normalization)", 'SCDE (Read Counts)', 'SCDE (New size normalization)', 
                        'MAST (mc no normalization)', 'MAST (FPKM no normalization)', 'MAST (readcount no normalization)'))
df3_HSMM_bulk$data_type = c("Read counts", "FPKM", "MC transcripts", "MC transcripts", "Read counts", "MC transcripts", "Read counts", 

"MC transcripts", "Read counts", "MC transcripts", "Read counts", 

"MC transcripts", "Read counts", "FPKM")

df3_HSMM_bulk$class = '3relative'


# only show new size normalization: 
df3.1_HSMM_bulk <- df3_HSMM_bulk
df3.1_HSMM_bulk[, 'Type'] <- c('MAST', 'MAST', 'MAST', 'SCDE', 'SCDE', 'DESeq1', 'DESeq1', 'DESeq2', 'DESeq2', 'edgeR', 'edgeR', 'Monocle', 'Monocle', 'Monocle') # geom_bar(stat = 'identity', position = 'dodge') 

tmp <- data.frame(Type = c('SCDE', 'DESeq1', 'DESeq2', 'edgeR'), 
                  data_type = c('FPKM', 'FPKM', 'FPKM', 'FPKM'),
                  class = '3relative', 
                  pre = NA, rec = NA, f1 = NA)
df_res <- rbind(df3.1_HSMM_bulk, tmp) 
df_res <-  melt(df_res)

pdf('./supplementary_figures/fig2a_si_HSMM_bulk.pdf', width = 3, height = 3)
qplot(factor(Type), value, stat = "identity", geom = 'bar', position = 'dodge', fill = data_type, data = subset(df_res, Type %in% c('Monocle', 'edgeR', 'DESeq2', 'SCDE'))) + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + monocle_theme_opts() + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) + nm_theme()
dev.off()

pdf('./supplementary_figures/fig2a_si_HSMM_bulk_helper.pdf', width = 6, height = 2)
qplot(factor(Type), value, stat = "identity", geom = 'bar', position = 'dodge', fill = data_type, data = df_res) + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + monocle_theme_opts() + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) 
dev.off()

HSMM_bulk_pval_df <- data.frame(monocle_p = monocle_p_HSMM_bulk, 
                                        monocle_p_readcount = monocle_p_readcount_HSMM_bulk,
                                        mc_mode_size_norm_permutate_ratio_by_geometric_mean = new_mc_size_norm_monocle_p_ratio_by_geometric_mean_HSMM_bulk, 
                                        default_edgeR_p = default_edgeR_p_HSMM_bulk[names(monocle_p_HSMM_bulk)], 
                                        mc_default_edgeR_p = mc_default_edgeR_p[names(monocle_p_HSMM_bulk)],         
                                        default_deseq2_p = default_deseq2_p_HSMM_bulk[names(monocle_p_HSMM_bulk)], 
                                        mc_default_deseq2_p = mc_default_deseq2_p_HSMM_bulk[names(monocle_p_HSMM_bulk)], 
                                        default_deseq_p = default_deseq_p_HSMM_bulk[names(monocle_p_HSMM_bulk)], 
                                        mc_default_deseq_p = mc_default_deseq_p_HSMM_bulk[names(monocle_p_HSMM_bulk)], 
                                        scde_p = scde_p_HSMM_bulk[names(monocle_p_HSMM_bulk)], 
                                        mc_scde_p = mc_scde_p_HSMM_bulk[names(monocle_p_HSMM_bulk)], 
                                        # mast_abs_pval_no_norm = mast_abs_pval_no_norm_HSMM_bulk, 
                                        mast_mc_pval_no_norm = mast_mc_pval_no_norm_HSMM_bulk, 
                                        mast_std_pval_no_norm = mast_std_pval_no_norm_HSMM_bulk, 
                                        mast_count_pval_no_norm = mast_count_pval_no_norm_HSMM_bulk)

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

roc_df_list <- apply(HSMM_bulk_pval_df[select_genes, ], 2, function(x) generate_roc_df(x, true_data[select_genes] != 1))

roc_df_list <- lapply(roc_df_list, function(x) {colnames(x) <- c('tpr', 'fpr', 'auc'); x} )
roc_df <- do.call(rbind, roc_df_list)
roc_df$method <-  str_split_fixed(row.names(roc_df), '\\.', 2)[, 1]

auc <- unique(roc_df[, c('method', 'auc')])
row.names(auc) <- auc$method
hmcols <- blue2green2red(nrow(auc))

test <- c(paste(sort(colnames(HSMM_bulk_pval_df)), ':', " AUC = ", auc[order(colnames(HSMM_bulk_pval_df)), 2]))
save(file = './roc_df', roc_df)
roc_df$software <- revalue(roc_df$method, c("monocle_p" = 'Monocle', "monocle_p_readcount" = 'Monocle', "mc_mode_size_norm_permutate_ratio_by_geometric_mean" = 'Monocle', 
                                            "default_edgeR_p" = 'edgeR', "mc_default_edgeR_p" = 'edgeR',
                                            "default_deseq2_p" = 'DESeq2', "mc_default_deseq2_p" = 'DESeq2',
                                            "default_deseq_p" = 'DESeq', "mc_default_deseq_p" = 'DESeq',
                                            "scde_p" = 'SCDE', "mc_scde_p" = 'SCDE',
                                            "mast_mc_pval_no_norm" = 'MAST', 
                                            "mast_std_pval_no_norm" = "MAST", 
                                            "mast_count_pval_no_norm" = "MAST"))

roc_df$Type <- revalue(roc_df$method, c("monocle_p" = 'FPKM', "monocle_p_readcount" = 'read counts', "mc_mode_size_norm_permutate_ratio_by_geometric_mean" = 'estimated transcript counts', 
                                            "default_edgeR_p" = 'read counts', "mc_default_edgeR_p" = 'estimated transcript counts',
                                            "default_deseq2_p" = 'read counts', "mc_default_deseq2_p" = 'estimated transcript counts',
                                            "default_deseq_p" = 'read counts', "mc_default_deseq_p" = 'estimated transcript counts',
                                            "scde_p" = 'read counts', "mc_scde_p" = 'estimated transcript counts',
                                            "mast_mc_pval_no_norm" = 'estimated transcript counts', 
                                            "mast_std_pval_no_norm" = "FPKM", 
                                            "mast_count_pval_no_norm" = "read counts"))

cols <- c("FPKM" = "#F2756D", "read.counts" = "#6F94CC", "transcript.counts" = "#26B24B", 'read counts' = "#6F94CC", "estimated transcript counts" = "#7BAE41")
pdf('./supplementary_figures/roc.pdf', height = 2.5, width = 3)
qplot(fpr, tpr, data= subset(roc_df, software %in% c('Monocle', 'edgeR', 'DESeq2', 'SCDE')), geom="step", color = Type) + #, linetype = Type
    xlab("False positive rate") +
    ylab("True positive rate") +
      ylim(c(0, 1.0)) + facet_wrap(~software) + 
  xlim(c(0, 1.0)) + monocle_theme_opts() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
   scale_color_manual(values = cols, name = "Type") + nm_theme()
dev.off()

pdf('./supplementary_figures/roc_helper.pdf', height = 13, width = 14)
qplot(fpr, tpr, data=roc_df, geom="step", linetype = Type, color = Type) + 
    xlab("False positive rate") +
    ylab("True positive rate") +
      ylim(c(0, 1.0)) + facet_wrap(~software) + 
  xlim(c(0, 1.0)) + monocle_theme_opts() +
   scale_color_manual(values = cols, name = "Type") #+ nm_theme()
dev.off()

# #ensure the color is the same: 
# #"tpr"      "fpr"      "auc"      "method"   "software" "Type"
cols <- c("FPKM" = "#F2756D", "read.counts" = "#6F94CC", "transcript.counts" = "#26B24B")
roc_df$Type <- revalue(roc_df$Type, c("estimated transcript counts" =  "transcript.counts", "FPKM" = "FPKM", "read counts" = "read.counts"))
tmp <- data.frame(tpr = NA, fpr = NA, auc = NA, method = NA, 
                  software = c('DESeq', 'DESeq2', 'edgeR', 'MAST', 'SCDE'), 
                  Type = c('FPKM', 'FPKM', 'FPKM', 'FPKM', 'FPKM'))
roc_df <- rbind(roc_df, tmp) 
#roc_df <-  melt(df_res)
#                              software
# Type                          DESeq DESeq2 edgeR MAST Monocle SCDE
#   estimated transcript counts  1139   1239  1234 2484    1249 1207
#   FPKM                         1229      0     0 1242     996    0
#   read counts                     0   1289  1424    0    1313 1275

pdf('./supplementary_figures/roc_auc_bar.pdf', height = 3, width = 3)
ggplot(aes(software, auc), data = subset(roc_df, software %in% c('Monocle', 'edgeR', 'DESeq2', 'SCDE'))) + geom_bar(position = 'dodge', stat = 'identity', aes(fill=Type)) + 
    xlab("") +
    # ylim(0.5, 1.0) + 
    # monocle_theme_opts() +  theme(axis.text.x=element_text(angle=30, hjust=1)) + 
    scale_fill_manual(values = cols, name = "Software", label = test)  + nm_theme()
dev.off()

pdf('./supplementary_figures/roc_auc_bar.pdf', height = 3, width = 3)
qplot(method, auc, stat = 'identity', data=roc_df, fill=method, geom="bar") + 
    xlab("") +
    ylim(c(0, 1.0)) + 
    monocle_theme_opts() +  theme(axis.text.x=element_text(angle=30, hjust=1)) + 
    scale_fill_manual(values = cols, name = "Software", label = test)  + nm_theme()
dev.off()

pdf('./supplementary_figures/roc_helper.pdf', height = 3, width = 14)
qplot(fpr, tpr, data=roc_df, color=method, geom="step") + 
    xlab("False positive rate") +
    ylab("True positive rate") +
      ylim(c(0, 1.0)) + 
  xlim(c(0, 1.0)) + monocle_theme_opts() +
   scale_color_manual(values = hmcols, name = "Test", label = test)  # + nm_theme()
dev.off()

#show the values of auc 
unique(roc_df[, c('method','auc')])
#                                                             auc
# monocle_p.1                                           0.9365970
# monocle_p_readcount.1                                 0.8637075
# mc_mode_size_norm_permutate_ratio_by_geometric_mean.1 0.8878245
# default_edgeR_p.1                                     0.8411780
# mc_default_edgeR_p.1                                  0.9112060
# default_deseq2_p.1                                    0.9222525
# mc_default_deseq2_p.1                                 0.9217450
# default_deseq_p.1                                     0.8520790
# mc_default_deseq_p.1                                  0.9043890
# scde_p.1                                              0.9296265
# mc_scde_p.1                                           0.9438740
# mast_mc_pval_no_norm.1                                0.8229900
# mast_std_pval_no_norm.1                               0.8265460
# mast_count_pval_no_norm.1                             0.8265460

save.image('./RData/deg_benchmark_analysis_HSMM_bulk.RData')


