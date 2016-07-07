######################################################
#update the edgeRTest: 
load_all_libraries <- function ()
{
    # library(monocle)
    # library(xacHelper)
    library(lmtest)
    library(MASS)
    library(DESeq)
    library(reshape2)
    library(DESeq2)
    library(fitdistrplus)
    library(VGAM)
    library(pscl)
    library(zoo)
    library(VennDiagram)
    library(stringr)
    library(limma)
    library(reshape2)
    library(plyr)
    library(modeest)
    library(MASS)
    library(mixsmsn)
    library(doMC)
    library(data.table) 
    library(boot)
    library(testthat)
    require(pheatmap)
    library(gplots)
    library(RColorBrewer)
    library(piano)
    library(ggdendro)
    library(stringr)
    library(grid)
    library(reshape2)
    library(snow)
    library(princurve)
    library(matrixStats)
    library(colorRamps)
    library(edgeR)
    library(plyr)
    library(pheatmap)
    library(stringr)
    library(piano)
    library(grid)
    library(colorRamps)
    library(RColorBrewer)
    library(R.utils)
    library(sp)
    library(raster)
    library(venneuler)
    library(colorRamps)
    library(scde)
    library(stringr)
    library(plyr)
}
edgeR_test <- function (counts = exprs(count_cds), sf = sizeFactors(count_cds), group = Time_ori, glm = F, normalize = T) {
    print(paste('size factors are ', sf))
    y <- DGEList(counts = round(counts), group = group)
    if (!glm) {
        
        if(normalize){
          y$sample$norm.factors <- sf
        }
        else
            y <- calcNormFactors(y)

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
    contrast = conditions, fullModelFormulaStr = "count~condition",
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

######################################################

load('./RData/prepare_lung_data.RData')
# library(devtools)
# load_all('~/Projects/monocle-dev')
library(monocle)

library(xacHelper)

library(scde)

#this function enables to select cells from particular time groups: 
if(!exists("conditions"))
    conditions <- c('E14.5', 'E18.5')

#######################run the following scripts in remote sever#############################

#generate the pvals from the statistical test (permutation based or from software: monocle/DESeq/SCDE)

#double check the cells selected are the same as we used in the paper:  identical(colnames(new_abs_cds_14_18), cells_used)
new_abs_cds_14_18 <- absolute_cds[, intersect(valid_cells, row.names(subset(pData(absolute_cds), Time %in% conditions)))]
new_abs_cds_14_18 <- estimateSizeFactors(new_abs_cds_14_18) #calculate the size factor for performing the relative absolute expression tests
 
table(pData(new_abs_cds_14_18)$Time) #check the number of used cells at each Time Point

new_mc_cds_14_18 <- mc_adj_cds[, colnames(new_abs_cds_14_18)]
new_mc_cds_14_18 <- estimateSizeFactors(new_mc_cds_14_18)

new_std_cds_14_18 <- standard_cds[, colnames(new_abs_cds_14_18)]
# new_std_cds_14_18 <- estimateSizeFactors(new_std_cds_14_18) #FPKM values are already on relative scale

#prepare the readcount data for DESeq / SCDE: 
read_countdata <- read_countdata[row.names(new_abs_cds_14_18), colnames(new_abs_cds_14_18)] 

#create a cds for readcount data to perform the default DEG tests for DESeq and SCDE : 
count_cds <- newCellDataSet(read_countdata[row.names(new_abs_cds_14_18), colnames(new_abs_cds_14_18)],
                            phenoData = new("AnnotatedDataFrame", data = pData(new_abs_cds_14_18)),
                            featureData = new("AnnotatedDataFrame", data = fData(new_abs_cds_14_18)),
                            expressionFamily = negbinomial(),
                            lowerDetectionLimit = 1)

count_cds <- estimateSizeFactors(count_cds)

#calculate the pval with the readcount with scde: (calculate the scde associate DEG test result LOCALLY, this part is super slow, parallel doesn't work) 
std_scde_res_list <- scde_DEG(dir = NULL, count_cds = count_cds, DEG_attribute = 'Time', contrast = conditions, n.cores = 1)
save(std_scde_res_list, file = paste(conditions[1], conditions[2], 'std_scde_res_list', sep = ''))

std_scde_res_list_no_normalize <- scde_DEG(dir = NULL, count_cds = count_cds, DEG_attribute = 'Time', contrast = conditions, n.cores = 1)
save(std_scde_res_list_no_normalize, file = paste(conditions[1], conditions[2], 'std_scde_res_list_no_normalize', sep = ''))

#calculate the pval with the normalized transcripts with scde: 
abs_scde_res_list <- scde_DEG(dir = NULL, count_cds = new_abs_cds_14_18, DEG_attribute = 'Time', contrast = conditions, n.cores = 1, normalize = T)
save(abs_scde_res_list, file = paste(conditions[1], conditions[2], 'abs_scde_res_list', sep = ''))

abs_scde_res_list_no_normalize <- scde_DEG(dir = NULL, count_cds = new_abs_cds_14_18, DEG_attribute = 'Time', contrast = conditions, n.cores = 1)
save(abs_scde_res_list_no_normalize, file = paste(conditions[1], conditions[2], 'abs_scde_res_list_no_normalize', sep = ''))

mc_scde_res_list_no_normalize <- scde_DEG(dir = NULL, count_cds = new_mc_cds_14_18, DEG_attribute = 'Time', contrast = conditions, n.cores = 1, normalize = F)
save(mc_scde_res_list_no_normalize, file = paste(conditions[1], conditions[2], 'mc_scde_res_list_no_normalize', sep = ''))

mc_scde_res_list <- scde_DEG(dir = NULL, count_cds = new_mc_cds_14_18, DEG_attribute = 'Time', contrast = conditions, n.cores = 1, normalize = T)
save(mc_scde_res_list, file = paste(conditions[1], conditions[2], 'mc_scde_res_list', sep = ''))

load(paste(conditions[1], conditions[2], 'std_scde_res_list', sep = ''))
load(paste(conditions[1], conditions[2], 'std_scde_res_list_no_normalize', sep = ''))
load(paste(conditions[1], conditions[2], 'abs_scde_res_list', sep = ''))
load(paste(conditions[1], conditions[2], 'abs_scde_res_list_no_normalize', sep = ''))
load(paste(conditions[1], conditions[2], 'mc_scde_res_list_no_normalize', sep = ''))
load(paste(conditions[1], conditions[2], 'mc_scde_res_list', sep = ''))

# load all other necessary packages: 
load_all_libraries()

# #perform  the stastical tests on the data: 
new_std_diff_test_res <- differentialGeneTest(new_std_cds_14_18[1:transcript_num, ], 
                                              fullModelFormulaStr = "~Time", 
                                              reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)  

new_std_diff_test_res_relative <- differentialGeneTest(new_std_cds_14_18[1:transcript_num, ], 
                                              fullModelFormulaStr = "~Time", 
                                              reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)  

new_size_norm_abs_diff_test_res <- differentialGeneTest(new_abs_cds_14_18[1:transcript_num, ], 
                                                        fullModelFormulaStr = "~Time", 
                                                        reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

new_size_norm_mc_diff_test_res <- differentialGeneTest(new_mc_cds_14_18[1:transcript_num, ], 
                                                       fullModelFormulaStr = "~Time", 
                                                       reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

new_cds_diff_test_res <- differentialGeneTest(count_cds[1:transcript_num, ], 
                                              fullModelFormulaStr = "~Time", 
                                              reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

new_cds_diff_test_res_no_relative <- differentialGeneTest(count_cds[1:transcript_num, ], 
                                                          fullModelFormulaStr = "~Time", 
                                                          reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)

Time_ori <- pData(new_abs_cds_14_18)$Time #Time as input for newCountDataSet

#calculate the pval with the readcount with DESeq: 
readcount_sf_mat <- matrix(rep(sizeFactors(count_cds), nrow(count_cds)), nrow = nrow(count_cds), byrow = T, dimnames = dimnames(count_cds))
abs_sf_mat <- matrix(rep(sizeFactors(new_abs_cds_14_18), nrow(new_abs_cds_14_18)), nrow = nrow(new_abs_cds_14_18), byrow = T)
mc_sf_mat <- matrix(rep(sizeFactors(new_mc_cds_14_18), nrow(new_mc_cds_14_18)), nrow = nrow(new_mc_cds_14_18), byrow = T)

read_count_d <- newCountDataSet(round(read_countdata[, ]), Time_ori) #/ readcount_sf_mat[, colnames(read_countdata)]
std_dtable_pool_max_nbinomTest <- DESeq1_test(read_count_d, pd = pData(count_cds), disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomTest', scale = T) 
row.names(std_dtable_pool_max_nbinomTest$dtalbe) <- std_dtable_pool_max_nbinomTest$dtalbe$id

#calculate the pval with the normalized transcripts with DESeq: 
abs_count_d <- newCountDataSet(round(t(t(exprs(new_abs_cds_14_18)) )), (Time_ori)) #Don't normalized the data by size factor
abs_dtable_pool_max_nbinomTest <- DESeq1_test(abs_count_d, pd = pData(new_abs_cds_14_18), disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomTest', scale = T) 
row.names(abs_dtable_pool_max_nbinomTest$dtalbe) <- abs_dtable_pool_max_nbinomTest$dtalbe$id

#DESeq glm: (GLM tests are more relevant to our software)
std_dtable_pool_max_nbinomGLMTest <- DESeq1_test(read_count_d[1:transcript_num, ], pd = pData(count_cds), disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomGLMTest', scale = T) 
row.names(std_dtable_pool_max_nbinomGLMTest$dtalbe) <-  row.names(read_count_d[1:transcript_num, ])
abs_dtable_pool_max_nbinomGLMTest <- DESeq1_test(abs_count_d[1:transcript_num, ], pd = pData(new_abs_cds_14_18), disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomGLMTest', scale = T) 
row.names(abs_dtable_pool_max_nbinomGLMTest$dtalbe) <- row.names(new_abs_cds_14_18[1:transcript_num, ])

#permutation results for two group test on the FPKM value (gold standard for monocle on FPKM values with tobit model): 
Time_order <- order(Time_ori) #order the E14.5 cell at the begining (ensure E14.5 cells number is 43 while E18.5 is 74)
std_split_cds <- split(t(exprs(new_std_cds_14_18[1:transcript_num, Time_order])), col(t(exprs(new_std_cds_14_18[1:transcript_num, Time_order])), as.factor = T))
std_fc <- esApply(new_std_cds_14_18[1:transcript_num, ], 1, mean_fc, grp0 = conditions[1], grp1 = conditions[2])
std_split_fc <- split(t(std_fc), col(t(std_fc), as.factor = T))

permuation_pval <- function (x, fc, alpha = table(pData(new_abs_cds_14_18)$Time)[conditions[1]], beta = table(pData(new_abs_cds_14_18)$Time)[conditions[2]], permutate_num = 10000,
    return_fc = F)
{
    fc_vec <- rep(0, permutate_num)
    if (is.na(fc))
        return(1)
    for (i in 1:permutate_num) {
        x <- sample(x, length(x))
        x_0 <- (mean(x[1:alpha]))
        x_1 <- (mean(x[(alpha + 1):(alpha + beta)]))
        mean_fc <- log2(x_1/x_0)
        fc_vec[i] <- mean_fc
    }
    if (fc > 0)
        pval <- sum(fc <= fc_vec)/length(fc_vec)
    else pval <- sum(fc >= fc_vec)/length(fc_vec)
    if (return_fc)
        return(list(pval = pval, fc_vec = fc_vec))
    else return(pval)
}
closeAllConnections()
save.image('./RData/deg_benchmark_analysis_tmp.RData')
std_permutate_pval <- mcmapply(permuation_pval, std_split_cds, std_split_fc, mc.cores = detectCores()) #multiple cores 

closeAllConnections()
# permutation results for two group test on the readcounts value (gold standard for DESEeq / SCDE on readcount data): 
# no size factor normalization: 
readcount_split_cds <- split(t(exprs(count_cds[1:transcript_num, Time_order])), col(t(exprs(count_cds[1:transcript_num, Time_order])), as.factor = T))
readcount_fc <- esApply(count_cds[1:transcript_num, ], 1, mean_fc, grp0 = conditions[1], grp1 = conditions[2]) #valid_gene_id

readcount_split_cds <- split(t(round(exprs(count_cds[1:transcript_num, Time_order]) / readcount_sf_mat[1:transcript_num, Time_order])), col(t(exprs(count_cds[1:transcript_num, Time_order])), as.factor = T))
readcount_fc <- apply(round(exprs(count_cds)[1:transcript_num, Time_order] / readcount_sf_mat[1:transcript_num, Time_order]), 1, mean_fc, grp0 = conditions[1], grp1 = conditions[2], grp = pData(count_cds[1:transcript_num, Time_order])$Time) #valid_gene_id
readcount_split_fc <- split(t(readcount_fc), col(t(readcount_fc), as.factor = T))
readcount_permutate_pval <- mcmapply(permuation_pval, readcount_split_cds, readcount_split_fc, mc.cores = detectCores()) #multiple core
closeAllConnections()

# no size factortwo-group permutation tests (the same as above)
# permutation results for two group test on the normalized absolute value (gold standard for DESEeq / SCDE on readcount data): 
# readcount_permutate_pval <- permu_two_group_gen(count_cds[, Time_order], alpha = 43, beta = 60, grp0 = conditions[1], grp1 = conditions[2], group = pData(count_cds)$Time)
#mode_size_norm_permutate_ratio_by_geometric_mean <- cal_perm_pval_size_norm(new_abs_cds_14_18[1:transcript_num, Time_order])
#mc_mode_size_norm_permutate_ratio_by_geometric_mean <- cal_perm_pval_size_norm(new_mc_cds_14_18[1:transcript_num, Time_order])

# new implementation with size factor normalization: 
abs_split_cds <- split(round(t(exprs(new_abs_cds_14_18[1:transcript_num, Time_order]) / abs_sf_mat[1:transcript_num, Time_order])), col(t(exprs(new_abs_cds_14_18[1:transcript_num, Time_order])), as.factor = T))
abs_fc <- apply(round(exprs(new_abs_cds_14_18)[1:transcript_num, Time_order] / abs_sf_mat[1:transcript_num, Time_order]), 1, mean_fc, grp0 = conditions[1], grp1 = conditions[2], grp = pData(new_abs_cds_14_18[1:transcript_num, Time_order])$Time) #valid_gene_id
abs_split_fc <- split(t(abs_fc), col(t(abs_fc), as.factor = T))
mode_size_norm_permutate_ratio_by_geometric_mean <- mcmapply(permuation_pval, abs_split_cds, abs_split_fc, mc.cores = detectCores()) #multiple core
closeAllConnections()

mc_split_cds <- split(round(t(exprs(new_mc_cds_14_18[1:transcript_num, Time_order]) / mc_sf_mat[1:transcript_num, Time_order])), col(t(exprs(new_mc_cds_14_18[1:transcript_num, Time_order])), as.factor = T))
mc_fc <- apply(round(exprs(new_mc_cds_14_18)[1:transcript_num, Time_order] / mc_sf_mat[1:transcript_num, Time_order]), 1, mean_fc, grp0 = conditions[1], grp1 = conditions[2], grp = pData(new_mc_cds_14_18[1:transcript_num, Time_order])$Time) #valid_gene_id
mc_split_fc <- split(t(mc_fc), col(t(mc_fc), as.factor = T))
mc_mode_size_norm_permutate_ratio_by_geometric_mean <- mcmapply(permuation_pval, mc_split_cds, mc_split_fc, mc.cores = detectCores()) #multiple core
closeAllConnections()

# # # add the DEG tests using edgeR / DESeq2: 
edgeR_res <- edgeR_test(counts = exprs(count_cds)[1:transcript_num, ], sf = sizeFactors(count_cds), glm = F)
abs_edgeR_res <- edgeR_test(exprs(new_abs_cds_14_18[1:transcript_num, ]), sf = sizeFactors(new_abs_cds_14_18), group = Time_ori, glm = F)

edgeR_res_glm <- edgeR_test(counts = exprs(count_cds)[1:transcript_num, ], sf = sizeFactors(count_cds), glm = T)
abs_edgeR_res_glm <- edgeR_test(exprs(new_abs_cds_14_18[1:transcript_num, ]), sf = sizeFactors(new_abs_cds_14_18), group = Time_ori, glm = T)
  
deseq2_res <- DESeq2_deg(dir = NULL, count_cds[1:transcript_num, ], Time = Time_ori, pd = pData(count_cds))
abs_deseq2_res <- DESeq2_deg(dir = NULL, new_abs_cds_14_18[1:transcript_num, ], Time = Time_ori, pd = pData(new_abs_cds_14_18))

#add the analysis for the mc counts: 
new_mc_cds_14_18_d <- newCountDataSet(round(exprs(new_mc_cds_14_18))[, ], Time_ori) 
mc_dtable_pool_max_nbinomGLMTest <- DESeq1_test(new_mc_cds_14_18_d[1:transcript_num, ], pd = pData(new_mc_cds_14_18), disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomGLMTest', scale = T) 
row.names(mc_dtable_pool_max_nbinomGLMTest$dtalbe) <-  row.names(new_mc_cds_14_18[1:transcript_num, ])

mc_dtable_pool_max_nbinomTest <- DESeq1_test(new_mc_cds_14_18_d, pd = pData(new_mc_cds_14_18), disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomTest', scale = T) 
row.names(mc_dtable_pool_max_nbinomTest$dtalbe) <- mc_dtable_pool_max_nbinomTest$dtalbe$id

mc_edgeR_res_glm <- edgeR_test(exprs(new_mc_cds_14_18[1:transcript_num, ]), sf = sizeFactors(new_mc_cds_14_18), group = Time_ori, glm = T)
mc_edgeR_res <- edgeR_test(exprs(new_mc_cds_14_18[1:transcript_num, ]), sf = sizeFactors(new_mc_cds_14_18), group = Time_ori, glm = F)
mc_deseq2_res <- DESeq2_deg(dir = NULL, new_mc_cds_14_18[1:transcript_num, ], Time = Time_ori, pd = pData(new_mc_cds_14_18))

# prepare to generate the data for create the precision/recall/F1 score data.frame: 
#new monocle_p: 
monocle_p <- new_std_diff_test_res[, 'pval'] 
names(monocle_p) <- row.names(new_std_diff_test_res)
#use cds: 
monocle_p_readcount <- new_cds_diff_test_res[, 'pval'] 
names(monocle_p_readcount) <- row.names(new_cds_diff_test_res)

#abs
new_abs_size_norm_monocle_p_ratio_by_geometric_mean <- new_size_norm_abs_diff_test_res[, 'pval']
names(new_abs_size_norm_monocle_p_ratio_by_geometric_mean) <- row.names(new_size_norm_abs_diff_test_res)
#mc 
new_mc_size_norm_monocle_p_ratio_by_geometric_mean <- new_size_norm_mc_diff_test_res[, 'pval']
names(new_mc_size_norm_monocle_p_ratio_by_geometric_mean) <- row.names(new_size_norm_mc_diff_test_res)

#deseq
default_deseq_p <- std_dtable_pool_max_nbinomGLMTest$dtalbe[, 'pval'] #std_dtable_pool_max_nbinomTest
names(default_deseq_p) <- row.names(std_dtable_pool_max_nbinomGLMTest$dtalbe) #std_dtable_pool_max_nbinomGLMTest
abs_default_deseq_p <- abs_dtable_pool_max_nbinomGLMTest$dtalbe[, 'pval'] #abs_dtable_pool_max_nbinomTest
names(abs_default_deseq_p) <- row.names(abs_dtable_pool_max_nbinomGLMTest$dtalbe)
mc_default_deseq_p <- mc_dtable_pool_max_nbinomGLMTest$dtalbe[, 'pval'] #abs_dtable_pool_max_nbinomTest
names(mc_default_deseq_p) <- row.names(mc_dtable_pool_max_nbinomGLMTest$dtalbe)

#deseq2
default_deseq2_p = results(deseq2_res)$pvalue
names(default_deseq2_p) <- row.names(results(deseq2_res))
abs_default_deseq2_p = results(abs_deseq2_res)$pvalue 
names(abs_default_deseq2_p) <- row.names(results(abs_deseq2_res))
mc_default_deseq2_p = results(mc_deseq2_res)$pvalue 
names(mc_default_deseq2_p) <- row.names(results(mc_deseq2_res))

# scde
# save.image('tmp_benchmark_analysis.RData')
scde_p <- std_scde_res_list$pval[names(monocle_p)] 
abs_scde_p <- abs_scde_res_list$pval
mc_scde_p <- mc_scde_res_list$pval

#prepare the data for the DESeq2 / edgeR: 
default_edgeR_p_glm = edgeR_res_glm$lrt$table$PValue
names(default_edgeR_p_glm) <- row.names(edgeR_res_glm$lrt$table)
abs_default_edgeR_p_glm = abs_edgeR_res_glm$lrt$table$PValue
names(abs_default_edgeR_p_glm) <- row.names(abs_edgeR_res_glm$lrt$table)
mc_edgeR_p_glm <- mc_edgeR_res_glm$lrt$table$PValue
names(mc_edgeR_p_glm) <- row.names(mc_edgeR_res_glm$lrt$table)

#non_glm: 
default_edgeR_p = edgeR_res$et$table$PValue
names(default_edgeR_p) <- row.names(edgeR_res$et$table)
abs_default_edgeR_p = abs_edgeR_res$et$table$PValue
names(abs_default_edgeR_p) <- row.names(abs_edgeR_res$et$table)
mc_default_edgeR_p = mc_edgeR_res$et$table$PValue
names(mc_default_edgeR_p) <- row.names(mc_edgeR_res$et$table)

# generate df3.1 for making the comparision: 

#select only high expressed genes: 
high_gene_list <- esApply(absolute_cds[, colnames(AT12_cds_subset_all_gene)], 1, function(x) sum(x > 1) > 15)
#overlap genes: 
permutation_sets <- list(FPKM = names(std_permutate_pval[which(std_permutate_pval < .01)]), 
                         # TPM = names(TPM_permutate_pval[TPM_permutate_pval < 0.01]),
                         # no_sf_TPM_permutate_pval = names(no_sf_TPM_permutate_pval[no_sf_TPM_permutate_pval < 0.01]),
                         transcript = names(mode_size_norm_permutate_ratio_by_geometric_mean[mode_size_norm_permutate_ratio_by_geometric_mean < 0.01]),
                         estimate_transcript = names(mc_mode_size_norm_permutate_ratio_by_geometric_mean[mc_mode_size_norm_permutate_ratio_by_geometric_mean < 0.01]),
                         readcount = names(readcount_permutate_pval[readcount_permutate_pval < 0.01]))
overlap_genes <- Reduce(intersect,  permutation_sets) 

TP_permutation_sets <- list(FPKM = length(std_permutate_pval[which(std_permutate_pval < .01)]), 
                         # TPM = length(TPM_permutate_pval[TPM_permutate_pval < 0.01]),
                         # no_sf_TPM_permutate_pval = length(no_sf_TPM_permutate_pval[no_sf_TPM_permutate_pval < 0.01]),
                         transcript = length(mode_size_norm_permutate_ratio_by_geometric_mean[mode_size_norm_permutate_ratio_by_geometric_mean < 0.01]),
                         estimate_transcript = length(mc_mode_size_norm_permutate_ratio_by_geometric_mean[mc_mode_size_norm_permutate_ratio_by_geometric_mean < 0.01]),
                         readcount = length(readcount_permutate_pval[readcount_permutate_pval < 0.01]))

#select genes for benchmarking the performance: (this should matach with the roc_auc plot) 
select_genes <- row.names(new_std_cds_14_18[1:transcript_num])[esApply(new_std_cds_14_18[1:transcript_num], 1, function(x) sum(x > 0.5) > 10)]

#generate the dataframe for making the benchmarking plots: 
df3 <- plot_pre_rec_f1(test_p_list = list(monocle_p = monocle_p, 
                                        monocle_p_readcount = monocle_p_readcount,
                                        mode_size_norm_permutate_ratio_by_geometric_mean = new_abs_size_norm_monocle_p_ratio_by_geometric_mean,
                                        mc_mode_size_norm_permutate_ratio_by_geometric_mean = new_mc_size_norm_monocle_p_ratio_by_geometric_mean, 
                                        default_edgeR_p = default_edgeR_p_glm, 
                                        abs_default_edgeR_p = abs_default_edgeR_p_glm,         
                                        default_deseq2_p = default_deseq2_p, 
                                        abs_default_deseq2_p = abs_default_deseq2_p, 
                                        default_deseq_p = default_deseq_p, 
                                        abs_default_deseq_p = abs_default_deseq_p, 
                                        scde_p = scde_p, 
                                        abs_scde_p = abs_scde_p,

                                        # tpm_count_monocle_p = tpm_count_monocle_p, 
                                        # tpm_count_edgeR_p_glm = tpm_count_edgeR_p_glm, 
                                        # tpm_count_deseq2_p = tpm_count_deseq2_p, 
                                        # tpm_count_deseq_p = tpm_count_deseq_p, 
                                        # tpm_count_scde_p = tpm_count_scde_p, 

                                        # mc_count_monocle_p = new_mc_size_norm_monocle_p_ratio_by_geometric_mean, 
                                        mc_count_edgeR_p_glm = mc_edgeR_p_glm, 
                                        mc_count_deseq2_p = mc_default_deseq2_p, 
                                        mc_count_deseq_p = mc_default_deseq_p, 
                                        mc_count_scde_p = mc_scde_p
                                        ),
                     permutate_pval = list(monocle_p = mode_size_norm_permutate_ratio_by_geometric_mean,#std_permutate_pval, #readcount_permutate_pval, #std_permutate_pval, 
                                           monocle_p_readcount = mode_size_norm_permutate_ratio_by_geometric_mean,#readcount_permutate_pval, 
                                           mode_size_norm_permutate_ratio_by_geometric_mean = mode_size_norm_permutate_ratio_by_geometric_mean,
                                           mc_mode_size_norm_permutate_ratio_by_geometric_mean = mode_size_norm_permutate_ratio_by_geometric_mean,#mc_mode_size_norm_permutate_ratio_by_geometric_mean,
                                           default_edgeR_p = mode_size_norm_permutate_ratio_by_geometric_mean,#readcount_permutate_pval, #readcount_permutate_pval use readcount instead of std_permutate_pval
                                           abs_default_edgeR_p = mode_size_norm_permutate_ratio_by_geometric_mean,#mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           default_deseq2_p = mode_size_norm_permutate_ratio_by_geometric_mean,#readcount_permutate_pval, #readcount_permutate_pval use readcount instead of std_permutate_pval
                                           abs_default_deseq2_p = mode_size_norm_permutate_ratio_by_geometric_mean,#mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           default_deseq_p = mode_size_norm_permutate_ratio_by_geometric_mean,#readcount_permutate_pval, #readcount_permutate_pval use readcount instead of std_permutate_pval
                                           abs_default_deseq_p = mode_size_norm_permutate_ratio_by_geometric_mean,#mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           # abs_default_deseq_p_new_norm = mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           scde_p = mode_size_norm_permutate_ratio_by_geometric_mean,#readcount_permutate_pval, 
                                           abs_scde_p = mode_size_norm_permutate_ratio_by_geometric_mean,#mode_size_norm_permutate_ratio_by_geometric_mean,

                                           # tpm_count_monocle_p = TPM_permutate_pval, #readcount_permutate_pval, #std_permutate_pval, 
                                           # tpm_count_edgeR_p_glm = TPM_permutate_pval, 
                                           # tpm_count_deseq2_p = TPM_permutate_pval,
                                           # tpm_count_deseq_p = TPM_permutate_pval,
                                           # tpm_count_scde_p = TPM_permutate_pval, 

                                           mc_count_edgeR_p_glm = mode_size_norm_permutate_ratio_by_geometric_mean,#mc_mode_size_norm_permutate_ratio_by_geometric_mean, 
                                           mc_count_deseq2_p = mode_size_norm_permutate_ratio_by_geometric_mean,#mc_mode_size_norm_permutate_ratio_by_geometric_mean,
                                           mc_count_deseq_p = mode_size_norm_permutate_ratio_by_geometric_mean,#mc_mode_size_norm_permutate_ratio_by_geometric_mean,
                                           mc_count_scde_p = mode_size_norm_permutate_ratio_by_geometric_mean#mc_mode_size_norm_permutate_ratio_by_geometric_mean
                                           ),
                     select_genes, #gene_list, overlap_genes, high_gene_list
                     return_df = T, #na.rm = T, 
                     p_thrsld = 0.05 #0.05
                     # rownames = c('monocle (FPKM)', 'monocle (readcount)', 'monocle (New size normalization)', 'monocle (New size normalization, Estimate transcript)', 
                     #    'edgeR (edgeR size normalization)', 'edgeR (New Size normalization)', 'DESeq2 (DESeq2 size normalization)', 'DESeq2 (New Size normalization)',
                     #    'DESeq (DESeq size normalization)', "DESeq (New Size normalization)", 'SCDE (Read Counts)', 'SCDE (New size normalization)', 

                     #     'monocle (TPM counts)', 'edgeR (TPM counts)', 'DESeq2 (TPM counts)','DESeq (TPM counts)', 'SCDE (TPM counts)'
                        # )
                     )
df3$data_type = c("MC transcripts", "MC transcripts", "MC transcripts", "MC transcripts", 
    # "TPM counts", "TPM counts", "TPM counts", "TPM counts", "TPM counts", 

"Spikein transcripts", "Read counts", "Spikein transcripts", "Read counts", 

"Spikein transcripts", "Read counts", "Spikein transcripts", "Read counts", 

"MC transcripts", "Spikein transcripts", "Read counts", "FPKM")

df3$class = '3relative'


# only show new size normalization: 
df3.1 <- df3
df3.1[, 'Type'] <- c('SCDE', 'DESeq1', 'DESeq2', 'edgeR', 
 # 'SCDE', 'DESeq1', 'DESeq2', 'edgeR', 'Monocle', 
'SCDE', 'SCDE', 'DESeq1', 'DESeq1', 'DESeq2', 'DESeq2', 'edgeR', 'edgeR', 'Monocle', 'Monocle', 'Monocle', 'Monocle') # geom_bar(stat = 'identity', position = 'dodge') 

pdf('./supplementary_figures/fig2a_si.pdf', width = 3, height = 2)
qplot(factor(Type), value, stat = "identity", geom = 'bar', position = 'dodge', fill = data_type, data = melt(df3.1)) + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + monocle_theme_opts() + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) + nm_theme()
dev.off()

pdf('./supplementary_figures/fig2a_si_helper.pdf', width = 10, height = 4)
qplot(factor(Type), value, stat = "identity", geom = 'bar', position = 'dodge', fill = data_type, data = melt(df3.1)) + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + monocle_theme_opts() + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) #+ nm_theme()
dev.off()

# distribution of pval vs the mean expression: 
pval_mean <- data.frame(std = std_permutate_pval, abs = mode_size_norm_permutate_ratio_by_geometric_mean, 
                        mc = mc_mode_size_norm_permutate_ratio_by_geometric_mean, read = readcount_permutate_pval, 
                        fpkm_mean = apply(new_std_cds_14_18[1:transcript_num, ], 1, mean), 
                        spike_mean = apply(new_abs_cds_14_18[1:transcript_num, ], 1, mean), 
                        mc_mean = apply(new_mc_cds_14_18[1:transcript_num, ], 1, mean), 
                        read_mean = apply(count_cds[1:transcript_num, ], 1, mean) )


# # fig 3h: 
# # show only the spike-in / mc algorithm test: 
mc_spikein_df <- plot_pre_rec_f1(test_p_list = list(mode_size_norm_permutate_ratio_by_geometric_mean = new_abs_size_norm_monocle_p_ratio_by_geometric_mean,
                                    mc_mode_size_norm_permutate_ratio_by_geometric_mean = new_mc_size_norm_monocle_p_ratio_by_geometric_mean),
                 permutate_pval = list(mode_size_norm_permutate_ratio_by_geometric_mean = mode_size_norm_permutate_ratio_by_geometric_mean,
                                       mc_mode_size_norm_permutate_ratio_by_geometric_mean = mc_mode_size_norm_permutate_ratio_by_geometric_mean),
                 row.names(absolute_cds), #gene_list, overlap_genes, high_gene_list
                 return_df = T, #na.rm = T, 
                 p_thrsld = 0.01, #0.05
                 rownames = c('monocle (New size normalization)', 'monocle (New size normalization, Estimate transcript)'))
mc_spikein_df$data_type = c("Spikein transcripts", "estimated transcripts")

mc_spikein_df[, 'Type'] <- c('Monocle', 'Monocle') # geom_bar(stat = 'identity', position = 'dodge') 
colnames(mc_spikein_df)[1:3] <- c('Precision', 'Recall', 'F1 score')

pdf('./main_figures/fig3h.pdf', width = 1.7, height = 1.9)
ggplot(aes(factor(Type), value,  fill = data_type), data = melt(mc_spikein_df)) + geom_bar(position = position_dodge(), stat = 'identity') + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) + nm_theme() + xlab('') + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

pdf('./tmp/fig3h_helper.pdf', width = 3, height = 2)
ggplot(aes(factor(Type), value,  fill = data_type), data = melt(mc_spikein_df)) + geom_bar(position = position_dodge(), stat = 'identity') + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank())
dev.off()


# save.image('./RData/deg_benchmark_analysis.RData')
save.image(paste('./RData/deg_benchmark_analysis', conditions[1], conditions[2], '.RData', sep = ''))

###############################################################################################################################################################################################
#benchmark with TPM times median total mRNAs: 

#create a cds with TPM times the a median expression level: 
lung_TPM_mat_14_18 <- exprs(TPM_cds[row.names(new_abs_cds_14_18), colnames(new_abs_cds_14_18)]) * median(pData(new_abs_cds_14_18)$Total_mRNAs) / 10^6

lung_TPM_count_cds_14_18 <- newCellDataSet(lung_TPM_mat_14_18,
                            phenoData = new("AnnotatedDataFrame", data = pData(new_abs_cds_14_18)),
                            featureData = new("AnnotatedDataFrame", data = fData(new_abs_cds_14_18)),
                            expressionFamily = negbinomial(),
                            lowerDetectionLimit = 1)
lung_TPM_count_cds_14_18 <- estimateSizeFactors(lung_TPM_count_cds_14_18)
lung_TPM_count_cds_14_18 <- estimateDispersions(lung_TPM_count_cds_14_18)

lung_TPM_count_scde_res_list <- scde_DEG(dir = NULL, count_cds = lung_TPM_count_cds_14_18, DEG_attribute = 'Time', contrast = conditions, n.cores = 1)
lung_TPM_count_diff_test_res <- differentialGeneTest(lung_TPM_count_cds_14_18[1:transcript_num, ], 
                                                        fullModelFormulaStr = "~Time", 
                                                        reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)

Time_order <- order(Time_ori) #order the E14.5 cell at the begining (ensure E14.5 cells number is 43 while E18.5 is 74)

lung_TPM_count_d <- newCountDataSet(round(lung_TPM_mat_14_18)[, ], Time_ori) 
lung_TPM_count_dtable_pool_max_nbinomGLMTest <- DESeq1_test(lung_TPM_count_d[1:transcript_num, ], pd = pData(lung_TPM_mat_14_18), disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomGLMTest', scale = T) 
row.names(lung_TPM_count_dtable_pool_max_nbinomGLMTest$dtalbe) <-  row.names(lung_TPM_count_d[1:transcript_num, ])

permutation results for two group test on the FPKM value (gold standard for monocle on FPKM values with tobit model): 
TPM_split_cds <- split(t(exprs(lung_TPM_count_cds_14_18[1:transcript_num, Time_order])), col(t(exprs(lung_TPM_count_cds_14_18[1:transcript_num, Time_order])), as.factor = T))
TPM_fc <- esApply(lung_TPM_count_cds_14_18[1:transcript_num, Time_order], 1, mean_fc, grp0 = conditions[1], grp1 = conditions[2], grp = pData(lung_TPM_count_cds_14_18[1:transcript_num, Time_order])$Time)
TPM_split_fc <- split(t(TPM_fc), col(t(TPM_fc), as.factor = T))
TPM_permutate_pval <- mcmapply(permuation_pval, TPM_split_cds, TPM_split_fc, mc.cores = detectCores()) #multiple cores 

TPM_sf_mat <- matrix(rep(sizeFactors(lung_TPM_count_cds_14_18), nrow(lung_TPM_count_cds_14_18)), nrow = nrow(lung_TPM_count_cds_14_18), byrow = T)
TPM_split_cds <- split(t(round(exprs(lung_TPM_count_cds_14_18[1:transcript_num, Time_order]) / TPM_sf_mat[1:transcript_num, Time_order])), col(t(exprs(lung_TPM_count_cds_14_18[1:transcript_num, Time_order])), as.factor = T))
TPM_round_fc <- apply(round(exprs(lung_TPM_count_cds_14_18[1:transcript_num, Time_order]) / TPM_sf_mat[1:transcript_num, Time_order]), col(t(exprs(lung_TPM_count_cds_14_18[1:transcript_num, Time_order]))), 1, mean_fc, grp0 = conditions[1], grp1 = conditions[2], pData(lung_TPM_count_cds_14_18[1:transcript_num, Time_order])$Time) #valid_gene_id
TPM_round_split_fc <- split(t(TPM_round_fc), col(t(TPM_round_fc), as.factor = T))
TPM_round_permutate_pval <- mcmapply(permuation_pval, TPM_split_cds, TPM_round_split_fc, mc.cores = detectCores()) #multiple core

lung_TPM_count_edgeR_res <- edgeR_test(counts = exprs(lung_TPM_count_cds_14_18)[1:transcript_num, ], glm = F)
lung_TPM_count_deseq2_res <- DESeq2_deg(dir = NULL, lung_TPM_count_cds_14_18[1:transcript_num, ], Time = Time_ori, pd = pData(count_cds))

###############################################################################################################################################################################################
#make the new figure: 
#new monocle_p: 
tpm_count_monocle_p <- lung_TPM_count_diff_test_res[, 'pval'] 
names(tpm_count_monocle_p) <- row.names(lung_TPM_count_diff_test_res)

#deseq
tpm_count_deseq_p <- lung_TPM_count_dtable_pool_max_nbinomGLMTest$dtalbe[, 'pval'] #std_dtable_pool_max_nbinomTest
names(tpm_count_deseq_p) <- row.names(lung_TPM_count_dtable_pool_max_nbinomGLMTest$dtalbe) #std_dtable_pool_max_nbinomGLMTest

# scde
tpm_count_scde_p <- lung_TPM_count_scde_res_list$pval 
tpm_count_scde_p <- tpm_count_scde_p[names(tpm_count_monocle_p)]
#prepare the data for the DESeq2 / edgeR: 
tpm_count_edgeR_p_glm = lung_TPM_count_edgeR_res$et$table$PValue
names(tpm_count_edgeR_p_glm) <- row.names(lung_TPM_count_edgeR_res$et$table)

tpm_count_deseq2_p = results(lung_TPM_count_deseq2_res)$pvalue
names(tpm_count_deseq2_p) <- row.names(results(lung_TPM_count_deseq2_res))
 
#generate the dataframe for making the benchmarking plots: 
df3 <- plot_pre_rec_f1(test_p_list = list(tpm_count_monocle_p = tpm_count_monocle_p, 
                                        tpm_count_edgeR_p_glm = tpm_count_edgeR_p_glm, 
                                        tpm_count_deseq2_p = tpm_count_deseq2_p, 
                                        tpm_count_deseq_p = tpm_count_deseq_p, 
                                        tpm_count_scde_p = tpm_count_scde_p),
                     permutate_pval = list(tpm_count_monocle_p = TPM_permutate_pval, #readcount_permutate_pval, #std_permutate_pval, 
                                           tpm_count_edgeR_p_glm = TPM_permutate_pval, 
                                           tpm_count_deseq2_p = TPM_permutate_pval,
                                           tpm_count_deseq_p = TPM_permutate_pval,
                                           tpm_count_scde_p = TPM_permutate_pval),
                     row.names(absolute_cds[, ]), #1:transcript_num gene_list, overlap_genes, high_gene_list
                     return_df = T, #na.rm = T, 
                     p_thrsld = 0.05, #0.05
                     rownames = c('monocle (TPM counts)', 
                        'edgeR (TPM counts)', 'DESeq2 (TPM counts)',
                        'DESeq (TPM counts)', 'SCDE (TPM counts)'))
df3$data_type = c("TPM counts", "TPM counts", "TPM counts", "TPM counts", "TPM counts")

df3$class = '3relative'


# only show new size normalization: 
df3.1 <- df3
df3.1[, 'Type'] <- c('SCDE', 'DESeq1', 'DESeq2', 'edgeR', 'Monocle') # geom_bar(stat = 'identity', position = 'dodge') 

pdf('./supplementary_figures/fig2a_si_tpm_counts.pdf', width = 3, height = 2)
qplot(factor(Type), value, stat = "identity", geom = 'bar', position = 'dodge', fill = data_type, data = melt(df3.1)) + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + monocle_theme_opts() + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) + nm_theme()
dev.off()

pdf('./supplementary_figures/fig2a_si_tpm_counts_helper.pdf', width = 3, height = 2)
qplot(factor(Type), value, stat = "identity", geom = 'bar', position = 'dodge', fill = data_type, data = melt(df3.1)) + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + monocle_theme_opts() + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) #+ nm_theme()
dev.off()

save.image('./RData/TPM_median_total_mRNA_analysis.RData')

# ###############################################################################################################################################################################################
#add the test without size factor normalization: 
no_sf_mode_size_norm_permutate_ratio_by_geometric_mean <- cal_perm_pval_size_norm(new_abs_cds_14_18[1:transcript_num, Time_order])
no_sf_mc_mode_size_norm_permutate_ratio_by_geometric_mean <- cal_perm_pval_size_norm(new_mc_cds_14_18[1:transcript_num, Time_order])

mc_sf_mat <- matrix(rep(sizeFactors(new_mc_cds_14_18), nrow(new_mc_cds_14_18)), nrow = nrow(new_mc_cds_14_18), byrow = T)
mc_split_cds <- split(t(exprs(new_mc_cds_14_18[1:transcript_num, Time_order]) / mc_sf_mat[1:transcript_num, Time_order]), col(t(exprs(new_mc_cds_14_18[1:transcript_num, Time_order])), as.factor = T))
mc_fc <- apply(exprs(new_mc_cds_14_18)[1:transcript_num, Time_order] / mc_sf_mat[1:transcript_num, Time_order], 1, mean_fc, grp0 = conditions[1], grp1 = conditions[2], grp = pData(new_mc_cds_14_18[1:transcript_num, Time_order])$Time) #valid_gene_id
mc_split_fc <- split(t(mc_fc), col(t(mc_fc), as.factor = T))
mc_mode_size_norm_permutate_ratio_by_geometric_mean <- mcmapply(permuation_pval, mc_split_cds, mc_split_fc, mc.cores = detectCores()) #multiple core
closeAllConnections()

no_sf_readcount_split_cds <- split(t(exprs(count_cds[1:transcript_num, Time_order])), col(t(exprs(count_cds[1:transcript_num, Time_order])), as.factor = T))
no_sf_readcount_fc <- apply(exprs(count_cds[1:transcript_num, Time_order]), 1, mean_fc, grp0 = conditions[1], grp1 = conditions[2], grp = pData(count_cds[, Time_order])$Time) #valid_gene_id
no_sf_readcount_split_fc <- split(t(readcount_fc), col(t(readcount_fc), as.factor = T))
no_sf_readcount_permutate_pval <- mcmapply(permuation_pval, no_sf_readcount_split_cds, no_sf_readcount_split_fc, mc.cores = detectCores()) #multiple core


TPM_cds_14_18 <- TPM_cds[, colnames(new_mc_cds_14_18)]
#do the analysis for the TPM values: 
no_sf_TPM_split_cds <- split(t(exprs(TPM_cds_14_18[1:transcript_num, Time_order])), col(t(exprs(TPM_cds_14_18[1:transcript_num, Time_order])), as.factor = T))
no_sf_TPM_fc <- apply(exprs(TPM_cds_14_18[1:transcript_num, Time_order]), 1, mean_fc, grp0 = conditions[1], grp1 = conditions[2], grp = pData(TPM_cds_14_18[, Time_order])$Time) #valid_gene_id
no_sf_TPM_split_fc <- split(t(no_sf_TPM_fc), col(t(no_sf_TPM_fc), as.factor = T))

no_sf_TPM_permutate_pval <- mcmapply(permuation_pval, no_sf_TPM_split_cds, no_sf_TPM_split_fc, mc.cores = detectCores()) #multiple core

#do the anlaysis for the TPM with size factor: 
TPM_split_cds <- split(t(exprs(TPM_cds_14_18[1:transcript_num, Time_order]) / readcount_sf_mat[1:transcript_num, Time_order]), col(t(exprs(TPM_cds_14_18[1:transcript_num, Time_order])), as.factor = T))
TPM_fc <- apply(exprs(TPM_cds_14_18[1:transcript_num, Time_order] / readcount_sf_mat[1:transcript_num, Time_order]), 1, mean_fc, grp0 = conditions[1], grp1 = conditions[2], grp = pData(TPM_cds_14_18[, Time_order])$Time) #valid_gene_id
TPM_split_fc <- split(t(readcount_fc), col(t(readcount_fc), as.factor = T))

TPM_permutate_pval <- mcmapply(permuation_pval, TPM_split_cds, TPM_split_fc, mc.cores = detectCores()) #multiple core

save(file = './no_size_permutation.RData', 
  no_sf_mode_size_norm_permutate_ratio_by_geometric_mean, no_sf_mc_mode_size_norm_permutate_ratio_by_geometric_mean, no_sf_readcount_split_cds, no_sf_readcount_fc, no_sf_readcount_split_fc, no_sf_readcount_permutate_pval)
#analysis from Census counts: 
mc_scde_res_list <- scde_DEG(dir = NULL, count_cds = new_mc_cds_14_18, DEG_attribute = 'Time', contrast = conditions, n.cores = 1, normalize = T)
save(file = './RData/mc_scde_res_list', mc_scde_res_list)

mc_scde_res_list_no_normalize <- scde_DEG(dir = NULL, count_cds = new_mc_cds_14_18, DEG_attribute = 'Time', contrast = conditions, n.cores = 1, normalize = F)
save(file = './RData/mc_scde_res_list_no_normalize', mc_scde_res_list_no_normalize)

new_mc_cds_14_18_d <- newCountDataSet(round(exprs(new_mc_cds_14_18))[, ], Time_ori) 
mc_dtable_pool_max_nbinomGLMTest <- DESeq1_test(new_mc_cds_14_18_d[1:transcript_num, ], pd = pData(new_mc_cds_14_18), disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomGLMTest', scale = T) 
row.names(mc_dtable_pool_max_nbinomGLMTest$dtalbe) <-  row.names(new_mc_cds_14_18[1:transcript_num, ])
mc_edgeR_res_glm <- edgeR_test(exprs(new_mc_cds_14_18[1:transcript_num, ]), group = Time_ori, glm = T)
mc_edgeR_res <- edgeR_test(exprs(new_mc_cds_14_18[1:transcript_num, ]), group = Time_ori, glm = F)
mc_deseq2_res <- DESeq2_deg(dir = NULL, new_mc_cds_14_18[1:transcript_num, ], Time = Time_ori, pd = pData(new_mc_cds_14_18))

TPM_scde_res_list_no_normalize <- scde_DEG(dir = NULL, count_cds = TPM_cds_14_18, DEG_attribute = 'Time', contrast = conditions, n.cores = 1, normalize = F)
new_tpm_cds_14_18_d <- newCountDataSet(round(exprs(TPM_cds_14_18))[, ], Time_ori) 
tpm_dtable_pool_max_nbinomGLMTest <- DESeq1_test(TPM_cds_14_18[1:transcript_num, ], pd = pData(TPM_cds_14_18), disp_method = 'pooled', sharing = 'maximum', test_type = 'nbinomGLMTest', scale = T) 
row.names(tpm_dtable_pool_max_nbinomGLMTest$dtalbe) <-  row.names(TPM_cds_14_18[1:transcript_num, ])
tpm_edgeR_res_glm <- edgeR_test(exprs(TPM_cds_14_18[1:transcript_num, ]), group = Time_ori, glm = T)
tpm_deseq2_res <- DESeq2_deg(dir = NULL, TPM_cds_14_18[1:transcript_num, ], Time = Time_ori, pd = pData(TPM_cds_14_18))

#pseudotime test on lineage 2/3:    
std_pseudotime_test_lineage2_res <- differentialGeneTest(std_AT12_cds_subset_all_gene[1:transcript_num, pData(std_AT12_cds_subset_all_gene)$State %in% c(1, 2)], 
                                              fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                              reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)  
std_pseudotime_test_lineage3_res <- differentialGeneTest(std_AT12_cds_subset_all_gene[1:transcript_num, pData(std_AT12_cds_subset_all_gene)$State %in% c(1, 3)], 
                                                      fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                                      reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)
abs_pseudotime_test_lineage2_res <- differentialGeneTest(abs_AT12_cds_subset_all_gene[1:transcript_num, pData(abs_AT12_cds_subset_all_gene)$State %in% c(1, 2)], 
                                                fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                                reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)  
abs_pseudotime_test_lineage3_res <- differentialGeneTest(abs_AT12_cds_subset_all_gene[1:transcript_num, pData(abs_AT12_cds_subset_all_gene)$State %in% c(1, 3)], 
                                                      fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                                      reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)
mc_pseudotime_test_lineage2_res <- differentialGeneTest(mc_AT12_cds_subset_all_gene[1:transcript_num, pData(mc_AT12_cds_subset_all_gene)$State %in% c(1, 2)], 
                                                fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                                reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)    
mc_pseudotime_test_lineage3_res <- differentialGeneTest(mc_AT12_cds_subset_all_gene[1:transcript_num, pData(mc_AT12_cds_subset_all_gene)$State %in% c(1, 3)], 
                                                      fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                                      reducedModelFormulaStr = "~1", cores = detectCores(), relative = T)









