  library(monocle)
  library(xacHelper)
  elife_directory = "./"

  load_all_libraries()
  
#all functions for the supplementary files: 

#Supplementary files:

#########################################################################################################
  #test this: 
  abs_gd_fit_res <- mcesApply(absolute_cds[ ], 1, gd_fit_pval, cores = detectCores(), required_packages = c('VGAM', 'fitdistrplus', 'MASS', 'pscl'), exprs_thrsld = 10, pseudo_cnt = 0.01)
  closeAllConnections()
  std_gd_fit_res <- mcesApply(standard_cds[, ], 1, gd_fit_pval, cores = detectCores(), required_packages = c('VGAM', 'fitdistrplus', 'MASS', 'pscl'), exprs_thrsld = 10, pseudo_cnt = 0.01)
  closeAllConnections()
  tpm_gd_fit_res <- mcesApply(TPM_cds[, ], 1, gd_fit_pval, cores = detectCores(), required_packages = c('VGAM', 'fitdistrplus', 'MASS', 'pscl'), exprs_thrsld = 10, pseudo_cnt = 0.01)
  closeAllConnections()
  read_gd_fit_res <- mcesApply(count_cds[, ], 1, gd_fit_pval, cores = detectCores(), required_packages = c('VGAM', 'fitdistrplus', 'MASS', 'pscl'), exprs_thrsld = 10, pseudo_cnt = 0.01)

  abs_gd_fit_res <- unlist(abs_gd_fit_res) 
  read_gd_fit_res <- unlist(abs_gd_fit_res) 

  abs_gd_fit_df <- matrix(abs_gd_fit_res, nrow(absolute_cds), ncol = 11, byrow = T)
  dimnames(abs_gd_fit_df) <- list(row.names(absolute_cds), c("ln_pvalue", "nb_pvalue", "ln_pvalue.glm.link", "ln_pvalue.glm.log", "ln_pvalue.chisq", "nb_pvalue.glm", "nb_pvalue.chisq", "zinb_pvalue.chisq", "zanb_pvalue.chisq", "zinb_pvalue", "zanb_pvalue"))
  read_gd_fit_df <- matrix(read_gd_fit_res, nrow(absolute_cds), ncol = 11, byrow = T)
  dimnames(read_gd_fit_df) <- list(row.names(absolute_cds), c("ln_pvalue", "nb_pvalue", "ln_pvalue.glm.link", "ln_pvalue.glm.log", "ln_pvalue.chisq", "nb_pvalue.glm", "nb_pvalue.chisq", "zinb_pvalue.chisq", "zanb_pvalue.chisq", "zinb_pvalue", "zanb_pvalue"))

  #select only nb and zinb and calculate the number of genes pass goodness of fit and number of genes can be fitted: 
  abs_gd_fit_res <- cal_gd_statistics(abs_gd_fit_df[, c('nb_pvalue', 'zinb_pvalue')], percentage = F, type = 'absolute')
  readcount_gd_fit_res <- cal_gd_statistics(read_gd_fit_df[, c('nb_pvalue', 'zinb_pvalue')], percentage = F,  type = 'readcount')
  gd_fit_res <- rbind(abs_gd_fit_res, readcount_gd_fit_res)
  gd_fit_res <- cbind(gd_fit_res, data_type = row.names(gd_fit_res))
  row.names(gd_fit_res) <- NULL
  gd_fit_res <- as.data.frame(gd_fit_res)
  
  gd_fit_res_num <- subset(gd_fit_res, data_type == 'gd_fit_num')
  gd_fit_res_success_num <- subset(gd_fit_res, data_type == 'success_fit_num')

