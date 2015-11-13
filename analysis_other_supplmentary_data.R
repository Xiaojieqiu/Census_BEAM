  # load('prepare_lung_data.RData')

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
  # std_gd_fit_res <- mcesApply(standard_cds[, ], 1, gd_fit_pval, cores = detectCores(), required_packages = c('VGAM', 'fitdistrplus', 'MASS', 'pscl'), exprs_thrsld = 10, pseudo_cnt = 0.01)
  # closeAllConnections()
  # tpm_gd_fit_res <- mcesApply(TPM_cds[, ], 1, gd_fit_pval, cores = detectCores(), required_packages = c('VGAM', 'fitdistrplus', 'MASS', 'pscl'), exprs_thrsld = 10, pseudo_cnt = 0.01)
  # closeAllConnections()
  read_gd_fit_res <- mcesApply(read_countdata_cds[, ], 1, gd_fit_pval, cores = detectCores(), required_packages = c('VGAM', 'fitdistrplus', 'MASS', 'pscl'), exprs_thrsld = 10, pseudo_cnt = 0.01)

  abs_gd_fit_res <- unlist(abs_gd_fit_res) 
  read_gd_fit_res <- unlist(read_gd_fit_res) 

  abs_gd_fit_df <- matrix(abs_gd_fit_res, nrow(absolute_cds), ncol = 11, byrow = T)
  dimnames(abs_gd_fit_df) <- list(row.names(absolute_cds), c("ln_pvalue", "nb_pvalue", "ln_pvalue.glm.link", "ln_pvalue.glm.log", "ln_pvalue.chisq", "nb_pvalue.glm", "nb_pvalue.chisq", "zinb_pvalue.chisq", "zanb_pvalue.chisq", "zinb_pvalue", "zanb_pvalue"))
  read_gd_fit_df <- matrix(read_gd_fit_res, nrow(absolute_cds), ncol = 11, byrow = T)
  dimnames(read_gd_fit_df) <- list(row.names(absolute_cds), c("ln_pvalue", "nb_pvalue", "ln_pvalue.glm.link", "ln_pvalue.glm.log", "ln_pvalue.chisq", "nb_pvalue.glm", "nb_pvalue.chisq", "zinb_pvalue.chisq", "zanb_pvalue.chisq", "zinb_pvalue", "zanb_pvalue"))

save('analysis_other_supplementary_data.RData')