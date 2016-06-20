#load the data 
load('./RData/analysis_UMI_data.RData')

library(monocle)
# library(devtools)
# load_all('~/Projects/monocle-dev')
library(xacHelper)
load_all_libraries()

pData(UMI_cds)$Total_mRNAs <- esApply(UMI_cds, 2, sum)
pData(UMI_cds)$endogenous_RNA <- esApply(UMI_cds[-ERCC_ids, ], 2, sum)

#1. test the normalization algorithm by using the TPM dataset converted from UMI data

#convert UMI data to TPM (a relative abundance measure) for the following analysis: 
UMI_TPM <- esApply(UMI_cds, 2, function(x) x / sum(x) * 10^6)
UMI_TPM <- UMI_TPM[, esApply(UMI_cds, 2, sum) > 5e4] #remove cells have terrible recovery as well as TPM distribution 
UMI_tpm_align <- melt(UMI_TPM)
UMI_tpm_align <- UMI_tpm_align[UMI_tpm_align$value > 0.1, ]
  
pdf('./tmp/UMI_tpm_align_distr.pdf', width = 20, height = 20) 
  qplot(value, log = 'x', geom = 'density', color = Var2, data = UMI_tpm_align) + facet_wrap(~Var2) + theme_bw() + theme(legend.position = 'none') + xlab('TPM')
dev.off()

UMI_TPM_cds <- newCellDataSet(UMI_TPM, 
                        phenoData = new("AnnotatedDataFrame", data = pData(UMI_cds[, colnames(UMI_TPM)])),
                        featureData = new("AnnotatedDataFrame", data = fData(UMI_cds)),
                        expressionFamily=tobit(), 
                        lowerDetectionLimit=1)
pData(UMI_TPM_cds)$Total_mRNAs <- esApply(UMI_TPM_cds, 2, sum)
pData(UMI_TPM_cds)$endogenous_RNA <- esApply(UMI_TPM_cds[-ERCC_ids, ], 2, sum)

UMI_TPM_ercc_controls <- newCellDataSet(exprs(UMI_TPM_cds[ERCC_ids, ]), 
                                  phenoData = new("AnnotatedDataFrame", data = pData(UMI_cds[, colnames(UMI_TPM)])),
                                  featureData = new("AnnotatedDataFrame", data = fData(UMI_cds[ERCC_ids, ])), 
                                  expressionFamily=tobit(), 
                                  lowerDetectionLimit=1)

#We use the UMI counts as the ladder instead of ERCC spike-in transcripts: 
#1. Only 50 ERCC species detected by their experiment (92 species was added)
#2. Cell-seq with UMI is lower in recovery efficiency 
#3. Authors from the paper mentioned low recovery efficiency based on spike-in regression (~ 0.025) and potential ERCC degradation

#subset dataframe for ERCC annotation: 
ERCC_ids_names <- row.names(UMI_cds)[ERCC_ids]
subset_spike_df <- spike_df[ERCC_ids_names, ] #spike_df from the monocle package 
identical(ERCC_ids_names, row.names(subset_spike_df))

#there are variations between total UMI counts of the detected ERCC species for each cell, so median is used 
subset_spike_df$numMolecules <- as.numeric(esApply(UMI_cds[ERCC_ids, ], 1, median)) + 1

#regression model between TPM and UMI counts 
UMI_molModels <- esApply(UMI_TPM_ercc_controls, 2, function(cell_exprs, input.ERCC.annotation) {
    #print (cell_exprs)
    spike_df <- input.ERCC.annotation 
    spike_df <- cbind(spike_df, cell_exprs[row.names(spike_df)]) 
    colnames(spike_df)[length(colnames(spike_df))] <- "FPKM"
    #spike_df$numMolecules <- spike_df$conc_attomoles_ul_Mix1*(1000*10^(-3)*1/2500000*10^(-18)*6.02214179*10^(23))
    #spike_df$rounded_numMolecules <- round(spike_df$conc_attomoles_ul_Mix1*(1000*10^(-3)*1/2500000*10^(-18)*6.02214179*10^(23)))
    # print (mean(spike_df$rounded_numMolecules))
    
    spike_df <- subset(spike_df, FPKM >= 1e-10)
    spike_df$log_fpkm <- log10(spike_df$FPKM) 
    spike_df$log_numMolecules <- log10(spike_df$numMolecules)
    
    # print (paste("geometric mean of log numMol ", mean(spike_df$log_numMolecules), "geometric mean of log FPKM ", mean(spike_df$log_fpkm)))
           
   molModel <- tryCatch({
     #molModel <- vgam (rounded_numMolecules ~ sm.ns(log_fpkm, df=3), data=spike_df, family=negbinomial(zero=NULL))
     molModel <- rlm(log_numMolecules ~ log_fpkm, data=spike_df)
     
     #
     #        qp <- qplot(FPKM, numMolecules, data=spike_df, log="xy") +
     #      geom_abline(color="green") +
     # geom_smooth(aes(y=10^log_numMolecules), method="lm") +
     # geom_smooth(aes(y=10^log_numMolecules), method="rlm", color="red") +
     # geom_line(aes(y=10^predict(molModel,type="response")))
     # print(qp)
     molModel
   }, 
   error = function(e) { print(e); NULL })
   molModel
  }, subset_spike_df)

# #calculate the corresponding transcript counts of the mode with the model learnt from the UMI data (from TPM to UMI)
t_estimate <- estimate_t(UMI_TPM_cds)
UMI_mode_absolute_counts <- mapply(function(cell_exprs, molModel) {
tryCatch({
  norm_df <- data.frame(log_fpkm=log10(as.numeric(cell_exprs)))
  res <- 10^predict(molModel, type="response", newdata=norm_df)
}, 
error = function(e) {
  rep(NA, length(cell_exprs))
})
}, 
split(as.vector(t_estimate), names(t_estimate)), 
UMI_molModels) 

pdf('./supplementary_figures/UMI_mode_dist.pdf', width = 2, height = 2)
qplot(UMI_mode_absolute_counts) + xlab('Transcript count for most frequent log10(TPM)') + nm_theme()
dev.off()

#obtain the regression parameters from the UMI_molModels
UMI_models_kd_df <- do.call(rbind.data.frame, lapply(UMI_molModels, function(x) data.frame(b = coef(x)[1], k = coef(x)[2])))

pdf('./supplementary_figures/UMI_kb_line.pdf', width = 2, height = 2)
qplot(k, b, data = UMI_models_kd_df) + geom_smooth(method = 'rlm') + 
    xlab('Slope from TPM vs. Median UMI counts \n for ERCC spike-in transcripts') + ylab('Intercept from TPM vs. Median UMI counts \n for ERCC spike-in transcripts') +  nm_theme()
dev.off()

rlm(b ~ k, data = UMI_models_kd_df) 

UMI_norm_recovery_all_correct_mc <- relative2abs(UMI_TPM_cds, expected_total_mRNAs = pData(UMI_cds[, colnames(UMI_TPM_cds)])$Total_mRNAs, 
                                                 t_estimate = estimate_t(exprs(UMI_TPM_cds)),  kb_slope = -1.9266285 , kb_slope_rng = c(-2.1, -1.9), 
                                                 kb_intercept = 0.8270289, kb_intercept_rng = c(0.8270289, 0.8270289), 
                                                 return_all = T, cores = 1)

qplot(k, b, data = as.data.frame(t(UMI_norm_recovery_all_correct_mc$k_b_solution))) + geom_point(aes(x = k, y = b, color = 'red'), data = UMI_models_kd_df)
rlm(b ~k , data = as.data.frame(t(UMI_norm_recovery_all_correct_mc$k_b_solution))) 

pdf('./tmp/UMI_algorithm_total_10e5.pdf', width = 2, height = 2)
qplot(apply(UMI_norm_recovery_all_correct_mc$norm_cds[-ERCC_ids, ], 2, sum), pData(UMI_cds[, colnames(UMI_TPM_cds)])$endogenous_RNA, log = 'xy') + 
    geom_smooth(method = 'rlm') + geom_abline(color = 'red') + xlab('Total endogenous mRNA (algorithm)') + ylab('Total endogenous mRNA (spike-in regression)') + 
    nm_theme()
dev.off()

pdf('./tmp/UMI_algorithm_total.pdf', width = 2, height = 2)
qplot(apply(UMI_norm_recovery_all_correct_mc$norm_cds[-ERCC_ids, ], 2, sum), pData(UMI_cds[, colnames(UMI_TPM_cds)])$endogenous_RNA, log = 'xy') + 
    geom_smooth(method = 'rlm') + geom_abline(color = 'red') + expand_limits(x = c(1, 1e6), y = c(1, 1e06)) + xlab('Total endogenous mRNA (algorithm)') + ylab('Total endogenous mRNA (spike-in regression)') + 
    nm_theme()
dev.off()

pdf('./supplementary_figures/UMI_algorithm_mean.pdf', width = 2, height = 2)
qplot(apply(UMI_norm_recovery_all_correct_mc$norm_cds[-ERCC_ids, ], 1, mean), esApply(UMI_cds[-ERCC_ids,  colnames(UMI_TPM_cds)], 1, mean), log = 'xy') + 
    geom_smooth(method = 'rlm') + xlab('Mean endogenous mRNA count \n (algorithm)') + ylab('Mean endogenous mRNA count (regression)') + 
    nm_theme()
dev.off()

################################################################################################################################################################################################################################
#recover without setting range for m
UMI_norm_recovery_default_c <- relative2abs(UMI_TPM_cds, expected_total_mRNAs = median(pData(UMI_cds[, colnames(UMI_TPM_cds)])$Total_mRNAs), 
                                                 t_estimate = estimate_t(exprs(UMI_TPM_cds)), # m = -1.9266285 , m_rng = c(-2.1, -1.9), 
                                                 kb_intercept = 0.8270289, kb_intercept_rng = c(0.8270289, 0.8270289), #fix_c
                                                 return_all = T, cores = detectCores())

qplot(k, b, data = as.data.frame(t(UMI_norm_recovery_default_c$k_b_solution))) + geom_point(aes(x = k, y = b, color = 'red'), data = UMI_models_kd_df)
rlm(b ~k , data = as.data.frame(t(UMI_norm_recovery_default_c$k_b_solution))) 

pdf('./tmp/UMI_algorithm_total_10e5_optim_m.pdf', width = 2, height = 2)
qplot(apply(UMI_norm_recovery_default_c$norm_cds[-ERCC_ids, ], 2, sum), pData(UMI_cds[, colnames(UMI_TPM_cds)])$endogenous_RNA, log = 'xy') + 
    geom_smooth(method = 'rlm') + geom_abline(color = 'red') + xlab('Total endogenous mRNA (algorithm)') + ylab('Total endogenous mRNA (spike-in regression)') + 
    nm_theme()
dev.off()

pdf('./tmp/UMI_algorithm_total_optim_m.pdf', width = 2, height = 2)
qplot(apply(UMI_norm_recovery_default_c$norm_cds[-ERCC_ids, ], 2, sum), pData(UMI_cds[, colnames(UMI_TPM_cds)])$endogenous_RNA, log = 'xy') + 
    geom_smooth(method = 'rlm') + geom_abline(color = 'red') + expand_limits(x = c(1, 1e6), y = c(1, 1e06)) + xlab('Total endogenous mRNA (algorithm)') + ylab('Total endogenous mRNA (spike-in regression)') + 
    nm_theme()
dev.off()

pdf('./supplementary_figures/UMI_algorithm_mean_optim_m.pdf', width = 2, height = 2)
qplot(apply(UMI_norm_recovery_default_c$norm_cds[-ERCC_ids, ], 1, mean), esApply(UMI_cds[-ERCC_ids,  colnames(UMI_TPM_cds)], 1, mean), log = 'xy') + 
    geom_smooth(method = 'rlm') + xlab('Mean endogenous mRNA count \n (algorithm)') + ylab('Mean endogenous mRNA count (regression)') + 
    nm_theme()
dev.off()

pdf('./supplementary_figures/ERCC_UMI_algorithm_mean_optim_m.pdf', width = 2, height = 2)
qplot(apply(UMI_norm_recovery_default_c$norm_cds[ERCC_ids, ], 1, mean), esApply(UMI_cds[ERCC_ids,  colnames(UMI_TPM_cds)], 1, mean), log = 'xy') + 
    geom_smooth(method = 'rlm') + xlab('Mean endogenous mRNA count \n (algorithm)') + ylab('Mean endogenous mRNA count (regression)') + 
    nm_theme()
dev.off()
################################################################################################################################################################################################################################

save.image('./RData/umi_normalization.RData')

