#test the UMI data using the census: 
#load the data 
load('./RData/analysis_UMI_data.RData')

# library(monocle)
library(devtools)
load_all('~/Projects/monocle-dev')

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
subset_spike_df$numMolecules <- input.ERCC.annotation[row.names(subset_spike_df), 'numMolecules'] #as.numeric(esApply(UMI_cds[ERCC_ids, ], 1, median)) + 1

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
qplot(UMI_mode_absolute_counts) + xlab('Transcript count for most frequent log10(TPM)') + nm_theme() + ylab('Cells')
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
UMI_norm_recovery_default_c <- relative2abs(UMI_TPM_cds, expected_total_mRNAs = median(apply(UMI_absolute_counts, 2, sum)), 
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

#recovery by setting correct capture rate: 
UMI_norm_recovery_capture_rate <- relative2abs(UMI_TPM_cds, expected_total_mRNAs = median(apply(UMI_absolute_counts, 2, sum)), 
                                            t_estimate = estimate_t(exprs(UMI_TPM_cds)), expected_capture_rate = 0.03529354,
                                            return_all = T, cores = detectCores())

qplot(k, b, data = as.data.frame(t(UMI_norm_recovery_capture_rate$k_b_solution))) + geom_point(aes(x = k, y = b, color = 'red'), data = UMI_models_kd_df)
rlm(b ~k , data = as.data.frame(t(UMI_norm_recovery_capture_rate$k_b_solution))) 

pdf('./tmp/UMI_algorithm_total_capture_rate.pdf', width = 2, height = 2)
qplot(apply(UMI_norm_recovery_capture_rate$norm_cds[-ERCC_ids, ], 2, sum), pData(UMI_cds[, colnames(UMI_TPM_cds)])$endogenous_RNA, log = 'xy') + 
  geom_smooth(method = 'rlm') + geom_abline(color = 'red') + xlab('Total endogenous mRNA (algorithm)') + ylab('Total endogenous mRNA (spike-in regression)') + 
  nm_theme()
dev.off()

pdf('./tmp/UMI_algorithm_total_capture_rate.pdf', width = 2, height = 2)
qplot(apply(UMI_norm_recovery_capture_rate$norm_cds[-ERCC_ids, ], 2, sum), pData(UMI_cds[, colnames(UMI_TPM_cds)])$endogenous_RNA, log = 'xy') + 
  geom_smooth(method = 'rlm') + geom_abline(color = 'red') + expand_limits(x = c(1, 1e6), y = c(1, 1e06)) + xlab('Total endogenous mRNA (algorithm)') + ylab('Total endogenous mRNA (spike-in regression)') + 
  nm_theme()
dev.off()

pdf('./supplementary_figures/UMI_algorithm_mean_capture_rate.pdf', width = 2, height = 2)
qplot(apply(UMI_absolute_counts[-ERCC_ids,  ], 1, mean), apply(UMI_norm_recovery_capture_rate$norm_cds[-ERCC_ids, ], 1, mean)) +
  xlim(c(0.1, 1000)) + ylim(c(0.1, 1000)) + scale_x_log10() + scale_y_log10() + geom_abline(color ='red') + 
  geom_smooth(method = 'rlm') + xlab('Mean endogenous mRNA count (regression)') + ylab('Mean endogenous mRNA count \n (algorithm)') + 
  nm_theme()
dev.off()

pdf('./supplementary_figures/ERCC_UMI_algorithm_mean_capture_rate.pdf', width = 2, height = 2)
qplot(apply(UMI_norm_recovery_capture_rate$norm_cds[ERCC_ids, ], 1, mean), apply(UMI_absolute_counts[ERCC_ids,  colnames(UMI_TPM_cds)], 1, mean), log = 'xy') + 
  geom_smooth(method = 'rlm') + xlab('Mean endogenous mRNA count \n (algorithm)') + ylab('Mean endogenous mRNA count (regression)') + 
  nm_theme()
dev.off()

#calculated c from the formula: 
UMI_absolute_counts <- mapply(function(cell_exprs, molModel) {
  tryCatch({
    norm_df <- data.frame(log_fpkm=log10(as.numeric(cell_exprs)))
    res <- 10^predict(molModel, type="response", newdata=norm_df)
  }, 
  error = function(e) {
    rep(NA, length(cell_exprs))
  })
}, 
split(exprs(UMI_TPM_cds), rep(1:ncol(exprs(UMI_TPM_cds)), each = nrow(exprs(UMI_TPM_cds)))), 
UMI_molModels) 

c_mean_y_ij <- mean(log10(subset_spike_df[row.names(UMI_TPM_ercc_controls)[esApply(UMI_TPM_ercc_controls, 1, function(x) all(x >0)) > 0], 'numMolecules']))
UMI_norm_recovery_default_c <- relative2abs(UMI_TPM_cds, expected_total_mRNAs = median(apply(UMI_absolute_counts, 2, sum)), 
                                            t_estimate = estimate_t(exprs(UMI_TPM_cds)), # m = -1.9266285 , m_rng = c(-2.1, -1.9), 
                                            kb_intercept = c_mean_y_ij, kb_intercept_rng = c(c_mean_y_ij, c_mean_y_ij), #fix_c
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
qplot(apply(UMI_absolute_counts[-ERCC_ids,  ], 1, mean), apply(UMI_norm_recovery_default_c$norm_cds[-ERCC_ids, ], 1, mean), log = 'xy') + 
  geom_abline(color = 'red') + 
  geom_smooth(method = 'rlm') + xlab('Mean endogenous mRNA count (regression)') + ylab('Mean endogenous mRNA count \n (algorithm)') + 
  nm_theme()
dev.off()

pdf('./supplementary_figures/ERCC_UMI_algorithm_mean_optim_m.pdf', width = 2, height = 2)
qplot(apply(UMI_norm_recovery_default_c$norm_cds[ERCC_ids, ], 1, mean), apply(UMI_absolute_counts[ERCC_ids,  ], 1, mean), log = 'xy') + 
  geom_smooth(method = 'rlm') + xlab('Mean endogenous mRNA count \n (algorithm)') + ylab('Mean endogenous mRNA count (regression)') + 
  nm_theme()
dev.off()

################################################################################################################################################################################################################################
round(spike_df$conc_attomoles_ul_Mix1*(10*10^(-3)*1/40000*10^(-18)*6.02214179*10^(23)))
input.ERCC.annotation$numMolecules <- input.ERCC.annotation$conc_attomoles_ul_Mix1*(1000*10^(-3)*1/2500000*10^(-18)*6.02214179*10^(23))

#calculate the capture rate: 
#calculate the capture rate: 
#generate the figure for calculating capture efficiency: 
ercc_exprs <- exprs(UMI_cds[ERCC_ids, ])

#fit the data: 
detected_ercc_spikein <- apply(ercc_exprs, 2, function(x) as.numeric(x > 10e-4))

#total number of detected spikein for each cell
num_time_spike_in_detect_df <- data.frame(Time = 1, 
                                          rounded_numMolecules = unique(round(subset_spike_df$numMolecules)),
                                          detected = 0, tau = 0)

#calculate number of observed cases for a particular spikein in a particular time point,
#also calculate tau, the probability for a particular spikein to be detected 
subset_spike_df$numMolecules <- round(subset_spike_df$numMolecules)

for(t in unique(1)){
  for(round_numMolecules in sort(subset_spike_df$numMolecules)){
    print(t)
    print(round_numMolecules)
    ind <- num_time_spike_in_detect_df$Time == t & num_time_spike_in_detect_df$rounded_numMolecules == round_numMolecules
    num_time_spike_in_detect_df[ind, 3] <- 
      sum(detected_ercc_spikein[subset_spike_df$numMolecules == round_numMolecules])
    
    num_spikein <- sum(subset_spike_df$numMolecules == round_numMolecules)
    num_time_spike_in_detect_df[ind, 4] <- num_time_spike_in_detect_df[ind, 3] / (ncol(ercc_exprs) * num_spikein)
  }
}

#write the objective function: 
optim_p_from_tau <- function(x, tau, rounded_numMolecules) {
  res <- 0
  for(i in 1:length(tau)){
    tmp <- (tau[i] - (1 - (1 - x)^rounded_numMolecules[i]))^2
    res <- res + tmp
  }
  res
}

optim_p <- function(num_time_spike_in_detect_df, Time, rounded_numMolecules) {
  tau <- num_time_spike_in_detect_df[num_time_spike_in_detect_df$Time == Time, 4]
  message('initial guess for tau, ', tau)
  optim_res <- optim(par = c(x = tau[1]), optim_p_from_tau,
                     gr = NULL, tau = tau,
                     rounded_numMolecules = rounded_numMolecules,
                     method = c("Brent"), 
                     lower = c(0), 
                     upper = c(1),
                     hessian = FALSE)
  optim_p_val <- optim_res$par
  predict_tau <- rep(0, length(tau))
  optim_predict_tau <- rep(0, length(tau))
  empirical_p <- rep(0, length(tau))
  
  for(i in 1:length(tau)) {
    optim_predict_tau[i] <- (1 - (1 - optim_p_val)^rounded_numMolecules[i])
  }
  
  for(i in 1:length(tau)) {
    predict_tau[i] <- (1 - (1 - tau[1])^rounded_numMolecules[i])
  }
  
  empirical_p <- num_time_spike_in_detect_df[order(num_time_spike_in_detect_df$rounded_numMolecules), 'optim_predicted_tau']
  
  list(optim_p_val = optim_p_val, optim_predict_tau = optim_predict_tau, predict_tau = predict_tau, empirical_p = empirical_p)
}

#show the results: 
num_time_spike_in_detect_df <- subset(num_time_spike_in_detect_df, rounded_numMolecules > 0 & rounded_numMolecules < 1000)
optim_res <- optim_p(num_time_spike_in_detect_df, 1, unique(num_time_spike_in_detect_df$rounded_numMolecules))

num_time_spike_in_detect_df$optim_predicted_tau <- optim_res$optim_predict_tau
num_time_spike_in_detect_df$predicted_tau <- optim_res$predicted_tau
num_time_spike_in_detect_df$optim_p_val <- optim_res$optim_p_val
num_time_spike_in_detect_df$empirical_p <- optim_res$empirical_p

#make the plots: 
ggplot(aes(rounded_numMolecules, tau), data = num_time_spike_in_detect_df) + 
  geom_line(aes(rounded_numMolecules, optim_predicted_tau)) + scale_x_log10() 
ggplot(aes(rounded_numMolecules, tau), data = num_time_spike_in_detect_df) + 
  geom_line(aes(rounded_numMolecules, predicted_tau)) + scale_x_log10() 

pdf('./supplementary_figures//umi_sequencing_efficiency.pdf', width = 2, height = 2)
ggplot(aes(rounded_numMolecules, tau), data = num_time_spike_in_detect_df) + geom_point(aes(color = Time), alpha = 0.5, size = 1) + 
  geom_line(aes(rounded_numMolecules, optim_predicted_tau, color = Time), size = 0.2) + scale_x_log10() + nm_theme() + xlab('Spike-in molecules') + ylab('Probability of observation (P)')
dev.off()

pdf('./supplementary_figures//umi_sequencing_efficiency_helper.pdf', width = 2, height = 2)
ggplot(aes(rounded_numMolecules, tau), data = num_time_spike_in_detect_df) + geom_point(aes(color = Time), size = 1) + 
  geom_line(aes(rounded_numMolecules, optim_predicted_tau, color = Time), size = 0.2) + scale_x_log10() 
dev.off()

################################################################################################################################################################################################################################
save.image('./RData/umi_normalization.RData')