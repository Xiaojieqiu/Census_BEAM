# library(monocle)
library(devtools)
load_all('~/Projects/monocle-dev')

library(DDRTree)
library(R.matlab)
library(igraph)
library(MASS)
library(reshape2)
library(xacHelper)

rpkm <- read.table("./data/Embryo_trilineage/rpkm.txt", header = T, row.names = 1)
ercc.rpkm <- read.table("./data/Embryo_trilineage/ercc.rpkm.txt", header = T, row.names = 1)
ercc.amount <- read.table("./data/Embryo_trilineage/E-MTAB-3929.additional.1/ERCCamounts.txt", header = T, row.names = 1, sep = '\t')
#EMTAB_idf <-  read.table("/Users/xqiu/Dropbox (Personal)/Single-cell RNA-seq/Embryo_trilineage/E-MTAB-3929.idf.txt", header = T, row.names = 1, sep = '\t')
EMTAB_sdrf <-  read.table("./data/Embryo_trilineage/E-MTAB-3929.sdrf.txt", header = T, row.names = 1, sep = '\t')

# writeMat('rpkm.Mat', rpkm = as.matrix(rpkm))

# #make a new cds: 
f_df <- data.frame(gene_short_name = row.names(rpkm), row.names = row.names(rpkm))
fd <- new("AnnotatedDataFrame", data = f_df)
pd <- new("AnnotatedDataFrame", data = EMTAB_sdrf)

# Now, make a new CellDataSet using the RNA counts
embryo_rpkm <- newCellDataSet(as.matrix(rpkm[row.names(f_df), row.names(EMTAB_sdrf)]),
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit=1,
                              expressionFamily=tobit())

embryo_rpkm <- detectGenes(embryo_rpkm, min_expr = 0)
qplot(pData(embryo_rpkm)$num_genes_expressed, color = pData(embryo_rpkm)$`Characteristics.developmental.stage.`, geom = 'density', log = 'x')

#convert to absolute transcript counts:
ercc.amount_all <- merge(ercc.amount, spike_df, by.x = 'row.names', by.y = 'row.names')
row.names(ercc.amount_all) <- ercc.amount_all$Row.names

# Fit a robust linear model of spike concentration vs. FPKM for each cell
identical(row.names(ercc.rpkm), row.names(spike_df)) #test the row.names are the same
ESC_molModels <- apply(ercc.rpkm, 2, function(cell_exprs, spike_df) {
  cell_exprs <- as.numeric(cell_exprs)
  spike_df$fpkm <- cell_exprs
  spike_df <- spike_df[cell_exprs > 1e-10, ]
  spike_df$log_fpkm <- log10(cell_exprs[cell_exprs > 1e-10]) 
  spike_df$log_numMolecules <- log10(as.numeric(spike_df$Number.of.molecules))
  spike_df <- subset(spike_df, Number.of.molecules > 0)
  
  # print (paste("geometric mean of log numMol ", mean(spike_df$log_numMolecules), "geometric mean of log FPKM ", mean(spike_df$log_fpkm)))
  
  molModel <- tryCatch({
    #molModel <- vgam (rounded_numMolecules ~ sm.ns(log_fpkm, df=3), data=spike_df, family=negbinomial(zero=NULL))
    molModel <- rlm(log_numMolecules ~ log_fpkm, data=spike_df)
    
    molModel
  }, 
  error = function(e) { print(e); NULL })
  molModel
}, ercc.amount_all)

#test the relationship between k/b: 
ESC_kb_df <- t(rbind.data.frame(lapply(ESC_molModels, function(x) {
  if(is.null(x)) 
    c(0, 0)
  else
    c(b = coef(x)[1], k = coef(x)[2])
} )))
colnames(ESC_kb_df) <- c('b', 'k')

# Now use the per-cell linear models to produce a matrix of absolute transcript abundances 
# for each gene in the genome, in each cell
norm_fpkms <- mapply(function(cell_exprs, molModel) {
  tryCatch({
    norm_df <- data.frame(log_fpkm=log10(cell_exprs))
    res <- 10^predict(molModel, type="response", newdata=norm_df)
  }, 
  error = function(e) {
    rep(NA, length(cell_exprs))
  })
}, 
split(exprs(embryo_rpkm), rep(1:ncol(exprs(embryo_rpkm)), each = nrow(exprs(embryo_rpkm)))), 
ESC_molModels)

row.names(norm_fpkms) <- row.names(exprs(embryo_rpkm))
colnames(norm_fpkms) <- colnames(exprs(embryo_rpkm))

norm_fpkms_melted <- melt(norm_fpkms)

# Now let's generate a new CellDataSet that uses the absolute transcript counts
fpkm_matrix_abs <- norm_fpkms
colnames(fpkm_matrix_abs) <- colnames(embryo_rpkm)
row.names(fpkm_matrix_abs) <- row.names(embryo_rpkm)
#fpkm_matrix <- fpkm_matrix[,-1]

absolute_embryo_cds <- newCellDataSet(as.matrix(fpkm_matrix_abs), 
                                      phenoData = pd, 
                                      featureData = fd, 
                                      expressionFamily=negbinomial(), 
                                      lowerDetectionLimit=1)

pData(absolute_embryo_cds)$Total_mRNAs <- esApply(absolute_embryo_cds, 2, sum)
pData(absolute_embryo_cds)$endogenous_RNA <- esApply(absolute_embryo_cds, 2, function(x) sum(x))

valid_cells <- intersect(row.names(subset(as.data.frame(ESC_kb_df), k > 0 & k < 2)), row.names(subset(pData(absolute_embryo_cds), Total_mRNAs < 1e8)))
qplot(k, b, data = as.data.frame(ESC_kb_df[valid_cells, ]), alpha = 0.5, size = residuals(rlm(b ~ k, data = as.data.frame(ESC_kb_df[valid_cells, ]))) ) + geom_smooth(method = 'rlm') +
  xlab('Slope from RPKM vs. \n for ERCC spike-in transcripts') + ylab('Intercept from RPKM vs. \n for ERCC spike-in transcripts') +  nm_theme()

valid_cells <- row.names(ESC_kb_df[valid_cells[residuals(rlm(b ~ k, data = as.data.frame(ESC_kb_df[valid_cells, ]))) < 0.6], ])
qplot(k, b, data = as.data.frame(ESC_kb_df[valid_cells, ]), alpha = 0.5, size = residuals(rlm(b ~ k, data = as.data.frame(ESC_kb_df[valid_cells, ]))) ) + geom_smooth(method = 'rlm') +
  xlab('Slope from RPKM vs. \n for ERCC spike-in transcripts') + ylab('Intercept from RPKM vs. \n for ERCC spike-in transcripts') +  nm_theme()

rlm(b ~ k, data = as.data.frame(ESC_kb_df[valid_cells, ]))
valid_absolute_embryo_cds <- absolute_embryo_cds[, valid_cells]

# valid_absolute_embryo_cds <- estimateSizeFactors(valid_absolute_embryo_cds)
# valid_absolute_embryo_cds <- estimateDispersions(valid_absolute_embryo_cds)
# valid_absolute_embryo_cds <- detectGenes(valid_absolute_embryo_cds, min_expr = 0)

# qplot(pData(valid_absolute_embryo_cds)$Total_mRNAs, color = pData(valid_absolute_embryo_cds)$`Characteristics.developmental.stage.`, geom = 'density', log = 'x')
# qplot(pData(valid_absolute_embryo_cds)$num_genes_expressed, color = pData(valid_absolute_embryo_cds)$`Characteristics.developmental.stage.`, geom = 'density', log = 'x')

pdf('./supplementary_figures/ESC_kb_line.pdf', width = 2, height = 2)
qplot(k, b, data = as.data.frame(ESC_kb_df[valid_cells, ]), alpha = 0.5) + geom_smooth(method = 'rlm') +
  xlab('Slope from RPKM vs. \n for ERCC spike-in transcripts') + ylab('Intercept from RPKM vs. \n for ERCC spike-in transcripts') +  nm_theme()
dev.off()
#calculate the capture rate: 
#generate the figure for calculating capture efficiency: 
ercc_exprs <- ercc.rpkm[row.names(ercc.amount), valid_cells]

#fit the data: 
detected_ercc_spikein <- apply(ercc_exprs, 2, function(x) as.numeric(x > 10e-4))

#total number of detected spikein for each cell
num_time_spike_in_detect_df <- data.frame(Time = 1, 
                                          rounded_numMolecules = unique(ercc.amount$Number.of.molecules),
                                          detected = 0, tau = 0)

#calculate number of observed cases for a particular spikein in a particular time point,
#also calculate tau, the probability for a particular spikein to be detected 
for(t in unique(1)){
  for(round_numMolecules in sort(ercc.amount$Number.of.molecules)){
    print(t)
    print(round_numMolecules)
    ind <- num_time_spike_in_detect_df$Time == t & num_time_spike_in_detect_df$rounded_numMolecules == round_numMolecules
    num_time_spike_in_detect_df[ind, 3] <- 
      sum(detected_ercc_spikein[ercc.amount$Number.of.molecules == round_numMolecules])
    
    num_spikein <- sum(ercc.amount$Number.of.molecules == round_numMolecules)
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

pdf('./supplementary_figures/ESC_sequencing_efficiency.pdf', width = 2, height = 2)
ggplot(aes(rounded_numMolecules, tau), data = num_time_spike_in_detect_df) + geom_point(aes(color = Time), alpha = 0.5, size = 1) + 
  geom_line(aes(rounded_numMolecules, optim_predicted_tau, color = Time), size = 0.2) + scale_x_log10() + nm_theme() + xlab('Spike-in molecules') + ylab('Probability of observation (P)”')
dev.off()

pdf('./supplementary_figures/ESC_sequencing_efficiency_helper.pdf', width = 2, height = 2)
ggplot(aes(rounded_numMolecules, tau), data = num_time_spike_in_detect_df) + geom_point(aes(color = Time), size = 1) + 
  geom_line(aes(rounded_numMolecules, optim_predicted_tau, color = Time), size = 0.2) + scale_x_log10() 
dev.off()

t_estimate <- estimate_t(embryo_rpkm[, valid_cells])
ESC_mode_absolute_counts <- mapply(function(cell_exprs, molModel) {
  tryCatch({
    norm_df <- data.frame(log_fpkm=log10(as.numeric(cell_exprs)))
    res <- 10^predict(molModel, type="response", newdata=norm_df)
  }, 
  error = function(e) {
    rep(NA, length(cell_exprs))
  })
}, 
split(as.vector(t_estimate), names(t_estimate)), 
ESC_molModels[valid_cells]) 

ESC_mode_df <- data.frame(Lysate = ESC_mode_absolute_counts, cDNA = ESC_mode_absolute_counts * optim_res$optim_p_val)
mlt_ESC_mode_df <- melt(ESC_mode_df)
mlt_ESC_mode_df$variable <- as.character(mlt_ESC_mode_df$variable)
mlt_ESC_mode_df[mlt_ESC_mode_df$variable == 'Lysate', 'variable'] <- 'Lysate RNA'
# cell lysate
# cDNA #4682B4

pdf('./supplementary_figures/ESC_mode_dist.pdf', width = 2, height = 2)
qplot(value, data = mlt_ESC_mode_df, fill = variable, geom = 'density', log = 'x', alpha = 0.7) + 
  facet_wrap(~variable, scale = 'free_y', ncol = 1)  + xlab('Transcript count for most frequent \nlog10(TPM)') + 
  nm_theme() + ylab('Cells') + scale_fill_manual(values = alpha(c("#ED1F24", "#4682B4"), 1))
dev.off()

# abs_embryo_mat <- relative2abs(embryo_rpkm[, valid_cells], 
#                                    t_estimate = estimate_t(embryo_rpkm[, valid_cells]), 
#                                    cores=detectCores(), 
#                                    verbose = T)
# c_mean_y_ij <- mean(log10(ercc.amount[ercc.amount_all$conc_attomoles_ul_Mix1 > 80, 'Number.of.molecules']))
# ESC_norm_recovery_default_c  <- relative2abs(embryo_rpkm[, valid_cells],  
#                                             t_estimate = estimate_t(exprs(embryo_rpkm)[, valid_cells_3]), # m = -1.9266285 , m_rng = c(-2.1, -1.9), 
#                                             kb_intercept = c_mean_y_ij, kb_intercept_rng = c(c_mean_y_ij,  c_mean_y_ij), #fix_c
#                                             return_all = T, cores = detectCores())
c_mean_y_ij <- mean(log10(ercc.amount[ercc.amount_all$Number.of.molecules > 1, 'Number.of.molecules']))

ESC_norm_recovery_capture_rate <- relative2abs(embryo_rpkm[, valid_cells], expected_total_mRNAs = median(apply(absolute_embryo_cds[, valid_cells], 2, sum)), 
                                               calibrate_total_mRNA = F, 
                                               kb_intercept = 2.010404, kb_intercept_rng = c(2.010404,  2.010404), #fix_c
                                               expected_capture_rate = optim_res$optim_p_val,
                                               verbose = T, 
                                               return_all = T, cores = detectCores())

qplot(k, b, data = as.data.frame(t(ESC_norm_recovery_capture_rate$k_b_solution))) + geom_point(aes(x = k, y = b, color = 'red'), data = as.data.frame(ESC_kb_df[valid_cells, ]))
rlm(b ~k , data = as.data.frame(t(ESC_norm_recovery_capture_rate$k_b_solution))) 
qplot(ESC_mode_absolute_counts, ESC_norm_recovery_capture_rate$calibrated_mode$hypothetical_mode, log = 'xy')

pData(absolute_embryo_cds)$Total_mRNAs <- esApply(absolute_embryo_cds, 2, sum)
pdf('./supplementary_figures/ESC_algorithm_total_10e5_optim_m.pdf', width = 2, height = 2)
qplot(esApply(absolute_embryo_cds[, valid_cells], 2, sum), apply(ESC_norm_recovery_capture_rate$norm_cds, 2, sum), log = 'xy') + 
  geom_smooth(method = 'rlm') + geom_abline(color = 'red') + xlab('Total endogenous mRNA (algorithm)') + ylab('Total endogenous mRNA (spike-in regression)') + 
  nm_theme()
dev.off()

pdf('./supplementary_figures/ESC_algorithm_total_optim_m.pdf', width = 2, height = 2)
qplot( esApply(absolute_embryo_cds[, valid_cells], 2, sum), apply(ESC_norm_recovery_capture_rate$norm_cds[, ], 2, sum),log = 'xy') + 
  geom_smooth(method = 'rlm') + geom_abline(color = 'red') + expand_limits(x = c(1, 1e6), y = c(1, 1e06)) + xlab('Total endogenous mRNA (algorithm)') + ylab('Total endogenous mRNA (spike-in regression)') + 
  nm_theme()
dev.off()

ESC_mean <- apply(ESC_norm_recovery_capture_rate$norm_cds[,  ], 1, mean)
recovery_mean <- apply(absolute_embryo_cds[, valid_cells], 1, mean)

pdf('./supplementary_figures/ESC_algorithm_mean_optim_m.pdf', width = 2, height = 2)
qplot(ESC_mean[ESC_mean > 0 & recovery_mean > 0], recovery_mean[ESC_mean > 0 & recovery_mean > 0], log = 'xy') + 
  geom_abline(color = 'red') + 
  geom_smooth(method = 'rlm') + xlab('Mean endogenous mRNA count (regression)') + ylab('Mean endogenous mRNA count \n (algorithm)') + 
  nm_theme()
dev.off()

############################################################################################################################################################
#save the data

save.image('./RData/analysis_ESC_trilineage_every_spike.RData')
