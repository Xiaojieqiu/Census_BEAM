#use mse, correlation: 
library(tidyr)
library(grid)
library(gridExtra)
library(monocle) 
# library(devtools)
# load_all('~/Projects/monocle-dev')
library(xacHelper)
library(plyr)
library(data.table)
library(MASS)
library(modeest)
library(dplyr)
library(matrixStats)


####################################################
# Load prepared downsampling data and required code
####################################################
load("RData/prepare_downsampling_data_cell.RData")
# source("monocle_helper_functions.R")

####################################################
# Helper functions for analysis
####################################################

# New metrics for showing the performance of the recovery algorithm: 
calculate_mse_cor <- function(expression_cds, reference_cds) {
    expression_matrix <- exprs(expression_cds)[1:transcript_num, ]
    reference_matrix <- exprs(reference_cds)[1:transcript_num, colnames(expression_matrix)]
    mse <- mean((expression_matrix - reference_matrix)^2)
    cor_true_estimated <- mean(unlist(lapply(1:ncol(expression_matrix), function(x) cor(expression_matrix[, x], reference_matrix[, x]))))
    total_difference <- mean(abs(colSums(expression_matrix) - colSums(reference_matrix)))

    return(data.frame(mse = mse, cor = cor_true_estimated, total_difference = total_difference))
}

####################################################
# Calculate CDF curves of proportions of genes quantified within
# threshold of error of original measurement for quantiles of gene 
# expression for each of the downsampled depths.
####################################################

# Prepare dataframes for plotting
performance_transcript_counts = do.call(rbind.data.frame, lapply(cds_to_compare_transcript_counts, calculate_mse_cor, absolute_cds))
performance_transcript_counts$fraction <- names(cds_to_compare_transcript_counts)
performance_transcript_counts$Type <- "Census"
performance_transcript_counts_regression = do.call(rbind.data.frame, lapply(cds_to_compare_transcript_counts_regression, calculate_mse_cor, absolute_cds))
performance_transcript_counts_regression$fraction <- names(cds_to_compare_transcript_counts_regression)
performance_transcript_counts_regression$Type <- "spike-in regression"

# Combine all results and generate plot
performance_df = rbind(performance_transcript_counts, performance_transcript_counts_regression)

pdf("./supplementary_figures/censusperformance_rmse.pdf", width=1.5, height=1.5)
ggplot(performance_transcript_counts, aes(fraction, mse + 1, color=Type)) + 
    geom_point() +
    scale_y_log10() + 
    scale_color_brewer(palette="Set1") +
    ylab("proportion of cells sampled") +
    xlab("mse of expression to the full dataset") +
    nm_theme() 
dev.off()

pdf("./supplementary_figures/performance_rmse.pdf", width=1.5, height=1.5)
ggplot(performance_df, aes(fraction, mse + 1, color=Type)) + 
    geom_point() +
    scale_y_log10() + 
    scale_color_brewer(palette="Set1") +
    ylab("proportion of cells sampled") +
    xlab("mse of expression to the full dataset") +
    nm_theme() 
dev.off()

pdf("./supplementary_figures/performance_cor.pdf", width=1.7, height=1.5)
ggplot(performance_df, aes(fraction, cor, color=Type)) +
    geom_point() +
    geom_line() +
    ylim(c(0, 1)) + 
    scale_color_brewer(palette="Set1") +
    ylab("proportion of cells sampled") +
    xlab("correlation of expression to the full dataset") +
    theme(axis.text.x = element_text(angle = 30)) + 
    nm_theme() 
dev.off()

pdf("./supplementary_figures/performance_total_difference.pdf", width=1.5, height=1.5)
ggplot(performance_df, aes(fraction, total_difference, color=Type)) +
    geom_point() +
    scale_y_log10() + 
    scale_color_brewer(palette="Set1") +
    ylab("proportion of cells sampled") +
    xlab("mean difference of total mRNA to the full dataset") +
    nm_theme() 
dev.off()
