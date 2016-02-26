library(tidyr)
library(grid)
library(gridExtra)
library(monocle)
library(DevTree)
library(plyr)
library(data.table)
library(MASS)
library(modeest)
library(dplyr)
library(matrixStats)


####################################################
# Load prepared downsampling data and required code
####################################################
load("RData/prepare_downsampling_data.RData")
source("monocle_helper_functions.R")

####################################################
# Helper functions for analysis
####################################################

# Helper function to take expression matrix, reference matrix, a threshold, and set of quartile assignments and return
# Proportion of genes quantified within error_threshold of the original value in reference matrix for each quartile of gene expression
calculate_fraction_of_genes_detected_within_error = function(expression_matrix, reference_matrix, error_threshold, quartiles) {
    expression_matrix = exprs(expression_matrix)[names(quartiles), ]
    reference_matrix = exprs(reference_matrix)[names(quartiles), ]

    # Split dataframe into quantiles
    expression_matrix_by_quartile = split(as.data.frame(expression_matrix), quartiles)
    reference_matrix_by_quartile = split(as.data.frame(reference_matrix), quartiles)

    get_performance = function(expression_matrix, reference_matrix) {
        error = abs(expression_matrix - reference_matrix) / (reference_matrix)

        proportion_detected = colSums(error < error_threshold, na.rm=T) / colSums(reference_matrix > 0)
        return ( list(mean_value=mean(proportion_detected, na.rm=T), lower_bound=quantile(proportion_detected, na.rm=T)[[2]], upper_bound=quantile(proportion_detected, na.rm=T)[[4]]))
    }

    performance_by_quartile = as.data.frame(t(mapply(get_performance, expression_matrix_by_quartile, reference_matrix_by_quartile)))
    performance_by_quartile$quartile = row.names(performance_by_quartile)
    row.names(performance_by_quartile) = NULL
    return(performance_by_quartile)
}

# Helper function to massage results obtained with lapply calls to calculate_fraction_of_genes_detected_within_error
# into a more usable format
process_performance_dataframe = function(performance_results, quantification_type) {
    performance_results = as.data.frame(apply(performance_results, 2, unlist))  # lapply strategy has columns as lists... need to undo
    performance_results$depth = as.numeric(stringr::str_replace(row.names(performance_results), '\\..+', ''))
    performance_results$quantification_type = quantification_type
    performance_results$quartile = performance_results$quartile
    return(performance_results)
}

####################################################
# Calculate CDF curves of proportions of genes quantified within
# threshold of error of original measurement for quantiles of gene 
# expression for each of the downsampled depths.
####################################################

# Define quantiles of gene expression based on full depth FPKM data
valid_genes = row.names(original_depth_cds[rowSums(exprs(original_depth_cds)) > 0, ])
total_expression = apply(exprs(original_depth_cds)[valid_genes, ], 1, mean)
quartiles = cut(total_expression, breaks=quantile(total_expression, probs=seq(0,1, by=0.25)), include.lowest=TRUE)
names(quartiles) = names(total_expression)

# Calculate performance per quartile
ERROR_THRESHOLD = 0.20
performance_transcript_counts = as.data.frame(do.call(rbind, lapply(cds_to_compare_transcript_counts, calculate_fraction_of_genes_detected_within_error, reference_matrix=original_depth_cds_transcript_counts, error_threshold=ERROR_THRESHOLD, quartiles=quartiles)))
performance_fpkm = as.data.frame(do.call(rbind, lapply(cds_to_compare, calculate_fraction_of_genes_detected_within_error, reference_matrix=original_depth_cds, error_threshold=ERROR_THRESHOLD, quartiles=quartiles)))
performance_transcript_counts_regression = as.data.frame(do.call(rbind, lapply(cds_to_compare_transcript_counts_regression, calculate_fraction_of_genes_detected_within_error, reference_matrix=original_depth_cds_transcript_counts_regression, error_threshold=ERROR_THRESHOLD, quartiles=quartiles)))

# Prepare dataframes for plotting
performance_fpkm = process_performance_dataframe(performance_fpkm, "FPKM")
performance_transcript_counts = process_performance_dataframe(performance_transcript_counts, "spike-free transcript counts")
performance_transcript_counts_regression = process_performance_dataframe(performance_transcript_counts_regression, "spike-in regression transcript counts")

# Generate a single matrix with methods (FPKM, normalization, and spike-in regression)
performance = rbind(performance_transcript_counts_regression, performance_transcript_counts, performance_fpkm)

# Reorder by quartile magnitude
performance$quartile_start = stringr::str_match(performance$quartile, "([0-9e\\-\\.]+),")[, 2]
performance = performance %>% dplyr::arrange(as.numeric(quartile_start))
performance$quartile = factor(performance$quartile, levels=unique(performance$quartile))

# Make plot
pdf("fig7c_si.pdf", height=3, width=3)

ggplot(performance, aes(as.numeric(depth), as.numeric(mean_value), color=quantification_type)) +
    geom_point() +
    geom_line() +
    geom_hline(y=1, linetype="longdash", color="black", alpha=0.75) + 
    ylim(c(0, 1)) +
    facet_wrap(~ quartile, nrow=2) +
    illustrator_theme() +
    xlab("max depth per cell (read pairs)") +
    ylab("median proportion of genes quantified within 20% of original per cell") +
    scale_color_brewer(palette="Set1") +
    scale_fill_brewer(palette="Set1")
dev.off()
    

####################################################
# Quantify how normalization-based recovery algorithm performance
# changes with depth
####################################################

# Get median transcript count value across all cells for both regression and recovery algorithm
regression = melt(do.call(rbind, lapply(cds_to_compare_transcript_counts_regression, function(x) { rowMedians(exprs(x)[valid_genes, ]) })))
colnames(regression) = c("depth", "gene", "transcript_counts_regression")

recovery = melt(do.call(rbind, lapply(cds_to_compare_transcript_counts, function(x) { rowMedians(exprs(x)[valid_genes, ]) })))
colnames(recovery) = c("depth", "gene", "transcript_counts_recovery")

# Combine results into a single dataframe for plotting
combined_recovery_and_regression_by_depth = regression %>% left_join(recovery, by=c("depth", "gene"))

# Generate scatterplots at different depths to show how recovery algorithm performance changes with depth 
pdf("fig7a_si.pdf", width=3, height=3)

ggplot(combined_recovery_and_regression_by_depth, aes(log10(transcript_counts_regression + 0.1), log10(transcript_counts_recovery + 0.1))) +
    geom_point(size=0.3) + 
    geom_abline(color="red", alpha=0.75) + 
    facet_wrap(~ depth, ncol=3) + 
    xlab("log10(spike-in regression transcript counts + 0.1)") + 
    ylab("log10(spike-free transcript counts + 0.1)") + 
    illustrator_theme()

dev.off()

####################################################
# Quantify median error in regression vs. normalization
# plots from figure 2
####################################################
## Filter out transcripts not detected in either the spike-in regression or normalization
combined_recovery_and_regression_by_depth = subset(combined_recovery_and_regression_by_depth, ! (transcript_counts_regression == 0 & transcript_counts_recovery == 0))

median_proportion_of_error_per_depth = combined_recovery_and_regression_by_depth %>% group_by(depth) %>% summarize(median=median((transcript_counts_recovery - transcript_counts_regression) / sum(transcript_counts_regression)), quartile_2=quantile((transcript_counts_recovery - transcript_counts_regression) / sum(transcript_counts_regression))[[2]], quartile_3=quantile((transcript_counts_recovery - transcript_counts_regression) / sum(transcript_counts_regression))[[4]])

pdf("fig7b_si.pdf", width=3, height=3)
ggplot(median_proportion_of_error_per_depth, aes(depth, median)) +
    geom_hline(y=0, linetype="longdash", color="red", alpha=0.75) + 
    geom_errorbar(aes(ymax=quartile_3, ymin=quartile_2), color="#d3d3d3") +
    geom_point(size=1) +
    geom_line() +
    xlab("max depth per cell (read pairs)") +
    ylab("median normalized difference between recovery methods") + 
    illustrator_theme()
dev.off()


####################################################
# Verify that m and c parameters in normalization
# approach are relatively robust to depth
####################################################
# Plot of M and C parameter stability
recovery_algorithm_stability_df = data.frame(m=unlist(cds_to_compare_conversion_m_values), c=unlist(cds_to_compare_conversion_c_values), depth=as.numeric(names(cds_to_compare_conversion_c_values)))
recovery_algorithm_stability_df = reshape2::melt(recovery_algorithm_stability_df, id="depth")

# Temporary hack to fix depth of original dataset (not 10,000,000)
recovery_algorithm_stability_df$depth[recovery_algorithm_stability_df$depth == 10000000] = 5000000

pdf("fig7d_si.pdf", width=2.5, height=2.5)

ggplot(recovery_algorithm_stability_df, aes(depth, value, color=variable)) +
    geom_point() +
    geom_line() +
    xlab("max depth per cell (read pairs)") +
    ylab("spike-free method parameter value") +
    facet_grid(variable ~ ., scales="free_y") + 
    scale_color_brewer(palette="Set1") + 
    illustrator_theme()

dev.off()


####################################################
# Examine the mode of the transcript counts and log10(FPKM) distribution
# over the dowsampled datasets  
####################################################
get_proportion_modes_correct_per_depth = function(cds_list, original_depth_cds, error_threshold = 0.2) {
    # Calculate modes for each depth and join with metadata
    modes = as.data.frame(do.call(rbind, lapply(cds_list, function(x) estimate_t(x))))
    modes$depth = as.numeric(row.names(modes))
    modes_melted = melt(modes, id="depth")
    modes_melted = merge(modes_melted, pData(original_depth_cds), by.x="variable", by.y="row.names")

    # Calculate proportion of cells with mode within 20% of original at each depth
    proportion_correct = modes_melted %>%
        group_by(variable) %>%
            mutate(percent_error=(value - value[depth == max(depth)]) / value[depth == max(depth)]) %>%
            ungroup() %>%
        group_by(depth) %>%
            summarize(proportion_correct = sum(percent_error <= error_threshold) / n())

    return(proportion_correct)
}

names(input_data$isoform_matrices_to_compare) = downsampled_depths

# Calculate proportion of cells with mode within 20% of the original value
proportion_correct_fpkm = get_proportion_modes_correct_per_depth(input_data$isoform_matrices_to_compare, original_depth_cds)
proportion_correct_fpkm$method = "fpkm"

proportion_correct_transcript_counts_regression = get_proportion_modes_correct_per_depth(cds_to_compare_transcript_counts_regression, original_depth_cds)
proportion_correct_transcript_counts_regression$method = "spike-in regression"

proportion_correct_transcript_counts = get_proportion_modes_correct_per_depth(cds_to_compare_transcript_counts, original_depth_cds)
proportion_correct_transcript_counts$method = "spike-in free"

# Combine all results and generate plot
proportion_correct = rbind(proportion_correct_fpkm, proportion_correct_transcript_counts_regression, proportion_correct_transcript_counts)

ggplot(proportion_correct, aes(depth, proportion_correct, color=method)) +
    geom_point() +
    geom_line() +
    ylim(c(0, 1)) + 
    scale_color_brewer(palette="Set1") +
    ylab("proportion of cells within 20% of original value") +
    xlab("max depth per cell (total aligned PE reads)") +
    theme_bw() +
    ggsave("temp.png", height=7.5, width=7.5)