library(argparse)
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
        return ( list(mean_value=median(proportion_detected, na.rm=T), lower_bound=quantile(proportion_detected, na.rm=T)[[2]], upper_bound=quantile(proportion_detected, na.rm=T)[[4]]))
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
total_expression = apply(exprs(original_depth_cds)[valid_genes, ], 1, max)
quartiles = cut(total_expression, breaks=quantile(total_expression, probs=seq(0,1, by=0.2)), include.lowest=TRUE)
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
performance$quartile_start = stringr::str_match(performance$quartile, "([0-9e\\-]+),")[, 2]
performance = performance %>% dplyr::arrange(as.numeric(quartile_start))
performance$quartile = factor(performance$quartile, levels=unique(performance$quartile))

# Make plot
pdf("figure1_genes_quantified_close_to_original_value_over_depth.pdf", height=8, width=2)

ggplot(performance, aes(as.numeric(depth), as.numeric(mean_value), color=quantification_type)) +
    geom_errorbar(aes(ymin=as.numeric(lower_bound), ymax=as.numeric(upper_bound))) +  
    geom_line() +
    geom_point() + 
    facet_grid(quartile ~ .) +
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
regression = melt(do.call(rbind, lapply(cds_to_compare_transcript_counts_regression, function(x) { rowMedians(exprs(x)) })))
colnames(regression) = c("depth", "gene", "transcript_counts_regression")

recovery = melt(do.call(rbind, lapply(cds_to_compare_transcript_counts, function(x) { rowMedians(exprs(x)) })))
colnames(recovery) = c("depth", "gene", "transcript_counts_recovery")

# Combine results into a single dataframe for plotting
combined_recovery_and_regression_by_depth = regression %>% left_join(recovery, by=c("depth", "gene"))

# Generate scatterplots at different depths to show how recovery algorithm performance changes with depth 
pdf("figure2_normalization_vs_spike_regression_by_depth.pdf", width=3, height=3)

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

pdf("figure3_median_error_between_recovery_methods_by_depth.pdf", width=1.75, height=1.75)

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

pdf("figure4_stability_of_m_and_c_with_depth.pdf", width=2.5, height=2.5)

ggplot(recovery_algorithm_stability_df, aes(depth, value, color=variable)) +
    geom_point() +
    geom_line() +
    xlab("max depth per cell (read pairs)") +
    ylab("spike-free method parameter value") +
    facet_grid(variable ~ ., scales="free_y") + 
    scale_color_brewer(palette="Set1") + 
    illustrator_theme()

dev.off()    


#########################################
# Make plots of constructed trees
#########################################
plot_monocle_spanning_tree_vectorized = function(cds) {
    monocle::plot_spanning_tree(cds, cell_size=1, color_by="Time")
        #scale_fill_manual(values=shalek_custom_color_scale)
}

# Figure 8: Pairwise correlation between the original trajectory and each downsampled one
plot_pseudotime_correlations <- function (original, downsampled)
{
    # Extract a dataframe with cell maturation times
    maturation_df <- data.frame(cell = c(colnames(original), colnames(downsampled)), 
        maturation_level = 100 * c(pData(original)$Pseudotime/max(pData(original)$Pseudotime),
                                   pData(downsampled)$Pseudotime/max(pData(downsampled)$Pseudotime)),
        Type = c(rep("Original", dim(original)[2]), rep("Sampled", dim(downsampled)[2])), rownames = colnames(downsampled))

    # Calculate a correlation coefficient between the two CDS maturation times
    cor.coeff <- cor(pData(downsampled)$Pseudotime, pData(original)$Pseudotime, method = "spearman")

    # Calculate rank differences of cells between methods
    maturation_df = maturation_df %>% group_by(Type) %>% mutate(rank=rank(maturation_level))
    maturation_df = maturation_df %>% group_by(cell) %>% mutate(rank_diff=abs(mean(rank) - rank) * 2)

    # Define custom color scheme
    highlighting = c('TRUE'='black', 'FALSE'="#f6f6f6")

    # Generate plot
    p <- ggplot(aes(x = maturation_level, y = Type, group = cell), data = maturation_df) + 
            geom_line(aes(color = rank_diff > 20), alpha=0.5, size=0.25) +
            geom_point(size = 1) + 
            scale_color_manual(values=c(highlighting)) + 
            xlab("Pseudotime") + 
            annotate("text", x = 80, y = 2.2, label = paste("corr:", round(cor.coeff, 3)), size=3) +
            illustrator_theme()

    return(p)
}

## Generate plot objects for transcript counts from normalization and regression
original_depth_cds_transcript_counts_trajectory = plot_monocle_spanning_tree_vectorized(original_depth_cds_transcript_counts_ordered)
cds_to_compare_transcript_counts_trajectories = lapply(cds_to_compare_transcript_counts_ordered, plot_monocle_spanning_tree_vectorized)

original_depth_cds_transcript_counts_regression_trajectory = plot_monocle_spanning_tree_vectorized(original_depth_cds_transcript_counts_regression_ordered)
cds_to_compare_transcript_counts_regression_trajectories = lapply(cds_to_compare_transcript_counts_regression_ordered, plot_monocle_spanning_tree_vectorized)


## Layout plots of trees obtained at each depth for each method
pdf("figure5_normalized_transcript_counts_tree_with_depth.pdf", width=6, height=5)
grob_args <- c(cds_to_compare_transcript_counts_trajectories, 2)
names(grob_args) <- c(names(cds_to_compare_transcript_counts_trajectories), "nrow")
plot_grob = do.call(grid.arrange, grob_args)
dev.off()

pdf("figure6_regression_transcript_counts_tree_with_depth.pdf", width=6, height=5)
grob_args <- c(cds_to_compare_transcript_counts_regression_trajectories, 2)
names(grob_args) <- c(names(cds_to_compare_transcript_counts_regression_trajectories), "nrow")
plot_grob = do.call(grid.arrange, grob_args)
dev.off()

# Now plot correlations between original pseudotime and downsampled Pseudotime and plot for each method
normalization_pairwise_correlation_plots = lapply(cds_to_compare_transcript_counts_ordered, function(cds) { plot_pseudotime_correlations(original_depth_cds_transcript_counts_ordered, cds) + illustrator_theme() })
regression_pairwise_correlation_plots = lapply(cds_to_compare_transcript_counts_regression_ordered, function(cds) { plot_pseudotime_correlations(original_depth_cds_transcript_counts_regression_ordered, cds) + illustrator_theme() })

pdf("figure7_normalized_transcript_counts_pseudotime_correlation.pdf", width=8, height=3)
grob_args <- c(normalization_pairwise_correlation_plots, 2)
names(grob_args) <- c(names(normalization_pairwise_correlation_plots), "nrow")
plot_grob = do.call(grid.arrange, grob_args)
dev.off()

pdf("figure8_regression_transcript_counts_pseudotime_correlation.pdf", width=8, height=3)
grob_args <- c(regression_pairwise_correlation_plots, 2)
names(grob_args) <- c(names(regression_pairwise_correlation_plots), "nrow")
plot_grob = do.call(grid.arrange, grob_args)
dev.off()