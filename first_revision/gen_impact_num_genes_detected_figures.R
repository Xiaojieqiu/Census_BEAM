####################################
# Impact of number of genes detected on analysis
####################################
library(monocle)
library(ggplot2)
library(matrixStats)
library(grid)
library(gridExtra)
library(xacHelper)

# Helper function to get r squared as text between x and y
get_r_squared_text = function(y, x) {
    model = summary(lm(y ~ x))
    r_squared = format(round(model$r.squared, 4), nsmall = 2)
    paste("r-squared:", r_squared)
}

plot_smooth_scatter = function(x, y) {
    r_squared_text = get_r_squared_text(x, y)

    plotting_dataframe = data.frame(x=x, y=y)

    ggplot(plotting_dataframe, aes(x, log10(abs(y) + 0.1))) +
        stat_density2d(geom="tile", aes(fill=..density..^0.25, alpha=1), contour=FALSE) + 
        geom_point(size=0.5) +
        stat_density2d(geom="tile", aes(fill=..density..^0.25,     alpha=ifelse(..density..^0.25<0.4,0,1)), contour=FALSE) + 
        scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256)) + 
        annotate("text",  x=Inf, y = Inf, label = r_squared_text, vjust=1, hjust=1, size=4) +
        illustrator_theme()
}

load("RData/prepare_downsampling_data.RData")
load("RData/analysis_lung_data.RData")
load("RData/analysis_shalek_data.RData")

valid_genes_per_cell = exprs(original_depth_cds_transcript_counts) > 0 & exprs(original_depth_cds_transcript_counts_regression) > 0
valid_genes_per_cell[valid_genes_per_cell == FALSE] = NA
median_error_per_cell = colMedians(exprs(original_depth_cds_transcript_counts) * valid_genes_per_cell - exprs(original_depth_cds_transcript_counts_regression) * valid_genes_per_cell, na.rm=T)

number_of_genes_detected = colSums(exprs(original_depth_cds_transcript_counts_regression) > 0)
fpkm_mode = estimate_t(input_data$original_isoform_matrix)

# Plot of the mode of the FPKM distribution vs. number of genes detected for Lung data
pdf("number_of_genes_detected_impact_normalization_approach_a.pdf", width=2.5, height=2.5)
    plot_smooth_scatter(number_of_genes_detected, log10(abs(median_error_per_cell) + 0.1)) +
    xlab("number of genes detected") + 
    ylab("log10 median difference between regression\nand normalization transcript counts")
dev.off()

# Plot of error between transcript counts from regression and recovery vs. number of genes detected for Lung data
pdf("number_of_genes_detected_impact_normalization_approach_b.pdf", width=2.5, height=2.5)
    plot_smooth_scatter(number_of_genes_detected, log10(abs(fpkm_mode) + 0.1)) +
    xlab("number of genes detected") + 
    ylab("log10 FPKM mode per cell")
dev.off()

# Plot of the mode of the FPKM distribution vs. number of genes detected for Shalek data (only cells used for our analysis)
shalek_abs_cells_used = Shalek_abs[, pData(Shalek_abs)$experiment_name %in% c('Ifnar1_KO_LPS', 'Stat1_KO_LPS', "LPS_GolgiPlug", "LPS", "Unstimulated_Replicate")]
shalek_isoforms_cells_used = Shalek_isoform_fpkm_matrix[, pData(Shalek_abs)$experiment_name %in% c('Ifnar1_KO_LPS', 'Stat1_KO_LPS', "LPS_GolgiPlug", "LPS", "Unstimulated_Replicate")]

fpkm_mode_shalek = estimate_t(shalek_isoforms_cells_used)
number_of_genes_detected_shalek = colSums(exprs(shalek_abs_cells_used) > 0)

pdf("number_of_genes_detected_impact_normalization_approach_panel_c.pdf", width=2.5, height=2.5)
plot_smooth_scatter(number_of_genes_detected_shalek, log10(abs(fpkm_mode_shalek) + 0.1)) +
    xlab("number of genes detected") + 
    ylab("log10 FPKM mode per cell")
dev.off()
