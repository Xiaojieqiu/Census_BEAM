#abs_mat_1$norm_cds <- abs_mat_1_new_default$norm_cds

new_endogenous_RNA <- colSums(abs_mat_1$norm_cds[1:transcript_num, ])
outlier <- which(new_endogenous_RNA == max(new_endogenous_RNA))
new_endogenous_RNA <- new_endogenous_RNA#[-outlier]

absolute_cds_subset <- absolute_cds#[, -outlier]
pdf('./main_figures/fig3f_new.pdf', width = 2, height = 1.7)
pData(absolute_cds_subset)$endogenous_RNA <- esApply(absolute_cds_subset, 2, function(x) sum(x[1:transcript_num]))
qplot(pData(absolute_cds_subset)$endogenous_RNA,
      new_endogenous_RNA, log="xy", color=pData(absolute_cds_subset)$Time, size = I(4)) +
  geom_smooth(method="lm", color="black", size = .1) + geom_abline(color="red") +
  xlab("Total endogenous mRNA \n (spike-in)") +
  ylab("Total endogenous mRNA \n (spike-in free algorithm)") + #scale_size(range = c(0.25, 0.25)) +
  scale_color_discrete(name = "Time points") + nm_theme()
dev.off()

valid_Time <- pData(absolute_cds_subset)$Time
E14.5_cell <- rowMeans(exprs(absolute_cds_subset[, which(valid_Time == 'E14.5')])) #1:transcript_num
E16.5_cell <- rowMeans(exprs(absolute_cds_subset[, which(valid_Time == 'E16.5')]))
E18.5_cell <- rowMeans(exprs(absolute_cds_subset[, which(valid_Time == 'E18.5')]))
Adult_cell <- rowMeans(exprs(absolute_cds_subset[, which(valid_Time == 'Adult')]))
norm_cds <- abs_mat_1$norm_cds#[, -outlier]
mc_E14.5_cell <- rowMeans(norm_cds[, colnames(absolute_cds_subset)[which(valid_Time == 'E14.5')]])
# mc_E14.5_cell <- rowMeans(abs_mat_time_14$norm_cds[, colnames(abs_AT12_cds_subset_all_gene)[which(Time == 'E14.5')]])http://127.0.0.1:23551/graphics/plot_zoom_png?width=562&height=202
mc_E16.5_cell <- rowMeans(norm_cds[, colnames(absolute_cds_subset)[which(valid_Time == 'E16.5')]])
#mc_E16.5_cell <- rowMeans(abs_mat_time_16$norm_cds[, colnames(abs_AT12_cds_subset_all_gene)[which(Time == 'E16.5')]])
mc_E18.5_cell <- rowMeans(norm_cds[, colnames(absolute_cds_subset)[which(valid_Time == 'E18.5')]])
mc_Adult_cell <- rowMeans(norm_cds[, colnames(absolute_cds_subset)[which(valid_Time == 'Adult')]])
#mc_Adult_cell <- rowMeans(abs_mat_time_adult$norm_cds[, colnames(abs_AT12_cds_subset_all_gene)[which(Time == 'Adult')]])
mc_abs_exprs_df <- data.frame(spikein = as.vector(c(E14.5_cell, E16.5_cell, E18.5_cell, Adult_cell)),
                              mc_algorithm = as.vector(c(mc_E14.5_cell, mc_E16.5_cell, mc_E18.5_cell, mc_Adult_cell)),
                              cell = rep(c("E14.5_cell", "E16.5_cell", "E18.5_cell", "Adult_cell"), each = nrow(mc_adj_cds[, ])))

pdf('./main_figures/fig3g2_new.pdf', width = 3, height = 2)
qplot(spikein + 1, mc_algorithm + 1, log = 'xy',
      color = cell, data = mc_abs_exprs_df, size  = 1.5) + facet_wrap(~cell, scales = 'free', ncol = 2) + #  geom_smooth(method = 'rlm', aes(group = 199), size = .1) +
  scale_size(range = c(1.5, 1)) +  geom_abline() + xlab('Transcript counts (Spike-in)') + scale_size(range = c(0.25, 0.25)) +
  ylab('Transcript counts (Recovery algorithm)') + nm_theme()
dev.off()

