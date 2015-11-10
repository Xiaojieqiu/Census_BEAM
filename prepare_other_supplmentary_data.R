#all functions for the supplementary files: 

#Supplementary files:

#########################################################################################################
  #goodness of fit: 
  #generate the result of goodness of fit for each gene: 
  colnames(gd_fit_res_num)[1:2] <- c('NB', 'ZINB')
  test <- melt(gd_fit_res_num[, 1:3], id.vars = 'type')
  p1 <- qplot(as.factor(variable), as.numeric(value), geom = 'bar', stat = 'identity', data = test, fill = type) + facet_wrap('type') + nm_theme() + 
    theme(legend.position = 'none') + xlab('Fit types') + ylab('number of genes') + theme(strip.background = element_blank(),
         strip.text.x = element_blank()) + theme(axis.text.x = element_text(angle = 30, hjust = .9))
  p1 + xlab('')
  ggsave(paste(elife_directory, 'goodness_fit.pdf', sep = ''), height = 2, width = 3)
  
  #fig 3 SI: 
  quake_all_modes <- estimate_t(exprs(isoform_count_cds), return_all = T)

  cell_nanmes <- c("SRR1033974_0", "SRR1033922_0", "SRR1033866_0")
  cell_id <- which(colnames(isoform_count_cds) %in% cell_nanmes)
  three_cell_iso_df <- data.frame(Cell_id = rep(row.names(quake_all_modes)[cell_id], each = nrow(isoform_count_cds)), 
                  log10_FPKM = log10(c(exprs(isoform_count_cds)[, cell_id[1]], exprs(isoform_count_cds)[, cell_id[2]], exprs(isoform_count_cds)[, cell_id[3]])), 
                  Cell_mode = rep(log10(quake_all_modes[cell_id, 1]), each = nrow(isoform_count_cds)))

  three_cell_iso_df <- data.frame(Cell_id = rep(row.names(quake_all_modes)[which(quake_all_modes$best_cov_dmode <= 2)], each = nrow(isoform_count_cds)), 
                  log10_FPKM = log10(c(exprs(isoform_count_cds)[, which(quake_all_modes$best_cov_dmode <= 2)])), 
                  Cell_mode = rep(log10(quake_all_modes[which(quake_all_modes$best_cov_dmode <= 2), 1]), each = nrow(isoform_count_cds)))


  qplot(x = log10_FPKM, geom = 'histogram', data = three_cell_iso_df[, ], binwidth = .05, color = I('red'))  +
    geom_vline(aes(xintercept=log10(Cell_mode)), color = 'blue') + facet_wrap(~Cell_id) + xlim(-3, 5) + monocle_theme_opts() + xlab('log10 FPKM') + ylab('Isoform counts') + nm_theme()
  ggsave(filename = paste(elife_directory, 'SI/eLife_fig4_SI.pdf', sep = ''), width = 3, height = 1.2)

  10^mapply(function(cell_dmode, model) {
      predict(model, newdata = data.frame(log_fpkm = cell_dmode), type = 'response')
  }, as.list(unique(three_cell_iso_df$Cell_mode)), molModels_select[c(1,9,14)])


