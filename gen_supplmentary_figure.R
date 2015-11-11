  library(monocle)
  library(xacHelper)

  load_all_libraries()

  load('analysis_other_supplementary_data.RData')
  load('analysis_lung_data.RData')

  # figure b: 
  #concordance between alernative analysis improves when absolute copy numbers are used: 
  #(overlapping plot between absolute copy / read counts)
  #read counts: 
  default_deseq_p[is.na(default_deseq_p)] <- 1
  default_deseq2_p[is.na(default_deseq2_p)] <- 1
  default_edgeR_p[is.na(default_edgeR_p)] <- 1

  element_all <- c(
                  names(default_edgeR_p[default_edgeR_p < 0.01]), 
                  names(default_deseq2_p[default_deseq2_p < 0.01]), 
                  names(readcount_permutate_pval[which(readcount_permutate_pval < .01)]), 
                  names(default_deseq_p[default_deseq_p < 0.01]), 
                  names(monocle_p_readcount[monocle_p_readcount < 0.01]), 
                  names(scde_p[scde_p < 0.01]))
  sets_all <- c(
                rep(paste('edgeR', sep = ''), sum(default_edgeR_p < 0.01, na.rm = T)), 
                rep(paste('DESeq2', sep = ''), sum(default_deseq2_p < 0.01, na.rm = T)), 

                rep(paste('Permutation test', sep = ''), length(which(readcount_permutate_pval < .01))),
                rep(paste('DESeq', sep = ''), length(default_deseq_p[default_deseq_p < 0.01])), 
                rep(paste('Monocle', sep = ''), length(monocle_p_readcount[monocle_p_readcount < 0.01])), 
                rep(paste('SCDE', sep = ''), length(scde_p[scde_p < 0.01])))

  pdf(paste(elife_directory, 'eLife_fig2c.1.pdf', sep = ''))
  venneuler_venn(element_all, sets_all)
  table(sets_all)

  dev.off()

  #absolute transcript counts: 
  abs_default_deseq_p[is.na(abs_default_deseq_p)] <- 1
  abs_default_deseq2_p[is.na(abs_default_deseq2_p)] <- 1
  abs_default_edgeR_p[is.na(abs_default_edgeR_p)] <- 1

  abs_element_all <- c(
                   names(abs_default_edgeR_p[abs_default_edgeR_p < 0.01]), 
                   names(abs_default_deseq2_p[abs_default_deseq2_p < 0.01]), 
                   names(new_abs_size_norm_monocle_p_ratio_by_geometric_mean[which(new_abs_size_norm_monocle_p_ratio_by_geometric_mean < .01)]), 
                   names(abs_default_deseq_p[abs_default_deseq_p < 0.01]), 
                   names(abs_scde_p[abs_scde_p < 0.01]), 
                   names(mode_size_norm_permutate_ratio_by_geometric_mean[which(mode_size_norm_permutate_ratio_by_geometric_mean < 0.01)]))

  abs_sets_all <- c(
                rep(paste('edgeR', sep = ''), sum(abs_default_edgeR_p < 0.01, na.rm = T)), 
                rep(paste('DESeq2', sep = ''), sum(abs_default_deseq2_p < 0.01, na.rm = T)), 
                rep(paste('Monocle', sep = ''), length(new_abs_size_norm_monocle_p_ratio_by_geometric_mean[new_abs_size_norm_monocle_p_ratio_by_geometric_mean < 0.01])), 
                rep(paste('DESeq', sep = ''), sum(abs_default_deseq_p < 0.01, na.rm = T)), 
                rep(paste('SCDE', sep = ''), length(abs_scde_p[abs_scde_p < 0.01])),
                rep(paste('Permutation test', sep = ''), length(which(mode_size_norm_permutate_ratio_by_geometric_mean < 0.01))))

  pdf(paste(elife_directory, 'eLife_fig2c.2.pdf', sep = ''))
  venneuler_venn(abs_element_all, abs_sets_all)
  dev.off()

  table(sets_all)

  #add the barplot for the overlapping genes: 
  element_all_list <- list(
                  names(default_edgeR_p[default_edgeR_p < 0.01]), 
                  names(default_deseq2_p[default_deseq2_p < 0.01]), 
                  names(readcount_permutate_pval[which(readcount_permutate_pval < .01)]), 
                  names(default_deseq_p[default_deseq_p < 0.01]), 
                  names(monocle_p_readcount[monocle_p_readcount < 0.01]))

  abs_element_all_list <- list(
                   names(abs_default_edgeR_p[abs_default_edgeR_p < 0.01]), 
                   names(abs_default_deseq2_p[abs_default_deseq2_p < 0.01]), 
                   names(new_abs_size_norm_monocle_p_ratio_by_geometric_mean[which(new_abs_size_norm_monocle_p_ratio_by_geometric_mean < .01)]), 
                   names(abs_default_deseq_p[abs_default_deseq_p < 0.01]), 
                   names(mode_size_norm_permutate_ratio_by_geometric_mean[which(mode_size_norm_permutate_ratio_by_geometric_mean < 0.01)]))

  readcount_overlap <- Reduce(intersect, element_all_list)
  readcount_union <- Reduce(union, element_all_list)
  abs_overlap <- Reduce(intersect, abs_element_all_list)
  abs_union <- Reduce(union, abs_element_all_list)
  
  overlap_df <- data.frame(read_counts = length(readcount_overlap), transcript_counts = length(abs_overlap)) 
  union_df <- data.frame(read_counts = length(readcount_union), transcript_counts = length(abs_union))  

  qplot(variable, value, geom = 'bar', stat = 'identity', fill = variable, data = melt(overlap_df)) + xlab('') + ylab('number') + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = .9))
  ggsave(paste(elife_directory, '/SI/eLife_fig_SI_DEG_overlapping.pdf', sep = ''), width = 1, height = 1.1)

  qplot(variable, value, geom = 'bar', stat = 'identity', fill = variable, data = melt(union_df)) + xlab('') + ylab('number') + nm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = .9))
  ggsave(paste(elife_directory, '/SI/eLife_fig_SI_DEG_union.pdf', sep = ''), width = 1, height = 1.1)


  #

  #test the cell cycle: 
  #cyclin E, CDK2, Cyclin A, CDK1, CDK2, Cyclin B, CDK1
  #CCNE1, CCNE2; 
  #CDK2
  #CCNA1, CCNA2
  #CCNB1, CCNB2
    
  cc_markers <- which(fData(abs_AT12_cds_subset_all_gene)$gene_short_name %in% c('Ccne1', 'Ccne2', 'Cdk2', 'Ccna1', 'Ccna2', 'Ccnb1', 'Ccnb2', 'Cdk1'))
  colour_cell <- rep(0, length(new_cds$Lineage))
  names(colour_cell) <- as.character(new_cds$Time)
  colour_cell[names(colour_cell) == 'E14.5'] <- "#7CAF42"
  colour_cell[names(colour_cell) == 'E16.5'] <- "#00BCC3"
  colour_cell[names(colour_cell) == 'E18.5'] <- "#A680B9"
  colour_cell[names(colour_cell) == 'Adult'] <- "#F3756C"

  colour <- rep(0, length(new_cds$Lineage))
  names(colour) <- as.character(new_cds$Lineage)
  colour[names(colour) == 'AT1'] <- AT1_Lineage
  colour[names(colour) ==  'AT2'] <- AT2_Lineage

  plot_genes_branched_pseudotime2(abs_AT12_cds_subset_all_gene[cc_markers, ], cell_color_by = "Time", trajectory_color_by = "Lineage", fullModelFormulaStr = '~sm.ns(Pseudotime, df = 3)*Lineage', normalize = F, stretch = T, lineage_labels = c('AT1', 'AT2'), cell_size = 1, ncol = 2) + ylab('Transcript counts') ##nm_theme()
  ggsave(paste(elife_directory, 'Cell_cycle.pdf', sep = ''), height = 2, width = 3)

  #Shalek. show all the cells on the same graph:  
  Shalek_abs_subset <- Shalek_abs[,pData(Shalek_abs)$experiment_name %in% c('Ifnar1_KO_LPS', 'Stat1_KO_LPS',  "LPS", "LPS_GolgiPlug", "Unstimulated_Replicate")]

  pData(Shalek_abs_subset)$Total_mRNAs <- colSums(exprs(Shalek_abs_subset))
  Shalek_abs_subset <- Shalek_abs_subset[, pData(Shalek_abs_subset)$Total_mRNAs < 75000]

  DEG_union <- c(row.names(subset(Shalek_LPS_subset_DEG_res, qval < qval_thrsld)))

  Shalek_abs_subset <- setOrderingFilter(Shalek_abs_subset, DEG_union)
  Shalek_abs_subset <- reduceDimension(Shalek_abs_subset, use_vst = T, use_irlba=F, pseudo_expr = 0, covariates = as.vector(pData(Shalek_abs_subset)$num_genes_expressed) )
  Shalek_abs_subset <- orderCells(Shalek_abs_subset, num_path = 5)

  monocle::plot_spanning_tree(Shalek_abs_subset, color_by="interaction(experiment_name, time)", cell_size=5) + 
    scale_color_manual(values=shalek_custom_color_scale_plus_states) + 
      ggsave('figure_SI_all_cells_tree.pdf', , height=12, width = 7)

  plot_spanning_tree(Shalek_abs_subset, color_by="interaction(experiment_name, time)", cell_size=5) + #x = 1, y = 2, 
      scale_color_manual(values=shalek_custom_color_scale_plus_states) + 
      ggsave(filename = ggplot_name, height = 20, width = 20) 

#Explaining why E16.5d cell has bad correspondence: 
#fraction of clusters in each biotype: 
valid_class <- c("lincRNA", "processed_transcript", "protein_coding", "pseudogene", 'spike', "rRNA")
gene_class <- fData(read_countdata_cds)$biotype

read_countdata_cds_biotype <- apply(read_countdata_cds, 2, function(x) {
    sum_x <- sum(x)
    c(lincRNA = sum(x[gene_class == valid_class[1]]) / sum_x, 
        processed_transcript = sum(x[gene_class == valid_class[2]]) / sum_x, 
        protein_coding = sum(x[gene_class == valid_class[3]]) / sum_x, 
        pseudogene = sum(x[gene_class == valid_class[4]]) / sum_x, 
        MT_RNA = sum(x[gene_class %in% c('Mt_rRNA', 'Mt_tRNA')]) / sum_x, 
        spike = sum(x[gene_class == valid_class[5]]) / sum_x,         
        rRNA = sum(x[gene_class == valid_class[6]]) / sum_x,
        others = sum(x[!(gene_class %in% c(valid_class, 'Mt_rRNA', 'Mt_tRNA'))]) / sum_x
        ) 
    })
mlt_read_countdata_cds_biotype <- melt(read_countdata_cds_biotype)
mlt_read_countdata_cds_biotype$Time <- pData(standard_cds)[mlt_read_countdata_cds_biotype$Var2, 'Time']
qplot(Var2, value, geom = 'histogram', fill = Var1, data = mlt_read_countdata_cds_biotype, stat = 'identity', group = Time) + facet_wrap(~Time, scale  = 'free_x') + monocle_theme_opts()
ggsave(filename = paste(elife_directory, 'gene_type_percentage.pdf', sep = ''), width = 2, height = 3)

#number of read counts for spikein data: 
df <- data.frame(Time = pData(read_countdata_cds)$Time, sum_readcounts = esApply(read_countdata_cds[fData(read_countdata_cds)$biotype == 'spike', ], 2, sum))
# qplot(sum_readcounts, fill = Time, log = 'x') + facet_wrap(~Time)
qplot(sum_readcounts, fill = Time, log = 'x', data = df) + facet_wrap(~Time, ncol = 1, scales = 'free_y') + nm_theme() 
ggsave(filename = paste(elife_directory, '/SI/read_countdata_cds_sum_spikein.pdf', sep = ''), width = 2, height = 3)

#number of ERCC spike-in detected in each cell
ercc_controls_detected_df <- data.frame(loss = esApply(ercc_controls, 2, function(x) sum(x > 0)), Time = pData(absolute_cds[, colnames(loss_ercc_spikein)])$Time)
qplot(loss, fill = Time, data = ercc_controls_detected_df) + facet_wrap(~Time, ncol = 1) + nm_theme()
ggsave(filename = paste(elife_directory, '/SI/spikein_detected.pdf', sep = ''), width = 2, height = 3)

#readcount for the Shalek data: The Shalek data is great
#read the read count data for the genes: 
dir = "/net/trapnell/vol1/ajh24/proj/2015shalek_et_al_reanalysis/results/ahill/2015_05_07_input_files_for_monocle/cuffnorm_output_files/"
Shalek_sample_table <- read.delim(paste(dir, "/samples.table", sep = ''))
Shalek_norm_count <- read.delim(paste(dir, "/genes.count_table", sep = ''))
row.names(Shalek_norm_count) <- Shalek_norm_count$tracking_id
Shalek_norm_count <- Shalek_norm_count[, -1]

Shalek_read_countdata <- round(t(t(Shalek_norm_count) * Shalek_sample_table$internal_scale)) #convert back to the raw counts 
Shalek_read_countdata <- Shalek_read_countdata[row.names(Shalek_abs), paste(colnames(Shalek_abs), '_0', sep = '')]
colnames(Shalek_read_countdata) <- colnames(Shalek_abs)

Shalek_read_countdata_cds <- newCellDataSet(as.matrix(Shalek_read_countdata),
                       phenoData = new("AnnotatedDataFrame", data = pData(Shalek_abs)),
                       featureData = new("AnnotatedDataFrame", data = fData(Shalek_abs)),
                       expressionFamily = negbinomial(),
                       lowerDetectionLimit = 1)

pData(Shalek_read_countdata_cds)$Total_mRNAs <- esApply(Shalek_read_countdata_cds, 2, sum)
pData(Shalek_read_countdata_cds)$endogenous_RNA <- esApply(Shalek_read_countdata_cds, 2, sum)

Shalek_gene_df <- data.frame(experiment_name = pData(Shalek_read_countdata_cds[, c(colnames(Shalek_abs_subset_ko_LPS), colnames(Shalek_golgi_update))])$experiment_name, 
                         sum_readcounts = esApply(Shalek_read_countdata_cds[, c(colnames(Shalek_abs_subset_ko_LPS), colnames(Shalek_golgi_update))], 2, sum))

qplot(sum_readcounts, fill = experiment_name, log = 'x', data = Shalek_gene_df) + facet_wrap(~Time, ncol = 1, scales = 'free_y') + nm_theme() 
ggsave(filename = paste(elife_directory, 'Shalek_readcounts.pdf', sep = ''), width = 2, height = 3)

#   gene_names <- row.names(abs_branchTest_res_stretch[(abs_branchTest_res_stretch$qval < 0.05 & !is.na(abs_branchTest_res_stretch$ABCs)), ])[1:2881]
  

  #fit of distributions 
  #test this: 
  # abs_gd_fit_res <- mcesApply(absolute_cds[ ], 1, gd_fit_pval, cores = detectCores(), required_packages = c('VGAM', 'fitdistrplus', 'MASS', 'pscl'), exprs_thrsld = 10, pseudo_cnt = 0.01)
  # closeAllConnections()
  # std_gd_fit_res <- mcesApply(standard_cds[, ], 1, gd_fit_pval, cores = detectCores(), required_packages = c('VGAM', 'fitdistrplus', 'MASS', 'pscl'), exprs_thrsld = 10, pseudo_cnt = 0.01)
  # closeAllConnections()
  # tpm_gd_fit_res <- mcesApply(TPM_cds[, ], 1, gd_fit_pval, cores = detectCores(), required_packages = c('VGAM', 'fitdistrplus', 'MASS', 'pscl'), exprs_thrsld = 10, pseudo_cnt = 0.01)
  # closeAllConnections()
  # read_gd_fit_res <- mcesApply(count_cds[, ], 1, gd_fit_pval, cores = detectCores(), required_packages = c('VGAM', 'fitdistrplus', 'MASS', 'pscl'), exprs_thrsld = 10, pseudo_cnt = 0.01)

  # abs_gd_fit_res <- unlist(abs_gd_fit_res) 
  # read_gd_fit_res <- unlist(abs_gd_fit_res) 

  # abs_gd_fit_df <- matrix(abs_gd_fit_res, nrow(absolute_cds), ncol = 11, byrow = T)
  # dimnames(abs_gd_fit_df) <- list(row.names(absolute_cds), c("ln_pvalue", "nb_pvalue", "ln_pvalue.glm.link", "ln_pvalue.glm.log", "ln_pvalue.chisq", "nb_pvalue.glm", "nb_pvalue.chisq", "zinb_pvalue.chisq", "zanb_pvalue.chisq", "zinb_pvalue", "zanb_pvalue"))
  # read_gd_fit_df <- matrix(read_gd_fit_res, nrow(absolute_cds), ncol = 11, byrow = T)
  # dimnames(read_gd_fit_df) <- list(row.names(absolute_cds), c("ln_pvalue", "nb_pvalue", "ln_pvalue.glm.link", "ln_pvalue.glm.log", "ln_pvalue.chisq", "nb_pvalue.glm", "nb_pvalue.chisq", "zinb_pvalue.chisq", "zanb_pvalue.chisq", "zinb_pvalue", "zanb_pvalue"))

  # # 
  #select only nb and zinb and calculate the number of genes pass goodness of fit and number of genes can be fitted: 
  abs_gd_fit_res <- cal_gd_statistics(abs_gd_fit_df[, c('nb_pvalue', 'zinb_pvalue')], percentage = F, type = 'absolute', gene_list = valid_genes)
  readcount_gd_fit_res <- cal_gd_statistics(read_gd_fit_df[, c('nb_pvalue', 'zinb_pvalue')], percentage = F,  type = 'readcount', gene_list = valid_genes)
  gd_fit_res <- rbind(abs_gd_fit_res, readcount_gd_fit_res)
  gd_fit_res <- cbind(gd_fit_res, data_type = row.names(gd_fit_res))
  row.names(gd_fit_res) <- NULL
  gd_fit_res <- as.data.frame(gd_fit_res)
  
  gd_fit_res_num <- subset(gd_fit_res, data_type == 'gd_fit_num')
  gd_fit_res_success_num <- subset(gd_fit_res, data_type == 'success_fit_num')
  # 
  
    #generate the result of goodness of fit for each gene: 
  colnames(gd_fit_res_num)[1:2] <- c('NB', 'ZINB')
  test <- melt(gd_fit_res_num[, 1:3], id.vars = 'type')
  p1 <- qplot(as.factor(variable), as.numeric(value), geom = 'bar', stat = 'identity', data = test, fill = type) + facet_wrap('type') + nm_theme() + 
    theme(legend.position = 'none') + xlab('Fit types') + ylab('number of genes') + theme(strip.background = element_blank(),
         strip.text.x = element_blank()) + theme(axis.text.x = element_text(angle = 30, hjust = .9))
  p1 + xlab('')
  ggsave(paste(elife_directory, 'goodness_fit.pdf', sep = ''), height = 1.5, width = 1)

  colnames(gd_fit_res_success_num)[1:2] <- c('NB', 'ZINB')
  test <- melt(gd_fit_res_success_num[, 1:3], id.vars = 'type')

  p2 <- qplot(as.factor(variable), as.numeric(value), geom = 'bar', stat = 'identity', data = test, fill = type) + facet_wrap('type') + nm_theme() + 
     theme(legend.position = 'none') + xlab('Fit types') + ylab('number of genes') + theme(strip.background = element_blank(),
          strip.text.x = element_blank()) + theme(axis.text.x = element_text(angle = 30, hjust = .9))
  p2 + xlab('')
  ggsave(paste(elife_directory, 'goodness_fit2.pdf', sep = ''), height = 1.5, width = 1)

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

  #################### generate the figures for FigSC6: ###################
  
  #test on three other datasets for the differential gene expression: 
  #quake new data: 
  #NBt data: 
  #molecular cell data: 
  #several UMI data (use two Drop-seq): 
    #test mode, and the regression relationship between FPKM and UMI dataset (are the k/b also on a line)
    #differential gene expression test and comparing fitting of NB
    
  #check the influence of spike-in free estimation: 
  spike_free_standard_cds <- exprs(standard_cds)[1:transcript_num, ]

  pd <- new("AnnotatedDataFrame", data = pData(standard_cds)[colnames(spike_free_standard_cds),])
  fd <- new("AnnotatedDataFrame", data = fData(standard_cds)[rownames(spike_free_standard_cds),])

  spike_free_TPM <- newCellDataSet(apply(spike_free_standard_cds, 2, function(x) x / sum(x) * 10^6), 
                                        phenoData = pd, 
                                        featureData = fd, 
                                        expressionFamily=tobit(), 
                                        lowerDetectionLimit=1)
  pd <- new("AnnotatedDataFrame", data = pData(isoform_count_cds)[colnames(isoform_count_cds),])
  fd <- new("AnnotatedDataFrame", data = fData(isoform_count_cds)[rownames(isoform_count_cds)[1:(nrow(TPM_isoform_count_cds) - 97)],])
  
    spike_free_TPM_isoform_count_cds <- newCellDataSet(esApply(TPM_isoform_count_cds[1:(nrow(TPM_isoform_count_cds) - 97), ], 2, function(x) x / sum(x) * 10^6), 
                                            phenoData = pd, 
                                            featureData = fd, 
                                            expressionFamily = tobit(), 
                                            lowerDetectionLimit=1)
  
  #recover the transcript counts with the new algorithm (lower end ladder removed): 
  spike_free_Quake_norm_cds_optim_weight_fix_c <- relative2abs_optim_fix_c(relative_expr_matrix = exprs(spike_free_TPM), t_estimate = estimate_t(spike_free_TPM_isoform_count_cds, relative_expr_thresh = .1),                                                   
                                                                alpha_v = 1, total_RNAs = 50000, weight = 0.01, 
                                                                verbose = T, return_all = T, cores = 2, m =  -4.864207, c = mean(mean_m_c_select[1, ]))
  spike_free_optim_sum <- apply(Quake_norm_cds_optim_weight_fix_c$norm_cds[1:transcript_num, ], 2, sum)
  
  qplot(pData(absolute_cds)$endogenous_RNA[pData(absolute_cds)$endogenous_RNA > 1e3], 
      spike_free_optim_sum[pData(absolute_cds)$endogenous_RNA > 1e3], log="xy", color=pData(absolute_cds)$Time[pData(absolute_cds)$endogenous_RNA > 1e3], size = I(1)) + 
     geom_smooth(method="rlm", color="black", size = .1) + geom_abline(color="red") +  
    xlab("Total endogenous mRNA \n (spike-in)") +
    ylab("Total endogenous mRNA \n (spike-in free algorithm)") + #scale_size(range = c(0.25, 0.25)) + 
    scale_color_discrete(name = "Time points") + nm_theme()

  #benchmark the branching test (overlapping with group test as well as pseudotime tests): 
  #AT12_cds_subset_all_gene: remove the Cilia and Clara cells
  abs_group_test_res <- differentialGeneTest(abs_AT12_cds_subset_all_gene, 
                                                fullModelFormulaStr = "~Time", 
                                                reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)  
  abs_pseudotime_test_res <- differentialGeneTest(abs_AT12_cds_subset_all_gene, 
                                              fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)", 
                                              reducedModelFormulaStr = "~1", cores = detectCores(), relative = F)  

  #test whether or not the weight influence the results:
  abs_AT12_cds_subset_all_gene_res_no_weight <- branchTest(abs_AT12_cds_subset_all_gene[, ], cores = detectCores(), relative_expr = F, weighted = F)
  abs_AT12_cds_subset_all_gene_res_no_weight_relative <- branchTest(abs_AT12_cds_subset_all_gene[, ], cores = detectCores(), relative_expr = T, weighted = F)
  
  DEG_time_sets <- list(abs_group_test_res = row.names(abs_group_test_res[which(abs_group_test_res$qval < .01), ]), 
                        abs_pseudotime_test_res = row.names(abs_pseudotime_test_res[abs_pseudotime_test_res$qval < 0.01, ]),
                        abs_AT12_cds_subset_all_gene_res = row.names(abs_AT12_cds_subset_all_gene_res[abs_AT12_cds_subset_all_gene_res$qval < 0.01, ]),
                        abs_AT12_cds_subset_all_gene_res_no_weight = row.names(abs_AT12_cds_subset_all_gene_res_no_weight[abs_AT12_cds_subset_all_gene_res_no_weight$qval < 0.01, ]))
  
  overlap_genes <- Reduce(intersect,  DEG_time_sets) 
  pseudotime_element_all <- c(row.names(abs_group_test_res[which(abs_group_test_res$qval < .01), ]), 
                   row.names(abs_pseudotime_test_res[abs_pseudotime_test_res$qval < 0.01, ]),
                   row.names(abs_AT12_cds_subset_all_gene_res[abs_AT12_cds_subset_all_gene_res$qval < 0.01, ]), 
                   row.names(abs_AT12_cds_subset_all_gene_res_no_weight[abs_AT12_cds_subset_all_gene_res_no_weight$qval < 0.01, ]))
  pseudotime_sets_all <- c(rep(paste('Multiple timepoint test', sep = ''), length(which(abs_group_test_res$qval < .01))),
                rep(paste('Pseudotime test', sep = ''), length(which(abs_pseudotime_test_res$qval < 0.01))), 
                rep(paste('Branch test', sep = ''), length(which(abs_AT12_cds_subset_all_gene_res$qval < 0.01))), 
                rep(paste('Branch test (no weight)', sep = ''), length(which(abs_AT12_cds_subset_all_gene_res_no_weight$qval < 0.01))))

  venneuler_venn(pseudotime_element_all, pseudotime_sets_all)



#supplementary figures: 

branch_pseudotime_element_all <- c(
    #row.names(abs_group_test_res[which(abs_group_test_res$qval < .01), ]), 
    #           row.names(abs_pseudotime_test_res[abs_pseudotime_test_res$qval < 0.01, ]),
               row.names(relative_abs_AT12_cds_subset_all_gene[relative_abs_AT12_cds_subset_all_gene$qval < 0.01, ]), 
               # row.names(abs_AT12_cds_subset_all_gene_res_no_weight[abs_AT12_cds_subset_all_gene_res_no_weight$qval < 0.01, ]), 
    #           row.names(abs_AT12_cds_subset_all_gene_res_no_weight_relative[abs_AT12_cds_subset_all_gene_res_no_weight_relative$qval < 0.01, ]), 
               row.names(abs_pseudotime_test_lineage2_res[abs_pseudotime_test_lineage2_res$qval < 0.01, ]), 
               row.names(abs_pseudotime_test_lineage3_res[abs_pseudotime_test_lineage3_res$qval < 0.01, ])
               )
branch_pseudotime_sets_all <- c(
    #ep(paste('Multiple timepoint test', sep = ''), length(which(abs_group_test_res$qval < .01))),
    #        rep(paste('Pseudotime test', sep = ''), length(which(abs_pseudotime_test_res$qval < 0.01))), 
            rep(paste('Branch test', sep = ''), length(which(relative_abs_AT12_cds_subset_all_gene$qval < 0.01))), 
            # rep(paste('Branch test (no weight)', sep = ''), length(which(abs_AT12_cds_subset_all_gene_res_no_weight$qval < 0.01))), 
    #        rep(paste('Branch test (no weight, relative)', sep = ''), length(which(abs_AT12_cds_subset_all_gene_res_no_weight_relative$qval < 0.01))), 
            rep(paste('Pseudotime test (AT1 lineage)', sep = ''), length(which(abs_pseudotime_test_lineage2_res$qval < 0.01))), 
            rep(paste('Pseudotime test (AT2 lineage)', sep = ''), length(which(abs_pseudotime_test_lineage3_res$qval < 0.01)))
            )

# save(branch_pseudotime_element_all, branch_pseudotime_sets_all, file = 'branchTest_cmpr_subset')

# pdf(file = paste(elife_directory, 'eLife_fig_SI_branchTest_cmpr.pdf', sep = ''), height = 2, width = 3)
pdf(file = paste(elife_directory, 'eLife_fig_SI_branchTest_cmpr1.pdf', sep = ''))
venneuler_venn(branch_pseudotime_element_all, branch_pseudotime_sets_all)
dev.off()

#see the branch genes outside of AT1/2 pseudotime genes: 
branchGenes_example <- setdiff(row.names(subset(relative_abs_AT12_cds_subset_all_gene, qval < 0.01)),
   c(row.names(abs_pseudotime_test_lineage2_res[abs_pseudotime_test_lineage2_res$qval < 0.01, ]),
    row.names(abs_pseudotime_test_lineage3_res[abs_pseudotime_test_lineage3_res$qval < 0.01, ])))
plot_genes_branched_pseudotime2(abs_AT12_cds_subset_all_gene[branchGenes_example[6:10], ], color_by = "State", trajectory_color_by = 'Lineage', fullModelFormulaStr = '~sm.ns(Pseudotime, df = 3)*Lineage', normalize = T, stretch = T, lineage_labels = c('AT1', 'AT2'), cell_size = 1, ncol = 2, add_pval = T, reducedModelFormulaStr = '~sm.ns(Pseudotime, df = 3)') + nm_theme()+ ylab('Transcript counts') + xlab('Pseudotime')
plot_genes_branched_pseudotime2(abs_AT12_cds_subset_all_gene[branchGenes_example[6:10], ], color_by = "State", trajectory_color_by = 'Lineage', fullModelFormulaStr = '~sm.ns(Pseudotime, df = 3)*Lineage', normalize = T, stretch = T, lineage_labels = c('AT1', 'AT2'), cell_size = 1, ncol = 2, add_pval = T, reducedModelFormulaStr = '~sm.ns(Pseudotime, df = 3)') + nm_theme()+ ylab('Transcript counts') + xlab('Pseudotime')

pseudotime_element_all <- c(row.names(abs_group_test_res[which(abs_group_test_res$qval < .01), ]), 
                   row.names(abs_pseudotime_test_res[abs_pseudotime_test_res$qval < 0.01, ]),
                   row.names(abs_AT12_cds_subset_all_gene_res[abs_AT12_cds_subset_all_gene_res$qval < 0.01, ])
                   # row.names(abs_AT12_cds_subset_all_gene_res_no_weight[abs_AT12_cds_subset_all_gene_res_no_weight$qval < 0.01, ])
                   )
pseudotime_sets_all <- c(rep(paste('Multiple timepoint test', sep = ''), length(which(abs_group_test_res$qval < .01))),
                rep(paste('Pseudotime test', sep = ''), length(which(abs_pseudotime_test_res$qval < 0.01))), 
                rep(paste('Branch test', sep = ''), length(which(abs_AT12_cds_subset_all_gene_res$qval < 0.01)))
                # rep(paste('Branch test (no weight)', sep = ''), length(which(abs_AT12_cds_subset_all_gene_res_no_weight$qval < 0.01)))
                )

pdf(file = paste(elife_directory, 'eLife_fig_SI_branchTest_cmpr2.pdf', sep = ''))
venneuler_venn(pseudotime_element_all, pseudotime_sets_all)
dev.off()
    
# add a vertical line for the early / late lineage dependent genes to represent the bifurcation time points 

# ILR heatmap: don’t use blue / red color scheme for better representation of the idea of ILRs 

# tree branch plots with all cells on a lineage collapse to one branch 

# alpha - FDR plots for the two-group tests 

# alpha_fdr <- function(alpha_vec, est_pval, true_pval, est_pval_name, true_q_thrsld = 0.05, type = c('precision', 'recall', 'fdr')) {
#     qval <- p.adjust(est_pval, method = 'BH')
#     names(qval) <- est_pval_name
#     true_qval <- p.adjust(true_pval, method = 'BH')

#     P <- names(true_pval[true_qval <= true_q_thrsld])
#     N <- names(true_qval[true_qval > true_q_thrsld])

#     #FDR = v / (v + s)
#     unlist(lapply(alpha_vec, function(alpha, type, qval, P, N) {
#       fp <- setdiff(names(qval[qval <= alpha]), P) #false positive
#       tp <- intersect(names(qval[qval <= alpha]), P) #true positive
#       fn <- setdiff(names(qval[qval > alpha]), N)
      
#       if(type == 'precision') length(tp) / length(union(fp, tp)) #precision = tp / (tp + fp)
#       else if(type =='recall') length(tp) / length(union(fp, fn)) #recall = tp / (tp + fn)
#       else if(type == 'fdr')  length(fp) / length(union(fp, tp)) #fdr: fp / (fp + tp)
      
#     }, type = type, qval, P, N)) #/ sum(true_pval < alpha, na.rm = T)
# }

# alpha_fdr2 <- function() {
#     gene_list_true_data_list <- gene_list_true_data(p_thrsld = p_thrsld,
#         permutate_pval = permutate_pval[[ind]][gene_list],
#         na.rm = na.rm)
#     gene_list_new <- gene_list_true_data_list$gene_list
#     true_data <- gene_list_true_data_list$true_data
#     test_p_vec <- test_p_list[[ind]][gene_list_new]
#     TF_PN <- TF_PN_vec(true_data, test_p_vec)
# }

# alpha_vec <- seq(0, 1, length.out = 1000)
# true_pval <- permutate_pval #permutation_pval
# true_pval <- true_pval[gene_list] #gene_list
# #pval_df

# # fdr <- lapply(alpha_vec, function(x) alpha_fdr(pval_df[, 1], p.adjust(true_pval, method = 'BH'), x))
# # result <- mcmapply(alpha_fdr, split(t(alpha_vec), col(t(alpha_vec), as.factor = T)), split(as.matrix(pval_df), col(as.matrix(pval_df), as.factor = T)), true_pval, mc.cores = 8)
#   monocle_p <- new_std_diff_test_res[, 'pval'] 
#   names(monocle_p) <- row.names(new_std_diff_test_res)

# df3_pval_df <- data.frame(#monocle_p = monocle_p, 
#                       monocle_p_readcount = monocle_p_readcount,
#                       #mode_size_norm_permutate_ratio_by_geometric_mean = new_abs_size_norm_monocle_p_ratio_by_geometric_mean,
#                       #mc_mode_size_norm_permutate_ratio_by_geometric_mean = new_mc_size_norm_monocle_p_ratio_by_geometric_mean, 
#                       default_edgeR_p = default_edgeR_p, 
#                       #abs_default_edgeR_p = abs_default_edgeR_p,         
#                       default_deseq2_p = default_deseq2_p, 
#                       #abs_default_deseq2_p = abs_default_deseq2_p, 
#                       default_deseq_p = default_deseq_p, 
#                       #abs_default_deseq_p = abs_default_deseq_p, 
#                       scde_p = scde_p#, 
#                       #abs_scde_p = abs_scde_p
#                       )

# alpha_fdr_res <- apply(df3_pval_df, 2, function(x) alpha_fdr(alpha_vec, x, readcount_permutate_pval, row.names(df3_pval_df), type = 'fdr'))
# #alpha_fdr_res <- apply(pval_df, 2, function(x) alpha_fdr(alpha_vec, x, true_pval, row.names(pval_df), type = 'fdr'))

# p_alpha_fdr <- 
#     qplot(Var1 / 1000, value, geom = 'line', data = melt(alpha_fdr_res), color = Var2, linetype = Var2) + monocle_theme_opts() + 
#     geom_abline(color = 'red') + ggtitle('alpha VS fdr') + facet_wrap(~Var2, scale = 'free_y', ncol = round(sqrt(dim(alpha_fdr_res)))) + xlab('alpha') + ylab('fdr')

# abs_df3_pval_df <- data.frame(#monocle_p = monocle_p, 
#                       #monocle_p_readcount = monocle_p_readcount,
#                       mode_size_norm_permutate_ratio_by_geometric_mean = new_abs_size_norm_monocle_p_ratio_by_geometric_mean,
#                       mc_mode_size_norm_permutate_ratio_by_geometric_mean = new_mc_size_norm_monocle_p_ratio_by_geometric_mean, 
#                       #default_edgeR_p = default_edgeR_p, 
#                       abs_default_edgeR_p = abs_default_edgeR_p,         
#                       #default_deseq2_p = default_deseq2_p, 
#                       abs_default_deseq2_p = abs_default_deseq2_p, 
#                       #default_deseq_p = default_deseq_p, 
#                       abs_default_deseq_p = abs_default_deseq_p, 
#                       #scde_p = scde_p#, 
#                       abs_scde_p = abs_scde_p
#                       )

# abs_alpha_fdr_res <- apply(abs_df3_pval_df, 2, function(x) alpha_fdr(alpha_vec, x, mode_size_norm_permutate_ratio_by_geometric_mean, row.names(abs_df3_pval_df), type = 'fdr'))
# #alpha_fdr_res <- apply(pval_df, 2, function(x) alpha_fdr(alpha_vec, x, true_pval, row.names(pval_df), type = 'fdr'))

# p_abs_alpha_fdr <- 
#     qplot(Var1 / 1000, value, geom = 'line', data = melt(abs_alpha_fdr_res), color = Var2, linetype = Var2) + monocle_theme_opts() + 
#     geom_abline(color = 'red') + ggtitle('alpha VS fdr') + facet_wrap(~Var2, scale = 'free_y', ncol = round(sqrt(dim(abs_alpha_fdr_res)))) + xlab('alpha') + ylab('fdr')

# #find genes with expression goes up: 

# #or use the pval / qval from the global tests: 

######################################################################################################

# #fig b
# ABCs_df <- subset(ABCs_df, abs(ABCs) > 5)
# ABCs <- ABCs_df[, 'ABCs']
# names(ABCs) <- ABCs_df[, 'gene_short_name']
# pval <- abs_AT12_cds_subset_all_gene_res[ABCs_df[, 'gene_id'], 'pval']
# names(pval) <- names(ABCs)
# pval <- pval[!is.na(pval)]
# ABCs <- ABCs[names(pval)]

# gasRes <- auto_make_enrichment(gsaRes_go, 15, F, F, F, T, T)

# gasRes + nm_theme()
# ggsave(paste(elife_directory, 'eLife_fig3a.pdf', sep = ''), width = 6.5, height = 2.5)

# enrich_data_non_direction <- make_enrichment_df(std_bif_time_gsaRes_go, extract_num = 100,
#     custom_p_adjust = F, add_terms = F, direction = F)

# enrich_data_non_direction 

# enrich_data_non_direction <- enrich_data_non_direction[sort(as.vector(enrich_data_non_direction$"Stat (non-dir.)"),
#     index.return = T)$ix, ]
# qplot(x = 1:nrow(enrich_data_non_direction), y = abs(log(enrich_data_non_direction[, 'Stat (non-dir.)'])),
#             data = enrich_data_non_direction, geom = "bar", stat = "identity") +
#             coord_flip() + scale_x_discrete(limits = 1:nrow(enrich_data_non_direction),
#             labels = enrich_data_non_direction$Name) +
#             xlab("") + ylab("Normalized Enrichment Score")  + nm_theme()
# ggsave(paste(elife_directory, 'eLife_fig3a.pdf', sep = ''), height = 2, width = 4)

# #debug buildLineageBranchCellDataSet for weight_constant: 
# str_logfc_df_list <- calILRs(cds = std_AT12_cds_subset_all_gene[add_quake_gene_all_marker_ids, ], lineage_states = c(2, 3), stretch = T, cores = 1, 
#     trend_formula = "~sm.ns(Pseudotime, df = 3)*Lineage", ILRs_limit = 3, 
#     relative_expr = F, weighted = FALSE, label_by_short_name = F, 
#     useVST = FALSE, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "str_logfc_df", return_all = T)

# #calILRs for all genes is not feasible... 

# #fig c
# str_logfc_df <- t(str_logfc_df)
# str_logfc_df <- str_logfc_df[!is.na(str_logfc_df[, 1]), ]
# str_logfc_df <- str_logfc_df[abs(str_logfc_df[, 100]) > 1, ]
# plot_ILRs_heatmap(absolute_cds, str_logfc_df, abs_AT12_cds_subset_all_gene_ABCs, relative_abs_AT12_cds_subset_quake_gene, "ensemble_id", ABC_type = "all",
#     dist_method = "euclidean", hclust_method = "ward", ILRs_limit = 3,
#     cluster_num = 4)
# pdf(file = paste(elife_directory, 'eLife_fig3b.pdf', sep = ''))

# #fig 3e: 
# #test the the AT1/2 early late group with the proportity score: 
# c('Clic5', 'Muc1', 'S100g', 'Soat1')
# bifurcation_time[c('Clic5', 'Muc1', 'S100g', 'Soat1')]

# bifurcation_time <- detectBifurcationPoint(abs_AT12_cds_subset_all_gene_ILRs[, 27:100])
# valid_bifurcation_time  <- bifurcation_time[!is.na(bifurcation_time)]
# valid_bifurcation_time <- valid_bifurcation_time[unique(names(valid_bifurcation_time))]
# bif_time_gsaRes_go <- runGSA(valid_bifurcation_time, gsc = mouse_go_gsc, ncpus = 1) 



