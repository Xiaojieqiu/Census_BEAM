  use_select_algorithm = TRUE

  library(monocle)
  library(xacHelper)
  elife_directory = "./"

  load_all_libraries()
  # ########################################################prepare the Quake dataset########################################################
  
  #load some meta data from quake paper: 
  load('./Quake_data/SRR_cell_type') #cell type specification from Quake lung paper
  load('./Quake_data/quake_gene_name') #genes studied in Quake lung paper
  transcript_num <- 38919 #number of genes in the genome, excluding ERCC spikein 

  # Load up the Quake paper data as relative expression values (classic FPKM)
  
  sra_table <- read.delim("./Quake_data/quake_lung/SraRunTable.txt")
  sra_table$Run <- paste(sra_table$Run, "_0", sep="")
  
  fpkm_matrix <- read.delim("./Quake_data/quake_lung/standard_normalized_out/genes.fpkm_table", row.names="tracking_id")
  #fpkm_matrix <- fpkm_matrix[,-1]
  sample_sheet <- read.delim("./Quake_data/quake_lung/standard_normalized_out/samples.table", row.names="sample_id")
  sample_sheet <- merge(sample_sheet, sra_table, by.x="row.names", by.y="Run")
  row.names(sample_sheet) <- sample_sheet$Row.names 
  sample_sheet <- sample_sheet[,-1]
  sample_sheet$Time <- NULL
  sample_sheet$Time[sample_sheet$age == "Embryonic day 14.5"] <- "E14.5"
  sample_sheet$Time[sample_sheet$age == "Embryonic day 16.5"] <- "E16.5"
  sample_sheet$Time[sample_sheet$age == "Embryonic day 18.5"] <- "E18.5"
  sample_sheet$Time[sample_sheet$age == "post natal day 107"] <- "Adult"
  
  #fpkm_matrix <- fpkm_matrix * sample_sheet$internal_scale ############### this is wrong since each matrix row will multiple one value in the vector 
  # fpkm_matrix <- t(t(fpkm_matrix) * sample_sheet$internal_scale)
  gene_ann <- read.delim("./Quake_data/quake_lung/standard_normalized_out/genes.attr_table", row.names="tracking_id")
  
  gencode_biotypes <- read.delim("./Quake_data/quake_lung/gencode_biotypes.txt")
  gencode_biotypes <- gencode_biotypes[,c(1,2)]
  spike_names <- row.names(fpkm_matrix)[grepl("ERCC", row.names(fpkm_matrix))]
  spike_biotypes <- data.frame(gene_id=spike_names, biotype="spike")
  
  gencode_biotypes <- rbind(gencode_biotypes,spike_biotypes)
  
  gene_ann <- merge(gene_ann, gencode_biotypes, by.x="row.names", by.y="gene_id")
  row.names(gene_ann) <- gene_ann$Row.names
  gene_ann <- gene_ann[,c(-1)]
  
  pd <- new("AnnotatedDataFrame", data = sample_sheet)
  fd <- new("AnnotatedDataFrame", data = gene_ann)
  
  fpkm_matrix <- as.matrix(fpkm_matrix)
  fpkm_matrix <- fpkm_matrix[row.names(gene_ann),]
  fpkm_matrix <- fpkm_matrix[,row.names(sample_sheet)]
  
  standard_cds <- newCellDataSet(fpkm_matrix, 
                                 expressionFamily=tobit(), 
                                 phenoData = pd, 
                                 featureData = fd,
                                 lowerDetectionLimit=1)
  
  pData(standard_cds)$Total_mRNAs <- colSums(exprs(standard_cds))
  pData(standard_cds)$endogenous_RNA <- colSums(exprs(standard_cds[1:transcript_num, ]))
  
  valid_genes <- row.names(subset(fData(standard_cds), biotype %in% c("spike") == FALSE))
  
  # Load up the ERCC spike-in control metadata, which we will use to 
  # convert the standard FPKM values to absolute values 
  input.ERCC.annotation<-read.delim("./Quake_data/quake_lung/ERCC_specification.txt", header=T)
  colnames(input.ERCC.annotation)<-c("Resort_ID",
                                     "ERCC_ID",
                                     "subgroup",
                                     "conc_attomoles_ul_Mix1",
                                     "conc_attomoles_ul_Mix2",
                                     "exp_fch_ratio",
                                     "log2_Mix1_Mix2")
  
  # So we can index this data frame by ERCC transcript ID  							   
  rownames(input.ERCC.annotation)<-input.ERCC.annotation[,"ERCC_ID"]
  
  # Extract the spike-in FPKM values for each cell as a CellDataSet subset
  ercc_controls <- standard_cds[rownames(input.ERCC.annotation),]
  
  # Compute the mean spike FPKMs
  mean_spike_fpkms <- rowMeans(exprs(ercc_controls))
  spike_df <- input.ERCC.annotation 
  
  # Merge the means with the metadata so we can compare mean FPKM to spike-in concentration later
  spike_df <- cbind(spike_df, mean_spike_fpkms[row.names(spike_df)])
  colnames(spike_df)[length(colnames(spike_df))] <- "FPKM"
  spike_df$numMolecules <- spike_df$conc_attomoles_ul_Mix1*(10*10^(-3)*1/40000*10^(-18)*6.02214179*10^(23))
  pData(ercc_controls)$Cell <- row.names(pData(ercc_controls))
  
  # Fit a robust linear model of spike concentration vs. FPKM for each cell
  molModels <- esApply(ercc_controls, 2, function(cell_exprs, input.ERCC.annotation) {
    
    #print (cell_exprs)
    spike_df <- input.ERCC.annotation 
    spike_df <- cbind(spike_df, cell_exprs[row.names(spike_df)])
    colnames(spike_df)[length(colnames(spike_df))] <- "FPKM"
    spike_df$numMolecules <- spike_df$conc_attomoles_ul_Mix1*(10*10^(-3)*1/40000*10^(-18)*6.02214179*10^(23))
    spike_df$rounded_numMolecules <- round(spike_df$conc_attomoles_ul_Mix1*(10*10^(-3)*1/40000*10^(-18)*6.02214179*10^(23)))
    # print (mean(spike_df$rounded_numMolecules))
    
    spike_df <- subset(spike_df, FPKM >= 1e-10)
    spike_df$log_fpkm <- log10(spike_df$FPKM) 
    spike_df$log_numMolecules <- log10(spike_df$numMolecules)
    
    # print (paste("geometric mean of log numMol ", mean(spike_df$log_numMolecules), "geometric mean of log FPKM ", mean(spike_df$log_fpkm)))
           
   molModel <- tryCatch({
     #molModel <- vgam (rounded_numMolecules ~ sm.ns(log_fpkm, df=3), data=spike_df, family=negbinomial(zero=NULL))
     molModel <- rlm(log_numMolecules ~ log_fpkm, data=spike_df)
     
     #
     # 		    qp <- qplot(FPKM, numMolecules, data=spike_df, log="xy") +
     # 			geom_abline(color="green") +
     # geom_smooth(aes(y=10^log_numMolecules), method="lm") +
     # geom_smooth(aes(y=10^log_numMolecules), method="rlm", color="red") +
     # geom_line(aes(y=10^predict(molModel,type="response")))
     # print(qp)
     molModel
   }, 
   error = function(e) { print(e); NULL })
   molModel
  }, input.ERCC.annotation)
  
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
  split(exprs(standard_cds), rep(1:ncol(exprs(standard_cds)), each = nrow(exprs(standard_cds)))), 
  molModels)
  
  row.names(norm_fpkms) <- row.names(exprs(standard_cds))
  colnames(norm_fpkms) <- colnames(exprs(standard_cds))
  
  norm_fpkms_melted <- melt(norm_fpkms)
  total_RNA_df <- ddply(norm_fpkms_melted, .(Var2), function(x) { data.frame(total=sum(x$value), mode_mol=mlv(x$value[x$value > 0.0], method = "mfv")$M) }) #estimate the mode
  
  total_RNA_df <- merge(sample_sheet, total_RNA_df, by.x="row.names", by.y="Var2")
  
  # Now let's generate a new CellDataSet that uses the absolute transcript counts
  fpkm_matrix_abs <- norm_fpkms
  colnames(fpkm_matrix_abs) <- colnames(fpkm_matrix)
  row.names(fpkm_matrix_abs) <- row.names(fpkm_matrix)
  #fpkm_matrix <- fpkm_matrix[,-1]
  
  pd <- new("AnnotatedDataFrame", data = sample_sheet[colnames(fpkm_matrix_abs),])
  fd <- new("AnnotatedDataFrame", data = gene_ann[rownames(fpkm_matrix_abs),])
  
  absolute_cds <- newCellDataSet(fpkm_matrix_abs, 
                                 phenoData = pd, 
                                 featureData = fd, 
                                 expressionFamily=negbinomial(), 
                                 lowerDetectionLimit=1)
  
  pData(absolute_cds)$Total_mRNAs <- esApply(absolute_cds, 2, sum)
  pData(absolute_cds)$endogenous_RNA <- esApply(absolute_cds, 2, function(x) sum(x[1:transcript_num]))

  capture_efficiency <- 1 #we don't consider  the capture efficiency adjustment
  
  #generate the absolute transcript count cds recovered with ERCC-spikein transcripts after removing the lower end: 
  molModels_select <- esApply(ercc_controls, 2, function(cell_exprs, input.ERCC.annotation, volume = volume, dilution = dilution, select_above_thresh = T, capture_efficiency) {
    
    #print (cell_exprs)
    spike_df <- input.ERCC.annotation 
    spike_df <- cbind(spike_df, cell_exprs[row.names(spike_df)])
    colnames(spike_df)[length(colnames(spike_df))] <- "FPKM"
    spike_df$numMolecules <- spike_df$conc_attomoles_ul_Mix1*(volume*10^(-3)*dilution*10^(-18)*6.02214179*10^(23))
    if(!is.na(capture_efficiency))
      spike_df$numMolecules <- spike_df$numMolecules * capture_efficiency
    spike_df$rounded_numMolecules <- round(spike_df$conc_attomoles_ul_Mix1*(volume*10^(-3)*dilution*10^(-18)*6.02214179*10^(23)))
    # print (mean(spike_df$rounded_numMolecules))
    
    if(select_above_thresh){
      spike_df <- subset(spike_df, conc_attomoles_ul_Mix1 >= 800)
      spike_df <- subset(spike_df, FPKM >= 1e-10)    
    }
    spike_df$log_fpkm <- log10(spike_df$FPKM) 
    spike_df$log_numMolecules <- log10(spike_df$numMolecules)
    
    # print (paste("geometric mean of log numMol ", mean(spike_df$log_numMolecules), "geometric mean of log FPKM ", mean(spike_df$log_fpkm)))
    
    molModel <- tryCatch({
      #molModel <- vgam (rounded_numMolecules ~ sm.ns(log_fpkm, df=3), data=spike_df, family=negbinomial(zero=NULL))
      molModel <- rlm(log_numMolecules ~ log_fpkm, data=spike_df)
      
      #
      #     qp <- qplot(FPKM, numMolecules, data=spike_df, log="xy") +
      # geom_smooth(aes(y=10^log_numMolecules), method="lm") +
      # geom_smooth(aes(y=10^log_numMolecules), method="rlm", color="red") +
      # geom_line(aes(y=10^predict(molModel,type="response")))
      # print(qp)
      molModel
    }, 
    error = function(e) { print(e); NULL })
    molModel
  }, input.ERCC.annotation, volume = 10, dilution = 1/ 40000, capture_efficiency = capture_efficiency)
  
  #calculate m_mean_x_ij and c_mean_y_ij of all the ERCC spikein transcripts (j) for each cell i, lower end of ladder removed 
  #we don't use TPM_ercc_controls here because c_mean_y_ij will be the same (the absolute_cds_select will be the same no matter what we choose to recover the transcript counts)
  mean_m_c_select <- esApply(ercc_controls, 2, function(cell_exprs, input.ERCC.annotation, volume = volume, dilution = dilution, select_above_thresh = T, capture_efficiency) {
    
    #print (cell_exprs)
    spike_df <- input.ERCC.annotation 
    spike_df <- cbind(spike_df, cell_exprs[row.names(spike_df)])
    colnames(spike_df)[length(colnames(spike_df))] <- "FPKM"
    spike_df$numMolecules <- spike_df$conc_attomoles_ul_Mix1*(volume*10^(-3)*dilution*10^(-18)*6.02214179*10^(23))
    if(!is.na(capture_efficiency))
      spike_df$numMolecules <- spike_df$numMolecules * capture_efficiency
    spike_df$rounded_numMolecules <- round(spike_df$conc_attomoles_ul_Mix1*(volume*10^(-3)*dilution*10^(-18)*6.02214179*10^(23)))
    # print (mean(spike_df$rounded_numMolecules))
    
    if(select_above_thresh){
      spike_df <- subset(spike_df, conc_attomoles_ul_Mix1 >= 800)
      spike_df <- subset(spike_df, FPKM >= 1e-10)    
    }
    spike_df$log_fpkm <- log10(spike_df$FPKM) 
    spike_df$log_numMolecules <- log10(spike_df$numMolecules)
    
    #   print (paste("geometric mean of log numMol ", mean(spike_df$log_numMolecules), "geometric mean of log FPKM ", mean(spike_df$log_fpkm)))
    return(c(c_mean_y_ij = mean(spike_df$log_numMolecules), m_mean_x_ij = mean(spike_df$log_fpkm)))
    
    molModel <- tryCatch({
      #molModel <- vgam (rounded_numMolecules ~ sm.ns(log_fpkm, df=3), data=spike_df, family=negbinomial(zero=NULL))
      molModel <- rlm(log_numMolecules ~ log_fpkm, data=spike_df)
      
      #
      #     qp <- qplot(FPKM, numMolecules, data=spike_df, log="xy") +
      # geom_smooth(aes(y=10^log_numMolecules), method="lm") +
      # geom_smooth(aes(y=10^log_numMolecules), method="rlm", color="red") +
      # geom_line(aes(y=10^predict(molModel,type="response")))
      # print(qp)
      molModel
    }, 
    error = function(e) { print(e); NULL })
    molModel
  }, input.ERCC.annotation, volume = 10, dilution = 1/40000, capture_efficiency = capture_efficiency)
    
  # Load up the isoform-level read counts
  count_matrix <- read.delim("./Quake_data/quake_lung/standard_normalized_out/isoforms.fpkm_table", row.names="tracking_id")
  #fpkm_matrix <- fpkm_matrix[,-1]
  isoform_ann <- read.delim("./Quake_data/quake_lung/standard_normalized_out/isoforms.attr_table", row.names="tracking_id")
  isoform_ann$isoform_id <- row.names(isoform_ann)
  
  isoform_ann <- merge(isoform_ann, gene_ann[, c("gene_id", "biotype")], by.x="gene_id", by.y="row.names")
  #gencode_biotypes <- read.delim("gencode_biotypes.txt")
  row.names(isoform_ann) <- isoform_ann$isoform_id
  isoform_ann <- isoform_ann[,-1]
  
  pd <- new("AnnotatedDataFrame", data = sample_sheet)
  fd <- new("AnnotatedDataFrame", data = isoform_ann)
  
  count_matrix <- as.matrix(count_matrix)
  count_matrix <- count_matrix[row.names(isoform_ann),]
  count_matrix <- count_matrix[,row.names(sample_sheet)]
  
  isoform_count_cds <- newCellDataSet(count_matrix, 
                                      phenoData = pd, 
                                      featureData = fd, 
                                      expressionFamily=negbinomial(), 
                                      lowerDetectionLimit=1)
  pData(isoform_count_cds)$Total_mRNAs <- esApply(isoform_count_cds, 2, sum)
  pData(isoform_count_cds)$endogenous_RNA <- esApply(isoform_count_cds, 2, function(x) sum(x[1:transcript_num]))

  #convert the isoform fpkm value into absolute transcript counts: 
  iso_norm_fpkms <- mapply(function(cell_exprs, molModel) {
    tryCatch({
      norm_df <- data.frame(log_fpkm=log10(cell_exprs))
      res <- 10^predict(molModel, type="response", newdata=norm_df)
    }, 
    error = function(e) {
      rep(NA, length(cell_exprs))
    })
  }, 
  split(exprs(isoform_count_cds), rep(1:ncol(exprs(isoform_count_cds)), each = nrow(exprs(isoform_count_cds)))), 
  molModels_select)
  
  row.names(iso_norm_fpkms) <- row.names(exprs(isoform_count_cds))
  colnames(iso_norm_fpkms) <- colnames(exprs(isoform_count_cds))

  pd <- new("AnnotatedDataFrame", data = pData(isoform_count_cds))
  fd <- new("AnnotatedDataFrame", data = fData(isoform_count_cds))
  
  iso_absolute_cds <- newCellDataSet(iso_norm_fpkms, 
                                 phenoData = pd, 
                                 featureData = fd, 
                                 expressionFamily=negbinomial(), 
                                 lowerDetectionLimit=1)
 
  #use the algorithm to recover the transcript counts based on mode inferred from isoform data 
  norm_matrix <- relative2abs(standard_cds, t_estimate = estimate_t(exprs(isoform_count_cds)))
  mc_adj_cds <- newCellDataSet(as.matrix(norm_matrix),
                               phenoData = new("AnnotatedDataFrame", data = pData(standard_cds)),
                               featureData = new("AnnotatedDataFrame", data = fData(standard_cds)),
                               expressionFamily = negbinomial(),
                               lowerDetectionLimit = 1)
  pData(mc_adj_cds)$Total_mRNAs <- esApply(mc_adj_cds, 2, sum)
  pData(mc_adj_cds)$endogenous_RNA <- esApply(mc_adj_cds, 2, function(x) sum(x[1:transcript_num]))
  
  #read the read count data for the genes: 
  dir = "/net/trapnell/vol1/home/xqiu/Projects/BEAM/Xiaojie_reproduce/quake_lung/standard_normalized_out"
  sample_table <- read.delim(paste(dir, "/samples.table", sep = ''))
  norm_count <- read.delim(paste(dir, "/genes.count_table", sep = ''))
  row.names(norm_count) <- norm_count$tracking_id
  norm_count <- norm_count[, -1]
  
  read_countdata <- round(t(t(norm_count) * sample_table$internal_scale)) #convert back to the raw counts 
  read_countdata <- read_countdata[row.names(standard_cds), colnames(standard_cds)]

  read_countdata_cds <- newCellDataSet(as.matrix(read_countdata),
                               phenoData = new("AnnotatedDataFrame", data = pData(standard_cds)),
                               featureData = new("AnnotatedDataFrame", data = fData(standard_cds)),
                               expressionFamily = negbinomial(),
                               lowerDetectionLimit = 1)
  pData(read_countdata_cds)$Total_mRNAs <- esApply(read_countdata_cds, 2, sum)
  pData(read_countdata_cds)$endogenous_RNA <- esApply(read_countdata_cds, 2, function(x) sum(x[1:transcript_num]))
  
  #convert FPKM values to TPM (a relative abundance measure) for the following analysis: 
  TPM <- esApply(standard_cds, 2, function(x) x / sum(x) * 10^6)
  tpm_align <- melt(TPM)
  tpm_align <- tpm_align[tpm_align$value > 0.1, ]
  
  pdf('tpm_align_distr.pdf')
  qplot(value, log = 'x', geom = 'density', color = Var2, data = tpm_align) + theme_bw() + theme(legend.position = 'none') + xlab('TPM')
  dev.off()

  TPM_cds <- newCellDataSet(TPM, 
                            phenoData = new("AnnotatedDataFrame", data = pData(standard_cds)),
                            featureData = new("AnnotatedDataFrame", data = fData(standard_cds)),
                            expressionFamily=tobit(), 
                            lowerDetectionLimit=1)
  pData(TPM_cds)$Total_mRNAs <- esApply(TPM_cds, 2, sum)
  pData(TPM_cds)$endogenous_RNA <- esApply(TPM_cds[1:transcript_num, ], 2, sum)

  TPM_ercc_controls <- newCellDataSet(TPM[row.names(ercc_controls), colnames(ercc_controls)], 
                                      phenoData = new("AnnotatedDataFrame", data = pData(ercc_controls)),
                                      featureData = new("AnnotatedDataFrame", data = fData(ercc_controls)), 
                                      expressionFamily=tobit(), 
                                      lowerDetectionLimit=1)
  
  #generate all the parameters used in the recovery algorithm, we can make a gigantic scatter plot for those parameters too 
  tpm_test_selected <- all_algorithm_parameters(ercc_control = TPM_ercc_controls, select_above_thresh = T, capture_efficiency = capture_efficiency)
  tpm_test_selected <- do.call(rbind.data.frame, tpm_test_selected)
  tpm_test_selected$Cell <- str_split_fixed(row.names(tpm_test_selected), "[*.*]", 2)[, 1] 
  # tpm_test$Time <- pData(absolute_cds[, tpm_test$Cell])$Time
  fraction <- 1- pData(absolute_cds)[, 'endogenous_RNA'] / pData(absolute_cds)[, 'Total_mRNAs'] 
  names(fraction) <- str_split_fixed(colnames(absolute_cds), "[*.*]", 2)[, 1]
  tpm_test_selected$z_fraction <- log10(fraction[tpm_test_selected$Cell])
  tpm_test_selected$R_plus_z <- tpm_test_selected[, 'R_rho_ij'] + tpm_test_selected[, 'z_fraction']
  
  pdf('molModel_kb_relationship.pdf')
  qplot(unlist(lapply(molModels_select, function(x) coef(x)[2])), unlist(lapply(molModels_select, function(x) coef(x)[1]))) + 
    geom_abline(intercept = mean(tpm_test_selected$c_mean_y_ij), slope = -mean(tpm_test_selected$m_mean_x_ij))
  dev.off()

  #plot all the parameters
  tpm_data <- tpm_test_selected[, -which(colnames(tpm_test_selected) %in% c('Cell', 'Time'))]
  pdf('compare_all_algorithm_param.pdf')
  plotmatrix2(tpm_data, myLab= 'Parameters related to the recovery algorithm', read_count = NULL, useCounts = T, logMode = F, pseudocount = 0, alpha = 1)
  dev.off()

  norm_fpkms_select <- mapply(function(cell_exprs, molModel) {
    tryCatch({
      norm_df <- data.frame(log_fpkm=log10(cell_exprs))
      res <- 10^predict(molModel, type="response", newdata=norm_df)
    }, 
    error = function(e) {
      rep(NA, length(cell_exprs))
    })
  }, 
  split(exprs(standard_cds), rep(1:ncol(exprs(standard_cds)), each = nrow(exprs(standard_cds)))), 
  molModels_select)
  
  row.names(norm_fpkms_select) <- row.names(exprs(standard_cds))
  colnames(norm_fpkms_select) <- colnames(exprs(standard_cds))
  
  norm_fpkms_melted_select <- melt(norm_fpkms_select)
  # total_RNA_df <- ddply(norm_fpkms_melted, .(Var2), function(x) { data.frame(total=sum(x$value), mode_mol=mlv(x$value[x$value > 0.0], method = "mfv")$M) })
  
  # total_RNA_df <- merge(sample_sheet, total_RNA_df, by.x="row.names", by.y="Var2")
  
  # Now let's generate a new CellDataSet that uses the absolute transcript counts
  fpkm_matrix_abs_select <- norm_fpkms_select
  colnames(fpkm_matrix_abs_select) <- colnames(standard_cds)
  row.names(fpkm_matrix_abs_select) <- row.names(fpkm_matrix)
  #fpkm_matrix <- fpkm_matrix[,-1]
  
  pd <- new("AnnotatedDataFrame", data = pData(absolute_cds)[colnames(fpkm_matrix_abs),])
  fd <- new("AnnotatedDataFrame", data = fData(absolute_cds)[rownames(fpkm_matrix_abs),])
  
  absolute_cds_select <- newCellDataSet(fpkm_matrix_abs_select, 
                                        phenoData = pd, 
                                        featureData = fd, 
                                        expressionFamily=negbinomial(), 
                                        lowerDetectionLimit=1)
  pd <- new("AnnotatedDataFrame", data = pData(isoform_count_cds)[colnames(isoform_count_cds),])
  fd <- new("AnnotatedDataFrame", data = fData(isoform_count_cds)[rownames(isoform_count_cds),])
  
  TPM_isoform_count_cds <- newCellDataSet(esApply(isoform_count_cds, 2, function(x) x / sum(x) * 10^6), 
                                          phenoData = pd, 
                                          featureData = fd, 
                                          expressionFamily=tobit(), 
                                          lowerDetectionLimit=1)
  
  #recover the transcript counts with the new algorithm (lower end ladder removed): 
  Quake_norm_cds_optim_weight_fix_c <- relative2abs_optim_fix_c(relative_expr_matrix = exprs(TPM_cds), t_estimate = estimate_t(TPM_isoform_count_cds, relative_expr_thresh = .1),                                                   
                                                                alpha_v = 1, total_RNAs = 50000, weight = 0.01, 
                                                                verbose = T, return_all = T, cores = 2, m =  -4.864207, c = mean(mean_m_c_select[1, ]))
  optim_sum <- apply(Quake_norm_cds_optim_weight_fix_c$norm_cds[1:transcript_num, ], 2, sum)
  cmpr_Quake_norm_cds_optim_weight_fix_c <- relative2abs(TPM_cds, t_estimate = estimate_t(TPM_isoform_count_cds, relative_expr_thresh = .1),                                                   
                                                                alpha_v = 1, total_RNAs = 50000, weight = 0.01, 
                                                                verbose = T, return_all = T, cores = 2, m =  -4.864207, c = mean(mean_m_c_select[1, ]))
  pdf('opt_m_fix_c_algorithm.pdf')
  qplot(esApply(absolute_cds_select[1:transcript_num, ], 2, sum)[optim_sum > 10], apply(Quake_norm_cds_optim_weight_fix_c$norm_cds[1:transcript_num, ], 2, sum)[optim_sum > 10], log = 'xy') + geom_abline() + geom_smooth(method = 'rlm') + xlab('spikein total RNA counts') + 
    geom_smooth(method = 'rlm') + ylab('Optimized total RNA counts') + ggtitle('Quake lung data: Optimizing m with fixed c directly calculated')
  dev.off()

  Quake_norm_cds_optim_mc <- relative2abs_mc(relative_expr_matrix = exprs(TPM_cds), t_estimate = estimate_t(TPM_isoform_count_cds, relative_expr_thresh = .1),                                                   
                                             verbose = T, return_all = T, cores = 2, m =  -4.277778, c = 2.932929)

  pdf('true_mc_algorithm.pdf')
  qplot(esApply(absolute_cds_select[1:transcript_num, ], 2, sum)[optim_sum > 10], apply(Quake_norm_cds_optim_mc$norm_cds[1:transcript_num, ], 2, sum)[optim_sum > 10], log = 'xy') + geom_abline() + geom_smooth(method = 'rlm') + xlab('spikein total RNA counts') + 
    geom_smooth(method = 'rlm') + ylab('Optimized total RNA counts') + ggtitle('Quake lung data: Use the true m/c values')
  dev.off()

  pdf('lower_end_removed.pdf')
  qplot(esApply(absolute_cds[1:transcript_num, ], 2, sum), esApply(absolute_cds_select[1:transcript_num, ], 2, sum), log = 'xy') + geom_abline() + geom_smooth(method = 'rlm') + xlab('spikein total RNA counts (all expressed spikein transcripts)') + 
    geom_smooth(method = 'rlm') + ylab('spikein total RNA counts (lower end removed)') + ggtitle('Comparing the effect of removing lower end or not')
  dev.off()

  ################################## this long section is used to generate all the data for making the figure 1e: ##################################
  #generate the test p-val as well as the permutation pval
  
  #select cells for two group tests (E18.5 vs E14.5): 
  valid_cells <- colnames(absolute_cds)[colSums(round(exprs(absolute_cds[valid_genes,])))  > 9500] #find cells before change expression values to ensure the cells are consistent with previous analysis 
  
  #use the new algorithm recovered transcript counts for all downstream analysis (lower end ladder removed): 
  if(use_select_algorithm) {
    exprs(absolute_cds) <- exprs(absolute_cds_select)
    exprs(mc_adj_cds) <- Quake_norm_cds_optim_weight_fix_c$norm_cds
    pData(mc_adj_cds)$endogenous_RNA <- esApply(mc_adj_cds, 2, function(x) sum(x[1:transcript_num]))
    pData(absolute_cds)$endogenous_RNA <- esApply(absolute_cds, 2, function(x) sum(x[1:transcript_num]))
    fraction <- 1- pData(absolute_cds)[, 'endogenous_RNA'] / pData(absolute_cds)[, 'Total_mRNAs'] 
  }
  #################################################################USE marker genes from Quake paper to reduce dimension and then find the bifurcation trajectories#######################################################################
  AT2 <- c("S100g", "Slc34a2", "Lamp3", "Bex2", "Cd36", "Etv5", "Scd1", "Il33", "Fabp5", "Hc", "Sftpa1", "Dlk1", "Cxcl15", "Chsy1", "Lcn2", "Fasn", "Sftpb", "Egfl6", "Lyz1", "Sftpc", "Soat1", "Glrx", "Mlc1", "Retnla", "Mid1ip1", "Chi3l1", "Rab27a", "Fabp12", "Lgi3", "Lyz2")
  
  AT1 <- c("Clic5", "Pdpn", "Akap5", "Lmo7", "Ahnak", "Sdpr", "Emp2", "Col4a3", "Rtkn2", "Cav1", "Timp3", "Ager", "Msn", "Limch1", "Hopx", "Ptrf", "S100a6", "Dpysl2", "Agrn", "Akap2", "Pmp22", "Tinagl1", "Lgals3", "S100a14", "Malat1", "Vegfa", "Samhd1", "Qk", "Sema3a", "Aqp5")
  
  Clara <- c("Cyp2f2", "Scgb3a2", "Chad", "Hp", "Scgb1a1", "Krt15", "Cbr2", "Itm2a", "Upk3a", "Rassf9", "Col23a1", "Lypd2", "Nupr1", "Ppap2b", "Cd200", "Mtus1", "Mia1", "1810010H24Rik", "Dcxr", "Aldh1a7", "Osgin1", "Kdr", "Ccnd1", "Ccnd2", "Gsta3", "Tacstd2", "Nrarp", "Iyd", "Pir")
  
  Ciliated <- c("Ncs1", "Ccdc67", "1700007G11Rik", "Foxj1", "1110017D15Rik", "Ccdc19", "Stk33", "Spag17", "Ccdc39", "Lrriq1", "4930451C15Rik", "Efhc1", "Lrrc23", "Ttc18", "Mlf1", "6820408C15Rik", "Traf3ip1", "Kif9", "Wdr16", "Lrrc36", "Fank1", "Enkur", "1700009P17Rik", "Tekt4", "Ccdc40", "Hs6st2", "Ccdc113", "Dnph1", "Mcidas")
  
  
  Ciliated_novel <- c("1110017D15Rik", "1700007G11Rik", "1700009P17Rik", "4930451C15Rik", "6820408C15Rik", "Dnph1", "Ccdc113", "Ccdc19", "Ccdc39", "Ccdc40", "Ccdc67", "Efhc1", "Fank1", "Foxj1", "Hs6st2", "Kndc1", "Lrrc23", "Lrriq1", "Ncs1", "Stk33", "Tekt4", "Ttc18", "Wdr16", "Kif9", "Lrrc36", "Mlf1", "Spag17", "Traf3ip1", "Enkur", "Mcidas")
  
  Calra_novel <- c("Chad", "Cyp2f2", "Scgb3a2", "Hp", "Itm2a", "Krt15", "Col23a1", "Upk3a", "Cd200", "Rassf9", "Cbr2", "Nupr1", "Scgb1a1", "1810010H24Rik", "Lypd2", "Ppap2b", "Ccnd2", "Kdr", "Osgin1", "Cd24a", "Dcxr", "Ccnd1", "Mia1", "Aldh1a7", "Gsto1", "Iyd", "Aldh1a1", "C3", "Bcl6", "Tacstd2")
  
  AT1_novel <- c("Clic5", "Akap5", "Ahnak", "Lmo7", "Pdpn", "Emp2", "Timp3", "Malat1", "Cav1", "Dpysl2", "Sema3a", "Col4a3", "Ager", "Ptrf", "Rtkn2", "Lgals3", "Sdpr", "Cldn18", "Samhd1", "Agrn", "Pmp22", "Qk", "Spnb2", "S100a14", "Scnn1g", "Clic3", "Msn", "Hopx", "Tspan8", "Vegfa")
  
  AT2_novel  <-  c("Etv5", "Lamp3", "Slc34a2", "S100g", "Scd1", "Ppp1r14c", "Sftpa1", "Hc", "Egfl6", "Fabp5", "Glrx", "Cd36", "Dlk1", "Il33", "Lcn2", "Cxcl15",   "Sftpb", "Slc31a1", "Bex2", "Chsy1", "Lyz1", "Sftpc", "Fabp12", "Fasn", "Soat1", "Chi3l1", "Rab27a", "Mid1ip1", "Nucb2", "Rab27b")
  
  AT1_early <- c("Aqp5", "Pdpn", "Rtkn2", "Ager", "Emp2", "Cav1", "Clic5", "Lmo7", "S100a6", "Col4a3", "Akap5", "Cryab")
  AT1_late <- c("Sdpr", "S100a14")
  
  AT2_early <- c("Fabp5", "Lamp3", "Cd36", "Scd1", "Sftpb", "Slc34a2", "Abca3", "Sftpa1", "Egfl6", "Soat1", "Bex2", "Muc1", "Sftpc")
  AT2_late <- c("Lcn2", "Il33", "Hc", "Trf", "Lyz2", "S100g", "Lyz1")
  
  all_markers <- unique(c(AT1, AT2, Clara, Ciliated, Ciliated_novel, Calra_novel, AT1_novel, AT2_novel, AT1_early, AT1_late, AT2_early, AT2_late))
  AT1_marker <- unique(c(AT1, AT1_novel, AT1_early, AT1_late))
  AT2_marker <- unique(c(AT2, AT2_novel, AT2_early, AT2_late))
  AT1_id <- row.names(subset(fData(absolute_cds), 
                             fData(absolute_cds)$gene_short_name %in% AT1_marker)) 
  AT2_id <- row.names(subset(fData(absolute_cds), 
                             fData(absolute_cds)$gene_short_name %in% AT2_marker)) 
  
  
  add_quake_gene_all_marker <- unique(c(all_markers, as.character(quake_gene_name)))
  
  all_markers_ids <- row.names(subset(fData(absolute_cds), 
                                      fData(absolute_cds)$gene_short_name %in% all_markers))
  
  add_quake_gene_all_marker_ids <- row.names(subset(fData(absolute_cds), 
                                                    fData(absolute_cds)$gene_short_name %in% add_quake_gene_all_marker)) 
  
  quake_clustering <- c("1110017D15Rik", "1700007G11Rik", "1700009P17Rik", "4930451C15Rik", "6820408C15Rik", "BC048355", "Ccdc113", "Ccdc19", "Ccdc39", "Ccdc40", "Ccdc67", "Efhc1", "Fank1", "Foxj1", "Hs6st2", "Kndc1", "Lrrc23", "Lrriq1", "Ncs1", "Stk33", "Tekt4", "Ttc18", "Wdr16", "Kif9", "Lrrc36", "Mlf1", "Spag17", "Traf3ip1", "Enkur", "Gm6320", "Chad", "Cyp2f2", "Scgb3a2", "Hp", "Itm2a", "Krt15", "Col23a1", "Upk3a", "Cd200", "Rassf9", "Cbr2", "Nupr1", "Scgb1a1", "1810010H24Rik", "Lypd2", "Ppap2b", "Ccnd2", "Kdr", "Osgin1", "Cd24a", "Dcxr", "Ccnd1", "Mia1", "Aldh1a7", "Gsto1", "Iyd", "Aldh1a1", "C3", "Bcl6", "Tacstd2", "Clic5", "Akap5", "Ahnak", "Lmo7", "Pdpn", "Emp2", "Timp3", "Malat1", "Cav1", "Dpysl2", "Sema3a", "Col4a3", "Ager", "Ptrf", "Rtkn2", "Lgals3", "Sdpr", "Cldn18", "Samhd1", "Agrn", "Pmp22", "Qk", "Spnb2", "S100a14", "Scnn1g", "Clic3", "Msn", "Hopx", "Tspan8", "Vegfa", "Etv5", "Lamp3", "Slc34a2", "S100g", "Scd1", "Ppp1r14c", "Sftpa1", "Hc", "Egfl6", "Fabp5", "Glrx", "Cd36", "Dlk1", "Il33", "Lcn2", "Cxcl15", "Sftpb", "Slc31a1", "Bex2", "Chsy1", "Lyz1", "Sftpc", "Fabp12", "Fasn", "Soat1", "Chi3l1", "Rab27a", "Mid1ip1", "Nucb2", "Rab27b")
  known_markers <- c("Foxj1", "Scgb1a1", "Pdpn", "Ager", "Sftpb", "Sftpc")
  # figure 3 in quake paper
  quake_pca_loading <- c("Marcks", "Mdk", "Mia1", "Ranbp1", "Fam60a", "Cks1b", "Hmgn2", "Grb10", "Erdr1", "Gpc3", "Tubb5", "Myef2", "2810417H13Rik", "Kif20b", "Tspan6", "Khsrp", "Gm12070", "Slc24a5",      "Gm5506", "Cdkn1c", "Top2a", "Cldn6", "Tuba1a", "Col18a1", "Cdca3", "Cnn3", "Snrpf", "Smc2", "Otub1", "Gm6682", "Phgdh", "Ptprs", "Cdk1", "Tspan3", "Set", "H19", "Mcm3",           "Hist1h4d", "Cbx3", "Cbx5", "Clic5", "Msln", "Rtkn2", "Cav1", "Hopx", "Sdpr", "Gm11744", "Tmem213", "2200002D01Rik", "Tmem37", "S100a6", "Lgals3", "Ahnak", "Tns1", "Gprc5a",       "1190002H23Rik", "Lmo7", "Ager", "Pmp22", "Timp3", "Cryab", "Igfbp6", "Clic3", "Myo1c", "Vegfa", "AU021092", "Emp2", "Pxdn", "Meg3", "Adrb2", "Sec14l3", "Krt19", "Tinagl1",        "Aqp1", "Myl6", "Itgb6", "Agrn", "Scnn1a", "Akap5", "Cav2", "Pi4k2b", "Cbr2", "Scp2", "Cd36", "Ank3", "Pon3", "Elovl1", "Rnase4", "Idi1", "Slc34a2", "Gclc", "Rpl9", "Emb", "Lyz2", "Sftpd", "Lamp3", "S100g", "Sftpa1", "Rn4.5s", "Hc", "Tmem97", "Oat", "Lyz1", "0610007L01Rik", "Acly", "Nupr1", "Chi3l1", "Rpl23", "Baiap2l1", "Cxcl15", "Snhg11", "Car8", "H2-Aa", "H2-Ab1", "Cd200", "Scd1", "Il33", "Atp6v1c2", "Sat1", "Dbi")
  # figure 4 in quake paper
  quake_molecular_switch <- c("Lcn2", "Il33", "Hc", "Trf", "Lyz2", "S100g", "Lyz1", "Fabp5", "Lamp3", "Cd36", "Scd1", "Sftpb", "Slc34a2", "Abca3", "Sftpa1", "Egfl6", "Soat1", "Bex2", "Muc1", "Sftpc",            "Aqp5", "Pdpn", "Rtkn2", "Ager", "Emp2", "Cav1", "Clic5", "Lmo7", "S100a6", "Col4a3", "Akap5", "Cryab", "Sdpr", "S100a14")
  
  gene_list_all <- unique(c(quake_clustering, known_markers, quake_pca_loading, quake_molecular_switch)) #all genes discovered in the Quake lung paper
  gene_list_all_marker_ids <- row.names(subset(fData(absolute_cds), 
                                               fData(absolute_cds)$gene_short_name %in% gene_list_all))
  
  # use only the cells from E18.5 for the analysis to confirm the clustering for the E18.5 day in Quake paper (Fig .1): 
  absolute_cds <- estimateSizeFactors(absolute_cds)
#   absolute_cds <- estimateDispersions(absolute_cds, dispModelName='blind', cores=detectCores()) #fix the bug of estimateDispersions
  
  absolute_cds_subset <- absolute_cds[all_markers_ids, pData(absolute_cds)$Time == 'E18.5']
  absolute_cds_subset@expressionFamily <- tobit() #fix the bug of estimateDispersions
  absolute_cds_subset <- reduceDimension(absolute_cds_subset, use_irlba = F, use_vst = F) 
  absolute_cds_subset <- orderCells(absolute_cds_subset, num_paths = 2, reverse = F) #SRR1033962_0 
  row.names(SRR_cell_type) <- paste(SRR_cell_type[, 6], '0', sep = '_')
  pData(absolute_cds_subset)$Cell_type <- SRR_cell_type[colnames(absolute_cds_subset), 'putative_cell_type']
  absolute_cds_subset@expressionFamily <- negbinomial()

  pdf('quake_absolute_cds_subset.pdf')
  plot_spanning_tree(absolute_cds_subset, color_by="Cell_type", show_backbone=T, show_cell_names = F)
  dev.off()

  ########test the tree constructed with either absolute transcript counts or FPKM dataset after removing the ciliated and clara cells (use either add_quake_gene_all_marker_ids or all_AT12_markers as ordering genes)
  # remove the Clara, Ciliated cells for tree construction (use cell type assignment from Quake paper): 
  other_cell_names <- SRR_cell_type[SRR_cell_type[, 5] %in% c('ciliated', 'Clara'), 6]
  all_AT12_markers <- unique(c(AT1, AT2, AT1_novel, AT2_novel, AT1_early, AT1_late, AT2_early, AT2_late))
  all_AT12_markers_ids <- row.names(subset(fData(absolute_cds), 
                                           fData(absolute_cds)$gene_short_name %in% all_AT12_markers))
  AT12_cds_subset <- reduceDimension(standard_cds[all_markers_ids, 
                                                  !(colnames(standard_cds) %in% paste(other_cell_names, "_0", sep = ''))], use_irlba = F, use_vst = F) 
  
  # absolute_cds.quake_maker <- orderCells(absolute_cds.quake_maker, num_paths = 2, reverse = F)
  AT12_cds_subset <- orderCells(AT12_cds_subset, num_paths = 2, reverse = F) #SRR1033962_0 
  
  pdf('quake_AT12_cds_subset.pdf')
  plot_spanning_tree(AT12_cds_subset, color_by="Time", show_backbone=T, show_cell_names = F)
  dev.off()

  # use all marker genes and the quake_gene_list used before for performing the tree construction used for the later analysis: 
  AT12_cds_subset_all_gene <- reduceDimension(standard_cds[add_quake_gene_all_marker_ids, !(colnames(standard_cds) %in% paste(other_cell_names, "_0", sep = ''))], use_irlba = F, use_vst = F) 
  AT12_cds_subset_all_gene <- orderCells(AT12_cds_subset_all_gene, num_paths = 2, reverse = F) #SRR1033962_0 
  
  pData(AT12_cds_subset_all_gene)$Cell_type <- SRR_cell_type[colnames(AT12_cds_subset_all_gene), 'putative_cell_type']
  
  pdf('quake_AT12_cds_subset_all_gene.pdf')
  plot_spanning_tree(AT12_cds_subset_all_gene, color_by="Time", show_backbone=T, show_cell_names = F)
  dev.off()

  # use the tree construction information for the absolute_cds (we can also use use_for_ordering to generate the tree)
  #pData(absolute_cds) <- pData(AT12_cds_subset_all_gene)
  
  # the above tree construction approach can be better done with the use_for_ordering function: 
  test <- setOrderingFilter(standard_cds[, !(colnames(standard_cds) %in% paste(other_cell_names, "_0", sep = ''))], add_quake_gene_all_marker_ids)
  test <- reduceDimension(test[, !(colnames(test) %in% paste(other_cell_names, "_0", sep = ''))], use_irlba = F, use_vst = F) 
  test <- orderCells(test, num_paths = 2, reverse = F) #SRR1033962_0 
  plot_spanning_tree(test, color_by="Time", show_backbone=T, show_cell_names = F)
  
  # switch the state 2/3 to ensure that we get the exact plots in the paper: 
  State <- pData(AT12_cds_subset_all_gene)$State
  pData(AT12_cds_subset_all_gene)$State[State == 2] <- 3
  pData(AT12_cds_subset_all_gene)$State[State == 3] <- 2
  
  # do the same thing with transcript counts data: 
  absolute_cds@expressionFamily <- tobit()
  AT12_cds_subset_all_gene2 <- reduceDimension(absolute_cds[add_quake_gene_all_marker_ids, !(colnames(absolute_cds) %in% paste(other_cell_names, "_0", sep = ''))], max_components = 2, use_irlba=T, fun="exp", use_vst = F)  
  AT12_cds_subset_all_gene2 <- orderCells(AT12_cds_subset_all_gene2, num_paths = 2, reverse = T) #SRR1033962_0 
  
  plot_spanning_tree(AT12_cds_subset_all_gene2, color_by="Time", show_backbone=T, show_cell_names = F)
  absolute_cds@expressionFamily <- negbinomial()
  # calculate size factors as well as dispersion parameters for later analysis: 
  absolute_cds <- estimateSizeFactors(absolute_cds)
  #absolute_cds <- estimateDispersions(absolute_cds, cores=detectCores())
  rm(sizeFactors)
  AT12_cds_subset_all_gene <- estimateSizeFactors(AT12_cds_subset_all_gene)
  sizeFactors(AT12_cds_subset_all_gene) <- sizeFactors(AT12_cds_subset_all_gene)
  AT12_cds_subset_all_gene@dispFitInfo <- absolute_cds@dispFitInfo
  
  AT12_cds_subset_all_gene@dispFitInfo <- absolute_cds@dispFitInfo
  
  ##########Generate the objects for making the go enrichment barplot in figure 2 (the above reconstructed trees can be ignored?)e################
  standard_cds <- estimateSizeFactors(standard_cds, method = 'weighted-median')
  
  quake_id <- row.names(subset(fData(absolute_cds), 
                     fData(absolute_cds)$gene_short_name %in% quake_gene_name))
  quake_maker_cds <- setOrderingFilter(standard_cds, quake_id)
  quake_maker_cds <- setOrderingFilter(standard_cds, row.names(AT12_cds_subset_all_gene))
  
  quake_maker_cds <- reduceDimension(quake_maker_cds, use_irlba = F, use_vst = F) 
  quake_maker_cds <- orderCells(quake_maker_cds, num_paths = 2, reverse = F) #SRR1033962_0 
  quake_maker_cds <- orderCells(quake_maker_cds, num_paths = 2, reverse = T) #SRR1033962_0 

  plot_spanning_tree(quake_maker_cds, color_by="State", show_backbone=T, show_cell_names = F)
  
  #assign pseudotime associated data calculated with FPKM values to all datasets: USE THE ORIGINAL CELL ORDERING
  abs_AT12_cds_subset_all_gene <- absolute_cds[, colnames(AT12_cds_subset_all_gene)] #AT12_cds_subset_all_gene ONLY INCLUDES CELLS OTHER THAN THE CLARA OR CILIATED CELLS
  pData(abs_AT12_cds_subset_all_gene) <- pData(AT12_cds_subset_all_gene[, colnames(AT12_cds_subset_all_gene)])
  mc_AT12_cds_subset_all_gene <- mc_adj_cds[, colnames(AT12_cds_subset_all_gene)] #AT12_cds_subset_all_gene ONLY INCLUDES CELLS OTHER THAN THE CLARA OR CILIATED CELLS
  pData(mc_AT12_cds_subset_all_gene) <- pData(AT12_cds_subset_all_gene[, colnames(AT12_cds_subset_all_gene)])
  std_AT12_cds_subset_all_gene <- standard_cds[, colnames(AT12_cds_subset_all_gene)]
  pData(std_AT12_cds_subset_all_gene) <- pData(AT12_cds_subset_all_gene[, colnames(std_AT12_cds_subset_all_gene)])

  #estimate size and dispersion parameters
  absolute_cds <- estimateDispersions(absolute_cds)
  mc_adj_cds <- estimateSizeFactors(mc_adj_cds)
  mc_adj_cds <- estimateDispersions(mc_adj_cds)

  abs_AT12_cds_subset_all_gene <- estimateSizeFactors(abs_AT12_cds_subset_all_gene)
  abs_AT12_cds_subset_all_gene <- estimateDispersions(abs_AT12_cds_subset_all_gene)
  mc_AT12_cds_subset_all_gene <- estimateSizeFactors(mc_AT12_cds_subset_all_gene)
  mc_AT12_cds_subset_all_gene <- estimateDispersions(mc_AT12_cds_subset_all_gene)

  #save to the folder: 
#  save.image('prepare_lung_data.RData')

