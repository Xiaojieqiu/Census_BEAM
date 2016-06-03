#this script tries to perform the downsampling of cells for showing the recovery algorithm

library(argparse)
library(tidyr)
library(grid)
library(gridExtra)
# library(monocle)
library(devtools)
load_all('~/Projects/monocle-dev')
library(xacHelper)
library(plyr)
library(data.table)
library(MASS)
library(modeest)
library(dplyr)
library(matrixStats)


####################################################
# Load in data from initial lung analysis (mostly need marker genes and a few other sets of genes)
####################################################
load("RData/prepare_lung_data.RData")
valid_cells = c("SRR1033854_thout_0", "SRR1033855_thout_0", "SRR1033856_thout_0", "SRR1033859_thout_0", "SRR1033860_thout_0", "SRR1033861_thout_0", "SRR1033862_thout_0", "SRR1033863_thout_0", "SRR1033864_thout_0", "SRR1033867_thout_0", "SRR1033869_thout_0", "SRR1033871_thout_0", "SRR1033872_thout_0", "SRR1033874_thout_0", "SRR1033875_thout_0", "SRR1033876_thout_0", "SRR1033877_thout_0", "SRR1033878_thout_0", "SRR1033879_thout_0", "SRR1033880_thout_0", "SRR1033881_thout_0", "SRR1033882_thout_0", "SRR1033883_thout_0", "SRR1033884_thout_0", "SRR1033885_thout_0", "SRR1033886_thout_0", "SRR1033887_thout_0", "SRR1033888_thout_0", "SRR1033889_thout_0", "SRR1033890_thout_0", "SRR1033891_thout_0", "SRR1033892_thout_0", "SRR1033893_thout_0", "SRR1033894_thout_0", "SRR1033895_thout_0", "SRR1033896_thout_0", "SRR1033897_thout_0", "SRR1033898_thout_0", "SRR1033899_thout_0", "SRR1033900_thout_0", "SRR1033901_thout_0", "SRR1033902_thout_0", "SRR1033904_thout_0", "SRR1033905_thout_0", "SRR1033906_thout_0", "SRR1033907_thout_0", "SRR1033908_thout_0", "SRR1033909_thout_0", "SRR1033910_thout_0", "SRR1033911_thout_0", "SRR1033914_thout_0", "SRR1033915_thout_0", "SRR1033916_thout_0", "SRR1033917_thout_0", "SRR1033919_thout_0", "SRR1033921_thout_0", "SRR1033922_thout_0", "SRR1033923_thout_0", "SRR1033926_thout_0", "SRR1033927_thout_0", "SRR1033928_thout_0", "SRR1033929_thout_0", "SRR1033931_thout_0", "SRR1033932_thout_0", "SRR1033933_thout_0", "SRR1033934_thout_0", "SRR1033936_thout_0", "SRR1033937_thout_0", "SRR1033938_thout_0", "SRR1033939_thout_0", "SRR1033940_thout_0", "SRR1033941_thout_0", "SRR1033942_thout_0", "SRR1033943_thout_0", "SRR1033945_thout_0", "SRR1033946_thout_0", "SRR1033947_thout_0", "SRR1033948_thout_0", "SRR1033949_thout_0", "SRR1033950_thout_0", "SRR1033951_thout_0", "SRR1033952_thout_0", "SRR1033953_thout_0", "SRR1033954_thout_0", "SRR1033955_thout_0", "SRR1033956_thout_0", "SRR1033957_thout_0", "SRR1033958_thout_0", "SRR1033959_thout_0", "SRR1033960_thout_0", "SRR1033961_thout_0", "SRR1033962_thout_0", "SRR1033963_thout_0", "SRR1033964_thout_0", "SRR1033965_thout_0", "SRR1033966_thout_0", "SRR1033967_thout_0", "SRR1033968_thout_0", "SRR1033969_thout_0", "SRR1033970_thout_0", "SRR1033971_thout_0", "SRR1033972_thout_0", "SRR1033973_thout_0", "SRR1033974_thout_0", "SRR1033975_thout_0", "SRR1033976_thout_0", "SRR1033977_thout_0", "SRR1033978_thout_0", "SRR1033979_thout_0", "SRR1033980_thout_0", "SRR1033981_thout_0", "SRR1033982_thout_0", "SRR1033983_thout_0", "SRR1033984_thout_0", "SRR1033985_thout_0", "SRR1033986_thout_0", "SRR1033987_thout_0", "SRR1033988_thout_0", "SRR1033989_thout_0", "SRR1033990_thout_0", "SRR1033991_thout_0", "SRR1033992_thout_0", "SRR1033993_thout_0", "SRR1033994_thout_0", "SRR1033995_thout_0", "SRR1033996_thout_0", "SRR1033997_thout_0", "SRR1033998_thout_0", "SRR1033999_thout_0", "SRR1034000_thout_0", "SRR1034001_thout_0", "SRR1034002_thout_0", "SRR1034003_thout_0", "SRR1034004_thout_0", "SRR1034005_thout_0", "SRR1034006_thout_0", "SRR1034007_thout_0", "SRR1034008_thout_0", "SRR1034009_thout_0", "SRR1034010_thout_0", "SRR1034011_thout_0", "SRR1034012_thout_0", "SRR1034013_thout_0", "SRR1034014_thout_0", "SRR1034015_thout_0", "SRR1034016_thout_0", "SRR1034017_thout_0", "SRR1034018_thout_0", "SRR1034019_thout_0", "SRR1034020_thout_0", "SRR1034021_thout_0", "SRR1034022_thout_0", "SRR1034023_thout_0", "SRR1034024_thout_0", "SRR1034025_thout_0", "SRR1034026_thout_0", "SRR1034027_thout_0", "SRR1034028_thout_0", "SRR1034029_thout_0", "SRR1034030_thout_0", "SRR1034031_thout_0", "SRR1034032_thout_0", "SRR1034033_thout_0", "SRR1034034_thout_0", "SRR1034035_thout_0", "SRR1034036_thout_0", "SRR1034037_thout_0", "SRR1034038_thout_0", "SRR1034039_thout_0", "SRR1034040_thout_0", "SRR1034041_thout_0", "SRR1034042_thout_0", "SRR1034043_thout_0", "SRR1034044_thout_0", "SRR1034045_thout_0", "SRR1034046_thout_0", "SRR1034047_thout_0", "SRR1034048_thout_0", "SRR1034049_thout_0", "SRR1034050_thout_0", "SRR1034051_thout_0", "SRR1034052_thout_0", "SRR1034053_thout_0")

####################################################
# Helper functions for analysis
####################################################

# Wraps conversion to transcript counts while also returning the m and c statistics from each run
get_recovered_transcript_counts_with_stats = function(cds, kb_intercept = 2.62) {
    closeAllConnections()

    set.seed(1:(300*ncol(cds))) #set the seed so that we get the same result every time
    normalized_cds_stats = relative2abs(cds, estimate_t(exprs(cds)), cores=detectCores(), return_all=T) #kb_intercept = kb_intercept,  also get m and c from return="all"
    closeAllConnections()

    # Update the expression matrix
    exprs(cds) = as.matrix(normalized_cds_stats$norm_cds)

    # Update expression family and estimate size factors, etc.
    cds@expressionFamily = negbinomial()
    cds = estimateSizeFactors(cds)
    cds = estimateDispersions(cds)

    # Put the CDS in the list and return with m and c parameters
    normalized_cds_stats$norm_cds = cds
    return(normalized_cds_stats[c("norm_cds", "kb_intercept", "kb_slope")])
}

# Wraps code borrowed from lung analysis to do spike-in regression
convert_monocle_cds_to_transcript_counts_regression <- function(cds) {
  # Prepare spike DF and their expression values
  spike_names <- row.names(cds)[grepl("ERCC", row.names(cds))]
  
  ercc_controls <- cds[spike_names,]
  mean_spike_fpkms <- rowMeans(exprs(ercc_controls))
  
  spike_df <- cbind(spike_df, mean_spike_fpkms[row.names(spike_df)])
  pData(ercc_controls)$Cell = row.names(pData(ercc_controls))
  
  # Perform regression on spike DF
  molModels <- esApply(ercc_controls, 2, function(cell_exprs, input.ERCC.annotation, volume = volume, dilution = dilution, select_above_thresh = T, capture_efficiency) {
    
    spike_df <- input.ERCC.annotation 
    spike_df <- cbind(spike_df, cell_exprs[row.names(spike_df)])
    colnames(spike_df)[length(colnames(spike_df))] <- "FPKM"
    spike_df$numMolecules <- spike_df$conc_attomoles_ul_Mix1*(volume*10^(-3)*dilution*10^(-18)*6.02214179*10^(23))
    if(!is.na(capture_efficiency))
      spike_df$numMolecules <- spike_df$numMolecules * capture_efficiency
    spike_df$rounded_numMolecules <- round(spike_df$conc_attomoles_ul_Mix1*(volume*10^(-3)*dilution*10^(-18)*6.02214179*10^(23)))
    
    if(select_above_thresh){
      spike_df <- subset(spike_df, conc_attomoles_ul_Mix1 >= 800)
      spike_df <- subset(spike_df, FPKM >= 1e-10)    
    }
    spike_df$log_fpkm <- log10(spike_df$FPKM) 
    spike_df$log_numMolecules <- log10(spike_df$numMolecules)
    
    molModel <- tryCatch({
      molModel <- rlm(log_numMolecules ~ log_fpkm, data=spike_df)
      molModel
    }, 
    error = function(e) { print(e); NULL })
    molModel
  }, spike_df, volume = 10, dilution = 1/ 40000, capture_efficiency = 1)
  
  
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
  split(exprs(cds), rep(1:ncol(exprs(cds)), each = nrow(exprs(cds)))), 
  molModels)
  
  row.names(norm_fpkms) <- row.names(exprs(cds))
  colnames(norm_fpkms) <- colnames(exprs(cds))
  
  # Now generate a new CellDataSet that uses the absolute transcript counts
  fpkm_matrix_abs <- norm_fpkms
  colnames(fpkm_matrix_abs) <- colnames(exprs(cds))
  row.names(fpkm_matrix_abs) <- row.names(exprs(cds))
  
  pd <- new("AnnotatedDataFrame", data = pData(cds)[colnames(fpkm_matrix_abs),])
  fd <- new("AnnotatedDataFrame", data = fData(cds)[rownames(fpkm_matrix_abs),])
  
  absolute_cds <- newCellDataSet(fpkm_matrix_abs, 
                                 phenoData = pd, 
                                 featureData = fd, 
                                 expressionFamily=negbinomial(), 
                                 lowerDetectionLimit=1)
  
  return( absolute_cds )
}

set.seed(5)
MIN_PROPORTION = 0.1
MAX_PROPORTION = 1
STEP = 0.1 #0.1
REPS_PER = 3
EXTRA_PROPORTIONS = c(0.85, 0.95) 

downsampled_proportions = sort(rep(c(seq(MIN_PROPORTION, MAX_PROPORTION, by=STEP), EXTRA_PROPORTIONS) , REPS_PER))
names(downsampled_proportions) = downsampled_proportions  # will tie CDS objects to proportion for later
cds_downsampled_cells = lapply(downsampled_proportions, function(x) { standard_cds[, sample(ncol(standard_cds), round(ncol(standard_cds) * x))] })

####################################################
# Convert the CDS to transcript counts using normalization algorithm
####################################################

# Convert the dataset to transcript counts and keep info about m and c values
original_depth_cds_conversion_stats = get_recovered_transcript_counts_with_stats(standard_cds)
original_depth_cds_transcript_counts = original_depth_cds_conversion_stats$norm_cds

# #test the consistency
# original_depth_cds_conversion_stats2 = get_recovered_transcript_counts_with_stats(original_depth_cds, input_data$original_isoform_matrix)
# original_depth_cds_transcript_counts2 = original_depth_cds_conversion_stats2$norm_cds

closeAllConnections()
cds_to_compare_conversion_stats = mapply(get_recovered_transcript_counts_with_stats, cds_downsampled_cells)
closeAllConnections()
cds_to_compare_transcript_counts = cds_to_compare_conversion_stats[1, ]
cds_to_compare_conversion_m_values = cds_to_compare_conversion_stats[2, ]
cds_to_compare_conversion_c_values = cds_to_compare_conversion_stats[3, ]

# Also convert using spike-in regression for comparison (no stats need to be collected here)
original_depth_cds_transcript_counts_regression = convert_monocle_cds_to_transcript_counts_regression(standard_cds)
cds_to_compare_transcript_counts_regression = mapply(convert_monocle_cds_to_transcript_counts_regression, cds_downsampled_cells)

#########################################
# Save results for figure generation script
#########################################
save.image("RData/prepare_downsampling_data.RData")
