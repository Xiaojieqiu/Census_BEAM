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
# Load in data from initial lung analysis (mostly need marker genes and a few other sets of genes)
####################################################
load("RData/prepare_lung_data.RData")
valid_cells = c("SRR1033854_thout_0", "SRR1033855_thout_0", "SRR1033856_thout_0", "SRR1033859_thout_0", "SRR1033860_thout_0", "SRR1033861_thout_0", "SRR1033862_thout_0", "SRR1033863_thout_0", "SRR1033864_thout_0", "SRR1033867_thout_0", "SRR1033869_thout_0", "SRR1033871_thout_0", "SRR1033872_thout_0", "SRR1033874_thout_0", "SRR1033875_thout_0", "SRR1033876_thout_0", "SRR1033877_thout_0", "SRR1033878_thout_0", "SRR1033879_thout_0", "SRR1033880_thout_0", "SRR1033881_thout_0", "SRR1033882_thout_0", "SRR1033883_thout_0", "SRR1033884_thout_0", "SRR1033885_thout_0", "SRR1033886_thout_0", "SRR1033887_thout_0", "SRR1033888_thout_0", "SRR1033889_thout_0", "SRR1033890_thout_0", "SRR1033891_thout_0", "SRR1033892_thout_0", "SRR1033893_thout_0", "SRR1033894_thout_0", "SRR1033895_thout_0", "SRR1033896_thout_0", "SRR1033897_thout_0", "SRR1033898_thout_0", "SRR1033899_thout_0", "SRR1033900_thout_0", "SRR1033901_thout_0", "SRR1033902_thout_0", "SRR1033904_thout_0", "SRR1033905_thout_0", "SRR1033906_thout_0", "SRR1033907_thout_0", "SRR1033908_thout_0", "SRR1033909_thout_0", "SRR1033910_thout_0", "SRR1033911_thout_0", "SRR1033914_thout_0", "SRR1033915_thout_0", "SRR1033916_thout_0", "SRR1033917_thout_0", "SRR1033919_thout_0", "SRR1033921_thout_0", "SRR1033922_thout_0", "SRR1033923_thout_0", "SRR1033926_thout_0", "SRR1033927_thout_0", "SRR1033928_thout_0", "SRR1033929_thout_0", "SRR1033931_thout_0", "SRR1033932_thout_0", "SRR1033933_thout_0", "SRR1033934_thout_0", "SRR1033936_thout_0", "SRR1033937_thout_0", "SRR1033938_thout_0", "SRR1033939_thout_0", "SRR1033940_thout_0", "SRR1033941_thout_0", "SRR1033942_thout_0", "SRR1033943_thout_0", "SRR1033945_thout_0", "SRR1033946_thout_0", "SRR1033947_thout_0", "SRR1033948_thout_0", "SRR1033949_thout_0", "SRR1033950_thout_0", "SRR1033951_thout_0", "SRR1033952_thout_0", "SRR1033953_thout_0", "SRR1033954_thout_0", "SRR1033955_thout_0", "SRR1033956_thout_0", "SRR1033957_thout_0", "SRR1033958_thout_0", "SRR1033959_thout_0", "SRR1033960_thout_0", "SRR1033961_thout_0", "SRR1033962_thout_0", "SRR1033963_thout_0", "SRR1033964_thout_0", "SRR1033965_thout_0", "SRR1033966_thout_0", "SRR1033967_thout_0", "SRR1033968_thout_0", "SRR1033969_thout_0", "SRR1033970_thout_0", "SRR1033971_thout_0", "SRR1033972_thout_0", "SRR1033973_thout_0", "SRR1033974_thout_0", "SRR1033975_thout_0", "SRR1033976_thout_0", "SRR1033977_thout_0", "SRR1033978_thout_0", "SRR1033979_thout_0", "SRR1033980_thout_0", "SRR1033981_thout_0", "SRR1033982_thout_0", "SRR1033983_thout_0", "SRR1033984_thout_0", "SRR1033985_thout_0", "SRR1033986_thout_0", "SRR1033987_thout_0", "SRR1033988_thout_0", "SRR1033989_thout_0", "SRR1033990_thout_0", "SRR1033991_thout_0", "SRR1033992_thout_0", "SRR1033993_thout_0", "SRR1033994_thout_0", "SRR1033995_thout_0", "SRR1033996_thout_0", "SRR1033997_thout_0", "SRR1033998_thout_0", "SRR1033999_thout_0", "SRR1034000_thout_0", "SRR1034001_thout_0", "SRR1034002_thout_0", "SRR1034003_thout_0", "SRR1034004_thout_0", "SRR1034005_thout_0", "SRR1034006_thout_0", "SRR1034007_thout_0", "SRR1034008_thout_0", "SRR1034009_thout_0", "SRR1034010_thout_0", "SRR1034011_thout_0", "SRR1034012_thout_0", "SRR1034013_thout_0", "SRR1034014_thout_0", "SRR1034015_thout_0", "SRR1034016_thout_0", "SRR1034017_thout_0", "SRR1034018_thout_0", "SRR1034019_thout_0", "SRR1034020_thout_0", "SRR1034021_thout_0", "SRR1034022_thout_0", "SRR1034023_thout_0", "SRR1034024_thout_0", "SRR1034025_thout_0", "SRR1034026_thout_0", "SRR1034027_thout_0", "SRR1034028_thout_0", "SRR1034029_thout_0", "SRR1034030_thout_0", "SRR1034031_thout_0", "SRR1034032_thout_0", "SRR1034033_thout_0", "SRR1034034_thout_0", "SRR1034035_thout_0", "SRR1034036_thout_0", "SRR1034037_thout_0", "SRR1034038_thout_0", "SRR1034039_thout_0", "SRR1034040_thout_0", "SRR1034041_thout_0", "SRR1034042_thout_0", "SRR1034043_thout_0", "SRR1034044_thout_0", "SRR1034045_thout_0", "SRR1034046_thout_0", "SRR1034047_thout_0", "SRR1034048_thout_0", "SRR1034049_thout_0", "SRR1034050_thout_0", "SRR1034051_thout_0", "SRR1034052_thout_0", "SRR1034053_thout_0")

####################################################
# Helper functions for analysis
####################################################

# Wraps generation of CDS objects for use with lapply
get_monocle_cds <- function(fpkm_matrix, pd, fd, valid_cells) {

    row.names(fpkm_matrix) = fpkm_matrix[, 1]
    row.names(pd) = pd[, 1]
    row.names(fd) = fd[, 1]

    # Filter to a set of valid cells
    fpkm_matrix = fpkm_matrix[, valid_cells]
    pd = pd[valid_cells, ]

    # Normalize the gene short names to all uppercase
    fd$gene_short_name = toupper(fd$gene_short_name)

    # Make sure the order of each dataframe is the same
    fpkm_matrix = fpkm_matrix[, row.names(pd)]
    fd = fd[row.names(fpkm_matrix), ]

    # Initialize data structures for monocle and return HSMM
    pd = new("AnnotatedDataFrame", data = pd)
    fd = new("AnnotatedDataFrame", data = fd)
    cds = new("CellDataSet", exprs = as.matrix(fpkm_matrix), phenoData = pd, featureData = fd, expressionFamily=tobit(), lowerDetectionLimit=0.1)

    return(cds)
}

# Wraps conversion to transcript counts while also returning the m and c statistics from each run
get_recovered_transcript_counts_with_stats = function(cds, isoform_matrix) {
    closeAllConnections()

    normalized_cds_stats = monocle::relative2abs(cds, monocle::estimate_t(isoform_matrix), cores=detectCores(), return_all=T) # also get m and c from return="all"
    closeAllConnections()

    # Update the expression matrix
    exprs(cds) = as.matrix(normalized_cds_stats$norm_cds)

    # Update expression family and estimate size factors, etc.
    cds@expressionFamily = negbinomial()
    cds = estimateSizeFactors(cds)
    cds = estimateDispersions(cds)

    # Put the CDS in the list and return with m and c parameters
    normalized_cds_stats$norm_cds = cds
    return(normalized_cds_stats[c("norm_cds", "m", "c")])
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


####################################################
# Load and process sample metadata
####################################################
sra_table = read.delim("data/Quake_data/quake_lung/SraRunTable.txt")
sra_table$Run = paste(sra_table$Run, "_thout_0", sep="")  # match sample names in other files

sample_sheet = read.delim("data/Quake_data/quake_lung/standard_normalized_out/samples.table", row.names="sample_id")
row.names(sample_sheet) = stringr::str_replace(row.names(sample_sheet), "_0", "_thout_0") # match sample names in other files

# Merge sample data to generate sample sheet
sample_sheet = merge(sample_sheet, sra_table, by.x="row.names", by.y="Run")
row.names(sample_sheet) = sample_sheet$Row.names

# Switch sample name conventions to short names
sample_sheet$Time = NULL
sample_sheet$Time[sample_sheet$age == "Embryonic day 14.5"] <- "E14.5"
sample_sheet$Time[sample_sheet$age == "Embryonic day 16.5"] <- "E16.5"
sample_sheet$Time[sample_sheet$age == "Embryonic day 18.5"] <- "E18.5"
sample_sheet$Time[sample_sheet$age == "post natal day 107"] <- "Adult"


####################################################
# Define files to read in and reorder according to depth
####################################################
original_expression_matrix = "data/Quake_data/quake_lung/downsampled_normalized_out/7645000_cuffnorm_output/genes.fpkm_table"
original_isoform_matrix = "data/Quake_data/quake_lung/downsampled_normalized_out/7645000_cuffnorm_output/isoforms.fpkm_table"
gene_annotation = "data/Quake_data/quake_lung/downsampled_normalized_out/7645000_cuffnorm_output/genes.attr_table"

expression_matrices_to_compare = Sys.glob("data/Quake_data/quake_lung/downsampled_normalized_out/*_cuffnorm_output/genes.fpkm_table")
isoform_matrices_to_compare = Sys.glob("data/Quake_data/quake_lung/downsampled_normalized_out/*_cuffnorm_output/isoforms.fpkm_table")
downsampled_depths = stringr::str_match(expression_matrices_to_compare, "([0-9]+)_cuffnorm_output")[,2]

order_by_depth = order(as.numeric(as.character(downsampled_depths)))

downsampled_depths = downsampled_depths[order_by_depth]
expression_matrices_to_compare = expression_matrices_to_compare[order_by_depth]
isoform_matrices_to_compare = isoform_matrices_to_compare[order_by_depth]


####################################################
# Read in data files
####################################################
# Read in input files
input_files = c("original_expression_matrix"=original_expression_matrix, "original_isoform_matrix"=original_isoform_matrix, "gene_annotation"=gene_annotation)
input_data=lapply(input_files, function(i){read.table(i, sep="\t", header=TRUE)})
input_data$sample_metadata = sample_sheet  # add in the sample sheet

# Read in expression matrices not at original depth that you want to compare
names(expression_matrices_to_compare) = downsampled_depths
input_data$expression_matrices_to_compare = lapply(expression_matrices_to_compare, function(i){read.table(i, sep="\t", header=TRUE)})

# Also read in isoform expression matrices for transcript count estimations
input_data$isoform_matrices_to_compare = lapply(isoform_matrices_to_compare, function(i){read.table(i, sep="\t", header=TRUE)})


####################################################
# Generate CDS objects for data analysis and preprocess data (only including valid cells)
####################################################

# Generate CDS objects
original_depth_cds=get_monocle_cds(input_data$original_expression_matrix, input_data$sample_metadata, input_data$gene_annotation, valid_cells)
cds_to_compare=lapply(input_data$expression_matrices_to_compare, get_monocle_cds, input_data$sample_metadata, input_data$gene_annotation, valid_cells)

# Ensure that the original dataframe and new ones have same gene ordering of genes
cds_to_compare = lapply(cds_to_compare, function(cds) { cds[row.names(original_depth_cds), ]})

# Prepare the isoform matrices
## Set row names
row.names(input_data$original_isoform_matrix) = input_data$original_isoform_matrix$tracking_id
input_data$isoform_matrices_to_compare = lapply(input_data$isoform_matrices_to_compare, function(expression_matrix) { row.names(expression_matrix) = input_data$original_isoform_matrix$tracking_id; return(expression_matrix)})

## Make sure they have the same cells as the expression matrices
input_data$original_isoform_matrix = input_data$original_isoform_matrix[, colnames(original_depth_cds)]
input_data$isoform_matrices_to_compare = lapply(input_data$isoform_matrices_to_compare, function(expression_matrix) { return (expression_matrix[, colnames(original_depth_cds)])})


####################################################
# Convert the CDS to transcript counts using normalization algorithm
####################################################

# Convert the dataset to transcript counts and keep info about m and c values
original_depth_cds_conversion_stats = get_recovered_transcript_counts_with_stats(original_depth_cds, input_data$original_isoform_matrix)
original_depth_cds_transcript_counts = original_depth_cds_conversion_stats$norm_cds

closeAllConnections()
cds_to_compare_conversion_stats = mapply(get_recovered_transcript_counts_with_stats, cds_to_compare, input_data$isoform_matrices_to_compare)
closeAllConnections()
cds_to_compare_transcript_counts = cds_to_compare_conversion_stats[1, ]
cds_to_compare_conversion_m_values = cds_to_compare_conversion_stats[2, ]
cds_to_compare_conversion_c_values = cds_to_compare_conversion_stats[3, ]

# Also convert using spike-in regression for comparison (no stats need to be collected here)
original_depth_cds_transcript_counts_regression = convert_monocle_cds_to_transcript_counts_regression(original_depth_cds)
cds_to_compare_transcript_counts_regression = mapply(convert_monocle_cds_to_transcript_counts_regression, cds_to_compare)

#########################################
# Reduce dimension and order cells
#########################################
## use all marker genes and the quake_gene_list used before for performing the tree construction used for the later analysis: 
get_root_state <- function(cds, root_cell ) {
  pData(cds[, root_cell])$State
}

reduce_dimension_and_order_cells_transcript_counts = function(cds, root_cell) {
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  cds = reduceDimension(cds[add_quake_gene_all_marker_ids, ], use_irlba = F, use_vst=T, method="ICA", scaling=F, pseudo_expr=0) 
  cds = orderCells(cds, num_paths = 2, reverse = T)
  root_state = get_root_state(cds, root_cell)  # get state that was assigned to root cell
  cds = orderCells(cds, root_state, num_paths = 2, reverse = T)  # now set root cell state as root state
  return(cds)
}

# Reduce dimension and order all transcript counts matrices 
original_depth_cds_transcript_counts_ordered = reduce_dimension_and_order_cells_transcript_counts(original_depth_cds_transcript_counts, "SRR1033943_thout_0")
original_depth_cds_transcript_counts_regression_ordered = reduce_dimension_and_order_cells_transcript_counts(original_depth_cds_transcript_counts_regression, "SRR1033943_thout_0")

cds_to_compare_transcript_counts_ordered = lapply(cds_to_compare_transcript_counts, reduce_dimension_and_order_cells_transcript_counts, "SRR1033943_thout_0") 
cds_to_compare_transcript_counts_regression_ordered = lapply(cds_to_compare_transcript_counts_regression, reduce_dimension_and_order_cells_transcript_counts, "SRR1033943_thout_0") 

cor(pData(original_depth_cds_transcript_counts_regression_ordered)$Pseudotime, pData(cds_to_compare_transcript_counts_regression_ordered$"7645000")$Pseudotime)
cor(pData(original_depth_cds_transcript_counts_ordered)$Pseudotime, pData(cds_to_compare_transcript_counts_ordered$"7645000")$Pseudotime)

#########################################
# Save results for figure generation script
#########################################
save.image("RData/prepare_downsampling_data.RData")