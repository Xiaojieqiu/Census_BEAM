#########################################################
# run_all_analyses.sh
# submitted with manuscript, "Branched single-cell trajectories reveal regulators of cell fate decisions"
#
# Autmates all analysis and figure panel generation for the manuscript. Please see README.md for information
# about dependencies and platform requirements, paticularly if you only plan to run certain analyses.
#
# No arguments required, simply run:
# bash run_all_analyses.sh
#########################################################


#########################################################
# Pre-analysis setup
#########################################################

## Install all R packages
Rscript install_packages.R

## Download data required for the analysis and decompress 
wget http://www.gs.washington.edu/~xqiu/proj/BEAM_analysis_data.tar.gz
tar -zxvf BEAM_analysis_data.tar.gz 
rm BEAM_analysis_data.tar.gz

## Make required directories to store figures and data that will be generated
mkdir main_figures supplementary_figures supplementary_data tmp RData nbt_2nd_sub_reviewers

## Download the HSMMSingleCell package to local directory
wget http://bioconductor.org/packages/release/data/experiment/src/contrib/HSMMSingleCell_0.104.0.tar.gz
wget https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_1.0.1.tar.gz 

#########################################################
# Analysis (WARNING: see README.md for dependency
# information if running individual commands)
#
# The following section includes scripts that analyze data and save .RData files
# for figure generation using other scripts
#########################################################

## Pre-process lung data before downstream analysis
Rscript prepare_lung_data.R 

## Perform BEAM analysis for the lung dataset
Rscript analysis_lung_data.R 

## Sample m,c space to calculate the value for the optimization function used in the spike-in free recovery algorithm
Rscript spikein_free_algorithm_sampling.R

## Perform DEG analysis on DESeq, DESeq2, edgeR, SCDE and monocle to benchmark the DEG test performance
Rscript deg_benchmark_analysis.R 

## Perform all analysis for the HSMM dataset
Rscript analysis_HSMM_data.R 

# Perform all analysis for the UMI dataset 
Rscript analysis_UMI_data.R 

## Perform analysis for the Shalek dataset
Rscript analysis_shalek_data.R 

# Perform goodness of fit analysis on read counts or transcript counts data 
Rscript analysis_distribution_fitting.R

#########################################################
# Figure Generation (WARNING: see README.md for dependency
# information if running individual commands)
#
# The following section contains scripts that generate
# figure panels from the paper.
#########################################################

## Generate figures based on Lung data
Rscript gen_lung_figures.R

## Generate figures based on Shalek data
Rscript gen_shalek_figures.R

## Generate supplementary figures
Rscript gen_supplementary_figure.R 

###### The following are scripts developed for the resubmission #####
##
#########################################################
# New figures added in manuscript and rebuttal figures for resubmission 
# All the scripts below are depend on previous figures 
# excepting deg_benchmark_analysis_HSMM_bulk.R and roc_curves.R 
#
# The following section contains scripts that generates
# figure panels from the paper.
#########################################################

## Perform cell dowsampling for testing robustness of BEAM 
#(WARNING: this script takes about 3 days to finish on a cluster with 64 cores)
Rscript analysis_cell_downsampling.R  

## DEG bencmark for HSMM data using bulk RNA-seq as gold standard
Rscript deg_benchmark_analysis_HSMM_bulk.R       

## This R markdown file can be used to generate mode_of_isoform_rebuttal.html   
## Warning: run this script in Rstudio from R command line  
mode_of_isoform_rebuttal.Rmd      

## function for perform ROC analysis, required for deg_benchmark_analysis_HSMM_bulk  
Rscript roc_curves.R

## Perform the benchmark between beam test and two-group test
Rscript BEAM_test.R                   

## Perform dimension reduction using LLE, diffusion and tSNE 
Rscript dimension_reduction.R                    

## Perform branch time point analysis (for shalek data, in particular)
Rscript branchTimePoint.R             

## Perform spike-in free recovery using UMI data 
Rscript umi_normalization.R

## Assess the impact of number of genes detected on the dimension reduction 
Rscript gen_impact_num_genes_detected_figures.R  

## Perform read count downsampling 
Rscript prepare_lung_downsampling_data.R  

## Perform clustering analysis for the unstimulated cells 
Rscript subpopulation.R

## Compare with other existing software (MAST)
Rscript cmpr_three_packages.R         

## Making downsampled figure from lung data 
#(WARNING: this script takes about 2 days to finish on a cluster with 64 cores)
Rscript gen_lung_downsampling_figures.R          

## Script to generate other figures not in scripts above 
Rscript remaining_review_figures.R        



