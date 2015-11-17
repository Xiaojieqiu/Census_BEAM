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
mkdir main_figures supplementary_figures supplementary_data tmp RData


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
