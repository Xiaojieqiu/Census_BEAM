################################################################################################
#download the data necessary for the analysis: 
wget -r --no-parent --no-directories http://www.gs.washington.edu/~xqiu/proj/BEAM/

################################################################################################
#each of the following analysis can be run in parallel or separately: 

################################################################
#prepare the data for the lung dataset analysis
R CMD prepare_lung_data.R 
#perform the BEAM analysis for the lung dataset
R CMD analysis_Lung_data.R 
#perform the analysis for the spike-in free algorithm
R CMD spike-in_free_algorithm_benchmark.R 
#perform benchmark analysis
R CMD benchmark_analysis.R 
################################################################

################################################################
#perform the analysis for the HSMM dataset
R CMD analysis_HSMM_data.R 
################################################################

################################################################
#perform the analysis for the UMI dataset 
R CMD analysis_UMI_data.R 
################################################################

################################################################
#perform the analysis for the Shalek dataset
R CMD analysis_shalek_data.R 
################################################################

################################################################
#perform the analysis for making supplementary figures 
R CMD prepare_other_supplementary_data.R 
################################################################


################################################################################################
#when the objects are generated from the above run, the following scripts 
#can be used to generate figures in the manuscript
################################################################################################

################################################################################################
#the following script can be run in parallel or separately: 
#generate the figures: 
R CMD gen_lung_figures.R
R CMD gen_shalek_figures.R
R CMD gen_supplementary_figure.R 

################################################################################################
#done 