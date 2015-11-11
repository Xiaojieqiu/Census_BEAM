#download the data necessary for the analysis: 
wget -r --no-parent --no-directories http://www.gs.washington.edu/~xqiu/proj/BEAM/

#prepare the data for the lung analysis
R CMD prepare_lung_data.R 

#parallel the running of all the analysis: 
R CMD analysis_Lung_data.R 
R CMD analysis_HSMM_data.R 
R CMD analysis_UMI_data.R 
R CMD analysis_shalek_data.R 
R CMD prepare_other_supplementary_data.R 
R CMD spikein_free_algorithm_benchmark.R 
R CMD benchmark_analysis.R 

#generate the figures: 
R CMD gen_lung_figures.R
R CMD gen_shalek_figures.R
R CMD gen_supplementary_figure.R 

#done 