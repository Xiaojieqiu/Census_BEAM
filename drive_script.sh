################################################################################################
#on the cluster, go to a directory where monocle / xacHelper located:  
################################################################################################
#install the package: 
Rscript install_packages.R ##can this be done well automatically?  

################################################################################################
#download the data necessary for the analysis: 
wget -r --no-parent --no-directories http://www.gs.washington.edu/~xqiu/proj/BEAM_analysis_data.tar.gz
#untar the file and then a directory named data which included all data necessary for reproducing the BEAM analysis will be provided: 
tar -zxvf BEAM_analysis_data.tar.gz 

#make the directories to store the figures generatated in the script 
mkdir main_figures supplementary_figure supplementary_data tmp 

################################################################################################
#each of the following analysis can be run in parallel or separately: 

################################################################
#prepare the data for the lung dataset analysis
Rscript prepare_lung_data.R 
#perform the BEAM analysis for the lung dataset
Rscript analysis_Lung_data.R 
#perform the analysis for the spike-in free algorithm
Rscript spikein_free_algorithm_benchmark.R 
#perform benchmark analysis
Rscript benchmark_analysis.R 
################################################################

################################################################
#perform the analysis for the HSMM dataset
Rscript analysis_HSMM_data.R 
################################################################

################################################################
#perform the analysis for the UMI dataset 
Rscript analysis_UMI_data.R 
################################################################

################################################################
#perform the analysis for the Shalek dataset
Rscript analysis_shalek_data.R 
################################################################

################################################################
#perform the analysis for making supplementary figures 
Rscript prepare_other_supplementary_data.R 
################################################################


################################################################################################
#when the objects are generated from the above run, the following scripts 
#can be used to generate figures in the manuscript
################################################################################################

################################################################################################
#the following script can be run in parallel or separately: 
#generate the figures: 
Rscript gen_lung_figures.R
Rscript gen_shalek_figures.R
Rscript gen_supplementary_figure.R 

################################################################################################
#done 