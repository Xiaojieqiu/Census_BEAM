################################################################################################
#on the cluster, go to a directory where monocle / xacHelper located:  
################################################################################################
#install the package: 
##make sure you install all the packages listed in the script before you run all the downstream 
#analysis: 
Rscript install_packages.R ##can this be done well automatically?  

################################################################################################
#download the data necessary for the analysis: 
wget -r http://www.gs.washington.edu/~xqiu/proj/BEAM_analysis_data.tar.gz
#untar the file and then a directory named data which included all data necessary for reproducing 
#the BEAM analysis will be provided: 
tar -zxvf BEAM_analysis_data.tar.gz 

#make the directories to store the figures generatated in the script 
mkdir main_figures supplementary_figures supplementary_data tmp 

################################################################################################
#each of the following analysis can be run in parallel or separately: 

################################################################
#prepare the data for the lung dataset analysis
Rscript prepare_lung_data.R 

##notee that the following dependent on the prepare_lung_data.R but each of them are independent 
#and can be parallelized

#perform the BEAM analysis for the lung dataset
Rscript analysis_lung_data.R 
#perform the analysis for the spike-in free algorithm
Rscript spikein_free_algorithm_benchmark.R  ###test whether or not we can parallelize this? and make sure it will be avoided when testing on cluster
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
Rscript analysis_other_supplementary_data.R 
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