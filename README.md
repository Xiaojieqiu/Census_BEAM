#Analysis for "Branched single-cell trajectories reveal regulators of cell fate decisions"

## Introduction
This distribution automates generation of all analysis and figure generation performed for our manuscript, "Branched single-cell trajectories reveal regulators of cell fate decisions". In addition to the scripts that recreate the analysis, this distribution additionally includes:

1. An alpha version of Monocle 2.0 (file: "monocle_1.99.0.tar.gz"), which will be released publicly as the next major release of Monocle through both GitHub and Bioconductor as an R package.

2. A helper package (file: "xacHelper_0.0.0.9000.tar.gz") that provides helper functions for certain analyses in the paper, unrelated the the general functionality of Monocle or our proposed method, BEAM.

Since this is an extensive analysis across several different single-cell RNA-seq datasets, we have broken our analysis down into multiple scripts so that users can selectively run parts of the analysis or paralellize the entire analysis as needed.

##Requirements
- Memory: We have tested all scripts with 100G of memory (RAM), but the memory usage in reality is much lower for most of the analysis
- Operating system: We have tested on our Linux-based computing environment, but MacOSX and other Unix-based operating systems should also work
- Tools: 
	1. You must have R version 3.1 or higher and permissions to install additional packages
	2. One library used during analysis requires rJava, which relies on Java being available (we have tested with Java SE Development Kit 7u17)
- Internet connection: Used to download packages for installation and a tarball with data required for the analysis (too large to provide with this distribution)

##Recreating the analysis
The entire analysis can be run in a series simply by running:
bash run_all_analyses.sh

Note that this script will first install all necessary R packages. Once the analysis has run the following folders and files will be populated:
- data: Data to recreate the analysis (too large to provide with this distribution) are downloaded to this folder
- main_figures: PDFs of panels from all figures in the main text that contain analysis (does not include diagrams) are saved here
- supplementary_figures: PDFs of panels from all supplementary figures that contain analysis (does not include diagrams) are saved here
- supplementary_data: XLS files of supplementary data submitted with the manuscript are saved here
- temp: Temporary files (mostly RData files) for use across different scripts are saved here

###Running portions of the analysis
While running the provided script (run_all_analyses.sh) is the simplest option, in some cases, it may be desirable to run parts of the analysis in isolation. While we have tried to minimize dependencies between analyses, some dependencies must be considered if choosing this route.

Before running any scripts, you must:
- Download the required data and unzip:
wget http://www.gs.washington.edu/~xqiu/proj/BEAM_analysis_data.tar.gz
tar -zxvf BEAM_analysis_data.tar.gz
rm BEAM_analysis_data.tar.gz

- Make the relevant destination directories
mkdir main_figures supplementary_figures supplementary_data temp

- Install all R packages (this may require user interaction to confirm permission to update existing packages):
Rscript install_packages.R


Once this is complete, individual analyses may be run. The following is a list of all scripts and their dependencies (all dependencies must be run before running a given script):

Lung Data Analysis
- prepare_lung_data.R (no dependencies)
- analysis_lung_data.R (depends on results from prepare_lung_data.R)
- gen_lung_figures.R (depends on results from analysis_lung_data.R, spikein_free_algorithm_benchmark.R, and benchmark_analysis.R)

Benchmarking Analysis
- spikein_free_algorithm_benchmark.R (depends on results from prepare_lung_data.R)
- benchmark_analysis.R (depends on results from prepare_lung_data.R)

HSMM Data
- analysis_HSMM_data.R (no dependencies)

Shalek Data
analysis_shalek_data.R (no dependencies)
gen_shalek_figures.R (depdends on results from analysis_shalek_data.R)

Other:
- analysis_UMI_data.R (no dependencies)
- analysis_other_supplementary_data.R (depends on prepare_lung_data.R)  
- gen_supplementary_figure.R (depends on results from analysis_lung_data.R, spikein_free_algorithm_benchmark>R, benchmark_analysis.R, analysis_other_supplementary_data.R, and analysis_HSMM_data.R)

##Additional Information 
We have also provided a Jupyter R notebook HTML file with all the figures generated in the manuscript included. This html file can be viewed for the main text and supplementary figures respectively at: 
- Main Text Figures: http://www.gs.washington.edu/~xqiu/proj/BEAM/BEAM_jupyter_notebook/main_figures_notebook.html 
- Supplementary Figures: http://www.gs.washington.edu/~xqiu/proj/BEAM/BEAM_jupyter_notebook/supplementary_figures_notebook.html
