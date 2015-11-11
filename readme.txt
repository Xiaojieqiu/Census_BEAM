#README

Our analysis for the manuscript is mainly centered in the BEAM analysis offered from the monocle package (file: "monocle_1.99.0.tar.gz"). The xacHelper package (file: "xacHelper_0.0.0.9000.tar.gz") is used to facilitate the analysis. 

Since this is a huge analysis on four datasets, in order to reproduce all the figures exactly in our manuscript, we break this analysis down into several parts so that users can parallelize the reproduction. 

To reproduce the results in the manuscript, a cluster with at 100 G memory is required. In the command line, we need to first download data required for this analysis. Then we need to install the monocle and xacHelper packages as well as other dependent packages. After all the data are downloaded and the R packages are installed, we can perform the analysis for lung, HSMM, UMI and Shalek datasets individually or in parallel (for example, using the sge system). The prepare_other_supplementary_data.R script depends on the objects generated from prepare_lung_data.R script. For each script, there will be a RData file with all objects generated in the script saved with the names as the some as the script. 

In terms of time, the following are some rough estimation based on a cluster with 64 cores and 100 G memory used: 
Lung dataset analysis: (~ 2 day)
Prepare data for other supplementary file: (~ 6h)
HSMM dataset analysis: (~ 3 day)
Shalek datast analysis: (~ 1h)
UMI dataset analysis: (~ 2h)

For generating figures, each script will need less than 2 hrs. 

The instructions on how to reproduce this analysis is included in the "drive_script.sh".
 
Note that we also provided a Jupyter R notebook html file with all the figures generated in the manuscript included. The packages version for the session we run this analysis is included in the notebook too. Specially, for using the venneuler package which depends on rJava, Java SE Development Kit 7u17 is used in our case. 