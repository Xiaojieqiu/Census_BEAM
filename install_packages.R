 #install all packages: 
 packages = c("ggplot2", "VGAM", "igraph", "pRlyr", "combinat", "fastICA", "irlba", "matrixStats", "reshape2", "R.utils", "snow", "tsne", "lle", ###DDRTree package 
            "stringr", "modeest", "Hmisc", "boot", "doMC", "data.table", "fitdistrplus", "ggdendro", "gplots", "princurve", "sp", "devtools",
            "lmtest", "MASS", "mixsmsn", "pheatmap", "plyr", "pscl", "RColorBrewer", "VennDiagram", "zoo", "raster", "colorRamps", "grid")
 install.packages(packages, repo = 'http://cran.fhcrc.org/')

 bio_packages = c("Biobase", "BiocGenerics", "cummeRbund", "limma", "edgeR", "DESeq", "DESeq2", "piano")
 source("http://bioconductor.org/biocLite.R")
 biocLite(bio_packages)


 #note that the latest HSMMSingleCell version 0.104.0 is only available in R 3.2. You need to download it if you use it in R 3.1. 
 install.packages('./HSMMSingleCell_0.104.0.tar.gz', dependencies = TRUE)
 install.packages('./monocle_1.99.0.tar.gz', dependencies = TRUE) 
 install.packages('./xacHelper_0.0.0.9000.tar.gz', dependencies = TRUE)

 #install scde: 
 library(devtools)
 install_github('hms-dbmi/scde')

 	
