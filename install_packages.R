 #install all packages: 
 packages = c("ggplot2", "VGAM", "igraph", "pRlyr", "combinat", "fastICA", "irlba", "matrixStats", "reshape2", "R.utils", "snow", 
            "stringr", "modeest", "Hmisc", "boot", "doMC", "data.table", "fitdistrplus", "ggdendro", "gplots", "princurve", "sp", "devtools",
            "lmtest", "MASS", "mixsmsn", "pheatmap", "plyr", "pscl", "RColorBrewer", "VennDiagram", "zoo", "raster", "colorRamps", "grid")
 install.packages(packages, repo = 'http://cran.fhcrc.org/')

 bio_packages = c("Biobase", "BiocGenerics",  "limma", "edgeR", "DESeq", "DESeq2", "piano")
 source("http://bioconductor.org/biocLite.R")
 biocLite(bio_packages)

 # go to https://github.com/settings/tokens and generate personal tokens for install the private monocle / devtree package: 
 # install_github("cole-trapnell-lab/monocle-dev", auth_token = "2b5f9747e17c8512f1ecd2bf76f5df4730be21e2")
 # install_github("cole-trapnell-lab/branch-diff", auth_token = "2b5f9747e17c8512f1ecd2bf76f5df4730be21e2")

 install.packages('./xacHelper_0.0.0.9000.tar.gz', dependencies = TRUE)
 install.packages('./monocle_1.99.0.tar.gz', dependencies = TRUE)

 #install scde: 
 library(devtools)
 install_github('hms-dbmi/scde')

 	