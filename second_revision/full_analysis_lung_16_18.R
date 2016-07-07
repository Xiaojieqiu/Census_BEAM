conditions <- c('E16.5', 'E18.5')
script.dir <- dirname(sys.frame(1)$ofile)
source(paste(script.dir, '/deg_benchmark_analysis.R', sep = ''))
source(paste(script.dir, '/first_revision/cmpr_three_packages.R', sep = ''))
source(paste(script.dir, '/second_revision/deg_benchmark_lung_roc_auc.R', sep = ''))

