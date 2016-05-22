load('./RData/prepare_lung_data.RData')
library(monocle)
library(xacHelper)

load_all_libraries()

############################make the landscape heatmap: 
mc_select <- coef(rlm(unlist(lapply(molModels_select, function(x) coef(x)[1])) ~ unlist(lapply(molModels_select, function(x) coef(x)[2]))))
mc <- coef(rlm(unlist(lapply(molModels, function(x) coef(x)[1])) ~ unlist(lapply(molModels, function(x) coef(x)[2]))))

optim_mc_func_fix_c_test_optim(c(as.numeric(mc_select[2]), as.numeric(mc_select[1])))
x_list <- split(expand.grid(c(seq(-6, -1, length.out = 2), -4.403166,  as.numeric(mc_select[2])), 
               c(seq(0, 4, length.out = 2), 2.77514, as.numeric(mc_select[1]))), 1:16)

# test the function: whether or not it will run fine
optim_mc_func_fix_c_test_optim(x_list[[1]])

# # run in parallel (mclapply cannot deal with situations when NAs are returned?)
optimization_landscape_3d <- mclapply(X = split(expand.grid(c(seq(-6, -1, length.out = 100), -4.403166,  as.numeric(mc_select[2])), 
               c(seq(0, 4, length.out = 100), 2.77514, as.numeric(mc_select[1]))), 1:102^2), optim_mc_func_fix_c_test_optim, mc.cores = detectCores())

# optimization_landscape_3d <- lapply(X = split(expand.grid(c(seq(-6, -1, length.out = 100), -4.403166,  as.numeric(mc_select[2])), 
#                    c(seq(0, 4, length.out = 100), 2.77514, as.numeric(mc_select[1]))), 1:102^2), optim_mc_func_fix_c_test_optim)


#update the mc optimization landscape: 
optim_mc_func_fix_c(mc_select[2], kb_intercept = mc_select[1], t_estimate = estimate_t(standard_cds),
          relative_expr_matrix = exprs(standard_cds), split_relative_expr_matrix = split_relative_exprs,
          alpha = df$mode_transcript, total_RNAs = pData(absolute_cds)$endo,
          cores = 1, weight_mode=0.17, weight_relative_expr=0.50, weight_total_rna=0.33, verbose = F,  ...)

# # provide perfect parameters: 
optimization_landscape_3d_perfect_parameters <- mclapply(X = split(expand.grid(c(seq(-6, -1, length.out = 100), -4.403166,  as.numeric(mc[2])), 
               c(seq(0, 4, length.out = 100), 2.77514, as.numeric(mc[1]))), 1:102^2), optim_mc_func_fix_c_test_optim, t_estimate = estimate_t(standard_cds),
          	   relative_expr_matrix = exprs(standard_cds), split_relative_expr_matrix = split_relative_exprs,
          	   alpha = df$mode_transcript, total_RNAs = pData(absolute_cds)$endo,
          	   cores = 1, weight_mode=0.17, weight_relative_expr=0.50, weight_total_rna=0.33, mc.cores = detectCores())

# # with default parameters from the paper: 
optimization_landscape_3d_perfect_parameters <- mclapply(X = split(expand.grid(c(seq(-6, -1, length.out = 100), -4.403166,  as.numeric(mc[2])), 
               c(seq(0, 4, length.out = 100), 2.77514, as.numeric(mc[1]))), 1:102^2), optim_mc_func_fix_c_test_optim, mc.cores = detectCores())

save.image('./RData/spikein_free_algorithm_sampling.RData')
