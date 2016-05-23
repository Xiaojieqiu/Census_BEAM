load('./RData/prepare_lung_data.RData')
library(devtools)
load_all('~/Projects/monocle-dev')

# load_all_libraries()

############################make the landscape heatmap: 
mc_select <- coef(rlm(unlist(lapply(molModels_select, function(x) coef(x)[1])) ~ unlist(lapply(molModels_select, function(x) coef(x)[2]))))
mc <- coef(rlm(unlist(lapply(molModels, function(x) coef(x)[1])) ~ unlist(lapply(molModels, function(x) coef(x)[2]))))

# optim_mc_func_fix_c((mc_select[2]), as.numeric(mc_select[1]))
x_list <- split(expand.grid(c(seq(-6, -1, length.out = 2), -4.403166,  as.numeric(mc_select[2])), 
               c(seq(0, 4, length.out = 2), 2.77514, as.numeric(mc_select[1]))), 1:16)

# test the function: whether or not it will run fine
# optim_mc_func_fix_c(x_list[[1]])

# # run in parallel (mclapply cannot deal with situations when NAs are returned?)
# optimization_landscape_3d <- mclapply(X = split(expand.grid(c(seq(-6, -1, length.out = 100), -4.403166,  as.numeric(mc_select[2])), 
#                c(seq(0, 4, length.out = 100), 2.77514, as.numeric(mc_select[1]))), 1:102^2), optim_mc_func_fix_c, mc.cores = detectCores())

# optimization_landscape_3d <- lapply(X = split(expand.grid(c(seq(-6, -1, length.out = 100), -4.403166,  as.numeric(mc_select[2])), 
#                    c(seq(0, 4, length.out = 100), 2.77514, as.numeric(mc_select[1]))), 1:102^2), optim_mc_func_fix_c)


# split_relative_exprs <- split(exprs(standard_cds), rep(1:ncol(exprs(standard_cds)), each = nrow(exprs(standard_cds))))

# test <- mapply(function(cell_dmode, model) {
#   predict(model, newdata = data.frame(log_fpkm = log10(cell_dmode)), type = 'response')
# }, as.list(estimate_t(exprs(standard_cds)[1:transcript_num, ])), molModels) #molModels_select

# df <- pData(absolute_cds)
# df$mode_transcript <- 10^test
# df$estimate_mode <- estimate_t(exprs(standard_cds))

#update the mc optimization landscape: 
optim_mc_func_fix_c(mc_select[2], kb_intercept = mc_select[1], t_estimate = estimate_t(standard_cds),
          relative_expr_matrix = exprs(standard_cds), split_relative_expr_matrix = split_relative_exprs,
          alpha = df$mode_transcript, total_RNAs = pData(absolute_cds)$endo,
          cores = 1, weight_mode=0.17, weight_relative_expr=0.50, weight_total_rna=0.33, verbose = F)

# # provide perfect parameters: 
optim_mc_func_fix_c_simulation <- function (mc_list, t_estimate = estimate_t(TPM_isoform_count_cds),
          relative_expr_matrix = relative_expr_matrix, split_relative_expr_matrix = split_relative_exprs,
          alpha = rep(1, ncol(relative_expr_matrix)), total_RNAs = rep(150000, ncol(relative_expr_matrix)),
          cores = 1, weight_mode=0.17, weight_relative_expr=0.50, weight_total_rna=0.33, verbose = F,  ...) {
  data('spike_df') #add the spikein dataset

  kb_slope_val <- mc_list[1]
  kb_intercept_val <- mc_list[2]

  cell_num <- ncol(relative_expr_matrix)
  names(t_estimate) <- colnames(relative_expr_matrix)
  split_t <- split(Matrix::t(t_estimate), col(as.matrix(Matrix::t(t_estimate)), as.factor = T))

  total_rna_df <- data.frame(Cell = colnames(relative_expr_matrix), t_estimate = t_estimate, alpha_v = alpha)

  t_k_b_solution <- tryCatch({
    k_b_solution <- plyr::ddply(total_rna_df, .(Cell), function(x) {
      a_matrix <- matrix(as.numeric(c(log10(x[, "t_estimate"]), 1,
                           kb_slope_val, -1)), ncol = 2, nrow = 2, byrow = T)

      b_matrix <- matrix(c(log10(x[, "alpha_v"]), -kb_intercept_val), nrow = 2, byrow = T)
      k_b_solution <- Matrix::t(solve(a_matrix, b_matrix))
    })
    k_b_solution},
    error = function(e) {print(e); c(NA, NA)}
  )

  if(any(is.na(t_k_b_solution)))
    return(NA)

  cell_dmode <- tryCatch({
    if(cores > 1){
      cell_dmode <- unlist(mcmapply(opt_norm_t, split_t, split_relative_expr_matrix, alpha, kb_slope = kb_slope_val, kb_intercept = kb_intercept_val, pseudocnt = 0.01, mc.cores = cores))
      adj_est_std_cds <- unlist(mcmapply(opt_norm_t, split_t, split_relative_expr_matrix, alpha, kb_slope = kb_slope_val, kb_intercept = kb_intercept_val, pseudocnt = 0.01, return_norm = T, mc.cores = cores))
    }
    else {
      cell_dmode <- unlist(mapply(opt_norm_t, split_t, split_relative_expr_matrix, alpha, kb_slope = kb_slope_val, kb_intercept = kb_intercept_val, pseudocnt = 0.01))
      adj_est_std_cds <- unlist(mapply(opt_norm_t, split_t, split_relative_expr_matrix,  alpha, kb_slope = kb_slope_val, kb_intercept = kb_intercept_val, pseudocnt = 0.01, return_norm = T))
    }
    cell_dmode},
    error = function(e) {print(e); NA}
  )

  if(any(is.na(cell_dmode)))
    return(NA)

  sum_total_cells_rna <- colSums(adj_est_std_cds)

  #Jenn-Shannon divergence:
  p_df <- makeprobs(relative_expr_matrix) #relative expression
  p_list <- split(Matrix::t(p_df), 1:ncol(p_df))
  q_df <- makeprobs(adj_est_std_cds) #no rounding
  q_list <- split(Matrix::t(q_df), 1:ncol(q_df))

  dist_divergence <- mcmapply(function(x, y) {
    JSdistVec(x, y)
  }, p_list, q_list, mc.cores = cores)


  gm_dist_divergence <- exp(mean(log(dist_divergence)))

  #total RNA MSE
  sum_total_cells_rna_finite <- sum_total_cells_rna[is.finite(sum_total_cells_rna)]
  total_RNAs_finite <- total_RNAs[is.finite(sum_total_cells_rna)]
  total_rna_obj <- exp(mean(log(((sum_total_cells_rna_finite -  total_RNAs_finite)/total_RNAs_finite)^2))) #use geometric mean to avoid outlier cells

  #mode MSE
  mode_obj <- exp(mean(log(((cell_dmode[is.finite(cell_dmode)] - alpha[is.finite(cell_dmode)])/alpha[is.finite(cell_dmode)])^2)))

  relative_expr_obj <- gm_dist_divergence * 10 #give more weight on this

  #objective
  res <- weight_mode * mode_obj + weight_relative_expr * relative_expr_obj + weight_total_rna * total_rna_obj

  if(verbose){
    message('current m, c values are ', paste(kb_slope_val, kb_intercept_val, sep = ', '))
    message('objective is:')
    print(res)

    message('\n\ntotal_rna_obj is ', total_rna_obj)
    message('mode_obj is ', mode_obj)
    message('relative_expr_obj is ', relative_expr_obj)

    message('\n\nmean modes are:')
    print (mean(cell_dmode))
    message('mean target modes are:')
    print (mean(alpha))
    message('mean mode delta is:')
    print (mean(cell_dmode - alpha))
    message('mean total RNA delta is:')
    print (mean(sum_total_cells_rna_finite -  total_RNAs_finite))

  }
  #   return(list(m = m_val, c = c_val, dmode_rmse_weight_total = dmode_rmse_weight_total, gm_dist_divergence = gm_dist_divergence, dist_divergence_round = dist_divergence_round,
  #               cell_dmode = cell_dmode, t_k_b_solution = t_k_b_solution, sum_total_cells_rna = sum_total_cells_rna, optim_res = res))
  #
  if(is.finite(res))
    return(res)
  else
    return(1e10 * runif(1)) #Census should not run this part since non-finite values are removed
}

optimization_landscape_3d_perfect_parameters <- mclapply(X = split(expand.grid(c(seq(-6, -1, length.out = 2), -4.403166,  as.numeric(mc[2])), 
               c(seq(0, 4, length.out = 2), 2.77514, as.numeric(mc[1]))), 1:4^2), optim_mc_func_fix_c_simulation, t_estimate = estimate_t(standard_cds),
          	   relative_expr_matrix = exprs(standard_cds), split_relative_expr_matrix = split_relative_exprs,
          	   alpha = df$mode_transcript, total_RNAs = pData(absolute_cds)$endo,
          	   cores = 1, weight_mode=0.17, weight_relative_expr=0.50, weight_total_rna=0.33, mc.cores = detectCores())


optimization_landscape_3d_perfect_parameters <- mclapply(X = split(expand.grid(c(seq(-6, -1, length.out = 100), -4.403166,  as.numeric(mc[2])), 
               c(seq(0, 4, length.out = 100), 2.77514, as.numeric(mc[1]))), 1:102^2), optim_mc_func_fix_c_simulation, t_estimate = estimate_t(standard_cds),
          	   relative_expr_matrix = exprs(standard_cds), split_relative_expr_matrix = split_relative_exprs,
          	   alpha = df$mode_transcript, total_RNAs = pData(absolute_cds)$endo,
          	   cores = 1, weight_mode=0.17, weight_relative_expr=0.50, weight_total_rna=0.33, mc.cores = detectCores())

save.image('./RData/spikein_free_algorithm_sampling.RData')
