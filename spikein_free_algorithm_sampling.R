load('./RData/prepare_lung_data.RData')
# library(devtools)
# load_all('~/Projects/monocle-dev')
library(monocle)
library(MASS)
library(sp)
library(plyr)
library(raster)
library(xacHelper)
library(grid)
library(RColorBrewer)
load_all_libraries()

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


split_relative_exprs <- split(exprs(standard_cds), rep(1:ncol(exprs(standard_cds)), each = nrow(exprs(standard_cds))))

test <- mapply(function(cell_dmode, model) {
  predict(model, newdata = data.frame(log_fpkm = log10(cell_dmode)), type = 'response')
}, as.list(estimate_t(exprs(standard_cds)[1:transcript_num, ])), molModels) #molModels_select

df <- pData(absolute_cds)
df$mode_transcript <- 10^test
df$estimate_mode <- estimate_t(exprs(standard_cds))

#update the mc optimization landscape: 
optim_mc_func_fix_c(mc_select[2], kb_intercept = mc_select[1], t_estimate = estimate_t(standard_cds),
          relative_expr_matrix = exprs(standard_cds), split_relative_expr_matrix = split_relative_exprs,
          alpha = df$mode_transcript, total_RNAs = pData(absolute_cds)$endo,
          cores = 1, weight_mode=0.17, weight_relative_expr=0.50, weight_total_rna=0.33, verbose = F)

# # provide perfect parameters: 
optim_mc_func_fix_c_simulation <- function (mc_list, t_estimate = estimate_t(TPM_isoform_count_cds),
          relative_expr_matrix = relative_expr_matrix, split_relative_expr_matrix = split_relative_exprs,
          alpha = rep(1, ncol(relative_expr_matrix)), total_RNAs = rep(100000, ncol(relative_expr_matrix)),
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
    return(list(m = as.numeric(kb_slope_val), c = as.numeric(kb_intercept_val), optim_res = res))
  #
  # if(is.finite(res))
  #   return(res)
  # else
  #   return(1e10 * runif(1)) #Census should not run this part since non-finite values are removed
}

optimization_landscape_3d_perfect_parameters <- mclapply(X = split(expand.grid(c(seq(-6, -1, length.out = 2), -4.403166,  as.numeric(mc[2])), 
               c(seq(0, 4, length.out = 2), 2.77514, as.numeric(mc[1]))), 1:4^2), optim_mc_func_fix_c_simulation, t_estimate = estimate_t(standard_cds),
          	   relative_expr_matrix = exprs(standard_cds), split_relative_expr_matrix = split_relative_exprs,
          	   alpha = df$mode_transcript, total_RNAs = pData(absolute_cds)$endo,
          	   cores = 1, weight_mode=0.17, weight_relative_expr=0.50, weight_total_rna=0.33, mc.cores = detectCores())


optimization_landscape_3d <- mclapply(X = split(expand.grid(c(seq(-6, -1, length.out = 100), -4.403166,  as.numeric(mc[2])), 
               c(seq(0, 4, length.out = 100), 2.77514, as.numeric(mc[1]))), 1:102^2), optim_mc_func_fix_c_simulation, t_estimate = estimate_t(standard_cds),
          	   relative_expr_matrix = exprs(standard_cds), split_relative_expr_matrix = split_relative_exprs,
          	   alpha = df$mode_transcript, total_RNAs = pData(absolute_cds)$endo,
          	   cores = 1, weight_mode=0.17, weight_relative_expr=0.50, weight_total_rna=0.33, mc.cores = detectCores())


############################make the landscape heatmap: 
optimization_landscape_3d_trim <- lapply(optimization_landscape_3d, function(x) x[c('m', 'c', 'optim_res')])
optimization_matrix<- do.call(rbind.data.frame, optimization_landscape_3d_trim)

optimization_matrix_filt <- subset(optimization_matrix, is.nan(optim_res) == FALSE & is.finite(optim_res))
max_optim_score <- 3
optimization_matrix_filt$optim_res[optimization_matrix_filt$optim_res > max_optim_score] <- max_optim_score

spdf <- SpatialPointsDataFrame( data.frame( x = as.numeric(as.character(optimization_matrix_filt$m)), y = as.numeric(as.character(optimization_matrix_filt$c)) ) , data = data.frame( z = as.numeric(as.character(optimization_matrix_filt$optim_res)) ) )

# Make an evenly spaced raster, the same extent as original data
e <- extent( spdf )

# Determine ratio between x and y dimensions
ratio <- ( e@xmax - e@xmin ) / ( e@ymax - e@ymin )

# Create template raster to sample to
r <- raster( nrows = 10 , ncols = floor( 9 * ratio ) , ext = extent(spdf) ) 
rf <- rasterize( spdf , r , field = "z" , fun = mean )

# We can then plot this using `geom_tile()` or `geom_raster()`
rdf <- data.frame( rasterToPoints( rf ) )

optimal_solution <- head(arrange(optimization_matrix_filt, optim_res), 1)
pdf('./main_figures/fig3e.pdf', width = 1.38, height = 1.25)
ggplot( NULL ) + geom_raster( data = rdf , aes( x , y , fill = log10(layer) ) ) + 
    annotate("text", x = -4.1, y = 2.7, label = "True (m,c)", color="magenta", size=2) + 
    annotate("point", x = mc_select[2], y = mc_select[1], color="magenta", size = 1) + #-4.636364 3.434343
    #annotate("text", x = -3.7, y = 3.2, label = "True (m,c)") + 
    annotate("text", x = -3.85, y = 3.2, label = "Algorithm (m,c)", color="red", size=2) + 
    annotate("point", x = as.numeric(as.character(optimal_solution$m)), y = as.numeric(as.character(optimal_solution$c)), color="red", size=1) + 
    scale_fill_gradientn(guide=guide_legend(title=expression(paste(log[10](F)))), colours=brewer.pal(name="YlGnBu", n=7)) +
    xlab("m") + ylab("c") +
    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + scale_size(range = c(0.1, 2)) + 
    theme(panel.background = element_rect(fill='white')) + nm_theme()
dev.off()

#create the helper pdf file to annotate the figure: 
pdf('./tmp/fig3e_helper.pdf', width = 5, height = 1.5)
ggplot( NULL ) + geom_raster( data = rdf , aes( x , y , fill = log10(layer) ) ) + 
    annotate("text", x = -3.85, y = 3.2, label = "True (m,c)", color="magenta", size=2) + 
    annotate("point", x = mc_select[2], y = mc_select[1], color="magenta", size = 1) + #-4.636364 3.434343
    #annotate("text", x = -3.7, y = 3.2, label = "True (m,c)") + 
    annotate("text", x = -4.1, y = 2.7, label = "Algorithm (m,c)", color="red", size=2) + 
    annotate("point", x = as.numeric(as.character(optimal_solution$m)), y = as.numeric(as.character(optimal_solution$c)), color="red", size=1) + 
    scale_fill_gradientn(guide=guide_legend(title=expression(paste(log[10](F)))), colours=brewer.pal(name="YlGnBu", n=7)) +
    xlab("m") + ylab("c") +
    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + scale_size(range = c(0.1, 2)) + 
    theme(panel.background = element_rect(fill='white')) + nm_theme()
dev.off()

save.image('./RData/spikein_free_algorithm_sampling.RData')
