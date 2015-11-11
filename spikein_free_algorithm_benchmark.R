############################make the landscape heatmap: 
mc_select <- coef(rlm(unlist(lapply(molModels_select, function(x) coef(x)[1])) ~ unlist(lapply(molModels_select, function(x) coef(x)[2]))))

optim_mc_func_fix_c_test_optim(c(as.numeric(mc_select[2]), as.numeric(mc_select[1])))
x_list <- split(expand.grid(c(seq(-6, -1, length.out = 2), -4.403166,  as.numeric(mc_select[2])), 
                   c(seq(0, 4, length.out = 2), 2.77514, as.numeric(mc_select[1]))), 1:16)

# test the function: whether or not it will run fine
optim_mc_func_fix_c_test_optim(x_list[[1]])
# mclapply cannot deal with situations when NAs are returned 
# optimization_landscape_3d <- mclapply(X = split(expand.grid(c(seq(-6, -1, length.out = 100), -4.403166,  as.numeric(mc_select[2])), 
#                    c(seq(0, 4, length.out = 100), 2.77514, as.numeric(mc_select[1]))), 1:102^2), optim_mc_func_fix_c_test_optim, mc.cores = detectCores())

optimization_landscape_3d <- lapply(X = split(expand.grid(c(seq(-6, -1, length.out = 100), -4.403166,  as.numeric(mc_select[2])), 
                   c(seq(0, 4, length.out = 100), 2.77514, as.numeric(mc_select[1]))), 1:102^2), optim_mc_func_fix_c_test_optim)

#this doesn't work, need to use roster and sp: 
optimization_landscape_2d <- melt(optimization_landscape_3d)
qplot(Var1, Var2, fill=optim_res, geom="tile", data=optimization_landscape_2d) + scale_fill_gradientn(colours=rainbow(7)) #hmcols
ggsave(filename = paste(elife_directory, 'eLife_fig4E.pdf', sep = ''), width = 1.38, height = 1.25)

save.image()
