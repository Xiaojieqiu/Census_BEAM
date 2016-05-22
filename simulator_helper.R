return_spike_true_rpc_matrix <- function(libraries, total_ladder_transcripts, ladder, total_mRNA){
	t(laply(libraries, function(lib){
			E <- lib$E
			M <- lib$M
			S <- lib$S
			L <- lib$L
			T <- c(M, L)
			E[names(M)]
			#spike-in regression: 

			# Fit a robust linear model of spike concentration vs. FPKM for each cell
			input.ERCC.annotation <- data.frame(numMolecules = ladder)
			ladder_expr <- data.frame(prop = ladder / (total_ladder_transcripts + total_mRNA) * 10^6)
			endo_expr <- data.frame(prop = E / (total_ladder_transcripts + total_mRNA)  * 10^6)
			save(file = 'debug_simulation', endo_expr,ladder_expr,  input.ERCC.annotation)
			molModels <- apply(ladder_expr, 2, function(cell_exprs, input.ERCC.annotation) {
			  
			  #print (cell_exprs)
			  spike_df <- input.ERCC.annotation 
			  spike_df <- cbind(spike_df, FPKM = cell_exprs)
			  spike_df$rounded_numMolecules <- round(spike_df$numMolecules)
			  # print (mean(spike_df$rounded_numMolecules))
			  
			  spike_df <- subset(spike_df, FPKM >= 1e-10)
			  spike_df$log_fpkm <- log10(spike_df$FPKM) 
			  spike_df$log_numMolecules <- log10(spike_df$numMolecules)
			  
			  # print (paste("geometric mean of log numMol ", mean(spike_df$log_numMolecules), "geometric mean of log FPKM ", mean(spike_df$log_fpkm)))
			         
			 molModel <- tryCatch({
			   #molModel <- vgam (rounded_numMolecules ~ sm.ns(log_fpkm, df=3), data=spike_df, family=negbinomial(zero=NULL))
			   molModel <- rlm(log_numMolecules ~ log_fpkm, data=spike_df)
			   
			   molModel
			 }, 
			 error = function(e) { print(e); NULL })
			 molModel
			}, input.ERCC.annotation)

			# Now use the per-cell linear models to produce a matrix of absolute transcript abundances 
			# for each gene in the genome, in each cell
			split_relative_exprs <- split(endo_expr, col(endo_expr, as.factor = T))
			norm_fpkms <- mapply(function(cell_exprs, molModel) {
  			  norm_df <- data.frame(log_fpkm=log10(unlist(cell_exprs)))
			  tryCatch({
			    # print(head(norm_df))
			    res <- 10^predict(molModel, type="response", newdata=norm_df)
			  }, 
			  error = function(e) {
			    rep(NA, length(cell_exprs))
			  })
			}, 
			split_relative_exprs, 
			molModels)

			norm_fpkms
		}))
}

return_capture_stats_df <- function(libraries, total_ladder_transcripts, ladder, capture_rate) {
	ldply(libraries, function(lib){
			E <- lib$E
			M <- lib$M
			S <- lib$S
			L <- lib$L
			T <- c(M, L)

			hypothetical_ladder <- total_ladder_transcripts * capture_rate * (ladder / sum(ladder))

			hypothetical_ladder_tpm <- 1e6 * hypothetical_ladder / ((total_ladder_transcripts + 1e6) * capture_rate)

			
			no_endo_regression <- rlm (log10(hypothetical_ladder) ~ log10(hypothetical_ladder_tpm))

			captured_ladder <- S[names(L)]
			captured_ladder <- captured_ladder[captured_ladder>0]

			#ladder_high_end <- which(captured_ladder == max(captured_ladder))
			phi = ladder / sum(ladder)
			z_i = sum(captured_ladder) / sum(S)

			#true_fixed_c = mean(log10(captured_ladder[ladder_high_end]))

			#print (head(M))
			tpm_true_no_ladder <- as.matrix(1e6*M / sum(M))
			#print (dim(tpm_true_no_ladder))
			mode_log_tpm = tryCatch({estimate_t(tpm_true_no_ladder)}, error=function(e) {print (e); NA })
			#mode_log_tpm=NA
			mode_molecules_in_cell = 10^dmode(log10(M[M > 0]))
			mode_molecules_in_capture = 10^dmode(log10(S[S > 0]))
			data.frame(mode_log_tpm=mode_log_tpm, 
					   mode_molecules_in_cell=mode_molecules_in_cell, 
					   mode_molecules_in_capture=mode_molecules_in_capture)
			})
}

return_spike_rpc_df <- function(libraries, total_ladder_transcripts, ladder, true_rpc_matrix, library_rpc_matrix) {
	res <- ldply(libraries, function(lib){
			E <- lib$E
			M <- lib$M
			S <- lib$S
			L <- lib$L
			T <- c(M, L)
			read_sample <- lib$read_sample
			tpm_est <- 1e6*read_sample / sum(read_sample)
			mode_log_tpm_estimate = tryCatch({estimate_t(as.matrix(tpm_est))}, error=function(e) {print (e); NA })
			# print(length(S))
			# print(length(total_ladder_transcripts))
			# print(length(ladder))

			fit_df <- data.frame(log_S = log10(S[names(L)]), 
								 log_ladder=log10(total_ladder_transcripts * ladder/sum(ladder)), 
								 log_tpm=log10(tpm_est[names(L)]))
			#print (fit_df)
			fit_df <- subset(fit_df, is.finite(log_S) & is.finite(log_tpm))
			tryCatch({

				fit <- rlm(log_ladder ~ log_tpm, data=subset(fit_df, log_tpm > -5))
				spike_rpc <- 10^predict(fit, newdata=data.frame(log_tpm=log10(tpm_est[names(M)])))
				#print (head(spike_rpc))
				true_nzgenes = names(E)[which(E > 0)]
				library_nzgenes = names(S[names(M)])[which(S[names(M)] > 0)]
				ladder_nzgenes = names(S[names(L)])[which(S[names(L)] > 0)]
				median_ladder_library_rpc = median(S[names(L)])
				median_estimated_tpm = median(tpm_est[tpm_est>0])

				ladder_low_end <- which(ladder == min(ladder))
				ladder_low_end_tpm = mean(tpm_est[names(L)[ladder_low_end]])

				ladder_high_end <- which(ladder == max(ladder))
				ladder_high_end_tpm = mean(tpm_est[names(L)[ladder_high_end]])
				#print (tpm_est[names(L)[ladder_low_end]])
			    spike_mse_true_rpc <- mean((spike_rpc[true_nzgenes] - E[true_nzgenes])^2, na.rm=T)
			    spike_mse_library_rpc <- mean((spike_rpc[library_nzgenes] - S[library_nzgenes])^2, na.rm=T)

			    num_genes_converged_true_rpc <- sum((abs(spike_rpc[true_nzgenes] - E[true_nzgenes]) / E[true_nzgenes]) < 0.2) / length(true_nzgenes)
			    num_genes_converged_library_rpc <- sum((abs(spike_rpc[library_nzgenes] - S[library_nzgenes]) / S[library_nzgenes]) < 0.2) / length(library_nzgenes)

			    molecules_for_mode_log_tpm <- 10^predict(fit, newdata=data.frame(log_tpm=log10(mode_log_tpm_estimate)))[1]


			 #    qp <- qplot(log_tpm, log_ladder, data=subset(fit_df, log_tpm > -5)) + 
			 #    	geom_abline(intercept=coefficients(fit)[1], slope=coefficients(fit)[2]) +
			 #    	scale_x_continuous(limits=c(1.5, 5)) +
			 #    	scale_y_continuous(limits=c(0.8, 3)) 
				# print(qp) 

				data.frame(mse_true_rpc = spike_mse_true_rpc,
						   mse_library_rpc = spike_mse_library_rpc,
						   num_genes_converged_true_rpc=num_genes_converged_true_rpc,
						   num_genes_converged_library_rpc=num_genes_converged_library_rpc,
						   cor_true_estimated=NA,
						   cor_library_estimated=NA,
						   num_true_nzgenes=length(true_nzgenes),
						   num_library_nzgenes=length(library_nzgenes),
						   num_ladder_nzgenes=length(ladder_nzgenes),
						   median_ladder_library_rpc=median_ladder_library_rpc,
						   median_estimated_tpm=median_estimated_tpm,
						   ladder_low_end_tpm=ladder_low_end_tpm,
						   ladder_high_end_tpm=ladder_high_end_tpm,
						   regression_k=coefficients(fit)[2],
						   regression_b=coefficients(fit)[1],
						   molecules_for_mode_log_tpm=molecules_for_mode_log_tpm,
						   mode_log_tpm_estimate=mode_log_tpm_estimate)
			}, error = function(e) { print (e); data.frame(mse_true_rpc = NA,
						   								   mse_library_rpc = NA,
						   								   num_genes_converged_true_rpc=NA,
						   								   num_genes_converged_library_rpc=NA,
						   								   cor_true_estimated=NA,
						   								   cor_library_estimated=NA,
						   								   num_true_nzgenes=length(true_nzgenes),
						   								   num_library_nzgenes=length(library_nzgenes),
						   								   num_ladder_nzgenes=length(ladder_nzgenes),
						   								   median_ladder_library_rpc=median_ladder_library_rpc,
						   								   median_estimated_tpm=median_estimated_tpm,
						   								   ladder_low_end_tpm=ladder_low_end_tpm,
						   								   ladder_high_end_tpm=ladder_high_end_tpm,
						   					  			   regression_k=NA,
						   					  			   regression_b=NA,
						   					  			   molecules_for_mode_log_tpm=NA,
						   					  			   mode_log_tpm_estimate=mode_log_tpm_estimate) })

			
		})

	spike_rpc_matrix <- return_spike_rpc_matrix(libraries, total_ladder_transcripts, ladder)
	spike_rpc_matrix <- return_spike_rpc_matrix(libraries, total_ladder_transcripts, ladder)

	sd_vec <- apply(true_rpc_matrix, 1, sd)
	num_genes_converged_true_rpc <- abs(spike_rpc_matrix - true_rpc_matrix)
	num_genes_converged_true_rpc <- num_genes_converged_true_rpc < 2 * sd_vec
	num_genes_converged_true_rpc <- colSums(num_genes_converged_true_rpc) / nrow(true_rpc_matrix > 0)

	cor_true_estimated <- unlist(lapply(1:10, function(x) cor(spike_rpc_matrix[, x], true_rpc_matrix[, x])))

	sd_vec <- apply(library_rpc_matrix, 1, sd)
	num_genes_converged_library_rpc <- abs(spike_rpc_matrix - library_rpc_matrix)
	num_genes_converged_library_rpc <- num_genes_converged_library_rpc < 2 * sd_vec
	num_genes_converged_library_rpc <- colSums(num_genes_converged_library_rpc) / nrow(library_rpc_matrix > 0)

	cor_library_estimated <- unlist(lapply(1:10, function(x) cor(library_rpc_matrix[, x], true_rpc_matrix[, x])))

	res$num_genes_converged_true_rpc <- num_genes_converged_true_rpc
	res$num_genes_converged_library_rpc <- num_genes_converged_library_rpc
	res$cor_true_estimated <- cor_true_estimated
	res$cor_library_estimated <- cor_library_estimated
	res
}

return_spike_rpc_matrix <- function(libraries, total_ladder_transcripts, ladder){
	t(laply(libraries, function(lib){
			E <- lib$E
			M <- lib$M
			S <- lib$S
			L <- lib$L
			T <- c(M, L)
			read_sample <- lib$read_sample
			tpm_est <- 1e6*read_sample / sum(read_sample)
			mode_log_tpm_estimate = tryCatch({estimate_t(as.matrix(tpm_est))}, error=function(e) {print (e); NA })
			fit_df <- data.frame(log_S = log10(S[names(L)]), 
								 log_ladder=log10(total_ladder_transcripts*ladder/sum(ladder)), 
								 log_tpm=log10(tpm_est[names(L)]))
			#print (fit_df)
			fit_df <- subset(fit_df, is.finite(log_S) & is.finite(log_tpm))
			spike_rpc <- tryCatch({

				fit <- rlm(log_ladder ~ log_tpm, data=subset(fit_df, log_tpm > -5))
				spike_rpc <- 10^predict(fit, newdata=data.frame(log_tpm=log10(tpm_est[names(M)])))
				spike_rpc
			}, error = function(e) NULL)
			spike_rpc
		}))
}

return_spike_tpm_mode <- function(libraries, total_ladder_transcripts, ladder){

	t(laply(libraries, function(lib){
				E <- lib$E
				M <- lib$M
				S <- lib$S
				L <- lib$L
				T <- c(M, L)
				read_sample <- lib$read_sample
				tpm_est <- 1e6*read_sample / sum(read_sample)
				mode_log_tpm_estimate = tryCatch({estimate_t(as.matrix(tpm_est))}, error=function(e) {print (e); NA })
				fit_df <- data.frame(log_S = log10(S[names(L)]), 
									 log_ladder=log10(total_ladder_transcripts*ladder/sum(ladder)), 
									 log_tpm=log10(tpm_est[names(L)]))
				#print (fit_df)
				fit_df <- subset(fit_df, is.finite(log_S) & is.finite(log_tpm))
				spike_rpc <- tryCatch({

					fit <- rlm(log_ladder ~ log_tpm, data=subset(fit_df, log_tpm > -5))
					spike_rpc <- 10^predict(fit, newdata=data.frame(log_tpm=log10(mode_log_tpm_estimate)))
					spike_rpc
				}, error = function(e) NULL)
				spike_rpc
			}))
}





