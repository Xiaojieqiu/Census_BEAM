source('./simulator_helper.R')

simulate_sequencing_with_ladder <- function(original_fpkm_dist,
										#genes_expressed=seq(1000, 25000, by=500), 
										genes_expressed=10000,
										total_mRNAs=150000, #300000 
										total_ladder_transcripts=50000,
										ladder=rep(c(2^(seq(from=0,to=24, by=2))), 4),
										ladder_degradation_rates=0,#c(0.01, 0.05, 0.15, 0.25, seq(0.1, 1, by=0.1)),
										mRNA_degradation_rates=0,#c(0.01, 0.05, 0.15, 0.25, seq(0.1, 1, by=0.1)),
										capture_rates=0.1, #0.25
										reads_per_cell=1000000,
										cells=15,
										run_rel2abs=TRUE,
										return_matrices=FALSE,
										calibration_trials=10,
										... # Additional args passed to relative2abs
										#mRNA_range=seq(10000, 500000, by=10000),
										#seq_range=seq(10000, 2000000, by=10000)
										){

	true_rpc_matrix <- NA
	library_rpc_matrix <- NA 
	tpm_matrix <- NA

	full_res <- ldply(genes_expressed, function(num_genes){
	original_fpkm_dist <- original_fpkm_dist[original_fpkm_dist > 0]
	expr_profile <- original_fpkm_dist[sample(length(original_fpkm_dist), num_genes)]
	expr_profile <- expr_profile / sum(expr_profile)
	ldply(total_mRNAs, function(total_mRNA) {
	ldply(capture_rates, function(capture_rate){
	ldply(reads_per_cell, function(reads) {
	ldply(mRNA_degradation_rates, function(mRNA_degradation_rate){
	ldply(ladder_degradation_rates, function(ladder_degradation_rate){
		libraries <- llply(1:cells, function(i) {

			#the library 
			E <- rmultinom(1, total_mRNA, expr_profile)[,1]
			names(E) <- paste("endo", 1:length(E), sep="_")
			#trial_modes <- colHistMode(T)
			tpm_endo = 1e6*E / sum(E)

			#mRNA_degradation_rate = 0
			
			M <- rmultinom(1, (1-mRNA_degradation_rate)*total_mRNA, E)[,1]
			names(M) <- paste("endo", 1:length(M), sep="_")

			#ladder_degradation_rate = 0
			# print('hello')
			# print(ladder/sum(ladder))
			L <- rmultinom(1, (1-ladder_degradation_rate)*total_ladder_transcripts, ladder/sum(ladder))[,1]
			names(L) <- paste("ladder", 1:length(L), sep="_")
			T <- c(M, L)

			tpm_true <- 1e6*T / sum(T)

			S <- rmultinom(1, sum(T)*capture_rate, tpm_true)[,1]
			names(S) <- names(T)
			tpm_captured = 1e6*S / sum(S)

			read_sample <- rmultinom(1, reads, tpm_captured)[,1]
			names(read_sample) <- names(T)

			tpm_est = 1e6*read_sample / sum(read_sample)
			names(tpm_est) <- names(T)

			list(E=E, M=M, L=L, S=S, read_sample=read_sample)
		})


		#use spike-in to calculate the true RPC: 
		true_rpc_matrix <<- t(laply(libraries, function(lib){
			E <- lib$E
			M <- lib$M
			S <- lib$S
			L <- lib$L
			T <- c(M, L)
			E[names(M)]
			}))


		library_rpc_matrix <<- t(laply(libraries, function(lib){
			E <- lib$E
			M <- lib$M
			S <- lib$S
			L <- lib$L
			T <- c(M, L)
			S[names(M)]
			}))

		tpm_matrix <<- t(laply(libraries, function(lib){
			E <- lib$E
			M <- lib$M
			S <- lib$S
			L <- lib$L
			T <- c(M, L)
			read_sample <- lib$read_sample[names(M)]
			tpm_est <- 1e6*read_sample / sum(read_sample)

			}))
		row.names(tpm_matrix) <- c(names(libraries[1]$M), names(libraries[1]$L))
		colnames(tpm_matrix) <- names(libraries)

		spike_true_rpc_matrix <<- return_spike_true_rpc_matrix(libraries, total_ladder_transcripts, ladder, total_mRNA)

	 	split_relative_exprs <- split(tpm_matrix, col(tpm_matrix, as.factor = T))
		total_mRNA = rep(total_mRNA, length(split_relative_exprs)) 
       	capture_rate = rep(capture_rate, length(split_relative_exprs))
       	reads = rep(reads_per_cell, length(split_relative_exprs))
		save(file = 'debug_simulation', split_relative_exprs, 
               ladder, total_ladder_transcripts, total_mRNA, capture_rate, reads)

		calibrated_mc_modes <- lapply(1:length(split_relative_exprs), 
			   monocle:::calibrate_mode, 
			   tpm_distribution = split_relative_exprs, 
               ladder = ladder, 
               total_ladder_transcripts = total_ladder_transcripts,
               total_mRNA = total_mRNA, 
               capture_rate = capture_rate,
               reads = reads)

		calibrated_mc_modes_df <- do.call(rbind.data.frame, calibrated_mc_modes)
		mc_model <- rlm(b ~ k, data = calibrated_mc_modes_df)
		calibrated_mc <- list(m = coef(mc_model)[2], c = coef(mc_model)[1])
		calibrated_modes_cell <- calibrated_mc_modes_df[, 1]

		capture_stats_df <- return_capture_stats_df(libraries, total_ladder_transcripts, ladder, capture_rate)
		spike_rpc_df <- return_spike_rpc_df(libraries, total_ladder_transcripts, ladder, true_rpc_matrix, library_rpc_matrix)
		spike_rpc_matrix <- return_spike_rpc_matrix(libraries, total_ladder_transcripts, ladder)

		# return(spike_rpc_df)


		kb_lm <- rlm(regression_b ~ regression_k, data=spike_rpc_df)
		spike_m <- coefficients(kb_lm)[2]
		spike_c <- coefficients(kb_lm)[1]

		spike_tpm_mode <<- return_spike_tpm_mode(libraries, total_ladder_transcripts, ladder)
		# print('spike_tpm_mode: ')
		# print (spike_tpm_mode)

		#print (dim(tpm_matrix))
		colnames(tpm_matrix) <- paste("Cell_", seq(1,ncol(tpm_matrix)), sep="")
		TPM_cds <- newCellDataSet(tpm_matrix)

		TPM_thresh <- 0
		
		spike_res <- data.frame(total_mRNA=total_mRNA,
								num_genes=num_genes,
								capture_rate=capture_rate,
								reads=reads,
								mRNA_degradation_rate=mRNA_degradation_rate,
								ladder_degradation_rate=ladder_degradation_rate,
								method="spike", 
								mode_log_tpm = capture_stats_df$mode_log_tpm,
								mode_molecules_in_cell = capture_stats_df$mode_molecules_in_cell,
								mode_molecules_in_capture = capture_stats_df$mode_molecules_in_capture,
								cor_true_estimated=spike_rpc_df$cor_true_estimated,
						   		cor_library_estimated=spike_rpc_df$cor_library_estimated,
								mse_true_rpc=spike_rpc_df$mse_true_rpc,
								mse_library_rpc=spike_rpc_df$mse_library_rpc,
								num_genes_converged_true_rpc=spike_rpc_df$num_genes_converged_true_rpc,
						   		num_genes_converged_library_rpc=spike_rpc_df$num_genes_converged_library_rpc,
						   		#num_true_nzgenes=spike_rpc_df$num_true_nzgenes,
						   	    #num_library_nzgenes=spike_rpc_df$num_library_nzgenes,
						   	    #num_ladder_nzgenes=spike_rpc_df$num_ladder_nzgenes,
						   	    #median_ladder_library_rpc=spike_rpc_df$median_ladder_library_rpc,
						   	    #median_estimated_tpm=spike_rpc_df$median_estimated_tpm,
						   	    ladder_low_end_tpm=spike_rpc_df$ladder_low_end_tpm,
						   	    ladder_high_end_tpm=spike_rpc_df$ladder_high_end_tpm,
								regression_k=spike_rpc_df$regression_k,
								regression_b=spike_rpc_df$regression_b,
								spike_m=spike_m,
								spike_c=spike_c,
								calibrated_m=calibrated_mc[["m"]],
								calibrated_c=calibrated_mc[["c"]],
								calibrated_mode = mean(calibrated_modes_cell),
								mode_log_tpm_estimate=spike_rpc_df$mode_log_tpm_estimate,
								molecules_for_mode_log_tpm=spike_rpc_df$molecules_for_mode_log_tpm,
								t_estimate=estimate_t(exprs(TPM_cds), relative_expr_thresh=TPM_thresh))

		if (run_rel2abs){
			message("running relative2abs())")
			print ("****")
			# print(qplot(tpm_matrix[,1], log="x"))
			# print (calibrated_modes)
			estimated_rpc_matrix <- as.matrix(relative2abs(TPM_cds,
											 #total_RNAs=1000000, 
											 estimate_t(exprs(TPM_cds), relative_expr_thresh=0.1),
											 #alpha_v=spike_rpc_df$molecules_for_mode_log_tpm,
											 #expected_mRNA_mode=calibrated_modes_capture,
											 verbose=F,
											 #kb_slope=calibrated_mc[["m"]],
											 #kb_intercept=calibrated_mc[["c"]],
											 #kb_slope_rng=sort(c(calibrated_mc[["m"]]*0.8, calibrated_mc[["m"]]*1.2)),
											 #c_rng=sort(c(calibrated_mc[["c"]]*0.8, calibrated_mc[["c"]]*1.2)),
											 expected_total_mRNAs = total_mRNA  / capture_rate,
											 expected_capture_rate = capture_rate,
											 reads_per_cell=reads,
											 use_fixed_intercept=T,
											 ...))
			print ("****")

		    print(qplot(estimated_rpc_matrix[,1], true_rpc_matrix[,1], log="xy") + geom_abline(color="red") + 
		    		ggtitle(paste('capture_rate:', capture_rate)))


		    # exclude genes that aren't expressed at all - in this simulation,
		    # the methods won't detect expression for these genes, which means
		    # we will overestimate our accuracy if we leave these genes in:
		    true_nzgenes <- rowSums(true_rpc_matrix) > 0
			library_nzgenes <- rowSums(library_rpc_matrix) > 0

			spike_free_mse_true_rpc <- colMeans((estimated_rpc_matrix[true_nzgenes,] - true_rpc_matrix[true_nzgenes,])^2)
			spike_free_mse_library_rpc <- colMeans((estimated_rpc_matrix[library_nzgenes,] - library_rpc_matrix[library_nzgenes,])^2)
			#print(dim(estimated_rpc_matrix - true_rpc_matrix))
			#print (spike_free_mse_true_rpc)
			#print (spike_free_mse_library_rpc)

			#num_genes_converged_true_rpc <- colSums((abs(estimated_rpc_matrix[true_nzgenes,] - true_rpc_matrix[true_nzgenes,]) / true_rpc_matrix[true_nzgenes,]) < 0.2) / length(sum(true_nzgenes))
			#num_genes_converged_library_rpc <- colSums((abs(estimated_rpc_matrix[library_nzgenes,] - library_rpc_matrix[library_nzgenes,]) / library_rpc_matrix[library_nzgenes,]) < 0.2) / length(sum(library_nzgenes))
			
			#print (class(estimated_rpc_matrix))
			#print (class(true_rpc_matrix))
			save(file = 'test_simulation', estimated_rpc_matrix, spike_tpm_mode, library_rpc_matrix, spike_m, spike_c, true_rpc_matrix, TPM_cds, TPM_thresh, total_mRNA, capture_rate, reads)

			sd_vec <- apply(true_rpc_matrix, 1, sd)
			num_genes_converged_true_rpc <- abs(estimated_rpc_matrix - true_rpc_matrix)
			num_genes_converged_true_rpc <- num_genes_converged_true_rpc < 2 * sd_vec
			num_genes_converged_true_rpc <- colSums(num_genes_converged_true_rpc) / nrow(true_rpc_matrix > 0)

			cor_true_estimated <- mean(unlist(lapply(1:10, function(x) cor(estimated_rpc_matrix[, x], true_rpc_matrix[, x]))))

			sd_vec <- apply(library_rpc_matrix, 1, sd)
			num_genes_converged_library_rpc <- abs(estimated_rpc_matrix - library_rpc_matrix)
			num_genes_converged_library_rpc <- num_genes_converged_library_rpc < 2 * sd_vec
			num_genes_converged_library_rpc <- colSums(num_genes_converged_library_rpc) / nrow(library_rpc_matrix > 0)

			cor_library_estimated <- mean(unlist(lapply(1:10, function(x) cor(estimated_rpc_matrix[, x], library_rpc_matrix[, x]))))


			#print (colSums((abs(estimated_rpc_matrix[library_nzgenes,] - library_rpc_matrix[library_nzgenes,])/ library_rpc_matrix[library_nzgenes,])))
			spike_free_res <- data.frame(total_mRNA=total_mRNA,
								num_genes=num_genes,
								capture_rate=capture_rate,
								reads=reads,
								mRNA_degradation_rate=mRNA_degradation_rate,
								ladder_degradation_rate=ladder_degradation_rate,
								method="spike-free", 
								mode_log_tpm = capture_stats_df$mode_log_tpm,
								mode_molecules_in_cell = capture_stats_df$mode_molecules_in_cell,
								mode_molecules_in_capture = capture_stats_df$mode_molecules_in_capture,
								mse_true_rpc=spike_free_mse_true_rpc,
								mse_library_rpc=spike_free_mse_library_rpc,
								num_genes_converged_true_rpc=num_genes_converged_true_rpc, 
						   		num_genes_converged_library_rpc=num_genes_converged_library_rpc,
								cor_true_estimated=cor_true_estimated,
						   		cor_library_estimated=cor_library_estimated,
						   	    ladder_low_end_tpm=spike_rpc_df$ladder_low_end_tpm,
						   	    ladder_high_end_tpm=spike_rpc_df$ladder_high_end_tpm,
								regression_k=NA,
								regression_b=NA,
								spike_m=spike_m,
								spike_c=spike_c,
							    calibrated_m=calibrated_mc[["m"]],
								calibrated_c=calibrated_mc[["c"]],
								calibrated_mode = mean(calibrated_modes_cell),
								#true_fixed_ladder_c=mean(capture_stats_df$true_fixed_c),
								mode_log_tpm_estimate=NA,
								molecules_for_mode_log_tpm=NA,
								t_estimate=estimate_t(exprs(TPM_cds), relative_expr_thresh=TPM_thresh))

			res <- rbind(spike_res, spike_free_res)
			# res[, c("capture_rate", "num_genes_converged_true_rpc", "num_genes_converged_library_rpc", "method")]

			return(res)
		}else{
			res <- spike_res
		}

	})
	})
	})
	})
	})
	})

	if (return_matrices == FALSE){
		return (full_res)
	}else{
		return (list(true_rpc_matrix=true_rpc_matrix, spike_true_rpc_matrix = spike_true_rpc_matrix, library_rpc_matrix=library_rpc_matrix, tpm_matrix=tpm_matrix, spike_rpc_matrix=spike_rpc_matrix, estimated_rpc_matrix=estimated_rpc_matrix))
	}
}

