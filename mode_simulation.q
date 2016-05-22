library(VGAM)
library(plyr)
library(ggplot2)
library(monocle)
library(dplyr)
library(MASS)
theme_set(theme_bw())

library(devtools)
load_all('/Users/xqiu/Dropbox (Personal)/Projects/monocle-dev')

theme_set(theme_bw(base_size=6))

volume <- 10
dilution <- 40000
my_spike_df <- subset(spike_df, conc_attomoles_ul_Mix1 > 800) # only use spikes above our normal detection threshhold
ladder_molecules <- my_spike_df$conc_attomoles_ul_Mix1 * (volume * 10^(-3) * 1/dilution * 10^(-18) * 6.02214129 * 10^(23))

dmode <- function(x, breaks="Sturges") {
  den <- density(x, kernel=c("gaussian"))
  ( den$x[den$y==max(den$y)] )
}

#############  mode here is calculated based on the downsampled distribution when assuming exp_profile is true transcript count distribution
generate_multinomial_profiles <- function(expr_profile, mRNA_range=c(10000, 25000, 50000, 100000, 250000, 500000, 750000, 1000000)){
	names(mRNA_range) <- mRNA_range
	
	res <- ldply(mRNA_range, function(mRNA_yield) {
		
			T <- rmultinom(1, mRNA_yield, expr_profile)
			#trial_modes <- colHistMode(T)
			tpm = 1e6*T / sum(T)
			mode_log_tpm = estimate_t(tpm)
			mode_molecules = 10^dmode(log10(T[T > 0]))
			data.frame(mRNA_yield=mRNA_yield, tpm=tpm, mode_log_tpm=mode_log_tpm, mode_molecules=mode_molecules, expression=T)
		})
	#res <- cbind(res)
	#colnames(res) <- "mRNA_yield"
	# res$mRNA_yield <- as.numeric(res$mRNA_yield)
	res

}

#original_fpkm_dist <- IFM[,1]
#write.table(IFM[,1], "hypothetical_rna_seq_profile.txt", sep="\t", row.names=F, col.names=F, quote=F)

original_fpkm_dist <- as.vector(read.delim("/Users/xqiu/Dropbox (Personal)/Projects/tmp_BEAM_proj/hypothetical_rna_seq_profile.txt", header=F)[,1])

original_fpkm_dist <- original_fpkm_dist[original_fpkm_dist > 0.01]
original_fpkm_dist <- original_fpkm_dist / sum(original_fpkm_dist)

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
               total_mRNA = rep(total_mRNA, length(split_relative_exprs)), 
               capture_rate = rep(capture_rate, length(split_relative_exprs)),
               reads = rep(reads_per_cell, length(split_relative_exprs)))

		# calibrated_modes_cell <- unlist(lapply(1:length(split_relative_exprs), 
  #              monocle:::calibrate_mode, 
  #              tpm_distribution = split_relative_exprs, 
  #              ladder = ladder, 
  #              total_ladder_transcripts = total_ladder_transcripts,
  #              total_mRNA = rep(total_mRNA, length(split_relative_exprs)), 
  #              capture_rate = rep(capture_rate, length(split_relative_exprs)),
  #              reads = rep(reads_per_cell, length(split_relative_exprs))))

		calibrated_mc_modes_df <- do.call(rbind.data.frame, calibrated_mc_modes)
		mc_model <- rlm(b ~ k, data = calibrated_mc_modes_df)
		calibrated_mc <- list(m = coef(mc_model)[2], c = coef(mc_model)[1])
		calibrated_modes_cell <- calibrated_mc_modes_df[, 1]

		# calibrated_modes_capture <- as.vector(unlist(llply(libraries, function(lib){
		# 	read_sample <- lib$read_sample
		# 	tpm_est <- 1e6*read_sample / sum(read_sample)
		# 	monocle:::calibrate_mode(tpm_est, total_mRNA, capture_rate, reads_per_cell)
		# })))

		# print (calibrated_modes_cell / calibrated_modes_capture)

		capture_stats_df <- ldply(libraries, function(lib){
			E <- lib$E
			M <- lib$M
			S <- lib$S
			L <- lib$L
			T <- c(M, L)


			#true_fixed_c <- mean(log10(captured_ladder[captured_ladder > 0]))
			#ladder_tpm <-  1e6*captured_ladder / (sum(S)) 
			
			#read_sample <- lib$read_sample
			#ladder_tpm <- 1e6*S[names(captured_ladder)] / sum(S)

			hypothetical_ladder <- total_ladder_transcripts * capture_rate * (ladder / sum(ladder))

			hypothetical_ladder_tpm <- 1e6 * hypothetical_ladder / ((total_ladder_transcripts + 1e6) * capture_rate)

			
			no_endo_regression <- rlm (log10(hypothetical_ladder) ~ log10(hypothetical_ladder_tpm))

			#print(qplot(hypothetical_ladder_tpm, hypothetical_ladder, log="xy") + geom_abline(slope=true_fixed_m, intercept=true_fixed_c))
			#print (summary(no_endo_regression))
			#print (true_fixed_m)

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
		print(dim(capture_stats_df))
		spike_rpc_df <- ldply(libraries, function(lib){
			E <- lib$E
			M <- lib$M
			S <- lib$S
			L <- lib$L
			T <- c(M, L)
			read_sample <- lib$read_sample
			tpm_est <- 1e6*read_sample / sum(read_sample)
			mode_log_tpm_estimate = tryCatch({estimate_t(as.matrix(tpm_est))}, error=function(e) {print (e); NA })
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

		spike_rpc_matrix <<- t(laply(libraries, function(lib){
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

		kb_lm <- lm(regression_b ~ regression_k, data=spike_rpc_df)
		spike_m <- coefficients(kb_lm)[2]
		spike_c <- coefficients(kb_lm)[1]

		#print (spike_rpc_df)

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
								#mse_true_rpc=spike_rpc_df$mse_true_rpc,
								#mse_library_rpc=spike_rpc_df$mse_library_rpc,
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
			print(qplot(tpm_matrix[,1], log="x"))
			#print (calibrated_modes)
			estimated_rpc_matrix <- as.matrix(relative2abs(TPM_cds,
											 #total_RNAs=1000000, 
											 estimate_t(exprs(TPM_cds), relative_expr_thresh=TPM_thresh),
											 #alpha_v=spike_rpc_df$molecules_for_mode_log_tpm,
											 #expected_mRNA_mode=calibrated_modes_capture,
											 verbose=F,
											 #kb_slope=calibrated_mc[["m"]],
											 #kb_intercept=calibrated_mc[["c"]],
											 #kb_slope_rng=sort(c(calibrated_mc[["m"]]*0.8, calibrated_mc[["m"]]*1.2)),
											 #c_rng=sort(c(calibrated_mc[["c"]]*0.8, calibrated_mc[["c"]]*1.2)),
											 expected_total_mRNAs = total_mRNA,
											 expected_capture_rate = capture_rate,
											 reads_per_cell=reads,
											 use_fixed_intercept=FALSE,
											 ...))

		    print(qplot(estimated_rpc_matrix[,1], true_rpc_matrix[,1], log="xy") + geom_abline(color="red"))

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
			num_genes_converged_true_rpc <- abs(estimated_rpc_matrix - true_rpc_matrix)
			num_genes_converged_true_rpc <- num_genes_converged_true_rpc / true_rpc_matrix
			num_genes_converged_true_rpc[is.finite(num_genes_converged_true_rpc) == FALSE] <- NA
			num_genes_converged_true_rpc[is.nan(num_genes_converged_true_rpc)] <- NA
			num_genes_converged_true_rpc <- colSums(num_genes_converged_true_rpc < 0.2, na.rm=TRUE) / colSums(true_rpc_matrix > 0) 

			num_genes_converged_library_rpc <- abs(estimated_rpc_matrix - library_rpc_matrix)
			num_genes_converged_library_rpc <- num_genes_converged_library_rpc / library_rpc_matrix
			num_genes_converged_library_rpc[is.finite(num_genes_converged_library_rpc) == FALSE] <- NA
			num_genes_converged_library_rpc[is.nan(num_genes_converged_library_rpc)] <- NA
			num_genes_converged_library_rpc <- colSums(num_genes_converged_library_rpc < 0.2, na.rm=TRUE) / colSums(library_rpc_matrix > 0) 
			

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
								#mse_true_rpc=spike_free_mse_true_rpc,
								#mse_library_rpc=spike_free_mse_library_rpc,
								num_genes_converged_true_rpc=num_genes_converged_true_rpc,
						   		num_genes_converged_library_rpc=num_genes_converged_library_rpc,
						   		#num_true_nzgenes=spike_rpc_df$num_true_nzgenes,
						   	    #num_library_nzgenes=spike_rpc_df$num_library_nzgenes,
						   	    #num_ladder_nzgenes=spike_rpc_df$num_ladder_nzgenes,
						   	    #median_ladder_library_rpc=spike_rpc_df$median_ladder_library_rpc,
						   	    #median_estimated_tpm=spike_rpc_df$median_estimated_tpm,
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
		return (list(true_rpc_matrix=true_rpc_matrix, library_rpc_matrix=library_rpc_matrix, tpm_matrix=tpm_matrix, spike_rpc_matrix=spike_rpc_matrix, estimated_rpc_matrix=estimated_rpc_matrix))
	}
}



spike_sim_rel2abs_by_capture_rate <- simulate_sequencing_with_ladder(original_fpkm_dist, 
												genes_expressed=5000,
												total_mRNAs=125000,
												#capture_rates=c(0.01, 0.05, 0.15, 0.25, seq(0.1, 1, by=0.1)),
												capture_rate = 0.1,
												ladder=my_spike_df$conc_attomoles_ul_Mix1,
												total_ladder_transcripts=round(sum(ladder_molecules)),
											 	mRNA_degradation_rates=0.0, #c(0.0, 0.25, 0.5),
											 	reads_per_cell=1e6,
												#mRNA_degradation_rates=0,
												#mRNA_degradation_rates=0.1,
												#ladder_degradation_rate=0.15,
												ladder_degradation_rates=0,
												cells=10,
												run_rel2abs=TRUE,
												weight_mode=0.33,
												weight_total_rna=0.33,
												weight_relative_expr=0.33)

# calibrate_c(total_mRNA=300000,
# 			capture_rate=1, 
# 			ladder=my_spike_df$conc_attomoles_ul_Mix1,
# 			total_ladder_transcripts=round(sum(my_spike_df$conc_attomoles_ul_Mix1)))

demo_sim <- simulate_sequencing_with_ladder(original_fpkm_dist, return_matrices=TRUE)

pdf("original_fpkm_dist.pdf", width=1.5, height=1.5)
qplot(original_fpkm_dist * 1e6, fill=I("black")) +
	 xlab("Bulk relative abundance (TPM)") +
	 ylab("Genes") +
	 scale_x_log10() +
	 #scale_y_continuous(limits=c(0,500), breaks=c(0,100,200,300,400,500)) +
	 monocle:::monocle_theme_opts()
dev.off()

pdf("demo_true_RPC_matrix.pdf", width=1.5, height=1.5)
qplot(demo_sim$true_rpc_matrix[,1], fill=I("red")) +
	 xlab("Transcripts in cell") +
	 ylab("Genes") +
	 scale_x_log10(limits=c(1,10000), breaks=c(1,10,100, 1000)) +
	 scale_y_continuous(limits=c(0,500), breaks=c(0,100,200,300,400,500)) +
	 monocle:::monocle_theme_opts()
dev.off()

pdf("demo_library_RPC_matrix.pdf", width=1.5, height=1.5)
qplot(demo_sim$library_rpc_matrix[,1], fill=I("steelblue")) +
	 xlab("Transcripts captured in library") +
	 ylab("Genes") +
	 scale_x_log10(limits=c(1,10000), breaks=c(1,10,100, 1000)) +
	 scale_y_continuous(limits=c(0,500), breaks=c(0,100,200,300,400,500)) +
	 monocle:::monocle_theme_opts()
dev.off()

pdf("demo_TPM_matrix.pdf", width=1.5, height=1.5)
qplot(demo_sim$tpm_matrix[,1], log="x", fill=I("darkgrey")) +
	 xlab("Relative abundance (TPM)") +
	 ylab("Genes") +
	 monocle:::monocle_theme_opts()
dev.off()


mode_df <- simulate_sequencing_with_ladder(original_fpkm_dist,
										   genes_expressed=c(2500, 5000, 7500, 10000, 20000),
										   total_mRNAs=c(10000, 25000, 50000, 125000, 250000, 500000, 750000, 1000000),
										   capture_rates=c(0.01, 0.05, 0.15, 0.25, seq(0.1, 1, by=0.1)),
										   reads_per_cell=c(10000,50000,100000,500000,1000000,2000000),
										   ladder=my_spike_df$conc_attomoles_ul_Mix1,
										   total_ladder_transcripts=round(sum(ladder_molecules)),
										   run_rel2abs=FALSE)

#qplot(reads, mode_log_tpm - mode_log_tpm_estimate, color=as.factor(total_mRNAs), geom="line", data=mode_df, log="x")

#qplot(reads, mode_molecules - molecules_for_mode_log_tpm, color=as.factor(total_mRNAs), geom="line", data=mode_df, log="x")

# selected_df <- subset(mode_df, num_genes == 5000)
# qplot(capture_rate, mode_molecules - molecules_for_mode_log_tpm, data=selected_df, log="x") + 
# 	geom_smooth() + 
# 	geom_hline(yintercept=0) +
# 	facet_grid(reads~total_mRNAs, scale="free_y")


#selected_df <- subset(mode_df, reads==1e6)

num_genes_colors <- brewer.pal(length(unique(mode_df$num_genes)) + 3, "Blues")[-c(1:3)]

# This shows that the actual, true most frequent mRNA count in the cell 
# is one until the cell either expresses few genes or a LOT of mRNAs from them
selected_df <- subset(mode_df, reads==1e6 & capture_rate == 1.0)
pdf ("fig1d.pdf", width=2.25, height=1.5)
qplot(total_mRNA, mode_molecules_in_cell, size=I(0.15), color=as.factor(num_genes), data=selected_df, log="x") + 
	geom_smooth(se=F, size=I(0.25)) + scale_color_manual(values=num_genes_colors) +
	geom_hline(yintercept=1, size=I(0.25), color="orange") +
	theme(legend.key.size = unit(0.15, "in")) +
	ylab("Transcripts per cell (mode)") +
	xlab("Total transcripts in cell") +
	guides(color=guide_legend(title="Genes")) +
	theme(legend.key=element_blank())+
	monocle:::monocle_theme_opts()
dev.off()

# This shows when the most frequent mRNA count in the capture as a function of capture efficiency
selected_df <- subset(mode_df, reads==1e6 & total_mRNA == 125000)
pdf ("fig1e.pdf", width=2.25, height=1.5)
qplot(capture_rate, mode_molecules_in_capture, size=I(0.15), color=as.factor(num_genes), data=selected_df) + 
	geom_smooth(se=F, size=I(0.25)) + scale_color_manual(values=num_genes_colors) +
	geom_hline(yintercept=1, size=I(0.25), color="orange") +
	theme(legend.key.size = unit(0.15, "in")) +
	ylab("Transcripts per library (mode)") +
	xlab("Transcript capture rate") +
	scale_x_continuous(labels=scales::percent) + 
	guides(color=guide_legend(title="Genes")) +
	theme(legend.key=element_blank())+
	monocle:::monocle_theme_opts()
dev.off()


# This shows when the molecule count for the mode of log(TPM) distribution as a function of capture efficiency
selected_df <- subset(mode_df, reads==1e6 & total_mRNA == 125000)
pdf ("fig1f.pdf", width=2.35, height=1.5)
qplot(capture_rate, molecules_for_mode_log_tpm, size=I(0.15), color=as.factor(num_genes), data=selected_df) + 
	geom_smooth(se=F, size=I(0.25)) + scale_color_manual(values=num_genes_colors) +
	geom_hline(yintercept=1, size=I(0.25), color="orange") +
    theme(legend.key.size = unit(0.15, "in")) + 
	ylab("Transcript count for\nmost frequent log(FPKM)") +
	xlab("Transcript capture rate") +
	scale_x_continuous(labels=scales::percent) + 
	guides(color=guide_legend(title="Genes")) +
	theme(legend.key=element_blank())+
	monocle:::monocle_theme_opts()
dev.off()


# This shows when the molecule count for the mode of log(TPM) distribution as a function of total mRNAs
selected_df <- subset(mode_df, reads==1e6 & capture_rate == 0.1)
pdf ("molecules_for_mode_log_tpm_by_total_mRNAs.pdf")
qplot(total_mRNA, molecules_for_mode_log_tpm, size=I(0.15), color=as.factor(num_genes), data=selected_df) + 
	geom_smooth(se=F)
dev.off()

# ##### These figures argue explain why the k and b regression parameters lie on a line:
# selected_df <- subset(mode_df, reads==1e6 & num_genes == 5000 & total_mRNA == 250000 & capture_rate %in% c(0.01, 0.1, 0.25, 0.5, 1.0))
# capture_rate_colors <- brewer.pal(length(unique(selected_df$capture_rate)) + 3, "Greens")[-c(1:3)]
# pdf("fig1g.pdf", width=2.35, height=1.5)
# qplot(regression_k, regression_b, size=I(0.15), color=as.factor(capture_rate*100), data=selected_df) + 
# 	geom_smooth(method="lm", size=I(0.25), se=F) +
#     theme(legend.key.size = unit(0.15, "in")) +
#     scale_color_manual(values=capture_rate_colors) +
# 	ylab("Intercept from FPKM vs.\n ERCC transcript counts") +
# 	xlab("Slope from FPKM vs.\n ERCC transcript counts") +
# 	guides(color=guide_legend(title="Capture\n(percent)")) +
# 	theme(legend.key=element_blank())+
# 	monocle:::monocle_theme_opts()
# dev.off()



##### These figures argue explain why the k and b regression parameters lie on a line:
selected_df <- subset(mode_df, reads==1e6 & num_genes == 5000 & total_mRNA == 250000 & capture_rate %in% c(0.01, 0.1, 0.5, 1.0))
capture_rate_colors <- brewer.pal(length(unique(selected_df$capture_rate)) + 3, "Greens")[-c(1:3)]
pdf("fig1g.pdf", width=2.15, height=1.5)
qplot(regression_k, regression_b, size=I(0.15), color=as.factor(capture_rate*100), data=selected_df) + 
	geom_smooth(method="lm", size=I(0.25), se=F) +
    theme(legend.key.size = unit(0.15, "in")) +
    scale_color_manual(values=capture_rate_colors) +
	ylab("Intercept from FPKM vs.\n ERCC transcript counts") +
	xlab("Slope from FPKM vs.\n ERCC transcript counts") +
	guides(color=guide_legend(title="Capture\n(percent)")) +
	theme(legend.key=element_blank())+
	monocle:::monocle_theme_opts()
dev.off()

selected_df <- subset(mode_df, reads==1e6 & num_genes==5000 & capture_rate == 0.1 & total_mRNA %in% c(10000, 50000, 250000, 500000, 1000000))
#capture_rate_colors <- brewer.pal(length(unique(selected_df$capture_rate)) + 3, "Greens")[-c(1:3)]
total_mRNA_colors <- brewer.pal(length(unique(selected_df$total_mRNA)) + 3, "Reds")[-c(1:3)]
pdf("fig1i.pdf", width=2.15, height=1.5)
qplot(regression_k, regression_b, size=I(0.15), color=as.factor(total_mRNA), data=selected_df) + 
	geom_smooth(method="lm", size=I(0.25), se=F) +
    theme(legend.key.size = unit(0.15, "in")) +
    scale_color_manual(values=total_mRNA_colors) +
	ylab("Intercept from FPKM vs.\n ERCC transcript counts") +
	xlab("Slope from FPKM vs.\n ERCC transcript counts") +
	guides(color=guide_legend(title="Total\ntranscripts")) +
	theme(legend.key=element_blank())+
	monocle:::monocle_theme_opts()
dev.off()

# selected_df <- subset(mode_df, num_genes==5000 & total_mRNA == 250000 & capture_rate %in% c(0.01, 0.1, 0.5, 1.0))
# capture_rate_colors <- brewer.pal(length(unique(selected_df$capture_rate)) + 3, "Greens")[-c(1:3)]
# pdf("fig1h.pdf", width=2.15, height=1.35)
# qplot(reads, regression_b, color=as.factor(capture_rate * 100), size=I(0.15),  data=selected_df, log="x") + 
# 	geom_smooth(size=I(0.25), se=F) +
#     theme(legend.key.size = unit(0.15, "in")) +
#     scale_color_manual(values=capture_rate_colors) +
# 	ylab("Intercept from FPKM vs.\n ERCC transcript counts") +
# 	xlab("Reads") +
# 	guides(color=guide_legend(title="Capture\n(percent)")) +
# 	theme(legend.key=element_blank())+
# 	monocle:::monocle_theme_opts()
# dev.off()

selected_df <- subset(mode_df, num_genes==5000 & total_mRNA == 125000 & capture_rate %in% c(0.01, 0.1, 0.5, 1.0))
capture_rate_colors <- brewer.pal(length(unique(selected_df$capture_rate)) + 3, "Greens")[-c(1:3)]
pdf("fig1h.pdf", width=2.15, height=1.35)
qplot(reads, regression_b, color=as.factor(capture_rate * 100), size=I(0.15),  data=selected_df, log="x") + 
	geom_smooth(size=I(0.25), se=F) +
    theme(legend.key.size = unit(0.15, "in")) +
    scale_color_manual(values=capture_rate_colors) +
	ylab("Intercept from FPKM vs.\n ERCC transcript counts") +
	xlab("Reads") +
	guides(color=guide_legend(title="Capture\n(percent)")) +
	theme(legend.key=element_blank())+
	monocle:::monocle_theme_opts()
dev.off()



selected_df <- simulate_sequencing_with_ladder(original_fpkm_dist, 
												total_mRNAs=125000,
												capture_rate=c(0.1, 0.25, 0.5),
												#capture_rate = 0.25,
												ladder=my_spike_df$conc_attomoles_ul_Mix1,
												total_ladder_transcripts=round(sum(ladder_molecules)),
											 	mRNA_degradation_rates=0.0, #c(0.0, 0.25, 0.5),
											 	reads_per_cell=c(1e6),
												#mRNA_degradation_rates=0,
												#mRNA_degradation_rates=0.1,
												#ladder_degradation_rate=0.15,
												#ladder_degradation_rates=0.0,
												cells=100,
												run_rel2abs=FALSE)
#total_mRNA_colors <- brewer.pal(length(unique(selected_df$total_mRNA)) + 3, "Reds")[-c(1:3)]
pdf("fig1h.pdf", width=2.15, height=1.5)
qplot(regression_k, regression_b, size=I(0.2), color=as.numeric(ladder_low_end_tpm/ladder_high_end_tpm), data=selected_df) + 
	#geom_smooth(aes(group=as.factor(capture_rate)), color=I("black"), method="lm", size=I(0.25), se=F) +
    theme(legend.key.size = unit(0.15, "in")) +
    scale_color_distiller(palette="BuPu") +
	ylab("Intercept from FPKM vs.\n ERCC transcript counts") +
	xlab("Slope from FPKM vs.\n ERCC transcript counts") +
	#theme(legend.position="none") +
	guides(color=guide_legend(title="Spike TPM\nlow / high")) +
	theme(legend.key=element_blank())+
	monocle:::monocle_theme_opts()
dev.off()

selected_df <- subset(mode_df, reads==1e6 & num_genes == 5000 & total_mRNA == 125000 & capture_rate == 0.1)
capture_rate_colors <- brewer.pal(length(unique(selected_df$capture_rate)) + 3, "Greens")[-c(1:3)]
#pdf("fig1g.pdf", width=2.15, height=1.5)
qplot(regression_k, regression_b, color=mse_library_rpc, data=selected_df) + 
	#geom_smooth(method="lm", size=I(0.25), se=F) +
    theme(legend.key.size = unit(0.15, "in")) +
    #scale_color_manual(values=capture_rate_colors) +
	ylab("Intercept from FPKM vs.\n ERCC transcript counts") +
	xlab("Slope from FPKM vs.\n ERCC transcript counts") +
	guides(color=guide_legend(title="Capture\n(percent)")) +
	theme(legend.key=element_blank())+
	monocle:::monocle_theme_opts()
#dev.off()


# selected_df <- subset(mode_df, reads==2e6)
# pdf("k_vs_b_by_capture_rate_50k_reads.pdf")
# qplot(regression_k, regression_b, size=I(0.15), color=capture_rate, data=selected_df) + 
# 	facet_grid(total_mRNA~num_genes, scale="free") + geom_smooth(method="lm")
# dev.off()


# mc_df <- as.data.frame(mode_df %>% 
# 		#mutate(captured_mRNAs=total_mRNA * capture_rate) %>%
# 		select(reads, num_genes, total_mRNA, regression_b, regression_k) %>%
# 		group_by(reads, num_genes, total_mRNA) %>% 
# 		do(model=lm(regression_b ~ regression_k, data = .)) %>% 
# 		mutate(c=coefficients(model)[1], m=coefficients(model)[2]) %>%
# 		select(-model))
# qplot(m, c, color=log10(reads/(total_mRNA * capture_rate)), data=as.data.frame(mc_df))

# qplot(log(reads/captured_mRNAs), regression_b/regression_k, size=I(0.15), color=capture_rate, data=mode_df) + 
# 	facet_grid(~total_mRNAs, scale="free")

# qplot(capture_rate, regression_b, data=selected_df) 




# explore_kb_relationship <- function(original_fpkm_distribution,
# 										#genes_expressed=seq(1000, 25000, by=500), 
# 										genes_expressed=c(1000, 2500, 5000, 7500, 10000, 20000),
# 										total_mRNAs=c(10000, 25000, 50000, 100000, 250000, 500000, 750000, 1000000),
# 										replicates=5
# 										){
# 	ldply(1:replicates, function(i) {
# 		ldply(genes_expressed, function(num_genes){
# 			expr_profile <- original_fpkm_dist[sample(length(original_fpkm_dist), num_genes)]
# 			expr_profile <- expr_profile / sum(expr_profile)
# 			ldply(total_mRNAs, function(total_mRNA) {
# 				T <- rmultinom(1, total_mRNA, expr_profile)
# 				#trial_modes <- colHistMode(T)
# 				tpm_true = 1e6*T / sum(T)

# 				mode_log_tpm = estimate_t(tpm_true)
# 				mode_molecules_in_cell = 10^dmode(log10(T[T > 0]))

# 				fit_df <- data.frame(log_T = log10(T), log_tpm=log10(tpm_true))
# 				fit_df <- subset(fit_df, is.finite(log_T) & is.finite(log_tpm))
# 				fit <- lm(log_T ~ log_tpm, data=fit_df)
# 				print (summary(fit))
# 				b <- coefficients(fit)[1]
# 				k <- coefficients(fit)[2]
# 				#print(qplot(tpm_true, T, log="xy"))
# 				data.frame(rep_id=i,
# 						   num_genes=num_genes, 
# 						   total_mRNAs=total_mRNA, 
# 						   regression_k=k,
# 						   regression_b=b)
# 			})
# 		})
# 	})
# }

# kb_df <- explore_kb_relationship(original_fpkm_dist)


# qplot(regression_k, regression_b, size=I(0.15), data=kb_df)


# selected_df <- subset(mode_df, reads==1e6 & num_genes == 5000)
# qplot(capture_rate, mode_molecules - molecules_for_mode_log_tpm, size=I(0.15), data=selected_df) + 
# 	#geom_smooth() + 
# 	#geom_hline(yintercept=0, color="red") +
# 	facet_grid(num_genes~total_mRNAs, scale="free_y")

# qplot(capture_rate, mode_molecules, size=I(0.15), data=selected_df) + 
# 	#geom_smooth() + 
# 	#geom_hline(yintercept=0, color="red") +
# 	facet_grid(num_genes~total_mRNAs, scale="free_y")

# qplot(capture_rate, tpm_mse, size=I(0.15),  data=selected_df) + 
# 	scale_y_sqrt()+
# 	#geom_smooth() + 
# 	#geom_hline(yintercept=0, color="red") +
# 	facet_grid(num_genes~total_mRNAs, scale="free_y")


spike_sim_rel2abs_df <- simulate_sequencing_with_ladder(original_fpkm_dist, 
												total_mRNAs=300000,
												capture_rate=c(0.1, 0.25, 0.5),
												#capture_rate = 0.25,
												ladder=my_spike_df$conc_attomoles_ul_Mix1,
												total_ladder_transcripts=round(sum(ladder_molecules)),
											 	mRNA_degradation_rates=0.0, #c(0.0, 0.25, 0.5),
											 	reads_per_cell=c(10000,50000,100000,250000,500000,1000000,2000000),
												#mRNA_degradation_rates=0,
												#mRNA_degradation_rates=0.1,
												#ladder_degradation_rate=0.15,
												ladder_degradation_rates=c(0.0, 0.01, 0.1, 0.2, seq(0.4, 0.8, by=0.2)),
												cells=10,
												run_rel2abs=TRUE)


qplot(regression_k, regression_b, size=I(0.15), data=spike_sim_rel2abs_df) + facet_grid(mRNA_degradation_rate~capture_rate)

selected_df <- subset(spike_sim_rel2abs_df, reads==1e6 & capture_rate == 0.1)
pdf("Supplemental_Figure_X4a.pdf", width=2, height=2)
qplot(ladder_degradation_rate, num_genes_converged_library_rpc, data=selected_df, size=I(0.15), color=method, fill=method) + 
	xlab("Ladder degradation") +
	ylab("Genes correctly quantified") +
	scale_x_continuous(label=scales::percent) +
	scale_y_continuous(limits=c(0,1), label=scales::percent) +
	#facet_grid(~capture_rate) + 
	scale_color_brewer(palette="Set1") +
	scale_fill_brewer(palette="Set1") +
	stat_summary(fun.data = "mean_cl_boot", size=0.25) +
	stat_summary(fun.data = "mean_cl_boot", size=0.35, geom="line") +
    stat_summary(color="black", fun.data = "mean_cl_boot", size=0.75, geom="point") +

	theme(legend.position="none") +
	monocle:::monocle_theme_opts()
dev.off()


spike_sim_rel2abs_by_capture_rate <- simulate_sequencing_with_ladder(original_fpkm_dist, 
												total_mRNAs=125000,
												capture_rates=c(0.01, 0.05, 0.15, 0.25, seq(0.1, 1, by=0.1)),
												#capture_rates = 0.5,
												ladder=my_spike_df$conc_attomoles_ul_Mix1,
												total_ladder_transcripts=round(sum(ladder_molecules)),
											 	mRNA_degradation_rates=0.0, #c(0.0, 0.25, 0.5),
											 	reads_per_cell=1e6,
												#mRNA_degradation_rates=0,
												#mRNA_degradation_rates=0.1,
												#ladder_degradation_rate=0.15,
												ladder_degradation_rates=0,
												cells=10,
												run_rel2abs=TRUE)

#qplot(capture_rate, spike_c, data=spike_sim_rel2abs_by_capture_rate)
#qplot(regression_k, regression_b, color=capture_rate, data=spike_sim_rel2abs_by_capture_rate) + geom_abline(slope=-1, intercept=spike_sim_rel2abs_by_capture_rate$true_fixed_ladder_c)

#qplot(regression_k, regression_b, color=log10(estimated_ladder_tpm), data=spike_sim_rel2abs_by_capture_rate)

#selected_df <- subset(spike_sim_rel2abs_df)
pdf("Supplemental_Figure_X4b.pdf", width=2, height=2)
qplot(capture_rate, num_genes_converged_library_rpc, fill=method, color=method, size=I(0.15), data=spike_sim_rel2abs_by_capture_rate) +  
	xlab("Capture") +
	ylab("Genes correctly quantified") +
	scale_x_continuous(label=scales::percent) +
	scale_y_continuous(limits=c(0,1),label=scales::percent) +
	scale_color_brewer(palette="Set1") +
	scale_fill_brewer(palette="Set1") +
	#geom_smooth(size=0.25) +
    stat_summary(fun.data = "mean_cl_boot", size=0.25) +
	stat_summary(fun.data = "mean_cl_boot", size=0.35, geom="line") +
    stat_summary(color="black", fun.data = "mean_cl_boot", size=0.75, geom="point") +
	theme(legend.position="none") +
	monocle:::monocle_theme_opts()
dev.off()

spike_sim_rel2abs_by_reads <- simulate_sequencing_with_ladder(original_fpkm_dist, 
												total_mRNAs=300000,
												#capture_rates=c(0.01, 0.05, 0.15, 0.25, seq(0.1, 1, by=0.1)),
												capture_rates=0.1,
												#capture_rate = 0.25,
												ladder=my_spike_df$conc_attomoles_ul_Mix1,
												total_ladder_transcripts=round(sum(ladder_molecules)),
											 	mRNA_degradation_rates=0.0, #c(0.0, 0.25, 0.5),
											 	reads_per_cell=c(10000,50000,100000,250000,500000,1000000,2000000, 3e6, 8e6),

												#mRNA_degradation_rates=0,
												#mRNA_degradation_rates=0.1,
												#ladder_degradation_rate=0.15,
												ladder_degradation_rates=0,
												cells=10,
												run_rel2abs=TRUE)
#selected_df <- subset(spike_sim_rel2abs_df, capture_rate == 0.1 & ladder_degradation_rate== 0.0)
pdf("Supplemental_Figure_X4c.pdf", width=2, height=2)
qplot(reads, num_genes_converged_library_rpc, fill=method, color=method, size=I(0.15), data=spike_sim_rel2abs_by_reads) + 
	xlab("Reads") +
	ylab("Genes correctly quantified") +
	scale_x_log10(breaks=c(10000,100000, 1000000, 8e6)) +
	scale_y_continuous(limits=c(0,1),label=scales::percent) +
	scale_color_brewer(palette="Set1") +
	scale_fill_brewer(palette="Set1") +
    stat_summary(fun.data = "mean_cl_boot", size=0.25) +
	stat_summary(fun.data = "mean_cl_boot", size=0.35, geom="line") +
    stat_summary(color="black", fun.data = "mean_cl_boot", size=0.75, geom="point") +
	theme(legend.position="none") +
	monocle:::monocle_theme_opts()
dev.off()


#selected_df <- subset(spike_sim_rel2abs_by_reads, reads == 1e6)
qplot(regression_k, regression_b, color=estimated_ladder_tpm, data=spike_sim_rel2abs_by_reads) + 
	geom_smooth(method="lm", size=I(0.25), se=F) +
    theme(legend.key.size = unit(0.15, "in")) +
    #scale_color_manual(values=total_mRNA_colors) +
	ylab("Intercept from FPKM vs.\n ERCC transcript counts") +
	xlab("Slope from FPKM vs.\n ERCC transcript counts") +
	guides(color=guide_legend(title="Total\ntranscripts")) +
	theme(legend.key=element_blank())+
	monocle:::monocle_theme_opts()

#qplot(mRNA_degradation_rate, mse, data=spike_sim_rel2abs_df, color=method) + scale_y_sqrt()



spike_sim_rel2abs_by_reads <- simulate_sequencing_with_ladder(original_fpkm_dist, 
												total_mRNAs=300000,
												#capture_rates=c(0.01, 0.05, 0.15, 0.25, seq(0.1, 1, by=0.1)),
												capture_rates=1,
												#capture_rate = 0.25,
												ladder=my_spike_df$conc_attomoles_ul_Mix1,
												total_ladder_transcripts=round(sum(ladder_molecules)),
											 	mRNA_degradation_rates=0.9, #c(0.0, 0.25, 0.5),
											 	reads_per_cell=2e5,
												#mRNA_degradation_rates=0,
												#mRNA_degradation_rates=0.9,
												#ladder_degradation_rate=0.15,
												ladder_degradation_rates=0.1,
												cells=10,
												run_rel2abs=TRUE, return_matrices=TRUE)
qplot(rowMeans(spike_sim_rel2abs_by_reads$spike_rpc_matrix), rowMeans(spike_sim_rel2abs_by_reads$estimated_rpc_matrix), log="xy") + geom_abline()


# selected_df <- subset(mode_df, num_genes==5000 & reads==1e6 & capture_rate == 0.1 & total_mRNA == 250000)
# qplot(ladder_degradation_rate, mse, data=selected_df, color=method) + 
# 	facet_grid(~capture_rate) + 
# 	scale_y_log10() + 
# 	geom_smooth()


#####################################

# Pick around 5,000 transcripts to "express" from this hypothetical cell, consistent
# with what we typically see with current gen SC-RNA-Seq
reduced_fpkm_dist_5k <- original_fpkm_dist[sample(length(original_fpkm_dist), 5000)]
reduced_fpkm_dist_5k <- reduced_fpkm_dist_5k / sum(reduced_fpkm_dist_5k)

profiles <- generate_multinomial_profiles(reduced_fpkm_dist_5k)
#pdf("RPC_dist_by_yield.pdf", width=5, height=2)
qplot(expression, geom="histogram", data=profiles, log="x") + 
	facet_wrap(~mRNA_yield, ncol=4) + 
	geom_vline(aes(xintercept=mode_molecules))+
	scale_x_log10(breaks=c(1, 10, 100, 1000)) +
	 theme(axis.text.y=element_text(size=6)) +
	 theme(axis.text.x=element_text(size=6)) +
	 theme(axis.title.y=element_text(size=6)) +
	 theme(axis.title.x=element_text(size=6)) +
	 xlab("Expression (RPC)") +
	 ylab("Transcript isoforms") +
	 theme(panel.border = element_blank(), axis.line = element_line()) +
	 theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
	 theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
 	 theme(panel.background=element_blank()) +
	 theme(legend.position = "none") +
	 theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
 	 theme(strip.text = element_text(size=6))
#dev.off()

qplot(tpm, geom="histogram", data=profiles, log="x") + 
	facet_wrap(~mRNA_yield, ncol=4) + 
	geom_vline(aes(xintercept=mode_log_tpm))+
	scale_x_log10(breaks=c(1, 10, 100, 1000)) +
	 theme(axis.text.y=element_text(size=6)) +
	 theme(axis.text.x=element_text(size=6)) +
	 theme(axis.title.y=element_text(size=6)) +
	 theme(axis.title.x=element_text(size=6)) +
	 xlab("Expression (RPC)") +
	 ylab("Transcript isoforms") +
	 theme(panel.border = element_blank(), axis.line = element_line()) +
	 theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
	 theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
 	 theme(panel.background=element_blank()) +
	 theme(legend.position = "none") +
	 theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
 	 theme(strip.text = element_text(size=6))

pdf("RPC_dist_sampled.pdf", width=2, height=2)
qplot(reduced_fpkm_dist_5k, geom="histogram", log="x") +
theme(axis.text.y=element_text(size=6)) +
	 theme(axis.text.x=element_text(size=6)) +
	 theme(axis.title.y=element_text(size=6)) +
	 theme(axis.title.x=element_text(size=6)) +
	 xlab("Relative abundance") +
	 ylab("Genes") +
	 theme(panel.border = element_blank(), axis.line = element_line()) +
	 theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
	 theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
 	 theme(panel.background=element_blank()) +
	 theme(legend.position = "none") +
	 theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
 	 theme(strip.text = element_text(size=6))
dev.off()

pdf("RPC_dist_orig.pdf", width=2, height=2)
qplot(original_fpkm_dist, geom="histogram", log="x") +
theme(axis.text.y=element_text(size=6)) +
	 theme(axis.text.x=element_text(size=6)) +
	 theme(axis.title.y=element_text(size=6)) +
	 theme(axis.title.x=element_text(size=6)) +
	 xlab("Relative abundance") +
	 ylab("Genes") +
	 theme(panel.border = element_blank(), axis.line = element_line()) +
	 theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
	 theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
 	 theme(panel.background=element_blank()) +
	 theme(legend.position = "none") +
	 theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
 	 theme(strip.text = element_text(size=6))
dev.off()

fpkm_dist_df <- rbind(data.frame(type="Bulk cell population", expression=original_fpkm_dist),
					  data.frame(type="Hypothetical single cell", expression=reduced_fpkm_dist_5k))
pdf("RPC_dist_combined.pdf", width=2, height=2)
qplot(expression, fill=type, data=fpkm_dist_df, geom="histogram", log="x") +
facet_wrap(~type, ncol=1, scales="free_y")+
theme(axis.text.y=element_text(size=6)) +
	 theme(axis.text.x=element_text(size=6)) +
	 theme(axis.title.y=element_text(size=6)) +
	 theme(axis.title.x=element_text(size=6)) +
	 xlab("Relative abundance") +
	 ylab("Transcript isoforms") +
	 theme(panel.border = element_blank(), axis.line = element_line()) +
	 theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
	 theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
 	 theme(panel.background=element_blank()) +
	 theme(legend.position = "none") +
	 theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
 	 theme(strip.text = element_text(size=6)) + scale_fill_brewer(palette="Set1")
dev.off()


