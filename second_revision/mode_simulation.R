library(VGAM)
library(plyr)
library(ggplot2)
library(monocle)
library(dplyr)
library(MASS)
theme_set(theme_bw())
# use relative2abs from debug_monocle2
#library(devtools)
#load_all('/Users/xqiu/Dropbox (Personal)/Projects/monocle-dev') 
library(monocle)
source('./simulator.R')
theme_set(theme_bw(base_size=6))

volume <- 10
dilution <- 40000
#now we can estimate the transcript counts corresponding to the mode of TPM distribution, no subset is required 
my_spike_df <- spike_df
my_spike_df <- subset(my_spike_df, conc_attomoles_ul_Mix1 > 800) # only use spikes above our normal detection threshhold
ladder_molecules <- my_spike_df$conc_attomoles_ul_Mix1 * (volume * 10^(-3) * 1/dilution * 10^(-18) * 6.02214129 * 10^(23))

#original_fpkm_dist <- IFM[,1]
#write.table(IFM[,1], "hypothetical_rna_seq_profile.txt", sep="\t", row.names=F, col.names=F, quote=F)

#read the bulk RNA-seq distribution
#use a true single-cell RNA-seq dataset (maybe pooled single-cell RNA-seq dataset?) 
original_fpkm_dist <- as.vector(read.delim("/Users/xqiu/Dropbox (Personal)/Projects/tmp_BEAM_proj/hypothetical_rna_seq_profile.txt", header=F)[,1])

original_fpkm_dist <- original_fpkm_dist[original_fpkm_dist > 0.01]
original_fpkm_dist <- original_fpkm_dist / sum(original_fpkm_dist) * 10^6

dmode <- function(x, breaks="Sturges") {
  if (length(x) < 2) return (0);
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
	res

}

demo_sim <- simulate_sequencing_with_ladder(original_fpkm_dist, return_matrices=TRUE)

pdf("./tmp/original_fpkm_dist.pdf", width=1.5, height=1.5)
qplot(original_fpkm_dist * 1e6, fill=I("black")) +
	 xlab("Bulk relative abundance (TPM)") +
	 ylab("Genes") +
	 scale_x_log10() +
	 #scale_y_continuous(limits=c(0,500), breaks=c(0,100,200,300,400,500)) +
	 monocle:::monocle_theme_opts()
dev.off()

pdf("./tmp/demo_true_RPC_matrix.pdf", width=1.5, height=1.5)
qplot(demo_sim$true_rpc_matrix[,1], fill=I("red")) +
	 xlab("Transcripts in cell") +
	 ylab("Genes") +
	 scale_x_log10(limits=c(1,10000), breaks=c(1,10,100, 1000)) +
	 scale_y_continuous(limits=c(0,500), breaks=c(0,100,200,300,400,500)) +
	 monocle:::monocle_theme_opts()
dev.off()

pdf("./tmp/demo_library_RPC_matrix.pdf", width=1.5, height=1.5)
qplot(demo_sim$library_rpc_matrix[,1], fill=I("steelblue")) +
	 xlab("Transcripts captured in library") +
	 ylab("Genes") +
	 scale_x_log10(limits=c(1,10000), breaks=c(1,10,100, 1000)) +
	 scale_y_continuous(limits=c(0,500), breaks=c(0,100,200,300,400,500)) +
	 monocle:::monocle_theme_opts()
dev.off()

pdf("./tmp/demo_TPM_matrix.pdf", width=1.5, height=1.5)
qplot(demo_sim$tpm_matrix[,1], log="x", fill=I("darkgrey")) +
	 xlab("Relative abundance (TPM)") +
	 ylab("Genes") +
	 monocle:::monocle_theme_opts()
dev.off()


mode_df <- simulate_sequencing_with_ladder(original_fpkm_dist,
                                genes_expressed=c(2500, 5000, 7500, 10000, 20000),
                                total_mRNAs=c(10000, 25000, 37500, 50000, 100000, 150000, 250000, 500000, 750000, 1000000),
                                capture_rates=c(0.01, 0.05, 0.15, 0.25, seq(0.1, 1, by=0.1)),
                                reads_per_cell=c(500000,1000000,2000000), #10000,50000,100000,
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
selected_df <- subset(mode_df, reads==1e6 & total_mRNA == 100000)
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
selected_df <- subset(mode_df, reads==1e6 & total_mRNA == 100000)
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
selected_df <- subset(mode_df, reads==1e6 & capture_rate == 0.25)
pdf ("molecules_for_mode_log_tpm_by_total_mRNAs.pdf")
qplot(total_mRNA, molecules_for_mode_log_tpm, size=I(0.15), color=as.factor(num_genes), data=selected_df) + 
	geom_smooth(se=F)
dev.off()

##### These figures argue explain why the k and b regression parameters lie on a line:
selected_df <- subset(mode_df, reads==1e6 & num_genes == 5000 & total_mRNA == 100000 & capture_rate %in% c(0.01, 0.1, 0.25, 0.5, 1.0))
capture_rate_colors <- brewer.pal(length(unique(selected_df$capture_rate)) + 3, "Greens")[-c(1:3)]
pdf("./main_figures/fig1g.pdf", width=2.35, height=1.5)
qplot(regression_k, regression_b, size=I(0.5), color=as.factor(capture_rate*100), data=selected_df) + 
	geom_smooth(method="lm", size=I(0.25), se=F) +
    theme(legend.key.size = unit(0.15, "in")) +
    scale_color_manual(values=capture_rate_colors) +
	ylab("Intercept from FPKM vs.\n ERCC transcript counts") +
	xlab("Slope from FPKM vs.\n ERCC transcript counts") +
	guides(color=guide_legend(title="Capture\n(percent)")) +
	theme(legend.key=element_blank())+
	monocle:::monocle_theme_opts() #+ 
dev.off()

selected_df <- subset(mode_df, reads==1e6 & num_genes==5000 & capture_rate == 0.25 & total_mRNA %in% c(10000, 37500, 50000, 100000, 150000, 500000, 1000000))
#capture_rate_colors <- brewer.pal(length(unique(selected_df$capture_rate)) + 3, "Greens")[-c(1:3)]
total_mRNA_colors <- brewer.pal(length(unique(selected_df$total_mRNA)) + 3, "Reds")[-c(1:2)]
pdf("./main_figures/fig1i.pdf", width=2.15, height=1.5)
qplot(regression_k, regression_b, size=I(0.5), color=as.factor(total_mRNA), data=selected_df) + 
	geom_smooth(method="lm", size=I(0.25), se=F) +
    theme(legend.key.size = unit(0.15, "in")) +
    scale_color_manual(values=total_mRNA_colors) +
	ylab("Intercept from FPKM vs.\n ERCC transcript counts") +
	xlab("Slope from FPKM vs.\n ERCC transcript counts") +
	guides(color=guide_legend(title="Total\ntranscripts")) +
	theme(legend.key=element_blank())+
	monocle:::monocle_theme_opts()
dev.off()

selected_df$fraction_spike_in <- sum(ladder_molecules) / (sum(ladder_molecules) + selected_df$total_mRNA)
selected_df$fraction_spike_in <- round(selected_df$fraction_spike_in, 2)
pdf("./main_figures/fig1i_fraction.pdf", width=2.15, height=1.5)
qplot(regression_k, regression_b, size=I(0.15), color=as.factor(fraction_spike_in), data=selected_df) + 
	geom_smooth(method="lm", size=I(0.25), se=F) +
    theme(legend.key.size = unit(0.15, "in")) +
    scale_color_manual(values=total_mRNA_colors) +
	ylab("Intercept from FPKM vs.\n ERCC transcript counts") +
	xlab("Slope from FPKM vs.\n ERCC transcript counts") +
	guides(color=guide_legend(title="Fraction of spike-in")) +
	theme(legend.key=element_blank())+
	monocle:::monocle_theme_opts()
dev.off()

selected_df <- subset(mode_df, num_genes==5000 & total_mRNA == 250000 & capture_rate %in% c(0.01, 0.1, 0.5, 1.0))
capture_rate_colors <- brewer.pal(length(unique(selected_df$capture_rate)) + 3, "Greens")[-c(1:3)]
pdf("./main_figures/fig1h_intercept.pdf", width=2.15, height=1.35)
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

pdf("./main_figures/fig1h_slope.pdf", width=2.15, height=1.35)
qplot(reads, regression_k, color=as.factor(capture_rate * 100), size=I(0.15),  data=selected_df, log="x") + 
	geom_smooth(size=I(0.25), se=F) +
    theme(legend.key.size = unit(0.15, "in")) +
    scale_color_manual(values=capture_rate_colors) +
	ylab("Slope from FPKM vs.\n ERCC transcript counts") +
	xlab("Reads") +
	guides(color=guide_legend(title="Capture\n(percent)")) +
	theme(legend.key=element_blank())+
	monocle:::monocle_theme_opts()
dev.off()

selected_df <- subset(mode_df, reads==1e6 & total_mRNA == 100000 & capture_rate %in% 0.25)
num_genes_colors <- brewer.pal(length(unique(selected_df$num_genes)) + 3, "Greens")[-c(1:3)]
pdf("./main_figures/fig1_num_genes.pdf", width=2.35, height=1.5)
qplot(regression_k, regression_b, size=I(0.5), color=as.factor(num_genes), data=selected_df) + 
  geom_smooth(method="lm", size=I(0.25), se=F) +
  theme(legend.key.size = unit(0.15, "in")) +
  scale_color_manual(values=num_genes_colors) +
  ylab("Intercept from FPKM vs.\n ERCC transcript counts") +
  xlab("Slope from FPKM vs.\n ERCC transcript counts") +
  guides(color=guide_legend(title="Number of \n expressed genes")) +
  theme(legend.key=element_blank())+
  monocle:::monocle_theme_opts()
dev.off()

selected_df <- simulate_sequencing_with_ladder(original_fpkm_dist, 
												total_mRNAs=100000,
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
pdf("./main_figures/fig1h.pdf", width=2.15, height=1.5)
qplot(regression_k, regression_b, size=I(0.5), color=as.numeric(ladder_low_end_tpm/ladder_high_end_tpm), data=selected_df) + 
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
#pdf("./main_figures/fig1g.pdf", width=2.15, height=1.5)
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


selected_df <- subset(mode_df, reads==2e6)
pdf("./main_figures/k_vs_b_by_capture_rate_50k_reads.pdf")
qplot(regression_k, regression_b, size=I(0.15), color=capture_rate, data=selected_df) + 
	facet_grid(total_mRNA~num_genes, scale="free") + geom_smooth(method="lm")
dev.off()

total_mRNAs <- 100000 #* 0.05

spike_sim_rel2abs_df <- simulate_sequencing_with_ladder(original_fpkm_dist, 
												total_mRNAs=total_mRNAs, #
												capture_rate=c(0.1, 0.25, 0.5),
												#capture_rate = 0.25,
												ladder=spike_df$conc_attomoles_ul_Mix1,
												total_ladder_transcripts=round(sum(spike_df$rounded_numMolecules)),
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
pdf("./supplementary_figures/Supplemental_Figure_X4a.pdf", width=2, height=2)
qplot(ladder_degradation_rate, num_genes_converged_true_rpc, data=selected_df, size=I(0.15), color=method, fill=method) + 
	xlab("Ladder degradation") +
	ylab("Genes correctly quantified") + nm_theme() + 
	scale_x_continuous(label=scales::percent) +
	scale_y_continuous(limits=c(0,1), label=scales::percent) +
	#facet_grid(~capture_rate) + 
	scale_color_brewer(palette="Set1") +
	scale_fill_brewer(palette="Set1") +
	stat_summary(fun.data = "mean_cl_boot", size=0.25) +
	stat_summary(fun.data = "mean_cl_boot", size=0.35, geom="line") +
    stat_summary(color="black", fun.data = "mean_cl_boot", size=0.75, geom="point") 

	# theme(legend.position="none") +
	# monocle:::monocle_theme_opts()
dev.off()

pdf("./supplementary_figures/Supplemental_Figure_X_deg_cor.pdf", width=2, height=2)
qplot(ladder_degradation_rate, cor_true_estimated, fill=method, color=method, size=I(0.15), data=selected_df) +  
	xlab("Ladder degradation") +
	ylab("Correlation between estimated transcript counts \n and true transcript count") +
	scale_x_continuous(label=scales::percent) +
	scale_y_continuous(label=scales::percent) + #limits=c(0,1),
	scale_color_brewer(palette="Set1") +
	scale_fill_brewer(palette="Set1") +
	#geom_smooth(size=0.25) +
    stat_summary(fun.data = "mean_cl_boot", size=0.25) +
	stat_summary(fun.data = "mean_cl_boot", size=0.35, geom="line") +
    stat_summary(color="black", fun.data = "mean_cl_boot", size=0.75, geom="point") + nm_theme()
dev.off()

spike_sim_rel2abs_by_capture_rate <- simulate_sequencing_with_ladder(original_fpkm_dist, 
												total_mRNAs=total_mRNAs, #  * 0.25,
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
pdf("./supplementary_figures/Supplemental_Figure_X4b.pdf", width=2, height=2)
qplot(capture_rate, num_genes_converged_true_rpc, fill=method, color=method, size=I(0.15), data=spike_sim_rel2abs_by_capture_rate) +  
	xlab("Capture") +
	ylab("Genes correctly quantified") +
	scale_x_continuous(label=scales::percent) +
	scale_y_continuous(limits=c(0,1),label=scales::percent) + #
	scale_color_brewer(palette="Set1") +
	scale_fill_brewer(palette="Set1") +
	#geom_smooth(size=0.25) +
    stat_summary(fun.data = "mean_cl_boot", size=0.25) +
	stat_summary(fun.data = "mean_cl_boot", size=0.35, geom="line") +
    stat_summary(color="black", fun.data = "mean_cl_boot", size=0.75, geom="point") +
	nm_theme()
dev.off()

pdf("./supplementary_figures/Supplemental_Figure_X4b_2.pdf", width=2, height=2)
qplot(capture_rate, cor_true_estimated, fill=method, color=method, size=I(0.15), data=spike_sim_rel2abs_by_capture_rate) +  
	xlab("Capture") +
	ylab("Correlation between estimated transcript counts \n and true transcript count") +
	scale_x_continuous(label=scales::percent) +
	scale_y_continuous(label=scales::percent) + #limits=c(0,1),
	scale_color_brewer(palette="Set1") +
	scale_fill_brewer(palette="Set1") +
	#geom_smooth(size=0.25) +
    stat_summary(fun.data = "mean_cl_boot", size=0.25) +
	stat_summary(fun.data = "mean_cl_boot", size=0.35, geom="line") +
    stat_summary(color="black", fun.data = "mean_cl_boot", size=0.75, geom="point") +
	nm_theme()
dev.off()

pdf("./supplementary_figures/Supplemental_Figure_X4b_3.pdf", width=2, height=2)
qplot(capture_rate, mse_true_rpc, fill=method, color=method, size=I(0.15), data=spike_sim_rel2abs_by_capture_rate, log = 'y') +  
	xlab("Capture") +
	ylab("MSE between estimated transcript counts \n and true transcript count") +
	scale_x_continuous(label=scales::percent) +
	# scale_y_continuous(label=scales::percent) + #limits=c(0,1),
	scale_color_brewer(palette="Set1") +
	scale_fill_brewer(palette="Set1") +
	#geom_smooth(size=0.25) +
    stat_summary(fun.data = "mean_cl_boot", size=0.25) +
	stat_summary(fun.data = "mean_cl_boot", size=0.35, geom="line") +
    stat_summary(color="black", fun.data = "mean_cl_boot", size=0.75, geom="point") +
	nm_theme()
dev.off()

pdf("./supplementary_figures/Supplemental_Figure_X4b_3.pdf", width=2, height=2)
qplot(capture_rate, mse_true_rpc, fill=method, color=method, size=I(0.15), data=spike_sim_rel2abs_by_capture_rate, log = 'y') +  
	xlab("Capture") +
	ylab("MSE between estimated transcript counts \n and true transcript count") +
	scale_x_continuous(label=scales::percent) +
	# scale_y_continuous(label=scales::percent) + #limits=c(0,1),
	scale_color_brewer(palette="Set1") +
	scale_fill_brewer(palette="Set1") +
	#geom_smooth(size=0.25) +
    stat_summary(fun.data = "mean_cl_boot", size=0.25) +
	stat_summary(fun.data = "mean_cl_boot", size=0.35, geom="line") +
    stat_summary(color="black", fun.data = "mean_cl_boot", size=0.75, geom="point") +
	nm_theme()
dev.off()

pdf("./supplementary_figures/Supplemental_Figure_X4b_4.pdf", width=2, height=2)
qplot(capture_rate, round(molecules_for_mode_log_tpm), fill=method, color=method, size=I(0.15), data=spike_sim_rel2abs_by_capture_rate, log = 'y') +  
	xlab("Capture") +
	ylab("Transcript counts \n corresponding to TPM mode") +
	scale_x_continuous(label=scales::percent) +
	# scale_y_continuous(label=scales::percent) + #limits=c(0,1),
	scale_color_brewer(palette="Set1") +
	scale_fill_brewer(palette="Set1") +
	#geom_smooth(size=0.25) +
    stat_summary(fun.data = "mean_cl_boot", size=0.25) +
	stat_summary(fun.data = "mean_cl_boot", size=0.35, geom="line") +
    stat_summary(color="black", fun.data = "mean_cl_boot", size=0.75, geom="point") +
	nm_theme()
dev.off()

subset_spike_sim_rel2abs_by_capture_rate <- spike_sim_rel2abs_by_capture_rate[, c('capture_rate', 'method', 'molecules_for_mode_log_tpm', 'mode_molecules_in_capture')]
mlt_df <- melt(subset_spike_sim_rel2abs_by_capture_rate, id.vars = c('capture_rate', 'method'), measure.vars = c('molecules_for_mode_log_tpm', 'mode_molecules_in_capture'))

pdf("./supplementary_figures/Supplemental_Figure_X4b_4_all.pdf", width=2, height=2)
qplot(capture_rate, round(value), fill=method, color=method, linetype = variable, size=I(0.15), data=mlt_df, log = 'y') +  
  xlab("Capture") +
  ylab("Transcript counts \n corresponding to TPM mode") +
  scale_x_continuous(label=scales::percent) +
  # scale_y_continuous(label=scales::percent) + #limits=c(0,1),
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  #geom_smooth(size=0.25) +
  stat_summary(fun.data = "mean_cl_boot", size=0.25) +
  stat_summary(fun.data = "mean_cl_boot", size=0.35, geom="line") +
  stat_summary(color="black", fun.data = "mean_cl_boot", size=0.75, geom="point") +
  nm_theme()
dev.off()

pdf("./supplementary_figures/Supplemental_Figure_X4b_4_all.pdf", width=2, height=2)
qplot(capture_rate, round(value), fill=method, color=method, linetype = variable, size=I(0.15), data=mlt_df, log = 'y') +  
  xlab("Capture") +
  ylab("Transcript counts \n corresponding to TPM mode") +
  scale_x_continuous(label=scales::percent) +
  # scale_y_continuous(label=scales::percent) + #limits=c(0,1),
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  #geom_smooth(size=0.25) +
  stat_summary(fun.data = "mean_cl_boot", size=0.25) +
  stat_summary(fun.data = "mean_cl_boot", size=0.35, geom="line") +
  stat_summary(color="black", fun.data = "mean_cl_boot", size=0.75, geom="point") #  nm_theme()
dev.off()

spike_sim_rel2abs_by_reads <- simulate_sequencing_with_ladder(original_fpkm_dist, 
												total_mRNAs=total_mRNAs ,
												#capture_rates=c(0.01, 0.05, 0.15, 0.25, seq(0.1, 1, by=0.1)),
												capture_rates=0.25,
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
pdf("./supplementary_figures/Supplemental_Figure_X4c.pdf", width=2, height=2)
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
	nm_theme()
dev.off()


# selected_df <- subset(spike_sim_rel2abs_by_reads, reads == 1e6)
# qplot(regression_k, regression_b, color=estimated_ladder_tpm, data=spike_sim_rel2abs_by_reads) + 
# 	geom_smooth(method="lm", size=I(0.25), se=F) +
#     theme(legend.key.size = unit(0.15, "in")) +
#     #scale_color_manual(values=total_mRNA_colors) +
# 	ylab("Intercept from FPKM vs.\n ERCC transcript counts") +
# 	xlab("Slope from FPKM vs.\n ERCC transcript counts") +
# 	guides(color=guide_legend(title="Total\ntranscripts")) +
# 	theme(legend.key=element_blank())+
# 	monocle:::monocle_theme_opts()

#qplot(mRNA_degradation_rate, mse, data=spike_sim_rel2abs_df, color=method) + scale_y_sqrt()



# spike_sim_rel2abs_by_reads <- simulate_sequencing_with_ladder(original_fpkm_dist, 
# 												total_mRNAs=total_mRNAs,
# 												#capture_rates=c(0.01, 0.05, 0.15, 0.25, seq(0.1, 1, by=0.1)),
# 												capture_rates=1,
# 												#capture_rate = 0.25,
# 												ladder=my_spike_df$conc_attomoles_ul_Mix1,
# 												total_ladder_transcripts=round(sum(ladder_molecules)),
# 											 	mRNA_degradation_rates=0.9, #c(0.0, 0.25, 0.5),
# 											 	reads_per_cell=2e5,
# 												#mRNA_degradation_rates=0,
# 												#mRNA_degradation_rates=0.9,
# 												#ladder_degradation_rate=0.15,
# 												ladder_degradation_rates=0.1,
# 												cells=10,
# 												run_rel2abs=TRUE, return_matrices=TRUE)
# qplot(rowMeans(spike_sim_rel2abs_by_reads$spike_rpc_matrix), rowMeans(spike_sim_rel2abs_by_reads$estimated_rpc_matrix), log="xy") + geom_abline()
# 

# selected_df <- subset(mode_df, num_genes==5000 & reads==1e6 & capture_rate == 0.1 & total_mRNA == 250000)
# qplot(ladder_degradation_rate, mse, data=selected_df, color=method) + 
# 	facet_grid(~capture_rate) + 
# 	scale_y_log10() + 
# 	geom_smooth()


# #####################################
# 
# # Pick around 5,000 transcripts to "express" from this hypothetical cell, consistent
# with what we typically see with current gen SC-RNA-Seq
reduced_fpkm_dist_5k <- original_fpkm_dist[sample(length(original_fpkm_dist), 5000)]
reduced_fpkm_dist_5k <- reduced_fpkm_dist_5k / sum(reduced_fpkm_dist_5k)

# # profiles <- generate_multinomial_profiles(reduced_fpkm_dist_5k)
# # pdf("./supplementary_figures/RPC_dist_by_yield.pdf", width=5, height=2)
# # qplot(expression, geom="histogram", data=profiles, log="x") + 
# # 	facet_wrap(~mRNA_yield, ncol=4) + 
# # 	geom_vline(aes(xintercept=mode_molecules))+
# # 	scale_x_log10(breaks=c(1, 10, 100, 1000)) +
# # 	 theme(axis.text.y=element_text(size=6)) +
# # 	 theme(axis.text.x=element_text(size=6)) +
# # 	 theme(axis.title.y=element_text(size=6)) +
# # 	 theme(axis.title.x=element_text(size=6)) +
# # 	 xlab("Expression (RPC)") +
# # 	 ylab("Transcript isoforms") +
# # 	 theme(panel.border = element_blank(), axis.line = element_line()) +
# # 	 theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
# # 	 theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
# #  	 theme(panel.background=element_blank()) +
# # 	 theme(legend.position = "none") +
# # 	 theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
# #  	 theme(strip.text = element_text(size=6))
# # dev.off()
# # 
# # qplot(tpm, geom="histogram", data=profiles, log="x") + 
# # 	facet_wrap(~mRNA_yield, ncol=4) + 
# # 	geom_vline(aes(xintercept=mode_log_tpm))+
# # 	scale_x_log10(breaks=c(1, 10, 100, 1000)) +
# # 	 theme(axis.text.y=element_text(size=6)) +
# # 	 theme(axis.text.x=element_text(size=6)) +
# # 	 theme(axis.title.y=element_text(size=6)) +
# # 	 theme(axis.title.x=element_text(size=6)) +
# # 	 xlab("Expression (RPC)") +
# # 	 ylab("Transcript isoforms") +
# # 	 theme(panel.border = element_blank(), axis.line = element_line()) +
# # 	 theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
# # 	 theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
# #  	 theme(panel.background=element_blank()) +
# # 	 theme(legend.position = "none") +
# # 	 theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
# #  	 theme(strip.text = element_text(size=6))
# # 
# # pdf("./supplementary_figures/RPC_dist_sampled.pdf", width=2, height=2)
# # qplot(reduced_fpkm_dist_5k, geom="histogram", log="x") +
# # theme(axis.text.y=element_text(size=6)) +
# # 	 theme(axis.text.x=element_text(size=6)) +
# # 	 theme(axis.title.y=element_text(size=6)) +
# # 	 theme(axis.title.x=element_text(size=6)) +
# # 	 xlab("Relative abundance") +
# # 	 ylab("Genes") +
# # 	 theme(panel.border = element_blank(), axis.line = element_line()) +
# # 	 theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
# # 	 theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
# #  	 theme(panel.background=element_blank()) +
# # 	 theme(legend.position = "none") +
# # 	 theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
# #  	 theme(strip.text = element_text(size=6))
# # dev.off()
# # 
# # pdf("./supplementary_figures/RPC_dist_orig.pdf", width=2, height=2)
# # qplot(original_fpkm_dist, geom="histogram", log="x") +
# # theme(axis.text.y=element_text(size=6)) +
# # 	 theme(axis.text.x=element_text(size=6)) +
# # 	 theme(axis.title.y=element_text(size=6)) +
# # 	 theme(axis.title.x=element_text(size=6)) +
# # 	 xlab("Relative abundance") +
# # 	 ylab("Genes") +
# # 	 theme(panel.border = element_blank(), axis.line = element_line()) +
# # 	 theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
# # 	 theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
# #  	 theme(panel.background=element_blank()) +
# # 	 theme(legend.position = "none") +
# # 	 theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
# #  	 theme(strip.text = element_text(size=6))
# # dev.off()
# # 
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
# 

#The following code will generate figures for supplementary figures 1 / 2: 

#supplementary figure 1: 
# mode_log_tpm_estimate=spike_rpc_df$mode_log_tpm_estimate,
# molecules_for_mode_log_tpm=estimated_rpc_all$calibrated_modes_df[, 1],
#supplementary figure 2: 
##### These figures argue explain why the k and b regression parameters lie on a line:
selected_df <- subset(mode_df, reads==1e6 & num_genes == 5000 & total_mRNA == 100000 & capture_rate %in% c(0.25, seq(0.3, 1, by = 0.1)))
capture_rate_colors <- brewer.pal(10, "PiYG")[-1]
pdf("./supplementary_figures/fig2g_simulation.pdf", width=2.15, height=1.5)
qplot(mode_log_tpm_estimate, molecules_for_mode_log_tpm, size=I(0.5), color=as.factor(capture_rate*100), data=selected_df) + 
	geom_smooth(method="lm", size=I(0.25), se=F) +
    theme(legend.key.size = unit(0.15, "in")) +
    scale_color_manual(values=capture_rate_colors) +
	xlab("Mode of TPM") +
	ylab("Mode of transcript counts") +
	guides(color=guide_legend(title="Capture\n(percent)")) +
	theme(legend.key=element_blank())+
	monocle:::monocle_theme_opts()
dev.off()


selected_df <- subset(mode_df, reads==1e6 & num_genes==5000 & capture_rate == 0.25 & total_mRNA %in% c(50000, 100000, 150000, 500000, 1000000))
#capture_rate_colors <- brewer.pal(length(unique(selected_df$capture_rate)) + 3, "Greens")[-c(1:3)]
total_mRNA_colors <- brewer.pal(length(unique(selected_df$total_mRNA)) + 3, "Reds")[-c(1:3)]
pdf("./supplementary_figures/fig2i_total_mRNAs_simulation.pdf", width=2.15, height=1.5)
qplot(mode_log_tpm_estimate, molecules_for_mode_log_tpm, size=I(0.5), color=as.factor(total_mRNA), data=selected_df) + 
	geom_smooth(method="lm", size=I(0.25), se=F) +
    theme(legend.key.size = unit(0.15, "in")) +
    scale_color_manual(values=total_mRNA_colors) +
  xlab("Mode of TPM") +
  ylab("Mode of transcript counts") +
  #theme(legend.position="none") +
	guides(color=guide_legend(title="Total\ntranscripts")) +
	theme(legend.key=element_blank())+
	monocle:::monocle_theme_opts()
dev.off()

selected_df <- subset(mode_df, reads==1e6  & total_mRNA == 100000 & capture_rate > 0.2)
num_genes_colors <- brewer.pal(length(unique(selected_df$num_genes)) + 3, "PiYG")[-c(1:3)]
pdf("./supplementary_figures/fig2g_simulation_num_genes.pdf", width=2.15, height=1.5)
qplot(mode_log_tpm_estimate, molecules_for_mode_log_tpm, size=I(0.5), color=as.factor(num_genes), data=subset(selected_df, num_genes > 2500)) + 
  geom_smooth(method="lm", size=I(0.25), se=F) +
  theme(legend.key.size = unit(0.15, "in")) +
  scale_color_manual(values=num_genes_colors) +
  xlab("Mode of TPM") +
  ylab("Mode of transcript counts") +
  guides(color=guide_legend(title="Number of \n expressed genes")) +
  theme(legend.key=element_blank())+
  monocle:::monocle_theme_opts()
dev.off()

selected_df <- subset(mode_df, num_genes==5000 & total_mRNA == 100000 & capture_rate %in% c(0.25, seq(0.3, 1, by = 0.1)))
capture_rate_colors <- brewer.pal(9, "Greens")
pdf("./supplementary_figures/fig2h_slope_simulation.pdf", width=2.15, height=1.35)
qplot(mode_log_tpm_estimate, molecules_for_mode_log_tpm, size=I(0.5), color=as.numeric(ladder_low_end_tpm/ladder_high_end_tpm), data=selected_df) + #
  geom_smooth(aes(group=as.factor(capture_rate)), color = 'black', method="lm", size=I(0.25), se=F) +
  theme(legend.key.size = unit(0.15, "in")) +
  scale_color_distiller(palette="BuPu") +
  xlab("Mode of TPM") +
  ylab("Mode of transcript counts") +
  #theme(legend.position="none") +
  guides(color=guide_legend(title="Spike TPM\nlow / high")) +
  theme(legend.key=element_blank())+
  monocle:::monocle_theme_opts()
dev.off()


# #simulate lung experiment: 
# load('/Users/xqiu/Dropbox (Personal)/Projects/BEAM/lung_endo_total')
# 
# simulate_lung_mode_df <- simulate_sequencing_with_ladder(original_fpkm_dist,
# 										   genes_expressed=pData(detectGenes(standard_cds))$num_genes_expressed,
# 										   total_mRNAs=lung_endo_total,
# 										   capture_rates=c(0.2326, 0.2987, 0.3093, 0.09174),
# 										   reads_per_cell=colSums(read_countdata),
# 										   ladder=spike_df$conc_attomoles_ul_Mix1,
# 										   total_ladder_transcripts=round(sum(spike_df$numMolecules)),
# 										   cells=10,
# 										   run_rel2abs=F)
# 
# #find the index for each capture rate: 
# expected_capture_rate <- c(rep(0.2326, sum(Time == 'E18.5')), rep(0.2987, sum(Time == 'E14.5')), rep(0.3093, sum(Time == 'Adult')), rep(0.09174, sum(Time == 'E16.5')))
# index <- rep(0, 183)
# for(i in 1:183) {
# 	index[i] <- which(simulate_lung_mode_df$capture_rate == expected_capture_rate[i] & 
# 		simulate_lung_mode_df$total_mRNA == lung_endo_total[i])
# }
# 
# qplot(regression_k, regression_b, data=simulate_lung_mode_df[index, ], color = total_mRNA) + 
# geom_smooth(method = 'rlm') + #xlim(0.6, 1.1) + ylim(-2, 0.5) + 
# theme(axis.text.y=element_text(size=6)) +
# 	 theme(axis.text.x=element_text(size=6)) +
# 	 theme(axis.title.y=element_text(size=6)) +
# 	 theme(axis.title.x=element_text(size=6)) +
# 	 xlab("Relative abundance") +
# 	 ylab("Transcript isoforms") +
# 	 theme(panel.border = element_blank(), axis.line = element_line()) +
# 	 theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
# 	 theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
#  	 theme(panel.background=element_blank()) +
# 	 theme(legend.position = "none") +
# 	 theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
#  	 theme(strip.text = element_text(size=6)) + scale_fill_brewer(palette="Set1")
# 
# subset_simulate_lung_mode_df <- simulate_lung_mode_df[index, ]
# qplot(mode_molecules_in_cell, data=subset_simulate_lung_mode_df, color = Time, geom = 'density', log = 'x') + geom_wrap(~Time)
# 
# 
# qplot(mode_molecules_in_cell, data=simulate_lung_mode_df[index, ], color = Time, geom = 'density', log = 'x') + 
# geom_smooth(method = 'rlm') + #xlim(0.6, 1.1) + ylim(-2, 0.5) + 
# theme(axis.text.y=element_text(size=6)) +
# 	 theme(axis.text.x=element_text(size=6)) +
# 	 theme(axis.title.y=element_text(size=6)) +
# 	 theme(axis.title.x=element_text(size=6)) +
# 	 xlab("Relative abundance") +
# 	 ylab("Transcript isoforms") +
# 	 theme(panel.border = element_blank(), axis.line = element_line()) +
# 	 theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
# 	 theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
#  	 theme(panel.background=element_blank()) +
# 	 theme(legend.position = "none") +
# 	 theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
#  	 theme(strip.text = element_text(size=6)) + scale_fill_brewer(palette="Set1")

save.image('./RData/mode_simulation.RData')

