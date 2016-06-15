utosome_sim_res_counts<- sim_monoallelic_expr(sandberg_count[autosomal_gene_records,], 
                                               formulaStr="~embryo + Spike_factor", 
                                               allele_thresh=0.1,
                                               dispersion=0,
                                               cores=detectCores() / 2)
autosome_sim_res_counts$Region <- "Autosomal"

rpc_expr_df <- autosome_sim_res %>% select(gene_id, mean_expression) %>% rename(mean_rpc = mean_expression)
autosome_sim_res_counts <- inner_join(autosome_sim_res_counts, rpc_expr_df)

autosome_sim_res_counts$expr_interval <- cut_interval(log(autosome_sim_res_counts$mean_rpc), 25)
expected_autosome_sim_res_counts <- subset(autosome_sim_res_counts, variable=="expected_monoallelic_calls" & is.na(value) == FALSE)
obs_autosome_sim_res_counts <- subset(autosome_sim_res_counts, variable=="obs_monoallelic_calls" & is.na(value) == FALSE)



exp_med_df <- expected_autosome_sim_res_counts %>% 
  group_by(expr_interval) %>% 
  mutate(interval_expression=median(mean_rpc),
         med_expr=median(value), 
         up = quantile(value, probs = 0.97, na.rm=TRUE),
         lo = quantile(value, probs = 0.05, na.rm=TRUE)) %>% 
  select(interval_expression, expr_interval, med_expr, up, lo) %>% distinct()


obs_med_df <- obs_autosome_sim_res_counts %>% 
  group_by(expr_interval) %>% 
  mutate(interval_expression=median(mean_rpc),
         med_expr=median(value), 
         up = quantile(value, probs = 0.95, na.rm=TRUE),
         lo = quantile(value, probs = 0.05, na.rm=TRUE)) %>% 
  select(interval_expression, expr_interval, med_expr, up, lo) %>% distinct()

#qplot(x=interval_expression, y=value, geom="boxplot", data=exp_med_df) + geom_line(aes(x=interval_expression, y=exp_med_df, group=1), color="red", data=obs_med_df) 
genes_interval_test_df <- inner_join(obs_autosome_sim_res_counts, exp_med_df)
genes_above_interval <- subset(genes_interval_test_df, value > up)

pdf("fig6f.pdf",  width=1.5, height=1)
ggplot(aes(x=interval_expression, y=med_expr, ymax=up, ymin=lo), data=exp_med_df) + 
  scale_x_continuous(trans="log2", limits=c(1,32)) + 
  scale_y_continuous(labels=scales::percent) + 
  geom_linerange() + 
  xlab("Transcript counts") + 
  ylab("Monoallelic genes\n(read counts)") + 
  geom_point(aes(x=interval_expression, y=value, group=1), data=genes_above_interval, color="red", position="jitter", size=0.01, alpha=0.05) + 
  geom_line(aes(x=interval_expression, y=med_expr, group=1), color="steelblue") +
  geom_line(aes(x=interval_expression, y=med_expr, group=1), color="red", data=obs_med_df) +
  monocle:::monocle_theme_opts()
dev.off()

