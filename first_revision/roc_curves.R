recall<-function(alpha_val, dd)
{
    num_diff_exp <- sum(dd$diff_exp == "True", na.rm=T) 
    num_true_positives <- sum(dd$diff_exp == "True" & dd$q_value <= alpha_val & dd$correct_direction, na.rm=T)
    if (num_diff_exp > 0)
        return (num_true_positives / num_diff_exp)
    else
        return (1.0)
}
​
precision<-function(alpha_val, dd)
{
    num_positives <- sum(dd$diff_exp != "No call" & dd$q_value <= alpha_val, na.rm=T) 
    num_true_positives <- sum(dd$diff_exp == "True" & dd$q_value <= alpha_val & dd$correct_direction, na.rm=T)
    if (num_positives > 0)
        return (num_true_positives / num_positives)
    else
        return (1.0)
}
​
true_positives<-function(alpha_val, dd)
{
    num_true_positives <- sum(dd$diff_exp == "True" & dd$q_value <= alpha_val & dd$correct_direction, na.rm=T)
    return(num_true_positives)
}
​
false_positives<-function(alpha_val, dd)
{
    num_positives <- sum(dd$diff_exp != "No call" & dd$q_value <= alpha_val, na.rm=T) 
    num_true_positives <- sum(dd$diff_exp == "True" & dd$q_value <= alpha_val & dd$correct_direction, na.rm=T)
    return(num_positives - num_true_positives)
}
​
generate_precision_recall_df <-function(dd)
{
    df <- rbind(data.frame(precision = precision(alpha, dd),
                           recall = recall(alpha, dd),
						   true_positives = true_positives(alpha, dd),
				           false_positives = false_positives(alpha, dd)))
    return (df)
}
​
​
count_psig<-function(alpha_val, dd)
{
    sum(dd$diff_exp != "No call" & dd$p_value <= alpha_val & dd$correct_direction, na.rm=T)
}
​
count_qsig<-function(alpha_val, dd)
{
    sum(dd$diff_exp != "No call" & dd$q_value <= alpha_val & dd$correct_direction, na.rm=T)
}
​
reduce_df <- function(dd)
{
	dd[,c("read_depth", "num_C1_replicates", "num_C2_replicates", "tracking_id", "method", "p_value", "q_value", "diff_exp", "correct_direction")]
}
​
generate_roc_df <-function(dd)
{
​
	calls <- subset(dd, diff_exp != "No call")
	
	p_vals <- calls$p_value
	classification <- calls$diff_exp == "True"
	
	pred_p_value <- prediction(calls$p_value, classification)
	perf_tpr_fpr <- performance(pred_p_value, "tpr", "fpr")
	
	pred_q_value <- prediction(calls$q_value, classification)
	perf_fdr <- performance(pred_q_value, "tpr", "tpr")
	
	total_tps <- sum(calls$diff_exp == "True")
	total_fps <- nrow(calls) - total_tps
	
	tp <- perf_tpr_fpr@x.values * total_tps
	fp <- perf_tpr_fpr@x.values * total_fps
	
	td <- perf_fdr@x.values * total_tps
	fd <- perf_fdr@x.values * total_fps
	
	roc_df <- data.frame(
        fpr = perf_tpr_fpr@y.value,
        tpr = perf_tpr_fpr@x.value
    )
    
    return (roc_df)
}

plot_basic_roc<-function(exp_df)
{
    reduced_dd <- reduce_df(exp_df)
	roc_df <- ddply(reduced_dd, .(read_depth, num_C1_replicates, num_C2_replicates, method), generate_roc_df)
    
    
    roc_plot <- qplot(fpr, tpr, data=roc_df, color=method, geom="line") + 
        xlab("False positive rate") +
        ylab("True positive rate") +
        #opts(plot.margin = unit(c(pmrg,pmrg,pmrg,pmrg), "lines"))+
        ylim(c(0, 1.0)) + 
		xlim(c(0, 1.0))
    roc_plot <- roc_plot + scale_color_manual(values = color_map)
    print(roc_plot)
}

pred <- prediction(predictions, labels)
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf, col=rainbow(10))

generate_roc_df <-function(p_value, classification, type = 'fpr') {
	pred_p_value <- prediction(p_value, classification)
	perf_tpr_fpr <- performance(pred_p_value, "tpr", "fpr")
	
    fpr = perf_tpr_fpr@y.values

    tpr = perf_tpr_fpr@x.values
    
    data.frame(tpr = tpr, fpr = fpr)
}

generate_roc_df(HSMM_bulk_pval_df[select_genes, 2], true_data != 1)






