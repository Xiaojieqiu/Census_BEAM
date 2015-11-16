library(monocle)
library(xacHelper)

load_all_libraries()
library(sp)
library(igraph)
library(grid)

##load the data: 
load('analysis_lung_data.RData')
load('spikein_free_algorithm_benchmark.RData')
#color scheme: 
prog_cell_state = "#979797"
AT1_cell_state = "#F05662" 
AT2_cell_state = "#7990C8" 
AT1_Lineage = "#BD1C7C"
AT2_Lineage = "#337DB9" 

#########################################################################################################
#debug the functions (see runme_submission_test.R): 

#########################################################################################################
# generate the figures for submission paper: 
# generate the figure again with package from Develop branch: 

#########################################################################################################
# figure 1a: 
# panel a:
# submission_directory <- "/Users/xqiu/Dropbox (Cole Trapnell's Lab)/manuscript/Manuscript_submission/AI/"
submission_directory <- "./"

#marker from Quake paper AT1 (Pdpn, Ager) and AT2 (Sftpc, Sftpb) cells) with examples of cell cycle genes
markers <- c('Pdpn', 'Sftpb', 'Ccnb2', 'Cdk1') #Ccnb2, Cdk1
abs_AT12_cds_subset_all_gene@reducedDimS <- AT12_cds_subset_all_gene@reducedDimS
abs_AT12_cds_subset_all_gene@reducedDimA <- AT12_cds_subset_all_gene@reducedDimA
abs_AT12_cds_subset_all_gene@reducedDimK <- AT12_cds_subset_all_gene@reducedDimK
abs_AT12_cds_subset_all_gene@reducedDimW <- AT12_cds_subset_all_gene@reducedDimW
abs_AT12_cds_subset_all_gene@minSpanningTree <- AT12_cds_subset_all_gene@minSpanningTree

pdf('./main_figures/fig1a.pdf', height = 2, width = 2.5)
plot_spanning_tree(abs_AT12_cds_subset_all_gene, color_by="Time", show_backbone=T, backbone_color = 'black',
    markers=markers, show_cell_names = F, show_all_lineages = F, cell_link_size = 0.2) + 
        scale_size(range = c(0.1, 2.5)) + nm_theme()
dev.off()

plot_spanning_tree(abs_AT12_cds_subset_all_gene, color_by="State", show_backbone=T, backbone_color = 'black',
    markers=NULL, show_cell_names = F, show_all_lineages = F, cell_link_size = 0.2) + 
        scale_size(range = c(0.1, 2.5)) + nm_theme()

#figure 1b: 
# markers <- c('Soat1', 'S100g', 'Clic5', 'Muc1')
example_ids <- row.names(subset(fData(absolute_cds), gene_short_name %in% 
                              markers))
# E14.5: #7CAF42
# E16.5: #00BCC3
# E18.5: #A680B9
# Adult: #F3756C

# "E18.5" "E14.5" "E16.5" "Adult"
new_cds <- buildLineageBranchCellDataSet(abs_AT12_cds_subset_all_gene[1:10, ], lineage_labels = c('AT1', 'AT2'))

colour_cell <- rep(0, length(new_cds$Lineage))
names(colour_cell) <- as.character(new_cds$Time)
colour_cell[names(colour_cell) == 'E14.5'] <- "#7CAF42"
colour_cell[names(colour_cell) == 'E16.5'] <- "#00BCC3"
colour_cell[names(colour_cell) == 'E18.5'] <- "#A680B9"
colour_cell[names(colour_cell) == 'Adult'] <- "#F3756C"

colour <- rep(0, length(new_cds$Lineage))
names(colour) <- as.character(new_cds$Lineage)
colour[names(colour) == 'AT1'] <- AT1_Lineage
colour[names(colour) ==  'AT2'] <- AT2_Lineage

pdf('./main_figures/fig1b.pdf', height = 2, width = 3)
plot_genes_branched_pseudotime2(abs_AT12_cds_subset_all_gene[example_ids, ], cell_color_by = "Time", trajectory_color_by = "Lineage", fullModelFormulaStr = '~sm.ns(Pseudotime, df = 3)*Lineage', normalize = F, stretch = T, lineage_labels = c('AT1', 'AT2'), cell_size = 1, ncol = 2) + nm_theme()+ ylab('Transcript counts') + xlab('Pseudotime')
dev.off()

#########################################################################################################
#figure 2: 
#fig 2b: 
fig2_genes <- c("Pdpn", "Sftpb", 'Hprt', 'Pgk1')#, 'Ubc', 'Rpl5', 'Puf60', 'Nucb2') #Qiange housekeeping genes: Hprt, Pgk1

fig2_genes_ids <- row.names(subset(fData(absolute_cds), gene_short_name %in% 
                              fig2_genes))

colour_cell <- rep(0, length(new_cds$Lineage))
names(colour_cell) <- as.character(new_cds$State)
colour_cell[names(colour_cell) == '1'] <- prog_cell_state
colour_cell[names(colour_cell) == '2'] <- AT1_Lineage
colour_cell[names(colour_cell) == '3'] <- AT2_Lineage

colour <- rep(0, length(new_cds$Lineage))
names(colour) <- as.character(new_cds$Lineage)
colour[names(colour) == 'AT1'] <- AT1_Lineage
colour[names(colour) ==  'AT2'] <- AT2_Lineage

abs_house_keeping_marker_branchTest_res <- branchTest(abs_AT12_cds_subset_all_gene[fig2_genes_ids, ], cores = 1, relative_expr = F, weighted = T)
abs_house_keeping_marker_branchTest_res[fig2_genes_ids, 'pval']

pdf('./main_figures/fig2b.pdf', width = 2.6, height = 1.75)
plot_genes_branched_pseudotime2(abs_AT12_cds_subset_all_gene[fig2_genes_ids, ], color_by = "State", panel_order = fig2_genes, 
    trajectory_color_by = "Lineage", trend_formula = '~sm.ns(Pseudotime, df = 3)*Lineage', reducedModelFormulaStr = '~sm.ns(Pseudotime, df = 3)', 
     normalize = T, stretch = T, lineage_labels = c('AT1', 'AT2'), cell_size = 1, ncol = 2, add_pval = T) + nm_theme()+ ylab('Transcript counts') + xlab('Pseudotime')
dev.off()

#fig 2c: 
markers <- c('Rtkn2', 'Sdpr', 'Egfl6', 'Hc') #Rtkn2: AT1 early; Sdpr: AT1 late; Egfl6: AT2 early; Hc: AT2 late
example_ids <- row.names(subset(fData(absolute_cds), gene_short_name %in% markers))

rMax <- function(df) {apply(df, 1, function(x) if(all(is.na(x))) NA else max(abs(x), na.rm = T))} #calculate row max
abs_AT12_cds_subset_all_gene_ILRs_list$str_raw_div_df <- abs_AT12_cds_subset_all_gene_ILRs_list$str_branchA_expression_curve_matrix - abs_AT12_cds_subset_all_gene_ILRs_list$str_branchB_expression_curve_matrix
abs_AT12_cds_subset_all_gene_ILRs_list$str_norm_div_df <- abs_AT12_cds_subset_all_gene_ILRs_list$str_raw_div_df / rMax(abs_AT12_cds_subset_all_gene_ILRs_list$str_raw_div_df) #calculate normalized divergence
abs_AT12_cds_subset_all_gene_ILRs_list$log_str_raw_div_df <- log2((abs_AT12_cds_subset_all_gene_ILRs_list$str_branchA_expression_curve_matrix + .1)/(abs_AT12_cds_subset_all_gene_ILRs_list$str_branchB_expression_curve_matrix + .1))
abs_AT12_cds_subset_all_gene_ILRs_list$norm_str_logfc_df <- abs_AT12_cds_subset_all_gene_ILRs_list$str_logfc_df / rMax(abs_AT12_cds_subset_all_gene_ILRs_list$log_str_raw_div_df)

all_abs_bifurcation_time <- detectBifurcationPoint(abs_AT12_cds_subset_all_gene_ILRs_list$norm_str_logfc_df[, 27:100])

bif_time <- all_abs_bifurcation_time[markers]

# State 1: #F3756C
# State 2: #A2A738
# State 3: #29B67A

colour_cell <- rep(0, length(new_cds$Lineage))
names(colour_cell) <- as.character(new_cds$State)
colour_cell[names(colour_cell) == '1'] <- prog_cell_state
colour_cell[names(colour_cell) == '2'] <- AT1_cell_state
colour_cell[names(colour_cell) == '3'] <- AT2_cell_state

colour <- rep(0, length(new_cds$Lineage))
names(colour) <- as.character(new_cds$Lineage)
colour[names(colour) == 'AT1'] <- AT1_Lineage
colour[names(colour) ==  'AT2'] <- AT2_Lineage

pdf('./tmp/submission_fig2b_time.pdf', height = 2, width = 3)
plot_genes_branched_pseudotime2(abs_AT12_cds_subset_all_gene[example_ids, ], color_by = "State", panel_order = markers, trajectory_color_by = 'Lineage', 
    fullModelFormulaStr = '~sm.ns(Pseudotime, df = 3)*Lineage', reducedModelFormulaStr = '~sm.ns(Pseudotime, df = 3)', normalize = T, stretch = T,
    lineage_labels = c('AT1', 'AT2'), cell_size = 1, ncol = 2, bifurcation_time  = abs(bif_time)) + nm_theme()+ ylab('Transcript counts') + xlab('Pseudotime')
dev.off()

abs_str_logfc_df_list <- calILRs(cds = abs_AT12_cds_subset_all_gene[add_quake_gene_all_marker_ids, ], lineage_states = c(2, 3), stretch = T,cores = 1, 
    trend_formula = "~sm.ns(Pseudotime, df = 3) * Lineage", ILRs_limit = 3, 
    relative_expr = T, weighted = T, label_by_short_name = F, 
    useVST = T, round_exprs = FALSE, pseudocount = 0, output_type = "all", file = "str_logfc_df", return_all = T)

#make the plot for the comparing of the bifurcation timing: 
abs_bifurcation_time <- detectBifurcationPoint(abs_str_logfc_df_list$str_norm_div_df[, ], ILRs_threshold = 0.3)
names(abs_bifurcation_time) <- fData(abs_AT12_cds_subset_all_gene[row.names(abs_str_logfc_df_list$str_norm_div_df), ])$gene_short_name
abs_valid_bifurcation_time  <- abs_bifurcation_time[!is.na(abs_bifurcation_time)]
abs_valid_bifurcation_time <- abs_valid_bifurcation_time[unique(names(abs_valid_bifurcation_time))]

abs_bif_df <- data.frame(bifurcation_time_point = abs_bifurcation_time[c(AT1_early, AT1_late, AT2_early, AT2_late)], 
                                    type = c(rep("AT1 early", length(AT1_early)), 
                                             rep("AT1 late", length(AT1_late)), 
                                             rep("AT2 early", length(AT2_early)), 
                                             rep("AT2 late", length(AT2_late))))
timing_example_ids <- row.names(subset(fData(absolute_cds), gene_short_name %in% 
                              c(AT1_early, AT1_late, AT2_early, AT2_late)))

abs_marker_genes_branch_pval <- branchTest(abs_AT12_cds_subset_all_gene[add_quake_gene_all_marker_ids, ], cores = 1, relative_expr = T, weighted = F)

valid_timing_id <- timing_example_ids[abs_marker_genes_branch_pval[timing_example_ids, 'qval'] < 0.01]

abs_bif_df <- abs_bif_df[!is.na(abs_bif_df$bifurcation_time_point), ]
data <- subset(abs_bif_df, abs(bifurcation_time_point) > 27)
data <- data[as.character(fData(absolute_cds[valid_timing_id, ])$gene_short_name), ]
data <- subset(data, !is.na(type))

pdf('./main_figures/fig2c.pdf', width = 2.25, height = 1.25)
qplot(type, abs(bifurcation_time_point), color = type, geom = c('jitter', 'boxplot'), data = data, alpha = I(0.7)) + 
    xlab('') + ylab('bifurcation time point') + #geom_boxplot(stat = "identity", aes(ymin = `0%`, lower = `25%`, middle = `50%`, upper = `75%`, ymax = `100%`)) 
    nm_theme()
dev.off()

#figure 2d: 
add_quake_gene_all_marker_ids_branchTest_res <- weihgted_relative_abs_AT12_cds_subset_all_gene[add_quake_gene_all_marker_ids, ]
valid_add_quake_gene_all_marker_ids <- row.names(subset(add_quake_gene_all_marker_ids_branchTest_res, qval < 0.1))

jet.colors <- colorRampPalette(c("#4F64AD", "#F6F7FB", "#F2991F"))
bk <- seq(-3.1,3.1, by=0.1)
hmcols <- jet.colors(length(bk) - 1)

#Quake figure in the paper: 
valid_add_quake_gene_all_marker_ids_res_no_fit_log <- plot_genes_branched_heatmap(abs_AT12_cds_subset_all_gene[add_quake_gene_all_marker_ids, ], num_clusters=4, norm_method = "log", use_fitting_curves = F, scaling = F, hmcols = hmcols, file_name = paste(submission_directory, 'tmp/submission_fig2d.pdf', sep = ''))
valid_add_quake_gene_all_marker_ids_res_FIT_log <- plot_genes_branched_heatmap(abs_AT12_cds_subset_all_gene[add_quake_gene_all_marker_ids, ], num_clusters=4, norm_method = "log", use_fitting_curves = T, scaling = F, hmcols = hmcols, file_name = paste(submission_directory, 'tmp/submission_fig2d.1.pdf', sep = ''))

#add the AT1/2 early/late annotation: 
#add cell cycle genes
cell_cycle_markers <- c('Ccnb2', 'Cdk1', 'Ccna2', 'Ccna1', 'Ccne1', 'Ccne2')
cell_cycle_markers_id <- row.names(subset(fData(absolute_cds), gene_short_name %in% cell_cycle_markers))

cell_cycle_timing_example_ids_df <- weihgted_relative_abs_AT12_cds_subset_all_gene[c(timing_example_ids, cell_cycle_markers_id), c('gene_short_name', 'pval')]

add_annotation_row <- data.frame(significance = as.numeric(cell_cycle_timing_example_ids_df$pval < 0.05), row.names = row.names(cell_cycle_timing_example_ids_df), 
                                Type = c(rep("Markers", length(timing_example_ids)), rep("Cell cycle", length(cell_cycle_markers))))

time_annotated_heatmap <- plot_genes_branched_heatmap(abs_AT12_cds_subset_all_gene[c(timing_example_ids, cell_cycle_markers_id), ], 
    num_clusters=4, norm_method = "log", use_fitting_curves = F, scaling = F, hmcols = hmcols, use_gene_short_name = T,
    add_annotation_row = add_annotation_row, file_name = paste(submission_directory, 'main_figures/fig2d.pdf', sep = ''), show_rownames = T)

#########################################################################################################
#figure 3: 
#fig 3a: 
transcript_num <- 38919
test <- exprs(iso_absolute_cds)[1:38919, ] #filter out the spike-in trnascripts
#iso_absolute_cds: it will be exactly one 
mode <-apply(test, 2, function(x) mlv(round(x)[round(x) > .1], method = "mfv")$M) #calculate the mode of transcript counts

pdf('./main_figures/fig3a.pdf', width = 1.38, height = 1.25)
ggplot(data = data.frame(mode = mode), aes(x = mode)) + geom_bar(fill = I('red'), size = .5) + #geom_vline(x = 1, linetype = 'longdash', color = I('blue'), size = .2) + 
    xlab('Mode of transcript counts') + ylab('Cells') + scale_x_continuous(breaks = 1:10) + monocle_theme_opts() + nm_theme() #make figure 1.a
dev.off()

# mode_df <- data.frame(mode = estimate_t(absolute_cds), Time = pData(absolute_cds)$Time)
# ggplot(aes(x = mode), data = mode_df) + geom_bar(fill = I('red'), size = .5) + 
# geom_vline(x = 1, linetype = 'longdash', color = I('blue'), size = .2) + facet_wrap(~ Time, scales = 'free') + 
#     xlab('Mode of transcript counts') + ylab('Cells') + scale_x_continuous(breaks = 1:10) + monocle_theme_opts() + nm_theme() #make figure 1.a

#fig 3b: 
test <- mapply(function(cell_dmode, model) {
    predict(model, newdata = data.frame(log_fpkm = log10(cell_dmode)), type = 'response')
}, as.list(estimate_t(exprs(isoform_count_cds)[1:119469, ])), molModels_select)

df <- pData(absolute_cds)
df$mode_transcript <- 10^test
df$estimate_mode <- estimate_t(exprs(isoform_count_cds))

#make figure 3b
pdf('./main_figures/fig3b.pdf', width = 2.2, height = 1.4)
qplot(ceiling(mode_transcript), fill = I('red'), data = df)  + xlab('Transcript count for most frequent log10(FPKM)') + ylab('Cells') + nm_theme() #+ geom_vline(x = 1, linetype = 'longdash', color = I('blue'), size = .1)
dev.off()

#fig 3c: 
Time <- pData(absolute_cds)$Time
kb_df <- t(rbind.data.frame(lapply(molModels, function(x) c(b = coef(x)[1], k = coef(x)[2]))))
colnames(kb_df) <- c('b', 'k')

t <- -kb_df[, 'b'] / kb_df[, 'k']
pdf('./main_figures/fig3c.pdf', width = 2, height = 2)
qplot(k, b, data = as.data.frame(kb_df), color = Time) + scale_size(range = c(0.1, 2.5)) + nm_theme() 
dev.off()

pdf('./tmp/fig3c_helper.pdf', width = 2, height = 2)
qplot(k, b, data = as.data.frame(kb_df), size = t, color = Time) + scale_size(range = c(0.1, 2.5)) #+ nm_theme() 
dev.off()

#fig 4e: 

############################make the landscape heatmap: 
optimization_landscape_3d_trim <- lapply(optimization_landscape_3d, function(x) x[c('m', 'c', 'optim_res')])
optimization_matrix<- do.call(rbind.data.frame, optimization_landscape_3d_trim)

optimization_matrix_filt <- subset(optimization_matrix, is.nan(optim_res) == FALSE & is.finite(optim_res))
max_optim_score <- 3
optimization_matrix_filt$optim_res[optimization_matrix_filt$optim_res > max_optim_score] <- max_optim_score

spdf <- SpatialPointsDataFrame( data.frame( x = optimization_matrix_filt$m , y = optimization_matrix_filt$c ) , data = data.frame( z = optimization_matrix_filt$optim_res ) )

# Make an evenly spaced raster, the same extent as original data
e <- extent( spdf )

# Determine ratio between x and y dimensions
ratio <- ( e@xmax - e@xmin ) / ( e@ymax - e@ymin )

# Create template raster to sample to
r <- raster( nrows = 56 , ncols = floor( 56 * ratio ) , ext = extent(spdf) )
rf <- rasterize( spdf , r , field = "z" , fun = mean )

# We can then plot this using `geom_tile()` or `geom_raster()`
rdf <- data.frame( rasterToPoints( rf ) )

optimal_solution <- head(arrange(optimization_matrix_filt, optim_res), 1)
pdf('./main_figures/fig3e.pdf', width = 1.38, height = 1.25)
ggplot( NULL ) + geom_raster( data = rdf , aes( x , y , fill = log10(layer) ) ) + 
    annotate("text", x = -3.85, y = 3.1, label = "True (m,c)", color="magenta", size=2) + 
    annotate("point", x = -4.277778, y = 2.932929, color="magenta", size = 1) + 
    #annotate("text", x = -3.7, y = 3.2, label = "True (m,c)") + 
    annotate("text", x = -5.1, y = 2.7, label = "Algorithm (m,c)", color="red", size=2) + 
    annotate("point", x = optimal_solution$m, y = optimal_solution$c, color="red", size=1) + 
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
    annotate("text", x = -3.85, y = 3.1, label = "True (m,c)", color="magenta", size=2) + 
    annotate("point", x = -4.277778, y = 2.932929, color="magenta", size = 1) + 
    #annotate("text", x = -3.7, y = 3.2, label = "True (m,c)") + 
    annotate("text", x = -5.1, y = 2.7, label = "Algorithm (m,c)", color="red", size=2) + 
    annotate("point", x = optimal_solution$m, y = optimal_solution$c, color="red", size=1) + 
    scale_fill_gradientn(guide=guide_legend(title=expression(paste(log[10](F)))), colours=brewer.pal(name="YlGnBu", n=7)) +
    xlab("m") + ylab("c") +
    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + scale_size(range = c(0.1, 2)) + 
    theme(panel.background = element_rect(fill='white'))
dev.off()

#fig 3f: 
pdf('./main_figures/fig3f.pdf', width = 2, height = 1.7)
qplot(pData(absolute_cds)$endogenous_RNA[pData(absolute_cds)$endogenous_RNA > 1e3], 
      pData(mc_adj_cds)$endogenous_RNA[pData(absolute_cds)$endogenous_RNA > 1e3], log="xy", color=pData(absolute_cds)$Time[pData(absolute_cds)$endogenous_RNA > 1e3], size = I(1)) + 
     geom_smooth(method="lm", color="black", size = .1) + geom_abline(color="red") +  
    xlab("Total endogenous mRNA \n (spike-in)") +
    ylab("Total endogenous mRNA \n (spike-in free algorithm)") + #scale_size(range = c(0.25, 0.25)) + 
    scale_color_discrete(name = "Time points") + nm_theme()
dev.off()

#fig 3g:
Time <- pData(abs_AT12_cds_subset_all_gene)$Time
E14.5_cell <- rowMeans(exprs(abs_AT12_cds_subset_all_gene[, which(Time == 'E14.5')]))
E16.5_cell <- rowMeans(exprs(abs_AT12_cds_subset_all_gene[, which(Time == 'E16.5')]))
E18.5_cell <- rowMeans(exprs(abs_AT12_cds_subset_all_gene[, which(Time == 'E18.5')]))
Adult_cell <- rowMeans(exprs(abs_AT12_cds_subset_all_gene[, which(Time == 'Adult')]))

mc_E14.5_cell <- rowMeans(exprs(mc_adj_cds[, colnames(abs_AT12_cds_subset_all_gene)[which(Time == 'E14.5')]]))
mc_E16.5_cell <- rowMeans(exprs(mc_adj_cds[, colnames(abs_AT12_cds_subset_all_gene)[which(Time == 'E16.5')]]))
mc_E18.5_cell <- rowMeans(exprs(mc_adj_cds[, colnames(abs_AT12_cds_subset_all_gene)[which(Time == 'E18.5')]]))
mc_Adult_cell <- rowMeans(exprs(mc_adj_cds[, colnames(abs_AT12_cds_subset_all_gene)[which(Time == 'Adult')]]))

mc_abs_exprs_df <- data.frame(spikein = as.vector(c(E14.5_cell, E16.5_cell, E18.5_cell, Adult_cell)), 
    mc_algorithm = as.vector(c(mc_E14.5_cell, mc_E16.5_cell, mc_E18.5_cell, mc_Adult_cell)),
    cell = rep(c("E14.5_cell", "E16.5_cell", "E18.5_cell", "Adult_cell"), each = nrow(mc_adj_cds)))

# pdf('./main_figures/fig3g2.pdf', width = 3, height = 2)
# qplot(spikein + 1, mc_algorithm + 1, log = 'xy', 
#     color = cell, data = mc_abs_exprs_df, size  = 1.5) + facet_wrap(~cell, scales = 'free', ncol = 2) + #  geom_smooth(method = 'rlm', aes(group = 199), size = .1) + 
#     scale_size(range = c(1.5, 1)) +  geom_abline() + xlab('Transcript counts (Spike-in)') + scale_size(range = c(0.25, 0.25)) + 
#         ylab('Transcript counts (Recovery algorithm)') + nm_theme()
# dev.off()

pdf('./main_figures/fig3g.pdf', width = 3, height = 2)
ggplot(mc_abs_exprs_df) + aes(x=spikein + 1, y= mc_algorithm + 1) + scale_x_log10() + scale_y_log10() + facet_wrap(~cell, scales = 'free', ncol = 2) + 
xlab('Transcript counts (Spike-in)') + #scale_size(range = c(0.25, 0.25)) + 
        ylab('Transcript counts (Recovery algorithm)')  + 
  #stat_density2d(geom="tile", aes(fill=..density..^1, alpha=1), contour=FALSE) + 
  geom_point(size=0.5, aes(color = cell)) + geom_abline() + nm_theme()
  # stat_density2d(geom="tile", aes(fill=..density..^1, alpha=ifelse(..density..^1<0.4,0,1)), contour=FALSE) 
dev.off()
  
# # fig 3h: 
# # show only the spike-in / mc algorithm test: 
mc_spikein_df <- plot_pre_rec_f1(test_p_list = list(mode_size_norm_permutate_ratio_by_geometric_mean = new_abs_size_norm_monocle_p_ratio_by_geometric_mean,
                                      mc_mode_size_norm_permutate_ratio_by_geometric_mean = new_mc_size_norm_monocle_p_ratio_by_geometric_mean),
                   permutate_pval = list(mode_size_norm_permutate_ratio_by_geometric_mean = mode_size_norm_permutate_ratio_by_geometric_mean,
                                         mc_mode_size_norm_permutate_ratio_by_geometric_mean = mc_mode_size_norm_permutate_ratio_by_geometric_mean),
                   row.names(absolute_cds), #gene_list, overlap_genes, high_gene_list
                   return_df = T, #na.rm = T, 
                   p_thrsld = 0.01, #0.05
                   rownames = c('monocle (New size normalization)', 'monocle (New size normalization, Estimate transcript)'))
mc_spikein_df$data_type = c("Spikein transcripts", "estimated transcripts")

mc_spikein_df[, 'Type'] <- c('Monocle', 'Monocle') # geom_bar(stat = 'identity', position = 'dodge') 
colnames(mc_spikein_df)[1:3] <- c('Precision', 'Recall', 'F1 score')

pdf('./main_figures/fig_3h.pdf', width = 5, height = 1.5)
ggplot(aes(factor(Type), value,  fill = data_type), data = melt(mc_spikein_df)) + geom_bar(position = position_dodge(), stat = 'identity') + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank()) + nm_theme() + xlab('') + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

pdf('./tmp/fig_3h_helper.pdf', width = 3, height = 2)
ggplot(aes(factor(Type), value,  fill = data_type), data = melt(mc_spikein_df)) + geom_bar(position = position_dodge(), stat = 'identity') + #facet_wrap(~variable) + 
ggtitle(title) + scale_fill_discrete('Type') + xlab('Type') + ylab('') + facet_wrap(~variable, scales = 'free_x') +  theme(axis.text.x = element_text(angle = 30, hjust = .9)) + 
ggtitle('') + theme(strip.text.x = element_blank(), strip.text.y = element_blank()) + theme(strip.background = element_blank())
dev.off()

########################################################################################################
# figure 4: 
# panel a: 

# branch_gene_ABCs <- abs_AT12_cds_subset_all_gene_ABCs[row.names(weihgted_relative_abs_AT12_cds_subset_all_gene[weihgted_relative_abs_AT12_cds_subset_all_gene[, 'qval'] <= 0.01, ]), ]

# AT1_gene_num <- sum(branch_gene_ABCs[, 1] > 0, na.rm = T)
# AT2_gene_num <- sum(branch_gene_ABCs[, 1] < 0, na.rm = T)

# Lineage_gene_num <- c(AT1_gene_num, AT2_gene_num)
# names(Lineage_gene_num) <- c("AT1", "AT2")
# colour[names(colour) == 'AT1'] <- AT1_Lineage
# colour[names(colour) ==  'AT2'] <- AT2_Lineage

# pdf('./main_figure/submission_fig4a.pdf', width = 1, height = 1.2)
# qplot(c('AT1', 'AT2'), Lineage_gene_num, stat = 'identity', geom = 'bar', fill = c(AT1_Lineage, AT2_Lineage)) + 
#     xlab('Lineages') + ylab('Gene number') +  scale_fill_manual(name = "", values = c(AT1_Lineage, AT2_Lineage), labels = c("AT1", "AT2")) + 
#     nm_theme()
# dev.off()

# panel b: 
branch_gene_str_norm_div_df <- abs_AT12_cds_subset_all_gene_ILRs_list$str_logfc_df[row.names(weihgted_relative_abs_AT12_cds_subset_all_gene[weihgted_relative_abs_AT12_cds_subset_all_gene[, 'qval'] <= 0.05, ]), ]
quake_branch_genes <- row.names(subset(weihgted_relative_abs_AT12_cds_subset_all_gene, qval <= 0.05))

fData(abs_AT12_cds_subset_all_gene)$num_cell_expressed <- esApply(abs_AT12_cds_subset_all_gene[, ], 1, function(x) sum(round(x) > 0))
valid_expressed_genes <- row.names(subset(fData(abs_AT12_cds_subset_all_gene), num_cell_expressed > 5))

all_AT12_heatmap <- plot_genes_branched_heatmap(abs_AT12_cds_subset_all_gene[quake_branch_genes[quake_branch_genes %in% valid_expressed_genes], ], 
    stretch = T, file_name = paste(submission_directory, "./main_figures/fig4b.pdf", sep = ''))

clusters <- as.numeric(all_AT12_heatmap$annotation_row$Cluster); 
names(clusters) <- fData(abs_AT12_cds_subset_all_gene[row.names(all_AT12_heatmap$annotation_row), ])$gene_short_name

branchGenes_gsa_results_branch_genes <- collect_gsa_hyper_results(abs_AT12_cds_subset_all_gene[, ], mouse_go_gsc, clusters)

pdf(paste(submission_directory, "./main_figures/fig4b_go_enrichment.pdf", sep = ''), height=100, width=15)
plot_gsa_hyper_heatmap(abs_AT12_cds_subset_all_gene, branchGenes_gsa_results_branch_genes, significance = 1e-2)
dev.off()

#save the hyper_df into a xls table: 
save_hyper_df(branchGenes_gsa_results_branch_genes, './supplementary_data/lung_hyper_df.xls') 

#panel c: 
#find the mutual inhibition gene pairs: 
motif_TFs <- c("ELF1", "EWSR1", "FLI1", "RAR", "BACH1", "STAT5A", "PPARG", "STAT5B", "SOX2", "TCFE2A", "RXR", "ARNT", "GATA1", "NR1H2", "BATF", "VDR", "MAX", "TAL1", "FOS", "HAND1", "POU5F1", "TCF3", "MYC", "TLX1", "NR1H3", "CEBPA", "MAF", "NFE2", "RXRA", "SMAD4", "SMAD3", "SMAD2", "MAFK", "STAT1", "DDIT3", "STAT2", "HIF1A", "JUN", "NFIC",  "HSF1", "RFX1", "SPI1", "MEF2C", "MYOD1", "MEF2A", "CDX2", "FOXO1", "TBP", "RORA", "REST", "NR2E3", "FOXO3", "GATA1", "GATA2", "ATOH1", "HOXC9", "FLI1", "GATA3", "GATA4", "FOXF2", "RREB1", "YY1", "RELA", "RXRA", "GABPA", "HNF4G", "MECOM", "JUNB", "HNF4A", "JUN", "TFAP2A", "NFE2L2", "PRDM1", "TFAP2C", "SOX3", "SOX2", "ELK1", "SOX6", "SOX9", "SRF", "MEIS1", "NR2C2", "FOXH1", "SRY", "T", "FOXQ1", "HOXA5", "ELK4", "LHX3", "JUND", "NKX3-2", "HOXA9", "NKX3-1", "BHLHE40", "RUNX1", "DUX4", "NKX2-5", "RUNX2", "TCF3", "PLAG1", "SREBF1", "KLF5", "MAFF", "ESRRA", "ESRRB", "RFX5", "MAFB", "FOXA1", "NR4A2", "TEAD1", "EN1", "MAFK", "USF1", "FOXP1", "BRCA1", "USF2", "FOXP2", "SREBF2", "NRF1", "ETS1", "EBF1", "RFX1", "MZF1", "RFX2", "KLF1", "FOXI1", "TCF12", "KLF4", "E2F1", "HLF", "ZBTB33", "HNF1B", "E2F3", "ELF1", "FOSL2", "HNF1A", "E2F4", "FOXA2", "E2F6", "ELF5", "PPARG", "PAX6", "SPI1", "PAX5", "TP63", "PAX4", "NFKB1", "CTCF", "ZEB1", "PAX2", "FEV", "CRX", "FOS", "MAX", "GFI1B", "HSF1", "NOBOX", "SPIB", "SOX17", "NFIL3", "MYB", "FOSL1", "MYC", "NR2F1", "EGR1", "ZFP423", "AR", "EGR2", "ZFX", "ESR1", "TP53", "ESR2", "ZNF143", "HLTF", "MYCN", "FOXC1", "SPZ1", "NFYB", "EHF", "NFYA", "NR3C1", "TCF7L2", "TCFCP2L1", "STAT6", "STAT4", "REL", "POU2F2", "HINFP", "BCL6", "MYOG", "THAP1", "GFI1", "NFATC2", "FOXD1", "FOXD3", "CEBPA", "INSM1", "ZNF263", "ERG", "CEBPB", "FOXL1", "CREB1", "STAT1", "STAT3", "SP1", "SP2", "IRF1", "IRF2", "PBX1", "NR5A2", "NHLH1")
motif_TFs_add <- c("RARA", "E2A", "TRP63", "TRP53", "ZFP143", "TFCP2L1")

motif_Tfs_id <- row.names(subset(fData(absolute_cds[, ]), toupper(gene_short_name) %in% c(motif_TFs, motif_TFs_add)))

#genes are not shown in the gene_short_name: 
setdiff(motif_TFs, toupper(fData(absolute_cds)$gene_short_name))
# "RAR:RXR"      "TCFE2A"      "DUX4"     "TP63"     "TP53"     "ZNF143" "TCFCP2L1" "ZNF263"
#RARA E2A - ? TRP63 TRP53 ZFP143 TFCP2L1 ?
#update the gene names: 

TF_5k_enrichment_gsc <- loadGSCSafe("./data/Lung_JASPAR_5kb_hits_olap.gmt", encoding="latin1") 

valid_gene_id_cell_new <- row.names(absolute_cds[which(rowSums(exprs(standard_cds) >= 1) > 4), ])
valid_cell_genes_TF_enrichment_results_5k <- collect_gsa_hyper_results(abs_AT12_cds_subset_all_gene[valid_gene_id_cell_new, ], TF_5k_enrichment_gsc, clusters)
hyper_df <- plot_gsa_hyper_heatmap(abs_AT12_cds_subset_all_gene[valid_gene_id_cell_new, ], valid_cell_genes_TF_enrichment_results_5k, significance = 1e-2)

load('hyper_df')

hyper_df[, 'first'] <- str_split_fixed(hyper_df[, 'gene_set'], "::", 2)[,1]
hyper_df[, 'second'] <- str_split_fixed(hyper_df[, 'gene_set'], "::", 2)[,2]

valid_gene_id_20_cell <- row.names(absolute_cds[which(rowSums(exprs(standard_cds) >= 1) > 40), ])
valid_hyper_df <- subset(hyper_df, (toupper(first) %in% toupper(fData(abs_AT12_cds_subset_all_gene[valid_gene_id_20_cell, ])$gene_short_name)) | 
    (toupper(second) %in% toupper(fData(abs_AT12_cds_subset_all_gene[valid_gene_id_20_cell, ])$gene_short_name))
    ) 

valid_hyper_df$sig <- T
cluster_color <- c('#E52027' = '#E52027', '#357EB9' = '#357EB9', '#4BAE49' = '#4BAE49', '#974F9F' = '#974F9F', '#F47D20' = '#F47D20', '#F6EE38' = '#F6EE38')
valid_hyper_df$cluster_color <- cluster_color[as.numeric(valid_hyper_df$cluster_id)]

pdf('./main_figures/fig4c.pdf', width = 2, height = 6)
qplot(cluster_id, gene_set, fill=cluster_color, geom="tile", data=valid_hyper_df) + nm_theme() + scale_fill_manual(values=cluster_color)  + 
    theme(axis.text.x = element_text(angle = 0, hjust = 1)) + xlab('') + ylab('')
dev.off()

#Figure 4d: 
colour_cell <- rep(0, length(new_cds$Lineage))
names(colour_cell) <- as.character(new_cds$State)
colour_cell[names(colour_cell) == '1'] <- prog_cell_state
colour_cell[names(colour_cell) == '2'] <- AT1_cell_state
colour_cell[names(colour_cell) == '3'] <- AT2_cell_state

colour <- rep(0, length(new_cds$Lineage))
names(colour) <- as.character(pData(new_cds)$Lineage)
colour[names(colour) == 'AT1'] <- AT1_Lineage
colour[names(colour) ==  'AT2'] <- AT2_Lineage

gene_grn_list <- infer_branch_gene_grn(TF_enrichment_gsc = TF_5k_enrichment_gsc, file = 'gene_regulatory_net_up5k', p_thrsld = 0.01, branchTest_res = weihgted_relative_abs_AT12_cds_subset_all_gene)
branch_motif_Tfs <- gene_grn_list$branch_tfs[toupper(gene_grn_list$branch_tfs) %in% toupper(valid_hyper_df$first) | 
              toupper(gene_grn_list$branch_tfs) %in% toupper(valid_hyper_df$second)]
branch_motif_Tfs_id <- row.names(subset(fData(abs_AT12_cds_subset_all_gene), toupper(gene_short_name) %in% branch_motif_Tfs)) 

pdf('./main_figures/fig4d.pdf', width = 5, height = 3)
# pdf('fig4d.pdf', width = 5, height = 3)
plot_genes_branched_pseudotime2(abs_AT12_cds_subset_all_gene[branch_motif_Tfs_id, ], cell_color_by = "State", 
                trajectory_color_by = "Lineage", fullModelFormulaStr = '~sm.ns(Pseudotime, df = 3)*Lineage', normalize = F, stretch = T,
                lineage_labels = c('AT1', 'AT2'), cell_size = 1, ncol = 4, reducedModelFormulaStr = "~sm.ns(Pseudotime, df=3)", add_pval = T) + 
                ylab('Transcript counts') + nm_theme() + xlab('Pseudotime')
dev.off()

# #Tcf12 is not the branching genes anymore
# plot_genes_branched_pseudotime2(abs_AT12_cds_subset_all_gene['ENSMUSG00000032228.10', ], cell_color_by = "State", 
#   trajectory_color_by = "Lineage", fullModelFormulaStr = '~sm.ns(Pseudotime, df = 3)*Lineage', normalize = F, stretch = T,
#   lineage_labels = c('AT1', 'AT2'), cell_size = 1, ncol = 7, reducedModelFormulaStr = "~ sm.ns(Pseudotime, df=3)", add_pval = T) + 
#   nm_theme()+ ylab('Transcript counts') + nm_theme() + xlab('Pseudotime')

#########################################################################################################

save.image('gen_lung_figures.RData')
