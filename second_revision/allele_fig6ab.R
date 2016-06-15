so_switch_test_res <- allelic_bias_test(sandberg[row.names(subset(fData(sandberg), gene_id %in% sig_genes)),], fullModelFormulaStr="~stage", cores=detectCores() / 2)
iso_switch_test_res <- merge(distinct(fData(sandberg)[,c("gene_id", "symbol", "chr")]), iso_switch_test_res, by.x="gene_id", by.y="row.names")
row.names(iso_switch_test_res) <- iso_switch_test_res$gene_id
iso_switch_test_res$qval <- 1.0
iso_switch_test_res$qval[iso_switch_test_res$status == "OK"] <- p.adjust(subset(iso_switch_test_res, status == "OK")$pval, method="fdr")

allele_switch_genes <- subset(iso_switch_test_res, qval < 0.01 & chr == "X")
allele_switch_records <- row.names(subset(fData(sandberg), symbol %in% allele_switch_genes$symbol))
plot_genes_jitter(sandberg[allele_switch_records[1:16],], 
                  grouping="stage", color_by="allele", plot_trend=TRUE) + 
  facet_wrap( ~ feature_label, scales="free_y", ncol=1)

allele_switch_records <- row.names(subset(fData(sandberg), symbol %in% c("Nanog")))

plot_genes_jitter(sandberg[allele_switch_records,], 
                  grouping="stage", color_by="allele", plot_trend=TRUE, min_expr=0.1) + 
  facet_wrap( ~ feature_label, scales="free_y", ncol=1)

allele_switch_records <- row.names(subset(fData(sandberg), symbol %in% c("Xist", "Tsix", "Jpx", "Firre")))
#iso_switch_test_res <- allelic_bias_test(sandberg[allele_switch_records,], fullModelFormulaStr="~stage", cores=2)
plot_genes_jitter(sandberg[allele_switch_records,is.na(pData(sandberg)$sex) == FALSE & pData(sandberg)$sex == "F"], 
                  grouping="stage", color_by="allele", plot_trend=TRUE, min_expr=0.1) + 
  facet_wrap( ~ feature_label, scales="free_y", ncol=1)

allele_switch_records <- row.names(subset(fData(sandberg), symbol %in% c("Gm44279", "Mrpl33", "Cbx7", "Dock5")))
plot_genes_jitter(sandberg[allele_switch_records,], 
                  grouping="stage", color_by="allele", plot_trend=TRUE, min_expr=0.1) + 
  facet_wrap( ~ feature_label, scales="free_y", ncol=1)

#########??????? no sandberg_female_X ####
## X_inact_allele_matrix <- calculate_allele_ratio_matrix(sandberg_female_X, 
##                                                        subset(fData(sandberg_female_X), num_cells_expressed >= 1)$symbol,
##                                                         cores=detectCores() / 2)

allele_matrix <- calculate_allele_ratio_matrix(sandberg, 
                                               c("Xist", "Tsix", "Jpx", "Firre"),
                                               cores=detectCores() / 2)

allele_switch_genes <- subset(iso_switch_test_res, qval < 1e-4)
allele_matrix <- calculate_allele_ratio_matrix(sandberg, 
                                               allele_switch_genes$symbol,
                                               cores=detectCores() / 2)

num_allele_switch_clusters <- 8

mean_exprs <- rowMeans(exprs(sandberg[row.names(subset(fData(sandberg), gene %in% row.names(allele_switch_genes))), ]))
valid_genes <- unique(as.character(fData(sandberg[names(mean_exprs)[mean_exprs > 1], ])$gene))

pdf ("fig6a.pdf", height=4, width=3)
#iso_switch_heatmap <- 
plot_allele_ratio_heatmap(sandberg, 
                          allele_matrix[intersect(valid_genes, row.names(allele_matrix)) , ],
                          num_clusters=num_allele_switch_clusters)
dev.off()

Deng_sexed_cells <- sandberg[,pData(sandberg)$stage %in% c("4-cell", "16-cell", "early-blast")]
Deng_sexed_cells <- detectGenes(Deng_sexed_cells, min_expr=1)

Deng_Y_markers <- row.names(subset(fData(sandberg), symbol %in% c("Eif2s3y", "Zfy2", "Uba1y", "Kdm5d")))
all_Y_genes <- row.names(subset(fData(sandberg), chr == "Y" & num_cells_expressed >= 5))

Y_marker_exprs <- limma:::removeBatchEffect(vstExprs(Deng_sexed_cells[all_Y_genes,], round_vals=F), batch=pData(Deng_sexed_cells)$stage)
sex_cluster <- hclust(dist(t(Y_marker_exprs)), method="ward.D2")
plot(sex_cluster)
sex_cluster_cut <- cutree(sex_cluster, k=2)
#Y_marker_exprs <- colMeans(round(Y_marker_exprs))
sex_cluster_cut[sex_cluster_cut == 1] <- "Male"
sex_cluster_cut[sex_cluster_cut == 2] <- "Female"

pData(Deng_sexed_cells)$Sex_cluster <- sex_cluster_cut[row.names(pData(Deng_sexed_cells))]
sex_embryo_table <- table(pData(Deng_sexed_cells)$Sex_cluster, pData(Deng_sexed_cells)$embryo)
sex_embryo_table <- as.matrix(sex_embryo_table)
embryo_sexes <- sex_embryo_table["Female",] >  sex_embryo_table["Male",]
male_embryos <- names(embryo_sexes[embryo_sexes == FALSE])

Deng_sexed_cells$Deng_Sex <- factor(rep("Female", ncol(Deng_sexed_cells)), levels=c("Male", "Female"))
Deng_sexed_cells$Deng_Sex[pData(Deng_sexed_cells)$embryo %in% male_embryos] <- "Male"
#pData(sandberg)$Deng_Sex[Y_marker_exprs > 0] <- "Male"

gene_metadata <- read.delim("genes.attr_table")
gene_metadata <- subset(gene_metadata, grepl("ERCC", tracking_id) == FALSE)

locus_cols <- str_split_fixed(gene_metadata$locus, ":|-", 3) 
gene_metadata$gene_id <- str_split_fixed(gene_metadata$gene_id, "\\.", 2)[,1]

locus_df <- data.frame(gene_id=gene_metadata$gene_id, start=as.integer(locus_cols[,2]), stop=as.integer(locus_cols[,3]))


allele_switch_records <- row.names(subset(fData(Deng_sexed_cells), chr == "X" & is.allele.informative == 1))

X_inact_allele_matrix <- calculate_allele_ratio_matrix(Deng_sexed_cells[allele_switch_records, pData(Deng_sexed_cells)$Deng_Sex == "Female"], 
                                                       subset(fData(Deng_sexed_cells), num_cells_expressed >= 1)$symbol,
                                                       cores=detectCores() / 2)

X_inact_allele_matrix_melted <- melt(X_inact_allele_matrix)
colnames(X_inact_allele_matrix_melted) <- c("gene_id", "stage", "expression")
X_inact_allele_matrix_melted <- merge(X_inact_allele_matrix_melted, locus_df, by="gene_id")
X_inact_allele_matrix_melted <- subset(X_inact_allele_matrix_melted, stage %in% c("4-cell", "16-cell", "early-blast"))

X_inact_df <- X_inact_allele_matrix_melted %>%
  group_by(stage) %>%
  arrange(stage, start) %>%
  mutate(moving_avg_expr = rollmean(x = expression, 10, align = "right", fill = NA))

pdf("fig6b.pdf", width=1.35, height=1.5)
qplot(as.numeric(start) / 1e6, moving_avg_expr, data=X_inact_df, color=stage, geom="step", size=I(0.25)) + 
  xlab("X chr position (Mb)") + 
  ylab("Maternal expression") + 
  theme(legend.position = "none") + 
  geom_hline(yintercept=0.5, linetype=1) + 
  geom_vline(xintercept=103.414466, linetype=2) + 
  scale_color_brewer(palette="Set1") + 
  scale_y_continuous(limits=c(0,1), label=scales::percent) + monocle:::monocle_theme_opts()
dev.off()

############
X_inact_targets <- Deng_sexed_cells[row.names(subset(fData(Deng_sexed_cells), gene_short_name %in%   c("Xist", "Rlim", "Kif4", "Huwe1", "Atrx"))),
                                    pData(Deng_sexed_cells)$stage == "early-blast" &
                                      pData(Deng_sexed_cells)$Deng_Sex == "Female" ]
sim_escape(X_inact_targets, formulaStr="~1", allele_thresh=0.1, dispersion=rpc_dispersion, nsims=1000, cores=detectCores() / 2)



