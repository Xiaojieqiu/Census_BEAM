cape_candidate_gene_ids <- unique(as.character(subset(fData(Deng_sexed_cells), chr == "X" & is.allele.informative == 1)$gene_id))

#x_escape_candidate_gene_ids <- intersect(x_escape_candidate_gene_ids, subset(gencode_gene_metadata, biotype %in% c("protein_coding", "lincRNA", "antisense"))$gene_id)
x_escape_candidate_paternal_allele_ids <- row.names(subset(fData(Deng_sexed_cells), gene_id %in% x_escape_candidate_gene_ids & allele == "P"))

#Deng_sexed_cells_paternal_alleles <- Deng_sexed_cells[x_escape_candidate_paternal_allele_ids,]
#Deng_sexed_cells_paternal_alleles <- detectGenes(Deng_sexed_cells_paternal_alleles, min_exp=0.1)
#x_escape_candidate_gene_ids <- unique(as.character(subset(fData(Deng_sexed_cells_paternal_alleles), num_cells_expressed >= 1)$gene_id))


sandberg_escape_cds <- Deng_sexed_cells[row.names(subset(fData(Deng_sexed_cells), gene_id %in% x_escape_candidate_gene_ids)),
                                        pData(Deng_sexed_cells)$stage == "early-blast" &
                                          pData(Deng_sexed_cells)$Deng_Sex == "Female" ]
sandberg_escape_sim_res_early_blast<- sim_escape(sandberg_escape_cds, formulaStr="~1", allele_thresh=0.1, dispersion=rpc_dispersion, nsims=1000, cores=detectCores() / 2)
sandberg_escape_sim_res_early_blast$stage= "early-blast"

#subset(fData(sandberg_escape_cds), gene_id %in% subset(sandberg_escape_sim_res, paternal_expr_qval < 0.05)$gene_id)

sandberg_escape_cds <- Deng_sexed_cells[row.names(subset(fData(Deng_sexed_cells), gene_id %in% x_escape_candidate_gene_ids)),
                                        pData(Deng_sexed_cells)$stage == "16-cell" &
                                          pData(Deng_sexed_cells)$Deng_Sex == "Female"]
sandberg_escape_sim_res_16_cell<- sim_escape(sandberg_escape_cds, formulaStr="~1", allele_thresh=0.1, dispersion=rpc_dispersion, nsims=1000, cores=detectCores() / 2)
sandberg_escape_sim_res_16_cell$stage= "16-cell"

sandberg_escape_cds <- Deng_sexed_cells[row.names(subset(fData(Deng_sexed_cells), gene_id %in% x_escape_candidate_gene_ids)),
                                        pData(Deng_sexed_cells)$stage == "4-cell" &
                                          pData(Deng_sexed_cells)$Deng_Sex == "Female"]
sandberg_escape_sim_res_4_cell<- sim_escape(sandberg_escape_cds, formulaStr="~1", allele_thresh=0.1, dispersion=rpc_dispersion, nsims=1000, cores=detectCores() / 2)
sandberg_escape_sim_res_4_cell$stage= "4-cell"

escape_sim_res <- rbind(sandberg_escape_sim_res_4_cell, sandberg_escape_sim_res_16_cell, sandberg_escape_sim_res_early_blast)

escape_sim_res$stage <- factor(escape_sim_res$stage, levels=c("4-cell", "16-cell", "early-blast"))
escape_sim_res$gene_id <- row.names(escape_sim_res)

allelic_bias_genes <- subset(escape_sim_res, allelic_bias_qval < 0.05)$gene_id
subset(fData(Deng_sexed_cells), gene_id %in% allelic_bias_genes)

escape_sim_res <- escape_sim_res[,c("gene_id", "stage", "maternal_expr_qval", "paternal_expr_qval")]
#Xist results:
subset(escape_sim_res, gene_id == "ENSMUSG00000086503")

#
#escape_sim_res <- subset(escape_sim_res, monoallele_qval )

escape_sim_res <- melt(escape_sim_res)

escape_df <- escape_sim_res %>% group_by(stage) %>%
  mutate(n_detected = sum(value < 0.01)) %>% 
  group_by(stage, variable) %>% 
  mutate(percent_allele = sum(value < 0.01) / n_detected) %>%
  select(stage, variable, percent_allele) %>%
  distinct()


#escape_sim_res <- merge(escape_sim_res, fData(Deng_sexed_cells)[,c("gene_id", "gene_short_name", "num_cells_expressed")], by.x="row.names", by.y="gene_id")
pdf("fig6c.pdf", width=1.35, height=1)
ggplot(aes(stage, percent_allele, fill=variable), data=escape_df) + 
  geom_bar(position="dodge", stat="identity") + scale_fill_manual(values=c("#C51B7D", "#4D9221")) + 
  ylab("Detectable alleles")+
  scale_y_continuous(label=scales::percent) +
  theme(axis.title.x=element_blank()) +
  theme(legend.position="none") + 
  #scale_fill_brewer(palette="Set1") + 
  monocle:::monocle_theme_opts()
dev.off()

early_blast_escape_genes <- subset(sandberg_escape_sim_res_early_blast, paternal_expr_qval < 0.05)$gene_id
subset(fData(Deng_sexed_cells), gene_id %in% early_blast_escape_genes)

known_escape_gene_names <- c("Xist", "Mid1", "Firre", "Kdm6a", "Pbdc1", "Eif2s3x", "Ddx3x", "Shroom4", "BC022960", "Kdm5c", "Car5b", "Bgn", "1810030O07Rik")
known_escape_genes <- unique(subset(fData(Deng_sexed_cells), symbol %in% c(known_escape_gene_names))$gene_id)
missing_escapees <- setdiff(known_escape_genes, early_blast_escape_genes)
missing_escapees <- subset(fData(Deng_sexed_cells), gene_id %in% missing_escapees)

early_blast_maternal <- subset(sandberg_escape_sim_res_early_blast, maternal_expr_qval < 0.01)$gene_id #maternal_expression????????
subset(fData(Deng_sexed_cells), gene_id %in% early_blast_maternal)

female_cds <- Deng_sexed_cells[, pData(Deng_sexed_cells)$Deng_Sex == "Female" ]

plot_genes_jitter(female_cds[row.names(subset(fData(Deng_sexed_cells), gene_short_name %in%  c("Xist", "Jpx", "Eif2s3x",  "Ddx3x", "Kdm6a", "Pbdc1",  "Kdm5c", "Car5b"))),], 
                  grouping="stage", color_by="allele", plot_trend=TRUE, min_expr=0.5) 

plot_genes_jitter(female_cds[row.names(subset(fData(Deng_sexed_cells), gene_short_name %in%  c("Pgk1", "Gprasp1", "Klf8", "Mid1", "Fam3a", "Cetn2", "Dlg3"))),], 
                  grouping="stage", color_by="allele", plot_trend=TRUE, min_expr=0.5) 

early_blast_inactivated_genes <- subset(sandberg_escape_sim_res_early_blast, paternal_expr_qval > 0.05 & maternal_expr_qval < 0.05)$gene_id
subset(fData(Deng_sexed_cells), gene_id %in% early_blast_inactivated_genes)

subset(sandberg_escape_sim_res_early_blast, gene_id %in% subset(fData(Deng_sexed_cells), 
                                                                gene_short_name %in%  c("Xist", "Rlim", "Kif4", "Huwe1", "Atrx"))$gene)

fData(female_cds[row.names(subset(fData(Deng_sexed_cells), 
                                  gene_short_name %in%  c("Xist", "Rlim", "Kif4", "Huwe1", "Atrx"))),])

X_inact_targets <- Deng_sexed_cells[row.names(subset(fData(Deng_sexed_cells), gene_short_name %in%  c("Xist", "Rlim", "Kif4", "Huwe1", "Atrx"))),
                                    pData(Deng_sexed_cells)$stage == "early-blast" &
                                      pData(Deng_sexed_cells)$Deng_Sex == "Female" ]
sim_escape(X_inact_targets, formulaStr="~num_genes_expressed", allele_thresh=0.1, dispersion=0, nsims=1000, cores=detectCores() / 2)


subset_female_cds <- female_cds[row.names(subset(fData(Deng_sexed_cells), 
                            gene_short_name %in%  c("Xist", "Rlim", "Kif4", "Huwe1", "Atrx"))),]

# subset_female_cds <- subset_female_cds[esApply(subset_female_cds, 1, mean) > 0.1, ]

pdf ("fig6d.pdf", width=1.5, height=4.25)
plot_genes_jitter(subset_female_cds,
                  grouping="stage", color_by="allele", plot_trend=FALSE, min_expr=0, cell_size= 0.15) + 
  stat_summary(aes_string(x="stage", y="expression", color="allele", group="allele"), fun.data = "mean_cl_boot", size=0.15, geom="line") +
  stat_summary(aes_string(group="allele"), color="black", fun.data = "mean_cl_boot", size=0.15, geom="point") +
  stat_summary(aes_string(color="allele"), fun.data = "mean_cl_boot", size=0.15, geom="point") + scale_size(range = c(0.1, 0.15)) + 
  scale_color_manual(values=c("#C51B7D", "#4D9221")) +
  # scale_y_log10() +  
  theme(legend.position="none") + 
  nm_theme()
dev.off()

