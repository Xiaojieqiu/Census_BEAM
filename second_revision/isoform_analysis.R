library(stringr)
library(scales)
# library(monocle)
library(devtools)
load_all('/Users/xqiu/Dropbox (Personal)/Projects/monocle-dev') 

library(reshape2)
library(plyr)
library(dplyr)
library(piano)
library(corrplot)
library(matrixStats)

#########################
theme_set(theme_bw())

# Set global ggplot2 properties for making print-scaled PDF panels
theme_set(theme_bw(base_size=6))

time_colors = c("#D7191C", "#FDAE61", "#ABD9E9", "#2C7BB6")

# Some auxiliarry code:
source ("support_functions.q")


#########################


####################

#### Load up the HSMMs

HSMM_fpkm_matrix <- read.delim("~/Dropbox (Cole Trapnell's Lab)/Shared Data/Muscle/HSMM/sc-RNA-Seq/HSMM_cuffnorm_out/genes.fpkm_table")
row.names(HSMM_fpkm_matrix) <- HSMM_fpkm_matrix$tracking_id
HSMM_fpkm_matrix <- HSMM_fpkm_matrix[,-1]

HSMM_isoform_fpkm_matrix <- read.delim("~/Dropbox (Cole Trapnell's Lab)/Shared Data/Muscle/HSMM/sc-RNA-Seq/HSMM_cuffnorm_out/isoforms.fpkm_table")
row.names(HSMM_isoform_fpkm_matrix) <- HSMM_isoform_fpkm_matrix$tracking_id
HSMM_isoform_fpkm_matrix <- HSMM_isoform_fpkm_matrix[,-1]

sample_sheet <- read.delim("~/Dropbox (Cole Trapnell's Lab)/Shared Data/Muscle/HSMM/sc-RNA-Seq/sample_sheet.txt")
sample_sheet$cell_id <- paste(sample_sheet$cell_id, "_0", sep="")
row.names(sample_sheet) <- sample_sheet$cell_id
sample_sheet <- sample_sheet[colnames(HSMM_fpkm_matrix),]

cell_is_valid_singleton <- row.names(subset(sample_sheet, Control == FALSE & Unusual.Shape == FALSE & Debris == FALSE & Clump == FALSE & Cells.in.Well == 1)) 
sample_sheet <- sample_sheet[cell_is_valid_singleton,]

HSMM_fpkm_matrix <- HSMM_fpkm_matrix[,row.names(sample_sheet)]
HSMM_isoform_fpkm_matrix <- HSMM_isoform_fpkm_matrix[,row.names(sample_sheet)]

gene_ann <- read.delim("~/Dropbox (Cole Trapnell's Lab)/Shared Data/Muscle/HSMM/sc-RNA-Seq/HSMM_cuffnorm_out/genes.attr_table")
isoform_ann <- read.delim("~/Dropbox (Cole Trapnell's Lab)/Shared Data/Muscle/HSMM/sc-RNA-Seq/HSMM_cuffnorm_out/isoforms.attr_table")

gencode_biotypes <- read.delim("~/Dropbox (Cole Trapnell's Lab)/Shared Data/Muscle/HSMM/sc-RNA-Seq/gencode_biotypes.txt")

gene_ann <- merge(gene_ann, gencode_biotypes, by = "gene_id")
row.names(gene_ann) <- gene_ann$tracking_id
gene_ann <- gene_ann[,c("gene_id", "gene_short_name", "biotype")]

isoform_ann <- merge(isoform_ann, gencode_biotypes, by = "gene_id")
row.names(isoform_ann) <- isoform_ann$tracking_id
isoform_ann <- isoform_ann[,c("gene_id", "gene_short_name", "biotype", "tss_id")]


gene_ann <- gene_ann[row.names(HSMM_fpkm_matrix),]
isoform_ann <- isoform_ann[row.names(HSMM_isoform_fpkm_matrix),]


pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_ann)
HSMM <- newCellDataSet(as.matrix(HSMM_fpkm_matrix), 
                                     phenoData = pd, 
                                     featureData = fd)
HSMM_fpkm_matrix_adj <- relative2abs(HSMM, estimate_t(HSMM_isoform_fpkm_matrix), cores=detectCores() / 2)

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_ann)
HSMM <-  newCellDataSet(as.matrix(HSMM_fpkm_matrix_adj), 
               phenoData = pd, 
               featureData = fd, 
               expressionFamily=negbinomial(), 
               lowerDetectionLimit=1)
pData(HSMM)$Total_mRNAs <- colSums(exprs(HSMM))


pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = isoform_ann)
HSMM_isoform <- newCellDataSet(as.matrix(HSMM_isoform_fpkm_matrix), 
                                     phenoData = pd, 
                                     featureData = fd)
HSMM_isoform <- detectGenes(HSMM_isoform, min_expr = 1)


HSMM_isoform_fpkm_matrix_adj <- relative2abs(HSMM_isoform, estimate_t(HSMM_isoform_fpkm_matrix), cores=detectCores() / 2)


pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = isoform_ann)
HSMM_isoform <- newCellDataSet(as.matrix(HSMM_isoform_fpkm_matrix_adj), 
						   phenoData = pd, 
						   featureData = fd, 
						   expressionFamily=negbinomial(), 
						   lowerDetectionLimit=1)

pData(HSMM_isoform)$Total_mRNAs <- colSums(exprs(HSMM_isoform))

HSMM <- detectGenes(HSMM, min_expr = 0.1)
HSMM_isoform <- detectGenes(HSMM_isoform, min_expr = 0.1)

PDGFRA_expr <- colSums(exprs(HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("PDGFRA"))),]))
SPHK1_expr <- colSums(exprs(HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("SPHK1"))),]))
ANPEP_expr <- colSums(exprs(HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("ANPEP"))),]))
MEF2C_expr <- colSums(exprs(HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("MEF2C"))),]))
MYF5_expr <- colSums(exprs(HSMM[row.names(subset(fData(HSMM), gene_short_name %in% c("MYF5"))),]))


pData(HSMM)$PDGFRA <- as.vector(PDGFRA_expr)
pData(HSMM)$SPHK1 <- as.vector(SPHK1_expr)
pData(HSMM)$ANPEP <- as.vector(ANPEP_expr)
pData(HSMM)$MEF2C <- as.vector(MEF2C_expr)
pData(HSMM)$MYF5 <- as.vector(MYF5_expr)


SPHK1_thresh <- 1
PDGFRA_thresh <- 1
ANPEP_thresh <- Inf #0.1
MEF2C_thresh <- 5 #0.1
MYF5_thresh <- 1 #0.1

#cell_is_hsmm <- (pData(HSMM)$MEF2C > MEF2C_thresh | pData(HSMM)$MYF5 > MYF5_thresh) 
cell_is_hsmm <- (pData(HSMM)$MEF2C > MEF2C_thresh | pData(HSMM)$MYF5 > MYF5_thresh) 

names(cell_is_hsmm) <- row.names(pData(HSMM))
pData(HSMM)$CellType <- "Myoblast"
pData(HSMM)$CellType[cell_is_hsmm == FALSE] <-  "Fibroblast"

pData(HSMM_isoform)$CellType <- "Myoblast"
pData(HSMM_isoform)$CellType[cell_is_hsmm == FALSE] <-  "Fibroblast"


HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) + 2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) - 2*sd(log10(pData(HSMM)$Total_mRNAs)))
HSMM_myo <- HSMM[,pData(HSMM)$CellType == "Myoblast" & pData(HSMM)$Total_mRNAs > lower_bound & pData(HSMM)$Total_mRNAs < upper_bound]								  
HSMM_myo_isoform <- HSMM_isoform[,row.names(pData(HSMM_myo))]                 

qplot(Total_mRNAs, data=pData(HSMM_myo), color=Time, geom="density")

#######

HSMM_expressed_genes <- row.names(subset(fData(HSMM_myo), 
									num_cells_expressed >= 15 & 
									biotype %in% c("protein_coding", "lincRNA")))

HSMM_myo <- estimateSizeFactors(HSMM_myo)

HSMM_myo <- estimateDispersions(HSMM_myo, cores=detectCores() / 2)

diff_test_res <- differentialGeneTest(HSMM_myo[HSMM_expressed_genes,], 
                                      fullModelFormulaStr="~Time", cores=detectCores() / 2)

high_disp_genes <- row.names(subset(dispersionTable(HSMM_myo), mean_expression >= 1 & dispersion_empirical >= 2.5 * dispersion_fit))
# sig_gene_names <- row.names(subset(diff_test_res, qval < 1e-5))

HSMM_myo <- setOrderingFilter(HSMM_myo, high_disp_genes)

HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, method="ICA", use_vst=TRUE, pseudo_expr=0)

HSMM_myo <- orderCells(HSMM_myo, reverse=F, num_paths=1)

marker <- c("MYH3", "CDK1", "MEF2C", "H19")
marker_id <- row.names(subset(fData(HSMM_myo), gene_short_name %in% marker))
plot_spanning_tree(HSMM_myo, marker=marker, color_by = 'Time', cell_size = 4)
plot_genes_in_pseudotime(HSMM_myo[marker_id, ], cell_size = 4, color_by = 'Time')


pData(HSMM_myo_isoform)$Pseudotime <- pData(HSMM_myo)$Pseudotime



############

HSMM_myo_isoform <- detectGenes(HSMM_myo_isoform)
HSMM_myo_isoform <- estimateSizeFactors(HSMM_myo_isoform)
expr_matrix <- as.data.frame(t(t(exprs(HSMM_myo_isoform) / sizeFactors(HSMM_myo_isoform))))
expr_matrix <- cbind(expr_matrix, fData(HSMM_myo_isoform)$gene_id)
colnames(expr_matrix)[length(colnames(expr_matrix))] <- "gene_id"
expr_matrix <- as.data.frame(expr_matrix)

expr_matrix$gene_id <- as.factor(expr_matrix$gene_id)
expr_matrix$isoform_id <- row.names(expr_matrix)

# target_gene <- "MYH3"

# target_genes <- row.names(subset(fData(HSMM_myo_isoform), gene_short_name %in% c("MYH3", "MYOG", "MEF2C", "TPM2")))
# #target_genes <- 1:1500

# iso_switch_test_res <- iso_switch_test(HSMM_myo_isoform[target_genes,], cores=2)
# plot_genes_in_pseudotime(HSMM_myo_isoform[row.names(subset(fData(HSMM_myo_isoform), gene_short_name %in% c("TPM2"))),], color_by="Time", label_by_short_name=FALSE)
# plot_genes_in_pseudotime(HSMM_myo_isoform[row.names(subset(fData(HSMM_myo_isoform), gene_short_name %in% c("MEF2C"))),], color_by="Time", label_by_short_name=FALSE)

# plot_genes_in_pseudotime(HSMM_myo_isoform[row.names(subset(fData(HSMM_myo_isoform), gene_short_name %in% c("ACTA1"))),], color_by="Time", label_by_short_name=FALSE)


# plot_genes_in_pseudotime(HSMM_myo_isoform[row.names(subset(fData(HSMM_myo_isoform), 
#                          gene_short_name %in% c("TPM1") & num_cells_expressed >= 15)),], 
#                          color_by="Time", 
#                          label_by_short_name=FALSE,
#                          min_expr=0.1)

# splicing_factors <- HSMM_myo[row.names(subset(fData(HSMM_myo), 
#                              gene_short_name %in% c("PTBP1", "PTBP2", "MBNL1", "MBNL2", "MBNL3", "RAVER1", "SRSF1", "SRSF2"))),]
# splicing_factor_diff_test <- differentialGeneTest(splicing_factors)
# plot_genes_in_pseudotime(splicing_factors, 
#                          color_by="Time", 
#                          label_by_short_name=TRUE,
#                          min_expr=0.5)

# plot_genes_in_pseudotime(HSMM_myo_isoform[row.names(subset(fData(HSMM_myo_isoform), gene_short_name %in% c("PTBP1") & num_cells_expressed >=15 )),], color_by="Time", label_by_short_name=FALSE, min_expr=0.5)

# muscle_ids <- row.names(subset(fData(HSMM_myo), gene_short_name %in% c("MYOG", "MYOD1", "MYH3")))
# subset_vst <- vstExprs(HSMM_myo[muscle_ids,], round_vals=FALSE)
# vst_subset <- vstExprs(HSMM_myo, round_vals=FALSE)[muscle_ids,]


############


diff_gene_test_res <- differentialGeneTest(HSMM_myo[HSMM_expressed_genes,], cores=detectCores() / 2)

target_genes <- row.names(subset(diff_gene_test_res, qval < 0.05))
target_isoforms <- row.names(subset(fData(HSMM_myo_isoform), gene_id %in% target_genes))

iso_switch_test_res <- iso_switch_test(HSMM_myo_isoform[target_isoforms,], cores=detectCores() / 2)
iso_switch_test_res <- merge(fData(HSMM_myo), iso_switch_test_res, by="row.names")
iso_switch_test_res$qval <- 1.0
iso_switch_test_res$qval[iso_switch_test_res$status == "OK"] <- p.adjust(subset(iso_switch_test_res, status == "OK")$pval, method="fdr")

# Load some gene sets:
# Gene Ontology
#go_bp_gsc <- loadGSCSafe("~/Dropbox (Cole Trapnell's Lab)/Shared Data/GMT/v4/c5.bp.v4.0.symbols.gmt")
#go_cc_gsc <- loadGSCSafe("~/Dropbox (Cole Trapnell's Lab)/Shared Data/GMT/v4/c5.cc.v4.0.symbols.gmt")
#go_mf_gsc <- loadGSCSafe("~/Dropbox (Cole Trapnell's Lab)/Shared Data/GMT/v4/c5.mf.v4.0.symbols.gmt")

go_mf_gsc <- loadGSCSafe("~/Dropbox (Cole Trapnell's Lab)/Shared Data/GMT/EM_pathways/Human/symbol/GO/Human_GO_mf_with_GO_iea_symbol.gmt")
names(go_mf_gsc$gsc) <- str_split_fixed(names(go_mf_gsc$gsc), "%", 2)[,1]

go_cc_gsc <- loadGSCSafe("~/Dropbox (Cole Trapnell's Lab)/Shared Data/GMT/EM_pathways/Human/symbol/GO/Human_GO_cc_with_GO_iea_symbol.gmt")
names(go_cc_gsc$gsc) <- str_split_fixed(names(go_cc_gsc$gsc), "%", 2)[,1]

go_bp_gsc <- loadGSCSafe("~/Dropbox (Cole Trapnell's Lab)/Shared Data/GMT/EM_pathways/Human/symbol/GO/Human_GO_bp_with_GO_iea_symbol.gmt")
names(go_bp_gsc$gsc) <- str_split_fixed(names(go_bp_gsc$gsc), "%", 2)[,1]


# Reactome pathways
human_reactome_gsc <- loadGSCSafe("~/Dropbox (Cole Trapnell's Lab)/Shared Data/GMT/EM_pathways/Human/symbol/Pathways/Human_Reactome_June_20_2014_symbol.gmt", encoding="latin1")
names(human_reactome_gsc$gsc) <- str_split_fixed(names(human_reactome_gsc$gsc), "%", 2)[,1]


gsa_gene_ids = target_genes

#plot_spanning_tree(CD8_T_cells_ordered, color_by="Condition", cell_size=I(0.15), markers=c("APOC1"))
#plot_genes_branched_pseudotime(CD8_T_cells_ordered[c("APOC1", "MALAT1"),], color_by="Condition", branch_point=1, min_expr=0.1, ncol=2, add_pval=TRUE)

iso_switching_genes <- subset(iso_switch_test_res, qval < 0.1)$Row.names
iso_switching_clusters <- rep(1, length(iso_switching_genes))
names(iso_switching_clusters) <- iso_switching_genes
iso_switch_reactome_gsa_hyper_heatmap <- collect_gsa_hyper_results(HSMM_myo[HSMM_expressed_genes,], human_reactome_gsc, iso_switching_clusters)
#pdf ("bp5_reactome_gsa_hyper_heatmap.pdf", height=48, width=16)
plot_gsa_hyper_heatmap(HSMM_myo,  iso_switch_reactome_gsa_hyper_heatmap, significance=0.1)
#dev.off()

write.table(sort(subset(iso_switch_test_res, qval < 0.1)$gene_short_name), "isoform_switching_gene_names.txt", quote=F, sep="\t", col.names=F, row.names=F)

HSMM_myo_isoform <- estimateDispersions(HSMM_myo_isoform)

# plot_isoform_covariance <- function(isoform_cds, gene_name){
#   gene_cds <- isoform_cds[row.names(subset(fData(isoform_cds), gene_short_name %in% gene_name & num_cells_expressed >= 15)),]
#   isoform_exprs <- vstExprs(gene_cds)

#   isoform_expr_cor <- cor(t(isoform_exprs))
#   #invalid_isos <- rowSum(is.na(isoform_expr_cor)) != 0
#   isoform_expr_cor[is.na(isoform_expr_cor)] = 0
#   diag(isoform_expr_cor) <- 0
#   #print (isoform_expr_cor)
#   #print (dim(isoform_expr_cor))
#   corrplot(isoform_expr_cor, "ellipse", order="hclust")

# }
# plot_isoform_covariance(HSMM_myo_isoform, "MYL6")
# plot_isoform_covariance(HSMM_myo_isoform, "TNNT2")
# plot_isoform_covariance(HSMM_myo_isoform, "TPM1")
# plot_isoform_covariance(HSMM_myo_isoform, "TPM2")
# plot_isoform_covariance(HSMM_myo_isoform, "FN1")
# plot_isoform_covariance(HSMM_myo_isoform, "IL32")
# plot_isoform_covariance(HSMM_myo_isoform, "ACTA1")


plot_isoform_expression_heatmap <- function(isoform_cds, gene_names, num_clusters=6, add_annotation_row=NULL){
  gene_cds <- isoform_cds[row.names(subset(fData(isoform_cds), gene_short_name %in% gene_names & num_cells_expressed >= 15)),]
  isoform_exprs <- vstExprs(gene_cds)
  isoform_exprs <- t(scale(t(isoform_exprs)))
  isoform_exprs[isoform_exprs > 3] <- 3
  isoform_exprs[isoform_exprs < -3] <- -3
  isoform_exprs[is.na(isoform_exprs)] <- 0
  #pheatmap(isoform_exprs, cluster_cols=F, show_rownames=F, show_colnames=F)
  ph_res <- plot_pseudotime_heatmap(gene_cds, 
                                    use_gene_short_name=FALSE, 
                                    num_clusters=num_clusters, 
                                    add_annotation_row=add_annotation_row,
                                    return_heatmap=TRUE)
  return(ph_res)
}

num_iso_switch_clusters <- 6

#pdf ("iso_switch_heatmap.pdf", height=6, width=6)
iso_switch_heatmap <- plot_isoform_expression_heatmap(HSMM_myo_isoform, 
                                                     #c("TPM1", "TPM2"),
                                                     subset(iso_switch_test_res, qval < 0.01)$gene_short_name, 
                                                     num_clusters=num_iso_switch_clusters)
#dev.off()

tree_row <- iso_switch_heatmap$tree_row
tree_row_gene_clusters <- cutree(tree_row, num_iso_switch_clusters)
fData_iso <- fData(HSMM_myo_isoform[names(tree_row_gene_clusters),])
fData_iso$isoform_switching_cluster <- tree_row_gene_clusters

cluster_summaries <- fData_iso %>% 
    mutate(total_genes = length(unique(gene_id))) %>% 
    group_by(isoform_switching_cluster) %>% 
    summarise(num_isoforms = n(),
              num_genes = length(unique(gene_id)), 
              fraction_genes=length(unique(gene_id))/unique(total_genes)[1])
cluster_fraction_genes <- cluster_summaries$fraction_genes
names(cluster_fraction_genes) <- cluster_summaries$isoform_switching_cluster

additional_row_annotation <- data.frame(cluster_fraction_genes = cluster_fraction_genes[fData_iso$isoform_switching_cluster])
row.names(additional_row_annotation) <- row.names(fData_iso)

pdf ("fig5a_1.pdf", height=6, width=6)
iso_switch_heatmap <- plot_isoform_expression_heatmap(HSMM_myo_isoform, 
                                                     #c("TPM1", "TPM2"),
                                                     subset(iso_switch_test_res, qval < 0.01)$gene_short_name, 
                                                     num_clusters=num_iso_switch_clusters,
                                                     add_annotation_row=additional_row_annotation)
dev.off()


gene_summaries <- fData_iso %>% add_rownames() %>%
    group_by(gene_id) %>% 
    summarise(num_isoforms=length(unique(rowname)), 
              num_clusters=length(unique(isoform_switching_cluster)))


pdf("fig5b.pdf", width=1.25, height=1)
qplot(num_clusters, data=gene_summaries) + 
  xlab("Isoform clusters per gene") + 
  ylab("Isoform-switching\ngenes")+
  monocle:::monocle_theme_opts() + nm_theme()
dev.off()

actins <- HSMM_myo_isoform[row.names(subset(fData(HSMM_myo_isoform), 
                             gene_short_name %in% c("ACTA1", "ACTA2", "ACTB", "ACTG1", "ACTG2", "ACTC1") & num_cells_expressed >= 15)),]
#actins <- differentialGeneTest(actins)
pdf ("fig5c.pdf", width=1.25, height=4.5)
plot_genes_in_pseudotime(actins, 
                         color_by="gene_short_name", 
                         label_by_short_name=FALSE,
                         cell_size=1,
                         min_expr=0.5) + 
  theme(legend.position="none") + 
  scale_color_brewer(palette="Set1") + nm_theme()
dev.off()

isoform_universe <- row.names(subset(fData(HSMM_myo_isoform), gene_id %in% HSMM_expressed_genes))
iso_switch_gsa_hyper_heatmap <- collect_gsa_hyper_results(HSMM_myo_isoform[isoform_universe,], human_reactome_gsc, tree_row_gene_clusters)
pdf ("iso_switch_gsa_hyper_heatmap_reactome_gsc.pdf", height=6, width=6)
plot_gsa_hyper_heatmap(HSMM_myo_isoform,  iso_switch_gsa_hyper_heatmap, significance=0.1)
dev.off()

isoform_universe <- row.names(subset(fData(HSMM_myo_isoform), gene_id %in% HSMM_expressed_genes))
iso_switch_gsa_hyper_heatmap <- collect_gsa_hyper_results(HSMM_myo_isoform[isoform_universe,], go_bp_gsc, tree_row_gene_clusters)
pdf ("iso_switch_gsa_hyper_heatmap_go_bp.pdf", height=16, width=6)
plot_gsa_hyper_heatmap(HSMM_myo_isoform,  iso_switch_gsa_hyper_heatmap, significance=0.1)
dev.off()

iso_switch_gsa_hyper_heatmap <- collect_gsa_hyper_results(HSMM_myo_isoform[isoform_universe,], go_mf_gsc, tree_row_gene_clusters)
pdf ("iso_switch_gsa_hyper_heatmap_go_mf.pdf", height=6, width=6)
plot_gsa_hyper_heatmap(HSMM_myo_isoform,  iso_switch_gsa_hyper_heatmap, significance=0.1)
dev.off()

pdf("diff_gene_heatmap.pdf")
plot_pseudotime_heatmap(HSMM_myo[target_genes,])
dev.off()

num_de_gene_clusters <- 5
de_genes <- unique(row.names(subset(diff_gene_test_res, qval < 0.1)))
de_genes <- de_genes[duplicated(fData(HSMM_myo[de_genes,])$gene_short_name) == FALSE]

pdf("diff_gene_heatmap.pdf")
de_heatmap <- plot_pseudotime_heatmap(HSMM_myo[de_genes,], return_heatmap=TRUE, num_clusters=num_de_gene_clusters)
dev.off()

tree_row <- de_heatmap$tree_row
tree_row_gene_clusters <- cutree(tree_row, num_de_gene_clusters)

gene_universe <- row.names(subset(fData(HSMM_myo), gene_id %in% HSMM_expressed_genes))
de_gsa_hyper_heatmap <- collect_gsa_hyper_results(HSMM_myo[gene_universe,], go_cc_gsc, tree_row_gene_clusters, gsSizeLim=c(10, Inf))
pdf ("de_gsa_hyper_heatmap_go_cc.pdf", height=12, width=6)
plot_gsa_hyper_heatmap(HSMM_myo_isoform,  de_gsa_hyper_heatmap, significance=0.1)
dev.off()

gene_universe <- row.names(subset(fData(HSMM_myo), gene_id %in% HSMM_expressed_genes))
de_gsa_hyper_heatmap <- collect_gsa_hyper_results(HSMM_myo[gene_universe,], go_mf_gsc, tree_row_gene_clusters, gsSizeLim=c(10, Inf))
pdf ("de_gsa_hyper_heatmap_go_mf.pdf", height=6, width=6)
plot_gsa_hyper_heatmap(HSMM_myo_isoform,  de_gsa_hyper_heatmap, significance=0.1)
dev.off()


gene_universe <- row.names(subset(fData(HSMM_myo), gene_id %in% HSMM_expressed_genes))
de_gsa_hyper_heatmap <- collect_gsa_hyper_results(HSMM_myo[gene_universe,], go_bp_gsc, tree_row_gene_clusters, gsSizeLim=c(15, Inf))
pdf ("de_gsa_hyper_heatmap_go_bp.pdf", height=48, width=6)
plot_gsa_hyper_heatmap(HSMM_myo_isoform,  de_gsa_hyper_heatmap, significance=0.1)
dev.off()

####

library(biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host="www.ensembl.org") #useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl",mart)
listFilters(mart)

listFilters(mart)[70,]

library(Gviz)
biomTrack <- BiomartGeneRegionTrack(genome = "hg19", name = "ENSEMBL", symbol = "ABCB5", transcriptAnnotation = "symbol", filters=list(hgnc_symbol="ABCB5"))
plotTracks(biomTrack)


HSMM_myo_isoform <- detectGenes(HSMM_myo_isoform, min_expr=1)

TPM1_isoform_ids <- row.names(subset(fData(HSMM_myo_isoform), 
                         gene_short_name %in% c("TPM1") & num_cells_expressed >= 15))
TPM1_isoform_ids <- setdiff(TPM1_isoform_ids, c("ENST00000559108.1", "ENST00000560131.1"))
TPM1_isoforms <- HSMM_myo_isoform[TPM1_isoform_ids,]

# plot_isoform_covariance(TPM1_isoforms, "TPM1")

# plot_genes_in_pseudotime(TPM1_isoforms, 
#                          color_by="Time", 
#                          label_by_short_name=FALSE,
#                          min_expr=0.1)

library(stringr)
transcript_ids <- row.names(TPM1_isoforms)
transcript_ids <- str_split_fixed(transcript_ids, "\\.", 2)[,1]
biomTrack <- BiomartGeneRegionTrack(genome = "hg19", name = "GENCODE", symbol = "TPM1", filters=list(ensembl_transcript_id=transcript_ids))

pdf("TPM1_gviz.pdf", width=12, height=2)
plotTracks(biomTrack, transcriptAnnotation = "transcript", shape = "box", )
dev.off()

pdf("TPM1_collapsed.pdf", width=4, height=1)
plotTracks(biomTrack, collapseTranscripts = "meta", shape = "box", transcriptAnnotation = "symbol")
dev.off()

all_ids <- c("ENST00000357980.4", 
                "ENST00000267996.7",
                "ENST00000288398.6",
                "ENST00000358278.3",
                "ENST00000560445.1",
                "ENST00000403994.3",
                "ENST00000559556.1",
                "ENST00000404484.4")

# exon1a_ids <- c("ENST00000357980.4", "ENST00000267996.7")
# exon1b_ids <- c("ENST00000288398.6",
#                 "ENST00000358278.3",
#                 "ENST00000560445.1",
#                 "ENST00000403994.3",
#                 "ENST00000559556.1",
#                 "ENST00000404484.4")



exon2a_ids <- c("ENST00000357980.4", 
                "ENST00000267996.7")
exon2b_ids <- c("ENST00000357980.4", 
                "ENST00000288398.6",
                "ENST00000358278.3",
                "ENST00000560445.1",
                "ENST00000403994.3",
                "ENST00000559556.1",
                "ENST00000404484.4")

exon2a_expr <- signif(colSums(exprs(TPM1_isoforms[exon2a_ids,])))
exon2b_expr <- signif(colSums(exprs(TPM1_isoforms[exon2b_ids,])))

exon2_total_expr <- signif(colSums(exprs(TPM1_isoforms[unique(c(exon2a_ids, exon2b_ids)),])))

exon2 <- rbind(data.frame(Pseudotime=pData(TPM1_isoforms)$Pseudotime, mRNAs=exon2a_expr, PSI=exon2a_expr/exon2_total_expr, event="Exon 2a"),
               data.frame(Pseudotime=pData(TPM1_isoforms)$Pseudotime, mRNAs=exon2b_expr, PSI=exon2b_expr/exon2_total_expr, event="Exon 2b"))

pdf("fig5e.pdf", width=1.25, height=1.5)
qplot(Pseudotime, PSI, data=exon2, color=event, size=I(0.5)) + stat_smooth(method="gam", color="black", formula=y~sm.ns(x, df=5), size=I(0.5), se=F) + 
  facet_wrap(~event, ncol=1) + 
  theme(legend.position="none") + 
  scale_color_manual(values=list("Exon 2a"="#00A550", "Exon 2b"="#00ADEF"))+ nm_theme() + 
  monocle:::monocle_theme_opts() + scale_y_continuous(limits=c(0,1), labels = scales::percent)
dev.off()     
      

exon6a_ids <- c("ENST00000288398.6",
                "ENST00000358278.3",
                "ENST00000404484.4")
exon6b_ids <- c("ENST00000357980.4", 
                "ENST00000267996.7",
                "ENST00000560445.1",
                "ENST00000403994.3",
                "ENST00000559556.1")

exon6a_expr <- signif(colSums(exprs(TPM1_isoforms[exon6a_ids,])))
exon6b_expr <- signif(colSums(exprs(TPM1_isoforms[exon6b_ids,])))

exon6_total_expr <- signif(colSums(exprs(TPM1_isoforms[unique(c(exon6a_ids, exon6b_ids)),])))

exon6 <- rbind(data.frame(Pseudotime=pData(TPM1_isoforms)$Pseudotime, mRNAs=exon6a_expr, PSI=exon6a_expr/exon6_total_expr, event="Exon 6a"),
               data.frame(Pseudotime=pData(TPM1_isoforms)$Pseudotime, mRNAs=exon6b_expr, PSI=exon6b_expr/exon6_total_expr, event="Exon 6b"))

pdf("fig5f.pdf", width=1.25, height=1.5)
qplot(Pseudotime, PSI, data=exon6, color=event, size=I(0.5)) + stat_smooth(method="gam", color="black", formula=y~sm.ns(x, df=5), size=I(0.5), se=F) + 
  facet_wrap(~event, ncol=1) + 
  scale_color_manual(values=list("Exon 6a"="#ED1C24", "Exon 6b"="#2E3092"))+
  theme(legend.position="none") + nm_theme() + 
  monocle:::monocle_theme_opts() + scale_y_continuous(limits=c(0,1), labels = scales::percent)
dev.off()     
      


exon9a_ids <- c("ENST00000288398.6",
                "ENST00000403994.3"
                )
exon9b_ids <- c("ENST00000288398.6",
                "ENST00000403994.3",
                "ENST00000560445.1")

exon9d_ids <- c("ENST00000357980.4", 
                "ENST00000267996.7",
                #"ENST00000560445.1",
                "ENST00000358278.3",
                "ENST00000404484.4",
                "ENST00000559556.1")

exon9a_expr <- signif(colSums(exprs(TPM1_isoforms[exon9a_ids,])))
exon9b_expr <- signif(colSums(exprs(TPM1_isoforms[exon9b_ids,])))
exon9d_expr <- signif(colSums(exprs(TPM1_isoforms[exon9d_ids,])))
exon9_total_expr <- signif(colSums(exprs(TPM1_isoforms[unique(c(exon9a_ids, exon9b_ids, exon9d_ids)),])))

exon9 <- rbind(data.frame(Pseudotime=pData(TPM1_isoforms)$Pseudotime, mRNAs=exon9a_expr, PSI=exon9a_expr/exon9_total_expr, event="Exon 9a"),
               data.frame(Pseudotime=pData(TPM1_isoforms)$Pseudotime, mRNAs=exon9b_expr, PSI=exon9b_expr/exon9_total_expr, event="Exon 9b"),
               data.frame(Pseudotime=pData(TPM1_isoforms)$Pseudotime, mRNAs=exon9d_expr, PSI=exon9d_expr/exon9_total_expr, event="Exon 9d"))

pdf("fig5g.pdf", width=1.25, height=1.5)
qplot(Pseudotime, PSI, data=exon9, color=event, size=I(0.5)) + stat_smooth(method="gam", color="black", formula=y~sm.ns(x, df=5), size=I(0.5), se=F) + 
  facet_wrap(~event, ncol=1) + 
  theme(legend.position="none") + 
  scale_color_manual(values=list("Exon 9a"="#00A79D", "Exon 9b"="#1B75BB", "Exon 9d"="#A87B4F"))+ nm_theme() + 
  monocle:::monocle_theme_opts() + scale_y_continuous(limits=c(0,1), breaks = c(0, 0.5, 1), labels = scales::percent)
dev.off()     

####
   

promoter1_ids <- c("ENST00000357980.4", 
                "ENST00000267996.7",
                "ENST00000288398.6",
                "ENST00000358278.3",
                "ENST00000560445.1",
                "ENST00000403994.3",
                "ENST00000559556.1")
promoter2_ids <- c(
                "ENST00000404484.4")

promoter1_expr <- signif(colSums(exprs(TPM1_isoforms[promoter1_ids,])))
promoter2_expr <- signif(colSums(exprs(TPM1_isoforms[promoter2_ids,])))

tpm_total_expr <- signif(colSums(exprs(TPM1_isoforms[unique(c(promoter1_ids, promoter2_ids)),])))

promoter_df <- rbind(data.frame(Pseudotime=pData(TPM1_isoforms)$Pseudotime, mRNAs=promoter1_expr, PSI=promoter1_expr/tpm_total_expr, event="Promoter 1"),
               data.frame(Pseudotime=pData(TPM1_isoforms)$Pseudotime, mRNAs=promoter2_expr, PSI=promoter2_expr/tpm_total_expr, event="Promoter 2"))

#pdf("exon9.pdf", width=1.25, height=1.5)
qplot(Pseudotime, PSI, data=promoter_df, color=event, size=I(0.15)) + stat_smooth(method="gam", color="black", formula=y~sm.ns(x, df=5), size=I(0.5), se=F) + 
  facet_wrap(~event, ncol=1) + 
  theme(legend.position="none") + 
  monocle:::monocle_theme_opts() + scale_y_continuous(limits=c(0,1), labels = scales::percent)
#dev.off()     

#################
# Now repeat the isoform switching test with read counts:


HSMM_isoform_count_matrix <- read.delim("~/Dropbox (Cole Trapnell's Lab)/Shared Data/Muscle/HSMM/sc-RNA-Seq/HSMM_cuffnorm_out/isoforms.count_table")
row.names(HSMM_isoform_count_matrix) <- HSMM_isoform_count_matrix$tracking_id
HSMM_isoform_count_matrix <- HSMM_isoform_count_matrix[,-1]
HSMM_isoform_count_matrix <- HSMM_isoform_count_matrix[,row.names(sample_sheet)]


pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = isoform_ann)
HSMM_isoform_counts<- newCellDataSet(as.matrix(HSMM_isoform_count_matrix), 
                                     phenoData = pd, 
                                     featureData = fd,
                                     expressionFamily=negbinomial.size(), 
                                     lowerDetectionLimit=1)
HSMM_isoform_counts <- detectGenes(HSMM_isoform_counts, min_expr = 1)

HSMM_isoform_counts <- HSMM_isoform_counts[,row.names(pData(HSMM_myo))]                 

HSMM_isoform_counts <-estimateSizeFactors(HSMM_isoform_counts)
HSMM_isoform_counts <-estimateDispersions(HSMM_isoform_counts)

pData(HSMM_isoform_counts)$Pseudotime <- pData(HSMM_myo_isoform)$Pseudotime

count_iso_switch_test_res <- iso_switch_test(HSMM_isoform_counts[target_isoforms,], cores=detectCores() / 2)
count_iso_switch_test_res <- merge(fData(HSMM_myo), count_iso_switch_test_res, by="row.names")
count_iso_switch_test_res$qval <- 1.0
count_iso_switch_test_res$qval[count_iso_switch_test_res$status == "OK"] <- p.adjust(subset(count_iso_switch_test_res, status == "OK")$pval, method="fdr")

#plot_genes_in_pseudotime(HSMM_isoform_counts[row.names(subset(fData(HSMM_isoform_counts), gene_short_name %in% c("TPM1") & num_cells_expressed >=15 )),], color_by="Time", label_by_short_name=FALSE, min_expr=0.5)

setdiff(subset(count_iso_switch_test_res, qval < 0.01)$gene_short_name, subset(iso_switch_test_res, qval < 0.01)$gene_short_name)

setdiff(subset(iso_switch_test_res, qval < 0.01)$gene_short_name, subset(count_iso_switch_test_res, qval < 0.01)$gene_short_name)



pdf ("iso_switch_heatmap.pdf", height=6, width=6)
iso_switch_heatmap_counts <- plot_isoform_expression_heatmap(HSMM_isoform_counts, 
                                                     #c("TPM1", "TPM2"),
                                                     subset(count_iso_switch_test_res, qval < 0.01)$gene_short_name, 
                                                     num_clusters=num_iso_switch_clusters)
dev.off()
### 
tree_row <- iso_switch_heatmap_counts$tree_row
tree_row_gene_clusters <- cutree(tree_row, num_iso_switch_clusters)
fData_iso <- fData(HSMM_isoform_counts[names(tree_row_gene_clusters),])
fData_iso$isoform_switching_cluster <- tree_row_gene_clusters

cluster_summaries <- fData_iso %>% 
    mutate(total_genes = length(unique(gene_id))) %>% 
    group_by(isoform_switching_cluster) %>% 
    summarise(num_isoforms = n(),
              num_genes = length(unique(gene_id)), 
              fraction_genes=length(unique(gene_id))/unique(total_genes)[1])
cluster_fraction_genes <- cluster_summaries$fraction_genes
names(cluster_fraction_genes) <- cluster_summaries$isoform_switching_cluster

additional_row_annotation <- data.frame(cluster_fraction_genes = cluster_fraction_genes[fData_iso$isoform_switching_cluster])
row.names(additional_row_annotation) <- row.names(fData_iso)

pdf ("fig5a_2.pdf", height=6, width=6)
iso_switch_heatmap <- plot_isoform_expression_heatmap(HSMM_myo_isoform, 
                                                     #c("TPM1", "TPM2"),
                                                     subset(iso_switch_test_res, qval < 0.01)$gene_short_name, 
                                                     num_clusters=num_iso_switch_clusters,
                                                     add_annotation_row=additional_row_annotation)
dev.off()



###


test_gene <- HSMM_myo_isoform[row.names(subset(fData(HSMM_myo_isoform), 
                             gene_short_name %in% c("VASH2"))),]
#actins <- differentialGeneTest(actins)
pdf ("test_gene_RPC.pdf", width=1.25, height=8)
plot_genes_in_pseudotime(test_gene, 
                         color_by="Time", 
                         label_by_short_name=FALSE,
                         cell_size=0.15,
                         min_expr=0.5) + 
  theme(legend.position="none")  +
  monocle:::monocle_theme_opts()
  #scale_color_brewer(palette="Set1")
dev.off()

test_gene <- HSMM_isoform_counts[row.names(subset(fData(HSMM_isoform_counts), 
                             gene_short_name %in% c("VASH2"))),]
#actins <- differentialGeneTest(actins)
pdf ("test_gene_read_counts.pdf", width=1.25, height=8)
plot_genes_in_pseudotime(test_gene, 
                         color_by="Time", 
                         label_by_short_name=FALSE,
                         cell_size=0.15,
                         min_expr=0.5) + 
  theme(legend.position="none")  +
  monocle:::monocle_theme_opts()
  #scale_color_brewer(palette="Set1")
dev.off()


####


TPM2 <- HSMM_myo_isoform[row.names(subset(fData(HSMM_myo_isoform), 
                             gene_short_name %in% c("TPM2") & num_cells_expressed >= 15)),]
#actins <- differentialGeneTest(actins)
pdf ("TPM2_RPC.pdf", width=1.25, height=4.5)
plot_genes_in_pseudotime(TPM2, 
                         color_by="Time", 
                         label_by_short_name=FALSE,
                         cell_size=0.15,
                         min_expr=0.5) + 
  theme(legend.position="none")  +
  monocle:::monocle_theme_opts()
  #scale_color_brewer(palette="Set1")
dev.off()


TNNT2 <- HSMM_myo_isoform[row.names(subset(fData(HSMM_myo_isoform), 
                             gene_short_name %in% c("TNNT2") & num_cells_expressed >= 15)),]
#actins <- differentialGeneTest(actins)
pdf ("TNNT2_RPC.pdf", width=1.25, height=6)
plot_genes_in_pseudotime(TNNT2, 
                         color_by="Time", 
                         label_by_short_name=FALSE,
                         cell_size=0.15,
                         min_expr=0.5) + 
  theme(legend.position="none")  +
  monocle:::monocle_theme_opts()
  #scale_color_brewer(palette="Set1")
dev.off()


MYL6 <- HSMM_myo_isoform[row.names(subset(fData(HSMM_myo_isoform), 
                             gene_short_name %in% c("MYL6") & num_cells_expressed >= 15)),]
#actins <- differentialGeneTest(actins)
pdf ("MYL6_RPC.pdf", width=1.25, height=6)
plot_genes_in_pseudotime(MYL6, 
                         color_by="Time", 
                         label_by_short_name=FALSE,
                         cell_size=0.15,
                         min_expr=0.5) + 
  theme(legend.position="none")  +
  monocle:::monocle_theme_opts()
  #scale_color_brewer(palette="Set1")
dev.off()

###

TPM2 <- HSMM_isoform_counts[row.names(subset(fData(HSMM_myo_isoform), 
                             gene_short_name %in% c("TPM2") & num_cells_expressed >= 15)),]
#actins <- differentialGeneTest(actins)
pdf ("TPM2_read_count.pdf", width=1.25, height=4.5)
plot_genes_in_pseudotime(TPM2, 
                         color_by="Time", 
                         label_by_short_name=FALSE,
                         cell_size=0.15,
                         min_expr=0.5) + 
  theme(legend.position="none")  +
  monocle:::monocle_theme_opts()
  #scale_color_brewer(palette="Set1")
dev.off()


TNNT2 <- HSMM_isoform_counts[row.names(subset(fData(HSMM_myo_isoform), 
                             gene_short_name %in% c("TNNT2") & num_cells_expressed >= 15)),]
#actins <- differentialGeneTest(actins)
pdf ("TNNT2_read_count.pdf", width=1.25, height=6)
plot_genes_in_pseudotime(TNNT2, 
                         color_by="Time", 
                         label_by_short_name=FALSE,
                         cell_size=0.15,
                         min_expr=0.5) + 
  theme(legend.position="none")  +
  monocle:::monocle_theme_opts()
  #scale_color_brewer(palette="Set1")
dev.off()


MYL6 <- HSMM_isoform_counts[row.names(subset(fData(HSMM_myo_isoform), 
                             gene_short_name %in% c("MYL6") & num_cells_expressed >= 15)),]
#actins <- differentialGeneTest(actins)
pdf ("MYL6_read_count.pdf", width=1.25, height=6)
plot_genes_in_pseudotime(MYL6, 
                         color_by="Time", 
                         label_by_short_name=FALSE,
                         cell_size=0.15,
                         min_expr=0.5) + 
  theme(legend.position="none")  +
  monocle:::monocle_theme_opts()
  #scale_color_brewer(palette="Set1")
dev.off()


# save.image('./RData/analysis_hsmm_isoform.RData')


