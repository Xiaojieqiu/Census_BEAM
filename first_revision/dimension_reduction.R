# Test different dimension reduction approaches
# load('./RData/prepare_lung_data.RData')
# load('./RData/analysis_shalek_data.RData')

library(lle)
library(tsne)
library(monocle)
library(xacHelper) #diffusion_maps, diffusion_maps_bw, tSNE, LLE
load_all_libraries()

reduceDimension_custom <- function (cds, max_components = 2, method = c("ICA", "DDRTree"),
    pseudo_expr = NULL, residualModelFormulaStr = NULL, use_vst = NULL,
    verbose = FALSE, use_irlba = NULL, ...)
{
    FM <- exprs(cds)
    if (is.null(use_irlba)) {
        message("Warning: argument 'use_irlba' is deprecated and will be removed in a future release")
    }
    if (is.null(use_vst) && cds@expressionFamily@vfamily == "negbinomial") {
        use_vst = TRUE
        pseudo_expr = 0
    }
    if (cds@expressionFamily@vfamily == "negbinomial") {
        if (is.null(use_vst))
            use_vst = TRUE
        if (is.null(pseudo_expr))
            pseudo_expr = 0
    }
    else {
        if (is.null(use_vst))
            use_vst = FALSE
        if (is.null(pseudo_expr))
            pseudo_expr = 1
    }
    if (use_vst == FALSE && cds@expressionFamily@vfamily == "negbinomial") {
        checkSizeFactors(cds)
        size_factors <- sizeFactors(cds)
        FM <- Matrix::t(Matrix::t(FM)/size_factors)
    }
    if (is.null(fData(cds)$use_for_ordering) == FALSE && nrow(subset(fData(cds),
        use_for_ordering == TRUE)) > 0)
        FM <- FM[fData(cds)$use_for_ordering, ]
    if (cds@expressionFamily@vfamily == "binomialff") {
        ncounts <- FM > 0
        ncounts[ncounts != 0] <- 1
        FM <- Matrix::t(Matrix::t(ncounts) * log(1 + ncol(ncounts)/rowSums(ncounts)))
    }
    if (cds@expressionFamily@vfamily != "binomialff") {
        FM <- FM + pseudo_expr
    }
    FM <- FM[apply(FM, 1, sd) > 0, ]
    if (cds@expressionFamily@vfamily != "binomialff") {
        if (use_vst) {
            VST_FM <- vstExprs(cds, expr_matrix = FM, round_vals = FALSE)
            if (is.null(VST_FM) == FALSE) {
                FM <- VST_FM
            }
            else {
                stop("Error: set the variance-stabilized value matrix with vstExprs(cds) <- computeVarianceStabilizedValues() before calling this function with use_vst=TRUE")
            }
        }
        else {
            FM <- log2(FM)
        }
    }
    if (is.null(residualModelFormulaStr) == FALSE) {
        if (verbose)
            message("Removing batch effects")
        X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr),
            data = pData(cds), drop.unused.levels = TRUE)
        fit <- limma::lmFit(FM, X.model_mat, ...)
        beta <- fit$coefficients[, -1, drop = FALSE]
        beta[is.na(beta)] <- 0
        FM <- as.matrix(FM) - beta %*% t(X.model_mat[, -1])
        if (cds@expressionFamily@vfamily != "binomialff") {
            if (use_vst == FALSE) {
                FM <- 2^FM
            }
        }
    }
    if (verbose)
        message("Reducing to independent components")
    FM <- Matrix::t(scale(Matrix::t(FM)))
    FM <- FM[apply(FM, 1, sd) > 0, ]
    if (nrow(FM) == 0) {
        stop("Error: all rows have standard deviation zero")
    }
    if (is.function(method)) {        
        reducedDim <- method(FM, ...)
        return(reducedDim)
    }
}

#test the LLE/tSNE function: 
# LLE_dm <- LLE(log2(exprs(AT12_cds_subset) + 1 ))

# pdf('./nbt_2nd_sub_reviewers/LLE_dim_reduction.pdf', width = 3, height = 2)
# qplot(LLE_dm[1, ], LLE_dm[2, ], color = pData(AT12_cds_subset)$Time) + monocle_theme_opts()
# dev.off()

# tSNE_dm <- tSNE(log2(exprs(AT12_cds_subset) + 1))
# pdf('./nbt_2nd_sub_reviewers/LLE_dim_reduction.pdf', width = 3, height = 2)
# qplot(tSNE_dm[1, ], tSNE_dm[2, ], color = pData(AT12_cds_subset)$Time) + monocle_theme_opts()
# dev.off()

## test passing LLE/tSNE into the reduceDimension function 

#################################### Lung data ####################################
# pd <- new("AnnotatedDataFrame", data = pData(AT12_cds_subset))
# fd <- new("AnnotatedDataFrame", data = fData(AT12_cds_subset))
# test_cds <-  newCellDataSet(exprs(AT12_cds_subset), 
#                                    phenoData = pd, 
#                                    featureData = fd, 
#                                    expressionFamily=tobit(), 
#                                    lowerDetectionLimit=1)

lung_LLE_res <- reduceDimension_custom(AT12_cds_subset, method = LLE)
pdf('./nbt_2nd_sub_reviewers/lung_LLE_dim_reduction.pdf', width = 3, height = 2)
qplot(lung_LLE_res[1, ], lung_LLE_res[2, ], color = pData(AT12_cds_subset)$Time) + 
    xlab('LLE Reduce Dimension 1') + ylab('LLE Reduce Dimension 2') + nm_theme() 
dev.off()

lung_tsne_res <- reduceDimension_custom(AT12_cds_subset, method = tSNE)
pdf('./nbt_2nd_sub_reviewers/lung_tsne_dim_reduction.pdf', width = 3, height = 2)
qplot(lung_tsne_res[1, ], lung_tsne_res[2, ], color = pData(AT12_cds_subset)$Time) + 
    xlab('tSNE Reduce Dimension 1') + ylab('tSNE Reduce Dimension 2') + nm_theme() 
dev.off()

# diffusion_maps_dm <- diffusion_maps(log2(exprs(AT12_cds_subset) + 1 ))
# qplot(diffusion_maps_dm[1, ], diffusion_maps_dm[2, ], color = pData(AT12_cds_subset)$Time)


#test on passing diffusion map to reduceDimension_custom function: 
#note that diffusion map function is implemented in xacHelper helper package 
#which is based on code from Fabian J. Theis1, etc. group. Destiny package can be also used 

lung_diffusion_maps <- reduceDimension_custom(AT12_cds_subset, method = diffusion_maps)
pdf('./nbt_2nd_sub_reviewers/lung_dm_dim_reduction.pdf', width = 3, height = 2)
qplot(lung_diffusion_maps[1, ], lung_diffusion_maps[2, ], color = pData(AT12_cds_subset)$Time) + 
    xlab('LLE Reduce Dimension 1') + ylab('LLE Reduce Dimension 2') + nm_theme() 
dev.off()

#################################### Shalek data ####################################
# pd <- new("AnnotatedDataFrame", data = pData(AT12_cds_subset))
# fd <- new("AnnotatedDataFrame", data = fData(AT12_cds_subset))
# test_cds <-  newCellDataSet(exprs(AT12_cds_subset), 
#                                    phenoData = pd, 
#                                    featureData = fd, 
#                                    expressionFamily=tobit(), 
#                                    lowerDetectionLimit=1)

ko_LLE_res <- reduceDimension_custom(Shalek_abs_subset_ko_LPS, method = LLE)
pdf('./nbt_2nd_sub_reviewers/Shalek_ko_LLE_dim_reduction.pdf', width = 3, height = 2)
qplot(ko_LLE_res[1, ], ko_LLE_res[2, ], color = pData(Shalek_abs_subset_ko_LPS)$stim_time) + ##check
    xlab('LLE Reduce Dimension 1') + ylab('LLE Reduce Dimension 2') + nm_theme() 
dev.off()

ko_tsne_res <- reduceDimension_custom(Shalek_abs_subset_ko_LPS, method = tSNE)
pdf('./nbt_2nd_sub_reviewers/Shalek_ko_tsne_dim_reduction.pdf', width = 3, height = 2)
qplot(ko_tsne_res[1, ], ko_tsne_res[2, ], color = pData(Shalek_abs_subset_ko_LPS)$stim_time) + 
    xlab('tSNE Reduce Dimension 1') + ylab('tSNE Reduce Dimension 2') + nm_theme() 
dev.off()

ko_diffusion_maps <- reduceDimension_custom(Shalek_abs_subset_ko_LPS, method = diffusion_maps)
pdf('./nbt_2nd_sub_reviewers/Shalek_ko_dm_dim_reduction.pdf', width = 3, height = 2)
qplot(ko_diffusion_maps[1, ], ko_diffusion_maps[2, ], color = pData(Shalek_abs_subset_ko_LPS)$stim_time) + 
    xlab('DM Reduce Dimension 1') + ylab('DM Reduce Dimension 2') + nm_theme() 
dev.off()

Golgi_LLE_res <- reduceDimension_custom(Shalek_golgi_update, method = LLE)
pdf('./nbt_2nd_sub_reviewers/Shalek_golgi_LLE_dim_reduction.pdf', width = 3, height = 2)
qplot(Golgi_LLE_res[1, ], Golgi_LLE_res[2, ], color = pData(Shalek_golgi_update)$stim_time) + 
    xlab('LLE Reduce Dimension 1') + ylab('LLE Reduce Dimension 2') + nm_theme() 
dev.off()

Golgi_tsne_res <- reduceDimension_custom(Shalek_golgi_update, method = tSNE)
pdf('./nbt_2nd_sub_reviewers/Shalek_golgi_tsne_dim_reduction.pdf', width = 3, height = 2)
qplot(Golgi_tsne_res[1, ], Golgi_tsne_res[2, ], color = pData(Shalek_golgi_update)$stim_time) + 
    xlab('tSNE Reduce Dimension 1') + ylab('tSNE Reduce Dimension 2') + nm_theme() 
dev.off()

Golgi_diffusion_maps <- reduceDimension_custom(Shalek_golgi_update, method = diffusion_maps)
pdf('./nbt_2nd_sub_reviewers/Shalek_golgi_dm_dim_reduction.pdf', width = 3, height = 2)
qplot(Golgi_diffusion_maps[1, ], Golgi_diffusion_maps[2, ], color = pData(Shalek_golgi_update)$stim_time) + 
    xlab('DM Reduce Dimension 1') + ylab('DM Reduce Dimension 2') + nm_theme() 
dev.off()

save.image('./RData/dimension_reduction.RData')















