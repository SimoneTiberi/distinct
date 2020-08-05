#' Subset from the `Kang18_8vs8()` object of the \code{muscData} package.
#' 
#' @rdname Kang_subset
#' @name Kang_subset
#' @aliases Kang_subset
#' 
#' @param Kang_subset contains a \code{\linkS4class{SingleCellExperiment}} object, 
#' representing a subset of 6 samples (3 individuals observed across 2 conditions) and 100 genes selected from the `Kang18_8vs8()` object of the \code{muscData} package.
#' Below the code used to subset the original dataset.
#
#' @examples
#' #################### 
#' # Object 'Kang_subset' is generated as follows:
#' ####################
#' # library(muscData)
#' # sce = Kang18_8vs8()
#' # 
#' # library(scater)
#' # sce = computeLibraryFactors(sce)
#' # sce = logNormCounts(sce)
#' # 
#' # Select genes with at least 1000 non-zero cells:
#' # sce = sce[ rowSums(assays(sce)$counts > 0) >= 1000, ]
#' # 
#' # randomly select 100 of these genes:
#' # set.seed(61217)
#' # sel = sample( rownames(sce), size = 100)
#' # sce = sce[ rownames(sce) %in% sel, ]
#' # 
#' # select 3 individuals only:
#' # ind_selected = levels(factor(colData(sce)$ind))[1:3]
#' # sce = sce[, sce$ind %in% ind_selected]
#' # 
#' # make a sample_id column:
#' # colData(sce)$sample_id = factor(paste(colData(sce)$stim, colData(sce)$ind, sep = "_"))
#' # 
#' # create an experiment_info object containing sample-group information:
#' # experiment_info = unique(data.frame(sample_id = colData(sce)$sample_id, 
#' #                                     stim = colData(sce)$stim) )
#' # metadata(sce)$experiment_info = data.frame(experiment_info, row.names = NULL)
#' # 
#' # remove unnecessary information to reduce storage space:
#' # sce$cluster = NULL;
#' # sce$multiplets = NULL;
#' # rowData(sce) = NULL;
#' # colnames(sce) = NULL;
#' # reducedDim(sce) = NULL
#' #
#' # Kang_subset = sce
#' # save(Kang_subset, file = "Kang_subset.RData")
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{distinct_test}}
NULL