#' Results from \code{\link{distinct_test}} function
#' 
#' @rdname res
#' @name res
#' @aliases res
#' 
#' @param res contains a \code{\linkS4class{data.frame}} object,
#' with the results obtained applying \code{\link{distinct_test}} function to \code{\link{Kang_subset}} dataset.
#' Below the code used to obtain `res`.
#
#' @examples
#' # load the input data:
#' # data("Kang_subset", package = "distinct")
#' # Kang_subset
#' # 
#' # create the design of the study:
#' # samples = Kang_subset@metadata$experiment_info$sample_id
#' # group = Kang_subset@metadata$experiment_info$stim
#' # design = model.matrix(~group)
#' # rownames of the design must indicate sample ids:
#' # rownames(design) = samples
#' # design
#' # 
#' # Note that the sample names in `colData(x)$name_sample` have to be the same ones as those in `rownames(design)`.
#' # rownames(design)
#' # unique(SingleCellExperiment::colData(Kang_subset)$sample_id)
#' # 
#' # In order to obtain a finer ranking for the most significant genes, if computational resources are available, we encourage users to increase P_4 (i.e., the number of permutations when a raw p-value is < 0.001) and set P_4 = 20,000 (by default P_4 = 10,000).
#' # 
#' # The group we would like to test for is in the second column of the design, therefore we will specify: column_to_test = 2
#' # 
#' # set.seed(61217)
#' # res = distinct_test(
#' #   x = Kang_subset, 
#' #   name_assays_expression = "logcounts",
#' #   name_cluster = "cell",
#' #   design = design,
#' #   column_to_test = 2,
#' #   min_non_zero_cells = 20,
#' #   n_cores = 2)
#' #   
#' #   save(res, file = "res.RData")
#' #   saveRDS(res, file = "res.rds")
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{distinct_test}}
NULL