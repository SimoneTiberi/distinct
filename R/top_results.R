#' Test for differential state between two groups of samples, based on scRNA-seq data.
#'
#' \code{top_results} returns the significant results obtained via \code{\link{distinct_test}}.
#' 
#' @param res a \code{\linkS4class{data.frame}} with results as returned from \code{\link{distinct_test}}.
#' @param cluster a character indicating the cluster(s) whose results have to be returned. 
#' Results from all clusters are returned by default ("all").
#' @param significance numeric, results with adjusted p-value < significance will be returned.
#' @param global logical indicating whether to filter results according to p_adj.glb (when TRUE), or p_adj.loc (when FALSE).
#' @return A \code{\linkS4class{data.frame}} object.
#' Columns `gene` and `cluster_id` contain the gene and cell-cluster name, while `p_val`, `p_adj.loc` and `p_adj.glb` report the raw p-values, locally and globally adjusted p-values, via Benjamini and Hochberg (BH) correction.
#' In locally adjusted p-values (`p_adj.loc`) BH correction is applied in each cluster separately, while in globally adjusted p-values (`p_adj.glb`) BH correction is performed to the results from all clusters.
#' @examples
#' # load the input data:
#' data("Kang_subset", package = "distinct")
#' Kang_subset
#' 
#' # create the design of the study:
#' samples = Kang_subset@metadata$experiment_info$sample_id
#' group = Kang_subset@metadata$experiment_info$stim
#' design = model.matrix(~group)
#' # rownames of the design must indicate sample ids:
#' rownames(design) = samples
#' design
#' 
#' # The group we would like to test for is in the second column of the design, therefore we will specify: column_to_test = 2
#' 
#' set.seed(61217)
#' res = distinct_test(
#'   x = Kang_subset, 
#'   name_assays_expression = "logcounts",
#'   name_cluster = "cell",
#'   design = design,
#'   column_to_test = 2,
#'   min_non_zero_cells = 20,
#'   n_cores = 2)
#' 
#' # Visualize significant results:
#' top_results(res)
#' 
#' # Visualize significant results from a specified cluster of cells:
#' top_results(res, cluster = "Dendritic cells")
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{plot_cdfs}}, \code{\link{plot_densities}}
#' 
#' @export
top_results = function(res, 
                       cluster = "all",
                       significance = 0.01,
                       global = TRUE){
  stopifnot(
    is.data.frame(res),
    is.character(cluster), length(cluster) > 0L,
    is.numeric(significance), length(significance) == 1L,
    is.logical(global), length(global) == 1L
  )
  
  if(global){
    sel_rows = res$p_adj.glb < significance
  }else{
    sel_rows = res$p_adj.loc < significance
  }
  
  # if cluster is not "all", then I must filter clusters too:
  if(cluster != "all"){
    # check that "cluster" is among cluster names in cluster_id slot:
    if( !any(cluster %in% res$cluster_id) ){
      message("'cluster' not found in 'res$cluster_id' make sure cluster is typed correctly.")
      return(NULL)
    }
    
    sel_cluster = res$cluster %in% cluster
    sel_rows = sel_rows & sel_cluster
  }
  
  res[sel_rows,]
}
