#' Test for differential state between two groups of samples, based on scRNA-seq data.
#'
#' \code{distinct_test} tests for differential state between two groups of samples.
#' 
#' @param x a \code{\linkS4class{SummarizedExperiment}} or a \code{\linkS4class{SingleCellExperiment}} object.
#' @param name_assays_expression a character ("counts" by default), 
#' indicating the name of the assays(x) element which stores the expression data (i.e., assays(x)$name_assays_expression).
#' @param name_cluster a character ("cluster_id" by default), 
#' indicating the name of the colData(x) element which stores the cluster id of each cell (i.e., colData(x)$name_colData_cluster).
#' @param name_sample a character ("sample_id" by default), 
#' indicating the name of the colData(x) element which stores the sample id of each cell (i.e., colData(x)$name_colData_sample).
#' @param name_group a character ("group_id" by default), 
#' indicating the name of the colData(x) element which stores the group id of each cell (i.e., colData(x)$name_colData_group).
#' @param logarithm a boolean, indicating whether to use raw counts (if FALSE, by default) or log(counts + 1) (if TRUE).
#' @param P  the number of permutations to use.
#' @param N_breaks the number of breaks at which to evaluate the comulative density function.
#' @param min_non_zero_cells the minimum number of non-zero cells in each cluster for a gene to be evaluated.
#' @return A \code{\linkS4class{data.frame}} object.
#' Columns `gene` and `cluster_id` contain the gene and cell-cluster name, while `p_val`, `p_adj.loc` and `p_adj.glb` report the raw p-values, locally and globally adjusted p-values, via Benjamini and Hochberg (BH) correction.
#' In locally adjusted p-values (`p_adj.loc`) BH correction is applied in each cluster separately, while in globally adjusted p-values (`p_adj.glb`) BH correction is performed to the results from all clusters.
#' @examples
#' data("Kang_subset", package = "distinct")
#' Kang_subset
#' 
#' set.seed(61217)
#' res = distinct_test(
#'   x = Kang_subset, 
#'   name_assays_expression = "counts",
#'   name_cluster = "cell",
#'   name_group = "stim",
#'   name_sample = "sample_id",
#'   P = 10^3, 
#'   min_non_zero_cells = 20)
#' 
#' # Visualize results:
#' head(res)
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @export
distinct_test = function(x, 
                         name_assays_expression = "counts",
                         name_cluster = "cluster_id",
                         name_sample = "sample_id",
                         name_group = "group_id",
                         logarithm = FALSE, 
                         P = 10^3, 
                         N_breaks = 10, 
                         min_non_zero_cells = 20){
  stopifnot(
    ( is(x, "SummarizedExperiment") | is(x, "SingleCellExperiment") ),
    is.character(name_assays_expression), length(name_assays_expression) == 1L,
    is.character(name_cluster), length(name_cluster) == 1L,
    is.character(name_sample), length(name_sample) == 1L,
    is.character(name_group), length(name_group) == 1L,
    is.logical(logarithm), length(logarithm) == 1L,
    is.numeric(P), length(P) == 1L,
    is.numeric(N_breaks), length(N_breaks) == 1L,
    is.numeric(min_non_zero_cells), length(min_non_zero_cells) == 1L
  )
  
  # lower-bound for min_non_zero_cells:
  if(min_non_zero_cells < 1){
    message("'min_non_zero_cells' must be at least 1.")
    return(NULL)
  }
  
  # count matrix:
  sel = which(names(assays(x)) == name_assays_expression)
  if( length(sel) == 0 ){
    message("'name_assays_expression' not found in names(assays(x))")
    return(NULL)
  }
  if( length(sel) > 1 ){
    message("'name_assays_expression' found multiple times in names(assays(x))")
    return(NULL)
  }
  counts = assays(x)[[sel]]
  counts = as.matrix(counts)
  # remove rows with 0 counts:
  counts = counts[ rowSums(counts > 0) > 0, ]
  
  # cluster ids:
  sel = which(names(colData(x)) == name_cluster)
  if( length(sel) == 0 ){
    message("'name_cluster' not found in names(colData(x))")
    return(NULL)
  }
  if( length(sel) > 1 ){
    message("'name_cluster' found multiple times in names(colData(x))")
    return(NULL)
  }
  cluster_ids = factor(colData(x)[[sel]])
  n_clusters = nlevels(cluster_ids)
  cluster_ids_num = as.numeric(cluster_ids)-1
  
  # sample ids:
  sel = which(names(colData(x)) == name_sample)
  if( length(sel) == 0 ){
    message("'name_sample' not found in names(colData(x))")
    return(NULL)
  }
  if( length(sel) > 1 ){
    message("'name_sample' found multiple times in names(colData(x))")
    return(NULL)
  }
  sample_ids = factor(colData(x)[[sel]])
  
  # group ids (1 per cell)
  sel = which(names(colData(x)) == name_group)
  if( length(sel) == 0 ){
    message("'name_group' not found in names(colData(x))")
    return(NULL)
  }
  if( length(sel) > 1 ){
    message("'name_group' found multiple times in names(colData(x))")
    return(NULL)
  }
  group_ids = factor(colData(x)[[sel]])
  group_ids = as.numeric(group_ids)
  n_groups = nlevels(factor(group_ids))
  
  # select experimental info:
  experiment_info = unique(data.frame(sample_id = sample_ids, 
                                      group_id = group_ids, 
                                      stringsAsFactors = FALSE) )
  
  # sample ids from experiment_info:
  sample_ids = factor(sample_ids, levels = experiment_info$sample_id)
  n_samples = nlevels(sample_ids)
  sample_ids_num = as.numeric(sample_ids)-1
  
  # group ids from experiment_info:
  group_ids_of_samples = factor(experiment_info$group_id)
  group_ids_of_samples = as.numeric(group_ids_of_samples)
  
  groups = unique(group_ids_of_samples)
  n_samples_per_group = vapply( groups, function(g) sum(group_ids_of_samples == g), FUN.VALUE = numeric(1) )
  # n_samples_per_group contains the samples of each group (e.g., 3 2)
  n_samples_per_group_per_sample = n_samples_per_group[match(group_ids_of_samples,groups)]
  # n_samples_per_group_per_sample contains the samples of each group that samples belong to, (e.g., 3 3 3 2 2) 
  
  if(n_groups < 2){
    message("At least 2 groups should be provided")
    return(NULL)
  }
  
  # REMOVE WHEN implementing comparisons btw more than 2 groups:
  if(n_groups > 2){
    message("At most 2 groups should be provided: comparisons between more than 2 groups will be implemented (soon) in future releases.")
    return(NULL)
  }
  
  message("Data loaded, starting differential testing")
  
  p_val = .Call(`_distinct_perm_test`,
                logarithm, # if TRUE, use log2(counts + 1); if FALSE, use counts
                P, # number of permutations
                N_breaks, # number of breaks at which to evaluate the cdf
                cluster_ids_num, # ids of clusters (cell-population) for every cell
                n_clusters, # total number of clusters
                sample_ids_num, # ids of samples for every cell
                n_samples, # total number of samples
                group_ids_of_samples, # ids of groups (1 or 2) for every sample
                min_non_zero_cells, # min number of cells with > 0 expression in each group
                counts,
                1)[[1]] # [[1]]: results returned as a 1 element list
  
  message("Differential testing completed, returning results")
  
  # set -1s to NA, so that we don't use these elements when adjusting p-values:
  p_val[ p_val == -1 ] = NA
  
  # locally adjusted p-values:
  res_adjusted_locally = apply(p_val, 2, p.adjust, method = "BH")
  
  p_val = c(p_val)
  res_adjusted_locally = c(res_adjusted_locally)
  # globally adjusted p-values:
  res_adjusted_globally = p.adjust(p_val, method = "BH")
  
  gene_names = rownames(counts)
  if(is.null(gene_names)){
    gene_names = seq_len(nrow(counts))
  }
  
  res = data.frame(
    gene = rep(gene_names, times = n_clusters),
    cluster_id = rep( levels(cluster_ids), each = nrow(counts) ),
    p_val = p_val,
    p_adj.loc = res_adjusted_locally,
    p_adj.glb = res_adjusted_globally
  )
  
  # set to 1 pvals (and adjusted pvals) which were NA (not analyzed:)
  res$p_val[is.na(res$p_val)] = 1
  res$p_adj.loc[is.na(res$p_adj.loc)] = 1
  res$p_adj.glb[is.na(res$p_adj.glb)] = 1
  
  # TODO: return an object that can be accesses via "topTags", "gene" or similar(see PB objects):
  res
}
