#' Test for differential state between two groups of samples, based on scRNA-seq data.
#'
#' \code{distinct_test} tests for differential state between two groups of samples, based on scRNA-seq data.
#' 
#' @param x a \code{\linkS4class{SummarizedExperiment}} object.
#' @param name_assays_expression a character ("counts" by default), 
#' indicating the name of the assays(x) element which stores the expression data (i.e., assays(x)$name_assays_expression).
#' @param name_cluster a character ("cluster_id" by default), 
#' indicating the name of the colData(x) element which stores the cluster id of each cell (i.e., colData(x)$name_colData_cluster).
#' @param name_sample a character ("sample_id" by default), 
#' indicating the name of the colData(x) element which stores the sample id of each cell (i.e., colData(x)$name_colData_sample).
#' @param name_group a character ("group_id" by default), 
#' indicating the name of the colData(x) element which stores the group id of each cell (i.e., colData(x)$name_colData_group).
#' @param name_metadata_experiment_info a character ("experiment_info" by default), 
#' indicating the name of the metadata(x) element which stores the group id of each cell (i.e., metadata(x)$name_metadata_experiment_info).
#' @param logarithm a boolean, indicating whether to use raw counts (if FALSE, by default) or log(counts + 1) (if TRUE).
#' @param P  the number of permutations to use.
#' @param N_breaks the number of breaks at which to evaluate the comulative density function.
#' @param min_non_zero_cells_per_group the minimum number of non-zero cells in each group for a gene to be evaluated.
#' @return A \code{\linkS4class{data.frame}} object.
#' @examples
#' library(HDCytoData)
#' main = Weber_BCR_XL_sim_main_SE()
#' main
#' 
#' library(diffcyt)
#' # normalize data:
#' main <- transformData(main, cofactor = 5)
#' 
#' # select all State markers:
#' # sel_cols = colData(main)$marker_class == "state"
#' # select 1 marker only:
#' sel_cols = colData(main)$marker_name == "pNFkB"
#' 
#' library(SummarizedExperiment)
#' # create a SummarizedExperiment object with seleceted columns, and inverting row-column structure:
#' sce <- SummarizedExperiment(
#'   assays = list(exprs = t(assays(main)$exprs[, sel_cols]) ), 
#'   colData = rowData(main), 
#'   rowData = colnames(main)[sel_cols], 
#'   metadata = list(experiment_info =  metadata(main)$experiment_info ))
#'   sce
#'   
#' # Perform differential analyses, within each cluster of cells, between conditions:
#' res = distinct_test(
#'   sce, 
#'   name_assays_expression = "exprs",
#'   name_cluster = "population_id",
#'   name_group = "group_id",
#'   name_sample = "sample_id",
#'   P = 10^2, 
#'   min_non_zero_cells_per_group = 0)
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
                         name_metadata_experiment_info = "experiment_info",
                         logarithm = FALSE, 
                         P = 10^3, 
                         N_breaks = 10, 
                         min_non_zero_cells_per_group = 20){
  # x is a single cell experiment
  
  # check counts, clusters, samples, etc... are all available.
  # make it general, maybe also working without a SummarizedExperiment structure?
  
  # TODO: remove genes with only 0 counts from 'counts'
  
  # count matrix:
  sel = which(names(assays(x)) == name_assays_expression)
  if( length(sel) == 0 ){
    return("'name_assays_expression' not found in names(assays(x))")
  }
  if( length(sel) > 1 ){
    return("'name_assays_expression' found multiple times in names(assays(x))")
  }
  counts = assays(x)[[sel]]
  
  # cluster ids:
  sel = which(names(colData(x)) == name_cluster)
  if( length(sel) == 0 ){
    return("'name_cluster' not found in names(colData(x))")
  }
  if( length(sel) > 1 ){
    return("'name_cluster' found multiple times in names(colData(x))")
  }
  cluster_ids = factor(colData(x)[[sel]])
  n_clusters <- nlevels(cluster_ids)
  cluster_ids_num = as.numeric(cluster_ids)-1
  
  # sample ids:
  sel = which(names(colData(x)) == name_sample)
  if( length(sel) == 0 ){
    return("'name_sample' not found in names(colData(x))")
  }
  if( length(sel) > 1 ){
    return("'name_sample' found multiple times in names(colData(x))")
  }
  sample_ids = factor(colData(x)[[sel]])
  
  # select experimental info:
  sel = which(names(metadata(x)) == name_metadata_experiment_info)
  if( length(sel) == 0 ){
    return("'name_metadata_experiment_info' not found in names(metadata(x))")
  }
  if( length(sel) > 1 ){
    return("'name_metadata_experiment_info' found multiple times in names(metadata(x))")
  }
  experimental_info =  metadata(x)[[sel]]
  
  # sample ids from experiment_info:
  sel_sample = which(colnames(experimental_info) == name_sample)
  if( length(sel_sample) == 0 ){
    return("'name_sample' not found in colnames(metadata(x)$name_metadata_experiment_info)")
  }
  if( length(sel_sample) > 1 ){
    return("'name_sample' found multiple times in colnames(metadata(x)$name_metadata_experiment_info)")
  }
  
  levels(sample_ids) = factor(experimental_info[,sel_sample])
  n_samples <- nlevels(sample_ids)
  sample_ids_num = as.numeric(sample_ids)-1
  
  # group ids from experiment_info:
  sel_group = which(colnames(experimental_info) == name_group)
  if( length(sel_group) == 0 ){
    return("'name_group' not found in colnames(metadata(x)$name_metadata_experiment_info)")
  }
  if( length(sel_group) > 1 ){
    return("'name_group' found multiple times in colnames(metadata(x)$name_metadata_experiment_info)")
  }
  
  group_ids_of_samples = factor(experimental_info[,sel_group])
  group_ids_of_samples = as.numeric(group_ids_of_samples)
  
  # group ids (1 per cell)
  sel = which(names(colData(x)) == name_group)
  if( length(sel) == 0 ){
    return("'name_group' not found in names(colData(x))")
  }
  if( length(sel) > 1 ){
    return("'name_group' found multiple times in names(colData(x))")
  }
  group_ids = factor(colData(x)[[sel]])
  group_ids = as.numeric(group_ids)
  n_groups = nlevels(factor(group_ids))
  
  groups = unique(group_ids_of_samples)
  n_samples_per_group = vapply( groups, function(g) sum(group_ids_of_samples == g), FUN.VALUE = numeric(1) )
  # n_samples_per_group contains the samples of each group (e.g., 3 2)
  n_samples_per_group_per_sample = n_samples_per_group[match(group_ids_of_samples,groups)]
  # n_samples_per_group_per_sample contains the samples of each group that samples belong to, (e.g., 3 3 3 2 2) 
  
  if(n_groups < 2){
    return("At least 2 groups should be provided")
  }
  
  # REMOVE WHEN implementing comparisons btw more than 2 groups:
  if(n_groups > 2){
    return("At most 2 groups should be provided: comparisons between more than 2 groups will be implemented (soon) in future releases.")
  }
  
  # TODO: print something before running C++ functions
  
  p_val = .Call(`_distinct_perm_test`,
                logarithm, # if TRUE, use log2(counts + 1); if FALSE, use counts
                P, # number of permutations
                N_breaks, # number of breaks at which to evaluate the cdf
                cluster_ids_num, # ids of clusters (cell-population) for every cell
                n_clusters, # total number of clusters
                sample_ids_num, # ids of samples for every cell
                n_samples, # total number of samples
                group_ids_of_samples, # ids of groups (1 or 2) for every sample
                min_non_zero_cells_per_group, # min number of cells with > 0 expression in each group
                group_ids, # ids of groups (cell-population) for every cell
                counts,
                1)[[1]] # [[1]]: results returned as a 1 element list
  
  # TODO: print something after running C++ functions
  
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
