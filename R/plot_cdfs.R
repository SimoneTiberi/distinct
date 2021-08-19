#' Plot sample-specific CDFs.
#'
#' \code{plot_cdfs} returns a \code{\link{ggplot}} object
#' with the sample-specific cumulative density function (CDF) estimates, for a specified cluster and gene.
#' 
#' @param x a \code{\linkS4class{SummarizedExperiment}} or a \code{\linkS4class{SingleCellExperiment}} object.
#' @param name_assays_expression a character ("logcounts" by default), 
#' indicating the name of the assays(x) element which stores the expression data (i.e., assays(x)$name_assays_expression).
#' We strongly encourage using normalized data, such as counts per million (CPM) or log-CPM.
#' @param name_cluster a character ("cluster_id" by default), 
#' indicating the name of the colData(x) element which stores the cluster id of each cell (i.e., colData(x)$name_colData_cluster).
#' @param name_sample a character ("sample_id" by default), 
#' indicating the name of the colData(x) element which stores the sample id of each cell (i.e., colData(x)$name_colData_sample).
#' @param name_group a character ("group_id" by default), 
#' indicating the name of the colData(x) element which stores the group id of each cell (i.e., colData(x)$name_colData_group).
#' @param cluster a character, indicating the name of the cluster to plot.
#' @param gene a character, indicating the name of the gene to plot.
#' @param group_level a logical, indicating whether to plot group-level (if TRUE) or sample-level curves (if FALSE).
#' @param pad a logical element indicating whether to plot the lines of the CDF when 0 and 1 (TRUE) or not (FALSE).
#' @param size a numeric argument defining the width of lines, passed to \code{\link{stat_ecdf}}.
#' @return A \code{\link{ggplot}} object.
#' @examples
#' data("Kang_subset", package = "distinct")
#' Kang_subset
#' 
#' plot_cdfs(x = Kang_subset,
#'           gene = "ISG15",
#'           cluster = "Dendritic cells",
#'           name_assays_expression = "logcounts",
#'           name_cluster = "cell",
#'           name_sample = "sample_id",
#'           name_group = "stim")
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{distinct_test}}, \code{\link{plot_densities}}
#' 
#' @export
plot_cdfs = function(x, 
                     name_assays_expression = "logcounts",
                     name_cluster = "cluster_id",
                     name_sample = "sample_id",
                     name_group = "group_id",
                     cluster,
                     gene,
                     group_level = FALSE,
                     pad = TRUE,
                     size = 0.75){
  
  stopifnot(
    ( is(x, "SummarizedExperiment") | is(x, "SingleCellExperiment") ),
    is.character(name_assays_expression), length(name_assays_expression) == 1L,
    is.character(name_cluster), length(name_cluster) == 1L,
    is.character(name_sample), length(name_sample) == 1L,
    is.character(name_group), length(name_group) == 1L,
    is.character(cluster), length(cluster) == 1L,
    is.character(gene), length(gene) == 1L,
    is.logical(group_level), length(group_level) == 1L
  )
  
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
  cluster_ids = colData(x)[[sel]]
  
  # sample ids:
  if(!group_level){
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
  }
  
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
  
  # sel gene:
  sel_gene = rownames(x) == gene
  sel_gene[is.na(sel_gene)] = FALSE
  
  if( all(!sel_gene) ){ # if all sel_gene FALSE:
    message("'gene' not found in `rownames(x)`")
    return(NULL)
  }
  
  # sel cluster:
  sel_cluster = cluster_ids == cluster
  sel_cluster[is.na(sel_cluster)] = FALSE
  
  if( all(!sel_cluster) ){ # if all sel_gene FALSE:
    message("'cluster' not found in `colData(x)[[name_cluster]]`")
    return(NULL)
  }
  
  
  if(group_level){
    DF = data.frame(x = counts[sel_gene, sel_cluster], 
                    group = group_ids[sel_cluster])
    gg = ggplot(DF, aes(x=x, group=group, colour = group))
  }else{
    DF = data.frame(x = counts[sel_gene, sel_cluster], 
                    group = group_ids[sel_cluster],
                    sample = sample_ids[sel_cluster])
    gg = ggplot(DF, aes(x=x, group=sample, colour = group))
  }
  
  # Density plot:
  gg +
    geom_hline(yintercept=0, linetype="dashed", color = "grey") +
    geom_hline(yintercept=1, linetype="dashed", color = "grey") +
    stat_ecdf(pad = pad, size = size) +
    theme_bw() + 
    theme(panel.grid = element_blank()) +
    labs(title = paste(cluster, "-", gene),
         x = name_assays_expression,
         y = "CDF")
}
