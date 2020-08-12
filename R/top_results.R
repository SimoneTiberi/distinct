#' Filter significant results.
#'
#' \code{top_results} returns the significant results obtained via \code{\link{distinct_test}}.
#' 
#' @param res a \code{\linkS4class{data.frame}} with results as returned from \code{\link{distinct_test}}.
#' @param cluster a character indicating the cluster(s) whose results have to be returned. 
#' Results from all clusters are returned by default ("all").
#' @param significance numeric, results with adjusted p-value < significance will be returned.
#' @param global logical indicating whether to filter results according to p_adj.glb (when TRUE), or p_adj.loc (when FALSE).
#' @param up_down a character indicating whether to return: all results ("both" or "BOTH"), only up-regulated results ("up" or "UP") or down-regulated results ("down" or "DOWN").
#' In `res`, a FC > 1 (or log2FC > 0) indicates up-regulation of group1 (compared to group2); while a FC < 1 (or log2FC < 0) indicates down-regulation of group1 (compared to group2).
#' #' @param up_down a character indicating how results should be sorted.
#' Results can be sorted by globally adjusted p-value ("p_adj.glb", default choice), locally adjusted p-value ("p_adj.loc"), raw p-value ("p_val") or (log2)fold-change ("FC" or "log2FC").
#' @return A \code{\linkS4class{data.frame}} object.
#' Columns `gene` and `cluster_id` contain the gene and cell-cluster name, while `p_val`, `p_adj.loc` and `p_adj.glb` report the raw p-values, locally and globally adjusted p-values, via Benjamini and Hochberg (BH) correction.
#' In locally adjusted p-values (`p_adj.loc`) BH correction is applied in each cluster separately, while in globally adjusted p-values (`p_adj.glb`) BH correction is performed to the results from all clusters.
#' @examples
#' # Check the examples of 'distinct_test'
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{distinct_test}}, \code{\link{log2_FC}}
#' 
#' @export
top_results = function(res, 
                       cluster = "all",
                       significance = 0.01,
                       global = TRUE,
                       up_down = "both",
                       sort_by = "p_adj.glb"){
  stopifnot(
    is.data.frame(res),
    is.character(cluster), length(cluster) > 0L,
    is.numeric(significance), length(significance) == 1L,
    is.logical(global), length(global) == 1L,
    is.character(up_down), length(up_down) == 1L,
    is.character(sort_by), length(sort_by) == 1L
  )
  
  if( !(up_down %in% c("both", "BOTH", "up", "UP", "down", "DOWN")) ){
    message("'up_down' should be one of: 'both', 'BOTH', 'up', 'UP', 'down', 'DOWN'")
    return(NULL)
  }
  
  if( !(sort_by %in% c("p_adj.glb", "p_adj.loc", "p_val", "FC", "log2FC")) ){
    message("'sort_by' should be one of: 'p_adj.glb', 'p_adj.loc', 'p_val', 'FC', 'log2FC'")
    return(NULL)
  }
  
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
  
  # filter by gene (and cluster)
  res = res[sel_rows,]
  
  # filter by UP-down regulation if required:
  if(up_down %in% c("up", "UP", "down", "DOWN")){
    sel_log2FC = grep("log2FC", colnames(res))
    if(length(sel_log2FC) == 0){
      message("To filter up- or down-regulated results, first compute 'FC' and 'log2FC' via 'log2_FC' function; otherwise avoid filtering and set up_down = 'both'.")
      return(NULL)
    }
    
    if(up_down %in% c("up", "UP")){ # filter UP-regulated results only:
      sel_genes = res[,sel_log2FC] >= 0
      res = res[ sel_genes, ]
    }else{
      sel_genes = res[,sel_log2FC] < 0
      res = res[ sel_genes, ]
    }
  }
  
  # order results:
  if( sort_by == "p_adj.glb"){
    res = res[ order(res$p_adj.glb), ]
  }
  if( sort_by == "p_adj.loc"){
    res = res[ order(res$p_adj.loc), ]
  }
  if( sort_by == "p_val"){
    res = res[ order(res$p_val), ]
  }
  
  if( sort_by %in% c("FC", "log2FC") ){
    sel_log2FC = grep("log2FC", colnames(res))
    if(length(sel_log2FC) == 0){
      message("To sort by FC or log2FC, first compute 'FC' and 'log2FC' via 'log2_FC' function; otherwise sort by (raw or adjusted) p-value.")
      return(NULL)
    }
    res = res[ order( abs(res[ ,sel_log2FC ]), decreasing = TRUE), ]
  }
  
  res
}
