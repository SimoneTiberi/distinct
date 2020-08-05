# this is a R wrapper, to parallelize computation of distinct test via foreach:
perm_test_parallel_covariates_R = function(P, P_2, P_3, P_4, 
                                N_breaks, cluster_ids, sample_ids, 
                                n_samples, group_ids_of_samples, min_non_zero_cells, counts,
                                n_cores, design_batches){
  suppressWarnings({
    cl <- makeCluster(n_cores, setup_strategy = "sequential")
  })
  registerDoParallel(cl, n_cores);
  
  # change the final_order to sort clusters by number of cells:
  final_order = order(table(cluster_ids), decreasing = TRUE) - 1
  # -1 because index goes from 0 to n-1 in C++!
  
  p_values_ALL = foreach(cl_id = final_order,
                         .packages=c("distinct", "Matrix"), # Matrix package needed to treat the sparce matrix
                         .errorhandling = "stop") %dorng%{
                           # select elements belonging to the cl_id-th cluster:
                           sel = cluster_ids == cl_id
                           
                           .Call(`_distinct_perm_test_parallel_covariates`,
                                 P,
                                 P_2,
                                 P_3,
                                 P_4,# number of permutations
                                 N_breaks, # number of breaks at which to evaluate the cdf
                                 sample_ids[sel], # ids of samples for every cell
                                 n_samples, # total number of samples
                                 group_ids_of_samples, # ids of groups (1 or 2) for every sample
                                 min_non_zero_cells, # min number of cells with > 0 expression in each group
                                 counts[ sel, ],
                                 design_batches)[[1]]
                         }
  stopCluster(cl)
  stopImplicitCluster()
  
  # order(final_order) will revert the original oder:
  do.call(cbind, p_values_ALL[order(final_order)])
}
