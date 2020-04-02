// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// Compare each group to a "baseline": a mean group

// [[Rcpp::export]]
Rcpp::NumericMatrix
  perm_test_MG_baseline(bool const& logarithm, // if TRUE, use log2(counts + 1); if FALSE, use counts
                        unsigned int const& P, // number of permutations
                        unsigned int const& N_breaks, // number of breaks at which to evaluate the cdf
                        Rcpp::IntegerVector const& cluster_ids, // ids of clusters (cell-population) for every cell
                        unsigned int const& n_clusters, // total number of clusters
                        Rcpp::IntegerVector const& sample_ids, // ids of samples for every cell
                        unsigned int const& n_samples, // total number of samples
                        Rcpp::IntegerVector& group_ids_of_samples, // ids of groups (1 or 2) for every sample
                        double const& min_counts_per_group, // min number of cells with > 0 expression in each group
                        Rcpp::IntegerVector const& n_samples_per_group_per_sample,
                        Rcpp::NumericMatrix const& counts, // count matrix (rows = genes, cols = cells)
                        unsigned int const& n_groups) // number of groups
  {
    unsigned int n_genes = counts.nrow();
    double T_obs, T_tmp;
    
    Rcpp::NumericVector seq_len_N = as<NumericVector>(wrap(seq(1, N_breaks)));
    
    unsigned int n_cells; // number of cells PER group
    bool keep_gene; // logical , indicating whether to keep genes or not (if filtered due to low expression)
    Rcpp::IntegerVector sample_ids_one_cl, sel_both; // sel_A and sel_B define the indeces of the counts matrix referring to group A and B
    
    // create arma vectors, for group and cluster ids
    //arma::vec group_ids_arma = Rcpp::as<arma::vec>(group_ids);
    arma::vec cluster_ids_arma = Rcpp::as<arma::vec>(cluster_ids);
    
    // Logical vector to check if samples have at least 1 cell
    Rcpp::LogicalVector sample_has_cells(n_samples);
    
    // turn sample_id into numeric:
    group_ids_of_samples = group_ids_of_samples-1; // -1 to make the index start from 0
    //arma::vec group_ids_of_samples_arma = Rcpp::as<arma::vec>(group_ids_of_samples);
    //arma::uvec samples_in_group_A_arma = arma::find( group_ids_of_samples_arma == 0);
    //arma::uvec samples_in_group_B_arma = arma::find( group_ids_of_samples_arma == 1);
    //Rcpp::IntegerVector samples_in_group_A = Rcpp::wrap(samples_in_group_A_arma);
    //Rcpp::IntegerVector samples_in_group_B = Rcpp::wrap(samples_in_group_B_arma);
    
    //unsigned int n_samples_A = samples_in_group_A.length(); // number of samples in group A
    //unsigned int n_samples_B = samples_in_group_B.length(); // number of samples in group B
    // CHECK n_samples_A + n_samples_B == n_samples !!!
    arma::uvec ids;
    
    Rcpp::NumericVector counts_one_gene_permuted, T_perm(P), mean_cdf_A(N_breaks), mean_cdf_baseline(N_breaks), cdfs_one_break, counts_one_gene, counts_one_gene_A, counts_one_gene_B, breaks, counts_one_gene_one_sample;
    Rcpp::NumericMatrix mean_cdf(N_breaks, n_groups), res(n_genes, n_clusters); // res: final p-values;
    std::fill(res.begin(), res.end(), -1);
    
    for (unsigned int cl_id = 0; cl_id < n_clusters; ++cl_id) { // loop over clusters
      // select the cells belonging to the cluster 'cl_id':
      ids = arma::find( cluster_ids_arma == cl_id );
      sel_both = Rcpp::wrap(ids);
      
      sample_ids_one_cl = sample_ids[sel_both];
      
      // calculate the total number of cells in the cluster 'cl_id' (across both groups)
      n_cells = sel_both.length(); 
      
      // compute P permutations of the n_cells cells (sample without replacement).
      Rcpp::IntegerMatrix PERMUTATIONS(n_cells,P);
      Rcpp::IntegerVector permutation = Rcpp::seq_len(n_cells)-1; // define the vector for each permutation, -1 because Rcpp indexes start from 0 !
      for (unsigned int p = 0; p < P; ++p) {
        // permutation = Rcpp::seq_len(n_cells)-1;
        std::random_shuffle(permutation.begin(), permutation.end(), randWrapper);
        PERMUTATIONS(_,p) = permutation;
        // using Rcpp samplers will allow to control the seed from R.
      }
      
      //Rcpp::IntegerMatrix sel(n_groups, n_cells);
      //for (unsigned int group = 0; group < n_groups; ++group) {
      // select data belonging to the first and second group:
      //  ids = arma::find( group_ids_arma == (group+1) && cluster_ids_arma == cl_id );
      //  sel(group,_) = Rcpp::wrap(ids);
      //}
      
      // CHECK that every sample has at least 1 cell, independent of the gene (to be checked 1 for every cluster):
      for (unsigned int sample = 0; sample < n_samples; ++sample) {
        sample_has_cells[sample] = is_true(Rcpp::any(sample_ids_one_cl == sample));
      }
      // SPEED UP: COMPUTE "sample_ids_one_cl == sample" ONCE FOR EVERY CLUSTER ONLY HERE!!!
      
      // Rcpp::sugar provides a comparator with one value for vectors but not for matrices
      for (unsigned int gene = 0; gene < n_genes; ++gene) {
        counts_one_gene = counts(gene,_); // SELECT cl_id too!!!
        
        // counts_one_gene_A = counts_one_gene[sel_A];
        // counts_one_gene_B = counts_one_gene[sel_B];
        // TO DO: sub-select groups with > min counts and test differences between those only ?
        
        counts_one_gene = counts_one_gene[sel_both];
        
        keep_gene = (sum(counts_one_gene > 0) >= (n_groups * min_counts_per_group));
        // keep_gene is a filter to remove lowly abundant genes: in each group, we want at least `min_counts_per_group` cells to have non-zero expression.
        
        if(keep_gene){
          if(logarithm){ // log data if `logarithm` is true
            for (unsigned int cell = 0; cell < counts_one_gene.length(); ++cell) {
              counts_one_gene[cell] = log2(counts_one_gene[cell] + 1);
            }
          }
          
          // compute the grid at which we ought to calculate the cdf:
          double min_ = min(counts_one_gene);
          double len_out = ( max(counts_one_gene) - min_ )/(N_breaks+1);
          breaks = min_ + len_out * seq_len_N; 
          
          // loop over samples in each group (already computed above, mayube in R)
          // do mean_cdf for each group already / n_samples_group
          std::fill(mean_cdf.begin(), mean_cdf.end(), 0);
          for (unsigned int sample = 0; sample < n_samples; ++sample) {
            // if no sample_ids_one_cl == sample ???
            // put a contraint here!
            if(sample_has_cells[sample]){
              counts_one_gene_one_sample = counts_one_gene[sample_ids_one_cl == sample];
              // compute the cdf of each sample in counts_one_gene_list[sample]:
              for (unsigned int b = 0; b < N_breaks ; ++b) {
                mean_cdf(b, group_ids_of_samples[sample]) += Rcpp::mean( breaks[b] > counts_one_gene_one_sample )/n_samples_per_group_per_sample[sample];
                // put counts_one_gene[sample_ids_one_cl == sample] instead of counts_one_gene_one_sample
              }
            }
          }
          
          // create a baseline cdf to compare against (mean_cdf_baseline) based on ALL cells from ALL samples in ALL groups:
          for (unsigned int b = 0; b < N_breaks ; ++b) {
            mean_cdf_baseline[b] = Rcpp::mean( breaks[b] > counts_one_gene );
          }
          
          // loop over group g1 = 1:(n_groups-1)
          // loop over group g2 = (g1+1):n_groups
          // add all differences btw mean cdfs
          T_obs = 0;
          for (unsigned int group_1 = 0; group_1 < n_groups; ++group_1) {
            mean_cdf_A = mean_cdf(_, group_1);
            T_obs += Rcpp::sum(Rcpp::abs(mean_cdf_A - mean_cdf_baseline));
          }
          
          //  test permuted data:
          std::fill(T_perm.begin(), T_perm.end(), 0);
          for (unsigned int p = 0; p < P; ++p) {
            // obtain the permuted data:
            counts_one_gene_permuted = counts_one_gene[PERMUTATIONS(_,p)];
            
            // loop over samples in each group (already computed above, mayube in R)
            // do mean_cdf for each group already / n_samples_group
            std::fill(mean_cdf.begin(), mean_cdf.end(), 0);
            for (unsigned int sample = 0; sample < n_samples; ++sample) {
              // if no sample_ids_one_cl == sample ???
              // put a contraint here!
              if(sample_has_cells[sample]){
                counts_one_gene_one_sample = counts_one_gene_permuted[sample_ids_one_cl == sample];
                // compute the cdf of each sample in counts_one_gene_list[sample]:
                for (unsigned int b = 0; b < N_breaks ; ++b) {
                  mean_cdf(b, group_ids_of_samples[sample]) += Rcpp::mean( breaks[b] > counts_one_gene_one_sample )/n_samples_per_group_per_sample[sample];
                  // put counts_one_gene_permuted[sample_ids_one_cl == sample] instead of counts_one_gene_one_sample
                }
              }
            }
            
            // create a baseline cdf to compare against (mean_cdf_baseline) based on ALL cells from ALL samples in ALL groups:
            //for (unsigned int b = 0; b < N_breaks ; ++b) {
            //  mean_cdf_baseline[b] = Rcpp::mean( breaks[b] > counts_one_gene_permuted );
            //}
            // mean_cdf_baseline unaltered: it includes ALL cells, regardless of sample allocation.
            
            // loop over group g1 = 1:(n_groups-1)
            // loop over group g2 = (g1+1):n_groups
            // add all differences btw mean cdfs
            T_tmp = 0;
            for (unsigned int group_1 = 0; group_1 < n_groups; ++group_1) {
              mean_cdf_A = mean_cdf(_, group_1);
              T_tmp += Rcpp::sum(Rcpp::abs(mean_cdf_A - mean_cdf_baseline));
            }
            T_perm[p] = T_tmp;
          }
          
          // Rcout << "T_perm = " << std::endl << T_perm << std::endl;
          
          res(gene,cl_id) = ( P * Rcpp::mean(T_perm >= T_obs) + 1 )/(P+1);
          //  Rcpp::sum(T_perm >= T_obs)  OK, but /P gives 0...
        }
      }
    }
    return res;
  }
