// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
Rcpp::NumericMatrix
  perm_samples_test(bool const& logarithm, // if TRUE, use log2(counts + 1); if FALSE, use counts
            unsigned int const& P, // number of permutations
            unsigned int const& N_breaks, // number of breaks at which to evaluate the cdf
            Rcpp::IntegerVector const& cluster_ids, // ids of clusters (cell-population) for every cell
            unsigned int const& n_clusters, // total number of clusters
            Rcpp::IntegerVector const& sample_ids, // ids of samples for every cell
            unsigned int const& n_samples, // total number of samples
            Rcpp::IntegerVector const& group_ids_of_samples, // ids of groups (1 or 2) for every sample
            double const& min_counts_per_group, // min number of cells with > 0 expression in each group
            Rcpp::IntegerVector const& group_ids, // ids of groups (cell-population) for every cell
            Rcpp::NumericMatrix const& counts) // count matrix (rows = genes, cols = cells)
  {
    unsigned int n_genes = counts.nrow(), sample_p;
    double T_obs;
    
    Rcpp::NumericVector seq_len_N = as<NumericVector>(wrap(seq(1, N_breaks)));
    
    unsigned int n_cells; // number of cells PER group
    bool keep_gene; // logical , indicating whether to keep genes or not (if filtered due to low expression)
    Rcpp::IntegerVector sample_ids_one_cl, sel_A, sel_B, sel_both; // sel_A and sel_B define the indeces of the counts matrix referring to group A and B
    
    // create arma vectors, for group and cluster ids
    arma::vec group_ids_arma = Rcpp::as<arma::vec>(group_ids);
    arma::vec cluster_ids_arma = Rcpp::as<arma::vec>(cluster_ids);
    
    // turn sample_id into numeric:
    arma::vec group_ids_of_samples_arma = Rcpp::as<arma::vec>(group_ids_of_samples);
    arma::uvec samples_in_group_A_arma = arma::find( group_ids_of_samples_arma == 1);
    arma::uvec samples_in_group_B_arma = arma::find( group_ids_of_samples_arma == 2);
    Rcpp::IntegerVector samples_in_group_A = Rcpp::wrap(samples_in_group_A_arma);
    Rcpp::IntegerVector samples_in_group_B = Rcpp::wrap(samples_in_group_B_arma);
    
    // Logical vector to check if samples have at least 1 cell
    Rcpp::LogicalVector sample_has_cells(n_samples);
    
    unsigned int n_samples_A = samples_in_group_A.length(); // number of samples in group A
    unsigned int n_samples_B = samples_in_group_B.length(); // number of samples in group B
    // CHECK n_samples_A + n_samples_B == n_samples !!!
    arma::uvec ids_A, ids_B, ids_both;
    
    Rcpp::NumericVector mean_cdf_A(N_breaks), mean_cdf_B(N_breaks), T_perm(P), cdfs_one_break_A, cdfs_one_break_B, cdfs_one_break, counts_one_gene, counts_one_gene_A, counts_one_gene_B, breaks, counts_one_gene_one_sample;
    Rcpp::NumericMatrix CDF(n_samples, N_breaks), res(n_genes, n_clusters); // res: final p-values;
    std::fill(res.begin(), res.end(), -1);
    
    for (unsigned int cl_id = 0; cl_id < n_clusters; ++cl_id) { // loop over clusters
      // select data belonging to the first and second group:
      ids_A = arma::find( group_ids_arma == 1 && cluster_ids_arma == cl_id );
      ids_B = arma::find( group_ids_arma == 2 && cluster_ids_arma == cl_id );
      ids_both = arma::find( cluster_ids_arma == cl_id );
      
      sel_A = Rcpp::wrap(ids_A);
      sel_B = Rcpp::wrap(ids_B);
      sel_both = Rcpp::wrap(ids_both);
      
      sample_ids_one_cl = sample_ids[sel_both];
      
      // calculate the total number of cells in the cluster 'cl_id' (across both groups)
      n_cells = sel_both.length(); 
      
      // compute P permutations of the n_samples cells (sample without replacement).
      Rcpp::IntegerMatrix PERMUTATIONS(n_samples,P);
      Rcpp::IntegerVector permutation = Rcpp::seq_len(n_samples)-1; // define the vector for each permutation, -1 because Rcpp indexes start from 0 !
      for (unsigned int p = 0; p < P; ++p) {
        std::random_shuffle(permutation.begin(), permutation.end(), randWrapper);
        PERMUTATIONS(_,p) = permutation;
        // using Rcpp samplers will allow to control the seed from R.
      }
      
      // CHECK that every sample has at least 1 cell, independent of the gene (to be checked 1 for every cluster):
      for (unsigned int sample = 0; sample < n_samples; ++sample) {
        sample_has_cells[sample] = is_true(Rcpp::any(sample_ids_one_cl == sample));
      }
      // SPEED UP: COMPUTE "sample_ids_one_cl == sample" ONCE FOR EVERY CLUSTER ONLY HERE!!!
      
      // Rcpp::sugar provides a comparator with one value for vectors but not for matrices
      for (unsigned int gene = 0; gene < n_genes; ++gene) {
        counts_one_gene = counts(gene,_); // SELECT cl_id too!!!
        counts_one_gene_A = counts_one_gene[sel_A];
        counts_one_gene_B = counts_one_gene[sel_B];
        
        keep_gene = (sum(counts_one_gene_A > 0) >= min_counts_per_group) & (sum(counts_one_gene_B > 0) >= min_counts_per_group);
        // keep_gene is a filter to remove lowly abundant genes: in each group, we want at least `min_counts_per_group` cells to have non-zero expression.
        
        if(keep_gene){
          counts_one_gene = counts_one_gene[sel_both];
          
          if(logarithm){ // log data if `logarithm` is true
            for (unsigned int cell = 0; cell < counts_one_gene.length(); ++cell) {
              counts_one_gene[cell] = log2(counts_one_gene[cell] + 1);
            }
          }
          
          // compute the grid at which we ought to calculate the cdf:
          double min_ = min(counts_one_gene);
          double len_out = ( max(counts_one_gene) - min_ )/(N_breaks+1);
          breaks = min_ + len_out * seq_len_N; 
          
          // compute the cdf of each sample:
          for (unsigned int sample = 0; sample < n_samples; ++sample) {
            // if no sample_ids_one_cl == sample ???
            // put a contraint here!
            if(sample_has_cells[sample]){
              counts_one_gene_one_sample = counts_one_gene[sample_ids_one_cl == sample];
              // compute the cdf of each sample in counts_one_gene_list[sample]:
              for (unsigned int b = 0; b < N_breaks ; ++b) {
                CDF(sample, b) = Rcpp::mean( breaks[b] > counts_one_gene_one_sample );
                // put counts_one_gene[sample_ids_one_cl == sample] instead of counts_one_gene_one_sample
              }
            }
          }
          
          // compute the MEAN cdf for group A:
          std::fill(mean_cdf_A.begin(), mean_cdf_A.end(), 0);
          for (unsigned int sample = 0; sample < n_samples_A; ++sample) {
            // if no sample_ids_one_cl == sample ???
            // put a contraint here!
            if(sample_has_cells[sample]){
              // compute the cdf of each sample in counts_one_gene_list[sample]:
              for (unsigned int b = 0; b < N_breaks ; ++b) {
                mean_cdf_A[b] += CDF(sample, b)/n_samples_A;
                // put counts_one_gene[sample_ids_one_cl == sample] instead of counts_one_gene_one_sample
              }
            }
          }
          
          // compute the MEAN cdf for group B:
          std::fill(mean_cdf_B.begin(), mean_cdf_B.end(), 0);
          for (unsigned int sample = n_samples_A; sample < n_samples; ++sample) {
            // if no sample_ids_one_cl == sample ???
            // put a contraint here!
            if(sample_has_cells[sample]){
              // compute the cdf of each sample in counts_one_gene_list[sample]:
              for (unsigned int b = 0; b < N_breaks ; ++b) {
                mean_cdf_B[b] += CDF(sample, b)/n_samples_B;
                // put counts_one_gene[sample_ids_one_cl == sample] instead of counts_one_gene_one_sample
              }
            }
          }
          
          // test statistic for the observed data:
          T_obs = Rcpp::sum(Rcpp::abs(mean_cdf_A - mean_cdf_B));
          
          //Rcout << "A = " << std::endl << mean_cdf_A << std::endl;
          //Rcout << "B = " << std::endl << mean_cdf_B << std::endl;
          //Rcout << "T_obs = " << std::endl << T_obs << std::endl;
          
          //  test permuted data:
          for (unsigned int p = 0; p < P; ++p) {
            // obtain the permuted data:
            //counts_one_gene_permuted = counts_one_gene[PERMUTATIONS(_,p)];
            
            // split data, 1 per sample:
            std::fill(mean_cdf_A.begin(), mean_cdf_A.end(), 0);
            for (unsigned int sample = 0; sample < n_samples_A; ++sample) {
              // if no sample_ids_one_cl == sample ???
              // put a contraint here!
              if(sample_has_cells[sample]){
                counts_one_gene_one_sample = counts_one_gene[sample_ids_one_cl == PERMUTATIONS(sample,p) ];
                // compute the cdf of each sample in counts_one_gene_list[sample]:
                for (unsigned int b = 0; b < N_breaks ; ++b) {
                  mean_cdf_A[b] += Rcpp::mean( breaks[b] > counts_one_gene_one_sample )/n_samples_A;
                  // put counts_one_gene[sample_ids_one_cl == sample] instead of counts_one_gene_one_sample
                }
              }
            }
            
            std::fill(mean_cdf_B.begin(), mean_cdf_B.end(), 0);
            for (unsigned int sample = n_samples_A; sample < n_samples; ++sample) {
              // if no sample_ids_one_cl == sample ???
              // put a contraint here!
              if(sample_has_cells[sample]){
                counts_one_gene_one_sample = counts_one_gene[sample_ids_one_cl == PERMUTATIONS(sample,p) ];
                // compute the cdf of each sample in counts_one_gene_list[sample]:
                for (unsigned int b = 0; b < N_breaks ; ++b) {
                  mean_cdf_B[b] += Rcpp::mean( breaks[b] > counts_one_gene_one_sample )/n_samples_B;
                  // put counts_one_gene[sample_ids_one_cl == sample] instead of counts_one_gene_one_sample
                }
              }
            }
            
            // compute the MEAN cdf for group A:
            std::fill(mean_cdf_A.begin(), mean_cdf_A.end(), 0);
            for (unsigned int sample = 0; sample < n_samples_A; ++sample) {
              // take index of the permuted sample:
              sample_p = PERMUTATIONS(sample,p);
              if(sample_has_cells[sample_p]){
                for (unsigned int b = 0; b < N_breaks ; ++b) {
                  mean_cdf_A[b] += CDF(sample_p, b)/n_samples_A;
                  // put counts_one_gene[sample_ids_one_cl == sample] instead of counts_one_gene_one_sample
                }
              }
            }
            
            // compute the MEAN cdf for group B:
            std::fill(mean_cdf_B.begin(), mean_cdf_B.end(), 0);
            for (unsigned int sample = n_samples_A; sample < n_samples; ++sample) {
              // take index of the permuted sample:
              sample_p = PERMUTATIONS(sample,p);
              if(sample_has_cells[sample_p]){
                for (unsigned int b = 0; b < N_breaks ; ++b) {
                  mean_cdf_B[b] += CDF(sample_p, b)/n_samples_B;
                  // put counts_one_gene[sample_ids_one_cl == sample] instead of counts_one_gene_one_sample
                }
              }
            }
            
            // test statistic for the permuted data:
            T_perm[p] = Rcpp::sum(Rcpp::abs(mean_cdf_A - mean_cdf_B));
          }
          
          res(gene,cl_id) = ( P * Rcpp::mean(T_perm >= T_obs) + 1 )/(P+1);
          //  Rcpp::sum(T_perm >= T_obs)  OK, but /P gives 0...
        }
      }
    }
    return res;
  }
