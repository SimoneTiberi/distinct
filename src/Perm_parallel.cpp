#include <RcppArmadillo.h>
//#include <Rcpp.h>
#include <cmath>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

unsigned int perm_one_gene_parallel(arma::vec const& breaks, 
                           arma::mat& cdf_A, 
                           arma::mat& cdf_B,
                           unsigned int const& n_samples,
                           arma::uvec const& sample_has_cells,
                           arma::vec const& sample_ids_one_cl,
                           arma::vec const& counts_one_gene,
                           arma::vec const& group_ids_of_samples,
                           unsigned int const& N_breaks,
                           unsigned int const& n_samples_A,
                           unsigned int const& n_samples_B,
                           arma::vec& T_perm,
                           unsigned int const& P,
                           unsigned int const& n_both,
                           arma::umat const& PERM,
                           double T_obs)
{
  unsigned int res = 0;
  unsigned int tmp2, cell, b, sample, p;
  
  /* Permutation testing */
  T_perm.fill(0.0);
  arma::vec counts_one_gene_permuted(n_both);
  
  for (p=0 ; p<P ; p++) {
    
    for (b=0 ; b<n_both ; b++) {
      counts_one_gene_permuted(b) = counts_one_gene( PERM(b,p) );
    }
    
    cdf_A.fill(0.0);
    cdf_B.fill(0.0);
    for (sample=0 ; sample<n_samples ; sample++ ) {
      if ( (sample_has_cells(sample))==1 ) { // check that samples have at least 1 cell
        arma::uvec tmp1 = find( sample_ids_one_cl == sample ); // find the cells associated to the sample `sample`
        tmp2 = tmp1.n_elem; // counts number of cells
        arma::vec counts_one_gene_one_sample(tmp2); // create a vector to store these cells
        
        for (b=0 ; b<tmp2 ; b++) { // store cells of sample `sample`
          counts_one_gene_one_sample(b) = counts_one_gene_permuted( tmp1(b) );
        }
        
        // sort counts in ascending order
        std::sort(counts_one_gene_one_sample.begin(), counts_one_gene_one_sample.end());
        
        cell = 0; // cell index
        
        if ( (group_ids_of_samples(sample)) == 1.0) { // if sample belongs to group A, add to cdf_A
          
          for (b=0 ; b<N_breaks ; b++) { // loop over the number of breaks of the CDF
            // loop until cell index is smaller than tmp2 (length of counts_one_gene_one_sample)
            while(cell < tmp2){
              if (breaks(b)>counts_one_gene_one_sample(cell) ) { // if break is > cell, keep counting and go to next cell
                cell += 1;
              }
              else{ // if break is < cell, stop loop and save cdf for break b
                cdf_A(b) += 1.0 * cell/(tmp2*n_samples_A);
                break;
              }
            }
            
            if(cell == tmp2){
              cdf_A(b) += 1.0/n_samples_A;
            }
            
          } //end of b for loop
        } // end of GROUP 1
        else { // if sample belongs to group B, add to cdf_B
          
          for (b=0 ; b<N_breaks ; b++) { // loop over the number of breaks of the CDF
            // loop until cell index is smaller than tmp2 (length of counts_one_gene_one_sample)
            while(cell < tmp2){
              if (breaks(b)>counts_one_gene_one_sample(cell) ) { // if break is > cell, keep counting and go to next cell
                cell += 1;
              }
              else{ // if break is < cell, stop loop and save cdf for break b
                cdf_B(b) += 1.0 * cell/(tmp2*n_samples_B);
                break;
              }
            }
            
            if(cell == tmp2){
              cdf_B(b) += 1.0/n_samples_B;
            }
            
          } //end of b for loop
        } // end of GROUP 2
        
      }
    }
    
    /* Permuted test statistic */
    for (b=0 ; b<N_breaks ; b++) {
      T_perm(p) += fabs( cdf_A(b) - cdf_B(b) );
    }
    
  }   // End of loop for permutations
  
  for (p=0 ; p<P ; ++p) {
    if (T_perm(p)>=T_obs) {
      res += 1;
    }
  }
  
  return res;
}

// [[Rcpp::export]]
List perm_test_parallel(unsigned int const& P,                             // number of permutations for all gene-cluster combinations
                        unsigned int const& P_2,                           // number of permutations when p < 0.1
                        unsigned int const& P_3,                           // number of permutations when p < 0.01
                        unsigned int const& P_4,                           // number of permutations when p < 0.001
                        unsigned int const& N_breaks,                      // number of breaks at which to evaluate the cdf
                        arma::vec const& sample_ids,                             // ids of samples for every cell
                        unsigned int const& n_samples,                     // total number of samples
                        arma::vec const& group_ids_of_samples,                   // ids of groups (1 or 2) for every sample
                        double const& min_non_zero_cells,                  // min number of cells with > 0 expression in each cluster to consider the gene for testing
                        arma::sp_mat const& counts)
{
  /* Set the number of cores used */
  arma::vec T_perm(P), T_perm_2(P_2 - P), T_perm_3(P_3 - P_2), T_perm_4(P_4 - P_3);
  
  /* Global variables */
  unsigned int tmp2, cell, b, sample, p, n_cells;
  
  double T_obs, p_val;
  
  unsigned int n_genes = counts.n_cols, tmp_res;
  arma::uvec samples_in_group_A = arma::find( group_ids_of_samples == 1);
  arma::uvec samples_in_group_B = arma::find( group_ids_of_samples == 2);
  double s, mn, len_out;
  unsigned int n_samples_A = samples_in_group_A.n_elem;                        // number of samples in group A
  unsigned int n_samples_B = samples_in_group_B.n_elem;                        // number of samples in group B
  arma::vec res(n_genes); res.fill(-1);                             // array of p-values
  //arma::mat chunk = core_load(n_genes,nCores);
  arma::vec seq0(N_breaks), breaks(N_breaks);
  for (b=0 ; b<N_breaks ; b++) {
    seq0(b) = b + 1.0;
  }
  
  arma::vec cdf_A(N_breaks);
  arma::vec cdf_B(N_breaks); 
  
  /* Start of cluster loop - declare the variables whose size won't change */
  arma::uvec sample_has_cells(n_samples);
  // select data belonging to both groups:
  //arma::uvec ids_both = arma::find( cluster_ids == cl_id );
  n_cells = counts.n_rows; //ids_both.n_elem;
  arma::vec counts_one_gene(n_cells);
  
  // vector with data for one gene (all cells of all clusters):
  //arma::vec counts_gene(n_cells_overall);
  
  //arma::uvec sample_ids_one_cl(n_cells); 
  //for (b=0 ; b<n_cells ; b++) {
  //  sample_ids_one_cl(b) = sample_ids( ids_both(b) );
  //}
  
  /* permutation martix */ 
  arma::umat PERM(n_cells,P);
  Rcpp::IntegerMatrix PERMUTATIONS(n_cells,P);
  Rcpp::IntegerVector permutation = Rcpp::seq_len(n_cells)-1; // define the vector for each permutation, -1 because Rcpp indexes start from 0 !
  for (p=0 ; p<P ; p++) {
    // PERM.col(p) = randperm(n_cells);
    std::random_shuffle(permutation.begin(), permutation.end(), randWrapper);
    PERMUTATIONS(_,p) = permutation;
    // PERM.col(p) = Rcpp::as<arma::vec>(permutation);
  }
  PERM = Rcpp::as<arma::umat>(PERMUTATIONS);
  
  /* SECOND permutation martix */ 
  arma::umat PERM_2(n_cells,P_2-P);
  Rcpp::IntegerMatrix PERMUTATIONS_2(n_cells,P_2-P);
  for (p=0 ; p<(P_2-P) ; p++) {
    // PERM.col(p) = randperm(n_cells);
    std::random_shuffle(permutation.begin(), permutation.end(), randWrapper);
    PERMUTATIONS_2(_,p) = permutation;
    // PERM.col(p) = Rcpp::as<arma::vec>(permutation);
  }
  PERM_2 = Rcpp::as<arma::umat>(PERMUTATIONS_2);
  
  /* THIRD permutation martix */ 
  arma::umat PERM_3(n_cells,P_3-P_2);
  Rcpp::IntegerMatrix PERMUTATIONS_3(n_cells,P_3-P_2);
  for (p=0 ; p<(P_3-P_2) ; p++) {
    // PERM.col(p) = randperm(n_cells);
    std::random_shuffle(permutation.begin(), permutation.end(), randWrapper);
    PERMUTATIONS_3(_,p) = permutation;
    // PERM.col(p) = Rcpp::as<arma::vec>(permutation);
  }
  PERM_3 = Rcpp::as<arma::umat>(PERMUTATIONS_3);
  
  /* THIRD permutation martix */ 
  arma::umat PERM_4(n_cells,P_4-P_3);
  Rcpp::IntegerMatrix PERMUTATIONS_4(n_cells,P_4-P_3);
  for (p=0 ; p<(P_4-P_3) ; p++) {
    // PERM.col(p) = randperm(n_cells);
    std::random_shuffle(permutation.begin(), permutation.end(), randWrapper);
    PERMUTATIONS_4(_,p) = permutation;
    // PERM.col(p) = Rcpp::as<arma::vec>(permutation);
  }
  PERM_4 = Rcpp::as<arma::umat>(PERMUTATIONS_4);
  
  sample_has_cells.fill(0);
  for (sample=0 ; sample<n_samples ; sample++) {
    sample_has_cells(sample) = any(sample_ids==sample);
  }
  
  /* Flag to interrupt the code from R: */
  // Rcpp::checkUserInterrupt();
  
  /* gene loop */ 
  for (unsigned int gene = 0; gene < n_genes; ++gene) {
    /* Determine if there are sufficient counts */
    //counts_gene = counts.col(gene);
    // arma::sp_mat counts_gene(counts.col(gene));
    
    // TODO: turn counts_one_gene into arma::sp_mat  (currently arma::vec):
    counts_one_gene = counts.col(gene);
    
    s = sum(counts_one_gene > 0);
    
    //s = 0.0;
    //for (b=0 ; b<n_cells ; b++) {
    //  // counts_one_gene(b) = counts_gene( ids_both(b) );
    //  if ( counts_one_gene(b) > 0 ) {
    //    s += 1.0;
    //  }
    //} 
    
    /* */
    if (s >= min_non_zero_cells) {
      // if keep_gene is false (0), res(gene,cl_id) will be -1 (filled above).
      // if keep_gene is true (1),res(gene,cl_id) will be set at 1 (the observed value) and values added below (and then devided by P+1).
      
      /* CDF grid */ 
      mn = counts_one_gene.min();
      len_out = ( counts_one_gene.max() - mn ) / (N_breaks + 1.0);
      for (b=0 ; b<N_breaks ; b++) {
        breaks(b) = mn + len_out*seq0(b) ;
      }
      
      /* compute T_obs */ 
      cdf_A.fill(0.0);
      cdf_B.fill(0.0);
      for ( sample=0 ; sample<n_samples ; sample++ ) {
        if ( (sample_has_cells(sample))==1 ) { // check that samples have at least 1 cell
          arma::uvec tmp1 = find( sample_ids == sample ); // find the cells associated to the sample `sample`
          tmp2 = tmp1.n_elem; // counts number of cells
          arma::vec counts_one_gene_one_sample(tmp2); // create a vector to store these cells
          for ( b=0 ; b<tmp2 ; b++) { // store cells of sample `sample`
            counts_one_gene_one_sample(b) = counts_one_gene( tmp1(b) );
          }
          
          // sort counts in ascending order
          std::sort(counts_one_gene_one_sample.begin(), counts_one_gene_one_sample.end());
          
          cell = 0; // cell index
          
          if ( (group_ids_of_samples(sample)) == 1.0) { // if sample belongs to group A, add to cdf_A
            
            for ( b=0 ; b<N_breaks ; b++) { // loop over the number of breaks of the CDF
              // loop until cell index is smaller than tmp2 (length of counts_one_gene_one_sample)
              while(cell < tmp2){
                if (breaks(b)>counts_one_gene_one_sample(cell) ) { // if break is > cell, keep counting and go to next cell
                  cell += 1;
                }
                else{ // if break is < cell, stop loop and save cdf for break b
                  cdf_A(b) += 1.0 * cell/(tmp2*n_samples_A);
                  break;
                }
              }
              
              if(cell == tmp2){
                cdf_A(b) += 1.0/n_samples_A;
              }
              
            } //end of b for loop
          } // end of GROUP 1
          else { // if sample belongs to group B, add to cdf_B
            
            for ( b=0 ; b<N_breaks ; b++) { // loop over the number of breaks of the CDF
              // loop until cell index is smaller than tmp2 (length of counts_one_gene_one_sample)
              while(cell < tmp2){
                if (breaks(b)>counts_one_gene_one_sample(cell) ) { // if break is > cell, keep counting and go to next cell
                  cell += 1;
                }
                else{ // if break is < cell, stop loop and save cdf for break b
                  cdf_B(b) += 1.0 * cell/(tmp2*n_samples_B);
                  break;
                }
              }
              
              if(cell == tmp2){
                cdf_B(b) += 1.0/n_samples_B;
              }
              
            } //end of b for loop
          } // end of GROUP 2
          
        }
      }
      
      /* Obseved test statistic */
      T_obs = 0.0;
      for ( b=0 ; b<N_breaks ; b++) {
        T_obs += fabs( cdf_A(b) - cdf_B(b) )  ;
      }
      
      // run 100 permutations on all gene-cluster combinations
      tmp_res = 1 + perm_one_gene_parallel(breaks, 
                                  cdf_A, 
                                  cdf_B,
                                  n_samples,
                                  sample_has_cells,
                                  sample_ids,
                                  counts_one_gene,
                                  group_ids_of_samples,
                                  N_breaks,
                                  n_samples_A,
                                  n_samples_B,
                                  T_perm,
                                  P,
                                  n_cells,
                                  PERM,
                                  T_obs);
      
      p_val = (1.0 * tmp_res)/( P + 1.0 );
      
      if( p_val <= 0.1 ){ // if p_val < 0.1, use 10 * perm
        if(P_2 > P){
          tmp_res += perm_one_gene_parallel(breaks, 
                                   cdf_A, 
                                   cdf_B,
                                   n_samples,
                                   sample_has_cells,
                                   sample_ids,
                                   counts_one_gene,
                                   group_ids_of_samples,
                                   N_breaks,
                                   n_samples_A,
                                   n_samples_B,
                                   T_perm_2, // 2nd T_perm_2
                                   P_2 - P,  // 2nd P = (P_2 - P)
                                   n_cells,
                                   PERM_2,    // 2nd PERM
                                   T_obs);
          
          p_val = (1.0 * tmp_res)/( P_2 + 1.0 );
        }
        
        if( p_val <= 0.01 ){ // if p_val < 0.01, use 100 * perm
          if(P_3 > P_2){
            tmp_res += perm_one_gene_parallel(breaks, 
                                     cdf_A, 
                                     cdf_B,
                                     n_samples,
                                     sample_has_cells,
                                     sample_ids,
                                     counts_one_gene,
                                     group_ids_of_samples,
                                     N_breaks,
                                     n_samples_A,
                                     n_samples_B,
                                     T_perm_3,  // 3rd T_perm_3
                                     P_3 - P_2, // 3rd P = P_3-(P_2)
                                     n_cells,
                                     PERM_3,    // 2nd PERM
                                     T_obs);
            
            p_val = (1.0 * tmp_res)/( P_3 + 1.0 );
          }
          
          if( p_val <= 0.001 ){ // if p_val < 0.001, use 1,000 * perm
            if(P_4 > P_3){
              tmp_res += perm_one_gene_parallel(breaks, 
                                       cdf_A, 
                                       cdf_B,
                                       n_samples,
                                       sample_has_cells,
                                       sample_ids,
                                       counts_one_gene,
                                       group_ids_of_samples,
                                       N_breaks,
                                       n_samples_A,
                                       n_samples_B,
                                       T_perm_4,  // 3rd T_perm_3
                                       P_4 - P_3, // 3rd P = P_3-(P_2)
                                       n_cells,
                                       PERM_4,    // 2nd PERM
                                       T_obs);
              
              p_val = (1.0 * tmp_res)/( P_4 + 1.0 );
            }
          }
        }
      }
      
      res(gene) = p_val;
    }
  }       // End of loop for genes
  
  List output = List::create( res );
  return output;
}
