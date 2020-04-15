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

// [[Rcpp::export]]
List perm_test(bool const logarithm,                            // if TRUE, use log2(counts + 1); if FALSE, use counts
                unsigned int const P,                             // number of permutations
                unsigned int const N_breaks,                      // number of breaks at which to evaluate the cdf
                arma::vec cluster_ids,                            // ids of clusters (cell-population) for every cell
                unsigned int const n_clusters,                    // total number of clusters
                arma::vec sample_ids,                             // ids of samples for every cell
                unsigned int const n_samples,                     // total number of samples
                arma::vec group_ids_of_samples,                   // ids of groups (1 or 2) for every sample
                double const min_non_zero_cells,                  // min number of cells with > 0 expression in each cluster to consider the gene for testing
                arma::vec group_ids,                              // ids of groups (cell-population) for every cell
                arma::mat& counts,                                // count matrix (rows = genes, cols = cells)
                unsigned int const nCores )                       // total number of cores used
{
  
  /* Set the number of cores used */
  //omp_set_num_threads(nCores);
  
  /* Global variables */
  unsigned int n_genes = counts.n_rows, n_cells_overall = counts.n_cols;
  arma::uvec samples_in_group_A = arma::find( group_ids_of_samples == 1);
  arma::uvec samples_in_group_B = arma::find( group_ids_of_samples == 2);
  double s, mn, len_out;
  unsigned int n_samples_A = samples_in_group_A.n_elem;                        // number of samples in group A
  unsigned int n_samples_B = samples_in_group_B.n_elem;                        // number of samples in group B
  arma::mat res(n_genes,n_clusters); res.fill(-1);                             // array of p-values
  //arma::mat chunk = core_load(n_genes,nCores);
  arma::vec seq0(N_breaks), breaks(N_breaks);
  for (unsigned int b=0 ; b<N_breaks ; b++) {
    seq0(b) = b + 1.0;
  }
  
  unsigned int tmp2, p_tmp;
  double tmp3, T_obs;
  arma::vec cdf_A(N_breaks);
  arma::vec cdf_B(N_breaks); 
  
  arma::vec T_perm(P);
  
  /* Start of cluster loop - declare the variables whose size won't change */
  unsigned int b, p, n_cells,  sample; //, core;
  arma::uvec sample_has_cells(n_samples);
  for (unsigned int cl_id = 0; cl_id < n_clusters; ++cl_id) { // loop over clusters
    // select data belonging to both groups:
    arma::uvec ids_both = arma::find( cluster_ids == cl_id );
    unsigned int n_both = ids_both.n_elem;
    arma::vec counts_one_gene(n_both);
    arma::vec counts_one_gene_permuted(n_both);
    
    // vector with data for one gene (all cells of all clusters):
    arma::vec counts_gene(n_cells_overall);
    
    n_cells = ids_both.n_elem;
    arma::uvec sample_ids_one_cl(n_cells); 
    for (b=0 ; b<n_cells ; b++) {
      sample_ids_one_cl(b) = sample_ids( ids_both(b) );
    }
    
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
    
    sample_has_cells.fill(0);
    for (sample=0 ; sample<n_samples ; sample++) {
      sample_has_cells(sample) = any(sample_ids_one_cl==sample);
    }
    
    /* Parallel code, blocked for now */ 
    //    # pragma omp parallel for shared( cl_id, ids_A, ids_B, ids_both, chunk, counts, n_A, n_B, n_both, seq0, sample_has_cells, n_samples_A, n_samples_B, sample_ids_one_cl, PERM, res ) private(b, p, sample)
    //    for (core=0 ; core<nCores ; core++) {
    //      unsigned int start_gene = chunk(core,0);
    //      unsigned int finish_gene = chunk(core,1);
    //      unsigned int gene;
    //      for (gene=start_gene ; gene<=finish_gene ; gene++) {
    
    /* Flag to interrupt the code from R: */
    // Rcpp::checkUserInterrupt();
    
    /* gene loop */ 
    for (unsigned int gene = 0; gene < n_genes; ++gene) {
      /* Determine if there are sufficient counts */
      counts_gene = counts.row(gene).t();
      
      s = 0.0;
      for (b=0 ; b<n_both ; b++) {
        counts_one_gene(b) = counts_gene( ids_both(b) );
        if ( counts_one_gene(b) > 0 ) {
          s += 1.0;
        }
      } 
      
      /* */
      if (s >= min_non_zero_cells) {
        //res(gene,cl_id) = 0.0;
        res(gene,cl_id) = 1.0;
        // if keep_gene is false (0), res(gene,cl_id) will be -1 (filled above).
        // if keep_gene is true (1),res(gene,cl_id) will be set at 1 (the observed value) and values added below (and then devided by P+1).
        
        if (logarithm==1) {
          b = counts_one_gene.n_elem;
          for (p=0 ; p<b ; p++) {
            counts_one_gene(p) = log2( counts_one_gene(p) + 1.0);
          }
        }
        
        /* CDF grid */ 
        mn      = counts_one_gene.min();
        len_out = ( counts_one_gene.max() - mn ) / (N_breaks + 1.0);
        for (b=0 ; b<N_breaks ; b++) {
          breaks(b) = mn + len_out*seq0(b) ;
        }
        
        cdf_A.fill(0.0);
        for (sample=0 ; sample<n_samples_A ; sample++ ) {
          if ( (sample_has_cells(sample))==1 ) {
            arma::uvec tmp1 = find( sample_ids_one_cl == sample );
            tmp2 = tmp1.n_elem;
            arma::vec counts_one_gene_one_sample(tmp2);
            for (b=0 ; b<tmp2 ; b++) {
              counts_one_gene_one_sample(b) = counts_one_gene( tmp1(b) );
            }
            for (b=0 ; b<N_breaks ; b++) {
              tmp3 =0.0; 
              for (p=0 ; p<tmp2 ; p++) {
                if (breaks(b)>counts_one_gene_one_sample(p)) {
                  tmp3 += 1.0;
                }
              }
              tmp3 /= tmp2*n_samples_A;
              cdf_A(b) += tmp3;
            }
          }
        }
        
        cdf_B.fill(0.0);
        for (sample = n_samples_A; sample < n_samples; sample++) {
          if ( (sample_has_cells(sample))==1 ) {
            arma::uvec tmp1 = find( sample_ids_one_cl == sample );
            tmp2 = tmp1.n_elem;
            arma::vec counts_one_gene_one_sample(tmp2);
            for (b=0 ; b<tmp2 ; b++) {
              counts_one_gene_one_sample(b) = counts_one_gene( tmp1(b) );
            }
            for (b=0 ; b<N_breaks ; b++) {
              tmp3 = 0.0;
              for (p=0 ; p<tmp2 ; p++) {
                if (breaks(b)>counts_one_gene_one_sample(p)) {
                  tmp3 += 1.0;
                }
              }
              tmp3/= tmp2*n_samples_B;
              cdf_B(b) += tmp3;
            }
          }
        }
        
        /* Obseved test statistic */
        T_obs = 0.0;
        for (b=0 ; b<N_breaks ; b++) {
          T_obs += fabs( cdf_A(b) - cdf_B(b) )  ;
        }
        
        
        /* Permutation testing */
        T_perm.fill(0.0);
        for (p=0 ; p<P ; p++) {
          
          for (b=0 ; b<n_both ; b++) {
            counts_one_gene_permuted(b) = counts_one_gene( PERM(b,p) );
          }
          
          cdf_A.fill(0.0);
          for (sample=0 ; sample<n_samples_A ; sample++ ) {
            if ( (sample_has_cells(sample))==1 ) {
              arma::uvec tmp1 = find( sample_ids_one_cl == sample );
              tmp2 = tmp1.n_elem;
              arma::vec counts_one_gene_one_sample(tmp2);
              for (b=0 ; b<tmp2 ; b++) {
                counts_one_gene_one_sample(b) = counts_one_gene_permuted( tmp1(b) );
              }
              for (b=0 ; b<N_breaks ; b++) {
                tmp3 = 0.0;
                for (p_tmp=0 ; p_tmp<tmp2 ; p_tmp++) {
                  if (breaks(b)>counts_one_gene_one_sample(p_tmp)) {
                    tmp3 += 1.0;
                  }
                }
                tmp3 /= tmp2*n_samples_A;
                cdf_A(b) += tmp3;
              }
            }
          }
          
          cdf_B.fill(0.0);
          for (sample = n_samples_A; sample < n_samples; sample++) {
            if ( (sample_has_cells(sample))==1 ) {
              arma::uvec tmp1 = find( sample_ids_one_cl == sample );
              tmp2 = tmp1.n_elem;
              arma::vec counts_one_gene_one_sample(tmp2);
              for (b=0 ; b<tmp2 ; b++) {
                counts_one_gene_one_sample(b) = counts_one_gene_permuted( tmp1(b) );
              }
              for (b=0 ; b<N_breaks ; b++) {
                tmp3 = 0.0;
                for (p_tmp=0 ; p_tmp<tmp2 ; p_tmp++) {
                  if (breaks(b)>counts_one_gene_one_sample(p_tmp)) {
                    tmp3 += 1.0;
                  }
                }
                tmp3 /= tmp2*n_samples_B;
                cdf_B(b) += tmp3;
              }
            }
          }
          
          /* Permuted test statistic */
          for (b=0 ; b<N_breaks ; b++) {
            T_perm(p) += fabs( cdf_A(b) - cdf_B(b) );
          }
          
        }   // End of loop for permutations
        for (p=0 ; p<P ; ++p) {
          if (T_perm(p)>=T_obs) {
            res(gene,cl_id) += 1.0;
          }
        }
        res(gene,cl_id) /= ( P + 1.0 );
      }     // End of if-keep-gene
    }       // End of loop for genes
    //    }         // End of loop for cores
  }           // End of loop for clusters
  
  List output = List::create( res );
  return output;
}
