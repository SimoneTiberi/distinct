#include <RcppArmadillo.h>
#include "Rfast.h"
#include <cmath>
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rfast)]]
using namespace arma;
using namespace Rcpp;
using namespace Rfast;

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
// inline int randWrapper(const int n) { return floor(unif_rand()*n); }

unsigned int perm_one_gene_covariates(arma::vec const& breaks, 
                                      arma::mat& cdf_A, 
                                      arma::mat& cdf_B,
                                      unsigned int const& n_samples,
                                      arma::uvec const& sample_has_cells,
                                      arma::uvec const& sample_ids_one_cl,
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

// Author of the function below: Stefanos Fafalios, taken from Rfast CRAN package
double calc_f(vec nix, double n, vec ni2hi2, double S, double x, int size){
  double sum1 = 0.0, sum2 = 0.0;
  
  for(int i = 0; i < size; i++){
    sum1+=log1p(nix[i]);
    sum2+=ni2hi2[i]/(1+nix[i]);
  }
  
  return sum1+n*log(S-x*sum2);
}

// Author of the function below: Stefanos Fafalios, taken from Rfast CRAN package
vec gold_rat3(double n, vec ni, vec ni2, double S, vec hi2,const int size, const double tol=1e-07){
  double a = 0, b = 50;
  const double ratio=0.618033988749895;
  double x1=b-ratio*b, x2=ratio*b;
  vec nix1 = ni*x1, nix2 = ni*x2, ni2hi2 = ni2%hi2;
  
  double f1 = calc_f(nix1, n, ni2hi2, S, x1, size);
  double f2 = calc_f(nix2, n, ni2hi2, S, x2, size);
  double bmina = b - a;
  while (abs(bmina)>tol){
    if(f2>f1){
      b=x2;
      bmina = b - a;
      x2=x1;
      f2=f1;
      x1=b - ratio * (bmina);
      nix1 = ni*x1;
      f1 = calc_f(nix1, n, ni2hi2, S, x1, size);
    }
    else {
      a=x1;
      bmina = b - a;
      x1=x2;
      f1=f2;
      x2=a + ratio * (bmina);
      nix2 = ni*x2;
      f2 = calc_f(nix2, n, ni2hi2, S, x2, size);
    }
  }
  vec ret(2);
  ret(0) = 0.5*(x1+x2);
  ret(1) = (f1+f2)/2;
  
  return ret;
}

// Author of the function below: Stefanos Fafalios, taken from Rfast CRAN package
vec my_rint_reg(arma::mat const& x, 
                arma::vec const& y,
                Rcpp::IntegerVector id, 
                unsigned int const& n,
                unsigned int const& p,
                double const& tol,
                unsigned int const& maxiters){
  // do these bits once only in every cluster of cells, not n_genes times!
  int idmx,idmn;
  maximum<int>(id.begin(),id.end(),idmx);
  minimum<int>(id.begin(),id.end(),idmn);
  mat xx(p,p),sx(idmx,p),sxy(p,1),mx(idmx,p);
  vec my(idmx);
  
  vec ni=Tabulate<vec,IntegerVector>(id,idmx);
  
  xx = cross_x<mat,mat>(x);
  for(unsigned int i=0;i<p;i++)
    sx.col(i) = group_sum_helper<vec,vec,IntegerVector>(x.col(i), id, &idmn,&idmx);
  sxy = cross_x_y<mat,mat,vec>(x,y);
  colvec sy = group_sum_helper<colvec,vec,IntegerVector>(y, id, &idmn,&idmx);
  mx = sx.each_col()/ni;
  my = sy/ni;
  
  mat b1 = solve(xx,sxy,solve_opts::fast);
  vec tmp = y - x*b1;
  double S = sum_with<square2<double>, vec>(tmp);
  vec tmp2 = my-mx*b1;
  vec hi2 = tmp2%tmp2;
  vec ni2 = ni%ni;
  
  vec d(2);
  d = gold_rat3(n, ni, ni2, S, hi2,idmx, tol);
  vec oneplnid = 1+ni*d(0);
  vec b2 = solve(xx - d(0)* cross_x_y<mat,mat,vec>(sx.each_col()/oneplnid, sx), sxy -
    d(0) * cross_x_y<mat,mat,vec>(sx, sy/oneplnid),solve_opts::fast);
  int i = 2;
  
  while(i++<maxiters && sum(abs(b2-b1.col(0))) > tol) {
    b1.col(0) = b2;
    
    tmp = y - x*b1;
    S = sum_with<square2<double>, vec>(tmp);
    tmp2 = my-mx*b1;
    hi2 = tmp2%tmp2;
    
    d = gold_rat3(n, ni, ni2, S, hi2,idmx, tol);
    oneplnid = 1+ni*d(0);
    b2 = solve(xx - d(0) * cross_x_y<mat,mat,vec>(sx.each_col()/oneplnid, sx), sxy -
      d(0) * cross_x_y<mat,mat,vec>(sx, sy/oneplnid),solve_opts::fast);
  }
  
  return b2;
}

arma::vec my_residuals(arma::mat const& X, 
                 arma::vec const& Y,
                 Rcpp::IntegerVector const& id, 
                 unsigned int const& n, 
                 unsigned int const& p, 
                 double const& tol,
                 unsigned int const& maxiters){
  vec coef = my_rint_reg( X, Y, id, n, p, tol, maxiters);
  
  vec res;
  res = Y - X * coef;
  return res;
} 

// [[Rcpp::export]]
List perm_test_covariates(unsigned int const& P,                             // number of permutations for all gene-cluster combinations
                          unsigned int const& P_2,                           // number of permutations when p < 0.1
                          unsigned int const& P_3,                           // number of permutations when p < 0.01
                          unsigned int const& P_4,                           // number of permutations when p < 0.001
                          unsigned int const& N_breaks,                      // number of breaks at which to evaluate the cdf
                          arma::vec const& cluster_ids,                            // ids of clusters (cell-population) for every cell
                          unsigned int const& n_clusters,                    // total number of clusters
                          arma::vec const& sample_ids,                             // ids of samples for every cell
                          unsigned int const& n_samples,                     // total number of samples
                          arma::vec const& group_ids_of_samples,                   // ids of groups (1 or 2) for every sample
                          double const& min_non_zero_cells,                  // min number of cells with > 0 expression in each cluster to consider the gene for testing
                          arma::sp_mat const& counts,                                // count matrix (rows = genes, cols = cells)
                          arma::mat& design) // design, excluding group
{
  arma::vec T_perm(P), T_perm_2(P_2 - P), T_perm_3(P_3 - P_2), T_perm_4(P_4 - P_3);
  
  /* Global variables */
  unsigned int tmp2, cell, b, sample, p, n_cells, k, K = design.n_cols;
  
  double T_obs, p_val;
  
  unsigned int n_genes = counts.n_cols, tmp_res; //, n_cells_overall = counts.n_rows;
  arma::uvec samples_in_group_A = arma::find( group_ids_of_samples == 1);
  arma::uvec samples_in_group_B = arma::find( group_ids_of_samples == 2);
  double s, mn, len_out;
  unsigned int n_samples_A = samples_in_group_A.n_elem;                        // number of samples in group A
  unsigned int n_samples_B = samples_in_group_B.n_elem;                        // number of samples in group B
  arma::mat res(n_genes,n_clusters); res.fill(-1);                             // array of p-values
  arma::vec seq0(N_breaks), breaks(N_breaks);
  for (b=0 ; b<N_breaks ; b++) {
    seq0(b) = b + 1.0;
  }
  
  arma::vec cdf_A(N_breaks);
  arma::vec cdf_B(N_breaks); 
  
  /* Start of cluster loop - declare the variables whose size won't change */
  arma::uvec sample_has_cells(n_samples);
  for (unsigned int cl_id = 0; cl_id < n_clusters; ++cl_id) { // loop over clusters
    // select data belonging to both groups:
    arma::uvec ids_both = arma::find( cluster_ids == cl_id );
    unsigned int n_both = ids_both.n_elem;
    arma::vec counts_one_gene(n_both);
    
    // vector with data for one gene (all cells of all clusters):
    //arma::vec counts_gene(n_cells_overall);
    
    n_cells = ids_both.n_elem;
    arma::uvec sample_ids_one_cl(n_cells); 
    for (b=0 ; b<n_cells ; b++) {
      sample_ids_one_cl(b) = sample_ids( ids_both(b) );
    }
    
    Rcpp::IntegerVector sample_ids_one_cl_vec(n_cells);
    sample_ids_one_cl_vec = sample_ids_one_cl +1;
    
    /* permutation martix */ 
    arma::umat PERM(n_cells,P);
    Rcpp::IntegerMatrix PERMUTATIONS(n_cells,P);
    Rcpp::IntegerVector permutation = Rcpp::seq_len(n_cells)-1; // define the vector for each permutation, -1 because Rcpp indexes start from 0 !
    for (p=0 ; p<P ; p++) {
      PERM.col(p) = randperm(n_cells);
      // std::random_shuffle(permutation.begin(), permutation.end(), randWrapper);
      // PERMUTATIONS(_,p) = permutation;
      // PERM.col(p) = Rcpp::as<arma::vec>(permutation);
    }
    // PERM = Rcpp::as<arma::umat>(PERMUTATIONS);
    
    /* SECOND permutation martix */ 
    arma::umat PERM_2(n_cells,P_2-P);
    Rcpp::IntegerMatrix PERMUTATIONS_2(n_cells,P_2-P);
    for (p=0 ; p<(P_2-P) ; p++) {
      PERM_2.col(p) = randperm(n_cells);
      // std::random_shuffle(permutation.begin(), permutation.end(), randWrapper);
      // PERMUTATIONS_2(_,p) = permutation;
      // PERM.col(p) = Rcpp::as<arma::vec>(permutation);
    }
    // PERM_2 = Rcpp::as<arma::umat>(PERMUTATIONS_2);
    
    /* THIRD permutation martix */ 
    arma::umat PERM_3(n_cells,P_3-P_2);
    Rcpp::IntegerMatrix PERMUTATIONS_3(n_cells,P_3-P_2);
    for (p=0 ; p<(P_3-P_2) ; p++) {
      PERM_3.col(p) = randperm(n_cells);
      // std::random_shuffle(permutation.begin(), permutation.end(), randWrapper);
      // PERMUTATIONS_3(_,p) = permutation;
      // PERM.col(p) = Rcpp::as<arma::vec>(permutation);
    }
    // PERM_3 = Rcpp::as<arma::umat>(PERMUTATIONS_3);
    
    /* THIRD permutation martix */ 
    arma::umat PERM_4(n_cells,P_4-P_3);
    Rcpp::IntegerMatrix PERMUTATIONS_4(n_cells,P_4-P_3);
    for (p=0 ; p<(P_4-P_3) ; p++) {
      PERM_4.col(p) = randperm(n_cells);
      // std::random_shuffle(permutation.begin(), permutation.end(), randWrapper);
      // PERMUTATIONS_4(_,p) = permutation;
      // PERM.col(p) = Rcpp::as<arma::vec>(permutation);
    }
    // PERM_4 = Rcpp::as<arma::umat>(PERMUTATIONS_4);
    
    sample_has_cells.fill(0);
    for (sample=0 ; sample<n_samples ; sample++) {
      sample_has_cells(sample) = any(sample_ids_one_cl==sample);
    }
    
    // only when NOT all samples have cells, rescale the id:
    unsigned int id = 1;
    bool cond = all(sample_has_cells);
    if( cond == false ){
      for (sample=0 ; sample<n_samples ; sample++) {
        if ( (sample_has_cells(sample))==1 ) { // check that samples have at least 1 cell
          for (unsigned int cell = 0; cell < n_cells; ++cell) {
            if(sample_ids(cell) == sample){
              sample_ids_one_cl_vec(cell) = id;
            }
          }
          id += 1;
        }
      }
    }
    
    // design has k cols (number of covariates) and n_sample rows
    // create an X_design matrix with k cols and N_cells rows
    arma::mat X_design(n_cells, K);
    // compute design matrix: select for every cell the corresponding sample's design row via sample_ids_one_cl[cell]
    for (unsigned int cell = 0; cell < n_cells; ++cell) {
      for(k = 0; k < K; k++){
        X_design(cell, k) = design(sample_ids_one_cl(cell), k);
      }
    }
    
    /* Flag to interrupt the code from R: */
    Rcpp::checkUserInterrupt();
    
    /* gene loop */ 
    for (unsigned int gene = 0; gene < n_genes; ++gene) {
      //Rcpp::checkUserInterrupt();
      
      /* Determine if there are sufficient counts */
      //counts_gene = counts.col(gene);
      arma::sp_mat counts_gene(counts.col(gene));
      
      s = 0.0;
      for (b=0 ; b<n_both ; b++) {
        counts_one_gene(b) = counts_gene( ids_both(b) );
        if ( counts_one_gene(b) > 0 ) {
          s += 1.0;
        }
      } 
      
      /* */
      if (s >= min_non_zero_cells) {
        // residuals:
        counts_one_gene = my_residuals(X_design, counts_one_gene, sample_ids_one_cl_vec, n_cells, K, 1e-08, 100);
        
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
            arma::uvec tmp1 = find( sample_ids_one_cl == sample ); // find the cells associated to the sample `sample`
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
        tmp_res = 1 + perm_one_gene_covariates(breaks, 
                                               cdf_A, 
                                               cdf_B,
                                               n_samples,
                                               sample_has_cells,
                                               sample_ids_one_cl,
                                               counts_one_gene,
                                               group_ids_of_samples,
                                               N_breaks,
                                               n_samples_A,
                                               n_samples_B,
                                               T_perm,
                                               P,
                                               n_both,
                                               PERM,
                                               T_obs);
        
        p_val = (1.0 * tmp_res)/( P + 1.0 );
        
        if( p_val <= 0.1 ){ // if p_val < 0.1, use 10 * perm
          if(P_2 > P){
            tmp_res += perm_one_gene_covariates(breaks, 
                                                cdf_A, 
                                                cdf_B,
                                                n_samples,
                                                sample_has_cells,
                                                sample_ids_one_cl,
                                                counts_one_gene,
                                                group_ids_of_samples,
                                                N_breaks,
                                                n_samples_A,
                                                n_samples_B,
                                                T_perm_2, // 2nd T_perm_2
                                                P_2 - P,  // 2nd P = (P_2 - P)
                                                n_both,
                                                PERM_2,    // 2nd PERM
                                                T_obs);
            
            p_val = (1.0 * tmp_res)/( P_2 + 1.0 );
          }
          
          if( p_val <= 0.01 ){ // if p_val < 0.01, use 100 * perm
            if(P_3 > P_2){
              tmp_res += perm_one_gene_covariates(breaks, 
                                                  cdf_A, 
                                                  cdf_B,
                                                  n_samples,
                                                  sample_has_cells,
                                                  sample_ids_one_cl,
                                                  counts_one_gene,
                                                  group_ids_of_samples,
                                                  N_breaks,
                                                  n_samples_A,
                                                  n_samples_B,
                                                  T_perm_3,  // 3rd T_perm_3
                                                  P_3 - P_2, // 3rd P = P_3-(P_2)
                                                  n_both,
                                                  PERM_3,    // 2nd PERM
                                                  T_obs);
              
              p_val = (1.0 * tmp_res)/( P_3 + 1.0 );
            }
            
            if( p_val <= 0.001 ){ // if p_val < 0.001, use 1,000 * perm
              if(P_4 > P_3){
                tmp_res += perm_one_gene_covariates(breaks, 
                                                    cdf_A, 
                                                    cdf_B,
                                                    n_samples,
                                                    sample_has_cells,
                                                    sample_ids_one_cl,
                                                    counts_one_gene,
                                                    group_ids_of_samples,
                                                    N_breaks,
                                                    n_samples_A,
                                                    n_samples_B,
                                                    T_perm_4,  // 3rd T_perm_3
                                                    P_4 - P_3, // 3rd P = P_3-(P_2)
                                                    n_both,
                                                    PERM_4,    // 2nd PERM
                                                    T_obs);
                
                p_val = (1.0 * tmp_res)/( P_4 + 1.0 );
              }
            }
          }
        }
        
        res(gene,cl_id) = p_val;
      }
    }       // End of loop for genes
  }         // End of loop for clusters
  
  List output = List::create( res );
  return output;
}
