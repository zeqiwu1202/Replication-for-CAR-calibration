#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;









// [[Rcpp::export]]
Rcpp::List generate_objective_func(const int str_num, const vec& str_indicator, const vec& A, 
                                      const vec& pai, const mat& hX_1, const mat& hX_0) {
  
  int n = pai.n_elem;
  int J = hX_1.n_cols;
  mat mat_A = zeros<mat>(2 * str_num * J, n);
  mat hX_0_demean = hX_0;
  mat hX_1_demean = hX_1;
  vec B = zeros<vec>(2 * str_num * J);
  rowvec col_means0 = mean(hX_0, 0);
  rowvec col_means1 = mean(hX_1, 0);
  mat mask = zeros<mat>(n, J);
  hX_0_demean = hX_0.each_row() - col_means0;
  hX_1_demean = hX_1.each_row() - col_means1;
  
  //rowvec mm = mean(hX_0_demean, 0);
  //std::cout<<mm<<std::endl;
  
  
  mat transformed_hX_0 = hX_0_demean.each_col() % (A-pai);  
  mat transformed_hX_1 = hX_1_demean.each_col() % (A-pai);
  
  for(int k=1;k<str_num+1;k++){
    
    
    mat_A.rows((k-1) * J, k * J-1) = transformed_hX_0.t();
    mat_A.rows((k-1) * J + str_num * J, k * J - 1 + str_num * J) = transformed_hX_1.t();
    
    
  }
  
  mat abs_mat_A = abs(mat_A);
  
  vec row_max = max(abs_mat_A, 1);
  row_max = clamp(row_max, 0.001, datum::inf);
  mat mat_A_normalized = mat_A.each_col() / row_max;
  // mat_A_normalized.row(2 * str_num * J) = ones<vec>(n).t();
  // B(2 * str_num * J) = 1;
  
  return Rcpp::List::create(
    Rcpp::Named("mat_A") = mat_A_normalized,
    Rcpp::Named("B") = B
  );
  
}









