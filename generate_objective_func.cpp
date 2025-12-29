#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;


// // [[Rcpp::export]]
// Rcpp::List generate_objective_func(const int str_num, const vec& str_indicator, const vec& A, 
//                       const vec& pai, const mat& hX_1, const mat& hX_0) {
//   
//   int n = pai.n_elem;
//   int J = hX_1.n_cols;
//   mat mat_A = zeros<mat>(2 * str_num * J, n);
//   mat mask = zeros<mat>(n, J);
//   vec B = zeros<vec>(2 * str_num * J);
//   
//   mat transformed_hX_0 = hX_0.each_col() % ((1-A) / (1-pai));  
//   mat transformed_hX_1 = hX_1.each_col() % (A / pai);
// 
//   for(int k=1;k<str_num+1;k++){
//     uvec row_indeces = find(str_indicator == k);
//     int str_size = row_indeces.n_elem;
//     if (str_size == 0) {continue;}
//     mask.rows(row_indeces) = ones(str_size, J);
//     
//     mat masked_transformed_hX_0 = mask % transformed_hX_0;
//     mat masked_transformed_hX_1 = mask % transformed_hX_1;
//     mat_A.rows((k-1) * J, k * J-1) = masked_transformed_hX_0.t();
//     mat_A.rows((k-1) * J + str_num * J, k * J - 1 + str_num * J) = masked_transformed_hX_1.t();
//     B.subvec((k-1) * J, k * J-1) = mean(mask % hX_0, 0).t();
//     B.subvec((k-1) * J + str_num * J, k * J - 1 + str_num * J) = mean(mask % hX_1, 0).t();
//     mask = zeros<mat>(n, J);
//   }
//   
//   // mat_A.row(2 * str_num * J) = ones<vec>(n).t();
//   // B(2 * str_num * J) = 1;
//   return Rcpp::List::create(
//     Rcpp::Named("mat_A") = mat_A,
//     Rcpp::Named("B") = B
//   );
//     
// }


// 
// // [[Rcpp::export]]
// Rcpp::List generate_objective_func(const int str_num, const vec& str_indicator, const vec& A, 
//                                    const vec& pai, const mat& hX_1, const mat& hX_0) {
//   
//   int n = pai.n_elem;
//   int J = hX_1.n_cols;
//   mat mat_A = zeros<mat>(2 * str_num * J, n);
//   mat mask = zeros<mat>(n, J);
//   vec B = zeros<vec>(2 * str_num * J);
//   
//   mat transformed_hX_0 = hX_0.each_col() % ((1-A) / (1-pai));  
//   mat transformed_hX_1 = hX_1.each_col() % (A / pai);
//   
//   for(int k=1;k<str_num+1;k++){
//     uvec row_indeces = find(str_indicator == k);
//     int str_size = row_indeces.n_elem;
//     if (str_size == 0) {continue;}
//     mask.rows(row_indeces) = ones(str_size, J);
//     
//     mat masked_transformed_hX_0 = mask % transformed_hX_0;
//     mat masked_transformed_hX_1 = mask % transformed_hX_1;
//     mat_A.rows((k-1) * J, k * J-1) = masked_transformed_hX_0.t();
//     mat_A.rows((k-1) * J + str_num * J, k * J - 1 + str_num * J) = masked_transformed_hX_1.t();
//     B.subvec((k-1) * J, k * J-1) = mean(mask % hX_0, 0).t();
//     B.subvec((k-1) * J + str_num * J, k * J - 1 + str_num * J) = mean(mask % hX_1, 0).t();
//     mask = zeros<mat>(n, J);
//   }
//   
//   mat abs_mat_A = abs(mat_A);
//   
//   vec row_max = max(abs_mat_A, 1);
//   
//   
//   mat mat_A_normalized = mat_A.each_col() / row_max;
//   B = B / row_max;
//   
//   
//   // mat_A_normalized.row(2 * str_num * J) = ones<vec>(n).t();
//   // B(2 * str_num * J) = 1;
//   return Rcpp::List::create(
//     // Rcpp::Named("mat_A") = mat_A_normalized,
//     // Rcpp::Named("B") = ones<vec>(2 * str_num * J)
//     Rcpp::Named("mat_A") = mat_A_normalized,
//     Rcpp::Named("B") = B
//   );
//   
// }






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





// [[Rcpp::export]]
Rcpp::List generate_objective_func_EL(const int str_num, const vec& str_indicator, const vec& A, 
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
  mat mat_A_normalized = mat_A.each_col() / row_max;
  // mat_A_normalized.row(2 * str_num * J) = ones<vec>(n).t();
  // B(2 * str_num * J) = 1;
  
  return Rcpp::List::create(
    Rcpp::Named("mat_A") = mat_A_normalized,
    Rcpp::Named("B") = B
  );
  
}





// [[Rcpp::export]]
Rcpp::List generate_objective_func_non_standard(const int str_num, const vec& str_indicator, const vec& A, 
                                   const vec& pai, const mat& hX_1, const mat& hX_0) {
  
  int n = pai.n_elem;
  int J = hX_1.n_cols;
  mat mat_A = zeros<mat>(2 * str_num * J + 2, n);
  mat mask = zeros<mat>(n, J);
  vec B = zeros<vec>(2 * str_num * J + 2);
  
  mat transformed_hX_0 = hX_0.each_col() % (1-A);  
  mat transformed_hX_1 = hX_1.each_col() % A;
  
  for(int k=1;k<str_num+1;k++){
    uvec row_indeces = find(str_indicator == k);
    int str_size = row_indeces.n_elem;
    if (str_size == 0) {continue;}
    mask.rows(row_indeces) = ones(str_size, J);
    
    mat masked_transformed_hX_0 = mask % transformed_hX_0;
    mat masked_transformed_hX_1 = mask % transformed_hX_1;
    mat_A.rows((k-1) * J, k * J-1) = masked_transformed_hX_0.t();
    mat_A.rows((k-1) * J + str_num * J, k * J - 1 + str_num * J) = masked_transformed_hX_1.t();
    B.subvec((k-1) * J, k * J-1) = mean(mask % hX_0, 0).t();
    B.subvec((k-1) * J + str_num * J, k * J - 1 + str_num * J) = mean(mask % hX_1, 0).t();
    mask = zeros<mat>(n, J);
  }
  mat mat_A_normalized = mat_A.each_col() / B;
  
  
  mat_A_normalized.row(2 * str_num * J) = A.t();
  mat_A_normalized.row(2 * str_num * J) = (1-A).t();
  B(2 * str_num * J) = 1;
  return Rcpp::List::create(
    Rcpp::Named("mat_A") = mat_A_normalized,
    Rcpp::Named("B") = ones<vec>(2 * str_num * J + 2)
  );
  
}




// [[Rcpp::export]]
Rcpp::List generate_objective_func_agg(const vec& A, const vec& pai, const mat& hX_1, const mat& hX_0) {
  
  int n = pai.n_elem;

  mat mat_A = zeros<mat>(2, n);

  vec B = zeros<vec>(2);
  
  vec transformed_hX_0 = hX_0.col(0) % ((1-A) / (1-pai));  
  vec transformed_hX_1 = hX_1.col(0) % (A / pai);
  B(0) = mean(hX_0.col(0));
  B(1) = mean(hX_1.col(0));
  mat_A.row(0) = transformed_hX_0.t() / B(0);
  mat_A.row(1) = transformed_hX_1.t() / B(1);


  
  


  return Rcpp::List::create(
    Rcpp::Named("mat_A") = mat_A,
    Rcpp::Named("B") = ones<vec>(2)
  );
  
}