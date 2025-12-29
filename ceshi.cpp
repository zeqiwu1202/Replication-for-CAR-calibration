#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
// [[Rcpp::export]]
arma::mat multiplyRowsByVectorElements(arma::mat& mat, const arma::vec& vec) {
  if (mat.n_rows != vec.n_elem) {
    Rcpp::stop("Matrix row count must match vector length");
  }
  
  return mat.each_col() % vec;  
}
