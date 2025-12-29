#include "RcppArmadillo.h"
#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <cmath>
#include <iostream>
#include <typeinfo>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

Rcpp::NumericVector Csample(int n, int num,
                            bool replace,
                            arma::vec proba) {
  Rcpp::NumericVector a(n);
  for(int i =0;i<n;i++){
    a(i) = i + 1;
  }
  Rcpp::NumericVector ret = Rcpp::RcppArmadillo::sample(a, num, replace, proba);
  return ret;
}

bool same(arma::vec a, arma::vec b){
  bool re = TRUE;
  if(a.n_elem == b.n_elem){
    for(unsigned int i = 0; i < a.n_elem; i++){
      if(a[i]!=b[i]){
        re = FALSE;break;
      }
    }
  }
  else{re = FALSE;}
  return re;
}

//[[Rcpp::export]]
arma::mat PStrR(arma::mat M){
  unsigned n = M.n_cols;
  arma::vec ind(n, arma::fill::zeros);
  unsigned int i, j;
  for(i = 0; i < n; i++){
    for(j = i + 1; j < n; j++){
      if(same(M.col(j), M.col(i))){
        ind(i) = 1; break;
      }
    }
  }
  return M.cols(find(ind == 0));
}

arma::uvec ReturnCol(arma::mat M, arma::vec V){
  unsigned int leng = V.n_elem;
  arma::uvec Ind = find(M.row(0) == V(0));
  for(unsigned int i = 1; i < leng; i++){
    arma::uvec a(1); a(0) = i;
    arma::uvec ind = find(M.submat(a, Ind) == V(i));
    Ind = Ind(ind);
  }
  return Ind + 1;
}

//[[Rcpp::export]]
arma::vec sqchoose(std::string method, double pi, int strt_num){
  arma::vec sq(strt_num,arma::fill::zeros);
  if(method == "SRS"){
    sq.fill(pi*(1.0 - pi));
  }
  else if(method == "WEI"){
    sq.fill(1.0/12);
  }
  return sq;
}

arma::mat HHunequal(arma::mat P, arma::mat D, arma::vec cov_profile, unsigned int cov_num,
                    arma::vec level_num, arma::vec omega, double p = 0.85, double pi = 0.5){
  int strt_num = P.n_cols; 
  
  arma::vec brid(2); brid(0) = 1.0/pi; brid(1) = -1.0/(1.0-pi);
  arma::uvec rvec = ReturnCol(P, cov_profile);
  int r = rvec(0);
  arma::uvec sub(2 + cov_num); sub(0) = 0; sub(1) = r;
  int temp = 0;
  for(unsigned int j = 0; j < cov_num; j++){
    sub(j + 2) = strt_num + temp + cov_profile(j);
    temp += level_num(j);
  }
  arma::uvec a(1); a(0) = 0;
  arma::vec arg = trans(omega) * D.submat(sub, a);
  double argval = arg(0, 0);
  
  arma::vec T_One(1); 
  if(argval > -0.000001 && argval < 0.000001){
    T_One = arma::randu<arma::vec>(1); 
    D.submat(sub, a) = D.submat(sub, a) + brid(sum(T_One > 0.5));
    D(D.n_rows - 1, 0) = -sum(T_One > 0.5) + 1; 
    return D;
  }
  if(argval >= 0.000001){
    T_One = arma::randu<arma::vec>(1);
    D.submat(sub, a) = D.submat(sub, a) + brid(sum(T_One > 1 - p));
    D(D.n_rows - 1, 0) = -sum(T_One > 1 - p)+ 1; 
    return D;
  }
  else{
    T_One = arma::randu<arma::vec>(1);
    D.submat(sub, a) = D.submat(sub, a) + brid(sum(T_One > p));
    D(D.n_rows - 1, 0) = -sum(T_One > p) + 1;
    return D;
  }
}

//[[Rcpp::export]]
arma::vec SRS(int n, double pi){
  arma::vec proba(2);
  proba(0) = 1 - pi;
  proba(1) = pi;
  Rcpp::NumericVector a = Csample(2,n,TRUE,proba);
  arma::vec s = Rcpp::as<arma::vec>(a);
  return s - 1;
}

//[[Rcpp::export]]
arma::vec BCD(arma::mat profiles, double pi, double lambda){
  int cov_num = profiles.n_rows;
  arma::vec omeganew(2 + cov_num, arma::fill::zeros); omeganew(1) = 1;
  int n = profiles.n_cols;
  arma::vec level_num = max(profiles, 1);
  arma::mat P = PStrR(profiles);
  int strt_num = P.n_cols;
  arma::vec D(2 + strt_num + sum(level_num),arma::fill::zeros);
  arma::vec Assign(n);
  for(int i = 0; i< n; i++){
    D = HHunequal(P, D, profiles.col(i), cov_num, level_num, omeganew, lambda, pi);
    Assign(i) = D(D.n_rows-1,0);
  }
  return Assign;
}

double AdaptFUN(double d, double nk, double pi){
  if(nk == 0 || d == 0){
    return pi;
  }
  else if(d > 0){
    return pi*(1.0 - d/nk);
  }
  else{
    return (pi-1.0)*d/nk+pi;
  }
}

arma::field<arma::vec> AdaptOne(arma::vec D,arma::vec nk,arma::mat PStr,arma::vec cov_profile,double pi){
  arma::field<arma::vec> result(3,1);
  arma::vec T_One = arma::randu<arma::vec>(1);
  int r = ReturnCol(PStr,cov_profile)(0);
  double p = AdaptFUN(D(r - 1),nk(r-1),pi);
  double T = sum(T_One < p);
  result(0,0) = T;
  arma::vec brid(2);
  brid(1) = 1.0/pi;
  brid(0) = -1.0/(1.0-pi);
  D(r - 1) = D(r - 1) + brid(T);
  nk(r-1) = nk(r-1) + abs(brid(T));
  result(1,0) = D;
  result(2,0) = nk;
  return result;
}

//[[Rcpp::export]]
arma::vec WEI(arma::mat profiles,double pi){
  arma::mat P = PStrR(profiles);
  int strt_num = P.n_cols;
  arma::vec D(strt_num,arma::fill::zeros);
  arma::vec nk(strt_num,arma::fill::zeros);
  double n = profiles.n_cols;
  arma::field<arma::vec> result(2,1);
  arma::vec Assign(n);
  for(int i = 0; i < n; i++){
    result = AdaptOne(D,nk,P,profiles.col(i),pi);
    Assign(i) = result(0,0)(0);
    D = result(1,0);
    nk = result(2,0);
  }
  return Assign;
}

arma::mat permutations(int n){
  arma::mat Z(n, arma::prod(arma::linspace<arma::vec>(1, n, n))); 
  Z(0, 0) = 1;
  for(int i = 1; i < n; i++){
    int nip = arma::prod(arma::linspace<arma::vec>(1, i, i)); 
    Z.cols(0, nip - 1).row(i).fill(i + 1);
    arma::mat xtemp = Z.cols(0, nip - 1).rows(0, i);
    arma::uvec ind(2 * i + 1); 
    ind.subvec(0, i) = arma::linspace<arma::uvec>(0, i, i + 1);
    ind.subvec(i + 1, i + i) = arma::linspace<arma::uvec>(0, i - 1, i); 
    for(int j = 0; j < i; j++){
      Z.cols((j + 1) * nip, (j + 2) * nip - 1).rows(0, i) = xtemp.rows(ind.subvec(j + 1, j + 1 + i));
    }
  }
  return Z;
}

// [[Rcpp::export]]
arma::mat Bpert(int bsize, double pi){
  arma::mat BP = permutations(bsize);
  int i = 0;
  for(i = 1; i < (floor(bsize*(1.0-pi))+1); i++){
    BP.replace(i,1);
  }
  for(i = (floor(bsize*(1.0-pi))+1); i < (bsize + 1); i++){
    BP.replace(i,2);
  }
  arma::mat BPT = PStrR(BP);
  return BPT;
}

// [[Rcpp::export]]
arma::field<arma::mat> StrROne(arma::mat D, arma::mat PS, arma::vec cov_profile,
                               unsigned int cov_num, arma::vec level_num, int bsize, 
                               arma::mat B, arma::mat BG, arma::vec strp, double pi){
  arma::field<arma::mat> Res(4, 1); 
  
  arma::vec brid(2); brid(0) = 1.0/pi; brid(1) = -1.0/(1.0-pi);
  int Psize = PS.n_cols, Bsize = B.n_cols; 
  
  int r = ReturnCol(PS, cov_profile)(0);
  arma::uvec sub(2 + cov_num); sub(0) = 0; sub(1) = r;
  int temp = 0; 
  for(unsigned int j = 0; j < cov_num; j++){
    sub(2 + j) = Psize + temp + cov_profile(j);
    temp += level_num(j);
  }
  strp(r - 1) = strp(r - 1) + 1;
  arma::uvec ustrp = arma::conv_to<arma::uvec>::from(strp); 
  int pos = ustrp(r - 1) % bsize; 
  int T_STR;
  if(pos == 0){
    T_STR = BG(bsize - 1, r - 1); 
    BG.col(r - 1) = B.col(arma::randi<arma::vec>(Bsize, arma::distr_param(1, Bsize))(0) - 1);
  }
  else{
    T_STR = BG(pos - 1, r - 1);
  }
  D.rows(sub) = D.rows(sub) + brid(T_STR - 1);
  
  Res(0, 0) = strp; 
  Res(1, 0) = BG; 
  Res(2, 0) = T_STR; 
  Res(3, 0) = D; 
  return Res;
}

//[[Rcpp::export]]
arma::vec SBR(arma::mat profiles, double pi, int bsize = 6){
  unsigned int cov_num = profiles.n_rows;
  arma::vec level_num = max(profiles,1);
  arma::mat P = PStrR(profiles);
  int strt_num = P.n_cols;
  int n = profiles.n_cols;
  arma::vec assignew(n,arma::fill::zeros);
  arma::vec D(1 + strt_num + sum(level_num),arma::fill::zeros);
  arma::mat B = Bpert(bsize,pi);
  int Bsize = B.n_cols;
  arma::uvec shuffle = arma::randi<arma::uvec>(n,arma::distr_param(0,Bsize-1));
  arma::mat BG = B.cols(shuffle);
  arma::vec strp(strt_num,arma::fill::zeros);
  arma::field<arma::mat> Res(4,1);
  for(int i = 0; i < n; i++){
    Res = StrROne(D,P,profiles.col(i),cov_num,level_num,bsize,B,BG,strp,pi);
    strp = Res(0,0);
    BG = Res(1,0);
    assignew(i) = Res(2,0)(0,0);
    D = Res(3,0);
  }
  return assignew - 1;
}

//[[Rcpp::export]]
arma::vec PocSimue(arma::mat profiles, double pi, arma::vec weight, double p = 0.85){
  unsigned int cov_num = profiles.n_rows;
  arma::vec level_num = max(profiles,1);
  arma::mat P = PStrR(profiles);
  int strt_num = P.n_cols;
  int n = profiles.n_cols;
  arma::vec Assign(n);
  arma::vec D(2+strt_num+sum(level_num),arma::fill::zeros);
  arma::vec omega(2+cov_num);
  omega.subvec(2,1+cov_num) = abs(weight)/sum(abs(weight));
  for(int i = 0; i< n; i++){
    D = HHunequal(P, D, profiles.col(i), cov_num, level_num, omega, p, pi);
    Assign(i) = D(D.n_rows-1,0);
  }
  return Assign;
}

//[[Rcpp::export]]
arma::vec HHue(arma::mat profiles, double pi, arma::vec omega, double p = 0.85){
  unsigned int cov_num = profiles.n_rows;
  arma::vec level_num = max(profiles,1);
  arma::mat P = PStrR(profiles);
  int strt_num = P.n_cols;
  int n = profiles.n_cols;
  arma::vec Assign(n);
  arma::vec D(2+strt_num+sum(level_num),arma::fill::zeros);
  for(int i = 0; i< n; i++){
    D = HHunequal(P, D, profiles.col(i), cov_num, level_num, omega, p, pi);
    Assign(i) = D(D.n_rows-1,0);
  }
  return Assign;
}
