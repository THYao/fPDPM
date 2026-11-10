#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
double log_dmvn_wdb(NumericVector& xVec_, NumericVector& mu_, arma::mat& W_, double& s2_) {
  int L_ = W_.n_rows;
  int d_ = W_.n_cols;
  // arma::mat I_L(L_,L_, arma::fill::eye);
  arma::mat I_d(d_,d_, arma::fill::eye);
  
  double detVal;
  double detSign;
  bool detOK = log_det(detVal, detSign, s2_ * I_d + W_.t() * W_);
  detVal = detVal + (L_-d_) * log(s2_);
  // Rcout << "detVal: " << detVal << std::endl;
  
  arma::mat M = W_.t() * W_ + s2_ * I_d;
  arma::mat M_inv = inv(M);
  arma::mat WM_inv = W_ * M_inv;
  
  arma::vec aVec(xVec_ - mu_);
  arma::mat expNumMat = aVec.t() * aVec - aVec.t() * WM_inv * W_.t() * aVec;
  double expNum = expNumMat(0,0)/s2_;
  // Rcout << "expNum: " << expNum << std::endl;
  
  double res = detVal + expNum + L_*log(2*arma::datum::pi);
  res = res * (-0.5);
  return res;
}


