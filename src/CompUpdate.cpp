#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericVector CompUpdate(NumericMatrix X) {
  NumericVector d;
  NumericMatrix U;
  NumericMatrix V;

  arma::svd(U,d,V,X);
  return d;
}
