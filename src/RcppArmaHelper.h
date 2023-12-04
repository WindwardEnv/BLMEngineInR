#ifndef __RcppArmaHelper_H_INCLUDED__
#define __RcppArmaHelper_H_INCLUDED__

#include <RcppArmadillo.h>
#include <vector>

// Rcpp to armadillo conversion functions
arma::vec RcppVectorToVector(Rcpp::NumericVector XRcppVec);
arma::mat RcppVectorToMatrix(Rcpp::NumericVector XRcppVec);
arma::mat RcppMatrixToMatrix(Rcpp::NumericMatrix XRcppMat);
arma::mat RcppMatrixToMatrix(Rcpp::IntegerMatrix XRcppMat);

//armadillo to Rcpp conversion functions
Rcpp::NumericVector VectorToRcppVector(arma::vec Xvec);
Rcpp::NumericMatrix MatrixToRcppMatrix(arma::mat Xmat);
Rcpp::NumericVector MatrixToRcppVector(arma::mat Xmat);
Rcpp::NumericMatrix RcppMatMult(Rcpp::NumericMatrix A, Rcpp::NumericMatrix B);
Rcpp::NumericVector RcppMatMult(Rcpp::NumericMatrix A, Rcpp::NumericVector B);
Rcpp::NumericVector RcppMatMult(Rcpp::NumericVector A, Rcpp::NumericMatrix B);

// specific mathy functions
arma::mat SvdInverse(arma::mat X);
Rcpp::NumericMatrix SvdInverse(Rcpp::NumericMatrix XRcppMat);


#endif // __RcppArmaHelper_H_INCLUDED__
