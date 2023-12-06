#ifndef __RCPPARMAHELPER_H_INCLUDED__
#define __RCPPARMAHELPER_H_INCLUDED__

//should these be #included in the header file??
#ifndef RcppArmadillo__RcppArmadillo__h
  #include <RcppArmadillo.h>
#endif
#include <vector>

// Rcpp to armadillo conversion functions
// arma::vec RcppVectorToVector(Rcpp::NumericVector XRcppVec);
arma::mat RcppVectorToMatrix(Rcpp::NumericVector XRcppVec);//used in CalcResidual.cpp and RcppArmaHelper.cpp
arma::mat RcppMatrixToMatrix(Rcpp::NumericMatrix XRcppMat);//used in RcppArmaHelper.cpp
arma::mat RcppIntMatrixToMatrix(Rcpp::IntegerMatrix XRcppMat);//wants to be used in CalcResidual.cpp

//armadillo to Rcpp conversion functions
// Rcpp::NumericVector VectorToRcppVector(arma::vec Xvec);
Rcpp::NumericVector MatrixToRcppVector(arma::mat Xmat);//used in RcppArmaHelper.cpp
Rcpp::NumericMatrix MatrixToRcppMatrix(arma::mat Xmat);//used in RcppArmaHelper.cpp

//matrix mathy functions using arma classes
arma::mat SvdInverse(arma::mat X); //used in RcppArmaHelper.cpp (for Rcpp class version of function)

//matrix mathy functions using Rcpp classes
// Rcpp::NumericMatrix RcppMatMult(Rcpp::NumericMatrix A, Rcpp::NumericMatrix B);
Rcpp::NumericVector RcppMatMult(Rcpp::NumericMatrix A, Rcpp::NumericVector B);//used in CalcStep.cpp
// Rcpp::NumericVector RcppMatMult(Rcpp::NumericVector A, Rcpp::NumericMatrix B);
Rcpp::NumericMatrix SvdInverse(Rcpp::NumericMatrix XRcppMat);//used in CalcStep.cpp


#endif // __RCPPARMAHELPER_H_INCLUDED__
