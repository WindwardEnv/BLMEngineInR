#include <RcppArmadillo.h>
#include <vector>
#include "RcppArmaHelper.h"

// Rcpp to armadillo conversion functions
arma::vec RcppVectorToVector(Rcpp::NumericVector XRcppVec){
  std::vector<double> Xstd = Rcpp::as<std::vector<double>>(XRcppVec);//Rcpp::NumericVector -> std::vector<double>
  arma::vec Xvec(Xstd); //std::vector -> arma::vec
  return Xvec;
}

arma::mat RcppVectorToMatrix(Rcpp::NumericVector XRcppVec){
  std::vector<double> Xstd = Rcpp::as<std::vector<double>>(XRcppVec);//Rcpp::NumericVector -> std::vector<double>
  arma::vec Xvec(Xstd); //std::vector -> arma::vec
  arma::mat Xmat = arma::reshape(Xvec, XRcppVec.size(), 1);//arma::vec -> arma::mat
  return Xmat;
}

arma::mat RcppMatrixToMatrix(Rcpp::NumericMatrix XRcppMat){
  int Xrows = XRcppMat.nrow();
  int Xcols = XRcppMat.ncol();
  std::vector<double> Xstd = Rcpp::as<std::vector<double>>(XRcppMat);//Rcpp::NumericMatrix -> std::vector<double>
  arma::vec Xvec(Xstd); //std::vector -> arma::vec
  arma::mat Xmat = arma::reshape(Xvec, Xrows, Xcols);//arma::vec -> arma::mat
  return Xmat;
}
arma::mat RcppMatrixToMatrix(Rcpp::IntegerMatrix XRcppMat){
  int Xrows = XRcppMat.nrow();
  int Xcols = XRcppMat.ncol();
  std::vector<int> Xstd = Rcpp::as<std::vector<int>>(XRcppMat);//Rcpp::NumericMatrix -> std::vector<double>
  arma::vec Xvec(Xstd); //std::vector -> arma::vec
  arma::mat Xmat = arma::reshape(Xvec, Xrows, Xcols);//arma::vec -> arma::mat
  return Xmat;
}


//armadillo to Rcpp conversion functions
Rcpp::NumericVector VectorToRcppVector(arma::vec Xvec){
  std::vector<double> Xstd = arma::conv_to< std::vector<double> >::from(Xvec);//arma::vec -> std::vector
  Rcpp::NumericVector XRcppvec = Rcpp::wrap(Xstd);//std::vector -> Rcpp::NumericVector
  return XRcppvec;
}

Rcpp::NumericMatrix MatrixToRcppMatrix(arma::mat Xmat){
  int Xrows = Xmat.n_rows;
  int Xcols = Xmat.n_cols;
  std::vector<double> Xstd = arma::conv_to< std::vector<double> >::from(Xmat.as_col());//arma::mat -> std::vector
  Rcpp::NumericVector XRcppvec = Rcpp::wrap(Xstd);//std::vector -> Rcpp::NumericVector
  XRcppvec.attr("dim") = Rcpp::Dimension(Xrows, Xcols);//set dimensions
  Rcpp::NumericMatrix XRcppmat = Rcpp::as<Rcpp::NumericMatrix>(XRcppvec);//Rcpp::NumericVector -> Rcpp::NumericMatrix
  return XRcppmat;
}

Rcpp::NumericVector MatrixToRcppVector(arma::mat Xmat){
  std::vector<double> Xstd = arma::conv_to< std::vector<double> >::from(Xmat.as_col());//arma::mat -> std::vector
  Rcpp::NumericVector XRcppvec = Rcpp::wrap(Xstd);//std::vector -> Rcpp::NumericVector
  return XRcppvec;
}

Rcpp::NumericMatrix RcppMatMult(Rcpp::NumericMatrix A, Rcpp::NumericMatrix B){
  arma::mat Amat = RcppMatrixToMatrix(A);
  arma::mat Bmat = RcppMatrixToMatrix(B);
  arma::mat Cmat = Amat * Bmat;
  Rcpp::NumericMatrix C = MatrixToRcppMatrix(Cmat);
  return C;
}

Rcpp::NumericVector RcppMatMult(Rcpp::NumericMatrix A, Rcpp::NumericVector B){
  arma::mat Amat = RcppMatrixToMatrix(A);
  arma::mat Bmat = RcppVectorToMatrix(B);
  arma::mat Cmat = Amat * Bmat;
  Rcpp::NumericVector C = MatrixToRcppVector(Cmat);
  return C;
}

Rcpp::NumericVector RcppMatMult(Rcpp::NumericVector A, Rcpp::NumericMatrix B){
  arma::mat Amat = RcppVectorToMatrix(A);
  arma::mat Bmat = RcppMatrixToMatrix(B);
  arma::mat Cmat = Amat * Bmat;
  Rcpp::NumericVector C = MatrixToRcppVector(Cmat);
  return C;
}


// specific mathy functions
arma::mat SvdInverse(arma::mat X){

  arma::mat Xinv;
  arma::mat U, UT;
  arma::vec s;
  arma::mat V;

  arma::svd(U,s,V,X);
  arma::mat dinv = arma::diagmat(1/s);
  UT = trans(U);
  Xinv = V * dinv * UT;

  return Xinv;
}

Rcpp::NumericMatrix SvdInverse(Rcpp::NumericMatrix XRcppMat){

  arma::mat X = RcppMatrixToMatrix(XRcppMat);

  arma::mat Xinv = SvdInverse(X);

  return MatrixToRcppMatrix(Xinv);
}
