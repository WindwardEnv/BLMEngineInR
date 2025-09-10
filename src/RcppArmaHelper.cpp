// Copyright 2024 Windward Environmental LLC
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <vector>
#include <cfloat>
#include "RcppArmaHelper.h"

/*
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

arma::mat RcppIntMatrixToMatrix(Rcpp::IntegerMatrix XRcppMat){
  int Xrows = XRcppMat.nrow();
  int Xcols = XRcppMat.ncol();
  std::vector<double> Xstd = Rcpp::as<std::vector<double>>(XRcppMat);//Rcpp::NumericMatrix -> std::vector<double>
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
*/

Rcpp::NumericVector MatrixToRcppVector(arma::mat Xmat){
  std::vector<double> Xstd = arma::conv_to< std::vector<double> >::from(Xmat.as_col());//arma::mat -> std::vector
  Rcpp::NumericVector XRcppvec = Rcpp::wrap(Xstd);//std::vector -> Rcpp::NumericVector
  return XRcppvec;
}

/*
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
arma::mat SvdInverse(arma::mat X) {

  arma::mat Xinv;
  arma::mat U, UT;
  arma::vec s;
  arma::mat V;

  arma::svd(U,s,V,X);
  arma::mat dinv = arma::diagmat(1/s);

  bool tmp = false;  
  if (tmp) {
    Rcpp::NumericMatrix RcppU = MatrixToRcppMatrix(U);
    //Rcpp::Rcout << RcppU << std::endl;
  }  
  if (tmp) {
    Rcpp::NumericMatrix RcppV = MatrixToRcppMatrix(V);
    //Rcpp::Rcout << RcppV << std::endl;
  }  
  if (tmp) {
    Rcpp::NumericMatrix RcppX = MatrixToRcppMatrix(X);
    //Rcpp::Rcout << RcppX << std::endl;
  }
  if (tmp) {
    Rcpp::NumericMatrix Rcppdinv = MatrixToRcppMatrix(s);
    //Rcpp::Rcout << 1 / Rcppdinv << std::endl;
    //Rcpp::Rcout << Rcppdinv << std::endl;
  }  
  
  double max_dinv = arma::max(1 / s);
  double min_dinv = arma::min(1 / s);
  double condition = max_dinv / min_dinv;
  const double condition_limit = 4.5036e+15;//707945784384137.38;//100000000000000000000.0;//pow(10, DBL_DIG * 0.99);
  //                             012345678901234567890123456789012
  int n_zero = 0;
  if (condition > condition_limit) {
    double near_zero_limit = min_dinv * condition_limit;
    for (unsigned int i = 0; i < s.size(); i++) {
      if (dinv(i,i) > near_zero_limit) {
        dinv(i,i) = 0.0;
        n_zero++;
      }
    }
  }

  //Rcpp::Rcout << "     " << condition << ", " << max_dinv << ", " << min_dinv << ", " << n_zero << std::endl;

  UT = trans(U);
  Xinv = V * dinv * UT;

  if (tmp) {
    Rcpp::NumericMatrix RcppXinv = MatrixToRcppMatrix(Xinv);
    //Rcpp::Rcout << RcppXinv << std::endl;
  }

  return Xinv;
}

Rcpp::NumericMatrix SvdInverse(Rcpp::NumericMatrix XRcppMat){

  arma::mat X = RcppMatrixToMatrix(XRcppMat);

  arma::mat Xinv = SvdInverse(X);

  return MatrixToRcppMatrix(Xinv);
}
*/
