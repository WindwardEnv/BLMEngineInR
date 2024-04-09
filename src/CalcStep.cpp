#include <RcppArmadillo.h>
#include "RcppArmaHelper.h"
#include "CHESSFunctions.h"
#include <armadillo>

//' @title Calculate the Newton-Raphson step
//'
//' @param JacobianMatrix numeric matrix (NComp x NComp), the Jacobian matrix
//'   ("Z")
//' @param Resid numeric vector (NComp), the residuals = calculated totals -
//'   known totals
//' @param NComp integer, the number of components
//' @param CompType character vector (NComp), the type of component. It should be
//'   a fixed set of values (MassBal, FixedConc, Substituted, ChargeBal, SurfPot)
//' @param CompName character vector (NComp), the names of the components
//'
//' @return numeric vector (NComp), the N-R step to take for each component ("X"
//'   in C(i+1) = C(i) - X)
//'
Rcpp::NumericVector CalcStep(Rcpp::NumericMatrix JacobianMatrix,
                             Rcpp::NumericVector Resid,
                             Rcpp::NumericVector CompConc,
                             Rcpp::NumericVector TotMoles,
                             Rcpp::NumericVector CalcTotMoles,
                             int NComp,
                             Rcpp::CharacterVector CompType,
                             Rcpp::CharacterVector CompName){
  /* output */
  Rcpp::NumericVector CompConcStep(NComp);
  CompConcStep.names() = CompName;

  /* variables */
  int iComp, iComp2;
  int i, j;
  int NSolve;
  NSolve = 0;
  for (iComp = 0; iComp < NComp; iComp++){
    if ((CompType(iComp) != "FixedConc") && (CompType(iComp) != "FixedAct")) { 
      NSolve++; 
    }
  }
  
  Rcpp::NumericVector CompConcStepSolve(NSolve);
  Rcpp::NumericVector ResidSolve(NSolve);  
  Rcpp::NumericMatrix JacobianMatrixSolve(NSolve, NSolve);
  Rcpp::NumericMatrix JacobianMatrixInv(NSolve, NSolve);

  /*arma::mat ArmaCompConcStepSolve(NSolve, 1);
  arma::mat ArmaResidSolve(NSolve, 1);
  arma::mat ArmaJacobianMatSolve(NSolve, NSolve);
  arma::mat ArmaJacobianMatInv(NSolve, NSolve);*/
  
  //Pull out sub-set that should be solved
  i = 0;
  for (iComp = 0; iComp < NComp; iComp++){
    if ((CompType(iComp) != "FixedConc") && (CompType(iComp) != "FixedAct")) {
      ResidSolve(i) = Resid(iComp);
      //ArmaResidSolve(i, 0) = Resid(iComp);
      j = 0;
      for (iComp2 = 0; iComp2 < NComp; iComp2++){
        if ((CompType(iComp2) != "FixedConc") && 
            (CompType(iComp2) != "FixedAct")) {
          JacobianMatrixSolve(i, j) = JacobianMatrix(iComp, iComp2);
          //ArmaJacobianMatSolve(i, j) = JacobianMatrix(iComp, iComp2);
          j++;
        }
      }
      i++;
    }
  }


  /*Rcpp::Rcout << "ResidSolve = [" << ResidSolve << "]" << std::endl;
  std::cout << "ArmaResidSolve = [";
  for (i = 0; i < NSolve; i++) {
    if (i != 0) { std::cout << " "; }
    std::cout << ArmaResidSolve(i, 0);
  }
  std::cout << "]" << std::endl;

  Rcpp::Rcout << "JacobianMatrixSolve = [" << JacobianMatrixSolve << "]" << std::endl;
  std::cout << "ArmaJacobianMatSolve = [";
  for (i = 0; i < NSolve; i++) {
    if (i != 0) { std::cout << std::endl; }
    for (j = 0; j < NSolve; j++) {
      if (j != 0) { std::cout << " "; }
      std::cout << ArmaJacobianMatSolve(i, j);
    }
  }
  std::cout << "]" << std::endl; */

 
  try {
    // find the matrix inverse of JacobianMatrix by SVD
    JacobianMatrixInv = SvdInverse(JacobianMatrixSolve);
    if (JacobianMatrixInv.nrow() != NSolve) {
      throw 10;
    }
    CompConcStepSolve = RcppMatMult(JacobianMatrixInv, ResidSolve);
    
    /*// Doing the calculations with Armadillo classes
    ArmaJacobianMatInv = SvdInverse(ArmaJacobianMatSolve);
    if (ArmaJacobianMatInv.n_rows != NSolve) {
      // If we wanted to not do SVD...
      ArmaJacobianMatInv = arma::inv(ArmaJacobianMatSolve);
      throw 10;
    }
    ArmaCompConcStepSolve = ArmaJacobianMatInv * ArmaResidSolve;  */
  }
  catch (int e) {
    if (e == 10) {
      Rcpp::Rcout << "Singluar Matrix?" << std::endl;
      i = 0;
      for (iComp = 0; iComp < NComp; iComp++){
        if ((CompType(iComp) != "FixedConc") && (CompType(iComp) != "FixedAct")) {
          //CompConcStepSolve(i) = CompConc(iComp) * (1 - (TotMoles(iComp) / CalcTotMoles(iComp) + 1) / 5);          
          CompConcStepSolve(i) = CompConc(iComp) * (-1) * (TotMoles(iComp) / CalcTotMoles(iComp));
          i++;
        }
      }

    } else {
      Rcpp::Rcout << "Unknown exception occurred" << std::endl;
    }
  }
  /*Rcpp::Rcout << "CompConcStepSolve = [" << CompConcStepSolve << "]" << std::endl;
  std::cout << "ArmaCompConcStepSolve = [";
  for (i = 0; i < NSolve; i++) {
    if (i != 0) { std::cout << " "; }
    std::cout << ArmaCompConcStepSolve(i, 0);
  }
  std::cout << "]" << std::endl;*/
  /*std::cout << "ArmaCompConcStepSolve1N = [";
  for (i = 0; i < NSolve; i++) {
    if (i != 0) { std::cout << " "; }
    std::cout << ArmaCompConcStepSolve1N(0, i);
  }
  std::cout << "]" << std::endl;*/

  i = 0;
  for (iComp = 0; iComp < NComp; iComp++){
    if ((CompType(iComp) == "FixedConc") || (CompType(iComp) == "FixedAct")) {
      CompConcStep(iComp) = 0;
    } else {
      CompConcStep(iComp) = CompConcStepSolve(i);      
      //CompConcStep(iComp) = ArmaCompConcStepSolve(i, 0);
      i++;
    }
  }
  return CompConcStep;
}

