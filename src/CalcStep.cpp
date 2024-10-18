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
  Rcpp::LogicalVector CompSolve(NComp);
  NSolve = 0;
  for (iComp = 0; iComp < NComp; iComp++){
    if ((CompType(iComp) != CTYPE_FIXEDCONC) && (CompType(iComp) != CTYPE_FIXEDACT)
        //&& (CompType(iComp) != CTYPE_DONNANHA) && (CompType(iComp) != CTYPE_DONNANFA)
       ) { 
      CompSolve[iComp] = true;
      NSolve++; 
    } else {
      CompSolve[iComp] = false;
    }
  }
  
  Rcpp::NumericVector CompConcStepSolve(NSolve);
  /*Rcpp::NumericVector ResidSolve(NSolve);  
  Rcpp::NumericMatrix JacobianMatrixSolve(NSolve, NSolve);
  Rcpp::NumericMatrix JacobianMatrixInv(NSolve, NSolve);*/

  arma::mat ArmaCompConcStepSolve(NSolve, 1);
  arma::mat ArmaResidSolve(NSolve, 1);
  arma::mat ArmaJacobianMatSolve(NSolve, NSolve);
  arma::mat ArmaJacobianMatInv(NSolve, NSolve);
  
  //Pull out sub-set that should be solved
  i = 0;
  for (iComp = 0; iComp < NComp; iComp++){
    if (CompSolve[iComp]) {
      //ResidSolve(i) = Resid(iComp);
      ArmaResidSolve(i, 0) = Resid(iComp);
      j = 0;
      for (iComp2 = 0; iComp2 < NComp; iComp2++){
        if (CompSolve[iComp2]) {
          //JacobianMatrixSolve(i, j) = JacobianMatrix(iComp, iComp2);
          ArmaJacobianMatSolve(i, j) = JacobianMatrix(iComp, iComp2);
          j++;
        }
      }
      i++;
    }
  }


  try {
    
    /*// find the matrix inverse of JacobianMatrix by SVD
    JacobianMatrixInv = SvdInverse(JacobianMatrixSolve);

    if (JacobianMatrixInv.nrow() != NSolve) {
      throw ERROR_MATRIX_INVERSION;
    }
    CompConcStepSolve = RcppMatMult(JacobianMatrixInv, ResidSolve);*/
    
    // Doing the calculations with Armadillo classes
    bool InverseSuccess;
    // Start with the general inverse
    // (https://arma.sourceforge.net/docs.html#inv)
    InverseSuccess = arma::inv(ArmaJacobianMatInv, ArmaJacobianMatSolve, 
                               arma::inv_opts::allow_approx);    
    
    if (!InverseSuccess) {
      // Try to calculate the Moore-Penrose pseudo-invserse (uses SVD)
      // (https://arma.sourceforge.net/docs.html#pinv)
      InverseSuccess = arma::pinv(ArmaJacobianMatInv, ArmaJacobianMatSolve);
    }
    if (!InverseSuccess) { 

      /*Rcpp::Rcout << "ResidSolve = [" << ResidSolve << "]" << std::endl;
      Rcpp::Rcout << "ArmaResidSolve = [";
      for (i = 0; i < NSolve; i++) {
        if (i != 0) { Rcpp::Rcout << " "; }
        Rcpp::Rcout << ArmaResidSolve(i, 0);
      }
      Rcpp::Rcout << "]" << std::endl;

      Rcpp::Rcout << "JacobianMatrixSolve = [" << JacobianMatrixSolve << "]" << std::endl;
      Rcpp::Rcout << "ArmaJacobianMatSolve = [";
      for (i = 0; i < NSolve; i++) {
        if (i != 0) { Rcpp::Rcout << std::endl; }
        for (j = 0; j < NSolve; j++) {
          if (j != 0) { Rcpp::Rcout << " "; }
          Rcpp::Rcout << ArmaJacobianMatSolve(i, j);
        }
      }
      Rcpp::Rcout << "]" << std::endl;*/

      throw ERROR_MATRIX_INVERSION; 
    }
    /*ArmaJacobianMatInv = SvdInverse(ArmaJacobianMatSolve);
    if (ArmaJacobianMatInv.n_rows != NSolve) {
      // If we wanted to not do SVD...
      ArmaJacobianMatInv = arma::inv(ArmaJacobianMatSolve);
      throw ERROR_SINGULAR_MATRIX;
    }*/
    ArmaCompConcStepSolve = ArmaJacobianMatInv * ArmaResidSolve;
    CompConcStepSolve = MatrixToRcppVector(ArmaCompConcStepSolve);
  }
  catch (int e) {
    /*if (e == ERROR_SINGULAR_MATRIX) {
      Rcpp::Rcout << "Singluar Matrix" << std::endl;
    } else if (e == ERROR_MATRIX_INVERSION) {
      Rcpp::Rcout << "Matrix inversion failed." << std::endl;
    } else {
      Rcpp::Rcout << "Unknown exception occurred" << std::endl;
    }*/
    
    // default to a brute-force step value
    i = 0;
    for (iComp = 0; iComp < NComp; iComp++){
      if (CompSolve[iComp]) {
        //CompConcStepSolve(i) = CompConc(iComp) * (1 - (TotMoles(iComp) / CalcTotMoles(iComp) + 1) / 5);          
        CompConcStepSolve(i) = CompConc(iComp) * (-1) * (TotMoles(iComp) / CalcTotMoles(iComp));
        i++;
      }
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
    if (CompSolve[iComp]) {
      CompConcStep(iComp) = CompConcStepSolve(i);      
      //CompConcStep(iComp) = ArmaCompConcStepSolve(i, 0);
      i++;
    } else {
      CompConcStep(iComp) = 0;
    }
  }
  return CompConcStep;
}

Rcpp::NumericVector CalcStepBrute(int NComp,
                                  Rcpp::CharacterVector CompName, 
                                  Rcpp::CharacterVector CompType, 
                                  Rcpp::NumericVector CompConc, 
                                  Rcpp::NumericVector TotConc, 
                                  Rcpp::NumericVector CalcTotConc,
                                  bool DoTox,
                                  int MetalComp,
                                  double CACalculated,
                                  double CATarget) {
  /*outputs*/
  Rcpp::NumericVector CompConcStep(NComp);
    CompConcStep.names() = CompName;
  
  /*varaibles*/
  int iComp;

  for (iComp = 0; iComp < NComp; iComp++){
    if ((CompType(iComp) != CTYPE_FIXEDCONC) && (CompType(iComp) != CTYPE_FIXEDACT)) { 
      CompConcStep(iComp) = CompConc(iComp) * (1 - TotConc(iComp) / CalcTotConc(iComp));
      if (CompConcStep(iComp) < (CompConc(iComp) - TotConc(iComp))) {
        CompConcStep(iComp) = (CompConc(iComp) - TotConc(iComp));
      }
    } else if (DoTox && (iComp == MetalComp)) {
      CompConcStep(iComp) = CompConc(iComp) * (1 - CATarget / CACalculated);
    } else {
      CompConcStep(iComp) = 0.0;
    }
  }

  return CompConcStep;
}