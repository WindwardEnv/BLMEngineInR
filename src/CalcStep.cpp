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
      throw ERROR_MATRIX_INVERSION; 
    }
    ArmaCompConcStepSolve = ArmaJacobianMatInv * ArmaResidSolve;
    CompConcStepSolve = MatrixToRcppVector(ArmaCompConcStepSolve);
  }
  catch (int e) {
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
