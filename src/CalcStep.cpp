#include <RcppArmadillo.h>
#include "RcppArmaHelper.h"
#include "CHESSFunctions.h"

//' @title Calculate the Newton-Raphson step
//'
//' @param JacobianMatrix numeric matrix (NComp x NComp), the Jacobian matrix
//'   ("Z")
//' @param Resid numeric vector (NComp), the residuals = calculated totals -
//'   known totals
//' @param NComp integer, the number of components
//' @param CompType character vector (NComp), the type of component. It should be
//'   a fixed set of values (MassBal, FixedAct, Substituted, ChargeBal, SurfPot)
//' @param CompName character vector (NComp), the names of the components
//'
//' @return numeric vector (NComp), the N-R step to take for each component ("X"
//'   in C(i+1) = C(i) - X)
//'
 Rcpp::NumericVector CalcStep(Rcpp::NumericMatrix JacobianMatrix,
                              Rcpp::NumericVector Resid,
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
     if (CompType(iComp) != "FixedAct"){ NSolve++; }
   }
   Rcpp::NumericVector CompConcStepSolve(NSolve);
   Rcpp::NumericVector ResidSolve(NSolve);
   Rcpp::NumericMatrix JacobianMatrixSolve(NSolve, NSolve);
   Rcpp::NumericMatrix JacobianMatrixInv(NSolve, NSolve);

   //Pull out sub-set that should be solved
   i = 0;
   for (iComp = 0; iComp < NComp; iComp++){
     if (CompType(iComp) != "FixedAct"){
       ResidSolve(i) = Resid(iComp);
       j = 0;
       for (iComp2 = 0; iComp2 < NComp; iComp2++){
         if (CompType(iComp2) != "FixedAct"){
           JacobianMatrixSolve(i, j) = JacobianMatrix(iComp, iComp2);
           j++;
         }
       }
       i++;
     }
   }

   // find the matrix inverse of JacobianMatrix by SVD
   JacobianMatrixInv = SvdInverse(JacobianMatrixSolve);

   // // If we wanted to not do SVD...
   // JacobianMatrixInv = solve(JacobianMatrix.solve)

   CompConcStepSolve = RcppMatMult(JacobianMatrixInv, ResidSolve);

   i = 0;
   for (iComp = 0; iComp < NComp; iComp++){
     if (CompType(iComp) == "FixedAct"){
       CompConcStep(iComp) = 0;
     } else {
       CompConcStep(iComp) = CompConcStepSolve(i);
       i++;
     }
   }

   return CompConcStep;

 }

