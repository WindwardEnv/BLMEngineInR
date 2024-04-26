#include <Rcpp.h>
#include "CHESSFunctions.h"

//' @title Iterative step improvement in component concentrations
//'
//' @description CompUpdate calculates an iterative improvement on the component
//' concentrations based on the Newton-Raphson solution from the current
//' iteration.
//'
//' @details If the iteration would cause the adjusted component concentrations
//' to be less than zero, then the component concentration is simply divided by
//' 10 for this iteration.
//'
//' @param NComp integer, the number of components
//' @param CompConcStep numeric vector (NComp) of adjustments to the component
//'   concentrations
//' @param CompName numeric vector (NComp) with the names of the components
//' @param CompConc (INPUT & OUTPUT) numeric vector (NComp) of component
//'   concentrations, input values are from this iteration
//'
//' @return  numeric cector CompConc (NComp) modified for the next iteration
//'
void CompUpdate(int NComp, 
                Rcpp::NumericVector CompConcStep,
                Rcpp::CharacterVector CompType,
                Rcpp::NumericVector &CompConc){

  /* variables */
  Rcpp::NumericVector oldCompConc(NComp); //oldCompConc = CompConc;
  // Rcpp::LogicalVector ltzero(NComp);
  int iComp;


  // newCompConc = CompConc - CompConcStep;
  // newCompConc[ltzero] = CompConc[ltzero] / 10;

  for (iComp = 0; iComp < NComp; iComp++){
    oldCompConc[iComp] = CompConc[iComp];
    if (CompConcStep[iComp] >= oldCompConc[iComp]) {
      CompConc[iComp] = oldCompConc[iComp] / 10;
    } else {
      CompConc[iComp] = oldCompConc[iComp] - CompConcStep[iComp];
    }
    //if (((CompType[iComp] == TYPE_DONNANHA) || (CompType[iComp] == TYPE_DONNANFA))) {
    //  if (CompConc[iComp] < 1.0) { CompConc[iComp] = 1.0; }
    //}
  }//NEXT iComp

}
