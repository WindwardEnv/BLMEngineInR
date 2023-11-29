#include <Rcpp.h>
using namespace Rcpp;

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
//' @param CompConc numeric vector (NComp) of component concentrations, input values
//'   are from this iteration
//' @param CompName numeric vector (NComp) with the names of the components
//'
//' @return  numeric cector CompConc (NComp) modified for the next iteration
//'

// [[Rcpp::export]]
Rcpp::NumericVector CompUpdate(unsigned int NComp,
                               Rcpp::NumericVector CompConcStep,
                               Rcpp::NumericVector CompConc,
                               Rcpp::CharacterVector CompName){

  /* outputs */
  Rcpp::NumericVector newCompConc(NComp);
  newCompConc.names() = CompName;

  /* variables */
  // Rcpp::NumericVector newCompConc(NComp); //oldCompConc = CompConc;
  // Rcpp::LogicalVector ltzero(NComp);
  unsigned int iComp;


  // newCompConc = CompConc - CompConcStep;
  // newCompConc[ltzero] = CompConc[ltzero] / 10;

  for (iComp = 0; iComp < NComp; iComp++){
    if (CompConcStep(iComp) >= CompConc(iComp)) {
      newCompConc(iComp) = CompConc(iComp) / 10;
    } else {
      newCompConc(iComp) = CompConc(iComp) - CompConcStep(iComp);
    }
  }//NEXT iComp

  return newCompConc;

}
