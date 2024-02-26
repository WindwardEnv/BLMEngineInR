#include <Rcpp.h>
#include "CHESSFunctions.h"

Rcpp::NumericVector CalcWHAMSpecCharge(unsigned int NSpec, 
                                       Rcpp::CharacterVector SpecActCorr,
                                       Rcpp::NumericVector SpecMoles,
                                       Rcpp::IntegerVector SpecCharge,
                                       Rcpp::IntegerVector SpecMC,
                                       int AqueousMC) {
  
  /* output */
  Rcpp::NumericVector WHAMSpecCharge(2);
  
  /* variables */
  unsigned int iSpec;
  
  WHAMSpecCharge(iHA) = 0;
  WHAMSpecCharge(iFA) = 0;
  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    if ((SpecMC(iSpec) == AqueousMC) && (SpecCharge(iSpec) != 0)) {
      if (SpecActCorr(iSpec) == "WHAMHA") {
        WHAMSpecCharge(iHA) += SpecMoles(iSpec) * SpecCharge(iSpec);
      } else if (SpecActCorr(iSpec) == "WHAMFA") {
        WHAMSpecCharge(iFA) += SpecMoles(iSpec) * SpecCharge(iSpec);
      }
    }
  }

  return WHAMSpecCharge;

}