#include <Rcpp.h>
#include "CHESSFunctions.h"

Rcpp::NumericVector CalcWHAMSpecCharge(int NSpec, 
                                       Rcpp::CharacterVector SpecActCorr,
                                       Rcpp::NumericVector SpecConc,
                                       Rcpp::IntegerVector SpecCharge,
                                       Rcpp::IntegerVector SpecMC,
                                       int AqueousMC,
                                       Rcpp::NumericVector HumicSubstGramsPerLiter) {
  
  /* output */
  Rcpp::NumericVector WHAMSpecCharge(2);
  
  /* variables */
  int iSpec;
  
  WHAMSpecCharge(iHA) = 0;
  WHAMSpecCharge(iFA) = 0;
  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    if (SpecActCorr(iSpec) == "WHAMHA") {
      WHAMSpecCharge(iHA) += SpecConc(iSpec) / HumicSubstGramsPerLiter[iHA] * SpecCharge(iSpec);
    } else if (SpecActCorr(iSpec) == "WHAMFA") {
      WHAMSpecCharge(iFA) += SpecConc(iSpec) / HumicSubstGramsPerLiter[iFA] * SpecCharge(iSpec);
    }
  }

  //WHAMSpecCharge = WHAMSpecCharge / HumicSubstGramsPerLiter;
  //charge / g HS = (charge / L) * (L / g HS)

  return WHAMSpecCharge;

}