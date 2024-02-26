#include <Rcpp.h>
#include "CHESSFunctions.h"

Rcpp::NumericVector CalcIonicStrengthEffects(double IonicStrength,
                                             Rcpp::NumericVector WHAMSpecCharge,
                                             unsigned int NSpec,
                                             Rcpp::IntegerVector SpecCharge,
                                             Rcpp::NumericVector SpecK,
                                             Rcpp::CharacterVector SpecActCorr,
                                             Rcpp::NumericVector wP) {
  /* output */
  Rcpp::NumericVector SpecKISAdj(NSpec);
    SpecKISAdj.names() = SpecK.names();

  /* variables */
  Rcpp::NumericVector W = wP * log10(IonicStrength);
  Rcpp::NumericVector WZ2 = 2 * W * WHAMSpecCharge;
  //unsigned int iHA = 0;
  //unsigned int iFA = 1;
  unsigned int iSpec;
  
  //check that WHAMSpecCharge is negative for both -- do we need to??
  
  //adjust the intrinsic K's for WHAM species based on the charge
  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    SpecKISAdj(iSpec) = SpecK(iSpec);
    if (SpecCharge(iSpec) != 0) {
      if (SpecActCorr(iSpec) == "WHAMHA") {
        SpecKISAdj(iSpec) = SpecK(iSpec) * exp(WZ2(iHA) * abs(SpecCharge(iSpec)));
      } else if (SpecActCorr(iSpec) == "WHAMFA") {
        SpecKISAdj(iSpec) = SpecK(iSpec) * exp(WZ2(iFA) * abs(SpecCharge(iSpec)));
      }
    }    
  }

  return SpecKISAdj;

}