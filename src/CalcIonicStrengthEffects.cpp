#include <Rcpp.h>
#include "CHESSFunctions.h"

Rcpp::NumericVector CalcIonicStrengthEffects(double IonicStrength,
                                             Rcpp::NumericVector WHAMSpecCharge,
                                             int NSpec,
                                             Rcpp::IntegerVector SpecCharge,
                                             Rcpp::NumericVector SpecK,
                                             Rcpp::CharacterVector SpecActCorr,
                                             Rcpp::NumericVector wP) {
  /* output */
  Rcpp::NumericVector SpecKISAdj(NSpec);
    SpecKISAdj.names() = SpecK.names();

  /* variables */
  Rcpp::NumericVector W = wP * log10(IonicStrength);
  Rcpp::NumericVector WZ2 = -2 * W * WHAMSpecCharge;
  int iSpec;
  
  //adjust the intrinsic K's for WHAM species based on the charge
  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    SpecKISAdj(iSpec) = SpecK(iSpec);
    if (SpecCharge(iSpec) != 0) {
      if (SpecActCorr(iSpec) == ACTYPE_WHAMHA) {
        SpecKISAdj(iSpec) = SpecK(iSpec) * exp(WZ2(iHA) * SpecCharge(iSpec));
      } else if (SpecActCorr(iSpec) == ACTYPE_WHAMFA) {
        SpecKISAdj(iSpec) = SpecK(iSpec) * exp(WZ2(iFA) * SpecCharge(iSpec));
      }
    }    
  }

  return SpecKISAdj;

}