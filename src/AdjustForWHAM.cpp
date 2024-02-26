#include <Rcpp.h>
#include "CHESSFunctions.h"

void AdjustForWHAM(unsigned int NComp,
                   unsigned int NSpec,
                   Rcpp::CharacterVector CompName,
                   Rcpp::CharacterVector SpecActCorr,
                   Rcpp::IntegerVector SpecCharge,
                   Rcpp::NumericVector WHAMSpecCharge,
                   Rcpp::NumericVector SpecCtoM,
                   Rcpp::NumericVector &SpecConc,
                   Rcpp::NumericVector &TotConc,
                   Rcpp::NumericVector &TotMoles) {

  /* variables */
  unsigned int iComp;
  unsigned int iSpec;

  // Set -Z_HS to the "known" total for the Donnan component
  for (iComp = 0; iComp < NComp; iComp++) {
    if (CompName(iComp) == "DonnanHA") {
      TotMoles(iComp) = abs(WHAMSpecCharge(iHA));
      TotConc(iComp) = abs(WHAMSpecCharge(iHA)) / SpecCtoM(iComp);
    } else if (CompName(iComp) == "DonnanFA") {
      TotMoles(iComp) = abs(WHAMSpecCharge(iFA));
      TotConc(iComp) = abs(WHAMSpecCharge(iFA)) / SpecCtoM(iComp);
    }
  }

  /* The cation species concentrations should be canceled out if Z_HS is 
     positive, and anion species concentrations should be canceled out if Z_HS
     is negative */
  for (iSpec = NComp; iSpec < NSpec; iSpec++) {
    if ((SpecActCorr(iSpec) == "DonnanHA") && 
        (((WHAMSpecCharge(iHA) < 0) && (SpecCharge(iSpec) < 0)) || 
         ((WHAMSpecCharge(iHA) > 0) && (SpecCharge(iSpec) > 0)))) {
      SpecConc(iSpec) = 0;
    } else if ((SpecActCorr(iSpec) == "DonnanFA") &&
               (((WHAMSpecCharge(iFA) < 0) && (SpecCharge(iSpec) < 0)) || 
                ((WHAMSpecCharge(iFA) > 0) && (SpecCharge(iSpec) > 0)))) {
      SpecConc(iSpec) = 0;
    }
  }

}