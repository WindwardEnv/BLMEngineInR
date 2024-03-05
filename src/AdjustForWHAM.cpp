#include <Rcpp.h>
#include "CHESSFunctions.h"

void AdjustForWHAM(int NComp,
                   int NSpec,
                   Rcpp::CharacterVector CompName,
                   Rcpp::CharacterVector SpecActCorr,
                   Rcpp::IntegerVector SpecCharge,
                   Rcpp::NumericVector WHAMSpecCharge,
                   Rcpp::NumericVector SpecCtoM,
                   Rcpp::NumericVector &SpecConc,
                   Rcpp::NumericVector &TotConc,
                   Rcpp::NumericVector &TotMoles) {

  /* variables */
  int iComp;
  int iSpec;

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

void WHAMAdjustments(
  int NMass,   
  Rcpp::NumericVector MassAmt,
  Rcpp::NumericVector &MassAmtAdj,
  int NComp,
  Rcpp::CharacterVector CompType,
  Rcpp::NumericVector &TotConc,
  Rcpp::NumericVector &TotMoles,
  int NSpec,
  Rcpp::NumericVector &SpecConc,
  Rcpp::IntegerVector SpecMC,
  Rcpp::CharacterVector SpecActCorr,
  Rcpp::IntegerVector SpecCharge,
  Rcpp::NumericVector SpecKTempAdj,
  Rcpp::NumericVector &SpecKISTempAdj,
  Rcpp::NumericVector &SpecCtoMAdj,  
  double IonicStrength,
  Rcpp::NumericVector &WHAMSpecCharge,
  int AqueousMC,
  Rcpp::IntegerVector WHAMDonnanMC,
  Rcpp::NumericVector SolHS,
  Rcpp::NumericVector wMolWt,
  Rcpp::NumericVector wRadius,
  Rcpp::NumericVector wP,
  double wDLF,
  double wKZED
) {

  /* variables */
  int iComp;
  int iSpec;

  //Calculate the charge on the organic matter
  WHAMSpecCharge = CalcWHAMSpecCharge(NSpec, SpecActCorr,
                                      SpecConc * SpecCtoMAdj, 
                                      SpecCharge, SpecMC, AqueousMC);

  //Adjust organic matter specific binding based on ionic strength
  SpecKISTempAdj = CalcIonicStrengthEffects(IonicStrength, WHAMSpecCharge, 
                                            NSpec, SpecCharge, SpecKTempAdj, 
                                            SpecActCorr, wP);

  //Calculate the portion of the solution that's in the diffuse layer
  MassAmtAdj = CalcDonnanLayerVolume(NMass, NSpec, IonicStrength, MassAmt, 
                                     AqueousMC, WHAMDonnanMC, wMolWt, 
                                     wRadius, wDLF, wKZED, WHAMSpecCharge, 
                                     SolHS);
  SpecCtoMAdj = MassAmtAdj[SpecMC];

  // Set -Z_HS to the "known" total for the Donnan component
  for (iComp = 0; iComp < NComp; iComp++) {
    if (CompType(iComp) == "DonnanHA") {
      TotMoles(iComp) = abs(WHAMSpecCharge(iHA));
      TotConc(iComp) = abs(WHAMSpecCharge(iHA)) / SpecCtoMAdj(iComp);
    } else if (CompType(iComp) == "DonnanFA") {
      TotMoles(iComp) = abs(WHAMSpecCharge(iFA));
      TotConc(iComp) = abs(WHAMSpecCharge(iFA)) / SpecCtoMAdj(iComp);
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