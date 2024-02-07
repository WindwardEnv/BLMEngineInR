#include <cmath>
#include <Rcpp.h>
#include "CHESSFunctions.h"

//' @title Calculate the Diffuse Layer Volumes and Adjust CtoM
//' 
//' @description This function calculates a diffuse layer volume for the
//'   non-specific binding of ions to charged organic matter, then adjusts the
//'   SpecCtoM for aqueous inorganic species and Donnan species 
//' 
//' @details {text}
//' 
//' @references 
//'   Tipping E. (1994). WHAM--A chemical equilibrium model and computer code
//'     for waters, sediments, and soils incorporating a discrete 
//'     site/electrostatic model of ion-binding by humic substances. Computers & 
//'     Geosciences, vol. 20, iss. 6, pp. 973-1023.
//' 
//' @author Kelly Croteau (kellyc@windwardenv.com)
//' 
//' @param NSpec integer, number of species
//' @param IonicStrength double, {text}
//' @param SpecCtoM numeric vector, {text}
//' @param SpecActCorr character vector, {text}
//' @param SpecMC integer vector, {text}
//' @param AqueousMC integer, {text}
//' @param wMolWt numeric vector, {text}
//' @param wRadius numeric vector, {text}
//' @param wDLF double, {text}
//' @param wKZED double, {text}
//' @param WHAMSpecCharge numeric vector, {text}
//' @param SolHS numeric vector, {text}
//' 
//' @return Rcpp::NumericVector - a modified SpecCtoM
//' 
Rcpp::NumericVector CalcDonnanLayerVolume(unsigned int NSpec,
                                          double IonicStrength,
                                          Rcpp::NumericVector SpecCtoM,
                                          Rcpp::CharacterVector SpecActCorr,
                                          Rcpp::IntegerVector SpecMC,
                                          int AqueousMC,
                                          Rcpp::NumericVector wMolWt,
                                          Rcpp::NumericVector wRadius,
                                          double wDLF,
                                          double wKZED,
                                          Rcpp::NumericVector WHAMSpecCharge,
                                          Rcpp::NumericVector SolHS) {

  /* outputs */
  Rcpp::NumericVector DonnanLayerSpecCtoM = SpecCtoM;
  
  /* constants */
  const double Avogadro = 6.022E+23;
  const double pi = 3.14159265358979323846;
  
  /* variables */
  unsigned int iSpec, iHS;  
  double IKappa = 0.000000000304 / pow(IonicStrength, 0.5);
  double M, r, Zmod, VDmax, Tmp, Total_VD_CtoM;
  Rcpp::NumericVector VDP(2);
  Rcpp::NumericVector VD_CtoM(2);
  
  // calculate the max diffuse layer volumes
  for (iHS = 0; iHS < 2; iHS++) {
    r = wRadius(iHS);
    M = wMolWt(iHS);
    Zmod = abs(WHAMSpecCharge(iHS));  
    VDmax = 1000 * (Avogadro / M) * (4 * pi / 3) * (pow(r + IKappa, 3) - pow(r, 3));
    VDP(iHS) = (VDmax * wKZED * Zmod) / (1 + wKZED * Zmod);
    VDP(iHS) *= SolHS(iHS) * wMolWt(iHS);
  }
  
  // adjust for diffuse layer overlap with wDLF
  Tmp = wDLF / (wDLF + VDP(iHA) + VDP(iFA));
  VD_CtoM = Tmp * VDP;
  
  // Volume should never be 0. We're setting a minimum value to avoid
  // numerical issues.
  if (VD_CtoM(iHA) < 0.0) { VD_CtoM(iHA) = 0.0; }
  if (VD_CtoM(iFA) < 0.0) { VD_CtoM(iFA) = 0.0; }

  Total_VD_CtoM = VD_CtoM(iHA) + VD_CtoM(iFA);
  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    if (SpecActCorr(iSpec) == "DonnanHA") {
      DonnanLayerSpecCtoM(iSpec) = VD_CtoM(iHA);
    } else if (SpecActCorr(iSpec) == "DonnanFA") {
      DonnanLayerSpecCtoM(iSpec) = VD_CtoM(iFA);
    } else if (SpecMC(iSpec) == AqueousMC) {
      /* QUESTION:
      Should we be checking that this is an inorganic species with 
      ((SpecActCorr(iSpec) != "WHAMHA") && (SpecActCorr(iSpec) != "WHAMFA"))?
      */
      DonnanLayerSpecCtoM(iSpec) = SpecCtoM(iSpec) - Total_VD_CtoM;
    }
  }

  return DonnanLayerSpecCtoM;

}


