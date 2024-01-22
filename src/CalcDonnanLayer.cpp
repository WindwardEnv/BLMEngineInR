#include <Rcpp.h>
#include <cmath>
#include "CHESSFunctions.h"

Rcpp::NumericVector CalcDonnanLayer (unsigned int NSpec,
                                     double IonicStrength,
                                     Rcpp::NumericVector SpecCtoM,
                                     Rcpp::IntegerVector SpecMC,
                                     int AqueousMC,
                                     Rcpp::IntegerVector DonnanMC,
                                     Rcpp::NumericVector wMolWt,
                                     Rcpp::NumericVector wRadius,
                                     double wDLF,
                                     double wKZED,
                                     Rcpp::NumericVector WHAMSpecCharge,
                                     Rcpp::NumericVector SolHS) {

  /* outputs */
  Rcpp::NumericVector DonnanLayerSpecCtoM = SpecCtoM;
  
  /* variables */
  unsigned int iSpec;
  double AvagadrosNumber = 6.022E+23;
  double IKappa = 0.000000000304 / pow(IonicStrength, 0.5);
  //double VDP_HA; //Donnan Layer volume (L) per gram of HA
  //double VD_HA_CToM;
  Rcpp::NumericVector VDP(2);
  Rcpp::NumericVector VD_CtoM(2);
  
  /*----------------------------------------------------------------*
  |  Calculate the diffuse layer volume, VD_HA/FA, adjust CToM,    |
  |  and set the equilibrium constant for DL                       |
  *----------------------------------------------------------------*/

  // humic acid
  VDP(0) = 1000 * (AvagadrosNumber / wMolWt(0)) * 4 * (3.14159 / 3);
  VDP(0) *= (pow((wRadius(0) + IKappa), 3) - pow(wRadius(0), 3));
  VDP(0) *= wKZED * abs(WHAMSpecCharge(0)) / (1 + wKZED * abs(WHAMSpecCharge(0)));
  VDP(0) *= SolHS(0) * wMolWt(0);  

  // fulvic acid
  VDP(1) = 1000 * (AvagadrosNumber / wMolWt(1)) * 4 * (3.14159 / 3);
  VDP(1) *= (pow((wRadius(1) + IKappa), 3) - pow(wRadius(1), 3));
  VDP(1) *= wKZED * abs(WHAMSpecCharge(1)) / (1 + wKZED * abs(WHAMSpecCharge(1)));
  VDP(1) *= SolHS(1) * wMolWt(1);
  
  VD_CtoM(0) = VDP(0) / (1 + VDP(0) / wDLF + VDP(1) / wDLF);
  VD_CtoM(1) = VDP(1) / (1 + VDP(0) / wDLF + VDP(1) / wDLF);

  // Volume should never be 0. We're setting a minimum value to avoid
  // numerical issues.
  if (VD_CtoM(0) < 0.0) { VD_CtoM(0) = 0.0; }
  if (VD_CtoM(1) < 0.0) { VD_CtoM(1) = 0.0; }

  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    if (SpecMC(iSpec) == AqueousMC) {
      DonnanLayerSpecCtoM(iSpec) = DonnanLayerSpecCtoM(iSpec) - VD_CtoM(0) - VD_CtoM(1);
    } else if (SpecMC(iSpec) == DonnanMC(0)) {
      DonnanLayerSpecCtoM(iSpec) = VD_CtoM(0);
    } else if (SpecMC(iSpec) == DonnanMC(1)) {
      DonnanLayerSpecCtoM(iSpec) = VD_CtoM(1);
    }
  }

  return DonnanLayerSpecCtoM;

}


