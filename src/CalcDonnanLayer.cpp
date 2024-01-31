#include <cmath>
#include <Rcpp.h>
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
  
  /* constants */
  const double Avogadro = 6.022E+23;
  const double pi = 3.14159265358979323846;
  const unsigned int iHA = 0; // humic acid = 0
  const unsigned int iFA = 1; // fulvic acid = 1
  
  /* variables */
  unsigned int iSpec, iHS;  
  double IKappa = 0.000000000304 / pow(IonicStrength, 0.5);
  double M, r, Zmod, VDmax, Tmp, Total_VD_CtoM;
  //double VDP_HA; //Donnan Layer volume (L) per gram of HA
  //double VD_HA_CToM;
  Rcpp::NumericVector VDP(2);
  Rcpp::NumericVector VD_CtoM(2);
  
  /*----------------------------------------------------------------*
  |  Calculate the diffuse layer volume, VD_HA/FA, adjust CToM,    |
  |  and set the equilibrium constant for DL                       |
  *----------------------------------------------------------------*/

  /*
  From Tipping 1993, the maximum diffuse layer volume, VDmax is given by 
    VDmax = 1000 * (N / M) * (4 * pi / 3) * [(r + 1 / kappa) ^ 3 - r ^ 3]
  where: N is Avogadro's number
         M is the molecular weight
         r is the radius of the organic molecule
         kappa is ... the Debye-Huckel characteristic distance?
  VDmax should be modified to prevent the situation where the counterion charge
  predicted to be lower than in the bulk solution when the net humic charge is
  low:
    VDmaxbar = (VDmax * KZED * Zmod) / (1 + KZED * Zmod)
  where KZED is a non-optimized value, usually set to 1000
        Zmod is the modulus of Z, the net humic charge, i.e., abs(Z)
  The actual diffuse layer volume is calculated as:
    VD(FA) = [fDL * VDmaxbar(FA)] / [fDL + VDmaxbar(FA) + VDmaxbar(HA)]
           = VDmaxbar(FA) / [1 + VDmaxbar(FA) / fDL + VDmaxbar(HA) / fDL]
    VD(HA) = [fDL * VDmaxbar(HA)] / [fDL + VDmaxbar(FA) + VDmaxbar(HA)]
  where fDL is set to a number between 1 and 0, that is the asymptotic max total
  diffuse layer volume...so wDLF...
  */


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
  //VD_CtoM(0) = Tmp * VDP(0);
  //VD_CtoM(1) = Tmp * VDP(1);
  
  // Volume should never be 0. We're setting a minimum value to avoid
  // numerical issues.
  if (VD_CtoM(iHA) < 0.0) { VD_CtoM(iHA) = 0.0; }
  if (VD_CtoM(iFA) < 0.0) { VD_CtoM(iFA) = 0.0; }

  Total_VD_CtoM = VD_CtoM(iHA) + VD_CtoM(iFA);
  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    if (SpecMC(iSpec) == AqueousMC) {
      DonnanLayerSpecCtoM(iSpec) = DonnanLayerSpecCtoM(iSpec) - Total_VD_CtoM;
    } else if (SpecMC(iSpec) == DonnanMC(iHA)) {
      DonnanLayerSpecCtoM(iSpec) = VD_CtoM(iHA);
    } else if (SpecMC(iSpec) == DonnanMC(iFA)) {
      DonnanLayerSpecCtoM(iSpec) = VD_CtoM(iFA);
    }
  }

  return DonnanLayerSpecCtoM;

}


