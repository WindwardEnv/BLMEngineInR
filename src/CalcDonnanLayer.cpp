// Copyright 2024 Windward Environmental LLC
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <Rcpp.h>
#include "CHESSFunctions.h"

void CalcDonnanLayerParams(int NSpec,
                           double WHAMIonicStrength,
                           Rcpp::NumericVector wMolWt,
                           Rcpp::NumericVector wRadius,
                           double wDLF,
                           double wKZED,
                           Rcpp::NumericVector WHAMSpecCharge,
                           Rcpp::NumericVector HumicSubstGramsPerLiter,
                           Rcpp::NumericVector &MaxVolDiffusePerGramHS,
                           Rcpp::NumericVector &MaxVolDiffuse,
                           Rcpp::NumericVector &VolDiffuse) {
  
  /* outputs */
  
  /* constants */
  const double Avogadro = 6E+23;//6.022E+23;
  const double pi = 3.1425;//3.14159265358979323846;
  
  /* variables */
  Rcpp::NumericVector VTerm1(2), VTerm2(2), VTerm3(2), ZTerm(2);//intermediates
  double Denom;
  
  // calculate the max diffuse layer volumes  
  VTerm1 = wRadius + (3.04E-10 / sqrt(WHAMIonicStrength));
  VTerm2 = pow(VTerm1, 3) - pow(wRadius, 3);
  VTerm3 = (4 * pi / 3) * VTerm2;
  MaxVolDiffusePerGramHS = Avogadro * VTerm3 * (1000 / wMolWt);
  MaxVolDiffuse = MaxVolDiffusePerGramHS * HumicSubstGramsPerLiter;
  ZTerm = wKZED * Rcpp::abs(WHAMSpecCharge);
  MaxVolDiffuse = MaxVolDiffuse * ZTerm / (1 + ZTerm);
  
  // adjust for diffuse layer overlap with wDLF
  Denom = 1 + ((MaxVolDiffuse(0) + MaxVolDiffuse(1)) / wDLF);
  VolDiffuse = MaxVolDiffuse / Denom;
  
  // Volume should never be 0. We're setting a minimum value to avoid
  // numerical issues.
  if (VolDiffuse(iHA) < 0.0) { VolDiffuse(iHA) = 0.0; }
  if (VolDiffuse(iFA) < 0.0) { VolDiffuse(iFA) = 0.0; }

}


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
//' @param WHAMIonicStrength double, {text}
//' @param SpecCtoM numeric vector, {text}
//' @param SpecMC integer vector, {text}
//' @param AqueousMC integer, {text}
//' @param wMolWt numeric vector, {text}
//' @param wRadius numeric vector, {text}
//' @param wDLF double, {text}
//' @param wKZED double, {text}
//' @param WHAMSpecCharge numeric vector, {text}
//' @param HumicSubstGramsPerLiter NumericVector, length of 2, grams per liter 
//'   of each organic matter component (HA and FA) in solution
//' 
//' @return Rcpp::NumericVector - a modified SpecCtoM
//' 
Rcpp::NumericVector CalcDonnanLayerVolume(int NMass,
                                          int NSpec,
                                          double WHAMIonicStrength,
                                          Rcpp::NumericVector MassAmt,
                                          int AqueousMC,
                                          Rcpp::IntegerVector WHAMDonnanMC,
                                          Rcpp::NumericVector wMolWt,
                                          Rcpp::NumericVector wRadius,
                                          double wDLF,
                                          double wKZED,
                                          Rcpp::NumericVector WHAMSpecCharge,
                                          Rcpp::NumericVector HumicSubstGramsPerLiter) {

  /* outputs */
  Rcpp::NumericVector MassAmtAdj(NMass);// = SpecCtoM;
    MassAmtAdj.names() = MassAmt.names();
  
  /* variables */
  int iMass;
  Rcpp::NumericVector MaxVolDiffusePerGramHS(2);
  Rcpp::NumericVector MaxVolDiffuse(2);
  Rcpp::NumericVector VolDiffuse(2);
  double VolSolution;  

  CalcDonnanLayerParams(NSpec, WHAMIonicStrength, wMolWt, wRadius, wDLF, wKZED, 
                        WHAMSpecCharge, HumicSubstGramsPerLiter, 
                        MaxVolDiffusePerGramHS, MaxVolDiffuse, VolDiffuse);
  
  // Update the Aqueous and Donnan mass compartments with the new volumes
  VolSolution = MassAmt[AqueousMC] - VolDiffuse(iHA) - VolDiffuse(iFA);

  /*Rcpp::Rcout << "VolSolution = " << VolSolution 
              << ", VolDiffuse[iHA] = " << VolDiffuse[iHA]              
              << ", VolDiffuse[iFA] = " << VolDiffuse[iFA]
              << std::endl;*/

  for (iMass = 0; iMass < NMass; iMass++) {
    if (iMass == AqueousMC) {
      MassAmtAdj[iMass] = VolSolution;
    } else if (iMass == WHAMDonnanMC[iHA]) {
      MassAmtAdj[iMass] = VolDiffuse[iHA];
    } else if (iMass == WHAMDonnanMC[iFA]) {
      MassAmtAdj[iMass] = VolDiffuse[iFA];
    } else {
      MassAmtAdj[iMass] = MassAmt[iMass];
    }
  }
  /* QUESTION:
      Should we be checking that this is an inorganic species with 
      ((SpecType(iSpec) != STYPE_WHAMHA) && (SpecType(iSpec) != STYPE_WHAMFA))?

      (KEC, 2024-03-07):
      I looked into the Tipping code, and it seems that VolSol is only used when
      dealing with calculating the totals of inorganic species, not the DOC
      species. However, I think it's useful to have the MassAmt reflect the
      different volumes, so I changed ExpandWHAM so that DOC species are added
      to a Water_Bulk mass compartment, which should be a copy of the Water MC
      that doesn't get changes to MassAmt.
      */

  return MassAmtAdj;

}

