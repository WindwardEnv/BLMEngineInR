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

#include <math.h>
#include <Rcpp.h>
#include "CHESSFunctions.h"

//' @title AdjustForWHAMBeforeCalcSpecies
//' 
//' @description Make any necessary adjustments to the simulation that involve
//'   WHAM organic matter. The adjustments in this functionhappen before 
//'   CalcSpecConc is called
//' 
//' @details 
//'   This function does several tasks necessary when simulating WHAM
//'   organic matter. It will:
//'     1. Adjust the organic matter binding constants for ionic strength
//'        effects,
//'     2. Calculate the volume of the Donnan layer(s) (i.e., diffuse binding
//'        layers),
//' 
//' @author Kelly Croteau (kellyc@windwardenv.com)
//' 
//' @param NMass int, the number of mass compartments
//' @param MassAmt NumericVector, The amount of each mass compartment.
//' @param MassAmtAdj (OUTPUT), NumericVector, The amount of each mass 
//'   compartment, adjusted for the Donnan Layer volumes (aqueous compartment 
//'   will be its original value minus the total volume of the Donnan Layers).
//' @param NSpec int, the number of chemical species for which we have formation
//'   reactions in the simulation
//' @param SpecMC IntegerVector, length NSpec, the mass compartment of the
//'   chemical species for which we have formation reactions
//' @param SpecActCorr CharacterVector, length NSpec, the activity correction
//'   method of the chemical species for which we have formation reactions
//' @param SpecCharge IntegerVector, length NSpec, the charge of the chemical
//'   species for which we have formation reactions
//' @param SpecKTempAdj NumericVector, length NSpec, the equilibrium coefficient
//'   of the formation reactions, adjusted for temperature only
//' @param SpecKISTempAdj (OUTPUT) NumericVector, length NSpec, the equilibrium
//'   coefficient of the formation reactions, adjusted for temperature and ionic
//'   strength effects
//' @param SpecCtoMAdj (OUTPUT) NumericVector, length NSpec, the
//'   concentration to mass conversion factor of the chemical species for which
//'   we have formation reactions, adjusted to take into account diffuse binding
//'   to organic matter.
//' @param WHAMIonicStrength double, the ionic strength of the solution
//' @param WHAMSpecCharge NumericVector, length 2, the net humic charge for
//'   humic (0) and fulvic (1) substances.
//' @param AqueousMC int, the mass compartment corresponding to water/aqueous
//' @param WHAMDonnanMC IntegerVector, length 2, the mass compartments
//'   corresponding to the humic acid (0) and fulvic acid (1) Donnan layers.
//' @param HumicSubstGramsPerLiter NumericVector, length of 2, grams per liter 
//'   of each organic matter component (HA and FA) in solution
//' @param wMolWt NumericVector, WHAM's molecular weight parameter for organic
//'   matter
//' @param wRadius NumericVector, WHAM's molecular radius parameter for organic
//'   matter
//' @param wP NumericVector, WHAM's P parameter...
//' @param wDLF double, WHAM's Double layer overlap factor
//' @param wKZED double, WHAM's Constant to control DDL at low ZED
//' 
void AdjustForWHAMBeforeCalcSpecies(
  int NMass,   
  Rcpp::NumericVector MassAmt,
  Rcpp::NumericVector &MassAmtAdj,
  int NSpec,
  Rcpp::CharacterVector SpecType,
  Rcpp::IntegerVector SpecMC,
  Rcpp::IntegerVector SpecCharge,
  Rcpp::NumericVector SpecKTempAdj,
  Rcpp::NumericVector &SpecKISTempAdj,
  Rcpp::NumericVector &SpecCtoMAdj,  
  double WHAMIonicStrength,
  Rcpp::NumericVector WHAMSpecCharge,
  int AqueousMC,
  Rcpp::IntegerVector WHAMDonnanMC,
  Rcpp::NumericVector HumicSubstGramsPerLiter,
  Rcpp::NumericVector wMolWt,
  Rcpp::NumericVector wRadius,
  Rcpp::NumericVector wP,
  double wDLF,
  double wKZED
) {

  int iSpec;

  //Adjust organic matter specific binding based on ionic strength
  SpecKISTempAdj = CalcIonicStrengthEffects(WHAMIonicStrength, WHAMSpecCharge, 
                                            NSpec, SpecCharge, SpecKTempAdj, 
                                            SpecType, wP);

  //Calculate the portion of the solution that's in the diffuse layer
  MassAmtAdj = CalcDonnanLayerVolume(NMass, NSpec, WHAMIonicStrength, MassAmt, 
                                     AqueousMC, WHAMDonnanMC, wMolWt, 
                                     wRadius, wDLF, wKZED, WHAMSpecCharge, 
                                     HumicSubstGramsPerLiter);
  //SpecCtoMAdj = MassAmtAdj[SpecMC];
  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    if ((SpecType[iSpec] == STYPE_WHAMHA) || (SpecType[iSpec] == STYPE_WHAMFA)) {
      SpecCtoMAdj[iSpec] = MassAmt[AqueousMC];
    } else {
      SpecCtoMAdj[iSpec] = MassAmtAdj[SpecMC[iSpec]];
    }
  }

  // Calculate Donnan Layer concentrations


}

//' @title AdjustForWHAMAfterCalcSpecies
//' 
//' @description Make any necessary adjustments to the simulation that involve
//'   WHAM organic matter. The adjustments in this function should happen after
//'   CalcSpecConc is called.
//' 
//' @details 
//'   This function does several tasks necessary when simulating WHAM
//'   organic matter. It will:
//'     1. Set the species concentrations of ions with the same charge as the
//'        net humic charge that are binding diffusely equal to zero (only 
//'        counterions should bind diffusely).
//'     2. Calculate the charge on the WHAM organic matter species,
//'     3. Make the totals of the Donnan components equal to the net humic
//'        charge to be balanced (so that Newton-Rapshon will work on it), and
//' 
//' @author Kelly Croteau (kellyc@windwardenv.com)
//' 
//' @param NComp int, the total number of components in the simulation,
//'   including the input components, defined components (and including the
//'   defined components that get added by ExpandWHAM)
//' @param CompType CharacterVector, length of NComp, the type of each component
//'   in the simulation
//' @param TotConc (INPUT & OUTPUT) NumericVector, length of NComp, the total
//'   concentrations of each component in the simulation (units of e.g., mol/L
//'   and mol/kg)
//' @param TotMoles (INPUT & OUTPUT) NumericVector, the total moles of each
//'   component in the simulation (units of mol)
//' @param NSpec int, the number of chemical species for which we have formation
//'   reactions in the simulation
//' @param SpecConc (INPUT & OUTPUT) NumericVector, length NSpec, the
//'   concentrations of each species for which we have formation reactions
//' @param SpecMC IntegerVector, length NSpec, the mass compartment of the
//'   chemical species for which we have formation reactions
//' @param SpecActCorr CharacterVector, length NSpec, the activity correction
//'   method of the chemical species for which we have formation reactions
//' @param SpecCharge IntegerVector, length NSpec, the charge of the chemical
//'   species for which we have formation reactions
//' @param SpecCtoMAdj NumericVector, length NSpec, the concentration to mass
//'   conversion factor of the chemical species for which we have formation
//'   reactions, adjusted to take into account diffuse binding to organic
//'   matter.
//' @param WHAMSpecCharge (OUTPUT) NumericVector, length 2, the net humic charge
//'   for humic (0) and fulvic (1) substances.
//' @param AqueousMC int, the mass compartment corresponding to water/aqueous
//' @param HumicSubstGramsPerLiter NumericVector, length of 2, grams per liter 
//'   of each organic matter component (HA and FA) in solution
//' 
void AdjustForWHAMAfterCalcSpecies(int NComp,
                                   Rcpp::CharacterVector CompType,
                                   Rcpp::NumericVector &TotConc,
                                   Rcpp::NumericVector &TotMoles,
                                   int NSpec,
                                   Rcpp::CharacterVector SpecName,
                                   Rcpp::CharacterVector SpecType,
                                   Rcpp::NumericVector &SpecConc,
                                   Rcpp::NumericVector SpecKISTempAdj,
                                   Rcpp::IntegerMatrix SpecStoich,
                                   Rcpp::NumericVector SpecActivityCoef,
                                   Rcpp::IntegerVector SpecMC,
                                   Rcpp::IntegerVector SpecCharge,
                                   Rcpp::NumericVector SpecCtoMAdj,  
                                   Rcpp::NumericVector &WHAMSpecCharge,
                                   int AqueousMC,
                                   Rcpp::NumericVector HumicSubstGramsPerLiter,
                                   bool UpdateZED,
                                   bool DoWHAMSimpleAdjust,
                                   bool DoDonnanSimpleAdjust,
                                   double ConvCrit,
                                   int MaxIter) {

  /* variables */
  //const double ConvCrit = 0.000001;
  //const int MaxIter = 100;
  int iComp;
  Rcpp::NumericVector NewWHAMSpecCharge(2);
  Rcpp::NumericVector CompConc(NComp);
  for (iComp = 0; iComp < NComp; iComp++) { CompConc(iComp) = SpecConc(iComp); }
  
  if (DoWHAMSimpleAdjust) {
    // Update the WHAM component concentrations
    for (iComp = 0; iComp < NComp; iComp++) {
      if ((CompType(iComp) == CTYPE_WHAMHA) || (CompType(iComp) == CTYPE_WHAMFA)) {
        SimpleAdjustComp(iComp, ConvCrit, MaxIter, TotMoles(iComp), 
                        NComp, CompConc,
                        NSpec, SpecConc, SpecKISTempAdj, SpecStoich, SpecName,
                        SpecType, SpecActivityCoef, SpecCtoMAdj, SpecCharge,
                        WHAMSpecCharge, true);
      }
    }
  }

  //Calculate the charge on the organic matter
  if (UpdateZED) {
    NewWHAMSpecCharge = CalcWHAMSpecCharge(NSpec, SpecType, SpecConc, 
                                           SpecCharge, SpecMC, AqueousMC, 
                                           HumicSubstGramsPerLiter);
    WHAMSpecCharge(iHA) = WHAMSpecCharge(iHA) + ((NewWHAMSpecCharge(iHA) - WHAMSpecCharge(iHA)) / 5);
    WHAMSpecCharge(iFA) = WHAMSpecCharge(iFA) + ((NewWHAMSpecCharge(iFA) - WHAMSpecCharge(iFA)) / 5);
  }

  // Set -Z_HS*Th to the "known" total for the Donnan component
  for (iComp = 0; iComp < NComp; iComp++) {
    if (CompType(iComp) == CTYPE_DONNANHA) {
      TotMoles[iComp] = std::fabs(WHAMSpecCharge[iHA]) * HumicSubstGramsPerLiter[iHA];
      TotConc[iComp] = TotMoles[iComp] / SpecCtoMAdj[iComp];
    } else if (CompType(iComp) == CTYPE_DONNANFA) {
      TotMoles[iComp] = std::fabs(WHAMSpecCharge[iFA]) * HumicSubstGramsPerLiter[iFA];
      TotConc[iComp] = TotMoles[iComp] / SpecCtoMAdj[iComp];
    }
  }

  if (DoDonnanSimpleAdjust){
    // Update the Donnan layers the way WHAM does it
    for (iComp = 0; iComp < NComp; iComp++) {
      if ((CompType(iComp) == CTYPE_DONNANHA) || (CompType(iComp) == CTYPE_DONNANFA)) {
        SimpleAdjustComp(iComp, ConvCrit, MaxIter, TotMoles(iComp), 
                        NComp, CompConc,
                        NSpec, SpecConc, SpecKISTempAdj, SpecStoich, SpecName,
                        SpecType, SpecActivityCoef, SpecCtoMAdj, SpecCharge,
                        WHAMSpecCharge, true);
      }
    }
  }
  

}

