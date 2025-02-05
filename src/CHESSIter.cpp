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

//' @title CHESSIter
//' 
//' @description Perform calculations for one iteration of CHESS
//' 
//' @details
//'   This function will take a CompConcStep input and will update all of the
//'   necessary calculations to return the residual and error information that
//'   results from that step. The updated variables and arrays are also
//'   returned. This function is meant to consolidate the calculations needed
//'   for a CHESS iteration into one function, so if other methods that are not
//'   strictly Newton-Raphson will be used, they can be without repeating (and
//'   potentially messing up) the code. This also helps to make it clearer what
//'   variables are updated with each iteration so that they won't be
//'   overwritten in the main routine unless specifically asked to. 
//' 
//' @author Kelly Croteau (kellyc@windwardenv.com)
//' 
//' @param CompConcStep NumericVector, length of NComp, the number to substract
//'   from CompConc in this iteration
//' @param NMass,MassAmt mass compartment information
//' @param NComp,CompName,CompType,CompPosInSpec component information
//' @param NSpec,SpecName,SpecMC,SpecActCorr,SpecStoich species info
//' @param SpecCharge,SpecKTempAdj more species information
//' @param DoWHAM boolean, whether WHAM organic matter is in this simulation
//' @param AqueousMC,WHAMDonnanMC which mass compartments correspond to the
//'   aqueous solution and Donnan layers, respectively
//' @param HumicSubstGramsPerLiter,wMolWt,wRadius,wP,wDLF,wKZED WHAM parameters
//' @param SysTempKelvin double, the temperature of the solution, in Kelvin
//' @param DoTox boolean, is this a toxicity run?
//' @param MetalName,MetalCompNBLMetal,BLMetalSpecs,CATarget params for tox mode
//' @param MassAmtAdj (OUTPUT), NumericVector, The amount of each mass 
//'   compartment, adjusted for the Donnan Layer volumes (aqueous compartment 
//'   will be its original value minus the total volume of the Donnan Layers).
//' @param TotConc (INPUT & OUTPUT) NumericVector, length of NComp, the total
//'   concentrations of each component in the simulation (units of e.g., mol/L
//'   and mol/kg)
//' @param SpecKISTempAdj (INPUT & OUTPUT) NumericVector, length NSpec, the
//'   equilibrium coefficient of the formation reactions, adjusted for
//'   temperature and ionic strength effects
//' @param SpecCtoMAdj (INPUT & OUTPUT) NumericVector, length NSpec, the
//'   concentration to mass conversion factor of the chemical species for which
//'   we have formation reactions, adjusted to take into account diffuse binding
//'   to organic matter.
//' @param SpecConc (INPUT & OUTPUT) NumericVector, length NSpec, the
//'   concentrations of each species for which we have formation reactions
//' @param CalcTotMoles (OUTPUT) NumericVector, length of NComp, the calculated
//'   total moles of each component in the simulation
//' @param WhichMax (OUTPUT) integer, the position in the component
//'   vectors of the component with the highest absolute error
//' @param IonicStrength (OUTPUT) double, the ionic strength of the solution
//' @param WHAMIonicStrength (OUTPUT) double, the ionic strength of the 
//'   solution, excluding organic matter contribtions
//' @param Resid (OUTPUT) numeric vector (NComp), the residuals = 
//'   calculated totals - known totals
//' @param CompError (OUTPUT) numeric vector (NComp), the absolute error
//'   fraction for each component in this iteration = abs(Resid / TotMoles)
//' 
//' @return MaxError, double, the highest absolute error fraction in this
//'   iteration =max(abs(Resid / TotMoles))
//' 
double CHESSIter(
  Rcpp::NumericVector CompConcStep,
  int NMass,
  Rcpp::NumericVector MassAmt,
  int NComp,
  Rcpp::CharacterVector CompName,
  Rcpp::CharacterVector CompType,
  Rcpp::IntegerVector CompPosInSpec,
  int NSpec,
  Rcpp::CharacterVector SpecName,
  Rcpp::CharacterVector SpecType,
  Rcpp::IntegerVector SpecMC,
  Rcpp::CharacterVector SpecActCorr,
  Rcpp::IntegerMatrix SpecStoich,
  Rcpp::IntegerVector SpecCharge,
  Rcpp::NumericVector SpecKTempAdj,
  bool DoWHAM,
  bool UpdateZED,
  int AqueousMC,
  Rcpp::IntegerVector WHAMDonnanMC,
  Rcpp::NumericVector HumicSubstGramsPerLiter,
  Rcpp::NumericVector wMolWt,
  Rcpp::NumericVector wRadius,
  Rcpp::NumericVector wP,
  double wDLF,
  double wKZED,
  double SysTempKelvin,
  bool DoTox,
  Rcpp::String MetalName,
  int MetalComp,
  int NBLMetal,
  Rcpp::IntegerVector BLMetalSpecs,
  double CATarget,
  Rcpp::NumericVector &MassAmtAdj,
  Rcpp::NumericVector &TotConc,
  Rcpp::NumericVector &TotMoles,
  Rcpp::NumericVector &SpecKISTempAdj,
  Rcpp::NumericVector &SpecCtoMAdj,
  Rcpp::NumericVector &SpecConc,
  Rcpp::NumericVector &SpecActivityCoef,
  Rcpp::NumericVector &CalcTotMoles,
  Rcpp::NumericVector &WHAMSpecCharge,
  int &WhichMax,
  double &IonicStrength,
  double &WHAMIonicStrength,
  double &ChargeBalance,
  Rcpp::NumericVector &Resid,
  Rcpp::NumericVector &CompError,
  bool DoWHAMSimpleAdjust,
  bool DoDonnanSimpleAdjust,
  double ConvCrit,
  int MaxIter
) {
  
  /*outputs*/
  double MaxError;
  Rcpp::NumericVector CalcTotConc(NComp);
    CalcTotConc.names() = CompName;
  Rcpp::NumericVector SpecMoles(NSpec);
    SpecMoles.names() = SpecName;
      
  /*variables*/
  Rcpp::NumericVector CompConc(NComp);
  Rcpp::NumericVector CompCtoMAdj(NComp);
  
  // Update the component free concentrations
  CompConc = SpecConc[CompPosInSpec];
  CompUpdate(NComp, CompConcStep, CompType, CompConc);
  SpecConc[CompPosInSpec] = clone(CompConc);
  UpdateFixedComps(NComp, CompType, TotConc, SpecActivityCoef, 
                  SpecConc, CompConc);
  
  // Calculate the ionic strength and activity coefficients
  SpecMoles = SpecConc * SpecCtoMAdj;
  ChargeBalance = CalcChargeBalance(NSpec, SpecMoles, SpecCharge, 
                                    SpecMC, AqueousMC);
  IonicStrength = CalcIonicStrength(NSpec, SpecConc, SpecCharge, 
                                    SpecMC, AqueousMC, SpecType, false);
  WHAMIonicStrength = CalcIonicStrength(NSpec, SpecConc, SpecCharge, 
                                    SpecMC, AqueousMC, SpecType, true);
  SpecActivityCoef = CalcActivityCoef(NSpec, SpecName, SpecActCorr, SpecCharge, 
                                      WHAMIonicStrength, SysTempKelvin);
  UpdateFixedComps(NComp, CompType, TotConc, SpecActivityCoef, 
                  SpecConc, CompConc);
  SpecMoles = SpecConc * SpecCtoMAdj;

  if (DoWHAM) {
    AdjustForWHAMBeforeCalcSpecies(NMass, MassAmt, MassAmtAdj, NSpec, SpecType,
      SpecMC, SpecCharge, SpecKTempAdj, SpecKISTempAdj, SpecCtoMAdj,
      WHAMIonicStrength, WHAMSpecCharge, AqueousMC, WHAMDonnanMC, 
      HumicSubstGramsPerLiter, wMolWt, wRadius, wP, wDLF, wKZED);
    //CompCtoMAdj = SpecCtoMAdj[CompPosInSpec];
    //TotMoles = TotConc * CompCtoM;Adj;
  }

  // Calculate the species concentrations
  SpecConc = CalcSpecConc(NComp, NSpec, CompConc, SpecKISTempAdj, SpecStoich, 
                          SpecName, SpecType, SpecActivityCoef, DoWHAM, 
                          SpecCharge, WHAMSpecCharge);

  UpdateFixedComps(NComp, CompType, TotConc, SpecActivityCoef, 
                  SpecConc, CompConc);
  
  if (DoWHAM) {
    
    AdjustForWHAMAfterCalcSpecies(NComp, CompType, TotConc, TotMoles, NSpec, 
      SpecName, SpecType, SpecConc, SpecKISTempAdj, SpecStoich, SpecActivityCoef, SpecMC,
      SpecCharge, SpecCtoMAdj, WHAMSpecCharge, AqueousMC, 
      HumicSubstGramsPerLiter, UpdateZED,
      DoWHAMSimpleAdjust, DoDonnanSimpleAdjust,
      ConvCrit, MaxIter);

  }       

  // Calculate the total moles & conc from species concentrations
  SpecMoles = SpecConc * SpecCtoMAdj;
  CalcTotMoles = CalcIterationTotalMoles(NComp, NSpec, SpecMoles, SpecStoich);
  CalcIterationTotals(NComp, NSpec, SpecConc, SpecCtoMAdj, SpecStoich,
                      CalcTotMoles, CalcTotConc);

  // Calculate the residuals and error fraction for each component
  CalcResidAndError(NComp, CalcTotMoles, TotMoles, CompType, 
                    SpecType, Resid, CompError);

  // Adjust Resid and CompError for toxicity mode
  if (DoTox) {
    AdjustForToxMode(NBLMetal, BLMetalSpecs, MetalComp, CATarget, SpecConc,
                      Resid, CompError);
  }

  // Determine which component has the highest error fraction
  MaxError = MaxCompError(NComp, CompError, WhichMax);

  return MaxError;

}
