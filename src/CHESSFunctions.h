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

#ifndef __CHESSFUNCTIONS_H__
#define __CHESSFUNCTIONS_H__

//should these be #included in the header file??
#ifndef Rcpp__Rcpp__h
  #include <Rcpp.h>
#endif
#include <string>

void CalcIterationTotals(int NComp,
                         int NSpec,
                         Rcpp::NumericVector SpecConc,
                         Rcpp::NumericVector SpecCtoM,
                         Rcpp::IntegerMatrix SpecStoich,
                         Rcpp::NumericVector &CalcTotMoles,
                         Rcpp::NumericVector &CalcTotConc);

Rcpp::NumericVector CalcIterationTotalMoles(int NComp,
                                            int NSpec,
                                            Rcpp::NumericVector SpecMoles,
                                            Rcpp::IntegerMatrix SpecStoich);

double CalcCompTotalMoles(int iComp,
                          int NSpec,
                          Rcpp::NumericVector SpecMoles,
                          Rcpp::IntegerMatrix SpecStoich);

void CalcResidAndError(int NComp,
                       Rcpp::NumericVector CalcTotMoles,
                       Rcpp::NumericVector TotMoles,
                       Rcpp::CharacterVector CompType,
                       Rcpp::CharacterVector SpecType,
                       Rcpp::NumericVector &Resid,
                       Rcpp::NumericVector &CompError);

void AdjustForToxMode(int NBLMetal, 
                      Rcpp::IntegerVector BLMetalSpecs, 
                      int MetalComp,
                      double CATarget,
                      Rcpp::NumericVector SpecConc,
                      Rcpp::NumericVector &Resid,
                      Rcpp::NumericVector &CompError);
                   
double CalcCA(int NBLMetal, 
              Rcpp::IntegerVector BLMetalSpecs, 
              Rcpp::NumericVector SpecConc);

double MaxCompError(int NComp, Rcpp::NumericVector CompError, 
                    int &WhichMax);
                    
Rcpp::NumericVector CalcSpecConc(int NComp,
                                 int NSpec,
                                 Rcpp::NumericVector CompConc,
                                 Rcpp::NumericVector SpecK,
                                 Rcpp::IntegerMatrix SpecStoich,
                                 Rcpp::CharacterVector SpecName,
                                 Rcpp::CharacterVector SpecType,
                                 Rcpp::NumericVector SpecActivityCoef,
                                 bool DoWHAM,
                                 Rcpp::IntegerVector SpecCharge,
                                 Rcpp::NumericVector WHAMSpecCharge);

Rcpp::NumericVector CalcStep(Rcpp::NumericMatrix JacobianMatrix,
                             Rcpp::NumericVector Resid,
                             Rcpp::NumericVector CompConc,
                             Rcpp::NumericVector TotMoles,
                             Rcpp::NumericVector CalcTotMoles,
                             int NComp,
                             Rcpp::CharacterVector CompType,
                             Rcpp::CharacterVector CompName);

void CompUpdate(int NComp, 
                Rcpp::NumericVector CompConcStep,
                Rcpp::CharacterVector CompType,
                Rcpp::NumericVector &CompConc);

void SimpleAdjustComp(int iComp,
                      double ConvCrit,
                      int MaxIter,
                      double TotMolesi,
                      int NComp,
                      Rcpp::NumericVector &CompConc,
                      int NSpec,
                      Rcpp::NumericVector &SpecConc,
                      Rcpp::NumericVector SpecKISTempAdj,
                      Rcpp::IntegerMatrix SpecStoich,
                      Rcpp::CharacterVector SpecName,
                      Rcpp::CharacterVector SpecType,
                      Rcpp::NumericVector SpecActivityCoef,
                      Rcpp::NumericVector SpecCtoMAdj,
                      Rcpp::IntegerVector SpecCharge,
                      Rcpp::NumericVector WHAMSpecCharge,
                      bool DoWHAM);

Rcpp::NumericVector InitialGuess(Rcpp::NumericVector &TotConc,
                                 Rcpp::NumericVector SpecCtoM, 
                                 Rcpp::CharacterVector CompType,
									               Rcpp::NumericVector SpecK,
                                 Rcpp::IntegerMatrix SpecStoich,
                                 Rcpp::CharacterVector SpecName,
                                 int NComp,
                                 int NSpec, 
                                 bool DoTox,
                                 int NBLMetal, 
                                 Rcpp::IntegerVector BLMetalSpecs, 
                                 int MetalComp,
                                 double CATarget);

Rcpp::NumericMatrix Jacobian(int NComp, //number of components
                              int NSpec, //number of species
                              Rcpp::CharacterVector CompType,
                              Rcpp::IntegerMatrix SpecStoich, //formation reaction stoichiometry (NSpec x NComp)
                              Rcpp::NumericVector SpecConc, //species concentrations
                              Rcpp::IntegerVector SpecMC,
                              Rcpp::NumericVector SpecCtoM, //concentration to mass conversion for each species
                              Rcpp::CharacterVector SpecType,
                              Rcpp::IntegerVector SpecCharge,
                              Rcpp::NumericVector SpecK,
                              Rcpp::NumericVector SpecActivityCoef,
                              double WHAMIonicStrength,
                              bool DoWHAM,
                              Rcpp::NumericVector HumicSubstGramsPerLiter,
                              Rcpp::NumericVector WHAMSpecCharge,
                              Rcpp::NumericVector wP,
                              Rcpp::NumericVector wMolWt,
                              Rcpp::NumericVector wRadius,
                              double wDLF,
                              double wKZED,
                              int AqueousMC,
                              int MetalComp, //position of the metal component
                              int BLComp,//position of BL component
                              int NBLMetal, //number of BL-Metal species
                              Rcpp::IntegerVector BLMetalSpecs, //positions of BL-metal species
                              bool DoTox,
                 bool DodVidCj,
                 bool DodVidCjDonnan,
                 bool DodKidCj,
                 bool DoGammai,
                 bool DoJacDonnan,
                 bool DoJacWHAM);

void UpdateFixedComps(int NComp, 
                      Rcpp::CharacterVector CompType, 
                      Rcpp::NumericVector InCompConc, 
                      Rcpp::NumericVector SpecActivityCoef,
                      Rcpp::NumericVector &SpecConc, 
                      Rcpp::NumericVector &CompConc);

Rcpp::NumericVector TempCorrection(double SysTempKelvin,
                                   int NSpec,
                                   Rcpp::NumericVector SpecK,
                                   Rcpp::NumericVector SpecTempKelvin,
                                   Rcpp::NumericVector SpecDeltaH);

double CalcIonicStrength(int NSpec,
                         Rcpp::NumericVector SpecConc,
                         Rcpp::IntegerVector SpecCharge,
                         Rcpp::IntegerVector SpecMC,
                         int AqueousMC, 
                         Rcpp::CharacterVector SpecType,
                         bool ExcludeOrgMatter);

double CalcChargeBalance(int NSpec,
                         Rcpp::NumericVector SpecMoles,
                         Rcpp::IntegerVector SpecCharge,
                         Rcpp::IntegerVector SpecMC,
                         int AqueousMC);

Rcpp::NumericVector CalcActivityCoef(int NSpec,
                                     Rcpp::CharacterVector SpecName, 
                                     Rcpp::CharacterVector SpecActCorr,
                                     Rcpp::IntegerVector SpecCharge,
                                     double IonicStrength,
                                     double SysTempKelvin);

Rcpp::NumericVector CalcWHAMSpecCharge(int NSpec, 
                                       Rcpp::CharacterVector SpecType,
                                       Rcpp::NumericVector SpecConc,
                                       Rcpp::IntegerVector SpecCharge,
                                       Rcpp::IntegerVector SpecMC,
                                       int AqueousMC,
                                       Rcpp::NumericVector HumicSubstGramsPerLiter);

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
                           Rcpp::NumericVector &VolDiffuse);

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
                                          Rcpp::NumericVector HumicSubstGramsPerLiter);

Rcpp::NumericVector CalcIonicStrengthEffects(double WHAMIonicStrength,
                                             Rcpp::NumericVector WHAMSpecCharge,
                                             int NSpec,
                                             Rcpp::IntegerVector SpecCharge,
                                             Rcpp::NumericVector SpecK,
                                             Rcpp::CharacterVector SpecType,
                                             Rcpp::NumericVector wP);

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
);

void AdjustForWHAMAfterCalcSpecies(
  int NComp,
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
  int MaxIter
);

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
);

// humic substance position indices
const int iHA = 0; // humic acid = 0
const int iFA = 1; // fulvic acid = 1

// console feedback flags
const Rcpp::String FLAG_DEBUG = "Debug";
const Rcpp::String FLAG_QUIET = "Quiet";

//Component Types
const Rcpp::String CTYPE_FIXEDCONC = "FixedConc";
const Rcpp::String CTYPE_FIXEDACT = "FixedAct";
const Rcpp::String CTYPE_MASSBAL = "MassBal";
const Rcpp::String CTYPE_DONNANHA = "DonnanHA";
const Rcpp::String CTYPE_DONNANFA = "DonnanFA";
const Rcpp::String CTYPE_WHAMHA = "WHAMHA";
const Rcpp::String CTYPE_WHAMFA = "WHAMFA";

//Species types
const Rcpp::String STYPE_NORMAL = "Normal";
const Rcpp::String STYPE_DONNANHA = "DonnanHA";
const Rcpp::String STYPE_DONNANFA = "DonnanFA";
const Rcpp::String STYPE_WHAMHA = "WHAMHA";
const Rcpp::String STYPE_WHAMFA = "WHAMFA";

//Activity correction methods
const Rcpp::String ACTYPE_NONE = "None";
const Rcpp::String ACTYPE_DEBYE = "Debye";
const Rcpp::String ACTYPE_DAVIES = "Davies";

// error type constants
const int ERROR_MATRIX_INVERSION = 10;
const int ERROR_SINGULAR_MATRIX = 11;
const int ERROR_JACOBIAN_NAN = 20;

// status messages
const Rcpp::String STATUS_SPEC_DNC = "CHESS did not converge.";
const Rcpp::String STATUS_JAC_ERR = "Jacobian matrix error.";
const Rcpp::String STATUS_MAT_ERR = "Matrix inversion error.";

#endif //__CHESSFUNCTIONS_H__