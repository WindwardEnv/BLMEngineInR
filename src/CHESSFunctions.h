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
                       Rcpp::CharacterVector SpecActCorr,
                       Rcpp::NumericVector &Resid,
                       Rcpp::NumericVector &CompError);

void AdjustForToxMode(int NBLMetal, 
                      Rcpp::IntegerVector BLMetalSpecs, 
                      int MetalComp,
                      double CATarget,
                      Rcpp::NumericVector SpecConc,
                      Rcpp::NumericVector &Resid,
                      Rcpp::NumericVector &CompError);

void CalcToxError (int NBLMetal, 
                   Rcpp::IntegerVector BLMetalSpecs, 
                   int MetalComp,
                   double CATarget,
                   Rcpp::NumericVector SpecConc,
                   double &ToxResid,
                   double &ToxError) ;

double MaxCompError(int NComp, Rcpp::NumericVector CompError, 
                    int &WhichMax);
                    
Rcpp::NumericVector CalcSpecConc(int NComp,
                                 int NSpec,
                                 Rcpp::NumericVector CompConc,
                                 Rcpp::NumericVector SpecK,
                                 Rcpp::IntegerMatrix SpecStoich,
                                 Rcpp::CharacterVector SpecName,
                                 Rcpp::CharacterVector SpecActCorr,
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

Rcpp::NumericVector CalcStepBrute(int NComp, 
                                  Rcpp::CharacterVector CompName, 
                                  Rcpp::CharacterVector CompType, 
                                  Rcpp::NumericVector CompConc, 
                                  Rcpp::NumericVector TotMoles, 
                                  Rcpp::NumericVector CalcTotMoles);

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
                      Rcpp::CharacterVector SpecActCorr,
                      Rcpp::NumericVector SpecActivityCoef,
                      Rcpp::NumericVector SpecCtoMAdj,
                      Rcpp::IntegerVector SpecCharge,
                      Rcpp::NumericVector WHAMSpecCharge);

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

Rcpp::NumericMatrix Jacobian (int NComp, //number of components
                              int NSpec, //number of species
                              Rcpp::CharacterVector CompName, //names of components
                              Rcpp::NumericVector TotConc,
                              Rcpp::IntegerMatrix SpecStoich, //formation reaction stoichiometry (NSpec x NComp)
                              Rcpp::NumericVector SpecConc, //species concentrations
                              Rcpp::IntegerVector SpecMC,
                              Rcpp::NumericVector SpecCtoM, //concentration to mass conversion for each species
                              Rcpp::CharacterVector SpecActCorr,
                              Rcpp::IntegerVector SpecCharge,
                              Rcpp::NumericVector SpecK,
                              double IonicStrength,
                              bool DoWHAM,
                              Rcpp::NumericVector HumicSubstGramsPerLiter,
                              Rcpp::NumericVector WHAMSpecCharge,
                              Rcpp::NumericVector wP,
                              Rcpp::NumericVector wMolWt,
                              Rcpp::NumericVector wRadius,
                              double wDLF,
                              double wKZED,
                              Rcpp::NumericVector MassAmtAdj,
                              int AqueousMC,
                              Rcpp::IntegerVector WHAMDonnanMC,
                              int MetalComp, //position of the metal component
                              int NBLMetal, //number of BL-Metal species
                              Rcpp::IntegerVector BLMetalSpecs, //positions of BL-metal species
                              bool DoTox);

Rcpp::NumericMatrix NumericalJacobian(
  int NMass,
  Rcpp::NumericVector MassAmt,
  int NComp,
  Rcpp::CharacterVector CompName,
  Rcpp::CharacterVector CompType,
  Rcpp::IntegerVector CompPosInSpec,
  int NSpec,
  Rcpp::CharacterVector SpecName,
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
  Rcpp::NumericVector MassAmtAdj,
  Rcpp::NumericVector TotConc,
  Rcpp::NumericVector TotMoles,
  Rcpp::NumericVector SpecKISTempAdj,
  Rcpp::NumericVector SpecCtoMAdj,
  Rcpp::NumericVector SpecConc,
  Rcpp::NumericVector SpecActivityCoef,
  Rcpp::NumericVector WHAMSpecCharge,
  double IonicStrength,
  Rcpp::NumericVector Resid) ;

void UpdateTotals(int NComp, 
                  int NSpec, 
                  bool DoTox, 
                  Rcpp::CharacterVector CompType, 
                  Rcpp::CharacterVector CompName, 
                  Rcpp::String MetalName, 
                  Rcpp::IntegerMatrix SpecStoich, 
                  Rcpp::NumericVector SpecConc, 
                  Rcpp::NumericVector CompCtoM, 
                  Rcpp::NumericVector &TotMoles, 
                  Rcpp::NumericVector &TotConc);

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
                         Rcpp::NumericVector SpecMoles,
                         Rcpp::IntegerVector SpecCharge,
                         Rcpp::IntegerVector SpecMC,
                         int AqueousMC, 
                         Rcpp::CharacterVector SpecActCorr,
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
                                       Rcpp::CharacterVector SpecActCorr,
                                       Rcpp::NumericVector SpecConc,
                                       Rcpp::IntegerVector SpecCharge,
                                       Rcpp::IntegerVector SpecMC,
                                       int AqueousMC,
                                       Rcpp::NumericVector HumicSubstGramsPerLiter);

void CalcDonnanLayerParams(int NSpec,
                           double IonicStrength,
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
                                          double IonicStrength,
                                          Rcpp::NumericVector MassAmt,
                                          int AqueousMC,
                                          Rcpp::IntegerVector WHAMDonnanMC,
                                          Rcpp::NumericVector wMolWt,
                                          Rcpp::NumericVector wRadius,
                                          double wDLF,
                                          double wKZED,
                                          Rcpp::NumericVector WHAMSpecCharge,
                                          Rcpp::NumericVector HumicSubstGramsPerLiter);

Rcpp::NumericVector CalcIonicStrengthEffects(double IonicStrength,
                                             Rcpp::NumericVector WHAMSpecCharge,
                                             int NSpec,
                                             Rcpp::IntegerVector SpecCharge,
                                             Rcpp::NumericVector SpecK,
                                             Rcpp::CharacterVector SpecActCorr,
                                             Rcpp::NumericVector wP);

void AdjustForWHAM(
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
  Rcpp::NumericVector HumicSubstGramsPerLiter,
  Rcpp::NumericVector wMolWt,
  Rcpp::NumericVector wRadius,
  Rcpp::NumericVector wP,
  double wDLF,
  double wKZED
);

void AdjustForWHAMBeforeCalcSpecies(
  int NMass,   
  Rcpp::NumericVector MassAmt,
  Rcpp::NumericVector &MassAmtAdj,
  int NSpec,
  Rcpp::IntegerVector SpecMC,
  Rcpp::CharacterVector SpecActCorr,
  Rcpp::IntegerVector SpecCharge,
  Rcpp::NumericVector SpecKTempAdj,
  Rcpp::NumericVector &SpecKISTempAdj,
  Rcpp::NumericVector &SpecCtoMAdj,  
  double IonicStrength,
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

void AdjustForWHAMAfterCalcSpecies(int NComp,
                                   Rcpp::CharacterVector CompType,
                                   Rcpp::NumericVector &TotConc,
                                   Rcpp::NumericVector &TotMoles,
                                   int NSpec,
                                   Rcpp::CharacterVector SpecName,
                                   Rcpp::NumericVector &SpecConc,
                                   Rcpp::NumericVector SpecKISTempAdj,
                                   Rcpp::IntegerMatrix SpecStoich,
                                   Rcpp::NumericVector SpecActivityCoef,
                                   Rcpp::IntegerVector SpecMC,
                                   Rcpp::CharacterVector SpecActCorr,
                                   Rcpp::IntegerVector SpecCharge,
                                   Rcpp::NumericVector SpecCtoMAdj,  
                                   Rcpp::NumericVector &WHAMSpecCharge,
                                   int AqueousMC,
                                   Rcpp::NumericVector HumicSubstGramsPerLiter,
                                   bool UpdateZED);

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
  Rcpp::NumericVector &Resid,
  Rcpp::NumericVector &CompError
);


const int iHA = 0; // humic acid = 0
const int iFA = 1; // fulvic acid = 1

const Rcpp::String FLAG_DEBUG = "Debug";
const Rcpp::String FLAG_QUIET = "Quiet";

const Rcpp::String TYPE_FIXEDCONC = "FixedConc";
const Rcpp::String TYPE_FIXEDACT = "FixedAct";
const Rcpp::String TYPE_DONNANHA = "DonnanHA";
const Rcpp::String TYPE_DONNANFA = "DonnanFA";
const Rcpp::String TYPE_MASSBAL = "MassBal";

const Rcpp::String ACTYPE_NONE = "None";
const Rcpp::String ACTYPE_DEBYE = "Debye";
const Rcpp::String ACTYPE_DAVIES = "Davies";
const Rcpp::String ACTYPE_DONNANHA = "DonnanHA";
const Rcpp::String ACTYPE_DONNANFA = "DonnanFA";
const Rcpp::String ACTYPE_WHAMHA = "WHAMHA";
const Rcpp::String ACTYPE_WHAMFA = "WHAMFA";

const int ERROR_MATRIX_INVERSION = 10;
const int ERROR_SINGULAR_MATRIX = 11;
const int ERROR_JACOBIAN_NAN = 20;


#endif //__CHESSFUNCTIONS_H__