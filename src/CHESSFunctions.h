#ifndef __CHESSFUNCTIONS_H__
#define __CHESSFUNCTIONS_H__

//should these be #included in the header file??
#ifndef Rcpp__Rcpp__h
  #include <Rcpp.h>
#endif
#include <string>

Rcpp::List CalcResidualList (int NComp,
                             int NSpec,
                             Rcpp::NumericVector SpecConc,
                             Rcpp::IntegerMatrix SpecStoich,
                             Rcpp::NumericVector TotMoles,
                             Rcpp::NumericVector SpecCtoM,
                             Rcpp::CharacterVector CompName,
                             Rcpp::CharacterVector CompType,
                             int MetalComp,
                             int NBLMetal,
                             Rcpp::IntegerVector BLMetalSpecs,
                             double CATarget,
                             bool DoTox);

void CalcIterationTotals(int NComp,
                         int NSpec,
                         Rcpp::NumericVector SpecConc,
                         Rcpp::NumericVector SpecCtoM,
                         Rcpp::IntegerMatrix SpecStoich,
                         Rcpp::NumericVector &CalcTotMoles,
                         Rcpp::NumericVector &CalcTotConc);

Rcpp::NumericVector CalcResidualsOnly(int NComp,
                                      Rcpp::NumericVector CalcTotMoles,
                                      Rcpp::NumericVector TotMoles,
                                      Rcpp::CharacterVector CompType);

void CalcResidAndError(int NComp,
                       Rcpp::NumericVector CalcTotMoles,
                       Rcpp::NumericVector TotMoles,
                       Rcpp::CharacterVector CompType,
                       Rcpp::NumericVector &Resid,
                       Rcpp::NumericVector &CompError);

void AdjustForToxMode(int NBLMetal, 
                      Rcpp::IntegerVector BLMetalSpecs, 
                      int MetalComp,
                      double CATarget,
                      Rcpp::NumericVector SpecConc,
                      Rcpp::NumericVector &Resid,
                      Rcpp::NumericVector &CompError);

double MaxCompError(int NComp, Rcpp::NumericVector CompError, 
                    int &WhichMax);
                    
Rcpp::NumericVector CalcSpecConc(int NComp,
                                 int NSpec,
                                 Rcpp::NumericVector CompConc,
                                 Rcpp::NumericVector SpecK,
                                 Rcpp::IntegerMatrix SpecStoich,
                                 Rcpp::CharacterVector SpecName,
                                 Rcpp::CharacterVector SpecActCorr,
                                 Rcpp::NumericVector SpecActivityCoef);

Rcpp::NumericVector CalcStep(Rcpp::NumericMatrix JacobianMatrix,
                              Rcpp::NumericVector Resid,
                              int NComp,
                              Rcpp::CharacterVector CompType,
                              Rcpp::CharacterVector CompName);

void CompUpdate(int NComp, 
                Rcpp::NumericVector CompConcStep,
                Rcpp::NumericVector &CompConc);

Rcpp::NumericVector InitialGuess(Rcpp::NumericVector TotConc,
                                 Rcpp::NumericVector SpecCtoM, 
                                 Rcpp::CharacterVector CompType,
									               Rcpp::NumericVector SpecK,
                                 Rcpp::IntegerMatrix SpecStoich,
                                 Rcpp::CharacterVector SpecName,
                                 int NComp,
                                 int NSpec);

Rcpp::NumericMatrix Jacobian (int NComp, //number of components
                              int NSpec, //number of species
                              Rcpp::IntegerMatrix SpecStoich, //formation reaction stoichiometry (NSpec x NComp)
                              Rcpp::NumericVector SpecConc, //species concentrations
                              Rcpp::NumericVector SpecCtoM, //concentration to mass conversion for each species
                              Rcpp::CharacterVector CompName, //names of components
                              int MetalComp, //position of the metal component
                              int NBLMetal, //number of BL-Metal species
                              Rcpp::IntegerVector BLMetalSpecs, //positions of BL-metal species
                              bool DoTox);

void UpdateTotals(int NComp, 
                  int NSpec, 
                  bool DoTox, 
                  Rcpp::CharacterVector CompType, 
                  Rcpp::CharacterVector CompName, 
                  Rcpp::String MetalName, 
                  Rcpp::IntegerMatrix SpecStoich, 
                  Rcpp::NumericVector SpecMoles, 
                  Rcpp::NumericVector CompCtoM, 
                  Rcpp::NumericVector &TotMoles, 
                  Rcpp::NumericVector &TotConc);

Rcpp::List UpdateTotalsList(int NComp, 
                            int NSpec, 
                            bool DoTox, 
                            Rcpp::CharacterVector CompType, 
                            Rcpp::CharacterVector CompName, 
                            Rcpp::String MetalName, 
                            Rcpp::IntegerMatrix SpecStoich, 
                            Rcpp::NumericVector SpecMoles, 
                            Rcpp::NumericVector CompCtoM, 
                            Rcpp::NumericVector TotMoles, 
                            Rcpp::NumericVector TotConc);

Rcpp::NumericVector TempCorrection(double SysTempKelvin,
                                   int NSpec,
                                   Rcpp::NumericVector SpecK,
                                   Rcpp::NumericVector SpecTempKelvin,
                                   Rcpp::NumericVector SpecDeltaH);

double CalcIonicStrength(int NSpec,
                         Rcpp::NumericVector SpecMoles,
                         Rcpp::IntegerVector SpecCharge,
                         Rcpp::IntegerVector SpecMC,
                         int AqueousMC);

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
                                       Rcpp::NumericVector SpecMoles,
                                       Rcpp::IntegerVector SpecCharge,
                                       Rcpp::IntegerVector SpecMC,
                                       int AqueousMC);

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
                                          Rcpp::NumericVector SolHS);

Rcpp::NumericVector CalcIonicStrengthEffects(double IonicStrength,
                                             Rcpp::NumericVector WHAMSpecCharge,
                                             int NSpec,
                                             Rcpp::IntegerVector SpecCharge,
                                             Rcpp::NumericVector SpecK,
                                             Rcpp::CharacterVector SpecActCorr,
                                             Rcpp::NumericVector wP);

void AdjustForWHAM(int NComp,
                   int NSpec,
                   Rcpp::CharacterVector CompName,
                   Rcpp::CharacterVector SpecActCorr,
                   Rcpp::IntegerVector SpecCharge,
                   Rcpp::NumericVector WHAMSpecCharge,
                   Rcpp::NumericVector SpecCtoM,
                   Rcpp::NumericVector &SpecConc,
                   Rcpp::NumericVector &TotConc,
                   Rcpp::NumericVector &TotMoles);


const int iHA = 0; // humic acid = 0
const int iFA = 1; // fulvic acid = 1

#endif //__CHESSFUNCTIONS_H__