#ifndef __CHESSFUNCTIONS_H__
#define __CHESSFUNCTIONS_H__

//should these be #included in the header file??
#ifndef Rcpp__Rcpp__h
  #include <Rcpp.h>
#endif
#include <string>

Rcpp::List CalcResidualList (unsigned int NComp,
                             unsigned int NSpec,
                             Rcpp::NumericVector SpecConc,
                             Rcpp::IntegerMatrix SpecStoich,
                             Rcpp::NumericVector TotMoles,
                             Rcpp::NumericVector SpecCtoM,
                             Rcpp::CharacterVector CompName,
                             Rcpp::CharacterVector CompType,
                             unsigned int MetalComp,
                             unsigned int NBLMetal,
                             Rcpp::IntegerVector BLMetalSpecs,
                             double CATarget,
                             bool DoTox);

void CalcIterationTotals(unsigned int NComp,
                         unsigned int NSpec,
                         Rcpp::NumericVector SpecConc,
                         Rcpp::NumericVector SpecCtoM,
                         Rcpp::IntegerMatrix SpecStoich,
                         Rcpp::NumericVector &CalcTotMoles,
                         Rcpp::NumericVector &CalcTotConc);

Rcpp::NumericVector CalcResidualsOnly(unsigned int NComp,
                                      Rcpp::NumericVector CalcTotMoles,
                                      Rcpp::NumericVector TotMoles,
                                      Rcpp::CharacterVector CompType);

void CalcResidAndError(unsigned int NComp,
                       Rcpp::NumericVector CalcTotMoles,
                       Rcpp::NumericVector TotMoles,
                       Rcpp::CharacterVector CompType,
                       Rcpp::NumericVector &Resid,
                       Rcpp::NumericVector &CompError);

void AdjustForToxMode(unsigned int NBLMetal, 
                      Rcpp::IntegerVector BLMetalSpecs, 
                      unsigned int MetalComp,
                      double CATarget,
                      Rcpp::NumericVector SpecConc,
                      Rcpp::NumericVector &Resid,
                      Rcpp::NumericVector &CompError);

double MaxCompError(unsigned int NComp, Rcpp::NumericVector CompError, 
                    unsigned int &WhichMax);
                    
Rcpp::NumericVector CalcSpecConc(unsigned int NComp,
                                 unsigned int NSpec,
                                 Rcpp::NumericVector CompConc,
                                 Rcpp::NumericVector SpecK,
                                 Rcpp::IntegerMatrix SpecStoich,
                                 Rcpp::CharacterVector SpecName,
                                 Rcpp::NumericVector SpecActivityCoef);

Rcpp::NumericVector CalcStep(Rcpp::NumericMatrix JacobianMatrix,
                              Rcpp::NumericVector Resid,
                              unsigned int NComp,
                              Rcpp::CharacterVector CompType,
                              Rcpp::CharacterVector CompName);

void CompUpdate(unsigned int NComp, 
                Rcpp::NumericVector CompConcStep,
                Rcpp::NumericVector &CompConc);

Rcpp::NumericVector InitialGuess(Rcpp::NumericVector TotConc,
                                    Rcpp::CharacterVector CompType,
									                  Rcpp::NumericVector SpecK,
                                    Rcpp::IntegerMatrix SpecStoich,
                                    Rcpp::CharacterVector SpecName,
                                    unsigned int NComp,
                                    unsigned int NSpec);

Rcpp::NumericMatrix Jacobian (unsigned int NComp, //number of components
                              unsigned int NSpec, //number of species
                              Rcpp::IntegerMatrix SpecStoich, //formation reaction stoichiometry (NSpec x NComp)
                              Rcpp::NumericVector SpecConc, //species concentrations
                              Rcpp::NumericVector SpecCtoM, //concentration to mass conversion for each species
                              Rcpp::CharacterVector CompName, //names of components
                              unsigned int MetalComp, //position of the metal component
                              unsigned int NBLMetal, //number of BL-Metal species
                              Rcpp::IntegerVector BLMetalSpecs, //positions of BL-metal species
                              bool DoTox);

void UpdateTotals(unsigned int NComp, 
                  unsigned int NSpec, 
                  bool DoTox, 
                  Rcpp::CharacterVector CompType, 
                  Rcpp::CharacterVector CompName, 
                  Rcpp::String MetalName, 
                  Rcpp::IntegerMatrix SpecStoich, 
                  Rcpp::NumericVector SpecMoles, 
                  Rcpp::NumericVector CompCtoM, 
                  Rcpp::NumericVector &TotMoles, 
                  Rcpp::NumericVector &TotConc);

Rcpp::List UpdateTotalsList(unsigned int NComp, 
                            unsigned int NSpec, 
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
                                   unsigned int NSpec,
                                   Rcpp::NumericVector SpecK,
                                   Rcpp::NumericVector SpecTempKelvin,
                                   Rcpp::NumericVector SpecDeltaH);

double CalcIonicStrength(unsigned int NSpec,
                         Rcpp::NumericVector SpecMoles,
                         Rcpp::IntegerVector SpecCharge,
                         Rcpp::IntegerVector SpecMC);

double CalcChargeBalance(unsigned int NSpec,
                         Rcpp::NumericVector SpecMoles,
                         Rcpp::IntegerVector SpecCharge,
                         Rcpp::IntegerVector SpecMC);

Rcpp::NumericVector CalcActivityCoef(unsigned int NSpec,
                                     Rcpp::CharacterVector SpecActCorr,
                                     Rcpp::IntegerVector SpecCharge,
                                     double IonicStrength,
                                     double SysTempKelvin);

Rcpp::NumericVector CalcWHAMSpecCharge(unsigned int NSpec, 
                                       Rcpp::CharacterVector SpecActCorr,
                                       Rcpp::NumericVector SpecMoles,
                                       Rcpp::IntegerVector SpecCharge,
                                       Rcpp::IntegerVector SpecMC,
                                       int AqueousMC);

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
                                          Rcpp::NumericVector SolHS);

Rcpp::NumericVector CalcIonicStrengthEffects(double IonicStrength,
                                             Rcpp::NumericVector WHAMSpecCharge,
                                             unsigned int NSpec,
                                             Rcpp::IntegerVector SpecCharge,
                                             Rcpp::NumericVector SpecK,
                                             Rcpp::CharacterVector SpecActCorr,
                                             Rcpp::NumericVector wP);

const unsigned int iHA = 0; // humic acid = 0
const unsigned int iFA = 1; // fulvic acid = 1

#endif //__CHESSFUNCTIONS_H__