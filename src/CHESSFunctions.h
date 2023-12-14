#ifndef __CHESSFUNCTIONS_H__
#define __CHESSFUNCTIONS_H__

#include <string>
#include <Rcpp.h>

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
                    
Rcpp::NumericVector CalcSpecConc(Rcpp::NumericVector CompConc,
                                 Rcpp::NumericVector SpecK,
                                 Rcpp::IntegerMatrix SpecStoich,
                                 Rcpp::CharacterVector SpecName,
                                 unsigned int NComp,
                                 unsigned int NSpec);

Rcpp::NumericVector CalcLogSpecConc(Rcpp::NumericVector LogCompConc,
                                    Rcpp::NumericVector SpecLogK,
                                    Rcpp::IntegerMatrix SpecStoich,
                                    Rcpp::CharacterVector SpecName,
                                    unsigned int NComp,
                                    unsigned int NSpec);

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
                  Rcpp::CharacterVector CompType, 
                  Rcpp::CharacterVector CompName, 
                  Rcpp::String MetalName, 
                  Rcpp::NumericVector &TotMoles, 
                  Rcpp::IntegerMatrix SpecStoich, 
                  Rcpp::NumericVector SpecMoles, 
                  Rcpp::NumericVector &TotConc, 
                  Rcpp::NumericVector SpecCtoM, 
                  bool DoTox);

Rcpp::List UpdateTotalsList(unsigned int NComp, 
                            unsigned int NSpec, 
                            Rcpp::CharacterVector CompType, 
                            Rcpp::CharacterVector CompName, 
                            Rcpp::String MetalName, 
                            Rcpp::NumericVector TotMoles, 
                            Rcpp::IntegerMatrix SpecStoich, 
                            Rcpp::NumericVector SpecMoles, 
                            Rcpp::NumericVector TotConc, 
                            Rcpp::NumericVector SpecCtoM, 
                            bool DoTox);

Rcpp::NumericVector TempCorrection(double SysTemp,
                                   unsigned int NSpec,
                                   Rcpp::NumericVector SpecK,
                                   Rcpp::NumericVector SpecTemp,
                                   Rcpp::NumericVector SpecDeltaH);
#endif //__CHESSFUNCTIONS_H__