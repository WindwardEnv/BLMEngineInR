#ifndef __CHESSFUNCTIONS_H__
#define __CHESSFUNCTIONS_H__

#include <string>
#include <Rcpp.h>

Rcpp::List CalcResidual (unsigned int NComp,
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

Rcpp::NumericVector CompUpdate(unsigned int NComp,
                               Rcpp::NumericVector CompConcStep,
                               Rcpp::NumericVector CompConc,
                               Rcpp::CharacterVector CompName);

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
                  std::string MetalName, 
                  Rcpp::NumericVector &TotMoles, 
                  Rcpp::IntegerMatrix SpecStoich, 
                  Rcpp::NumericVector SpecMoles, 
                  Rcpp::NumericVector &TotConc, 
                  Rcpp::NumericVector SpecCtoM, 
                  bool DoTox);

#endif //__CHESSFUNCTIONS_H__