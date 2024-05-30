j.#ifndef __EXTERNALCHESSFUNCTIONS_H__
#define __EXTERNALCHESSFUNCTIONS_H__

  #include <Rcpp.h>
#endif

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

#endif //__EXTERNALCHESSFUNCTIONS_H__