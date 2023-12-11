#include <Rcpp.h>
#include "CHESSFunctions.h"

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
                  bool DoTox) {
    
    unsigned int iComp, iSpec;

            
    for (iComp = 0; iComp < NComp; iComp++) {
        if ((CompType(iComp) == "FixedAct") | (DoTox & (CompName(iComp) == MetalName))) {
            TotMoles(iComp) = 0;
            for (iSpec = 0; iSpec < NSpec; iSpec++) {
                TotMoles(iComp) += SpecStoich(iSpec, iComp) * SpecMoles(iSpec);
            }
            TotConc(iComp) = TotMoles(iComp) * SpecCtoM(iComp);
        }
    }
}

// [[Rcpp::export]]
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
                            bool DoTox) {
    
    unsigned int iComp, iSpec;

    for (iComp = 0; iComp < NComp; iComp++) {
        if ((CompType(iComp) == "FixedAct") || (DoTox & (CompName(iComp) == MetalName))) {
            TotMoles(iComp) = 0;
            for (iSpec = 0; iSpec < NSpec; iSpec++) {
                TotMoles(iComp) += SpecStoich(iSpec, iComp) * SpecMoles(iSpec);
            }
            TotConc(iComp) = TotMoles(iComp) * SpecCtoM(iComp);
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("TotMoles") = TotMoles,
        Rcpp::Named("TotConc") = TotConc
    );
}