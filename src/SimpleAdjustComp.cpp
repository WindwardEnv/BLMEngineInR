#include <Rcpp.h>
#include "CHESSFunctions.h"

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
                      Rcpp::NumericVector WHAMSpecCharge) {
      
  /* variables */
  double CalcTotMolesi;
  double CompErrori;
  int Iter = 0;
  
  // calculate the total moles for just this component
  CalcTotMolesi = 
      CalcCompTotalMoles(iComp, NSpec, SpecConc * SpecCtoMAdj, SpecStoich);
  
  // calculate a starting CompError value for iComp
  CompErrori = abs(TotMolesi - CalcTotMolesi) / TotMolesi;

  while ((CompErrori > ConvCrit) && (Iter < MaxIter)) {
    Iter++;
    
    // adjust the component concentration
    CompConc(iComp) = CompConc(iComp) * TotMolesi / CalcTotMolesi;

    // make sure it's not unreasonable
    if (CompConc(iComp) <= 0.0) {
      CompConc(iComp) = 1.0E-20;
      break;
    }

    // calculate new species concentrations
    SpecConc = CalcSpecConc(NComp, NSpec, CompConc, SpecKISTempAdj, SpecStoich, 
                            SpecName, SpecActCorr, SpecActivityCoef, true, 
                            SpecCharge, WHAMSpecCharge);

    // calculate the total moles for just this component
    CalcTotMolesi = 
      CalcCompTotalMoles(iComp, NSpec, SpecConc * SpecCtoMAdj, SpecStoich);
    
    // update the CompError
    CompErrori = abs(TotMolesi - CalcTotMolesi) / TotMolesi;
            
  }

}