#include <Rcpp.h>
#include <cmath>
//#include "RcppArmaHelper.h"
#include "CHESSFunctions.h"

//' @title Calculate the Resid and CompError
//'
//' @description Calculate the residuals of the speciation problem, and the
//'   component Error.
//'
//' @details The residuals of the speciation problem are the calculated total
//'   moles for each component minus their known concentrations. This function
//'   does not account for the different toxic metal residual needed in a 
//'   toxicity run.
//'
//' @param NComp integer, the number of components
//' @param TotMoles numeric vector (NComp), the total moles of each component
//' @param CalcTotMoles numeric vector (NComp), the calculated total moles of 
//'   each component = sum(SpecConc * SpecCtoM * SpecStoich[,j]) for each 
//'   component j
//' @param CompType character vector (NComp), the type of component. It should
//'   be a fixed set of values (MassBal, FixedConc, Substituted, ChargeBal,
//'   SurfPot)
//' @param Resid (return value) numeric vector (NComp), the residuals = 
//'   calculated totals - known totals
//' @param CompError (return value) numeric vector (NComp), the absolute error
//'   fraction for each component in this iteration = abs(Resid / TotMoles)
//'
//' return void
//'
//' @name CalcResidualAndError
//' @usage CalcResidAndError(NComp, CalcTotMoles, TotMoles, CompType, Resid, CompError);
void CalcResidAndError(int NComp,
                       Rcpp::NumericVector CalcTotMoles,
                       Rcpp::NumericVector TotMoles,
                       Rcpp::CharacterVector CompType,
                       Rcpp::NumericVector &Resid,
                       Rcpp::NumericVector &CompError) {

  /* variables: */
  int iComp; // loop counters

  // Calculate the residuals
  Resid = CalcTotMoles - TotMoles;
  for (iComp = 0; iComp < NComp; iComp++) {
    if ((CompType(iComp) == "FixedConc") || (CompType(iComp) == "FixedAct")) {
      Resid(iComp) = 0.0;
    }
  }
  
  // Calculate the error fraction for each component
  CompError = abs(Resid / TotMoles);

}



