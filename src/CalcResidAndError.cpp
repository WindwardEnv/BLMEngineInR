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

#include <Rcpp.h>
#include <math.h>
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
                       Rcpp::CharacterVector SpecType,
                       Rcpp::NumericVector &Resid,
                       Rcpp::NumericVector &CompError) {

  /* variables: */
  int iComp; // loop counters
  /*Rcpp::NumericVector TotHSMoles(2);
  
  TotHSMoles(iHA) = 0.0;
  TotHSMoles(iFA) = 0.0;
  for (iComp = 0; iComp < NComp; iComp++) {
    if (SpecActCorr(iComp) == ACTYPE_WHAMHA) {
      TotHSMoles(iHA) += TotMoles(iComp);
    } else if (SpecActCorr(iComp) == ACTYPE_WHAMFA) {
      TotHSMoles(iFA) = TotMoles(iComp);
    }
  }*/

  // Calculate the residuals
  Resid = CalcTotMoles - TotMoles;
  
  // Calculate the error fraction for each component
  CompError = Rcpp::abs(Resid / TotMoles);
  
  for (iComp = 0; iComp < NComp; iComp++) {
    if ((CompType(iComp) == CTYPE_FIXEDCONC) || (CompType(iComp) == CTYPE_FIXEDACT)) {
      Resid(iComp) = 0.0;
      CompError(iComp) = 0.0;
    //} else if (SpecType(iComp) == STYPE_WHAMHA) {
    //  CompError(iComp) = Rcpp::abs(Resid(iComp) / TotHSMoles[iHA]);
    //} else if (SpecType(iComp) == STYPE_WHAMFA) {
    //  CompError(iComp) = Rcpp::abs(Resid(iComp) / TotHSMoles[iFA]);
    }
  }  
  
}



