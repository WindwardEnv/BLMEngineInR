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
#include "CHESSFunctions.h"

//' @title Calculate the Moles and Concentrations For This Iteration
//' 
//' @description Sum the relevant species concentrations to calculate the 
//'   resulting total concentrations and moles for each component.
//'
//' @param NComp integer, the number of components
//' @param NSpec integer, the number of species
//' @param SpecConc numeric vector (NSpec), the concentration of chemical
//'   species
//' @param SpecStoich integer matrix (NSpec x NComp), the stoichiometry of
//'   species formatin reactions
//' @param SpecCtoM numeric vector (NSpec), the concentration-to-mass conversion
//'   factor for each chemical species
//' @param CalcTotConc (return parameter) numeric vector (NComp), the 
//'   calculated total concentration of each component = CalcTotMoles / CtoM
//' @param CalcTotMoles (return parameter) numeric vector (NComp), the 
//'   calculated total moles of each component = sum(SpecConc * SpecCtoM * 
//'   SpecStoich[,j]) for each component j
//' 
//' @return void
//' 
//' @name CalcIterationTotals
//' @usage CalcIterationTotals(NComp, NSpec, SpecConc, SpecCtoM, SpecStoich, 
//'   CalcTotMoles, CalcTotConc);
void CalcIterationTotals(int NComp,
                         int NSpec,
                         Rcpp::NumericVector SpecConc,
                         Rcpp::NumericVector SpecCtoM,
                         Rcpp::IntegerMatrix SpecStoich,
                         Rcpp::NumericVector &CalcTotMoles,
                         Rcpp::NumericVector &CalcTotConc) {

  /* variables: */
  int iComp, iSpec; // loop counters
  Rcpp::NumericVector SpecMoles; // moles of each species
  
  // Calculate the total moles and concentrations from species concentrations
  SpecMoles = SpecConc * SpecCtoM;
  for (iComp = 0; iComp < NComp; iComp++){
    CalcTotMoles(iComp) = 0;
    for (iSpec = 0; iSpec < NSpec; iSpec++){
      if (SpecStoich(iSpec, iComp) != 0){
        CalcTotMoles(iComp) += (SpecMoles(iSpec) * SpecStoich(iSpec, iComp));
      }
    }
    CalcTotConc(iComp) = CalcTotMoles(iComp) / SpecCtoM(iComp); 
  }

}

