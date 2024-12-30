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

//' @title Calculate the Moles For This Iteration
//' 
//' @description Sum the relevant species moles to calculate the 
//'   resulting total moles for each component.
//'
//' @param NComp integer, the number of components
//' @param NSpec integer, the number of species
//' @param SpecStoich integer matrix (NSpec x NComp), the stoichiometry of
//'   species formatin reactions
//' 
//' @return CalcTotMoles (return parameter) numeric vector (NComp), the 
//'   calculated total moles of each component = sum(SpecConc * SpecCtoM * 
//'   SpecStoich[,j]) for each component j
//' 
//' @name CalcIterationTotalMoles
//' @usage CalcIterationTotalMoles(NComp, NSpec, SpecMoles, SpecStoich, 
//'   CalcTotMoles);
Rcpp::NumericVector CalcIterationTotalMoles(int NComp,
                                            int NSpec,
                                            Rcpp::NumericVector SpecMoles,
                                            Rcpp::IntegerMatrix SpecStoich) {

  /* output */
  Rcpp::NumericVector CalcTotMoles(NComp);

  /* variables: */
  int iComp, iSpec; // loop counters
  
  // Calculate the total moles and concentrations from species concentrations
  for (iComp = 0; iComp < NComp; iComp++){
    CalcTotMoles(iComp) = 0;
    for (iSpec = 0; iSpec < NSpec; iSpec++){
      if (SpecStoich(iSpec, iComp) != 0) {
        CalcTotMoles(iComp) += SpecMoles(iSpec) * SpecStoich(iSpec, iComp);
      }
    }
  }

  return CalcTotMoles;

}

double CalcCompTotalMoles(int iComp,
                          int NSpec,
                          Rcpp::NumericVector SpecMoles,
                          Rcpp::IntegerMatrix SpecStoich) {

  /* output */
  double CalcTotMolesi;

  /* variables: */
  int iSpec; // loop counters
  
  // Calculate the total moles and concentrations from species concentrations
  CalcTotMolesi = 0.0;
  for (iSpec = 0; iSpec < NSpec; iSpec++){
    if (SpecStoich(iSpec, iComp) != 0) {
      CalcTotMolesi += SpecMoles(iSpec) * SpecStoich(iSpec, iComp);
    }
  }

  return CalcTotMolesi;

}
