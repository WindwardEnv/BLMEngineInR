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

//' @title Iterative step improvement in component concentrations
//'
//' @description CompUpdate calculates an iterative improvement on the component
//' concentrations based on the Newton-Raphson solution from the current
//' iteration.
//'
//' @details If the iteration would cause the adjusted component concentrations
//' to be less than zero, then the component concentration is simply divided by
//' 10 for this iteration.
//'
//' @param NComp integer, the number of components
//' @param CompConcStep numeric vector (NComp) of adjustments to the component
//'   concentrations
//' @param CompType character vector (NComp) with the types of each component
//' @param CompConc (INPUT & OUTPUT) numeric vector (NComp) of component
//'   concentrations, input values are from this iteration
//'
void CompUpdate(int NComp, 
                Rcpp::NumericVector CompConcStep,
                Rcpp::CharacterVector CompType,
                Rcpp::NumericVector &CompConc){

  /* variables */
  Rcpp::NumericVector oldCompConc(NComp);
  int iComp;

  for (iComp = 0; iComp < NComp; iComp++) {
    oldCompConc[iComp] = CompConc[iComp];
    if (CompConcStep[iComp] >= oldCompConc[iComp]) {
      CompConc[iComp] = oldCompConc[iComp] / 10;
    } else {
      CompConc[iComp] = oldCompConc[iComp] - CompConcStep[iComp];
    }
    if (((CompType[iComp] == CTYPE_DONNANHA) || (CompType[iComp] == CTYPE_DONNANFA))) {
      if (CompConc[iComp] < 1.0) { CompConc[iComp] = 1.0; }
    }
  }//NEXT iComp

}
