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

//' @title Update fixed activity/concentration components
//' 
//' @description Components which are fixed activity or fixed concentration 
//'   should maintain an activity or concentration equal to their input 
//'   concentration. This function resets the SpecConc and CompConc vectors to
//'   those values.
//' 
//' @author Kelly Croteau (kellyc@windwardenv.com)
//' 
//' @param NComp integer, the combined number of components in the simulation,
//'   including the input components, defined components (and including the
//'   defined components that get added by ExpandWHAM)
//' @param CompType character vector (NComp), the type of each component in the
//'   simulation
//' @param InCompConc numeric vector (NComp) of the input value of each
//'   component's concentration, whether it's on the basis of total
//'   concentration, fixed activity, or fixed concentration
//' @param SpecActivityCoef {type}, {text}
//' @param SpecConc numeric vector (NSpec), the concentrations of each species 
//'   for which we have formation reactions
//' @param CompConc numeric vector (NComp) of component concentrations
//' 
//' 
void UpdateFixedComps(int NComp, 
                      Rcpp::CharacterVector CompType, 
                      Rcpp::NumericVector InCompConc, 
                      Rcpp::NumericVector SpecActivityCoef,
                      Rcpp::NumericVector &SpecConc, 
                      Rcpp::NumericVector &CompConc) {
  /* variables */
  int iComp;

  for(iComp = 0; iComp < NComp; iComp++) {
    if (CompType[iComp] == CTYPE_FIXEDACT) {
      SpecConc[iComp] = InCompConc[iComp] / SpecActivityCoef[iComp];
      CompConc[iComp] = SpecConc[iComp];
    } else if (CompType[iComp] == CTYPE_FIXEDCONC){
      SpecConc[iComp] = InCompConc[iComp];
      CompConc[iComp] = SpecConc[iComp];
    }
  }
  
}
