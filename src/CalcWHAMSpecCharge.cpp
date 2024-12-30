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

Rcpp::NumericVector CalcWHAMSpecCharge(int NSpec, 
                                       Rcpp::CharacterVector SpecType,
                                       Rcpp::NumericVector SpecConc,
                                       Rcpp::IntegerVector SpecCharge,
                                       Rcpp::IntegerVector SpecMC,
                                       int AqueousMC,
                                       Rcpp::NumericVector HumicSubstGramsPerLiter) {
  
  /* output */
  Rcpp::NumericVector WHAMSpecCharge(2);
  
  /* variables */
  int iSpec;
  
  WHAMSpecCharge(iHA) = 0;
  WHAMSpecCharge(iFA) = 0;
  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    if (SpecType(iSpec) == STYPE_WHAMHA) {
      WHAMSpecCharge(iHA) += SpecConc(iSpec) / HumicSubstGramsPerLiter[iHA] * SpecCharge(iSpec);
    } else if (SpecType(iSpec) == STYPE_WHAMFA) {
      WHAMSpecCharge(iFA) += SpecConc(iSpec) / HumicSubstGramsPerLiter[iFA] * SpecCharge(iSpec);
    }
  }

  //WHAMSpecCharge = WHAMSpecCharge / HumicSubstGramsPerLiter;
  //charge / g HS = (charge / L) * (L / g HS)

  return WHAMSpecCharge;

}
