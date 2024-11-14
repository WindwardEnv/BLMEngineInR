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

Rcpp::NumericVector CalcIonicStrengthEffects(double WHAMIonicStrength,
                                             Rcpp::NumericVector WHAMSpecCharge,
                                             int NSpec,
                                             Rcpp::IntegerVector SpecCharge,
                                             Rcpp::NumericVector SpecK,
                                             Rcpp::CharacterVector SpecType,
                                             Rcpp::NumericVector wP) {
  /* output */
  Rcpp::NumericVector SpecKISAdj(NSpec);
    SpecKISAdj.names() = SpecK.names();

  /* variables */
  Rcpp::NumericVector W = wP * log10(WHAMIonicStrength);
  Rcpp::NumericVector WZ2 = -2 * W * WHAMSpecCharge;
  int iSpec;
  
  //adjust the intrinsic K's for WHAM species based on the charge
  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    SpecKISAdj(iSpec) = SpecK(iSpec);
    if (SpecCharge(iSpec) != 0) {
      if (SpecType(iSpec) == STYPE_WHAMHA) {
        SpecKISAdj(iSpec) = SpecK(iSpec) * exp(WZ2(iHA) * SpecCharge(iSpec));
      } else if (SpecType(iSpec) == STYPE_WHAMFA) {
        SpecKISAdj(iSpec) = SpecK(iSpec) * exp(WZ2(iFA) * SpecCharge(iSpec));
      }
    }    
  }

  return SpecKISAdj;

}
