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

//' @title Find the MaxError and WhichMax
//'
//' @description Determine which component has the maximum absolute error 
//'   fraction.
//'
//' @param NComp integer, the number of components
//' @param CompError numeric vector (NComp), the absolute error fraction for
//'   each component in this iteration =abs(Resid / TotMoles)
//' @param MaxError (return value) numeric, the highest absolute error fraction
//'   in this iteration =max(abs(Resid / TotMoles))
//' @param WhichMax (return value) integer, the position in the component
//'   vectors of the component with the highest absolute error
//' 
//' @return void
//'
//' @name MaxCompError
//' @usage MaxError = MaxCompError(NComp, CompError, WhichMax)
double MaxCompError(int NComp, Rcpp::NumericVector CompError, 
                    int &WhichMax) {

  /* output */
  double MaxError;

  /* variable */
  int iComp;

  // Determine which component has the highest error fraction
  MaxError = CompError(0);
  WhichMax = 0;
  for (iComp = 1; iComp < NComp; iComp++){
    if (CompError(iComp) > MaxError){
      MaxError = CompError(iComp);
      WhichMax = iComp;
    }
  }

  return MaxError;
}
