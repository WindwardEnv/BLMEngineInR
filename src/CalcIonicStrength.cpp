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

//' @title Calculate the Ionic Strength
//'
//' @description Calculate the ionic strength of the solution given the species
//'   concentrations and charges.
//'
//' @author Kelly Croteau (kellyc@windwardenv.com)
//'
//' @param NSpec integer, the number of chemical species for which we have
//'   formation reactions in the simulation
//' @param SpecMoles numeric vector (NSpec), the concentrations of each species
//'   for which we have formation reactions
//' @param SpecCharge signed integer vector (NSpec), the charge of the chemical
//'   species for which we have formation reactions
//' @param SpecMC integer vector (NSpec), the mass compartment of the chemical
//'   species for which we have formation reactions
//'
//' @return double, the ionic strength
//'
double CalcIonicStrength(int NSpec,
                         Rcpp::NumericVector SpecMoles,
                         Rcpp::IntegerVector SpecCharge,
                         Rcpp::IntegerVector SpecMC,
                         int AqueousMC, 
                         Rcpp::CharacterVector SpecType,
                         bool ExcludeOrgMatter) {
  /* output */
  double IonicStrength = 0;

  /* variables */
  int iSpec;
  bool IncludeSpec;

  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    IncludeSpec = (SpecMC(iSpec) == AqueousMC);      
    //IncludeSpec = IncludeSpec || (SpecType(iSpec) == STYPE_DONNANHA);
    //IncludeSpec = IncludeSpec || (SpecType(iSpec) == STYPE_DONNANFA);
    if (ExcludeOrgMatter) {
      IncludeSpec &= (SpecType(iSpec) != STYPE_WHAMHA);
      IncludeSpec &= (SpecType(iSpec) != STYPE_WHAMFA);    
    }
    if (IncludeSpec){ 
      IonicStrength += pow(SpecCharge(iSpec), 2) * SpecMoles(iSpec);
    } 
  }
  IonicStrength *= 0.5;  

  if (IonicStrength > 100.0) { IonicStrength = 100.0; }
  if (IonicStrength <= 0.0) { IonicStrength = 0.1; }

  return IonicStrength;
}
