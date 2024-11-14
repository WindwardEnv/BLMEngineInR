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

#include <math.h>
#include <Rcpp.h>
#include "CHESSFunctions.h"

///' @title SimpleAdjustComp
//' 
//' @description Change the component concentration of the selected component to
//'   match the supplied known total. The component concentration is updated 
//'   based on a simple ratio of known total to calculated total moles.
//' 
//' @author Kelly Croteau (kellyc@windwardenv.com)
//' 
//' @param iComp int, the component to update the concentration of
//' @param ConvCrit double, what CompError (=abs(TotMolesi - 
//'   CalcTotMoles(iComp)) / TotMolesi) is considered converged?
//' @param MaxIter int, how many iterations of updates should be done before 
//'   giving up?
//' @param TotMolesi double, the known total moles of component iComp
//' @param NComp int, the number of components
//' @param CompConc NumericVector, (updated and returned) the vector of 
//'   component concentrations
//' @param NSpec int, the number of formation reaction species
//' @param SpecConc NumericVector, (updated and returned) the vector of species 
//'   concentration
//' @param SpecKISTempAdj NumericVector, the equilibrium coefficients for the 
//'   formation reactions, adjusted for temperature and ionic strength
//' @param SpecStoich IntegerMatrix, matrix of the stochiometric coefficients 
//'   for the formation reactions
//' @param SpecName CharacterVector, vector of species names
//' @param SpecType CharacterVector, vector of the STYPE_* for each species
//' @param SpecActivityCoef NumericVector, the activity coefficient of each 
//'   species
//' @param SpecCtoMAdj NumericVector, coefficient to convert from concentrations
//'   to moles for each species
//' @param SpecCharge IntegerVector, the charge of the ion formed in each 
//'   species formation reaction
//' @param WHAMSpecCharge NumericVector, vector of length 2 - the charge of each
//'   WHAM humic species 
//' @param DoWHAM bool, does this simulation include WHAM calculations?
//' 
void SimpleAdjustComp(int iComp,
                      double ConvCrit,
                      int MaxIter,
                      double TotMolesi,
                      int NComp,
                      Rcpp::NumericVector &CompConc,
                      int NSpec,
                      Rcpp::NumericVector &SpecConc,
                      Rcpp::NumericVector SpecKISTempAdj,
                      Rcpp::IntegerMatrix SpecStoich,
                      Rcpp::CharacterVector SpecName,
                      Rcpp::CharacterVector SpecType,
                      Rcpp::NumericVector SpecActivityCoef,
                      Rcpp::NumericVector SpecCtoMAdj,
                      Rcpp::IntegerVector SpecCharge,
                      Rcpp::NumericVector WHAMSpecCharge) {
      
  /* variables */
  double CalcTotMolesi;
  double CompErrori;
  int Iter = 0;
  
  // calculate the total moles for just this component
  CalcTotMolesi = 
      CalcCompTotalMoles(iComp, NSpec, SpecConc * SpecCtoMAdj, SpecStoich);
  
  // calculate a starting CompError value for iComp
  CompErrori = std::fabs(TotMolesi - CalcTotMolesi) / TotMolesi;

  while ((CompErrori > ConvCrit) && (Iter < MaxIter)) {
    Iter++;
    
    // adjust the component concentration
    CompConc(iComp) = CompConc(iComp) * TotMolesi / CalcTotMolesi;

    // make sure it's not unreasonable
    if (CompConc(iComp) <= 0.0) {
      CompConc(iComp) = 1.0E-20;
      break;
    }

    // calculate new species concentrations
    SpecConc = CalcSpecConc(NComp, NSpec, CompConc, SpecKISTempAdj, SpecStoich, 
                            SpecName, SpecType, SpecActivityCoef, true, 
                            SpecCharge, WHAMSpecCharge);

    // calculate the total moles for just this component
    CalcTotMolesi = 
      CalcCompTotalMoles(iComp, NSpec, SpecConc * SpecCtoMAdj, SpecStoich);
    
    // update the CompError
    CompErrori = std::fabs(TotMolesi - CalcTotMolesi) / TotMolesi;
            
  }

}
