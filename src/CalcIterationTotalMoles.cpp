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
      if (SpecStoich(iSpec, iComp) != 0){
        CalcTotMoles(iComp) += SpecMoles(iSpec) * SpecStoich(iSpec, iComp);
      }
    }
  }

  return CalcTotMoles;

}
