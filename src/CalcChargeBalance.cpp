#include <Rcpp.h>
#include "CHESSFunctions.h"

//' @title Calculate the Charge Balance
//'
//' @description Calculate the charge balance of the solution given the species
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
//' @return double, the net charge balance in solution
//'
double CalcChargeBalance(unsigned int NSpec,
                         Rcpp::NumericVector SpecMoles,
                         Rcpp::IntegerVector SpecCharge,
                         Rcpp::IntegerVector SpecMC) {
  /* output */
  double ChargeBal;
  
  /* variables */
  unsigned int iSpec;
  
  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    if (SpecMC(iSpec) == 1L) {
      ChargeBal += SpecCharge(iSpec) * SpecMoles(iSpec);
    }
  }
  
  return ChargeBal;

}


Rcpp::NumericVector CalcWHAMSpecCharge(unsigned int NSpec, 
                                       Rcpp::CharacterVector SpecActCorr,
                                       Rcpp::NumericVector SpecMoles,
                                       Rcpp::IntegerVector SpecCharge,
                                       Rcpp::IntegerVector SpecMC,
                                       int AqueousMC) {
  
  /* output */
  Rcpp::NumericVector WHAMSpecCharge(2);
  
  /* variables */
  unsigned int iSpec;
  unsigned int iHA = 0;
  unsigned int iFA = 1;

  WHAMSpecCharge(iHA) = 0;
  WHAMSpecCharge(iFA) = 0;
  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    if (SpecMC(iSpec) == AqueousMC) {
      if (SpecActCorr(iSpec) == "WHAMHA") {
        WHAMSpecCharge(iHA) += SpecMoles(iSpec) * SpecCharge(iSpec);
      } else if (SpecActCorr(iSpec) == "WHAMFA") {
        WHAMSpecCharge(iFA) += SpecMoles(iSpec) * SpecCharge(iSpec);
      }
    }
  }

  return WHAMSpecCharge;

}