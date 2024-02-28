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
double CalcChargeBalance(int NSpec,
                         Rcpp::NumericVector SpecMoles,
                         Rcpp::IntegerVector SpecCharge,
                         Rcpp::IntegerVector SpecMC,
                         int AqueousMC) {
  /* output */
  double ChargeBal;
  
  /* variables */
  int iSpec;
  
  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    if (SpecMC(iSpec) == AqueousMC) {
      ChargeBal += SpecCharge(iSpec) * SpecMoles(iSpec);
    }
  }
  
  return ChargeBal;

}

