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
                         int AqueousMC) {
  /* output */
  double IonicStrength = 0;

  /* variables */
  int iSpec;

  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    if (SpecMC(iSpec) == AqueousMC) {
      IonicStrength += SpecCharge(iSpec) * SpecCharge(iSpec) * SpecMoles(iSpec);
    }
  }
  IonicStrength *= 0.5;

  return IonicStrength;
}