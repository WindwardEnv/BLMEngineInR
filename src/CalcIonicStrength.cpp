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
                         Rcpp::CharacterVector SpecActCorr,
                         bool ExcludeOrgMatter) {
  /* output */
  double IonicStrength = 0;

  /* variables */
  int iSpec;
  bool IncludeSpec;

  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    IncludeSpec = (SpecMC(iSpec) == AqueousMC);
    if (ExcludeOrgMatter) {
      IncludeSpec &= (SpecActCorr(iSpec) != ACTYPE_WHAMHA);
      IncludeSpec &= (SpecActCorr(iSpec) != ACTYPE_WHAMFA);
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