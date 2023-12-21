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
double CalcIonicStrength(unsigned int NSpec,
                         Rcpp::NumericVector SpecMoles,
                         Rcpp::IntegerVector SpecCharge,
                         Rcpp::IntegerVector SpecMC) {
  /* output */
  double IonicStrength;

  /* variables */
  unsigned int iSpec;

  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    if (SpecMC(iSpec) == 1L) {
      IonicStrength += SpecCharge(iSpec) * SpecCharge(iSpec) * SpecMoles(iSpec);
    }
  }
  IonicStrength *= 0.5;

  return IonicStrength;
}

Rcpp::NumericVector ExtDebyeHuckel(double IonicStrength, 
                                   int MaxCharge, 
                                   double SysTempKelvin) {

  /* outputs */
  Rcpp::NumericVector ActivityCoefDebye(MaxCharge + 1);

  /* variables */
  int ZZ;
  double A, B, C, D;
  int E, F;

  A = 0.27 + (0.0008 * SysTempKelvin);
  B = 0.33;
  C = -A * sqrt(IonicStrength);
  D = B * sqrt(IonicStrength);

  for (ZZ = 0; ZZ <= MaxCharge; ZZ++) {
    E = abs(ZZ);
    F = ZZ * ZZ;
    ActivityCoefDebye(ZZ) = pow(10, F * C / (1 + (3 * E * D)));
  }
  
/* This is what's done in the PB code...Is this even necessary?? Doesn't the 
   math work out to these same equations if you just use the last one? 
  switch(E) {
    case 0: ActivityCoef = 1;
    case 1: ActivityCoef = pow(10, C / (1 + (3 * D)));
    case 2: ActivityCoef = pow(10, 4 * C / (1 + (6 * D)));
    case 3: ActivityCoef = pow(10, 9 * C / (1 + (9 * D)));
    case 4: ActivityCoef = pow(10, 16 * C / (1 + (12 * D)));
    default : ActivityCoef = pow(10, F * C / (1 + (3 * E * D)));
  }
*/
  return ActivityCoefDebye;
}

Rcpp::NumericVector CalcActivityCoef(unsigned int NSpec,
                                     Rcpp::CharacterVector SpecActCorr,
                                     Rcpp::IntegerVector SpecCharge,
                                     double IonicStrength,
                                     double SysTempKelvin) {

  /* outputs */
  Rcpp::NumericVector ActivityCoef(NSpec);

  /* variables */
  unsigned int iSpec;
  int MaxCharge = Rcpp::max(Rcpp::abs(SpecCharge));
  Rcpp::NumericVector ActivityCoefDebye;
  
  ActivityCoefDebye = ExtDebyeHuckel(IonicStrength, MaxCharge, SysTempKelvin);

  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    if (SpecActCorr(iSpec) == "None") {
      ActivityCoef(iSpec) = 1.0;
    //} else if (SpecActCorr(iSpec) == "Davies") {
    //  ActivityCoef(iSpec) = Davies(IonicStrength, SpecCharge(iSpec));
    } else if (SpecActCorr(iSpec) == "Debye") {
      ActivityCoef(iSpec) = ActivityCoefDebye(abs(SpecCharge(iSpec)));
    }
  }

  return ActivityCoef;

}
