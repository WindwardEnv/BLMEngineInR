#include <Rcpp.h>
#include "CHESSFunctions.h"

//' Correct the Equilibrium Coefficients for Temperature
//' 
//' The units for `SysTempKelvin` and `SpecTempKelvin` should be Kelvin.
//' 
//' @param SysTempKelvin double, the system temperature, in Kelvin
//' @param NSpec integer, the number of chemical species for which we have
//'   formation reactions in the simulation
//' @param SpecK numeric vector (NSpec), the equilibrium coefficient of the 
//'   formation reactions
//' @param SpecDeltaH numeric vector (NSpec), the enthalpy change of the
//'   formation reactions (J/mol)
//' @param SpecTempKelvin numeric vector (NSpec), the temperature associated with 
//'   K/logK and DeltaH of the formation reactions, in Kelvin
//' 
//' @return SpecKTempAdj numeric vector (NSpec), the the enthalpy corrected 
//'   log Ks for species
//' 
Rcpp::NumericVector TempCorrection(double SysTempKelvin,
                                   int NSpec,
                                   Rcpp::NumericVector SpecK,
                                   Rcpp::NumericVector SpecTempKelvin,
                                   Rcpp::NumericVector SpecDeltaH) {
  /* Outputs: */
	Rcpp::NumericVector SpecKTempAdj(NSpec);
	
  /* Local variables: */
	int iSpec;
	double Rcon = 8.314;//universal gas constant (Rcon = 8.314 J / (mol*K))
  double T0, T1, T2;//temparary temperature variables 
  
  /*
	*  Check to see if temperature correction can be performed
	*/
	T1 = 1 / SysTempKelvin;
	for (iSpec = 0; iSpec < NSpec; iSpec++) {
	  if ((SpecTempKelvin(iSpec) != 0) & (SpecDeltaH(iSpec) != 0)) {
		  T0 = 1 / SpecTempKelvin(iSpec);
			T2 = SpecDeltaH(iSpec) * (T0 - T1) / Rcon;
			SpecKTempAdj(iSpec) = SpecK(iSpec) * std::exp(T2); 
    } else {
        SpecKTempAdj(iSpec) = SpecK(iSpec);
    }
  }

  return SpecKTempAdj;

}
