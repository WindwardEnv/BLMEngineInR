#include <Rcpp.h>
#include "CHESSFunctions.h"

Rcpp::NumericVector TempCorrection(double SysTemp,
                                   unsigned int NSpec,
                                   Rcpp::NumericVector SpecK,
                                   Rcpp::NumericVector SpecTemp,
                                   Rcpp::NumericVector SpecDeltaH) {
  /*Outputs:*/
	Rcpp::NumericVector SpecKTempAdj(NSpec);//the enthalpy corrected log Ks for species
	
  /*Local variables:*/
	unsigned int iSpec;
	double Rcon = 8.314;//universal gas constant (Rcon = 8.314)	
  double T0, T1, T2;//temparary temperature variables 
  
  /*
	*  Check to see if temperature correction can be performed
	*/
	for (iSpec = 0; iSpec < NSpec; iSpec++) {
	  if ((SpecTemp(iSpec) != 0) & (SpecDeltaH(iSpec) != 0)) {
		  T0 = 1 / SysTemp;
			T1 = 1 / SpecTemp(iSpec);
			T2 = SpecDeltaH(iSpec) * (T0 - T1) / Rcon;
			SpecKTempAdj(iSpec) = SpecK(iSpec) * std::exp(T2); 
    } else {
        SpecKTempAdj(iSpec) = SpecK(iSpec);
    }
  }
  return SpecKTempAdj;
}
