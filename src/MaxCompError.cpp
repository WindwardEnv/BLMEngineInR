#include <Rcpp.h>
#include "CHESSFunctions.h"

//' @title Find the MaxError and WhichMax
//'
//' @description Determine which component has the maximum absolute error 
//'   fraction.
//'
//' @param NComp integer, the number of components
//' @param CompError numeric vector (NComp), the absolute error fraction for
//'   each component in this iteration =abs(Resid / TotMoles)
//' @param MaxError (return value) numeric, the highest absolute error fraction
//'   in this iteration =max(abs(Resid / TotMoles))
//' @param WhichMax (return value) integer, the position in the component
//'   vectors of the component with the highest absolute error
//' 
//' @return void
//'
//' @name MaxCompError
//' @usage MaxError = MaxCompError(NComp, CompError, WhichMax)
double MaxCompError(int NComp, Rcpp::NumericVector CompError, 
                    int &WhichMax) {

  /* output */
  double MaxError;

  /* variable */
  int iComp;

  // Determine which component has the highest error fraction
  MaxError = CompError(0);
  WhichMax = 0;
  for (iComp = 1; iComp < NComp; iComp++){
    if (CompError(iComp) > MaxError){
      MaxError = CompError(iComp);
      WhichMax = iComp;
    }
  }

  return MaxError;
}
