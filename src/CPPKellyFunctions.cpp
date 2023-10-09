#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> CppCalcSpecConc(std::vector<double> CConc,
                                    std::vector<double> K,
                                    Rcpp::IntegerMatrix Stoich,
                                    unsigned int NComp,
                                    unsigned int NSpec){
  /*
   * Purpose:
   * Inputs:
   *   CConc
   *   K
   *   Stoich
   *   NComp
   *   NSpec
   * Outputs:
   *   SConc(NSpec) vector
   */

  /* Variable definitions */
  std::vector<double> SConc(NSpec);//species concentrations
  double Tmp;//temporary variable

  for (unsigned int iSpec = 0; iSpec < NSpec; iSpec ++){
    Tmp = 1;
    for (unsigned int iComp = 0; iComp < NComp; iComp ++){
      Tmp *= std::pow(CConc[iComp], Stoich(iSpec, iComp));
    }
    SConc[iSpec] = Tmp * K[iSpec];
  }
  return(SConc);
}
