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
  std::vector<double> SConc(NSpec);
  double Tmp;
  for (unsigned int iSpec = 0; iSpec < NSpec; iSpec ++){
    Tmp = 1;
    for (unsigned int iComp = 0; iComp < NComp; iComp ++){
      Tmp *= std::pow(CConc[iComp], Stoich(iSpec, iComp));
    }
    SConc[iSpec] = Tmp * K[iSpec];
  }
  return(SConc);
}
