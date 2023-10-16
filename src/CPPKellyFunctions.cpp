#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

//' Calculate species concentrations
//'
//' @param CConc vector of
//' @param K
// [[Rcpp::export]]
std::vector<double> CppCalcSpecConc(std::vector<double> CConc,
                                   std::vector<double> K,
                                   Rcpp::IntegerMatrix Stoich,
                                   unsigned int NComp,
                                   unsigned int NSpec){
 /* Variable definitions */
 std::vector<double> SConc(NSpec);//species concentrations

 for (unsigned int iSpec = 0; iSpec < NSpec; iSpec ++){
   SConc[iSpec] = K[iSpec];
   for (unsigned int iComp = 0; iComp < NComp; iComp ++){
     SConc[iSpec] *= std::pow(CConc[iComp], Stoich(iSpec, iComp));
   }
 }
 return(SConc);
}

//' Calculate species concentrations
//'
//' @param LogCConc vector of log10-transformed component concentrations
//' @param logK vector of log10-transformed equil constants
// [[Rcpp::export]]
std::vector<double> CppCalcLogSpecConc(std::vector<double> LogCConc,
                                       std::vector<double> LogK,
                                       Rcpp::IntegerMatrix Stoich,
                                       unsigned int NComp,
                                       unsigned int NSpec){
  /* Variable definitions */
  std::vector<double> LogSConc(NSpec);//species concentrations

  for (unsigned int iSpec = 0; iSpec < NSpec; iSpec ++){
    LogSConc[iSpec] = LogK[iSpec];
    for (unsigned int iComp = 0; iComp < NComp; iComp ++){
      LogSConc[iSpec] += (LogCConc[iComp] * Stoich(iSpec, iComp));
    }
  }
  return(LogSConc);
}
