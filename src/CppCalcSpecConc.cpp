#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

//' Calculate Species Concentrations
//'
//' Use `CppCalcSpecConc` to calculate species concentrations from a known set of
//' free component ion concentrations.
//'
//' This is an internal function that will, for each species `i` in `NSpec`
//' calculate the equilibrium concentration, which in its most basic form would
//' be calculated as \eqn{SConc_{i} = K_{i} *
//' \prod_{j=1}^{n}(CConc^SpecStoich_{i,j})}. Further modifications to these
//' calculations would include temperature correction, ionic strength
//' corrections, and diffuse double layer calculations for organic matter
//' binding. As of 19-Oct-2023, none of these corrections have been implemented.
//'
//' This is the C++ version of this function.
//'
//' @param CConc A vector of component concentrations for each of `NComp` components.
//' @param K A vector of reaction equilibrium constants for each of `NSpec` reactions.
//' @param SpecStoich A matrix of reaction stoichiometry, with `NSpec` rows and `NComp` columns.
//' @param NComp The number of components in the equilibrium system.
//' @param NSpec The number of species (reactions) in the equilibrium system.
//'
//' @returns A vector of `NSpec` species concentrations.
//'
//' @examples
//' # data(TestDataFreeConc, TestDataK, TestDataStoich)
//' # CppCalcSpecConc(CConc = TestDataFreeConc[1:2], K = TestDataK[3:4], SpecStoich = TestDataStoich[3:4,])
//'
//' @noRd
// [[Rcpp::export]]
std::vector<double> CppCalcSpecConc(std::vector<double> CConc,
                                   std::vector<double> K,
                                   Rcpp::IntegerMatrix SpecStoich,
                                   unsigned int NComp,
                                   unsigned int NSpec){
 /* Variable definitions */
 std::vector<double> SConc(NSpec);//species concentrations

 for (unsigned int iSpec = 0; iSpec < NSpec; iSpec ++){
   SConc[iSpec] = K[iSpec];
   for (unsigned int iComp = 0; iComp < NComp; iComp ++){
     SConc[iSpec] *= std::pow(CConc[iComp], SpecStoich(iSpec, iComp));
   }
 }
 return(SConc);
}

//' Calculate Species Concentrations
//'
//' Use `CppCalcLogSpecConc` to calculate species concentrations from a known set of
//' free component ion concentrations.
//'
//' This is an internal function that will, for each species `i` in `NSpec`
//' calculate the equilibrium concentration, which in its most basic form would
//' be calculated as \eqn{SConc_{i} = K_{i} *
//' \prod_{j=1}^{n}(CConc^SpecStoich_{i,j})}. Further modifications to these
//' calculations would include temperature correction, ionic strength
//' corrections, and diffuse double layer calculations for organic matter
//' binding. As of 19-Oct-2023, none of these corrections have been implemented.
//'
//' This is the C++ version of this function.
//'
//' @param LogCConc A vector of log10-transformed component concentrations for each of `NComp` components.
//' @param LogK A vector of log10-transformed reaction equilibrium constants for each of `NSpec` reactions.
//' @param SpecStoich A matrix of reaction stoichiometry, with `NSpec` rows and `NComp` columns.
//' @param NComp The number of components in the equilibrium system.
//' @param NSpec The number of species (reactions) in the equilibrium system.
//'
//' @returns A vector of `NSpec` log10-transformed species concentrations.
//'
//' @examples
//' # data(TestDataFreeConc, TestDataK, TestDataStoich)
//' # Log_TestDataFreeConc = log10(TestDataFreeConc)
//' # Log_TestDataK = log10(TestDataK)
//' # 10^CppCalcLogSpecConc(LogCConc = Log_TestDataFreeConc[1:2], LogK = Log_TestDataK[3:4], SpecStoich = TestDataStoich[3:4,])
//'
//' @noRd
// [[Rcpp::export]]
std::vector<double> CppCalcLogSpecConc(std::vector<double> LogCConc,
                                       std::vector<double> SpecLogK,
                                       Rcpp::IntegerMatrix SpecStoich,
                                       unsigned int NComp,
                                       unsigned int NSpec){
  /* Variable definitions */
  std::vector<double> LogSConc(NSpec);//species concentrations

  for (unsigned int iSpec = 0; iSpec < NSpec; iSpec ++){
    LogSConc[iSpec] = SpecLogK[iSpec];
    for (unsigned int iComp = 0; iComp < NComp; iComp ++){
      LogSConc[iSpec] += (LogCConc[iComp] * SpecStoich(iSpec, iComp));
    }
  }
  return(LogSConc);
}
