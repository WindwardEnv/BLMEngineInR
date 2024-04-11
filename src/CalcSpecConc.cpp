#include <cmath>
#include <Rcpp.h>
#include "CHESSFunctions.h"

//' Calculate Species Concentrations
//'
//' Use `CalcSpecConc` to calculate species concentrations from a known set of
//' free component ion concentrations.
//'
//' This is an internal function that will, for each species `i` in `NSpec`
//' calculate the equilibrium concentration, which in its most basic form would
//' be calculated as \eqn{SpecConc_{i} = K_{i} *
//' \prod_{j=1}^{n}(CompConc^SpecStoich_{i,j})}. Further modifications to these
//' calculations would include temperature correction, ionic strength
//' corrections, and diffuse double layer calculations for organic matter
//' binding. As of 19-Oct-2023, none of these corrections have been implemented.
//'
//' This is the C++ version of this function.
//'
//' @param CompConc A vector of component concentrations for each of `NComp` components.
//' @param SpecK A vector of reaction equilibrium constants for each of `NSpec` reactions.
//' @param SpecStoich A matrix of reaction stoichiometry, with `NSpec` rows and `NComp` columns.
//' @param SpecName character vector (NSpec), the names of the chemical species
//' @param NComp The number of components in the equilibrium system.
//' @param NSpec The number of species (reactions) in the equilibrium system.
//'
//' @returns A vector of `NSpec` species concentrations.
//'
//' @noRd
Rcpp::NumericVector CalcSpecConc(int NComp,
                                 int NSpec,
                                 Rcpp::NumericVector CompConc,
                                 Rcpp::NumericVector SpecK,
                                 Rcpp::IntegerMatrix SpecStoich,
                                 Rcpp::CharacterVector SpecName,
                                 Rcpp::CharacterVector SpecActCorr,
                                 Rcpp::NumericVector SpecActivityCoef,
                                 bool DoWHAM,
                                 Rcpp::IntegerVector SpecCharge,
                                 Rcpp::NumericVector WHAMSpecCharge) {
  /* outputs */
  Rcpp::NumericVector SpecConc(NSpec);//species concentrations
    SpecConc.names() = SpecName;

  /* Variables */
  int iSpec, iComp; //loop counters
  Rcpp::NumericVector SpecActivity(NSpec);
  Rcpp::NumericVector CompActivity(NComp);

  for (iComp = 0; iComp < NComp; iComp++) {
    CompActivity(iComp) = CompConc(iComp) * SpecActivityCoef(iComp);
    SpecActivity(iComp) = CompActivity(iComp);
  }

  for (iSpec = NComp; iSpec < NSpec; iSpec++) {
    SpecActivity(iSpec) = SpecK(iSpec);
    for (iComp = 0; iComp < NComp; iComp++) {
      if (SpecStoich(iSpec, iComp) != 0) {
          /*if ((SpecActCorr(iSpec) == "DonnanHA") || 
            (SpecActCorr(iSpec) == "DonnanFA")) {
              Rcpp::Rcout << SpecActCorr(iSpec) << std::endl;
          }*/
          /*if ((SpecActCorr(iComp) == "DonnanHA") || 
              (SpecActCorr(iComp) == "DonnanFA")) {
            SpecActivity(iSpec) *= std::pow(CompConc(iComp), 
                                            SpecStoich(iSpec, iComp));
          } else {
            SpecActivity(iSpec) *= std::pow(CompActivity(iComp), 
                                            SpecStoich(iSpec, iComp));
          }
          
        } else {*/
          SpecActivity(iSpec) *= std::pow(CompActivity(iComp),
                                        SpecStoich(iSpec, iComp));
        //}
      }
    }
  }

  if (DoWHAM) {
    for (iSpec = NComp; iSpec < NSpec; iSpec++) {
      if ((SpecActCorr(iSpec) == "DonnanHA") && 
          (((WHAMSpecCharge(iHA) < 0) && (SpecCharge(iSpec) < 0)) || 
          ((WHAMSpecCharge(iHA) > 0) && (SpecCharge(iSpec) > 0)))) {
        SpecActivity[iSpec] = 0.0;
      } else if ((SpecActCorr(iSpec) == "DonnanFA") &&
                (((WHAMSpecCharge(iFA) < 0) && (SpecCharge(iSpec) < 0)) || 
                  ((WHAMSpecCharge(iFA) > 0) && (SpecCharge(iSpec) > 0)))) {
        SpecActivity[iSpec] = 0.0;
      }
    }
  }
    

  SpecConc = SpecActivity / SpecActivityCoef;
  return(SpecConc);
}
