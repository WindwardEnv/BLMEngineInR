#include <Rcpp.h>
#include <cmath>
#include "CHESSFunctions.h"
//' Set the initial guess for component concentrations
//'
//' The initial guess is calculated by first assuming that component
//' concentrations are equal to the total concentrations.  Then species are
//' calculated based on these components, and the ratio of component totals to
//' actual totals to adjust the estimates of the free component concentrations.
//' Since there are a lot of inter-dependencies between the components, a few
//' iterations are needed to develop good intitial guesses.
//'
//' @param TotConc numeric vector (`NComp`); the total concentrations of components
//' @param CompType character vector (`NComp`); the component type
//' @param SpecK A vector of reaction equilibrium constants for each of `NSpec` reactions.
//' @param SpecStoich A matrix of reaction stoichiometry, with `NSpec` rows and `NComp` columns.
//' @param SpecName character vector (NSpec), the names of the chemical species
//' @param NComp The number of components in the equilibrium system.
//' @param NSpec The number of species (reactions) in the equilibrium system.#' @param NComp integer; the number of components
//'
//' @return CompConc The initial guesses for the free ion concentrations of components
//'
//' @keywords internal
//'
//' @noRd
//[[Rcpp::export]]
Rcpp::NumericVector InitialGuess(Rcpp::NumericVector TotConc,
                                    Rcpp::CharacterVector CompType,
									                  Rcpp::NumericVector SpecK,
                                    Rcpp::IntegerMatrix SpecStoich,
                                    Rcpp::CharacterVector SpecName,
                                    unsigned int NComp,
                                    unsigned int NSpec){
  /* outputs */
  Rcpp::NumericVector CompConc(NComp);//component concentrations

  /* Variables */
  unsigned int iSpec, iComp, iRound; //loop counters
  Rcpp::NumericVector SpecConc(NSpec);//species concentrations
  Rcpp::NumericVector CalcTot(NComp);//calculated total component concentrations

  /* Seed component concentrations with total concentrations */
  for (iComp = 0; iComp < NComp; iComp++){
	  CompConc(iComp) = TotConc(iComp);
  }

  /* Perform three rounds of adjustments */
  for (iRound =1; iRound <= 3; iRound++){
	  /* Calc Species */
	  for (iSpec = 0; iSpec < NSpec; iSpec ++){
		SpecConc(iSpec) = SpecK(iSpec);
		for (iComp = 0; iComp < NComp; iComp ++){
		   SpecConc(iSpec) = SpecConc(iSpec) * pow(CompConc(iComp), SpecStoich(iSpec, iComp));
		}
	  }
      /* Calc totals */
      for (iComp = 0; iComp < NComp; iComp++){
	    CalcTot(iComp) = 0;
	    for (iSpec = 0; iSpec < NSpec; iSpec++){
	      if (SpecStoich(iSpec, iComp) != 0){
             CalcTot(iComp) += SpecConc(iSpec) * SpecStoich(iSpec, iComp);
          }
        }
      }
	  /* Adjust component concentrations */
      for (iComp = 0; iComp < NComp; iComp++){
        if (CompType(iComp) == "MassBal"){
          CompConc(iComp) = CompConc(iComp)*(TotConc(iComp)/CalcTot(iComp));
        }
      }
  }
  return(CompConc);
}
