#include <Rcpp.h>
#include "ExternalCHESSFunctions.h"


//' @title Update Fixed Activity and Metal Totals
//'
//' @description Update the totals that change with each iteration - the totals
//'   for fixed activity components, and the metal component during tox mode.
//'   (returns a list)
//'
//' @param NComp integer, the combined number of components in the simulation,
//'   including the input components, defined components (and including the
//'   defined components that get added by ExpandWHAM)
//' @param NSpec integer, the number of chemical species for which we have
//'   formation reactions in the simulation
//' @param DoTox logical, TRUE for toxicity mode where the MetalName component
//'   concentration is adjusted to try to match the CATarget with BLMetalSpecs
//' @param CompType character vector (NComp), the type of each component in the
//'   simulation
//' @param CompName character vector (NComp), the name of each component in the
//'   simulation
//' @param MetalName The name of the component that corresponds to the metal
//'   associated with toxic effects.
//' @param SpecStoich signed integer matrix (NSpec x NComp), the reaction
//'   stoichiometry of the formation reactions.
//' @param SpecConc numeric vector (NSpec), the concentrations of each species 
//'   for which we have formation reactions
//' @param CompCtoM numeric vector (NSpec), the concentration to mass conversion
//'   factor of the components
//' @param TotMoles numeric vector (NComp), the total moles of each component in
//'   the simulation (units of mol)
//' @param TotConc numeric vector (NComp), the total concentrations of each
//'   component in the simulation (units of e.g., mol/L and mol/kg)
//'
//' @return Rcpp::List (Rcpp::NumericVector TotMoles,
//'   Rcpp::NumericVector TotConc)
//'
Rcpp::List UpdateTotalsList(int NComp,
                            int NSpec,
                            bool DoTox,
                            Rcpp::CharacterVector CompType,
                            Rcpp::CharacterVector CompName,
                            Rcpp::String MetalName,
                            Rcpp::IntegerMatrix SpecStoich,
                            Rcpp::NumericVector SpecConc,
                            Rcpp::NumericVector CompCtoM,
                            Rcpp::NumericVector TotMoles,
                            Rcpp::NumericVector TotConc) {

  int iComp, iSpec;

  for (iComp = 0; iComp < NComp; iComp++) {
    if ((CompType(iComp) == "FixedConc") || 
        (CompType(iComp) == "FixedAct") || 
        (DoTox & (CompName(iComp) == MetalName))) {
      TotConc(iComp) = 0;
      for (iSpec = 0; iSpec < NSpec; iSpec++) {
        TotConc(iComp) += SpecStoich(iSpec, iComp) * SpecConc(iSpec);
      }
      TotMoles(iComp) = TotConc(iComp) * CompCtoM(iComp);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("TotMoles") = TotMoles,
    Rcpp::Named("TotConc") = TotConc
  );
}

void UpdateFixedComps(int NComp, 
                      Rcpp::CharacterVector CompType, 
                      Rcpp::NumericVector InCompConc, 
                      Rcpp::NumericVector SpecActivityCoef,
                      Rcpp::NumericVector &SpecConc, 
                      Rcpp::NumericVector &CompConc) {
  /* variables */
  int iComp;

  for(iComp = 0; iComp < NComp; iComp++) {
    if (CompType[iComp] == "FixedAct") {
      SpecConc[iComp] = InCompConc[iComp] / SpecActivityCoef[iComp];
      CompConc[iComp] = SpecConc[iComp];
    } else if (CompType[iComp] == "FixedConc"){
      SpecConc[iComp] = InCompConc[iComp];
      CompConc[iComp] = SpecConc[iComp];
    }
  }
  
}