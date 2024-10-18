#include <Rcpp.h>
#include <math.h>
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
//' @param TotConc numeric vector (`NComp`); the total concentrations of
//'   components
//' @param SpecCtoM, numericVector (NSpec); the concentration to mass conversion
//'   factor for each species
//' @param CompType character vector (`NComp`); the component type
//' @param SpecK A vector of reaction equilibrium constants for each of `NSpec`
//'   reactions.
//' @param SpecStoich A matrix of reaction stoichiometry, with `NSpec` rows and
//'   `NComp` columns.
//' @param SpecName character vector (NSpec), the names of the chemical species
//' @param NComp The number of components in the equilibrium system.
//' @param NSpec The number of species (reactions) in the equilibrium system.
//'
//' @return CompConc The initial guesses for the free ion concentrations of
//'   components
//'
//' @keywords internal
//'
//' @noRd
//'
Rcpp::NumericVector InitialGuess(Rcpp::NumericVector &TotConc,
                                 Rcpp::NumericVector SpecCtoM,
                                 Rcpp::CharacterVector CompType,
									               Rcpp::NumericVector SpecK,
                                 Rcpp::IntegerMatrix SpecStoich,
                                 Rcpp::CharacterVector SpecName,
                                 int NComp,
                                 int NSpec, 
                                 bool DoTox,
                                 int NBLMetal, 
                                 Rcpp::IntegerVector BLMetalSpecs, 
                                 int MetalComp,
                                 double CATarget) {
  
  /* outputs */
  Rcpp::NumericVector CompConc(NComp);//component concentrations

  /* Variables */
  int iSpec, iComp, iRound; //loop counters
  Rcpp::NumericVector SpecConc(NSpec);//species concentrations
  //Rcpp::NumericVector SpecMoles(NSpec);
  Rcpp::NumericVector CalcTotConc(NComp);//calculated total component concentrations
  Rcpp::NumericVector CalcTotMoles(NComp);//calculated total component moles 
  Rcpp::NumericVector TotMoles(NComp);
  double CACalculated;
  int i;
  /*int CompTotStoich;
  double CompStoichMoles;
  Rcpp::NumericVector MaxSpecConc(NSpec);

  for (iSpec = 0; iSpec < NSpec; iSpec++) {    
    MaxSpecConc(iSpec) = INFINITY;
    for (iComp = 0; iComp < NComp; iComp++) {
      if ((((CompType(iComp) == CTYPE_MASSBAL) || 
            (CompType(iComp) == CTYPE_WHAMHA) ||
            (CompType(iComp) == CTYPE_WHAMFA))) && 
          (TotConc(iComp) < MaxSpecConc(iSpec))) {
        MaxSpecConc(iSpec) = TotConc(iComp);
      }
    }
  }*/

  /* Seed component concentrations with total concentrations */
  for (iComp = 0; iComp < NComp; iComp++) {    
    /*if ((CompType[iComp] == CTYPE_FIXEDACT) ||
        (CompType[iComp] == CTYPE_FIXEDCONC)) {
      CompConc(iComp) = TotConc(iComp);
    } else if ((CompType(iComp) == CTYPE_MASSBAL) || 
                 (CompType(iComp) == CTYPE_WHAMHA) ||
                 (CompType(iComp) == CTYPE_WHAMFA)) {
      CompTotStoich = 0;
      CompStoichMoles = 0.0;
      for (iSpec = 0; iSpec < NSpec; iSpec++) {
        CompTotStoich += SpecStoich(iSpec, iComp);
        CompStoichMoles += SpecStoich(iSpec, iComp) * SpecCtoM(iSpec);
      }
      CompConc(iComp) = TotConc(iComp) / CompStoichMoles;
    } else */if ((CompType(iComp) == CTYPE_DONNANHA) ||
               (CompType(iComp) == CTYPE_DONNANFA)) {
      CompConc(iComp) = 10.0;
    } else {
      CompConc(iComp) = TotConc(iComp);
    }
    TotMoles(iComp) = TotConc(iComp) * SpecCtoM(iComp);
  }

  /* Perform three rounds of adjustments */
  for (iRound = 1; iRound <= 3; iRound++) {
	  /* Calc Species */
    //CalcSpecConc(NComp, NSpec, CompConc, SpecK, SpecStoich, SpecName, 
    //             SpecActCorr, SpecActivityCoef, DoWHAM, SpecCharge, WHAMSpecCharge);
	  for (iSpec = 0; iSpec < NSpec; iSpec ++) {
      SpecConc(iSpec) = SpecK(iSpec);
      for (iComp = 0; iComp < NComp; iComp ++) {
        SpecConc(iSpec) = SpecConc(iSpec) * pow(CompConc(iComp), SpecStoich(iSpec, iComp));
      }
	  }
    /* Calc totals */
    CalcIterationTotals(NComp, NSpec, SpecConc, SpecCtoM, SpecStoich, 
                        CalcTotMoles, CalcTotConc);
    /*for (iComp = 0; iComp < NComp; iComp++) {
      CalcTotConc(iComp) = 0;
      for (iSpec = 0; iSpec < NSpec; iSpec++) {
        if (SpecStoich(iSpec, iComp) != 0){
            CalcTotConc(iComp) += SpecConc(iSpec) * SpecStoich(iSpec, iComp);
        }
      }
    }*/

    if (DoTox) {
      // Sum toxic BL-bound metal species
      CACalculated = 0;
      for (i = 0; i < NBLMetal; i++){
        iSpec = BLMetalSpecs[i];
        CACalculated += SpecConc[iSpec];
      }
      CompConc(MetalComp) = CompConc(MetalComp) * (CATarget / CACalculated);
      TotConc(MetalComp) = TotConc(MetalComp) * (CATarget / CACalculated);
      TotMoles(MetalComp) = TotMoles(MetalComp) * (CATarget / CACalculated);
    }

	  /* Adjust component concentrations */
    for (iComp = 0; iComp < NComp; iComp++) {
      /*if ((iRound == 1) && ((CompType(iComp) == CTYPE_DONNANHA) ||
                                   (CompType(iComp) == CTYPE_DONNANFA))) {
        CompConc(iComp) = 10.0;//CompConc(iComp) * (TotMoles(iComp) / CalcTotMoles(iComp) + 1) / 2;//
      } else */if ((CompType(iComp) == CTYPE_MASSBAL) || 
                 (CompType(iComp) == CTYPE_WHAMHA) ||
                 (CompType(iComp) == CTYPE_WHAMFA)) {
        //CompConc(iComp) = CompConc(iComp) * (TotConc(iComp) / CalcTotConc(iComp));
        CompConc(iComp) = CompConc(iComp) * (TotMoles(iComp) / CalcTotMoles(iComp));
      }
    }
  }
  return(CompConc);
}
