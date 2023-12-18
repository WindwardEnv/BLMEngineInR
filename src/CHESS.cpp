#include <Rcpp.h>
#include "CHESSFunctions.h"

//' @title CHemical Equilibria in Soils and Solutions
//'
//' @description Given a chemical system, equilibria equations, and total
//'   concentrations of components, calculate the species concentrations of each
//'   chemical product in the system.
//'
//' @param QuietFlag character, one of "Very Quiet" (only print out when run is
//'   done), "Quiet" (print out Obs=iObs), or "Debug" (print out lots of info)
//' @param ConvergenceCriteria numeric, the maximum value of MaxError that counts
//'   as convergence by the Newton-Raphson root-finding algorithm
//' @param MaxIter integer, the maximum number of iterations the Newton-Raphson
//'   root-finding algorithm should do before giving up
//' @param NComp integer, number of components
//' @param NSpec integer, number of species reactions
//' @param NBLMetal integer, the number of biotic ligand-bound metal species that
//'   are associated with toxic effects.
//' @param SpecK numeric vector (NSpec), the equilibrium coefficient of the
//'   formation reactions
//' @param SpecStoich signed integer matrix (NSpec x NComp), the reaction
//'   stoichiometry of the formation reactions
//' @param SpecCtoM numeric vector (NSpec), the concentration to mass conversion
//'   factor of the chemical species for which we have formation reactions
//' @param SpecName character vector (NSpec), the name of the chemical species
//'   for which we have formation reactions
//' @param CompType character vector (NComp), the type of each component in the
//'   simulation
//' @param CompName character vector (NComp), the name of each component in the
//'   simulation
//' @param TotMoles numeric vector (NComp), the total moles of each component in
//'   the simulation (units of mol)
//' @param TotConc numeric vector (NComp), the total concentrations of each
//'   component in the simulation (units of e.g., mol/L and mol/kg)
//' @param DoTox logical, TRUE for toxicity mode where the MetalName component
//'   concentration is adjusted to try to match the CATarget with BLMetalSpecs
//' @param MetalName character string, the name of the toxic metal
//' @param MetalComp integer, the position of the metal in the component arrays
//'   (i.e., which is the toxic metal component) Note: this are base-1 indexed.
//' @param BLMetalSpecs integer vector, the positions of the species in the
//'   arrays which contribute to toxicity (i.e., which species are the toxic
//'   metal bound to the relevant biotic ligand) Note: these are base-1 indexed.
//' @param CATarget numeric, the target critical accumulation in units of mol /
//'   kg (only used when DoTox == TRUE)
//'
//' @return list with the following elements:
//'   \desribe{
//'     \item{SpecConc}{numeric vector (NSpec), the concentrations of each
//'       species for which we have formation reactions}
//'     \item{Iter}{integer, the number of Newton-Raphson iterations that we
//'       needed to reach convergence}
//'     \item{MaxError}{numeric, the highest final absolute error fraction
//'       =max(abs(Resid / TotMoles))}
//'     \item{CalcTotConc}{numeric vector (NComp), the calculated total
//'       concentrations of each component in the simulation (units of e.g.,
//'       mol/L and mol/kg)}
//'   }
//' @export
//'
//[[Rcpp::export]]
Rcpp::List CHESS(Rcpp::String QuietFlag,
                 double ConvergenceCriteria,
                 unsigned int MaxIter,
                 unsigned int NComp,
                 unsigned int NSpec,
                 unsigned int NBLMetal,
                 Rcpp::NumericVector SpecK,
                 Rcpp::NumericVector SpecTemp,
                 Rcpp::NumericVector SpecDeltaH,
                 Rcpp::IntegerMatrix SpecStoich,
                 Rcpp::NumericVector SpecCtoM,
                 Rcpp::CharacterVector SpecName,
                 Rcpp::CharacterVector CompType,
                 Rcpp::CharacterVector CompName,
                 Rcpp::NumericVector TotMoles,
                 Rcpp::NumericVector TotConc,
                 double SysTemp,
                 bool DoTox,
                 Rcpp::String MetalName,
                 unsigned int MetalComp,
                 Rcpp::IntegerVector BLMetalSpecs,
                 double CATarget) {

  /*outputs*/
  Rcpp::NumericVector SpecConc(NSpec); // species concentrations after optimization
  unsigned int Iter = 0;
  double MaxError;
  Rcpp::NumericVector CalcTotConc(NComp); //the calculated total concentrations of each component in the simulation (units of e.g., mol/L and mol/kg)}

  /*variables*/
  Rcpp::NumericMatrix JacobianMatrix(NComp);
  Rcpp::NumericVector CompConcStep(NComp);
  Rcpp::NumericVector CompConc(NComp);
  Rcpp::NumericVector SpecKTempAdj(NSpec);
  Rcpp::List ResidResults;
  unsigned int WhichMax;
  Rcpp::NumericVector Resid (NComp);
  Rcpp::NumericVector CompError(NComp); // the absolute ratios of residuals to totals
  Rcpp::NumericVector CalcTotMoles(NComp); // the calculated total concentrations (mol)
  Resid.names() = CompName;
  CompError.names() = CompName;
  CalcTotConc.names() = CompName;
  CalcTotMoles.names() = CompName;

   SpecKTempAdj = TempCorrection(SysTemp, NSpec, SpecK, SpecTemp, SpecDeltaH);

   // Get initial values for component concentrations
   CompConc = InitialGuess(TotConc, CompType, SpecKTempAdj, SpecStoich,
                           SpecName, NComp, NSpec);

   // Initialize Species Concentrations
   SpecConc = CalcSpecConc(CompConc, SpecKTempAdj, SpecStoich, SpecName,
                           NComp, NSpec);

   // Update Total Concentrations for Fixed Activity & Metal
   UpdateTotals(NComp, NSpec, CompType, CompName, MetalName, TotMoles, 
                SpecStoich, (SpecConc * SpecCtoM), TotConc, SpecCtoM, DoTox);
   
  // Calculate Residuals for the first time...

  // Calculate the total moles & conc from species concentrations
  CalcIterationTotals(NComp, NSpec, SpecConc, SpecCtoM, SpecStoich, 
                      CalcTotMoles, CalcTotConc);

  // Calculate the residuals and error fraction for each component
  CalcResidAndError(NComp, CalcTotMoles, TotMoles, CompType, Resid, CompError);

  // Adjust Resid and CompError for toxicity mode
  if (DoTox) {
    AdjustForToxMode(NBLMetal, BLMetalSpecs, MetalComp, CATarget, SpecConc, 
                     Resid, CompError);
  }

  // Determine which component has the highest error fraction
  MaxError = MaxCompError(NComp, CompError, WhichMax);

  // Begin iterating
  Iter = 0;
  while ((MaxError > ConvergenceCriteria) & (Iter <= MaxIter)) {

    // update the iteration counter
    Iter++;

    // Calculate the Jacobian Matrix
    JacobianMatrix = Jacobian(NComp, NSpec, SpecStoich, SpecConc, SpecCtoM, 
                              CompName, MetalComp, NBLMetal, BLMetalSpecs, DoTox);

    // Calculate the Newton-Raphson step
    CompConcStep = CalcStep(JacobianMatrix = JacobianMatrix, Resid = Resid,
                            NComp = NComp, CompType = CompType, CompName=CompName);

    //Update the component free concentrations
    CompConc = SpecConc.import(SpecConc.begin(), SpecConc.begin() + NComp);
    CompUpdate(NComp, CompConcStep, CompConc);

    // Calculate the species concentrations
    SpecConc = CalcSpecConc(CompConc, SpecKTempAdj, SpecStoich, SpecName, NComp, 
                            NSpec);

    // Update Total Concentrations for Fixed Activity & Metal
    UpdateTotals(NComp, NSpec, CompType, CompName, MetalName, TotMoles,
                SpecStoich, (SpecConc * SpecCtoM), TotConc, SpecCtoM, DoTox);
  
    // Calculate the total moles & conc from species concentrations
    CalcIterationTotals(NComp, NSpec, SpecConc, SpecCtoM, SpecStoich, 
                        CalcTotMoles, CalcTotConc);

    // Calculate the residuals and error fraction for each component
    CalcResidAndError(NComp, CalcTotMoles, TotMoles, CompType, Resid, CompError);

    // Adjust Resid and CompError for toxicity mode
    if (DoTox) {
      AdjustForToxMode(NBLMetal, BLMetalSpecs, MetalComp, CATarget, SpecConc, 
                      Resid, CompError);
    }

    // Determine which component has the highest error fraction
    MaxError = MaxCompError(NComp, CompError, WhichMax);

    if(QuietFlag == "Debug"){
      Rcpp::Rcout << "Iter=" << Iter << 
                  ", WhichMax=" << CompName(WhichMax) <<
                  ", SpecConc[WhichMax]= " << SpecConc(WhichMax) <<
                  ", CalcTotConc[WhichMax]= " << CalcTotConc(WhichMax) <<
                  ", MaxError=" << MaxError << std::endl;
    }

  }//while ((MaxError > ConvergenceCriteria) & (Iter <= MaxIter))

  return Rcpp::List::create(
      Rcpp::Named("SpecConc") = SpecConc,
      Rcpp::Named("Iter") = Iter,
      Rcpp::Named("MaxError") = MaxError,
      Rcpp::Named("CalcTotConc") = CalcTotConc
  );

}


