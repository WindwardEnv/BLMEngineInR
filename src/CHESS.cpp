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
                 Rcpp::IntegerMatrix SpecStoich,
                 Rcpp::NumericVector SpecCtoM,
                 Rcpp::CharacterVector SpecName,
                 Rcpp::CharacterVector CompType,
                 Rcpp::CharacterVector CompName,
                 Rcpp::NumericVector TotMoles,
                 Rcpp::NumericVector TotConc,
                 bool DoTox,
                 Rcpp::String MetalName,
                 unsigned int MetalComp,
                 Rcpp::IntegerVector BLMetalSpecs,
                 double CATarget) {

   /*outputs*/
   Rcpp::NumericVector SpecConc(NSpec); // species concentrations after optimization
   unsigned int Iter = 0;
   double MaxError;
   Rcpp::NumericVector CalcTotConc; //the calculated total concentrations of each component in the simulation (units of e.g., mol/L and mol/kg)}

   /*variables*/
   Rcpp::NumericMatrix JacobianMatrix(NComp);
   Rcpp::NumericVector CompConcStep(NComp);
   Rcpp::NumericVector CompConc(NComp);
   Rcpp::List ResidResults;
   Rcpp::NumericVector Resid (NComp);
   unsigned int WhichMax;
   Rcpp::NumericVector CalcTotMoles(NComp);

   // Get initial values for component concentrations
   CompConc = InitialGuess(TotConc = TotConc,
                          CompType = CompType,
                          SpecK = SpecK,
                          SpecStoich = SpecStoich,
                          SpecName = SpecName,
                          NComp = NComp,
                          NSpec = NSpec);

   // Initialize Species Concentrations
   SpecConc = CalcSpecConc(CompConc = CompConc,
                           SpecK = SpecK,
                           SpecStoich = SpecStoich,
                           SpecName = SpecName,
                           NComp = NComp,
                           NSpec = NSpec);

   // Update Total Concentrations for Fixed Activity & Metal
   UpdateTotals(NComp, NSpec, CompType, CompName, MetalName, TotMoles, 
                SpecStoich, (SpecConc * SpecCtoM), TotConc, SpecCtoM, DoTox);
   
   // Calculate Residuals for the first time
   ResidResults = CalcResidual(
     NComp = NComp,
     NSpec = NSpec,
     SpecConc = SpecConc,
     SpecStoich = SpecStoich,
     TotMoles = TotMoles,
     SpecCtoM = SpecCtoM,
     CompName = CompName,
     CompType = CompType,
     MetalComp = MetalComp,
     NBLMetal = NBLMetal,
     BLMetalSpecs = BLMetalSpecs,
     CATarget = CATarget,
     DoTox = DoTox
   );
    Resid = ResidResults["Resid"];
    MaxError = ResidResults["MaxError"];
    WhichMax = ResidResults["WhichMax"];
    CalcTotConc = ResidResults["CalcTotConc"];
    CalcTotMoles = ResidResults["CalcTotMoles"];

   // Begin iterating
   Iter = 0;
   while ((MaxError > ConvergenceCriteria) & (Iter <= MaxIter)){

    Iter++;

    JacobianMatrix = Jacobian(NComp = NComp, NSpec = NSpec, SpecStoich = SpecStoich,
                              SpecConc = SpecConc, SpecCtoM = SpecCtoM, CompName = CompName,
                              MetalComp = MetalComp, NBLMetal = NBLMetal,
                              BLMetalSpecs = BLMetalSpecs, DoTox = DoTox);

    CompConcStep = CalcStep(JacobianMatrix = JacobianMatrix, Resid = Resid,
                            NComp = NComp, CompType = CompType, CompName=CompName);

    CompConc = CompUpdate(
      NComp = NComp,
      CompConcStep = CompConcStep,
      CompConc = SpecConc.import(SpecConc.begin(), SpecConc.begin() + NComp),
      CompName = CompName);

    SpecConc = CalcSpecConc(CompConc = CompConc,
                            SpecK = SpecK,
                            SpecStoich = SpecStoich,
                            SpecName = SpecName,
                            NComp = NComp,
                            NSpec = NSpec);
    // Update Total Concentrations for Fixed Activity & Metal
    UpdateTotals(NComp, NSpec, CompType, CompName, MetalName, TotMoles,
                 SpecStoich, (SpecConc * SpecCtoM), TotConc, SpecCtoM, DoTox);

    ResidResults = CalcResidual(
        NComp = NComp,
        NSpec = NSpec,
        SpecConc = SpecConc,
        SpecStoich = SpecStoich,
        TotMoles = TotMoles,
        SpecCtoM = SpecCtoM,
        CompName = CompName,
        CompType = CompType,
        MetalComp = MetalComp,
        NBLMetal = NBLMetal,
        BLMetalSpecs = BLMetalSpecs,
        CATarget = CATarget,
        DoTox = DoTox
    );
    Resid = ResidResults["Resid"];
    MaxError = ResidResults["MaxError"];
    WhichMax = ResidResults["WhichMax"];
    CalcTotConc = ResidResults["CalcTotConc"];
    CalcTotMoles = ResidResults["CalcTotMoles"];

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


