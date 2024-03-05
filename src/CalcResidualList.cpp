#include <Rcpp.h>
#include "ExternalCHESSFunctions.h"

//' @title Calculate the Residual
 //'
 //' @description
 //'   Calculate the residuals of the speciation problem - returns a list.
 //'
 //' @details The residuals of the speciation problem are the calculated total
 //'   moles for each component minus their known concentrations. In the case of
 //'   a toxicity run (DoTox = TRUE), the residual for the `MetalComp` component
 //'   is instead the sum of the `BLMetalSpecs` concentrations minus the
 //'   `CATarget`. This, paired with a modified set of derivatives for the
 //'   `MetalComp` in the Jacobian matrix, results in Newton-Rapshon changing the
 //'   metal concentration in order to produce chemical conditions that would
 //'   result in the critical accumulation known to cause a toxic effect.
 //'
 //' @param NComp integer, the number of components
 //' @param NSpec integer, the number of species
 //' @param SpecConc numeric vector (NSpec), the concentration of chemical
 //'   species
 //' @param SpecStoich integer matrix (NSpec x NComp), the stoichiometry of
 //'   species formatin reactions
 //' @param TotMoles numeric vector (NComp), the total moles of each component
 //' @param SpecCtoM numeric vector (NSpec), the concentration-to-mass conversion
 //'   factor for each chemical species
 //' @param CompName character vector (NComp), the names of the components
 //' @param CompType character vector (NComp), the type of component. It should
 //'   be a fixed set of values (MassBal, FixedConc, Substituted, ChargeBal,
 //'   SurfPot)
 //' @param MetalComp integer, the position in component vectors of the toxic
 //'   metal component
 //' @param BLMetalSpecs integer vector, the position in the species vectors of
 //'   the metal-biotic ligand species associated with toxicity
 //' @param CATarget numeric, the target critical accumulation value for a
 //'   toxicity run, in mol/kg
 //' @param DoTox logical, if TRUE = toxicity run, FALSE = speciation run
 //'
 //' @return A `list` object with the following components:
 //' \describe{
 //'  \item{\code{MaxError}}{numeric, the highest absolute error fraction in this
 //'   iteration =max(abs(Resid / TotMoles))}
 //'  \item{\code{WhichMax}}{integer, the position in the component vectors of the
 //'    component with the highest absolute error}
 //'  \item{\code{Resid}}{numeric vector (NComp), the residuals =
 //'    calculated totals - known totals}
 //'  \item{\code{CompError}}{numeric vector (NComp), the absolute error fraction
 //'    for each component in this iteration =abs(Resid / TotMoles)}
 //'  \item{\code{CalcTotConc}}{numeric vector (NComp), the calculated total
 //'    concentration of each component = CalcTotMoles / CtoM}
 //'  \item{\code{CalcTotMoles}}{numeric vector (NComp), the calculated total
 //'    moles of each component = sum(SpecConc * SpecCtoM * SpecStoich[,j]) for
 //'    each component j}
 //' }
 //'
 Rcpp::List CalcResidualList (int NComp,
                              int NSpec,
                              Rcpp::NumericVector SpecConc,
                              Rcpp::IntegerMatrix SpecStoich,
                              Rcpp::NumericVector TotMoles,
                              Rcpp::NumericVector SpecCtoM,
                              Rcpp::CharacterVector CompName,
                              Rcpp::CharacterVector CompType,
                              int MetalComp,
                              int NBLMetal,
                              Rcpp::IntegerVector BLMetalSpecs,
                              double CATarget,
                              bool DoTox) {
   /* outputs */
   double MaxError;// maximum of absolute ratios of residuals to totals
   int WhichMax;// which component has the highest absolute error
   Rcpp::NumericVector Resid(NComp);// residuals
   Rcpp::NumericVector CompError(NComp); // the absolute ratios of residuals to totals
   Rcpp::NumericVector CalcTotConc(NComp); // the calculated total concentrations (mol/L)
   Rcpp::NumericVector CalcTotMoles(NComp); // the calculated total concentrations (mol)
   Resid.names() = CompName;
   CompError.names() = CompName;
   CalcTotConc.names() = CompName;
   CalcTotMoles.names() = CompName;

   /* variables: */
   int iComp, iSpec, i;
   Rcpp::NumericVector SpecMoles;
   double CalcCA;
   Rcpp::NumericVector CompCtoM(NComp);

   // Rcpp::Rcout << "Calculate CalcTotMoles" << std::endl;
   SpecMoles = SpecConc * SpecCtoM;
   for (iComp = 0; iComp < NComp; iComp++) {
     CalcTotMoles(iComp) = 0;
     for (iSpec = 0; iSpec < NSpec; iSpec++) {
       if (SpecStoich(iSpec, iComp) != 0){
         CalcTotMoles(iComp) += SpecMoles(iSpec) * SpecStoich(iSpec, iComp);
       }
     }
   }

   // arma::mat SpecMolesMat = RcppVectorToMatrix(SpecMoles);
   // arma::mat SpecStoichMat = RcppIntMatrixToMatrix(SpecStoich);
   // arma::mat CalcTotMolesMat = SpecMolesMat * SpecStoichMat;
   // CalcTotMoles = MatrixToRcppMatrix(CalcTotMolesMat);

   // Rcpp::Rcout << "Calculate Resid" << std::endl;
   Resid = CalcTotMoles - TotMoles;
   for (iComp = 0; iComp < NComp; iComp++) {
     if ((CompType(iComp) == "FixedConc") || (CompType(iComp) == "FixedAct")) {
       Resid(iComp) = 0.0;
     }
   }

   // Rcpp::Rcout << "Calculate CompError" << std::endl;
   CompError = abs(Resid / TotMoles);
   if (DoTox) {
     // Rcpp::Rcout << "DoTox" << std::endl;
     // CalcCA = sum(SpecConc(BLMetalSpecs - 1));
     CalcCA = 0;
     for (i = 0; i < NBLMetal; i++) {
       iSpec = BLMetalSpecs(i) - 1;
       CalcCA += SpecConc(iSpec);
     }
     Resid(MetalComp - 1) = CalcCA - CATarget;
     // sum( (mol/L or mol/kgww) ) - mol = sum(mol) - mol = mol
     CompError(MetalComp - 1) = abs(Resid(MetalComp - 1) / CATarget);
   }


   // Rcpp::Rcout << "Calculate MaxError" << std::endl;
   MaxError = CompError(0);
   WhichMax = 0;
   for (iComp = 1; iComp < NComp; iComp++){
     if (CompError(iComp) > MaxError){
       MaxError = CompError(iComp);
       WhichMax = iComp;
     }
   }


   // Rcpp::Rcout << "Calculate CalcTotConc" << std::endl;
   for (iComp = 0; iComp < NComp; iComp++){
     CompCtoM(iComp) = SpecCtoM(iComp);
   }
   CalcTotConc = CalcTotMoles / CompCtoM;

   return Rcpp::List::create(
     Rcpp::Named("Resid") = Resid,
     Rcpp::Named("MaxError") = MaxError,
     Rcpp::Named("WhichMax") = WhichMax,
     Rcpp::Named("CompError") = CompError,
     Rcpp::Named("CalcTotConc") = CalcTotConc,
     Rcpp::Named("CalcTotMoles") = CalcTotMoles
   );
 }
