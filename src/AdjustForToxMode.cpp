#include <math.h>
#include <Rcpp.h>
#include "CHESSFunctions.h"

//' @title Adjust the residual and error for toxicity mode
//'
//' @description Calculate a new residual and error value for the toxic metal
//'   component, based instead on the difference between the target and 
//'   calculated critical accumulations.
//'
//' @details This should only be run when DoTox = TRUE. The residual for the
//'   `MetalComp` component is instead the sum of the `BLMetalSpecs`
//'   concentrations minus the `CATarget`. This, paired with a modified set of 
//'   derivatives for the `MetalComp` in the Jacobian matrix, results in 
//'   Newton-Rapshon changing the metal concentration in order to produce
//'   chemical conditions that would result in the critical accumulation known
//'   to cause a toxic effect.
//'
//' @param NBLMetal integer, the number of biotic ligand-bound metal species
//' @param BLMetalSpecs integer vector, the position in the species vectors of
//'   the metal-biotic ligand species associated with toxicity
//' @param MetalComp integer, the position in component vectors of the toxic
//'   metal component
//' @param CATarget numeric, the target critical accumulation value for a
//'   toxicity run, in mol/kg
//' @param SpecConc numeric vector (NSpec), the concentration of chemical
//'   species
//' @param Resid (return value) numeric vector (NComp), the residuals =
//'   calculated totals - known totals, modified for toxicity mode upon return
//' @param CompError (return value) numeric vector (NComp), the absolute error
//'   fraction for each component in this iteration = abs(Resid / TotMoles),
//'   modified for toxicity mode upon return
//'
//' @name AdjustForToxMode
//' @usage AdjustForToxMode(NBLMetal, BLMetalSpecs, MetalComp, CATarget, 
//'   SpecConc, Resid, CompError)
void AdjustForToxMode(int NBLMetal, 
                      Rcpp::IntegerVector BLMetalSpecs, 
                      int MetalComp,
                      double CATarget,
                      Rcpp::NumericVector SpecConc,
                      Rcpp::NumericVector &Resid,
                      Rcpp::NumericVector &CompError) {

  /* variables */
  double CACalc;
  
  // Adjust Resid and CompError for toxicity mode
  CACalc = CalcCA(NBLMetal, BLMetalSpecs, SpecConc);

  // Resid and CompError for the metal component are based on these species
  Resid[MetalComp] = CACalc - CATarget;
  CompError[MetalComp] = std::fabs(Resid[MetalComp] / CATarget);

}

void CalcToxError (int NBLMetal, 
                   Rcpp::IntegerVector BLMetalSpecs, 
                   double CATarget,
                   Rcpp::NumericVector SpecConc,
                   double &ToxResid,
                   double &ToxError) {
  /* variables */
  double CACalc;
  
  CACalc = CalcCA(NBLMetal, BLMetalSpecs, SpecConc);
  ToxResid = CACalc - CATarget;
  ToxError = std::fabs(ToxResid) / CATarget;
}

double CalcCA(int NBLMetal, 
              Rcpp::IntegerVector BLMetalSpecs, 
              Rcpp::NumericVector SpecConc) {

  /* variables */
  double CACalc;
  int i, iSpec;

  // Sum toxic BL-bound metal species
  CACalc = 0;
  for (i = 0; i < NBLMetal; i++) {
    iSpec = BLMetalSpecs[i];
    CACalc += SpecConc[iSpec];
  }

  return CACalc;
}