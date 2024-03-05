#include <Rcpp.h>
#include "CHESSFunctions.h"

double CHESSIter(
  Rcpp::NumericVector CompConcStep,
  int NMass,
  Rcpp::NumericVector MassAmt,
  int NComp,
  Rcpp::CharacterVector CompName,
  Rcpp::CharacterVector CompType,
  Rcpp::IntegerVector CompPosInSpec,
  int NSpec,
  Rcpp::CharacterVector SpecName,
  Rcpp::IntegerVector SpecMC,
  Rcpp::CharacterVector SpecActCorr,
  Rcpp::IntegerMatrix SpecStoich,
  Rcpp::IntegerVector SpecCharge,
  Rcpp::NumericVector SpecKTempAdj,
  bool DoWHAM,
  int AqueousMC,
  Rcpp::IntegerVector WHAMDonnanMC,
  Rcpp::NumericVector SolHS,
  Rcpp::NumericVector wMolWt,
  Rcpp::NumericVector wRadius,
  Rcpp::NumericVector wP,
  double wDLF,
  double wKZED,
  double SysTempKelvin,
  bool DoTox,
  Rcpp::String MetalName,
  int MetalComp,
  int NBLMetal,
  Rcpp::IntegerVector BLMetalSpecs,
  double CATarget,
  Rcpp::NumericVector &MassAmtAdj,
  Rcpp::NumericVector &TotConc,
  Rcpp::NumericVector &SpecKISTempAdj,
  Rcpp::NumericVector &SpecCtoMAdj,
  Rcpp::NumericVector &SpecConc,
  Rcpp::NumericVector &CalcTotMoles,
  int &WhichMax,
  double &IonicStrength,
  Rcpp::NumericVector &Resid,
  Rcpp::NumericVector &CompError
) {
  
  /*outputs*/
  double MaxError;
  Rcpp::NumericVector CalcTotConc(NComp);
    CalcTotConc.names() = CompName;
  Rcpp::NumericVector SpecMoles(NSpec);
    SpecMoles.names() = SpecName;
      
  /*variables*/
  //Rcpp::NumericVector MassAmtAdj(NMass);
  //  MassAmtAdj.names() = MassAmt.names();
  Rcpp::NumericVector CompConc(NComp);
  Rcpp::NumericVector TotMoles(NComp);
    TotMoles.names() = CompName;
  Rcpp::NumericVector CompCtoMAdj(NComp);
  //Rcpp::NumericVector SpecKISTempAdj(NSpec);
  //Rcpp::NumericVector SpecCtoMAdj(NSpec);
  Rcpp::NumericVector SpecActivityCoef(NSpec);
  Rcpp::NumericVector WHAMSpecCharge(2);
  
  // Update the component free concentrations
  CompConc = SpecConc[CompPosInSpec];
  CompUpdate(NComp, CompConcStep, CompType, CompConc);
  SpecConc[CompPosInSpec] = clone(CompConc);

  // Calculate the ionic strength and activity coefficients
  IonicStrength = CalcIonicStrength(NSpec, SpecConc * SpecCtoMAdj, SpecCharge, 
                                    SpecMC, AqueousMC);
  SpecActivityCoef = CalcActivityCoef(NSpec, SpecName, SpecActCorr, SpecCharge, 
                                      IonicStrength, SysTempKelvin);
  
  // Calculate the species concentrations
  SpecConc = CalcSpecConc(NComp, NSpec, CompConc, SpecKISTempAdj, SpecStoich, 
                          SpecName, SpecActCorr, SpecActivityCoef);
  
  if (DoWHAM) {

    WHAMAdjustments(NMass, MassAmt, MassAmtAdj, 
                    NComp, CompType, TotConc, TotMoles,
                    NSpec, SpecConc, SpecMC, SpecActCorr, SpecCharge, 
                    SpecKTempAdj, SpecKISTempAdj, SpecCtoMAdj,  
                    IonicStrength, WHAMSpecCharge, AqueousMC, WHAMDonnanMC,
                    SolHS, wMolWt, wRadius, wP, wDLF, wKZED);
    
  }       
  SpecMoles = SpecConc * SpecCtoMAdj;
  CompCtoMAdj = SpecCtoMAdj[CompPosInSpec];
  TotMoles = TotConc * CompCtoMAdj;

  // Update Total Concentrations for Fixed Activity & Metal
  UpdateTotals(NComp, NSpec, DoTox, CompType, CompName, MetalName,
                SpecStoich, (SpecConc * SpecCtoMAdj), CompCtoMAdj, 
                TotMoles, TotConc);

  // Calculate the total moles & conc from species concentrations
  CalcTotMoles = CalcIterationTotalMoles(NComp, NSpec, SpecConc * SpecCtoMAdj, 
                                         SpecStoich);

  // Calculate the residuals and error fraction for each component
  CalcResidAndError(NComp, CalcTotMoles, TotMoles, CompType, Resid, 
                    CompError);

  // Adjust Resid and CompError for toxicity mode
  if (DoTox) {
    AdjustForToxMode(NBLMetal, BLMetalSpecs, MetalComp, CATarget, SpecConc,
                      Resid, CompError);
  }

  // Determine which component has the highest error fraction
  MaxError = MaxCompError(NComp, CompError, WhichMax);

  return MaxError;

}