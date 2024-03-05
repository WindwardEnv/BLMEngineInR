#include <Rcpp.h>
#include "CHESSFunctions.h"

//' @title CHESSIter
//' 
//' @description Perform calculations for one iteration of CHESS
//' 
//' @details
//'   This function will take a CompConcStep input and will update all of the
//'   necessary calculations to return the residual and error information that
//'   results from that step. The updated variables and arrays are also
//'   returned. This function is meant to consolidate the calculations needed
//'   for a CHESS iteration into one function, so if other methods that are not
//'   strictly Newton-Raphson will be used, they can be without repeating (and
//'   potentially messing up) the code. This also helps to make it clearer what
//'   variables are updated with each iteration so that they won't be
//'   overwritten in the main routine unless specifically asked to. 
//' 
//' @author Kelly Croteau (kellyc@windwardenv.com)
//' 
//' @param CompConcStep NumericVector, length of NComp, the number to substract
//'   from CompConc in this iteration
//' @param NMass,MassAmt mass compartment information
//' @param NComp,CompName,CompType,CompPosInSpec component information
//' @param NSpec,SpecName,SpecMC,SpecActCorr,SpecStoich species info
//' @param SpecCharge,SpecKTempAdj more species information
//' @param DoWHAM boolean, whether WHAM organic matter is in this simulation
//' @param AqueousMC,WHAMDonnanMC which mass compartments correspond to the
//'   aqueous solution and Donnan layers, respectively
//' @param SolHS,wMolWt,wRadius,wP,wDLF,wKZED WHAM parameters
//' @param SysTempKelvin double, the temperature of the solution, in Kelvin
//' @param DoTox boolean, is this a toxicity run?
//' @param MetalName,MetalCompNBLMetal,BLMetalSpecs,CATarget params for tox mode
//' @param MassAmtAdj (OUTPUT), NumericVector, The amount of each mass 
//'   compartment, adjusted for the Donnan Layer volumes (aqueous compartment 
//'   will be its original value minus the total volume of the Donnan Layers).
//' @param TotConc (INPUT & OUTPUT) NumericVector, length of NComp, the total
//'   concentrations of each component in the simulation (units of e.g., mol/L
//'   and mol/kg)
//' @param SpecKISTempAdj (INPUT & OUTPUT) NumericVector, length NSpec, the
//'   equilibrium coefficient of the formation reactions, adjusted for
//'   temperature and ionic strength effects
//' @param SpecCtoMAdj (INPUT & OUTPUT) NumericVector, length NSpec, the
//'   concentration to mass conversion factor of the chemical species for which
//'   we have formation reactions, adjusted to take into account diffuse binding
//'   to organic matter.
//' @param SpecConc (INPUT & OUTPUT) NumericVector, length NSpec, the
//'   concentrations of each species for which we have formation reactions
//' @param CalcTotMoles (OUTPUT) NumericVector, length of NComp, the calculated
//'   total moles of each component in the simulation
//' @param WhichMax (OUTPUT) integer, the position in the component
//'   vectors of the component with the highest absolute error
//' @param IonicStrength (OUTPU) double, the ionic strength of the solution
//' @param Resid (OUTPUT) numeric vector (NComp), the residuals = 
//'   calculated totals - known totals
//' @param CompError (OUTPUT) numeric vector (NComp), the absolute error
//'   fraction for each component in this iteration = abs(Resid / TotMoles)
//' 
//' @return MaxError, double, the highest absolute error fraction in this
//'   iteration =max(abs(Resid / TotMoles))
//' 
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
  Rcpp::NumericVector &SpecActivityCoef,
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
  //Rcpp::NumericVector SpecActivityCoef(NSpec);
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
  
  for(int iComp = 0; iComp < NComp; iComp++) {
    if (CompType[iComp] == "FixedAct") {
      SpecConc[iComp] = TotConc[iComp] / SpecActivityCoef[iComp];
      CompConc[iComp] = SpecConc[iComp];
    } else if (CompType[iComp] == "FixedConc"){
      SpecConc[iComp] = TotConc[iComp];
      CompConc[iComp] = SpecConc[iComp];
    }
  }

  // Calculate the species concentrations
  SpecConc = CalcSpecConc(NComp, NSpec, CompConc, SpecKISTempAdj, SpecStoich, 
                          SpecName, SpecActCorr, SpecActivityCoef);
  
  if (DoWHAM) {
    AdjustForWHAM(NMass, MassAmt, MassAmtAdj, 
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
  /*UpdateTotals(NComp, NSpec, DoTox, CompType, CompName, MetalName,
                SpecStoich, SpecConc, CompCtoMAdj, 
                TotMoles, TotConc);*/
  
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