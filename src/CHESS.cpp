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
//' @param ConvergenceCriteria numeric, the maximum value of MaxError that 
//'   counts as convergence by the Newton-Raphson root-finding algorithm
//' @param MaxIter integer, the maximum number of iterations the Newton-Raphson
//'   root-finding algorithm should do before giving up
//' @param NComp integer, number of components
//' @param NSpec integer, number of species reactions
//' @param NBLMetal integer, the number of biotic ligand-bound metal species 
//'   that are associated with toxic effects.
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
//' \describe{
//'   \item{SpecConc}{numeric vector (NSpec), the concentrations of each species
//'     for which we have formation reactions}
//'   \item{FinalIter}{integer, the number of Newton-Raphson iterations that we
//'     needed to reach convergence}
//'   \item{FinalMaxError}{numeric, the highest final absolute error fraction  
//'     =max(abs(Resid / TotMoles))}
//'   \item{CalcTotConc}{numeric vector (NComp), the calculated total
//'     concentrations of each component in the simulation (units of e.g., mol/L
//'     and mol/kg)}
//' }
//' @export
//'
//[[Rcpp::export]]
Rcpp::List CHESS(Rcpp::String QuietFlag,
                 double ConvergenceCriteria,
                 unsigned int MaxIter,
                 unsigned int NMass,
                 unsigned int NComp,
                 unsigned int NSpec,
                 unsigned int NBLMetal,
                 Rcpp::IntegerVector SpecMC,
                 Rcpp::NumericVector SpecK,
                 Rcpp::NumericVector SpecTempKelvin,
                 Rcpp::NumericVector SpecDeltaH,
                 Rcpp::IntegerMatrix SpecStoich,
                 Rcpp::IntegerVector SpecCharge,
                 Rcpp::CharacterVector SpecActCorr,
                 Rcpp::NumericVector SpecCtoM,
                 Rcpp::CharacterVector SpecName,
                 Rcpp::CharacterVector CompType,
                 Rcpp::CharacterVector CompName,
                 Rcpp::NumericVector TotConc,
                 bool DoWHAM,
                 int AqueousMC,
                 Rcpp::NumericVector SolHS,
                 Rcpp::NumericVector wMolWt,
                 Rcpp::NumericVector wRadius,
                 Rcpp::NumericVector wP,
                 double wDLF,
                 double wKZED,
                 double SysTempKelvin,
                 bool DoTox,
                 Rcpp::String MetalName,
                 unsigned int MetalComp,
                 Rcpp::IntegerVector BLMetalSpecs,
                 double CATarget) {

  /*outputs*/
  Rcpp::NumericVector SpecConc(NSpec); // species concentrations after optimization
  Rcpp::NumericVector SpecAct(NSpec);
  unsigned int Iter = 0;  
  double MaxError;
  Rcpp::NumericVector CalcTotConc(NComp); // the calculated total concentrations of each component in the simulation (units of e.g., mol/L and mol/kg)}

  /*variables*/
  Rcpp::NumericMatrix JacobianMatrix(NComp);
  Rcpp::NumericVector CompConcStep(NComp);
  Rcpp::NumericVector CompConc(NComp);
  Rcpp::NumericVector SpecKTempAdj(NSpec);
  Rcpp::NumericVector SpecKISTempAdj(NSpec);
  Rcpp::NumericVector SpecCtoMAdj = clone(SpecCtoM);
  Rcpp::NumericVector SpecMoles(NSpec);
  Rcpp::NumericVector TotMoles(NComp);
  Rcpp::NumericVector CompCtoM(NComp);
  Rcpp::NumericVector CompCtoMAdj(NComp);
  Rcpp::List ResidResults;
  unsigned int WhichMax;
  Rcpp::NumericVector Resid(NComp);
  Rcpp::NumericVector CompError(NComp);
  Rcpp::NumericVector CalcTotMoles(NComp);
  double IonicStrength;
  Rcpp::NumericVector SpecActivityCoef(NSpec);
  Rcpp::IntegerVector CompPosInSpec(NComp);
  Rcpp::NumericVector WHAMSpecCharge(2);
  unsigned int iComp;
  unsigned int iSpec;
  TotMoles.names() = CompName;
  Resid.names() = CompName;
  CompError.names() = CompName;
  CalcTotConc.names() = CompName;
  CalcTotMoles.names() = CompName;
  SpecAct.names() = SpecName;

  for (iComp = 0; iComp < NComp; iComp++){
    CompPosInSpec(iComp) = iComp;
  }

  CompCtoM = SpecCtoM[CompPosInSpec];
  TotMoles = TotConc * CompCtoM;

  SpecKTempAdj = TempCorrection(SysTempKelvin, NSpec, SpecK, SpecTempKelvin, 
                                SpecDeltaH);
  /*Rcpp::NumericVector TmpVector = SpecKTempAdj / SpecK;
  Rcpp::Rcout << TmpVector << std::endl;*/

  // Get initial values for component concentrations
  CompConc = InitialGuess(TotConc, SpecCtoMAdj, CompType, SpecKTempAdj, SpecStoich,
                          SpecName, NComp, NSpec);
  SpecConc[CompPosInSpec] = clone(CompConc);

  // Calculate the ionic strength and activity coefficients
  IonicStrength = CalcIonicStrength(NSpec, SpecConc * SpecCtoMAdj, SpecCharge, 
                                    SpecMC);
  SpecActivityCoef = CalcActivityCoef(NSpec, SpecName, SpecActCorr, SpecCharge, 
                                      IonicStrength, SysTempKelvin);

  // Initialize Species Concentrations
  SpecConc = CalcSpecConc(NComp, NSpec, CompConc, SpecKTempAdj, SpecStoich, 
                          SpecName, SpecActCorr, SpecActivityCoef);

  if (DoWHAM) {
    //Calculate the charge on the organic matter
    WHAMSpecCharge = CalcWHAMSpecCharge(NSpec, SpecActCorr, SpecConc, 
                                        SpecCharge, SpecMC, AqueousMC);

    //Adjust organic matter specific binding based on ionic strength
    IonicStrength = CalcIonicStrength(NSpec, SpecConc * SpecCtoMAdj, SpecCharge, 
                                      SpecMC);
    SpecKISTempAdj = CalcIonicStrengthEffects(IonicStrength, WHAMSpecCharge, 
                                              NSpec, SpecCharge, SpecKTempAdj, 
                                              SpecActCorr, wP);
    
    //Adjust diffuse binding 
    SpecCtoMAdj = CalcDonnanLayerVolume(NSpec, IonicStrength, SpecCtoM, 
                                        SpecActCorr, SpecMC, AqueousMC, wMolWt, 
                                        wRadius, wDLF, wKZED, WHAMSpecCharge, 
                                        SolHS);

    AdjustForWHAM(NComp, NSpec, CompName, SpecActCorr, SpecCharge, 
                  WHAMSpecCharge, SpecCtoMAdj, SpecConc, TotConc, TotMoles);

  } else {
    SpecKISTempAdj = SpecKTempAdj;
    SpecCtoMAdj = SpecCtoM;
  }
  SpecMoles = SpecConc * SpecCtoMAdj;
  CompCtoMAdj = SpecCtoMAdj[CompPosInSpec];
  TotMoles = TotConc * CompCtoMAdj;

  // Update Total Concentrations for Fixed Activity & Metal
  UpdateTotals(NComp, NSpec, DoTox, CompType, CompName, MetalName,
               SpecStoich, (SpecConc * SpecCtoMAdj), CompCtoMAdj, 
               TotMoles, TotConc);

  // Calculate Residuals for the first time...

  // Calculate the total moles & conc from species concentrations
  CalcIterationTotals(NComp, NSpec, SpecConc, SpecCtoMAdj, SpecStoich,
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

    if (QuietFlag == "Debug") {
      Rcpp::Rcout << "Iter=" << Iter << 
        ", WhichMax=" << CompName(WhichMax) << 
        ", MaxError=" << MaxError << 
        ", Resid=" << Resid <<
        ", WHAMSpecCharge=" << WHAMSpecCharge <<
        ", TotMoles=" << TotMoles <<
        ", CalcTotConc=" << CalcTotConc <<
        ", CalcTotMoles=" << CalcTotMoles <<
        ", SpecConc=" << SpecConc << 
        ", SpecCtoMAdj=" << SpecCtoMAdj << 
        std::endl;
    }

    // update the iteration counter
    Iter++;

    // Calculate the Jacobian Matrix
    JacobianMatrix = Jacobian(NComp, NSpec, SpecStoich, SpecConc, SpecCtoMAdj,
                              CompName, MetalComp, NBLMetal, BLMetalSpecs, 
                              DoTox);

    // Calculate the Newton-Raphson step
    CompConcStep = CalcStep(JacobianMatrix, Resid, NComp, CompType, CompName);

    if (FALSE) {//(DoWHAM) {
      // Update the Donnan layers the way WHAM does it
      for (iComp = 0; iComp < NComp; iComp++) {
        if ((CompType(iComp) == "DonnanHA") || 
            (CompType(iComp) == "DonnanFA")) {
          CompConcStep(iComp) = CompConc(iComp) * 
            (1 - (TotMoles(iComp) / CalcTotMoles(iComp))) / 2;
        }
      }
    }

    // Update the component free concentrations
    CompConc = SpecConc[CompPosInSpec];
    CompUpdate(NComp, CompConcStep, CompConc);
    SpecConc[CompPosInSpec] = clone(CompConc);
    
    // Calculate the ionic strength and activity coefficients
    IonicStrength = CalcIonicStrength(NSpec, SpecConc * SpecCtoMAdj, SpecCharge, 
                                      SpecMC);
    SpecActivityCoef = CalcActivityCoef(NSpec, SpecName, SpecActCorr, SpecCharge, 
                                        IonicStrength, SysTempKelvin);

    // Calculate the species concentrations
    SpecConc = CalcSpecConc(NComp, NSpec, CompConc, SpecKISTempAdj, SpecStoich, 
                            SpecName, SpecActCorr, SpecActivityCoef);

    if (DoWHAM) {
      //Calculate the charge on the organic matter
      WHAMSpecCharge = CalcWHAMSpecCharge(NSpec, SpecActCorr,
                                          SpecConc * SpecCtoMAdj, 
                                          SpecCharge, SpecMC, AqueousMC);
      
      //Adjust organic matter specific binding based on ionic strength
      SpecKISTempAdj = CalcIonicStrengthEffects(IonicStrength, WHAMSpecCharge, 
                                                NSpec, SpecCharge, SpecKTempAdj, 
                                                SpecActCorr, wP);
      
      //Calculate the portion of the solution that's in the diffuse layer
      SpecCtoMAdj = CalcDonnanLayerVolume(NSpec, IonicStrength, SpecCtoM, 
                                          SpecActCorr, SpecMC, AqueousMC, wMolWt, 
                                          wRadius, wDLF, wKZED, WHAMSpecCharge, 
                                          SolHS);
      
      // Adjust species and component totals
      AdjustForWHAM(NComp, NSpec, CompName, SpecActCorr, SpecCharge, 
                    WHAMSpecCharge, SpecCtoMAdj, SpecConc, TotConc, TotMoles);
    
    } else {
      SpecKISTempAdj = SpecKTempAdj;
      SpecCtoMAdj = SpecCtoM;
    }       
    SpecMoles = SpecConc * SpecCtoMAdj;
    CompCtoMAdj = SpecCtoMAdj[CompPosInSpec];
    TotMoles = TotConc * CompCtoMAdj;

    // Update Total Concentrations for Fixed Activity & Metal
    UpdateTotals(NComp, NSpec, DoTox, CompType, CompName, MetalName,
                 SpecStoich, (SpecConc * SpecCtoMAdj), CompCtoMAdj, 
                 TotMoles, TotConc);

    // Calculate the total moles & conc from species concentrations
    CalcIterationTotals(NComp, NSpec, SpecConc, SpecCtoMAdj, SpecStoich,
                        CalcTotMoles, CalcTotConc);

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

  } // while ((MaxError > ConvergenceCriteria) & (Iter <= MaxIter))

  SpecMoles = SpecConc * SpecCtoMAdj;
  SpecAct = SpecConc * SpecActivityCoef;

  if (QuietFlag == "Debug") {
      Rcpp::Rcout << "Iter=" << Iter << 
        ", WhichMax=" << CompName(WhichMax) << 
        ", MaxError=" << MaxError << 
        ", WHAMSpecCharge=" << WHAMSpecCharge <<
        ", Resid=" << Resid <<
        ", SpecConc=" << SpecConc << 
        ", CalcTotConc=" << CalcTotConc << std::endl;
    }

  return Rcpp::List::create(
      Rcpp::Named("SpecConc") = SpecConc,
      Rcpp::Named("SpecAct") = SpecAct,
      Rcpp::Named("SpecMoles") = SpecMoles,
      Rcpp::Named("FinalIter") = Iter,
      Rcpp::Named("FinalMaxError") = MaxError,
      Rcpp::Named("CalcTotConc") = CalcTotConc);
}
