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
//' @param DoPartialStepsAlways boolean, Should CHESS do strict Newton-Raphson
//'   iterations (FALSE), or try to improve the simulation with partial N-R
//'   steps (trying to prevent oscillations).
//' @param NMass integer, number of mass compartments
//' @param MassName CharacterVector (NMass), the names of the mass compartments
//' @param MassAmt NumericVector (NMass), The amount of each mass compartment. 
//' @param NComp integer, number of components
//' @param CompName character vector (NComp), the name of each component in the
//'   simulation
//' @param CompType character vector (NComp), the type of each component in the
//'   simulation
//' @param TotConc numeric vector (NComp), the total concentrations of each
//'   component in the simulation (units of e.g., mol/L and mol/kg)
//' @param NSpec integer, number of species reactions
//' @param SpecName character vector (NSpec), the name of the chemical species
//'   for which we have formation reactions
//' @param SpecMC IntegerVector (NSpec), the mass compartment of the chemical
//'   species for which we have formation reactions
//' @param SpecK numeric vector (NSpec), the equilibrium coefficient of the
//'   formation reactions.
//' @param SpecTempKelvin NumericVector (NSpec), the temperature associated with
//'   K/logK and DeltaH of the formation reactions 
//' @param SpecDeltaH numeric vector (NSpec), the enthalpy change of the
//'   formation reactions
//' @param SpecStoich signed integer matrix (NSpec x NComp), the reaction
//'   stoichiometry of the formation reactions
//' @param SpecCharge signed integer vector (NSpec), the charge of the chemical
//'   species for which we have formation reactions
//' @param SpecActCorr character vector (NSpec), the activity correction method
//'   of the chemical species for which we have formation reactions
//' @param DoWHAM boolean, true=there are WHAM species, false=no WHAM species
//' @param AqueousMC integer, the (1-based) position of the water/aqueous mass
//'   compartment. (transformed to 0-based at the beginning of the function)
//' @param WHAMDonnanMC the mass compartments corresponding to the
//'   humic acid (0) and fulvic acid (1) Donnan layers. (transformed to 0-based
//'   at the beginning of the function)
//' @param SolHS numeric (2), moles of each organic matter component in solution
//' @param wMolWt numeric (2), WHAM's molecular weight parameter for organic
//'   matter
//' @param wRadius numeric (2), WHAM's molecular radius parameter for organic
//'   matter
//' @param wP numeric (2), WHAM's P parameter...
//' @param wDLF numeric (2), WHAM's Double layer overlap factor
//' @param wKZED numeric (2), WHAM's Constant to control DDL at low ZED
//' @param SysTempKelvin double; input temperature for the current observation,
//'   in Kelvin
//' @param DoTox logical, TRUE for toxicity mode where the MetalName component
//'   concentration is adjusted to try to match the CATarget with BLMetalSpecs
//' @param MetalName character string, the name of the toxic metal
//' @param MetalComp integer, the position of the metal in the component arrays
//'   (i.e., which is the toxic metal component) Note: this are base-1 indexed.
//' @param NBLMetal integer, the number of biotic ligand-bound metal species 
//'   that are associated with toxic effects.
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
                 int MaxIter,
                 bool DoPartialStepsAlways,
                 int NMass,
                 Rcpp::CharacterVector MassName,
                 Rcpp::NumericVector MassAmt,
                 int NComp,
                 Rcpp::CharacterVector CompName,
                 Rcpp::CharacterVector CompType,
                 Rcpp::NumericVector TotConc,
                 int NSpec,
                 Rcpp::CharacterVector SpecName,
                 Rcpp::IntegerVector SpecMC,
                 Rcpp::NumericVector SpecK,
                 Rcpp::NumericVector SpecTempKelvin,
                 Rcpp::NumericVector SpecDeltaH,
                 Rcpp::IntegerMatrix SpecStoich,
                 Rcpp::IntegerVector SpecCharge,
                 Rcpp::CharacterVector SpecActCorr,
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
                 double CATarget) {

  /*outputs*/
  int Iter = 0;  
  double MaxError;
  Rcpp::NumericVector SpecConc(NSpec); // species concentrations after optimization
  Rcpp::NumericVector SpecAct(NSpec);
    SpecAct.names() = SpecName;
  Rcpp::NumericVector CalcTotConc(NComp); // the calculated total concentrations of each component in the simulation (units of e.g., mol/L and mol/kg)}
    CalcTotConc.names() = CompName;
  Rcpp::NumericVector SpecMoles(NSpec);
    SpecMoles.names() = SpecName;

  /*variables*/
  int iComp;
  Rcpp::NumericVector MassAmtAdj(NMass);
    MassAmtAdj.names() = MassName;
  Rcpp::NumericVector CompConcStep(NComp);
  Rcpp::NumericVector CompConc(NComp);
  Rcpp::NumericVector TotMoles(NComp);
    TotMoles.names() = CompName;
  Rcpp::NumericVector CompCtoM(NComp);
  Rcpp::NumericVector CompCtoMAdj(NComp);
  Rcpp::NumericVector CalcTotMoles(NComp);
    CalcTotMoles.names() = CompName;
  Rcpp::IntegerVector CompPosInSpec(NComp);
  Rcpp::NumericVector SpecKTempAdj(NSpec);
  Rcpp::NumericVector SpecKISTempAdj(NSpec);
  Rcpp::NumericVector SpecCtoM(NSpec);
  Rcpp::NumericVector SpecCtoMAdj(NSpec);
  Rcpp::NumericVector SpecActivityCoef(NSpec);
  Rcpp::NumericVector WHAMSpecCharge(2);
  int WhichMax;
  double IonicStrength;
  Rcpp::NumericVector Resid(NComp);
    Resid.names() = CompName;
  Rcpp::NumericVector CompError(NComp);
    CompError.names() = CompName;
  Rcpp::NumericMatrix JacobianMatrix(NComp);

  bool DoPartialSteps;

  double MaxErrorFull;
  Rcpp::NumericVector CompConcStepFull(NComp);
  Rcpp::NumericVector MassAmtAdjFull(NMass);
  Rcpp::NumericVector TotConcFull(NComp);
  Rcpp::NumericVector SpecKISTempAdjFull(NSpec);
  Rcpp::NumericVector SpecCtoMAdjFull(NSpec);
  Rcpp::NumericVector SpecConcFull(NSpec);
  Rcpp::NumericVector CalcTotMolesFull(NComp);
  int WhichMaxFull;
  double IonicStrengthFull;
  Rcpp::NumericVector ResidFull(NComp);
  Rcpp::NumericVector CompErrorFull(NComp);

  double StepSizeAlt;
  double MaxErrorAlt;
  Rcpp::NumericVector CompConcStepAlt(NComp);
  Rcpp::NumericVector MassAmtAdjAlt(NMass);
  Rcpp::NumericVector TotConcAlt(NComp);
  Rcpp::NumericVector SpecKISTempAdjAlt(NSpec);
  Rcpp::NumericVector SpecCtoMAdjAlt(NSpec);
  Rcpp::NumericVector SpecConcAlt(NSpec);
  Rcpp::NumericVector CalcTotMolesAlt(NComp);
  int WhichMaxAlt;
  double IonicStrengthAlt;
  Rcpp::NumericVector ResidAlt(NComp);
  Rcpp::NumericVector CompErrorAlt(NComp);

  /*double StepSizeInterp;
  double MaxErrorInterp;
  Rcpp::NumericVector CompConcStepInterp(NComp);
  Rcpp::NumericVector MassAmtAdjInterp(NMass);
  Rcpp::NumericVector TotConcInterp(NComp);
  Rcpp::NumericVector SpecKISTempAdjInterp(NSpec);
  Rcpp::NumericVector SpecCtoMAdjInterp(NSpec);
  Rcpp::NumericVector SpecConcInterp(NSpec);
  Rcpp::NumericVector CalcTotMolesInterp(NComp);
  int WhichMaxInterp;
  double IonicStrengthInterp;
  Rcpp::NumericVector ResidInterp(NComp);
  Rcpp::NumericVector CompErrorInterp(NComp);*/

  int Counter;
  
  // Initialize some variables
  for (iComp = 0; iComp < NComp; iComp++) {
    CompPosInSpec(iComp) = iComp;
  }
  SpecMC = SpecMC - 1;
  AqueousMC = AqueousMC - 1;
  WHAMDonnanMC = WHAMDonnanMC - 1;
  MassAmtAdj = clone(MassAmt);
  SpecCtoM = MassAmtAdj[SpecMC];
  SpecCtoMAdj = clone(SpecCtoM);
  CompCtoM = SpecCtoM[CompPosInSpec];
  TotMoles = TotConc * CompCtoM;

  SpecKTempAdj = TempCorrection(SysTempKelvin, NSpec, SpecK, SpecTempKelvin, 
                                SpecDeltaH);
  SpecKISTempAdj = clone(SpecKTempAdj);
  /*Rcpp::NumericVector TmpVector = SpecKTempAdj / SpecK;
  Rcpp::Rcout << TmpVector << std::endl;*/

  // Get initial values for component concentrations
  CompConc = InitialGuess(TotConc, SpecCtoMAdj, CompType, SpecKISTempAdj, 
                          SpecStoich, SpecName, NComp, NSpec);
  SpecConc[CompPosInSpec] = clone(CompConc);

  // Calculate the ionic strength and activity coefficients
  IonicStrength = CalcIonicStrength(NSpec, SpecConc * SpecCtoMAdj, SpecCharge, 
                                    SpecMC, AqueousMC);
  SpecActivityCoef = CalcActivityCoef(NSpec, SpecName, SpecActCorr, SpecCharge, 
                                      IonicStrength, SysTempKelvin);

  // Initialize Species Concentrations
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
      Rcpp::Rcout << "Iter=" << Iter 
        << ", WhichMax=" << CompName(WhichMax) 
        << ", MaxError=" << MaxError 
        << ", WHAMSpecCharge=" << WHAMSpecCharge 
        << ", Resid=" << Resid 
        << ", TotMoles=" << TotMoles 
        << ", CalcTotConc=" << CalcTotConc 
        << ", CalcTotMoles=" << CalcTotMoles 
        << ", SpecConc=" << SpecConc 
        << ", SpecMoles=" << SpecMoles
        << ", MassAmtAdj=" << MassAmtAdj
        << std::endl;
    }

    // update the iteration counter
    Iter++;

    // Calculate the Jacobian Matrix
    JacobianMatrix = Jacobian(NComp, NSpec, SpecStoich, SpecConc, SpecCtoMAdj,
                              CompName, MetalComp, NBLMetal, BLMetalSpecs, 
                              DoTox);

    // Calculate the Newton-Raphson step
    CompConcStepFull = CalcStep(JacobianMatrix, Resid, NComp, CompType, CompName);

    // Do a full N-R step
    SpecKISTempAdjFull = clone(SpecKISTempAdj);
    SpecCtoMAdjFull = clone(SpecCtoMAdj);
    SpecConcFull = clone(SpecConc);
    TotConcFull = clone(TotConc);
    MaxErrorFull = CHESSIter(CompConcStepFull, NMass, MassAmt, NComp, CompName, 
                        CompType, CompPosInSpec, NSpec, SpecName, SpecMC, 
                        SpecActCorr, SpecStoich, SpecCharge, SpecKTempAdj, 
                        DoWHAM, AqueousMC, WHAMDonnanMC, SolHS, wMolWt, 
                        wRadius, wP, wDLF, wKZED, SysTempKelvin, DoTox, 
                        MetalName, MetalComp, NBLMetal, BLMetalSpecs, CATarget, 
                        MassAmtAdjFull, TotConcFull, SpecKISTempAdjFull, 
                        SpecCtoMAdjFull, SpecConcFull, CalcTotMolesFull, 
                        WhichMaxFull, IonicStrengthFull, ResidFull, 
                        CompErrorFull);

    if (MaxErrorFull > MaxError) {
      // Do a half-step
      StepSizeAlt = std::pow(0.5, 0.5);//this just ensures it starts at 0.5
      DoPartialSteps = true;
    } else if (MaxErrorFull > (MaxError * 0.99)) {
      // Do a double-step
      StepSizeAlt = std::pow(1.5, 0.5);
      DoPartialSteps = true;
    } else {
      DoPartialSteps = false;
    }
    
    DoPartialSteps = DoPartialSteps && DoPartialStepsAlways;
    //MaxErrorInterp = MaxError + 999;
    MaxErrorAlt = MaxError + 999;

    if (DoPartialSteps) {      
      if (QuietFlag == "Debug") {
        Rcpp::Rcout << 
        "Counter StepSizeAlt MaxErrorAlt"
        << std::endl;
      }

      Counter = 0;
      while((Counter < 3) && (MaxErrorAlt > MaxError)) {
        Counter++;
        StepSizeAlt = std::pow(StepSizeAlt, 2);
        SpecKISTempAdjAlt = clone(SpecKISTempAdj);
        SpecCtoMAdjAlt = clone(SpecCtoMAdj);
        SpecConcAlt = clone(SpecConc);
        TotConcAlt = clone(TotConc);
        CompConcStepAlt = CompConcStepFull * StepSizeAlt;
        MaxErrorAlt = CHESSIter(CompConcStepAlt, NMass, MassAmt, NComp, CompName, 
                            CompType, CompPosInSpec, NSpec, SpecName, SpecMC, 
                            SpecActCorr, SpecStoich, SpecCharge, SpecKTempAdj, 
                            DoWHAM, AqueousMC, WHAMDonnanMC, SolHS, wMolWt, 
                            wRadius, wP, wDLF, wKZED, SysTempKelvin, DoTox, 
                            MetalName, MetalComp, NBLMetal, BLMetalSpecs, CATarget, 
                            MassAmtAdjAlt, TotConcAlt, SpecKISTempAdjAlt, 
                            SpecCtoMAdjAlt, SpecConcAlt, CalcTotMolesAlt, 
                            WhichMaxAlt, IonicStrengthAlt, ResidAlt, 
                            CompErrorAlt);
        if (QuietFlag == "Debug") {
          Rcpp::Rcout << 
          Counter << " " << 
          StepSizeAlt << " " << 
          MaxErrorAlt
          << std::endl;
        }
      }

      /*if (MaxErrorAlt < MaxError) {
       
        // linearly interpolate for a potential best step
        StepSizeInterp = 1 - MaxErrorFull * (1 - StepSizeAlt) /
          (MaxErrorFull - MaxErrorAlt);
        SpecKISTempAdjInterp = clone(SpecKISTempAdj);
        SpecCtoMAdjInterp = clone(SpecCtoMAdj);
        SpecConcInterp = clone(SpecConc);
        TotConcInterp = clone(TotConc);
        CompConcStepInterp = CompConcStepFull * StepSizeInterp;
        MaxErrorInterp = CHESSIter(CompConcStepInterp, NMass, MassAmt, NComp, CompName, 
                            CompType, CompPosInSpec, NSpec, SpecName, SpecMC, 
                            SpecActCorr, SpecStoich, SpecCharge, SpecKTempAdj, 
                            DoWHAM, AqueousMC, WHAMDonnanMC, SolHS, wMolWt, 
                            wRadius, wP, wDLF, wKZED, SysTempKelvin, DoTox, 
                            MetalName, MetalComp, NBLMetal, BLMetalSpecs, CATarget, 
                            MassAmtAdjInterp, TotConcInterp, SpecKISTempAdjInterp, 
                            SpecCtoMAdjInterp, SpecConcInterp, CalcTotMolesInterp, 
                            WhichMaxInterp, IonicStrengthInterp, ResidInterp, 
                            CompErrorInterp);
      }

      if (MaxErrorInterp > MaxErrorAlt){
        // interpolation didn't help us
        MaxErrorInterp = MaxError + 999;
      }*/

    }
    
    /*if (MaxErrorInterp < MaxError) {
      // Store the best partial step in the main variables
      MaxError = MaxErrorInterp;
      WhichMax = WhichMaxInterp;
      IonicStrength = IonicStrengthInterp;
      Resid = clone(ResidInterp);
      CompError = clone(CompErrorInterp);
      MassAmtAdj = clone(MassAmtAdjInterp);
      TotConc = clone(TotConcInterp);
      SpecKISTempAdj = clone(SpecKISTempAdjInterp);
      SpecCtoMAdj = clone(SpecCtoMAdjInterp);
      SpecConc = clone(SpecConcInterp);
      CalcTotMoles = clone(CalcTotMolesInterp);
    } else */if (MaxErrorAlt < MaxError) {
      //Store Alt in the main veriables
      MaxError = MaxErrorAlt;
      WhichMax = WhichMaxAlt;
      IonicStrength = IonicStrengthAlt;
      Resid = clone(ResidAlt);
      CompError = clone(CompErrorAlt);
      MassAmtAdj = clone(MassAmtAdjAlt);
      TotConc = clone(TotConcAlt);
      SpecKISTempAdj = clone(SpecKISTempAdjAlt);
      SpecCtoMAdj = clone(SpecCtoMAdjAlt);
      SpecConc = clone(SpecConcAlt);
      CalcTotMoles = clone(CalcTotMolesAlt);
    } else {
      // Partial Steps has failed, just revert to full
      MaxError = MaxErrorFull;
      WhichMax = WhichMaxFull;
      IonicStrength = IonicStrengthFull;
      Resid = clone(ResidFull);
      CompError = clone(CompErrorFull);
      MassAmtAdj = clone(MassAmtAdjFull);
      TotConc = clone(TotConcFull);
      SpecKISTempAdj = clone(SpecKISTempAdjFull);
      SpecCtoMAdj = clone(SpecCtoMAdjFull);
      SpecConc = clone(SpecConcFull);
      CalcTotMoles = clone(CalcTotMolesFull);
    }
    
    SpecMoles = SpecConc * SpecCtoMAdj;
    CompCtoMAdj = SpecCtoMAdj[CompPosInSpec];
    TotMoles = TotConc * CompCtoMAdj;



    /*
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
    */

  } // while ((MaxError > ConvergenceCriteria) & (Iter <= MaxIter))

  SpecMoles = SpecConc * SpecCtoMAdj;
  SpecAct = SpecConc * SpecActivityCoef;

  if (QuietFlag == "Debug") {
      Rcpp::Rcout << "Iter=" << Iter 
        << ", WhichMax=" << CompName(WhichMax) 
        << ", MaxError=" << MaxError 
        << ", WHAMSpecCharge=" << WHAMSpecCharge 
        << ", Resid=" << Resid 
        << ", TotMoles=" << TotMoles 
        << ", CalcTotConc=" << CalcTotConc 
        << ", CalcTotMoles=" << CalcTotMoles 
        << ", SpecConc=" << SpecConc 
        << ", SpecMoles=" << SpecMoles
        << ", MassAmtAdj=" << MassAmtAdj
        << std::endl;
    }

  return Rcpp::List::create(
      Rcpp::Named("FinalIter") = Iter,
      Rcpp::Named("FinalMaxError") = MaxError,
      Rcpp::Named("SpecConc") = SpecConc,
      Rcpp::Named("SpecAct") = SpecAct,
      Rcpp::Named("SpecMoles") = SpecMoles,
      Rcpp::Named("CalcTotConc") = CalcTotConc,
      Rcpp::Named("MassAmt") = MassAmtAdj);
}
