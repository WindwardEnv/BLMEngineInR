#include <math.h>
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
//' @param SpecType character vector (NSpec), the type or category of the 
//'   chemical species for which we have formation reactions.
//' @param SpecMCR IntegerVector (NSpec), the mass compartment of the chemical
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
//' @param AqueousMCR integer, the (1-based) position of the water/aqueous mass
//'   compartment. (transformed to 0-based at the beginning of the function)
//' @param WHAMDonnanMCR the mass compartments corresponding to the
//'   humic acid (0) and fulvic acid (1) Donnan layers. (transformed to 0-based
//'   at the beginning of the function)
//' @param HumicSubstGramsPerLiter NumericVector, length of 2, grams per liter 
//'   of each organic matter component (HA and FA) in solution
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
//' @param MetalCompR integer, the position of the metal in the component arrays
//'   (i.e., which is the toxic metal component) Note: this is base-1 indexed on
//'   input then converted.
//' @param NBLMetal integer, the number of biotic ligand-bound metal species 
//'   that are associated with toxic effects.
//' @param BLMetalSpecsR integer vector, the positions of the species in the
//'   arrays which contribute to toxicity (i.e., which species are the toxic
//'   metal bound to the relevant biotic ligand) Note: these are base-1 indexed
//'   on input then converted.
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
                 int NMass,
                 Rcpp::CharacterVector MassName,
                 Rcpp::NumericVector MassAmt,
                 int NComp,
                 Rcpp::CharacterVector CompName,
                 Rcpp::CharacterVector CompType,
                 Rcpp::NumericVector TotConc,
                 int NSpec,
                 Rcpp::CharacterVector SpecName,
                 Rcpp::CharacterVector SpecType,
                 Rcpp::IntegerVector SpecMCR,
                 Rcpp::NumericVector SpecK,
                 Rcpp::NumericVector SpecTempKelvin,
                 Rcpp::NumericVector SpecDeltaH,
                 Rcpp::IntegerMatrix SpecStoich,
                 Rcpp::IntegerVector SpecCharge,
                 Rcpp::CharacterVector SpecActCorr,
                 bool DoWHAM,
                 int AqueousMCR,
                 Rcpp::IntegerVector WHAMDonnanMCR,
                 Rcpp::NumericVector HumicSubstGramsPerLiter,
                 Rcpp::NumericVector wMolWt,
                 Rcpp::NumericVector wRadius,
                 Rcpp::NumericVector wP,
                 double wDLF,
                 double wKZED,
                 double SysTempKelvin,
                 bool DoTox,
                 Rcpp::String MetalName,
                 int MetalCompR,
                 int NBLMetal,
                 Rcpp::IntegerVector BLMetalSpecsR,
                 double CATarget) {

  /*outputs*/
  int Iter = 0;  
  double MaxError;
  Rcpp::NumericVector SpecConc(NSpec);
  Rcpp::NumericVector SpecAct(NSpec);
    SpecAct.names() = SpecName;
  Rcpp::NumericVector CalcTotConc(NComp);
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
  Rcpp::NumericVector SpecActivityCoef(NSpec, 1.0);
    //for (int iSpec = 0; iSpec < NSpec; iSpec++) {SpecActivityCoef[iSpec] = 1.0;}
  Rcpp::NumericVector WHAMSpecCharge(2, -0.0001);
    //WHAMSpecCharge[0] = -0.0001;
    //WHAMSpecCharge[1] = -0.0001;
  int WhichMax;
  double IonicStrength = 0.1;
  //double ChargeBalance = 0.0;
  Rcpp::NumericVector Resid(NComp);
    Resid.names() = CompName;
  Rcpp::NumericVector CompError(NComp);
    CompError.names() = CompName;
  Rcpp::NumericMatrix JacobianMatrix(NComp);
  Rcpp::NumericMatrix NumericalJacobianMatrix(NComp);

  Rcpp::String StatusMessage = "";

  // Initialize some variables
  for (iComp = 0; iComp < NComp; iComp++) { CompPosInSpec(iComp) = iComp; }
  Rcpp::IntegerVector SpecMC = clone(SpecMCR) - 1;
  int AqueousMC = AqueousMCR - 1;
  Rcpp::IntegerVector WHAMDonnanMC = clone(WHAMDonnanMCR) - 1;
  int MetalComp = MetalCompR - 1;
  Rcpp::IntegerVector BLMetalSpecs = clone(BLMetalSpecsR) - 1;
  MassAmtAdj = clone(MassAmt);
  SpecCtoM = MassAmtAdj[SpecMC];
  SpecCtoMAdj = clone(SpecCtoM);
  CompCtoM = SpecCtoM[CompPosInSpec];

  bool DebugpH4NoES = false;
  if (DebugpH4NoES) { SysTempKelvin = 288; }
  // Do the temperature adjustments on the binding constants
  SpecKTempAdj = TempCorrection(SysTempKelvin, NSpec, SpecK, SpecTempKelvin, 
                                SpecDeltaH);
  SpecKISTempAdj = clone(SpecKTempAdj);

  // Get initial values for component concentrations
  CompConc = InitialGuess(TotConc, SpecCtoMAdj, CompType, SpecKISTempAdj, 
                          SpecStoich, SpecName, NComp, NSpec, 
                          DoTox, NBLMetal, BLMetalSpecs, MetalComp, CATarget);
  SpecConc[CompPosInSpec] = clone(CompConc);  
  TotMoles = TotConc * CompCtoM;

    /* BEGIN SPECIATION ITERATIONS */

    // Run through CHESS calculations with initial values
    MaxError = CHESSIter(CompConcStep, NMass, MassAmt, NComp, CompName, CompType,
                        CompPosInSpec, NSpec, SpecName, SpecType, SpecMC, SpecActCorr,
                        SpecStoich, SpecCharge, SpecKTempAdj, DoWHAM, false, AqueousMC, 
                        WHAMDonnanMC, HumicSubstGramsPerLiter, wMolWt, wRadius, 
                        wP, wDLF, wKZED, SysTempKelvin, DoTox, MetalName,
                        MetalComp, NBLMetal, BLMetalSpecs, CATarget, MassAmtAdj,
                        TotConc, TotMoles, SpecKISTempAdj, SpecCtoMAdj, SpecConc,
                        SpecActivityCoef, CalcTotMoles, WHAMSpecCharge,
                        WhichMax, IonicStrength, Resid, CompError);

    // Begin iterating
    Iter = 0;
    bool UpdateZED = true;
    while ((Iter == 0) || 
           ((MaxError > ConvergenceCriteria) & (Iter < MaxIter))) {

      /*if (QuietFlag == FLAG_DEBUG) {
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
      }*/

      // update the iteration counter
      Iter++;
      if (Iter <= 6) {
        //Rcpp::NumericVector CompConcStepBrute(NComp);
        CompConcStep = CalcStepBrute(NComp, CompName, CompType, CompConc, 
                                         TotMoles, CalcTotMoles);
      } else {
        //UpdateZED = true;//!UpdateZED;//             
        try {
          // Calculate the Jacobian Matrix
          JacobianMatrix = Jacobian(NComp, NSpec, CompName, TotConc, SpecStoich, SpecConc, 
                                  SpecMC, SpecCtoMAdj, SpecType, SpecCharge, 
                                  SpecKISTempAdj, IonicStrength, DoWHAM, 
                                  HumicSubstGramsPerLiter, WHAMSpecCharge, 
                                  wP, wMolWt, wRadius, wDLF, wKZED, MassAmtAdj, 
                                  AqueousMC, WHAMDonnanMC,  MetalComp, NBLMetal, BLMetalSpecs, 
                                  DoTox);
        
          CompConcStep = CalcStep(JacobianMatrix, Resid, 
                                        CompConc, TotMoles, CalcTotMoles, 
                                        NComp, CompType, CompName);
        }
        catch (int e) {
          if (e == ERROR_JACOBIAN_NAN) {
            StatusMessage += STATUS_JAC_ERR;
            break;
          }          
        }
      }

      
      
      /*double StepBrute;
      double StepNR;
      for (iComp = 0; iComp < NComp; iComp++) {
        StepBrute = CompConcStepBrute(iComp);
        StepNR = CompConcStep(iComp);
        if (StepBrute == StepNR) { continue; }
        if ((std::fabs(StepNR) > CompConc(iComp)) ||
            (std::signbit(StepBrute) != std::signbit(StepNR)) ||
            std::fabs((StepBrute - StepNR) / (StepBrute + StepNR)) > 1.0) {
          //CompConcStep(iComp) = StepBrute;
        }
      }*/

      if (QuietFlag == FLAG_DEBUG) {
        //Rcpp::Rcout << "JacobianMatrix = [" <<std::endl << JacobianMatrix << "]" << std::endl;
        //Rcpp::Rcout << "NumericalJacobianMatrix = [" <<std::endl << NumericalJacobianMatrix << "]" << std::endl;
        Rcpp::Rcout << "iComp\tSpecName\tSpecConc\tResid\tError\tStep" << std::endl;
        for (iComp = 0; iComp < NComp; iComp++) {
          Rcpp::Rcout << iComp << "\t" 
                      << SpecName[iComp] << "\t" 
                      << SpecConc[iComp] << "   " 
                      << Resid[iComp] << "   " 
                      << CompError[iComp] << "   "
                      << CompConcStep[iComp]
                      << std::endl;
          //Rcpp::Rcout << "Resid[" << SpecName[iComp] << "]=" << 
          //  Resid[iComp] << std::endl;
          //Rcpp::Rcout << "CompConcStep[" << SpecName[iComp] << "]=" << 
          //  CompConcStep[iComp] << std::endl;
        }
      }   

      // Do a N-R step
      MaxError = CHESSIter(CompConcStep, NMass, MassAmt, NComp, CompName, 
                          CompType, CompPosInSpec, NSpec, SpecName, SpecType, SpecMC, 
                          SpecActCorr, SpecStoich, SpecCharge, SpecKTempAdj, 
                          DoWHAM, UpdateZED, AqueousMC, WHAMDonnanMC, HumicSubstGramsPerLiter, wMolWt, 
                          wRadius, wP, wDLF, wKZED, SysTempKelvin, DoTox, 
                          MetalName, MetalComp, NBLMetal, BLMetalSpecs, CATarget, 
                          //parameters that are modified and returned:
                          MassAmtAdj, TotConc, TotMoles, SpecKISTempAdj, 
                          SpecCtoMAdj, SpecConc, SpecActivityCoef,
                          CalcTotMoles, WHAMSpecCharge, WhichMax, 
                          IonicStrength, Resid, CompError);

      SpecMoles = SpecConc * SpecCtoMAdj;
      //CompCtoMAdj = SpecCtoMAdj[CompPosInSpec];
      //TotMoles = TotConc * CompCtoMAdj;
      CalcTotConc = CalcTotMoles / CompCtoM;//Adj;

    } // while ((MaxError > ConvergenceCriteria) & (Iter <= MaxIter))

    /* END SPECIATION ITERATIONS */

  }

  if ((MaxError > ConvergenceCriteria)) { 
    StatusMessage += STATUS_SPEC_DNC;    
  }

  SpecMoles = SpecConc * SpecCtoMAdj;
  SpecAct = SpecConc * SpecActivityCoef;

  if (QuietFlag == FLAG_DEBUG) {
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
      Rcpp::Named("StatusMessage") = StatusMessage,
      Rcpp::Named("FinalIter") = Iter,
      Rcpp::Named("FinalToxIter") = ToxIter,
      Rcpp::Named("FinalMaxError") = MaxError,
      Rcpp::Named("IonicStrength") = IonicStrength,
      Rcpp::Named("SpecConc") = SpecConc,
      Rcpp::Named("SpecAct") = SpecAct,
      Rcpp::Named("SpecMoles") = SpecMoles,
      Rcpp::Named("CalcTotConc") = CalcTotConc,
      Rcpp::Named("CalcTotMoles") = CalcTotMoles,
      Rcpp::Named("WHAMSpecCharge") = WHAMSpecCharge,
      Rcpp::Named("MassAmt") = MassAmtAdj);
}
