// Copyright 2024 Windward Environmental LLC
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <math.h>
#include <iostream>
#include <fstream>
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
//' @param WHAMMolWt numeric (2), WHAM's molecular weight parameter for organic
//'   matter
//' @param WHAMRadius numeric (2), WHAM's molecular radius parameter for organic
//'   matter
//' @param WHAMP numeric (2), WHAM's P parameter...
//' @param WHAMDLF numeric (2), WHAM's Double layer overlap factor
//' @param WHAMKZED numeric (2), WHAM's Constant to control DDL at low ZED
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
                 Rcpp::NumericVector WHAMMolWt,
                 Rcpp::NumericVector WHAMRadius,
                 Rcpp::NumericVector WHAMP,
                 double WHAMDLF,
                 double WHAMKZED,
                 double SysTempKelvin,
                 bool DoTox,
                 Rcpp::String MetalName,
                 int MetalCompR,
                 int BLCompR,
                 int NBLMetal,
                 Rcpp::IntegerVector BLMetalSpecsR,
                 double CATarget,
                 bool DodVidCj,
                 bool DodVidCjDonnan,
                 bool DodKidCj,
                 bool DoGammai,
                 bool DoJacDonnan,
                 bool DoJacWHAM,
                 bool DoWHAMSimpleAdjust,
                 bool DoDonnanSimpleAdjust) {

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
  Rcpp::NumericVector CompConcStep(NComp, 0.0);
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
  double WHAMIonicStrength = 0.1;
  double ChargeBalance = 0.0;
  Rcpp::NumericVector Resid(NComp);
    Resid.names() = CompName;
  Rcpp::NumericVector CompError(NComp);
    CompError.names() = CompName;
  Rcpp::NumericMatrix JacobianMatrix(NComp);

  Rcpp::String StatusMessage = "";

  // Don't Be Hasty Variables
  bool Skip = false;
  int BadCount = 0;
  int DontBeHasty = 0;
  int GoodCount = 0;
  int NothingCount = 0;
  int SkipCount = 0;
  double Relax = 1.0;
  bool RelaxOn = false;
  double OldMaxError;
  Rcpp::NumericVector OldCompConc(NComp, 0.0);
  int SlowCount = 0;

  // Initialize some variables
  for (iComp = 0; iComp < NComp; iComp++) { CompPosInSpec(iComp) = iComp; }
  Rcpp::IntegerVector SpecMC = clone(SpecMCR) - 1;
  int AqueousMC = AqueousMCR - 1;
  Rcpp::IntegerVector WHAMDonnanMC = clone(WHAMDonnanMCR) - 1;
  int MetalComp = MetalCompR - 1;
  int BLComp = BLCompR - 1;
  Rcpp::IntegerVector BLMetalSpecs = clone(BLMetalSpecsR) - 1;
  MassAmtAdj = clone(MassAmt);
  SpecCtoM = MassAmtAdj[SpecMC];
  SpecCtoMAdj = clone(SpecCtoM);
  CompCtoM = SpecCtoM[CompPosInSpec];

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
                      WHAMDonnanMC, HumicSubstGramsPerLiter, WHAMMolWt, WHAMRadius, 
                      WHAMP, WHAMDLF, WHAMKZED, SysTempKelvin, DoTox, MetalName,
                      MetalComp, NBLMetal, BLMetalSpecs, CATarget, MassAmtAdj,
                      TotConc, TotMoles, SpecKISTempAdj, SpecCtoMAdj, SpecConc,
                      SpecActivityCoef, CalcTotMoles, WHAMSpecCharge,
                      WhichMax, IonicStrength, WHAMIonicStrength, ChargeBalance, Resid, CompError,
                      DoWHAMSimpleAdjust, DoDonnanSimpleAdjust,
                      ConvergenceCriteria, MaxIter);
  CompConc = SpecConc[CompPosInSpec];
  
  if (QuietFlag == FLAG_DEBUG) {
    Rcpp::Rcout << "Iter\tMaxError\tWhichMax" << std::endl;
  }

  // Begin iterating
  Iter = 0;
  bool UpdateZED = true;
  OldMaxError = MaxError; // needed for Don't Be Hasty Routine
  while ((Iter == 0) || 
          ((MaxError > ConvergenceCriteria) & (Iter < MaxIter))) {

    if (QuietFlag == FLAG_DEBUG) {
      Rcpp::Rcout << Iter << "\t" << MaxError << "\t" << CompName(WhichMax) << std::endl;
    }

    // update the iteration counter
    Iter++;
      //UpdateZED = true;//!UpdateZED;//             
      try {
        // Calculate the Jacobian Matrix
        JacobianMatrix = Jacobian(NComp, NSpec, CompType, SpecStoich, SpecConc, 
                                SpecMC, SpecCtoMAdj, SpecType, SpecCharge, 
                                SpecKISTempAdj, SpecActivityCoef, WHAMIonicStrength, DoWHAM, 
                                HumicSubstGramsPerLiter, WHAMSpecCharge, 
                                WHAMP, WHAMMolWt, WHAMRadius, WHAMDLF, WHAMKZED, 
                                AqueousMC, MetalComp, BLComp, NBLMetal, BLMetalSpecs, 
                                DoTox, 
                                DodVidCj, DodVidCjDonnan, DodKidCj, DoGammai, 
                                DoJacDonnan, DoJacWHAM);

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
    //}

      // Don't Be Hasty Routine Start
      if (Iter > 20) {
        RelaxOn = true;
      }
      if (Skip) { 
        Skip = false;
        continue;
      }
      if (RelaxOn) {
        if (SlowCount > 2) {
          SimpleAdjustComp(WhichMax, ConvergenceCriteria, MaxIter, TotMoles(WhichMax), 
                NComp, CompConc,
                NSpec, SpecConc, SpecKISTempAdj, SpecStoich, SpecName,
                SpecType, SpecActivityCoef, SpecCtoMAdj, SpecCharge,
                WHAMSpecCharge);
          CompConcStep(WhichMax) = 0.0;
          //SlowCount = 0;
          GoodCount = 0;
          NothingCount = 0;
          BadCount = 0;
        } else if ((MaxError > (1.3 * OldMaxError)) || (NothingCount > 10)) {
          GoodCount = 0;
          NothingCount = 0;
          SlowCount = 0;
          BadCount++;
          if (BadCount > 5) {
            SkipCount++;
            if (SkipCount > 5) { BadCount = 1;}
          } else {
            Relax = Relax / 1.5;
            if (Relax == 0) { Relax = 0.05; }
            Skip = true;
          }
        } else if (MaxError < OldMaxError) { 
          BadCount = 0;
          GoodCount++;
          if ((GoodCount > (3 + Iter / 3)) &&
              (GoodCount > DontBeHasty) && 
              (Relax < 1)) { 
            NothingCount = 0;
            DontBeHasty = GoodCount / 5;
            Relax = Relax * 1.5;
            if (Relax > 1) { 
              Relax = 1.0; 
            }
          }
          //if (MaxError > (0.95 * OldMaxError)) {
          if ((MaxError - (Iter - MaxIter) * (MaxError - OldMaxError)) > ConvergenceCriteria) {
            SlowCount++;
          } else {
            SlowCount = 0;
          }
        } else {
          SlowCount++;
          NothingCount++;
        }
      } else {
        Relax = 1.0;
      }
      if (!Skip) {
        OldCompConc = clone(CompConc);          
      } else {
        CompConc = clone(OldCompConc);
      }
      for (iComp = 0; iComp < NComp; iComp++){
        CompConcStep(iComp) = CompConcStep(iComp) * Relax;   
      }  
      OldMaxError = MaxError;
      // Don't Be Hasty Routine End    

    if (false) {//(Iter >= 6) {
      double BestMult = 1.0;
      double BestMaxError = MaxError;
      Rcpp::NumericVector CompConcStepSave(NComp);
      Rcpp::NumericVector CompConcSave(NComp);
      Rcpp::NumericVector MassAmtAdjSave(NMass);
      Rcpp::NumericVector TotConcSave(NComp);
      Rcpp::NumericVector TotMolesSave(NComp);
      Rcpp::NumericVector SpecKISTempAdjSave(NSpec);
      Rcpp::NumericVector SpecCtoMAdjSave(NSpec);
      Rcpp::NumericVector SpecConcSave(NSpec);
      Rcpp::NumericVector SpecActivityCoefSave(NSpec);
      Rcpp::NumericVector CalcTotMolesSave(NComp);
      Rcpp::NumericVector WHAMSpecChargeSave(2);
      int WhichMaxSave;
      double IonicStrengthSave;
      double WHAMIonicStrengthSave;
      double ChargeBalanceSave; 
      Rcpp::NumericVector ResidSave(NComp);
      Rcpp::NumericVector CompErrorSave(NComp);
      CompConcStepSave = clone(CompConcStep);
      CompConcSave = clone(CompConc);
      MassAmtAdjSave = clone(MassAmtAdj);
      TotConcSave = clone(TotConc);
      TotMolesSave = clone(TotMoles);
      SpecKISTempAdjSave = clone(SpecKISTempAdj);
      SpecCtoMAdjSave = clone(SpecCtoMAdj);
      SpecConcSave = clone(SpecConc);
      SpecActivityCoefSave = clone(SpecActivityCoef);
      CalcTotMolesSave = clone(CalcTotMoles);
      WHAMSpecChargeSave = clone(WHAMSpecCharge);
      WhichMaxSave = WhichMax;
      IonicStrengthSave = IonicStrength;
      WHAMIonicStrengthSave = WHAMIonicStrength;
      ChargeBalanceSave = ChargeBalance;
      ResidSave = clone(Resid);
      CompErrorSave = clone(CompError);

      std::ofstream myfile;
      myfile.open("C:/Users/kellyc/Documents/BLM Development/engine/BLMEngineInR_scrap/CompError and spec Output.txt");

      myfile << "StepMult";
      for (iComp = 0; iComp < NComp; iComp++){
        myfile << "\t" << CompName[iComp];
      }
      for (int iSpec = 0; iSpec < NSpec; iSpec++) {
        myfile << "\t" << "[" << SpecName[iSpec] << "]";
      }
      myfile << std::endl;

      double StepMult = 0.01;
      while (StepMult <= 100) {     
        for (iComp = 0; iComp < NComp; iComp++){
          CompConcStep(iComp) = CompConcStepSave(iComp) * StepMult;   
        }      
        MaxError = CHESSIter(CompConcStep, NMass, MassAmt, NComp, CompName, 
                         CompType, CompPosInSpec, NSpec, SpecName, SpecType, SpecMC, 
                         SpecActCorr, SpecStoich, SpecCharge, SpecKTempAdj, 
                         DoWHAM, UpdateZED, AqueousMC, WHAMDonnanMC, HumicSubstGramsPerLiter, WHAMMolWt, 
                         WHAMRadius, WHAMP, WHAMDLF, WHAMKZED, SysTempKelvin, DoTox, 
                         MetalName, MetalComp, NBLMetal, BLMetalSpecs, CATarget, 
                         //parameters that are modified and returned:
                         MassAmtAdj, TotConc, TotMoles, SpecKISTempAdj, 
                         SpecCtoMAdj, SpecConc, SpecActivityCoef,
                         CalcTotMoles, WHAMSpecCharge, WhichMax, 
                         IonicStrength, WHAMIonicStrength, ChargeBalance, Resid, CompError,
                         DoWHAMSimpleAdjust, DoDonnanSimpleAdjust,
                         ConvergenceCriteria, MaxIter);
        myfile << StepMult;
        for (iComp = 0; iComp < NComp; iComp++) {
          myfile << "\t" << CompError[iComp];
        }
        for (int iSpec = 0; iSpec < NSpec; iSpec++) {
          myfile << "\t" << SpecConc[iSpec];
        }
        myfile << std::endl;   

        if (MaxError < BestMaxError) {
          BestMult = StepMult;
          BestMaxError = MaxError;
        }

        CompConcStep = clone(CompConcStepSave);
        CompConc = clone(CompConcSave);
        MassAmtAdj = clone(MassAmtAdjSave);
        TotConc = clone(TotConcSave);
        TotMoles = clone(TotMolesSave);
        SpecKISTempAdj = clone(SpecKISTempAdjSave);
        SpecCtoMAdj = clone(SpecCtoMAdjSave);
        SpecConc = clone(SpecConcSave);
        SpecActivityCoef = clone(SpecActivityCoefSave);
        CalcTotMoles = clone(CalcTotMolesSave);
        WHAMSpecCharge = clone(WHAMSpecChargeSave);
        WhichMax = WhichMaxSave;
        IonicStrength = IonicStrengthSave;
        WHAMIonicStrength = WHAMIonicStrengthSave;
        ChargeBalance = ChargeBalanceSave;
        Resid = clone(ResidSave);
        CompError = clone(CompErrorSave);
        StepMult *= pow(10, (std::log10(100) - std::log10(0.01)) / 22);

      }

      myfile.close();
      //StepMult = 1.0;
      for (iComp = 0; iComp < NComp; iComp++){
        CompConcStep(iComp) = CompConcStepSave(iComp) * BestMult;   
      }
      
    }
    
    // Do a N-R step
    MaxError = CHESSIter(CompConcStep, NMass, MassAmt, NComp, CompName, 
                         CompType, CompPosInSpec, NSpec, SpecName, SpecType, SpecMC, 
                         SpecActCorr, SpecStoich, SpecCharge, SpecKTempAdj, 
                         DoWHAM, UpdateZED, AqueousMC, WHAMDonnanMC, HumicSubstGramsPerLiter, WHAMMolWt, 
                         WHAMRadius, WHAMP, WHAMDLF, WHAMKZED, SysTempKelvin, DoTox, 
                         MetalName, MetalComp, NBLMetal, BLMetalSpecs, CATarget, 
                         //parameters that are modified and returned:
                         MassAmtAdj, TotConc, TotMoles, SpecKISTempAdj, 
                         SpecCtoMAdj, SpecConc, SpecActivityCoef,
                         CalcTotMoles, WHAMSpecCharge, WhichMax, 
                         IonicStrength, WHAMIonicStrength, ChargeBalance, Resid, CompError,
                         DoWHAMSimpleAdjust, DoDonnanSimpleAdjust,
                         ConvergenceCriteria, MaxIter);
    CompConc = SpecConc[CompPosInSpec];
    SpecMoles = SpecConc * SpecCtoMAdj;
    //CompCtoMAdj = SpecCtoMAdj[CompPosInSpec];
    //TotMoles = TotConc * CompCtoMAdj;
    CalcTotConc = CalcTotMoles / CompCtoM;//Adj;

  } // while ((MaxError > ConvergenceCriteria) & (Iter <= MaxIter))

  /* END SPECIATION ITERATIONS */


  if ((MaxError > ConvergenceCriteria)) { 
    StatusMessage += STATUS_SPEC_DNC;    
  }

  SpecMoles = SpecConc * SpecCtoMAdj;
  SpecAct = SpecConc * SpecActivityCoef;

  if (QuietFlag == FLAG_DEBUG) {
    Rcpp::Rcout << Iter << "\t" << MaxError << "\t" << CompName(WhichMax) << std::endl;
  }

  return Rcpp::List::create(
      Rcpp::Named("StatusMessage") = StatusMessage,
      Rcpp::Named("FinalIter") = Iter,
      Rcpp::Named("FinalMaxError") = MaxError,
      Rcpp::Named("IonicStrength") = IonicStrength,
      Rcpp::Named("WHAMIonicStrength") = WHAMIonicStrength,
      Rcpp::Named("ChargeBalance") = ChargeBalance,
      Rcpp::Named("SpecConc") = SpecConc,
      Rcpp::Named("SpecAct") = SpecAct,
      Rcpp::Named("SpecMoles") = SpecMoles,
      Rcpp::Named("CalcTotConc") = CalcTotConc,
      Rcpp::Named("CalcTotMoles") = CalcTotMoles,
      Rcpp::Named("WHAMSpecCharge") = WHAMSpecCharge,
      Rcpp::Named("MassAmt") = MassAmtAdj);
}
