#include <Rcpp.h>
#include "CHESSFunctions.h"

//' @title AdjustForWHAM
//' 
//' @description Make any necessary adjustments to the simulation that involve
//'   WHAM organic matter.
//' 
//' @details 
//'   This function does several tasks necessary when simulating WHAM
//'   organic matter. It will:
//'     1. Calculate the charge on the WHAM organic matter species,
//'     2. Adjust the organic matter binding constants for ionic strength
//'        effects,
//'     3. Calculate the volume of the Donnan layer(s) (i.e., diffuse binding
//'        layers),
//'     4. Make the totals of the Donnan components equal to the net humic
//'        charge to be balanced (so that Newton-Rapshon will work on it), and
//'     5. Set the species concentrations of ions with the same charge as the
//'        net humic charge that are binding diffusely equal to zero (only 
//'        counterions should bind diffusely).
//' 
//' @author Kelly Croteau (kellyc@windwardenv.com)
//' 
//' @param NMass int, the number of mass compartments
//' @param MassAmt NumericVector, The amount of each mass compartment.
//' @param MassAmtAdj (OUTPUT), NumericVector, The amount of each mass 
//'   compartment, adjusted for the Donnan Layer volumes (aqueous compartment 
//'   will be its original value minus the total volume of the Donnan Layers).
//' @param NComp int, the total number of components in the simulation,
//'   including the input components, defined components (and including the
//'   defined components that get added by ExpandWHAM)
//' @param CompType CharacterVector, length of NComp, the type of each component
//'   in the simulation
//' @param TotConc (INPUT & OUTPUT) NumericVector, length of NComp, the total
//'   concentrations of each component in the simulation (units of e.g., mol/L
//'   and mol/kg)
//' @param TotMoles (INPUT & OUTPUT) NumericVector, the total moles of each
//'   component in the simulation (units of mol)
//' @param NSpec int, the number of chemical species for which we have formation
//'   reactions in the simulation
//' @param SpecConc (INPUT & OUTPUT) NumericVector, length NSpec, the
//'   concentrations of each species for which we have formation reactions
//' @param SpecMC IntegerVector, length NSpec, the mass compartment of the
//'   chemical species for which we have formation reactions
//' @param SpecActCorr CharacterVector, length NSpec, the activity correction
//'   method of the chemical species for which we have formation reactions
//' @param SpecCharge IntegerVector, length NSpec, the charge of the chemical
//'   species for which we have formation reactions
//' @param SpecKTempAdj NumericVector, length NSpec, the equilibrium coefficient
//'   of the formation reactions, adjusted for temperature only
//' @param SpecKISTempAdj (OUTPUT) NumericVector, length NSpec, the equilibrium
//'   coefficient of the formation reactions, adjusted for temperature and ionic
//'   strength effects
//' @param SpecCtoMAdj (INPUT & OUTPUT) NumericVector, length NSpec, the
//'   concentration to mass conversion factor of the chemical species for which
//'   we have formation reactions, adjusted to take into account diffuse binding
//'   to organic matter.
//' @param IonicStrength double, the ionic strength of the solution
//' @param WHAMSpecCharge (OUTPUT) NumericVector, length 2, the net humic charge
//'   for humic (0) and fulvic (1) substances.
//' @param AqueousMC int, the mass compartment corresponding to water/aqueous
//' @param WHAMDonnanMC IntegerVector, length 2, the mass compartments
//'   corresponding to the humic acid (0) and fulvic acid (1) Donnan layers.
//' @param HumicSubstGramsPerLiter NumericVector, length of 2, grams per liter 
//'   of each organic matter component (HA and FA) in solution
//' @param wMolWt NumericVector, WHAM's molecular weight parameter for organic
//'   matter
//' @param wRadius NumericVector, WHAM's molecular radius parameter for organic
//'   matter
//' @param wP NumericVector, WHAM's P parameter...
//' @param wDLF double, WHAM's Double layer overlap factor
//' @param wKZED double, WHAM's Constant to control DDL at low ZED
//' 
void AdjustForWHAM(
  int NMass,   
  Rcpp::NumericVector MassAmt,
  Rcpp::NumericVector &MassAmtAdj,
  int NComp,
  Rcpp::CharacterVector CompType,
  Rcpp::NumericVector &TotConc,
  Rcpp::NumericVector &TotMoles,
  int NSpec,
  Rcpp::NumericVector &SpecConc,
  Rcpp::IntegerVector SpecMC,
  Rcpp::CharacterVector SpecActCorr,
  Rcpp::IntegerVector SpecCharge,
  Rcpp::NumericVector SpecKTempAdj,
  Rcpp::NumericVector &SpecKISTempAdj,
  Rcpp::NumericVector &SpecCtoMAdj,  
  double IonicStrength,
  Rcpp::NumericVector &WHAMSpecCharge,
  int AqueousMC,
  Rcpp::IntegerVector WHAMDonnanMC,
  Rcpp::NumericVector HumicSubstGramsPerLiter,
  Rcpp::NumericVector wMolWt,
  Rcpp::NumericVector wRadius,
  Rcpp::NumericVector wP,
  double wDLF,
  double wKZED
) {

  /* variables */
  int iComp;
  int iSpec;

  //Calculate the charge on the organic matter
  WHAMSpecCharge = CalcWHAMSpecCharge(NSpec, SpecActCorr, SpecConc, 
                                      SpecCharge, SpecMC, AqueousMC, 
                                      HumicSubstGramsPerLiter);
  WHAMSpecCharge = WHAMSpecCharge / HumicSubstGramsPerLiter;

  //Adjust organic matter specific binding based on ionic strength
  SpecKISTempAdj = CalcIonicStrengthEffects(IonicStrength, WHAMSpecCharge, 
                                            NSpec, SpecCharge, SpecKTempAdj, 
                                            SpecActCorr, wP);

  //Calculate the portion of the solution that's in the diffuse layer
  MassAmtAdj = CalcDonnanLayerVolume(NMass, NSpec, IonicStrength, MassAmt, 
                                     AqueousMC, WHAMDonnanMC, wMolWt, 
                                     wRadius, wDLF, wKZED, WHAMSpecCharge, 
                                     HumicSubstGramsPerLiter);
  SpecCtoMAdj = MassAmtAdj[SpecMC];

  // Set -Z_HS to the "known" total for the Donnan component
  for (iComp = 0; iComp < NComp; iComp++) {
    if (CompType(iComp) == TYPE_DONNANHA) {
      TotMoles(iComp) = abs(WHAMSpecCharge(iHA)) * HumicSubstGramsPerLiter(iHA);
      TotConc(iComp) = TotMoles(iComp) / SpecCtoMAdj(iComp);
    } else if (CompType(iComp) == TYPE_DONNANFA) {
      TotMoles(iComp) = abs(WHAMSpecCharge(iFA)) * HumicSubstGramsPerLiter(iFA);
      TotConc(iComp) = TotMoles(iComp) / SpecCtoMAdj(iComp);
    }
  }

  /* The cation species concentrations should be canceled out if Z_HS is 
     positive, and anion species concentrations should be canceled out if Z_HS
     is negative */
  for (iSpec = NComp; iSpec < NSpec; iSpec++) {
    if ((SpecActCorr(iSpec) == ACTYPE_DONNANHA) && 
        (((WHAMSpecCharge(iHA) < 0) && (SpecCharge(iSpec) < 0)) || 
         ((WHAMSpecCharge(iHA) > 0) && (SpecCharge(iSpec) > 0)))) {
      SpecConc(iSpec) = 0;
    } else if ((SpecActCorr(iSpec) == ACTYPE_DONNANFA) &&
               (((WHAMSpecCharge(iFA) < 0) && (SpecCharge(iSpec) < 0)) || 
                ((WHAMSpecCharge(iFA) > 0) && (SpecCharge(iSpec) > 0)))) {
      SpecConc(iSpec) = 0;
    }
  }

}

//' @title AdjustForWHAMBeforeCalcSpecies
//' 
//' @description Make any necessary adjustments to the simulation that involve
//'   WHAM organic matter. The adjustments in this functionhappen before 
//'   CalcSpecConc is called
//' 
//' @details 
//'   This function does several tasks necessary when simulating WHAM
//'   organic matter. It will:
//'     1. Adjust the organic matter binding constants for ionic strength
//'        effects,
//'     2. Calculate the volume of the Donnan layer(s) (i.e., diffuse binding
//'        layers),
//' 
//' @author Kelly Croteau (kellyc@windwardenv.com)
//' 
//' @param NMass int, the number of mass compartments
//' @param MassAmt NumericVector, The amount of each mass compartment.
//' @param MassAmtAdj (OUTPUT), NumericVector, The amount of each mass 
//'   compartment, adjusted for the Donnan Layer volumes (aqueous compartment 
//'   will be its original value minus the total volume of the Donnan Layers).
//' @param NSpec int, the number of chemical species for which we have formation
//'   reactions in the simulation
//' @param SpecMC IntegerVector, length NSpec, the mass compartment of the
//'   chemical species for which we have formation reactions
//' @param SpecActCorr CharacterVector, length NSpec, the activity correction
//'   method of the chemical species for which we have formation reactions
//' @param SpecCharge IntegerVector, length NSpec, the charge of the chemical
//'   species for which we have formation reactions
//' @param SpecKTempAdj NumericVector, length NSpec, the equilibrium coefficient
//'   of the formation reactions, adjusted for temperature only
//' @param SpecKISTempAdj (OUTPUT) NumericVector, length NSpec, the equilibrium
//'   coefficient of the formation reactions, adjusted for temperature and ionic
//'   strength effects
//' @param SpecCtoMAdj (OUTPUT) NumericVector, length NSpec, the
//'   concentration to mass conversion factor of the chemical species for which
//'   we have formation reactions, adjusted to take into account diffuse binding
//'   to organic matter.
//' @param IonicStrength double, the ionic strength of the solution
//' @param WHAMSpecCharge NumericVector, length 2, the net humic charge for
//'   humic (0) and fulvic (1) substances.
//' @param AqueousMC int, the mass compartment corresponding to water/aqueous
//' @param WHAMDonnanMC IntegerVector, length 2, the mass compartments
//'   corresponding to the humic acid (0) and fulvic acid (1) Donnan layers.
//' @param HumicSubstGramsPerLiter NumericVector, length of 2, grams per liter 
//'   of each organic matter component (HA and FA) in solution
//' @param wMolWt NumericVector, WHAM's molecular weight parameter for organic
//'   matter
//' @param wRadius NumericVector, WHAM's molecular radius parameter for organic
//'   matter
//' @param wP NumericVector, WHAM's P parameter...
//' @param wDLF double, WHAM's Double layer overlap factor
//' @param wKZED double, WHAM's Constant to control DDL at low ZED
//' 
void AdjustForWHAMBeforeCalcSpecies(
  int NMass,   
  Rcpp::NumericVector MassAmt,
  Rcpp::NumericVector &MassAmtAdj,
  int NSpec,
  Rcpp::IntegerVector SpecMC,
  Rcpp::CharacterVector SpecActCorr,
  Rcpp::IntegerVector SpecCharge,
  Rcpp::NumericVector SpecKTempAdj,
  Rcpp::NumericVector &SpecKISTempAdj,
  Rcpp::NumericVector &SpecCtoMAdj,  
  double IonicStrength,
  Rcpp::NumericVector WHAMSpecCharge,
  int AqueousMC,
  Rcpp::IntegerVector WHAMDonnanMC,
  Rcpp::NumericVector HumicSubstGramsPerLiter,
  Rcpp::NumericVector wMolWt,
  Rcpp::NumericVector wRadius,
  Rcpp::NumericVector wP,
  double wDLF,
  double wKZED
) {

  int iSpec;

  //Adjust organic matter specific binding based on ionic strength
  SpecKISTempAdj = CalcIonicStrengthEffects(IonicStrength, WHAMSpecCharge, 
                                            NSpec, SpecCharge, SpecKTempAdj, 
                                            SpecActCorr, wP);

  //Calculate the portion of the solution that's in the diffuse layer
  MassAmtAdj = CalcDonnanLayerVolume(NMass, NSpec, IonicStrength, MassAmt, 
                                     AqueousMC, WHAMDonnanMC, wMolWt, 
                                     wRadius, wDLF, wKZED, WHAMSpecCharge, 
                                     HumicSubstGramsPerLiter);
  //SpecCtoMAdj = MassAmtAdj[SpecMC];
  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    if ((SpecActCorr(iSpec) == ACTYPE_WHAMHA) || (SpecActCorr(iSpec) == ACTYPE_WHAMFA)) {
      SpecCtoMAdj[iSpec] = MassAmt[AqueousMC];
    } else {
      SpecCtoMAdj[iSpec] = MassAmtAdj[SpecMC[iSpec]];
    }
  }

}

//' @title AdjustForWHAMAfterCalcSpecies
//' 
//' @description Make any necessary adjustments to the simulation that involve
//'   WHAM organic matter. The adjustments in this function should happen after
//'   CalcSpecConc is called.
//' 
//' @details 
//'   This function does several tasks necessary when simulating WHAM
//'   organic matter. It will:
//'     1. Set the species concentrations of ions with the same charge as the
//'        net humic charge that are binding diffusely equal to zero (only 
//'        counterions should bind diffusely).
//'     2. Calculate the charge on the WHAM organic matter species,
//'     3. Make the totals of the Donnan components equal to the net humic
//'        charge to be balanced (so that Newton-Rapshon will work on it), and
//' 
//' @author Kelly Croteau (kellyc@windwardenv.com)
//' 
//' @param NComp int, the total number of components in the simulation,
//'   including the input components, defined components (and including the
//'   defined components that get added by ExpandWHAM)
//' @param CompType CharacterVector, length of NComp, the type of each component
//'   in the simulation
//' @param TotConc (INPUT & OUTPUT) NumericVector, length of NComp, the total
//'   concentrations of each component in the simulation (units of e.g., mol/L
//'   and mol/kg)
//' @param TotMoles (INPUT & OUTPUT) NumericVector, the total moles of each
//'   component in the simulation (units of mol)
//' @param NSpec int, the number of chemical species for which we have formation
//'   reactions in the simulation
//' @param SpecConc (INPUT & OUTPUT) NumericVector, length NSpec, the
//'   concentrations of each species for which we have formation reactions
//' @param SpecMC IntegerVector, length NSpec, the mass compartment of the
//'   chemical species for which we have formation reactions
//' @param SpecActCorr CharacterVector, length NSpec, the activity correction
//'   method of the chemical species for which we have formation reactions
//' @param SpecCharge IntegerVector, length NSpec, the charge of the chemical
//'   species for which we have formation reactions
//' @param SpecCtoMAdj NumericVector, length NSpec, the concentration to mass
//'   conversion factor of the chemical species for which we have formation
//'   reactions, adjusted to take into account diffuse binding to organic
//'   matter.
//' @param WHAMSpecCharge (OUTPUT) NumericVector, length 2, the net humic charge
//'   for humic (0) and fulvic (1) substances.
//' @param AqueousMC int, the mass compartment corresponding to water/aqueous
//' @param HumicSubstGramsPerLiter NumericVector, length of 2, grams per liter 
//'   of each organic matter component (HA and FA) in solution
//' 
void AdjustForWHAMAfterCalcSpecies(int NComp,
                                   Rcpp::CharacterVector CompType,
                                   Rcpp::NumericVector &TotConc,
                                   Rcpp::NumericVector &TotMoles,
                                   int NSpec,
                                   Rcpp::CharacterVector SpecName,
                                   Rcpp::NumericVector &SpecConc,
                                   Rcpp::NumericVector SpecKISTempAdj,
                                   Rcpp::IntegerMatrix SpecStoich,
                                   Rcpp::NumericVector SpecActivityCoef,
                                   Rcpp::IntegerVector SpecMC,
                                   Rcpp::CharacterVector SpecActCorr,
                                   Rcpp::IntegerVector SpecCharge,
                                   Rcpp::NumericVector SpecCtoMAdj,  
                                   Rcpp::NumericVector &WHAMSpecCharge,
                                   int AqueousMC,
                                   Rcpp::NumericVector HumicSubstGramsPerLiter,
                                   bool UpdateZED) {

  /* variables */
  const double ConvCrit = 0.000001;
  const int MaxIter = 10;
  int iComp;
  Rcpp::NumericVector NewWHAMSpecCharge(2);
  Rcpp::NumericVector CompConc(NComp);
  for (iComp = 0; iComp < NComp; iComp++) { CompConc(iComp) = SpecConc(iComp); }
  
  // Update the WHAM component concentrations
  for (iComp = 0; iComp < NComp; iComp++) {
    if ((SpecActCorr(iComp) == ACTYPE_WHAMHA) || (SpecActCorr(iComp) == ACTYPE_WHAMFA)) {

      SimpleAdjustComp(iComp, ConvCrit, MaxIter, TotMoles(iComp), 
                       NComp, CompConc,
                       NSpec, SpecConc, SpecKISTempAdj, SpecStoich, SpecName,
                       SpecActCorr, SpecActivityCoef, SpecCtoMAdj, SpecCharge,
                       WHAMSpecCharge);

    }
  }

  //Calculate the charge on the organic matter
  if (UpdateZED) {
    NewWHAMSpecCharge = CalcWHAMSpecCharge(NSpec, SpecActCorr, SpecConc, 
                                           SpecCharge, SpecMC, AqueousMC, 
                                           HumicSubstGramsPerLiter);
    WHAMSpecCharge(iHA) = WHAMSpecCharge(iHA) + ((NewWHAMSpecCharge(iHA) - WHAMSpecCharge(iHA)) / 5);
    WHAMSpecCharge(iFA) = WHAMSpecCharge(iFA) + ((NewWHAMSpecCharge(iFA) - WHAMSpecCharge(iFA)) / 5);
  }

  // Set -Z_HS*Th to the "known" total for the Donnan component
  for (iComp = 0; iComp < NComp; iComp++) {
    if (CompType(iComp) == TYPE_DONNANHA) {
      TotMoles[iComp] = abs(WHAMSpecCharge[iHA]) * HumicSubstGramsPerLiter[iHA];
      TotConc[iComp] = TotMoles[iComp] / SpecCtoMAdj[iComp];
    } else if (CompType(iComp) == TYPE_DONNANFA) {
      TotMoles[iComp] = abs(WHAMSpecCharge[iFA]) * HumicSubstGramsPerLiter[iFA];
      TotConc[iComp] = TotMoles[iComp] / SpecCtoMAdj[iComp];
    }
  }

  // Update the Donnan layers the way WHAM does it
  for (iComp = 0; iComp < NComp; iComp++) {
    if ((CompType(iComp) == TYPE_DONNANHA) || (CompType(iComp) == TYPE_DONNANFA)) {

      SimpleAdjustComp(iComp, ConvCrit, MaxIter, TotMoles(iComp), 
                       NComp, CompConc,
                       NSpec, SpecConc, SpecKISTempAdj, SpecStoich, SpecName,
                       SpecActCorr, SpecActivityCoef, SpecCtoMAdj, SpecCharge,
                       WHAMSpecCharge);
    }
  }

}


