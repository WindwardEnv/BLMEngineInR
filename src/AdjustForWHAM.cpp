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
//' @param SolHS NumericVector, length of 2, moles of each organic matter 
//'   component in solution
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
  Rcpp::NumericVector SolHS,
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
  WHAMSpecCharge = CalcWHAMSpecCharge(NSpec, SpecActCorr,
                                      SpecConc * SpecCtoMAdj, 
                                      SpecCharge, SpecMC, AqueousMC);

  //Adjust organic matter specific binding based on ionic strength
  SpecKISTempAdj = CalcIonicStrengthEffects(IonicStrength, WHAMSpecCharge, 
                                            NSpec, SpecCharge, SpecKTempAdj, 
                                            SpecActCorr, wP);

  //Calculate the portion of the solution that's in the diffuse layer
  MassAmtAdj = CalcDonnanLayerVolume(NMass, NSpec, IonicStrength, MassAmt, 
                                     AqueousMC, WHAMDonnanMC, wMolWt, 
                                     wRadius, wDLF, wKZED, WHAMSpecCharge, 
                                     SolHS);
  SpecCtoMAdj = MassAmtAdj[SpecMC];

  // Set -Z_HS to the "known" total for the Donnan component
  for (iComp = 0; iComp < NComp; iComp++) {
    if (CompType(iComp) == "DonnanHA") {
      TotMoles(iComp) = abs(WHAMSpecCharge(iHA));
      TotConc(iComp) = abs(WHAMSpecCharge(iHA)) / SpecCtoMAdj(iComp);
    } else if (CompType(iComp) == "DonnanFA") {
      TotMoles(iComp) = abs(WHAMSpecCharge(iFA));
      TotConc(iComp) = abs(WHAMSpecCharge(iFA)) / SpecCtoMAdj(iComp);
    }
  }

  /* The cation species concentrations should be canceled out if Z_HS is 
     positive, and anion species concentrations should be canceled out if Z_HS
     is negative */
  for (iSpec = NComp; iSpec < NSpec; iSpec++) {
    if ((SpecActCorr(iSpec) == "DonnanHA") && 
        (((WHAMSpecCharge(iHA) < 0) && (SpecCharge(iSpec) < 0)) || 
         ((WHAMSpecCharge(iHA) > 0) && (SpecCharge(iSpec) > 0)))) {
      SpecConc(iSpec) = 0;
    } else if ((SpecActCorr(iSpec) == "DonnanFA") &&
               (((WHAMSpecCharge(iFA) < 0) && (SpecCharge(iSpec) < 0)) || 
                ((WHAMSpecCharge(iFA) > 0) && (SpecCharge(iSpec) > 0)))) {
      SpecConc(iSpec) = 0;
    }
  }

}