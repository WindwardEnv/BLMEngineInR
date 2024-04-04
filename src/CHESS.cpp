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
//' @param MetalComp integer, the position of the metal in the component arrays
//'   (i.e., which is the toxic metal component) Note: this is base-1 indexed on
//'   input then converted.
//' @param NBLMetal integer, the number of biotic ligand-bound metal species 
//'   that are associated with toxic effects.
//' @param BLMetalSpecs integer vector, the positions of the species in the
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
                 Rcpp::NumericVector HumicSubstGramsPerLiter,
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
  Rcpp::NumericVector SpecActivityCoef(NSpec);
    for (int iSpec = 0; iSpec < NSpec; iSpec++) {SpecActivityCoef[iSpec] = 1.0;}
  Rcpp::NumericVector WHAMSpecCharge(2);
    WHAMSpecCharge[0] = -0.0001;
    WHAMSpecCharge[1] = -0.0001;
  int WhichMax;
  double IonicStrength = 0.1;
  //double ChargeBalance = 0.0;
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
  Rcpp::NumericVector TotMolesFull(NComp);
  Rcpp::NumericVector SpecKISTempAdjFull(NSpec);
  Rcpp::NumericVector SpecCtoMAdjFull(NSpec);
  Rcpp::NumericVector SpecConcFull(NSpec);
  Rcpp::NumericVector SpecActivityCoefFull(NSpec);
  Rcpp::NumericVector CalcTotMolesFull(NComp);
  int WhichMaxFull;
  double IonicStrengthFull;
  Rcpp::NumericVector ResidFull(NComp);
  Rcpp::NumericVector CompErrorFull(NComp);
  Rcpp::NumericVector WHAMSpecChargeFull(2);

  double StepSizeAlt;
  double MaxErrorAlt;
  Rcpp::NumericVector CompConcStepAlt(NComp);
  Rcpp::NumericVector MassAmtAdjAlt(NMass);
  Rcpp::NumericVector TotConcAlt(NComp);
  Rcpp::NumericVector TotMolesAlt(NComp);
  Rcpp::NumericVector SpecKISTempAdjAlt(NSpec);
  Rcpp::NumericVector SpecCtoMAdjAlt(NSpec);
  Rcpp::NumericVector SpecConcAlt(NSpec);
  Rcpp::NumericVector SpecActivityCoefAlt(NSpec);
  Rcpp::NumericVector CalcTotMolesAlt(NComp);
  int WhichMaxAlt;
  double IonicStrengthAlt;
  Rcpp::NumericVector ResidAlt(NComp);
  Rcpp::NumericVector CompErrorAlt(NComp);
  Rcpp::NumericVector WHAMSpecChargeAlt(2);

  /*double StepSizeInterp;
  double MaxErrorInterp;
  Rcpp::NumericVector CompConcStepInterp(NComp);
  Rcpp::NumericVector MassAmtAdjInterp(NMass);
  Rcpp::NumericVector TotConcInterp(NComp);
  Rcpp::NumericVector TotMolesInterp(NComp);
  Rcpp::NumericVector SpecKISTempAdjInterp(NSpec);
  Rcpp::NumericVector SpecCtoMAdjInterp(NSpec);
  Rcpp::NumericVector SpecConcInterp(NSpec);
  Rcpp::NumericVector SpecActivityCoefInterp(NSpec);
  Rcpp::NumericVector CalcTotMolesInterp(NComp);
  int WhichMaxInterp;
  double IonicStrengthInterp;
  Rcpp::NumericVector ResidInterp(NComp);
  Rcpp::NumericVector CompErrorInterp(NComp);
  Rcpp::NumericVector WHAMSpecChargeInterp(2);*/

  int Counter;

  int CompToChange = 0;
  double CompChangeAmt = 0.0;
  
  // Initialize some variables
  for (iComp = 0; iComp < NComp; iComp++) {
    CompPosInSpec(iComp) = iComp;
  }
  SpecMC = SpecMC - 1;
  AqueousMC = AqueousMC - 1;
  WHAMDonnanMC = WHAMDonnanMC - 1;
  MetalComp = MetalComp - 1;
  BLMetalSpecs = BLMetalSpecs - 1;
  MassAmtAdj = clone(MassAmt);
  SpecCtoM = MassAmtAdj[SpecMC];
  SpecCtoMAdj = clone(SpecCtoM);
  CompCtoM = SpecCtoM[CompPosInSpec];
  TotMoles = TotConc * CompCtoM;

  // Do the temperature adjustments on the binding constants
  //SysTempKelvin = 288;
  SpecKTempAdj = TempCorrection(SysTempKelvin, NSpec, SpecK, SpecTempKelvin, 
                                SpecDeltaH);
  SpecKISTempAdj = clone(SpecKTempAdj);

  //for (int iSpec = 0; iSpec < NSpec; iSpec++){
  //  Rcpp::Rcout << SpecName[iSpec] << " logK = " << log10(SpecKTempAdj[iSpec]) << std::endl;
  //}
  

  // Get initial values for component concentrations
  CompConc = InitialGuess(TotConc, SpecCtoMAdj, CompType, SpecKISTempAdj, 
                          SpecStoich, SpecName, NComp, NSpec);
  SpecConc[CompPosInSpec] = clone(CompConc);

  if (false) {
    //LAST WHAM ITER FROM PB BLM
    
    WHAMSpecCharge[0] = -0.001055782;
    WHAMSpecCharge[1] = -0.002757042;
    TotMoles[10] = abs(WHAMSpecCharge[0]) * HumicSubstGramsPerLiter[0];
    TotMoles[11] = abs(WHAMSpecCharge[1]) * HumicSubstGramsPerLiter[1];
    MassAmtAdj[AqueousMC] =  0.999939401;
    MassAmtAdj[WHAMDonnanMC[iHA]] = 5.65E-07;
    MassAmtAdj[WHAMDonnanMC[iFA]] = 6.00E-05;
    for (int iSpec = 0; iSpec < NSpec; iSpec++) {    
      if ((SpecActCorr[iSpec] == "WHAMHA") || (SpecActCorr[iSpec] == "WHAMFA")) {
        SpecCtoMAdj[iSpec] = MassAmt[AqueousMC];
      } else {
        SpecCtoMAdj[iSpec] = MassAmtAdj[SpecMC[iSpec]];
      }
    }	
    TotConc[10] = TotMoles[10] / SpecCtoMAdj[10];
    TotConc[11] = TotMoles[11] / SpecCtoMAdj[11];   

    CompConc[0] = 9.59993363565727E-08;
    CompConc[1] = 3.68278171517698E-05;
    CompConc[2] = 0.000061754236008998;
    CompConc[3] = 0.000137555074214308;
    CompConc[4] = 6.71355101983406E-06;
    CompConc[5] = 9.85195915618222E-05;
    CompConc[6] = 6.69941662102341E-06;
    CompConc[7] = 1.94743882686021E-07;
    CompConc[8] = 2.76256091289306E-08;
    CompConc[10] = 4.001807;
    CompConc[11] = 1.712652;
    CompConc[12] = 3.90672800784937E-13;//HA monodentate start
    CompConc[13] = 1.47731944186428E-12;
    CompConc[14] = 5.08513478543844E-12;
    CompConc[15] = 1.34857790600327E-11;
    CompConc[16] = 1.45907868350776E-10;
    CompConc[17] = 1.68683394297017E-10;
    CompConc[18] = 1.70597781345935E-10;
    CompConc[19] = 1.70737081609003E-10;
    CompConc[20] = 1.62318730585003E-16;//HA bidentate start
    CompConc[21] = 8.72449447745709E-16;
    CompConc[22] = 1.62691821629889E-15;
    CompConc[23] = 1.62770370103955E-15;
    CompConc[24] = 8.77842067578099E-16;
    CompConc[25] = 1.6450367927837E-15;
    CompConc[26] = 1.64789379713423E-15;
    CompConc[27] = 1.23358746610706E-15;
    CompConc[28] = 1.65308852570421E-15;
    CompConc[29] = 1.65314162960448E-15;
    CompConc[30] = 1.65428833978588E-15;
    CompConc[31] = 1.65448109374049E-15;
    CompConc[32] = 8.56767594077975E-14;//FA monodentate start
    CompConc[33] = 1.09999978532507E-12;
    CompConc[34] = 1.24955490848289E-11;
    CompConc[35] = 6.18627774290784E-11;
    CompConc[36] = 3.58950110306434E-09;
    CompConc[37] = 6.22238094341529E-09;
    CompConc[38] = 6.28905895093181E-09;
    CompConc[39] = 6.29003322009559E-09;
    CompConc[40] = 1.23989844969647E-18;//FA bidentate start
    CompConc[41] = 1.15544203776794E-16;
    CompConc[42] = 1.29419496326124E-14;
    CompConc[43] = 1.30726020153824E-14;
    CompConc[44] = 1.16407613144101E-16;
    CompConc[45] = 6.23586088186863E-14;
    CompConc[46] = 8.38894475929258E-14;
    CompConc[47] = 2.59685374121812E-16;
    CompConc[48] = 1.43888965697356E-13;
    CompConc[49] = 1.43984008436532E-13;
    CompConc[50] = 1.51818397868478E-13;
    CompConc[51] = 1.52391684761759E-13;
    
    //for (int iComp = 10; iComp < NComp; iComp++) {
    //  if (iComp != 10) {CompType[iComp] = "FixedConc";}
    //}
    
    SpecConc[CompPosInSpec] = clone(CompConc);
    IonicStrength = 0.000534679;
    SpecActivityCoef = CalcActivityCoef(NSpec, SpecName, SpecActCorr, SpecCharge, 
                                      IonicStrength, SysTempKelvin);
    SpecConc = CalcSpecConc(NComp, NSpec, CompConc, SpecKTempAdj, SpecStoich, 
                            SpecName, SpecActCorr, SpecActivityCoef);
  }

  // Run through CHESS calculations with initial values
  MaxError = CHESSIter(CompConcStep, NMass, MassAmt, NComp, CompName, CompType,
                       CompPosInSpec, NSpec, SpecName, SpecMC, SpecActCorr,
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
  while ((MaxError > ConvergenceCriteria) & (Iter <= MaxIter)) {

    /*if (QuietFlag == "Debug") {
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
    UpdateZED = !UpdateZED;
    
    // Calculate the Jacobian Matrix
    JacobianMatrix = Jacobian(NComp, NSpec, CompName, TotConc, SpecStoich, SpecConc, 
                              SpecMC, SpecCtoMAdj, SpecActCorr, SpecCharge, 
                              SpecKISTempAdj, IonicStrength, DoWHAM, 
                              HumicSubstGramsPerLiter, WHAMSpecCharge, 
                              wP, wMolWt, wRadius, wDLF, wKZED, MassAmtAdj, 
                              AqueousMC, WHAMDonnanMC,  MetalComp, NBLMetal, BLMetalSpecs, 
                              DoTox);

    if (false) {

      Rcpp::NumericVector dZHAdC(NComp), dZFAdC(NComp);
      for (int j = 0; j < NComp; j++) {
        dZHAdC[j] = 0;        
        dZFAdC[j] = 0;
        for (int k = 0; k < NSpec; k++) {
          if (SpecActCorr[k] == "WHAMHA") {
            dZHAdC[j] += (SpecConc[k] * SpecCharge[k] * SpecStoich(k, j));
          } else if (SpecActCorr[k] == "WHAMFA") {
            dZFAdC[j] += (SpecConc[k] * SpecCharge[k] * SpecStoich(k, j));
          }
        }
        dZHAdC[j] = dZHAdC[j] / (HumicSubstGramsPerLiter[iHA] * SpecConc[j]);        
        dZFAdC[j] = dZFAdC[j] / (HumicSubstGramsPerLiter[iFA] * SpecConc[j]);
      }
      Rcpp::Rcout << "dZHAdC = " << dZHAdC << std::endl;
      Rcpp::Rcout << "dZFAdC = " << dZFAdC << std::endl;

      double L_HA, L_FA;
      L_HA = (1000 * (6.022E+23) / wMolWt[iHA]) * (4 * 3.14159 / 3) * 
             (pow(wRadius[iHA] - (3.04E-10) / sqrt(IonicStrength), 3) - 
              pow(wRadius[iHA], 3)) * 
             HumicSubstGramsPerLiter[iHA];
      L_FA = (1000 * (6.022E+23) / wMolWt[iFA]) * (4 * 3.14159 / 3) * 
             (pow(wRadius[iFA] - (3.04E-10) / sqrt(IonicStrength), 3) - 
              pow(wRadius[iFA], 3)) * 
             HumicSubstGramsPerLiter[iFA];
      
      Rcpp::NumericVector dVmaxHAdC(NComp), dVmaxFAdC(NComp);
      Rcpp::NumericVector dVHAdC(NComp), dVFAdC(NComp);
      for (int j = 0; j < NComp; j++) {
        dVmaxHAdC[j] = (dZHAdC[j] * L_HA * WHAMSpecCharge[iHA] * wKZED) / 
          (abs(WHAMSpecCharge[iHA]) * pow(1 + wKZED * abs(WHAMSpecCharge[iHA]), 2));
        dVmaxFAdC[j] = (dZFAdC[j] * L_FA * WHAMSpecCharge[iFA] * wKZED) / 
          (abs(WHAMSpecCharge[iFA]) * pow(1 + wKZED * abs(WHAMSpecCharge[iFA]), 2));

        dVHAdC[j] = (dVmaxHAdC[j] * wDLF * (wDLF + MassAmtAdj[WHAMDonnanMC[iHA]] + MassAmtAdj[WHAMDonnanMC[iFA]]) - 
          MassAmtAdj[WHAMDonnanMC[iHA]] * wDLF * (dVmaxHAdC[j] + dVmaxFAdC[j])) / 
          pow(wDLF + MassAmtAdj[WHAMDonnanMC[iHA]] + MassAmtAdj[WHAMDonnanMC[iHA]], 2);

        dVFAdC[j] = (dVmaxFAdC[j] * wDLF * (wDLF + MassAmtAdj[WHAMDonnanMC[iHA]] + MassAmtAdj[WHAMDonnanMC[iFA]]) - 
          MassAmtAdj[WHAMDonnanMC[iFA]] * wDLF * (dVmaxHAdC[j] + dVmaxFAdC[j])) / 
          pow(wDLF + MassAmtAdj[WHAMDonnanMC[iHA]] + MassAmtAdj[WHAMDonnanMC[iHA]], 2);

      }
      Rcpp::Rcout << "dVHAdC = " << dVHAdC << std::endl;
      Rcpp::Rcout << "dVFAdC = " << dVFAdC << std::endl;

      Rcpp::NumericVector dTDLHAdC_term1(NComp), dTDLFAdC_term1(NComp);
      Rcpp::NumericVector dTDLHAdC_term2(NComp), dTDLFAdC_term2(NComp);
      for (int j = 0; j < NComp; j++) {
        dTDLHAdC_term1[j] = 0;
        dTDLHAdC_term2[j] = 0;
        for (int i = 0; i < NSpec; i++) {
          dTDLHAdC_term1[j] += (SpecStoich(i, 10) * SpecStoich(i, j) * SpecConc[i]);
          dTDLHAdC_term2[j] += (SpecStoich(i, 10) * SpecConc[i]);
        }
        dTDLHAdC_term1[j] = dTDLHAdC_term1[j] * MassAmtAdj[WHAMDonnanMC[iHA]] / SpecConc[j];
        dTDLHAdC_term2[j] = dTDLHAdC_term2[j] * dVHAdC[j];

        dTDLFAdC_term1[j] = 0;
        dTDLFAdC_term2[j] = 0;
        for (int i = 0; i < NSpec; i++) {
          dTDLFAdC_term1[j] += (SpecStoich(i, 11) * SpecStoich(i, j) * SpecConc[i]);
          dTDLFAdC_term2[j] += (SpecStoich(i, 11) * SpecConc[i]);
        }
        dTDLFAdC_term1[j] = dTDLFAdC_term1[j] * MassAmtAdj[WHAMDonnanMC[iFA]] / SpecConc[j];
        dTDLFAdC_term2[j] = dTDLFAdC_term2[j] * dVFAdC[j];
      }
      Rcpp::Rcout << "dTDLHAdC_term1 = " << dTDLHAdC_term1 << std::endl;
      Rcpp::Rcout << "dTDLFAdC_term1 = " << dTDLFAdC_term1 << std::endl;
      Rcpp::Rcout << "dTDLHAdC_term2 = " << dTDLHAdC_term2 << std::endl;
      Rcpp::Rcout << "dTDLFAdC_term2 = " << dTDLFAdC_term2 << std::endl;

      Rcpp::NumericVector JacobianDLHACheck(NComp), JacobianDLFACheck(NComp);
      for (int j = 0; j < NComp; j++) {
        JacobianDLHACheck[j] = dTDLHAdC_term1[j] + dTDLHAdC_term2[j] - dZHAdC[j] * WHAMSpecCharge[iHA] * HumicSubstGramsPerLiter[iHA] / abs(WHAMSpecCharge[iHA]);
        JacobianDLFACheck[j] = dTDLFAdC_term1[j] + dTDLFAdC_term2[j] - dZFAdC[j] * WHAMSpecCharge[iFA] * HumicSubstGramsPerLiter[iFA] / abs(WHAMSpecCharge[iFA]);
      }
      Rcpp::Rcout << "JacobianDLHACheck = " << JacobianDLHACheck << std::endl;
      Rcpp::Rcout << "JacobianDLFACheck = " << JacobianDLFACheck << std::endl;
      
      Rcpp::NumericVector dRDOCHA1dC_term1(NComp), dRDOCFA1dC_term1(NComp);
      Rcpp::NumericVector dRDOCHA12dC_term1(NComp), dRDOCFA12dC_term1(NComp);
      Rcpp::NumericVector dRDOCHA1dC_term2(NComp), dRDOCFA1dC_term2(NComp);
      Rcpp::NumericVector dRDOCHA12dC_term2(NComp), dRDOCFA12dC_term2(NComp);
      Rcpp::NumericVector JacobianDOCHA1Check(NComp), JacobianDOCFA1Check(NComp);
      Rcpp::NumericVector JacobianDOCHA12Check(NComp), JacobianDOCFA12Check(NComp);
      for (int j = 0; j < NComp; j++) {
        dRDOCHA1dC_term1[j] = 0;
        dRDOCFA1dC_term1[j] = 0;
        dRDOCHA12dC_term1[j] = 0;
        dRDOCFA12dC_term1[j] = 0;
        dRDOCHA1dC_term2[j] = 0;
        dRDOCFA1dC_term2[j] = 0;
        dRDOCHA12dC_term2[j] = 0;
        dRDOCFA12dC_term2[j] = 0;
        for (int i = 0; i < NSpec; i++) {
          dRDOCHA1dC_term1[j] += (SpecStoich(i, 12) * SpecStoich(i, j) * SpecConc[i] / SpecConc[j]);
          dRDOCFA1dC_term1[j] += (SpecStoich(i, 32) * SpecStoich(i, j) * SpecConc[i] / SpecConc[j]);
          dRDOCHA12dC_term1[j] += (SpecStoich(i, 20) * SpecStoich(i, j) * SpecConc[i] / SpecConc[j]);
          dRDOCFA12dC_term1[j] += (SpecStoich(i, 40) * SpecStoich(i, j) * SpecConc[i] / SpecConc[j]);

          dRDOCHA1dC_term2[j] += (-2 * wP[iHA] * log10(IonicStrength) * SpecCharge[i] * SpecStoich(i, 12) * SpecConc[i] * dZHAdC[j]);
          dRDOCFA1dC_term2[j] += (-2 * wP[iFA] * log10(IonicStrength) * SpecCharge[i] * SpecStoich(i, 32) * SpecConc[i] * dZFAdC[j]);
          dRDOCHA12dC_term2[j] += (-2 * wP[iHA] * log10(IonicStrength) * SpecCharge[i] * SpecStoich(i, 20) * SpecConc[i] * dZHAdC[j]);
          dRDOCFA12dC_term2[j] += (-2 * wP[iFA] * log10(IonicStrength) * SpecCharge[i] * SpecStoich(i, 40) * SpecConc[i] * dZFAdC[j]);

          /*JacobianDOCHA1Check[j] += (SpecStoich(i, 12) * SpecConc[i] * 
            (dZHAdC[j] * (-2 * SpecCharge[i] * wP[iHA] * log10(IonicStrength)) + 
              SpecStoich(i, j) / SpecConc[j]));
          JacobianDOCFA1Check[j] += (SpecStoich(i, 32) * SpecConc[i] * 
            (dZFAdC[j] * (-2 * SpecCharge[i] * wP[iFA] * log10(IonicStrength)) + 
              SpecStoich(i, j) / SpecConc[j]));
          JacobianDOCHA12Check[j] += (SpecStoich(i, 20) * SpecConc[i] * 
            (dZHAdC[j] * (-2 * SpecCharge[i] * wP[iHA] * log10(IonicStrength)) + 
              SpecStoich(i, j) / SpecConc[j]));
          JacobianDOCFA12Check[j] += (SpecStoich(i, 40) * SpecConc[i] * 
            (dZFAdC[j] * (-2 * SpecCharge[i] * wP[iFA] * log10(IonicStrength)) + 
              SpecStoich(i, j) / SpecConc[j]));*/
        }    
        JacobianDOCHA1Check[j] = dRDOCHA1dC_term1[j] + dRDOCHA1dC_term2[j];
        JacobianDOCFA1Check[j] = dRDOCFA1dC_term1[j] + dRDOCFA1dC_term2[j];
        JacobianDOCHA12Check[j] = dRDOCHA12dC_term1[j] + dRDOCHA12dC_term2[j];
        JacobianDOCFA12Check[j] = dRDOCFA12dC_term1[j] + dRDOCFA12dC_term2[j];
      }
      Rcpp::Rcout << "dRDOCHA1dC_term1 = " << dRDOCHA1dC_term1 << std::endl;
      Rcpp::Rcout << "dRDOCFA1dC_term1 = " << dRDOCFA1dC_term1 << std::endl;
      Rcpp::Rcout << "dRDOCHA12dC_term1 = " << dRDOCHA12dC_term1 << std::endl;
      Rcpp::Rcout << "dRDOCFA12dC_term1 = " << dRDOCFA12dC_term1 << std::endl;
      Rcpp::Rcout << "dRDOCHA1dC_term2 = " << dRDOCHA1dC_term2 << std::endl;
      Rcpp::Rcout << "dRDOCFA1dC_term2 = " << dRDOCFA1dC_term2 << std::endl;
      Rcpp::Rcout << "dRDOCHA12dC_term2 = " << dRDOCHA12dC_term2 << std::endl;
      Rcpp::Rcout << "dRDOCFA12dC_term2 = " << dRDOCFA12dC_term2 << std::endl;
      Rcpp::Rcout << "JacobianDOCHA1Check = " << JacobianDOCHA1Check << std::endl;
      Rcpp::Rcout << "JacobianDOCFA1Check = " << JacobianDOCFA1Check << std::endl;
      Rcpp::Rcout << "JacobianDOCHA12Check = " << JacobianDOCHA12Check << std::endl;
      Rcpp::Rcout << "JacobianDOCFA12Check = " << JacobianDOCFA12Check << std::endl;
      
      
      double CalcTotCheck, TotCheck, ResidCheck;
      for (int iSpec = 0; iSpec < NComp; iSpec++) {
        Rcpp::Rcout << SpecConc[iSpec] << ",";
      }
      for (int iComp = 0; iComp < NComp; iComp++) {
        if ((SpecActCorr[iComp] == "WHAMHA") || (SpecActCorr[iComp] == "WHAMFA")) {
          CalcTotCheck = 0;
          for (int iSpec = 0; iSpec < NSpec; iSpec++) {
            CalcTotCheck += SpecStoich(iSpec, iComp) * SpecConc[iSpec];
          }
          TotCheck = TotConc[iComp];
        } else {
          CalcTotCheck = 0;
          for (int iSpec = 0; iSpec < NSpec; iSpec++) {
            CalcTotCheck += SpecStoich(iSpec, iComp) * SpecConc[iSpec] * SpecCtoMAdj[iSpec];
          }
          if ((SpecActCorr[iComp] == "DonnanHA")) {
            TotCheck = 0;
            for (int iSpec = 0; iSpec < NSpec; iSpec++) {
              if (SpecActCorr[iSpec] == "WHAMHA") {
                TotCheck += SpecConc[iSpec] * SpecCharge[iSpec];
              }            
            }
            TotCheck = abs(TotCheck);
          } else if ((SpecActCorr[iComp] == "DonnanFA")) {
            TotCheck = 0;
            for (int iSpec = 0; iSpec < NSpec; iSpec++) {
              if (SpecActCorr[iSpec] == "WHAMFA") {
                TotCheck += SpecConc[iSpec] * SpecCharge[iSpec];
              }            
            }
            TotCheck = abs(TotCheck);
          } else {
            TotCheck = TotConc[iComp] * SpecCtoMAdj[iComp];
          }
        }
        ResidCheck = CalcTotCheck - TotCheck;
        Rcpp::Rcout << Resid[iComp] << ",";
      }
      
      Rcpp::Rcout 
        << WHAMSpecCharge[0] << "," 
        << WHAMSpecCharge[1] << "," 
        << MassAmtAdj[AqueousMC] << ","
        << MassAmtAdj[WHAMDonnanMC[0]] << ","
        << MassAmtAdj[WHAMDonnanMC[1]] << ",";
      for (int iComp1 = 0; iComp1 < NComp; iComp1++) {
        if ((iComp1 == 10) || (iComp1 == 11) || (iComp1 == 12) || (iComp1 == 20)){
          for (int iComp2 = 0; iComp2 < NComp; iComp2++) {
            if ((iComp2 == 0) || (iComp2 == 3) || (iComp2 == 10) || (iComp2 == 11) || (iComp2 == 12) || (iComp2 == 20)){
              Rcpp::Rcout << JacobianMatrix(iComp1, iComp2) << ",";              
            }
          }
        }        
      }
      Rcpp::Rcout << std::endl << std::endl;

      if (false) {
          //LAST WHAM ITER FROM PB BLM
          
          WHAMSpecCharge[0] = -0.001055782;
          WHAMSpecCharge[1] = -0.002757042;
          MassAmtAdj[AqueousMC] =  0.999939401;
          MassAmtAdj[WHAMDonnanMC[iHA]] = 5.65E-07;
          MassAmtAdj[WHAMDonnanMC[iFA]] = 6.00E-05;
          for (int iSpec = 0; iSpec < NSpec; iSpec++) {    
            if ((SpecActCorr[iSpec] == "WHAMHA") || (SpecActCorr[iSpec] == "WHAMFA")) {
              SpecCtoMAdj[iSpec] = MassAmt[AqueousMC];
            } else {
              SpecCtoMAdj[iSpec] = MassAmtAdj[SpecMC[iSpec]];
            }
          }	

          CompConc[0] = 9.59993363565727E-08;
          CompConc[1] = 3.68278171517698E-05;
          CompConc[2] = 0.000061754236008998;
          CompConc[3] = 0.000137555074214308;
          CompConc[4] = 6.71355101983406E-06;
          CompConc[5] = 9.85195915618222E-05;
          CompConc[6] = 6.69941662102341E-06;
          CompConc[7] = 1.94743882686021E-07;
          CompConc[8] = 2.76256091289306E-08;
          CompConc[10] = 4.001807;
          CompConc[11] = 1.712652;
          CompConc[12] = 3.90672800784937E-13;
          CompConc[13] = 1.47731944186428E-12;
          CompConc[14] = 5.08513478543844E-12;
          CompConc[15] = 1.34857790600327E-11;
          CompConc[16] = 1.45907868350776E-10;
          CompConc[17] = 1.68683394297017E-10;
          CompConc[18] = 1.70597781345935E-10;
          CompConc[19] = 1.70737081609003E-10;
          CompConc[20] = 1.62318730585003E-16;
          CompConc[21] = 8.72449447745709E-16;
          CompConc[22] = 1.62691821629889E-15;
          CompConc[23] = 1.62770370103955E-15;
          CompConc[24] = 8.77842067578099E-16;
          CompConc[25] = 1.6450367927837E-15;
          CompConc[26] = 1.64789379713423E-15;
          CompConc[27] = 1.23358746610706E-15;
          CompConc[28] = 1.65308852570421E-15;
          CompConc[29] = 1.65314162960448E-15;
          CompConc[30] = 1.65428833978588E-15;
          CompConc[31] = 1.65448109374049E-15;
          CompConc[32] = 8.56767594077975E-14;
          CompConc[33] = 1.09999978532507E-12;
          CompConc[34] = 1.24955490848289E-11;
          CompConc[35] = 6.18627774290784E-11;
          CompConc[36] = 3.58950110306434E-09;
          CompConc[37] = 6.22238094341529E-09;
          CompConc[38] = 6.28905895093181E-09;
          CompConc[39] = 6.29003322009559E-09;
          CompConc[40] = 1.23989844969647E-18;
          CompConc[41] = 1.15544203776794E-16;
          CompConc[42] = 1.29419496326124E-14;
          CompConc[43] = 1.30726020153824E-14;
          CompConc[44] = 1.16407613144101E-16;
          CompConc[45] = 6.23586088186863E-14;
          CompConc[46] = 8.38894475929258E-14;
          CompConc[47] = 2.59685374121812E-16;
          CompConc[48] = 1.43888965697356E-13;
          CompConc[49] = 1.43984008436532E-13;
          CompConc[50] = 1.51818397868478E-13;
          CompConc[51] = 1.52391684761759E-13;
          
          SpecConc[CompPosInSpec] = clone(CompConc);
          IonicStrength = 0.000534679;
          SpecActivityCoef = CalcActivityCoef(NSpec, SpecName, SpecActCorr, SpecCharge, 
                                            IonicStrength, SysTempKelvin);
          SpecConc = CalcSpecConc(NComp, NSpec, CompConc, SpecKTempAdj, SpecStoich, 
                                  SpecName, SpecActCorr, SpecActivityCoef);
        }

      
      CompToChange = 3; 
      for (int iComp = 0; iComp < NComp; iComp++) {
        CompConcStepFull[iComp] = 0.0;
      }    
      CompChangeAmt -= 0.1 * SpecConc[CompToChange];
      if (CompChangeAmt > SpecConc[CompToChange]) {
        CompChangeAmt = 0.9 * SpecConc[CompToChange];
      }
      CompConcStepFull[CompToChange] = CompChangeAmt;

    } else {
      // Calculate the Newton-Raphson step
      CompConcStepFull = CalcStep(JacobianMatrix, Resid, NComp, CompType, CompName);
    }


    if (QuietFlag == "Debug") {
      Rcpp::Rcout << "JacobianMatrix = [" <<std::endl << JacobianMatrix << "]" << std::endl;
      Rcpp::Rcout << "iComp\tSpecName\tSpecConc\tResid\tError\tStep" << std::endl;
      for (iComp = 0; iComp < NComp; iComp++) {
        Rcpp::Rcout << iComp << "\t" 
                    << SpecName[iComp] << "\t" 
                    << SpecConc[iComp] << "   " 
                    << Resid[iComp] << "   " 
                    << CompError[iComp] << "   "
                    << CompConcStepFull[iComp]
                    << std::endl;
        //Rcpp::Rcout << "Resid[" << SpecName[iComp] << "]=" << 
        //  Resid[iComp] << std::endl;
        //Rcpp::Rcout << "CompConcStepFull[" << SpecName[iComp] << "]=" << 
        //  CompConcStepFull[iComp] << std::endl;
      }
    }   

    // Do a full N-R step
    SpecKISTempAdjFull = clone(SpecKISTempAdj);
    SpecCtoMAdjFull = clone(SpecCtoMAdj);
    SpecConcFull = clone(SpecConc);
    SpecActivityCoefFull = clone(SpecActivityCoef);
    TotConcFull = clone(TotConc);
    TotMolesFull = clone(TotMoles);
    WHAMSpecChargeFull = clone(WHAMSpecCharge);
    MassAmtAdjFull = clone(MassAmtAdj);
    MaxErrorFull = CHESSIter(CompConcStepFull, NMass, MassAmt, NComp, CompName, 
                        CompType, CompPosInSpec, NSpec, SpecName, SpecMC, 
                        SpecActCorr, SpecStoich, SpecCharge, SpecKTempAdj, 
                        DoWHAM, UpdateZED, AqueousMC, WHAMDonnanMC, HumicSubstGramsPerLiter, wMolWt, 
                        wRadius, wP, wDLF, wKZED, SysTempKelvin, DoTox, 
                        MetalName, MetalComp, NBLMetal, BLMetalSpecs, CATarget, 
                        MassAmtAdjFull, TotConcFull, TotMolesFull, SpecKISTempAdjFull, 
                        SpecCtoMAdjFull, SpecConcFull, SpecActivityCoefFull,
                        CalcTotMolesFull, WHAMSpecChargeFull, WhichMaxFull, 
                        IonicStrengthFull, ResidFull, CompErrorFull);

    DoPartialSteps = false;
    if (DoPartialStepsAlways) {
      if (MaxErrorFull > MaxError) {
        // Do a half-step
        StepSizeAlt = std::pow(0.5, 0.5);//this just ensures it starts at 0.5
        DoPartialSteps = true;
        if (QuietFlag == "Debug") { Rcpp::Rcout << "Little Step." << std::endl; }
      } else if (MaxErrorFull > (MaxError * 0.9)) {
        // Do a double-step
        StepSizeAlt = std::pow(1.5, 0.5);
        DoPartialSteps = true;
        if (QuietFlag == "Debug") { Rcpp::Rcout << "Big Step." << std::endl; }
      }
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
        SpecActivityCoef = clone(SpecActivityCoef);
        TotConcAlt = clone(TotConc);
        TotMolesAlt = clone(TotMoles);
        WHAMSpecChargeAlt = clone(WHAMSpecCharge);
        MassAmtAdjAlt = clone(MassAmtAdj);
        CompConcStepAlt = CompConcStepFull * StepSizeAlt;
        MaxErrorAlt = CHESSIter(CompConcStepAlt, NMass, MassAmt, NComp, CompName, 
                            CompType, CompPosInSpec, NSpec, SpecName, SpecMC, 
                            SpecActCorr, SpecStoich, SpecCharge, SpecKTempAdj, 
                            DoWHAM, UpdateZED, AqueousMC, WHAMDonnanMC, HumicSubstGramsPerLiter, wMolWt, 
                            wRadius, wP, wDLF, wKZED, SysTempKelvin, DoTox, 
                            MetalName, MetalComp, NBLMetal, BLMetalSpecs, CATarget, 
                            MassAmtAdjAlt, TotConcAlt, TotMolesAlt, SpecKISTempAdjAlt, 
                            SpecCtoMAdjAlt, SpecConcAlt, SpecActivityCoefAlt, 
                            CalcTotMolesAlt, WHAMSpecChargeAlt, WhichMaxAlt, 
                            IonicStrengthAlt, ResidAlt, CompErrorAlt);
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
        TotMolesInterp = clone(TotMoles);
        WHAMSpecChargeInterp = clone(WHAMSpecCharge);
        MassAmtAdjInterp = clone(MassAmtAdj);
        CompConcStepInterp = CompConcStepFull * StepSizeInterp;
        MaxErrorInterp = CHESSIter(CompConcStepInterp, NMass, MassAmt, NComp, CompName, 
                            CompType, CompPosInSpec, NSpec, SpecName, SpecMC, 
                            SpecActCorr, SpecStoich, SpecCharge, SpecKTempAdj, 
                            DoWHAM, WHAMSpecChargeInterp, AqueousMC, WHAMDonnanMC, HumicSubstGramsPerLiter, wMolWt, 
                            wRadius, wP, wDLF, wKZED, SysTempKelvin, DoTox, 
                            MetalName, MetalComp, NBLMetal, BLMetalSpecs, CATarget, 
                            MassAmtAdjInterp, TotConcInterp, TotMolesInterp, SpecKISTempAdjInterp, 
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
      TotMoles = clone(TotMolesInterp);
      SpecKISTempAdj = clone(SpecKISTempAdjInterp);
      SpecCtoMAdj = clone(SpecCtoMAdjInterp);
      SpecConc = clone(SpecConcInterp);
      SpecActivityCoef = clone(SpecActivityCoefInterp);
      CalcTotMoles = clone(CalcTotMolesInterp);
      WHAMSpecCharge = clone(WHAMSpecChargeInterp);
    } else */if (MaxErrorAlt < MaxError) {
      //Store Alt in the main veriables
      MaxError = MaxErrorAlt;
      WhichMax = WhichMaxAlt;
      IonicStrength = IonicStrengthAlt;
      Resid = clone(ResidAlt);
      CompError = clone(CompErrorAlt);
      MassAmtAdj = clone(MassAmtAdjAlt);
      TotConc = clone(TotConcAlt);
      TotMoles = clone(TotMolesAlt);
      SpecKISTempAdj = clone(SpecKISTempAdjAlt);
      SpecCtoMAdj = clone(SpecCtoMAdjAlt);
      SpecConc = clone(SpecConcAlt);
      SpecActivityCoef = clone(SpecActivityCoefAlt);
      CalcTotMoles = clone(CalcTotMolesAlt);
      WHAMSpecCharge = clone(WHAMSpecChargeAlt);
    } else {
      // Partial Steps has failed, just revert to full
      MaxError = MaxErrorFull;
      WhichMax = WhichMaxFull;
      IonicStrength = IonicStrengthFull;
      Resid = clone(ResidFull);
      CompError = clone(CompErrorFull);
      MassAmtAdj = clone(MassAmtAdjFull);
      TotConc = clone(TotConcFull);
      TotMoles = clone(TotMolesFull);
      SpecKISTempAdj = clone(SpecKISTempAdjFull);
      SpecCtoMAdj = clone(SpecCtoMAdjFull);
      SpecConc = clone(SpecConcFull);
      SpecActivityCoef = clone(SpecActivityCoefFull);
      CalcTotMoles = clone(CalcTotMolesFull);
      WHAMSpecCharge = clone(WHAMSpecChargeFull);
    }
    
    SpecMoles = SpecConc * SpecCtoMAdj;
    CompCtoMAdj = SpecCtoMAdj[CompPosInSpec];
    //TotMoles = TotConc * CompCtoMAdj;
    CalcTotConc = CalcTotMoles / CompCtoMAdj;

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
      Rcpp::Named("IonicStrength") = IonicStrength,
      Rcpp::Named("SpecConc") = SpecConc,
      Rcpp::Named("SpecAct") = SpecAct,
      Rcpp::Named("SpecMoles") = SpecMoles,
      Rcpp::Named("CalcTotConc") = CalcTotConc,
      Rcpp::Named("CalcTotMoles") = CalcTotMoles,
      Rcpp::Named("WHAMSpecCharge") = WHAMSpecCharge,
      Rcpp::Named("MassAmt") = MassAmtAdj);
}
