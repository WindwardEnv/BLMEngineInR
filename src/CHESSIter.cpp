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
//' @param HumicSubstGramsPerLiter,wMolWt,wRadius,wP,wDLF,wKZED WHAM parameters
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
  bool UpdateZED,
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
  double CATarget,
  Rcpp::NumericVector &MassAmtAdj,
  Rcpp::NumericVector &TotConc,
  Rcpp::NumericVector &TotMoles,
  Rcpp::NumericVector &SpecKISTempAdj,
  Rcpp::NumericVector &SpecCtoMAdj,
  Rcpp::NumericVector &SpecConc,
  Rcpp::NumericVector &SpecActivityCoef,
  Rcpp::NumericVector &CalcTotMoles,
  Rcpp::NumericVector &WHAMSpecCharge,
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
  //Rcpp::NumericVector TotMoles(NComp);
  //  TotMoles.names() = CompName;
  Rcpp::NumericVector CompCtoMAdj(NComp);
  //Rcpp::NumericVector SpecKISTempAdj(NSpec);
  //Rcpp::NumericVector SpecCtoMAdj(NSpec);
  //Rcpp::NumericVector SpecActivityCoef(NSpec);
  //Rcpp::NumericVector WHAMSpecCharge(2);
  //double ChargeBalance;
  //double WHAMIonicStrength;
  
  // Update the component free concentrations
  CompConc = SpecConc[CompPosInSpec];
  CompUpdate(NComp, CompConcStep, CompType, CompConc);
  SpecConc[CompPosInSpec] = clone(CompConc);
  
  if (false) {
    //FIRST/INITAL WHAM ITER FROM PB BLM
    IonicStrength = 0.1;
    WHAMSpecCharge[0] = -0.0001;
    WHAMSpecCharge[1] = -0.0001;
    //Inorganics
    CompConc[0] = 9.1253E-07;
    CompConc[1] = 8.1101E-05;
    CompConc[2] = 1.3596E-04;
    CompConc[3] = 1.7606E-04;
    CompConc[4] = 8.5928E-06;
    CompConc[5] = 2.1655E-04;
    CompConc[6] = 8.5738E-06;
    CompConc[7] = 3.3690E-10;
    CompConc[8] = 3.5523E-08;
    //Donnan ratios
    CompConc[10] = 10.0;
    CompConc[11] = 10.0;
    //DOC-HA-#H monodentate protonated sites
    CompConc[12] = 3.2178E-14;//1
    CompConc[13] = 1.2610E-13;//2
    CompConc[14] = 4.9356E-13;//3
    CompConc[15] = 1.9227E-12;//4
    CompConc[16] = 6.1071E-11;//5
    CompConc[17] = 2.0008E-10;//6
    CompConc[18] = 2.3923E-10;//7
    CompConc[19] = 2.4265E-10;//8
    //DOC-HA-##H bidentate protonated sites
    CompConc[20] = 1.22811E-18;//1
    CompConc[21] = 1.76377E-17;//2
    CompConc[22] = 2.60392E-16;//3
    CompConc[23] = 2.63466E-16;//4
    CompConc[24] = 1.76647E-17;//5
    CompConc[25] = 2.58895E-16;//6
    CompConc[26] = 2.70011E-16;//7
    CompConc[27] = 1.37116E-16;//8
    CompConc[28] = 2.71588E-16;//9
    CompConc[29] = 2.71803E-16;//10
    CompConc[30] = 2.71480E-16;//11
    CompConc[31] = 2.72246E-16;//12
    //DOC-FA_#H monodentate protonated sites
    CompConc[32] = 1.3651E-14;//1
    CompConc[33] = 1.7721E-13;//2
    CompConc[34] = 2.3002E-12;//3
    CompConc[35] = 2.9796E-11;//4
    CompConc[36] = 1.1013E-09;//5
    CompConc[37] = 5.9714E-09;//6
    CompConc[38] = 6.3791E-09;//7
    CompConc[39] = 6.3854E-09;//8
    //DOC-FA-##H bidentate protonated sites
    CompConc[40] = 6.3143E-20;//1
    CompConc[41] = 1.0617E-17;//2
    CompConc[42] = 4.2553E-15;//3
    CompConc[43] = 4.5503E-15;//4
    CompConc[44] = 1.0639E-17;//5
    CompConc[45] = 1.0188E-14;//6
    CompConc[46] = 5.9012E-14;//7
    CompConc[47] = 1.7888E-15;//8
    CompConc[48] = 7.1700E-13;//9
    CompConc[49] = 7.6671E-13;//10
    CompConc[50] = 1.7130E-12;//11
    CompConc[51] = 9.9220E-12;//12
  } else if (false) {
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

  // Calculate the ionic strength and activity coefficients
  //ChargeBalance = CalcChargeBalance(NSpec, SpecConc * SpecCtoMAdj, SpecCharge, 
  //                                  SpecMC, AqueousMC);
  IonicStrength = CalcIonicStrength(NSpec, SpecConc * SpecCtoMAdj, SpecCharge, 
                                    SpecMC, AqueousMC, SpecActCorr, true);
  SpecActivityCoef = CalcActivityCoef(NSpec, SpecName, SpecActCorr, SpecCharge, 
                                      IonicStrength, SysTempKelvin);
  
  if (DoWHAM) {
    //WHAMIonicStrength = CalcIonicStrength(NSpec, SpecConc * SpecCtoMAdj, SpecCharge, 
    //                                SpecMC, AqueousMC, SpecActCorr, true);
    AdjustForWHAMBeforeCalcSpecies(NMass, MassAmt, MassAmtAdj, NSpec, SpecMC,
      SpecActCorr, SpecCharge, SpecKTempAdj, SpecKISTempAdj, SpecCtoMAdj,
      IonicStrength, WHAMSpecCharge, AqueousMC, WHAMDonnanMC, 
      HumicSubstGramsPerLiter, wMolWt, wRadius, wP, wDLF, wKZED);
    CompCtoMAdj = SpecCtoMAdj[CompPosInSpec];
    TotMoles = TotConc * CompCtoMAdj;
  }

  //for (int iSpec = 0; iSpec < NSpec; iSpec++){
  //  Rcpp::Rcout << "SpecLogK[" << SpecName[iSpec] << "]=" << 
  //    log10(SpecKISTempAdj[iSpec]) << std::endl;
  //}

  // Calculate the species concentrations
  SpecConc = CalcSpecConc(NComp, NSpec, CompConc, SpecKISTempAdj, SpecStoich, 
                          SpecName, SpecActCorr, SpecActivityCoef);

  UpdateFixedComps(NComp, CompType, TotConc, SpecActivityCoef, 
                  SpecConc, CompConc);
  
  if (DoWHAM) {
    /*AdjustForWHAM(NMass, MassAmt, MassAmtAdj, 
                    NComp, CompType, TotConc, TotMoles,
                    NSpec, SpecConc, SpecMC, SpecActCorr, SpecCharge, 
                    SpecKTempAdj, SpecKISTempAdj, SpecCtoMAdj,  
                    IonicStrength, WHAMSpecCharge, AqueousMC, WHAMDonnanMC,
                    HumicSubstGramsPerLiter, wMolWt, wRadius, wP, wDLF, wKZED);*/
    AdjustForWHAMAfterCalcSpecies(NComp, CompType, TotConc, TotMoles, NSpec, 
      SpecConc, SpecActivityCoef, SpecMC, SpecActCorr, SpecCharge, SpecCtoMAdj, 
      WHAMSpecCharge, AqueousMC, HumicSubstGramsPerLiter, UpdateZED);
    
    AdjustDonnanRatio(NComp, NSpec, CompConc, CompType, TotMoles, SpecKISTempAdj,
                      SpecStoich, SpecName, SpecActCorr, SpecActivityCoef, SpecCtoMAdj);
    SpecConc = CalcSpecConc(NComp, NSpec, CompConc, SpecKISTempAdj, SpecStoich, 
                          SpecName, SpecActCorr, SpecActivityCoef);

  }       

  /*for (int iSpec = 0; iSpec < NSpec; iSpec++){
    Rcpp::Rcout << "SpecConc[" << SpecName[iSpec] << "]=" << 
      SpecConc[iSpec] << std::endl;
  }*/
  

  // Calculate the total moles & conc from species concentrations
  CalcTotMoles = CalcIterationTotalMoles(NComp, NSpec, SpecConc * SpecCtoMAdj, 
                                         SpecStoich);
  CalcIterationTotals(NComp, NSpec, SpecConc, SpecCtoMAdj, SpecStoich,
                      CalcTotMoles, CalcTotConc);

  double CalcTotDLCheck = 0;
  double ZCheck = 0;
  if (false) {
    for (int iSpec = 0; iSpec < NSpec; iSpec++) {
      if (((iSpec >= 12) && (iSpec <= 31)) || 
          ((iSpec >= 124) && (iSpec <= 247))) {
        ZCheck += ((SpecConc[iSpec] / HumicSubstGramsPerLiter[0]) * SpecCharge[iSpec]);
      }    
    }
    //ZCheck = ZCheck * HumicSubstGramsPerLiter[0];

    CalcTotDLCheck += SpecConc[74] * 2; //Cu
    CalcTotDLCheck += SpecConc[75] * 2; //Ca
    CalcTotDLCheck += SpecConc[76] * 2; //Mg
    CalcTotDLCheck += SpecConc[77] * 1; //Na
    CalcTotDLCheck += SpecConc[78] * 1; //K
    CalcTotDLCheck += SpecConc[81] * 2; //CO3
    CalcTotDLCheck += SpecConc[82] * 1;//H
    CalcTotDLCheck += SpecConc[84] * 1;//HCO3
    CalcTotDLCheck += SpecConc[86] * 1;//MgHCO3
    CalcTotDLCheck += SpecConc[89] * 1;//CaHCO3
    CalcTotDLCheck += SpecConc[92] * 1;//CuOH
    CalcTotDLCheck += SpecConc[95] * 1;//CuCl
    CalcTotDLCheck += SpecConc[98] * 1;//CuHCO3
    CalcTotDLCheck *= MassAmtAdj[WHAMDonnanMC[iHA]];
    
    //Tmp += SpecConc[85] * 0;//H2CO3
    //Tmp += SpecConc[87] * 0;//MgCO3
    //Tmp += SpecConc[88] * 0;//MgSO4
    //Tmp += SpecConc[90] * 0;//CaCO3
    //Tmp += SpecConc[91] * 0;//CaSO4
    //Tmp += SpecConc[93] * 0;//Cu(OH)2
    //Tmp += SpecConc[94] * 0;//CuSO4
    //Tmp += SpecConc[96] * 0;//CuCO3

    //Tmp = 0;
    //Tmp += SpecConc[83] * -1;//OH
    //Tmp += SpecConc[79] * -2; //SO4
    //Tmp += SpecConc[80] * -1; //Cl
    //Tmp += SpecConc[97] * -2;//Cu(CO3)2
    //Tmp = Tmp * MassAmtAdj[WHAMDonnanMC[iHA]];
  }

  // Calculate the residuals and error fraction for each component
  CalcResidAndError(NComp, CalcTotMoles, TotMoles, CompType, 
                    SpecActCorr,  
                    Resid, CompError);

  if (false){
    Rcpp::Rcout 
      << Resid[0] << "," 
      << Resid[10] << "," 
      << Resid[11] << "," 
      << Resid[12] << "," 
      << Resid[13] << "," 
      << Resid[17] << "," 
      << Resid[20] << "," 
      << Resid[22] << "," 
      << Resid[32] << "," 
      << Resid[33] << "," 
      << Resid[37] << "," 
      << Resid[40] << "," 
      << Resid[42] << "," 
      << WHAMSpecCharge[0] << "," 
      << WHAMSpecCharge[1] << "," 
      << MassAmtAdj[AqueousMC] << ","
      << MassAmtAdj[WHAMDonnanMC[0]] << ","
      << MassAmtAdj[WHAMDonnanMC[1]] << ","
      << std::endl;
  }
  

  // Adjust Resid and CompError for toxicity mode
  if (DoTox) {
    AdjustForToxMode(NBLMetal, BLMetalSpecs, MetalComp, CATarget, SpecConc,
                      Resid, CompError);
  }

  // Determine which component has the highest error fraction
  MaxError = MaxCompError(NComp, CompError, WhichMax);

  return MaxError;

}