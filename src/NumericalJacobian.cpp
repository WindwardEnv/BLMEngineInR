#include <Rcpp.h>
#include "CHESSFunctions.h"

Rcpp::NumericMatrix NumericalJacobian(
    int NMass,
    Rcpp::NumericVector MassAmt,
    int NComp,
    Rcpp::CharacterVector CompName,
    Rcpp::CharacterVector CompType,
    Rcpp::IntegerVector CompPosInSpec,
    int NSpec,
    Rcpp::CharacterVector SpecName,
    Rcpp::CharacterVector SpecType,
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
    Rcpp::NumericVector MassAmtAdj,
    Rcpp::NumericVector TotConc,
    Rcpp::NumericVector TotMoles,
    Rcpp::NumericVector SpecKISTempAdj,
    Rcpp::NumericVector SpecCtoMAdj,
    Rcpp::NumericVector SpecConc,
    Rcpp::NumericVector SpecActivityCoef,
    Rcpp::NumericVector WHAMSpecCharge,
    double IonicStrength,
    Rcpp::NumericVector Resid) {

  /* outputs: */
  Rcpp::NumericMatrix JacobianMatrix(NComp);
  Rcpp::rownames(JacobianMatrix) = CompName;
  Rcpp::colnames(JacobianMatrix) = CompName;

  /* variables: */
  int iComp1;//the row of the Jacobian = dResid(iComp1)
  int iComp2;//the column of the Jacobian = dCompConc(iComp2)
  int i;
  double dResidComp1;
  double dConcComp2;
  //double MaxErrorMod;
  Rcpp::NumericVector CompConcStepMod(NComp);
  Rcpp::NumericVector MassAmtAdjMod(NMass);
  Rcpp::NumericVector TotConcMod(NComp);
  Rcpp::NumericVector TotMolesMod(NComp);
  Rcpp::NumericVector SpecKISTempAdjMod(NSpec);
  Rcpp::NumericVector SpecCtoMAdjMod(NSpec);
  Rcpp::NumericVector SpecConcMod(NSpec);
  Rcpp::NumericVector SpecActivityCoefMod(NSpec);
  Rcpp::NumericVector CalcTotMolesMod(NComp);
  int WhichMaxMod;
  double IonicStrengthMod;
  Rcpp::NumericVector ResidMod(NComp);
  Rcpp::NumericVector CompErrorMod(NComp);
  Rcpp::NumericVector WHAMSpecChargeMod(2);

  for (iComp2 = 0; iComp2 < NComp; iComp2++) {
    // change the concentration of iComp2 by a tiny amount
    // Calculate the modified residual
    SpecKISTempAdjMod = clone(SpecKISTempAdj);
    SpecCtoMAdjMod = clone(SpecCtoMAdj);
    SpecConcMod = clone(SpecConc);
    SpecActivityCoef = clone(SpecActivityCoef);
    TotConcMod = clone(TotConc);
    TotMolesMod = clone(TotMoles);
    WHAMSpecChargeMod = clone(WHAMSpecCharge);
    MassAmtAdjMod = clone(MassAmtAdj);
    for (i = 0; i < NComp; i++) {
      if (i == iComp2) {
        CompConcStepMod[i] = -0.00001 * SpecConc[iComp2];
      } else {
        CompConcStepMod[i] = 0;
      }
    }

    CHESSIter(CompConcStepMod, NMass, MassAmt, NComp, CompName,
              CompType, CompPosInSpec, NSpec, SpecName, SpecType, SpecMC,
              SpecActCorr, SpecStoich, SpecCharge, SpecKTempAdj,
              DoWHAM, true, AqueousMC, WHAMDonnanMC,
              HumicSubstGramsPerLiter, wMolWt,
              wRadius, wP, wDLF, wKZED, SysTempKelvin, DoTox,
              MetalName, MetalComp, NBLMetal, BLMetalSpecs, CATarget,
              //variables that will be modified and returned again:
              MassAmtAdjMod, TotConcMod, TotMolesMod, SpecKISTempAdjMod,
              SpecCtoMAdjMod, SpecConcMod, SpecActivityCoefMod,
              CalcTotMolesMod, WHAMSpecChargeMod, WhichMaxMod,
              IonicStrengthMod, ResidMod, CompErrorMod);

    dConcComp2 = (SpecConcMod[iComp2] - SpecConc[iComp2]);
    
    for (iComp1 = 0; iComp1 < NComp; iComp1++) {
      dResidComp1 = (ResidMod[iComp1] - Resid[iComp1]);
      if (dConcComp2 == 0) {
        JacobianMatrix(iComp1, iComp2) = 0;
      } else {
        JacobianMatrix(iComp1, iComp2) = dResidComp1 / dConcComp2;
      }
    };//NEXT iComp1

  };//NEXT iComp2

  return JacobianMatrix;
}
