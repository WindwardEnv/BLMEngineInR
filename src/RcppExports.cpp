// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// CHESS
Rcpp::List CHESS(Rcpp::String QuietFlag, double ConvergenceCriteria, int MaxIter, bool DoPartialStepsAlways, int NMass, Rcpp::CharacterVector MassName, Rcpp::NumericVector MassAmt, int NComp, Rcpp::CharacterVector CompName, Rcpp::CharacterVector CompType, Rcpp::NumericVector TotConc, int NSpec, Rcpp::CharacterVector SpecName, Rcpp::IntegerVector SpecMC, Rcpp::NumericVector SpecK, Rcpp::NumericVector SpecTempKelvin, Rcpp::NumericVector SpecDeltaH, Rcpp::IntegerMatrix SpecStoich, Rcpp::IntegerVector SpecCharge, Rcpp::CharacterVector SpecActCorr, bool DoWHAM, int AqueousMC, Rcpp::IntegerVector WHAMDonnanMC, Rcpp::NumericVector SolHS, Rcpp::NumericVector wMolWt, Rcpp::NumericVector wRadius, Rcpp::NumericVector wP, double wDLF, double wKZED, double SysTempKelvin, bool DoTox, Rcpp::String MetalName, int MetalComp, int NBLMetal, Rcpp::IntegerVector BLMetalSpecs, double CATarget);
RcppExport SEXP _BLMEngineInR_CHESS(SEXP QuietFlagSEXP, SEXP ConvergenceCriteriaSEXP, SEXP MaxIterSEXP, SEXP DoPartialStepsAlwaysSEXP, SEXP NMassSEXP, SEXP MassNameSEXP, SEXP MassAmtSEXP, SEXP NCompSEXP, SEXP CompNameSEXP, SEXP CompTypeSEXP, SEXP TotConcSEXP, SEXP NSpecSEXP, SEXP SpecNameSEXP, SEXP SpecMCSEXP, SEXP SpecKSEXP, SEXP SpecTempKelvinSEXP, SEXP SpecDeltaHSEXP, SEXP SpecStoichSEXP, SEXP SpecChargeSEXP, SEXP SpecActCorrSEXP, SEXP DoWHAMSEXP, SEXP AqueousMCSEXP, SEXP WHAMDonnanMCSEXP, SEXP SolHSSEXP, SEXP wMolWtSEXP, SEXP wRadiusSEXP, SEXP wPSEXP, SEXP wDLFSEXP, SEXP wKZEDSEXP, SEXP SysTempKelvinSEXP, SEXP DoToxSEXP, SEXP MetalNameSEXP, SEXP MetalCompSEXP, SEXP NBLMetalSEXP, SEXP BLMetalSpecsSEXP, SEXP CATargetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type QuietFlag(QuietFlagSEXP);
    Rcpp::traits::input_parameter< double >::type ConvergenceCriteria(ConvergenceCriteriaSEXP);
    Rcpp::traits::input_parameter< int >::type MaxIter(MaxIterSEXP);
    Rcpp::traits::input_parameter< bool >::type DoPartialStepsAlways(DoPartialStepsAlwaysSEXP);
    Rcpp::traits::input_parameter< int >::type NMass(NMassSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type MassName(MassNameSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type MassAmt(MassAmtSEXP);
    Rcpp::traits::input_parameter< int >::type NComp(NCompSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type CompName(CompNameSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type CompType(CompTypeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type TotConc(TotConcSEXP);
    Rcpp::traits::input_parameter< int >::type NSpec(NSpecSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type SpecName(SpecNameSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type SpecMC(SpecMCSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SpecK(SpecKSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SpecTempKelvin(SpecTempKelvinSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SpecDeltaH(SpecDeltaHSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type SpecStoich(SpecStoichSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type SpecCharge(SpecChargeSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type SpecActCorr(SpecActCorrSEXP);
    Rcpp::traits::input_parameter< bool >::type DoWHAM(DoWHAMSEXP);
    Rcpp::traits::input_parameter< int >::type AqueousMC(AqueousMCSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type WHAMDonnanMC(WHAMDonnanMCSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SolHS(SolHSSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type wMolWt(wMolWtSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type wRadius(wRadiusSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type wP(wPSEXP);
    Rcpp::traits::input_parameter< double >::type wDLF(wDLFSEXP);
    Rcpp::traits::input_parameter< double >::type wKZED(wKZEDSEXP);
    Rcpp::traits::input_parameter< double >::type SysTempKelvin(SysTempKelvinSEXP);
    Rcpp::traits::input_parameter< bool >::type DoTox(DoToxSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type MetalName(MetalNameSEXP);
    Rcpp::traits::input_parameter< int >::type MetalComp(MetalCompSEXP);
    Rcpp::traits::input_parameter< int >::type NBLMetal(NBLMetalSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type BLMetalSpecs(BLMetalSpecsSEXP);
    Rcpp::traits::input_parameter< double >::type CATarget(CATargetSEXP);
    rcpp_result_gen = Rcpp::wrap(CHESS(QuietFlag, ConvergenceCriteria, MaxIter, DoPartialStepsAlways, NMass, MassName, MassAmt, NComp, CompName, CompType, TotConc, NSpec, SpecName, SpecMC, SpecK, SpecTempKelvin, SpecDeltaH, SpecStoich, SpecCharge, SpecActCorr, DoWHAM, AqueousMC, WHAMDonnanMC, SolHS, wMolWt, wRadius, wP, wDLF, wKZED, SysTempKelvin, DoTox, MetalName, MetalComp, NBLMetal, BLMetalSpecs, CATarget));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BLMEngineInR_CHESS", (DL_FUNC) &_BLMEngineInR_CHESS, 36},
    {NULL, NULL, 0}
};

RcppExport void R_init_BLMEngineInR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
