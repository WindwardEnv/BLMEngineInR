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
Rcpp::List CHESS(Rcpp::String QuietFlag, double ConvergenceCriteria, int MaxIter, int NMass, Rcpp::CharacterVector MassName, Rcpp::NumericVector MassAmt, int NComp, Rcpp::CharacterVector CompName, Rcpp::CharacterVector CompType, Rcpp::NumericVector TotConc, int NSpec, Rcpp::CharacterVector SpecName, Rcpp::CharacterVector SpecType, Rcpp::IntegerVector SpecMCR, Rcpp::NumericVector SpecK, Rcpp::NumericVector SpecTempKelvin, Rcpp::NumericVector SpecDeltaH, Rcpp::IntegerMatrix SpecStoich, Rcpp::IntegerVector SpecCharge, Rcpp::CharacterVector SpecActCorr, bool DoWHAM, int AqueousMCR, Rcpp::IntegerVector WHAMDonnanMCR, Rcpp::NumericVector HumicSubstGramsPerLiter, Rcpp::NumericVector WHAMMolWt, Rcpp::NumericVector WHAMRadius, Rcpp::NumericVector WHAMP, double WHAMDLF, double WHAMKZED, double SysTempKelvin, bool DoTox, Rcpp::String MetalName, int MetalCompR, int BLCompR, int NBLMetal, Rcpp::IntegerVector BLMetalSpecsR, double CATarget, bool DodVidCj, bool DodVidCjDonnan, bool DodKidCj, bool DoGammai, bool DoJacDonnan, bool DoJacWHAM, bool DoWHAMSimpleAdjust, bool DoDonnanSimpleAdjust);
RcppExport SEXP _BLMEngineInR_CHESS(SEXP QuietFlagSEXP, SEXP ConvergenceCriteriaSEXP, SEXP MaxIterSEXP, SEXP NMassSEXP, SEXP MassNameSEXP, SEXP MassAmtSEXP, SEXP NCompSEXP, SEXP CompNameSEXP, SEXP CompTypeSEXP, SEXP TotConcSEXP, SEXP NSpecSEXP, SEXP SpecNameSEXP, SEXP SpecTypeSEXP, SEXP SpecMCRSEXP, SEXP SpecKSEXP, SEXP SpecTempKelvinSEXP, SEXP SpecDeltaHSEXP, SEXP SpecStoichSEXP, SEXP SpecChargeSEXP, SEXP SpecActCorrSEXP, SEXP DoWHAMSEXP, SEXP AqueousMCRSEXP, SEXP WHAMDonnanMCRSEXP, SEXP HumicSubstGramsPerLiterSEXP, SEXP WHAMMolWtSEXP, SEXP WHAMRadiusSEXP, SEXP WHAMPSEXP, SEXP WHAMDLFSEXP, SEXP WHAMKZEDSEXP, SEXP SysTempKelvinSEXP, SEXP DoToxSEXP, SEXP MetalNameSEXP, SEXP MetalCompRSEXP, SEXP BLCompRSEXP, SEXP NBLMetalSEXP, SEXP BLMetalSpecsRSEXP, SEXP CATargetSEXP, SEXP DodVidCjSEXP, SEXP DodVidCjDonnanSEXP, SEXP DodKidCjSEXP, SEXP DoGammaiSEXP, SEXP DoJacDonnanSEXP, SEXP DoJacWHAMSEXP, SEXP DoWHAMSimpleAdjustSEXP, SEXP DoDonnanSimpleAdjustSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type QuietFlag(QuietFlagSEXP);
    Rcpp::traits::input_parameter< double >::type ConvergenceCriteria(ConvergenceCriteriaSEXP);
    Rcpp::traits::input_parameter< int >::type MaxIter(MaxIterSEXP);
    Rcpp::traits::input_parameter< int >::type NMass(NMassSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type MassName(MassNameSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type MassAmt(MassAmtSEXP);
    Rcpp::traits::input_parameter< int >::type NComp(NCompSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type CompName(CompNameSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type CompType(CompTypeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type TotConc(TotConcSEXP);
    Rcpp::traits::input_parameter< int >::type NSpec(NSpecSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type SpecName(SpecNameSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type SpecType(SpecTypeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type SpecMCR(SpecMCRSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SpecK(SpecKSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SpecTempKelvin(SpecTempKelvinSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SpecDeltaH(SpecDeltaHSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type SpecStoich(SpecStoichSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type SpecCharge(SpecChargeSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type SpecActCorr(SpecActCorrSEXP);
    Rcpp::traits::input_parameter< bool >::type DoWHAM(DoWHAMSEXP);
    Rcpp::traits::input_parameter< int >::type AqueousMCR(AqueousMCRSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type WHAMDonnanMCR(WHAMDonnanMCRSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type HumicSubstGramsPerLiter(HumicSubstGramsPerLiterSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type WHAMMolWt(WHAMMolWtSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type WHAMRadius(WHAMRadiusSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type WHAMP(WHAMPSEXP);
    Rcpp::traits::input_parameter< double >::type WHAMDLF(WHAMDLFSEXP);
    Rcpp::traits::input_parameter< double >::type WHAMKZED(WHAMKZEDSEXP);
    Rcpp::traits::input_parameter< double >::type SysTempKelvin(SysTempKelvinSEXP);
    Rcpp::traits::input_parameter< bool >::type DoTox(DoToxSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type MetalName(MetalNameSEXP);
    Rcpp::traits::input_parameter< int >::type MetalCompR(MetalCompRSEXP);
    Rcpp::traits::input_parameter< int >::type BLCompR(BLCompRSEXP);
    Rcpp::traits::input_parameter< int >::type NBLMetal(NBLMetalSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type BLMetalSpecsR(BLMetalSpecsRSEXP);
    Rcpp::traits::input_parameter< double >::type CATarget(CATargetSEXP);
    Rcpp::traits::input_parameter< bool >::type DodVidCj(DodVidCjSEXP);
    Rcpp::traits::input_parameter< bool >::type DodVidCjDonnan(DodVidCjDonnanSEXP);
    Rcpp::traits::input_parameter< bool >::type DodKidCj(DodKidCjSEXP);
    Rcpp::traits::input_parameter< bool >::type DoGammai(DoGammaiSEXP);
    Rcpp::traits::input_parameter< bool >::type DoJacDonnan(DoJacDonnanSEXP);
    Rcpp::traits::input_parameter< bool >::type DoJacWHAM(DoJacWHAMSEXP);
    Rcpp::traits::input_parameter< bool >::type DoWHAMSimpleAdjust(DoWHAMSimpleAdjustSEXP);
    Rcpp::traits::input_parameter< bool >::type DoDonnanSimpleAdjust(DoDonnanSimpleAdjustSEXP);
    rcpp_result_gen = Rcpp::wrap(CHESS(QuietFlag, ConvergenceCriteria, MaxIter, NMass, MassName, MassAmt, NComp, CompName, CompType, TotConc, NSpec, SpecName, SpecType, SpecMCR, SpecK, SpecTempKelvin, SpecDeltaH, SpecStoich, SpecCharge, SpecActCorr, DoWHAM, AqueousMCR, WHAMDonnanMCR, HumicSubstGramsPerLiter, WHAMMolWt, WHAMRadius, WHAMP, WHAMDLF, WHAMKZED, SysTempKelvin, DoTox, MetalName, MetalCompR, BLCompR, NBLMetal, BLMetalSpecsR, CATarget, DodVidCj, DodVidCjDonnan, DodKidCj, DoGammai, DoJacDonnan, DoJacWHAM, DoWHAMSimpleAdjust, DoDonnanSimpleAdjust));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BLMEngineInR_CHESS", (DL_FUNC) &_BLMEngineInR_CHESS, 45},
    {NULL, NULL, 0}
};

RcppExport void R_init_BLMEngineInR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
