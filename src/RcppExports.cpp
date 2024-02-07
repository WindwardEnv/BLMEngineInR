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
Rcpp::List CHESS(Rcpp::String QuietFlag, double ConvergenceCriteria, unsigned int MaxIter, unsigned int NMass, unsigned int NComp, unsigned int NSpec, unsigned int NBLMetal, Rcpp::IntegerVector SpecMC, Rcpp::NumericVector SpecK, Rcpp::NumericVector SpecTempKelvin, Rcpp::NumericVector SpecDeltaH, Rcpp::IntegerMatrix SpecStoich, Rcpp::IntegerVector SpecCharge, Rcpp::CharacterVector SpecActCorr, Rcpp::NumericVector SpecCtoM, Rcpp::CharacterVector SpecName, Rcpp::CharacterVector CompType, Rcpp::CharacterVector CompName, Rcpp::NumericVector TotConc, bool DoWHAM, int AqueousMC, //Rcpp::IntegerVector DonnanMC, Rcpp::NumericVector SolHS, Rcpp::NumericVector wMolWt, Rcpp::NumericVector wRadius, Rcpp::NumericVector wP, double wDLF, double wKZED, double SysTempKelvin, bool DoTox, Rcpp::String MetalName, unsigned int MetalComp, Rcpp::IntegerVector BLMetalSpecs, double CATarget);
RcppExport SEXP _BLMEngineInR_CHESS(SEXP QuietFlagSEXP, SEXP ConvergenceCriteriaSEXP, SEXP MaxIterSEXP, SEXP NMassSEXP, SEXP NCompSEXP, SEXP NSpecSEXP, SEXP NBLMetalSEXP, SEXP SpecMCSEXP, SEXP SpecKSEXP, SEXP SpecTempKelvinSEXP, SEXP SpecDeltaHSEXP, SEXP SpecStoichSEXP, SEXP SpecChargeSEXP, SEXP SpecActCorrSEXP, SEXP SpecCtoMSEXP, SEXP SpecNameSEXP, SEXP CompTypeSEXP, SEXP CompNameSEXP, SEXP TotConcSEXP, SEXP DoWHAMSEXP, SEXP AqueousMCSEXP, SEXP DonnanMCSEXP, SEXP SolHSSEXP, SEXP wMolWtSEXP, SEXP wRadiusSEXP, SEXP wPSEXP, SEXP wDLFSEXP, SEXP wKZEDSEXP, SEXP SysTempKelvinSEXP, SEXP DoToxSEXP, SEXP MetalNameSEXP, SEXP MetalCompSEXP, SEXP BLMetalSpecsSEXP, SEXP CATargetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type QuietFlag(QuietFlagSEXP);
    Rcpp::traits::input_parameter< double >::type ConvergenceCriteria(ConvergenceCriteriaSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type MaxIter(MaxIterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NMass(NMassSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NComp(NCompSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NSpec(NSpecSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NBLMetal(NBLMetalSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type SpecMC(SpecMCSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SpecK(SpecKSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SpecTempKelvin(SpecTempKelvinSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SpecDeltaH(SpecDeltaHSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type SpecStoich(SpecStoichSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type SpecCharge(SpecChargeSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type SpecActCorr(SpecActCorrSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SpecCtoM(SpecCtoMSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type SpecName(SpecNameSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type CompType(CompTypeSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type CompName(CompNameSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type TotConc(TotConcSEXP);
    Rcpp::traits::input_parameter< bool >::type DoWHAM(DoWHAMSEXP);
    Rcpp::traits::input_parameter< int >::type AqueousMC(AqueousMCSEXP);
    Rcpp::traits::input_parameter< //Rcpp::IntegerVector >::type DonnanMC(DonnanMCSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SolHS(SolHSSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type wMolWt(wMolWtSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type wRadius(wRadiusSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type wP(wPSEXP);
    Rcpp::traits::input_parameter< double >::type wDLF(wDLFSEXP);
    Rcpp::traits::input_parameter< double >::type wKZED(wKZEDSEXP);
    Rcpp::traits::input_parameter< double >::type SysTempKelvin(SysTempKelvinSEXP);
    Rcpp::traits::input_parameter< bool >::type DoTox(DoToxSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type MetalName(MetalNameSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type MetalComp(MetalCompSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type BLMetalSpecs(BLMetalSpecsSEXP);
    Rcpp::traits::input_parameter< double >::type CATarget(CATargetSEXP);
    rcpp_result_gen = Rcpp::wrap(CHESS(QuietFlag, ConvergenceCriteria, MaxIter, NMass, NComp, NSpec, NBLMetal, SpecMC, SpecK, SpecTempKelvin, SpecDeltaH, SpecStoich, SpecCharge, SpecActCorr, SpecCtoM, SpecName, CompType, CompName, TotConc, DoWHAM, AqueousMC, DonnanMC, SolHS, wMolWt, wRadius, wP, wDLF, wKZED, SysTempKelvin, DoTox, MetalName, MetalComp, BLMetalSpecs, CATarget));
    return rcpp_result_gen;
END_RCPP
}
// CalcResidualList
Rcpp::List CalcResidualList(unsigned int NComp, unsigned int NSpec, Rcpp::NumericVector SpecConc, Rcpp::IntegerMatrix SpecStoich, Rcpp::NumericVector TotMoles, Rcpp::NumericVector SpecCtoM, Rcpp::CharacterVector CompName, Rcpp::CharacterVector CompType, unsigned int MetalComp, unsigned int NBLMetal, Rcpp::IntegerVector BLMetalSpecs, double CATarget, bool DoTox);
RcppExport SEXP _BLMEngineInR_CalcResidualList(SEXP NCompSEXP, SEXP NSpecSEXP, SEXP SpecConcSEXP, SEXP SpecStoichSEXP, SEXP TotMolesSEXP, SEXP SpecCtoMSEXP, SEXP CompNameSEXP, SEXP CompTypeSEXP, SEXP MetalCompSEXP, SEXP NBLMetalSEXP, SEXP BLMetalSpecsSEXP, SEXP CATargetSEXP, SEXP DoToxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type NComp(NCompSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NSpec(NSpecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SpecConc(SpecConcSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type SpecStoich(SpecStoichSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type TotMoles(TotMolesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SpecCtoM(SpecCtoMSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type CompName(CompNameSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type CompType(CompTypeSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type MetalComp(MetalCompSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NBLMetal(NBLMetalSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type BLMetalSpecs(BLMetalSpecsSEXP);
    Rcpp::traits::input_parameter< double >::type CATarget(CATargetSEXP);
    Rcpp::traits::input_parameter< bool >::type DoTox(DoToxSEXP);
    rcpp_result_gen = Rcpp::wrap(CalcResidualList(NComp, NSpec, SpecConc, SpecStoich, TotMoles, SpecCtoM, CompName, CompType, MetalComp, NBLMetal, BLMetalSpecs, CATarget, DoTox));
    return rcpp_result_gen;
END_RCPP
}
// CalcStep
Rcpp::NumericVector CalcStep(Rcpp::NumericMatrix JacobianMatrix, Rcpp::NumericVector Resid, unsigned int NComp, Rcpp::CharacterVector CompType, Rcpp::CharacterVector CompName);
RcppExport SEXP _BLMEngineInR_CalcStep(SEXP JacobianMatrixSEXP, SEXP ResidSEXP, SEXP NCompSEXP, SEXP CompTypeSEXP, SEXP CompNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type JacobianMatrix(JacobianMatrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Resid(ResidSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NComp(NCompSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type CompType(CompTypeSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type CompName(CompNameSEXP);
    rcpp_result_gen = Rcpp::wrap(CalcStep(JacobianMatrix, Resid, NComp, CompType, CompName));
    return rcpp_result_gen;
END_RCPP
}
// InitialGuess
Rcpp::NumericVector InitialGuess(Rcpp::NumericVector TotConc, Rcpp::CharacterVector CompType, Rcpp::NumericVector SpecK, Rcpp::IntegerMatrix SpecStoich, Rcpp::CharacterVector SpecName, unsigned int NComp, unsigned int NSpec);
RcppExport SEXP _BLMEngineInR_InitialGuess(SEXP TotConcSEXP, SEXP CompTypeSEXP, SEXP SpecKSEXP, SEXP SpecStoichSEXP, SEXP SpecNameSEXP, SEXP NCompSEXP, SEXP NSpecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type TotConc(TotConcSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type CompType(CompTypeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SpecK(SpecKSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type SpecStoich(SpecStoichSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type SpecName(SpecNameSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NComp(NCompSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NSpec(NSpecSEXP);
    rcpp_result_gen = Rcpp::wrap(InitialGuess(TotConc, CompType, SpecK, SpecStoich, SpecName, NComp, NSpec));
    return rcpp_result_gen;
END_RCPP
}
// Jacobian
Rcpp::NumericMatrix Jacobian(unsigned int NComp, unsigned int NSpec, Rcpp::IntegerMatrix SpecStoich, Rcpp::NumericVector SpecConc, Rcpp::NumericVector SpecCtoM, Rcpp::CharacterVector CompName, unsigned int MetalComp, unsigned int NBLMetal, Rcpp::IntegerVector BLMetalSpecs, bool DoTox);
RcppExport SEXP _BLMEngineInR_Jacobian(SEXP NCompSEXP, SEXP NSpecSEXP, SEXP SpecStoichSEXP, SEXP SpecConcSEXP, SEXP SpecCtoMSEXP, SEXP CompNameSEXP, SEXP MetalCompSEXP, SEXP NBLMetalSEXP, SEXP BLMetalSpecsSEXP, SEXP DoToxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type NComp(NCompSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NSpec(NSpecSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type SpecStoich(SpecStoichSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SpecConc(SpecConcSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SpecCtoM(SpecCtoMSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type CompName(CompNameSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type MetalComp(MetalCompSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NBLMetal(NBLMetalSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type BLMetalSpecs(BLMetalSpecsSEXP);
    Rcpp::traits::input_parameter< bool >::type DoTox(DoToxSEXP);
    rcpp_result_gen = Rcpp::wrap(Jacobian(NComp, NSpec, SpecStoich, SpecConc, SpecCtoM, CompName, MetalComp, NBLMetal, BLMetalSpecs, DoTox));
    return rcpp_result_gen;
END_RCPP
}
// UpdateTotalsList
Rcpp::List UpdateTotalsList(unsigned int NComp, unsigned int NSpec, bool DoTox, Rcpp::CharacterVector CompType, Rcpp::CharacterVector CompName, Rcpp::String MetalName, Rcpp::IntegerMatrix SpecStoich, Rcpp::NumericVector SpecMoles, Rcpp::NumericVector CompCtoM, Rcpp::NumericVector TotMoles, Rcpp::NumericVector TotConc);
RcppExport SEXP _BLMEngineInR_UpdateTotalsList(SEXP NCompSEXP, SEXP NSpecSEXP, SEXP DoToxSEXP, SEXP CompTypeSEXP, SEXP CompNameSEXP, SEXP MetalNameSEXP, SEXP SpecStoichSEXP, SEXP SpecMolesSEXP, SEXP CompCtoMSEXP, SEXP TotMolesSEXP, SEXP TotConcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type NComp(NCompSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NSpec(NSpecSEXP);
    Rcpp::traits::input_parameter< bool >::type DoTox(DoToxSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type CompType(CompTypeSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type CompName(CompNameSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type MetalName(MetalNameSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type SpecStoich(SpecStoichSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SpecMoles(SpecMolesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type CompCtoM(CompCtoMSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type TotMoles(TotMolesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type TotConc(TotConcSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateTotalsList(NComp, NSpec, DoTox, CompType, CompName, MetalName, SpecStoich, SpecMoles, CompCtoM, TotMoles, TotConc));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BLMEngineInR_CHESS", (DL_FUNC) &_BLMEngineInR_CHESS, 34},
    {"_BLMEngineInR_CalcResidualList", (DL_FUNC) &_BLMEngineInR_CalcResidualList, 13},
    {"_BLMEngineInR_CalcStep", (DL_FUNC) &_BLMEngineInR_CalcStep, 5},
    {"_BLMEngineInR_InitialGuess", (DL_FUNC) &_BLMEngineInR_InitialGuess, 7},
    {"_BLMEngineInR_Jacobian", (DL_FUNC) &_BLMEngineInR_Jacobian, 10},
    {"_BLMEngineInR_UpdateTotalsList", (DL_FUNC) &_BLMEngineInR_UpdateTotalsList, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_BLMEngineInR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
