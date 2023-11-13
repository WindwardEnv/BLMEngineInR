// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// CppCalcSpecConc
std::vector<double> CppCalcSpecConc(std::vector<double> CompConc, std::vector<double> SpecK, Rcpp::IntegerMatrix SpecStoich, unsigned int NComp, unsigned int NSpec);
RcppExport SEXP _BLMEngineInR_CppCalcSpecConc(SEXP CompConcSEXP, SEXP SpecKSEXP, SEXP SpecStoichSEXP, SEXP NCompSEXP, SEXP NSpecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type CompConc(CompConcSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type SpecK(SpecKSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type SpecStoich(SpecStoichSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NComp(NCompSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NSpec(NSpecSEXP);
    rcpp_result_gen = Rcpp::wrap(CppCalcSpecConc(CompConc, SpecK, SpecStoich, NComp, NSpec));
    return rcpp_result_gen;
END_RCPP
}
// CppCalcLogSpecConc
Rcpp::NumericVector CppCalcLogSpecConc(Rcpp::NumericVector LogCompConc, Rcpp::NumericVector SpecLogK, Rcpp::IntegerMatrix SpecStoich, unsigned int NComp, unsigned int NSpec);
RcppExport SEXP _BLMEngineInR_CppCalcLogSpecConc(SEXP LogCompConcSEXP, SEXP SpecLogKSEXP, SEXP SpecStoichSEXP, SEXP NCompSEXP, SEXP NSpecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type LogCompConc(LogCompConcSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SpecLogK(SpecLogKSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type SpecStoich(SpecStoichSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NComp(NCompSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NSpec(NSpecSEXP);
    rcpp_result_gen = Rcpp::wrap(CppCalcLogSpecConc(LogCompConc, SpecLogK, SpecStoich, NComp, NSpec));
    return rcpp_result_gen;
END_RCPP
}
// Jacobian
Rcpp::NumericMatrix Jacobian(unsigned int NComp, unsigned int NSpec, Rcpp::IntegerMatrix SpecStoich, Rcpp::NumericVector SpecConc, Rcpp::NumericVector CtoM);
RcppExport SEXP _BLMEngineInR_Jacobian(SEXP NCompSEXP, SEXP NSpecSEXP, SEXP SpecStoichSEXP, SEXP SpecConcSEXP, SEXP CtoMSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type NComp(NCompSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NSpec(NSpecSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type SpecStoich(SpecStoichSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type SpecConc(SpecConcSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type CtoM(CtoMSEXP);
    rcpp_result_gen = Rcpp::wrap(Jacobian(NComp, NSpec, SpecStoich, SpecConc, CtoM));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BLMEngineInR_CppCalcSpecConc", (DL_FUNC) &_BLMEngineInR_CppCalcSpecConc, 5},
    {"_BLMEngineInR_CppCalcLogSpecConc", (DL_FUNC) &_BLMEngineInR_CppCalcLogSpecConc, 5},
    {"_BLMEngineInR_Jacobian", (DL_FUNC) &_BLMEngineInR_Jacobian, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_BLMEngineInR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
