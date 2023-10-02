#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::depends(Rcpp)]]
extern "C" {
  FNCalcSpecConc_(int NComp, int NSpec, double *CConc, double *K, int **Stoich, double *SConc);
}

// [[Rcpp::export]]
Rcpp::NumericVector FCalcSpecConc(Rcpp::NumericVector CConc,
                                  Rcpp::NumericVector K,
                                  Rcpp::IntegerMatrix Stoich,
                                  int NComp, int NSpec) {

  Rcpp::NumericVector SConc(NSpec);

  FNCalcSpecConc_(NComp, NSpec, CConc, K, Stoich, SConc)

  return SConc;
}
