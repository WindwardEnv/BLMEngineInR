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

// [[Rcpp::export]]
Rcpp::NumericMatrix Jacobian (unsigned int NComp, //number of components
                              unsigned int NSpec, //number of species
                              Rcpp::IntegerMatrix SpecStoich, //formation reaction stoichiometry (NSpec x NComp)
                              Rcpp::NumericVector SpecConc, //species concentrations
                              Rcpp::NumericVector CtoM) {//concentration to mass conversion for each species
  // //inputs:
  // //  NComp
  // //  NSpec
  // //  SpecStoich
  // //  SpecConc
  // //  CtoM
  // //outputs:
  Rcpp::NumericMatrix Z(NComp, NComp);
  // //variables:
  // double Sum;
  // unsigned int iComp1, iComp2, iSpec;
  //
  // for (iComp1 = 0; iComp1 < NComp; iComp1++) {
  //   for (iComp2 = 0; iComp2 < NComp; iComp2++) {
  //     Sum = 0;
  //     for (iSpec = 0; iSpec < NSpec; iSpec++) {
  //       Sum =+ SpecStoich(iSpec, iComp2) * SpecStoich(iSpec, iComp1) * SpecConc[iSpec] * CtoM[iSpec];
  //     };//NEXT iSpec
  //     if (SpecConc(iComp2) != 0) {
  //       Z(iComp1, iComp2) = Sum / SpecConc[iComp2];
  //     }
  //   };//NEXT iComp2
  // };//NEXT iComp1
  //
  return Z;
}


// // You can include R code blocks in C++ files processed with sourceCpp
// // (useful for testing and development). The R code will be automatically
// // run after the compilation.
// //
//
// /*** R
// timesTwo(42)
// */
