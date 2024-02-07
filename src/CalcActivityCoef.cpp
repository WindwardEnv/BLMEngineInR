#include <Rcpp.h>
#include "CHESSFunctions.h"

Rcpp::NumericVector ExtDebyeHuckel(double IonicStrength, 
                                   int MaxCharge, 
                                   double SysTempKelvin) {

  /* outputs */
  Rcpp::NumericVector ActivityCoefDebye(MaxCharge + 1);

  /* variables */
  int AbsCharge;
  double A, B, C, D;
  int F;

  A = 0.27 + (0.0008 * SysTempKelvin);
  B = 0.33;
  C = -A * sqrt(IonicStrength);
  D = B * sqrt(IonicStrength);

  for (AbsCharge = 0; AbsCharge <= MaxCharge; AbsCharge++) {
    //E = abs(AbsCharge);//redundant
    F = pow(AbsCharge, 2);
    ActivityCoefDebye(AbsCharge) = pow(10, F * C / (1 + (3 * AbsCharge * D)));
  }
  
/* This is what's done in the PB code...Is this even necessary?? Doesn't the 
   math work out to these same equations if you just use the last one? 
  switch(E) {
    case 0: SpecActivityCoef = 1;
    case 1: SpecActivityCoef = pow(10, C / (1 + (3 * D)));
    case 2: SpecActivityCoef = pow(10, 4 * C / (1 + (6 * D)));
    case 3: SpecActivityCoef = pow(10, 9 * C / (1 + (9 * D)));
    case 4: SpecActivityCoef = pow(10, 16 * C / (1 + (12 * D)));
    default : SpecActivityCoef = pow(10, F * C / (1 + (3 * E * D)));
  }
*/
  return ActivityCoefDebye;
}

Rcpp::NumericVector Davies(double IonicStrength, int MaxCharge) {
  
  /* outputs */
  Rcpp::NumericVector ActivityCoefDavies(MaxCharge + 1);
  
  /* variables */
  int AbsCharge;
  double SqrtIonicStrength = sqrt(IonicStrength);
  double DaviesConst;
  double Tmp;

  // This part is the same for every charge...
  DaviesConst = -0.5 * ((SqrtIonicStrength / (1 + SqrtIonicStrength)) 
                        - 0.2 * IonicStrength);

  for (AbsCharge = 0; AbsCharge <= MaxCharge; AbsCharge++) {
    if (AbsCharge == 0) { 
      ActivityCoefDavies(AbsCharge) = 1;
    } else {
      Tmp = pow(AbsCharge, 2) * DaviesConst;
      if (Tmp > 3.0) {
        ActivityCoefDavies(AbsCharge) = 1;
      } else {
        ActivityCoefDavies(AbsCharge) = pow(10, Tmp);
      }
    }
  }
  
  return ActivityCoefDavies;

}


Rcpp::NumericVector CalcActivityCoef(unsigned int NSpec,
                                     Rcpp::CharacterVector SpecActCorr,
                                     Rcpp::IntegerVector SpecCharge,
                                     double IonicStrength,
                                     double SysTempKelvin) {

  /* outputs */
  Rcpp::NumericVector SpecActivityCoef(NSpec);

  /* variables */
  unsigned int iSpec;
  int MaxCharge = Rcpp::max(Rcpp::abs(SpecCharge));
  Rcpp::NumericVector ActivityCoefDebye;
  Rcpp::NumericVector ActivityCoefDavies;
  
  ActivityCoefDebye = ExtDebyeHuckel(IonicStrength, MaxCharge, SysTempKelvin);
  //Rcpp::Rcout << "ActivityCoefDebye = " << ActivityCoefDebye << std::endl;
  ActivityCoefDavies = Davies(IonicStrength, MaxCharge);
  
  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    if (SpecActCorr(iSpec) == "Debye") {
      SpecActivityCoef(iSpec) = ActivityCoefDebye(abs(SpecCharge(iSpec)));
    } else if (SpecActCorr(iSpec) == "Davies") {
      SpecActivityCoef(iSpec) = ActivityCoefDavies(abs(SpecCharge(iSpec)));
    } else {//if (SpecActCorr(iSpec) == "None") {
      SpecActivityCoef(iSpec) = 1.0;
    }
  }

  return SpecActivityCoef;

}
