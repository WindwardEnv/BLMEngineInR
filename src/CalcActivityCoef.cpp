#include <Rcpp.h>
#include "CHESSFunctions.h"

//' @title ExtDebyeHuckel
//' 
//' @description Calculates the Extended Debye-Huckel Activity Coefficients
//' 
//' @details This function will calculate the EDH activitiy coefficients for 
//'   ions with charges from 0 to MaxCharge (input). 
//' 
//' @author Kelly Croteau (kellyc@windwardenv.com)
//' 
//' @param IonicStrength double, the ionic strength of the solution
//' @param MaxCharge int, the maximum absolute charge to calculate activity
//'   coefficients for (e.g., if the ions in solutions have charges of -3, -2,
//'   -1, 0, +1, +2, then MaxCharge should be 3).
//' @param SysTempKelvin double, the temperature of the solution, in Kelvin
//' 
//' @return Rcpp::NumericVector of length MaxCharge + 1, with the activity 
//'   coefficient for a 0-charge ion in the 0th element, +/-1 charge ion in the
//'   1st element, etc.
//' 
//' @keywords internal
//' 
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

//' @title Davies
//' 
//' @description Calculates the Davies Activity Coefficients
//' 
//' @details This function will calculate the Davies activitiy coefficients for 
//'   ions with charges from 0 to MaxCharge (input).
//' 
//' @author Kelly Croteau (kellyc@windwardenv.com)
//' 
//' @param IonicStrength double, the ionic strength of the solution
//' @param MaxCharge int, the maximum absolute charge to calculate activity
//'   coefficients for (e.g., if the ions in solutions have charges of -3, -2,
//'   -1, 0, +1, +2, then MaxCharge should be 3).
//' 
//' @return Rcpp::NumericVector of length MaxCharge + 1, with the activity 
//'   coefficient for a 0-charge ion in the 0th element, +/-1 charge ion in the
//'   1st element, etc.
//' 
//' @keywords internal
//' 
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


//' @title CalcActivityCoef
//' 
//' @description {text}
//' 
//' @details {text}
//' 
//' @author Kelly Croteau (kellyc@windwardenv.com)
//' 
//' @param NSpec int, the number of chemical species for which we have
//'   formation reactions in the simulation
//' @param SpecName Rcpp::CharacterVector, the name of the chemical species
//'   for which we have formation reactions
//' @param SpecActCorr Rcpp::CharacterVector, (should be length NSpec) the
//'   activity correction method of the chemical species for which we have
//'   formation reactions
//' @param SpecCharge Rcpp::IntegerVector, (should be length NSpec) the ionic 
//'   charges of each species
//' @param IonicStrength double, the ionic strength of the solution
//' @param SysTempKelvin double, the temperature of the solution, in Kelvin
//' 
//' @return Rcpp::NumericVector of length NSpec with the activity coefficient
//'   for each chemical species.
//'
//' @family CHESSFunctions
//' 
Rcpp::NumericVector CalcActivityCoef(int NSpec,
                                     Rcpp::CharacterVector SpecName,
                                     Rcpp::CharacterVector SpecActCorr,
                                     Rcpp::IntegerVector SpecCharge,
                                     double IonicStrength,
                                     double SysTempKelvin) {

  /* outputs */
  Rcpp::NumericVector SpecActivityCoef(NSpec);
    SpecActivityCoef.names() = SpecName;

  /* variables */
  int iSpec;
  int MaxCharge = Rcpp::max(Rcpp::abs(SpecCharge));
  Rcpp::NumericVector ActivityCoefDebye;
  Rcpp::NumericVector ActivityCoefDavies;
  Rcpp::String tmp;
  
  ActivityCoefDebye = ExtDebyeHuckel(IonicStrength, MaxCharge, SysTempKelvin);
  //Rcpp::Rcout << "ActivityCoefDebye = " << ActivityCoefDebye << std::endl;
  ActivityCoefDavies = Davies(IonicStrength, MaxCharge);
  
  for (iSpec = 0; iSpec < NSpec; iSpec++) {
    if (SpecActCorr(iSpec) == ACTYPE_DEBYE) {
      SpecActivityCoef(iSpec) = ActivityCoefDebye(abs(SpecCharge(iSpec)));
    } else if (SpecActCorr(iSpec) == ACTYPE_DAVIES) {
      SpecActivityCoef(iSpec) = ActivityCoefDavies(abs(SpecCharge(iSpec)));
    } else {//if (SpecActCorr(iSpec) == ACTYPE_NONE) {
      SpecActivityCoef(iSpec) = 1.0;
    }
  }

  return SpecActivityCoef;

}
