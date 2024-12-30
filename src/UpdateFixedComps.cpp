#include <Rcpp.h>
#include "CHESSFunctions.h"

//' @title Update fixed activity/concentration components
//' 
//' @description Components which are fixed activity or fixed concentration 
//'   should maintain an activity or concentration equal to their input 
//'   concentration. This function resets the SpecConc and CompConc vectors to
//'   those values.
//' 
//' @author Kelly Croteau (kellyc@windwardenv.com)
//' 
//' @param NComp integer, the combined number of components in the simulation,
//'   including the input components, defined components (and including the
//'   defined components that get added by ExpandWHAM)
//' @param CompType character vector (NComp), the type of each component in the
//'   simulation
//' @param InCompConc numeric vector (NComp) of the input value of each
//'   component's concentration, whether it's on the basis of total
//'   concentration, fixed activity, or fixed concentration
//' @param SpecActivityCoef {type}, {text}
//' @param SpecConc numeric vector (NSpec), the concentrations of each species 
//'   for which we have formation reactions
//' @param CompConc numeric vector (NComp) of component concentrations
//' 
//' 
void UpdateFixedComps(int NComp, 
                      Rcpp::CharacterVector CompType, 
                      Rcpp::NumericVector InCompConc, 
                      Rcpp::NumericVector SpecActivityCoef,
                      Rcpp::NumericVector &SpecConc, 
                      Rcpp::NumericVector &CompConc) {
  /* variables */
  int iComp;

  for(iComp = 0; iComp < NComp; iComp++) {
    if (CompType[iComp] == CTYPE_FIXEDACT) {
      SpecConc[iComp] = InCompConc[iComp] / SpecActivityCoef[iComp];
      CompConc[iComp] = SpecConc[iComp];
    } else if (CompType[iComp] == CTYPE_FIXEDCONC){
      SpecConc[iComp] = InCompConc[iComp];
      CompConc[iComp] = SpecConc[iComp];
    }
  }
  
}
