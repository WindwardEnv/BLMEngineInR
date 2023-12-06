#include <Rcpp.h>
#include "CHESSFunctions.h"

//' @title Calculate the Jacobian matrix
//' @description Calculate the Jacobian matrix of the speciation problem
//'
//' @details This function calculates the Jacobian matrix of the speciation
//'   problem. The Jacobian matrix is an (`NComp` x `NComp`) matrix of
//'   derivatives, where each row gives the derivative of the residual each total
//'   component concentration, dR[i], with respect to, in each column, the
//'   component free concentration, dX[j]. (i.e., \eqn{Z_{i,j} = \delta Resid_{i}
//'   / \delta CompConc_{j}}).
//'
//' @param NComp integer, the number of components
//' @param NSpec integer, the number of species
//' @param SpecStoich integer matrix (NSpec x NComp), the stoichiometry of
//'   each reaction
//' @param SpecConc numeric vector (NSpec), the concentrations of each species
//' @param SpecCtoM numeric vector (NSpec), the factor to apply to
//'   concentrations (i.e., units of mol/L or mol/kg) to convert to masses (i.e.,
//'   units of mol).
//' @param SpecName character vector (NSpec), the species names
//' @param MetalComp integer, the position in component vectors of the toxic
//'   metal component
//' @param BLMetalSpecs integer vector, the position in the species vectors of
//'   the metal-biotic ligand species associated with toxicity
//' @param DoTox logical, if TRUE = toxicity run, FALSE = speciation run
//'
//' @return numeric matrix (NComp x NComp), the jacobian matrix (see details)
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix Jacobian (unsigned int NComp, //number of components
                              unsigned int NSpec, //number of species
                              Rcpp::IntegerMatrix SpecStoich, //formation reaction stoichiometry (NSpec x NComp)
                              Rcpp::NumericVector SpecConc, //species concentrations
                              Rcpp::NumericVector SpecCtoM, //concentration to mass conversion for each species
                              Rcpp::CharacterVector CompName, //names of components
                              unsigned int MetalComp, //position of the metal component
                              unsigned int NBLMetal, //number of BL-Metal species
                              Rcpp::IntegerVector BLMetalSpecs, //positions of BL-metal species
                              bool DoTox) {//TRUE = do toxicity, FALSE = speciation only
  /* outputs: */
  Rcpp::NumericMatrix JacobianMatrix(NComp);
  Rcpp::rownames(JacobianMatrix) = CompName;
  Rcpp::colnames(JacobianMatrix) = CompName;

  /* variables: */
  double Sum;
  unsigned int iComp1, iComp2, iSpec, i;

  /* Loop through the Jacobian */
  for (iComp1 = 0; iComp1 < NComp; iComp1++) {
    for (iComp2 = 0; iComp2 < NComp; iComp2++) {
      Sum = 0;
      if (DoTox && (iComp1 == (MetalComp - 1))){//subtract 1 because C++ is 0-based
        /* Toxicity mode the metal's derivatives are relative to the CA error */
        for (i = 0; i < NBLMetal; i++) {
          iSpec = BLMetalSpecs(i) - 1;
          Sum += (SpecStoich(iSpec, iComp2) * SpecStoich(iSpec, iComp1) *
            SpecConc(iSpec) * SpecCtoM(iSpec));
        };//NEXT iSpec
      } else {
        /* All others are based on Resid = CalcTotMoles - TotMoles */
        for (iSpec = 0; iSpec < NSpec; iSpec++) {
          Sum += (SpecStoich(iSpec, iComp2) * SpecStoich(iSpec, iComp1) *
            SpecConc(iSpec) * SpecCtoM(iSpec));
        };//NEXT iSpec
      }
      if (SpecConc(iComp2) != 0) {
        JacobianMatrix(iComp1, iComp2) = Sum / (SpecConc(iComp2) * SpecCtoM(iComp2));
      }
    };//NEXT iComp2
  };//NEXT iComp1

  return JacobianMatrix;
}
