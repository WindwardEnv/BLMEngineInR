#' @title Calculate the Jacobian matrix
#' @description Calculate the Jacobian matrix of the speciation problem
#'
#' @details This function calculates the Jacobian matrix of the speciation
#'   problem. The Jacobian matrix is an (`NComp` x `NComp`) matrix of
#'   derivatives, where each row gives the derivative of the residual each total
#'   component concentration, dR[i], with respect to, in each column, the
#'   component free concentration, dX[j]. (i.e., \eqn{Z_{i,j} = \delta Resid_{i}
#'   / \delta CompConc_{j}}).
#'
#' @param NComp integer, the number of components
#' @param NSpec integer, the number of species
#' @param SpecStoich integer matrix (NSpec x NComp), the stoichiometry of each
#'   reaction
#' @param SpecConc numeric vector (NSpec), the concentrations of each species
#' @param SpecCtoM numeric vector (NSpec), the factor to apply to concentrations
#'   (i.e., units of mol/L or mol/kg) to convert to masses (i.e., units of mol).
#' @param CompName character vector (NComp), the component names
#' @param MetalComp integer, the position in component vectors of the toxic
#'   metal component
#' @param NBLMetal inteher, the number of BL-Metal species
#' @param BLMetalSpecs integer vector (NBLMetal), the position in the species
#'   vectors of the metal-biotic ligand species associated with toxicity
#' @param DoTox logical, if TRUE = toxicity run, FALSE = speciation run
#'
#' @return numeric matrix (NComp x NComp), the jacobian matrix (see details)
#' @export
RJacobian = function(NComp, NSpec, SpecStoich, SpecConc, SpecCtoM, CompName,
                     MetalComp, NBLMetal, BLMetalSpecs, DoTox) {

  # output:
  # variables:
  #   Sum
  #   iComp1, iComp2, iSpec


  # Doing it the R way...
  SpecName = rownames(SpecStoich)
  JacobianMatrix2 = matrix(data = 0, nrow = NComp, ncol = NComp,
                           dimnames = list(CompName, CompName))
  SpecMoles = SpecConc * SpecCtoM
  SpecMolesByComp = matrix(SpecMoles, nrow = NSpec, ncol = NComp, byrow = F,
                           dimnames = list(SpecName, CompName)) * SpecStoich
  CompStoichDivMoles = SpecStoich *
    matrix(1 / SpecMoles[1:NComp], nrow = NSpec, ncol = NComp, byrow = T,
           dimnames = list(SpecName, CompName))
  JacobianMatrix2 = t(SpecMolesByComp) %*% CompStoichDivMoles
  if (DoTox){
    iComp1 = MetalComp
    for (iComp2 in 1:NComp) {
      Sum = 0
      for (iSpec in BLMetalSpecs){
        Sum = Sum + (SpecStoich[iSpec, iComp2] * SpecStoich[iSpec, iComp1] *
                       SpecConc[iSpec] * SpecCtoM[iSpec])
      }
      if (SpecConc[iComp2] != 0) {
        JacobianMatrix2[iComp1, iComp2] = Sum / (SpecConc[iComp2] * SpecCtoM[iComp2])
      }
    }

  }
  return(JacobianMatrix2)

  # # Doing it the loopy way
  # JacobianMatrix = matrix(data = 0, nrow = NComp, ncol = NComp,
  #                         dimnames = list(CompName, CompName))
  # for (iComp1 in 1:NComp) {
  #   for (iComp2 in 1:NComp) {
  #     Sum = 0
  #     if (DoTox && (iComp1 == MetalComp)){
  #       for (iSpec in BLMetalSpecs){
  #         Sum = Sum + (SpecStoich[iSpec, iComp2] * SpecStoich[iSpec, iComp1] *
  #                        SpecConc[iSpec] * SpecCtoM[iSpec])
  #       }
  #     } else {
  #       for (iSpec in 1:NSpec) {
  #         Sum = Sum + (SpecStoich[iSpec, iComp2] * SpecStoich[iSpec, iComp1] *
  #                        SpecConc[iSpec] * SpecCtoM[iSpec])
  #       }
  #     }
  #     if (SpecConc[iComp2] != 0) {
  #       JacobianMatrix[iComp1, iComp2] = Sum / (SpecConc[iComp2] * SpecCtoM[iComp2])
  #     }
  #   }
  # }
  # return(JacobianMatrix)
  #
  # # Comparing the two
  # JacobianMatrix2 - JacobianMatrix
  # JacobianMatrix2 / JacobianMatrix

}
