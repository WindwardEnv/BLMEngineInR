#' @title Calculate the Residual
#'
#' @description Calculate the residuals of the speciation problem.
#'
#' @details The residuals of the speciation problem are the calculated total concentrations \eqn{\Sum()}
#'
#'
#' @param NComp
#' @param NSpec
#' @param SpecConc
#' @param SpecStoich
#' @param TotMoles
#' @param SpecCtoM
#' @param CompName
#' @param CompType
#' @param MetalComp
#' @param BLMetalSpecs
#' @param CATarget
#' @param DoTox
#'
#' @return
#' @export
#'
#' @examples
RCalcResidual = function(NComp, NSpec, SpecConc, SpecStoich, TotMoles, SpecCtoM,
                         CompName, CompType, #CompSiteDens,
                         MetalComp, BLMetalSpecs, CATarget, DoTox){
  # inputs:
  #   NComp - number of components
  #   NSpec - number of species
  #   SpecConc - species free concentrations (mol/L or mol/gww)
  #   SpecStoich - species stoichiometry
  #   TotConc - component total concentrations (mol)
  #   SpecCtoM - mass compartment concentration to mass conversion for each species (L or gww)
  #   CompName - component names
  #   CompType - component types
  #   MetalComp - integer, position of the metal component in the component list
  #   BLMetalSpecs - integer vector, positon of the metal-bound biotic ligand in the species list
  #   CATarget - the target critical accumulation for toxicity runs
  #   DoTox - boolean, TRUE means this is a toxicity run, FALSE means speciation run
  # outputs
  #   Resid - vector(NComp) (maybe...useful for debugging)
  Resid = array(dim=NComp, dimnames = list(CompName))
  #   MaxError - double - maximum of absolute ratios of residuals to totals
  #   WhichMax - which component has the highest absolute error
  #   CalcTotConc - the calculated total concentrations (mol)
  # variables:
  #   double CalcTotConc
  #   double ThisError

  CalcTotMoles = array((SpecConc * SpecCtoM) %*% SpecStoich, dim = NComp, dimnames = list(CompName))
  # CalcTotConc = as.numeric(matrix(SpecConc, nrow = 1, ncol = NSpec) %*% SpecStoich)
  Resid = CalcTotMoles - TotMoles # * SpecCtoM[1:NComp]
  Resid[CompType == "FixedAct"] = 0.0
  ThisError = abs(Resid / TotMoles)
  if(DoTox){
    Resid[MetalComp] = sum(SpecConc[BLMetalSpecs]) - CATarget
    # sum( (mol/L or mol/kgww) ) - mol = sum(mol) - mol = mol
    ThisError[MetalComp] = abs(Resid[MetalComp] / CATarget)
  }
  MaxError = max(ThisError)
  WhichMax = which.max(ThisError)

  CalcTotConc = CalcTotMoles / SpecCtoM[1:NComp]

  return(list(
    MaxError = MaxError,
    WhichMax = WhichMax,
    Resid = Resid,
    CompError = ThisError,
    CalcTotConc = CalcTotConc,
    CalcTotMoles = CalcTotMoles
  ))
}
