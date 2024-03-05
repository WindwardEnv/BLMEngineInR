#' @title Calculate the Residual
#'
#' @description Calculate the residuals of the speciation problem.
#'
#' @details The residuals of the speciation problem are the calculated total
#'   moles for each component minus their known concentrations. In the case of a
#'   toxicity run (DoTox = TRUE), the residual for the `MetalComp` component is
#'   instead the sum of the `BLMetalSpecs` concentrations minus the `CATarget`.
#'   This, paired with a modified set of derivatives for the `MetalComp` in the
#'   Jacobian matrix, results in Newton-Rapshon changing the metal concentration
#'   in order to produce chemical conditions that would result in the critical
#'   accumulation known to cause a toxic effect.
#'
#' @param NComp integer, the number of components
#' @param NSpec integer, the number of species
#' @param SpecConc numeric vector (NSpec), the concentration of chemical species
#' @param SpecStoich integer matrix (NSpec x NComp), the stoichiometry of
#'   species formatin reactions
#' @param TotMoles numeric vector (NComp), the total moles of each component
#' @param SpecCtoM numeric vector (NSpec), the concentration-to-mass conversion
#'   factor for each chemical species
#' @param CompName character vector (NComp), the names of the components
#' @param CompType character vector (NComp), the type of component. It should be
#'   a fixed set of values (MassBal, FixedConc, Substituted, ChargeBal, SurfPot)
#' @param MetalComp integer, the position in component vectors of the toxic
#'   metal component
#' @param BLMetalSpecs integer vector, the position in the species vectors of
#'   the metal-biotic ligand species associated with toxicity
#' @param CATarget numeric, the target critical accumulation value for a
#'   toxicity run, in mol/kg
#' @param DoTox logical, if TRUE = toxicity run, FALSE = speciation run
#'
#' @return A `list` object with the following components:
#' \describe{
#'  \item{\code{MaxError}}{numeric, the highest absolute error fraction in this
#'   iteration =max(abs(Resid / TotMoles))}
#'  \item{\code{WhichMax}}{integer, the position in the component vectors of the
#'    component with the highest absolute error}
#'  \item{\code{Resid}}{numeric vector (NComp), the residuals =
#'    calculated totals - known totals}
#'  \item{\code{CompError}}{numeric vector (NComp), the absolute error fraction
#'    for each component in this iteration =abs(Resid / TotMoles)}
#'  \item{\code{CalcTotConc}}{numeric vector (NComp), the calculated total
#'    concentration of each component = CalcTotMoles / CtoM}
#'  \item{\code{CalcTotMoles}}{numeric vector (NComp), the calculated total
#'    moles of each component = sum(SpecConc * SpecCtoM * SpecStoich[,j]) for
#'    each component j}
#' }
#'
#' @export
RCalcResidual = function(NComp, NSpec, SpecConc, SpecStoich, TotMoles, SpecCtoM,
                         CompName, CompType,
                         MetalComp, BLMetalSpecs, CATarget, DoTox){
  # outputs
  #   MaxError - double - maximum of absolute ratios of residuals to totals
  #   WhichMax - which component has the highest absolute error
  #   Resid - vector(NComp) (maybe...useful for debugging)
  Resid = array(dim=NComp, dimnames = list(CompName))
  #   CompError - the absolute ratios of residuals to totals
  #   CalcTotConc - the calculated total concentrations (mol/L)
  #   CalcTotMoles - the calculated total concentrations (mol)
  # variables:
  #   double ThisError

  CalcTotMoles = array((SpecConc * SpecCtoM) %*% SpecStoich, dim = NComp, dimnames = list(CompName))
  # CalcTotConc = as.numeric(matrix(SpecConc, nrow = 1, ncol = NSpec) %*% SpecStoich)
  Resid = CalcTotMoles - TotMoles # * SpecCtoM[1:NComp]
  Resid[CompType == "FixedConc"] = 0.0
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
    MaxError = MaxError,#needed - double
    WhichMax = WhichMax,#could be int or string
    Resid = Resid,#needed - vector of double
    CompError = ThisError,
    CalcTotConc = CalcTotConc,#needed - vector of double
    CalcTotMoles = CalcTotMoles
  ))
}
