#' Calculate Species Concentrations
#'
#' Use `RCalcSpecConc` to calculate species concentrations from a known set of
#' free component ion concentrations.
#'
#' This is an internal function that will, for each species `i` in `NSpec`
#' calculate the equilibrium concentration, which in its most basic form would
#' be calculated as \eqn{SpecConc_{i} = K_{i} *
#' \prod_{j=1}^{n}(CompConc^SpecStoich_{i,j})}. Further modifications to these
#' calculations would include temperature correction, ionic strength
#' corrections, and diffuse double layer calculations for organic matter
#' binding. As of 19-Oct-2023, none of these corrections have been implemented.
#'
#' This is the base-R version of this function.
#'
#' @param CompConc A vector of component concentrations for each of `NComp`
#'   components.
#' @param SpecK A vector of reaction equilibrium constants for each of `NSpec`
#'   reactions.
#' @param SpecStoich A matrix of reaction stoichiometry, with `NSpec` rows and
#'   `NComp` columns.
#' @param NComp The number of components in the equilibrium system.
#' @param NSpec The number of species (reactions) in the equilibrium system.
#' @param SpecName character vector (NSpec), the name of the chemical species
#'   for which we have formation reactions
#'
#' @returns A vector of `NSpec` species concentrations.
#'
RCalcSpecConc = function(CompConc, SpecK, SpecStoich, NComp = length(CompConc),
                         NSpec = length(SpecK), SpecName = names(SpecK)){
  SpecConc = rep(1, NSpec)
  for (iSpec in 1:NSpec){
    CompConcStep = CompConc ^ SpecStoich[iSpec,]
    SpecConc[iSpec] = prod(CompConcStep) * SpecK[iSpec]
  }
  names(SpecConc) = SpecName
  return(SpecConc)
}

