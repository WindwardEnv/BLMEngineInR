#' Calculate Species Concentrations
#'
#' Use `RCalcSpecConc` to calculate species concentrations from a known set of
#' free component ion concentrations.
#'
#' This is an internal function that will, for each species `i` in `NSpec`
#' calculate the equilibrium concentration, which in its most basic form would
#' be calculated as \eqn{SConc_{i} = K_{i} *
#' \prod_{j=1}^{n}(CConc^Stoich_{i,j})}. Further modifications to these
#' calculations would include temperature correction, ionic strength
#' corrections, and diffuse double layer calculations for organic matter
#' binding. As of 19-Oct-2023, none of these corrections have been implemented.
#'
#' This is the base-R version of this function.
#'
#' @param CConc A vector of component concentrations for each of `NComp`
#'   components.
#' @param K A vector of reaction equilibrium constants for each of `NSpec`
#'   reactions.
#' @param Stoich A matrix of reaction stoichiometry, with `NSpec` rows and
#'   `NComp` columns.
#' @param NComp The number of components in the equilibrium system.
#' @param NSpec The number of species (reactions) in the equilibrium system.
#'
#' @returns A vector of `NSpec` species concentrations.
#'
#' @examples
#' # data(TestDataFreeConc, TestDataK, TestDataStoich)
#' # RCalcSpecConc(CConc = TestDataFreeConc[1:2], K = TestDataK[3:4], Stoich = TestDataStoich[3:4,])
#'
#' @noRd
RCalcSpecConc = function(CConc, K, Stoich, NComp = length(CConc),
                         NSpec = length(K)){#method=2

  # if(method == 1){
  #   # Method 1
  #   SConc = rep(1, NSpec)
  #   for (iSpec in 1:NSpec){
  #     Tmp = 1
  #     for (iComp in 1:NComp){
  #       Tmp = Tmp * CConc[iComp] ^ Stoich[iSpec, iComp]
  #     }
  #     SConc[iSpec] = Tmp * K[iSpec]
  #   }
  # } else if (method == 2){
  #   # Method 2
  SConc = rep(1, NSpec)
  for (iSpec in 1:NSpec){
    X = CConc ^ Stoich[iSpec,]
    SConc[iSpec] = prod(X) * K[iSpec]
  }
  # } else if (method == 3){
  #   # Method 3
  #   SConc = apply(CConc^t(Stoich), MAR = 2, FUN = prod) * K
  # }

  return(SConc)
}

