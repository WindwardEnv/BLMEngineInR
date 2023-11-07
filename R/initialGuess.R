#' Set the initial guess for component concentrations
#'
#' @param NComp integer(1); the number of components
#' @param CompName character vector (`NComp`); the names of the components
#' @param TotConcObs numeric vector (`NComp`); the total concentrations of components
#'
#' @return The initial guesses for the free ion concentrations of components
#'
#' @keywords internal
#'
#' @noRd
initialGuess = function(NComp, CompName, TotConcObs){

  CConc = TotConcObs
  names(CConc) = CompName

  out = CConc

  return(out)
}
