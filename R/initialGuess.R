#' Set the initial guess for component concentrations
#'
#' @param NComp integer(1); the number of components
#' @param CompName character vector (`NComp`); the names of the components
#' @param TotConc numeric vector (`NComp`); the total concentrations of components
#' @param CompType character vector (`NComp`); the component type
#'
#' @return The initial guesses for the free ion concentrations of components
#'
#' @keywords internal
#'
#' @noRd
initialGuess = function(NComp, CompName, TotConc, CompType){


  CompConc = TotConc / 100
  names(CompConc) = CompName

  CompConc[CompType == "FixedAct"] = TotConc[CompType == "FixedAct"]

  out = CompConc

  return(out)
}
