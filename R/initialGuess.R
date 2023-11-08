#' Set the initial guess for component concentrations
#'
#' @param NComp integer(1); the number of components
#' @param CompName character vector (`NComp`); the names of the components
#' @param TotConc numeric vector (`NComp`); the total concentrations of components
#'
#' @return The initial guesses for the free ion concentrations of components
#'
#' @keywords internal
#'
#' @noRd
initialGuess = function(NComp, CompName, TotConc){

  CompConc = TotConc
  names(CompConc) = CompName

  out = CompConc

  return(out)
}
