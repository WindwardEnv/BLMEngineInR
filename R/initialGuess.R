#' Set the initial guess for component concentrations
#'
#' @param NComp integer(1); the number of components
#' @param CompNames character vector (`NComp`); the names of the components
#' @param totConcObs numeric vector (`NComp`); the total concentrations of components
#'
#' @return The initial guesses for the free ion concentrations of components
#'
#' @keywords internal
#'
#' @noRd
initialGuess = function(NComp, CompNames, totConcObs){

  # # For now, we're going to use test data, setting the initial "guess" to the
  # # actual component free ion concentrations
  # if (inputFile == "Test") {
  #   data("TestDataFreeConc")
  #   CConc = TestDataFreeConc[1:NComp]
  # } else if (inputFile == "Full_Inorg"){
  #   data("Full_InorgDataFreeConc")
  #   CConc = Full_InorgDataFreeConc[1:NComp]
  # }
  CConc = totConcObs
  names(CConc) = CompNames

  out = CConc

  return(out)
}
