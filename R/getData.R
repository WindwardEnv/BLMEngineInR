#' Get data from the input file
#'
#' `getData` reads in the input file and prepares it for input to CHESS.
#'
#' @param inputFile character(1); the path and file name to a BLM input file
#' @param NComp integer(1); the number of components
#' @param CompNames character vector (`NComp`); the names of the components
#'
#' @return Returns a `list` object with the following components:
#' \describe{
#'  \item{\code{obsLabels}}{character matrix (`NObs` x 2); the site and sample label information for each observation}
#'  \item{\cdoe{totConcObs}}{numeric matrix (`NObs` x `NComp`); the total concentrations of each component}
#'  \item{\code{NObs}}{integer; the number of species}
#' }
#'
#' @keywords internal
#'
#' @noRd
getData = function(inputFile, NComp, CompNames){

  # for now, using the test data
  NObs = as.integer(1)
  obsLabels = matrix(
    c("A", "a"),
    nrow = NObs,
    ncol = 2,
    dimnames = list(Obs = 1:NObs, c("Site Label", "Sample Label"))
  )
  if (inputFile == "Test") {
    data("TestDataTotalConc", envir = environment())
    totConcObs = matrix(BLMEngineInR::TestDataTotalConc,
                        nrow = NObs, ncol = NComp,
                        dimnames = list(Obs = 1, Comps = CompNames))
  } else if (inputFile == "Full_Inorg"){
    data("Full_InorgDataTotalConc", envir = environment())
    totConcObs = matrix(BLMEngineInR::Full_InorgDataTotalConc,
                        nrow = NObs, ncol = NComp,
                        dimnames = list(Obs = 1, Comps = CompNames))
  }


  # read in input file
  # -get number of observations
  # -get labels
  # -get temperatures
  # -get table of input concentrations

  out = list(
    NObs = NObs,
    obsLabels = obsLabels,
    totConcObs = totConcObs
  )

  return(out)
}
