#' Get data from the input file
#'
#' `getData` reads in the input file and prepares it for input to CHESS.
#'
#' @param inputFile character(1); the path and file name to a BLM input file
#' @param NComp integer(1); the number of components
#'
#' @return Returns a `list` object with the following components:
#' \describe{
#'  \item{`obsLabels`}{character matrix (`NObs` x 2); the site and sample label information for each observation}
#'  \item{`totConcObs`}{numeric matrix (`NObs` x `NComp`); the total concentrations of each component}
#'  \item{`NObs`}{integer; the number of species}
#' }
#'
#' @noRd
getData = function(inputFile, NComp){

  # for now, using the test data
  data("TestDataTotalConc")
  NObs = 1
  obsLabels = matrix(
    c("A", "a"),
    nrow = NObs,
    ncol = 2,
    dimnames = list(Obs = 1:NObs, c("Site Label", "Sample Label"))
  )
  totConcObs = matrix(TestDataTotalConc,
                      nrow = NObs,
                      dimnames = list(Obs = 1, Comps = NULL))

  # read in input file
  # -get number of observations
  # -get labels
  # -get temperatures
  # -get table of input concentrations

  out = list(
    obsLabels = obsLabels,
    totConcObs = totConcObs,
    NObs = NObs
  )

  return(out)
}
