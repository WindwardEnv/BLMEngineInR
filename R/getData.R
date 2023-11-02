#' Get data from the input file
#'
#' `getData` reads in the input file and prepares it for input to CHESS.
#'
#' @param inputFile character(1); the path and file name to a BLM input file
#' @param NComp integer(1); the number of components
#' @param CompNames character vector (`NComp`); the names of the components
#' @param NInVar integer(1); the number of non-component input variables
#'   included in the input file
#' @param InVarName character vector (`NInVar`); the names of the non-component
#'   input variables included in the input file
#'
#' @return Returns a `list` object with the following components:
#' \describe{
#'  \item{\code{NObs}}{integer; the number of species}
#'  \item{\code{SiteLabels}}{character vector (\code{NObs}); the site label information for each observation}
#'  \item{\code{SampleLabels}}{character vector (\code{NObs}); the sample label information for each observation}
#'  \item{\code{Temp}}{(optional) numeric vector (\code{NObs}); input temperatures, in Kelvin(?)}
#'  \item{\code{HA}}{(optional) numeric vector (\code{NObs}); input humic acid fractions}
#'  \item{\code{AFA}}{(optional) numeric vector (\code{NObs}); input active fulvic acid fractions}
#'  \item{\cdoe{TotConcObs}}{numeric matrix (\code{NObs} x \code{NComp}); the total concentrations of each component}
#' }
#'
#' @keywords internal
#'
#' @noRd
getData = function(inputFile, NComp, CompNames, NInVar, InVarName){

  # **for now, using the test data**

  # read in input file
  # -get number of observations
  NObs = as.integer(1)
  out = list(NObs = NObs)

  # -get labels
  if (InVarFlag[InVarName == "Site"]){
    out$SiteLabels = array("A", dim = NObs)
  }
  if (InVarFlag[InVarName == "Sample"]){
    out$SampleLabels = array("a", dim = NObs)
  }

  # -get temperatures
  if (InVarFlag[InVarName == "Temp"]){
    out$Temp = array(25 + 273.15, dim = NObs)
  } else {
    # use standard temperature
    out$Temp = array(25 + 273.15, dim = NObs)
  }

  # -get humic acid fraction
  if (InVarFlag[InVarName == "HA"]){
    out$HA = array(0.1, dim = NObs)
  } else {
    # use default HA of 10%
    out$HA = array(0.1, dim = NObs)
  }

  # -get active fulvic acid fraction
  if (InVarFlag[InVarName == "AFA"]){
    out$AFA = array(1.0, dim = NObs)
  } else {
    out$AFA = array(1.0, dim = NObs)
  }

  # -get table of input concentrations
  if (inputFile == "Test") {
    data("TestDataTotalConc", envir = environment())
    TotConcObs = matrix(BLMEngineInR::TestDataTotalConc,
                        nrow = NObs, ncol = NComp,
                        dimnames = list(Obs = 1, Comps = CompNames))
  } else if (inputFile == "Full_Inorg"){
    data("Full_InorgDataTotalConc", envir = environment())
    TotConcObs = matrix(BLMEngineInR::Full_InorgDataTotalConc,
                        nrow = NObs, ncol = NComp,
                        dimnames = list(Obs = 1, Comps = CompNames))
  }
  out$TotConcObs = TotConcObs

  return(out)
}
