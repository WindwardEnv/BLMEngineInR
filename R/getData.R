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
getData = function(inputFile, NInLab, NInVar, InVarName, InVarType,
                   NInComp, InCompName){

  # **for now, using the test data**

  # read in input file
  # -get number of observations
  NObs = as.integer(1)
  out = list(NObs = NObs)

  # -get labels
  out$Labels = array("", dim = c(NObs, NInLab))
  # to do: read in labels

  # -get temperatures
  out$SysTemp = array(25 + 273.15, dim = NObs)
  # to do: read in temperature

  # -get pH
  if (any(InVarType == "pH")){
    out$pH = array(7.0, dim = NObs)
    # to do: read in pH, convert to H
  }

  # -get organic matter and parse into components

  if (any(InVarType == "PercHA")){
    out$HA = array(0.1, dim = NObs)
    # to do: read in HA...might have more than one
  }

  # -get active fulvic acid fraction
  if (InVarFlag[InVarName == "AFA"]){
    out$AFA = array(1.0, dim = NObs)
  }

  # -get table of input concentrations
  if (inputFile == "Test") {
    data("TestDataTotalConc", envir = environment())
    TotConcObs = matrix(BLMEngineInR::TestDataTotalConc,
                        nrow = NObs, ncol = NInComp,
                        dimnames = list(Obs = 1, Comps = InCompName))
  } else if (inputFile == "Full_Inorg"){
    data("Full_InorgDataTotalConc", envir = environment())
    TotConcObs = matrix(BLMEngineInR::Full_InorgDataTotalConc,
                        nrow = NObs, ncol = NInComp,
                        dimnames = list(Obs = 1, Comps = InCompName))
  }
  out$TotConcObs = TotConcObs

  return(out)
}
