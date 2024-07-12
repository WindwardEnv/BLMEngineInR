#' List Critical Accumulation Table
#'
#' List out the critical accumulation table for the user to allow them to pick
#' which CAT number they should specify for a toxicity run where the critical
#' value is coming from the table in the parameter file.
#'
#' @param ParamFile character string; the file name and path of the parameter
#'   file.
#'
#' @return A \code{data.frame} object with the CAT table in the given parameter
#'   file. Columns include:
#' \describe{
#'    \item{\code{Num}}{the number or index in the table}
#'    \item{\code{`CA (nmol/gw)`}}{the critical accumulation in units of
#'      nmol/gw}
#'    \item{\code{Species}}{species name or CA significance, such as HC5 or
#'      FAV}
#'    \item{\code{`Test Type`}}{acute or chronic}
#'    \item{\code{Duration}}{test duration (e.g., 48 h)}
#'    \item{\code{Lifestage}}{age or size of the organisms}
#'    \item{\code{Endpoint}}{toxicity endpoint (e.g., mortality, reproduction)}
#'    \item{\code{Quantifier}}{endpoint quantifier or effect level (e.g., LC50,
#'      EC10, NOEC)}
#'    \item{\code{References}}{citations of sources with the toxicity data that
#'      went into calculating the CA, or the citation of the HC5 or FAV}
#'    \item{\code{Miscellanous}}{other notes or comments (e.g., number of data
#'      points or methods of calculating)}
#' }
#'
#' @export
#'
#' @examples
#' ## Not Run
#' # ListCAT("my_parameter_file.dat")
#' ## End Not Run
ListCAT = function(ParamFile){

  # error catching
  stopifnot(file.exists(ParamFile))

  # read the dimensions of the various elements of the reaction list
  skipRows = 2
  tmp = read.csv(file = ParamFile, header = FALSE, skip = skipRows,
                 nrows = 9, strip.white = T)
  NMass = tmp[1, 1]
  NInLab = tmp[2, 1]
  NInVar = tmp[3, 1]
  NInComp = tmp[4, 1]
  NDefComp = tmp[5, 1]
  NSpec = tmp[6, 1]
  NPhase = tmp[7, 1]
  NSpecialDef = tmp[8, 1]
  NCAT = tmp[9, 1]

  # skip past other information
  skipRows = skipRows + 9 + 2
  skipRows = skipRows + NMass + 3
  skipRows = skipRows + NInLab + 2
  skipRows = skipRows + NInVar + 3
  skipRows = skipRows + NInComp + 3
  skipRows = skipRows + NDefComp + 4
  skipRows = skipRows + NSpec + 3
  skipRows = skipRows + NPhase + 2
  skipRows = skipRows + NSpecialDef + 3

  # Read in CAT table from parameter file
  CATab = read.csv(file = ParamFile, header = TRUE, skip = skipRows,
                   nrows = NCAT, strip.white = T)
  colnames(CATab) = c("Num","CA (nmol/gw)","Species","Test Type","Duration",
                      "Lifestage","Endpoint","Quantifier","References",
                      "Miscellaneous")

  # Return that table for the user
  return(CATab)

}
