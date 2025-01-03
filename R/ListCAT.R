# Copyright 2024 Windward Environmental LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
#' mypfile = system.file(file.path("extdata", "ParameterFiles",
#'                                 "Cu_full_organic.dat4"),
#'                       package = "BLMEngineInR",
#'                       mustWork = TRUE)
#' ListCAT(ParamFile = mypfile)
ListCAT = function(ParamFile) {

  # error catching
  stopifnot(file.exists(ParamFile))

  # read the dimensions of the various elements of the reaction list
  SkipRows = 2
  Tmp = read.csv(file = ParamFile, header = FALSE, skip = SkipRows,
                 nrows = 9, strip.white = TRUE)
  NMass = as.integer(Tmp[1, 1])
  NInLab = as.integer(Tmp[2, 1])
  NInVar = as.integer(Tmp[3, 1])
  NInComp = as.integer(Tmp[4, 1])
  NDefComp = as.integer(Tmp[5, 1])
  NSpec = as.integer(Tmp[6, 1])
  NPhase = as.integer(Tmp[7, 1])
  NSpecialDef = as.integer(Tmp[8, 1])
  NCAT = as.integer(Tmp[9, 1])

  if (NCAT > 0) {
    # skip past other information
    SkipRows = SkipRows + 9 + 2
    SkipRows = SkipRows + NMass + 3
    SkipRows = SkipRows + NInLab + 2
    SkipRows = SkipRows + NInVar + 3
    SkipRows = SkipRows + NInComp + 3
    SkipRows = SkipRows + NDefComp + 4
    SkipRows = SkipRows + NSpec + 3
    SkipRows = SkipRows + NPhase + 2
    SkipRows = SkipRows + NSpecialDef + 3

    # Read in CAT table from parameter file
    CATab = read.csv(file = ParamFile, header = TRUE, skip = SkipRows,
                     nrows = NCAT, strip.white = TRUE)
    colnames(CATab) = c("Num", "CA (nmol/gw)", "Species", "Test Type",
                        "Duration", "Lifestage", "Endpoint", "Quantifier",
                        "References", "Miscellaneous")
  } else {
    CATab = BlankProblem()$CATab
  }

  # Return that table for the user
  return(CATab)

}
