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

#' @name CriticalValues
#'
#' @title Edit Critical Values Table
#'
#' @description `AddCriticalValues` will add one or more rows to the critical
#'   accumulation table (CAT or CATab), while `RemoveCriticalValues` will remove
#'   one or more rows,
#'
#' @param ThisProblem A list object with a structure like that returned by
#'   `BlankProblem()`.
#' @param CATab a data.frame object with, at a minimum, columns named
#'   `CA`/`CA (nmol/gw)`, `Species`, `Test.Type`/`Test Type`, `Duration`,
#'   `Lifestage`, `Endpoint`, `Quantifier`, `References`, `Miscellaneous`.
#'   See optional parameter descriptions for further descriptions of each of
#'   those columns.
#' @param CA (optional) a numeric vector of the critical accumulation value(s)
#'   in nmol/gw.
#' @param Species (optional) a character vector of the species names to include
#'   for the corresponding `CA` value.
#' @param Test.Type (optional) a character vector of the test type (e.g.,
#'   "Acute" or "Chronic") to include for the corresponding `CA` value.
#' @param Duration (optional) a character vector of the Duration to include for
#'   the corresponding `CA` value. Can also be `"DIV=#.##"` for FAV, FCV, WQS,
#'   and HC5 critical values.
#' @param Lifestage (optional) a character vector of the organism Lifestage to
#'   include for the corresponding `CA` value. Can also be `"ACR=#.##"` for FAV,
#'   WQS, and HC5 critical values.
#' @param Endpoint (optional) a character vector of the Endpoint to include for
#'   the corresponding `CA` value. This can also be either `"FAV"`, `"FCV"`,
#'   `"HC5"`, `"WQS"`, `"CMC"`, or `"CCC"` to indicate this critical value
#'   calculates one of those water quality standards.
#' @param Quantifier (optional) a character vector of the Quantifier (e.g.,
#'   EC50, NOEC, ...) to include for the corresponding `CA` value. May also be
#'   `NA` if this is a WQS value.
#' @param References (optional) a character vector of the list of References to
#'   include for the corresponding `CA` value. Each `CA` value requires a single
#'   character string with no line breaks.
#' @param Miscellaneous (optional) a character vector of the miscellaneous
#'   information (e.g., how the value was calculated, test conditions not
#'   covered by other columns, etc.) to include for the corresponding `CA`
#'   value.
#' @param CAToRemove an integer vector - the indices/row numbers of the critical
#'   values to remove from the table.
#' @param DoCheck A logical value indicating whether checks should be performed
#'   on the incoming and outgoing problem objects. Defaults to `TRUE`, as you
#'   usually want to make sure something isn't awry, but the value is often set
#'   to `FALSE` when used internally (like in DefineProblem) so the problem is
#'   only checked once at the end.
#'
#'
#' @return The edited version of `ThisProblem`.
#'
#' @family problem manipulation functions
#'
#' @examples
#' my_new_problem = carbonate_system_problem
#'
#' my_new_problem = AddCriticalValues(
#'   ThisProblem = my_new_problem,
#'   CA = 12345,
#'   Species = "A. species",
#'   Test.Type = "Acute",
#'   Duration = "24h",
#'   Lifestage = "adult",
#'   Endpoint = "survival",
#'   Quantifier = "LC50",
#'   References = "thin air",
#'   Miscellaneous = "individual data point"
#' )
#'
#' lots_of_data = data.frame(CA = runif(26),
#'                           Species = paste0(LETTERS,". species"),
#'                           Test.Type = "Acute",
#'                           Duration = "24h",
#'                           Lifestage = "adult",
#'                           Endpoint = "survival",
#'                           Quantifier = "LC50",
#'                           References = "thin air")
#' my_new_problem = AddCriticalValues(
#'   ThisProblem = my_new_problem,
#'   CATab = lots_of_data
#' )
#'
#' my_new_problem = RemoveCriticalValues(
#'   ThisProblem = my_new_problem,
#'   CAToRemove = which((my_new_problem$CATab$Species == "A. species") &
#'                        is.na(my_new_problem$CATab$Miscellaneous))
#' )
#'
#' print(my_new_problem$CATab)
#'
NULL


#' @rdname CriticalValues
#' @export
AddCriticalValues = function(ThisProblem, CATab = data.frame(),
                             CA = CATab[, which(colnames(CATab) %in% c("CA","CA (nmol/gw)"))[1]],
                             Species = CATab$Species,
                             Test.Type = CATab[, which(colnames(CATab) %in% c("Test.Type","Test Type"))[1]],
                             Duration = CATab$Duration,
                             Lifestage = CATab$Lifestage,
                             Endpoint = CATab$Endpoint,
                             Quantifier = CATab$Quantifier,
                             References = CATab$References,
                             Miscellaneous = CATab$Miscellaneous,
                             DoCheck = TRUE) {

  if (DoCheck) {
    CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  }
  NewProblem = ThisProblem

  if ((NewProblem$ParamFile != "") &&
      !grepl("[(]modified[)]$", NewProblem$ParamFile)) {
    NewProblem$ParamFile = paste0(NewProblem$ParamFile, " (modified)")
  }

  NumNewCA = length(CA)
  if (is.null(Species)) { Species = rep(NA_character_, NumNewCA) }
  if (is.null(Test.Type)) { Test.Type = rep(NA_character_, NumNewCA) }
  if (is.null(Duration)) { Duration = rep(NA_character_, NumNewCA) }
  if (is.null(Lifestage)) { Lifestage = rep(NA_character_, NumNewCA) }
  if (is.null(Endpoint)) { Endpoint = rep(NA_character_, NumNewCA) }
  if (is.null(Quantifier)) { Quantifier = rep(NA_character_, NumNewCA) }
  if (is.null(References)) { References = rep(NA_character_, NumNewCA) }
  if (is.null(Miscellaneous)) { Miscellaneous = rep(NA_character_, NumNewCA) }

  NoIDInfo = is.na(Species) & is.na(Test.Type) & is.na(Duration) &
    is.na(Lifestage) & is.na(Endpoint) & is.na(Quantifier) &
    is.na(References) & is.na(Miscellaneous)
  if (any(NoIDInfo)) {
    warning("Some of the new critical values do not have any identifying information.")
  }

  CA = as.numeric(CA)
  Species = as.character(Species)
  Test.Type = as.character(Test.Type)
  Duration = as.character(Duration)
  Lifestage = as.character(Lifestage)
  Endpoint = as.character(Endpoint)
  Quantifier = as.character(Quantifier)
  References = as.character(References)
  Miscellaneous = as.character(Miscellaneous)

  NewCATab = data.frame(
    Num = NA_integer_,
    CA = CA,
    Species = Species,
    Test.Type = Test.Type,
    Duration = Duration,
    Lifestage = Lifestage,
    Endpoint = Endpoint,
    Quantifier = Quantifier,
    References = References,
    Miscellaneous = Miscellaneous
  )
  NewCATab$Num = ThisProblem$N["CAT"] + seq(1L, NumNewCA)
  NewProblem$N["CAT"] = NewProblem$N["CAT"] + NumNewCA

  NewProblem$CATab = rbind(
    ThisProblem$CATab,
    NewCATab
  )

  if (DoCheck) {
    CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  }

  return(NewProblem)

}


#' @rdname CriticalValues
#' @export
RemoveCriticalValues = function(ThisProblem, CAToRemove, DoCheck = TRUE) {

  if (DoCheck) {
    CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  }
  NewProblem = ThisProblem

  if ((NewProblem$ParamFile != "") &&
      !grepl("[(]modified[)]$", NewProblem$ParamFile)) {
    NewProblem$ParamFile = paste0(NewProblem$ParamFile, " (modified)")
  }

  CAToRemove = unique(CAToRemove)
  if (any(CAToRemove <= 0L)) {
    stop("Invalid index in CAToRemove (",
         CAToRemove[CAToRemove <= 0L],").")
  }
  if (any(CAToRemove > ThisProblem$N["CAT"])) {
    stop(paste0("There are ", ThisProblem$N["CAT"], " CAs, ",
                "trying to remove the #(",
                paste(CAToRemove[CAToRemove > ThisProblem$N["CAT"]],
                      collapse = ", "),
                ") element(s)."))
  }

  NewProblem$CATab = ThisProblem$CATab[-CAToRemove, ]
  NewProblem$N["CAT"] = ThisProblem$N["CAT"] - length(CAToRemove)
  if (NewProblem$N["CAT"] > 0) {
    NewProblem$CATab$Num = seq(1, NewProblem$N["CAT"])
  }

  if (DoCheck) {
    CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  }

  return(NewProblem)

}
