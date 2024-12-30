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

#' @name InLabs
#'
#' @title Add or remove input labels in a problem
#'
#' @param ThisProblem A list object with a structure like that returned by
#'   `BlankProblem()`.
#' @param InLabName A character vector with the name(s) of the new input
#'   label(s).
#' @param InLabToRemove A character vector with names or indices of the input
#'   label(s) to remove from `ThisProblem`.
#' @param DoCheck A logical value indicating whether checks should be performed
#'   on the incoming and outgoing problem objects. Defaults to `TRUE`, as you
#'   usually want to make sure something isn't awry, but the value is often set
#'   to `FALSE` when used internally (like in DefineProblem) so the problem is
#'   only checked once at the end.
#'
#' @return `ThisProblem`, with the edited input labels.
#'
#' @family problem manipulation functions
#'
#' @examples
#' my_new_problem = carbonate_system_problem
#' print(carbonate_system_problem$InLabName) # ID only
#'
#' my_new_problem = AddInLabs(ThisProblem = my_new_problem, InLabName = "ID2")
#' my_new_problem = RemoveInLabs(ThisProblem = my_new_problem, InLabToRemove = "ID")
#'
#' print(my_new_problem$InLabName) # ID2 only
#'
NULL

#' @rdname InLabs
#' @export
AddInLabs = function(ThisProblem, InLabName, DoCheck = TRUE) {

  if (DoCheck) {
    CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  }
  NewProblem = ThisProblem

  if ((NewProblem$ParamFile != "") &&
      !grepl("[(]modified[)]$", NewProblem$ParamFile)) {
    NewProblem$ParamFile = paste0(NewProblem$ParamFile, " (modified)")
  }

  # error checking
  if (any((InLabName %in% ThisProblem$InLabName))) {
    stop(paste0(
      "Input Label(s) (",
      paste(paste0("\"", InLabName[InLabName %in% ThisProblem$InLabName], "\""), collapse = ", "),
      ") already exist."
    ))
  }
  if (any(is.na(InLabName))) {
    stop("NA arguments not allowed.")
  }

  # add input labels
  NewProblem$InLabName = c(ThisProblem$InLabName,
                           trimws(as.character(InLabName)))
  NewProblem$N["InLab"] = ThisProblem$N["InLab"] + length(InLabName)

  if (DoCheck) {
    CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  }
  return(NewProblem)

}

#' @rdname InLabs
#' @export
RemoveInLabs = function(ThisProblem, InLabToRemove, DoCheck = TRUE) {

  if (DoCheck) {
    CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  }
  NewProblem = ThisProblem

  if ((NewProblem$ParamFile != "") &&
      !grepl("[(]modified[)]$", NewProblem$ParamFile)) {
    NewProblem$ParamFile = paste0(NewProblem$ParamFile, " (modified)")
  }

  InLabToRemoveOrig = InLabToRemove
  if (is.character(InLabToRemove)) {
    InLabToRemove = which(ThisProblem$InLabName %in% InLabToRemove)
  }
  InLabToRemove = as.integer(InLabToRemove)
  if (length(InLabToRemove) < 1) {
    stop(paste0("Input Label(s) \"", InLabToRemoveOrig, "\" does not exist."))
  }
  if (any(InLabToRemove <= 0L)) {
    stop("Invalid index in InLabToRemove (",
         InLabToRemove[InLabToRemove <= 0L],").")
  }
  if (any(InLabToRemove > ThisProblem$N["InLab"])) {
    stop("There are ", ThisProblem$N["InLab"], " input labels, ",
         "trying to remove the #(",
         InLabToRemove[InLabToRemove > ThisProblem$N["InLab"]],
         ") element(s).")
  }

  # Remove input labels
  NewProblem$InLabName = NewProblem$InLabName[-InLabToRemove]
  NewProblem$N["InLab"] = ThisProblem$N["InLab"] - length(InLabToRemove)

  if (DoCheck) {
    CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  }
  return(NewProblem)

}

