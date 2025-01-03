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

#' @name InVars
#'
#' @title Add or remove a input variables in a problem
#'
#' @param ThisProblem A list object with a structure like that returned by
#'   `BlankProblem()`.
#' @param InVarName A character vector with the name(s) of the input input
#'   variable(s).
#' @param InVarMCName A character vector with the name(s) of the mass
#'   compartments the new input variables are associated with. Does not need to
#'   be specified if `InVarMCR` is specified instead.
#' @param InVarType A character vector with the types of the new input
#'   variables. Must be one of "Temp", "pH", "WHAM-HA", "WHAM-FA", "WHAM-HAFA",
#'   "PercHA", "PercAFA", and "Misc".
#' @param InVarMCR (optional) A character vector with the indices of the mass
#'   compartments the new input variables are associated with. Only needs to be
#'   specified if `InVarMCName` is not specified.
#' @param InVarToRemove A character vector with names or indices of the input
#'   variable(s) to remove from `ThisProblem`.
#' @param DoCheck A logical value indicating whether checks should be performed
#'   on the incoming and outgoing problem objects. Defaults to `TRUE`, as you
#'   usually want to make sure something isn't awry, but the value is often set
#'   to `FALSE` when used internally (like in DefineProblem) so the problem is
#'   only checked once at the end.
#'
#' @return `ThisProblem`, with the changed input variables. If the input
#'   variable being added is pH, "H" and "OH" components will also be added as
#'   fixed activity components.
#'
#' @family problem manipulation functions
#'
#' @examples
#' print(carbonate_system_problem$InVar)
#' my_new_problem = carbonate_system_problem
#' my_new_problem = AddInVars(ThisProblem = my_new_problem,
#'                            InVarName = c("Humics", "Fulvics"),
#'                            InVarMCName = "Water",
#'                            InVarType = c("WHAM-HA","WHAM-FA"))
#' my_new_problem = RemoveInVars(ThisProblem = my_new_problem,
#'                               InVarToRemove = "Humics")
#' print(my_new_problem$InVar)
#'
NULL

#' @rdname InVars
#' @export
AddInVars = function(
  ThisProblem, InVarName, InVarMCName = NULL,
  InVarType = c("Temperature", "pH", "WHAM-FA", "WHAM-HA",
                "WHAM-HAFA", "PercHA", "PercAFA", "Misc"),
  InVarMCR = match(InVarMCName, ThisProblem$Mass$Name, nomatch = -1L),
  DoCheck = TRUE
) {

  if (DoCheck) {
    CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  }
  NewProblem = ThisProblem

  if ((NewProblem$ParamFile != "") &&
        !grepl("[(]modified[)]$", NewProblem$ParamFile)) {
    NewProblem$ParamFile = paste0(NewProblem$ParamFile, " (modified)")
  }

  # error checking
  if (any(is.na(c(InVarName, InVarType, InVarMCName, InVarMCR)))) {
    stop("NA arguments not allowed.")
  }
  if (!all(InVarMCName %in% ThisProblem$Mass$Name) ||
        !all(InVarMCR <= ThisProblem$N["Mass"])) {
    stop(
      "Mass compartment(s) specified in InVarMCName or InVarMCR does not exist."
    )
  }
  if (any((InVarName %in% ThisProblem$InVar$Name))) {
    stop(
      "Input Variable(s) (",
      paste(
        paste0("\"", InVarName[InVarName %in% ThisProblem$InVar$Name], "\""),
        collapse = ", "
      ),
      ") already exist."
    )
  }
  if (is.null(InVarMCName)) {
    InVarMCName = ThisProblem$Mass$Name[InVarMCR]
  }
  if (!all(InVarMCName %in% ThisProblem$Mass$Name[InVarMCR])) {
    stop("InVarMCName does not match InVarMCR.")
  }
  if (!all(InVarType %in%
             c("Temperature", "pH", "WHAM-FA", "WHAM-HA",
               "WHAM-HAFA", "PercHA", "PercAFA", "Misc"))) {
    stop(
      "Invalid InVarType specified (",
      InVarType[InVarType %in% c("Temperature", "pH", "WHAM-FA",
                                 "WHAM-HA", "WHAM-HAFA", "PercHA",
                                 "PercAFA", "Misc") == FALSE],
      ")"
    )
  }
  NInVarAdd = length(InVarName)
  if (length(InVarMCName) == 1) { InVarMCName = rep(InVarMCName, NInVarAdd)}
  if (length(InVarMCR) == 1) { InVarMCR = rep(InVarMCR, NInVarAdd)}
  if (length(InVarType) == 1) { InVarType = rep(InVarType, NInVarAdd)}

  # Add input variable
  NewProblem$InVar = rbind(
    ThisProblem$InVar,
    data.frame(Name = trimws(as.character(InVarName)),
               MCName = trimws(as.character(InVarMCName)),
               MCR = as.integer(InVarMCR),
               Type = trimws(as.character(InVarType)))
  )
  NewProblem$N["InVar"] = ThisProblem$N["InVar"] + NInVarAdd

  for (i in 1:NInVarAdd) {
    if (InVarType[i] == "pH") {
      NewProblem = AddDefComps(ThisProblem = NewProblem,
                               DefCompName = c("H", "OH"),
                               DefCompFromVar = c(InVarName[i], "KW/H"),
                               DefCompFromNum = NA,
                               DefCompCharge = c(1L, -1L),
                               DefCompMCName = InVarMCName[i],
                               DefCompType = "FixedAct",
                               DefCompActCorr = "Debye",
                               DefCompSiteDens = 1.0,
                               DefCompMCR = InVarMCR[i],
                               InDefComp = TRUE,
                               DoCheck = DoCheck)
      # warning(paste("Defined component added for pH as a fixed activity",
      #               "component. Change manually if you wish to represent",
      #               "pH as a fixed concentration."))
    }
    if (grepl("WHAM", InVarType[i]) &&
          (!is.na(NewProblem$WHAM$Ver) || !is.na(NewProblem$WHAM$File))) {
      NewProblem = ExpandWHAM(ThisProblem = NewProblem)
    }
  }

  if (DoCheck) {
    CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  }
  return(NewProblem)

}


#' @rdname InVars
#' @export
RemoveInVars = function(ThisProblem, InVarToRemove, DoCheck = TRUE) {

  if (DoCheck) {
    CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  }
  NewProblem = ThisProblem

  if ((NewProblem$ParamFile != "") &&
        !grepl("[(]modified[)]$", NewProblem$ParamFile)) {
    NewProblem$ParamFile = paste0(NewProblem$ParamFile, " (modified)")
  }

  InVarToRemoveOrig = InVarToRemove
  if (is.character(InVarToRemove)) {
    InVarToRemove = match(InVarToRemove, ThisProblem$InVar$Name)
  }
  if (any(is.na(InVarToRemove))) {
    stop("Input Variable(s) (",
         paste(InVarToRemoveOrig[is.na(InVarToRemove)], collapse = ", "),
         ") do not exist.")
  }
  if (any(InVarToRemove <= 0L)) {
    stop("Invalid index in InVarToRemove (",
         InVarToRemove[InVarToRemove <= 0L], ").")
  }
  if (any(InVarToRemove > ThisProblem$N["InVar"])) {
    stop("There are ", ThisProblem$N["InVar"], " input variables, ",
         "trying to remove the #(",
         InVarToRemove[InVarToRemove > ThisProblem$N["InVar"]],
         ") element(s).")
  }

  # Remove DefComps that depend on input variables
  DefCompToRemove = which(NewProblem$DefComp$FromVar %in%
                            NewProblem$InVar$Name[InVarToRemove])
  if (length(DefCompToRemove) >= 1) {
    NewProblem = RemoveDefComps(NewProblem, DefCompToRemove, DoCheck = DoCheck)
  }

  # Remove input variables
  NewProblem$InVar = ThisProblem$InVar[-InVarToRemove, , drop = FALSE]
  rownames(NewProblem$InVar) = NULL
  NewProblem$N["InVar"] = ThisProblem$N["InVar"] - length(InVarToRemove)

  if (DoCheck) {
    CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  }
  return(NewProblem)

}
