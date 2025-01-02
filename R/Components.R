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

#' @name Components
#'
#' @title Add or remove components in the problem
#'
#' @description A component should be added either as an input component (with
#'   `AddInComps`) or a defined component (with `AddDefComps`). Both of those
#'   functions will call the `AddComponents` function, but using either
#'   `AddInComps` and `AddDefComps` ensures that it's very clear where the
#'   inputs come from.
#'
#' @param ThisProblem A list object with a structure like that returned by
#'   `BlankProblem()`.
#' @param CompName,InCompName,DefCompName A character vector with the name(s) of
#'   the components to be added.
#' @param DefCompFromNum A numeric vector with the numeric values used to derive
#'   the component. Specify `NA` if the defined component uses a variable to
#'   define it.
#' @param DefCompFromVar A character vector with the variable names used to
#'   derive the component. Specify `NA` if the defined component uses a number
#'   to define it.
#' @param CompCharge,InCompCharge,DefCompCharge An integer vector with the
#'   charge(s) of the components to be added.
#' @param CompMCName,InCompMCName,DefCompMCName A character vector with the
#'   name(s) of the mass compartments the new components are associated with.
#'   Does not need to be specified if `CompMCR`/`InCompMCR`/`DefCompMCR` is
#'   specified instead.
#' @param CompType,InCompType,DefCompType A character vector with the types of
#'   the new input variables. Must be one of "MassBal", "FixedAct", "FixedConc",
#'   "DonnanHA", or "DonnanFA".
#' @param CompActCorr,InCompActCorr,DefCompActCorr A character vector with the
#'   activity correction method(s) of the new components. Must be one of "None",
#'   "Debye", "Davies", "DonnanHA", "DonnanFA", "WHAMHA", or "WHAMFA".
#'   Generally, "DonnanHA", "DonnanFA", "WHAMHA", and "WHAMFA" will only be used
#'   internally.
#' @param CompSiteDens,DefCompSiteDens A numeric vector with the binding site
#'   densities of the new components. `AddInComps` assumes a site density of
#'   1.0.
#' @param CompMCR,InCompMCR,DefCompMCR (optional) A character vector with the
#'   indices of the mass compartments the new components are associated with.
#'   Only needs to be specified if `CompMCName`/`InCompMCName`/`DefCompMCName`
#'   is not specified.
#' @param InDefComp A logical value indicating if this is a defined component
#'   from the parameter file (`TRUE`) or was added from another process, such as
#'   `ExpandWHAM` (`FALSE`).
#' @param ComponentToRemove,InCompToRemove,DefCompToRemove A character vector
#'   with names or indices of the component(s) to remove from `ThisProblem`. It
#'   is safer to use a name, since the index of the component may be different
#'   within `ThisProblem$Comp$Name` versus `ThisProblem$InCompName` versus
#'   `ThisProblem$DefComp$Name`.
#' @param DoCheck A logical value indicating whether checks should be performed
#'   on the incoming and outgoing problem objects. Defaults to `TRUE`, as you
#'   usually want to make sure something isn't awry, but the value is often set
#'   to `FALSE` when used internally (like in DefineProblem) so the problem is
#'   only checked once at the end.
#'
#' @return `ThisProblem`, with the edits done to the component list, including
#'   "trickle-down" changes, such as removing formation reactions that used a
#'   now-removed component.
#'
#' @seealso The in-depth example in [BlankProblem()] will show all problem
#'   manipulation functions in use.
#'
#' @family problem manipulation functions
#'
#' @examples
#' print(carbonate_system_problem$Comp)
#' my_new_problem = carbonate_system_problem
#' my_new_problem = AddInComps(ThisProblem = my_new_problem,
#'                             InCompName = "Ca",
#'                             InCompCharge = 2,
#'                             InCompMCName = "Water",
#'                             InCompType = "MassBal",
#'                             InCompActCorr = "Debye")
#' print(my_new_problem$Comp)
NULL


#' @rdname Components
#' @export
AddComponents = function(ThisProblem, CompName,  CompCharge, CompMCName = NULL,
                         CompType, CompActCorr, CompSiteDens = 1.0,
                         CompMCR = match(CompMCName, ThisProblem$Mass$Name, nomatch = -1L),
                         DoCheck = TRUE) {

  if (DoCheck) {
    CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  }
  NewProblem = ThisProblem

  if ((NewProblem$ParamFile != "") &&
      !grepl("[(]modified[)]$", NewProblem$ParamFile)) {
    NewProblem$ParamFile = paste0(NewProblem$ParamFile, " (modified)")
  }

  # coerce to correct types
  CompName = trimws(as.character(CompName))
  CompCharge = as.integer(CompCharge)
  if(!is.null(CompMCName)) {CompMCName = trimws(as.character(CompMCName))}
  CompMCR = as.integer(CompMCR)
  CompType = trimws(as.character(CompType))
  CompActCorr = trimws(as.character(CompActCorr))
  CompSiteDens = as.numeric(CompSiteDens)

  # error checking
  if (any(CompName %in% ThisProblem$Spec$Name)) {
    stop(paste0(
      "Input Component(s) (",
      paste(paste0("\"", CompName[CompName %in% ThisProblem$Spec$Name], "\""),
            collapse = ", "),
      ") already exist as a component or species."
    ))
  }
  if (!all(CompMCName %in% ThisProblem$Mass$Name) ||
      !all(CompMCR <= ThisProblem$N["Mass"])) {
    stop("Mass compartment(s) specified in CompMCName or CompMCR does not exist.")
  }

  if (any(is.na(c(CompName, CompCharge, CompMCName,
                  CompType, CompActCorr, CompMCR, CompSiteDens)))) {
    stop("NA arguments not allowed.")
  }
  if (!all(CompType %in% c("FixedConc","FixedAct","MassBal","DonnanHA",
                           "DonnanFA","WHAMHA","WHAMFA"))) {
    stop("CompType values must be FixedConc, FixedAct, MassBal, DonnanHA, DonnanFA, WHAMHA, WHAMFA.")
  }
  if(is.null(CompMCName)) {
    CompMCName = ThisProblem$Mass$Name[CompMCR]
  }
  if (!all(CompMCName %in% ThisProblem$Mass$Name[CompMCR])) {
    stop("CompMCName does not match CompMCR.")
  }
  NCompAdd = length(CompName)


  # add Components
  if (length(CompType) == 1){
    CompType = rep(CompType, NCompAdd)
  }
  NewProblem$Comp = rbind(ThisProblem$Comp,
                          data.frame(
                            Name = trimws(as.character(CompName)),
                            Charge = as.integer(CompCharge),
                            MCName = trimws(as.character(CompMCName)),
                            MCR = as.integer(CompMCR),
                            Type = trimws(as.character(CompType)),
                            ActCorr = trimws(as.character(CompActCorr)),
                            SiteDens = as.numeric(CompSiteDens)
                          ))
  NewProblem$SpecStoich = cbind(
    NewProblem$SpecStoich,
    matrix(data = 0L,
           nrow = ThisProblem$N["Spec"],
           ncol = NCompAdd,
           dimnames = list(NewProblem$Spec$Name, CompName))
  )
  NewProblem$PhaseStoich = cbind(
    NewProblem$PhaseStoich,
    matrix(data = 0L,
           nrow = ThisProblem$N["Phase"],
           ncol = NCompAdd,
           dimnames = list(NewProblem$Phase$Name, CompName))
  )
  NewProblem$N["Comp"] = ThisProblem$N["Comp"] + NCompAdd

  # Add components to species list
  SpecType = ifelse(CompType %in% c("FixedAct","FixedConc","MassBal"),
                    "Normal", CompType)
  NewProblem = AddSpecies(ThisProblem = NewProblem,
                          SpecEquation = paste(
                            CompName,"=",
                            ifelse(SpecType %in% c("DonnanHA","DonnanFA"),
                                   "",
                                   paste0("1 * ", CompName))),
                          SpecMCName = CompMCName,
                          SpecType = SpecType,
                          SpecActCorr = CompActCorr,
                          SpecLogK = 0,
                          SpecDeltaH = 0,
                          SpecTempKelvin = 0,
                          SpecMCR = CompMCR,
                          InSpec = FALSE,
                          DoCheck = FALSE)

  # Update Special Definitions
  NewProblem$BLMetal$SpecsR = match(NewProblem$BLMetal$Name, NewProblem$Spec$Name)

  if (DoCheck) {
    CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  }

  return(NewProblem)

}


#' @rdname Components
#' @export
RemoveComponents = function(ThisProblem, ComponentToRemove, DoCheck = TRUE) {

  if (DoCheck) {
    CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  }
  NewProblem = ThisProblem

  if ((NewProblem$ParamFile != "") &&
      !grepl("[(]modified[)]$", NewProblem$ParamFile)) {
    NewProblem$ParamFile = paste0(NewProblem$ParamFile, " (modified)")
  }

  ComponentToRemoveOrig = ComponentToRemove
  if (is.character(ComponentToRemove)) {
    ComponentToRemove = match(ComponentToRemove, ThisProblem$Comp$Name)
  }
  if (any(is.na(ComponentToRemove))) {
    stop(paste0("Component(s) (",
                paste(ComponentToRemoveOrig[is.na(ComponentToRemove)],
                      collapse = ", "),
                ") does not exist."))
  }
  if (any(ComponentToRemove <= 0L)) {
    stop("Invalid index in ComponentToRemove (",
         ComponentToRemove[ComponentToRemove <= 0L],").")
  }
  if (any(ComponentToRemove > ThisProblem$N["Comp"])) {
    stop(paste0("There are ", ThisProblem$N["Comp"], " Components, ",
                "trying to remove the #(",
                paste(ComponentToRemove[ComponentToRemove > ThisProblem$N["Comp"]],
                      collapse = ", "),
                ") element(s)."))
  }
  # if (ThisProblem$Comp$Name[ComponentToRemove] %in% ThisProblem$DefComp$Name) {
  #   stop("Trying to remove a defined component with RemoveComponents(). ",
  #        "Please use RemoveDefComps() instead.")
  # }
  # if (ThisProblem$Comp$Name[ComponentToRemove] %in% ThisProblem$InCompName) {
  #   stop("Trying to remove an input component with RemoveComponents(). ",
  #        "Please use RemoveInComps() instead.")
  # }

  # Remove Species that depend on the component
  SpeciesToRemove = unique(c(
    ComponentToRemove,
    which(apply(NewProblem$SpecStoich[, ComponentToRemove, drop = FALSE],
                MARGIN = 1, FUN = function(X){any(X != 0L)}))
  ))
  NewProblem = RemoveSpecies(ThisProblem = NewProblem,
                             SpeciesToRemove = SpeciesToRemove,
                             DoCheck = FALSE)

  # Remove Phases that depend on the component
  if (ThisProblem$N["Phase"] > 0) {
    PhasesToRemove = which(apply(NewProblem$PhaseStoich[, ComponentToRemove, drop = FALSE],
                                 MARGIN = 1, FUN = function(X){any(X != 0L)}))
    if (length(PhasesToRemove) >= 1) {
      NewProblem = RemovePhases(NewProblem, PhasesToRemove, DoCheck = FALSE)
    }
  }

  # Remove DefComps that depend on the component
  DefCompToRemove = which(NewProblem$DefComp$FromVar %in%
                            ThisProblem$Comp$Name[ComponentToRemove])
  if (length(DefCompToRemove) >= 1) {
    NewProblem = RemoveDefComps(NewProblem, DefCompToRemove, DoCheck = FALSE)
  }

  # Remove InComps that are the component to remove
  InCompToRemove = which(NewProblem$InCompName %in%
                           ThisProblem$Comp$Name[ComponentToRemove])
  if (length(InCompToRemove) >= 1) {
    NewProblem$InCompName = NewProblem$InCompName[-InCompToRemove]
    NewProblem$N["InComp"] = ThisProblem$N["InComp"] - length(InCompToRemove)
  }

  # Remove Components
  NewProblem$Comp = ThisProblem$Comp[-ComponentToRemove, , drop = FALSE]
  rownames(NewProblem$Comp) = NULL
  NewProblem$N["Comp"] = ThisProblem$N["Comp"] - length(ComponentToRemove)
  NewProblem$SpecStoich = NewProblem$SpecStoich[, -ComponentToRemove, drop = FALSE]
  NewProblem$PhaseStoich = NewProblem$PhaseStoich[, -ComponentToRemove, drop = FALSE]

  # Update Special Definitions
  NewProblem$BL = ThisProblem$BL[ThisProblem$BL$Name %in% ThisProblem$Comp$Name[ComponentToRemove] == FALSE, , drop = FALSE]
  NewProblem$BL$CompR = match(NewProblem$BL$Name, NewProblem$Comp$Name)
  NewProblem$N["BL"] = length(NewProblem$BL$Name)

  NewProblem$Metal = ThisProblem$Metal[ThisProblem$Metal$Name %in% ThisProblem$Comp$Name[ComponentToRemove] == FALSE, , drop = FALSE]
  NewProblem$Metal$CompR = match(NewProblem$Metal$Name, NewProblem$Comp$Name)
  NewProblem$N["Metal"] = length(NewProblem$Metal$Name)

  # Update SpecCompList
  # NewProblem$Spec$NC = rowSums(NewProblem$SpecStoich != 0L)
  NewProblem$SpecCompList = matrix(
    apply(
      NewProblem$SpecStoich,
      MARGIN = 1,
      FUN = function(X) {
        Tmp = sort(which(X != 0L))
        if (length(Tmp) < max(NewProblem$Spec$NC)) {
          Tmp = c(Tmp, rep(0L, max(NewProblem$Spec$NC) - length(Tmp)))
        }
        return(Tmp)
      }
    ),
    nrow = NewProblem$N["Spec"],
    ncol = max(NewProblem$Spec$NC),
    byrow = TRUE,
    dimnames = list(NewProblem$Spec$Name, NULL)
  )

  # Update PhaseCompList
  # NewProblem$Phase$NC = rowSums(NewProblem$PhaseStoich != 0L)
  if (NewProblem$N["Phase"] > 0) {
    NewProblem$PhaseCompList = matrix(
      apply(
        NewProblem$PhaseStoich,
        MARGIN = 1,
        FUN = function(X) {
          Tmp = sort(which(X != 0L))
          if (length(Tmp) < max(NewProblem$Phase$NC)) {
            Tmp = c(Tmp, rep(0L, max(NewProblem$Phase$NC) - length(Tmp)))
          }
          return(Tmp)
        }
      ),
      nrow = NewProblem$N["Phase"],
      ncol = max(NewProblem$Phase$NC),
      byrow = TRUE,
      dimnames = list(NewProblem$Phase$Name, NULL)
    )
  } else {
    NewProblem$PhaseCompList = matrix(data = 0L, nrow = 0, ncol = 0)
  }

  if (DoCheck) {
    CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  }

  return(NewProblem)
}


#' @rdname Components
#' @export
AddInComps = function(ThisProblem, InCompName, InCompCharge, InCompMCName = NULL,
                      InCompType, InCompActCorr,
                      InCompMCR = match(InCompMCName, ThisProblem$Mass$Name, nomatch = -1L),
                      DoCheck = TRUE) {


  if (DoCheck) {
    CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  }

  # Add the component
  NewProblem = AddComponents(ThisProblem = ThisProblem,
                             CompName = InCompName,
                             CompCharge = InCompCharge,
                             CompMCName = InCompMCName,
                             CompType = InCompType,
                             CompActCorr = InCompActCorr,
                             CompSiteDens = 1.0,
                             CompMCR = InCompMCR,
                             DoCheck = FALSE)

  # Add the input component
  NewProblem$N["InComp"] = NewProblem$N["InComp"] + length(InCompName)
  NewProblem$InCompName = c(NewProblem$InCompName, trimws(as.character(InCompName)))

  if (DoCheck) {
    CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  }

  return(NewProblem)

}


#' @rdname Components
#' @export
RemoveInComps = function(ThisProblem, InCompToRemove, DoCheck = TRUE) {

  if (DoCheck) {
    CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  }

  InCompToRemoveOrig = InCompToRemove
  if (is.character(InCompToRemove)) {
    InCompToRemove = match(InCompToRemove, ThisProblem$InCompName)
  }
  if (any(is.na(InCompToRemove))) {
    stop(paste0("Input component(s) (",
                paste(InCompToRemoveOrig[is.na(InCompToRemove)],
                      collapse = ", "),
                ") does not exist."))
  }

  NewProblem = ThisProblem

  NewProblem = RemoveComponents(ThisProblem = NewProblem,
                                ComponentToRemove = InCompToRemoveOrig,
                                DoCheck = FALSE)

  # NewProblem$InCompName = NewProblem$InCompName[-InCompToRemove]
  # RemoveComponents should already do this...

  if (DoCheck) {
    CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  }

  return(NewProblem)

}


#' @rdname Components
#' @export
AddDefComps = function(ThisProblem, DefCompName, DefCompFromNum = NULL,
                       DefCompFromVar = NULL, DefCompCharge, DefCompMCName = NULL,
                       DefCompType, DefCompActCorr, DefCompSiteDens,
                       DefCompMCR = match(DefCompMCName, ThisProblem$Mass$Name, nomatch = -1L),
                       InDefComp = TRUE, DoCheck = TRUE) {

  if (DoCheck) {
    CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  }
  NewProblem = ThisProblem

  # coerce to correct types
  DefCompName = trimws(as.character(DefCompName))
  DefCompFromNum = as.numeric(DefCompFromNum)
  DefCompFromVar = trimws(as.character(DefCompFromVar))
  DefCompCharge = as.integer(DefCompCharge)
  DefCompMCName = trimws(as.character(DefCompMCName))
  DefCompMCR = as.integer(DefCompMCR)
  DefCompType = trimws(as.character(DefCompType))
  DefCompActCorr = trimws(as.character(DefCompActCorr))
  DefCompSiteDens = as.numeric(DefCompSiteDens)

  # error checking
  if (any(DefCompName %in% ThisProblem$Spec$Name)) {
    stop(
      "Input Component(s) (",
      paste(paste0("\"", DefCompName[DefCompName %in% ThisProblem$Spec$Name], "\""),
            collapse = ", "),
      ") already exist as a component or species."
    )
  }
  if (any(is.na(c(DefCompName, DefCompCharge, DefCompMCName,
                  DefCompType, DefCompActCorr, DefCompMCR, DefCompSiteDens)))) {
    stop("NA arguments not allowed, except for DefCompFromNum and DefCompFromVar.")
  }
  if (!all(DefCompMCName %in% ThisProblem$Mass$Name) ||
      !all(DefCompMCR <= ThisProblem$N["Mass"])) {
    stop("Mass compartment(s) specified in DefCompMCName or DefCompMCR does not exist.")
  }
  if (length(DefCompFromVar) == 0){
    DefCompFromVar = array(NA, dim = length(DefCompName),
                           dimnames = list(DefCompName))
  }
  if (length(DefCompFromNum) == 0){
    DefCompFromNum = array(NA, dim = length(DefCompName),
                           dimnames = list(DefCompName))
  }
  if (all(is.na(DefCompFromVar) & is.na(DefCompFromNum))) {
    stop("Either DefCompFromNum or DefCompFromVar must be specified. ",
         "Use DefCompFromNum = 1 to specify the total concentration = DefCompSiteDens.")
  }
  if (any(DefCompFromVar %in%
          c(ThisProblem$Comp$Name, ThisProblem$InVar$Name, NA, "KW/H",
            paste0(ThisProblem$InVar$Name[ThisProblem$InVar$Type %in% c(
              "WHAM-HA","WHAM-FA","WHAM-HAFA")], "-",c("HA","FA"),"_")) == FALSE)) {
    stop("DefCompFromVar must be an input variable or component.")
  }
  if(is.null(DefCompMCName)) {
    DefCompMCName = ThisProblem$Mass$Name[DefCompMCR]
  }
  if (!all(DefCompMCName %in% ThisProblem$Mass$Name[DefCompMCR])) {
    stop("DefCompMCName does not match DefCompMCR.")
  }
  NDefCompAdd = length(DefCompName)

  # Add the defined component
  NewProblem$N["DefComp"] = ThisProblem$N["DefComp"] + NDefCompAdd
  NewProblem$DefComp = rbind(ThisProblem$DefComp,
                             data.frame(
                               Name = trimws(as.character(DefCompName)),
                               FromNum = as.numeric(DefCompFromNum),
                               FromVar = trimws(as.character(DefCompFromVar)),
                               Charge = as.integer(DefCompCharge),
                               MCName = trimws(as.character(DefCompMCName)),
                               MCR = as.integer(DefCompMCR),
                               Type = trimws(as.character(DefCompType)),
                               ActCorr = trimws(as.character(DefCompActCorr)),
                               SiteDens = as.numeric(DefCompSiteDens)
                             ))
  if (any(InDefComp)) {
    if (length(InDefComp) == 1) { InDefComp = rep(InDefComp, NDefCompAdd) }
    NewProblem$InDefCompName = c(NewProblem$InDefCompName, DefCompName[InDefComp])
    NewProblem$N["InDefComp"] = NewProblem$N["InDefComp"] + as.integer(sum(InDefComp))
  }

  # Add the component
  NewProblem = AddComponents(ThisProblem = NewProblem,
                             CompName = DefCompName,
                             CompCharge = DefCompCharge,
                             CompMCName = DefCompMCName,
                             CompType = DefCompType,
                             CompActCorr = DefCompActCorr,
                             CompSiteDens = DefCompSiteDens,
                             CompMCR = DefCompMCR,
                             DoCheck = FALSE)

  if (DoCheck) {
    CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  }

  return(NewProblem)

}


#' @rdname Components
#' @export
RemoveDefComps = function(ThisProblem, DefCompToRemove, DoCheck = TRUE) {

  if (DoCheck) {
    CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  }
  NewProblem = ThisProblem

  DefCompToRemoveOrig = DefCompToRemove
  if (is.character(DefCompToRemove)) {
    DefCompToRemove = which(ThisProblem$DefComp$Name %in% DefCompToRemove)
  }
  if (length(DefCompToRemove) < 1) {
    stop(paste0("Defined Component \"", DefCompToRemoveOrig, "\" does not exist."))
  }
  if (any(DefCompToRemove <= 0L)) {
    stop("Invalid index in DefCompToRemove (",
         DefCompToRemove[DefCompToRemove <= 0L],").")
  }
  if (any(DefCompToRemove > ThisProblem$N["DefComp"])) {
    stop(paste0("There are ", ThisProblem$N["DefComp"], " Defined Components, ",
                "trying to remove the #(",
                paste(DefCompToRemove[DefCompToRemove > ThisProblem$N["DefComp"]],
                      collapse = ", "),
                ") element(s)."))
  }

  # Remove components that are the DefComp
  ComponentToRemove = which(ThisProblem$Comp$Name %in%
                              ThisProblem$DefComp$Name[DefCompToRemove])
  NewProblem = RemoveComponents(ThisProblem = NewProblem,
                                ComponentToRemove = ComponentToRemove,
                                DoCheck = FALSE)

  # Remove InDefComp's that are the defined component to remove
  InDefCompToRemove = which(NewProblem$InDefCompName %in%
                              ThisProblem$DefComp$Name[DefCompToRemove])
  if (length(InDefCompToRemove) >= 1) {
    NewProblem$InDefCompName = NewProblem$InDefCompName[-InDefCompToRemove]
    NewProblem$N["InDefComp"] = ThisProblem$N["InDefComp"] - length(InDefCompToRemove)
  }


  # Remove DefComps
  NewProblem$DefComp = ThisProblem$DefComp[-DefCompToRemove, , drop= FALSE]
  NewProblem$N["DefComp"] = ThisProblem$N["DefComp"] - length(DefCompToRemove)

  if (DoCheck) {
    CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  }

  return(NewProblem)

}
