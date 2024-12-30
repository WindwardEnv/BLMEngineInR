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

#' @name Species
#'
#' @title Add or remove a species reactions in a problem
#'
#' @description Functions to add or remove species formation reactions.
#'
#' @param ThisProblem A list object with a structure like that returned by
#'   `BlankProblem()`.
#' @param SpecEquation A character vector giving the chemical equation for a
#'   formation reaction. This must include the stoichiometric coefficients for
#'   each reactant, even if it's 1. (e.g., the equation for the formation of
#'   calcium chloride would be `"CaCl2 = 1 * Ca + 2 * Cl"`). If `SpecName` is
#'   also supplied, then a partial equation with just the right hand side
#'   (reactants) can be supplied (i.e., `"= 1 * Ca + 2 * Cl"`). Can be omitted
#'   if either `SpecStoich` or both `SpecCompNames` and `SpecCompStoichs` are
#'   supplied.
#' @param SpecName A character vector with the name(s) of the species to add
#'   formation reactions for. Can be omitted if `SpecEquation` indicates the
#'   species name.
#' @param SpecMCName A character vector with the name(s) of the mass
#'   compartments the new species are associated with. Does not need to be
#'   specified if `SpecMCR` is specified instead.
#' @param SpecType A character vector with the species type. SpecType values
#'   must be either `"Normal"`, `"DonnanHA"`, `"DonnanFA"`, `"WHAMHA"`,
#'   `"WHAMFA"`. The default value is `"Normal"`, while the others are usually
#'   only needed for indicating species that are added from `ExpandWHAM`.
#' @param SpecActCorr A character vector with the activity correction method(s)
#'   of the new species. Must be one of "None", "Debye", "Davies", "DonnanHA",
#'   "DonnanFA", "WHAMHA", or "WHAMFA".
#' @param SpecCompNames A list where each element is a character vector of the
#'   component names used to form each species. See examples for clarification.
#'   Can be omitted if `SpecEquation` or `SpecStoich` is supplied.
#' @param SpecCompStoichs A list where each element is an integer vector of the
#'   stoichiometric coefficients of each component used to form each species.
#'   See examples for clarification. Can be omitted if `SpecEquation` or
#'   `SpecStoich` is supplied.
#' @param SpecStoich A matrix of stoichiometric coefficients, where each row
#'   corresponds to a chemical species and each column corresponds to a
#'   component. The columns should match `ThisProblem$Comp$Name` exactly. Can be
#'   omitted if either `SpecEquation` or both `SpecCompNames` and
#'   `SpecCompStoichs` are supplied.
#' @param SpecLogK A numeric vector with the log10-transformed equilibrium
#'   coefficients of the species formation reactions.
#' @param SpecDeltaH A numeric vector with the change in enthalpy of the species
#'   formation reactions.
#' @param SpecTempKelvin A numeric vector with the temperatures (in Kelvin)
#'   corresponding to `SpecDeltaH` values of the species formation reactions.
#' @param SpecMCR (optional) A character vector with the indices of the mass
#'   compartments the new species are associated with. Only needs to be
#'   specified if `SpecMCName` is not specified.
#' @param InSpec A logical value indicating if this is a species formation
#'   reaction indicated from the parameter file (`TRUE`, the default) or a
#'   reaction that was added from another process such as `ExpandWHAM`
#'   (`FALSE`). This should usually only be `FALSE` when another function is
#'   calling this function, such as `ExpandWHAM`.
#' @param SpeciesToRemove A character or integer vector indicating the names or
#'   indices (respectively) of the species formation reactions to remove.
#' @param DoCheck A logical value indicating whether checks should be performed
#'   on the incoming and outgoing problem objects. Defaults to `TRUE`, as you
#'   usually want to make sure something isn't awry, but the value is often set
#'   to `FALSE` when used internally (like in DefineProblem) so the problem is
#'   only checked once at the end.
#'
#' @return `ThisProblem`, with the species reaction(s) changed.
#'
#' @family problem manipulation functions
#'
#' @examples
#' print(carbonate_system_problem$Spec)
#' my_new_problem = carbonate_system_problem
#' my_new_problem = AddInComps(ThisProblem = my_new_problem,
#'                             InCompName = "Ca",
#'                             InCompCharge = 2,
#'                             InCompMCName = "Water",
#'                             InCompType = "MassBal",
#'                             InCompActCorr = "Debye")
#' my_new_problem = AddSpecies(ThisProblem = my_new_problem,
#'                             SpecEquation = c("CaCO3 = 1 * Ca + 1 * CO3",
#'                                              "CaHCO3 = 1 * Ca + 1 * H + 1 * CO3"),
#'                             SpecMCName = "Water",
#'                             SpecActCorr = "Debye",
#'                             SpecLogK = c(3.22, 11.44),
#'                             SpecDeltaH = c(14951, -3664),
#'                             SpecTempKelvin = 298.15)
#' print(my_new_problem$Spec)
#' my_new_problem = RemoveSpecies(ThisProblem = my_new_problem,
#'                                SpeciesToRemove = "CaCO3")
#' print(my_new_problem$Spec)
#'
NULL

#' @rdname Species
#' @export
AddSpecies = function(ThisProblem,
                      SpecEquation = character(),
                      SpecName = character(),
                      SpecMCName = NULL, SpecType = "Normal", SpecActCorr,
                      SpecCompNames = list(), SpecCompStoichs = list(),
                      SpecStoich = NULL,
                      SpecLogK, SpecDeltaH, SpecTempKelvin,
                      SpecMCR = match(SpecMCName, ThisProblem$Mass$Name, nomatch = -1L),
                      InSpec = TRUE, DoCheck = TRUE) {

  if (DoCheck) {
    CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  }
  NewProblem = ThisProblem

  if ((NewProblem$ParamFile != "") &&
      !grepl("[(]modified[)]$", NewProblem$ParamFile)) {
    NewProblem$ParamFile = paste0(NewProblem$ParamFile, " (modified)")
  }

  HasName = length(SpecName) > 0
  HasStoichCompsNames = (length(SpecCompNames) > 0) && (length(SpecCompStoichs) > 0)
  HasStoichMatrix = !is.null(SpecStoich)
  HasEquation = length(SpecEquation) > 0
  if (!((HasName && (HasStoichMatrix || HasStoichCompsNames)) || HasEquation)) {
    stop("Missing stoichiometry and/or species name information. ",
         "Must specify either: ",
         "(1) SpecEquation (and optionally SpecName); ",
         "(2) SpecName, SpecCompNames, and SpecCompStoichs; or ",
         "(3) SpecName and SpecStoich for each reaction.")
  }

  if (HasStoichMatrix) {
    mode(SpecStoich) = "integer"
    if (!is.null(colnames(SpecStoich))){
      SpecStoich = SpecStoich[, ThisProblem$Comp$Name, drop = FALSE]
      SpecStoich[is.na(SpecStoich)] = 0L
    }
  }
  if (HasStoichCompsNames) {
    tmp = StoichCompsToStoichMatrix(CompName = ThisProblem$Comp$Name,
                                    SpecName = SpecName,
                                    SpecCompNames = SpecCompNames,
                                    SpecCompStoichs = SpecCompStoichs)
    if (HasStoichMatrix) {
      SpecStoich = SpecStoich[, ThisProblem$Comp$Name, drop = FALSE]
      stopifnot(all(tmp == SpecStoich))
    } else {
      SpecStoich = tmp
      HasStoichMatrix = TRUE
    }
  }
  if (HasEquation) {

    if (HasName && (all(grepl("^=", trimws(SpecEquation)) ||
                        !grepl("=", SpecEquation)))) {
      SpecEquation = paste(SpecName, "=", trimws(gsub("=","",SpecEquation)))
    }

    tmp = EquationToStoich(SpecEquation, CompName = ThisProblem$Comp$Name)
    if (HasName) {
      stopifnot(SpecName == tmp$SpecName)
    } else {
      SpecName = tmp$SpecName
      HasName = TRUE
    }
    if (HasStoichCompsNames) {
      for (i in 1:length(SpecCompNames)) {
        SpecCompsAgg = stats::aggregate(SpecCompStoichs[[i]],
                                        by = list(SpecCompNames[[i]]), FUN = sum)
        NewOrder = stats::na.omit(match(ThisProblem$Comp$Name, SpecCompsAgg$Group.1))
        SpecCompNames[[i]] = SpecCompsAgg$Group.1[NewOrder]
        SpecCompStoichs[[i]] = SpecCompsAgg$x[NewOrder]
        stopifnot(SpecCompNames[[i]] == tmp$SpecCompNames[[i]],
                  SpecCompStoichs[[i]] == tmp$SpecCompStoichs[[i]])
      }
    } else {
      SpecCompNames = tmp$SpecCompNames
      SpecCompStoichs = tmp$SpecCompStoichs
      HasStoichCompsNames = TRUE
    }
    if (HasStoichMatrix) {
      stopifnot(SpecStoich == tmp$SpecStoich)
    } else {
      SpecStoich = tmp$SpecStoich
      HasStoichMatrix = TRUE
    }
    SpecEquation = StoichMatrixToEquation(SpecStoich = SpecStoich,
                                          SpecName = SpecName,
                                          CompName = ThisProblem$Comp$Name)
  } else {
    SpecEquation = StoichMatrixToEquation(SpecStoich = SpecStoich,
                                          SpecName = SpecName,
                                          CompName = ThisProblem$Comp$Name)
    HasEquation = TRUE
  }
  NSpecAdd = length(SpecName)
  SpecNC = as.integer(rowSums(SpecStoich != 0L))

  # error checking
  if (any(SpecName %in% ThisProblem$Spec$Name)) {
    stop(paste0(
      "Species (",
      paste(paste0("\"", SpecName[SpecName %in% ThisProblem$Spec$Name], "\""),
            collapse = ", "),
      ") already exist as a component or species."
    ))
  }
  if (!all(SpecMCName %in% ThisProblem$Mass$Name) ||
      !all(SpecMCR <= ThisProblem$N["Mass"])) {
    stop("Mass compartment(s) specified in SpecMCName or SpecMCR does not exist.")
  }
  if (any(is.na(c(SpecName, SpecMCName,
                  SpecLogK, SpecDeltaH, SpecTempKelvin,
                  SpecType, SpecActCorr, SpecMCR)))) {
    stop("NA arguments not allowed.")
  }
  if(is.null(SpecMCName)) {
    SpecMCName = ThisProblem$Mass$Name[SpecMCR]
  }
  if (!all(SpecActCorr %in% c("None","Debye","Davies"))) {
    stop("SpecActCorr values must be either None, Debye, or Davies.")
  }
  if (!all(SpecType %in% c("Normal","DonnanHA","DonnanFA","WHAMHA","WHAMFA"))) {
    stop("SpecType values must be either Normal, DonnanHA, DonnanFA, WHAMHA, WHAMFA.")
  }
  if (!all(SpecMCName %in% ThisProblem$Mass$Name[SpecMCR])) {
    stop("SpecMCName does not match SpecMCR.")
  }
  if (any(SpecNC != sapply(SpecCompNames, FUN = function(X){sum(X != "")}))) {
    stop("Stoichiometric inputs differ in length. Specify a matching set of names and values in SpecCompNames and SpecCompStoichs.")
  }

  # add Species
  NewProblem$Spec = rbind(ThisProblem$Spec,
                          data.frame(
                            Name = trimws(as.character(SpecName)),
                            Equation = SpecEquation,
                            Charge = as.integer(NA),
                            MCName = trimws(as.character(SpecMCName)),
                            MCR = as.integer(SpecMCR),
                            Type = trimws(as.character(SpecType)),
                            ActCorr = trimws(as.character(SpecActCorr)),
                            LogK = as.numeric(SpecLogK),
                            K = 10^as.numeric(SpecLogK),
                            DeltaH = as.numeric(SpecDeltaH),
                            TempKelvin = as.numeric(SpecTempKelvin),
                            NC = as.integer(SpecNC)
                          ))
  if (any(InSpec)) {
    if (length(InSpec) == 1) { InSpec = rep(InSpec, NSpecAdd) }
    NewProblem$InSpecName = c(NewProblem$InSpecName, SpecName[InSpec])
    NewProblem$N["InSpec"] = NewProblem$N["InSpec"] + as.integer(sum(InSpec))
  }


  SpecCompList = matrix(
    apply(
      SpecStoich,
      MARGIN = 1,
      FUN = function(X) {
        Tmp = sort(which(X != 0L))
        if (length(Tmp) < max(NewProblem$Spec$NC)) {
          Tmp = c(Tmp, rep(0L, max(NewProblem$Spec$NC) - length(Tmp)))
        }
        return(Tmp)
      }
    ),
    nrow = NSpecAdd,
    ncol = max(NewProblem$Spec$NC),
    byrow = TRUE,
    dimnames = list(SpecName, NULL)
  )
  if (ncol(SpecCompList) > ncol(NewProblem$SpecCompList)) {
    while (ncol(SpecCompList) > ncol(NewProblem$SpecCompList)) {
      NewProblem$SpecCompList = cbind(
        NewProblem$SpecCompList,
        matrix(0L, nrow = nrow(NewProblem$SpecCompList), ncol = 1)
      )
    }
  }
  NewProblem$SpecCompList = rbind(
    NewProblem$SpecCompList,
    SpecCompList
  )

  NewProblem$SpecStoich = rbind(
    ThisProblem$SpecStoich,
    SpecStoich
  )
  NewProblem$N["Spec"] = ThisProblem$N["Spec"] + NSpecAdd

  # These parts it's easier to just calculate again
  NewProblem$Spec$Charge = as.integer(drop(NewProblem$SpecStoich %*% NewProblem$Comp$Charge))

  # Put components at the beginning
  CompIndexes = match(NewProblem$Comp$Name, NewProblem$Spec$Name)
  SpecIndexes = which(NewProblem$Spec$Name %in% NewProblem$Comp$Name == FALSE)
  for (i in names(NewProblem)[grepl("^Spec", names(NewProblem))]) {
    NewProblem[[i]] = rbind(NewProblem[[i]][CompIndexes, , drop = FALSE],
                            NewProblem[[i]][SpecIndexes, , drop = FALSE])
  }
  rownames(NewProblem$Spec) = NULL

  if (DoCheck) {
    CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  }
  return(NewProblem)

}


#' @rdname Species
#' @export
RemoveSpecies = function(ThisProblem, SpeciesToRemove, DoCheck = TRUE) {

  if (DoCheck) {
    CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  }
  NewProblem = ThisProblem

  if ((NewProblem$ParamFile != "") &&
      !grepl("[(]modified[)]$", NewProblem$ParamFile)) {
    NewProblem$ParamFile = paste0(NewProblem$ParamFile, " (modified)")
  }

  SpeciesToRemoveOrig = SpeciesToRemove
  if (is.character(SpeciesToRemove)) {
    SpeciesToRemove = which(ThisProblem$Spec$Name %in% SpeciesToRemove)
  }
  if (length(SpeciesToRemove) < 1) {
    stop(paste0("Species \"", SpeciesToRemoveOrig, "\" does not exist."))
  }
  if (any(SpeciesToRemove <= 0L)) {
    stop("Invalid index in SpeciesToRemove (",
         SpeciesToRemove[SpeciesToRemove <= 0L],").")
  }
  if (any(SpeciesToRemove > ThisProblem$N["Spec"])) {
    stop(paste0("There are ", ThisProblem$N["Spec"], " Species, ",
                "trying to remove the #(",
                paste(SpeciesToRemove[SpeciesToRemove > ThisProblem$N["Spec"]],
                      collapse = ", "),
                ") element(s)."))
  }

  # Remove Species
  NewProblem$Spec = ThisProblem$Spec[-SpeciesToRemove, , drop = FALSE]
  rownames(NewProblem$Spec) = NULL
  NewProblem$SpecCompList = ThisProblem$SpecCompList[-SpeciesToRemove, , drop = FALSE]
  NewProblem$SpecStoich = ThisProblem$SpecStoich[-SpeciesToRemove, , drop = FALSE]
  NewProblem$N["Spec"] = ThisProblem$N["Spec"] - length(SpeciesToRemove)

  # Remove InSpec's that are the Species to remove
  InSpecToRemove = which(NewProblem$InSpecName %in%
                           ThisProblem$Spec$Name[SpeciesToRemove])
  if (length(InSpecToRemove) >= 1) {
    NewProblem$InSpecName = NewProblem$InSpecName[-InSpecToRemove]
    NewProblem$N["InSpec"] = ThisProblem$N["InSpec"] - length(InSpecToRemove)
  }

  # Update Special Definitions
  NewProblem$BLMetal = ThisProblem$BLMetal[ThisProblem$BLMetal$Name %in% ThisProblem$Spec$Name[SpeciesToRemove] == FALSE, , drop = FALSE]
  NewProblem$BLMetal$SpecsR = match(NewProblem$BLMetal$Name, NewProblem$Spec$Name)
  NewProblem$N["BLMetal"] = length(NewProblem$BLMetal$Name)

  if (DoCheck) {
    CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  }
  return(NewProblem)

}
