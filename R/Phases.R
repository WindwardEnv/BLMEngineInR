#' @name Phases
#'
#' @title Add or remove phase reactions in a problem
#'
#' @description PHASES ARE NOT CURRENTLY IMPLEMENTED. This function is here for
#'   as a placeholder since it will require much of the same support
#'   infrastructure once it is implemented, but no reactions are processed in
#'   CHESS.
#'
#' @param ThisProblem A list object with a structure like that returned by
#'   `BlankProblem()`.
#' @param PhaseEquation A character vector giving the chemical equation for a
#'   formation reaction. This must include the stoichiometric coefficients for
#'   each reactant, even if it's 1. (e.g., the equation for the formation of
#'   calcium chloride would be `"CaCl2 = 1 * Ca + 2 * Cl"`). If `PhaseName` is
#'   also supplied, then a partial equation with just the right hand side
#'   (reactants) can be supplied (i.e., `"= 1 * Ca + 2 * Cl"`). Can be omitted
#'   if either `PhaseStoich` or both `PhaseCompNames` and `PhaseCompStoichs` are
#'   supplied.
#' @param PhaseName A character vector with the name(s) of the species to add
#'   formation reactions for. Can be omitted if `SpecEquation` indicates the
#'   phase name.
#' @param PhaseCompNames A list where each element is a character vector of the
#'   component names used to form each phase. See examples for clarification.
#'   Can be omitted if `PhaseEquation` or `PhaseStoich` is supplied.
#' @param PhaseCompStoichs A list where each element is an integer vector of the
#'   stoichiometric coefficients of each component used to form each phase. See
#'   examples for clarification. Can be omitted if `PhaseEquation` or
#'   `PhaseStoich` is supplied.
#' @param PhaseStoich A matrix of stoichiometric coefficients, where each row
#'   corresponds to a phase reaction and each column corresponds to a component.
#'   The columns should match `ThisProblem$Comp$Name` exactly. Can be omitted if
#'   either `PhaseEquation` or both `PhaseCompNames` and `PhaseCompStoichs` are
#'   supplied.
#' @param PhaseLogK A numeric vector with the log10-transformed equilibrium
#'   coefficients of the phase formation reactions.
#' @param PhaseDeltaH A numeric vector with the change in enthalpy of the phase
#'   formation reactions.
#' @param PhaseTempKelvin A numeric vector with the temperatures (in Kelvin)
#'   corresponding to `PhaseDeltaH` values of the phase formation reactions.
#' @param PhaseMoles A numeric vector with the moles of the phase.
#' @param PhasesToRemove A character or integer vector indicating the names or
#'   indices (respectively) of the phase formation reactions to remove.
#' @param DoCheck A logical value indicating whether checks should be performed
#'   on the incoming and outgoing problem objects. Defaults to `TRUE`, as you
#'   usually want to make sure something isn't awry, but the value is often set
#'   to `FALSE` when used internally (like in DefineProblem) so the problem is
#'   only checked once at the end.
#'
#' @return `ThisProblem`, with the phase reaction(s) changed.
#'
#' @family problem manipulation functions
#'
#' @examples
#' print(carbonate_system_problem$Phase)
#' my_new_problem = carbonate_system_problem
#' my_new_problem = AddPhases(ThisProblem = my_new_problem,
#'                            PhaseEquation = "CO2(g) = 1 * CO3 + 2 * H",
#'                            PhaseLogK = -1.5,
#'                            PhaseDeltaH = 0,
#'                            PhaseTempKelvin = 0,
#'                            PhaseMoles = 10^-3.5)
#' print(my_new_problem$Phase)
#'
NULL


#' @rdname Phases
#' @export
AddPhases = function(ThisProblem,
                     PhaseEquation = character(),
                     PhaseName = character(),
                     PhaseCompNames = list(), PhaseCompStoichs = list(),
                     PhaseStoich = NULL,
                     PhaseLogK, PhaseDeltaH, PhaseTempKelvin, PhaseMoles,
                     DoCheck = TRUE) {

  if (DoCheck) {
    CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  }
  NewProblem = ThisProblem

  if ((NewProblem$ParamFile != "") &&
      !grepl("[(]modified[)]$", NewProblem$ParamFile)) {
    NewProblem$ParamFile = paste0(NewProblem$ParamFile, " (modified)")
  }

  HasName = length(PhaseName) > 0
  HasStoichCompsNames = (length(PhaseCompNames) > 0) && (length(PhaseCompStoichs) > 0)
  HasStoichMatrix = !is.null(PhaseStoich)
  HasEquation = length(PhaseEquation) > 0
  if (!((HasName && (HasStoichMatrix || HasStoichCompsNames)) | HasEquation)) {
    stop(c(
      "Inputs differ in length.  Must specify either: ",
      "  (1) PhaseEquation (and optionally PhaseName); ",
      "  (2) PhaseName, PhaseCompNames, and PhaseCompStoichs; or ",
      "  (3) PhaseName and PhaseStoich   for each reaction."
    ))
  }

  if (HasStoichCompsNames) {
    tmp = StoichCompsToStoichMatrix(CompName = ThisProblem$Comp$Name,
                                    SpecName = PhaseName,
                                    SpecCompNames = PhaseCompNames,
                                    SpecCompStoichs = PhaseCompStoichs)
    if (HasStoichMatrix) {
      stopifnot(all(tmp == PhaseStoich))
    } else {
      PhaseStoich = tmp
      HasStoichMatrix = TRUE
    }
  }
  if (HasEquation) {
    tmp = EquationToStoich(PhaseEquation, CompName = ThisProblem$Comp$Name)
    if (!HasName) {
      PhaseName = tmp$SpecName
      HasName = TRUE
    } else if (all(trimws(gsub("=.*", "", PhaseEquation)) == "")) {
      PhaseEquation = paste(PhaseName, "=", PhaseEquation)
    }
    if (HasStoichCompsNames) {
      for (i in 1:length(PhaseCompNames)) {
        stopifnot(PhaseCompNames[[i]] == tmp$SpecCompNames[[i]],
                  PhaseCompStoichs[[i]] == tmp$SpecCompStoichs[[i]])
      }
    } else {
      PhaseCompNames = tmp$SpecCompNames
      PhaseCompStoichs = tmp$SpecCompStoichs
      HasStoichCompsNames = TRUE
    }
    if (HasStoichMatrix) {
      stopifnot(PhaseStoich == tmp$SpecStoich)
    } else {
      PhaseStoich = tmp$SpecStoich
      HasStoichMatrix = TRUE
    }
  } else {
    PhaseEquation = StoichMatrixToEquation(SpecStoich = PhaseStoich,
                                     SpecName = PhaseName,
                                     CompName = ThisProblem$Comp$Name)
    HasEquation = TRUE
  }
  NPhaseAdd = length(PhaseName)
  PhaseNC = as.integer(rowSums(PhaseStoich != 0L))

  # error checking
  if (any(PhaseName %in% ThisProblem$Phase$Name)) {
    stop(paste0(
      "Phases (",
      paste(paste0("\"", PhaseName[PhaseName %in% ThisProblem$Phase$Name], "\""),
            collapse = ", "),
      ") already exist as a component or Phases."
    ))
  }
  if (any(is.na(c(PhaseName, PhaseLogK, PhaseDeltaH, PhaseTempKelvin)))) {
    stop("NA arguments not allowed.")
  }
  if (any(PhaseNC != sapply(PhaseCompStoichs, FUN = length))) {
    stop("Stoichiometric inputs differ in length. Specify a matching set of names and values in PhaseCompNames and PhaseCompStoichs.")
  }

  # add Phases
  NewProblem$Phase = rbind(ThisProblem$Phase,
                           data.frame(
                             Name = trimws(as.character(PhaseName)),
                             Equation = trimws(as.character(PhaseEquation)),
                             LogK = as.numeric(PhaseLogK),
                             K = 10^as.numeric(PhaseLogK),
                             DeltaH = as.numeric(PhaseDeltaH),
                             TempKelvin = as.numeric(PhaseTempKelvin),
                             NC = as.integer(PhaseNC),
                             Moles = as.numeric(PhaseMoles)
                           ))

  PhaseCompList = matrix(
    apply(
      PhaseStoich,
      MARGIN = 1,
      FUN = function(X) {
        Tmp = sort(which(X != 0L))
        if (length(Tmp) < max(NewProblem$Phase$NC)) {
          Tmp = c(Tmp, rep(0, max(NewProblem$Phase$NC) - length(Tmp)))
        }
        return(Tmp)
      }
    ),
    nrow = NPhaseAdd,
    ncol = max(NewProblem$Phase$NC),
    byrow = TRUE,
    dimnames = list(PhaseName, NULL)
  )
  if (ncol(PhaseCompList) > ncol(NewProblem$PhaseCompList)) {
    while (ncol(PhaseCompList) > ncol(NewProblem$PhaseCompList)) {
      NewProblem$PhaseCompList = cbind(
        NewProblem$PhaseCompList,
        matrix(0L, nrow = nrow(NewProblem$PhaseCompList), ncol = 1)
      )
    }
  }
  NewProblem$PhaseCompList = rbind(
    NewProblem$PhaseCompList,
    PhaseCompList
  )

  NewProblem$PhaseStoich = rbind(
    ThisProblem$PhaseStoich,
    PhaseStoich
  )
  NewProblem$N["Phase"] = ThisProblem$N["Phase"] + NPhaseAdd

  if (DoCheck) {
    CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  }
  return(NewProblem)

}

#' @rdname Phases
#' @export
RemovePhases = function(ThisProblem, PhasesToRemove, DoCheck = TRUE) {

  if (DoCheck) {
    CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  }
  NewProblem = ThisProblem

  if ((NewProblem$ParamFile != "") &&
      !grepl("[(]modified[)]$", NewProblem$ParamFile)) {
    NewProblem$ParamFile = paste0(NewProblem$ParamFile, " (modified)")
  }

  PhasesToRemoveOrig = PhasesToRemove
  if (is.character(PhasesToRemove)) {
    PhasesToRemove = which(ThisProblem$Phase$Name %in% PhasesToRemove)
  }
  if (length(PhasesToRemove) < 1) {
    stop(paste0("Phases \"", PhasesToRemoveOrig, "\" does not exist."))
  }
  if (any(PhasesToRemove <= 0L)) {
    stop("Invalid index in PhasesToRemove (",
         PhasesToRemove[PhasesToRemove <= 0L],").")
  }
  if (any(PhasesToRemove > ThisProblem$N["Phase"])) {
    stop(paste0("There are ", ThisProblem$N["Phase"], " Phases, ",
                "trying to remove the #(",
                paste(PhasesToRemove[PhasesToRemove > ThisProblem$N["Phase"]],
                      collapse = ", "),
                ") element(s)."))
  }

  # Remove Phases
  NewProblem$Phase = ThisProblem$Phase[-PhasesToRemove, , drop = FALSE]
  NewProblem$PhaseCompList = ThisProblem$PhaseCompList[-PhasesToRemove, , drop = FALSE]
  NewProblem$PhaseStoich = ThisProblem$PhaseStoich[-PhasesToRemove, , drop = FALSE]
  NewProblem$N["Phase"] = ThisProblem$N["Phase"] - length(PhasesToRemove)

  if (DoCheck) {
    CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  }
  return(NewProblem)

}
