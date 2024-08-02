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
#'   "PercHA", and "PercAFA".
#' @param InVarMCR (optional) A character vector with the indices of the mass
#'   compartments the new input variables are associated with. Only needs to be
#'   specified if `InVarMCName` is not specified.
#' @param InVarToRemove A character vector with names or indices of the input
#'   variable(s) to remove from `ThisProblem`.
#'
#' @return `ThisProblem`, with the changed input variables. If the input
#'   variable being added is pH, an "H" component will also be added as a fixed
#'   activity component.
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
AddInVars = function(ThisProblem, InVarName, InVarMCName = NULL,
                    InVarType,
                    InVarMCR = match(InVarMCName, ThisProblem$Mass$Name)) {

  CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  NewProblem = ThisProblem

  if ((NewProblem$ParamFile != "") &&
      !grepl("[(]modified[)]$", NewProblem$ParamFile)) {
    NewProblem$ParamFile = paste0(NewProblem$ParamFile, " (modified)")
  }

  # error checking
  if (any((InVarName %in% ThisProblem$InVar$Name))) {
    stop(paste0(
      "Input Variable(s) (",
      paste(paste0("\"", InVarName[InVarName %in% ThisProblem$InVar$Name], "\""), collapse = ", "),
      ") already exist."
    ))
  }
  if (!all(InVarMCName %in% ThisProblem$Mass$Name) ||
      !all(InVarMCR <= ThisProblem$N["Mass"])) {
    stop("Mass compartment(s) specified in InVarMCName or InVarMCR does not exist.")
  }
  if (any(is.na(c(InVarName, InVarType, InVarMCName, InVarMCR)))) {
    stop("NA arguments not allowed.")
  }
  if(is.null(InVarMCName)) {
    InVarMCName = ThisProblem$Mass$Name[InVarMCR]
  }
  if (!all(InVarMCName %in% ThisProblem$Mass$Name[InVarMCR])) {
    stop("InVarMCName does not match InVarMCR.")
  }
  stopifnot(all(InVarType %in%
                  c("Temperature", "pH", "WHAM-FA", "WHAM-HA",
                    "WHAM-HAFA", "PercHA", "PercAFA")))
  NInVarAdd = length(InVarName)
  if (length(InVarMCName) == 1) { InVarMCName = rep(InVarMCName, NInVarAdd)}
  if (length(InVarMCR) == 1) { InVarMCR = rep(InVarMCR, NInVarAdd)}
  if (length(InVarType) == 1) { InVarType = rep(InVarType, NInVarAdd)}

  # Add input variable
  NewProblem$InVar = rbind(ThisProblem$InVar,
                           data.frame(Name = trimws(as.character(InVarName)),
                                      MCName = trimws(as.character(InVarMCName)),
                                      MCR = as.integer(InVarMCR),
                                      Type = trimws(as.character(InVarType))))
  NewProblem$N["InVar"] = ThisProblem$N["InVar"] + NInVarAdd

  for (i in 1:NInVarAdd) {
    if (InVarType[i] == "pH") {
      NewProblem = AddDefComps(ThisProblem = NewProblem,
                              DefCompName = "H",
                              DefCompFromVar = InVarName[i],
                              DefCompFromNum = NA,
                              DefCompCharge = 1L,
                              DefCompMCName = InVarMCName[i],
                              DefCompType = "FixedAct",
                              DefCompActCorr = "Debye",
                              DefCompSiteDens = 1.0,
                              DefCompMCR = InVarMCR[i],
                              InDefComp = FALSE)
      # warning(paste("Defined component added for pH as a fixed activity",
      #               "component. Change manually if you wish to represent",
      #               "pH as a fixed concentration."))
    }
  }

  CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  return(NewProblem)

}


#' @rdname InVars
#' @export
RemoveInVars = function(ThisProblem, InVarToRemove) {

  CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  NewProblem = ThisProblem

  if ((NewProblem$ParamFile != "") &&
      !grepl("[(]modified[)]$", NewProblem$ParamFile)) {
    NewProblem$ParamFile = paste0(NewProblem$ParamFile, " (modified)")
  }

  if (length(InVarToRemove) < 1) {
    stop("Specify an input variable to add.")
  }

  InVarToRemoveOrig = InVarToRemove
  if (is.character(InVarToRemove)) {
    InVarToRemove = match(InVarToRemove, ThisProblem$InVar$Name)
  }
  if (any(is.na(InVarToRemove))) {
    stop(paste0("Input Variable(s) (",
                paste(InVarToRemoveOrig[is.na(InVarToRemove)], collapse = ", "),
                ") do not exist."))
  }
  if (any(InVarToRemove > ThisProblem$N["InVar"])) {
    stop(paste0("There are ", ThisProblem$N["InVar"], " input variables, ",
                "trying to remove the #(",
                InVarToRemove[InVarToRemove > ThisProblem$N["InVar"]],
                ") element(s)."))
  }

  # Remove DefComps that depend on input variables
  DefCompToRemove = which(NewProblem$DefComp$FromVar %in%
                            NewProblem$InVar$Name[InVarToRemove])
  if (length(DefCompToRemove) >= 1) {
    NewProblem = RemoveDefComps(NewProblem, DefCompToRemove)
  }

  # Remove input variables
  NewProblem$InVar = ThisProblem$InVar[-InVarToRemove, , drop = FALSE]
  rownames(NewProblem$InVar) = NULL
  NewProblem$N["InVar"] = ThisProblem$N["InVar"] - length(InVarToRemove)

  CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  return(NewProblem)

}
