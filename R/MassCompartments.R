#' @name MassCompartments
#'
#' @title Add or remove mass compartments in a problem
#'
#' @param ThisProblem A list object with a structure like that returned by
#'   `BlankProblem()`.
#' @param MassTable A `data.frame` object with, at a minimum, columns `Name`,
#'   `Amt`, and `Unit`, defining the characteristics of the mass compartment(s)
#'   to add.
#' @param MassName A character vector with the name(s) of the new mass
#'   compartment(s).
#' @param MassAmt A numeric vector with the mass compartment amount(s).
#' @param MassUnit A character vector with the units for the amount(s) of the
#'   mass compartment(s).
#' @param InMass A logical value or vector indicating if this mass compartment
#'   is in the parameter file (`TRUE`, default) or was created as a result of,
#'   e.g. the `ExpandWHAM` function (`FALSE`).
#' @param MCToRemove A character vector with names or indices of the mass
#'   compartment(s) to remove from `ThisProblem`.
#'
#' @return `ThisProblem`, with all the edited mass compartments, along with any
#'   components, input variables, etc. associated with those mass compartments
#'   edited.
#'
#' @inherit BlankProblem examples
#'
#' @family problem manipulation functions
NULL


#' @rdname MassCompartments
#' @export
AddMassCompartments = function(ThisProblem,
                               MassTable = data.frame(),
                               MassName = MassTable$Name,
                               MassAmt = MassTable$Amt,
                               MassUnit = MassTable$Unit,
                               InMass = TRUE) {

  NewProblem = ThisProblem

  # error checking
  if (any((MassName %in% ThisProblem$Mass$Name))) {
    stop(paste0(
      "Mass Compartment(s) (",
      paste(paste0("\"", MassName[MassName %in% ThisProblem$Mass$Name], "\""), collapse = ", "),
      ") already exist."
    ))
  }
  NMassAdd = length(MassName)
  if (any(is.na(c(MassName, MassAmt, MassUnit)))) {
    stop("NA arguments not allowed.")
  }

  # add mass compartment
  NewProblem$Mass = rbind(ThisProblem$Mass,
                          data.frame(Name = trimws(as.character(MassName)),
                                     Amt = as.numeric(MassAmt),
                                     Unit = trimws(as.character(MassUnit))))
  NewProblem$N["Mass"] = ThisProblem$N["Mass"] + NMassAdd
  if (any(InMass)) {
    if (length(InMass) == 1) { InMass = rep(InMass, NMassAdd) }
    NewProblem$InMassName = c(NewProblem$InMassName, MassName[InMass])
    NewProblem$N["InMass"] = NewProblem$N["InMass"] + sum(InMass)
  }

  # Update MCR indices
  NewProblem$Index$AqueousMCR = which(tolower(NewProblem$Mass$Name) %in% c("water", "aqueous"))
  NewProblem$Index$BioticLigMCR = which(
    grepl("BL", NewProblem$Mass$Name, ignore.case = TRUE) |
      grepl("gill", NewProblem$Mass$Name, ignore.case = TRUE)
  )
  if (any(NewProblem$Mass$Name %in% "Water_DonnanHA")) {
    NewProblem$Index$WHAMDonnanMCR[1] = which(NewProblem$Mass$Name == "Water_DonnanHA")
  } else {
    NewProblem$Index$WHAMDonnanMCR[1] = -1L
  }
  if (any(NewProblem$Mass$Name %in% "Water_DonnanFA")) {
    NewProblem$Index$WHAMDonnanMCR[2] = which(NewProblem$Mass$Name == "Water_DonnanFA")
  } else {
    NewProblem$Index$WHAMDonnanMCR[2] = -1L
  }

  return(NewProblem)

}


#' @rdname MassCompartments
#' @export
RemoveMassCompartments = function(ThisProblem, MCToRemove) {

  NewProblem = ThisProblem

  MCToRemoveOrig = MCToRemove
  if (is.character(MCToRemove)) {
    MCToRemove = match(MCToRemove, ThisProblem$Mass$Name)
  }
  if (any(is.na(MCToRemove))) {
    stop(paste0("Mass Compartment (", paste(MCToRemoveOrig[is.na(MCToRemove)], collapse = ", "), ") does not exist."))
  }
  if (any(MCToRemove > ThisProblem$N["Mass"])) {
    stop(paste0("There are ", ThisProblem$N["Mass"], " Mass Compartments, ",
                "trying to remove the #(",
                MCToRemove[MCToRemove > ThisProblem$N["Mass"]],
                ") element(s)."))
  }

  # Remove input variables that depend on the mass compartment
  InVarToRemove = which(NewProblem$InVar$MCR %in% MCToRemove)
  if (length(InVarToRemove) >= 1) {
    NewProblem = RemoveInVars(NewProblem, InVarToRemove)
  }

  # Remove DefComps that depend on the mass compartment
  DefCompToRemove = which(NewProblem$DefComp$MCR %in% MCToRemove)
  if (length(DefCompToRemove) >= 1) {
    NewProblem = RemoveDefComps(NewProblem, DefCompToRemove)
  }

  # Remove Components that depend on the mass compartment
  ComponentToRemove = which(NewProblem$Comp$MCR %in% MCToRemove)
  if (length(ComponentToRemove) >= 1) {
    NewProblem = RemoveComponents(NewProblem, ComponentToRemove)
  }

  # Remove Species that depend on the mass compartment
  SpeciesToRemove = which(NewProblem$Spec$MCR %in% MCToRemove)
  if (length(SpeciesToRemove) >= 1) {
    NewProblem = RemoveSpecies(NewProblem, SpeciesToRemove)
  }

  # remove mass compartment
  NewProblem$Mass = ThisProblem$Mass[-MCToRemove, , drop = FALSE]
  NewProblem$N["Mass"] = ThisProblem$N["Mass"] - length(MCToRemove)

  # Remove InMass's that are the Mass compartment to remove
  InMassToRemove = which(NewProblem$InMassName %in%
                           ThisProblem$Mass$Name[MCToRemove])
  if (length(InMassToRemove) >= 1) {
    NewProblem$InMassName = NewProblem$InMassName[-InMassToRemove]
    NewProblem$N["InMass"] = ThisProblem$N["InMass"] - length(InMassToRemove)
  }

  # Update mass compartment indices
  NewProblem$InVar$MCR = match(NewProblem$InVar$MCName, NewProblem$Mass$Name)
  NewProblem$Comp$MCR = match(NewProblem$Comp$MCName, NewProblem$Mass$Name)
  NewProblem$DefComp$MCR = match(NewProblem$DefComp$MCName, NewProblem$Mass$Name)
  NewProblem$Spec$MCR = match(NewProblem$Spec$MCName, NewProblem$Mass$Name)
  NewProblem$Index$AqueousMCR =
    which(tolower(NewProblem$Mass$Name) %in% c("water", "aqueous"))[1]
  NewProblem$Index$BioticLigMCR =
    which(grepl("BL", NewProblem$Mass$Name, ignore.case = TRUE) |
            grepl("gill", NewProblem$Mass$Name, ignore.case = TRUE))[1]
  if (any(NewProblem$Mass$Name %in% "Water_DonnanHA")) {
    NewProblem$Index$WHAMDonnanMCR[1] = which(NewProblem$Mass$Name == "Water_DonnanHA")
  } else {
    NewProblem$Index$WHAMDonnanMCR[1] = -1L
  }
  if (any(NewProblem$Mass$Name %in% "Water_DonnanFA")) {
    NewProblem$Index$WHAMDonnanMCR[2] = which(NewProblem$Mass$Name == "Water_DonnanFA")
  } else {
    NewProblem$Index$WHAMDonnanMCR[2] = -1L
  }

  return(NewProblem)

}
