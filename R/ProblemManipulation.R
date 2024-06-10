# BlankProblem = function() {
#
#   # assemble Output
#   Out = list(
#
#     # Counts
#     NMass = 0,
#     NInLab = 0,
#     NInVar = 0,
#     NInComp = 0,
#     NDefComp = 0,
#     NComp = 0,
#     NSpec = 0,
#     NPhase = 0,
#     NBL = 0,
#     NMetal = 0,
#     NBLMetal = 0,
#     NCAT = 0,
#
#     # Mass Compartment List
#     MassName = character(),
#     MassAmt = numeric(),
#     MassUnit = character(),
#     AqueousMCR = NA,
#     BioticLigMCR = NA,
#
#     # Input Labels
#     InLabName = character(),
#
#     # Input Variables
#     InVarName = character(),
#     InVarMCName = character(),
#     InVarMCR = integer(),
#     InVarType = character(),
#
#     # Input Components
#     InCompName = character(),
#     CompName = character(),
#     CompCharge = integer(),
#     CompMCName = character(),
#     CompMCR = integer(),
#     CompType = character(),
#     CompActCorr = character(),
#     CompSiteDens = numeric(),
#
#     # Defined Components
#     DefCompName = character(),
#     DefCompFromNum = numeric(),
#     DefCompFromVar = character(),
#     DefCompCharge = integer(),
#     DefCompMCName = character(),
#     DefCompMCR = integer(),
#     DefCompType = character(),
#     DefCompActCorr = character(),
#     DefCompSiteDens = numeric(),
#
#     # Formation Reactions
#     SpecName = character(),
#     SpecMCName = character(),
#     SpecMCR = integer(),
#     SpecActCorr = integer(),
#     SpecNC = integer(),
#     SpecCompList = matrix(data = 0, nrow = 0, ncol = 0),
#     SpecStoich = matrix(data = 0, nrow = 0, ncol = 0),
#     SpecLogK = numeric(),
#     SpecDeltaH = numeric(),
#     SpecTempKelvin = numeric(),
#     SpecCharge = integer(),
#     SpecK = numeric(),
#
#
#     # Phase List
#     PhaseName = character(),
#     PhaseNC = integer(),
#     PhaseCompList = matrix(data = 0, nrow = 0, ncol = 0),
#     PhaseStoich = matrix(data = 0, nrow = 0, ncol = 0),
#     PhaseLogK = numeric(),
#     PhaseDeltaH = numeric(),
#     PhaseTemp = numeric(),
#     PhaseMoles = numeric(),
#
#     # Special Definitions
#     BLName = character(),
#     BLCompR = integer(),
#     MetalName = character(),
#     MetalCompR = integer(),
#     BLMetalName = character(),
#     BLMetalSpecsR = integer(),
#     DoWHAM = FALSE,
#
#     # Critical Accumulation Table
#     CATab = data.frame(Num = integer(), CA = numeric(), Species = character(),
#                        Test.Type = character(), Duration = character(),
#                        Lifestage = character(), Endpoint = character(),
#                        Quantifier = character(), References = character(),
#                        Miscellaneous = character()),
#
#     # WHAM parameters
#     WHAMDonnanMCR = array(NA, dim = 2, dimnames = list(c("HA","FA"))),
#     wDLF = NA,
#     wKZED = NA,
#     wP = array(NA, dim = 2, dimnames = list(c("HA","FA"))),
#     wRadius = array(NA, dim = 2, dimnames = list(c("HA","FA"))),
#     wMolWt = array(NA, dim = 2, dimnames = list(c("HA","FA")))
#
#   )
#
#
#
#   return(Out)
#
# }
#
# RemoveMassCompartments = function(ThisProblem, MCToRemove) {
#
#   NewProblem = ThisProblem
#
#   MCToRemoveOrig = MCToRemove
#   if (is.character(MCToRemove)) {
#     MCToRemove = match(MCToRemove, ThisProblem$MassName)
#   }
#   if (any(is.na(MCToRemove))) {
#     stop(paste0("Mass Compartment (", paste(MCToRemoveOrig[is.na(MCToRemove)], collapse = ", "), ") does not exist."))
#   }
#   if (any(MCToRemove > ThisProblem$NMass)) {
#     stop(paste0("There are ", ThisProblem$NMass, " Mass Compartments, ",
#                 "trying to remove the #(",
#                 MCToRemove[MCToRemove > ThisProblem$NMass],
#                 ") element(s)."))
#   }
#
#   # Remove input variables that depend on the mass compartment
#   InVarToRemove = which(NewProblem$InVarMCR %in% MCToRemove)
#   if (length(InVarToRemove) >= 1) {
#     NewProblem = RemoveInVar(NewProblem, InVarToRemove)
#   }
#
#   # Remove DefComps that depend on the mass compartment
#   DefCompToRemove = which(NewProblem$DefCompMCR %in% MCToRemove)
#   if (length(DefCompToRemove) >= 1) {
#     NewProblem = RemoveDefComps(NewProblem, DefCompToRemove)
#   }
#
#   # Remove Components that depend on the mass compartment
#   ComponentToRemove = which(NewProblem$CompMCR %in% MCToRemove)
#   if (length(ComponentToRemove) >= 1) {
#     NewProblem = RemoveComponents(NewProblem, ComponentToRemove)
#   }
#
#   # Remove Species that depend on the mass compartment
#   SpeciesToRemove = which(NewProblem$SpecMCR %in% MCToRemove)
#   if (length(SpeciesToRemove) >= 1) {
#     NewProblem = RemoveSpecies(NewProblem, SpeciesToRemove)
#   }
#
#   # remove mass compartment
#   NewProblem$MassName = ThisProblem$MassName[-MCToRemove]
#   NewProblem$MassAmt = ThisProblem$MassAmt[-MCToRemove]
#   NewProblem$MassUnit = ThisProblem$MassUnit[-MCToRemove]
#   NewProblem$NMass = ThisProblem$NMass - length(MCToRemove)
#
#   # Update special definitions
#   NewProblem$InVarMCR = match(NewProblem$InVarMCName, NewProblem$MassName)
#   NewProblem$CompMCR = match(NewProblem$CompMCName, NewProblem$MassName)
#   NewProblem$DefCompMCR = match(NewProblem$DefCompMCName, NewProblem$MassName)
#   NewProblem$SpecMCR = match(NewProblem$SpecMCName, NewProblem$MassName)
#   NewProblem$AqueousMCR =
#     which(tolower(NewProblem$MassName) %in% c("water", "aqueous"))
#   NewProblem$BioticLigMCR =
#     which(grepl("BL", NewProblem$MassName, ignore.case = TRUE) |
#             grepl("gill", NewProblem$MassName, ignore.case = TRUE))
#   if (any(NewProblem$MassName %in% "Water_DonnanHA")) {
#     NewProblem$WHAMDonnanMCR[1] = which(NewProblem$MassName == "Water_DonnanHA")
#   } else {
#     NewProblem$WHAMDonnanMCR[1] = -1L
#   }
#   if (any(NewProblem$MassName %in% "Water_DonnanFA")) {
#     NewProblem$WHAMDonnanMCR[2] = which(NewProblem$MassName == "Water_DonnanFA")
#   } else {
#     NewProblem$WHAMDonnanMCR[2] = -1L
#   }
#
#   return(NewProblem)
#
# }
#
# RemoveInLab = function(ThisProblem, InLabToRemove) {
#
#   NewProblem = ThisProblem
#
#   InLabToRemoveOrig = InLabToRemove
#   if (is.character(InLabToRemove)) {
#     InLabToRemove = which(ThisProblem$InLabName %in% InLabToRemove)
#   }
#   if (length(InLabToRemove) < 1) {
#     stop(paste0("Mass Compartment \"", InLabToRemoveOrig, "\" does not exist."))
#   }
#   if (any(InLabToRemove > ThisProblem$NInLab)) {
#     stop(paste0("There are ", ThisProblem$NInLab, " input labels, ",
#                 "trying to remove the #(",
#                 InLabToRemove[InLabToRemove > ThisProblem$NInLab],
#                 ") element(s)."))
#   }
#
#   # Remove input labels
#   NewProblem$InLabName = NewProblem$InLabName[-InLabToRemove]
#   NewProblem$NInLab = ThisProblem$NInLab - length(InLabToRemove)
#
#   return(NewProblem)
#
# }
#
# RemoveInVar = function(ThisProblem, InVarToRemove) {
#
#   NewProblem = ThisProblem
#
#   if (length(InVarToRemove) < 1) {
#     stop("Specify an input variable to add.")
#   }
#
#   InVarToRemoveOrig = InVarToRemove
#   if (is.character(InVarToRemove)) {
#     InVarToRemove = match(InVarToRemove, ThisProblem$InVarName)
#   }
#   if (any(is.na(InVarToRemove))) {
#     stop(paste0("Input Variable(s) (",
#                 paste(InVarToRemoveOrig[is.na(InVarToRemove)], collapse = ", "),
#                 ") do not exist."))
#   }
#   if (any(InVarToRemove > ThisProblem$NInVar)) {
#     stop(paste0("There are ", ThisProblem$NInVar, " input variables, ",
#                 "trying to remove the #(",
#                 InVarToRemove[InVarToRemove > ThisProblem$NInVar],
#                 ") element(s)."))
#   }
#
#   # Remove DefComps that depend on input variables
#   DefCompToRemove = which(NewProblem$DefCompFromVar %in%
#                             NewProblem$InVarName[InVarToRemove])
#   if (length(DefCompToRemove) >= 1) {
#     NewProblem = RemoveDefComps(NewProblem, DefCompToRemove)
#   }
#
#   # Remove input variables
#   NewProblem$InVarName = NewProblem$InVarName[-InVarToRemove]
#   NewProblem$InVarMCName = NewProblem$InVarMCName[-InVarToRemove]
#   NewProblem$InVarMCR = NewProblem$InVarMCR[-InVarToRemove]
#   NewProblem$InVarType = NewProblem$InVarType[-InVarToRemove]
#   NewProblem$NInVar = ThisProblem$NInVar - length(InVarToRemove)
#
#   return(NewProblem)
#
# }
#
# RemoveComponents = function(ThisProblem, ComponentToRemove) {
#
#   NewProblem = ThisProblem
#
#   ComponentToRemoveOrig = ComponentToRemove
#   if (is.character(ComponentToRemove)) {
#     ComponentToRemove = match(ComponentToRemove, ThisProblem$CompName)
#   }
#   if (any(is.na(ComponentToRemove))) {
#     stop(paste0("Component(s) (",
#                 paste(ComponentToRemoveOrig[is.na(ComponentToRemove)],
#                       collapse = ", "),
#                 ") does not exist."))
#   }
#   if (any(ComponentToRemove > ThisProblem$NComp)) {
#     stop(paste0("There are ", ThisProblem$NComp, " Components, ",
#                 "trying to remove the #(",
#                 paste(ComponentToRemove[ComponentToRemove > ThisProblem$NComp],
#                       collapse = ", "),
#                 ") element(s)."))
#   }
#
#   # Remove Species that depend on the component
#   SpeciesToRemove = unique(c(
#     ComponentToRemove,
#     which(apply(NewProblem$SpecStoich[, ComponentToRemove, drop = FALSE],
#                 MARGIN = 1, FUN = function(X){any(X != 0L)}))
#   ))
#   NewProblem = RemoveSpecies(NewProblem, SpeciesToRemove)
#
#   # Remove DefComps that depend on the component
#   DefCompToRemove = which(NewProblem$DefCompFromVar %in%
#                             ThisProblem$CompName[ComponentToRemove])
#   if (length(DefCompToRemove) >= 1) {
#     NewProblem = RemoveDefComps(NewProblem, DefCompToRemove)
#   }
#
#   # Remove InComps that are the component to remove
#   InCompToRemove = which(NewProblem$InCompName %in%
#                            ThisProblem$CompName[ComponentToRemove])
#   if (length(InCompToRemove) >= 1) {
#     NewProblem$InCompName = NewProblem$InCompName[-InCompToRemove]
#     NewProblem$NInComp = ThisProblem$NInComp - length(InCompToRemove)
#   }
#
#   # Remove Components
#   NewProblem$CompName = ThisProblem$CompName[-ComponentToRemove]
#   NewProblem$CompCharge = ThisProblem$CompCharge[-ComponentToRemove]
#   NewProblem$CompMCName = ThisProblem$CompMCName[-ComponentToRemove]
#   NewProblem$CompMCR = ThisProblem$CompMCR[-ComponentToRemove]
#   NewProblem$CompType = ThisProblem$CompType[-ComponentToRemove]
#   NewProblem$CompActCorr = ThisProblem$CompActCorr[-ComponentToRemove]
#   NewProblem$CompSiteDens = ThisProblem$CompSiteDens[-ComponentToRemove]
#   NewProblem$NComp = ThisProblem$NComp - length(ComponentToRemove)
#   NewProblem$SpecStoich = NewProblem$SpecStoich[, -ComponentToRemove, drop = FALSE]
#
#   # Update Special Definitions
#   NewProblem$BLName = ThisProblem$BLName[
#     ThisProblem$BLName %in% ThisProblem$CompName[ComponentToRemove] == FALSE]
#   NewProblem$BLCompR = match(NewProblem$BLName, NewProblem$CompName)
#   NewProblem$NBL = length(NewProblem$BLName)
#
#   NewProblem$MetalName = ThisProblem$MetalName[
#     ThisProblem$MetalName %in% ThisProblem$CompName[ComponentToRemove] == FALSE]
#   NewProblem$MetalCompR = match(NewProblem$MetalName, NewProblem$CompName)
#   NewProblem$NMetal = length(NewProblem$MetalName)
#
#   # Update SpecCompList
#   NewProblem$SpecNC = rowSums(NewProblem$SpecStoich != 0L)
#   names(NewProblem$SpecNC) = NewProblem$SpecName
#   NewProblem$SpecCompList = t(apply(
#     NewProblem$SpecStoich,
#     MARGIN = 1,
#     FUN = function(X) {
#       Tmp = sort(which(X != 0L))
#       if (length(Tmp) < max(NewProblem$SpecNC)) {
#         Tmp = c(Tmp, rep(0, max(NewProblem$SpecNC) - length(Tmp)))
#       }
#       return(Tmp)
#     }
#   ))
#   rownames(NewProblem$SpecCompList) = NewProblem$SpecName
#
#   # Update PhaseCompList
#   NewProblem$PhaseNC = rowSums(NewProblem$PhaseStoich != 0L)
#   names(NewProblem$PhaseNC) = NewProblem$PhaseName
#   if (NewProblem$NPhase > 0) {
#     NewProblem$PhaseCompList = as.matrix(apply(
#       NewProblem$PhaseStoich,
#       MARGIN = 1,
#       FUN = function(X) {
#         Tmp = sort(which(X != 0L))
#         if (length(Tmp) < max(NewProblem$PhaseNC)) {
#           Tmp = c(Tmp, rep(0, max(NewProblem$PhaseNC) - length(Tmp)))
#         }
#         return(Tmp)
#       }
#     ))
#     rownames(NewProblem$PhaseCompList) = NewProblem$PhaseName
#   }
#
#   return(NewProblem)
# }
#
# RemoveSpecies = function(ThisProblem, SpeciesToRemove) {
#
#   NewProblem = ThisProblem
#
#   SpeciesToRemoveOrig = SpeciesToRemove
#   if (is.character(SpeciesToRemove)) {
#     SpeciesToRemove = which(ThisProblem$SpecName %in% SpeciesToRemove)
#   }
#   if (length(SpeciesToRemove) < 1) {
#     stop(paste0("Species \"", SpeciesToRemoveOrig, "\" does not exist."))
#   }
#   if (any(SpeciesToRemove > ThisProblem$NSpec)) {
#     stop(paste0("There are ", ThisProblem$NSpec, " Species, ",
#                 "trying to remove the #(",
#                 paste(SpeciesToRemove[SpeciesToRemove > ThisProblem$NSpec],
#                       collapse = ", "),
#                 ") element(s)."))
#   }
#
#   # Remove Species
#   NewProblem$SpecName = ThisProblem$SpecName[-SpeciesToRemove]
#   NewProblem$SpecMCName = ThisProblem$SpecMCName[-SpeciesToRemove]
#   NewProblem$SpecMCR = ThisProblem$SpecMCR[-SpeciesToRemove]
#   NewProblem$SpecActCorr = ThisProblem$SpecActCorr[-SpeciesToRemove]
#   NewProblem$SpecNC = ThisProblem$SpecNC[-SpeciesToRemove]
#   NewProblem$SpecCompList = ThisProblem$SpecCompList[-SpeciesToRemove, , drop = FALSE]
#   NewProblem$SpecLogK = ThisProblem$SpecLogK[-SpeciesToRemove]
#   NewProblem$SpecDeltaH = ThisProblem$SpecDeltaH[-SpeciesToRemove]
#   NewProblem$SpecTempKelvin = ThisProblem$SpecTempKelvin[-SpeciesToRemove]
#   NewProblem$SpecCharge = ThisProblem$SpecCharge[-SpeciesToRemove]
#   NewProblem$SpecK = ThisProblem$SpecK[-SpeciesToRemove]
#   NewProblem$SpecStoich = ThisProblem$SpecStoich[-SpeciesToRemove, , drop = FALSE]
#   NewProblem$NSpec = ThisProblem$NSpec - length(SpeciesToRemove)
#
#   # Update Special Definitions
#   NewProblem$BLMetalName = ThisProblem$BLMetalName[
#     ThisProblem$BLMetalName %in%
#       ThisProblem$SpecName[SpeciesToRemove] == FALSE]
#   NewProblem$BLMetalSpecsR = match(NewProblem$BLMetalName, NewProblem$SpecName)
#   NewProblem$NBLMetal = length(NewProblem$BLMetalName)
#
#   return(NewProblem)
#
# }
#
# RemoveDefComps = function(ThisProblem, DefCompToRemove) {
#
#   NewProblem = ThisProblem
#
#   DefCompToRemoveOrig = DefCompToRemove
#   if (is.character(DefCompToRemove)) {
#     DefCompToRemove = which(ThisProblem$DefCompName %in% DefCompToRemove)
#   }
#   if (length(DefCompToRemove) < 1) {
#     stop(paste0("Defined Component \"", DefCompToRemoveOrig, "\" does not exist."))
#   }
#   if (any(DefCompToRemove > ThisProblem$NDefComp)) {
#     stop(paste0("There are ", ThisProblem$NDefComp, " Defined Components, ",
#                 "trying to remove the #(",
#                 paste(DefCompToRemove[DefCompToRemove > ThisProblem$NDefComp],
#                       collapse = ", "),
#                 ") element(s)."))
#   }
#
#   # Remove components that are the DefComp
#   ComponentToRemove = which(ThisProblem$CompName %in%
#                               ThisProblem$DefCompName[DefCompToRemove])
#   NewProblem = RemoveComponents(NewProblem, ComponentToRemove)
#
#   # Remove DefComps
#   NewProblem$DefCompName = ThisProblem$DefCompName[-DefCompToRemove]
#   NewProblem$DefCompFromNum = ThisProblem$DefCompFromNum[-DefCompToRemove]
#   NewProblem$DefCompFromVar = ThisProblem$DefCompFromVar[-DefCompToRemove]
#   NewProblem$DefCompCharge = ThisProblem$DefCompCharge[-DefCompToRemove]
#   NewProblem$DefCompMCName = ThisProblem$DefCompMCName[-DefCompToRemove]
#   NewProblem$DefCompMCR = ThisProblem$DefCompMCR[-DefCompToRemove]
#   NewProblem$DefCompType = ThisProblem$DefCompType[-DefCompToRemove]
#   NewProblem$DefCompActCorr = ThisProblem$DefCompActCorr[-DefCompToRemove]
#   NewProblem$DefCompSiteDens = ThisProblem$DefCompSiteDens[-DefCompToRemove]
#   NewProblem$NDefComp = ThisProblem$NDefComp - length(DefCompToRemove)
#
#   return(NewProblem)
#
# }
#
# AddMassCompartments = function(ThisProblem, MassName, MassAmt, MassUnit) {
#
#   NewProblem = ThisProblem
#
#   # error checking
#   if (any((MassName %in% ThisProblem$MassName))) {
#     stop(paste0(
#       "Mass Compartment(s) (",
#       paste(paste0("\"", MassName[MassName %in% ThisProblem$MassName], "\""), collapse = ", "),
#       ") already exist."
#     ))
#   }
#   NMassAdd = length(MassName)
#   if ((NMassAdd != length(MassAmt)) || (NMassAdd != length(MassUnit))) {
#     stop("Inputs differ in length.  Be explicit for each mass compartment you're adding.")
#   }
#   if (any(is.na(c(MassName, MassAmt, MassUnit)))) {
#     stop("NA arguments not allowed.")
#   }
#
#   # add mass compartment
#   NewProblem$MassName = c(ThisProblem$MassName, MassName)
#   NewProblem$MassAmt = c(ThisProblem$MassAmt,
#                          array(MassAmt, dimnames = list(MassName)))
#   NewProblem$MassUnit = c(ThisProblem$MassUnit,
#                           array(MassUnit, dimnames = list(MassName)))
#   NewProblem$NMass = ThisProblem$NMass + NMassAdd
#
#   return(NewProblem)
#
# }
#
# AddInLabs = function(ThisProblem, InLabName) {
#
#   NewProblem = ThisProblem
#
#   # error checking
#   if (any((InLabName %in% ThisProblem$InLabName))) {
#     stop(paste0(
#       "Input Label(s) (",
#       paste(paste0("\"", InLabName[InLabName %in% ThisProblem$InLabName], "\""), collapse = ", "),
#       ") already exist."
#     ))
#   }
#   if (any(is.na(InLabName))) {
#     stop("NA arguments not allowed.")
#   }
#
#   # add input labels
#   NewProblem$InLabName = c(ThisProblem$InLabName, InLabName)
#   NewProblem$NInLab = ThisProblem$NInLab + length(InLabName)
#
#   return(NewProblem)
#
# }
#
# AddInVar = function(ThisProblem, InVarName, InVarMCName = NULL,
#                     InVarType,
#                     InVarMCR = match(InVarMCName, ThisProblem$MassName)) {
#
#   NewProblem = ThisProblem
#
#   # error checking
#   if (any((InVarName %in% ThisProblem$InVarName))) {
#     stop(paste0(
#       "Input Variable(s) (",
#       paste(paste0("\"", InVarName[InVarName %in% ThisProblem$InVarName], "\""), collapse = ", "),
#       ") already exist."
#     ))
#   }
#   if (!all(InVarMCName %in% ThisProblem$MassName) ||
#       !all(InVarMCR <= ThisProblem$NMass)) {
#     stop("Mass compartment(s) specified in InVarMCName or InVarMCR does not exist.")
#   }
#   if (any(is.na(c(InVarName, InVarType, InVarMCName, InVarMCR)))) {
#     stop("NA arguments not allowed.")
#   }
#   if(is.null(InVarMCName)) {
#     InVarMCName = ThisProblem$MassName[InVarMCR]
#   }
#   if (!all(InVarMCName %in% ThisProblem$MassName[InVarMCR])) {
#     stop("InVarMCName does not match InVarMCR.")
#   }
#   NInVarAdd = length(InVarName)
#   if ((NInVarAdd != length(InVarType)) ||
#       (NInVarAdd != length(InVarMCName)) ||
#       (NInVarAdd != length(InVarMCR))) {
#     stop("Inputs differ in length.  Be explicit for each input variable you're adding.")
#   }
#
#   # Add input variable
#   NewProblem$InVarName = c(ThisProblem$InVarName, InVarName)
#   NewProblem$InVarMCName = c(ThisProblem$InVarMCName,
#                              array(InVarMCName, dimnames = list(InVarName)))
#   NewProblem$InVarMCR = c(ThisProblem$InVarMCR,
#                           array(InVarMCR, dimnames = list(InVarName)))
#   NewProblem$InVarType = c(ThisProblem$InVarType,
#                            array(InVarType, dimnames = list(InVarName)))
#   NewProblem$NInVar = ThisProblem$NInVar + NInVarAdd
#
#   for (i in 1:NInVarAdd) {
#     if (InVarType[i] == "pH") {
#       NewProblem = AddDefComp(NewProblem,
#                               DefCompName = "H",
#                               DefCompFromVar = InVarName[i],
#                               DefCompCharge = 1L,
#                               DefCompMCName = InVarMCName[i],
#                               DefCompType = "FixedAct",
#                               DefCompActCorr = "Debye",
#                               DefCompSiteDens = 1.0)
#       warning(paste("Defined component added for pH as a fixed activity",
#                     "component. Change manually if you wish to represent",
#                     "pH as a fixed concentration."))
#     }
#   }
#
#   return(NewProblem)
#
# }
#
# AddComponents = function(ThisProblem, CompName,  CompCharge, CompMCName = NULL,
#                          CompType, CompActCorr, CompSiteDens,
#                          CompMCR = match(CompMCName, ThisProblem$MassName)) {
#
#   NewProblem = ThisProblem
#
#   # error checking
#   if (any(CompName %in% ThisProblem$SpecName)) {
#     stop(paste0(
#       "Input Component(s) (",
#       paste(paste0("\"", CompName[CompName %in% ThisProblem$SpecName], "\""),
#             collapse = ", "),
#       ") already exist as a component or species."
#     ))
#   }
#   if (!all(CompMCName %in% ThisProblem$MassName) ||
#       !all(CompMCR <= ThisProblem$NMass)) {
#     stop("Mass compartment(s) specified in CompMCName or CompMCR does not exist.")
#   }
#   if (any(is.na(c(CompName, CompCharge, CompMCName,
#                   CompType, CompActCorr, CompMCR, CompSiteDens)))) {
#     stop("NA arguments not allowed.")
#   }
#   if(is.null(CompMCName)) {
#     CompMCName = ThisProblem$MassName[CompMCR]
#   }
#   if (!all(CompMCName %in% ThisProblem$MassName[CompMCR])) {
#     stop("CompMCName does not match CompMCR.")
#   }
#   NCompAdd = length(CompName)
#   if ((NCompAdd != length(CompCharge)) ||
#       (NCompAdd != length(CompMCName)) ||
#       (NCompAdd != length(CompType)) ||
#       (NCompAdd != length(CompActCorr)) ||
#       (NCompAdd != length(CompMCR)) ||
#       (NCompAdd != length(CompSiteDens))) {
#     stop("Inputs differ in length.  Be explicit for each input component you're adding.")
#   }
#
#   # add Species for the component - before the reactions
#   if (ThisProblem$NComp > 0) {
#     CompsVec = 1:ThisProblem$NComp
#   } else {
#     CompsVec = integer()
#   }
#   if (ThisProblem$NSpec > ThisProblem$NComp) {
#     SpecRxnVec = (ThisProblem$NComp + 1):ThisProblem$NSpec
#   } else {
#     SpecRxnVec = integer()
#   }
#   NewProblem$SpecName = c(ThisProblem$CompName, CompName,
#                           ThisProblem$SpecName[SpecRxnVec])
#   # SpecNameRxns = setdiff(ThisProblem$SpecName, ThisProblem$CompName)
#
#   NewProblem$SpecMCName = c(ThisProblem$CompMCName,
#                             array(CompMCName, dimnames = list(CompName)),
#                             ThisProblem$SpecMCName[SpecRxnVec])
#   NewProblem$SpecMCR = c(ThisProblem$CompMCR,
#                          array(CompMCR, dimnames = list(CompName)),
#                          ThisProblem$SpecMCR[SpecRxnVec])
#   NewProblem$SpecActCorr = c(ThisProblem$CompActCorr,
#                              array(CompActCorr, dimnames = list(CompName)),
#                              ThisProblem$SpecActCorr[SpecRxnVec])
#   NewProblem$SpecLogK = c(ThisProblem$SpecLogK[CompsVec],
#                           array(0.0, dim = NCompAdd, dimnames = list(CompName)),
#                           ThisProblem$SpecLogK[SpecRxnVec])
#   NewProblem$SpecK = 10^(NewProblem$SpecLogK)
#   NewProblem$SpecDeltaH = c(ThisProblem$SpecDeltaH[CompsVec],
#                             array(0.0, dim = NCompAdd, dimnames = list(CompName)),
#                             ThisProblem$SpecDeltaH[SpecRxnVec])
#   NewProblem$SpecTempKelvin =  c(ThisProblem$SpecTempKelvin[CompsVec],
#                                  array(0.0, dim = NCompAdd, dimnames = list(CompName)),
#                                  ThisProblem$SpecTempKelvin[SpecRxnVec])
#   NewProblem$SpecNC = c(ThisProblem$SpecNC[CompsVec],
#                         array(1L, dim = NCompAdd, dimnames = list(CompName)),
#                         ThisProblem$SpecNC[SpecRxnVec])
#   if (ThisProblem$NSpec == 0) {
#     NewProblem$SpecCompList = matrix(0L,
#                                      nrow = NCompAdd,
#                                      ncol = 1L,
#                                      dimnames = list(CompName, NULL))
#   } else {
#     NewProblem$SpecCompList = rbind(
#       ThisProblem$SpecCompList[CompsVec, , drop = FALSE],
#       matrix(0L,
#              nrow = NCompAdd,
#              ncol = ncol(ThisProblem$SpecCompList),
#              dimnames = list(CompName, NULL)),
#       ThisProblem$SpecCompList[SpecRxnVec, , drop = FALSE]
#     )
#   }
#   NewProblem$SpecCompList[ThisProblem$NComp + (1:NCompAdd), 1] =
#     as.integer(ThisProblem$NComp + (1:NCompAdd))
#
#   NewProblem$SpecStoich = rbind(
#     ThisProblem$SpecStoich[CompsVec, , drop = FALSE],
#     matrix(0L, nrow = NCompAdd, ncol = ThisProblem$NComp,
#            dimnames = list(CompName, NULL)),
#     ThisProblem$SpecStoich[SpecRxnVec, , drop = F]
#   )
#   NewProblem$SpecStoich = cbind(
#     NewProblem$SpecStoich,
#     matrix(data = 0L,
#            nrow = ThisProblem$NSpec + NCompAdd,
#            ncol = NCompAdd,
#            dimnames = list(NULL, CompName))
#   )
#
#   for (i in ((1:NCompAdd) + ThisProblem$NComp)) {
#     NewProblem$SpecStoich[i, NewProblem$SpecCompList[i, 1]] = 1L
#   }
#   NewProblem$NSpec = ThisProblem$NSpec + NCompAdd
#
#   # add Components
#   NewProblem$CompName = c(ThisProblem$CompName, CompName)
#   NewProblem$CompCharge = c(ThisProblem$CompCharge,
#                             array(CompCharge, dimnames = list(CompName)))
#   NewProblem$CompMCName = c(ThisProblem$CompMCName,
#                             array(CompMCName, dimnames = list(CompName)))
#   NewProblem$CompMCR = c(ThisProblem$CompMCR,
#                          array(CompMCR, dimnames = list(CompName)))
#   NewProblem$CompType = c(ThisProblem$CompType,
#                           array(CompType, dimnames = list(CompName)))
#   NewProblem$CompActCorr = c(ThisProblem$CompActCorr,
#                              array(CompActCorr, dimnames = list(CompName)))
#   NewProblem$CompSiteDens = c(ThisProblem$CompSiteDens,
#                               array(1.0, dim = NCompAdd, dimnames = list(CompName)))
#   NewProblem$NComp = ThisProblem$NComp + NCompAdd
#
#   # Update SpecCompList
#   NewProblem$SpecNC = rowSums(NewProblem$SpecStoich != 0L)
#   names(NewProblem$SpecNC) = NewProblem$SpecName
#   NewProblem$SpecCompList = as.matrix(apply(
#     NewProblem$SpecStoich,
#     MARGIN = 1,
#     FUN = function(X) {
#       Tmp = sort(which(X != 0L))
#       if (length(Tmp) < max(NewProblem$SpecNC)) {
#         Tmp = c(Tmp, rep(0, max(NewProblem$SpecNC) - length(Tmp)))
#       }
#       return(Tmp)
#     }
#   ))
#   rownames(NewProblem$SpecCompList) = NewProblem$SpecName
#   NewProblem$SpecCharge = drop(NewProblem$SpecStoich %*% NewProblem$CompCharge)
#
#   # Update Special Definitions
#   NewProblem$BLMetalSpecsR = match(NewProblem$BLMetalName, NewProblem$SpecName)
#
#   return(NewProblem)
#
# }
#
# AddInComp = function(ThisProblem, InCompName, InCompCharge, InCompMCName = NULL,
#                      InCompType, InCompActCorr,
#                      InCompMCR = match(InCompMCName, ThisProblem$MassName)) {
#
#   NewProblem = ThisProblem
#
#   # error checking
#   if (any(InCompName %in% ThisProblem$SpecName)) {
#     stop(paste0(
#       "Input Component(s) (",
#       paste(paste0("\"", InCompName[InCompName %in% ThisProblem$SpecName], "\""),
#             collapse = ", "),
#       ") already exist as a component or species."
#     ))
#   }
#   if (!all(InCompMCName %in% ThisProblem$MassName) ||
#       !all(InCompMCR <= ThisProblem$NMass)) {
#     stop("Mass compartment(s) specified in InCompMCName or InCompMCR does not exist.")
#   }
#   if (any(is.na(c(InCompName, InCompCharge, InCompMCName,
#                   InCompType, InCompActCorr, InCompMCR)))) {
#     stop("NA arguments not allowed.")
#   }
#   if(is.null(InCompMCName)) {
#     InCompMCName = ThisProblem$MassName[InCompMCR]
#   }
#   if (!all(InCompMCName %in% ThisProblem$MassName[InCompMCR])) {
#     stop("InCompMCName does not match InCompMCR.")
#   }
#   NInCompAdd = length(InCompName)
#   if ((NInCompAdd != length(InCompCharge)) ||
#       (NInCompAdd != length(InCompMCName)) ||
#       (NInCompAdd != length(InCompType)) ||
#       (NInCompAdd != length(InCompActCorr)) ||
#       (NInCompAdd != length(InCompMCR))) {
#     stop("Inputs differ in length.  Be explicit for each input component you're adding.")
#   }
#
#   # Add the input component
#   NewProblem$NInComp = NewProblem$NInComp + NInCompAdd
#   NewProblem$InCompName = c(NewProblem$InCompName, InCompName)
#
#   # Add the component
#   NewProblem = AddComponents(ThisProblem = NewProblem,
#                              CompName = InCompName,
#                              CompCharge = InCompCharge,
#                              CompMCName = InCompMCName,
#                              CompType = InCompType,
#                              CompActCorr = InCompActCorr,
#                              CompSiteDens = rep(1.0, NInCompAdd))
#
#   return(NewProblem)
#
# }
#
# AddDefComp = function(ThisProblem, DefCompName, DefCompFromNum = NULL,
#                       DefCompFromVar = NULL, DefCompCharge, DefCompMCName = NULL,
#                       DefCompType, DefCompActCorr, DefCompSiteDens,
#                       DefCompMCR = match(DefCompMCName, ThisProblem$MassName)) {
#
#   NewProblem = ThisProblem
#
#   # error checking
#   if (any(DefCompName %in% ThisProblem$SpecName)) {
#     stop(paste0(
#       "Input Component(s) (",
#       paste(paste0("\"", DefCompName[DefCompName %in% ThisProblem$SpecName], "\""),
#             collapse = ", "),
#       ") already exist as a component or species."
#     ))
#   }
#   if (!all(DefCompMCName %in% ThisProblem$MassName) ||
#       !all(DefCompMCR <= ThisProblem$NMass)) {
#     stop("Mass compartment(s) specified in DefCompMCName or DefCompMCR does not exist.")
#   }
#   if (any(is.na(c(DefCompName, DefCompCharge, DefCompMCName,
#                   DefCompType, DefCompActCorr, DefCompMCR, DefCompSiteDens)))) {
#     stop("NA arguments not allowed, except for DefCompFromNum and DefCompFromVar.")
#   }
#   if (is.null(DefCompFromVar)){
#     DefCompFromVar = array(NA, dimnames = list(DefCompName))
#   }
#   if (is.null(DefCompFromNum)){
#     DefCompFromNum = array(NA, dimnames = list(DefCompName))
#   }
#   if(is.null(DefCompMCName)) {
#     DefCompMCName = ThisProblem$MassName[DefCompMCR]
#   }
#   if (!all(DefCompMCName %in% ThisProblem$MassName[DefCompMCR])) {
#     stop("DefCompMCName does not match DefCompMCR.")
#   }
#   NDefCompAdd = length(DefCompName)
#   if ((NDefCompAdd != length(DefCompCharge)) ||
#       (NDefCompAdd != length(DefCompMCName)) ||
#       (NDefCompAdd != length(DefCompType)) ||
#       (NDefCompAdd != length(DefCompActCorr)) ||
#       (NDefCompAdd != length(DefCompMCR)) ||
#       (NDefCompAdd != length(DefCompSiteDens))) {
#     stop("Inputs differ in length.  Be explicit for each input component you're adding.")
#   }
#
#   # Add the defined component
#   NewProblem$NDefComp = ThisProblem$NDefComp + NDefCompAdd
#   NewProblem$DefCompName = c(ThisProblem$DefCompName, DefCompName)
#   NewProblem$DefCompFromNum = c(ThisProblem$DefCompFromNum,
#                                 array(DefCompFromNum, dimnames = list(DefCompName)))
#   NewProblem$DefCompFromVar = c(ThisProblem$DefCompFromVar,
#                                 array(DefCompFromVar, dimnames = list(DefCompName)))
#   NewProblem$DefCompCharge = c(ThisProblem$DefCompCharge,
#                                array(DefCompCharge, dimnames = list(DefCompName)))
#   NewProblem$DefCompMCName = c(ThisProblem$DefCompMCName,
#                                array(DefCompMCName, dimnames = list(DefCompName)))
#   NewProblem$DefCompMCR = c(ThisProblem$DefCompMCR,
#                             array(DefCompMCR, dimnames = list(DefCompName)))
#   NewProblem$DefCompType = c(ThisProblem$DefCompType,
#                              array(DefCompType, dimnames = list(DefCompName)))
#   NewProblem$DefCompActCorr = c(ThisProblem$DefCompActCorr,
#                                 array(DefCompActCorr, dimnames = list(DefCompName)))
#   NewProblem$DefCompSiteDens = c(ThisProblem$DefCompSiteDens,
#                                  array(DefCompSiteDens, dimnames = list(DefCompName)))
#
#   # Add the component
#   NewProblem = AddComponents(ThisProblem = NewProblem,
#                              CompName = DefCompName,
#                              CompCharge = DefCompCharge,
#                              CompMCName = DefCompMCName,
#                              CompType = DefCompType,
#                              CompActCorr = DefCompActCorr,
#                              CompSiteDens = DefCompSiteDens)
#
#   return(NewProblem)
#
# }
#
# AddSpecies = function(ThisProblem, SpecName, SpecMCName = NULL, SpecActCorr,
#                       SpecCompNames = list(), SpecCompStoichs = list(),
#                       SpecLogK, SpecDeltaH, SpecTempKelvin,
#                       SpecMCR = match(SpecMCName, ThisProblem$MassName)) {
#
#   NewProblem = ThisProblem
#
#   # error checking
#   if (any(SpecName %in% ThisProblem$SpecName)) {
#     stop(paste0(
#       "Species (",
#       paste(paste0("\"", SpecName[SpecName %in% ThisProblem$SpecName], "\""),
#             collapse = ", "),
#       ") already exist as a component or species."
#     ))
#   }
#   if (!all(SpecMCName %in% ThisProblem$MassName) ||
#       !all(SpecMCR <= ThisProblem$NMass)) {
#     stop("Mass compartment(s) specified in SpecMCName or SpecMCR does not exist.")
#   }
#   if (any(is.na(c(SpecName, SpecMCName,
#                   SpecLogK, SpecDeltaH, SpecTempKelvin,
#                   SpecActCorr, SpecMCR)))) {
#     stop("NA arguments not allowed.")
#   }
#   if(is.null(SpecMCName)) {
#     SpecMCName = ThisProblem$MassName[SpecMCR]
#   }
#   if (!all(SpecMCName %in% ThisProblem$MassName[SpecMCR])) {
#     stop("SpecMCName does not match SpecMCR.")
#   }
#   NSpecAdd = length(SpecName)
#   if ((NSpecAdd != length(SpecMCName)) ||
#       (NSpecAdd != length(SpecActCorr)) ||
#       (NSpecAdd != length(SpecMCR)) ||
#       (NSpecAdd != length(SpecCompNames)) ||
#       (NSpecAdd != length(SpecCompStoichs))) {
#     stop("Inputs differ in length.  Be explicit for each input component you're adding.")
#   }
#   SpecNC = sapply(SpecCompNames, FUN = length)
#   if (any(SpecNC != sapply(SpecCompStoichs, FUN = length))) {
#     stop("Stoichiometric inputs differ in length.")
#   }
#
#   # add Species
#   NewProblem$SpecName = c(ThisProblem$SpecName, SpecName)
#   NewProblem$SpecMCName = c(ThisProblem$SpecMCName,
#                             array(SpecMCName, dimnames = list(SpecName)))
#   NewProblem$SpecMCR = c(ThisProblem$SpecMCR,
#                          array(SpecMCR, dimnames = list(SpecName)))
#   NewProblem$SpecActCorr = c(ThisProblem$SpecActCorr,
#                              array(SpecActCorr, dimnames = list(SpecName)))
#   NewProblem$SpecLogK = c(ThisProblem$SpecLogK,
#                           array(SpecLogK, dimnames = list(SpecName)))
#   NewProblem$SpecDeltaH = c(ThisProblem$SpecDeltaH,
#                             array(SpecDeltaH, dimnames = list(SpecName)))
#   NewProblem$SpecTempKelvin = c(ThisProblem$SpecTempKelvin,
#                                 array(SpecTempKelvin, dimnames = list(SpecName)))
#   NewProblem$SpecNC = c(ThisProblem$SpecNC,
#                         array(SpecNC, dimnames = list(SpecName)))
#   NewProblem$SpecCompList = rbind(
#     ThisProblem$SpecCompList,
#     matrix(0L, nrow = NSpecAdd, ncol = ncol(ThisProblem$SpecCompList),
#            dimnames = list(SpecName, NULL))
#   )
#   if (max(NewProblem$SpecNC) > ncol(ThisProblem$SpecCompList)) {
#     while (max(NewProblem$SpecNC) > ncol(NewProblem$SpecCompList)) {
#       NewProblem$SpecCompList = cbind(
#         NewProblem$SpecCompList,
#         matrix(0L, nrow = nrow(NewProblem$SpecCompList), ncol = 1)
#       )
#     }
#   }
#   NewProblem$SpecStoich = rbind(
#     ThisProblem$SpecStoich,
#     matrix(0L, nrow = NSpecAdd, ncol = ThisProblem$NComp,
#            dimnames = list(SpecName, ThisProblem$CompName))
#   )
#   for (i in 1:NSpecAdd) {
#     for (j in 1:SpecNC[i]) {
#       NewProblem$SpecCompList[i + ThisProblem$NSpec, j] =
#         match(SpecCompNames[[i]][j], ThisProblem$CompName)
#       NewProblem$SpecStoich[i + ThisProblem$NSpec,
#                             NewProblem$SpecCompList[i + ThisProblem$NSpec, j]] =
#         as.integer(SpecCompStoichs[[i]][j])
#     }
#   }
#   NewProblem$NSpec = ThisProblem$NSpec + NSpecAdd
#
#   # These parts it's easier to just calculate again
#   NewProblem$SpecK = 10^(NewProblem$SpecLogK)
#   NewProblem$SpecCharge = drop(NewProblem$SpecStoich %*% NewProblem$CompCharge)
#
#   return(NewProblem)
#
# }
#
#
