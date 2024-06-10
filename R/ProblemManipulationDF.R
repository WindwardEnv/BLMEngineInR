BlankProblemDF = function() {

  # assemble Output
  Out = list(

    # Counts
    N = c(
      Mass = 0L,
      InLab = 0L,
      InVar = 0L,
      InComp = 0L,
      DefComp = 0L,
      Comp = 0L,
      Spec = 0L,
      Phase = 0L,
      BL = 0L,
      Metal = 0L,
      BLMetal = 0L,
      CAT = 0L
    ),

    # Mass Compartment List
    Mass = data.frame(
      Name = character(),
      Amt = numeric(),
      Unit = character()
    ),


    # Input Labels
    InLabName = character(),

    # Input Variables
    InVar = data.frame(
      Name = character(),
      MCName = character(),
      MCR = integer(),
      Type = character()
    ),

    # Input Components
    InCompName = character(),

    # Defined Components
    DefComp = data.frame(
      Name = character(),
      FromNum = numeric(),
      FromVar = character(),
      Charge = integer(),
      MCName = character(),
      MCR = integer(),
      Type = character(),
      ActCorr = character(),
      SiteDens = numeric()
    ),

    # Components
    Comp = data.frame(
      Name = character(),
      Charge = integer(),
      MCName = character(),
      MCR = integer(),
      Type = character(),
      ActCorr = character(),
      SiteDens = numeric()
    ),

    # Formation Reactions
    Spec = data.frame(
      Name = character(),
      Equation = character(),
      Charge = integer(),
      MCName = character(),
      MCR = integer(),
      ActCorr = integer(),
      LogK = numeric(),
      K = numeric(),
      DeltaH = numeric(),
      TempKelvin = numeric(),
      NC = integer()
    ),
    SpecCompList = matrix(data = 0, nrow = 0, ncol = 0),
    SpecStoich = matrix(data = 0, nrow = 0, ncol = 0),

    # Phase List
    Phase = data.frame(
      Name = character(),
      Equation = character(),
      NC = integer(),
      LogK = numeric(),
      DeltaH = numeric(),
      Temp = numeric(),
      Moles = numeric()
    ),
    PhaseCompList = matrix(data = 0, nrow = 0, ncol = 0),
    PhaseStoich = matrix(data = 0, nrow = 0, ncol = 0),


    # Special Definitions
    BL = data.frame(
      Name = character(),
      CompR = integer()
    ),
    Metal = data.frame(
      Name = character(),
      CompR = integer()
    ),
    BLMetal = data.frame(
      Name = character(),
      SpecsR = integer()
    ),

    # Critical Accumulation Table
    CATab = data.frame(
      Num = integer(),
      CA = numeric(),
      Species = character(),
      Test.Type = character(),
      Duration = character(),
      Lifestage = character(),
      Endpoint = character(),
      Quantifier = character(),
      References = character(),
      Miscellaneous = character()
    ),

    # WHAM parameters
    WHAM = list(
      DoWHAM = FALSE,
      WHAMVer = NA,
      WdatFile = NA,
      wDLF = as.numeric(NA),
      wKZED = as.numeric(NA),
      wP = as.numeric(array(NA, dim = 2, dimnames = list(c("HA","FA")))),
      wRadius = as.numeric(array(NA, dim = 2, dimnames = list(c("HA","FA")))),
      wMolWt = as.numeric(array(NA, dim = 2, dimnames = list(c("HA","FA"))))
    ),

    Index = list(
      AqueousMCR = as.integer(NA),
      BioticLigMCR = as.integer(NA),
      WHAMDonnanMCR = as.integer(array(NA, dim = 2, dimnames = list(c("HA","FA"))))
    )

  )



  return(Out)

}

ConvertToList = function(ThisProblemDF) {

  ThisProblemList = list()

  CompositeNames = c("N", "Mass", "InLab", "InVar", "InComp", "DefComp", "Comp",
                     "Spec", "Phase", "BL", "Metal", "BLMetal")
  OrganizedNames = c("Index", "WHAM")
  AsIsNames = setdiff(names(ThisProblemDF), c(CompositeNames, OrganizedNames))

  for (i in CompositeNames) {
    if (is.data.frame(ThisProblemDF[[i]])) {
      for (j in names(ThisProblemDF[[i]])) {
        ThisProblemList[[paste0(i, j)]] = ThisProblemDF[[i]][, j]
        if (any(names(ThisProblemDF[[i]]) %in% "Name")) {
          names(ThisProblemList[[paste0(i, j)]]) = ThisProblemDF[[i]]$Name
        }
      }
    } else if (is.vector(ThisProblemDF[[i]]) || is.list(ThisProblemDF[[i]])) {
      for (j in names(ThisProblemDF[[i]])) {
        ThisProblemList[[paste0(i, j)]] = ThisProblemDF[[i]][j]
      }
    }
  }

  for (i in OrganizedNames) {
    ThisProblemList = c(ThisProblemList, ThisProblemDF[[i]])
  }

  for (i in AsIsNames) {
    ThisProblemList[[i]] = ThisProblemDF[[i]]
  }

  return(ThisProblemList)

}

ConvertToDF = function(ThisProblemList) {

  ThisProblemDF = list()

  ListNames = names(ThisProblemList)
  CompositeNames = c("N", "Mass", "InLab", "InVar", "InComp", "DefComp", "Comp",
                     "Spec", "Phase", "BL", "Metal", "BLMetal")
  OrganizedNames = list(
    Index = c("AqueousMCR", "BioticLigMCR", "WHAMDonnanMCR"),
    WHAM = c("DoWHAM", "WHAMVer", "WdatFile", "wDLF", "wKZED",
             "wP", "wRadius", "wMolWt")
  )

  ConvertedNames = c()

  MatchedByMultipleCompositeNames = ListNames[rowSums(sapply(
    CompositeNames,
    FUN = function(X){
      grepl(paste0("^", X), ListNames)
    })) > 1]

  for (i in CompositeNames) {

    MatchedByI = ListNames[grepl(paste0("^", i), ListNames)]
    MatchedByINotMultiple = setdiff(MatchedByI, MatchedByMultipleCompositeNames)
    if (length(MatchedByINotMultiple) > 0) {
      MemberNames = MatchedByINotMultiple
    } else {
      MemberNames = MatchedByI
    }
    VectorNames = MemberNames[sapply(ThisProblemList[MemberNames], FUN = is.vector)]
    OtherNames = setdiff(MemberNames, VectorNames)

    if (length(VectorNames) == 1) {
      ThisProblemDF[[VectorNames]] = ThisProblemList[[VectorNames]]
    } else {
      if (all(sapply(ThisProblemList[VectorNames], FUN = length) == 1)) {
        ThisProblemDF[[i]] = unlist(ThisProblemList[VectorNames])
      } else {
        ThisProblemDF[[i]] = as.data.frame(ThisProblemList[VectorNames])
      }
      names(ThisProblemDF[[i]]) = gsub(paste0("^", i), "", VectorNames)
    }

    for (j in OtherNames) {
      ThisProblemDF[[j]] = ThisProblemList[[j]]
    }

    ConvertedNames = c(ConvertedNames, MemberNames)

  }

  for (i in names(OrganizedNames)) {
    ThisProblemDF[[i]] = list()
    ThisProblemDF[[i]][OrganizedNames[[i]]] = ThisProblemList[OrganizedNames[[i]]]
    ConvertedNames = c(ConvertedNames, OrganizedNames[[i]])
  }


  AsIsNames = setdiff(ListNames, ConvertedNames)
  for (i in AsIsNames) {
    ThisProblemDF[[i]] = ThisProblemList[[i]]
  }

  ThisProblemDF$Spec$Equation = StoichToEquation(
    ThisProblemList$SpecStoich,
    ThisProblemList$SpecName,
    ThisProblemList$CompName)

  return(ThisProblemDF)

}

StoichToEquation = function(SpecStoich, SpecName, CompName) {
  paste(SpecName,
              apply(SpecStoich, MARGIN = 1, FUN = function(X){
                X.nonzero = X[X != 0]
                X.react = names(X.nonzero)
                return(gsub(" [+] -", " - ",
                            paste(paste(X.nonzero, X.react, sep = " * "),
                                  collapse = " + ")))
              }),
              sep = " = "
  )
}

AddMassCompartments = function(ThisProblem, MassName, MassAmt, MassUnit) {

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

AddInLabs = function(ThisProblem, InLabName) {

  NewProblem = ThisProblem

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

  return(NewProblem)

}

AddInVar = function(ThisProblem, InVarName, InVarMCName = NULL,
                    InVarType,
                    InVarMCR = match(InVarMCName, ThisProblem$Mass$Name)) {

  NewProblem = ThisProblem

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
      NewProblem = AddDefComp(ThisProblem = NewProblem,
                              DefCompName = "H",
                              DefCompFromVar = InVarName[i],
                              DefCompFromNum = NA,
                              DefCompCharge = 1L,
                              DefCompMCName = InVarMCName[i],
                              DefCompType = "FixedAct",
                              DefCompActCorr = "Debye",
                              DefCompSiteDens = 1.0,
                              DefCompMCR = InVarMCR[i])
      # warning(paste("Defined component added for pH as a fixed activity",
      #               "component. Change manually if you wish to represent",
      #               "pH as a fixed concentration."))
    }
  }

  return(NewProblem)

}

AddComponents = function(ThisProblem, CompName,  CompCharge, CompMCName = NULL,
                         CompType, CompActCorr, CompSiteDens,
                         CompMCR = match(CompMCName, ThisProblem$Mass$Name)) {

  NewProblem = ThisProblem

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
  if(is.null(CompMCName)) {
    CompMCName = ThisProblem$Mass$Name[CompMCR]
  }
  if (!all(CompMCName %in% ThisProblem$Mass$Name[CompMCR])) {
    stop("CompMCName does not match CompMCR.")
  }
  NCompAdd = length(CompName)

  # add Components
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
  NewProblem = AddSpecies(ThisProblem = NewProblem,
                          SpecName = CompName,
                          SpecMCName = CompMCName,
                          SpecActCorr = CompActCorr,
                          SpecCompNames = as.list(CompName),
                          SpecCompStoichs = as.list(rep(1, NCompAdd)),
                          SpecLogK = 0,
                          SpecDeltaH = 0,
                          SpecTempKelvin = 0,
                          SpecMCR = CompMCR)

  # Update Special Definitions
  NewProblem$BLMetal$SpecsR = match(NewProblem$BLMetal$Name, NewProblem$Spec$Name)

  return(NewProblem)

}

AddInComp = function(ThisProblem, InCompName, InCompCharge, InCompMCName = NULL,
                     InCompType, InCompActCorr,
                     InCompMCR = match(InCompMCName, ThisProblem$Mass$Name)) {

  NewProblem = ThisProblem

  # error checking
  if (any(InCompName %in% ThisProblem$Spec$Name)) {
    stop(paste0(
      "Input Component(s) (",
      paste(paste0("\"", InCompName[InCompName %in% ThisProblem$Spec$Name], "\""),
            collapse = ", "),
      ") already exist as a component or species."
    ))
  }
  if (!all(InCompMCName %in% ThisProblem$Mass$Name) ||
      !all(InCompMCR <= ThisProblem$N["Mass"])) {
    stop("Mass compartment(s) specified in InCompMCName or InCompMCR does not exist.")
  }
  if (any(is.na(c(InCompName, InCompCharge, InCompMCName,
                  InCompType, InCompActCorr, InCompMCR)))) {
    stop("NA arguments not allowed.")
  }
  if(is.null(InCompMCName)) {
    InCompMCName = ThisProblem$Mass$Name[InCompMCR]
  }
  if (!all(InCompMCName %in% ThisProblem$Mass$Name[InCompMCR])) {
    stop("InCompMCName does not match InCompMCR.")
  }
  NInCompAdd = length(InCompName)

  # Add the input component
  NewProblem$N["InComp"] = NewProblem$N["InComp"] + NInCompAdd
  NewProblem$InCompName = c(NewProblem$InCompName,
                            trimws(as.character(InCompName)))

  # Add the component
  NewProblem = AddComponents(ThisProblem = NewProblem,
                             CompName = InCompName,
                             CompCharge = InCompCharge,
                             CompMCName = InCompMCName,
                             CompType = InCompType,
                             CompActCorr = InCompActCorr,
                             CompSiteDens = 1.0)

  return(NewProblem)

}

AddDefComp = function(ThisProblem, DefCompName, DefCompFromNum = NULL,
                      DefCompFromVar = NULL, DefCompCharge, DefCompMCName = NULL,
                      DefCompType, DefCompActCorr, DefCompSiteDens,
                      DefCompMCR = match(DefCompMCName, ThisProblem$Mass$Name)) {

  NewProblem = ThisProblem

  # error checking
  if (any(DefCompName %in% ThisProblem$Spec$Name)) {
    stop(paste0(
      "Input Component(s) (",
      paste(paste0("\"", DefCompName[DefCompName %in% ThisProblem$Spec$Name], "\""),
            collapse = ", "),
      ") already exist as a component or species."
    ))
  }
  if (any(is.na(c(DefCompName, DefCompCharge, DefCompMCName,
                  DefCompType, DefCompActCorr, DefCompMCR, DefCompSiteDens)))) {
    stop("NA arguments not allowed, except for DefCompFromNum and DefCompFromVar.")
  }
  if (!all(DefCompMCName %in% ThisProblem$Mass$Name) ||
      !all(DefCompMCR <= ThisProblem$N["Mass"])) {
    stop("Mass compartment(s) specified in DefCompMCName or DefCompMCR does not exist.")
  }
  if (is.null(DefCompFromVar)){
    DefCompFromVar = array(NA, dimnames = list(DefCompName))
  }
  if (is.null(DefCompFromNum)){
    DefCompFromNum = array(NA, dimnames = list(DefCompName))
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

  # Add the component
  NewProblem = AddComponents(ThisProblem = NewProblem,
                             CompName = DefCompName,
                             CompCharge = DefCompCharge,
                             CompMCName = DefCompMCName,
                             CompType = DefCompType,
                             CompActCorr = DefCompActCorr,
                             CompSiteDens = DefCompSiteDens,
                             CompMCR = DefCompMCR)

  return(NewProblem)

}

AddSpecies = function(ThisProblem, SpecName, SpecMCName = NULL, SpecActCorr,
                      SpecCompNames = list(), SpecCompStoichs = list(),
                      SpecLogK, SpecDeltaH, SpecTempKelvin,
                      SpecMCR = match(SpecMCName, ThisProblem$Mass$Name)) {

  NewProblem = ThisProblem

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
                  SpecActCorr, SpecMCR)))) {
    stop("NA arguments not allowed.")
  }
  if(is.null(SpecMCName)) {
    SpecMCName = ThisProblem$Mass$Name[SpecMCR]
  }
  if (!all(SpecMCName %in% ThisProblem$Mass$Name[SpecMCR])) {
    stop("SpecMCName does not match SpecMCR.")
  }
  NSpecAdd = length(SpecName)
  if ((NSpecAdd != length(SpecCompNames)) ||
      (NSpecAdd != length(SpecCompStoichs))) {
    stop("Inputs differ in length.  Must specify SpecName, SpecCompNames, and SpecCompStoichs for each reaction.")
  }
  SpecNC = sapply(SpecCompNames, FUN = length)
  if (any(SpecNC != sapply(SpecCompStoichs, FUN = length))) {
    stop("Stoichiometric inputs differ in length. Specify a matching set of names and values in SpecCompNames and SpecCompStoichs.")
  }

  # add Species
  NewProblem$Spec = rbind(ThisProblem$Spec,
                          data.frame(
                            Name = trimws(as.character(SpecName)),
                            Equation = NA,
                            Charge = NA,
                            MCName = trimws(as.character(SpecMCName)),
                            MCR = as.integer(SpecMCR),
                            ActCorr = trimws(as.character(SpecActCorr)),
                            LogK = as.numeric(SpecLogK),
                            K = 10^as.numeric(SpecLogK),
                            DeltaH = as.numeric(SpecDeltaH),
                            TempKelvin = as.numeric(SpecTempKelvin),
                            NC = as.integer(SpecNC)
                          ))
  NewProblem$SpecCompList = rbind(
    ThisProblem$SpecCompList,
    matrix(0L, nrow = NSpecAdd, ncol = ncol(ThisProblem$SpecCompList),
           dimnames = list(SpecName, NULL))
  )
  if (max(NewProblem$Spec$NC) > ncol(ThisProblem$SpecCompList)) {
    while (max(NewProblem$Spec$NC) > ncol(NewProblem$SpecCompList)) {
      NewProblem$SpecCompList = cbind(
        NewProblem$SpecCompList,
        matrix(0L, nrow = nrow(NewProblem$SpecCompList), ncol = 1)
      )
    }
  }
  NewProblem$SpecStoich = rbind(
    ThisProblem$SpecStoich,
    matrix(0L, nrow = NSpecAdd, ncol = ThisProblem$N["Comp"],
           dimnames = list(SpecName, ThisProblem$Comp$Name))
  )
  for (i in 1:NSpecAdd) {
    ii = i + ThisProblem$N["Spec"]
    SpecCompIndexes = match(SpecCompNames[[i]], ThisProblem$Comp$Name)
    NewProblem$SpecCompList[ii, 1:SpecNC[i]] = sort(SpecCompIndexes)
    NewProblem$SpecStoich[ii, SpecCompIndexes] =
      as.integer(SpecCompStoichs[[i]])
  }
  NewProblem$N["Spec"] = ThisProblem$N["Spec"] + NSpecAdd

  # These parts it's easier to just calculate again
  NewProblem$Spec$Charge = drop(NewProblem$SpecStoich %*% NewProblem$Comp$Charge)
  NewProblem$Spec$Equation = StoichToEquation(SpecStoich = NewProblem$SpecStoich,
                                              SpecName = NewProblem$Spec$Name,
                                              CompName = NewProblem$Comp$Name)

  # Put components at the beginning
  CompIndexes = match(NewProblem$Comp$Name, NewProblem$Spec$Name)
  SpecIndexes = which(NewProblem$Spec$Name %in% NewProblem$Comp$Name == FALSE)
  for (i in names(NewProblem)[grepl("^Spec", names(NewProblem))]) {
    NewProblem[[i]] = rbind(NewProblem[[i]][CompIndexes, , drop = FALSE],
                            NewProblem[[i]][SpecIndexes, , drop = FALSE])
  }

  return(NewProblem)

}

AddPhases = function(ThisProblem, PhaseName, PhaseCompNames = list(),
                     PhaseCompStoichs = list(),
                     PhaseLogK, PhaseDeltaH, PhaseTempKelvin, PhaseMoles) {

  NewProblem = ThisProblem

  # error checking
  if (any(PhaseName %in% ThisProblem$Phase$Name)) {
    stop(paste0(
      "Phase(s) (",
      paste(paste0("\"", PhaseName[PhaseName %in% ThisProblem$Phase$Name], "\""),
            collapse = ", "),
      ") already exist as a phase."
    ))
  }
  if (any(PhaseName %in% ThisProblem$Spec$Name)) {
    stop(paste0(
      "Phase(s) (",
      paste(paste0("\"", PhaseName[PhaseName %in% ThisProblem$Spec$Name], "\""),
            collapse = ", "),
      ") already exist as a component or species"
    ))
  }
  if (any(is.na(c(PhaseName, PhaseLogK, PhaseDeltaH,
                  PhaseTempKelvin, PhaseMoles)))) {
    stop("NA arguments not allowed.")
  }
  NPhaseAdd = length(PhaseName)
  if ((NPhaseAdd != length(PhaseCompNames)) ||
      (NPhaseAdd != length(PhaseCompStoichs))) {
    stop("Inputs differ in length.  Must specify PhaseName, PhaseCompNames, and PhaseCompStoichs for each reaction.")
  }
  PhaseNC = sapply(PhaseCompNames, FUN = length)
  if (any(PhaseNC != sapply(PhaseCompStoichs, FUN = length))) {
    stop("Stoichiometric inputs differ in length. Specify a matching set of names and values in PhaseCompNames and PhaseCompStoichs.")
  }

  # add Phases
  NewProblem$Phase = rbind(ThisProblem$Phase,
                           data.frame(
                             Name = trimws(as.character(PhaseName)),
                             Equation = NA,
                             NC = as.integer(PhaseNC),
                             LogK = as.numeric(PhaseLogK),
                             K = 10^as.numeric(PhaseLogK),
                             DeltaH = as.numeric(PhaseDeltaH),
                             TempKelvin = as.numeric(PhaseTempKelvin),
                             Moles = as.numeric(PhaseMoles)
                           ))
    NewProblem$PhaseCompList = rbind(
    ThisProblem$PhaseCompList,
    matrix(0L, nrow = NPhaseAdd, ncol = ncol(ThisProblem$PhaseCompList),
           dimnames = list(PhaseName, NULL))
  )
  if (max(NewProblem$Phase$NC) > ncol(ThisProblem$PhaseCompList)) {
    while (max(NewProblem$Phase$NC) > ncol(NewProblem$PhaseCompList)) {
      NewProblem$PhaseCompList = cbind(
        NewProblem$PhaseCompList,
        matrix(0L, nrow = nrow(NewProblem$PhaseCompList), ncol = 1)
      )
    }
  }
  NewProblem$PhaseStoich = rbind(
    ThisProblem$PhaseStoich,
    matrix(0L, nrow = NPhaseAdd, ncol = ThisProblem$N["Comp"],
           dimnames = list(PhaseName, ThisProblem$Comp$Name))
  )
  for (i in 1:NPhaseAdd) {
    ii = i + ThisProblem$N["Phase"]
    PhaseCompIndexes = match(PhaseCompNames[[i]], ThisProblem$Comp$Name)
    NewProblem$PhaseCompList[ii, 1:PhaseNC[i]] = sort(PhaseCompIndexes)
    NewProblem$PhaseStoich[ii, PhaseCompIndexes] =
      as.integer(PhaseCompStoichs[[i]])
  }
  NewProblem$N["Phase"] = ThisProblem$N["Phase"] + NPhaseAdd

  NewProblem$Phase$Equation = StoichToEquation(SpecStoich = NewProblem$PhaseStoich,
                                               SpecName = NewProblem$Phase$Name,
                                               CompName = NewProblem$Comp$Name)

  return(NewProblem)

}

AddSpecicalDef = function(ThisProblem, Value, SpecialDef) {

  NewProblem = ThisProblem

  if (any(is.na(Value))) {
    stop("NA inputs not allowed.")
  }
  NSpecialDef = length(Value)
  if (NSpecialDef > length(SpecialDef)) {
    SpecialDef = rep(SpecialDef, NSpecialDef)
  }
  Value = trimws(as.character(Value))

  for (i in 1:NSpecialDef) {
    SpecialDef[i] = match.arg(SpecialDef[i], c("BL", "Metal", "BLMetal", "WHAM"))
    if (SpecialDef[i] == "WHAM") {
      Out = ConvertToList(NewProblem)
      #--Name of WHAM file or WHAM version
      if (any(grepl("WHAM", ThisProblem$InVar$Type))) {
        if (!is.na(ThisProblem$WHAM$WHAMVer) ||
            !is.na(ThisProblem$WHAM$WdatFile)) {
          stop("More than one WHAM version or file specified.")
        }
        Out$DoWHAM = TRUE
        NewProblem$DoWHAM = TRUE
        if (Value[i] %in% c("V", "VI", "VII")) {
          Out$WHAMVer = Value[i]
          Out$WdatFile = NA
          NewProblem$WHAMVer = Value[i]
          NewProblem$WdatFile = NA
        } else {
          Out$WHAMVer = NA
          Out$WdatFile = Value[i]
          NewProblem$WHAMVer = NA
          NewProblem$WdatFile = Value[i]
        }
      } else {
        stop("WHAM version or file specified without a WHAM input variable.")
      }

      # Add WHAM
      Out2 = do.call("ExpandWHAM", args = Out[formalArgs("ExpandWHAM")])
            Out[names(Out2)] = Out2
      NewProblem = ConvertToDF(Out)

    } else if (SpecialDef[i] %in% c("BL", "Metal")) {
      NewProblem[[SpecialDef[i]]] = rbind(NewProblem[[SpecialDef[i]]], data.frame(
        Name = Value[i],
        CompR = match(Value[i], ThisProblem$Comp$Name)
      ))
    } else if (SpecialDef[i] == "BLMetal") {
      NewProblem[[SpecialDef[i]]] = rbind(NewProblem[[SpecialDef[i]]], data.frame(
        Name = Value[i],
        SpecsR = match(Value[i], ThisProblem$Spec$Name)
      ))
    }
  }

  return(NewProblem)

}


RemoveMassCompartmentsDF = function(ThisProblem, MCToRemove) {

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
    NewProblem = RemoveInVar(NewProblem, InVarToRemove)
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
  NewProblem$Mass$Name = ThisProblem$Mass$Name[-MCToRemove]
  NewProblem$Mass$Amt = ThisProblem$Mass$Amt[-MCToRemove]
  NewProblem$Mass$Unit = ThisProblem$Mass$Unit[-MCToRemove]
  NewProblem$N["Mass"] = ThisProblem$N["Mass"] - length(MCToRemove)

  # Update mass compartment indices
  NewProblem$InVar$MCR = match(NewProblem$InVar$MCName, NewProblem$Mass$Name)
  NewProblem$Comp$MCR = match(NewProblem$Comp$MCName, NewProblem$Mass$Name)
  NewProblem$DefComp$MCR = match(NewProblem$DefComp$MCName, NewProblem$Mass$Name)
  NewProblem$Spec$MCR = match(NewProblem$Spec$MCName, NewProblem$Mass$Name)
  NewProblem$Index$AqueousMCR =
    which(tolower(NewProblem$Mass$Name) %in% c("water", "aqueous"))
  NewProblem$Index$BioticLigMCR =
    which(grepl("BL", NewProblem$Mass$Name, ignore.case = TRUE) |
            grepl("gill", NewProblem$Mass$Name, ignore.case = TRUE))
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

RemoveInLab = function(ThisProblem, InLabToRemove) {

  NewProblem = ThisProblem

  InLabToRemoveOrig = InLabToRemove
  if (is.character(InLabToRemove)) {
    InLabToRemove = which(ThisProblem$InLabName %in% InLabToRemove)
  }
  if (length(InLabToRemove) < 1) {
    stop(paste0("Mass Compartment \"", InLabToRemoveOrig, "\" does not exist."))
  }
  if (any(InLabToRemove > ThisProblem$N["InLab"])) {
    stop(paste0("There are ", ThisProblem$N["InLab"], " input labels, ",
                "trying to remove the #(",
                InLabToRemove[InLabToRemove > ThisProblem$N["InLab"]],
                ") element(s)."))
  }

  # Remove input labels
  NewProblem$InLabName = NewProblem$InLabName[-InLabToRemove]
  NewProblem$N["InLab"] = ThisProblem$N["InLab"] - length(InLabToRemove)

  return(NewProblem)

}

RemoveInVar = function(ThisProblem, InVarToRemove) {

  NewProblem = ThisProblem

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
  NewProblem$InVar$Name = NewProblem$InVar$Name[-InVarToRemove]
  NewProblem$InVar$MCName = NewProblem$InVar$MCName[-InVarToRemove]
  NewProblem$InVar$MCR = NewProblem$InVar$MCR[-InVarToRemove]
  NewProblem$InVar$Type = NewProblem$InVar$Type[-InVarToRemove]
  NewProblem$N["InVar"] = ThisProblem$N["InVar"] - length(InVarToRemove)

  return(NewProblem)

}

RemoveComponents = function(ThisProblem, ComponentToRemove) {

  NewProblem = ThisProblem

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
  if (any(ComponentToRemove > ThisProblem$N["Comp"])) {
    stop(paste0("There are ", ThisProblem$N["Comp"], " Components, ",
                "trying to remove the #(",
                paste(ComponentToRemove[ComponentToRemove > ThisProblem$N["Comp"]],
                      collapse = ", "),
                ") element(s)."))
  }

  # Remove Species that depend on the component
  SpeciesToRemove = unique(c(
    ComponentToRemove,
    which(apply(NewProblem$SpecStoich[, ComponentToRemove, drop = FALSE],
                MARGIN = 1, FUN = function(X){any(X != 0L)}))
  ))
  NewProblem = RemoveSpecies(NewProblem, SpeciesToRemove)

  # Remove DefComps that depend on the component
  DefCompToRemove = which(NewProblem$DefComp$FromVar %in%
                            ThisProblem$Comp$Name[ComponentToRemove])
  if (length(DefCompToRemove) >= 1) {
    NewProblem = RemoveDefComps(NewProblem, DefCompToRemove)
  }

  # Remove InComps that are the component to remove
  InCompToRemove = which(NewProblem$InCompName %in%
                           ThisProblem$Comp$Name[ComponentToRemove])
  if (length(InCompToRemove) >= 1) {
    NewProblem$InCompName = NewProblem$InCompName[-InCompToRemove]
    NewProblem$N["InComp"] = ThisProblem$N["InComp"] - length(InCompToRemove)
  }

  # Remove Components
  NewProblem$Comp$Name = ThisProblem$Comp$Name[-ComponentToRemove]
  NewProblem$Comp$Charge = ThisProblem$Comp$Charge[-ComponentToRemove]
  NewProblem$Comp$MCName = ThisProblem$Comp$MCName[-ComponentToRemove]
  NewProblem$Comp$MCR = ThisProblem$Comp$MCR[-ComponentToRemove]
  NewProblem$Comp$Type = ThisProblem$Comp$Type[-ComponentToRemove]
  NewProblem$Comp$ActCorr = ThisProblem$Comp$ActCorr[-ComponentToRemove]
  NewProblem$Comp$SiteDens = ThisProblem$Comp$SiteDens[-ComponentToRemove]
  NewProblem$N["Comp"] = ThisProblem$N["Comp"] - length(ComponentToRemove)
  NewProblem$SpecStoich = NewProblem$SpecStoich[, -ComponentToRemove, drop = FALSE]

  # Update Special Definitions
  NewProblem$BL$Name = ThisProblem$BL$Name[
    ThisProblem$BL$Name %in% ThisProblem$Comp$Name[ComponentToRemove] == FALSE]
  NewProblem$BL$CompR = match(NewProblem$BL$Name, NewProblem$Comp$Name)
  NewProblem$N["BL"] = length(NewProblem$BL$Name)

  NewProblem$Metal$Name = ThisProblem$Metal$Name[
    ThisProblem$Metal$Name %in% ThisProblem$Comp$Name[ComponentToRemove] == FALSE]
  NewProblem$Metal$CompR = match(NewProblem$Metal$Name, NewProblem$Comp$Name)
  NewProblem$N["Metal"] = length(NewProblem$Metal$Name)

  # Update SpecCompList
  NewProblem$Spec$NC = rowSums(NewProblem$SpecStoich != 0L)
  names(NewProblem$Spec$NC) = NewProblem$Spec$Name
  NewProblem$SpecCompList = t(apply(
    NewProblem$SpecStoich,
    MARGIN = 1,
    FUN = function(X) {
      Tmp = sort(which(X != 0L))
      if (length(Tmp) < max(NewProblem$Spec$NC)) {
        Tmp = c(Tmp, rep(0, max(NewProblem$Spec$NC) - length(Tmp)))
      }
      return(Tmp)
    }
  ))
  rownames(NewProblem$SpecCompList) = NewProblem$Spec$Name

  # Update PhaseCompList
  NewProblem$Phase$NC = rowSums(NewProblem$PhaseStoich != 0L)
  names(NewProblem$Phase$NC) = NewProblem$Phase$Name
  if (NewProblem$N["Phase"] > 0) {
    NewProblem$PhaseCompList = as.matrix(apply(
      NewProblem$PhaseStoich,
      MARGIN = 1,
      FUN = function(X) {
        Tmp = sort(which(X != 0L))
        if (length(Tmp) < max(NewProblem$Phase$NC)) {
          Tmp = c(Tmp, rep(0, max(NewProblem$Phase$NC) - length(Tmp)))
        }
        return(Tmp)
      }
    ))
    rownames(NewProblem$PhaseCompList) = NewProblem$Phase$Name
  }

  return(NewProblem)
}

RemoveSpecies = function(ThisProblem, SpeciesToRemove) {

  NewProblem = ThisProblem

  SpeciesToRemoveOrig = SpeciesToRemove
  if (is.character(SpeciesToRemove)) {
    SpeciesToRemove = which(ThisProblem$Spec$Name %in% SpeciesToRemove)
  }
  if (length(SpeciesToRemove) < 1) {
    stop(paste0("Species \"", SpeciesToRemoveOrig, "\" does not exist."))
  }
  if (any(SpeciesToRemove > ThisProblem$N["Spec"])) {
    stop(paste0("There are ", ThisProblem$N["Spec"], " Species, ",
                "trying to remove the #(",
                paste(SpeciesToRemove[SpeciesToRemove > ThisProblem$N["Spec"]],
                      collapse = ", "),
                ") element(s)."))
  }

  # Remove Species
  NewProblem$Spec$Name = ThisProblem$Spec$Name[-SpeciesToRemove]
  NewProblem$Spec$MCName = ThisProblem$Spec$MCName[-SpeciesToRemove]
  NewProblem$Spec$MCR = ThisProblem$Spec$MCR[-SpeciesToRemove]
  NewProblem$Spec$ActCorr = ThisProblem$Spec$ActCorr[-SpeciesToRemove]
  NewProblem$Spec$NC = ThisProblem$Spec$NC[-SpeciesToRemove]
  NewProblem$SpecCompList = ThisProblem$SpecCompList[-SpeciesToRemove, , drop = FALSE]
  NewProblem$Spec$LogK = ThisProblem$Spec$LogK[-SpeciesToRemove]
  NewProblem$Spec$DeltaH = ThisProblem$Spec$DeltaH[-SpeciesToRemove]
  NewProblem$Spec$TempKelvin = ThisProblem$Spec$TempKelvin[-SpeciesToRemove]
  NewProblem$Spec$Charge = ThisProblem$Spec$Charge[-SpeciesToRemove]
  NewProblem$Spec$K = ThisProblem$Spec$K[-SpeciesToRemove]
  NewProblem$SpecStoich = ThisProblem$SpecStoich[-SpeciesToRemove, , drop = FALSE]
  NewProblem$N["Spec"] = ThisProblem$N["Spec"] - length(SpeciesToRemove)

  # Update Special Definitions
  NewProblem$BLMetal$Name = ThisProblem$BLMetal$Name[
    ThisProblem$BLMetal$Name %in%
      ThisProblem$Spec$Name[SpeciesToRemove] == FALSE]
  NewProblem$BLMetal$SpecsR = match(NewProblem$BLMetal$Name, NewProblem$Spec$Name)
  NewProblem$N["BLMetal"] = length(NewProblem$BLMetal$Name)

  return(NewProblem)

}

RemoveDefComps = function(ThisProblem, DefCompToRemove) {

  NewProblem = ThisProblem

  DefCompToRemoveOrig = DefCompToRemove
  if (is.character(DefCompToRemove)) {
    DefCompToRemove = which(ThisProblem$DefComp$Name %in% DefCompToRemove)
  }
  if (length(DefCompToRemove) < 1) {
    stop(paste0("Defined Component \"", DefCompToRemoveOrig, "\" does not exist."))
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
  NewProblem = RemoveComponents(NewProblem, ComponentToRemove)

  # Remove DefComps
  NewProblem$DefComp$Name = ThisProblem$DefComp$Name[-DefCompToRemove]
  NewProblem$DefComp$FromNum = ThisProblem$DefComp$FromNum[-DefCompToRemove]
  NewProblem$DefComp$FromVar = ThisProblem$DefComp$FromVar[-DefCompToRemove]
  NewProblem$DefComp$Charge = ThisProblem$DefComp$Charge[-DefCompToRemove]
  NewProblem$DefComp$MCName = ThisProblem$DefComp$MCName[-DefCompToRemove]
  NewProblem$DefComp$MCR = ThisProblem$DefComp$MCR[-DefCompToRemove]
  NewProblem$DefComp$Type = ThisProblem$DefComp$Type[-DefCompToRemove]
  NewProblem$DefComp$ActCorr = ThisProblem$DefComp$ActCorr[-DefCompToRemove]
  NewProblem$DefComp$SiteDens = ThisProblem$DefComp$SiteDens[-DefCompToRemove]
  NewProblem$N["DefComp"] = ThisProblem$N["DefComp"] - length(DefCompToRemove)

  return(NewProblem)

}
