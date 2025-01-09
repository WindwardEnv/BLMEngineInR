ConvertWHAM7Files = function(DB7File, HAph7File, FAph7File,
                             RParamFile = NULL, RWHAMFile = NULL) {

  NewProblem = water_problem
  NewWHAM = BlankWHAM()

  # DB7File --------------------------------------------------------------------
  NHeaderLines =
    read.csv(file = DB7File, header = FALSE, skip = 1, nrows = 1)[1, 2]
  Header =
    read.csv(file = DB7File, header = TRUE, skip = 2, nrows = NHeaderLines - 1)

  NColMain = read.csv(
    file = DB7File,
    header = FALSE,
    skip = NHeaderLines + 3,
    nrows = 1
  )[1, 2]
  MainColHeader = as.character(read.csv(
    file = DB7File,
    header = FALSE,
    skip = NHeaderLines + 4,
    nrows = 1
  ))
  NMainLines = read.csv(
    file = DB7File,
    header = FALSE,
    skip = NHeaderLines + 6,
    nrows = 1
  )[1, 2]
  Main = read.csv(
    file = DB7File,
    header = FALSE,
    skip = NHeaderLines + 7,
    nrows = NMainLines,
    col.names = MainColHeader
  )
  MainComp = Main[Main[, 1] < 1000, ]
  NMainComp = nrow(MainComp)
  MainSpecColHeader = setdiff(MainColHeader, "MolWt")
  MainSpec = read.csv(
    file = DB7File,
    header = FALSE,
    skip = NHeaderLines + 7 + NMainComp,
    nrows = NMainLines - NMainComp,
    col.names = MainSpecColHeader
  )

  Tmp = paste(MainSpec$Name, "=")
  MainCompNames = MainComp$Name
  MainCompNum = MainComp$Number
  MaxNC = sum(grepl("Comp", MainSpecColHeader))
  for (i in 1:MaxNC) {
    CompColi = MainSpec[, paste0("Comp", i)]
    HasCompi = (CompColi != 0L)
    CompColi = CompColi[HasCompi]
    StoiColi = MainSpec[HasCompi, paste0("Stoi", i)]
    CompNamei = MainCompNames[match(CompColi, MainCompNum)]
    if (i != 1) {
      Tmp[HasCompi] = paste(Tmp[HasCompi], "+")
    }
    Tmp[HasCompi] = paste(Tmp[HasCompi], StoiColi, "*", CompNamei)
  }
  MainSpec$Equation = Tmp
  MainSpec$deltaH_J.mol = MainSpec$deltaH * 220 * 8.314 * log(10)

  # TODO: CO2, HCO3, and H2CO3 formation reactions have their enthalpy
  # correction done with an equation other than Van't Hoff: log10(K) = log10(K0)
  # + (A*(T^2) + B*T + C), where T is in degrees Celsius. If we want to match
  # this in BLM in R, we'll need to change how we specify it in the parameter
  # file and update the code to do this calculation. The below code is a
  # band-aid, and will ensure that the enthalpy correction is pretty close
  # between 15 and 25 degrees Celsius, but it starts to deviate pretty quickly
  # outside of that range.
  i = which(MainSpec$Name == "HCO3")
  Ai = as.numeric(Header[Header[, 1] == "A term in dH for HCO3", 3])
  Bi = as.numeric(Header[Header[, 1] == "B term in dH for HCO3", 3])
  Ci = as.numeric(Header[Header[, 1] == "C term in dH for HCO3", 3])
  MainSpec$deltaH_J.mol[i] =
    (Ai * (15^2) - Bi * 15 + Ci) * 8.341 * log(10) /
    (0.003354 - (1 / (15 + 273.15)))

  i = which(MainSpec$Name == "H2CO3")
  Ai = as.numeric(Header[Header[, 1] == "A term in dH for H2CO3", 3])
  Bi = as.numeric(Header[Header[, 1] == "B term in dH for H2CO3", 3])
  Ci = as.numeric(Header[Header[, 1] == "C term in dH for H2CO3", 3])
  MainSpec$deltaH_J.mol[i] =
    (Ai * (15^2) - Bi * 15 + Ci) * 8.341 * log(10) /
    (0.003354 - (1 / (15 + 273.15)))

  CompToAdd = which(MainComp$Name %in% c("H", "OH") == FALSE)
  NewProblem = AddInComps(
    ThisProblem = NewProblem,
    InCompName = MainComp$Name[CompToAdd],
    InCompCharge = MainComp$Charge[CompToAdd],
    InCompType = "MassBal",
    InCompMCName = "Water",
    InCompActCorr = "Debye"
  )
  NewProblem = AddSpecies(
    ThisProblem = NewProblem,
    SpecEquation = MainSpec$Equation,
    SpecMCName = "Water",
    SpecActCorr = "Debye",
    SpecLogK = MainSpec$logKf,
    SpecDeltaH = MainSpec$deltaH_J.mol,
    SpecTempKelvin = 1 / 0.003354
  )

  NewWHAM$DLF =
    as.numeric(Header[Header[, 1] == "Double layer overlap factor", 3])


  # PH7Files -------------------------------------------------------------------
  PH7Files = c(HAph7File, FAph7File)
  for (i in 1:2) {

    NHeaderLines =
      read.csv(file = PH7Files[i], header = FALSE, skip = 1, nrows = 1)[1, 2]
    Header = read.csv(
      file = PH7Files[i],
      header = TRUE,
      skip = 2,
      nrows = NHeaderLines - 1
    )

    NColMain = read.csv(
      file = PH7Files[i],
      header = FALSE,
      skip = NHeaderLines + 3,
      nrows = 1
    )[1, 2]
    MainColHeader = as.character(read.csv(
      file = PH7Files[i],
      header = FALSE,
      skip = NHeaderLines + 4,
      nrows = 1
    ))
    NMainLines = read.csv(
      file = PH7Files[i],
      header = FALSE,
      skip = NHeaderLines + 6,
      nrows = 1
    )[1, 2]
    Main = read.csv(
      file = PH7Files[i],
      header = FALSE,
      skip = NHeaderLines + 7,
      nrows = NMainLines,
      col.names = MainColHeader
    )

    WHAMLabel = c("nA (sites per g)", "pKHA", "pKHB", "dpKHA", "dpKHB", "P",
                  "proportion of sites in bidentate pairs",
                  "proportion of sites in tridentate triplets",
                  "molecular radius (m)", "molecular weight (daltons)")
    RLabel = c("nA", "pKA", "pKB", "dpKA", "dpKB", "P", "fprB", "fprT",
               "Radius", "MolWt")
    for (ii in 1:length(RLabel)) {
      NewWHAM[[RLabel[ii]]][i] =
        as.numeric(Header[Header[, 1] == WHAMLabel[ii], 3])
    }

    NewWHAM$dLK1A[i] = 0.0
    NewWHAM$dLK1B[i] = 0.0

    if (any(Main[, 1] %in% NewWHAM$MetalsTable$Metal == FALSE)) {
      NewWHAM$MetalsTable = rbind(
        NewWHAM$MetalsTable,
        data.frame(Metal = setdiff(Main[, 1], NewWHAM$MetalsTable$Metal),
                   pKMAHA = NA_real_, pKMAFA = NA_real_, dLK2 = NA_real_)
      )
    }
    RIndex = match(Main[, 1], NewWHAM$MetalsTable$Metal)
    NewWHAM$MetalsTable[RIndex, i + 1] = Main$logK_MA
    NewWHAM$MetalsTable$dLK2[RIndex] = Main$deltaLK2

    HasKsel = (Main$logKsel != 0)
    if (any(Main[HasKsel, 1] %in% NewWHAM$SpecKselTable$Spec == FALSE)) {
      NewWHAM$SpecKselTable = rbind(
        NewWHAM$SpecKselTable,
        data.frame(Spec = setdiff(Main[HasKsel, 1], NewWHAM$SpecKselTable$Spec),
                   KselHA = 1.0, KselFA = 1.0)
      )
    }
    RIndex = match(Main[, 1], NewWHAM$SpecKselTable$Spec)[HasKsel]
    NewWHAM$SpecKselTable[RIndex, i + 1] = 10 ^ Main$logKsel[HasKsel]
    # TODO: change the Ksel numbers and usage to logKsel instead

  }
  NewWHAM$KZED = as.numeric(
    Header[Header[, 1] == "constant to adjust DL volume at low charge", 3]
  )

  # Final adjustments ----------------------------------------------------------
  NewWHAM$File = NA_character_
  NewWHAM$Ver = "VII"
  WHAMVII = DefineWHAM(WHAMVer = "VII")
  NewWHAM$MonodentTable = WHAMVII$MonodentTable
  NewWHAM$BidentTable = WHAMVII$BidentTable
  NewWHAM$TridentTable = WHAMVII$TridentTable
  NewProblem = ExpandWHAM(ThisProblem = NewProblem, ThisWHAM = NewWHAM)

  CheckBLMObject(Object = NewProblem, Reference = BlankProblem())

  if (!is.null(RWHAMFile) || !is.null(RParamFile)) {
    if (!is.null(RWHAMFile)) {
      WriteWHAMFile(ThisWHAM = NewWHAM,
                    WHAMFile = RWHAMFile)
      NewProblem$WHAM$File = basename(RWHAMFile)
      NewWHAM$File = RWHAMFile
    }
    if (!is.null(RParamFile)) {
      NewProblem = WriteParamFile(ThisProblem = NewProblem,
                                  ParamFile = RParamFile)
    }
    return(invisible(NewProblem))
  } else {
    return(NewProblem)
  }

}
