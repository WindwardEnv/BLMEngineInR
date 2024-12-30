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

#' @title Convert From a Windows BLM Parameter File
#'
#' @param WindowsParamFile Character string with the file path of the
#'   Windows-format BLM parameter file. Typically will have the extension
#'   ".dat".
#' @param RParamFile (optional) Character string with the file path of the
#'   R-format BLM parameter file to save.
#' @param RWHAMFile (optional) Character string with the file path of the
#'   R-format WHAM parameter file to save.
#' @param MarineFile Boolean value - is this a marine file? If so, it uses a
#'   lower mass value. In the Windows BLM, this is equivalent to using the "/M"
#'   switch. Defaults to `FALSE`.
#'
#' @return The BLMEngineInR-compatible chemistry problem object. If RParamFile
#'   is provided, this will return invisibly.
#' @export
#'
ConvertWindowsParamFile = function(WindowsParamFile, RParamFile = NULL,
                                   RWHAMFile = NULL, MarineFile = FALSE) {

  # Read info from WindowsParamFile
  N = scan(
    file = WindowsParamFile,
    what = integer(),
    sep = ",",
    skip = 3,
    n = 4,
    quiet = TRUE
  )

  CompDF = utils::read.csv(file = WindowsParamFile, skip = 5, nrows = N[1], header = FALSE)
  colnames(CompDF) = c("Name", "Charge", "TypeNum", "ActivityNum", "SiteDen")
  CompDF$MCR = floor(CompDF$TypeNum / 10) + 1L
  CompDF$MCName = c("Water","BL")[CompDF$MCR]
  CompDF$Type = c("MassBal","FixedConc","Substituted","ChargeBal")[CompDF$TypeNum - 10*(CompDF$MCR - 1)]
  CompDF$Activity = c("None","Davies","Debye","MoleFraction","ChargeFraction","Debye")[CompDF$ActivityNum]
  CompDF$Name = gsub("Gill", "BL", CompDF$Name)

  if (N[2] > 0L) {
    SpecDF = utils::read.csv(file = WindowsParamFile,
                             skip = 7 + N[1],
                             nrows = N[2],
                             header = FALSE)
    if (ncol(SpecDF) == (6 + N[1])) {
      colnames(SpecDF) = c("Name", "TypeNum", "ActivityNum", "iTemp", "iMonte",
                           CompDF$Name, "LogK_VarLogK_DeltaH_Temp_Conc")
      SpecDF$LogK_VarLogK_DeltaH_Temp_Conc = trimws(SpecDF$LogK_VarLogK_DeltaH_Temp_Conc)
      while (any(nchar(SpecDF$LogK_VarLogK_DeltaH_Temp_Conc) !=
                nchar(gsub("  "," ", SpecDF$LogK_VarLogK_DeltaH_Temp_Conc)))) {
        SpecDF$LogK_VarLogK_DeltaH_Temp_Conc = gsub("  "," ", SpecDF$LogK_VarLogK_DeltaH_Temp_Conc)
      }
      LogK_VarLogK_DeltaH_Temp_Conc.mat =
        t(simplify2array(strsplit(SpecDF$LogK_VarLogK_DeltaH_Temp_Conc, split = " ")))
      SpecDF$LogK = as.numeric(LogK_VarLogK_DeltaH_Temp_Conc.mat[, 1])
      # SpecDF$VarLogK = as.numeric(LogK_VarLogK_DeltaH_Temp_Conc.mat[, 2])
      SpecDF$DeltaH = as.numeric(LogK_VarLogK_DeltaH_Temp_Conc.mat[, 3])
      SpecDF$TempKelvin = as.numeric(LogK_VarLogK_DeltaH_Temp_Conc.mat[, 4])
      # SpecDF$Conc = as.numeric(LogK_VarLogK_DeltaH_Temp_Conc.mat[, 5])
    } else {
      colnames(SpecDF) = c("Name", "TypeNum", "ActivityNum", "iTemp", "iMonte",
                           CompDF$Name, "LogK","VarLogK","DeltaH","TempKelvin","Conc")
    }
    SpecDF$MCR = floor(SpecDF$TypeNum / 10) + 1L
    SpecDF$MCName = c("Water","BL")[SpecDF$MCR]
    SpecDF$Activity = c("None","Davies","Debye","MoleFraction","ChargeFraction","Debye")[SpecDF$ActivityNum]
    SpecDF$Name = gsub("Gill", "BL", SpecDF$Name)
  }

  # Linked Lists
  if (N[4] > 0L) {
    stop("Linked Lists not implemented.")
  }
  # Phases
  if (N[3] > 0L) {
    PhasesDF = utils::read.csv(file = WindowsParamFile,
                               skip = 11 + N[1] + N[2] + N[4],
                               nrows = N[3],
                               header = FALSE)
    if (ncol(PhasesDF) == (4 + N[1])) {
      colnames(PhasesDF) = c("Name", "iTemp", "iMonte", CompDF$Name,
                             "LogKs_VarLogK_DeltaH_Temp_Moles")
      PhasesDF$LogKs_VarLogK_DeltaH_Temp_Moles = trimws(PhasesDF$LogKs_VarLogK_DeltaH_Temp_Moles)
      while(any(nchar(PhasesDF$LogKs_VarLogK_DeltaH_Temp_Moles) !=
                nchar(gsub("  "," ", PhasesDF$LogKs_VarLogK_DeltaH_Temp_Moles)))){
        PhasesDF$LogKs_VarLogK_DeltaH_Temp_Moles = gsub("  "," ", PhasesDF$LogKs_VarLogK_DeltaH_Temp_Moles)
      }
      LogKs_VarLogK_DeltaH_Temp_Moles.mat =
        t(simplify2array(strsplit(PhasesDF$LogKs_VarLogK_DeltaH_Temp_Moles, split = " ")))
      PhasesDF$LogK = as.numeric(LogKs_VarLogK_DeltaH_Temp_Moles.mat[, 1])
      # PhasesDF$VarLogK = as.numeric(LogKs_VarLogK_DeltaH_Temp_Moles.mat[, 2])
      PhasesDF$DeltaH = as.numeric(LogKs_VarLogK_DeltaH_Temp_Moles.mat[, 3])
      PhasesDF$TempKelvin = as.numeric(LogKs_VarLogK_DeltaH_Temp_Moles.mat[, 4])
      PhasesDF$Moles = as.numeric(LogKs_VarLogK_DeltaH_Temp_Moles.mat[, 5])
    } else {
      colnames(PhasesDF) = c("Name", "iTemp", "iMonte", CompDF$Name, "LogKs",
                             "VarLogK", "DeltaH", "TempKelvin", "Moles")
    }
  }

  # "User Notes" is used information and everything after the [END] flag is just
  # notes for humans to read.
  UserNotes = scan(
    file = WindowsParamFile,
    what = character(),
    skip = sum(N) + 41,
    sep = "\n",
    quiet = TRUE
  )
  TheEnd = which(grepl("[[]END[]]", UserNotes))[1]
  if (TheEnd == length(UserNotes)) {
    AfterTheEnd = NULL
  } else {
    AfterTheEnd = UserNotes[(TheEnd + 1):length(UserNotes)]
  }
  UserNotes = UserNotes[1:(TheEnd)]
  UserNotes = gsub("Gill", "BL", UserNotes)
  Tmp = which(grepl("[[]DOC[]]:", UserNotes))
  if (length(Tmp) > 0) {
    DOCComp = trimws(gsub("^[[]DOC[]]:( )?", "", UserNotes[Tmp]))
  } else {
    DOCComp = NULL
  }


  # Get thermodynamic database information
  if (any(grepl("[[]THERMO[]]: ?[[:digit:]a-zA-Z]", UserNotes))) {
    ThermoDBSName = toupper(trimws(
      gsub("[[]THERMO[]]:", "", UserNotes[which(grepl("THERMO", UserNotes))])
    ))
    # ideally, look in the same directory as the parameter file
    if (file.exists(file.path(dirname(WindowsParamFile), ThermoDBSName))) {
      NewProblem = ConvertWHAMVThermoFile(
        ThermoDBSName = file.path(dirname(WindowsParamFile), ThermoDBSName),
        RWHAMFile = RWHAMFile
      )
    } else if (ThermoDBSName == "WATER23.DBS") {
      NewProblem = BLMEngineInR::All_WATER23_reactions
    } else if (ThermoDBSName == "NIST_20170203.DBS") {
      NewProblem = BLMEngineInR::All_NIST20170203_reactions
    }
  } else {
    NewProblem = BLMEngineInR::water_problem
  }

  # Remove unnecessary components
  if ((NewProblem$N["InComp"] != 0L) &
      any(NewProblem$InCompName %in% CompDF$Name == FALSE)) {
    NewProblem = RemoveInComps(
      ThisProblem = NewProblem,
      InCompToRemove = NewProblem$InCompName[NewProblem$InCompName %in%
                                               CompDF$Name == FALSE]
    )
  }
  if ((NewProblem$N["InDefComp"] != 0L) &
      any(NewProblem$InDefCompName %in% c("H","OH") == FALSE)) {
    NewProblem = RemoveDefComps(
      ThisProblem = NewProblem,
      DefCompToRemove = NewProblem$InDefCompName[NewProblem$InDefCompName %in% c(
        "H","OH"
      ) == FALSE]
    )
  }
  if (any(SpecDF$Name == "OH")) {
    NewProblem = RemoveDefComps(
      ThisProblem = NewProblem,
      DefCompToRemove = "OH"
    )
  }

  # Add extra components not in thermodynamic database
  CompsToAdd = which((CompDF$Name %in% NewProblem$Comp$Name == FALSE) &
                       (CompDF$Name %in% c("DOC") == FALSE) &
                       !grepl("BL", CompDF$Name))
  if (length(CompsToAdd) > 0L) {
    NewProblem = AddComponents(
      ThisProblem = NewProblem,
      CompName = CompDF$Name[CompsToAdd],
      CompCharge = CompDF$Charge[CompsToAdd],
      CompMCName = CompDF$MCName[CompsToAdd],
      CompType = CompDF$Type[CompsToAdd],
      CompActCorr = CompDF$Activity[CompsToAdd]
    )
    if (any(CompDF$SiteDen[CompsToAdd] != 1)) {
      NewProblem$Comp$SiteDens[match(CompDF$Name[CompsToAdd],
                                     NewProblem$Comp$Name)] =
        CompDF$SiteDen[CompsToAdd]
      CheckBLMObject(Object = NewProblem, Reference = BlankProblem())
    }
  }

  # Add DefComps
  DefCompsToAdd = which(grepl("BL", CompDF$Name))
  if (length(DefCompsToAdd) > 0L) {
    if (!any(NewProblem$Mass$Name == "BL")) {
      NewProblem = AddMassCompartments(
        ThisProblem = NewProblem,
        MassName = "BL",
        MassAmt = 1,
        MassUnit = "kg wet"
      )
    }
    NewProblem = AddDefComps(
      ThisProblem = NewProblem,
      DefCompName = CompDF$Name[DefCompsToAdd],
      DefCompFromNum = if (MarineFile) { 1.78E-11 } else {1.78E-5},
      DefCompMCName = CompDF$MCName[DefCompsToAdd],
      DefCompCharge = CompDF$Charge[DefCompsToAdd],
      DefCompType = CompDF$Type[DefCompsToAdd],
      DefCompActCorr = CompDF$Activity[DefCompsToAdd],
      DefCompSiteDens = CompDF$SiteDen[DefCompsToAdd]
    )
  }

  # Species
  if (N[2] > 0L) {
    SpecStoich.mat = as.matrix(SpecDF[, CompDF$Name])
    rownames(SpecStoich.mat) = SpecDF$Name
    if(any(NewProblem$Comp$Name %in% CompDF$Name == FALSE)) {
      SpecStoich.mat = cbind(
        SpecStoich.mat,
        matrix(data = 0L, nrow = N[2],
               ncol = sum(NewProblem$Comp$Name %in% CompDF$Name == FALSE),
               dimnames = list(SpecDF$Name,
                               setdiff(NewProblem$Comp$Name, CompDF$Name)))
      )
    }
    SpecStoich.mat = SpecStoich.mat[, NewProblem$Comp$Name]
    NewProblem = AddSpecies(
      ThisProblem = NewProblem,
      SpecName = SpecDF$Name,
      SpecStoich = SpecStoich.mat,
      SpecMCName = SpecDF$MCName,
      SpecActCorr = SpecDF$Activity,
      SpecLogK = SpecDF$LogK,
      SpecDeltaH = SpecDF$DeltaH,
      SpecTempKelvin = SpecDF$TempKelvin
    )
  }

  # Phases
  if (N[3] > 0L) {
    PhasesStoich.mat = as.matrix(PhasesDF[, CompDF$Name])
    rownames(PhasesStoich.mat) = PhasesDF$Name
    if(any(NewProblem$Comp$Name %in% CompDF$Name == FALSE)) {
      PhasesStoich.mat = cbind(
        PhasesStoich.mat,
        matrix(data = 0L, nrow = N[3],
               ncol = sum(NewProblem$Comp$Name %in% CompDF$Name == FALSE),
               dimnames = list(PhasesDF$Name,
                               setdiff(NewProblem$Comp$Name, CompDF$Name)))
      )
    }
    PhasesStoich.mat = PhasesStoich.mat[, NewProblem$Comp$Name, drop = FALSE]
    NewProblem = AddPhases(
      ThisProblem = NewProblem,
      PhaseName = PhasesDF$Name,
      PhaseStoich = PhasesStoich.mat,
      PhaseLogK = PhasesDF$LogK,
      PhaseDeltaH = PhasesDF$DeltaH,
      PhaseTempKelvin = PhasesDF$TempKelvin,
      PhaseMoles = PhasesDF$Moles
    )
  }

  # DOC input variables
  if (!is.null(DOCComp)) {
    NewProblem = AddInVars(
      ThisProblem = NewProblem,
      InVarName = c(DOCComp, "HA"),
      InVarMCName = "Water",
      InVarType = c("WHAM-HAFA", "PercHA")
    )
    # NewProblem = ExpandWHAM(ThisProblem = NewProblem,
    #                         ThisWHAM = NewProblem$WHAM)
  }

  # Special Defs
  SpecialDefs.mat = trimws(simplify2array(strsplit(
    UserNotes[grepl("[[]Metal[]]:", UserNotes) |
                grepl("[[]BL[]]:", UserNotes) |
                grepl("[[]BL-Metal[]]:", UserNotes)], split = ":")))
  SpecialDefs.mat[1, ] = gsub("[[]", "", SpecialDefs.mat[1, ])
  SpecialDefs.mat[1, ] = gsub("[]]", "", SpecialDefs.mat[1, ])
  NewProblem = AddSpecialDefs(
    ThisProblem = NewProblem,
    Value = SpecialDefs.mat[2, ],
    SpecialDef = SpecialDefs.mat[1, ]
  )

  # Critical table
  CritFlags = which(grepl("CRITICAL", UserNotes) | (grepl("LA50", UserNotes)))
  if (length(CritFlags) == 1L) {
    CATab = data.frame(
      CA = as.numeric(gsub("[[]LA50[]]:", "",
                           gsub("[[]CRITICAL[]]:", "", UserNotes[CritFlags]))),
      Test.Type = "unknown"
    )
    if (any(grepl("HC5_LABEL", UserNotes))) {
      CATab$Endpoint = gsub("[[]HC5_LABEL[]]: ?", "",
                            UserNotes[grepl("HC5_LABEL", UserNotes)])
    } else if (any(grepl("ACUTE_DIV_LABEL", UserNotes))) {
      CATab$Endpoint = gsub("[[]ACUTE_DIV_LABEL[]]: ?", "",
                            UserNotes[grepl("ACUTE_DIV_LABEL", UserNotes)])
    } else if (any(grepl("CHRONIC_DIV_LABEL", UserNotes))) {
      CATab$Endpoint = gsub("[[]CHRONIC_DIV_LABEL[]]: ?", "",
                            UserNotes[grepl("CHRONIC_DIV_LABEL", UserNotes)])
    } else if (any(grepl("ACUTE", UserNotes) | grepl("CHRONIC", UserNotes))) {
      CATab$Endpoint = "WQS"
    }
    if (any(grepl("ACUTE_DIV", UserNotes))) {
      CATab$Test.Type = "Acute"
      CATab$Duration = gsub("[[]ACUTE_DIV[]]: ?", "DIV=",
                            UserNotes[grepl("[[]ACUTE_DIV[]]", UserNotes)])
    }
    if (any(grepl("CHRONIC_DIV", UserNotes))) {
      CATab$Test.Type = "Acute"
      CATab$Lifestage = gsub("[[]CHRONIC_DIV[]]: ?", "ACR=",
                             UserNotes[grepl("CHRONIC_DIV", UserNotes)])
    }
    CATab$Miscellaneous = paste0("Single critical value from \"",
                                 WindowsParamFile,
                                 "\" Windows BLM parameter file. See Notes for details")
  } else if (length(CritFlags) == 2L) {
    stopifnot(grepl("CRITICAL START", UserNotes[CritFlags[1]]) &
                grepl("CRITICAL END", UserNotes[CritFlags[2]]))
    CATab.text = UserNotes[(CritFlags[1] + 2):(CritFlags[2] - 1)]
    CATab = utils::read.csv(
      textConnection(CATab.text),
      header = FALSE,
      col.names = c(
        "CA",
        "Species",
        "Test.Type",
        "Lifestage",
        "Endpoint",
        "Quantifier",
        "References",
        "Miscellaneous"
      ),
    )
    CATab$Duration = NA
    i.bool = grepl("test duration:", CATab$Miscellaneous) & is.na(CATab$Duration)
    if (any(i.bool)) {
      CATab$Duration[i.bool] = gsub("; .+","", gsub("(.*);? ?test duration: (.+)", "\\2", CATab$Miscellaneous[i.bool]))
      CATab$Miscellaneous[i.bool] = gsub("test duration: .+", "", gsub("; test duration: .+$","",gsub("test duration: .+; ", "", CATab$Miscellaneous[i.bool])))
    }
    i.bool = grepl(" [(].+[)]", CATab$Test.Type) & is.na(CATab$Duration)
    if (any(i.bool)) {
      CATab$Duration[i.bool] = gsub("(.+) [(](.+)[)]", "\\2", CATab$Test.Type[i.bool])
      CATab$Test.Type[i.bool] = gsub("(.+) [(](.+)[)]", "\\1", CATab$Test.Type[i.bool])
    }
    i.bool = grepl("[,-] ", CATab$Test.Type) & is.na(CATab$Duration)
    if (any(i.bool)) {
      CATab$Duration[i.bool] = gsub("(.+) [,-](.+)", "\\2", CATab$Test.Type[i.bool])
      CATab$Test.Type[i.bool] = gsub("(.+) [,-](.+)", "\\1", CATab$Test.Type[i.bool])
    }
  } else {
    stop("Unknown critical value reporting.")
  }
  if (length(CritFlags) > 0L) {
    NewProblem = AddCriticalValues(
      ThisProblem = NewProblem,
      CATab = CATab
    )
  }

  # Labels
  NewProblem = RemoveInLabs(ThisProblem = NewProblem,
                            InLabToRemove = 1:NewProblem$N["InLab"])
  NewProblem = AddInLabs(
    ThisProblem = NewProblem,
    InLabName = c("Site Label", "Sample Label")
  )

  NewProblem$Notes = AfterTheEnd

  if (!is.null(RParamFile)) {
    NewProblem = WriteParamFile(ThisProblem = NewProblem,
                                ParamFile = RParamFile,
                                Notes = AfterTheEnd)
    return(invisible(NewProblem))
  } else {
    return(NewProblem)
  }

}
