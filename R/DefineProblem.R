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

#' Define the speciation problem
#'
#' `DefineProblem` reads in a parameter file, and sets up the required vectors
#' and matrices that will be needed to run the speciation calculations in CHESS.
#'
#' @param ParamFile the path and file name to a parameter file
#' @param WriteLog if TRUE, the CHESS.LOG file will be written, summarizing the
#'   current problem
#'
#' @returns Returns a `list` object with each list item named according to the
#'   template of BlankProblem
#'
#' @export
#'
#' @examples
#' mypfile = system.file(file.path("extdata","ParameterFiles","carbonate_system_only.dat4"),
#'                       package = "BLMEngineInR", mustWork = TRUE)
#' thisProblem = DefineProblem(mypfile)
#'
DefineProblem = function(ParamFile, WriteLog = FALSE) {

  # error catching
  stopifnot(file.exists(ParamFile))

  # read in parameter file
  # -get number of mass compartments
  # -get number of input variables
  # -get number of input components
  # -get number of defined components
  # -get number of species
  # -get number of phases
  # -get component information

  NewProblem = BlankProblem()
  NewProblem$ParamFile = ParamFile

  # read the dimensions of the various elements of the reaction list
  SkipRows = 2
  Tmp = read.csv(file = ParamFile, header = FALSE, skip = SkipRows,
                 nrows = 9, strip.white = TRUE)
  NMass = as.integer(Tmp[1, 1])
  NInLab = as.integer(Tmp[2, 1])
  NInVar = as.integer(Tmp[3, 1])
  NInComp = as.integer(Tmp[4, 1])
  NDefComp = as.integer(Tmp[5, 1])
  NSpec = as.integer(Tmp[6, 1])
  NPhase = as.integer(Tmp[7, 1])
  NSpecialDef = as.integer(Tmp[8, 1])
  NCAT = as.integer(Tmp[9, 1])
  stopifnot(NMass > 0, NInLab > 0, NInVar > 0, NInComp > 0, NSpec > 0)

  # read mass compartment list
  SkipRows = SkipRows + 9 + 2
  Tmp = read.csv(file = ParamFile, header = TRUE, skip = SkipRows,
                 nrows = NMass, strip.white = TRUE)
  NewProblem = AddMassCompartments(ThisProblem = NewProblem,
                                   MassName = as.character(trimws(Tmp[, 1])),
                                   MassAmt = as.numeric(Tmp[, 2]),
                                   MassUnit = as.character(trimws(Tmp[, 3])),
                                   DoCheck = FALSE)

  # read input Labels
  SkipRows = SkipRows + NMass + 3
  Tmp = scan(file = ParamFile, what = character(), nlines = NInLab, sep = "\n",
             skip = SkipRows, quiet = TRUE, strip.white = TRUE)
  NewProblem = AddInLabs(ThisProblem = NewProblem,
                         InLabName = as.character(Tmp),
                         DoCheck = FALSE)

  # read input variables -
  SkipRows = SkipRows + NInLab + 2
  Tmp = read.csv(file = ParamFile, header = TRUE, skip = SkipRows,
                 nrows = NInVar, strip.white = TRUE)
  NewProblem = AddInVars(ThisProblem = NewProblem,
                         InVarName = as.character(trimws(Tmp[, 1])),
                         InVarMCName = as.character(trimws(Tmp[, 2])),
                         InVarType = as.character(trimws(Tmp[, 3])),
                         DoCheck = FALSE)
  stopifnot("Temperature" %in% NewProblem$InVar$Type)
  # - Temperature = the temperature in degrees C
  # - pH = the -log[H]...you know, pH
  # - WHAM-HA, -FA, -HAFA = Windemere Humic Aqueous Model organic matter (input
  #   mg C/L), as all humic acid, all fulvic acid, or a mix of humics and
  #   fulvics, respectively.
  # - PercHA = optionally indicate the percent humic acid in a the WHAM-HAFA
  #   component for that compartment.
  # - PercAFA = optionally indicate the percent of active fulvic acid for the
  #   WHAM-FA or WHAM-HAFA component for that compartment

  # read component list and properties
  SkipRows = SkipRows + NInVar + 3
  Tmp = read.csv(file = ParamFile, header = TRUE, skip = SkipRows,
                 nrows = NInComp, strip.white = TRUE)
  NewProblem = AddInComps(ThisProblem = NewProblem,
                          InCompName = as.character(trimws(Tmp[, 1])),
                          InCompCharge = as.integer(Tmp[, 2]),
                          InCompMCName = as.character(trimws(Tmp[, 3])),
                          InCompType = as.character(trimws(Tmp[, 4])),
                          InCompActCorr = as.character(trimws(Tmp[, 5])),
                          DoCheck = FALSE)

  # read defined component list and properties
  SkipRows = SkipRows + NInComp + 3
  if (NDefComp > 0) {
    Tmp = read.csv(file = ParamFile, header = TRUE, skip = SkipRows,
                   nrows = NDefComp, strip.white = TRUE)
    DefCompName = as.character(trimws(Tmp[, 1]))
    NumericDefComp = grepl("^[[:digit:].]+[eE]?[+-]?[[:digit:]]?", Tmp[, 2])
    DefCompFromNum = as.numeric(array(NA, dim = NDefComp))
    DefCompFromNum[NumericDefComp] = as.numeric(Tmp[NumericDefComp, 2])
    DefCompFromVar = as.character(array(NA, dim = NDefComp))
    DefCompFromVar[!NumericDefComp] = trimws(Tmp[!NumericDefComp, 2])
    DefCompFromVar[!is.na(DefCompFromNum)] = NA
    DefCompFromNum[!is.na(DefCompFromVar)] = NA
    if (any((c("H", "OH") %in% NewProblem$Comp$Name) &
            (c("H", "OH") %in% DefCompName))) {
      iComp = match(c("H", "OH"), NewProblem$Comp$Name)
      iNewDefComp = match(c("H", "OH"), DefCompName)
      iOldDefComp = match(c("H", "OH"), NewProblem$DefComp$Name)
      NewProblem$Comp$Charge[iComp] = as.integer(trimws(Tmp$Charge[iNewDefComp]))
      NewProblem$Comp$MCName[iComp] = as.character(trimws(Tmp$Compartment[iNewDefComp]))
      NewProblem$Comp$MCR[iComp] = match(NewProblem$Comp$MCName[iComp], NewProblem$Mass$Name)
      NewProblem$Comp$Type[iComp] = as.character(trimws(Tmp$Type[iNewDefComp]))
      NewProblem$Comp$ActCorr[iComp] = as.character(trimws(Tmp$Activity[iNewDefComp]))
      NewProblem$Comp$SiteDens[iComp] = as.numeric(trimws(Tmp$Site.Den[iNewDefComp]))
      NewProblem$DefComp$FromNum[iOldDefComp] = DefCompFromNum[iNewDefComp]
      NewProblem$DefComp$FromVar[iOldDefComp] = DefCompFromVar[iNewDefComp]
      NewProblem$DefComp$Charge[iOldDefComp] = NewProblem$Comp$Charge[iComp]
      NewProblem$DefComp$MCName[iOldDefComp] = NewProblem$Comp$MCName[iComp]
      NewProblem$DefComp$MCR[iOldDefComp] = NewProblem$Comp$MCR[iComp]
      NewProblem$DefComp$Type[iOldDefComp] = NewProblem$Comp$Type[iComp]
      NewProblem$DefComp$ActCorr[iOldDefComp] = NewProblem$Comp$ActCorr[iComp]
      NewProblem$DefComp$SiteDens[iOldDefComp] = NewProblem$Comp$SiteDens[iComp]
      Tmp = Tmp[-iNewDefComp, ]
      DefCompName = DefCompName[-iNewDefComp]
      DefCompFromNum = DefCompFromNum[-iNewDefComp]
      DefCompFromVar = DefCompFromVar[-iNewDefComp]
    }
    if (length(DefCompName) > 0) {
      NewProblem = AddDefComps(
        ThisProblem = NewProblem,
        DefCompName = DefCompName,
        DefCompFromNum = DefCompFromNum,
        DefCompFromVar = DefCompFromVar,
        DefCompCharge = as.integer(Tmp[, 3]),
        DefCompMCName = as.character(trimws(Tmp[, 4])),
        DefCompType = as.character(trimws(Tmp[, 5])),
        DefCompActCorr = as.character(trimws(Tmp[, 6])),
        DefCompSiteDens = as.numeric(Tmp[, 7]),
        DoCheck = FALSE
      )
    }
  }

  # Check that pH is in the component list
  for (iMass in 1:NMass){#Must have either pH or H in non-BL compartments
    if (grepl("Water", NewProblem$Mass$Name[iMass], ignore.case = TRUE)) {
      if (!xor(("H" %in% NewProblem$InCompName),
                    ("pH" %in% NewProblem$InVar$Type[NewProblem$InVar$MCR == iMass]) &
                      ("H" %in% NewProblem$DefComp$Name[NewProblem$DefComp$MCR == iMass]))) {
        stop("Must specify either H as an input component, or pH with a",
             "corresponding H defined component.")
      }
    } else {
      if (("H" %in% NewProblem$Comp$Name[NewProblem$Comp$MCR == iMass]) ||
          ("pH" %in% NewProblem$InVar$Type[NewProblem$InVar$MCR == iMass])) {
        stop("pH/[H+] specified for non-water mass compartment.")
      }
    }
  }


  # read species information including stoichiometry, log Ks, etc.
  SkipRows = SkipRows + NDefComp + 4
  Tmp = scan(file = ParamFile, skip = SkipRows, sep = "\n", nlines = NSpec,
             what = "character", quiet = TRUE)
  TmpSplit = strsplit(Tmp, ",")
  for (i in 1:NSpec) {
    SpecNC_i = as.integer(trimws(TmpSplit[[i]][4]))
    NewProblem = AddSpecies(
      ThisProblem = NewProblem,
      SpecName = as.character(trimws(TmpSplit[[i]][1])),
      SpecMCName = as.character(trimws(TmpSplit[[i]][2])),
      SpecActCorr = as.character(trimws(TmpSplit[[i]][3])),
      SpecCompNames = list(trimws(as.character(TmpSplit[[i]][4 + seq(1, 2 * SpecNC_i, by = 2)]))),
      SpecCompStoichs = list(as.integer(trimws(TmpSplit[[i]][4 + seq(2, 2 * SpecNC_i, by = 2)]))),
      SpecLogK = as.numeric(trimws(TmpSplit[[i]][5 + SpecNC_i * 2])),
      SpecDeltaH = as.numeric(trimws(TmpSplit[[i]][6 + SpecNC_i * 2])),
      SpecTempKelvin = as.numeric(trimws(TmpSplit[[i]][7 + SpecNC_i * 2])),
      DoCheck = FALSE
    )
  }

  # Check that OH is either a component, defined component, or species
  for (iMass in 1:NMass) {
    if (grepl("Water", NewProblem$Mass$Name[iMass], ignore.case = TRUE)) {
      if (sum("OH" %in% c(NewProblem$InCompName, NewProblem$InSpecName,
                          NewProblem$InDefCompName)) != 1L) {
        stop("Must specify OH as either an input component, a defined",
             "component, or a formation reaction species.")
      }
    }
  }

  # -Get Phase information
  SkipRows = SkipRows + NSpec + 3
  if (NPhase > 0) {
    # read Phase information including stoichiometry, log Ks, etc.
    Tmp = scan(file = ParamFile, skip = SkipRows, sep = "\n", nlines = NPhase,
               what = "character", quiet = TRUE)
    TmpSplit = strsplit(Tmp, ",")
    for (i in 1:NPhase) {
      PhaseNC_i = as.integer(trimws(TmpSplit[[i]][2]))
      NewProblem = AddPhases(
        ThisProblem = NewProblem,
        PhaseName = as.character(trimws(TmpSplit[[i]][1])),
        PhaseCompNames = list(trimws(as.character(TmpSplit[[i]][2 + seq(1, 2 * PhaseNC_i, by = 2)]))),
        PhaseCompStoichs = list(as.integer(trimws(TmpSplit[[i]][2 + seq(2, 2 * PhaseNC_i, by = 2)]))),
        PhaseLogK = as.numeric(trimws(TmpSplit[[i]][3 + PhaseNC_i * 2])),
        PhaseDeltaH = as.numeric(trimws(TmpSplit[[i]][4 + PhaseNC_i * 2])),
        PhaseTempKelvin = as.numeric(trimws(TmpSplit[[i]][5 + PhaseNC_i * 2])),
        PhaseMoles = as.numeric(trimws(TmpSplit[[i]][6 + PhaseNC_i * 2])),
        DoCheck = FALSE
      )
    }
  }

  # -get Special definitions
  SkipRows = SkipRows + NPhase + 2
  if (NSpecialDef > 0) {
    Tmp = read.csv(file = ParamFile, header = TRUE, skip = SkipRows,
                   nrows = NSpecialDef, strip.white = TRUE)
    NewProblem = AddSpecialDefs(
      ThisProblem = NewProblem,
      Value = trimws(Tmp[, 2]),
      SpecialDef = trimws(Tmp[, 1]),
      DoCheck = FALSE
    )
  }

  # -get critical accumulation information
  # --> this part also needs to happen in ListCAT function
  SkipRows = SkipRows + NSpecialDef + 3
  if (NCAT > 0) {
    CATab = read.csv(file = ParamFile, header = TRUE, skip = SkipRows,
                     nrows = NCAT, strip.white = TRUE)
    colnames(CATab) = c("Num", "CA", "Species", "Test.Type", "Duration",
                        "Lifestage", "Endpoint", "Quantifier", "References",
                        "Miscellaneous")
    NewProblem = AddCriticalValues(ThisProblem = NewProblem, CATab = CATab, DoCheck = FALSE)

  }

  CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)

  if (WriteLog) {
    CHESSLog(ThisProblem = NewProblem)
  }

  return(NewProblem)
}
