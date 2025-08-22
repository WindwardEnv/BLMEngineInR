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
#' mypfile = system.file("extdata", "ParameterFiles",
#'                       "carbonate_system_only.dat4",
#'                       package = "BLMEngineInR", mustWork = TRUE)
#' thisProblem = DefineProblem(mypfile)
#'
DefineProblem = function(ParamFile, WriteLog = FALSE) {

  # error catching
  ParamFile = normalizePath(ParamFile, mustWork = TRUE)

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
  # - Misc = anything else you want available, for instance a variable you want
  #   to use for a defined component

  # read component list and properties
  SkipRows = SkipRows + NInVar + 3
  Tmp = read.csv(file = ParamFile, header = TRUE, skip = SkipRows,
                 nrows = NInComp, strip.white = TRUE)
  if (("OH" %in% trimws(Tmp[, 1])) && ("OH" %in% NewProblem$DefComp$Name)) {
    # if the "OH" defined component was added with pH, but is specified as an
    # input component, remove the defined component first.
    NewProblem = RemoveDefComps(
      ThisProblem = NewProblem,
      DefCompToRemove = intersect(intersect(c("OH"),
                                            NewProblem$DefComp$Name),
                                  trimws(Tmp[, 1]))
    )
  }
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
      # update the component and defined component specifications for "H", and
      # "OH, which were likely added with the input variables.
      IComp = match(c("H", "OH"), NewProblem$Comp$Name)
      INewDefComp = match(c("H", "OH"), DefCompName)
      IOldDefComp = match(c("H", "OH"), NewProblem$DefComp$Name)
      IBool = !is.na(INewDefComp) & !is.na(IOldDefComp) & !is.na(IComp)
      IComp = IComp[IBool]
      INewDefComp = INewDefComp[IBool]
      IOldDefComp = IOldDefComp[IBool]
      NewProblem$Comp$Charge[IComp] = as.integer(trimws(Tmp[INewDefComp, 3]))
      NewProblem$Comp$MCName[IComp] = as.character(trimws(Tmp[INewDefComp, 4]))
      NewProblem$Comp$MCR[IComp] =
        match(NewProblem$Comp$MCName[IComp], NewProblem$Mass$Name)
      NewProblem$Comp$Type[IComp] = as.character(trimws(Tmp[INewDefComp, 5]))
      NewProblem$Comp$ActCorr[IComp] = as.character(trimws(Tmp[INewDefComp, 6]))
      NewProblem$Comp$SiteDens[IComp] = as.numeric(trimws(Tmp[INewDefComp, 7]))
      NewProblem$DefComp$FromNum[IOldDefComp] = DefCompFromNum[INewDefComp]
      NewProblem$DefComp$FromVar[IOldDefComp] = DefCompFromVar[INewDefComp]
      NewProblem$DefComp$Charge[IOldDefComp] = NewProblem$Comp$Charge[IComp]
      NewProblem$DefComp$MCName[IOldDefComp] = NewProblem$Comp$MCName[IComp]
      NewProblem$DefComp$MCR[IOldDefComp] = NewProblem$Comp$MCR[IComp]
      NewProblem$DefComp$Type[IOldDefComp] = NewProblem$Comp$Type[IComp]
      NewProblem$DefComp$ActCorr[IOldDefComp] = NewProblem$Comp$ActCorr[IComp]
      NewProblem$DefComp$SiteDens[IOldDefComp] = NewProblem$Comp$SiteDens[IComp]
      Tmp = Tmp[-INewDefComp, ]
      DefCompName = DefCompName[-INewDefComp]
      DefCompFromNum = DefCompFromNum[-INewDefComp]
      DefCompFromVar = DefCompFromVar[-INewDefComp]
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

  # Check that pH is in the component list in the water compartment
  WaterMC = which(grepl("Water", NewProblem$Mass$Name, ignore.case = TRUE))
  if (length(WaterMC) == 0L) {
    stop("No water mass compartment.")
  }
  InCompH = ("H" %in% NewProblem$InCompName) &
    ("H" %in% NewProblem$Comp$Name[NewProblem$MCR == WaterMC])
  DefCompH =
    ("pH" %in% NewProblem$InVar$Type[NewProblem$InVar$MCR == WaterMC]) &
    ("H" %in% NewProblem$DefComp$Name[NewProblem$DefComp$MCR == WaterMC])
  if (!xor(InCompH, DefCompH)) {
    stop("Must specify either H as an input component, or pH with a",
         "corresponding H defined component in the water mass compartment.")
  }
  NonWaterMC = which(!grepl("Water", NewProblem$Mass$Name, ignore.case = TRUE))
  InCompH = ("H" %in% NewProblem$Comp$Name[NewProblem$Comp$MCR == NonWaterMC])
  DefCompH =
    ("pH" %in% NewProblem$InVar$Type[NewProblem$InVar$MCR == NonWaterMC])
  if (InCompH || DefCompH) {
    stop("pH or H specified for non-water mass compartment.")
  }

  # read species information including stoichiometry, log Ks, etc.
  SkipRows = SkipRows + NDefComp + 4
  Tmp = scan(file = ParamFile, skip = SkipRows, sep = "\n", nlines = NSpec,
             what = "character", quiet = TRUE)
  TmpSplit = strsplit(Tmp, ",")
  for (i in 1:NSpec) {

    if (("OH" %in% trimws(TmpSplit[[i]][1])) &&
          ("OH" %in% NewProblem$DefComp$Name)) {
      # if the "OH" defined component was added with pH, but is specified as a
      # species, remove the defined component first.
      NewProblem = RemoveDefComps(
        ThisProblem = NewProblem,
        DefCompToRemove = intersect(intersect(c("OH"),
                                              NewProblem$DefComp$Name),
                                    trimws(TmpSplit[[i]][1]))
      )
    }

    SpecNCTimes2 = as.integer(trimws(TmpSplit[[i]][4])) * 2
    NameInd = 4 + seq(1, SpecNCTimes2, by = 2)
    StoichInd = 4 + seq(2, SpecNCTimes2, by = 2)
    NewProblem = AddSpecies(
      ThisProblem = NewProblem,
      SpecName = as.character(trimws(TmpSplit[[i]][1])),
      SpecMCName = as.character(trimws(TmpSplit[[i]][2])),
      SpecActCorr = as.character(trimws(TmpSplit[[i]][3])),
      SpecCompNames = list(trimws(as.character(TmpSplit[[i]][NameInd]))),
      SpecCompStoichs = list(as.integer(trimws(TmpSplit[[i]][StoichInd]))),
      SpecLogK = as.numeric(trimws(TmpSplit[[i]][5 + SpecNCTimes2])),
      SpecDeltaH = as.numeric(trimws(TmpSplit[[i]][6 + SpecNCTimes2])),
      SpecTempKelvin = as.numeric(trimws(TmpSplit[[i]][7 + SpecNCTimes2])),
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
      PhaseNCTimes2 = as.integer(trimws(TmpSplit[[i]][2])) * 2
      NameInd = 2 + seq(1, PhaseNCTimes2, by = 2)
      StoichInd = 2 + seq(2, PhaseNCTimes2, by = 2)
      NewProblem = AddPhases(
        ThisProblem = NewProblem,
        PhaseName = trimws(TmpSplit[[i]][1]),
        PhaseCompNames = list(trimws(TmpSplit[[i]][NameInd])),
        PhaseCompStoichs = list(as.integer(trimws(TmpSplit[[i]][StoichInd]))),
        PhaseLogK = as.numeric(trimws(TmpSplit[[i]][3 + PhaseNCTimes2])),
        PhaseDeltaH = as.numeric(trimws(TmpSplit[[i]][4 + PhaseNCTimes2])),
        PhaseTempKelvin = as.numeric(trimws(TmpSplit[[i]][5 + PhaseNCTimes2])),
        PhaseMoles = as.numeric(trimws(TmpSplit[[i]][6 + PhaseNCTimes2])),
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
    colnames(CATab) = c("Num", "CA", "Species", "TestType", "Duration",
                        "Lifestage", "Endpoint", "Quantifier", "References",
                        "Miscellaneous")
    NewProblem = AddCriticalValues(
      ThisProblem = NewProblem,
      CATab = CATab,
      DoCheck = FALSE
    )

  }

  CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)

  if (WriteLog) {
    CHESSLog(ThisProblem = NewProblem)
  }

  return(NewProblem)
}
