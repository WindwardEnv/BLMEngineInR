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

#' @title Convert from a WHAM V thermodynamic database file
#'
#' @description This function will take a thermodynamic database file used by
#'   the Windemere Humic Aqueous Model (WHAM) V and the Windows BLM and convert
#'   it into a BLMEngineInR chemistry problem list.
#'
#' @param ThermoDBSName Character string with the file path of the WHAM
#'   thermodynamic database file to convert (typically ".dbs" file extension).
#' @param RWHAMFile (optional) Character string with the file path of the
#'   R-format WHAM parameter file to save (suggest file extension ".wdat").
#' @param RParamFile (optional) Character string with the file path of the
#'   R-format BLM parameter file to save (suggest file extension ".dat4").
#'
#' @return The BLMEngineInR-compatible chemistry problem object. If RWHAMFile or
#'   RParamFile are provided, this will return invisibly.
#' @export
ConvertWHAMVThermoFile = function(ThermoDBSName,
                                  RWHAMFile = NULL,
                                  RParamFile = NULL) {

  NewWHAM = DefineWHAM(WHAMVer = "V")

  WHAMHSParams = c("nA", "pKA", "pKB", "dpKA", "dpKB", "P", "fprB",
                   "Radius", "MolWt")
  Tmp = utils::read.csv(
    file = ThermoDBSName,
    header = FALSE,
    skip = 1L,
    nrows = 2L,
    col.names = c("Label", WHAMHSParams)
  )
  for (i in WHAMHSParams) {
    NewWHAM[[i]] = as.numeric(Tmp[, i])
    names(NewWHAM[[i]]) = c("HA", "FA")
  }
  NewWHAM$fprT = c(HA = 0.0, FA = 0.0)
  NewWHAM$dLK1A = NewWHAM$dpKA
  NewWHAM$dLK1B = NewWHAM$dpKB

  Tmp = utils::read.csv(file = ThermoDBSName, header = FALSE, skip = 3,
                        nrows = 2)
  NewWHAM$DLF = as.numeric(Tmp[1, 2])
  NewWHAM$KZED = as.numeric(Tmp[2, 2])

  NWHAMRxn = utils::read.csv(file = ThermoDBSName, header = FALSE,
                             skip = 5, nrows = 1)[, 2]
  WHAMRxnTab = utils::read.csv(
    file = ThermoDBSName, header = FALSE, skip = 6,
    nrows = NWHAMRxn,
    col.names = c("Num", "Name", "Charge", "C1", "C2", "C3", "S1", "S2", "S3",
                  "LogK", "DeltaH_kcal.mol", "pKMAHA", "pKMAFA")
  )
  NewWHAM$Notes = scan(file = ThermoDBSName, what = character(), sep = "\n",
                       skip = 7 + NWHAMRxn, quiet = TRUE)

  # eliminate duplicate species names
  if (any(duplicated(WHAMRxnTab$Name))) {

    NeedsRename = which(
      WHAMRxnTab$Name %in% WHAMRxnTab$Name[duplicated(WHAMRxnTab$Name)]
    )

    WHAMRxnSub = WHAMRxnTab[NeedsRename, ]
    WHAMRxnSub$OldName = WHAMRxnSub$Name

    WHAMRxnSub$C1_name = WHAMRxnTab$Name[match(WHAMRxnSub$C1, WHAMRxnTab$Num)]
    WHAMRxnSub$C2_name = WHAMRxnTab$Name[match(WHAMRxnSub$C2, WHAMRxnTab$Num)]
    WHAMRxnSub$C3_name = WHAMRxnTab$Name[match(WHAMRxnSub$C3, WHAMRxnTab$Num)]

    PatternTmp = "^[[:upper:]][[:upper:]]"
    WHAMRxnSub$C1_needs_paren = grepl(PatternTmp, WHAMRxnSub$C1_name)
    WHAMRxnSub$C2_needs_paren = grepl(PatternTmp, WHAMRxnSub$C2_name)
    WHAMRxnSub$C3_needs_paren = grepl(PatternTmp, WHAMRxnSub$C3_name)

    Tmp = rep("", nrow(WHAMRxnSub))
    IBool = (WHAMRxnSub$C1_needs_paren & (WHAMRxnSub$S1 > 1))
    Tmp[IBool] = paste0(Tmp[IBool], "(")
    IBool = (WHAMRxnSub$S1 > 0)
    Tmp[IBool] = paste0(Tmp, WHAMRxnSub$C1_name)[IBool]
    IBool = (WHAMRxnSub$C1_needs_paren & (WHAMRxnSub$S1 > 1))
    Tmp[IBool] = paste0(Tmp[IBool], ")")
    IBool = (WHAMRxnSub$S1 > 1)
    Tmp[IBool] = paste0(Tmp, WHAMRxnSub$S1)[IBool]

    IBool = (WHAMRxnSub$C2_needs_paren & (WHAMRxnSub$S2 > 1))
    Tmp[IBool] = paste0(Tmp[IBool], "(")
    IBool = (WHAMRxnSub$S2 > 0)
    Tmp[IBool] = paste0(Tmp, WHAMRxnSub$C2_name)[IBool]
    IBool = (WHAMRxnSub$C2_needs_paren & (WHAMRxnSub$S2 > 1))
    Tmp[IBool] = paste0(Tmp[IBool], ")")
    IBool = (WHAMRxnSub$S2 > 1)
    Tmp[IBool] = paste0(Tmp, WHAMRxnSub$S2)[IBool]

    IBool = (WHAMRxnSub$C3_needs_paren & (WHAMRxnSub$S3 > 1))
    Tmp[IBool] = paste0(Tmp[IBool], "(")
    IBool = (WHAMRxnSub$S3 > 0)
    Tmp[IBool] = paste0(Tmp, WHAMRxnSub$C3_name)[IBool]
    IBool = (WHAMRxnSub$C3_needs_paren & (WHAMRxnSub$S3 > 1))
    Tmp[IBool] = paste0(Tmp[IBool], ")")
    IBool = (WHAMRxnSub$S3 > 1)
    Tmp[IBool] = paste0(Tmp, WHAMRxnSub$S3)[IBool]

    WHAMRxnSub$Name = Tmp

    WHAMRxnTab$Name[match(WHAMRxnSub$Num, WHAMRxnTab$Num)] = WHAMRxnSub$Name
  }

  # Separate out the components and species
  WHAMCompTab = WHAMRxnTab[(WHAMRxnTab$Num < 100),
                           c("Num", "Name", "Charge", "pKMAHA", "pKMAFA")]
  WHAMSpecTab = WHAMRxnTab[(WHAMRxnTab$Num > 100), ]

  # Species reaction enthalpy changes need converting
  WHAMSpecTab$SpecTempKelvin = 1 / 0.003354
  WHAMSpecTab$DeltaH = WHAMSpecTab$DeltaH_kcal.mol * 220 * 8.314 * log(10)
  # Note: This is a conversion used by WHAM V. The WHAM deltaH values are in
  # units of kcal/mol while the Van't Hoff equation (which the BLM uses for
  # thermodynamic corrections) needs J/mol. Strictly speaking, you should
  # multiply by 4184 to convert from kcal to J, but Tipping was trying to do
  # this is log10-space and he wanted a simple equation. The Van't Hoff equation
  # is K2 = K1 * exp(deltaH / R * (1/T1 - 1/T2)). Tipping's version is log10K2 =
  # log10K1 + 220 * deltaH(kcal/mol) * (0.003354 - 1 / T2). Basically, the 220
  # is doing three things: 1) the conversion from kcal to J, 2) the conversion
  # from natural log to log10, and 3) dividing by the gas constant R = 8.314 J /
  # (mol * Kelvin). Doing these three actions with the official conversion: 4184
  # J/kcal / (8.314 J/(molK) * ln(10)) = 218.5576 mol*K / kcal. Tipping
  # apparently decided he only needed 2 significant digits, so he rounded it off
  # to 220. This makes his conversion from kcal to J, effectively, 4211.612
  # J/kcal.

  WHAMRxnTab$dLK2 = 0.13 # dLK2 is 0.13 for WHAM V, varies for WHAM VI/VII
  NewWHAM$MetalsTable = WHAMRxnTab[
    (WHAMRxnTab$pKMAHA != 999) & (WHAMRxnTab$pKMAFA != 999),
    c("Name", "pKMAHA", "pKMAFA", "dLK2")
  ]
  colnames(NewWHAM$MetalsTable)[1] = "Metal"
  NewWHAM$MetalsTable = NewWHAM$MetalsTable[order(NewWHAM$MetalsTable$Metal), ]
  rownames(NewWHAM$MetalsTable) = NULL

  # Now add the inorganic reactions to a BLM problem list
  NotHOHComps = (WHAMCompTab$Name %in% c("H", "OH") == FALSE)
  NewProblem = AddInComps(
    ThisProblem = BLMEngineInR::water_problem,
    InCompName = WHAMCompTab$Name[NotHOHComps],
    InCompCharge = WHAMCompTab$Charge[NotHOHComps],
    InCompType = "MassBal",
    InCompActCorr = "Debye",
    InCompMCName = "Water"
  )
  NewProblem = AddSpecies(
    ThisProblem = NewProblem,
    SpecName = WHAMSpecTab$Name,
    SpecCompNames = apply(
      WHAMSpecTab[, c("C1", "C2", "C3")],
      MARGIN = 1,
      FUN = function(X) {
        WHAMCompTab$Name[stats::na.omit(match(X, WHAMCompTab$Num))]
      }
    ),
    SpecCompStoichs = apply(
      WHAMSpecTab[, c("S1", "S2", "S3")],
      MARGIN = 1,
      FUN = function(X) { X[X != 0L] }
    ),
    SpecMCName = "Water",
    SpecActCorr = "Debye",
    SpecLogK = WHAMSpecTab$LogK,
    SpecDeltaH = WHAMSpecTab$DeltaH,
    SpecTempKelvin = WHAMSpecTab$SpecTempKelvin,
    DoCheck = TRUE
  )

  NewWHAM$File = NA_character_
  NewWHAM$Ver = "V"
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
