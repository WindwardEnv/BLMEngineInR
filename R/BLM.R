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

#' Run the Biotic Ligand Model
#'
#' `BLM` will run the Windward Environmental Biotic Ligand Model (BLM) with the
#' provided parameter file, input file, and options.
#'
#' @param ParamFile (optional) The path and file name of the parameter file
#' @param InputFile (optional) The path and file name of the chemistry input
#'   file
#' @param ThisProblem (optional) A problem list object, such as returned by
#'   `DefineProblem`.
#' @param AllInput (optional) An input chemistry list object, such as returned
#'   by `GetData`.
#' @param DoTox Should this be a speciation (TRUE) or toxicity (FALSE) run? In a
#'   speciation run, the total metal is input and the free metal and metal bound
#'   to the biotic ligand is calculated. In a toxicity run, the critical
#'   accumulation is input and the free and total metal concentrations that
#'   would result in that amount bound to the biotic ligand is calculated.
#' @param iCA (unnecessary unless DoTox = TRUE) Either the index of the critical
#'   accumulation value in the parameter file critical accumulation table, or
#'   the critical accumulation to use in nmol/gw. If this is a single value,
#'   then it will be applied to all observations. If it is a vector with the
#'   same length as the inputs, then each value given will be used for the
#'   corresponding observation.
#' @param QuietFlag Either "Quiet", "Very Quiet", or "Debug". With "Very Quiet",
#'   the simulation will run silently. With "Quiet", the simulation will print
#'   "Obs=1", "Obs=2", etc... to the console. With "Debug", intermediate
#'   information from the CHESS function will print to the console.
#' @param ConvergenceCriteria (numeric) The maximum allowed CompError in for the
#'   simulation to be considered complete. CompError = abs(CalcTotMoles -
#'   TotMoles) / TotMoles
#' @param MaxIter (integer) The maximum allowed CHESS iterations before the
#'   program should give up.
#'
#' @return A data frame with chemistry speciation information, including species
#'   concentrations, species activities, and total concentrations.
#'
#' @export
#'
#' @example tests/examples/examples-BLM.R
BLM = function(ParamFile = character(),
               InputFile = character(),
               ThisProblem = list(),
               AllInput = list(),
               DoTox = logical(),
               iCA = 1L,
               QuietFlag = c("Very Quiet", "Quiet", "Debug"),
               ConvergenceCriteria = 0.0001,
               MaxIter = 100L) {

  DodVidCj = TRUE
  DodVidCjDonnan = FALSE
  DodKidCj = FALSE
  DoGammai = TRUE
  DoJacDonnan = FALSE
  DoJacWHAM = TRUE
  DoWHAMSimpleAdjust = TRUE
  DoDonnanSimpleAdjust = TRUE

  StartTime = Sys.time()

  # error catching
  if ((length(ParamFile) == 0) && (length(ThisProblem) == 0)) {
    stop("Supply either ParamFile or ThisProblem arguments.")
  }
  if ((length(InputFile) == 0) && (length(AllInput) == 0)) {
    stop("Supply either InputFile or AllInput arguments.")
  }
  QuietFlag = match.arg(QuietFlag, c("Very Quiet", "Quiet", "Debug"))

  # 1. parse out parameter file in DefineProblem
  #   --> parameter file name
  #   <-- R variable that defines the problem for immediate use in CHESS
  if (length(ThisProblem) == 0) {
    if (!file.exists(ParamFile)) {
      stop("ParamFile ", ParamFile, " does not exist.")
    }
    ThisProblem = DefineProblem(ParamFile = ParamFile)
  }
  CheckBLMObject(Object = ThisProblem,
                 Reference = BlankProblem(),
                 BreakOnError = TRUE)

  # 2. Read InputFile
  #   --> input file name, component info from ParamFile
  #   <-- R variable with component concentrations (total/free dep on ParamFile)
  if (length(AllInput) == 0) {
    if (!file.exists(InputFile)) {
      stop("InputFile ", InputFile, " does not exist.")
    }
    AllInput = GetData(InputFile = InputFile, ThisProblem = ThisProblem)
  }
  CheckBLMObject(Object = AllInput,
                 Reference = BlankInputList(ThisProblem),
                 BreakOnError = TRUE)

  # Save some common variables for initializing arrays
  NObs = AllInput$NObs
  MassName = ThisProblem$Mass$Name
  NComp = ThisProblem$N["Comp"]
  CompName = ThisProblem$Comp$Name
  NSpec = ThisProblem$N["Spec"]
  SpecName = ThisProblem$Spec$Name
  BLComp = ThisProblem$BL$CompR
  if (DoTox) {

    if (any(ThisProblem$N[c("BL","Metal","BLMetal")]) < 1){
      stop("Cannot do tox mode without at least one each of BL, Metal, and BLMetal defined.")
    }

    if (is.integer(iCA) || all(as.integer(iCA) == iCA)) {
      CATargetDefault = ThisProblem$CATab$CA[iCA]
    } else if (is.double(iCA)) {
      CATargetDefault = iCA
    } else {
      stop("Unknown type of critical value iCA.")
    }
    if (length(iCA) == 1) {
      CATargetDefault = rep(CATargetDefault, NObs)
    } else if (length(iCA) != NObs) {
      stop("Critical value iCA should be either 1 or NObs in length.")
    }
    CATargetDefault = CATargetDefault * (10^-6) /
      ThisProblem$Comp$SiteDens[BLComp]
  } else {
    CATargetDefault = rep(NA, NObs)
  }

  # Initialize ThisInput as ThisProblem, with one observation's worth of
  # concentrations
  ThisInput = ConvertToList(ThisProblemDF = ThisProblem)
  ThisInput$InLab = array(character(ThisProblem$N["InLab"]),
                          dimnames = list(ThisProblem$InLabName))
  ThisInput$TotConc = array(numeric(NComp), dimnames = list(CompName))
  ThisInput$CompConc = array(numeric(NComp), dimnames = list(CompName))
  ThisInput$SpecConc = array(numeric(NSpec), dimnames = list(SpecName))
  ThisInput$HumicSubstGramsPerLiter =
    array(numeric(2), dimnames = list(c("HA", "FA")))

  FunctionInputs = list(
    QuietFlag = QuietFlag,
    ConvergenceCriteria = ConvergenceCriteria,
    MaxIter = MaxIter,
    DoTox = DoTox,
    CATarget = NA,
    DodVidCj = DodVidCj,
    DodVidCjDonnan = DodVidCjDonnan,
    DodKidCj = DodKidCj,
    DoGammai = DoGammai,
    DoJacDonnan = DoJacDonnan,
    DoJacWHAM = DoJacWHAM,
    DoWHAMSimpleAdjust = DoWHAMSimpleAdjust,
    DoDonnanSimpleAdjust = DoDonnanSimpleAdjust
  )
  ObsFunctionInputNames = formalArgs("CHESS")[
    formalArgs("CHESS") %in% names(FunctionInputs) == FALSE]

  # Initialize the output list
  OutputLabels = cbind(
    data.frame(Obs = 1:NObs),
    AllInput$InLabObs
  )
  MiscOutputCols = c("FinalIter", "FinalMaxError", "IonicStrength",
                     "WHAMIonicStrength", "ChargeBalance")
  ZCols = paste0("Z_", c("HA","FA"))
  if (ThisProblem$N["BLMetal"] > 0) {
    BLMetalCols = paste0(ThisProblem$BLMetal$Name, " (nmol/gw)")
    TotBLMetalCol = paste0("Total ", ThisProblem$BL$Name, "-",
                           ThisProblem$Metal$Name, " (nmol/gw)")
  } else {
    BLMetalCols = character()
    TotBLMetalCol = character()
  }
  MassAmtCols = paste0(MassName, " (",ThisProblem$Mass$Unit,")")
  SpecConcCols = paste0(SpecName," (mol/",ThisProblem$Mass$Unit[ThisProblem$Spec$MCR],")")
  TotConcCols = paste0("T.", CompName, " (mol/",ThisProblem$Mass$Unit[ThisProblem$Comp$MCR],")")
  SpecMolesCols = paste0(SpecName," (mol)")
  TotMolesCols = paste0("T.", CompName, " (mol)")
  SpecActCols = paste0("Act.", SpecName)
  SpecActCoefCols = paste0("Gamma.", SpecName)
  OutList = list(
    Inputs = cbind(
      OutputLabels,
      AllInput$InVarObs,
      AllInput$InCompObs
    ),
    Miscellaneous = cbind(
      OutputLabels,
      matrix(NA, nrow = NObs,
             ncol = length(MiscOutputCols) + length(ZCols) +
               length(MassAmtCols) + length(BLMetalCols) + 2,
             dimnames = list(NULL, c(MiscOutputCols, ZCols, MassAmtCols,
                                     BLMetalCols, "Status", "Message")))
    ),
    Concentrations = cbind(
      OutputLabels,
      matrix(NA, nrow = NObs, ncol = NSpec + NComp,
             dimnames = list(NULL, c(SpecConcCols, TotConcCols)))
    ),
    Moles = cbind(
      OutputLabels,
      matrix(NA, nrow = NObs, ncol = NSpec + NComp,
             dimnames = list(NULL, c(SpecMolesCols, TotMolesCols)))
    ),
    Activities = cbind(
      OutputLabels,
      matrix(NA, nrow = NObs, ncol = NSpec,
             dimnames = list(NULL, c(SpecActCols)))
    ),
    ActivityCoefficients = cbind(
      OutputLabels,
      matrix(NA, nrow = NObs, ncol = NSpec,
             dimnames = list(NULL, c(SpecActCoefCols)))
    )
  )
  if (DoTox) {
    OutList$Inputs = cbind(
      OutList$Inputs,
      data.frame(`Target Critical Accumulation (nmol/gw)` =
                   CATargetDefault / (10^-6) *
                   ThisProblem$Comp$SiteDens[BLComp])
    )
  }

  # For tox mode, fill in missing values and add WQC tab
  DoStandards = FALSE
  if (DoTox & (length(iCA) == 1)) {
    if (ThisProblem$CAT$Endpoint[iCA] %in% c("FAV","FCV","HC5","WQS","CMC","CCC")) {
      DoStandards = TRUE

      # If DIV is missing in Duration column, we'll assume the predicted metal
      # concentration should be used as-is. If the ACR is missing, we will only
      # report one toxic unit.
      DIV = NA
      ACR = NA
      NStandardsCols = 1
      NTUCols = 1
      StandardsCols = paste0(
        c(ThisProblem$CAT$Endpoint[iCA], ThisProblem$Metal$Name, "TU"),
        " (", c("\U00B5g/L", "\U00B5g/L", "unitless"), ")"
      )
      WQSCols = StandardsCols[1]
      TUCols = StandardsCols[3]
      if (grepl("^DIV=[[:digit:]]+", ThisProblem$CAT$Duration[iCA])) {
        NStandardsCols = NStandardsCols + 1
        DIV = as.numeric(gsub("^DIV=","", ThisProblem$CAT$Duration[iCA]))
        StandardsCols = c(
          StandardsCols[1],
          paste0(ThisProblem$CAT$Endpoint[iCA],"/DIV (\U00B5g/L)"),
          StandardsCols[2:3]
        )
        if (ThisProblem$CAT$Endpoint[iCA] == "FAV") {
          StandardsCols[2] = paste0("CMC=", StandardsCols[2])
        } else if (ThisProblem$CAT$Endpoint[iCA] == "FCV") {
          StandardsCols[2] = paste0("CCC=", StandardsCols[2])
        }
        WQSCols = StandardsCols[2]
      }

      if (grepl("^ACR=[[:digit:]]+", ThisProblem$CAT$Lifestage[iCA])) {
        ACR = as.numeric(gsub("^ACR=","", ThisProblem$CAT$Lifestage[iCA]))
        StandardsCols = c(
          StandardsCols[1:NStandardsCols],
          paste0(ThisProblem$CAT$Endpoint[iCA],"/ACR (\U00B5g/L)"),
          StandardsCols[NStandardsCols + 1],
          paste("Acute", utils::tail(StandardsCols, 1)),
          "Chronic TU (unitless)"
        )
        NTUCols = NTUCols + 1
        NStandardsCols = NStandardsCols + 1
        if (ThisProblem$CAT$Endpoint[iCA] == "FAV") {
          StandardsCols[NStandardsCols] = paste0("CCC=", StandardsCols[NStandardsCols])
        }
        WQSCols = c(WQSCols, StandardsCols[NStandardsCols])
        TUCols = utils::tail(StandardsCols, 2)
      }

      OutList$Standards = cbind(
        OutputLabels,
        matrix(NA, nrow = NObs, ncol = length(StandardsCols),
               dimnames = list(NULL, StandardsCols))
      )
    }
    EmptyOrInvalidMetalObs = is.na(AllInput$TotConcObs[, ThisProblem$Metal$Name]) |
      (AllInput$TotConcObs[, ThisProblem$Metal$Name] <= 0)
    AllInput$TotConcObs[EmptyOrInvalidMetalObs, ThisProblem$Metal$Name] = 1.0E-7
  }

  # Loop through each observation
  for (iObs in 1:NObs) {

    if (QuietFlag != "Very Quiet") {
      print(paste0("Obs=", iObs))
    }

    ThisInput$InLab = AllInput$InLabObs[iObs, ]
    ThisInput$SysTempKelvin = AllInput$SysTempKelvinObs[iObs]
    ThisInput$TotConc = AllInput$TotConcObs[iObs, CompName]
    ThisInput$HumicSubstGramsPerLiter = AllInput$HumicSubstGramsPerLiterObs[iObs, c("HA", "FA")]

    if (DoTox) {
      FunctionInputs$CATarget = CATargetDefault[iObs] * ThisInput$TotConc[BLComp]
    }

    # 3. Run the speciation problem
    #   --> R variable defining problem from step 1 and inputs from step 2
    #   <-- R variable with speciation outputs
    FunctionInputs[ObsFunctionInputNames] = ThisInput[ObsFunctionInputNames]

    if (is.na(ThisInput$SysTempKelvin) |
        any(is.na(ThisInput$TotConc)) |
        (ThisProblem$DoWHAM & any(is.na(ThisInput$HumicSubstGramsPerLiter[ThisProblem$Index$WHAMDonnanMCR > 0])))) {
      # Incomplete chemistry, so skip calling CHESS
      OutList$Miscellaneous$Status[iObs] = "Incomplete Chemistry"
    } else if ((ThisInput$SysTempKelvin <= 263) |
               any(ThisInput$TotConc <= 0) |
               (ThisInput$TotConc["H"] > 1) | (ThisInput$TotConc["H"] < 1.0E-14) |
               (any(ThisInput$HumicSubstGramsPerLiter[ThisProblem$Index$WHAMDonnanMCR > 0] <= 0))) {
      # Invalid chemistry, so skip calling CHESS
      OutList$Miscellaneous$Status[iObs] = "Invalid Chemistry"
    } else {
      # Complete and valid chemistry, so run CHESS
      RunCompleted = FALSE
      tryCatch({
        Tmp = do.call("CHESS", args = FunctionInputs)
        OutList$Miscellaneous[iObs, MiscOutputCols] = unlist(Tmp[MiscOutputCols])
        OutList$Miscellaneous[iObs, ZCols] = Tmp$WHAMSpecCharge
        OutList$Miscellaneous[iObs, MassAmtCols] = Tmp$MassAmt
        OutList$Concentrations[iObs, SpecConcCols] = Tmp$SpecConc
        OutList$Concentrations[iObs, TotConcCols] = Tmp$CalcTotConc
        OutList$Moles[iObs, SpecMolesCols] = Tmp$SpecMoles
        OutList$Moles[iObs, TotMolesCols] = Tmp$CalcTotMoles
        OutList$Activities[iObs, SpecActCols] = Tmp$SpecAct
        OutList$ActivityCoefficients[iObs, SpecActCoefCols] = Tmp$SpecActivityCoef

        if (is.na(OutList$Miscellaneous$FinalMaxError[iObs]) |
            (OutList$Miscellaneous$FinalMaxError[iObs] > ConvergenceCriteria)) {
          OutList$Miscellaneous$Status[iObs] = "Not Converged"
          OutList$Miscellaneous$Message[iObs] = Tmp$StatusMessage
        } else {
          OutList$Miscellaneous$Status[iObs] = "Okay"
        }
        RunCompleted = TRUE
      }, finally = {
        if (!RunCompleted) {
          OutList$Miscellaneous$Status[iObs] = "Run Error"
        }
      })

    }
  }

  # Make summary columns for organically-bound components
  if (ThisProblem$DoWHAM) {
    for (iComp in ThisProblem$InCompName){
      OrgCols = SpecMolesCols[grepl(iComp, SpecMolesCols) &
                                (grepl("DOC", SpecMolesCols) |
                                   grepl("Donnan", SpecMolesCols))]
      OutList$Concentrations[, paste0("TOrg.",iComp," (mol/L)")] =
        rowSums(OutList$Moles[, OrgCols, drop = FALSE])
    }
  }

  # Make summary columns for BL-Metal components
  if (ThisProblem$N["BLMetal"] > 0) {
    BLMetalCPPCols = paste0(ThisProblem$BLMetal$Name,
                            " (mol/",
                            ThisProblem$Mass$Unit[ThisProblem$Index$BioticLigMCR],
                            ")")
    OutList$Miscellaneous[, BLMetalCols] =
      OutList$Concentrations[, BLMetalCPPCols] *
      ThisProblem$Comp$SiteDens[BLComp] / (10 ^ -6) /
      AllInput$TotConcObs[, BLComp]
    OutList$Miscellaneous[, TotBLMetalCol] =
      rowSums(OutList$Miscellaneous[, BLMetalCols, drop = FALSE])
  }

  if (DoStandards) {
    PredMetalOutputCol = paste0(
      "T.", ThisProblem$Metal$Name, " (mol/",
      ThisProblem$Mass$Unit[ThisProblem$Comp$MCR[ThisProblem$Metal$CompR]], ")"
    )
    MetalMW = BLMEngineInR::MW[[ThisProblem$Metal$Name]]

    # Predicted Metal Column
    OutList$Standards[, StandardsCols[1]] = signif(
      x = OutList$Concentrations[, PredMetalOutputCol] * MetalMW * 10 ^ 6,
      digits = 4
    )

    # Input Metal Column
    OutList$Standards[, StandardsCols[NStandardsCols + 1]] =
      signif(OutList$Inputs[, ThisProblem$Metal$Name] * MetalMW * 10^6,
             digits = 4)

    # WQS column with DIV
    if (!is.na(DIV)) {
      OutList$Standards[, WQSCols[1]] =
        signif(OutList$Standards[, StandardsCols[1]] / DIV, digits = 4)
    }

    # Chronic WQS column by acute-to-chronic ratio
    if (!is.na(ACR)) {
      OutList$Standards[, WQSCols[2]] =
        signif(OutList$Standards[, StandardsCols[1]] / ACR, digits = 4)
    }

    # Calculate Toxic Units
    OutList$Standards[!EmptyOrInvalidMetalObs, TUCols] =
      signif(OutList$Standards[, StandardsCols[NStandardsCols + 1]] /
               OutList$Standards[, WQSCols, drop = FALSE],
             digits = 4)[!EmptyOrInvalidMetalObs, ]
  }

  TimeElapsed = Sys.time() - StartTime
  OutList$TimeElapsed = TimeElapsed
  if (QuietFlag != "Very Quiet") {
    print(TimeElapsed)
  }


  return(OutList)

}
