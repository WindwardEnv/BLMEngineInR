#' Run the Biotic Ligand Model
#'
#' `BLM` will run the Windward Environmental Biotic Ligand Model (BLM) with the
#' provided parameter file, input file, and options.
#'
#' @param ParamFile The path and file name of the parameter file
#' @param InputFile The path and file name of the chemistry input file
#' @param ThisProblem (optional) A problem list object, such as returned by
#'   `DefineProblem`.
#' @param AllInput (optional) An input chemistry list object, such as returned
#'   by `GetData`.
#' @param DoTox Should this be a speciation (TRUE) or toxicity (FALSE) run? In a
#'   speciation run, the total metal is input and the free metal and metal bound
#'   to the biotic ligand is calculated. In a toxicity run, the critical
#'   accumulation is input and the free and total metal concentrations that
#'   would result in that amount bound to the biotic ligand is calculated.
#' @param iCA (integer) The index of the critical accumulation value in the
#'   parameter file critical accumulation table.
#' @param QuietFlag Either "Quiet", "Very Quiet", or "Debug". With "Very Quiet",
#'   the simulation will run siliently. With "Quiet", the simulation will print
#'   "Obs=1", "Obs=2", etc... to the console. With "Debug", intermediate
#'   information from the CHESS function will print to the console.
#' @param ConvergenceCriteria (numeric) The maximum allowed CompError in for the
#'   simulation to be considered complete. CompError = abs(CalcTotMoles -
#'   TotMoles) / TotMoles
#' @param MaxIter (integer) The maximum allowed CHESS iterations before the
#'   program should give up.
#' @param DoPartialStepsAlways Should CHESS do strict Newton-Raphson iterations
#'   (FALSE), or try to improve the simulation with partial N-R steps (trying to
#'   prevent oscillations).
#'
#' @return A data frame with chemistry speciation information, including species
#'   concentrations, species activities, and total concentrations.
#'
#' @export
#'
#' @examples
#' ## Not run:
#' ## running the BLM function with a parameter file and input file:
#' # BLM(ParamFile = "path/mypfile.dat", InputFile = "path/myinputfile.blm")
#' #
#' ## running the BLM with parameter and input objects
#' # ThisProblem = DefineProblem(ParamFile = "path/mypfile.dat")
#' # AllInput = list(
#' #   NObs = 5,
#' #   InLabObs = as.matrix(data.frame(
#' #     row.names = "Obs",
#' #     Obs = 1:5,
#' #     ObsNum = as.character(1:5),
#' #     ID = "pH series",
#' #     ID2 = paste0("pH=",5:9)
#' #   )),
#' #   InVarObs = as.matrix(data.frame(
#' #     row.names = "Obs",
#' #     Obs = 1:5,
#' #     Temp = 15,
#' #     pH = 5:9,
#' #     DOC = 0.1,
#' #     HA = 10
#' #   )),
#' #   InCompObs = as.matrix(data.frame(
#' #     row.names = "Obs",
#' #     Obs = 1:5,
#' #     Zn = 1E-7,
#' #     Ca = 0.000299416,
#' #     Mg = 0.000501954,
#' #     Na = 0.00110049,
#' #     K = 5.37108e-05,
#' #     SO4 = 0.000799487,
#' #     Cl = 5.35921e-05,
#' #     CO3 = 0.00109987
#' #   ))
#' # )
#' # FunctionInputs = c(
#' #   ThisProblem[which(names(ThisProblem) %in% formalArgs("MatchInputsToProblem"))],
#' #   AllInput[which(names(AllInput) %in% formalArgs("MatchInputsToProblem"))])
#' # AllInput = c(AllInput, do.call("MatchInputsToProblem", args = FunctionInputs))
#' # ResultsTable_pH = BLM(ThisProblem = ThisProblem, AllInput = AllInput,
#' #                       DoTox = TRUE, iCA = 1)
#' ## End(Not run)
BLM = function(ParamFile = character(),
               InputFile = character(),
               ThisProblem = list(),
               AllInput = list(),
               DoTox = logical(),
               iCA = 1L,
               QuietFlag = c("Quiet", "Very Quiet", "Debug"),
               # writeOutputFile = F, outputFileName = NULL,
               # criticalSource = c("ParamFile","InputFile"),
               ConvergenceCriteria = 0.0001,
               MaxIter = 500L,
               DoPartialStepsAlways = FALSE) {

  StartTime = Sys.time()

  # error catching
  stopifnot(file.exists(ParamFile) || !is.null(ThisProblem))
  stopifnot(file.exists(InputFile) || !is.null(AllInput))
  QuietFlag = match.arg(QuietFlag)

  # 1. parse out parameter file in DefineProblem
  #   --> parameter file name
  #   <-- R variable that defines the problem for immediate use in CHESS
  if (length(ThisProblem) == 0) {
    ThisProblem = DefineProblem(ParamFile)
  }

  # 2. Read InputFile
  #   --> input file name, component info from ParamFile
  #   <-- R variable with component concentrations (total/free dep on ParamFile)
  if (length(AllInput) == 0) {
    FunctionInputs = ThisProblem[
      which(names(ThisProblem) %in% formalArgs("GetData"))]
    FunctionInputs$InputFile = InputFile
    AllInput = do.call("GetData", args = FunctionInputs)
  }

  # Inputs Error catching
  ReferenceProblemList = BlankProblem()
  if (typeof(ThisProblem) != "list") {
    stop("Invalid problem list - it's not a list.")
  }
  if(!all(names(ReferenceProblemList) %in% names(ThisProblem))) {
    print(setdiff(names(ReferenceProblemList), names(ThisProblem)))
    stop("Invalid problem list - missing elements.")
  }
  ThisProblem.types = sapply(ThisProblem, typeof)[match(names(ReferenceProblemList), names(ThisProblem))]
  Reference.types = sapply(ReferenceProblemList, typeof)
  if (!all(ThisProblem.types == Reference.types, na.rm = TRUE)) {
    print(
      data.frame(
        element.name = names(ReferenceProblemList),
        ThisProblem.type = ThisProblem.types,
        Reference.type = Reference.types
      )[ThisProblem.types != Reference.types, ]
    )
    stop("Invalid problem list - incorrect types.")
  }
  ReferenceInputList = list(NObs = 1L,
                            InLabObs = array("", dim = c(1, ThisProblem$NInLab)),
                            InVarObs = array(0.0, dim = c(1, ThisProblem$NInVar)),
                            InCompObs = array(0.0, dim = c(1, ThisProblem$NInComp)),
                            SysTempCelsiusObs = array(0.0, 1),
                            SysTempKelvinObs = array(0.0, 1),
                            pH = array(0.0, 1),
                            TotConcObs = array(0.0, dim = c(1, ThisProblem$NComp)),
                            HumicSubstGramsPerLiterObs = array(0.0, dim = c(1,2)))
  if (typeof(AllInput) != "list") {
    stop("Invalid inputs list - it's not a list.")
  }
  if(!all(names(ReferenceInputList) %in% names(AllInput))) {
    print(setdiff(names(ReferenceInputList), names(AllInput)))
    stop("Invalid inputs list - missing elements.")
  }
  AllInput.types = sapply(AllInput, typeof)[match(names(ReferenceInputList), names(AllInput))]
  Reference.types = sapply(ReferenceInputList, typeof)
  if (!all(AllInput.types == Reference.types, na.rm = TRUE)) {
    print(
      data.frame(
        element.name = names(ReferenceInputList),
        AllInput.type = AllInput.types,
        Reference.type = Reference.types
      )[AllInput.types != Reference.types, ]
    )
    stop("Invalid inputs list - incorrect types.")
  }


  # Save some common variables for initializing arrays
  NObs = AllInput$NObs
  MassName = ThisProblem$MassName
  NComp = ThisProblem$NComp
  CompName = ThisProblem$CompName
  NSpec = ThisProblem$NSpec
  SpecName = ThisProblem$SpecName
  BLComp = ThisProblem$BLComp
  if (DoTox) {
    CATargetDefault = ThisProblem$CATab$CA[iCA] * (10^-6) /
      ThisProblem$CompSiteDens[BLComp]
  } else {
    CATargetDefault = NA
  }

  # Initialize ThisInput as ThisProblem, with one observation's worth of
  # concentrations
  ThisInput = ThisProblem
  ThisInput$InLab = array(character(ThisProblem$NInLab),
                          dimnames = list(ThisProblem$InLabName))
  ThisInput$TotConc = array(numeric(NComp), dimnames = list(CompName))
  ThisInput$CompConc = array(numeric(NComp), dimnames = list(CompName))
  ThisInput$SpecConc = array(numeric(NSpec), dimnames = list(SpecName))
  ThisInput$HumicSubstGramsPerLiter = array(numeric(2), dimnames = list(c("HA", "FA")))

  FunctionInputs = list(
    QuietFlag = QuietFlag,
    ConvergenceCriteria = ConvergenceCriteria,
    MaxIter = MaxIter,
    DoPartialStepsAlways = DoPartialStepsAlways,
    DoTox = DoTox,
    # MetalComp = ThisProblem$MetalComp,
    CATarget = NA
  )
  ObsFunctionInputNames = formalArgs("CHESS")[
    formalArgs("CHESS") %in% names(FunctionInputs) == FALSE]

  # Initialize the output list
  OutputLabels = cbind(
    data.frame(Obs = 1:NObs),
    AllInput$InLabObs
  )
  # InputCols = c("Obs", ThisProblem$InLabName, ThisProblem$InVarName,
  #               paste0("Input.",ThisProblem$InCompName))
  MiscOutputCols = c("FinalIter", "FinalToxIter", "FinalMaxError", "IonicStrength")
  ZCols = paste0("Z_", c("HA","FA"))
  MassAmtCols = paste0(MassName, " (",ThisProblem$MassUnit,")")
  SpecConcCols = paste0(SpecName," (mol/",ThisProblem$MassUnit[ThisProblem$SpecMCR],")")
  TotConcCols = paste0("T.", CompName, " (mol/",ThisProblem$MassUnit[ThisProblem$CompMCR],")")
  SpecMolesCols = paste0(SpecName," (mol)")
  TotMolesCols = paste0("T.", CompName, " (mol)")
  SpecActCols = paste0("Act.", SpecName)
  OutList = list(
    Inputs = cbind(
      OutputLabels,
      AllInput$InVarObs,
      AllInput$InCompObs
    ),
    Miscellaneous = cbind(
      OutputLabels,
      matrix(NA, nrow = NObs, ncol = length(MiscOutputCols) + length(ZCols) + length(MassAmtCols) + 1,
             dimnames = list(NULL, c(MiscOutputCols, ZCols, MassAmtCols, "Status")))
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
    )
  )

  # For tox mode, fill in missing values
  if (DoTox) {
    EmptyOrInvalidMetalObs = is.na(AllInput$TotConcObs[, ThisProblem$MetalName]) |
      (AllInput$TotConcObs[, ThisProblem$MetalName] <= 0)
    AllInput$TotConcObs[EmptyOrInvalidMetalObs, ThisProblem$MetalName] = 1.0E-7
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
      FunctionInputs$CATarget = CATargetDefault * ThisInput$TotConc[BLComp]
    }

    # 3. Run the speciation problem
    #   --> R variable defining problem from step 1 and inputs from step 2
    #   <-- R variable with speciation outputs
    FunctionInputs[ObsFunctionInputNames] = ThisInput[ObsFunctionInputNames]

    if (is.na(ThisInput$SysTempKelvin) |
        any(is.na(ThisInput$TotConc)) |
        (ThisProblem$DoWHAM & any(is.na(ThisInput$HumicSubstGramsPerLiter[!is.na(ThisProblem$WHAMDonnanMCR)])))) {
      # Incomplete chemistry, so skip calling CHESS
      OutList$Miscellaneous$Status[iObs] = "Incomplete Chemistry"
    } else if ((ThisInput$SysTempKelvin <= 263) |
               any(ThisInput$TotConc <= 0) |
               (ThisInput$TotConc["H"] > 1) | (ThisInput$TotConc["H"] < 1.0E-14) |
               (any(ThisInput$HumicSubstGramsPerLiter[!is.na(ThisProblem$WHAMDonnanMCR)] <= 0))) {
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

        if (is.na(OutList$Miscellaneous$FinalMaxError[iObs]) |
            (OutList$Miscellaneous$FinalMaxError[iObs] > ConvergenceCriteria)) {
          OutList$Miscellaneous$Status[iObs] = "Not Converged"
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

  if (!DoTox) { OutList$Miscellaneous$FinalToxIter = NULL }

  print(Sys.time() - StartTime)


  return(OutList)

}
