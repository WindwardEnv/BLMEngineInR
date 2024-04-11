#' Run the Biotic Ligand Model
#'
#' `BLM` will run the Windward Environmental Biotic Ligand Model (BLM) with the
#' provided parameter file, input file, and options.
#'
#' @param ParamFile The path and file name of the parameter file
#' @param InputFile The path and file name of the chemistry input file
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
#' # BLM(ParamFile = "path/mypfile.dat", InputFile = "path/myinputfile.blm")
#' ## End(Not run)
BLM = function(ParamFile = character(),
               InputFile = character(),
               DoTox = logical(),
               iCA = 1L,
               QuietFlag = c("Quiet", "Very Quiet", "Debug"),
               # writeOutputFile = F, outputFileName = NULL,
               # criticalSource = c("ParamFile","InputFile"),
               ConvergenceCriteria = 0.001,
               MaxIter = 30L,
               DoPartialStepsAlways = FALSE) {

  StartTime = Sys.time()

  # error catching
  stopifnot(file.exists(ParamFile))
  stopifnot(file.exists(InputFile))
  QuietFlag = match.arg(QuietFlag)

  # 1. parse out parameter file in DefineProblem
  #   --> parameter file name
  #   <-- R variable that defines the problem for immediate use in CHESS
  ThisProblem = DefineProblem(ParamFile)

  # 2. Read InputFile
  #   --> input file name, component info from ParamFile
  #   <-- R variable with component concentrations (total/free dep on ParamFile)
  FunctionInputs = ThisProblem[
    which(names(ThisProblem) %in% formalArgs("GetData"))]
  FunctionInputs$InputFile = InputFile
  AllInput = do.call("GetData", args = FunctionInputs)

  # Save some common variables for initializing arrays
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

  # Initialize the output array
  MiscOutputCols = c("FinalIter", "FinalMaxError", "IonicStrength")
  SpecConcCols = paste0(SpecName," (mol/",ThisProblem$MassUnit[ThisProblem$SpecMCR],")")
  SpecMolesCols = paste0(SpecName," (mol)")
  SpecActCols = paste0("Act.", SpecName)
  TotConcCols = paste0("T.", CompName, " (mol/",ThisProblem$MassUnit[ThisProblem$CompMCR],")")
  TotMolesCols = paste0("T.", CompName, " (mol)")
  ZCols = paste0("Z_", c("HA","FA"))
  MassAmtCols = paste0(MassName, " (",ThisProblem$MassUnit,")")
  Out = data.frame(Obs = 1:AllInput$NObs)
  Out = cbind(Out, AllInput$InLabObs)
  Out = cbind(Out, AllInput$InVarObs)
  Out = cbind(Out, AllInput$InCompObs)
  InputCols = c("Obs", ThisProblem$InLabName, ThisProblem$InVarName,
                paste0("Input.",ThisProblem$InCompName))
  colnames(Out) = InputCols
  Out[, MiscOutputCols] = NA
  Out[, SpecConcCols] = NA
  Out[, SpecActCols] = NA
  Out[, SpecMolesCols] = NA
  Out[, TotConcCols] = NA
  Out[, TotMolesCols] = NA
  Out[, ZCols] = NA
  Out[, MassAmtCols] = NA
  # Out = array(
  #   NA,
  #   dim = c(AllInput$NObs, NSpec + NComp + length(MiscOutputCols)),
  #   dimnames = list(1:AllInput$NObs, c(MiscOutputCols, SpecName, TotConcCols))
  # )

  # Loop through each observation
  for (iObs in 1:AllInput$NObs) {

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

    Tmp = do.call("CHESS", args = FunctionInputs)
    Out[iObs, MiscOutputCols] = unlist(Tmp[MiscOutputCols])
    Out[iObs, SpecConcCols] = Tmp$SpecConc
    Out[iObs, SpecActCols] = Tmp$SpecAct
    Out[iObs, SpecMolesCols] = Tmp$SpecMoles
    Out[iObs, TotConcCols] = Tmp$CalcTotConc
    Out[iObs, TotMolesCols] = Tmp$CalcTotMoles
    Out[iObs, ZCols] = Tmp$WHAMSpecCharge
    Out[iObs, MassAmtCols] = Tmp$MassAmt

  }

  # Make summary columns for organically-bound components
  if (ThisProblem$DoWHAM) {
    for (iComp in ThisProblem$InCompName){
      OrgCols = SpecMolesCols[grepl(iComp, SpecMolesCols) &
                                (grepl("DOC", SpecMolesCols) |
                                   grepl("Donnan", SpecMolesCols))]
      # OrgCols = colnames(Out)[
      #   grepl(iComp, colnames(Out)) & grepl("[(]mol[)]", colnames(Out)) &
      #     (grepl("DOC", colnames(Out)) | grepl("Donnan", colnames(Out)))
      # ]
      Out[,paste0("TOrg.",iComp," (mol/L)")] = rowSums(Out[, OrgCols])
    }
  }

  print(Sys.time() - StartTime)

  return(Out)

}
