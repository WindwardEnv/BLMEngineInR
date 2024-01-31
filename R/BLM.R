#' Run the Biotic Ligand Model
#'
#' `BLM` will run the Windward Environmental Biotic Ligand Model (BLM) with the
#' provided parameter file, input file, and options.
#'
#' @param ParamFile the path and file name of the parameter file
#' @param InputFile the path and file name of the chemistry input file
# @param quiet logical. If `TRUE`, iteration information will be displayed in
#   the console.
# @param mode the mode to run the model in. Only values of `"speciation"` or
#   `"toxicity"` are supported, or partial matches to those character strings.
# @param writeOutputFile,outputFileName,criticalSource,convergenceCriteria
#   Other parameters that are not implemented, but expected to be needed.
#'
#' @return A data frame with chemistry speciation information, including total
#'   concentrations.
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
               MaxIter = 30L) {

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
  ThisInput$SolHS = array(numeric(2), dimnames = list(c("HA", "FA")))

  FunctionInputs = list(
    QuietFlag = QuietFlag,
    ConvergenceCriteria = ConvergenceCriteria,
    MaxIter = MaxIter,
    DoTox = DoTox,
    # MetalComp = ThisProblem$MetalComp,
    CATarget = NA
  )
  ObsFunctionInputNames = formalArgs("CHESS")[
    formalArgs("CHESS") %in% names(FunctionInputs) == FALSE]

  # Initialize the output array
  MiscOutputCols = c("FinalIter", "FinalMaxError")
  SpecActCols = paste0("Act.", SpecName)
  TotConcCols = paste0("T.", CompName)
  Out = data.frame(Obs = 1:AllInput$NObs)
  Out = cbind(Out, AllInput$InLabObs)
  Out = cbind(Out, AllInput$InVarObs)
  Out[, MiscOutputCols] = NA
  Out[, SpecName] = NA
  Out[, SpecActCols] = NA
  Out[, TotConcCols] = NA
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
    ThisInput$SolHS = AllInput$SolHSObs[iObs, c("HA", "FA")]

    if (DoTox) {
      FunctionInputs$CATarget = CATargetDefault * ThisInput$TotConc[BLComp]
    }

    # 3. Run the speciation problem
    #   --> R variable defining problem from step 1 and inputs from step 2
    #   <-- R variable with speciation outputs
    FunctionInputs[ObsFunctionInputNames] = ThisInput[ObsFunctionInputNames]

    Tmp = do.call("CHESS", args = FunctionInputs)
    Out[iObs, MiscOutputCols] = unlist(Tmp[MiscOutputCols])
    Out[iObs, SpecName] = Tmp$SpecConc[SpecName]
    Out[iObs, SpecActCols] = Tmp$SpecAct[1:NSpec]
    Out[iObs, TotConcCols] = Tmp$CalcTotConc[CompName]

  }

  return(Out)

}
