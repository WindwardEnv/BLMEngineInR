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
               QuietFlag = c("Quiet","Very Quiet","Debug"),
               # writeOutputFile = F, outputFileName = NULL,
               # criticalSource = c("ParamFile","InputFile"),
               ConvergenceCriteria = 0.001,
               MaxIter = 30L
               ){

  # error catching
  # stopifnot(file.exists(ParamFile))
  # stopifnot(file.exists(InputFile))
  # mode = match.arg(mode)

  # 1. parse out parameter file in DefineProblem
  #   --> parameter file name
  #   <-- R variable that defines the problem for immediate use in CHESS
  ThisProblem = DefineProblem(ParamFile)

  # 2. Read InputFile
  #   --> input file name
  #   <-- R variable with component concentrations (total or free dep on ParamFile)
  AllInput = do.call(GetData, args = c(list(InputFile=InputFile),
                                       ThisProblem[which(names(ThisProblem) %in% formalArgs(GetData))]))
  # globalVars = c(ThisProblem, ThisInput)

  # Save some common variables for initializing arrays
  NComp = ThisProblem$NComp
  CompName = ThisProblem$CompName
  NSpec = ThisProblem$NSpec
  SpecName = ThisProblem$SpecName

  # Initialize the output array
  out = array(numeric(1), dim = c(AllInput$NObs, NSpec),
              dimnames = list(1:AllInput$NObs, SpecName))

  # Initialize ThisInput as ThisProblem, with one observation's worth of
  # concentrations
  ThisInput = ThisProblem
  ThisInput$InLab = array(character(ThisProblem$NInLab),
                          dimnames = list(ThisProblem$InLabName))
  ThisInput$TotConc = array(numeric(NComp), dimnames = list(CompName))
  ThisInput$CompConc = array(numeric(NComp), dimnames = list(CompName))
  ThisInput$SpecConc = array(numeric(NSpec), dimnames = list(SpecName))

  # Loop through each observation
  for (iObs in 1:AllInput$NObs){

    if (QuietFlag != "Very Quiet"){ print(paste0("Obs=",iObs)) }

    ThisInput$InLab = AllInput$InLabObs[iObs,]
    ThisInput$TotConc = AllInput$TotConcObs[iObs,]

    # For now, we're going to use test data, setting the initial "guess" to the
    # actual component free ion concentrations
    ThisInput$CompConc = do.call(InitialGuess, args = ThisInput[formalArgs(InitialGuess)])
    ThisInput$LogCompConc = log10(ThisInput$CompConc)

    # 3. Run the speciation problem
    #   --> R variable defining problem from step 1 and inputs from step 2
    #   <-- R variable with speciation outputs
    out[iObs,] = do.call("RCHESS", args = ThisInput[formalArgs("RCHESS")])

  }

  return(out)

}
