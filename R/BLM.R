#' Run the Biotic Ligand Model
#'
#' `BLM` will run the Windward Environmental Biotic Ligand Model (BLM) with the
#' provided parameter file, input file, and options.
#'
#' @param paramFile the path and file name of the parameter file
#' @param inputFile the path and file name of the chemistry input file
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
#' # BLM(paramFile = "path/mypfile.dat", inputFile = "path/myinputfile.blm")
#' ## End(Not run)
BLM = function(paramFile, inputFile#, quiet = T, mode = c("speciation","toxicity"),
               # writeOutputFile = F, outputFileName = NULL,
               # criticalSource = c("paramFile","inputFile"),
               # convergenceCriteria = 0.001
               ){

  # error catching
  # stopifnot(file.exists(paramFile))
  # stopifnot(file.exists(inputFile))
  # mode = match.arg(mode)

  # 1. parse out parameter file in DefineProblem
  #   --> parameter file name
  #   <-- R variable that defines the problem for immediate use in CHESS
  thisProblem = defineProblem(paramFile)

  # 2. Read inputFile
  #   --> input file name
  #   <-- R variable with component concentrations (total or free dep on paramFile)
  allInput = do.call(getData, args = c(list(inputFile=inputFile),
                                       thisProblem[which(names(thisProblem) %in% formalArgs(getData))]))
  # globalVars = c(thisProblem, thisInput)

  # Save some common variables for initializing arrays
  NComp = thisProblem$NComp
  CompName = thisProblem$CompName
  NSpec = thisProblem$NSpec
  SpecName = thisProblem$SpecName

  # Initialize the output array
  out = array(numeric(1), dim = c(allInput$NObs, NSpec),
              dimnames = list(1:allInput$NObs, SpecName))

  # Initialize thisInput as thisProblem, with one observation's worth of
  # concentrations
  thisInput = thisProblem
  thisInput$obsLabels = array(character(thisProblem$NInLab),
                              dimnames = list(thisProblem$InLabName))
  thisInput$TotConcObs = array(numeric(NComp), dimnames = list(CompName))
  thisInput$CConc = array(numeric(NComp), dimnames = list(CompName))
  thisInput$SConc = array(numeric(NSpec), dimnames = list(SpecName))

  # Loop through each observation
  for (iObs in 1:allInput$NObs){
    thisInput$obsLabels = allInput$obsLabels[iObs,]
    thisInput$TotConcObs = allInput$TotConcObs[iObs,]

    # For now, we're going to use test data, setting the initial "guess" to the
    # actual component free ion concentrations
    thisInput$CConc = do.call(initialGuess, args = thisInput[formalArgs(initialGuess)])
    thisInput$LogCConc = log10(thisInput$CConc)

    # 3. Run the speciation problem
    #   --> R variable defining problem from step 1
    #   --> R variable with inputs from step 2
    #   <-- R variable with speciation outputs
    out[iObs,] = 10^do.call("CppCalcLogSpecConc", args = thisInput[formalArgs("CppCalcLogSpecConc")])

  }

  return(out)

}
