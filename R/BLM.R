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
  allInput = do.call(getData, args = c(list(inputFile=inputFile), thisProblem[formalArgs(getDAta)]))
  # globalVars = c(thisProblem, thisInput)

  # Save some common variables for inializing arrays
  NComp = thisProblem$NComp
  CompNames = thisProblem$CompNames
  NSpec = thisProblem$NSpec
  SpecNames = thisProblem$SpecNames

  # Initialize the output array
  out = array(numeric(1), dim = c(allInput$NObs, NSpec),
              dimnames = list(1:allInput$NObs, SpecNames))

  # Initialize thisInput as thisProblem, with one observation's worth of
  # concentrations
  thisInput = thisProblem
  thisInput$obsLabels = array(character(2), dimnames = list(c("Site Label", "Sample Label")))
  thisInput$totConcObs = array(numeric(NComp), dimnames = list(CompNames))
  thisInput$CConc = array(numeric(NComp), dimnames = list(CompNames))
  thisInput$SConc = array(numeric(NSpec), dimnames = list(SpecNames))

  # Loop through each observation
  for (iObs in 1:allInput$NObs){
    thisInput$obsLabels = allInput$obsLabels[iObs,]
    thisInput$totConcObs = allInput$totConcObs[iObs,]

    # For now, we're going to use test data, setting the initial "guess" to the
    # actual component free ion concentrations
    if (inputFile == "Test") {
      data("TestDataFreeConc", envir = environment())
      thisInput$CConc = BLMEngineInR::TestDataFreeConc[1:NComp]
    } else if (inputFile == "Full_Inorg"){
      data("Full_InorgDataFreeConc", envir = environment())
      thisInput$CConc = BLMEngineInR::Full_InorgDataFreeConc[1:NComp]
    } else {
      thisInput$CConc = do.call(initialGuess, args = thisInput[formalArgs(initialGuess)])
    }

    # 3. Run the speciation problem
    #   --> R variable defining problem from step 1
    #   --> R variable with inputs from step 2
    #   <-- R variable with speciation outputs
    out[iObs,] = do.call("CppCalcSpecConc", args = thisInput[formalArgs("CppCalcSpecConc")])

  }

  return(out)

}
