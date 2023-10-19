#' Run the Biotic Ligand Model
#'
#' `BLM` will run the Windward Environmental Biotic Ligand Model (BLM) with the
#' provided parameter file, input file, and options.
#'
#' @param paramFile the path and file name of the parameter file
#' @param inputFile the path and file name of the chemistry input file
#' @param quiet logical. If `TRUE`, iteration information will be displayed in
#'   the console.
#' @param mode the mode to run the model in. Only values of `"speciation"` or
#'   `"toxicity"` are supported, or partial matches to those character strings.
#' @param writeOutputFile,outputFileName,criticalSource,convergenceCriteria
#'   Other parameters that are not implemented, but expected to be needed.
#'
#' @return A data frame with chemistry speciation information, including total
#'   concentrations.
#' @export
#'
#' @examples
#' ## Not run:
#' BLM(paramFile = "path/mypfile.dat", inputFile = "path/myinputfile.blm")
#' ## End(Not run)
BLM = function(paramFile, inputFile, quiet = T, mode = c("speciation","toxicity"),
               writeOutputFile = F, outputFileName = NULL,
               criticalSource = c("paramFile","inputFile"),
               convergenceCriteria = 0.001){

  # error catching
  # stopifnot(file.exists(paramFile))
  # stopifnot(file.exists(inputFile))
  mode = match.arg(mode)

  # 1. parse out parameter file in DefineProblem
  #   --> parameter file name
  #   <-- R variable that defines the problem for immediate use in CHESS
  thisProblem = defineProblem(paramFile)

  # 2. Read inputFile
  #   --> input file name
  #   <-- R variable with component concentrations (total or free dep on paramFile)
  thisInput = getData(inputFile, thisProblem$NComp)

  # 3. Run the speciation problem
  #   --> R variable defining problem from step 1
  #   --> R variable with inputs from step 2
  #   <-- R variable with speciation outputs
  out = do.call(CppCalcSpecConc, args = c(thisProblem, thisInput))

  return(out)

}
