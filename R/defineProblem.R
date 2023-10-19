#' Define the speciation problem
#'
#' `defineProblem` reads in a parameter file, and sets up the required vectors
#' and matrices that will be needed to run the speciation calculations in CHESS.
#'
#' @param paramFile the path and file name to a parameter file
#'
#' @returns Returns whatever is needed for the speciation problem to be run. A
#' `list` object with the following components:
#' \describe{
#'  \item{`NComp`}{integer; the number of components}
#'  \item{`NSpec`}{integer; the number of species}
#'  \item{`logK`}{the log10-transformed equilibrium coefficients (numeric vector of length `NSpec`)}
#'  \item{`Stoich`}{the stoichiometry matrix of the reactions (integer matrix of size `NSpec` x `nComp`)}
#' }
#'
#' @noRd
defineProblem = function(paramFile){

  # # error catching
  # stopifnot(file.exists(paramFile))

  # for now, this function will be returning our test data
  data(TestDataFreeConc, TestDataK, TestDataStoich)
  NComp = length(TestDataFreeConc)
  NSpec = length(TestDataK)
  CConc = TestDataFreeConc[1:NComp]
  logK = log10(TestDataK)
  Stoich = TestDataStoich

  # read in parameter file
  # -get number of components
  # -get number of species
  # -get component information
  # --name (character string)
  # --charge (signed integer)
  # --component type
  # --activity correction type
  # --site density
  # -get species information
  # --name
  # --species type
  # --activity correction type
  # --stoichiometry
  # --equilibrium coefficient (logK)
  # --enthalpy change (deltaH)
  # -get special species/component information
  # --name of metal
  # --name of biotic ligand
  # --name of BL-metal complex(es)
  # --name of DOC component(s)
  # -get critical accumulation information --> this part also needs to happen in listCAT function
  # --number of critical accumulations in table
  # --critical accumulation table (CAT)
  # -Make WHAM species from DOC

  # assemble output
  out = list(
    NComp = NComp,
    NSpec = NSpec,
    logK = logK,
    Stoich = Stoich,
    CConc = CConc
  )
  return(out)
}
