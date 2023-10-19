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
#'  \item{`K`}{the equilibrium coefficients (numeric vector of length `NSpec`)}
#'  \item{`logK`}{the log10-transformed equilibrium coefficients (numeric vector of length `NSpec`)}
#'  \item{`Stoich`}{the stoichiometry matrix of the reactions (integer matrix of size `NSpec` x `nComp`)}
#'  \item{`CConc`}{the component free ion concentrations (numeric vector of length `nComp`)}
#' }
#'
#' @noRd
defineProblem = function(paramFile){

  # # error catching
  # stopifnot(file.exists(paramFile))

  # for now, this function will be returning our test data
  if (paramFile == "Test") {
    data("TestDataFreeConc", "TestDataK", "TestDataStoich")
    NComp = ncol(TestDataStoich)
    NSpec = nrow(TestDataStoich)
    CConc = TestDataFreeConc[1:NComp]
    K = TestDataK
    logK = log10(TestDataK)
    Stoich = TestDataStoich
  } else if (paramFile == "Full_Inorg"){
    data("Full_InorgDataFreeConc", "Full_InorgDataK", "Full_InorgDataStoich")
    NComp = ncol(Full_InorgDataStoich)
    NSpec = nrow(Full_InorgDataStoich)
    CConc = Full_InorgDataFreeConc[1:NComp]
    K = Full_InorgDataK
    logK = log10(Full_InorgDataK)
    Stoich = Full_InorgDataStoich
  }

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
    K = K,
    logK = logK,
    Stoich = Stoich,
    CConc = CConc
  )
  return(out)
}
