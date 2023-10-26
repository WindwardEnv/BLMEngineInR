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
#'  \item{\code{NComp}}{integer; the number of components}
#'  \item{\code{NSpec}}{integer; the number of species}
#'  \item{\code{CompNames}}{character vector of length `NComp`; component names}
#'  \item{\code{CompCharge}}{integer vector of length `NComp`; the charge of the components as free ions}
#'  \item{\code{CompType}}{integer vector of length `NComp`; the type of component. It should be a fixed set of values where 1 =...}
#'  \item{\code{CompActCorr}}{integer vector of length `NComp`; the method to use for activity corrections with this component, where 1 == }
#'  \item{\code{CompSiteDens}}{numeric vector of length `NComp`; the density of binding sites for this component - usually 1, except for DOC and BL components}
#'  \item{\code{SpecNames}}{character vector of length `NSpec`; species names}
#'  \item{\code{SpecType}}{integer vector of length `NSpec`; the type of chemical species, where 1 = ...}
#'  \item{\code{SpecActCorr}}{integer vector of length `NSpec`; the method to use for activity corrections with this speies where 1 = ...}
#'  \item{\code{Stoich}}{the stoichiometry matrix of the reactions (integer matrix of size `NSpec` x `nComp`)}
#'  \item{\code{K}}{the equilibrium coefficients (numeric vector of length `NSpec`)}
#'  \item{\code{logK}}{the log10-transformed equilibrium coefficients (numeric vector of length `NSpec`)}
#'  \item{\code{deltaH}}{the enthalpy change for each formation reaction (species)}
#' }
#'
#' @keywords internal
#'
#' @noRd
defineProblem = function(paramFile){

  # error catching
  stopifnot(file.exists(paramFile))

  # for now, this function will be returning our test data
  if (paramFile == "Test") {
    data("TestDataK", "TestDataStoich", envir = environment())
    NComp = ncol(BLMEngineInR::TestDataStoich)
    NSpec = nrow(BLMEngineInR::TestDataStoich)
    logK = log10(BLMEngineInR::TestDataK)
    Stoich = BLMEngineInR::TestDataStoich
  } else if (paramFile == "Full_Inorg"){
    data("Full_InorgDataK", "Full_InorgDataStoich", envir = environment())
    NComp = ncol(BLMEngineInR::Full_InorgDataStoich)
    NSpec = nrow(BLMEngineInR::Full_InorgDataStoich)
    logK = log10(BLMEngineInR::Full_InorgDataK)
    Stoich = BLMEngineInR::Full_InorgDataStoich
  }

  # read in parameter file
  # -get number of components
  # -get number of species
  # -get component information
  # --CompNames: name (character string)
  CompNames = colnames(Stoich)
  # --CompCharge: charge (signed integer)
  CompCharge = integer(NComp)
  # --CompType: component type
  CompType = integer(NComp)
  # --CompActCorr: activity correction type
  CompActCorr = integer(NComp)
  # --CompSiteDens: site density
  CompSiteDens = numeric(NComp)
  # -get species information
  # --SpecNames: name
  SpecNames = names(logK)
  # --SpecType: species type
  SpecType = integer(NSpec)
  # --SpecActCorr: activity correction type
  SpecActCorr = integer(NSpec)
  # --Stoich: stoichiometry
  # --logK: equilibrium coefficient
  # --detlaH: enthalpy change
  deltaH = numeric(NSpec)
  # -get special species/component information
  # --name of metal
  # --name of biotic ligand
  # --name of BL-metal complex(es)
  # --name of DOC component(s)
  # ---Make WHAM species from DOC
  # -get critical accumulation information --> this part also needs to happen in listCAT function
  # --number of critical accumulations in table
  # --critical accumulation table (CAT)

  # assemble output
  out = list(
    NComp = NComp,
    NSpec = NSpec,
    CompNames = CompNames,
    CompCharge = CompCharge,
    CompType = CompType,
    CompActCorr = CompActCorr,
    CompSiteDens = CompSiteDens,
    SpecNames = SpecNames,
    SpecType = SpecType,
    SpecActCorr = SpecActCorr,
    Stoich = Stoich,
    K = 10^logK,
    logK = logK,
    deltaH = deltaH
  )
  return(out)
}
