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
#'  \item{\code{CompName}}{character vector of length `NComp`; component names}
#'  \item{\code{CompCharge}}{integer vector of length `NComp`; the charge of the components as free ions}
#'  \item{\code{CompType}}{integer vector of length `NComp`; the type of component. It should be a fixed set of values where 1 =...}
#'  \item{\code{CompActCorr}}{integer vector of length `NComp`; the method to use for activity corrections with this component, where 1 == }
#'  \item{\code{CompSiteDens}}{numeric vector of length `NComp`; the density of binding sites for this component - usually 1, except for DOC and BL components}
#'  \item{\code{SpecName}}{character vector of length `NSpec`; species names}
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

  # read in parameter file
  # -get number of mass compartments
  # -get number of input variables
  # -get number of input components
  # -get number of defined components
  # -get number of species
  # -get number of phases
  # -get component information

  # read the dimensions of the various elements of the reaction list
  temp=read.delim(file=paramFile,header=FALSE,sep=",", skip=2,nrows=6)
  NMass    = temp[1,1]
  NInVar   = temp[2,1]
  NComp    = temp[3,1]
  NDefComp = temp[4,1]
  NSpecies = temp[5,1]
  NPhases  = temp[6,1]

#  MassName = character(NMass)
#  MassAmt = numeric(NMass)
#  MassUnit = character(NMass)
#  InVarName = character(NInVar)
#  InVarFlag = logical(NInVar)
#  CompName = character(NComp)
#  CompCharge = integer(NComp)
#  CompMC    = integer(NComp)
#  CompType  = integer(NComp)
#  CompActCorr = integer(NComp)
#  CompSiteDens = numeric(NComp)

  # read compartment list
  temp=read.delim(file=paramFile,header=FALSE,sep=",", skip=11,nrows=NMass)
  MassName = as.character(trimws(temp[,1]))
  MassAmt = as.numeric(temp[,2])
  MassUnit = as.character(trimws(temp[,3]))

  # read input variables
  temp=read.delim(file=paramFile,header=FALSE,sep=",", skip=14+NMass,nrows=NInVar)
  InVarName = as.character(trimws(temp[,1]))
  InVarFlag = as.logical(trimws(toupper(temp[,2])))

  # read component list and properties
  temp=read.delim(file=paramFile,header=FALSE,sep=",", skip=17+NMass+NInVar,nrows=NComp)
  CompName=as.character(trimws(temp[,1]))
  CompCharge= as.numeric(temp[,2])
  CompMC=match(trimws(temp[,3]),MassName)
  CompType= as.integer(temp[,4])
  CompActCorr= as.integer(temp[,5])
  CompSiteDens= as.numeric(temp[,6])  # do we still need this?


  # read defined component list and properties
  temp=read.delim(file=paramFile,header=FALSE,sep=",", skip=20+NMass+NInVar+NComp,nrows=NDefComp)
  DefCompName=as.character(trimws(temp[,1]))
  DefCompFrom= as.numeric(temp[,2])                 # we will eventually want to be able to read numbers or strings
  DefCompCharge= as.numeric(temp[,3])
  DefCompMC=match(trimws(temp[,4]),MassName)
  DefCompType= as.integer(temp[,5])
  DefCompActCorr= as.integer(temp[,6])
  DefCompSiteDens= as.numeric(temp[,7])  # do we still need this?

  #
  # concatenate DefComp and Comp
  CompName=c(CompName,DefCompName)
  CompCharge=c(CompCharge,DefCompCharge)
  CompMC=c(CompMC,DefCompMC)
  CompType=c(CompType,DefCompType)
  CompActCorr=c(CompActCorr,DefCompActCorr)
  CompSiteDens=c(CompSiteDens,DefCompSiteDens)
  #NComp=NComp+NDefComp   # don't increase the NComp dimension so that we can
                          # keep track of which portion of the component list
                          # will show up in the interface (NComp) and which
                          # does not (NDefComp)
                          # BUT . . . that means we need to consider the
                          # dimensions of both component lists in subsequent calcs

  # Create variables for species information
  SpeciesName=character(NSpecies)
  SpeciesMC=integer(NSpecies)
  SpeciesType=integer(NSpecies)
  SpeciesActCorr=integer(NSpecies)
  SpeciesNC=integer(NSpecies)                  # the number of components that form species(i)
  CompList = matrix(nrow=NSpecies,ncol=NComp+NDefComp)  # the list of components (by component number) that form species(i) for SpeciesNC(i) number of components
  SC=matrix(nrow=NSpecies,ncol=NComp+NDefComp)          # SC(i,j) = the amount of component(j) needed to form species(i)
  LogK=numeric(NSpecies)
  DeltaH=numeric(NSpecies)

  # read species information including stoichiometry, log Ks, etc.
  temp=scan(file=paramFile,skip=23+NMass+NInVar+NComp+NDefComp,sep="\n",nlines=NSpecies,what="character")
  for (i in 1: NSpecies) {
    SpeciesName[i]=as.character(trimws(strsplit(temp,",")[[i]][1]))
    SpeciesMC[i]=as.integer(match(trimws(strsplit(temp,",")[[i]][2]),MassName))
    SpeciesType[i]= as.integer(trimws(strsplit(temp,",")[[i]][3]))
    SpeciesActCorr[i]= as.integer(trimws(strsplit(temp,",")[[i]][4]))
    SpeciesNC[i]= as.integer(trimws(strsplit(temp,",")[[i]][5]))
    for (j in 1: SpeciesNC[i]){
      CompList[i,j]=match(trimws(strsplit(temp,",")[[i]][6+(j-1)*2]),CompName)
      SC[i,CompList[i,j]]=as.integer(trimws(strsplit(temp,",")[[i]][7+(j-1)*2]))
    }
    LogK[i]=as.numeric(trimws(strsplit(temp,",")[[i]][6+(SpeciesNC[i]-1)*2]))
    DeltaH[i]=as.numeric(trimws(strsplit(temp,",")[[i]][7+(SpeciesNC[i]-1)*2]))
  }


  # # read species list, tableau, and properties
  # temp=read.delim(file=paramFile,header=FALSE,sep=",", skip=23+NMass+NInVar+NComp+NDefComp,nrows=NSpecies)
  # SpeciesName=as.character(trimws(temp[,1]))
  # SpeciesMC=as.integer(match(trimws(strsplit(temp,",")[[i]][2]),MassName))
  # SpeciesType= as.integer(temp[,3])
  # SpeciesActCorr= as.integer(temp[,4])
  # SpeciesNC= as.numeric(temp[,5])
  # CompList = matrix(nrow=NSpecies,ncol=NComp)
  # for (i in 1: NSpecies) {
  #   for (j in 1: SpeciesNC[i]){
  #     CompList[i,j]=match(trimws(temp[i,5+j]),CompName)
  #   }
  # }

# CompNames = colnames(Stoich)
# # --CompCharge: charge (signed integer)
# CompCharge = integer(NComp)
# # --CompType: component type
# CompType = integer(NComp)
# # --CompActCorr: activity correction type
# CompActCorr = integer(NComp)
# # --CompSiteDens: site density
# CompSiteDens = numeric(NComp)

  # # -get species information
  # # --SpecNames: name
  #  SpecNames = names(logK)
  # # --SpecType: species type
  # SpecType = integer(NSpec)
  # # --SpecActCorr: activity correction type
  # SpecActCorr = integer(NSpec)
  # # --Stoich: stoichiometry
  # # --logK: equilibrium coefficient
  # # --detlaH: enthalpy change
  # deltaH = numeric(NSpec)
  # # -get special species/component information
  # # --name of metal
  # # --name of biotic ligand
  # # --name of BL-metal complex(es)
  # # --name of DOC component(s)
  # # ---Make WHAM species from DOC
  # # -get critical accumulation information --> this part also needs to happen in listCAT function
  # # --number of critical accumulations in table
  # # --critical accumulation table (CAT)

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
