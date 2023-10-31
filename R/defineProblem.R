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
  NSpec = temp[5,1]
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
  SpecName = character(NSpec)
  SpecMC = integer(NSpec)
  SpecType = integer(NSpec)
  SpecActCorr = integer(NSpec)
  SpecNC = integer(NSpec)                  # the number of components that form species(i)
  CompList = matrix(data = 0, nrow = NSpec, ncol = NComp + NDefComp)  # the list of components (by component number) that form species(i) for SpecNC(i) number of components
  Stoich = matrix(data = 0, nrow = NSpec, ncol = NComp + NDefComp)          # Stoich(i,j) = the amount of component(j) needed to form species(i)
  LogK = numeric(NSpec)
  DeltaH = numeric(NSpec)

  # read species information including stoichiometry, log Ks, etc.
  temp=scan(file=paramFile,skip=23+NMass+NInVar+NComp+NDefComp,sep="\n",nlines=NSpec,what="character")
  temp.split = strsplit(temp,",")
  for (i in 1:NSpec) {
    SpecName[i] = as.character(trimws(temp.split[[i]][1]))
    SpecMC[i] = as.integer(match(trimws(temp.split[[i]][2]), MassName))
    SpecType[i] = as.integer(trimws(temp.split[[i]][3]))
    SpecActCorr[i] = as.integer(trimws(temp.split[[i]][4]))
    SpecNC[i] = as.integer(trimws(temp.split[[i]][5]))
    for (j in 1:SpecNC[i]) {
      CompList[i, j] = match(trimws(temp.split[[i]][6 + (j - 1) * 2]), CompName)
      Stoich[i, CompList[i, j]] = as.integer(trimws(temp.split[[i]][7 + (j - 1) * 2]))
    }
    LogK[i] = as.numeric(trimws(temp.split[[i]][6 + SpecNC[i] * 2]))
    DeltaH[i] = as.numeric(trimws(temp.split[[i]][7 + SpecNC[i] * 2]))
  }

  # Add components to the species list
  SpecName = c(CompName, SpecName)
  SpecMC = c(CompMC, SpecMC)
  SpecType = c(CompType, SpecType)
  SpecActCorr = c(CompActCorr, SpecActCorr)
  SpecNC = c(array(1L, NComp + NDefComp), SpecNC)
  temp = matrix(0, nrow = NComp + NDefComp, ncol = NComp + NDefComp)
  temp[, 1] = 1:(NComp + NDefComp)
  CompList = rbind(temp, CompList)
  Stoich = rbind(diag(1, nrow = NComp + NDefComp, ncol = NComp + NDefComp), Stoich)
  LogK = c(rep(0,NComp + NDefComp), LogK)
  DeltaH = c(rep(0, NComp + NDefComp), DeltaH)

  # Trim down CompList...give it 2 extra columns in case any WHAM species need it
  CompList = CompList[,1:(max(which(apply(CompList, MARGIN = 2, FUN = sum)>0)) + 2)]



  # # read species list, tableau, and properties
  # temp=read.delim(file=paramFile,header=FALSE,sep=",", skip=23+NMass+NInVar+NComp+NDefComp,nrows=NSpec)
  # SpecName=as.character(trimws(temp[,1]))
  # SpecMC=as.integer(match(trimws(strsplit(temp,",")[[i]][2]),MassName))
  # SpecType= as.integer(temp[,3])
  # SpecActCorr= as.integer(temp[,4])
  # SpecNC= as.numeric(temp[,5])
  # CompList = matrix(nrow=NSpec,ncol=NComp)
  # for (i in 1: NSpec) {
  #   for (j in 1: SpecNC[i]){
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

  wdatFile = NULL
  WHAMVer = "V"
  # to do: Where is wdatFile coming from? Should it be passed within the
  # parameter file? Probably, since the WHAM calibration will definitely
  # change the results of any simulation.

  # # ---Make WHAM species from DOC
  # # -get critical accumulation information --> this part also needs to happen in listCAT function
  # # --number of critical accumulations in table
  # # --critical accumulation table (CAT)

  # assemble output
  out = list(
    NMass = NMass,
    NInVar = NInVar,
    NComp = NComp,
    NDefComp = NDefComp,
    NSpec = NSpec,
    NPhases = NPhases,

    InVarName = InVarName,
    InVarFlag = InVarFlag,

    CompName = CompName,
    CompCharge = CompCharge,
    CompMC = CompMC,
    CompType = CompType,
    CompActCorr = CompActCorr,
    CompSiteDens = CompSiteDens,

    DefCompName = DefCompName,
    DefCompFrom = DefCompFrom,
    DefCompCharge = DefCompCharge,
    DefCompMC = DefCompMC,
    DefCompType = DefCompType,
    DefCompActCorr = DefCompActCorr,
    DefCompSiteDens = DefCompSiteDens,

    SpecName = SpecName,
    SpecMC = SpecMC,
    SpecType = SpecType,
    SpecActCorr = SpecActCorr,
    SpecNC = SpecNC,
    CompList = CompList,
    Stoich = Stoich,
    LogK = LogK,
    DeltaH = DeltaH
  )

  # Expand WHAM components and species, if needed
  if (any((CompType == 6) | (CompType == "WHAM"))){
    out2 = do.call("expandWHAM",
                   args = c(out[which(names(out) %in% formalArgs("expandWHAM"))],
                            list(wdatFile=wdatFile, WHAMVer = WHAMVer)))
    out[names(out2)] = out2
  }

  return(out)
}
