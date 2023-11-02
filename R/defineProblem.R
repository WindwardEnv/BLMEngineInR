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
#'  \item{\code{CompType}}{character vector of length `NComp`; the type of component. It should be a fixed set of values (MassBal, FixedAct, Substituted, ChargeBal, SurfPot, WHAMV, WHAMVI, WHAMVII, WHAM=FileName.dat)
#'  \item{\code{CompActCorr}}{character vector of length `NComp`; the method to use for activity corrections with this component,  }
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
defineProblem = function(paramFile) {

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
  skipRows = 2
  tmp = read.csv(file = paramFile, header = FALSE, skip = skipRows, nrows = 9)
  NMass = tmp[1, 1]
  NInLab = tmp[2, 1]
  NInVar = tmp[3, 1]
  NInComp = tmp[4, 1]
  NDefComp = tmp[5, 1]
  NSpec = tmp[6, 1]
  NPhases  = tmp[7, 1]
  NBLD  = tmp[8, 1]
  NCAT  = tmp[9, 1]

  # read compartment list
  skipRows = skipRows + 9 + 2
  tmp = read.csv(file = paramFile, header = TRUE, skip = skipRows, nrows = NMass)
  MassName = as.character(trimws(tmp[, 1]))
  MassAmt = as.numeric(tmp[, 2])
  MassUnit = as.character(trimws(tmp[, 3]))

  # read input Labels
  skipRows = skipRows + NMass + 3
  tmp = scan(file = paramFile, what = character(), nlines = NInLab, sep = "\n",
             skip = skipRows, quiet = T)
  InLabName = as.character(trimws(tmp))

  # read input variables
  skipRows = skipRows + NInLab + 2
  tmp = read.csv(file = paramFile, header = TRUE, skip = skipRows, nrows = NInVar)
  InVarName = as.character(trimws(tmp[, 1]))
  InVarMC = match(trimws(tmp[, 2]), MassName)
  InVarType = as.character(trimws(tmp[, 3]))
  stopifnot(all(!is.na(InVarMC)))
  stopifnot(all((InVarType %in% c("Temperature","pH","WHAMV","WHAMVI","WHAMVII","PercHA","PercAFA")) |
                  grepl("WHAM=", InVarType)))
  # - Temperature = the temperature in degrees C
  # - pH = the -log[H]...you know, pH
  # - WHAM<#> = Windemere Humic Aqueous Model organic matter (input mg C/L),
  #   version V, VI, or VII....to be split out by expandWHAM. Assumes 100%
  #   fulvic acid if PercHA not given.
  # - WHAM=<filename> = WHAM component with a custom .wdat file.
  # - PercHA = optionally indicate the percent humic acid in the WHAM component
  #   for that compartment
  # - PercAFA = optionally indicate the percent of active fulvic acid for the
  #   WHAM component for that compartment

  # read component list and properties
  skipRows = skipRows + NInVar + 3
  tmp = read.csv(file = paramFile, header = TRUE, skip = skipRows, nrows = NInComp)
  InCompName = as.character(trimws(tmp[, 1]))
  CompCharge = as.integer(tmp[, 2])
  CompMC = match(trimws(tmp[, 3]), MassName)
  CompType = as.character(trimws(tmp[, 4]))
  CompActCorr = as.character(trimws(tmp[, 5]))

  # read defined component list and properties
  skipRows = skipRows + NInComp + 3
  tmp = read.csv(file = paramFile, header = TRUE, skip = skipRows, nrows = NDefComp)
  DefCompName = as.character(trimws(tmp[, 1]))
  DefCompFrom = as.numeric(tmp[, 2])                 # we will eventually want to be able to read numbers or strings
  DefCompCharge = as.integer(tmp[, 3])
  DefCompMC = match(trimws(tmp[, 4]), MassName)
  DefCompType = as.character(trimws(tmp[, 5]))
  DefCompActCorr = as.character(trimws(tmp[, 6]))

  # concatenate DefComp and Comp
  NComp = NInComp + NDefComp
  CompName = c(InCompName, DefCompName)
  CompCharge = c(CompCharge, DefCompCharge)
  CompMC = c(CompMC, DefCompMC)
  CompType = c(CompType, DefCompType)
  CompActCorr = c(CompActCorr, DefCompActCorr)
  # don't increase the NComp dimension so that we can
  # keep track of which portion of the component list will show up in the
  # interface (NComp) and which does not (NDefComp) BUT . . . that means we need
  # to consider the dimensions of both component lists in subsequent calcs KEC
  # Note: this is already becoming more complicated than it's worth (e.g., in
  # expandWHAM...let's instead increase NComp, then just subtract out NDefComp
  # if we ever need to refer specifically to the input components.

  # Add pH to the component list, if needed
  for (iMass in 1:NMass){#Must have either pH or H in non-BL compartments
    if (!grepl("BL", MassName[iMass])){
      stopifnot(xor(("H" %in% CompName[CompMC == iMass]),
                    ("pH" %in% InVarType[InVarMC == iMass])))
      if ("pH" %in% InVarType[InVarMC == iMass]) {
        NComp = NComp + 1
        CompName = c(CompName,"H")
        CompCharge = c(CompCharge, 1L)
        CompMC = c(CompMC, iMass)
        CompType = c(CompType, "FixedAct")
        CompActCorr = c(CompActCorr, "Debye")
      }
    }
  }

  # error catching
  stopifnot(all(!is.na(CompMC)))
  stopifnot(all(CompType %in% c("MassBal","FixedAct","Substituted","ChargeBal","SurfPot")))
  stopifnot(all(CompActCorr %in% c("None","Debye","Davies","MoleFraction","ChargeFraction","WHAM")))

  # Create variables for species information
  SpecName = character(NSpec)
  SpecMC = integer(NSpec)
  SpecActCorr = integer(NSpec)
  SpecNC = integer(NSpec)                  # the number of components that form species(i)
  CompList = matrix(data = 0, nrow = NSpec, ncol = NComp)  # the list of components (by component number) that form species(i) for SpecNC(i) number of components
  Stoich = matrix(data = 0, nrow = NSpec, ncol = NComp)          # Stoich(i,j) = the amount of component(j) needed to form species(i)
  LogK = numeric(NSpec)
  DeltaH = numeric(NSpec)

  # read species information including stoichiometry, log Ks, etc.
  skipRows = skipRows + NDefComp + 4
  tmp = scan(file = paramFile,skip = skipRows,sep = "\n",nlines = NSpec,what = "character", quiet = T)
  temp.split = strsplit(tmp, ",")
  for (i in 1:NSpec) {
    SpecName[i] = as.character(trimws(temp.split[[i]][1]))
    SpecMC[i] = as.integer(match(trimws(temp.split[[i]][2]), MassName))
    SpecActCorr[i] = as.character(trimws(temp.split[[i]][3]))
    SpecNC[i] = as.integer(trimws(temp.split[[i]][4]))
    for (j in 1:SpecNC[i]) {
      CompList[i, j] = match(trimws(temp.split[[i]][5 + (j - 1) * 2]), CompName)
      Stoich[i, CompList[i, j]] = as.integer(trimws(temp.split[[i]][6 + (j - 1) * 2]))
    }
    LogK[i] = as.numeric(trimws(temp.split[[i]][5 + SpecNC[i] * 2]))
    DeltaH[i] = as.numeric(trimws(temp.split[[i]][6 + SpecNC[i] * 2]))
  }

  # Add components to the species list
  NSpec = NComp + NSpec
  SpecName = c(CompName, SpecName)
  SpecMC = c(CompMC, SpecMC)
  SpecActCorr = c(CompActCorr, SpecActCorr)
  SpecNC = c(array(1L, NComp), SpecNC)
  tmp = matrix(0, nrow = NComp, ncol = NComp)
  tmp[, 1] = 1:(NComp)
  CompList = rbind(tmp, CompList)
  Stoich = rbind(diag(1, nrow = NComp, ncol = NComp), Stoich)
  dimnames(Stoich) = list(Spec = SpecName, Comp = CompName)
  LogK = c(rep(0, NComp), LogK)
  DeltaH = c(rep(0, NComp), DeltaH)

  # Trim down CompList...give it 2 extra columns in case any WHAM species need it
  CompList = CompList[, 1:(max(which(colSums(CompList)>0)) + 2)]

  # -get special species/component information
  # --name of metal
  # --name of biotic ligand
  # --name of BL-metal complex(es)
  # --name of DOC component(s)

  # --name of WHAM file
  if (any(grepl("WHAM",InVarType))){
    stopifnot(length(which(grepl("WHAM",InVarType))) == 1)
    tmp = gsub("WHAM","",InVarType[which(grepl("WHAM",InVarType))])
    if (tmp %in% c("V","VI","VII")){
      WHAMVer = tmp
      wdatFile = NULL
    } else {
      WHAMVer = NULL
      wdatFile = substr(tmp,2,nchar(tmp))
    }
  }
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
    NInLab = NInLab,
    NInVar = NInVar,
    NInComp = NInComp,
    NDefComp = NDefComp,
    NComp = NComp,
    NSpec = NSpec,
    NPhases = NPhases,
    NBLD = NBLD,
    NCAT = NCAT,

    MassName = MassName,
    MassAmt = MassAmt,
    MassUnit = MassUnit,

    InLabName = InLabName,

    InVarName = InVarName,
    InVarMC = InVarMC,
    InVarType = InVarType,

    InCompName = InCompName,

    CompName = CompName,
    CompCharge = CompCharge,
    CompMC = CompMC,
    CompType = CompType,
    CompActCorr = CompActCorr,
    CompSiteDens = CompSiteDens,

    # Defined Components
    DefCompName = DefCompName,
    DefCompFrom = DefCompFrom,
    DefCompCharge = DefCompCharge,
    DefCompMC = DefCompMC,
    DefCompType = DefCompType,
    DefCompActCorr = DefCompActCorr,
    DefCompSiteDens = DefCompSiteDens,

    #Formation Reactions
    SpecName = SpecName,
    SpecMC = SpecMC,
    SpecActCorr = SpecActCorr,
    SpecNC = SpecNC,
    CompList = CompList,
    Stoich = Stoich,
    LogK = LogK,
    DeltaH = DeltaH

    #Phase List
    #Biotic Ligand Definitions
    #Critical Accumultion Table
  )

  # Expand WHAM components and species, if needed
  if (any(InVarType == "WHAM")) {
    out2 = do.call("expandWHAM",
                   args = c(out[which(names(out) %in% formalArgs("expandWHAM"))],
                            list(
                              wdatFile = wdatFile, WHAMVer = WHAMVer
                            )))
    out[names(out2)] = out2
  }

  return(out)
}
