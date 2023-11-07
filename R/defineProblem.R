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
#'  \item{\code{CompType}}{character vector of length `NComp`; the type of component. It should be a fixed set of values (MassBal, FixedAct, Substituted, ChargeBal, SurfPot)}
#'  \item{\code{CompActCorr}}{character vector of length `NComp`; the method to use for activity corrections with this component,  }
#'  \item{\code{SpecName}}{character vector of length `NSpec`; species names}
#'  \item{\code{SpecActCorr}}{integer vector of length `NSpec`; the method to use for activity corrections with this speies where 1 = ...}
#'  \item{\code{Stoich}}{the stoichiometry matrix of the reactions (integer matrix of size `NSpec` x `nComp`)}
#'  \item{\code{K}}{the equilibrium coefficients (numeric vector of length `NSpec`)}
#'  \item{\code{logK}}{the log10-transformed equilibrium coefficients (numeric vector of length `NSpec`)}
#'  \item{\code{DeltaH}}{the enthalpy change for each formation reaction (species)}
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
  tmp = read.csv(file = paramFile, header = FALSE, skip = skipRows,
                 nrows = 9, strip.white = T)
  NMass = tmp[1, 1]
  NInLab = tmp[2, 1]
  NInVar = tmp[3, 1]
  NInComp = tmp[4, 1]
  NDefComp = tmp[5, 1]
  NSpec = tmp[6, 1]
  NPhase = tmp[7, 1]
  NSpecialDef = tmp[8, 1]
  NCAT = tmp[9, 1]
  stopifnot(NMass>0, NInLab>0, NInVar>0, NInComp>0, NSpec>0)

  # read compartment list
  skipRows = skipRows + 9 + 2
  tmp = read.csv(file = paramFile, header = TRUE, skip = skipRows,
                 nrows = NMass, strip.white = T)
  MassName = as.character(tmp[, 1])
  MassAmt = as.numeric(tmp[, 2])
  MassUnit = as.character(tmp[, 3])

  # read input Labels
  skipRows = skipRows + NMass + 3
  tmp = scan(file = paramFile, what = character(), nlines = NInLab, sep = "\n",
             skip = skipRows, quiet = T, strip.white = T)
  InLabName = as.character(tmp)

  # read input variables -
  skipRows = skipRows + NInLab + 2
  tmp = read.csv(file = paramFile, header = TRUE, skip = skipRows,
                 nrows = NInVar, strip.white = T)
  InVarName = as.character(tmp[, 1])
  InVarMC = match(tmp[, 2], MassName)
  InVarType = as.character(tmp[, 3])
  stopifnot(!any(duplicated(InVarName)))
  stopifnot(all(!is.na(InVarMC)))
  stopifnot(all(InVarType %in% c("Temperature","pH","WHAM-FA","WHAM-HA",
                                  "WHAM-HAFA","PercHA","PercAFA")))
  stopifnot("Temperature" %in% InVarType)
  # - Temperature = the temperature in degrees C
  # - pH = the -log[H]...you know, pH
  # - WHAM-HA, -FA, -HAFA = Windemere Humic Aqueous Model organic matter (input
  #   mg C/L), as all humic acid, all fulvic acid, or a mix of humics and
  #   fulvics, respectively.
  # - PercHA = optionally indicate the percent humic acid in a the WHAM-HAFA
  #   component for that compartment.
  # - PercAFA = optionally indicate the percent of active fulvic acid for the
  #   WHAM-FA or WHAM-HAFA component for that compartment

  # read component list and properties
  skipRows = skipRows + NInVar + 3
  tmp = read.csv(file = paramFile, header = TRUE, skip = skipRows,
                 nrows = NInComp, strip.white = T)
  NComp = NInComp
  InCompName = CompName = as.character(trimws(tmp[, 1]))
  CompCharge = as.integer(tmp[, 2])
  CompMC = match(trimws(tmp[, 3]), MassName)
  CompType = as.character(trimws(tmp[, 4]))
  CompActCorr = as.character(trimws(tmp[, 5]))
  CompSiteDens = array(1.0, dim = NInComp)

  # Add pH to the component list, if needed
  for (iMass in 1:NMass){#Must have either pH or H in non-BL compartments
    if (grepl("Water", MassName[iMass], ignore.case = T)){
      stopifnot(xor(("H" %in% CompName[CompMC == iMass]),
                    ("pH" %in% InVarType[InVarMC == iMass])))
      if ("pH" %in% InVarType[InVarMC == iMass]) {
        NComp = NComp + 1L
        CompName = c(CompName,"H")
        CompCharge = c(CompCharge, 1L)
        CompMC = c(CompMC, iMass)
        CompType = c(CompType, "FixedAct")
        CompActCorr = c(CompActCorr, "Debye")
        CompSiteDens = c(CompSiteDens, 1.0)
      }
    } else {
      if(("H" %in% CompName[CompMC == iMass]) | ("pH" %in% InVarType[InVarMC == iMass])){
        stop("pH/[H+] specified for non-water mass compartment.")
      }
    }
  }

  # read defined component list and properties
  skipRows = skipRows + NInComp + 3
  if (NDefComp > 0){
    tmp = read.csv(file = paramFile, header = TRUE, skip = skipRows,
                   nrows = NDefComp, strip.white = T)
    DefCompName = as.character(tmp[, 1])
    DefCompFromNum = as.numeric(tmp[, 2]) # we will eventually want to be able to read numbers or strings
    DefCompFromVar = as.character(tmp[, 2]) # we will eventually want to be able to read numbers or strings
    DefCompFromVar[!is.na(DefCompFromNum)] = NA
    DefCompFromNum[!is.na(DefCompFromVar)] = NA
    DefCompCharge = as.integer(tmp[, 3])
    DefCompMC = match(tmp[, 4], MassName)
    DefCompType = as.character(tmp[, 5])
    DefCompActCorr = as.character(tmp[, 6])
    DefCompSiteDens = as.numeric(tmp[, 7])
  } else {
    DefCompName = character()
    DefCompFromNum = numeric()
    DefCompFromVar = character()
    DefCompCharge = integer()
    DefCompMC = integer()
    DefCompType = character()
    DefCompActCorr = character()
    DefCompSiteDens = numeric()
  }

  # concatenate DefComp and Comp
  NComp = NComp + NDefComp
  CompName = c(CompName, DefCompName)
  CompCharge = c(CompCharge, DefCompCharge)
  CompMC = c(CompMC, DefCompMC)
  CompType = c(CompType, DefCompType)
  CompActCorr = c(CompActCorr, DefCompActCorr)
  CompSiteDens = c(CompSiteDens, DefCompSiteDens)

  # Create variables for species information
  SpecName = character(NSpec)
  SpecMC = integer(NSpec)
  SpecActCorr = integer(NSpec)
  SpecNC = integer(NSpec)                  # the number of components that form species(i)
  SpecCompList = matrix(data = 0, nrow = NSpec, ncol = NComp)  # the list of components (by component number) that form species(i) for SpecNC(i) number of components
  SpecStoich = matrix(data = 0, nrow = NSpec, ncol = NComp)          # Stoich(i,j) = the amount of component(j) needed to form species(i)
  SpecLogK = numeric(NSpec)
  SpecDeltaH = numeric(NSpec)
  SpecTemp = numeric(NSpec)

  # read species information including stoichiometry, log Ks, etc.
  skipRows = skipRows + NDefComp + 4
  tmp = scan(file = paramFile, skip = skipRows, sep = "\n", nlines = NSpec,
             what = "character", quiet = T)
  tmp.split = strsplit(tmp, ",")
  for (i in 1:NSpec) {
    SpecName[i] = as.character(trimws(tmp.split[[i]][1]))
    SpecMC[i] = as.integer(match(trimws(tmp.split[[i]][2]), MassName))
    SpecActCorr[i] = as.character(trimws(tmp.split[[i]][3]))
    SpecNC[i] = as.integer(trimws(tmp.split[[i]][4]))
    for (j in 1:SpecNC[i]) {
      SpecCompList[i, j] = match(trimws(tmp.split[[i]][5 + (j - 1) * 2]), CompName)
      SpecStoich[i, SpecCompList[i, j]] = as.integer(trimws(tmp.split[[i]][6 + (j - 1) * 2]))
    }
    SpecLogK[i] = as.numeric(trimws(tmp.split[[i]][5 + SpecNC[i] * 2]))
    SpecDeltaH[i] = as.numeric(trimws(tmp.split[[i]][6 + SpecNC[i] * 2]))
    SpecTemp[i] = as.numeric(trimws(tmp.split[[i]][7 + SpecNC[i] * 2]))
  }

  # -Get Phase information
  skipRows = skipRows + NSpec + 3
  # Create variables for Phase information
  PhaseName = character(NPhase)
  PhaseNC = integer(NPhase)                  # the number of components that form Phase(i)
  PhaseCompList = matrix(data = 0, nrow = NPhase, ncol = NComp)  # the list of components (by component number) that form Phase(i) for PhaseNC(i) number of components
  PhaseStoich = matrix(data = 0, nrow = NPhase, ncol = NComp)          # Stoich(i,j) = the amount of component(j) needed to form Phase(i)
  PhaseLogK = numeric(NPhase)
  PhaseDeltaH = numeric(NPhase)
  PhaseTemp = numeric(NPhase)
  PhaseMoles = numeric(NPhase)
  if (NPhase > 0){
    # read Phase information including stoichiometry, log Ks, etc.
    tmp = scan(file = paramFile, skip = skipRows, sep = "\n", nlines = NPhase,
               what = "character", quiet = T)
    tmp.split = strsplit(tmp, ",")
    for (i in 1:NPhase) {
      PhaseName[i] = as.character(trimws(tmp.split[[i]][1]))
      PhaseNC[i] = as.integer(trimws(tmp.split[[i]][2]))
      for (j in 1:PhaseNC[i]) {
        PhaseCompList[i, j] = match(trimws(tmp.split[[i]][3 + (j - 1) * 2]), CompName)
        PhaseStoich[i, PhaseCompList[i, j]] = as.integer(trimws(tmp.split[[i]][4 + (j - 1) * 2]))
      }
      PhaseLogK[i] = as.numeric(trimws(tmp.split[[i]][3 + PhaseNC[i] * 2]))
      PhaseDeltaH[i] = as.numeric(trimws(tmp.split[[i]][4 + PhaseNC[i] * 2]))
      PhaseTemp[i] = as.numeric(trimws(tmp.split[[i]][5 + PhaseNC[i] * 2]))
      PhaseMoles[i] = as.numeric(trimws(tmp.split[[i]][6 + PhaseNC[i] * 2]))
    }
  }

  # -get Special definitions
  skipRows = skipRows + NPhase + 2
  if (NSpecialDef > 0){
    tmp = read.csv(file = paramFile, header = TRUE, skip = skipRows,
                   nrows = NSpecialDef, strip.white = T)
  } else {
    tmp = data.frame(matrix(nrow = 0, ncol = 2))
  }
  # --name of metal
  NMetal = sum(tmp[,1] == "Metal")
  MetalName = as.character(trimws(tmp[tmp[,1] == "Metal",2]))
  stopifnot(all(MetalName %in% CompName))
  # --name of biotic ligand
  NBL = sum(tmp[,1] == "BL")
  BLName = as.character(trimws(tmp[tmp[,1] == "BL", 2]))
  stopifnot(all(BLName %in% CompName))
  # --name of BL-metal complex(es)
  NBLMetal = sum(tmp[,1] == "BL-Metal")
  BLMetalName = as.character(trimws(tmp[tmp[,1] == "BL-Metal", 2]))
  stopifnot(all(BLMetalName %in% SpecName))
  #--Name of WHAM file or WHAM version
  if (any(grepl("WHAM", InVarType))){
    if ("WHAM" %in% tmp[,1] == F){
      stop("WHAM input variable specified without WHAM version or file.")
    }
    if (sum("WHAM" %in% tmp[,1]) > 1){
      stop("More than one WHAM version or file specified.")
    }
    DoWHAM = TRUE
    WHAMString = tmp[tmp[,1] == "WHAM",2]
    if (WHAMString %in% c("V","VI","VII")){
      WHAMVer = WHAMString
      wdatFile = NULL
    } else {
      WHAMVer = NULL
      wdatFile = WHAMString
    }
  } else {
    DoWHAM = FALSE
  }


  # -get critical accumulation information --> this part also needs to happen in listCAT function
  skipRows = skipRows + NSpecialDef + 3
  if (NCAT > 0){
    CATab = read.csv(file = paramFile, header = TRUE, skip = skipRows,
                     nrows = NCAT, strip.white = T)
    colnames(CATab) = c("Num","CA","Species","Test.Type","Duration","Lifestage",
                        "Endpoint","Quantifier","References","Miscellaneous")
  } else {
    CATab = data.frame(Num=integer(), CA=numeric(), Species=character(),
                       Test.Type=character(), Duration=character(),
                       Lifestage=character(), Endpoint=character(),
                       Quantifier=character(), References=character(),
                       Miscellaneous=character())
  }

  # Add components to the species list
  NSpec = NComp + NSpec
  SpecName = c(CompName, SpecName)
  SpecMC = c(CompMC, SpecMC)
  SpecActCorr = c(CompActCorr, SpecActCorr)
  SpecNC = c(array(1L, NComp), SpecNC)
  tmp = matrix(0, nrow = NComp, ncol = NComp)
  tmp[, 1] = 1:(NComp)
  SpecCompList = rbind(tmp, SpecCompList)
  SpecStoich = rbind(diag(1, nrow = NComp, ncol = NComp), SpecStoich)
  dimnames(SpecStoich) = list(Spec = SpecName, Comp = CompName)
  SpecLogK = c(rep(0, NComp), SpecLogK)
  SpecDeltaH = c(rep(0, NComp), SpecDeltaH)
  SpecTemp = c(rep(0, NComp), SpecTemp)

  # Trim down CompList
  SpecCompList = SpecCompList[, 1:max(which(colSums(SpecCompList)>0))]

  # error catching
  stopifnot(
    all(!is.na(c(CompMC, SpecMC))),
    all(CompType %in% c("MassBal","FixedAct","Substituted","ChargeBal","SurfPot")),
    all(c(CompActCorr, SpecActCorr) %in% c("None","Debye","Davies","WHAM"))
  )



  # assemble output
  out = list(
    # Counts
    NMass = NMass,
    NInLab = NInLab,
    NInVar = NInVar,
    NInComp = NInComp,
    NDefComp = NDefComp,
    NComp = NComp,
    NSpec = NSpec,
    NPhase = NPhase,
    NSpecialDef = NSpecialDef,
    NBL = NBL,
    NMetal = NMetal,
    NBLMetal = NBLMetal,
    NCAT = NCAT,

    # Mass Compartment List
    MassName = MassName,
    MassAmt = MassAmt,
    MassUnit = MassUnit,

    # Input Labels
    InLabName = InLabName,

    # Input Variables
    InVarName = InVarName,
    InVarMC = InVarMC,
    InVarType = InVarType,

    # Input Components
    InCompName = InCompName,
    CompName = CompName,
    CompCharge = CompCharge,
    CompMC = CompMC,
    CompType = CompType,
    CompActCorr = CompActCorr,
    CompSiteDens = CompSiteDens,

    # Defined Components
    DefCompName = DefCompName,
    DefCompFromNum = DefCompFromNum,
    DefCompFromVar = DefCompFromVar,
    DefCompCharge = DefCompCharge,
    DefCompMC = DefCompMC,
    DefCompType = DefCompType,
    DefCompActCorr = DefCompActCorr,
    DefCompSiteDens = DefCompSiteDens,

    # Formation Reactions
    SpecName = SpecName,
    SpecMC = SpecMC,
    SpecActCorr = SpecActCorr,
    SpecNC = SpecNC,
    SpecCompList = SpecCompList,
    SpecStoich = SpecStoich,
    SpecLogK = SpecLogK,
    SpecDeltaH = SpecDeltaH,
    SpecTemp = SpecTemp,

    # Phase List
    PhaseName = PhaseName,
    PhaseNC = PhaseNC,
    PhaseCompList = PhaseCompList,
    PhaseStoich = PhaseStoich,
    PhaseLogK = PhaseLogK,
    PhaseDeltaH = PhaseDeltaH,
    PhaseTemp = PhaseTemp,
    PhaseMoles = PhaseMoles,

    # Special Definitions
    BLName = BLName,
    MetalName = MetalName,
    BLMetalName = BLMetalName,
    DoWHAM = DoWHAM,

    # Critical Accumulation Table
    CATab = CATab
  )

  # Expand WHAM components and species, if needed
  if (DoWHAM){
    out2 = do.call("expandWHAM",
                   args = c(out[which(names(out) %in% formalArgs("expandWHAM"))],
                            list(wdatFile = wdatFile, WHAMVer = WHAMVer)))
    out[names(out2)] = out2
  }

  stopifnot(!any(duplicated(c(InLabName, InVarName, SpecName))))

  SpecCharge = out$SpecStoich %*% out$CompCharge
  out$SpecCharge = SpecCharge

  return(out)
}
