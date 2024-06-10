#' Define the speciation problem
#'
#' `DefineProblem` reads in a parameter file, and sets up the required vectors
#' and matrices that will be needed to run the speciation calculations in CHESS.
#'
#' @param ParamFile the path and file name to a parameter file
#' @param WriteLog if TRUE, the CHESS.LOG file will be written, summarizing the
#'   current problem
#'
#' @returns Returns a `list` object with each list item named according to
#'
#' \describe{
#'  \item{\code{NComp}}{integer; the number of components}
#'  \item{\code{NSpec}}{integer; the number of species}
#'  \item{\code{CompName}}{character vector of length `NComp`; component names}
#'  \item{\code{CompCharge}}{integer vector of length `NComp`; the charge of
#'    the components as free ions}
#'  \item{\code{CompType}}{character vector of length `NComp`; the type of
#'    component. It should be a fixed set of values ("MassBal", "FixedAct",
#'    "FixedConc", "Substituted", "ChargeBal", "SurfPot", "DonnanChargeBal")}
#'  \item{\code{CompActCorr}}{character vector of length `NComp`; the method to
#'    use for activity corrections with this component,  }
#'  \item{\code{SpecName}}{character vector of length `NSpec`; species names}
#'  \item{\code{SpecActCorr}}{integer vector of length `NSpec`; the method to
#'    use for activity corrections with this speies where 1 = ...}
#'  \item{\code{SpecStoich}}{the stoichiometry matrix of the reactions
#'    (integer matrix of size `NSpec` x `nComp`)}
#'  \item{\code{SpecK}}{the equilibrium coefficients (numeric vector of length
#'    `NSpec`)}
#'  \item{\code{SpeclogK}}{the log10-transformed equilibrium coefficients
#'    (numeric vector of length `NSpec`)}
#'  \item{\code{SpecDeltaH}}{the enthalpy change for each formation reaction
#'    (species)}
#'
#'
#'  \item{\code{NMass, NInLab, NInVar, NInComp, NDefComp, NComp, NSpec, NPhase,
#'    NBL, NMetal, NBLMetal, NCAT}}{The counts of the mass
#'    compartments, input label fields, input variables, input components,
#'    defined components, total simulation components, species reactions,
#'    phases, special definitions, biotic ligand components associated with
#'    toxic effects, metal components associated with toxic effects, biotic
#'    ligand-bound metal species associated with toxic effects, and critical
#'    accumulations in the parameter file table.}
#'  \item{\code{MassName, MassAmt, MassUnit}}{The name, amount, and units for
#'    each mass compartment.}
#'  \item{\code{InLabName}}{The names of the input label fields.}
#'  \item{\code{InVarName, InVarMCR, InVarType}}{The name, mass compartment
#'    number, and type of each input variable.}
#'  \item{\code{InCompName}}{The names of the input components.}
#'  \item{\code{CompName, CompCharge, CompMCR, CompType, CompActCorr,
#'    CompSiteDens}}{The names, charges, mass compartments, types, activity
#'    correction methods, and site densities of each component.}
#'  \item{\code{DefCompName, DefCompFromNum, DefCompFromVar, DefCompCharge,
#'    DefCompMCR, DefCompType, DefCompActCorr, DefCompSiteDens}}{The ... of
#'    defined components.}
#'  \item{\code{SpecName, SpecMCR, SpecActCorr, SpecNC, SpecCompList, SpecStoich,
#'    SpecLogK, SpecDeltaH, SpecTempKelvin, SpecCtoM, SpecCharge}}{The...of species
#'    formation reactions}
#'  \item{\code{PhaseName, PhaseNC, PhaseCompList, PhaseStoich, PhaseLogK,
#'    PhaseDeltaH, PhaseTemp, PhaseMoles}}{The...of phases in the phase list.}
#'  \item{\code{BLName}}{The name of the component that corresponds to the
#'    biotic ligand associated with toxic effects.}
#'  \item{\code{MetalName}}{The name of the component that corresponds to the
#'    metal associated with toxic effects.}
#'  \item{\code{BLMetalName}}{The names of the species that are the
#'    biotic ligand-bound metal associated with toxic effects.}
#'  \item{\code{DoWHAM}}{logical, TRUE = there are WHAM species, FALSE = no WHAM
#'    species}
#'  \item{\code{CATab}}{data frame, the critical accumulation table from the
#'    parameter file.}
#' }
#'
#' @export
#'
#' @examples
#' ## Not Run
#' # thisProblem = DefineProblem("my_parameter_file.dat")
#' ## End Not Run
DefineProblem = function(ParamFile, WriteLog = FALSE) {

  # error catching
  stopifnot(file.exists(ParamFile))

  # read in parameter file
  # -get number of mass compartments
  # -get number of input variables
  # -get number of input components
  # -get number of defined components
  # -get number of species
  # -get number of phases
  # -get component information

  # read the dimensions of the various elements of the reaction list
  SkipRows = 2
  Tmp = read.csv(file = ParamFile, header = FALSE, skip = SkipRows,
                 nrows = 9, strip.white = TRUE)
  NMass = as.integer(Tmp[1, 1])
  NInLab = as.integer(Tmp[2, 1])
  NInVar = as.integer(Tmp[3, 1])
  NInComp = as.integer(Tmp[4, 1])
  NDefComp = as.integer(Tmp[5, 1])
  NSpec = as.integer(Tmp[6, 1])
  NPhase = as.integer(Tmp[7, 1])
  NSpecialDef = as.integer(Tmp[8, 1])
  NCAT = as.integer(Tmp[9, 1])
  stopifnot(NMass > 0, NInLab > 0, NInVar > 0, NInComp > 0, NSpec > 0)

  # read mass compartment list
  SkipRows = SkipRows + 9 + 2
  Tmp = read.csv(file = ParamFile, header = TRUE, skip = SkipRows,
                 nrows = NMass, strip.white = TRUE)
  MassName = as.character(trimws(Tmp[, 1]))
  MassAmt = as.numeric(Tmp[, 2])
  MassUnit = as.character(trimws(Tmp[, 3]))
  AqueousMCR = which(tolower(MassName) %in% c("water", "aqueous"))
  BioticLigMCR = which(grepl("BL", MassName, ignore.case = TRUE) |
                        grepl("gill", MassName, ignore.case = TRUE))

  # read input Labels
  SkipRows = SkipRows + NMass + 3
  Tmp = scan(file = ParamFile, what = character(), nlines = NInLab, sep = "\n",
             skip = SkipRows, quiet = TRUE, strip.white = TRUE)
  InLabName = as.character(Tmp)

  # read input variables -
  SkipRows = SkipRows + NInLab + 2
  Tmp = read.csv(file = ParamFile, header = TRUE, skip = SkipRows,
                 nrows = NInVar, strip.white = TRUE)
  InVarName = as.character(trimws(Tmp[, 1]))
  InVarMCName = as.character(trimws(Tmp[, 2]))
  InVarMCR = match(InVarMCName, MassName)
  InVarType = as.character(trimws(Tmp[, 3]))
  stopifnot(!any(duplicated(InVarName)))
  stopifnot(all(!is.na(InVarMCR)))
  stopifnot(all(InVarType %in% c("Temperature", "pH", "WHAM-FA", "WHAM-HA",
                                 "WHAM-HAFA", "PercHA", "PercAFA")))
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
  SkipRows = SkipRows + NInVar + 3
  Tmp = read.csv(file = ParamFile, header = TRUE, skip = SkipRows,
                 nrows = NInComp, strip.white = TRUE)
  NComp = NInComp
  InCompName = CompName = as.character(trimws(Tmp[, 1]))
  CompCharge = as.integer(Tmp[, 2])
  CompMCName = as.character(trimws(Tmp[, 3]))
  CompMCR = match(CompMCName, MassName)
  CompType = as.character(trimws(Tmp[, 4]))
  CompActCorr = as.character(trimws(Tmp[, 5]))
  CompSiteDens = array(1.0, dim = NInComp)

  # read defined component list and properties
  SkipRows = SkipRows + NInComp + 3
  if (NDefComp > 0) {
    Tmp = read.csv(file = ParamFile, header = TRUE, skip = SkipRows,
                   nrows = NDefComp, strip.white = TRUE)
    DefCompName = as.character(trimws(Tmp[, 1]))
    NumericDefComp = grepl("^[[:digit:].]+[eE]?[+-]?[[:digit:]]?", Tmp[, 2])
    DefCompFromNum = as.numeric(array(NA, dim = NDefComp))
    DefCompFromNum[NumericDefComp] = as.numeric(Tmp[NumericDefComp, 2])
    DefCompFromVar = as.character(array(NA, dim = NDefComp))
    DefCompFromVar[!NumericDefComp] = trimws(Tmp[!NumericDefComp, 2])
    DefCompFromVar[!is.na(DefCompFromNum)] = NA
    DefCompFromNum[!is.na(DefCompFromVar)] = NA
    DefCompCharge = as.integer(Tmp[, 3])
    DefCompMCName = as.character(trimws(Tmp[, 4]))
    DefCompMCR = match(DefCompMCName, MassName)
    DefCompType = as.character(trimws(Tmp[, 5]))
    DefCompActCorr = as.character(trimws(Tmp[, 6]))
    DefCompSiteDens = as.numeric(Tmp[, 7])
  } else {
    DefCompName = character()
    DefCompFromNum = numeric()
    DefCompFromVar = character()
    DefCompCharge = integer()
    DefCompMCName = character()
    DefCompMCR = integer()
    DefCompType = character()
    DefCompActCorr = character()
    DefCompSiteDens = numeric()
  }

  # Add pH to the component list, if needed
  for (iMass in 1:NMass){#Must have either pH or H in non-BL compartments
    if (grepl("Water", MassName[iMass], ignore.case = TRUE)) {
      if (!xor(("H" %in% CompName[CompMCR == iMass]),
                    ("pH" %in% InVarType[InVarMCR == iMass]) &
                      ("H" %in% DefCompName[DefCompMCR == iMass]))) {
        stop("Must specify either H as an input component, or pH with a
             corresponding H defined component.")
      }
      # if ("pH" %in% InVarType[InVarMCR == iMass]) {
      #   NComp = NComp + 1L
      #   CompName = c(CompName, "H")
      #   CompCharge = c(CompCharge, 1L)
      #   CompMCR = c(CompMCR, iMass)
      #   CompMCName = c(CompMCName, MassName[iMass])
      #   CompType = c(CompType, "FixedAct")
      #   CompActCorr = c(CompActCorr, "Debye")
      #   CompSiteDens = c(CompSiteDens, 1.0)
      # }
    } else {
      if (("H" %in% CompName[CompMCR == iMass]) ||
          ("pH" %in% InVarType[InVarMCR == iMass])) {
        stop("pH/[H+] specified for non-water mass compartment.")
      }
    }
  }


  # concatenate DefComp and Comp
  NComp = NComp + NDefComp
  CompName = c(CompName, DefCompName)
  CompCharge = c(CompCharge, DefCompCharge)
  CompMCR = c(CompMCR, DefCompMCR)
  CompMCName = c(CompMCName, DefCompMCName)
  CompType = c(CompType, DefCompType)
  CompActCorr = c(CompActCorr, DefCompActCorr)
  CompSiteDens = c(CompSiteDens, DefCompSiteDens)


  # Create variables for species information
  SpecName = character(NSpec)
  SpecMCR = integer(NSpec)
  SpecMCName = character(NSpec)
  SpecActCorr = integer(NSpec)
  SpecNC = integer(NSpec)
  SpecCompList = matrix(data = 0, nrow = NSpec, ncol = NComp)
  SpecStoich = matrix(data = 0, nrow = NSpec, ncol = NComp)
  SpecLogK = numeric(NSpec)
  SpecDeltaH = numeric(NSpec)
  SpecTempKelvin = numeric(NSpec)

  # read species information including stoichiometry, log Ks, etc.
  SkipRows = SkipRows + NDefComp + 4
  Tmp = scan(file = ParamFile, skip = SkipRows, sep = "\n", nlines = NSpec,
             what = "character", quiet = TRUE)
  TmpSplit = strsplit(Tmp, ",")
  for (i in 1:NSpec) {
    SpecName[i] = as.character(trimws(TmpSplit[[i]][1]))
    SpecMCName[i] = as.character(trimws(TmpSplit[[i]][2]))
    SpecMCR[i] = as.integer(match(SpecMCName[i], MassName))
    SpecActCorr[i] = as.character(trimws(TmpSplit[[i]][3]))
    SpecNC[i] = as.integer(trimws(TmpSplit[[i]][4]))
    for (j in 1:SpecNC[i]) {
      SpecCompList[i, j] =
        match(trimws(TmpSplit[[i]][5 + (j - 1) * 2]), CompName)
      SpecStoich[i, SpecCompList[i, j]] =
        as.integer(trimws(TmpSplit[[i]][6 + (j - 1) * 2]))
    }
    SpecLogK[i] = as.numeric(trimws(TmpSplit[[i]][5 + SpecNC[i] * 2]))
    SpecDeltaH[i] = as.numeric(trimws(TmpSplit[[i]][6 + SpecNC[i] * 2]))
    SpecTempKelvin[i] = as.numeric(trimws(TmpSplit[[i]][7 + SpecNC[i] * 2]))
  }

  # -Get Phase information
  SkipRows = SkipRows + NSpec + 3
  # Create variables for Phase information
  PhaseName = character(NPhase)
  PhaseNC = integer(NPhase)
  PhaseCompList = matrix(data = 0, nrow = NPhase, ncol = NComp)
  PhaseStoich = matrix(data = 0, nrow = NPhase, ncol = NComp)
  PhaseLogK = numeric(NPhase)
  PhaseDeltaH = numeric(NPhase)
  PhaseTemp = numeric(NPhase)
  PhaseMoles = numeric(NPhase)
  if (NPhase > 0) {
    # read Phase information including stoichiometry, log Ks, etc.
    Tmp = scan(file = ParamFile, skip = SkipRows, sep = "\n", nlines = NPhase,
               what = "character", quiet = TRUE)
    TmpSplit = strsplit(Tmp, ",")
    for (i in 1:NPhase) {
      PhaseName[i] = as.character(trimws(TmpSplit[[i]][1]))
      PhaseNC[i] = as.integer(trimws(TmpSplit[[i]][2]))
      for (j in 1:PhaseNC[i]) {
        PhaseCompList[i, j] =
          match(trimws(TmpSplit[[i]][3 + (j - 1) * 2]), CompName)
        PhaseStoich[i, PhaseCompList[i, j]] =
          as.integer(trimws(TmpSplit[[i]][4 + (j - 1) * 2]))
      }
      PhaseLogK[i] = as.numeric(trimws(TmpSplit[[i]][3 + PhaseNC[i] * 2]))
      PhaseDeltaH[i] = as.numeric(trimws(TmpSplit[[i]][4 + PhaseNC[i] * 2]))
      PhaseTemp[i] = as.numeric(trimws(TmpSplit[[i]][5 + PhaseNC[i] * 2]))
      PhaseMoles[i] = as.numeric(trimws(TmpSplit[[i]][6 + PhaseNC[i] * 2]))
    }
  }

  # -get Special definitions
  SkipRows = SkipRows + NPhase + 2
  if (NSpecialDef > 0) {
    Tmp = read.csv(file = ParamFile, header = TRUE, skip = SkipRows,
                   nrows = NSpecialDef, strip.white = TRUE)
  } else {
    Tmp = data.frame(matrix(nrow = 0, ncol = 2))
  }
  # --name of metal
  NMetal = sum(Tmp[, 1] == "Metal")
  MetalName = as.character(trimws(Tmp[Tmp[, 1] == "Metal", 2]))
  stopifnot(all(MetalName %in% CompName))

  # --name of biotic ligand
  NBL = sum(Tmp[, 1] == "BL")
  BLName = as.character(trimws(Tmp[Tmp[, 1] == "BL", 2]))
  stopifnot(all(BLName %in% CompName))

  # --name of BL-metal complex(es)
  NBLMetal = sum(Tmp[, 1] == "BL-Metal")
  BLMetalName = as.character(trimws(Tmp[Tmp[, 1] == "BL-Metal", 2]))
  stopifnot(all(BLMetalName %in% SpecName))

  #--Name of WHAM file or WHAM version
  if (any(grepl("WHAM", InVarType))) {
    if ("WHAM" %in% Tmp[, 1] == FALSE) {
      stop("WHAM input variable specified withOut WHAM version or file.")
    }
    if (sum("WHAM" %in% Tmp[, 1]) > 1) {
      stop("More than one WHAM version or file specified.")
    }
    DoWHAM = TRUE
    WHAMString = Tmp[Tmp[, 1] == "WHAM", 2]
    if (WHAMString %in% c("V", "VI", "VII")) {
      WHAMVer = WHAMString
      WdatFile = NA
    } else {
      WHAMVer = NA
      WdatFile = WHAMString
    }
  } else {
    DoWHAM = FALSE
  }


  # -get critical accumulation information
  # --> this part also needs to happen in ListCAT function
  SkipRows = SkipRows + NSpecialDef + 3
  if (NCAT > 0) {
    CATab = read.csv(file = ParamFile, header = TRUE, skip = SkipRows,
                     nrows = NCAT, strip.white = TRUE)
    colnames(CATab) = c("Num", "CA", "Species", "Test.Type", "Duration",
                        "Lifestage", "Endpoint", "Quantifier", "References",
                        "Miscellaneous")
  } else {
    CATab = data.frame(Num = integer(), CA = numeric(), Species = character(),
                       Test.Type = character(), Duration = character(),
                       Lifestage = character(), Endpoint = character(),
                       Quantifier = character(), References = character(),
                       Miscellaneous = character())
  }

  # Add components to the species list
  NSpec = NComp + NSpec
  SpecName = c(CompName, SpecName)
  SpecMCName = c(CompMCName, SpecMCName)
  SpecMCR = c(CompMCR, SpecMCR)
  SpecActCorr = c(CompActCorr, SpecActCorr)
  SpecNC = c(array(1L, NComp), SpecNC)
  Tmp = matrix(0, nrow = NComp, ncol = NComp)
  Tmp[, 1] = 1:(NComp)
  SpecCompList = rbind(Tmp, SpecCompList)
  SpecStoich = rbind(diag(1, nrow = NComp, ncol = NComp), SpecStoich)
  dimnames(SpecStoich) = list(Spec = SpecName, Comp = CompName)
  SpecLogK = c(rep(0, NComp), SpecLogK)
  SpecDeltaH = c(rep(0, NComp), SpecDeltaH)
  SpecTempKelvin = c(rep(0, NComp), SpecTempKelvin)

  # Trim down CompList
  SpecCompList = SpecCompList[, 1:max(which(colSums(SpecCompList) > 0))]

  # Final positions of special definition parameters
  MetalCompR = which(SpecName %in% MetalName)
  BLCompR = which(SpecName %in% BLName)

  # error catching
  stopifnot(
    all(!is.na(c(CompMCR, SpecMCR))),
    all(CompType %in% c("MassBal", "FixedAct", "FixedConc", "Substituted", "ChargeBal",
                        "SurfPot", "DonnanChargeBal")),
    all(c(CompActCorr, SpecActCorr) %in% c("None", "Debye", "Davies", "WHAM"))
  )

  # assemble Output
  Out = list(

    # Counts
    NMass = NMass,
    NInLab = NInLab,
    NInVar = NInVar,
    NInComp = NInComp,
    NDefComp = NDefComp,
    NComp = NComp,
    NSpec = NSpec,
    NPhase = NPhase,
    NBL = NBL,
    NMetal = NMetal,
    NBLMetal = NBLMetal,
    NCAT = NCAT,

    # Mass Compartment List
    MassName = MassName,
    MassAmt = MassAmt,
    MassUnit = MassUnit,
    AqueousMCR = AqueousMCR,
    BioticLigMCR = BioticLigMCR,

    # Input Labels
    InLabName = InLabName,

    # Input Variables
    InVarName = InVarName,
    InVarMCName = InVarMCName,
    InVarMCR = InVarMCR,
    InVarType = InVarType,

    # Input Components
    InCompName = InCompName,
    CompName = CompName,
    CompCharge = CompCharge,
    CompMCName = CompMCName,
    CompMCR = CompMCR,
    CompType = CompType,
    CompActCorr = CompActCorr,
    CompSiteDens = CompSiteDens,

    # Defined Components
    DefCompName = DefCompName,
    DefCompFromNum = DefCompFromNum,
    DefCompFromVar = DefCompFromVar,
    DefCompCharge = DefCompCharge,
    DefCompMCName = DefCompMCName,
    DefCompMCR = DefCompMCR,
    DefCompType = DefCompType,
    DefCompActCorr = DefCompActCorr,
    DefCompSiteDens = DefCompSiteDens,

    # Formation Reactions
    SpecName = SpecName,
    SpecMCName = SpecMCName,
    SpecMCR = SpecMCR,
    SpecActCorr = SpecActCorr,
    SpecNC = SpecNC,
    SpecCompList = SpecCompList,
    SpecStoich = SpecStoich,
    SpecLogK = SpecLogK,
    SpecDeltaH = SpecDeltaH,
    SpecTempKelvin = SpecTempKelvin,

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
    BLCompR = BLCompR,
    MetalName = MetalName,
    MetalCompR = MetalCompR,
    BLMetalName = BLMetalName,
    # BLMetalSpecsR = BLMetalSpecs,#this needs to be figured out after ExpandWHAM
    DoWHAM = DoWHAM,
    WHAMVer = WHAMVer,
    WdatFile = WdatFile,

    # Critical Accumulation Table
    CATab = CATab
  )

  # Expand WHAM components and species, if needed
  if (DoWHAM) {
    Out2 =
      do.call("ExpandWHAM",
              args = c(Out[which(names(Out) %in% formalArgs("ExpandWHAM"))],
                       list(WdatFile = WdatFile, WHAMVer = WHAMVer)))
    Out[names(Out2)] = Out2
  } else {
    Out$WHAMDonnanMCR = array(NA, dim = 2, dimnames = list(c("HA","FA")))
    Out$wDLF = NA
    Out$wKZED = NA
    Out$wP = array(NA, dim = 2, dimnames = list(c("HA","FA")))
    Out$wRadius = array(NA, dim = 2, dimnames = list(c("HA","FA")))
    Out$wMolWt = array(NA, dim = 2, dimnames = list(c("HA","FA")))
  }

  stopifnot(!any(duplicated(c(InLabName, InVarName, SpecName))))

  # Make SpecCtoM, SpecCharge, SpecK
  SpecCharge = drop(Out$SpecStoich %*% Out$CompCharge)
  Out$SpecCharge = SpecCharge
  Out$SpecK = 10 ^ Out$SpecLogK

  # Final positions of special definition parameters
  BLMetalSpecs = which(Out$SpecName %in% BLMetalName)
  Out$BLMetalSpecsR = BLMetalSpecs

  if (WriteLog) {
    CHESSLog(Out, ParamFile)
  }

  return(Out)
}
