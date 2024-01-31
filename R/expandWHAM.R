#' @title Expand the DOC component into WHAM components
#'
#' @param NMass integer; Number of mass compartments
#' @param MassName character vector of length `NMass`; Names of the mass
#'   compartments
#' @param NInVar integer; Number of input variables
#' @param InVarName character vector of length `NInVar`; Names of input
#'   variables
#' @param InVarMC integer vector of length `NInVar`;  Mass compartments of input
#'   variables
#' @param InVarType character vector of length `NInVar`; Types of input
#'   variables
#' @param NComp integer; Number of components (modified and returned)
#' @param CompName character vector of length `NComp`; component names (modified
#'   and returned)
#' @param CompCharge integer vector of length `NComp`; the charge of the
#'   components as free ions (modified and returned)
#' @param CompMC integer vector of length `NComp`; Which mass compartment the
#'   component belongs to (modified and returned)
#' @param CompType character vector of length `NComp`; the type of component
#'   (modified and returned)
#' @param CompActCorr character vector of length `NComp`; the method to use for
#'   activity corrections with this component (modified and returned)
#' @param CompSiteDens numeric vector of length `NComp`; the binding site
#'   density of each component (modified and returned)
#' @param NDefComp integer; Number of defined components (modified and returned)
#' @param DefCompName character vector of length `NDefComp`; defined component
#'   names (modified and returned)
#' @param DefCompFromNum numeric vector of length `NDefComp`; the number the
#'   defined component is formed from (modified and returned)
#' @param DefCompFromVar character vector of length `NDefComp`; the column used
#'   to form the defined component (modified and returned)
#' @param DefCompCharge integer vector of length `NDefComp`; the charge of the
#'   defined components as free ions (modified and returned)
#' @param DefCompMC integer vector of length `NDefComp`; Which mass compartment
#'   the defined component belongs to (modified and returned)
#' @param DefCompType character vector of length `NDefComp`; the type of defined
#'   component (modified and returned)
#' @param DefCompActCorr character vector of length `NDefComp`; the method to
#'   use for activity corrections with this defined component (modified and
#'   returned)
#' @param DefCompSiteDens numeric vector of length `NDefComp`; the binding site
#'   density of each defined component (modified and returned)
#' @param NSpec Number of species (modified and returned)
#' @param SpecName character vector of length `NSpec`; species names (modified
#'   and returned)
#' @param SpecMC integer vector of length `NSpec`; which mass compartment the
#'   speces belongs to (modified and returned)
#' @param SpecActCorr character vector of length `NSpec`; the method to use for
#'   activity corrections with this species (modified and returned)
#' @param SpecNC integer vector of length `NSpec`; the number of components used
#'   to create a given species (modified and returned)
#' @param SpecCompList integer matrix of `NSpec` rows and `max(SpecNC)+2`
#'   columns; the list of components used to create a given species (modified
#'   and returned)
#' @param SpecStoich integer matrix of `NSpec` rows and `NComp` columns; the
#'   stoichiometry matrix of the formation reactions (modified and returned)
#' @param SpecLogK numeric vector of length `NSpec`; the log10-transformed
#'   equilibrium coefficients (modified and returned)
#' @param SpecDeltaH numeric vector of length `NSpec`; the enthalpy change for
#'   each formation reaction/species (modified and returned)
#' @param SpecTempKelvin numeric vector of length `NSpec`; Tmperature at which
#'   the logK/deltaH were measured (modified and returned) in Kelvin
#' @param NPhase integer; Number of phases
#' @param PhaseCompList integer matrix of `NPhase` rows and `max(PhaseNC)+2`
#'   columns; the list of components used to create a given phase (modified and
#'   returned)
#' @param PhaseStoich integer matrix of `NPhase` rows and `NComp` columns; the
#'   stoichiometry matrix of the phase reactions (modified and returned)
#' @param WHAMVer a character string specifying the WHAM version to use, must be
#'   one of `"V"` (default), `"VI"`, or `"VII"`. Ignored if `WdatFile` is not
#'   `NULL`.
#' @param WdatFile (optional) a character string specifying the file path of a
#'   WHAM parameter file
#'
#' @keywords internal
#'
#' @noRd
ExpandWHAM = function(NMass,
                      MassName,
                      NInVar,
                      InVarName,
                      InVarMC,
                      InVarType,
                      NComp,
                      CompName,
                      CompCharge,
                      CompMC,
                      CompType,
                      CompActCorr,
                      CompSiteDens,
                      NDefComp,
                      DefCompName,
                      DefCompFromNum,
                      DefCompFromVar,
                      DefCompCharge,
                      DefCompMC,
                      DefCompType,
                      DefCompActCorr,
                      DefCompSiteDens,
                      NSpec,
                      SpecName,
                      SpecMC,
                      SpecActCorr,
                      SpecNC,
                      SpecCompList,
                      SpecStoich,
                      SpecLogK,
                      SpecDeltaH,
                      SpecTempKelvin,
                      NPhase,
                      PhaseCompList,
                      PhaseStoich,
                      WHAMVer = c("V", "VI", "VII"),
                      WdatFile = NULL) {

  # error catching and input cleanup
  if (is.null(WdatFile)) {
    WHAMVer = match.arg(WHAMVer)
    if (WHAMVer == "V") {
      WdatFile = system.file("extdata/WHAM_V.wdat",
                             package = "BLMEngineInR",
                             mustWork = TRUE)
    } else if (WHAMVer == "VI") {
      WdatFile = system.file("extdata/WHAM_VI.wdat",
                             package = "BLMEngineInR",
                             mustWork = TRUE)
    } else if (WHAMVer  == "VII") {
      WdatFile = system.file("extdata/WHAM_VII.wdat",
                             package = "BLMEngineInR",
                             mustWork = TRUE)
    }
  } else {
    WdatFile = normalizePath(WdatFile)
    stopifnot(file.exists(WdatFile))
  }


  # read WHAM data file-------------------------------------
  # header info
  SkipRows = 2
  Tmp = read.delim(
    file = WdatFile,
    header = FALSE,
    sep = ",",
    skip = SkipRows,
    nrows = 6
  )
  nMS = as.integer(Tmp[1, 2])#Number of monodentate sites#nolint: object_name_linter, line_length_linter.
  nBP = as.integer(Tmp[2, 2])#Number of bidentate pairs#nolint: object_name_linter, line_length_linter.
  nTG = as.integer(Tmp[3, 2])#Number of tridentate groups#nolint: object_name_linter, line_length_linter.
  nMP = as.integer(Tmp[4, 2])#Number of metals-OM parameters#nolint: object_name_linter, line_length_linter.
  wDLF = as.numeric(Tmp[5, 2])#Double layer overlap factor#nolint: object_name_linter, line_length_linter.
  wKZED = as.numeric(Tmp[6, 2])#Constant to control DDL at low ZED#nolint: object_name_linter, line_length_linter.

  # Parameters
  SkipRows = SkipRows + 6 + 1
  Tmp = read.delim(
    file = WdatFile,
    header = TRUE,
    sep = ",",
    skip = SkipRows,
    nrows = 12
  )
  nCOOH = array(as.numeric(Tmp[1, 3:4]), dimnames = list(c("HA", "FA"))) # nolint: object_name_linter, line_length_linter.
  pKHA = array(as.numeric(Tmp[2, 3:4]), dimnames = list(c("HA", "FA"))) # nolint: object_name_linter, line_length_linter.
  pKHB = array(as.numeric(Tmp[3, 3:4]), dimnames = list(c("HA", "FA"))) # nolint: object_name_linter, line_length_linter.
  dpKHA = array(as.numeric(Tmp[4, 3:4]), dimnames = list(c("HA", "FA"))) # nolint: object_name_linter, line_length_linter.
  dpKHB = array(as.numeric(Tmp[5, 3:4]), dimnames = list(c("HA", "FA"))) # nolint: object_name_linter, line_length_linter.
  fprB = array(as.numeric(Tmp[6, 3:4]), dimnames = list(c("HA", "FA"))) # nolint: object_name_linter, line_length_linter.
  fprT = array(as.numeric(Tmp[7, 3:4]), dimnames = list(c("HA", "FA"))) # nolint: object_name_linter, line_length_linter.
  dLK1A = array(as.numeric(Tmp[8, 3:4]), dimnames = list(c("HA", "FA"))) # nolint: object_name_linter, line_length_linter.
  dLK1B = array(as.numeric(Tmp[9, 3:4]), dimnames = list(c("HA", "FA"))) # nolint: object_name_linter, line_length_linter.
  wP = array(as.numeric(Tmp[10, 3:4]), dimnames = list(c("HA", "FA"))) # nolint: object_name_linter, line_length_linter.
  wRadius = array(as.numeric(Tmp[11, 3:4]), dimnames = list(c("HA", "FA"))) # nolint: object_name_linter, line_length_linter.
  wMolWt = array(as.numeric(Tmp[12, 3:4]), dimnames = list(c("HA", "FA"))) # nolint: object_name_linter, line_length_linter.

  # Monodentate Sites - these should always be the same, but we'll set things
  # up like this so we can add in this section if it's ever needed.
  # MonodentTable = data.frame(S=1:8, AbundDenom = c(rep(4,4),rep(8,4)))
  SkipRows = SkipRows + 12 + 3
  MonodentTable = read.delim(
    file = WdatFile,
    header = TRUE,
    sep = ",",
    skip = SkipRows,
    nrows = nMS
  )
  MonodentTable$FullyProt = paste0(MonodentTable$S, "H")
  MonodentTable$FullyDeprot = paste0(MonodentTable$S)
  MonodentTable$Strong1Weak2 = c(rep(1L, 4), rep(2L, 4))

  # Bidentate Pairs
  SkipRows = SkipRows + nMS + 3
  if (nBP > 0) {
    BidentTable = read.delim(
      file = WdatFile,
      header = TRUE,
      sep = ",",
      skip = SkipRows,
      nrows = nBP
    )
    names(BidentTable) = c("S1", "S2", "AbundDenom")
    BidentTable$FullyProt = paste0(BidentTable$S1, BidentTable$S2, "H")
    BidentTable$S1Deprot = paste0(BidentTable$S1, "-", BidentTable$S2, "H")
    BidentTable$S2Deprot = paste0(BidentTable$S2, "-", BidentTable$S1, "H")
    BidentTable$FullyDeprot = paste0(BidentTable$S1, BidentTable$S2)
    BidentTable$S1Strong1Weak2 = ifelse(BidentTable$S1 <= 4, 1, 2)
    BidentTable$S2Strong1Weak2 = ifelse(BidentTable$S2 <= 4, 1, 2)
  } else {
    BidentTable = data.frame()
  }


  # Tridentate Groups
  SkipRows = SkipRows + nBP + 3
  if (nTG > 0) {
    TridentTable = read.delim(
      file = WdatFile,
      header = TRUE,
      sep = ",",
      skip = SkipRows,
      nrows = nTG
    )
    names(TridentTable) = c("S1", "S2", "S3", "AbundDenom")
    TridentTable$FullyProt = paste0(TridentTable$S1, TridentTable$S2,
                                    TridentTable$S3, "H")
    TridentTable$S1Deprot = paste0(TridentTable$S1,
                                   "-",
                                   TridentTable$S2,
                                   TridentTable$S3,
                                   "H")
    TridentTable$S2Deprot = paste0(TridentTable$S2,
                                   "-",
                                   TridentTable$S1,
                                   TridentTable$S3,
                                   "H")
    TridentTable$S3Deprot = paste0(TridentTable$S3,
                                   "-",
                                   TridentTable$S1,
                                   TridentTable$S2,
                                   "H")
    TridentTable$S12Deprot = paste0(TridentTable$S1,
                                    TridentTable$S2,
                                    "-",
                                    TridentTable$S3,
                                    "H")
    TridentTable$S13Deprot = paste0(TridentTable$S1,
                                    TridentTable$S3,
                                    "-",
                                    TridentTable$S2,
                                    "H")
    TridentTable$S23Deprot = paste0(TridentTable$S2,
                                    TridentTable$S3,
                                    "-",
                                    TridentTable$S1,
                                    "H")
    TridentTable$FullyDeprot = paste0(TridentTable$S1, TridentTable$S2,
                                      TridentTable$S3)
    TridentTable$S1Strong1Weak2 = ifelse(TridentTable$S1 <= 4, 1, 2)
    TridentTable$S2Strong1Weak2 = ifelse(TridentTable$S2 <= 4, 1, 2)
    TridentTable$S3Strong1Weak2 = ifelse(TridentTable$S3 <= 4, 1, 2)
  } else {
    TridentTable = data.frame()
  }

  # Metals Parameters Table
  SkipRows = SkipRows + nTG + 3
  if (nMP > 0) {
    MetalsTable = read.delim(
      file = WdatFile,
      header = TRUE,
      sep = ",",
      skip = SkipRows,
      nrows = nMP
    )
    names(MetalsTable) = c("Metal", "pKMAHA", "pKMAFA", "dLK2", "KselHA", "KselFA")
    MetalsTable = MetalsTable[MetalsTable$Metal %in% c(CompName, SpecName), ]
    nMP = nrow(MetalsTable) # nolint: object_name_linter.
    MetalsTable$pKMBHA = 3 * MetalsTable$pKMAHA - 3
    MetalsTable$pKMBFA = 3.96 * MetalsTable$pKMAFA
  } else {
    MetalsTable = data.frame()
  }


  # Do the expansion ---------------------------------------

  # Initialize variables
  iH = which(CompName == "H") # nolint: object_name_linter.
  DonnanMC = array(NMass + c(1L, 2L), dim = 2, dimnames = list(c("HA", "FA")))
  SpecKsel = array(NA, dim = c(NSpec, 2), dimnames = list(SpecName, c("HA", "FA")))
  if (nMP > 0) {
    SpecKsel[match(MetalsTable$Metal, SpecName), ] =
      array(unlist(MetalsTable[, c("KselHA", "KselFA")]), dim = c(nMP, 2))
  }


  # Figure out the number of DOC components we're adding, and what fraction
  InVarWHAM = which(grepl("WHAM", InVarType))

  for (iInVar in InVarWHAM) {

    iMass = InVarMC[iInVar] # nolint: object_name_linter.

    if (InVarType[iInVar] == "WHAM-HA") {
      WHAMFracAdd = c("HA")
    } else if (InVarType[iInVar] == "WHAM-FA") {
      WHAMFracAdd = c("FA")
    } else if (InVarType[iInVar] == "WHAM-HAFA") {
      WHAMFracAdd = c("HA", "FA")
      if (!any(InVarType[InVarMC == InVarMC[iInVar]] %in% "PercHA")) {
        stop("Must have PercHA input variable in mass compartment if specifying a WHAM-HAFA input variable.") # nolint: line_length_linter.
      }
    }
    if ((InVarType[iInVar] %in% c("WHAM-FA", "WHAM-HA")) &&
          any(InVarType[InVarMC == InVarMC[iInVar]] %in% "PercHA")) {
      stop("PercHA input variable specified in mass compartment with WHAM-HA or WHAM-FA input variable.") # nolint: line_length_linter.
    }
    if ((InVarType[iInVar] %in% c("WHAM-HA")) &&
        any(InVarType[InVarMC == InVarMC[iInVar]] %in% "PercAFA")) {
      stop("PercAFA input variable specified in mass compartment with WHAM-HA input variable.") # nolint: line_length_linter.
    }
    NWHAMFracAdd = length(WHAMFracAdd)

    WHAMprefix = array(
      paste0(InVarName[iInVar], "-", WHAMFracAdd, "_"),
      dim = NWHAMFracAdd,
      dimnames = list(WHAMFracAdd)
    )



    # * Each component has the fully protonated species as the component
    # * Each component will have every possible combination of binding sites
    #   deprotonated
    # * Each component, when fully deprotonated, will bind to each metal nMP.
    #     monodentate example:
    #      component: FA1H
    #      species: FA1H, FA1, FA1-Mg, FA1-Ca, ...
    #     bidentate example:
    #      component: FA12H
    #      species: FA12H, FA1-2H, FA2-1H, FA12, FA12-Mg, FA12-Ca, ...
    #     tridentate example:
    #      component: FA123H
    #      species: FA123H, FA1-23H, FA2-13H, FA3-12H, FA23-1H, FA13-2H,
    #               FA12-3H, FA123, FA123-Mg, FA123-Ca, ...
    WNComp = (nMS + nBP + nTG) * NWHAMFracAdd
    StartComp = NComp + 1L + NWHAMFracAdd
    NComp = NComp + WNComp + NWHAMFracAdd
    WDonnanName = paste0("Donnan", WHAMFracAdd)
    WCompName = paste0(rep(WHAMprefix, each = WNComp / NWHAMFracAdd),
                       rep(
                         c(
                           MonodentTable$FullyProt,
                           BidentTable$FullyProt,
                           TridentTable$FullyProt
                         ),
                         times = NWHAMFracAdd
                       ))
    CompName = c(CompName, WDonnanName, WCompName)
    CompCharge = c(CompCharge,
                   array(0L, dim = NWHAMFracAdd, dimnames = list(WDonnanName)),
                   array(0L, dim = WNComp, dimnames = list(WCompName)))
    CompMC = c(CompMC,
               # array(NMass + match(WHAMFracAdd, c("HA", "FA")),
               array(DonnanMC[WHAMFracAdd],
                     dim = NWHAMFracAdd,
                     dimnames = list(WDonnanName)),
               array(iMass, dim = WNComp, dimnames = list(WCompName)))
    CompType = c(CompType,
                 array("FixedAct", NWHAMFracAdd, dimnames = list(WDonnanName)),
                 array("MassBal", WNComp, dimnames = list(WCompName)))
    CompActCorr = c(CompActCorr,
                    array("None", NWHAMFracAdd, dimnames = list(WDonnanName)),
                    array(rep(paste0("WHAM", WHAMFracAdd),
                              each = WNComp / NWHAMFracAdd),
                          dim = WNComp, dimnames = list(WCompName)))
    CompSiteDens = c(CompSiteDens,
                     array(1, NWHAMFracAdd, dimnames = list(WDonnanName)),
                     array(NA, WNComp, dimnames = list(WCompName)))

    StartDefComp = NDefComp + 1L + NWHAMFracAdd
    NDefComp = NDefComp + NWHAMFracAdd + WNComp
    DefCompName = c(DefCompName, WDonnanName, WCompName)
    DefCompFromNum = c(
      DefCompFromNum,
      array(NA, dim = NWHAMFracAdd, dimnames = list(WDonnanName)),
      array(NA, dim = WNComp, dimnames = list(WCompName)))
    DefCompFromVar = c(
      DefCompFromVar,
      array(NA, dim = NWHAMFracAdd, dimnames = list(WDonnanName)),
      array(rep(WHAMprefix, each = WNComp / NWHAMFracAdd),
            dim = WNComp, dimnames = list(WCompName)))
    DefCompCharge = c(
      DefCompCharge,
      array(0L, dim = NWHAMFracAdd, dimnames = list(WDonnanName)),
      array(0L, dim = WNComp, dimnames = list(WCompName)))
    DefCompMC = c(DefCompMC,
                  # array(NMass + match(WHAMFracAdd, c("HA", "FA")),
                  array(DonnanMC[WHAMFracAdd],
                        dim = NWHAMFracAdd,
                        dimnames = list(WDonnanName)),
                  array(iMass, dim = WNComp, dimnames = list(WCompName)))
    DefCompType = c(DefCompType,
                    array("FixedAct", dim = NWHAMFracAdd, dimnames = list(WDonnanName)),
                    array("MassBal", dim = WNComp, dimnames = list(WCompName)))
    DefCompActCorr = c(DefCompActCorr,
                       array(rep(paste0("WHAM", WHAMFracAdd),
                                 each = WNComp / NWHAMFracAdd),
                             dim = WNComp, dimnames = list(WCompName)))
    DefCompSiteDens = c(DefCompSiteDens,
                        array(1, dim = NWHAMFracAdd,
                              dimnames = list(WDonnanName)),
                        array(NA, dim = WNComp, dimnames = list(WCompName)))

    WNSpec = (nMS * (2L + nMP) +
                nBP * (4L + nMP) +
                nTG * (8L + nMP)) * NWHAMFracAdd
    WSpecName = array(paste0("newOCSpecies", 1:WNSpec), WNSpec)
    StartSpec = NSpec + + NWHAMFracAdd + 1L
    NSpec = NSpec + WNSpec + NWHAMFracAdd
    SpecName = c(SpecName, WDonnanName, WSpecName)
    SpecMC = c(SpecMC,
               # array(NMass + match(WHAMFracAdd, c("HA", "FA")),
               array(DonnanMC[WHAMFracAdd],
                     dim = NWHAMFracAdd,
                     dimnames = list(WDonnanName)),
               # this should always be water
               array(iMass, dim = WNSpec, dimnames = list(WSpecName)))
    SpecActCorr = c(SpecActCorr,
                    array("None", dim = NWHAMFracAdd,
                          dimnames = list(WDonnanName)),
                    array(rep(paste0("WHAM", WHAMFracAdd),
                              each = WNSpec / NWHAMFracAdd),
                          dim = WNSpec, dimnames = list(WSpecName)))
    SpecNC = c(SpecNC,
               array(1, dim = NWHAMFracAdd, dimnames = list(WDonnanName)),
               array(NA, dim = WNSpec, dimnames = list(WSpecName)))
    SpecCompList = rbind(SpecCompList,
                         matrix(0L, nrow = NWHAMFracAdd,
                                ncol = ncol(SpecCompList)),
                         matrix(0L, nrow = WNSpec, ncol = ncol(SpecCompList)))
    SpecStoich = rbind(
      cbind(
        SpecStoich,
        matrix(
          0L,
          nrow = NSpec - WNSpec - NWHAMFracAdd,
          ncol = NWHAMFracAdd + WNComp,
          dimnames = list(SpecName[1:(NSpec - WNSpec - NWHAMFracAdd)],
                          c(WDonnanName, WCompName))
        )
      ),
      matrix(
        0L,
        nrow = NWHAMFracAdd,
        ncol = NComp,
        dimnames = list(WDonnanName, CompName)
      ),
      matrix(
        0L,
        nrow = WNSpec,
        ncol = NComp,
        dimnames = list(WSpecName, CompName)
      )
    )
    SpecLogK = c(SpecLogK,
                 array(0.0, dim = NWHAMFracAdd, dimnames = list(WDonnanName)),
                 array(NA, dim = WNSpec, dimnames = list(WSpecName)))
    SpecDeltaH = c(SpecDeltaH,
                   array(0.0, dim = NWHAMFracAdd, dimnames = list(WDonnanName)),
                   array(0.0, dim = WNSpec, dimnames = list(WSpecName)))
    SpecTempKelvin = c(SpecTempKelvin,
                       array(298.15, dim = NWHAMFracAdd,
                             dimnames = list(WDonnanName)),
                       array(298.15, dim = WNSpec, dimnames = list(WSpecName)))
    SpecKsel = rbind(SpecKsel,
                     array(NA, dim = c(NWHAMFracAdd, 2), dimnames = list(WDonnanName)),
                     array(NA, dim = c(WNSpec, 2), dimnames = list(WSpecName)))

    MonodentpKH = numeric(nMS)
    MonodentAbundance = numeric(nMS)
    BidentAbundance = numeric(nBP)
    TridentAbundance = numeric(nTG)
    for (OMType in WHAMFracAdd) {

      ColspKM = paste0("pKM", c("A", "B"), OMType)

      # Monodentate sites
      NewCompNum = StartComp:(StartComp + nMS - 1)
      NewDefCompNum = StartDefComp:(StartDefComp + nMS - 1)
      MonodentpKH[1:4] = pKHA[OMType] + dpKHA[OMType] * (2 * MonodentTable$S[1:4] - 5) / 6
      MonodentpKH[5:8] = pKHB[OMType] + dpKHB[OMType] * (2 * MonodentTable$S[5:8] - 13) / 6
      MonodentAbundance = (1 - fprB[OMType] - fprT[OMType]) * nCOOH[OMType] / MonodentTable$AbundDenom
      CompSiteDens[NewCompNum] = MonodentAbundance * 1E-3 # the input is in milligrams, while nCOOH is mols/g
      DefCompSiteDens[NewDefCompNum] = MonodentAbundance * 1E-3 # the input is in milligrams, while nCOOH is mols/g

      # - fully protonated (components)
      NewSpecNum = StartSpec:(StartSpec + nMS - 1)
      SpecName[NewSpecNum] = paste0(WHAMprefix[OMType], MonodentTable$FullyProt)
      # SpecCharge[NewSpecNum] = 0L
      diag(SpecStoich[NewSpecNum, NewCompNum]) = 1L
      SpecLogK[NewSpecNum] = 0.0

      # - fully deprot
      NewSpecNum = NewSpecNum + nMS
      SpecName[NewSpecNum] = paste0(WHAMprefix[OMType], MonodentTable$FullyDeprot)
      # SpecCharge[NewSpecNum] = -1L
      diag(SpecStoich[NewSpecNum, NewCompNum]) = 1L
      SpecStoich[NewSpecNum, iH] = -1L
      SpecLogK[NewSpecNum] = -1 * MonodentpKH

      # bound to each metal
      for (iMetal in 1:nMP) {
        iMetalSpec = which(SpecName == MetalsTable$Metal[iMetal])#nolint: object_name_linter, line_length_linter.
        NewSpecNum = NewSpecNum + nMS
        # SpecCharge[NewSpecNum] = -1L + SpecCharge[iMetalSpec]
        SpecName[NewSpecNum] = paste0(WHAMprefix[OMType],
                                      MonodentTable$FullyDeprot,
                                      "-",
                                      MetalsTable$Metal[iMetal])
        SpecStoich[NewSpecNum, 1:NComp] = matrix(SpecStoich[iMetalSpec, ],
                                                 nrow = nMS,
                                                 ncol = NComp,
                                                 byrow = TRUE)
        diag(SpecStoich[NewSpecNum, NewCompNum]) = 1L
        SpecStoich[NewSpecNum, iH] = SpecStoich[NewSpecNum, iH] - 1L
        SpecLogK[NewSpecNum] = SpecLogK[iMetalSpec] -
          as.numeric(MetalsTable[iMetal, ColspKM[MonodentTable$Strong1Weak2]])
        SpecKsel[NewSpecNum, ] = SpecKsel[iMetalSpec, , drop = FALSE]
      }



      # Bidentate sites
      if (nBP > 0) {
        StartComp = max(NewCompNum) + 1
        StartDefComp = max(NewDefCompNum) + 1
        StartSpec = max(NewSpecNum) + 1
        NewCompNum = StartComp:(StartComp + nBP - 1)
        NewDefCompNum = StartDefComp:(StartDefComp + nBP - 1)
        BidentAbundance = fprB[OMType] * nCOOH[OMType] / BidentTable$AbundDenom
        CompSiteDens[NewCompNum] = BidentAbundance * 1E-3 # the input is in milligrams, while nCOOH is mols/g
        DefCompSiteDens[NewDefCompNum] = BidentAbundance * 1E-3 # the input is in milligrams, while nCOOH is mols/g

        # - fully protonated
        NewSpecNum = StartSpec:(StartSpec + nBP - 1)
        SpecName[NewSpecNum] = paste0(WHAMprefix[OMType], BidentTable$FullyProt)
        # SpecCharge[NewSpecNum] = 0L
        diag(SpecStoich[NewSpecNum, NewCompNum]) = 1L
        SpecLogK[NewSpecNum] = 0.0

        # - first site deprotonated
        NewSpecNum = NewSpecNum + nBP
        SpecName[NewSpecNum] = paste0(WHAMprefix[OMType], BidentTable$S1Deprot)
        # SpecCharge[NewSpecNum] = -1L
        diag(SpecStoich[NewSpecNum, NewCompNum]) = 1L
        SpecStoich[NewSpecNum, iH] = -1L
        SpecLogK[NewSpecNum] = -1 * MonodentpKH[BidentTable$S1]

        # - second site deprotonated
        NewSpecNum = NewSpecNum + nBP
        SpecName[NewSpecNum] = paste0(WHAMprefix[OMType], BidentTable$S2Deprot)
        # SpecCharge[NewSpecNum] = -1L
        diag(SpecStoich[NewSpecNum, NewCompNum]) = 1L
        SpecStoich[NewSpecNum, iH] = -1L
        SpecLogK[NewSpecNum] = -1 * MonodentpKH[BidentTable$S2]

        # - fully deprot
        NewSpecNum = NewSpecNum + nBP
        SpecName[NewSpecNum] = paste0(WHAMprefix[OMType], BidentTable$FullyDeprot)
        # SpecCharge[NewSpecNum] = -2L
        diag(SpecStoich[NewSpecNum, NewCompNum]) = 1L
        SpecStoich[NewSpecNum, iH] = -2L
        SpecLogK[NewSpecNum] = -1 * (MonodentpKH[BidentTable$S1] + MonodentpKH[BidentTable$S2])

        # bound to each metal
        for (iMetal in 1:nMP) {
          iMetalSpec = which(MetalsTable$Metal[iMetal] == SpecName)#nolint: object_name_linter, line_length_linter.
          NewSpecNum = NewSpecNum + nBP
          # SpecCharge[NewSpecNum] = -2L + SpecCharge[iMetalSpec]
          SpecName[NewSpecNum] = paste0(WHAMprefix[OMType],
                                        BidentTable$FullyDeprot,
                                        "-",
                                        MetalsTable$Metal[iMetal])
          SpecStoich[NewSpecNum, 1:NComp] = matrix(SpecStoich[iMetalSpec, ],
                                                   nrow = nBP,
                                                   ncol = NComp,
                                                   byrow = TRUE)
          diag(SpecStoich[NewSpecNum, NewCompNum]) = 1L
          SpecStoich[NewSpecNum, iH] = SpecStoich[NewSpecNum, iH] - 2L
          SpecLogK[NewSpecNum] = SpecLogK[iMetalSpec] -
            as.numeric(MetalsTable[iMetal, ColspKM[BidentTable$S1Strong1Weak2]] +
                         MetalsTable[iMetal, ColspKM[BidentTable$S2Strong1Weak2]])
          SpecKsel[NewSpecNum, ] = SpecKsel[iMetalSpec, , drop = FALSE]
        }

      }

      # Tridentate sites
      if (nTG > 0) {
        StartComp = max(NewCompNum) + 1
        StartDefComp = max(NewDefCompNum) + 1
        StartSpec = max(NewSpecNum) + 1
        NewCompNum = StartComp:(StartComp + nTG - 1)
        NewDefCompNum = StartDefComp:(StartDefComp + nTG - 1)
        TridentAbundance = fprT[OMType] * nCOOH[OMType] / TridentTable$AbundDenom
        CompSiteDens[NewCompNum] = TridentAbundance * 1E-3 # the input is in milligrams, while nCOOH is mols/g
        DefCompSiteDens[NewDefCompNum] = TridentAbundance * 1E-3 # the input is in milligrams, while nCOOH is mols/g

        # - fully protonated
        NewSpecNum = StartSpec:(StartSpec + nTG - 1)
        SpecName[NewSpecNum] = paste0(WHAMprefix[OMType], TridentTable$FullyProt)
        # SpecCharge[NewSpecNum] = 0L
        diag(SpecStoich[NewSpecNum, NewCompNum]) = 1L
        SpecLogK[NewSpecNum] = 0.0

        # - first site deprotonated
        NewSpecNum = NewSpecNum + nTG
        SpecName[NewSpecNum] = paste0(WHAMprefix[OMType], TridentTable$S1Deprot)
        # SpecCharge[NewSpecNum] = -1L
        diag(SpecStoich[NewSpecNum, NewCompNum]) = 1L
        SpecStoich[NewSpecNum, iH] = -1L
        SpecLogK[NewSpecNum] = -1 * MonodentpKH[TridentTable$S1]

        # - second site deprotonated
        NewSpecNum = NewSpecNum + nTG
        SpecName[NewSpecNum] = paste0(WHAMprefix[OMType], TridentTable$S2Deprot)
        # SpecCharge[NewSpecNum] = -1L
        diag(SpecStoich[NewSpecNum, NewCompNum]) = 1L
        SpecStoich[NewSpecNum, iH] = -1L
        SpecLogK[NewSpecNum] = -1 * MonodentpKH[TridentTable$S2]

        # - third site deprotonated
        NewSpecNum = NewSpecNum + nTG
        SpecName[NewSpecNum] = paste0(WHAMprefix[OMType], TridentTable$S3Deprot)
        # SpecCharge[NewSpecNum] = -1L
        diag(SpecStoich[NewSpecNum, NewCompNum]) = 1L
        SpecStoich[NewSpecNum, iH] = -1L
        SpecLogK[NewSpecNum] = -1 * MonodentpKH[TridentTable$S3]

        # - first & second sites deprotonated
        NewSpecNum = NewSpecNum + nTG
        SpecName[NewSpecNum] = paste0(WHAMprefix[OMType], TridentTable$S12Deprot)
        # SpecCharge[NewSpecNum] = -2L
        diag(SpecStoich[NewSpecNum, NewCompNum]) = 1L
        SpecStoich[NewSpecNum, iH] = -2L
        SpecLogK[NewSpecNum] = -1 * (MonodentpKH[TridentTable$S1] + MonodentpKH[TridentTable$S2])

        # - first & third sites deprotonated
        NewSpecNum = NewSpecNum + nTG
        SpecName[NewSpecNum] = paste0(WHAMprefix[OMType], TridentTable$S13Deprot)
        # SpecCharge[NewSpecNum] = -2L
        diag(SpecStoich[NewSpecNum, NewCompNum]) = 1L
        SpecStoich[NewSpecNum, iH] = -2L
        SpecLogK[NewSpecNum] = -1 * (MonodentpKH[TridentTable$S1] + MonodentpKH[TridentTable$S3])

        # - second & third sites deprotonated
        NewSpecNum = NewSpecNum + nTG
        SpecName[NewSpecNum] = paste0(WHAMprefix[OMType], TridentTable$S23Deprot)
        # SpecCharge[NewSpecNum] = -2L
        diag(SpecStoich[NewSpecNum, NewCompNum]) = 1L
        SpecStoich[NewSpecNum, iH] = -2L
        SpecLogK[NewSpecNum] = -1 * (MonodentpKH[TridentTable$S2] + MonodentpKH[TridentTable$S3])

        # - fully deprot
        NewSpecNum = NewSpecNum + nTG
        SpecName[NewSpecNum] = paste0(WHAMprefix[OMType], TridentTable$FullyDeprot)
        # SpecCharge[NewSpecNum] = -3L
        diag(SpecStoich[NewSpecNum, NewCompNum]) = 1L
        SpecStoich[NewSpecNum, iH] = -3L
        SpecLogK[NewSpecNum] = -1 * (MonodentpKH[TridentTable$S1] +
                                       MonodentpKH[TridentTable$S2] +
                                       MonodentpKH[TridentTable$S3])

        # bound to each metal
        for (iMetal in 1:nMP) {
          iMetalSpec = which(MetalsTable$Metal[iMetal] == SpecName) #nolint: object_name_linter, line_length_linter.
          NewSpecNum = NewSpecNum + nTG
          # SpecCharge[NewSpecNum] = -3L + SpecCharge[iMetalSpec]
          SpecName[NewSpecNum] = paste0(WHAMprefix[OMType],
                                        TridentTable$FullyDeprot,
                                        "-",
                                        MetalsTable$Metal[iMetal])
          SpecStoich[NewSpecNum, 1:NComp] = matrix(SpecStoich[iMetalSpec, ],
                                                   nrow = nTG,
                                                   ncol = NComp,
                                                   byrow = TRUE)
          diag(SpecStoich[NewSpecNum, NewCompNum]) = 1L
          SpecStoich[NewSpecNum, iH] = SpecStoich[NewSpecNum, iH] - 3L
          SpecLogK[NewSpecNum] = SpecLogK[iMetalSpec] -
            as.numeric(MetalsTable[iMetal, ColspKM[TridentTable$S1Strong1Weak2]] +
                         MetalsTable[iMetal, ColspKM[TridentTable$S2Strong1Weak2]] +
                         MetalsTable[iMetal, ColspKM[TridentTable$S3Strong1Weak2]])
          SpecKsel[NewSpecNum, ] = SpecKsel[iMetalSpec, , drop = FALSE]
        }

      }

      StartComp = max(NewCompNum) + 1
      StartDefComp = max(NewDefCompNum) + 1
      StartSpec = max(NewSpecNum) + 1

    }
  }

  # Cleanup
  names(CompSiteDens) = CompName
  names(DefCompSiteDens) = DefCompName
  # names(SpecCharge) = SpecName
  names(SpecMC) = SpecName
  # names(SpecType) = SpecName
  names(SpecActCorr) = SpecName
  rownames(SpecStoich) = SpecName
  names(SpecLogK) = SpecName
  names(SpecDeltaH) = SpecName
  names(SpecTempKelvin) = SpecName
  rownames(SpecKsel) = SpecName

  # Re-ordering species so components are in front
  Reorder = match(c(CompName, SpecName[SpecName %in% CompName == FALSE]), SpecName)
  SpecName = SpecName[Reorder]
  SpecMC = SpecMC[Reorder]
  SpecActCorr = SpecActCorr[Reorder]
  SpecNC = SpecNC[Reorder]
  SpecCompList = SpecCompList[Reorder, ]
  SpecStoich = SpecStoich[Reorder, ]
  SpecLogK = SpecLogK[Reorder]
  SpecDeltaH = SpecDeltaH[Reorder]
  SpecTempKelvin = SpecTempKelvin[Reorder]
  SpecKsel = SpecKsel[Reorder, ]

  SpecNC = rowSums(SpecStoich != 0L)
  names(SpecNC) = SpecName
  SpecCompList = t(apply(
    SpecStoich,
    MARGIN = 1,
    FUN = function(X) {
      Tmp = sort(which(X != 0L))
      if (length(Tmp) < max(SpecNC)) {
        Tmp = c(Tmp, rep(0, max(SpecNC) - length(Tmp)))
      }
      return(Tmp)
    }
  ))
  rownames(SpecCompList) = SpecName
  # CompNS = colSums(Stoich != 0L)
  # names(CompNS) = CompName
  # SpecList = t(apply(
  #   Stoich,
  #   MARGIN = 2,
  #   FUN = function(X) {
  #     Tmp = sort(which(X != 0L))
  #     if (length(Tmp) < max(CompNS)){
  #       Tmp = c(Tmp, rep(0,max(CompNS)-length(Tmp)))
  #     }
  #     return(Tmp)
  #   }))
  # rownames(SpecList) = CompName


  # Assemble output list
  return(list(

    # Components
    NComp = NComp,
    CompName = CompName,
    CompCharge = CompCharge,
    CompMC = CompMC,
    CompType = CompType,
    CompActCorr = CompActCorr,
    CompSiteDens = CompSiteDens,

    # Defined Components
    NDefComp = NDefComp,
    DefCompName = DefCompName,
    DefCompCharge = DefCompCharge,
    DefCompFromNum = DefCompFromNum,
    DefCompFromVar = DefCompFromVar,
    DefCompMC = DefCompMC,
    DefCompType = DefCompType,
    DefCompActCorr = DefCompActCorr,
    DefCompSiteDens = DefCompSiteDens,

    # Formation Reactions
    NSpec = NSpec,
    SpecName = SpecName,
    SpecMC = SpecMC,
    SpecActCorr = SpecActCorr,
    SpecNC = SpecNC,
    SpecCompList = SpecCompList,
    SpecStoich = SpecStoich,
    SpecLogK = SpecLogK,
    SpecDeltaH = SpecDeltaH,
    SpecTempKelvin = SpecTempKelvin,

    # Phase List
    PhaseCompList = PhaseCompList,
    PhaseStoich = PhaseStoich,

    # WHAM parameters - to be used later
    SpecKsel = SpecKsel,
    DonnanMC = DonnanMC,
    wDLF = wDLF,
    wKZED = wKZED,
    wP = wP,
    wRadius = wRadius,
    wMolWt = wMolWt

  ))

}
