#' @title Expand the DOC component into WHAM components
#'
#' @references
#'   Tipping E. (1994). WHAM--A chemical equilibrium model and
#'     computer code for waters, sediments, and soils incorporating a discrete
#'     site/electrostatic model of ion-binding by humic substances. Computers &
#'     Geosciences, vol. 20, iss. 6, pp. 973-1023.
#'
#' @param NMass integer; Number of mass compartments
#' @param MassName character vector of length `NMass`; Names of the mass
#'   compartments
#' @param NInVar integer; Number of input variables
#' @param InVarName character vector of length `NInVar`; Names of input
#'   variables
#' @param InVarMCR integer vector of length `NInVar`;  Mass compartments of input
#'   variables
#' @param InVarType character vector of length `NInVar`; Types of input
#'   variables
#' @param NComp integer; Number of components (modified and returned)
#' @param CompName character vector of length `NComp`; component names (modified
#'   and returned)
#' @param CompCharge integer vector of length `NComp`; the charge of the
#'   components as free ions (modified and returned)
#' @param CompMCR integer vector of length `NComp`; Which mass compartment the
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
#' @param DefCompMCR integer vector of length `NDefComp`; Which mass compartment
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
#' @param SpecMCR integer vector of length `NSpec`; which mass compartment the
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
                      MassAmt,
                      MassUnit,
                      NInVar,
                      InVarName,
                      InVarMCR,
                      InVarType,
                      NComp,
                      CompName,
                      CompCharge,
                      CompMCName,
                      CompMCR,
                      CompType,
                      CompActCorr,
                      CompSiteDens,
                      NDefComp,
                      DefCompName,
                      DefCompFromNum,
                      DefCompFromVar,
                      DefCompCharge,
                      DefCompMCName,
                      DefCompMCR,
                      DefCompType,
                      DefCompActCorr,
                      DefCompSiteDens,
                      NSpec,
                      SpecName,
                      SpecType,
                      SpecMCName,
                      SpecMCR,
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
      WdatFile = system.file("extdata/WHAM/WHAM_V.wdat",
                             package = "BLMEngineInR",
                             mustWork = TRUE)
    } else if (WHAMVer == "VI") {
      WdatFile = system.file("extdata/WHAM/WHAM_VI.wdat",
                             package = "BLMEngineInR",
                             mustWork = TRUE)
    } else if (WHAMVer  == "VII") {
      WdatFile = system.file("extdata/WHAM/WHAM_VII.wdat",
                             package = "BLMEngineInR",
                             mustWork = TRUE)
    }
  } else {
    WdatFile = normalizePath(WdatFile)
    stopifnot(file.exists(WdatFile))
  }


  # read WHAM data file-------------------------------------
  # header info
  SkipRows = 2L
  Tmp = read.delim(
    file = WdatFile,
    header = FALSE,
    sep = ",",
    skip = SkipRows,
    nrows = 7L
  )
  nMS = as.integer(Tmp[1, 2])#Number of monodentate sites#nolint: object_name_linter, line_length_linter.
  nBP = as.integer(Tmp[2, 2])#Number of bidentate pairs#nolint: object_name_linter, line_length_linter.
  nTG = as.integer(Tmp[3, 2])#Number of tridentate groups#nolint: object_name_linter, line_length_linter.
  nMP_file = as.integer(Tmp[4, 2])#Number of metals-OM parameters#nolint: object_name_linter, line_length_linter.
  nKsel_file = as.integer(Tmp[5, 2])#Number of non-standard selectivity coefficients#nolint: object_name_linter, line_length_linter.
  wDLF = as.numeric(Tmp[6, 2])#Double layer overlap factor#nolint: object_name_linter, line_length_linter.
  wKZED = as.numeric(Tmp[7, 2])#Constant to control DDL at low ZED#nolint: object_name_linter, line_length_linter.

  # Parameters
  SkipRows = SkipRows + 7L + 1L
  Tmp = read.delim(
    file = WdatFile,
    header = TRUE,
    sep = ",",
    skip = SkipRows,
    nrows = 12L
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
  SkipRows = SkipRows + 12L + 3L
  MonodentTable = read.delim(
    file = WdatFile,
    header = TRUE,
    sep = ",",
    skip = SkipRows,
    nrows = nMS
  )
  names(MonodentTable) = c("S","AbundDenom")
  MonodentTable$FullyProt = paste0(MonodentTable$S, "H")
  MonodentTable$FullyDeprot = paste0(MonodentTable$S)
  nStrong = nMS / 2
  MonodentTable$Strong1Weak2 = rep(c(1L, 2L), each = nStrong)#c(rep(1L, 4), rep(2L, 4))

  # Bidentate Pairs
  SkipRows = SkipRows + nMS + 3L
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
    BidentTable$S1Strong1Weak2 = ifelse(BidentTable$S1 <= nStrong, 1L, 2L)
    BidentTable$S2Strong1Weak2 = ifelse(BidentTable$S2 <= nStrong, 1L, 2L)
  } else {
    BidentTable = data.frame()
  }


  # Tridentate Groups
  SkipRows = SkipRows + nBP + 3L
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
    TridentTable$S1Strong1Weak2 = ifelse(TridentTable$S1 <= nStrong, 1L, 2L)
    TridentTable$S2Strong1Weak2 = ifelse(TridentTable$S2 <= nStrong, 1L, 2L)
    TridentTable$S3Strong1Weak2 = ifelse(TridentTable$S3 <= nStrong, 1L, 2L)
  } else {
    TridentTable = data.frame()
  }

  # Metals Parameters Table
  SkipRows = SkipRows + nTG + 3L
  if (nMP_file > 0) {
    MetalsTable = read.delim(
      file = WdatFile,
      header = TRUE,
      sep = ",",
      skip = SkipRows,
      nrows = nMP_file
    )
    names(MetalsTable) = c("Metal", "pKMAHA", "pKMAFA", "dLK2")
    MetalsTable = MetalsTable[MetalsTable$Metal %in% c(CompName, SpecName), ]
    nMP = nrow(MetalsTable) # nolint: object_name_linter.
    MetalsTable$pKMBHA = 3 * MetalsTable$pKMAHA - 3
    MetalsTable$pKMBFA = 3.96 * MetalsTable$pKMAFA
  } else {
    MetalsTable = data.frame()
    nMP = 0L
  }

  # Non-standard selectivity coefficients
  SkipRows = SkipRows + nMP_file + 3L
  if (nKsel_file > 0) {
    SpecKselTable = read.delim(
      file = WdatFile,
      header = TRUE,
      sep = ",",
      skip = SkipRows,
      nrows = nKsel_file
    )
    names(SpecKselTable) = c("Spec", "KselHA", "KselFA")
    SpecKselTable = SpecKselTable[SpecKselTable$Spec %in% SpecName, ]
    nKsel = nrow(SpecKselTable)
  } else {
    SpecKselTable = data.frame()
    nKsel = 0L
  }




  # Save original copies of arrays -------------------------
  NMass_orig = NMass
  MassName_orig = MassName
  MassAmt_orig = MassAmt
  MassUnit_orig = MassUnit

  NComp_orig = NComp
  CompName_orig = CompName
  CompCharge_orig = CompCharge
  CompMCName_orig = CompMCName
  CompMCR_orig = CompMCR
  CompType_orig = CompType
  CompActCorr_orig = CompActCorr
  CompSiteDens_orig = CompSiteDens

  NDefComp_orig = NDefComp
  DefCompName_orig = DefCompName
  DefCompFromNum_orig = DefCompFromNum
  DefCompFromVar_orig = DefCompFromVar
  DefCompCharge_orig = DefCompCharge
  DefCompMCName_orig = DefCompMCName
  DefCompMCR_orig = DefCompMCR
  DefCompType_orig = DefCompType
  DefCompActCorr_orig = DefCompActCorr
  DefCompSiteDens_orig = DefCompSiteDens

  NSpec_orig = NSpec
  SpecName_orig = SpecName
  SpecType_orig = SpecType
  SpecMCName_orig = SpecMCName
  SpecMCR_orig = SpecMCR
  SpecActCorr_orig = SpecActCorr
  SpecCharge_orig = SpecStoich %*% CompCharge
  SpecStoich_orig = SpecStoich
  SpecLogK_orig = SpecLogK
  SpecDeltaH_orig = SpecDeltaH
  SpecTempKelvin_orig = SpecTempKelvin


  # Do the expansion ---------------------------------------

  # Initialize variables
  iH = which(CompName == "H") # nolint: object_name_linter.
  SpecCharge = SpecStoich %*% CompCharge

  # Figure out the number of DOC components we're adding, and what fraction
  InVarWHAM = which(grepl("WHAM", InVarType))

  for (iInVar in InVarWHAM) {

    iMass = InVarMCR[iInVar] # nolint: object_name_linter.
    ChargedSpecName = SpecName_orig[(SpecCharge_orig != 0) &
                                      (SpecMCR_orig == iMass)]
    NChargedSpec = length(ChargedSpecName)
    SpecKsel = array(1, dim = c(NChargedSpec, 2),
                     dimnames = list(ChargedSpecName, c("HA", "FA")))
    if (nKsel > 0L) {
      SpecKsel[match(SpecKselTable$Spec, ChargedSpecName), ]  =
        array(unlist(SpecKselTable[, c("KselHA", "KselFA")]),
              dim = c(nKsel, 2))
    }
    ChargedSpecDonnanLogK = log10(SpecKsel) +
      SpecLogK_orig[match(ChargedSpecName, SpecName_orig)]

    if (InVarType[iInVar] == "WHAM-HA") {
      WHAMFracAdd = c("HA")
    } else if (InVarType[iInVar] == "WHAM-FA") {
      WHAMFracAdd = c("FA")
    } else if (InVarType[iInVar] == "WHAM-HAFA") {
      WHAMFracAdd = c("HA", "FA")
      if (!any(InVarType[InVarMCR == InVarMCR[iInVar]] %in% "PercHA")) {
        stop("Must have PercHA input variable in mass compartment if specifying a WHAM-HAFA input variable.") # nolint: line_length_linter.
      }
    }
    if ((InVarType[iInVar] %in% c("WHAM-FA", "WHAM-HA")) &&
          any(InVarType[InVarMCR == InVarMCR[iInVar]] %in% "PercHA")) {
      stop("PercHA input variable specified in mass compartment with WHAM-HA or WHAM-FA input variable.") # nolint: line_length_linter.
    }
    if ((InVarType[iInVar] %in% c("WHAM-HA")) &&
        any(InVarType[InVarMCR == InVarMCR[iInVar]] %in% "PercAFA")) {
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
    NDonnanComp = NWHAMFracAdd
    DonnanCompName = paste0("Donnan", WHAMFracAdd)
    DonnanMCName = paste(MassName_orig[iMass], DonnanCompName, sep = "_")
    DonnanMCR = array(NMass_orig + (1:NDonnanComp), dim = NDonnanComp,
                     dimnames = list(WHAMFracAdd))

    # BulkMassName = paste0(MassName_orig[iMass], "_Bulk")
    # BulkMCR = NMass_orig + NDonnanComp + 1

    NMass = NMass_orig + NDonnanComp #+ 1
    # MassName = c(MassName_orig, DonnanMCName, BulkMassName)
    MassName = c(MassName_orig, DonnanMCName)
    MassAmt = c(MassAmt_orig,
                array(1E-5, dim = NDonnanComp, dimnames = list(DonnanMCName)))
                # array(c(1E-5, MassAmt_orig[iMass]),
                #         dim = NDonnanComp + 1,
                #         dimnames = list(c(DonnanMCName, BulkMassName))))
    MassUnit = c(MassUnit_orig,
                 array(MassUnit_orig[iMass], dim = NDonnanComp,
                       dimnames = list(DonnanMCName)))
                 # array(MassUnit_orig[iMass],
                 #       dim = NDonnanComp + 1,
                 #       dimnames = list(c(DonnanMCName, BulkMassName))))

    NWHAMComp = (nMS + nBP + nTG) * NWHAMFracAdd
    StartComp = NComp_orig + 1L + NDonnanComp
    NComp = NComp_orig + NWHAMComp + NDonnanComp
    WHAMCompName = paste0(rep(WHAMprefix, each = NWHAMComp / NWHAMFracAdd),
                       rep(
                         c(
                           MonodentTable$FullyProt,
                           BidentTable$FullyProt,
                           TridentTable$FullyProt
                         ),
                         times = NWHAMFracAdd
                       ))
    CompName = c(CompName_orig, DonnanCompName, WHAMCompName)
    CompCharge = c(CompCharge_orig,
                   array(0L, dim = NDonnanComp,
                         dimnames = list(DonnanCompName)),
                   array(0L, dim = NWHAMComp, dimnames = list(WHAMCompName)))
    CompMCName = c(CompMCName_orig,
                   array(DonnanMCName,
                         dim = NDonnanComp,
                         dimnames = list(DonnanCompName)))
    CompMCR = c(CompMCR_orig,
               # array(NMass + match(WHAMFracAdd, c("HA", "FA")),
               array(DonnanMCR,#iMass, #DonnanMCR[WHAMFracAdd],
                     dim = NDonnanComp,
                     dimnames = list(DonnanCompName)),
               # array(BulkMCR, dim = NWHAMComp, dimnames = list(WHAMCompName)))
               array(iMass, dim = NWHAMComp, dimnames = list(WHAMCompName)))
    CompType = c(CompType_orig,
                 array(DonnanCompName, NDonnanComp,
                       dimnames = list(DonnanCompName)),
                 array("MassBal", NWHAMComp, dimnames = list(WHAMCompName)))
    CompActCorr = c(CompActCorr_orig,
                    array("None", dim = NDonnanComp,
                          dimnames = list(DonnanCompName)),
                    array("None", dim = NWHAMComp,
                          dimnames = list(WHAMCompName)))
    CompSiteDens = c(CompSiteDens_orig,
                     array(1, NDonnanComp, dimnames = list(DonnanCompName)),
                     array(NA, NWHAMComp, dimnames = list(WHAMCompName)))

    StartDefComp = NDefComp_orig + 1L + NDonnanComp
    NDefComp = NDefComp_orig + NDonnanComp + NWHAMComp
    DefCompName = c(DefCompName_orig, DonnanCompName, WHAMCompName)
    DefCompFromNum = c(
      DefCompFromNum_orig,
      array(NA, dim = NDonnanComp, dimnames = list(DonnanCompName)),
      array(NA, dim = NWHAMComp, dimnames = list(WHAMCompName)))
    DefCompFromVar = c(
      DefCompFromVar_orig,
      array(WHAMprefix, dim = NDonnanComp, dimnames = list(DonnanCompName)),
      array(rep(WHAMprefix, each = NWHAMComp / NWHAMFracAdd),
            dim = NWHAMComp, dimnames = list(WHAMCompName)))
    DefCompCharge = c(
      DefCompCharge_orig,
      array(0L, dim = NDonnanComp, dimnames = list(DonnanCompName)),
      array(0L, dim = NWHAMComp, dimnames = list(WHAMCompName)))
    DefCompMCName = c(DefCompMCName_orig,
                      array(DonnanMCName,
                            dim = NDonnanComp,
                            dimnames = list(DonnanCompName)),
                      array(MassName[iMass], dim = NWHAMComp,
                            dimnames = list(WHAMCompName)))
    DefCompMCR = c(DefCompMCR_orig,
                  # array(NMass + match(WHAMFracAdd, c("HA", "FA")),
                  array(DonnanMCR, #iMass, #DonnanMCR[WHAMFracAdd],
                        dim = NDonnanComp,
                        dimnames = list(DonnanCompName)),
                  # array(BulkMCR, dim = NWHAMComp, dimnames = list(WHAMCompName)))
                  array(iMass, dim = NWHAMComp, dimnames = list(WHAMCompName)))
    DefCompType = c(DefCompType_orig,
                    array(DonnanCompName, dim = NDonnanComp,
                          dimnames = list(DonnanCompName)),
                    array("MassBal", dim = NWHAMComp,
                          dimnames = list(WHAMCompName)))
    DefCompActCorr = c(DefCompActCorr_orig,
                       array("None", dim = NDonnanComp,
                             dimnames = list(DonnanCompName)),
                       array("None", dim = NWHAMComp,
                             dimnames = list(WHAMCompName)))
    DefCompSiteDens = c(DefCompSiteDens_orig,
                        array(1.0E-4, dim = NDonnanComp,
                              dimnames = list(DonnanCompName)),
                        array(NA, dim = NWHAMComp,
                              dimnames = list(WHAMCompName)))

    NWHAMSpec = (nMS * (2L + nMP) +
                   nBP * (4L + nMP) +
                   nTG * (8L + nMP)) * NWHAMFracAdd
    WHAMSpecName = array(paste0("newOCSpecies", 1:NWHAMSpec), NWHAMSpec)
    NDonnanSpec = (NChargedSpec + 1L) * NWHAMFracAdd
    DonnanSpecName = c(DonnanCompName,
                        paste0(rep(DonnanCompName, each = NChargedSpec), "-",
                               rep(ChargedSpecName, times = NWHAMFracAdd)))
    StartSpec = NSpec_orig + NDonnanSpec + 1L
    NSpec = NSpec_orig + NWHAMSpec + NDonnanSpec
    SpecName = c(SpecName_orig, DonnanSpecName, WHAMSpecName)
    SpecMCName = c(SpecMCName_orig,
                   array(c(DonnanMCName,
                           rep(DonnanMCName, each = NChargedSpec)),
                         dim = NDonnanSpec,
                         dimnames = list(DonnanSpecName)),
                   array(MassName[iMass], dim = NWHAMSpec,
                         dimnames = list(WHAMSpecName)))
    SpecMCR = c(SpecMCR_orig,
               array(c(DonnanMCR,
                       rep(DonnanMCR, each = NChargedSpec)), #iMass, #DonnanMCR[WHAMFracAdd],
                     dim = NDonnanSpec,
                     dimnames = list(DonnanSpecName)),
               # this should always be water
               # array(BulkMCR, dim = NWHAMSpec, dimnames = list(WHAMSpecName)))
               array(iMass, dim = NWHAMSpec, dimnames = list(WHAMSpecName)))
    SpecType = c(SpecType_orig,
                 array(c(DonnanCompName,
                         rep(DonnanCompName, each = NChargedSpec)),
                       dim = NDonnanSpec,
                       dimnames = list(DonnanSpecName)),
                 array(rep(paste0("WHAM", WHAMFracAdd),
                           each = NWHAMSpec / NWHAMFracAdd),
                       dim = NWHAMSpec, dimnames = list(WHAMSpecName)))
    SpecActCorr = c(SpecActCorr_orig,
                    array(c(rep("None", each = NWHAMFracAdd),
                            rep(SpecActCorr_orig[match(ChargedSpecName, SpecName_orig)],
                                times = NWHAMFracAdd)),
                          dim = NDonnanSpec,
                          dimnames = list(DonnanSpecName)),
                    array("None", dim = NWHAMSpec, dimnames = list(WHAMSpecName)))

    SpecStoich = rbind(
      cbind(
        SpecStoich_orig,
        matrix(
          0L,
          nrow = NSpec_orig,
          ncol = NWHAMFracAdd + NWHAMComp,
          dimnames = list(SpecName[1:NSpec_orig],
                          c(DonnanCompName, WHAMCompName))
        )
      ),
      matrix(
        0L,
        nrow = NDonnanSpec,
        ncol = NComp,
        dimnames = list(DonnanSpecName, CompName)
      ),
      matrix(
        0L,
        nrow = NWHAMSpec,
        ncol = NComp,
        dimnames = list(WHAMSpecName, CompName)
      )
    )
    SpecLogK = c(SpecLogK_orig,
                 array(NA, dim = NDonnanSpec, dimnames = list(DonnanSpecName)),
                 array(NA, dim = NWHAMSpec, dimnames = list(WHAMSpecName)))
    SpecDeltaH = c(SpecDeltaH_orig,
                   array(0.0, dim = NDonnanSpec,
                         dimnames = list(DonnanSpecName)),
                   array(0.0, dim = NWHAMSpec, dimnames = list(WHAMSpecName)))
    SpecTempKelvin = c(SpecTempKelvin_orig,
                       array(298.15, dim = NDonnanSpec,
                             dimnames = list(DonnanSpecName)),
                       array(298.15, dim = NWHAMSpec,
                             dimnames = list(WHAMSpecName)))

    MonodentpKH = numeric(nMS)
    MonodentAbundance = numeric(nMS)
    BidentAbundance = numeric(nBP)
    TridentAbundance = numeric(nTG)
    for (OMType in WHAMFracAdd) {

      # Donnan diffuse binding
      DonnanOMChargedSpecName = paste0("Donnan", OMType, "-", ChargedSpecName)
      SpecLogK[paste0("Donnan", OMType)] = 0.0
      SpecLogK[DonnanOMChargedSpecName] =
        ChargedSpecDonnanLogK[ChargedSpecName, OMType]
      SpecDeltaH[DonnanOMChargedSpecName] =
        SpecDeltaH_orig[match(ChargedSpecName, SpecName_orig)]
      SpecTempKelvin[DonnanOMChargedSpecName] =
        SpecTempKelvin_orig[match(ChargedSpecName, SpecName_orig)]
      SpecStoich[DonnanOMChargedSpecName, CompName_orig] =
        SpecStoich_orig[match(ChargedSpecName, SpecName_orig), CompName_orig]
      SpecStoich[DonnanOMChargedSpecName, paste0("Donnan", OMType)] =
        abs(SpecCharge_orig[match(ChargedSpecName, SpecName_orig)])


      ColspKM = paste0("pKM", c("A", "B"), OMType)

      # Monodentate sites
      NewCompNum = StartComp:(StartComp + nMS - 1)
      NewDefCompNum = StartDefComp:(StartDefComp + nMS - 1)
      MonodentpKH[1:nStrong] = pKHA[OMType] + dpKHA[OMType] *
        (2 * MonodentTable$S[1:nStrong] - 5) / 6
      MonodentpKH[(nStrong + 1):nMS] = pKHB[OMType] + dpKHB[OMType] *
        (2 * MonodentTable$S[(nStrong + 1):nMS] - 13) / 6
      MonodentAbundance = (1 - fprB[OMType] - fprT[OMType]) *
        nCOOH[OMType] / MonodentTable$AbundDenom
      CompSiteDens[NewCompNum] = MonodentAbundance * 2E-3 # the input is in mg C/L, while nCOOH is mols/g HS
      DefCompSiteDens[NewDefCompNum] = MonodentAbundance * 2E-3 # the input is in mg C/L, while nCOOH is mols/g HS

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
        SpecDeltaH[NewSpecNum] = SpecDeltaH[iMetalSpec]
        SpecTempKelvin[NewSpecNum] = SpecTempKelvin[iMetalSpec]
      }



      # Bidentate sites
      if (nBP > 0) {
        StartComp = max(NewCompNum) + 1
        StartDefComp = max(NewDefCompNum) + 1
        StartSpec = max(NewSpecNum) + 1
        NewCompNum = StartComp:(StartComp + nBP - 1)
        NewDefCompNum = StartDefComp:(StartDefComp + nBP - 1)
        BidentAbundance = fprB[OMType] * nCOOH[OMType] / BidentTable$AbundDenom
        CompSiteDens[NewCompNum] = BidentAbundance * 2E-3 # the input is in mg C/L, while nCOOH is mols/g HS
        DefCompSiteDens[NewDefCompNum] = BidentAbundance * 2E-3 # the input is in mg C/L, while nCOOH is mols/g HS

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
          SpecDeltaH[NewSpecNum] = SpecDeltaH[iMetalSpec]
          SpecTempKelvin[NewSpecNum] = SpecTempKelvin[iMetalSpec]
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
        CompSiteDens[NewCompNum] = TridentAbundance * 2E-3 # the input is in mg C/L, while nCOOH is mols/g HS
        DefCompSiteDens[NewDefCompNum] = TridentAbundance * 2E-3 # the input is in mg C/L, while nCOOH is mols/g HS

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
          SpecDeltaH[NewSpecNum] = SpecDeltaH[iMetalSpec]
          SpecTempKelvin[NewSpecNum] = SpecTempKelvin[iMetalSpec]
        }

      }

      StartComp = max(NewCompNum) + 1
      StartDefComp = max(NewDefCompNum) + 1
      StartSpec = max(NewSpecNum) + 1

    }
  }

  # Cleanup
  names(MassAmt) = MassName
  names(MassUnit) = MassName
  names(CompSiteDens) = CompName
  names(CompType) = CompName
  names(CompActCorr) = CompName
  names(DefCompSiteDens) = DefCompName
  names(SpecMCName) = SpecName
  names(SpecType) = SpecName
  names(SpecMCR) = SpecName
  names(SpecActCorr) = SpecName
  rownames(SpecStoich) = SpecName
  names(SpecLogK) = SpecName
  names(SpecDeltaH) = SpecName
  names(SpecTempKelvin) = SpecName
  WHAMDonnanMCR = array(DonnanMCR[c("HA","FA")], dim = 2,
                       dimnames = list(c("HA","FA")))

  # Re-ordering species so components are in front
  Reorder = match(c(CompName, SpecName[SpecName %in% CompName == FALSE]),
                  SpecName)
  SpecName = SpecName[Reorder]
  SpecMCName = SpecMCName[Reorder]
  SpecMCR = SpecMCR[Reorder]
  SpecActCorr = SpecActCorr[Reorder]
  SpecStoich = SpecStoich[Reorder, ]
  SpecLogK = SpecLogK[Reorder]
  SpecDeltaH = SpecDeltaH[Reorder]
  SpecTempKelvin = SpecTempKelvin[Reorder]

  SpecNC = as.integer(rowSums(SpecStoich != 0L))
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

    # Mass Compartment List
    NMass = NMass,
    MassName = MassName,
    MassAmt = MassAmt,
    MassUnit = MassUnit,
    WHAMDonnanMCR = WHAMDonnanMCR,

    # Components
    NComp = NComp,
    CompName = CompName,
    CompCharge = CompCharge,
    CompMCName = CompMCName,
    CompMCR = CompMCR,
    CompType = CompType,
    CompActCorr = CompActCorr,
    CompSiteDens = CompSiteDens,

    # Defined Components
    NDefComp = NDefComp,
    DefCompName = DefCompName,
    DefCompCharge = DefCompCharge,
    DefCompFromNum = DefCompFromNum,
    DefCompFromVar = DefCompFromVar,
    DefCompMCName = DefCompMCName,
    DefCompMCR = DefCompMCR,
    DefCompType = DefCompType,
    DefCompActCorr = DefCompActCorr,
    DefCompSiteDens = DefCompSiteDens,

    # Formation Reactions
    NSpec = NSpec,
    SpecName = SpecName,
    SpecType = SpecType,
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
    PhaseCompList = PhaseCompList,
    PhaseStoich = PhaseStoich,

    # WHAM parameters - to be used later
    wDLF = wDLF,
    wKZED = wKZED,
    wP = wP,
    wRadius = wRadius,
    wMolWt = wMolWt

  ))

}
