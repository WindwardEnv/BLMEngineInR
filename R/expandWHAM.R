#' @title Expand the DOC component into WHAM components
#'
#' @references
#'   Tipping E. (1994). WHAM--A chemical equilibrium model and
#'     computer code for waters, sediments, and soils incorporating a discrete
#'     site/electrostatic model of ion-binding by humic substances. Computers &
#'     Geosciences, vol. 20, iss. 6, pp. 973-1023.
#'
#' @param ThisProblem a list object following the template of BlankProblem
#' @param WHAMVer a character string specifying the WHAM version to use, must be
#'   one of `"V"` (default), `"VI"`, or `"VII"`. Ignored if `WdatFile` is not
#'   `NULL`.
#' @param WdatFile (optional) a character string specifying the file path of a
#'   WHAM parameter file
#'
#' @keywords internal
#'
#' @noRd
ExpandWHAM = function(ThisProblem,
                      WHAMVer = c("V", "VI", "VII"),
                      WdatFile = NA) {

  # error catching and input cleanup
  if (is.na(WdatFile)) {
    WHAMVer = match.arg(WHAMVer, choices =  c("V", "VI", "VII"))
    if (WHAMVer == "V") {
      WdatFile = system.file(file.path("extdata","WHAM","WHAM_V.wdat"),
                             package = "BLMEngineInR",
                             mustWork = TRUE)
    } else if (WHAMVer == "VI") {
      WdatFile = system.file(file.path("extdata","WHAM","WHAM_VI.wdat"),
                             package = "BLMEngineInR",
                             mustWork = TRUE)
    } else if (WHAMVer  == "VII") {
      WdatFile = system.file(file.path("extdata","WHAM","WHAM_VII.wdat"),
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
    MetalsTable = MetalsTable[MetalsTable$Metal %in% c(ThisProblem$Comp$Name, ThisProblem$Spec$Name), ]
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
    SpecKselTable = SpecKselTable[SpecKselTable$Spec %in% ThisProblem$Spec$Name, ]
    nKsel = nrow(SpecKselTable)
  } else {
    SpecKselTable = data.frame()
    nKsel = 0L
  }




  # Save original copies of arrays -------------------------
  NewProblem = ThisProblem

  # Do the expansion ---------------------------------------

  # Initialize variables
  iH = which(ThisProblem$Comp$Name == "H") # nolint: object_name_linter.

  # Figure out the number of DOC components we're adding, and what fraction
  InVarWHAM = which(grepl("WHAM", ThisProblem$InVar$Type))

  for (iInVar in InVarWHAM) {

    iMass = ThisProblem$InVar$MCR[iInVar] # nolint: object_name_linter.
    ChargedSpecName = ThisProblem$Spec$Name[(ThisProblem$Spec$Charge != 0) &
                                              (ThisProblem$Spec$MCR == iMass)]
    NChargedSpec = length(ChargedSpecName)
    SpecKsel = array(1, dim = c(NChargedSpec, 2),
                     dimnames = list(ChargedSpecName, c("HA", "FA")))
    if (nKsel > 0L) {
      SpecKsel[match(SpecKselTable$Spec, ChargedSpecName), ]  =
        array(unlist(SpecKselTable[, c("KselHA", "KselFA")]),
              dim = c(nKsel, 2))
    }
    ChargedSpecDonnanLogK = log10(SpecKsel) +
      ThisProblem$Spec$LogK[match(ChargedSpecName, ThisProblem$Spec$Name)]

    if (ThisProblem$InVar$Type[iInVar] == "WHAM-HA") {
      WHAMFracAdd = c("HA")
    } else if (ThisProblem$InVar$Type[iInVar] == "WHAM-FA") {
      WHAMFracAdd = c("FA")
    } else if (ThisProblem$InVar$Type[iInVar] == "WHAM-HAFA") {
      WHAMFracAdd = c("HA", "FA")
      if (!any(ThisProblem$InVar$Type[ThisProblem$InVar$MCR == ThisProblem$InVar$MCR[iInVar]] %in% "PercHA")) {
        stop("Must have PercHA input variable in mass compartment if specifying a WHAM-HAFA input variable.") # nolint: line_length_linter.
      }
    }
    if ((ThisProblem$InVar$Type[iInVar] %in% c("WHAM-FA", "WHAM-HA")) &&
        any(ThisProblem$InVar$Type[ThisProblem$InVar$MCR == ThisProblem$InVar$MCR[iInVar]] %in% "PercHA")) {
      stop("PercHA input variable specified in mass compartment with WHAM-HA or WHAM-FA input variable.") # nolint: line_length_linter.
    }
    if ((ThisProblem$InVar$Type[iInVar] %in% c("WHAM-HA")) &&
        any(ThisProblem$InVar$Type[ThisProblem$InVar$MCR == ThisProblem$InVar$MCR[iInVar]] %in% "PercAFA")) {
      stop("PercAFA input variable specified in mass compartment with WHAM-HA input variable.") # nolint: line_length_linter.
    }
    NWHAMFracAdd = length(WHAMFracAdd)

    WHAMprefix = array(
      paste0(ThisProblem$InVar$Name[iInVar], "-", WHAMFracAdd, "_"),
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

    DonnanCompName = paste0("Donnan", WHAMFracAdd)
    DonnanMCName = paste(ThisProblem$Mass$Name[iMass], DonnanCompName, sep = "_")
    NewProblem = AddMassCompartments(ThisProblem = NewProblem,
                                     MassName = DonnanMCName,
                                     MassAmt = 1E-5,
                                     MassUnit = ThisProblem$Mass$Unit[iMass],
                                     InMass = FALSE)
    NewProblem = AddDefComps(ThisProblem = NewProblem,
                             DefCompName = DonnanCompName,
                             DefCompFromVar = WHAMprefix,
                             DefCompCharge = 0L,
                             DefCompMCName = DonnanMCName,
                             DefCompType = DonnanCompName,
                             DefCompActCorr = "None",
                             DefCompSiteDens = 1.0E-4,
                             InDefComp = FALSE)

    MonodentpKH = numeric(nMS)
    MonodentAbundance = numeric(nMS)
    BidentAbundance = numeric(nBP)
    TridentAbundance = numeric(nTG)
    for (OMType in 1:NWHAMFracAdd) {

      ColspKM = paste0("pKM", c("A", "B"), WHAMFracAdd[OMType])
      OMSpecType = paste0("WHAM", WHAMFracAdd[OMType])

      # Donnan Species
      for (iSpec in match(ChargedSpecName, ThisProblem$Spec$Name)){
        NewProblem = AddSpecies(
          ThisProblem = NewProblem,
          SpecName = paste0(DonnanCompName[OMType],"-",ThisProblem$Spec$Name[iSpec]),
          SpecMCName = DonnanMCName[OMType],
          SpecType = DonnanCompName[OMType],
          SpecActCorr = ThisProblem$Spec$ActCorr[iSpec],
          SpecCompNames = list(c(DonnanCompName[OMType],
                                 ThisProblem$Comp$Name[ThisProblem$SpecCompList[iSpec, 1:ThisProblem$Spec$NC[iSpec]]])),
          SpecCompStoichs = list(c(abs(ThisProblem$Spec$Charge[iSpec]),
                                   ThisProblem$SpecStoich[iSpec, ThisProblem$SpecCompList[iSpec, 1:ThisProblem$Spec$NC[iSpec]]])),
          SpecLogK = SpecKsel[ThisProblem$Spec$Name[iSpec], OMType],
          SpecDeltaH = ThisProblem$Spec$DeltaH[iSpec],
          SpecTempKelvin = ThisProblem$Spec$TempKelvin[iSpec],
          InSpec = FALSE
        )
      }

      # Monodentate sites
      MonodentpKH[1:nStrong] = pKHA[OMType] + dpKHA[OMType] *
        (2 * MonodentTable$S[1:nStrong] - 5) / 6
      MonodentpKH[(nStrong + 1):nMS] = pKHB[OMType] + dpKHB[OMType] *
        (2 * MonodentTable$S[(nStrong + 1):nMS] - 13) / 6
      MonodentAbundance = (1 - fprB[OMType] - fprT[OMType]) *
        nCOOH[OMType] / MonodentTable$AbundDenom

      # Components - fully protonated
      NewProblem = AddDefComps(
        ThisProblem = NewProblem,
        DefCompName = paste0(WHAMprefix[OMType], MonodentTable$FullyProt),
        DefCompFromVar = WHAMprefix[OMType],
        DefCompCharge = 0L,
        DefCompMCName = ThisProblem$Mass$Name[iMass],
        DefCompType = OMSpecType,
        DefCompActCorr = "None",
        DefCompSiteDens = MonodentAbundance * 2E-3,
        InDefComp = FALSE
      )

      # - fully deprotonated
      NewProblem = AddSpecies(
        ThisProblem = NewProblem,
        SpecName = paste0(WHAMprefix[OMType], MonodentTable$FullyDeprot),
        SpecMCName = ThisProblem$Mass$Name[iMass],
        SpecType = OMSpecType,
        SpecActCorr = "None",
        SpecCompNames = as.list(as.data.frame(rbind(
          paste0(WHAMprefix[OMType], MonodentTable$FullyProt),
          rep(ThisProblem$Comp$Name[iH], nMS)
        ))),
        SpecCompStoichs = as.list(as.data.frame(rbind(
          rep(1, nMS),
          rep(-1, nMS)
        ))),
        SpecLogK = -1 * MonodentpKH,
        SpecDeltaH = 0,
        SpecTempKelvin = 298.15,
        InSpec = FALSE
      )

      # bound to each metal
      iDeprotSpec = which(NewProblem$Spec$Name %in% paste0(WHAMprefix[OMType], MonodentTable$FullyDeprot))
      for (iMetal in 1:nMP) {
        iMetalSpec = which(ThisProblem$Spec$Name == MetalsTable$Metal[iMetal])
        NewProblem = AddSpecies(
          ThisProblem = NewProblem,
          SpecEquation = paste0(
            WHAMprefix[OMType], MonodentTable$FullyDeprot, "-",
            ThisProblem$Spec$Equation[iMetalSpec],
            " -1 * H + 1 * ", WHAMprefix[OMType], MonodentTable$FullyProt),
          SpecMCName = ThisProblem$Mass$Name[iMass],
          SpecType = OMSpecType,
          SpecActCorr = "None",
          SpecLogK = ThisProblem$Spec$LogK[iMetalSpec] -
            as.numeric(MetalsTable[iMetal, ColspKM[MonodentTable$Strong1Weak2]]),
          SpecDeltaH = ThisProblem$Spec$DeltaH[iMetalSpec],
          SpecTempKelvin = ThisProblem$Spec$TempKelvin[iMetalSpec],
          InSpec = FALSE
        )
      }



      # Bidentate sites
      if (nBP > 0) {
        BidentAbundance = fprB[OMType] * nCOOH[OMType] / BidentTable$AbundDenom

        # Components - fully protonated
        NewProblem = AddDefComps(
          ThisProblem = NewProblem,
          DefCompName = paste0(WHAMprefix[OMType], BidentTable$FullyProt),
          DefCompFromVar = WHAMprefix[OMType],
          DefCompCharge = 0L,
          DefCompMCName = ThisProblem$Mass$Name[iMass],
          DefCompType = OMSpecType,
          DefCompActCorr = "None",
          DefCompSiteDens = BidentAbundance * 2E-3,
          InDefComp = FALSE
        )

        # - first site deprotonated
        NewProblem = AddSpecies(
          ThisProblem = NewProblem,
          SpecEquation = paste0(
            WHAMprefix[OMType], BidentTable$S1Deprot, " = ",
            "1 * ", WHAMprefix[OMType], BidentTable$FullyProt, " -1 * H"
          ),
          SpecMCName = ThisProblem$Mass$Name[iMass],
          SpecType = OMSpecType,
          SpecActCorr = "None",
          SpecLogK = -1 * MonodentpKH[BidentTable$S1],
          SpecDeltaH = 0,
          SpecTempKelvin = 298.15,
          InSpec = FALSE
        )

        # - second site deprotonated
        NewProblem = AddSpecies(
          ThisProblem = NewProblem,
          SpecEquation = paste0(
            WHAMprefix[OMType], BidentTable$S2Deprot, " = ",
            "1 * ", WHAMprefix[OMType], BidentTable$FullyProt, " -1 * H"
          ),
          SpecMCName = ThisProblem$Mass$Name[iMass],
          SpecType = OMSpecType,
          SpecActCorr = "None",
          SpecLogK = -1 * MonodentpKH[BidentTable$S2],
          SpecDeltaH = 0,
          SpecTempKelvin = 298.15,
          InSpec = FALSE
        )

        # - fully deprot
        NewProblem = AddSpecies(
          ThisProblem = NewProblem,
          SpecEquation = paste0(
            WHAMprefix[OMType], BidentTable$FullyDeprot, " = ",
            "1 * ", WHAMprefix[OMType], BidentTable$FullyProt, " -2 * H"
          ),
          SpecMCName = ThisProblem$Mass$Name[iMass],
          SpecType = OMSpecType,
          SpecActCorr = "None",
          SpecLogK = -1 * (MonodentpKH[BidentTable$S1] + MonodentpKH[BidentTable$S2]),
          SpecDeltaH = 0,
          SpecTempKelvin = 298.15,
          InSpec = FALSE
        )

        # bound to each metal
        iDeprotSpec = which(NewProblem$Spec$Name %in%
                              paste0(WHAMprefix[OMType], BidentTable$FullyDeprot))
        for (iMetal in 1:nMP) {
          iMetalSpec = which(ThisProblem$Spec$Name == MetalsTable$Metal[iMetal])
          NewProblem = AddSpecies(
            ThisProblem = NewProblem,
            SpecEquation = paste0(
              WHAMprefix[OMType], BidentTable$FullyDeprot, "-",
              ThisProblem$Spec$Equation[iMetalSpec],
              " -2 * H + 1 * ", WHAMprefix[OMType], BidentTable$FullyProt),
            SpecMCName = ThisProblem$Mass$Name[iMass],
            SpecType = OMSpecType,
            SpecActCorr = "None",
            SpecLogK = ThisProblem$Spec$LogK[iMetalSpec] -
              as.numeric(MetalsTable[iMetal, ColspKM[BidentTable$S1Strong1Weak2]] +
                           MetalsTable[iMetal, ColspKM[BidentTable$S2Strong1Weak2]]),
            SpecDeltaH = ThisProblem$Spec$DeltaH[iMetalSpec],
            SpecTempKelvin = ThisProblem$Spec$TempKelvin[iMetalSpec],
            InSpec = FALSE
          )

        }
      }

      # Tridentate sites
      if (nTG > 0) {

        TridentAbundance = fprT[OMType] * nCOOH[OMType] / TridentTable$AbundDenom

        # Components - fully protonated
        NewProblem = AddDefComps(
          ThisProblem = NewProblem,
          DefCompName = paste0(WHAMprefix[OMType], TridentTable$FullyProt),
          DefCompFromVar = WHAMprefix[OMType],
          DefCompCharge = 0L,
          DefCompMCName = ThisProblem$Mass$Name[iMass],
          DefCompType = OMSpecType,
          DefCompActCorr = "None",
          DefCompSiteDens = TridentAbundance * 2E-3,# the input is in mg C/L, while nCOOH is mols/g HS
          InDefComp = FALSE
        )

        # - first site deprotonated
        NewProblem = AddSpecies(
          ThisProblem = NewProblem,
          SpecEquation = paste0(
            WHAMprefix[OMType], TridentTable$S1Deprot, " = ",
            "1 * ", WHAMprefix[OMType], TridentTable$FullyProt, " -1 * H"
          ),
          SpecMCName = ThisProblem$Mass$Name[iMass],
          SpecType = OMSpecType,
          SpecActCorr = "None",
          SpecLogK = -1 * MonodentpKH[TridentTable$S1],
          SpecDeltaH = 0,
          SpecTempKelvin = 298.15,
          InSpec = FALSE
        )

        # - second site deprotonated
        NewProblem = AddSpecies(
          ThisProblem = NewProblem,
          SpecEquation = paste0(
            WHAMprefix[OMType], TridentTable$S2Deprot, " = ",
            "1 * ", WHAMprefix[OMType], TridentTable$FullyProt, " -1 * H"
          ),
          SpecMCName = ThisProblem$Mass$Name[iMass],
          SpecType = OMSpecType,
          SpecActCorr = "None",
          SpecLogK = -1 * MonodentpKH[TridentTable$S2],
          SpecDeltaH = 0,
          SpecTempKelvin = 298.15,
          InSpec = FALSE
        )

        # - third site deprotonated
        NewProblem = AddSpecies(
          ThisProblem = NewProblem,
          SpecEquation = paste0(
            WHAMprefix[OMType], TridentTable$S3Deprot, " = ",
            "1 * ", WHAMprefix[OMType], TridentTable$FullyProt, " -1 * H"
          ),
          SpecMCName = ThisProblem$Mass$Name[iMass],
          SpecType = OMSpecType,
          SpecActCorr = "None",
          SpecLogK = -1 * MonodentpKH[TridentTable$S3],
          SpecDeltaH = 0,
          SpecTempKelvin = 298.15,
          InSpec = FALSE
        )

        # - first & second sites deprotonated
        NewProblem = AddSpecies(
          ThisProblem = NewProblem,
          SpecEquation = paste0(
            WHAMprefix[OMType], TridentTable$S12Deprot, " = ",
            "1 * ", WHAMprefix[OMType], TridentTable$FullyProt, " -2 * H"
          ),
          SpecMCName = ThisProblem$Mass$Name[iMass],
          SpecType = OMSpecType,
          SpecActCorr = "None",
          SpecLogK = -1 * (MonodentpKH[TridentTable$S1] + MonodentpKH[TridentTable$S2]),
          SpecDeltaH = 0,
          SpecTempKelvin = 298.15,
          InSpec = FALSE
        )

        # - first & third sites deprotonated
        NewProblem = AddSpecies(
          ThisProblem = NewProblem,
          SpecEquation = paste0(
            WHAMprefix[OMType], TridentTable$S13Deprot, " = ",
            "1 * ", WHAMprefix[OMType], TridentTable$FullyProt, " -2 * H"
          ),
          SpecMCName = ThisProblem$Mass$Name[iMass],
          SpecType = OMSpecType,
          SpecActCorr = "None",
          SpecLogK = -1 * (MonodentpKH[TridentTable$S1] + MonodentpKH[TridentTable$S3]),
          SpecDeltaH = 0,
          SpecTempKelvin = 298.15,
          InSpec = FALSE
        )

        # - second & third sites deprotonated
        NewProblem = AddSpecies(
          ThisProblem = NewProblem,
          SpecEquation = paste0(
            WHAMprefix[OMType], TridentTable$S22Deprot, " = ",
            "1 * ", WHAMprefix[OMType], TridentTable$FullyProt, " -2 * H"
          ),
          SpecMCName = ThisProblem$Mass$Name[iMass],
          SpecType = OMSpecType,
          SpecActCorr = "None",
          SpecLogK = -1 * (MonodentpKH[TridentTable$S2] + MonodentpKH[TridentTable$S2]),
          SpecDeltaH = 0,
          SpecTempKelvin = 298.15,
          InSpec = FALSE
        )

        # - fully deprot
        NewProblem = AddSpecies(
          ThisProblem = NewProblem,
          SpecEquation = paste0(
            WHAMprefix[OMType], TridentTable$FullyDeprot, " = ",
            "1 * ", WHAMprefix[OMType], TridentTable$FullyProt, " -3 * H"
          ),
          SpecMCName = ThisProblem$Mass$Name[iMass],
          SpecType = OMSpecType,
          SpecActCorr = "None",
          SpecLogK = -1 * (MonodentpKH[TridentTable$S1] +
                             MonodentpKH[TridentTable$S2] +
                             MonodentpKH[TridentTable$S3]),
          SpecDeltaH = 0,
          SpecTempKelvin = 298.15,
          InSpec = FALSE
        )

        # bound to each metal
        iDeprotSpec = which(NewProblem$Spec$Name %in%
                              paste0(WHAMprefix[OMType], TridentTable$FullyDeprot))
        for (iMetal in 1:nMP) {
          iMetalSpec = which(ThisProblem$Spec$Name == MetalsTable$Metal[iMetal])
          NewProblem = AddSpecies(
            ThisProblem = NewProblem,
            SpecEquation = paste0(
              WHAMprefix[OMType], TridentTable$FullyDeprot, "-",
              ThisProblem$Spec$Equation[iMetalSpec],
              " -3 * H + 1 * ", WHAMprefix[OMType], TridentTable$FullyProt),
            SpecMCName = ThisProblem$Mass$Name[iMass],
            SpecType = OMSpecType,
            SpecActCorr = "None",
            SpecLogK = ThisProblem$Spec$LogK[iMetalSpec] -
              as.numeric(MetalsTable[iMetal, ColspKM[TridentTable$S1Strong1Weak2]] +
                           MetalsTable[iMetal, ColspKM[TridentTable$S2Strong1Weak2]] +
                           MetalsTable[iMetal, ColspKM[TridentTable$S3Strong1Weak2]]),
            SpecDeltaH = ThisProblem$Spec$DeltaH[iMetalSpec],
            SpecTempKelvin = ThisProblem$Spec$TempKelvin[iMetalSpec],
            InSpec = FALSE
          )

        }
      }

    }

  }

  # WHAM parameters - to be used later
  NewProblem$WHAM$DoWHAM = TRUE
  NewProblem$WHAM$WHAMVer = WHAMVer
  NewProblem$WHAM$WdatFile = WdatFile
  NewProblem$WHAM$wDLF = wDLF
  NewProblem$WHAM$wKZED = wKZED
  NewProblem$WHAM$wP = wP
  NewProblem$WHAM$wRadius = wRadius
  NewProblem$WHAM$wMolWt = wMolWt

  return(NewProblem)

}
