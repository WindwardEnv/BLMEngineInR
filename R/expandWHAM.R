#' Expand the DOC component into WHAM components
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
#' @param SpecTemp numeric vector of length `NSpec`; temperature at which the
#'   logK/deltaH were measured (modified and returned)
#' @param NPhase integer; Number of phases
#' @param PhaseCompList integer matrix of `NPhase` rows and `max(PhaseNC)+2`
#'   columns; the list of components used to create a given phase (modified and
#'   returned)
#' @param PhaseStoich integer matrix of `NPhase` rows and `NComp` columns; the
#'   stoichiometry matrix of the phase reactions (modified and returned)
#' @param WHAMVer a character string specifying the WHAM version to use, must be
#'   one of `"V"` (default), `"VI"`, or `"VII"`. Ignored if `wdatFile` is not
#'   `NULL`.
#' @param wdatFile (optional) a character string specifying the file path of a
#'   WHAM parameter file
#'
#' @keywords internal
#'
#' @noRd
expandWHAM = function(NMass, MassName,
                      NInVar, InVarName, InVarMC, InVarType,
                      NComp, CompName, CompCharge, CompMC, CompType, CompActCorr, CompSiteDens,
                      NDefComp, DefCompName, DefCompFromNum, DefCompFromVar, DefCompCharge, DefCompMC, DefCompType, DefCompActCorr, DefCompSiteDens,
                      NSpec, SpecName, SpecMC, SpecActCorr, SpecNC, SpecCompList, SpecStoich, SpecLogK, SpecDeltaH, SpecTemp,
                      NPhase, PhaseCompList, PhaseStoich,
                      WHAMVer = c("V", "VI", "VII"),
                      wdatFile = NULL) {


  # error catching and input cleanup
  if (is.null(wdatFile)){
    WHAMVer = match.arg(WHAMVer)
    if(WHAMVer == "V"){
      wdatFile = system.file("extdata/WHAM_V.wdat", package = "BLMEngineInR", mustWork = T)
    } else if (WHAMVer == "VI"){
      wdatFile = system.file("extdata/WHAM_VI.wdat", package = "BLMEngineInR", mustWork = T)
    } else if (WHAMVer  == "VII"){
      wdatFile = system.file("extdata/WHAM_VII.wdat", package = "BLMEngineInR", mustWork = T)
    }
  } else {
    wdatFile = normalizePath(wdatFile)
    stopifnot(file.exists(wdatFile))
  }


  # read WHAM data file
  {

    # header info
    skipRows = 2
    temp = read.delim(file = wdatFile, header = F, sep = ",",
                      skip = skipRows, nrows = 7)
    nMS = as.integer(temp[1,2])#Number of monodentate sites
    nBP = as.integer(temp[2,2])#Number of bidentate pairs
    nTG = as.integer(temp[3,2])#Number of tridentate groups
    nMP = as.integer(temp[4,2])#Number of metals-OM parameters
    wDLF = as.numeric(temp[5,2])#Double layer overlap factor
    wKZED = as.numeric(temp[6,2])#Constant to control DDL at low ZED
    wKsel = as.numeric(temp[7,2])#Selectivity coefficient Ksel

    # Parameters
    skipRows = skipRows + 7 + 1
    temp = read.delim(file = wdatFile, header = T, sep = ",",
                      skip = skipRows, nrows = 12)
    nCOOH = array(as.numeric(temp[1,3:4]), dimnames = list(c("HA","FA")))
    pKHA = array(as.numeric(temp[2,3:4]), dimnames = list(c("HA","FA")))
    pKHB = array(as.numeric(temp[3,3:4]), dimnames = list(c("HA","FA")))
    dpKHA = array(as.numeric(temp[4,3:4]), dimnames = list(c("HA","FA")))
    dpKHB = array(as.numeric(temp[5,3:4]), dimnames = list(c("HA","FA")))
    fprB = array(as.numeric(temp[6,3:4]), dimnames = list(c("HA","FA")))
    fprT = array(as.numeric(temp[7,3:4]), dimnames = list(c("HA","FA")))
    dLK1A = array(as.numeric(temp[8,3:4]), dimnames = list(c("HA","FA")))
    dLK1B = array(as.numeric(temp[9,3:4]), dimnames = list(c("HA","FA")))
    wP = array(as.numeric(temp[10,3:4]), dimnames = list(c("HA","FA")))
    wRadius = array(as.numeric(temp[11,3:4]), dimnames = list(c("HA","FA")))
    wMolWt = array(as.numeric(temp[12,3:4]), dimnames = list(c("HA","FA")))

    # Monodentate Sites - these should always be the same, but we'll set things
    # up like this so we can add in this section if it's ever needed.
    # MonodentTable = data.frame(S=1:8, AbundDenom = c(rep(4,4),rep(8,4)))
    skipRows = skipRows + 12 + 3
    MonodentTable = read.delim(file = wdatFile, header = T, sep = ",",
                               skip = skipRows, nrows = nMS)
    MonodentTable$FullyProt = paste0(MonodentTable$S,"H")
    MonodentTable$FullyDeprot = paste0(MonodentTable$S)
    MonodentTable$Strong1Weak2 = c(rep(1L,4),rep(2L,4))

    # Bidentate Pairs
    skipRows = skipRows + nMS + 3
    if (nBP > 0){
      BidentTable = read.delim(file = wdatFile, header = T, sep = ",",
                               skip = skipRows, nrows = nBP)
      names(BidentTable) = c("S1","S2","AbundDenom")
      BidentTable$FullyProt = paste0(BidentTable$S1,BidentTable$S2,"H")
      BidentTable$S1Deprot = paste0(BidentTable$S1,"-",BidentTable$S2,"H")
      BidentTable$S2Deprot = paste0(BidentTable$S2,"-",BidentTable$S1,"H")
      BidentTable$FullyDeprot = paste0(BidentTable$S1,BidentTable$S2)
      BidentTable$S1Strong1Weak2 = ifelse(BidentTable$S1<=4,1,2)
      BidentTable$S2Strong1Weak2 = ifelse(BidentTable$S2<=4,1,2)
    } else {
      BidentTable = data.frame()
    }


    # Tridentate Groups
    skipRows = skipRows + nBP + 3
    if (nTG > 0){
      TridentTable = read.delim(file = wdatFile, header = T, sep = ",",
                                skip = skipRows, nrows = nTG)
      names(TridentTable) = c("S1","S2","S3","AbundDenom")
      TridentTable$FullyProt = paste0(TridentTable$S1,TridentTable$S2,TridentTable$S3,"H")
      TridentTable$S1Deprot = paste0(TridentTable$S1,"-",TridentTable$S2,TridentTable$S3,"H")
      TridentTable$S2Deprot = paste0(TridentTable$S2,"-",TridentTable$S1,TridentTable$S3,"H")
      TridentTable$S3Deprot = paste0(TridentTable$S3,"-",TridentTable$S1,TridentTable$S2,"H")
      TridentTable$S12Deprot = paste0(TridentTable$S1,TridentTable$S2,"-",TridentTable$S3,"H")
      TridentTable$S13Deprot = paste0(TridentTable$S1,TridentTable$S3,"-",TridentTable$S2,"H")
      TridentTable$S23Deprot = paste0(TridentTable$S2,TridentTable$S3,"-",TridentTable$S1,"H")
      TridentTable$FullyDeprot = paste0(TridentTable$S1,TridentTable$S2,TridentTable$S3)
      TridentTable$S1Strong1Weak2 = ifelse(TridentTable$S1<=4,1,2)
      TridentTable$S2Strong1Weak2 = ifelse(TridentTable$S2<=4,1,2)
      TridentTable$S3Strong1Weak2 = ifelse(TridentTable$S3<=4,1,2)
    } else {
      TridentTable = data.frame()
    }

    # Metals Parameters Table
    skipRows = skipRows + nTG + 3
    if (nMP > 0){
      MetalsTable = read.delim(file = wdatFile, header = T, sep = ",",
                               skip = skipRows, nrows = nMP)
      names(MetalsTable) = c("Metal","pKMAHA","pKMAFA","dLK2")
      MetalsTable = MetalsTable[MetalsTable$Metal %in% c(CompName,SpecName),]
      nMP = nrow(MetalsTable)
      MetalsTable$pKMBHA = 3 * MetalsTable$pKMAHA - 3
      MetalsTable$pKMBFA = 3.96 * MetalsTable$pKMAFA
    } else {
      MetalsTable = data.frame()
    }

  }

  # Initialize variables
  iH = which(CompName == "H")

  # Figure out the number of DOC components we're adding, and what fraction
  InVarWHAM = which(grepl("WHAM",InVarType))

  for (iInVar in InVarWHAM){

    iMass = InVarMC[iInVar]

    if (InVarType[iInVar] == "WHAM-HA") {
      WHAMFracAdd = c("HA")
    } else if (InVarType[iInVar] == "WHAM-FA") {
      WHAMFracAdd = c("FA")
    } else if (InVarType[iInVar] == "WHAM-HAFA") {
      WHAMFracAdd = c("HA", "FA")
      if (!any(InVarType[InVarMC == InVarMC[iInVar]] %in% "PercHA")) {
        stop(
          "Must have PercHA input variable in mass compartment if specifying a WHAM-HAFA input variable."
        )
      }
    }
    if ((InVarType[iInVar] %in% c("WHAM-FA","WHAM-HA")) &
        any(InVarType[InVarMC == InVarMC[iInVar]] %in% "PercHA")) {
      stop(
        "PercHA input variable specified in mass compartment with WHAM-HA or WHAM-FA input variable."
      )
    }
    if ((InVarType[iInVar] %in% c("WHAM-HA")) &
        any(InVarType[InVarMC == InVarMC[iInVar]] %in% "PercAFA")) {
      stop(
        "PercAFA input variable specified in mass compartment with WHAM-HA input variable."
      )
    }
    nWHAMFracAdd = length(WHAMFracAdd)

    WHAMprefix = array(
      paste0(InVarName[iInVar], "-", WHAMFracAdd, "_"),
      dim = nWHAMFracAdd,
      dimnames = list(WHAMFracAdd)
    )



    # * Each component has the fully protonated species as the component
    # * Each component will have every possible combination of binding sites deprotonated
    # * Each component, when fully deprotonated, will bind to each metal nMP.
    #     monodentate example:
    #      component: FA1H
    #      species: FA1H, FA1, FA1-Mg, FA1-Ca, ...
    #     bidentate example:
    #      component: FA12H
    #      species: FA12H, FA1-2H, FA2-1H, FA12, FA12-Mg, FA12-Ca, ...
    #     tridentate example:
    #      component: FA123H
    #      species: FA123H, FA1-23H, FA2-13H, FA3-12H, FA23-1H, FA13-2H, FA12-3H, FA123, FA123-Mg, FA123-Ca, ...
    wNComp = (nMS + nBP + nTG) * nWHAMFracAdd
    startComp = NComp + 1L
    NComp = NComp + wNComp
    wCompName = paste0(rep(WHAMprefix, each = wNComp / nWHAMFracAdd),
                       rep(c(MonodentTable$FullyProt,
                             BidentTable$FullyProt,
                             TridentTable$FullyProt), times = nWHAMFracAdd))
    CompName = c(CompName, wCompName)
    CompCharge = c(CompCharge, array(0L, dim = wNComp, dimnames = list(wCompName)))
    CompMC = c(CompMC, array(iMass, dim = wNComp, dimnames = list(wCompName)))
    CompType = c(CompType, array("MassBal", wNComp, dimnames = list(wCompName)))
    CompActCorr = c(CompActCorr, array("WHAM", wNComp, dimnames = list(wCompName)))
    CompSiteDens = c(CompSiteDens, array(NA, wNComp, dimnames = list(wCompName)))

    startDefComp = NDefComp + 1L
    NDefComp = NDefComp + wNComp
    DefCompName = c(DefCompName, wCompName)
    DefCompFromNum = c(DefCompFromNum, array(NA, dim = wNComp, dimnames = list(wCompName)))
    # if (nWHAMFracAdd > 1){
      DefCompFromVar = c(DefCompFromVar, array(rep(WHAMprefix, each = wNComp / nWHAMFracAdd), dim=wNComp, dimnames = list(wCompName)))
    # } else {
    #   DefCompFromVar = c(DefCompFromVar, array(InVarName[iInVar], dim=wNComp, dimnames = list(wCompName)))
    # }
    DefCompCharge = c(DefCompCharge, array(0L, dim = wNComp, dimnames = list(wCompName)))
    DefCompMC = c(DefCompMC, array(iMass, dim = wNComp, dimnames = list(wCompName)))
    DefCompType = c(DefCompType, array("MassBal", wNComp, dimnames = list(wCompName)))
    DefCompActCorr = c(DefCompActCorr, array("WHAM", wNComp, dimnames = list(wCompName)))
    DefCompSiteDens = c(DefCompSiteDens, array(NA, wNComp, dimnames = list(wCompName)))

    wNSpec = (nMS * (2L + nMP) + nBP * (4L + nMP) + nTG * (8L + nMP)) * nWHAMFracAdd
    startSpec = NSpec + 1L
    NSpec = NSpec + wNSpec
    SpecName = c(SpecName, array(paste0("newOCSpecies",1:wNSpec),wNSpec))
    SpecMC = c(SpecMC, array(iMass, dim = wNSpec,dimnames=list(paste0("newOCSpecies",1:wNSpec))))#this should always be water
    # SpecType = c(SpecType, array(1L, wNSpec, dimnames=list(paste0("newOCSpecies",1:wNSpec))))
    SpecActCorr = c(SpecActCorr, array("WHAM", wNSpec, dimnames=list(paste0("newOCSpecies",1:wNSpec))))
    # SpecCharge = c(SpecCharge, array(NA, wNSpec, dimnames=list(paste0("newOCSpecies",1:wNSpec))))
    SpecNC = c(SpecNC, array(NA, wNSpec, dimnames=list(paste0("newOCSpecies",1:wNSpec))))
    SpecCompList = rbind(SpecCompList, matrix(0L, nrow = wNSpec, ncol = ncol(SpecCompList)))
    SpecStoich = rbind(cbind(SpecStoich,
                         matrix(0L,nrow=NSpec-wNSpec, ncol = wNComp,
                                dimnames = list(SpecName[1:(NSpec-wNSpec)],wCompName))),
                   matrix(0L, nrow = wNSpec, ncol = NComp,
                          dimnames = list(paste0("newOCSpecies",1:wNSpec),
                                          CompName)))
    SpecLogK = c(SpecLogK, array(NA,wNSpec,dimnames=list(paste0("newOCSpecies",1:wNSpec))))
    SpecDeltaH = c(SpecDeltaH, array(0.0,wNSpec,dimnames=list(paste0("newOCSpecies",1:wNSpec))))
    SpecTemp = c(SpecTemp, array(0.0,wNSpec,dimnames=list(paste0("newOCSpecies",1:wNSpec))))

    Monodent_pKH = numeric(nMS)
    Monodent_Abundance = numeric(nMS)
    Bident_Abundance = numeric(nBP)
    Trident_Abundance = numeric(nTG)
    for(OMType in WHAMFracAdd){

      pKM_cols = paste0("pKM",c("A","B"),OMType)

      # Monodentate sites
      newCompNum = startComp:(startComp + nMS - 1)
      newDefCompNum = startDefComp:(startDefComp + nMS - 1)
      Monodent_pKH[1:4] = pKHA[OMType] + dpKHA[OMType] * (2 * MonodentTable$S[1:4] - 5) / 6
      Monodent_pKH[5:8] = pKHB[OMType] + dpKHB[OMType] * (2 * MonodentTable$S[5:8] - 13) / 6
      Monodent_Abundance = (1 - fprB[OMType] - fprT[OMType]) * nCOOH[OMType] / MonodentTable$AbundDenom
      CompSiteDens[newCompNum] = Monodent_Abundance * 1E-3 # the input is in milligrams, while nCOOH is mols/g
      DefCompSiteDens[newDefCompNum] = Monodent_Abundance * 1E-3 # the input is in milligrams, while nCOOH is mols/g

      # - fully protonated (components)
      newSpecNum = startSpec:(startSpec + nMS - 1)
      SpecName[newSpecNum] = paste0(WHAMprefix[OMType], MonodentTable$FullyProt)
      # SpecCharge[newSpecNum] = 0L
      diag(SpecStoich[newSpecNum, newCompNum]) = 1L
      SpecLogK[newSpecNum] = 0.0

      # - fully deprot
      newSpecNum = newSpecNum + nMS
      SpecName[newSpecNum] = paste0(WHAMprefix[OMType], MonodentTable$FullyDeprot)
      # SpecCharge[newSpecNum] = -1L
      diag(SpecStoich[newSpecNum, newCompNum]) = 1L
      SpecStoich[newSpecNum,iH] = -1L
      SpecLogK[newSpecNum] = -1 * Monodent_pKH

      # bound to each metal
      for(iMetal in 1:nMP){
        iMetalSpec = which(MetalsTable$Metal[iMetal] == SpecName)
        newSpecNum = newSpecNum + nMS
        # SpecCharge[newSpecNum] = -1L + SpecCharge[iMetalSpec]
        SpecName[newSpecNum] = paste0(WHAMprefix[OMType], MonodentTable$FullyDeprot, "-",
                                      MetalsTable$Metal[iMetal])
        SpecStoich[newSpecNum, 1:NComp] = matrix(SpecStoich[iMetalSpec,], nrow = nMS, ncol = NComp, byrow = T)
        diag(SpecStoich[newSpecNum, newCompNum]) = 1L
        SpecStoich[newSpecNum,iH] = SpecStoich[newSpecNum,iH] -1L
        SpecLogK[newSpecNum] = -1 * as.numeric(MetalsTable[iMetal,pKM_cols[MonodentTable$Strong1Weak2]])
      }


      # Bidentate sites
      if (nBP > 0){
        startComp = max(newCompNum) + 1
        startDefComp = max(newDefCompNum) + 1
        startSpec = max(newSpecNum) + 1
        newCompNum = startComp:(startComp + nBP - 1)
        newDefCompNum = startDefComp:(startDefComp + nBP - 1)
        Bident_Abundance = fprB[OMType] * nCOOH[OMType] / BidentTable$AbundDenom
        CompSiteDens[newCompNum] = Bident_Abundance * 1E-3 # the input is in milligrams, while nCOOH is mols/g
        DefCompSiteDens[newDefCompNum] = Bident_Abundance * 1E-3 # the input is in milligrams, while nCOOH is mols/g

        # - fully protonated
        newSpecNum = startSpec:(startSpec + nBP - 1)
        SpecName[newSpecNum] = paste0(WHAMprefix[OMType], BidentTable$FullyProt)
        # SpecCharge[newSpecNum] = 0L
        diag(SpecStoich[newSpecNum, newCompNum]) = 1L
        SpecLogK[newSpecNum] = 0.0

        # - first site deprotonated
        newSpecNum = newSpecNum + nBP
        SpecName[newSpecNum] = paste0(WHAMprefix[OMType], BidentTable$S1Deprot)
        # SpecCharge[newSpecNum] = -1L
        diag(SpecStoich[newSpecNum, newCompNum]) = 1L
        SpecStoich[newSpecNum,iH] = -1L
        SpecLogK[newSpecNum] = -1 * Monodent_pKH[BidentTable$S1]

        # - second site deprotonated
        newSpecNum = newSpecNum + nBP
        SpecName[newSpecNum] = paste0(WHAMprefix[OMType], BidentTable$S2Deprot)
        # SpecCharge[newSpecNum] = -1L
        diag(SpecStoich[newSpecNum, newCompNum]) = 1L
        SpecStoich[newSpecNum,iH] = -1L
        SpecLogK[newSpecNum] = -1 * Monodent_pKH[BidentTable$S2]

        # - fully deprot
        newSpecNum = newSpecNum + nBP
        SpecName[newSpecNum] = paste0(WHAMprefix[OMType], BidentTable$FullyDeprot)
        # SpecCharge[newSpecNum] = -2L
        diag(SpecStoich[newSpecNum, newCompNum]) = 1L
        SpecStoich[newSpecNum,iH] = -2L
        SpecLogK[newSpecNum] = -1 * (Monodent_pKH[BidentTable$S1] +
                                   Monodent_pKH[BidentTable$S2])

        # bound to each metal
        for(iMetal in 1:nMP){
          iMetalSpec = which(MetalsTable$Metal[iMetal] == SpecName)
          newSpecNum = newSpecNum + nBP
          # SpecCharge[newSpecNum] = -2L + SpecCharge[iMetalSpec]
          SpecName[newSpecNum] = paste0(WHAMprefix[OMType], BidentTable$FullyDeprot, "-",
                                        MetalsTable$Metal[iMetal])
          SpecStoich[newSpecNum, 1:NComp] = matrix(SpecStoich[iMetalSpec,], nrow = nBP, ncol = NComp, byrow = T)
          diag(SpecStoich[newSpecNum, newCompNum]) = 1L
          SpecStoich[newSpecNum,iH] = SpecStoich[newSpecNum,iH] -2L
          SpecLogK[newSpecNum] = -1 * as.numeric(MetalsTable[iMetal,pKM_cols[BidentTable$S1Strong1Weak2]] +
                                               MetalsTable[iMetal,pKM_cols[BidentTable$S2Strong1Weak2]])
        }

      }

      # Tridentate sites
      if (nTG > 0){
        startComp = max(newCompNum) + 1
        startDefComp = max(newDefCompNum) + 1
        startSpec = max(newSpecNum) + 1
        newCompNum = startComp:(startComp + nTG - 1)
        newDefCompNum = startDefComp:(startDefComp + nTG - 1)
        Trident_Abundance = fprT[OMType] * nCOOH[OMType] / TridentTable$AbundDenom
        CompSiteDens[newCompNum] = Trident_Abundance * 1E-3 # the input is in milligrams, while nCOOH is mols/g
        DefCompSiteDens[newDefCompNum] = Trident_Abundance * 1E-3 # the input is in milligrams, while nCOOH is mols/g

        # - fully protonated
        newSpecNum = startSpec:(startSpec + nTG - 1)
        SpecName[newSpecNum] = paste0(WHAMprefix[OMType], TridentTable$FullyProt)
        # SpecCharge[newSpecNum] = 0L
        diag(SpecStoich[newSpecNum, newCompNum]) = 1L
        SpecLogK[newSpecNum] = 0.0

        # - first site deprotonated
        newSpecNum = newSpecNum + nTG
        SpecName[newSpecNum] = paste0(WHAMprefix[OMType], TridentTable$S1Deprot)
        # SpecCharge[newSpecNum] = -1L
        diag(SpecStoich[newSpecNum, newCompNum]) = 1L
        SpecStoich[newSpecNum,iH] = -1L
        SpecLogK[newSpecNum] = -1 * Monodent_pKH[TridentTable$S1]

        # - second site deprotonated
        newSpecNum = newSpecNum + nTG
        SpecName[newSpecNum] = paste0(WHAMprefix[OMType], TridentTable$S2Deprot)
        # SpecCharge[newSpecNum] = -1L
        diag(SpecStoich[newSpecNum, newCompNum]) = 1L
        SpecStoich[newSpecNum,iH] = -1L
        SpecLogK[newSpecNum] = -1 * Monodent_pKH[TridentTable$S2]

        # - third site deprotonated
        newSpecNum = newSpecNum + nTG
        SpecName[newSpecNum] = paste0(WHAMprefix[OMType], TridentTable$S3Deprot)
        # SpecCharge[newSpecNum] = -1L
        diag(SpecStoich[newSpecNum, newCompNum]) = 1L
        SpecStoich[newSpecNum,iH] = -1L
        SpecLogK[newSpecNum] = -1 * Monodent_pKH[TridentTable$S3]

        # - first & second sites deprotonated
        newSpecNum = newSpecNum + nTG
        SpecName[newSpecNum] = paste0(WHAMprefix[OMType], TridentTable$S12Deprot)
        # SpecCharge[newSpecNum] = -2L
        diag(SpecStoich[newSpecNum, newCompNum]) = 1L
        SpecStoich[newSpecNum,iH] = -2L
        SpecLogK[newSpecNum] = -1 * (Monodent_pKH[TridentTable$S1] + Monodent_pKH[TridentTable$S2])

        # - first & third sites deprotonated
        newSpecNum = newSpecNum + nTG
        SpecName[newSpecNum] = paste0(WHAMprefix[OMType], TridentTable$S13Deprot)
        # SpecCharge[newSpecNum] = -2L
        diag(SpecStoich[newSpecNum, newCompNum]) = 1L
        SpecStoich[newSpecNum,iH] = -2L
        SpecLogK[newSpecNum] = -1 * (Monodent_pKH[TridentTable$S1] + Monodent_pKH[TridentTable$S3])

        # - second & third sites deprotonated
        newSpecNum = newSpecNum + nTG
        SpecName[newSpecNum] = paste0(WHAMprefix[OMType], TridentTable$S23Deprot)
        # SpecCharge[newSpecNum] = -2L
        diag(SpecStoich[newSpecNum, newCompNum]) = 1L
        SpecStoich[newSpecNum,iH] = -2L
        SpecLogK[newSpecNum] = -1 * (Monodent_pKH[TridentTable$S2] + Monodent_pKH[TridentTable$S3])

        # - fully deprot
        newSpecNum = newSpecNum + nTG
        SpecName[newSpecNum] = paste0(WHAMprefix[OMType], TridentTable$FullyDeprot)
        # SpecCharge[newSpecNum] = -3L
        diag(SpecStoich[newSpecNum, newCompNum]) = 1L
        SpecStoich[newSpecNum,iH] = -3L
        SpecLogK[newSpecNum] = -1 *(Monodent_pKH[TridentTable$S1] +
                                  Monodent_pKH[TridentTable$S2] +
                                  Monodent_pKH[TridentTable$S3])

        # bound to each metal
        for(iMetal in 1:nMP){
          iMetalSpec = which(MetalsTable$Metal[iMetal] == SpecName)
          newSpecNum = newSpecNum + nTG
          # SpecCharge[newSpecNum] = -3L + SpecCharge[iMetalSpec]
          SpecName[newSpecNum] = paste0(WHAMprefix[OMType], TridentTable$FullyDeprot, "-",
                                        MetalsTable$Metal[iMetal])
          SpecStoich[newSpecNum, 1:NComp] = matrix(SpecStoich[iMetalSpec,], nrow = nTG, ncol = NComp, byrow = T)
          diag(SpecStoich[newSpecNum, newCompNum]) = 1L
          SpecStoich[newSpecNum,iH] = SpecStoich[newSpecNum,iH] -3L
          SpecLogK[newSpecNum] = -1 * as.numeric(MetalsTable[iMetal,pKM_cols[TridentTable$S1Strong1Weak2]] +
                                               MetalsTable[iMetal,pKM_cols[TridentTable$S2Strong1Weak2]] +
                                               MetalsTable[iMetal,pKM_cols[TridentTable$S3Strong1Weak2]])
        }

      }

      startComp = max(newCompNum) + 1
      startDefComp = max(newDefCompNum) + 1
      startSpec = max(newSpecNum) + 1

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
  names(SpecTemp) = SpecName

  # Re-ordering species so components are in front
  reorder = match(c(CompName, SpecName[SpecName %in% CompName == F]), SpecName)
  SpecName = SpecName[reorder]
  SpecMC = SpecMC[reorder]
  SpecActCorr = SpecActCorr[reorder]
  SpecNC = SpecNC[reorder]
  SpecCompList = SpecCompList[reorder,]
  SpecStoich = SpecStoich[reorder,]
  SpecLogK = SpecLogK[reorder]
  SpecDeltaH = SpecDeltaH[reorder]
  SpecTemp = SpecTemp[reorder]

  SpecNC = rowSums(SpecStoich != 0L)
  names(SpecNC) = SpecName
  SpecCompList = t(apply(
    SpecStoich,
    MARGIN = 1,
    FUN = function(X) {
      tmp = sort(which(X != 0L))
      if (length(tmp) < max(SpecNC)){
        tmp = c(tmp, rep(0,max(SpecNC) - length(tmp)))
      }
      return(tmp)
    }))
  rownames(SpecCompList) = SpecName
  # CompNS = colSums(Stoich != 0L)
  # names(CompNS) = CompName
  # SpecList = t(apply(
  #   Stoich,
  #   MARGIN = 2,
  #   FUN = function(X) {
  #     tmp = sort(which(X != 0L))
  #     if (length(tmp) < max(CompNS)){
  #       tmp = c(tmp, rep(0,max(CompNS)-length(tmp)))
  #     }
  #     return(tmp)
  #   }))
  # rownames(SpecList) = CompName

  # Assemble output list
  out = list(

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
    SpecTemp = SpecTemp,

    # Phase List
    PhaseCompList = PhaseCompList,
    PhaseStoich = PhaseStoich,

    # WHAM parameters - to be used later
    wDLF = wDLF,
    wKZED = wKZED,
    wKsel = wKsel,
    wP = wP,
    wRadius = wRadius,
    wMolWt = wMolWt

  )
  return(out)

}
