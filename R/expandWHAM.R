#' Expand the DOC component into WHAM components
#'
#' @param nMC Number of mass compartments
#' @param NComp Number of components
#' @param NSpec Number of species
#' @param MassName Names of the mass compartments
#' @param CompName character vector of length `NComp`; component names
#' @param CompCharge integer vector of length `NComp`; the charge of the
#'   components as free ions
#' @param CompMC integer vector of length `NComp`; Which mass compartment the
#'   component belongs to
#' @param CompType integer vector of length `NComp`; the type of component. It
#'   should be a fixed set of values where 1 =...
#' @param CompActCorr integer vector of length `NComp`; the method to use for
#'   activity corrections with this component, where 1 ==
#' @param CompSiteDens numeric vector of length `NComp`; the density of binding
#'   sites for this component - usually 1, except for DOC and BL components
#' @param SpecName character vector of length `NSpec`; species names
#' @param SpecMC integer vector of length `NSpec`; which mass compartment the
#'   speces belongs to
#' @param SpecType integer vector of length `NSpec`; the type of chemical
#'   species, where 1 = ...
#' @param SpecActCorr integer vector of length `NSpec`; the method to use for
#'   activity corrections with this speies where 1 = ...
#' @param SpecNC integer vector of length `NSpec`; the number of components used
#'   to create a given species
#' @param CompList the list of components used to create a given species
#' @param Stoich the stoichiometry matrix of the reactions (integer matrix of
#'   size `NSpec` x `nComp`)
#' @param LogK the log10-transformed equilibrium coefficients (numeric vector of
#'   length `NSpec`)
#' @param DeltaH the enthalpy change for each formation reaction (species)
#' @param WHAMVer a character string specifying the WHAM version to use, must be
#'   one of `"V"` (default), `"VI"`, or `"VII"`.
#' @param wdatFile (optional) a character string specifying the file path of a
#'   WHAM parameter file
#'
#' @keywords internal
#'
#' @noRd
expandWHAM = function(nMC, NComp, NSpec, MassName,
                      CompName, CompCharge, CompMC, CompType,
                      CompActCorr, CompSiteDens, SpecName, SpecMC,
                      SpecType,  SpecActCorr, SpecNC,
                      CompList, Stoich, LogK, DeltaH,
                      WHAMVer = c("V","VI","VII"), wdatFile=NULL) {

  # Debugging
  if(FALSE){

      Stoich = Full_InorgDataStoich
      K = Full_InorgDataK
      LogK = log10(K)
      nMC = 2
      MassName = c("Water","BL")
      NComp = ncol(Stoich)
      CompName = colnames(Stoich)
      CompCharge = array(as.integer(c(1,2,-1,2,2,1,1,-2,-1,-2,0)), dimnames = list(CompName))
      CompMC = array(c(rep(1L, NComp - 1), 2L), dimnames = list(CompName))
      CompType = array(c(2L, rep(1L,NComp - 2), 11L), dimnames = list(CompName))
      CompActCorr = array(c(2L,2L,1L,rep(2L,NComp-4),1L), dimnames = list(CompName))
      CompSiteDens = array(c(1.0,1.0,0.0006,rep(1.0, NComp-4),3E-5), dimnames = list(CompName))

      NSpec = nrow(Stoich)
      SpecName = rownames(Stoich)
      # SpecCharge = Full_InorgDataCharge#Stoich %*% CompCharge
      SpecMC = array(rep(1L, NSpec), dimnames = list(SpecName))
      SpecType = array(c(CompType, rep(1L,NSpec-NComp)), dimnames = list(SpecName))
      SpecActCorr = array(c(CompActCorr,rep(2L,NSpec-NComp)), dimnames = list(SpecName))
      SpecActCorr[grepl("BL",SpecName)] = 1L
      SpecNC = rowSums(Stoich != 0)
      CompList = t(apply(
        Stoich,
        MARGIN = 1,
        FUN = function(X) {
          tmp = sort(which(X != 0L))
          if (length(tmp) < max(SpecNC)+2) {
            tmp = c(tmp, rep(0, max(SpecNC) + 2 - length(tmp)))
          }
          return(tmp)
        }
      ))
      # CompNS = colSums(Stoich != 0L)
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
      DeltaH = array(rep(0.0, NSpec), dimnames = list(SpecName))
      WHAMVer = "V"

    }

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
  HAFA = c("HA","FA")
  iH = which(CompName == "H")
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
  wNComp = (nMS + nBP + nTG) * 2
  startComp = NComp + 1
  NComp = NComp + wNComp
  wCompName = paste0(rep(HAFA, each = wNComp / 2),
                     rep(c(MonodentTable$FullyProt,
                           BidentTable$FullyProt,
                           TridentTable$FullyProt), times = 2))
  CompName = c(CompName, wCompName)
  CompCharge = c(CompCharge, array(0L, dim = wNComp, dimnames = list(wCompName)))
  CompMC = c(CompMC, array(which(tolower(MassName) == "water"), dim = wNComp, dimnames = list(wCompName)))
  CompType = c(CompType, array(1L, wNComp, dimnames = list(wCompName)))
  CompActCorr = c(CompActCorr, array(0L, wNComp, dimnames = list(wCompName)))
  CompSiteDens = c(CompSiteDens, array(NA, wNComp, dimnames = list(wCompName)))

  wNSpec = (nMS * (2 + nMP) + nBP * (4 + nMP) + nTG * (8 + nMP)) * 2
  startSpec = NSpec + 1
  NSpec = NSpec + wNSpec
  SpecName = c(SpecName, array(paste0("newDOCSpecies",1:wNSpec),wNSpec))
  SpecMC = c(SpecMC, array(which(tolower(MassName) == "water"), dim = wNSpec,dimnames=list(paste0("newDOCSpecies",1:wNSpec))))#this should always be water
  SpecType = c(SpecType, array(1L, wNSpec, dimnames=list(paste0("newDOCSpecies",1:wNSpec))))
  SpecActCorr = c(SpecActCorr, array(0L, wNSpec, dimnames=list(paste0("newDOCSpecies",1:wNSpec))))
  # SpecCharge = c(SpecCharge, array(NA, wNSpec, dimnames=list(paste0("newDOCSpecies",1:wNSpec))))
  SpecNC = c(SpecNC, array(NA, wNSpec, dimnames=list(paste0("newDOCSpecies",1:wNSpec))))
  CompList = rbind(CompList, matrix(0L, nrow = wNSpec, ncol = ncol(CompList)))
  Stoich = rbind(cbind(Stoich,
                       matrix(0L,nrow=NSpec-wNSpec, ncol = wNComp,
                              dimnames = list(SpecName[1:(NSpec-wNSpec)],wCompName))),
                 matrix(0L, nrow = wNSpec, ncol = NComp,
                        dimnames = list(paste0("newDOCSpecies",1:wNSpec),
                                        CompName)))
  LogK = c(LogK, array(NA,wNSpec,dimnames=list(paste0("newDOCSpecies",1:wNSpec))))
  DeltaH = c(DeltaH, array(0,wNSpec,dimnames=list(paste0("newDOCSpecies",1:wNSpec))))

  Monodent_pKH = numeric(nMS)
  Monodent_Abundance = numeric(nMS)
  Bident_Abundance = numeric(nBP)
  Trident_Abundance = numeric(nTG)
  for(OMType in HAFA){

    pKM_cols = paste0("pKM",c("A","B"),OMType)

    # Monodentate sites
    newCompNum = startComp:(startComp + nMS - 1)
    Monodent_pKH[1:4] = pKHA[OMType] + dpKHA[OMType] * (2 * MonodentTable$S[1:4] - 5) / 6
    Monodent_pKH[5:8] = pKHB[OMType] + dpKHB[OMType] * (2 * MonodentTable$S[5:8] - 13) / 6
    Monodent_Abundance = (1 - fprB[OMType] - fprT[OMType]) * nCOOH[OMType] / MonodentTable$AbundDenom
    CompSiteDens[newCompNum] = Monodent_Abundance * 1E-3 # the input is in milligrams, while nCOOH is mols/g

    # - fully protonated (components)
    newSpecNum = startSpec:(startSpec + nMS - 1)
    SpecName[newSpecNum] = paste0(OMType, MonodentTable$FullyProt)
    # SpecCharge[newSpecNum] = 0L
    diag(Stoich[newSpecNum, newCompNum]) = 1L
    LogK[newSpecNum] = 0.0

    # - fully deprot
    newSpecNum = newSpecNum + nMS
    SpecName[newSpecNum] = paste0(OMType, MonodentTable$FullyDeprot)
    # SpecCharge[newSpecNum] = -1L
    diag(Stoich[newSpecNum, newCompNum]) = 1L
    Stoich[newSpecNum,iH] = -1L
    LogK[newSpecNum] = -1 * Monodent_pKH

    # bound to each metal
    for(iMetal in 1:nMP){
      iMetalSpec = which(MetalsTable$Metal[iMetal] == SpecName)
      newSpecNum = newSpecNum + nMS
      # SpecCharge[newSpecNum] = -1L + SpecCharge[iMetalSpec]
      SpecName[newSpecNum] = paste0(OMType, MonodentTable$FullyDeprot, "-",
                                     MetalsTable$Metal[iMetal])
      Stoich[newSpecNum, 1:NComp] = matrix(Stoich[iMetalSpec,], nrow = nMS, ncol = NComp, byrow = T)
      diag(Stoich[newSpecNum, newCompNum]) = 1L
      Stoich[newSpecNum,iH] = Stoich[newSpecNum,iH] -1L
      LogK[newSpecNum] = -1 * as.numeric(MetalsTable[iMetal,pKM_cols[MonodentTable$Strong1Weak2]])
    }


    # Bidentate sites
    if (nBP > 0){
      startComp = max(newCompNum) + 1
      startSpec = max(newSpecNum) + 1
      newCompNum = startComp:(startComp + nBP - 1)
      Bident_Abundance = fprB[OMType] * nCOOH[OMType] / BidentTable$AbundDenom
      CompSiteDens[newCompNum] = Bident_Abundance * 1E-3 # the input is in milligrams, while nCOOH is mols/g

      # - fully protonated
      newSpecNum = startSpec:(startSpec + nBP - 1)
      SpecName[newSpecNum] = paste0(OMType, BidentTable$FullyProt)
      # SpecCharge[newSpecNum] = 0L
      diag(Stoich[newSpecNum, newCompNum]) = 1L
      LogK[newSpecNum] = 0.0

      # - first site deprotonated
      newSpecNum = newSpecNum + nBP
      SpecName[newSpecNum] = paste0(OMType, BidentTable$S1Deprot)
      # SpecCharge[newSpecNum] = -1L
      diag(Stoich[newSpecNum, newCompNum]) = 1L
      Stoich[newSpecNum,iH] = -1L
      LogK[newSpecNum] = -1 * Monodent_pKH[BidentTable$S1]

      # - second site deprotonated
      newSpecNum = newSpecNum + nBP
      SpecName[newSpecNum] = paste0(OMType, BidentTable$S2Deprot)
      # SpecCharge[newSpecNum] = -1L
      diag(Stoich[newSpecNum, newCompNum]) = 1L
      Stoich[newSpecNum,iH] = -1L
      LogK[newSpecNum] = -1 * Monodent_pKH[BidentTable$S2]

      # - fully deprot
      newSpecNum = newSpecNum + nBP
      SpecName[newSpecNum] = paste0(OMType, BidentTable$FullyDeprot)
      # SpecCharge[newSpecNum] = -2L
      diag(Stoich[newSpecNum, newCompNum]) = 1L
      Stoich[newSpecNum,iH] = -2L
      LogK[newSpecNum] = -1 * (Monodent_pKH[BidentTable$S1] +
                                 Monodent_pKH[BidentTable$S2])

      # bound to each metal
      for(iMetal in 1:nMP){
        iMetalSpec = which(MetalsTable$Metal[iMetal] == SpecName)
        newSpecNum = newSpecNum + nBP
        # SpecCharge[newSpecNum] = -2L + SpecCharge[iMetalSpec]
        SpecName[newSpecNum] = paste0(OMType, BidentTable$FullyDeprot, "-",
                                      MetalsTable$Metal[iMetal])
        Stoich[newSpecNum, 1:NComp] = matrix(Stoich[iMetalSpec,], nrow = nBP, ncol = NComp, byrow = T)
        diag(Stoich[newSpecNum, newCompNum]) = 1L
        Stoich[newSpecNum,iH] = Stoich[newSpecNum,iH] -2L
        LogK[newSpecNum] = -1 * as.numeric(MetalsTable[iMetal,pKM_cols[BidentTable$S1Strong1Weak2]] +
                                             MetalsTable[iMetal,pKM_cols[BidentTable$S2Strong1Weak2]])
      }

    }

    # Tridentate sites
    if (nTG > 0){
      startComp = max(newCompNum) + 1
      startSpec = max(newSpecNum) + 1
      newCompNum = startComp:(startComp + nTG - 1)
      Trident_Abundance = fprT[OMType] * nCOOH[OMType] / TridentTable$AbundDenom
      CompSiteDens[newCompNum] = Trident_Abundance * 1E-3 # the input is in milligrams, while nCOOH is mols/g

      # - fully protonated
      newSpecNum = startSpec:(startSpec + nTG - 1)
      SpecName[newSpecNum] = paste0(OMType, TridentTable$FullyProt)
      # SpecCharge[newSpecNum] = 0L
      diag(Stoich[newSpecNum, newCompNum]) = 1L
      LogK[newSpecNum] = 0.0

      # - first site deprotonated
      newSpecNum = newSpecNum + nTG
      SpecName[newSpecNum] = paste0(OMType, TridentTable$S1Deprot)
      # SpecCharge[newSpecNum] = -1L
      diag(Stoich[newSpecNum, newCompNum]) = 1L
      Stoich[newSpecNum,iH] = -1L
      LogK[newSpecNum] = -1 * Monodent_pKH[TridentTable$S1]

      # - second site deprotonated
      newSpecNum = newSpecNum + nTG
      SpecName[newSpecNum] = paste0(OMType, TridentTable$S2Deprot)
      # SpecCharge[newSpecNum] = -1L
      diag(Stoich[newSpecNum, newCompNum]) = 1L
      Stoich[newSpecNum,iH] = -1L
      LogK[newSpecNum] = -1 * Monodent_pKH[TridentTable$S2]

      # - third site deprotonated
      newSpecNum = newSpecNum + nTG
      SpecName[newSpecNum] = paste0(OMType, TridentTable$S3Deprot)
      # SpecCharge[newSpecNum] = -1L
      diag(Stoich[newSpecNum, newCompNum]) = 1L
      Stoich[newSpecNum,iH] = -1L
      LogK[newSpecNum] = -1 * Monodent_pKH[TridentTable$S3]

      # - first & second sites deprotonated
      newSpecNum = newSpecNum + nTG
      SpecName[newSpecNum] = paste0(OMType, TridentTable$S12Deprot)
      # SpecCharge[newSpecNum] = -2L
      diag(Stoich[newSpecNum, newCompNum]) = 1L
      Stoich[newSpecNum,iH] = -2L
      LogK[newSpecNum] = -1 * (Monodent_pKH[TridentTable$S1] + Monodent_pKH[TridentTable$S2])

      # - first & third sites deprotonated
      newSpecNum = newSpecNum + nTG
      SpecName[newSpecNum] = paste0(OMType, TridentTable$S13Deprot)
      # SpecCharge[newSpecNum] = -2L
      diag(Stoich[newSpecNum, newCompNum]) = 1L
      Stoich[newSpecNum,iH] = -2L
      LogK[newSpecNum] = -1 * (Monodent_pKH[TridentTable$S1] + Monodent_pKH[TridentTable$S3])

      # - second & third sites deprotonated
      newSpecNum = newSpecNum + nTG
      SpecName[newSpecNum] = paste0(OMType, TridentTable$S23Deprot)
      # SpecCharge[newSpecNum] = -2L
      diag(Stoich[newSpecNum, newCompNum]) = 1L
      Stoich[newSpecNum,iH] = -2L
      LogK[newSpecNum] = -1 * (Monodent_pKH[TridentTable$S2] + Monodent_pKH[TridentTable$S3])

      # - fully deprot
      newSpecNum = newSpecNum + nTG
      SpecName[newSpecNum] = paste0(OMType, TridentTable$FullyDeprot)
      # SpecCharge[newSpecNum] = -3L
      diag(Stoich[newSpecNum, newCompNum]) = 1L
      Stoich[newSpecNum,iH] = -3L
      LogK[newSpecNum] = -1 *(Monodent_pKH[TridentTable$S1] +
                                Monodent_pKH[TridentTable$S2] +
                                Monodent_pKH[TridentTable$S3])

      # bound to each metal
      for(iMetal in 1:nMP){
        iMetalSpec = which(MetalsTable$Metal[iMetal] == SpecName)
        newSpecNum = newSpecNum + nTG
        # SpecCharge[newSpecNum] = -3L + SpecCharge[iMetalSpec]
        SpecName[newSpecNum] = paste0(OMType, TridentTable$FullyDeprot, "-",
                                      MetalsTable$Metal[iMetal])
        Stoich[newSpecNum, 1:NComp] = matrix(Stoich[iMetalSpec,], nrow = nTG, ncol = NComp, byrow = T)
        diag(Stoich[newSpecNum, newCompNum]) = 1L
        Stoich[newSpecNum,iH] = Stoich[newSpecNum,iH] -3L
        LogK[newSpecNum] = -1 * as.numeric(MetalsTable[iMetal,pKM_cols[TridentTable$S1Strong1Weak2]] +
                                             MetalsTable[iMetal,pKM_cols[TridentTable$S2Strong1Weak2]] +
                                             MetalsTable[iMetal,pKM_cols[TridentTable$S3Strong1Weak2]])
      }

    }

    startComp = max(newCompNum) + 1
    startSpec = max(newSpecNum) + 1

  }

  # Cleanup
  SpecNC = rowSums(Stoich != 0L)
  names(SpecNC) = SpecName
  CompList = t(apply(
    Stoich,
    MARGIN = 1,
    FUN = function(X) {
      tmp = sort(which(X != 0L))
      if (length(tmp) < max(SpecNC)){
        tmp = c(tmp, rep(0,max(SpecNC) - length(tmp)))
      }
      return(tmp)
    }))
  rownames(CompList) = SpecName
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
  # names(SpecCharge) = SpecName
  names(SpecMC) = SpecName
  names(SpecType) = SpecName
  names(SpecActCorr) = SpecName
  rownames(Stoich) = SpecName
  names(LogK) = SpecName
  names(DeltaH) = SpecName

  # every input should be returned, other than the WHAM type and other info
  out = list(
    NComp = NComp,
    NSpec = NSpec,
    CompName = CompName,
    CompCharge = CompCharge,
    CompMC = CompMC,
    CompType = CompType,
    CompActCorr = CompActCorr,
    CompSiteDens = CompSiteDens,
    SpecName = SpecName,
    # SpecCharge = SpecCharge,
    SpecMC = SpecMC,
    SpecType = SpecType,
    SpecActCorr = SpecActCorr,
    SpecNC = SpecNC,
    CompList = CompList,
    Stoich = Stoich,
    LogK = LogK,
    DeltaH = DeltaH,
    wDLF = wDLF,
    wKZED = wKZED,
    wKsel = wKsel,
    wP = wP,
    wRadius = wRadius,
    wMolWt = wMolWt
  )
  return(out)

}
