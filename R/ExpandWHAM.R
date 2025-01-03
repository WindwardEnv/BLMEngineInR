# Copyright 2024 Windward Environmental LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' @title Expand the DOC component into WHAM components
#'
#' @references
#'   Tipping E. (1994). WHAM--A chemical equilibrium model and
#'     computer code for waters, sediments, and soils incorporating a discrete
#'     site/electrostatic model of ion-binding by humic substances. Computers &
#'     Geosciences, vol. 20, iss. 6, pp. 973-1023.
#'
#' @param ThisProblem a list object following the template of `BlankProblem()`
#' @param ThisWHAM a list object following the template of `BlankWHAM()`
#'
#' @keywords internal
#'
#' @noRd
ExpandWHAM = function(ThisProblem,
                      ThisWHAM = ThisProblem$WHAM) {

  # Pull out variables from This WHAM --------------------------
  nMS = nrow(ThisWHAM$MonodentTable)
  nBP = nrow(ThisWHAM$BidentTable)
  nTG = nrow(ThisWHAM$TridentTable)
  nMPFile = nrow(ThisWHAM$MetalsTable)
  nKselFile = nrow(ThisWHAM$SpecKselTable)
  nCOOH = ThisWHAM$nA
  pKA = ThisWHAM$pKA
  pKB = ThisWHAM$pKB
  dpKA = ThisWHAM$dpKA
  dpKB = ThisWHAM$dpKB
  fprB = ThisWHAM$fprB
  fprT = ThisWHAM$fprT

  #Monodentate sites
  if (nMS > 0) {
    MonodentTable = ThisWHAM$MonodentTable
    MonodentTable$FullyProt = paste0(MonodentTable$S, "H")
    MonodentTable$FullyDeprot = paste0(MonodentTable$S)
    MonodentTable$Strong1Weak2 =
      match(x = tolower(substr(trimws(MonodentTable$StrongWeak), 1, 1)),
            table = c("s", "w"))
  } else {
    MonodentTable = data.frame()
  }

  # Bidentate Pairs
  if (nBP > 0) {
    BidentTable = ThisWHAM$BidentTable
    BidentTable$FullyProt = paste0(BidentTable$S1, BidentTable$S2, "H")
    BidentTable$S1Deprot = paste0(BidentTable$S1, "-", BidentTable$S2, "H")
    BidentTable$S2Deprot = paste0(BidentTable$S2, "-", BidentTable$S1, "H")
    BidentTable$FullyDeprot = paste0(BidentTable$S1, BidentTable$S2)
    BidentTable$S1Strong1Weak2 =
      MonodentTable$Strong1Weak2[match(BidentTable$S1, MonodentTable$S)]
    BidentTable$S2Strong1Weak2 =
      MonodentTable$Strong1Weak2[match(BidentTable$S2, MonodentTable$S)]
  } else {
    BidentTable = data.frame()
  }

  # Tridentate Groups
  if (nTG > 0) {
    TridentTable = ThisWHAM$TridentTable
    TridentTable$FullyProt = paste0(TridentTable$S1, TridentTable$S2,
                                    TridentTable$S3, "H")
    TridentTable$S1Deprot = paste0(TridentTable$S1, "-", TridentTable$S2,
                                   TridentTable$S3, "H")
    TridentTable$S2Deprot = paste0(TridentTable$S2, "-", TridentTable$S1,
                                   TridentTable$S3, "H")
    TridentTable$S3Deprot = paste0(TridentTable$S3, "-", TridentTable$S1,
                                   TridentTable$S2, "H")
    TridentTable$S12Deprot = paste0(TridentTable$S1, TridentTable$S2, "-",
                                    TridentTable$S3, "H")
    TridentTable$S13Deprot = paste0(TridentTable$S1, TridentTable$S3, "-",
                                    TridentTable$S2, "H")
    TridentTable$S23Deprot = paste0(TridentTable$S2, TridentTable$S3, "-",
                                    TridentTable$S1, "H")
    TridentTable$FullyDeprot = paste0(TridentTable$S1, TridentTable$S2,
                                      TridentTable$S3)
    TridentTable$S1Strong1Weak2 =
      MonodentTable$Strong1Weak2[match(TridentTable$S1, MonodentTable$S)]
    TridentTable$S2Strong1Weak2 =
      MonodentTable$Strong1Weak2[match(TridentTable$S2, MonodentTable$S)]
    TridentTable$S3Strong1Weak2 =
      MonodentTable$Strong1Weak2[match(TridentTable$S3, MonodentTable$S)]
  } else {
    TridentTable = data.frame()
  }

  # Metals Parameters Table
  if (nMPFile > 0) {
    MetalsTable = ThisWHAM$MetalsTable
    iBool = MetalsTable$Metal %in%
      c(ThisProblem$Comp$Name, ThisProblem$Spec$Name)
    MetalsTable = MetalsTable[iBool, ]
    nMP = nrow(MetalsTable)
    MetalsTable$pKMBHA = 3 * MetalsTable$pKMAHA - 3
    MetalsTable$pKMBFA = 3.96 * MetalsTable$pKMAFA
  } else {
    MetalsTable = data.frame()
    nMP = 0L
  }

  # Non-standard selectivity coefficients
  if (nKselFile > 0) {
    SpecKselTable = ThisWHAM$SpecKselTable
    SpecKselTable =
      SpecKselTable[SpecKselTable$Spec %in% ThisProblem$Spec$Name, ]
    nKsel = nrow(SpecKselTable)
  } else {
    SpecKselTable = data.frame()
    nKsel = 0L
  }

  # Save original copies of arrays -------------------------
  NewProblem = ThisProblem

  # Do the expansion ---------------------------------------

  # Initialize variables
  iH = which(ThisProblem$Comp$Name == "H")

  # Figure out the number of DOC components we're adding, and what fraction
  InVarWHAM = which(grepl("WHAM", ThisProblem$InVar$Type))

  for (iInVar in InVarWHAM) {

    iMass = ThisProblem$InVar$MCR[iInVar]
    PercHAInMC = any((ThisProblem$InVar$Type == "PercHA") &
                       (ThisProblem$InVar$MCR == iMass))
    PercAFAInMC = any((ThisProblem$InVar$Type == "PercAFA") &
                        (ThisProblem$InVar$MCR == iMass))


    # Charged species for Donnan Layer reactions
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

    if (ThisProblem$InVar$Type[iInVar] == "WHAM-HAFA") {
      WHAMFracAdd = c("HA", "FA")
      if (!PercHAInMC) {
        stop("Must have PercHA input variable in mass compartment if ",
             "specifying a WHAM-HAFA input variable.")
      }
    } else {
      if (PercHAInMC) {
        stop("PercHA input variable specified in mass compartment with ",
             "WHAM-HA or WHAM-FA input variable.")
      }
      if (ThisProblem$InVar$Type[iInVar] == "WHAM-HA") {
        WHAMFracAdd = c("HA")
        if (PercAFAInMC) {
          stop("PercAFA input variable specified in mass compartment with ",
               "WHAM-HA input variable.")
        }
      } else if (ThisProblem$InVar$Type[iInVar] == "WHAM-FA") {
        WHAMFracAdd = c("FA")
      }
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
    DonnanMCName = paste0(ThisProblem$Mass$Name[iMass], "_", DonnanCompName)
    NewProblem = AddMassCompartments(ThisProblem = NewProblem,
                                     MassName = DonnanMCName,
                                     MassAmt = 1E-5,
                                     MassUnit = ThisProblem$Mass$Unit[iMass],
                                     InMass = FALSE)
    NewProblem = AddDefComps(ThisProblem = NewProblem,
                             DefCompName = DonnanCompName,
                             # DefCompFromNum = 1.0,
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
    for (OMi in 1:NWHAMFracAdd) {

      OMType = WHAMFracAdd[OMi]
      ColspKM = paste0("pKM", c("A", "B"), OMType)
      OMSpecType = paste0("WHAM", OMType)

      # Donnan Species
      ChargedSpecDF =
        ThisProblem$Spec[match(ChargedSpecName, ThisProblem$Spec$Name), ]
      NewProblem = AddSpecies(
        ThisProblem = NewProblem,
        SpecEquation = paste0(
          DonnanCompName[OMi], "-",
          ChargedSpecDF$Equation, " + ",
          abs(ChargedSpecDF$Charge), " * ",
          DonnanCompName[OMi]
        ),
        SpecMCName = DonnanMCName[OMi],
        SpecType = DonnanCompName[OMi],
        SpecActCorr = ChargedSpecDF$ActCorr,
        SpecLogK = ChargedSpecDonnanLogK[ChargedSpecDF$Name, OMType],
        SpecDeltaH = ChargedSpecDF$DeltaH,
        SpecTempKelvin = ChargedSpecDF$TempKelvin,
        InSpec = FALSE
      )

      # Monodentate sites
      MonodentpKH[MonodentTable$Strong1Weak2 == 1L] =
        pKA[OMType] + dpKA[OMType] *
        (2 * MonodentTable$S[MonodentTable$Strong1Weak2 == 1L] - 5) / 6
      MonodentpKH[MonodentTable$Strong1Weak2 == 2L] =
        pKB[OMType] + dpKB[OMType] *
        (2 * MonodentTable$S[MonodentTable$Strong1Weak2 == 2L] - 13) / 6
      MonodentAbundance = (1 - fprB[OMType] - fprT[OMType]) *
        nCOOH[OMType] / MonodentTable$AbundDenom

      # Components - fully protonated
      NewProblem = AddDefComps(
        ThisProblem = NewProblem,
        DefCompName = paste0(WHAMprefix[OMi], MonodentTable$FullyProt),
        DefCompFromVar = WHAMprefix[OMi],
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
        SpecName = paste0(WHAMprefix[OMi], MonodentTable$FullyDeprot),
        SpecMCName = ThisProblem$Mass$Name[iMass],
        SpecType = OMSpecType,
        SpecActCorr = "None",
        SpecCompNames = as.list(as.data.frame(rbind(
          paste0(WHAMprefix[OMi], MonodentTable$FullyProt),
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
      for (iMetal in 1:nMP) {
        iMetalSpec = which(ThisProblem$Spec$Name == MetalsTable$Metal[iMetal])
        NewProblem = AddSpecies(
          ThisProblem = NewProblem,
          SpecEquation = paste0(
            WHAMprefix[OMi], MonodentTable$FullyDeprot, "-",
            ThisProblem$Spec$Equation[iMetalSpec],
            " -1 * H + 1 * ", WHAMprefix[OMi], MonodentTable$FullyProt
          ),
          SpecMCName = ThisProblem$Mass$Name[iMass],
          SpecType = OMSpecType,
          SpecActCorr = "None",
          SpecLogK = ThisProblem$Spec$LogK[iMetalSpec] -
            MetalsTable[iMetal, ColspKM[MonodentTable$Strong1Weak2]],
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
          DefCompName = paste0(WHAMprefix[OMi], BidentTable$FullyProt),
          DefCompFromVar = WHAMprefix[OMi],
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
            WHAMprefix[OMi], BidentTable$S1Deprot, " = ",
            "1 * ", WHAMprefix[OMi], BidentTable$FullyProt, " -1 * H"
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
            WHAMprefix[OMi], BidentTable$S2Deprot, " = ",
            "1 * ", WHAMprefix[OMi], BidentTable$FullyProt, " -1 * H"
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
            WHAMprefix[OMi], BidentTable$FullyDeprot, " = ",
            "1 * ", WHAMprefix[OMi], BidentTable$FullyProt, " -2 * H"
          ),
          SpecMCName = ThisProblem$Mass$Name[iMass],
          SpecType = OMSpecType,
          SpecActCorr = "None",
          SpecLogK = -1 * (MonodentpKH[BidentTable$S1] +
                             MonodentpKH[BidentTable$S2]),
          SpecDeltaH = 0,
          SpecTempKelvin = 298.15,
          InSpec = FALSE
        )

        # bound to each metal
        for (iMetal in 1:nMP) {
          iMetalSpec = which(ThisProblem$Spec$Name == MetalsTable$Metal[iMetal])
          NewProblem = AddSpecies(
            ThisProblem = NewProblem,
            SpecEquation = paste0(
              WHAMprefix[OMi], BidentTable$FullyDeprot, "-",
              ThisProblem$Spec$Equation[iMetalSpec],
              " -2 * H + 1 * ", WHAMprefix[OMi], BidentTable$FullyProt
            ),
            SpecMCName = ThisProblem$Mass$Name[iMass],
            SpecType = OMSpecType,
            SpecActCorr = "None",
            SpecLogK = ThisProblem$Spec$LogK[iMetalSpec] -
              (MetalsTable[iMetal, ColspKM[BidentTable$S1Strong1Weak2]] +
                 MetalsTable[iMetal, ColspKM[BidentTable$S2Strong1Weak2]]),
            SpecDeltaH = ThisProblem$Spec$DeltaH[iMetalSpec],
            SpecTempKelvin = ThisProblem$Spec$TempKelvin[iMetalSpec],
            InSpec = FALSE
          )

        }
      }

      # Tridentate sites
      if (nTG > 0) {

        # note: the input is in mg C/L, while nCOOH is mols/g HS
        TridentAbundance = fprT[OMType] * nCOOH[OMType] /
          TridentTable$AbundDenom * 2E-3

        # Components - fully protonated
        NewProblem = AddDefComps(
          ThisProblem = NewProblem,
          DefCompName = paste0(WHAMprefix[OMi], TridentTable$FullyProt),
          DefCompFromVar = WHAMprefix[OMi],
          DefCompCharge = 0L,
          DefCompMCName = ThisProblem$Mass$Name[iMass],
          DefCompType = OMSpecType,
          DefCompActCorr = "None",
          DefCompSiteDens = TridentAbundance,
          InDefComp = FALSE
        )

        # - first site deprotonated
        NewProblem = AddSpecies(
          ThisProblem = NewProblem,
          SpecEquation = paste0(
            WHAMprefix[OMi], TridentTable$S1Deprot, " = ",
            "1 * ", WHAMprefix[OMi], TridentTable$FullyProt, " -1 * H"
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
            WHAMprefix[OMi], TridentTable$S2Deprot, " = ",
            "1 * ", WHAMprefix[OMi], TridentTable$FullyProt, " -1 * H"
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
            WHAMprefix[OMi], TridentTable$S3Deprot, " = ",
            "1 * ", WHAMprefix[OMi], TridentTable$FullyProt, " -1 * H"
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
            WHAMprefix[OMi], TridentTable$S12Deprot, " = ",
            "1 * ", WHAMprefix[OMi], TridentTable$FullyProt, " -2 * H"
          ),
          SpecMCName = ThisProblem$Mass$Name[iMass],
          SpecType = OMSpecType,
          SpecActCorr = "None",
          SpecLogK = -1 * (MonodentpKH[TridentTable$S1] +
                             MonodentpKH[TridentTable$S2]),
          SpecDeltaH = 0,
          SpecTempKelvin = 298.15,
          InSpec = FALSE
        )

        # - first & third sites deprotonated
        NewProblem = AddSpecies(
          ThisProblem = NewProblem,
          SpecEquation = paste0(
            WHAMprefix[OMi], TridentTable$S13Deprot, " = ",
            "1 * ", WHAMprefix[OMi], TridentTable$FullyProt, " -2 * H"
          ),
          SpecMCName = ThisProblem$Mass$Name[iMass],
          SpecType = OMSpecType,
          SpecActCorr = "None",
          SpecLogK = -1 * (MonodentpKH[TridentTable$S1] +
                             MonodentpKH[TridentTable$S3]),
          SpecDeltaH = 0,
          SpecTempKelvin = 298.15,
          InSpec = FALSE
        )

        # - second & third sites deprotonated
        NewProblem = AddSpecies(
          ThisProblem = NewProblem,
          SpecEquation = paste0(
            WHAMprefix[OMi], TridentTable$S23Deprot, " = ",
            "1 * ", WHAMprefix[OMi], TridentTable$FullyProt, " -2 * H"
          ),
          SpecMCName = ThisProblem$Mass$Name[iMass],
          SpecType = OMSpecType,
          SpecActCorr = "None",
          SpecLogK = -1 * (MonodentpKH[TridentTable$S2] +
                             MonodentpKH[TridentTable$S2]),
          SpecDeltaH = 0,
          SpecTempKelvin = 298.15,
          InSpec = FALSE
        )

        # - fully deprot
        NewProblem = AddSpecies(
          ThisProblem = NewProblem,
          SpecEquation = paste0(
            WHAMprefix[OMi], TridentTable$FullyDeprot, " = ",
            "1 * ", WHAMprefix[OMi], TridentTable$FullyProt, " -3 * H"
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
        for (iMetal in 1:nMP) {
          iMetalSpec = which(ThisProblem$Spec$Name == MetalsTable$Metal[iMetal])
          NewProblem = AddSpecies(
            ThisProblem = NewProblem,
            SpecEquation = paste0(
              WHAMprefix[OMi], TridentTable$FullyDeprot, "-",
              ThisProblem$Spec$Equation[iMetalSpec],
              " -3 * H + 1 * ", WHAMprefix[OMi], TridentTable$FullyProt
            ),
            SpecMCName = ThisProblem$Mass$Name[iMass],
            SpecType = OMSpecType,
            SpecActCorr = "None",
            SpecLogK = ThisProblem$Spec$LogK[iMetalSpec] -
              (MetalsTable[iMetal, ColspKM[TridentTable$S1Strong1Weak2]] +
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
  NewProblem$DoWHAM = any(grepl("WHAM", NewProblem$Spec$Type))
  NewProblem$WHAM = ThisWHAM

  return(NewProblem)

}
