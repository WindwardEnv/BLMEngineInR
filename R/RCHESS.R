RCHESS = function(DoPartialSteps, QuietFlag, ConvergenceCriteria, MaxIter,
                  NComp, NSpec, NBLMetal,
                  SpecConc, SpecLogK, SpecStoich, SpecCtoM, SpecName,
                  CompType, CompName, TotMoles, TotConc,
                  DoTox, MetalName, MetalComp, BLMetalSpecs, CATarget) {
  # inputs:
  #   - DoPartialSteps: boolean, if TRUE, then will do partial Newton-Raphson steps
  #   - QuietFlag: character, one of "Very Quiet" (only print out when run is done), "Quiet" (print out Obs=iObs), or "Debug" (print out lots of info)
  #   - NComp: integer, number of components
  #   - NSpec: integer, number of species reactions
  #   - SpecConc: Species concentrations
  #   - SpecLogK: log-transformed equilibrium coefficients for each reaction
  #   - SpecStoich: (NSpec x NComp) stoichiometry matrix for each reaction
  #   - SpecCtoM: concentration to mass conversion for each species
  #   - SpecName: names of species
  #   - CompType: component types
  #   - CompName: names of components
  #   - TotMoles: the total moles of each component
  #   - DoTox: boolean
  #   - MetalName: character string, the name of the toxic metal
  #   - MetalComp: integer, the position of the metal in the component arrays (i.e., which is the toxic metal component)
  #   - BLMetalSpecs: integer vector, the positions of the species in the arrays which contribute to toxicity (i.e., which species are the toxic metal bound to the relevant biotic ligand)
  #   - CATarget
  # outputs:
  #   - SpecConc: species concentrations after optimization


  SpecK = 10 ^ SpecLogK

  # CompConc = initialGuess(NComp = NComp, CompName = CompName,
  #                         TotConc = TotConc, # mol / L or mol / kgw
  #                         CompType = CompType)
  CompConc = initialGuessV2(TotConc = TotConc,
                            CompType = CompType,
                            CompName = CompName,
                            SpecK = SpecK,
                            SpecStoich = SpecStoich,
                            NComp = NComp,
                            NSpec = NSpec,
                            SpecName = SpecName)

  LogCompConc = log10(CompConc)

  # Initialize Species Concentrations
  SpecConc = CalcSpecConc(CompConc = CompConc,
                          SpecK = SpecK,
                          SpecStoich = SpecStoich,
                          SpecName = SpecName,
                          NComp = NComp,
                          NSpec = NSpec)
  SpecConc[1:NComp] = CompConc

  # Update Total Concentrations for Fixed Activity & Metal
  for (iComp in which((CompType == "FixedAct") | (DoTox & (CompName == MetalName)))){
    TotMoles[iComp] = sum(SpecStoich[,iComp] * (SpecConc * SpecCtoM))
    TotConc[iComp] = TotMoles[iComp] * SpecCtoM[iComp]
  }

  # Calculate Residuals for the first time
  RR = CalcResidual(
    NComp = NComp,
    NSpec = NSpec,
    SpecConc = SpecConc,
    SpecStoich = SpecStoich,
    TotMoles = TotMoles,
    SpecCtoM = SpecCtoM,
    CompName = CompName,
    CompType = CompType,
    MetalComp = MetalComp,
    NBLMetal = NBLMetal,
    BLMetalSpecs = BLMetalSpecs,
    CATarget = CATarget,
    DoTox = DoTox
  )
  RRR = RCalcResidual(
    NComp = NComp,
    NSpec = NSpec,
    SpecConc = SpecConc,
    SpecStoich = SpecStoich,
    TotMoles = TotMoles,
    SpecCtoM = SpecCtoM,
    CompName = CompName,
    CompType = CompType,
    MetalComp = MetalComp,
    BLMetalSpecs = BLMetalSpecs,
    CATarget = CATarget,
    DoTox = DoTox
  )
  Resid = RR$Resid
  MaxError = RR$MaxError
  WhichMax = RR$WhichMax
  CalcTotConc = RR$CalcTotConc
  CalcTotMoles = RR$CalcTotMoles

  # Begin iterating
  Iter = 0
  while ((MaxError > ConvergenceCriteria) & (Iter <= MaxIter)){

    Iter = Iter + 1
    (MaxError_Last = MaxError)

    (JacobianMatrix = Jacobian(NComp = NComp, NSpec = NSpec, SpecStoich = SpecStoich,
                               SpecConc = SpecConc, SpecCtoM = SpecCtoM, CompName = CompName,
                               MetalComp = MetalComp, NBLMetal = NBLMetal,
                               BLMetalSpecs = BLMetalSpecs, DoTox = DoTox))
    # (RJacobianMatrix = RJacobian(NComp = NComp, NSpec = NSpec, SpecStoich = SpecStoich,
    #                              SpecConc = SpecConc, SpecCtoM = SpecCtoM, CompName = CompName,
    #                              MetalComp = MetalComp, NBLMetal = NBLMetal,
    #                              BLMetalSpecs = BLMetalSpecs, DoTox = DoTox))

    (CompConcStep = CalcStep(JacobianMatrix = JacobianMatrix, Resid = Resid,
                              NComp = NComp, CompType = CompType, CompName=CompName))
    # (RCompConcStep = RCalcStep(JacobianMatrix = JacobianMatrix, Resid = Resid,
    #                           NComp = NComp, CompType = CompType, CompName=CompName))

    (CompConc = CompUpdate(NComp = NComp, CompConcStep = CompConcStep, CompConc = SpecConc[1:NComp],
                            CompName = CompName))
    # (RCompConc = RCompUpdate(CompConcStep = CompConcStep, CompConc = SpecConc[1:NComp],
    #                         CompName = CompName))

    SpecConc = CalcSpecConc(CompConc = CompConc,
                             SpecK = SpecK,
                             SpecStoich = SpecStoich,
                             NComp = NComp,
                             NSpec = NSpec,
                             SpecName = SpecName)
    for (iComp in which(CompType == "FixedAct")){
      TotMoles[iComp] = sum(SpecStoich[,iComp] * (SpecConc * SpecCtoM))
      TotConc[iComp] = sum(SpecStoich[,iComp] * SpecConc)
    }

    RR = RCalcResidual(NComp = NComp, NSpec = NSpec,
                       SpecConc = SpecConc,
                       SpecStoich = SpecStoich,
                       TotMoles = TotMoles,
                       SpecCtoM = SpecCtoM,
                       CompName = CompName,
                       CompType = CompType,
                       MetalComp = MetalComp,
                       BLMetalSpecs = BLMetalSpecs,
                       CATarget = CATarget,
                       DoTox = DoTox)

    Resid = RR$Resid
    WhichMax = RR$WhichMax
    MaxError = RR$MaxError
    CalcTotConc = RR$CalcTotConc
    CalcTotMoles = RR$CalcTotMoles

    if (DoPartialSteps) {
      # Full Step
      (CompConc_Full = CompUpdate(NComp = NComp,
                                  CompConcStep = CompConcStep,
                                  CompConc = SpecConc[1:NComp],
                                  CompName = CompName))

      # SpecConc_Full = 10^CppCalcLogSpecConc(LogCompConc = log10(CompConc_Full),
      #                                       SpecLogK = SpecLogK,
      #                                       SpecStoich = SpecStoich,
      #                                       NComp = NComp,
      #                                       NSpec = NSpec)
      SpecConc_Full = CalcSpecConc(CompConc = CompConc_Full,
                                   SpecK = SpecK,
                                   SpecStoich = SpecStoich,
                                   NComp = NComp,
                                   NSpec = NSpec,
                                   SpecName = SpecName)
      TotConc_Full = TotConc
      TotMoles_Full = TotMoles
      for (iComp in which(CompType == "FixedAct")){
        TotMoles_Full[iComp] = sum(SpecStoich[,iComp] * (SpecConc_Full * SpecCtoM))
        TotConc_Full[iComp] = sum(SpecStoich[,iComp] * SpecConc_Full)
      }

      RR_Full = RCalcResidual(NComp = NComp, NSpec = NSpec,
                              SpecConc = SpecConc_Full,
                              SpecStoich = SpecStoich,
                              TotMoles = TotMoles_Full,
                              SpecCtoM = SpecCtoM,
                              CompName = CompName,
                              CompType = CompType,
                              # CompSiteDens = CompSiteDens,
                              MetalComp = MetalComp,
                              BLMetalSpecs = BLMetalSpecs,
                              CATarget = CATarget,
                              DoTox = DoTox)
      best_MaxError = 1L
      if (RR_Full$MaxError > MaxError_Last){
        # Half Step
        CompConc_Half = CompUpdate(NComp = NComp,
                                   CompConcStep = CompConcStep * 0.5,
                                   CompConc = SpecConc[1:NComp],
                                   CompName = CompName)

        # SpecConc_Half = 10^CppCalcLogSpecConc(LogCompConc = log10(CompConc_Half),
        #                                      SpecLogK = SpecLogK,
        #                                      SpecStoich = SpecStoich,
        #                                      NComp = NComp,
        #                                      NSpec = NSpec)
        SpecConc_Half = CalcSpecConc(CompConc = CompConc_Half,
                                     SpecK = SpecK,
                                     SpecStoich = SpecStoich,
                                     NComp = NComp,
                                     NSpec = NSpec,
                                     SpecName = SpecName)

        TotConc_Half = TotConc
        TotMoles_Half = TotMoles
        for (iComp in which(CompType == "FixedAct")){
          TotMoles_Half[iComp] = sum(SpecStoich[,iComp] * (SpecConc_Half * SpecCtoM))
          TotConc_Half[iComp] = sum(SpecStoich[,iComp] * SpecConc_Half)
        }

        RR_Half = RCalcResidual(NComp = NComp, NSpec = NSpec,
                                SpecConc = SpecConc_Half,
                                SpecStoich = SpecStoich,
                                TotMoles = TotMoles_Half,
                                SpecCtoM = SpecCtoM,
                                CompType = CompType,
                                CompName = CompName,
                                MetalComp = MetalComp,
                                BLMetalSpecs = BLMetalSpecs,
                                CATarget = CATarget,
                                DoTox = DoTox)

        # Best Step
        step_Best = 1 - RR_Full$MaxError * (1-0.5) / (RR_Full$MaxError - RR_Half$MaxError)
        CompConc_Best = CompUpdate(NComp = NComp,
                                   CompConcStep = CompConcStep * step_Best,
                                   CompName = CompName,
                                   CompConc = SpecConc[1:NComp])

        # SpecConc_Best = 10^CppCalcLogSpecConc(LogCompConc = log10(CompConc_Best),
        #                                       SpecLogK = SpecLogK,
        #                                       SpecStoich = SpecStoich,
        #                                       NComp = NComp,
        #                                       NSpec = NSpec)
        SpecConc_Best = CalcSpecConc(CompConc = CompConc_Best,
                                     SpecK = SpecK,
                                     SpecStoich = SpecStoich,
                                     NComp = NComp,
                                     NSpec = NSpec,
                                     SpecName = SpecName)

        TotConc_Best = TotConc
        TotMoles_Best = TotMoles
        for (iComp in which(CompType == "FixedAct")){
          TotConc_Best[iComp] = sum(SpecStoich[,iComp] * (SpecConc_Best * SpecCtoM))
          TotConc_Best[iComp] = sum(SpecStoich[,iComp] * SpecConc_Best)
        }

        RR_Best = RCalcResidual(NComp = NComp, NSpec = NSpec,
                                SpecConc = SpecConc_Best,
                                SpecStoich = SpecStoich,
                                TotMoles = TotMoles_Best,
                                SpecCtoM = SpecCtoM,
                                CompType = CompType,
                                CompName = CompName,
                                MetalComp = MetalComp,
                                BLMetalSpecs = BLMetalSpecs,
                                CATarget = CATarget,
                                DoTox = DoTox)

        best_MaxError = which.min(c(RR_Full$MaxError, RR_Half$MaxError, RR_Best$MaxError))

      } else{
        best_MaxError = 1L
      }

      if (best_MaxError == 1L){
        CompConc = CompConc_Full
        SpecConc = SpecConc_Full
        Resid = RR_Full$Resid
        WhichMax = RR_Full$WhichMax
        MaxError = RR_Full$MaxError
        CalcTotConc = RR_Full$CalcTotConc
        CalcTotMoles = RR_Full$CalcTotMoles
        # print("Full Step")
      } else if (best_MaxError == 2L){
        CompConc = CompConc_Half
        SpecConc = SpecConc_Half
        Resid = RR_Half$Resid
        WhichMax = RR_Half$WhichMax
        MaxError = RR_Half$MaxError
        CalcTotConc = RR_Half$CalcTotConc
        CalcTotMoles = RR_Half$CalcTotMoles
        # print("Half Step")
      } else if (best_MaxError == 3L){
        CompConc = CompConc_Best
        SpecConc = SpecConc_Best
        Resid = RR_Best$Resid
        WhichMax = RR_Best$WhichMax
        MaxError = RR_Best$MaxError
        CalcTotConc = RR_Best$CalcTotConc
        CalcTotMoles = RR_Best$CalcTotMoles
        print("Best Step")
      }
    }

    LogCompConc = log10(CompConc)

    if(QuietFlag %in% c("Very Quiet","Quiet") == F){
      print(paste0("Iter=",Iter, ", WhichMax=",CompName[WhichMax],
                   ", MaxError=",format(x = MaxError, scientific = T)))
    }

    # print(paste0("Iter=",Iter,
    #              ", WhichMax=",CompName[WhichMax],
    #              ", Resid[WhichMax]= ",signif(Resid[WhichMax], 6),
    #              ", SpecConc[WhichMax]= ",signif(SpecConc[WhichMax],6)))
    # print(paste0("Iter=",Iter, ", Resid[",CompOfIntName,"]= ",Resid[CompOfInt],
    #              ", CalcTotConc[",CompOfIntName,"]=",CalcTotConc[CompOfInt],
    #              ", CompConc[",CompOfIntName,"]=",CompConc[CompOfInt]))


  }

  out = list(SpecConc = SpecConc,
             Iter = Iter,
             MaxError = MaxError,
             CalcTotConc = CalcTotConc)

  return(out)
}

