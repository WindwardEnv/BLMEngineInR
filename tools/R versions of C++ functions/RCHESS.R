#' @title CHemical Equilibria in Soils and Solutions
#'
#' @description Given a chemical system, equilibria equations, and total
#'   concentrations of components, calculate the species concentrations of each
#'   chemical product in the system.
#'
#' @param DoPartialSteps logical, TRUE = do partial Newton-Raphson steps when a
#'   full step would send the MaxError higher; FALSE = only do full steps
#' @param QuietFlag character, one of "Very Quiet" (only print Out when run is
#'   done), "Quiet" (print Out Obs=iObs), or "Debug" (print Out lots of info)
#' @param ConvergenceCriteria numeric, the maximum value of MaxError that counts
#'   as convergence by the Newton-Raphson root-finding algorithm
#' @param MaxIter integer, the maximum number of iterations the Newton-Raphson
#'   root-finding algorithm should do before giving up
#' @param NComp integer, number of components
#' @param NSpec integer, number of species reactions
#' @param NBLMetal integer, the number of biotic ligand-bound metal species that
#'   are associated with toxic effects.
#' @param SpecK numeric vector (NSpec), the equilibrium coefficient of the
#'   formation reactions
#' @param SpecStoich signed integer matrix (NSpec x NComp), the reaction
#'   stoichiometry of the formation reactions
#' @param SpecCtoM numeric vector (NSpec), the concentration to mass conversion
#'   factor of the chemical species for which we have formation reactions
#' @param SpecName character vector (NSpec), the name of the chemical species
#'   for which we have formation reactions
#' @param CompType character vector (NComp), the type of each component in the
#'   simulation
#' @param CompName character vector (NComp), the name of each component in the
#'   simulation
#' @param TotMoles numeric vector (NComp), the total moles of each component in
#'   the simulation (units of mol)
#' @param TotConc numeric vector (NComp), the total concentrations of each
#'   component in the simulation (units of e.g., mol/L and mol/kg)
#' @param DoTox logical, TRUE for toxicity mode where the MetalName component
#'   concentration is adjusted to try to match the CATarget with BLMetalSpecs
#' @param MetalName character string, the name of the toxic metal
#' @param MetalComp integer, the position of the metal in the component arrays
#'   (i.e., which is the toxic metal component) Note: this are base-1 indexed.
#' @param BLMetalSpecs integer vector, the positions of the species in the
#'   arrays which contribute to toxicity (i.e., which species are the toxic
#'   metal bound to the relevant biotic ligand) Note: these are base-1 indexed.
#' @param CATarget numeric, the target critical accumulation in units of mol /
#'   kg (only used when DoTox == TRUE)
#'
#' @return list with the following elements:
#'   \desribe{
#'     \item{SpecConc}{numeric vector (NSpec), the concentrations of each
#'       species for which we have formation reactions}
#'     \item{Iter}{integer, the number of Newton-Raphson iterations that we
#'       needed to reach convergence}
#'     \item{MaxError}{numeric, the highest final absolute error fraction
#'       =max(abs(Resid / TotMoles))}
#'     \item{CalcTotConc}{numeric vector (NComp), the calculated total
#'       concentrations of each component in the simulation (units of e.g.,
#'       mol/L and mol/kg)}
#'   }
#' @export
#'
RCHESS = function(DoPartialSteps, QuietFlag, ConvergenceCriteria, MaxIter,
                  NComp, NSpec, NBLMetal,
                  SpecK, SpecStoich, SpecCtoM, SpecName,
                  CompType, CompName, TotMoles, TotConc,
                  DoTox, MetalName, MetalComp, BLMetalSpecs, CATarget) {
  # Outputs:
  #   - SpecConc: species concentrations after optimization

  CompConc = InitialGuess(TotConc = TotConc,
                          CompType = CompType,
                          # CompName = CompName,
                          SpecK = SpecK,
                          SpecStoich = SpecStoich,
                          SpecName = SpecName,
                          NComp = NComp,
                          NSpec = NSpec)

  # Initialize Species Concentrations
  SpecConc = CalcSpecConc(CompConc = CompConc,
                          SpecK = SpecK,
                          SpecStoich = SpecStoich,
                          SpecName = SpecName,
                          NComp = NComp,
                          NSpec = NSpec)
  SpecConc[1:NComp] = CompConc

  # Update Total Concentrations for Fixed Activity & Metal
  UpdateResult = UpdateTotalsList(NComp = NComp,
                                  NSpec = NSpec,
                                  CompType = CompType,
                                  CompName = CompName,
                                  MetalName = MetalName,
                                  TotMoles = TotMoles,
                                  SpecStoich = SpecStoich,
                                  SpecMoles = SpecConc * SpecCtoM,
                                  TotConc = TotConc,
                                  SpecCtoM = SpecCtoM,
                                  DoTox = DoTox)
  TotMoles = UpdateResult[["TotMoles"]]
  TotConc = UpdateResult[["TotConc"]]

  # Calculate Residuals for the first time
  ResidResults = CalcResidual(
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
  Resid = ResidResults$Resid
  MaxError = ResidResults$MaxError
  WhichMax = ResidResults$WhichMax
  CalcTotConc = ResidResults$CalcTotConc
  CalcTotMoles = ResidResults$CalcTotMoles

  # Begin iterating
  Iter = 0
  while ((MaxError > ConvergenceCriteria) && (Iter <= MaxIter)) {

    Iter = Iter + 1
    (MaxError_Last = MaxError)

    (JacobianMatrix = Jacobian(NComp = NComp, NSpec = NSpec,
                               SpecStoich = SpecStoich,
                               SpecConc = SpecConc, SpecCtoM = SpecCtoM,
                               CompName = CompName,
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
    for (iComp in which(CompType == "FixedConc")){
      TotMoles[iComp] = sum(SpecStoich[, iComp] * (SpecConc * SpecCtoM))
      TotConc[iComp] = sum(SpecStoich[, iComp] * SpecConc)
    }

    ResidResults = CalcResidual(NComp = NComp, NSpec = NSpec,
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
                                DoTox = DoTox)

    Resid = ResidResults$Resid
    WhichMax = ResidResults$WhichMax
    MaxError = ResidResults$MaxError
    CalcTotConc = ResidResults$CalcTotConc
    CalcTotMoles = ResidResults$CalcTotMoles

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
      for (iComp in which(CompType == "FixedConc")){
        TotMoles_Full[iComp] = sum(SpecStoich[,iComp] * (SpecConc_Full * SpecCtoM))
        TotConc_Full[iComp] = sum(SpecStoich[,iComp] * SpecConc_Full)
      }

      ResidResults_Full = CalcResidual(NComp = NComp, NSpec = NSpec,
                              SpecConc = SpecConc_Full,
                              SpecStoich = SpecStoich,
                              TotMoles = TotMoles_Full,
                              SpecCtoM = SpecCtoM,
                              CompName = CompName,
                              CompType = CompType,
                              # CompSiteDens = CompSiteDens,
                              MetalComp = MetalComp,
                              NBLMetal = NBLMetal,
                              BLMetalSpecs = BLMetalSpecs,
                              CATarget = CATarget,
                              DoTox = DoTox)
      best_MaxError = 1L
      if (ResidResults_Full$MaxError > MaxError_Last){
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
        for (iComp in which(CompType == "FixedConc")){
          TotMoles_Half[iComp] = sum(SpecStoich[,iComp] * (SpecConc_Half * SpecCtoM))
          TotConc_Half[iComp] = sum(SpecStoich[,iComp] * SpecConc_Half)
        }

        ResidResults_Half = CalcResidual(NComp = NComp, NSpec = NSpec,
                                SpecConc = SpecConc_Half,
                                SpecStoich = SpecStoich,
                                TotMoles = TotMoles_Half,
                                SpecCtoM = SpecCtoM,
                                CompType = CompType,
                                CompName = CompName,
                                MetalComp = MetalComp,
                                NBLMetal = NBLMetal,
                                BLMetalSpecs = BLMetalSpecs,
                                CATarget = CATarget,
                                DoTox = DoTox)

        # Best Step
        step_Best = 1 - ResidResults_Full$MaxError * (1-0.5) / (ResidResults_Full$MaxError - ResidResults_Half$MaxError)
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
        for (iComp in which(CompType == "FixedConc")){
          TotConc_Best[iComp] = sum(SpecStoich[,iComp] * (SpecConc_Best * SpecCtoM))
          TotConc_Best[iComp] = sum(SpecStoich[,iComp] * SpecConc_Best)
        }

        ResidResults_Best = CalcResidual(NComp = NComp, NSpec = NSpec,
                                SpecConc = SpecConc_Best,
                                SpecStoich = SpecStoich,
                                TotMoles = TotMoles_Best,
                                SpecCtoM = SpecCtoM,
                                CompType = CompType,
                                CompName = CompName,
                                MetalComp = MetalComp,
                                NBLMetal = NBLMetal,
                                BLMetalSpecs = BLMetalSpecs,
                                CATarget = CATarget,
                                DoTox = DoTox)

        best_MaxError = which.min(c(ResidResults_Full$MaxError, ResidResults_Half$MaxError, ResidResults_Best$MaxError))

      } else{
        best_MaxError = 1L
      }

      if (best_MaxError == 1L){
        CompConc = CompConc_Full
        SpecConc = SpecConc_Full
        Resid = ResidResults_Full$Resid
        WhichMax = ResidResults_Full$WhichMax
        MaxError = ResidResults_Full$MaxError
        CalcTotConc = ResidResults_Full$CalcTotConc
        CalcTotMoles = ResidResults_Full$CalcTotMoles
        # print("Full Step")
      } else if (best_MaxError == 2L){
        CompConc = CompConc_Half
        SpecConc = SpecConc_Half
        Resid = ResidResults_Half$Resid
        WhichMax = ResidResults_Half$WhichMax
        MaxError = ResidResults_Half$MaxError
        CalcTotConc = ResidResults_Half$CalcTotConc
        CalcTotMoles = ResidResults_Half$CalcTotMoles
        # print("Half Step")
      } else if (best_MaxError == 3L){
        CompConc = CompConc_Best
        SpecConc = SpecConc_Best
        Resid = ResidResults_Best$Resid
        WhichMax = ResidResults_Best$WhichMax
        MaxError = ResidResults_Best$MaxError
        CalcTotConc = ResidResults_Best$CalcTotConc
        CalcTotMoles = ResidResults_Best$CalcTotMoles
        print("Best Step")
      }
    }

    if(QuietFlag %in% c("Very Quiet","Quiet") == FALSE){
      print(paste0("Iter=",Iter, ", WhichMax=",CompName[WhichMax],
                   ", MaxError=",format(x = MaxError, scientific = TRUE)))
    }

    # print(paste0("Iter=",Iter,
    #              ", WhichMax=",CompName[WhichMax],
    #              ", Resid[WhichMax]= ",signif(Resid[WhichMax], 6),
    #              ", SpecConc[WhichMax]= ",signif(SpecConc[WhichMax],6)))
    # print(paste0("Iter=",Iter, ", Resid[",CompOfIntName,"]= ",Resid[CompOfInt],
    #              ", CalcTotConc[",CompOfIntName,"]=",CalcTotConc[CompOfInt],
    #              ", CompConc[",CompOfIntName,"]=",CompConc[CompOfInt]))


  }

  Out = list(SpecConc = SpecConc,
             Iter = Iter,
             MaxError = MaxError,
             CalcTotConc = CalcTotConc)

  return(Out)
}

