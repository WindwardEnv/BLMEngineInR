RCalcSpecConc = function(CompConc, SpecK, SpecStoich, NComp = length(CompConc),
                         NSpec = length(SpecK), SpecName = names(SpecK)){

  SpecConc = rep(1, NSpec)
  for (iSpec in 1:NSpec){
    CompConcStep = CompConc ^ SpecStoich[iSpec,]
    SpecConc[iSpec] = prod(CompConcStep) * SpecK[iSpec]
  }

  names(SpecConc) = SpecName
  return(SpecConc)
}



#' @title Iterative step improvement in component concentrations
#'
#' @description
#' RCompUpdate calculates an iterative improvement on the component concentrations based on the Newton-Raphson solution from the current iteration.
#'
#' @details
#' If the iteration would cause the adjusted component concentrations to be less than zero, then the component concentration is simply divided by 10 for this iteration.
#'
#'
#' @param CompConcStep : Vector (NComp) of adjustments to the component concentrations
#' @param CompConc : Vector (NComp) of component concentrations, input values are from this iteration
#' @param CompCtoM >>>>>>>>> I think we decided we don't need this <<<<<<<<<<<
#' @param CompName : Vector (NComp) with the names of the components
#'
#' @return  Vector CompConc (NComp) modified for the next iteration

RCompUpdate = function(CompConcStep, CompConc, CompCtoM, CompName){

  (oldCompConc = CompConc)
  # CompConc[CompConcStep < oldCompConc] = (oldCompConc - CompConcStep)[CompConcStep < oldCompConc]
  # CompConc = (oldCompConc - CompConcStep / CompCtoM)
  CompConc = (oldCompConc - CompConcStep)
  # CompConc[CompConcStep >= oldCompConc] = oldCompConc[CompConcStep >= oldCompConc] / 10
  CompConc[CompConc <= 0] = oldCompConc[CompConc <= 0] / 10


  # for (iComp in 1:NComp){
  #   if (CompConcStep[iComp] >= oldCompConc[iComp]) {
  #     CompConc[iComp] = oldCompConc[iComp] / 10
  #   } else {
  #     CompConc[iComp] = CompConc[iComp] - CompConcStep[iComp]
  #   }
  #   if (CompConc[iComp] <= 0.0) {
  #     CompConc[iComp] = 1E-10
  #   }
  # }#NEXT cc

  names(CompConc) = CompName
  return(CompConc)

}


RCHESSIter <- function(DoPartialSteps, QuietFlag,
                       NComp, NSpec,
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

  return()
}



start.time = Sys.time()
# Get the problem set up
paramFile = "scrap/parameter file format/full_organic.dat4"
inputFile = "scrap/parameter file format/full_organic.blm4"
QuietFlag = c("Very Quiet","Quiet","Debug")[3]
iCA = 1
DoTox = T
DoPartialSteps = F
ConvergenceCriteria = 0.001
MaxIter = 30L
{

  thisProblem = defineProblem(paramFile = paramFile)
  allInput = do.call("getData", args = c(
    thisProblem[names(thisProblem) %in% formalArgs("getData")],
    list(inputFile = inputFile)
  ))

  # Save some common variables for initializing arrays
  NComp = thisProblem$NComp
  CompName = thisProblem$CompName
  NSpec = thisProblem$NSpec
  SpecName = thisProblem$SpecName
  SpecStoich = thisProblem$SpecStoich
  SpecLogK = thisProblem$SpecLogK
  SpecK = array(10^SpecLogK, dim = NSpec, dimnames = list(SpecName))
  SpecCtoM = thisProblem$SpecCtoM
  CompType = thisProblem$CompType
  CompSiteDens = thisProblem$CompSiteDens

  BLName = thisProblem$BLName
  MetalName = thisProblem$MetalName
  BLMetalName = thisProblem$BLMetalName
  BLComps = which(SpecName %in% BLName)
  MetalComp = which(SpecName %in% MetalName)
  BLMetalSpecs = which(SpecName %in% BLMetalName)
  CATargetDefault = thisProblem$CATab$CA[iCA] * (10^-6) #thisProblem$DefCompFromNum[thisProblem$DefCompName == BLName]
  # 5.541E-8      =  0.05541 * 10^-6
  # (mol / kg) = (nmol / g) * (1 g-mol / 10^6 nmol-kg)

  TotConc = array(numeric(NComp), dimnames = list(CompName))
  TotMoles = array(numeric(NComp), dimnames = list(CompName))
  # CompConc = array(numeric(NComp), dimnames = list(CompName))
}

results.tab = as.data.frame(allInput$InLabObs)
for (iSpec in 1:NSpec){results.tab[,SpecName[iSpec]] = NA}
results.tab$Iter = NA
for (iComp in 1:NComp){results.tab[,paste0("T.",CompName[iComp])] = NA}
# results.tab$TotCu = NA
for (iObs in 1:allInput$NObs){

  if (QuietFlag != "Very Quiet"){print(paste0("Obs=",iObs))}

  TotConc = allInput$TotConcObs[iObs,]# mol/L
  TotMoles = TotConc * SpecCtoM[1:NComp] # mol/L * L = mol

  CATarget = CATargetDefault * TotConc[BLComps] / CompSiteDens[BLComps]
  # 9.86298E-13 = 5.541E-8 * 5.34E-10 / 3E-5
  # (mol / kg) = (mol / kg) * (mol / kg) * (kg / mol)
  #
  # 5.34E-10 = 1.78E-5 * 3E-5
  # (mol / kg) = (unitless) * (mol / kg)
  #
  # 9.86298E-13 = 5.541E-8 * 1.78E-5
  # (mol / kg) = (mol / kg) * (unitless)

  CompConc = initialGuess(NComp = NComp, CompName = CompName,
                          TotConc = TotConc, # mol / L or mol / kgw
                          CompType = CompType)

  LogCompConc = log10(CompConc)

  # Initialize Species Concentrations
  SpecConc = RCalcSpecConc(CompConc = CompConc,
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
  RR = RCalcResidual(
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

    (JacobianMatrix = RJacobian(NComp = NComp, NSpec = NSpec, SpecStoich = SpecStoich,
                                SpecConc = SpecConc, SpecCtoM = SpecCtoM, SpecName = SpecName,
                                MetalComp = MetalComp, BLMetalSpecs = BLMetalSpecs, DoTox = DoTox))

    (CompConcStep = RCalcStep(JacobianMatrix = JacobianMatrix, Resid = Resid,
                              NComp = NComp, CompType = CompType, CompName=CompName))

    (CompConc = RCompUpdate(CompConcStep = CompConcStep, CompConc = SpecConc[1:NComp],
                            CompCtoM = SpecCtoM[1:NComp], CompName = CompName))

    SpecConc = RCalcSpecConc(CompConc = CompConc,
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
                       # CompSiteDens = CompSiteDens,
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
      (CompConc_Full = RCompUpdate(CompConcStep = CompConcStep, CompConc = SpecConc[1:NComp],
                                   CompCtoM = SpecCtoM[1:NComp], CompName = CompName))

      # CompConc_Full[CompConc_Full > (TotConc / SpecCtoM[1:NComp])] =
      #   (TotConc/SpecCtoM[1:NComp])[CompConc_Full > (TotConc / SpecCtoM[1:NComp])]

      # SpecConc_Full = 10^CppCalcLogSpecConc(LogCompConc = log10(CompConc_Full),
      #                                       SpecLogK = SpecLogK,
      #                                       SpecStoich = SpecStoich,
      #                                       NComp = NComp,
      #                                       NSpec = NSpec)
      SpecConc_Full = RCalcSpecConc(CompConc = CompConc_Full,
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
        CompConc_Half = RCompUpdate(CompConcStep = CompConcStep * 0.5,
                                    CompConc = SpecConc[1:NComp],
                                    CompCtoM = SpecCtoM[1:NComp],
                                    CompName = CompName)

        # SpecConc_Half = 10^CppCalcLogSpecConc(LogCompConc = log10(CompConc_Half),
        #                                      SpecLogK = SpecLogK,
        #                                      SpecStoich = SpecStoich,
        #                                      NComp = NComp,
        #                                      NSpec = NSpec)
        SpecConc_Half = CppCalcSpecConc(CompConc = CompConc_Half,
                                        SpecK = SpecK,
                                        SpecStoich = SpecStoich,
                                        NComp = NComp,
                                        NSpec = NSpec)

        TotConc_Half = TotConc
        for (iComp in which(CompType == "FixedAct")){
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
        CompConc_Best = RCompUpdate(CompConcStep = CompConcStep * step_Best,
                                    CompName = CompName,
                                    CompConc = SpecConc[1:NComp], CompCtoM = SpecCtoM[1:NComp])

        # SpecConc_Best = 10^CppCalcLogSpecConc(LogCompConc = log10(CompConc_Best),
        #                                       SpecLogK = SpecLogK,
        #                                       SpecStoich = SpecStoich,
        #                                       NComp = NComp,
        #                                       NSpec = NSpec)
        SpecConc_Best = RCalcSpecConc(CompConc = CompConc_Best,
                                      SpecK = SpecK,
                                      SpecStoich = SpecStoich,
                                      NComp = NComp,
                                      NSpec = NSpec)

        TotConc_Best = TotConc
        for (iComp in which(CompType == "FixedAct")){
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

        # plot(x=NA, y= NA, xlab = "step", ylab = "MaxError", main = paste("Iter =",Iter),
        #      ylim = c(0,RR_Full$MaxError), xlim = c(0,1))
        # abline(h = MaxError_Last, col = "red")
        # segments(x0 = 1, x1 = 0.5, y0 = RR_Full$MaxError, y1 = RR_Half$MaxError, col = "red", lwd = 2)
        # segments(x0 = step_Best, x1 = 0.5, y0 = 0, y1 = RR_Half$MaxError, col = "red", lwd = 1, lty = 2)
        # segments(x0 = step_Best, x1 = 0.5, y0 = RR_Best$MaxError, y1 = RR_Half$MaxError, col = "red", lwd = 1)
        # for(stepi in seq(0,1,by=0.01)){
        #   CompConc_i = RCompUpdate(CompConcStep = CompConcStep * stepi, CompConc = SpecConc[1:NComp], CompCtoM = SpecCtoM[1:NComp])
        #
        #   SpecConc_i = CppCalcSpecConc(CompConc = CompConc_i,
        #                                SpecK = SpecK,
        #                                SpecStoich = SpecStoich,
        #                                NComp = NComp,
        #                                NSpec = NSpec)
        #   # SpecConc_i = 10^CppCalcLogSpecConc(LogCompConc = log10(CompConc_i),
        #   #                                      SpecLogK = SpecLogK,
        #   #                                      SpecStoich = SpecStoich,
        #   #                                      NComp = NComp,
        #   #                                      NSpec = NSpec)
        #
        #   TotConc_i = TotConc
        #   for (iComp in which(CompType == "FixedAct")){
        #     TotConc_i[iComp] = sum(SpecStoich[,iComp] * SpecConc_i)
        #   }
        #
        #   RR_i = RCalcResidual(NComp = NComp, NSpec = NSpec,
        #                          SpecConc = SpecConc_i,
        #                          SpecStoich = SpecStoich,
        #                          TotConc = TotConc_i,
        #                          SpecCtoM = SpecCtoM,
        #                          CompType = CompType)
        #   points(x = stepi, y = RR_i$MaxError)
        # }


        best_MaxError = which.min(c(RR_Full$MaxError, RR_Half$MaxError, RR_Best$MaxError))

      } else{
        best_MaxError = 1L
      }

      # plot(x=NA, y= NA, xlab = paste(CompOfIntName,"CompConc"),
      #      ylab = paste(CompOfIntName,"Resid"),
      #      log = "x",
      #      xlim = c(1E-12,1E-10),
      #      ylim = c(-1,1)*RR_Full$Resid[CompOfInt])
      # result.tab = data.frame(step3 = seq(0,1,by=0.1), CompConc_Cu = NA, Resid_Cu = NA)
      # for(i in 1:nrow(result.tab)){
      #   step3 = result.tab$step3[i]
      #   CompConc_Best = RCompUpdate(CompConcStep = CompConcStep * step3, CompConc = SpecConc[1:NComp],
      #                              CompType = CompType)
      #
      #   SpecConc_Best = 10^CppCalcLogSpecConc(LogCompConc = log10(CompConc_Best),
      #                                        SpecLogK = SpecLogK,
      #                                        SpecStoich = SpecStoich,
      #                                        NComp = NComp,
      #                                        NSpec = NSpec)
      #
      #   RR_Best = RCalcResidual(NComp = NComp, NSpec = NSpec,
      #                          SpecConc = SpecConc_Best,
      #                          SpecStoich = SpecStoich,
      #                          TotConc = TotConc,
      #                          SpecCtoM = SpecCtoM,
      #                          CompType = CompType)
      #   result.tab$CompConc_Cu[i] = CompConc_Best[CompOfInt]
      #   result.tab$Resid_Cu[i] = RR_Best$Resid[CompOfInt]
      #   # points(x = CompConc_Best[CompOfInt], y = RR_Best$Resid[CompOfInt])
      # }
      #
      # plot(x=result.tab$CompConc_Cu, y= result.tab$Resid_Cu,
      #      xlab = paste(CompOfIntName,"CompConc"),
      #      ylab = paste(CompOfIntName,"Resid"), type = "b")
      # text(x=result.tab$CompConc_Cu, y= result.tab$Resid_Cu,
      #      labels = result.tab$step3, pos = 2)


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

  SpecConc[all.BLMetal]

  results.tab[iObs, SpecName] = SpecConc
  results.tab$Iter[iObs] = Iter
  results.tab[iObs, paste0("T.",CompName)] = TotConc

  # results.tab$TotCu[iObs] = CalcTotConc[MetalComp]

}

end.time = Sys.time()
end.time - start.time


SpecConc_BLMetal_old = c(
  Cu = 7.43926E-11,
  BL1 = 5.29446E-10,
  `BL1-Cu` = 9.18887E-13,
  `BL1-CuOH` = 6.83595E-14
)
all.BLMetalName = c(MetalName, BLName, BLMetalName)
all.BLMetal = c(MetalComp, BLComps, BLMetalSpecs)

(results.tab[,all.BLMetalName] - SpecConc_BLMetal_old) / SpecConc_BLMetal_old

SpecConc_BLMetal_old
results.tab[,all.BLMetalName]

results.tab[,c("Cu","BL1","BL1-Cu","BL1-CuOH","T.Cu","T.BL1")]
results.tab[,c("BL1","BL1-Cu","BL1-CuOH","BL1-Ca","BL1-Mg","BL1-Na","BL1-H","T.BL1")]
sum(results.tab[,c("BL1-Cu","BL1-CuOH")]) * 10^6 / results.tab$T.BL1 * CompSiteDens[BLComps]
sum(results.tab[,c("BL1","BL1-Cu","BL1-CuOH","BL1-Ca","BL1-Mg","BL1-Na","BL1-H")])
