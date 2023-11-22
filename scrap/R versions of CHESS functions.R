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



RJacobian = function(NComp, NSpec, SpecStoich, SpecConc, SpecCtoM, SpecName,
                     iMetal, iBLMetal, iTox) {
  # inputs:
  #   NComp
  #   NSpec
  #   SpecStoich
  #   SpecConc
  #   SpecCtoM
  # outputs:
  JacobianMatrix2 = JacobianMatrix = matrix(data = 0, nrow = NComp, ncol = NComp,
             dimnames = list(SpecName[1:NComp], SpecName[1:NComp]))
  # variables:
  #   Sum
  #   iComp1, iComp2, iSpec


  SpecMoles = SpecConc * SpecCtoM
  SpecMolesByComp = matrix(SpecMoles, nrow = NSpec, ncol = NComp,
                           byrow = F,
                           dimnames = list(SpecName, SpecName[1:NComp])) * SpecStoich



  CompStoichDivMoles = SpecStoich *
    matrix(1/SpecMoles[1:NComp], nrow = NSpec, ncol = NComp, byrow = T,
           dimnames = list(SpecName, SpecName[1:NComp]))


  JacobianMatrix2 = t(SpecMolesByComp) %*% CompStoichDivMoles
  JacobianMatrix2 - JacobianMatrix
  JacobianMatrix2 / JacobianMatrix

  # for (iComp1 in 1:NComp) {
  #   for (iComp2 in 1:NComp) {
  #     Sum = sum(SpecStoich[,iComp2] * SpecStoich[,iComp1] * SpecMoles)
  #     if (SpecConc[iComp1] != 0) {
  #       JacobianMatrix2[iComp1, iComp2] = Sum / (SpecMoles[iComp2])
  #     }
  #   }
  #
  #   iComp1SpecMoles = matrix((SpecStoich[,iComp1] * SpecMoles), nrow = 1, ncol = NSpec)
  #   t(iComp1SpecMoles %*% SpecStoich) %*%
  #     matrix(1/SpecMoles[1:NComp], nrow = 1, ncol = NComp, dimnames = list(NULL, CompName))
  #
  #
  # }
  #
  # SpecStoich %*% SpecMoles
  #

  for (iComp1 in 1:NComp) {
    for (iComp2 in 1:NComp) {
      Sum = 0
      if (iTox && (iComp1 == iMetal)){
        for (iSpec in iBLMetal){
          Sum = Sum + (SpecStoich[iSpec, iComp2] * SpecStoich[iSpec, iComp1] *
                         SpecConc[iSpec] * SpecCtoM[iSpec])
        }
      } else {
        for (iSpec in 1:NSpec) {
          Sum = Sum + (SpecStoich[iSpec, iComp2] * SpecStoich[iSpec, iComp1] *
                         SpecConc[iSpec] * SpecCtoM[iSpec])
        }
      }
      if (SpecConc[iComp1] != 0) {
        JacobianMatrix[iComp1, iComp2] = Sum / (SpecConc[iComp2] * SpecCtoM[iComp2])
      }
    }
  }


  return(JacobianMatrix)
}

RCalcStep = function(JacobianMatrix, Resid, NComp, CompType, SpecName){

  i.solve = which(CompType != "FixedAct")
  n.solve = length(i.solve)

  (Resid.solve = matrix(Resid[i.solve], nrow = n.solve, ncol = 1,
                        dimnames = list(SpecName[i.solve],NULL)))
  (JacobianMatrix.solve = JacobianMatrix[i.solve, i.solve, drop = F])

  # find the matrix inverse of JacobianMatrix by SVD
  JacobianMatrixsvd = svd(JacobianMatrix.solve)
  dinv = diag(1/JacobianMatrixsvd$d, nrow = n.solve, ncol = n.solve)
  U = JacobianMatrixsvd$u
  V = JacobianMatrixsvd$v
  UT = t(U)
  (JacobianMatrixinv = V %*% dinv %*% UT)

  # JacobianMatrixinv = solve(JacobianMatrix.solve)
  (CompConcStep.solve = as.numeric(JacobianMatrixinv %*% Resid.solve))

  CompConcStep = array(0, dim = NComp, dimnames = list(SpecName[1:NComp]))
  CompConcStep[i.solve] = CompConcStep.solve

  # CompConcStep[CompType == "FixedAct"] = 0

  return(CompConcStep)

}

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

RCalcResidual = function(NComp, NSpec, SpecConc, SpecStoich, TotMoles, SpecCtoM,
                         CompName, CompType, #CompSiteDens,
                         iMetal, iBLMetal, CATarget, iTox){
  # inputs:
  #   NComp - number of components
  #   NSpec - number of species
  #   SpecConc - species free concentrations (mol/L or mol/gww)
  #   SpecStoich - species stoichiometry
  #   TotConc - component total concentrations (mol)
  #   SpecCtoM - mass compartment concentration to mass conversion for each species (L or gww)
  #   CompName - component names
  #   CompType - component types
  #   iMetal - integer, position of the metal component in the component list
  #   iBLMetal - integer vector, positon of the metal-bound biotic ligand in the species list
  #   CATarget - the target critical accumulation for toxicity runs
  #   iTox - boolean, TRUE means this is a toxicity run, FALSE means speciation run
  # outputs
  #   Resid - vector(NComp) (maybe...useful for debugging)
  Resid = array(dim=NComp, dimnames = list(CompName))
  #   MaxError - double - maximum of absolute ratios of residuals to totals
  #   WhichMax - which component has the highest absolute error
  #   CalcTotConc - the calculated total concentrations (mol)
  # variables:
  #   double CalcTotConc
  #   double ThisError

  CalcTotMoles = array((SpecConc * SpecCtoM) %*% SpecStoich, dim = NComp, dimnames = list(CompName))
  # CalcTotConc = as.numeric(matrix(SpecConc, nrow = 1, ncol = NSpec) %*% SpecStoich)
  Resid = CalcTotMoles - TotMoles # * SpecCtoM[1:NComp]
  Resid[CompType == "FixedAct"] = 0.0
  ThisError = abs(Resid / TotMoles)
  if(iTox){
    Resid[iMetal] = sum(SpecConc[iBLMetal]) - CATarget
    # sum( (mol/L or mol/gww) * (L or gww) ) - mol = sum(mol) - mol = mol
    ThisError[iMetal] = abs(Resid[iMetal] / CATarget)
  }
  MaxError = max(ThisError)
  WhichMax = which.max(ThisError)

  CalcTotConc = CalcTotMoles / SpecCtoM[1:NComp]

  return(list(
    MaxError = MaxError,
    WhichMax = WhichMax,
    Resid = Resid,
    CompError = ThisError,
    CalcTotConc = CalcTotConc,
    CalcTotMoles = CalcTotMoles
  ))
}

RCHESSIter <- function(DoPartialSteps, QuietFlag, NComp, NSpec,
                       CompConc, CompType, CompName,
                       SpecLogK, SpecStoich, SpecName,
                       SpecCtoM, TotMoles, TotConc,
                       iTox, MetalName, iMetal, iBLMetal, CATarget) {

  return()
}


SpecConc_BLMetal_old = c(
  Cu = 7.43926E-11,
  BL1 = 5.29446E-10,
  `BL1-Cu` = 9.18887E-13,
  `BL1-CuOH` = 6.83595E-14
)
all.BLMetalName = c(MetalName, BLName, BLMetalName)
all.BLMetal = c(iMetal, iBL, iBLMetal)

start.time = Sys.time()
# Get the problem set up
paramFile = "scrap/parameter file format/full_organic.dat4"
inputFile = "scrap/parameter file format/full_organic.blm4"
QuietFlag = c("Very Quiet","Quiet","Debug")[3]
iCA = 1
iTox = T
DoPartialSteps = F
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
  iBL = which(SpecName %in% BLName)
  iMetal = which(SpecName %in% MetalName)
  iBLMetal = which(SpecName %in% BLMetalName)
  CATargetDefault = thisProblem$CATab$CA[iCA] * (10^-6) #thisProblem$DefCompFromNum[thisProblem$DefCompName == BLName]
  # (mol bound / mol total) = (umol bound / umol total) * (1 mol / 10^6 umol)

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

  CATarget = CATargetDefault * TotConc[iBL] / CompSiteDens[iBL]
  # (mol bound / kg) = (mol bound / mol total) * (mol available to bind / kg) * (mol total / mol available to bind)

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
  for (iComp in which((CompType == "FixedAct") | (iTox & (CompName == MetalName)))){
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
    iMetal = iMetal,
    iBLMetal = iBLMetal,
    CATarget = CATarget,
    iTox = iTox
  )
  Resid = RR$Resid
  MaxError = RR$MaxError
  WhichMax = RR$WhichMax
  CalcTotConc = RR$CalcTotConc
  CalcTotMoles = RR$CalcTotMoles

  # Begin iterating
  Iter = 0
  while ((MaxError > 0.00001) & (Iter <= 30)){

    Iter = Iter + 1
    (MaxError_Last = MaxError)

    (JacobianMatrix = RJacobian(NComp = NComp, NSpec = NSpec, SpecStoich = SpecStoich,
                                SpecConc = SpecConc, SpecCtoM = SpecCtoM, SpecName = SpecName,
                                iMetal = iMetal, iBLMetal = iBLMetal, iTox = iTox))

    (CompConcStep = RCalcStep(JacobianMatrix = JacobianMatrix, Resid = Resid,
                              NComp = NComp, CompType = CompType, SpecName=SpecName))

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
                       iMetal = iMetal,
                       iBLMetal = iBLMetal,
                       CATarget = CATarget,
                       iTox = iTox)

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
                              iMetal = iMetal,
                              iBLMetal = iBLMetal,
                              CATarget = CATarget,
                              iTox = iTox)
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
                                iMetal = iMetal,
                                iBLMetal = iBLMetal,
                                CATarget = CATarget,
                                iTox = iTox)

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
                                iMetal = iMetal,
                                iBLMetal = iBLMetal,
                                CATarget = CATarget,
                                iTox = iTox)

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

  # results.tab$TotCu[iObs] = CalcTotConc[iMetal]

}



end.time = Sys.time()
end.time - start.time

(results.tab[,all.BLMetalName] - SpecConc_BLMetal_old) / SpecConc_BLMetal_old

SpecConc_BLMetal_old
results.tab[,all.BLMetalName]

results.tab[,c("Cu","BL1","BL1-Cu","BL1-CuOH","T.Cu","T.BL1")]
results.tab[,c("BL1","BL1-Cu","BL1-CuOH","BL1-Ca","BL1-Mg","BL1-Na","BL1-H","T.BL1")]
sum(results.tab[,c("BL1-Cu","BL1-CuOH")]) * 10^6 / results.tab$T.BL1 * CompSiteDens[iBL]
sum(results.tab[,c("BL1","BL1-Cu","BL1-CuOH","BL1-Ca","BL1-Mg","BL1-Na","BL1-H")])
