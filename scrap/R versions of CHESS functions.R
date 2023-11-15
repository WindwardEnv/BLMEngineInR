RCalcSpecConc = function(CompConc, SpecK, SpecStoich, NComp = length(CompConc),
                         NSpec = length(SpecK), SpecName = names(SpecK)){

  SpecConc = rep(1, NSpec)
  for (iSpec in 1:NSpec){
    X = CompConc ^ SpecStoich[iSpec,]
    SpecConc[iSpec] = prod(X) * SpecK[iSpec]
  }

  names(SpecConc) = SpecName
  return(SpecConc)
}



RJacobian = function(NComp, NSpec, SpecStoich, SpecConc, SpecCtoM, SpecName) {
  # inputs:
  #   NComp
  #   NSpec
  #   SpecStoich
  #   SpecConc
  #   SpecCtoM
  # outputs:
  Z = matrix(data = 0, nrow = NComp, ncol = NComp)
  # variables:
  #   Sum
  #   iComp1, iComp2, iSpec

  # for (iComp1 in 1:NComp){
  #   rowSums((SpecStoich[,iComp1]* SpecConc) %*% SpecStoich)
  #   for (iComp2 in 1:NComp){
  #     if (SpecConc[iComp2] != 0) {
  #       Sum = sum((SpecStoich[,iComp1]* SpecConc) * SpecStoich[,iComp2] )
  #       Z[iComp1, iComp2] = Sum / SpecConc[iComp2]
  #     }
  #   }
  # }

  for (iComp1 in 1:NComp) {
    for (iComp2 in 1:NComp) {
      Sum = 0
      for (iSpec in 1:NSpec) {
        Sum = Sum + (SpecStoich[iSpec, iComp2] * SpecStoich[iSpec, iComp1] *
          SpecConc[iSpec] * SpecCtoM[iSpec])
      }
      if (SpecConc[iComp1] != 0) {
        # Z[iComp1, iComp2] = Sum / (SpecConc[iComp1])
        Z[iComp1, iComp2] = Sum / (SpecConc[iComp2])
      }
    }
  }

  rownames(Z) = SpecName[1:NComp]
  colnames(Z) = SpecName[1:NComp]

  return(Z)
}

RCalcStep = function(Z, Resid, NComp, CompType, SpecName){

  i.solve = which(CompType != "FixedAct")
  n.solve = length(i.solve)

  (Resid.solve = matrix(Resid[i.solve], nrow = n.solve, ncol = 1))
  (Z.solve = Z[i.solve, i.solve, drop = F])

  # find the matrix inverse of Z by SVD
  Zsvd = svd(Z.solve)
  dinv = diag(1/Zsvd$d, nrow = n.solve, ncol = n.solve)
  U = Zsvd$u
  V = Zsvd$v
  UT = t(U)
  (Zinv = V %*% dinv %*% UT)

  # Zinv = solve(Z.solve)
  (X.solve = as.numeric(Zinv %*% Resid.solve))

  X = array(0, dim = NComp)
  X[i.solve] = X.solve

  # X[CompType == "FixedAct"] = 0

  names(X) = SpecName[1:NComp]

  return(X)

}

RCompUpdate = function(X, CompConc, CompCtoM, CompName){


  (oldCompConc = CompConc)
  # CompConc[X < oldCompConc] = (oldCompConc - X)[X < oldCompConc]
  # CompConc = (oldCompConc - X / CompCtoM)
  CompConc = (oldCompConc - X)
  # CompConc[X >= oldCompConc] = oldCompConc[X >= oldCompConc] / 10
  CompConc[CompConc <= 0] = oldCompConc[CompConc <= 0] / 10


  # for (iComp in 1:NComp){
  #   if (X[iComp] >= oldCompConc[iComp]) {
  #     CompConc[iComp] = oldCompConc[iComp] / 10
  #   } else {
  #     CompConc[iComp] = CompConc[iComp] - X[iComp]
  #   }
  #   if (CompConc[iComp] <= 0.0) {
  #     CompConc[iComp] = 1E-10
  #   }
  # }#NEXT cc

  names(CompConc) = CompName
  return(CompConc)

}

RCalcResidual = function(NComp, NSpec, SpecConc, SpecStoich, TotConc, SpecCtoM, CompType){
  # inputs:
  #   SpecConc - species free concentrations
  #   SpecStoich - species stoichiometry
  #   TotConc - component total concentrations
  #   SpecCtoM
  # outputs
  #   Resid - vector(NComp) (maybe...useful for debugging)
  Resid = array(dim=NComp)
  #   MaxError - double - maximum of absolute ratios of residuals to totals
  # variables:
  #   double CalcTotConc
  #   double ThisError

  CalcTotConc = as.numeric((SpecConc * SpecCtoM) %*% SpecStoich)
  # CalcTotConc = as.numeric(matrix(SpecConc, nrow = 1, ncol = NSpec) %*% SpecStoich)
  Resid = CalcTotConc - TotConc
  Resid[CompType == "FixedAct"] = 0.0
  ThisError = abs(Resid / TotConc)
  MaxError = max(ThisError)
  WhichMax = which.max(ThisError)

  # MaxError = 0.0
  # for (iComp in 1:NComp) {
  #   CalcTotConc = 0
  #   for (iSpec in 1:NSpec){
  #     CalcTotConc = CalcTotConc + (SpecConc[iSpec] * SpecStoich[iSpec, iComp] * SpecCtoM[iSpec])
  #   }
  #   Resid[iComp] = TotConc[iComp] - CalcTotConc
  #   ThisError = abs(Resid[iComp] / TotConc[iComp])
  #   if (ThisError > MaxError){
  #     MaxError = ThisError
  #   }
  # }

  return(list(MaxError = MaxError, Resid = Resid,
              WhichMax = WhichMax, CalcTotConc = CalcTotConc))
}

# Get the problem set up
paramFile = "scrap/parameter file format/full_organic.dat"
inputFile = "scrap/parameter file format/full_organic.blm4"
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
  SpecCtoM = thisProblem$SpecCtoM
  CompType = thisProblem$CompType

  TotConc = array(numeric(NComp), dimnames = list(CompName))
  # CompConc = array(numeric(NComp), dimnames = list(CompName))

  TotConc = allInput$TotConcObs[1,]
  TotConc = TotConc * thisProblem$SpecCtoM[1:NComp]

  # SpecConc[1:NComp] = do.call(initialGuess, args = thisInput[formalArgs(initialGuess)])
  # SpecConc[c(thisProblem$InCompName,"H","BL1")] =
  #   Full_OrganicDataFreeConc[c(thisProblem$InCompName,"H","BL1")]
  # SpecConc[thisProblem$DefCompName[grepl("DOC",thisProblem$DefCompName)]] =
  #   Full_OrganicDataFreeConc["DOC"] *
  #   thisProblem$DefCompSiteDens[grepl("DOC",thisProblem$DefCompName)]/1000
  CompConc = initialGuess(NComp = NComp, CompName = CompName,
                          TotConc = TotConc / thisProblem$SpecCtoM[1:NComp],
                          CompType = CompType)
  LogCompConc = log10(CompConc)
  SpecConc = RCalcSpecConc(CompConc = CompConc,
                           SpecK = 10^SpecLogK,
                           SpecStoich = SpecStoich,
                           NComp = NComp,
                           NSpec = NSpec)
  SpecConc[1:NComp] = CompConc
  for (iComp in which(CompType == "FixedAct")){
    TotConc[iComp] = sum(SpecStoich[,iComp] * (SpecConc * SpecCtoM))
  }

  RR = RCalcResidual(NComp = NComp, NSpec = NSpec,
                      SpecConc = SpecConc,
                      SpecStoich = SpecStoich,
                      TotConc = TotConc,
                      SpecCtoM = SpecCtoM,
                      CompType = CompType)
  Resid = RR$Resid
  MaxError = RR$MaxError

}


CompOfInt = which(CompName == "Cu")
CompOfIntName = CompName[CompOfInt]
# MaxError = 10^99
Iter = 0
while ((MaxError > 0.0000001) & (Iter <= 30)){
  Iter = Iter + 1
  (MaxError_Last = MaxError)

  (Z = RJacobian(NComp = NComp, NSpec = NSpec, SpecStoich = SpecStoich,
                SpecConc = SpecConc, SpecCtoM = SpecCtoM, SpecName = SpecName))
  (X = RCalcStep(Z = Z, Resid = Resid, NComp = NComp, CompType = CompType, SpecName=SpecName))

  # Full Step
  (CompConc_Full = RCompUpdate(X = X, CompConc = SpecConc[1:NComp],
                               CompCtoM = SpecCtoM[1:NComp], CompName = CompName))

  # CompConc_Full[CompConc_Full > (TotConc / SpecCtoM[1:NComp])] =
  #   (TotConc/SpecCtoM[1:NComp])[CompConc_Full > (TotConc / SpecCtoM[1:NComp])]

  # SpecConc_Full = 10^CppCalcLogSpecConc(LogCompConc = log10(CompConc_Full),
  #                                       SpecLogK = SpecLogK,
  #                                       SpecStoich = SpecStoich,
  #                                       NComp = NComp,
  #                                       NSpec = NSpec)
  SpecConc_Full = RCalcSpecConc(CompConc = CompConc_Full,
                                SpecK = 10^SpecLogK,
                                SpecStoich = SpecStoich,
                                NComp = NComp,
                                NSpec = NSpec,
                                SpecName = SpecName)
  TotConc_Full = TotConc
  for (iComp in which(CompType == "FixedAct")){
    TotConc_Full[iComp] = sum(SpecStoich[,iComp] * SpecConc_Full)
  }

  RR_Full = RCalcResidual(NComp = NComp, NSpec = NSpec,
                          SpecConc = SpecConc_Full,
                          SpecStoich = SpecStoich,
                          TotConc = TotConc_Full,
                          SpecCtoM = SpecCtoM,
                          CompType = CompType)
  best_MaxError = 1L
  # if (RR_Full$MaxError > MaxError_Last){
  #   # Half Step
  #   CompConc_Half = RCompUpdate(X = X * 0.5, CompConc = SpecConc[1:NComp], CompCtoM = SpecCtoM[1:NComp])
  #
  #   # SpecConc_Half = 10^CppCalcLogSpecConc(LogCompConc = log10(CompConc_Half),
  #   #                                      SpecLogK = SpecLogK,
  #   #                                      SpecStoich = SpecStoich,
  #   #                                      NComp = NComp,
  #   #                                      NSpec = NSpec)
  #   SpecConc_Half = CppCalcSpecConc(CompConc = CompConc_Half,
  #                                   SpecK = 10^SpecLogK,
  #                                   SpecStoich = SpecStoich,
  #                                   NComp = NComp,
  #                                   NSpec = NSpec)
  #
  #   TotConc_Half = TotConc
  #   for (iComp in which(CompType == "FixedAct")){
  #     TotConc_Half[iComp] = sum(SpecStoich[,iComp] * SpecConc_Half)
  #   }
  #
  #   RR_Half = RCalcResidual(NComp = NComp, NSpec = NSpec,
  #                          SpecConc = SpecConc_Half,
  #                          SpecStoich = SpecStoich,
  #                          TotConc = TotConc_Half,
  #                          SpecCtoM = SpecCtoM,
  #                          CompType = CompType)
  #
  #   # Best Step
  #   step_Best = 1 - RR_Full$MaxError * (1-0.5) / (RR_Full$MaxError - RR_Half$MaxError)
  #   CompConc_Best = RCompUpdate(X = X * step_Best, CompConc = SpecConc[1:NComp], CompCtoM = SpecCtoM[1:NComp])
  #
  #   # SpecConc_Best = 10^CppCalcLogSpecConc(LogCompConc = log10(CompConc_Best),
  #   #                                       SpecLogK = SpecLogK,
  #   #                                       SpecStoich = SpecStoich,
  #   #                                       NComp = NComp,
  #   #                                       NSpec = NSpec)
  #   SpecConc_Best = RCalcSpecConc(CompConc = CompConc_Best,
  #                                   SpecK = 10^SpecLogK,
  #                                   SpecStoich = SpecStoich,
  #                                   NComp = NComp,
  #                                   NSpec = NSpec)
  #
  #   TotConc_Best = TotConc
  #   for (iComp in which(CompType == "FixedAct")){
  #     TotConc_Best[iComp] = sum(SpecStoich[,iComp] * SpecConc_Best)
  #   }
  #
  #   RR_Best = RCalcResidual(NComp = NComp, NSpec = NSpec,
  #                           SpecConc = SpecConc_Best,
  #                           SpecStoich = SpecStoich,
  #                           TotConc = TotConc_Best,
  #                           SpecCtoM = SpecCtoM,
  #                           CompType = CompType)
  #
  #   # plot(x=NA, y= NA, xlab = "step", ylab = "MaxError", main = paste("Iter =",Iter),
  #   #      ylim = c(0,RR_Full$MaxError), xlim = c(0,1))
  #   # abline(h = MaxError_Last, col = "red")
  #   # segments(x0 = 1, x1 = 0.5, y0 = RR_Full$MaxError, y1 = RR_Half$MaxError, col = "red", lwd = 2)
  #   # segments(x0 = step_Best, x1 = 0.5, y0 = 0, y1 = RR_Half$MaxError, col = "red", lwd = 1, lty = 2)
  #   # segments(x0 = step_Best, x1 = 0.5, y0 = RR_Best$MaxError, y1 = RR_Half$MaxError, col = "red", lwd = 1)
  #   # for(stepi in seq(0,1,by=0.01)){
  #   #   CompConc_i = RCompUpdate(X = X * stepi, CompConc = SpecConc[1:NComp], CompCtoM = SpecCtoM[1:NComp])
  #   #
  #   #   SpecConc_i = CppCalcSpecConc(CompConc = CompConc_i,
  #   #                                SpecK = 10^SpecLogK,
  #   #                                SpecStoich = SpecStoich,
  #   #                                NComp = NComp,
  #   #                                NSpec = NSpec)
  #   #   # SpecConc_i = 10^CppCalcLogSpecConc(LogCompConc = log10(CompConc_i),
  #   #   #                                      SpecLogK = SpecLogK,
  #   #   #                                      SpecStoich = SpecStoich,
  #   #   #                                      NComp = NComp,
  #   #   #                                      NSpec = NSpec)
  #   #
  #   #   TotConc_i = TotConc
  #   #   for (iComp in which(CompType == "FixedAct")){
  #   #     TotConc_i[iComp] = sum(SpecStoich[,iComp] * SpecConc_i)
  #   #   }
  #   #
  #   #   RR_i = RCalcResidual(NComp = NComp, NSpec = NSpec,
  #   #                          SpecConc = SpecConc_i,
  #   #                          SpecStoich = SpecStoich,
  #   #                          TotConc = TotConc_i,
  #   #                          SpecCtoM = SpecCtoM,
  #   #                          CompType = CompType)
  #   #   points(x = stepi, y = RR_i$MaxError)
  #   # }
  #
  #
  #   best_MaxError = which.min(c(RR_Full$MaxError, RR_Half$MaxError, RR_Best$MaxError))
  #
  # } else{
  #   best_MaxError = 1L
  # }

  # plot(x=NA, y= NA, xlab = paste(CompOfIntName,"CompConc"),
  #      ylab = paste(CompOfIntName,"Resid"),
  #      log = "x",
  #      xlim = c(1E-12,1E-10),
  #      ylim = c(-1,1)*RR_Full$Resid[CompOfInt])
  # result.tab = data.frame(step3 = seq(0,1,by=0.1), CompConc_Cu = NA, Resid_Cu = NA)
  # for(i in 1:nrow(result.tab)){
  #   step3 = result.tab$step3[i]
  #   CompConc_Best = RCompUpdate(X = X * step3, CompConc = SpecConc[1:NComp],
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
    # print("Full Step")
  } else if (best_MaxError == 2L){
    CompConc = CompConc_Half
    SpecConc = SpecConc_Half
    Resid = RR_Half$Resid
    WhichMax = RR_Half$WhichMax
    MaxError = RR_Half$MaxError
    CalcTotConc = RR_Half$CalcTotConc
    # print("Half Step")
  } else if (best_MaxError == 3L){
    CompConc = CompConc_Best
    SpecConc = SpecConc_Best
    Resid = RR_Best$Resid
    WhichMax = RR_Best$WhichMax
    MaxError = RR_Best$MaxError
    CalcTotConc = RR_Best$CalcTotConc
    print("Best Step")
  }
  LogCompConc = log10(CompConc)


  print(paste0("Iter=",Iter, ", WhichMax=",CompName[WhichMax],
               ", MaxError=",format(x = MaxError, scientific = T)))
  # print(paste0("Iter=",Iter,
  #              ", WhichMax=",CompName[WhichMax],
  #              ", Resid[WhichMax]= ",signif(Resid[WhichMax], 6),
  #              ", SpecConc[WhichMax]= ",signif(SpecConc[WhichMax],6)))
  # print(paste0("Iter=",Iter, ", Resid[",CompOfIntName,"]= ",Resid[CompOfInt],
  #              ", CalcTotConc[",CompOfIntName,"]=",CalcTotConc[CompOfInt],
  #              ", CompConc[",CompOfIntName,"]=",CompConc[CompOfInt]))


}










abline(h = MaxError_Last, col = "red")
segments(x0 = 1, x1 = 0.5, y0 = RR_Full$MaxError, y1 = RR_Half$MaxError, col = "red", lwd = 2)
segments(x0 = step_Best, x1 = 0.5, y0 = 0, y1 = RR_Half$MaxError, col = "red", lwd = 1, lty = 2)
segments(x0 = step_Best, x1 = 0.5, y0 = RR_Best$MaxError, y1 = RR_Half$MaxError, col = "red", lwd = 1)


names(SpecConc) = SpecName
iBL =  which(SpecName == "BL1")
SpecConc_i = SpecConc
results.tab = data.frame(CompConc = seq(1E-10, 1.78E-5, length.out = 1000), Resid = NA)
for (i in 1:nrow(results.tab)){
  SpecConc_i[iBL] = results.tab$CompConc[i]
  SpecConc_i = CppCalcSpecConc(CompConc = SpecConc_i[1:NComp],
                               SpecK = 10^SpecLogK,
                               SpecStoich = SpecStoich,
                               NComp = NComp,
                               NSpec = NSpec)
  TotConc_i = TotConc
  for (iComp in which(CompType == "FixedAct")){
    TotConc_i[iComp] = sum(SpecStoich[,iComp] * SpecConc_i)
  }

  RR_i = RCalcResidual(NComp = NComp, NSpec = NSpec,
                       SpecConc = SpecConc_i,
                       SpecStoich = SpecStoich,
                       TotConc = TotConc_i,
                       SpecCtoM = SpecCtoM,
                       CompType = CompType)
  results.tab$Resid[i] = RR_i$Resid[iBL]
}
plot(x=results.tab$CompConc, y= results.tab$Resid,
     xlab = "CompConc[BL1]", ylab = "Resid[BL1]")
abline(v = TotConc[iBL], col = "red")

i = 999
SpecConc_i[iBL] = results.tab$CompConc[i]
SpecConc_i = CppCalcSpecConc(CompConc = SpecConc_i[1:NComp],
                             SpecK = 10^SpecLogK,
                             SpecStoich = SpecStoich,
                             NComp = NComp,
                             NSpec = NSpec)
TotConc_i = TotConc
for (iComp in which(CompType == "FixedAct")){
  TotConc_i[iComp] = sum(SpecStoich[,iComp] * SpecConc_i)
}

RR_i = RCalcResidual(NComp = NComp, NSpec = NSpec,
                     SpecConc = SpecConc_i,
                     SpecStoich = SpecStoich,
                     TotConc = TotConc_i,
                     SpecCtoM = SpecCtoM,
                     CompType = CompType)
results.tab$Resid[i] = RR_i$Resid[iBL]
points(x = SpecConc_i[iBL], y = RR_i$Resid[iBL], pch = 20, col = "green")

(Z_i = RJacobian(NComp = NComp, NSpec = NSpec, SpecStoich = SpecStoich,
                 SpecConc = SpecConc_i, SpecCtoM = SpecCtoM))
(X_i = RCalcStep(Z = Z_i, Resid = RR_i$Resid, NComp = NComp, CompType = CompType))

# Full Step
(SpecConc_i[1:NComp] = RCompUpdate(X = X_i, CompConc = SpecConc_i[1:NComp],
                          CompCtoM = SpecCtoM[1:NComp]))
