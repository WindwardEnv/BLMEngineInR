devtools::clean_dll()
devtools::load_all()

# test stuff
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

  NBLMetal = thisProblem$NBLMetal
  BLName = thisProblem$BLName
  MetalName = thisProblem$MetalName
  BLMetalName = thisProblem$BLMetalName
  BLComp = which(SpecName %in% BLName)
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

  CATarget = CATargetDefault * TotConc[BLComp] / CompSiteDens[BLComp]
  # 9.86298E-13 = 5.541E-8 * 5.34E-10 / 3E-5
  # (mol / kg) = (mol / kg) * (mol / kg) * (kg / mol)
  #
  # 5.34E-10 = 1.78E-5 * 3E-5
  # (mol / kg) = (unitless) * (mol / kg)
  #
  # 9.86298E-13 = 5.541E-8 * 1.78E-5
  # (mol / kg) = (mol / kg) * (unitless)

  # START CHESS ------------------------

  out = RCHESS(DoPartialSteps, QuietFlag, ConvergenceCriteria, MaxIter,
               NComp, NSpec, NBLMetal,
               SpecConc, SpecLogK, SpecStoich, SpecCtoM, SpecName,
               CompType, CompName, TotMoles, TotConc,
               DoTox, MetalName, MetalComp, BLMetalSpecs, CATarget)


  SpecConc = out$SpecConc

  # END CHESS ------------------------

  results.tab[iObs, SpecName] = SpecConc
  results.tab$Iter[iObs] = Iter
  results.tab[iObs, paste0("T.",CompName)] = TotConc

  # results.tab$TotCu[iObs] = CalcTotConc[MetalComp]

}

end.time = Sys.time()
end.time - start.time
