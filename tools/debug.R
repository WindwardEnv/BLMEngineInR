# devtools::clean_dll()
devtools::load_all()

# test stuff
start.time = Sys.time()
# Get the problem set up
ParamFile = "scrap/parameter file format/full_organic.dat4"
InputFile = "scrap/parameter file format/full_organic.blm4"
QuietFlag = c("Very Quiet","Quiet","Debug")[3]
iCA = 1
DoTox = T
DoPartialSteps = F
ConvergenceCriteria = 0.001
MaxIter = 30L
{

  ThisProblem = DefineProblem(ParamFile = ParamFile)
  AllInput = do.call("GetData", args = c(
    ThisProblem[names(ThisProblem) %in% formalArgs("GetData")],
    list(InputFile = InputFile)
  ))

  # Save some common variables for initializing arrays
  NComp = ThisProblem$NComp
  CompName = ThisProblem$CompName
  NSpec = ThisProblem$NSpec
  SpecName = ThisProblem$SpecName
  SpecStoich = ThisProblem$SpecStoich
  SpecLogK = ThisProblem$SpecLogK
  SpecK = array(10^SpecLogK, dim = NSpec, dimnames = list(SpecName))
  SpecCtoM = ThisProblem$SpecCtoM
  CompType = ThisProblem$CompType
  CompSiteDens = ThisProblem$CompSiteDens

  NBLMetal = ThisProblem$NBLMetal
  BLName = ThisProblem$BLName
  MetalName = ThisProblem$MetalName
  BLMetalName = ThisProblem$BLMetalName
  BLComp = which(SpecName %in% BLName)
  MetalComp = which(SpecName %in% MetalName)
  BLMetalSpecs = which(SpecName %in% BLMetalName)
  CATargetDefault = ThisProblem$CATab$CA[iCA] * (10^-6) #ThisProblem$DefCompFromNum[ThisProblem$DefCompName == BLName]
  # 5.541E-8      =  0.05541 * 10^-6
  # (mol / kg) = (nmol / g) * (1 g-mol / 10^6 nmol-kg)

  TotConc = array(numeric(NComp), dimnames = list(CompName))
  TotMoles = array(numeric(NComp), dimnames = list(CompName))
  # CompConc = array(numeric(NComp), dimnames = list(CompName))
}

ResultsTable = as.data.frame(AllInput$InLabObs)
for (iSpec in 1:NSpec){ResultsTable[,SpecName[iSpec]] = NA}
ResultsTable$Iter = NA
for (iComp in 1:NComp){ResultsTable[,paste0("T.",CompName[iComp])] = NA}
# ResultsTable$TotCu = NA
# for (iObs in 1:AllInput$NObs){
iObs = 1; {

  if (QuietFlag != "Very Quiet"){print(paste0("Obs=",iObs))}

  TotConc = AllInput$TotConcObs[iObs,]# mol/L
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

  ChessResults = CHESS(DoPartialSteps, QuietFlag, ConvergenceCriteria, MaxIter,
                       NComp, NSpec, NBLMetal,
                       SpecK, SpecStoich, SpecCtoM, SpecName,
                       CompType, CompName, TotMoles, TotConc,
                       DoTox, MetalName, MetalComp, BLMetalSpecs, CATarget)

  # END CHESS ------------------------

  ResultsTable[iObs, SpecName] = ChessResults$SpecConc
  ResultsTable$Iter[iObs] = ChessResults$Iter
  ResultsTable[iObs, paste0("T.",CompName)] = ChessResults$CalcTotConc

  # ResultsTable$TotCu[iObs] = CalcTotConc[MetalComp]

}

end.time = Sys.time()
end.time - start.time
