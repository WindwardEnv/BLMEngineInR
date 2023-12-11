ParamFile = "scrap/parameter file format/full_organic.dat"
InputFile = "scrap/parameter file format/full_organic.blm4"

thisProblem = defineProblem(ParamFile = ParamFile)
allInput = do.call("GetData", args = c(
  thisProblem[names(thisProblem) %in% formalArgs("GetData")],
  list(InputFile = InputFile)
))



# Save some common variables for initializing arrays
NComp = thisProblem$NComp
CompName = thisProblem$CompName
NSpec = thisProblem$NSpec
SpecName = thisProblem$SpecName

# Initialize the output array
OrgSummaryName = c("DOC","TOrg.H","TOrg.Cu","TOrg.Ca","TOrg.Mg","TOrg.Na","TOrg.K")
OrgSummarySpec = list(
  DOC = SpecName[grepl("^DOC-FA", SpecName) & !(grepl("H$",SpecName) |
                                                  grepl("Cu",SpecName) |
                                                  grepl("Ca",SpecName) |
                                                  grepl("Mg",SpecName) |
                                                  grepl("Na",SpecName) |
                                                  grepl("K",SpecName))],
  TOrg.H = SpecName[grepl("^DOC-FA", SpecName) & grepl("[[:digit:]]H$",SpecName)],
  TOrg.Cu = SpecName[grepl("^DOC-FA", SpecName) & grepl("Cu",SpecName)],#Cu2+ and CuOH+
  # TOrg.Cu = SpecName[grepl("^DOC-FA", SpecName) & grepl("Cu$",SpecName)],#Cu2+ onlt
  TOrg.Ca = SpecName[grepl("^DOC-FA", SpecName) & grepl("Ca",SpecName)],
  TOrg.Mg = SpecName[grepl("^DOC-FA", SpecName) & grepl("Mg",SpecName)],
  TOrg.Na = SpecName[grepl("^DOC-FA", SpecName) & grepl("Na",SpecName)],
  TOrg.K = SpecName[grepl("^DOC-FA", SpecName) & grepl("K",SpecName)]
)
NOrgSummary = length(OrgSummaryName)
out = array(numeric(1), dim = c(allInput$NObs, NSpec + NOrgSummary),
            dimnames = list(1:allInput$NObs, c(SpecName,OrgSummaryName)))

# Initialize thisInput as thisProblem, with one observation's worth of
# concentrations
thisInput = thisProblem
thisInput$InLab = array(character(thisProblem$NInLab),
                            dimnames = list(thisProblem$InLabName))
thisInput$TotConc = array(numeric(NComp), dimnames = list(CompName))
thisInput$CompConc = array(numeric(NComp), dimnames = list(CompName))
thisInput$SpecConc = array(numeric(NSpec), dimnames = list(SpecName))

# Loop through each observation
# for (iObs in 1:allInput$NObs){
iObs = 1; {
  thisInput$InLab = allInput$InLabObs[iObs,]
  thisInput$TotConc = allInput$TotConcObs[iObs,]

  # For now, we're going to use test data, setting the initial "guess" to the
  # actual component free ion concentrations
  # thisInput$CompConc = do.call(initialGuess, args = thisInput[formalArgs(initialGuess)])
  thisInput$CompConc[c(thisProblem$InCompName,"H","BL1")] =
    Full_OrganicDataFreeConc[c(thisProblem$InCompName,"H","BL1")]
  thisInput$CompConc[thisProblem$DefCompName[grepl("DOC",thisProblem$DefCompName)]] =
    Full_OrganicDataFreeConc["DOC"] *
    thisProblem$DefCompSiteDens[grepl("DOC",thisProblem$DefCompName)]/1000

  thisInput$LogCompConc = log10(thisInput$CompConc)

  # 3. Run the speciation problem
  #   --> R variable defining problem from step 1
  #   --> R variable with inputs from step 2
  #   <-- R variable with speciation outputs
  out[iObs,SpecName] = 10^do.call("CppCalcLogSpecConc",
                                  args = thisInput[formalArgs("CppCalcLogSpecConc")])
  #organics
  for (i in OrgSummaryName){
    out[iObs, i] = sum(out[iObs,OrgSummarySpec[[i]]])
  }

  tiff("scrap/test CalcSpecConc with Organic problem.tif")
  plot(x = Full_OrganicDataFreeConc, y = out[iObs, names(Full_OrganicDataFreeConc)],
       xlab = "Current BLM", ylab = "BLM In R Function",
       main = "reverse calculation (non-thermodynamic effects not included)",
       log = "xy", pch = 20, xlim = c(10^-15, 10^0), ylim = c(10^-15, 10^0))
  text(x = Full_OrganicDataFreeConc, y = out[iObs,names(Full_OrganicDataFreeConc)],
       labels = names(Full_OrganicDataFreeConc), pos = 4, cex = 0.75)
  abline(a = 0, b = 1)
  dev.off()

}


cbind(
  thisProblem$SpecStoich[c("Cu","H","BL1","DOC-FA_1H",
                           "OH","CuOH","BL1-Cu","BL1-CuOH","DOC-FA_1-Cu","DOC-FA_1-CuOH"),
                         c("Cu","H","BL1","DOC-FA_1H")],
  thisProblem$SpecLogK[c("Cu","H","BL1","DOC-FA_1H",
                         "OH","CuOH","BL1-Cu","BL1-CuOH","DOC-FA_1-Cu","DOC-FA_1-CuOH")]
)
