#' Set the initial guess for component concentrations
#'
#' @param NComp integer; the number of components
#' @param CompName character vector (`NComp`); the names of the components
#' @param TotConc numeric vector (`NComp`); the total concentrations of components
#' @param CompType character vector (`NComp`); the component type
#'
#' @return The initial guesses for the free ion concentrations of components
#'
#' @keywords internal
#'
#' @noRd
RinitialGuess = function(NComp, CompName, TotConc, CompType){


  CompConc = TotConc / 2
  names(CompConc) = CompName

  CompConc[CompType == "FixedConc"] = TotConc[CompType == "FixedConc"]

  return(CompConc)
}

# # Try getting closer to the answer with initial guess...
#
# # Initialize Species Concentrations
# SpecConc = RCalcSpecConc(CompConc = CompConc,
#                          SpecK = SpecK,
#                          SpecStoich = SpecStoich,
#                          SpecName = SpecName,
#                          NComp = NComp,
#                          NSpec = NSpec)
# SpecConc[1:NComp] = CompConc
#
# # Update Total Concentrations for Fixed Activity & Metal
# for (iComp in which((CompType == "FixedConc") | (iTox & (CompName == MetalName)))){
#   TotMoles[iComp] = sum(SpecStoich[,iComp] * (SpecConc * SpecCtoM))
#   TotConc[iComp] = TotMoles[iComp] * SpecCtoM[iComp]
# }
#
# # Calculate Residuals for the first time
# RR = RCalcResidual(
#   NComp = NComp,
#   NSpec = NSpec,
#   SpecConc = SpecConc,
#   SpecStoich = SpecStoich,
#   TotMoles = TotMoles,
#   SpecCtoM = SpecCtoM,
#   CompName = CompName,
#   CompType = CompType,
#   iMetal = iMetal,
#   iBLMetal = iBLMetal,
#   CATarget = CATarget,
#   iTox = iTox
# )
# Resid = RR$Resid
# MaxError = RR$MaxError
# WhichMax = RR$WhichMax
# CalcTotConc = RR$CalcTotConc
# CalcTotMoles = RR$CalcTotMoles
#
# if (MaxError > 1){
#   Iter = 0
#   while((MaxError > 1) & (Iter < 10)){
#     Iter = Iter + 1
#
#     if (iTox & (WhichMax == iMetal)){
#       CompConc[WhichMax] = CompConc[WhichMax] * (CATarget / sum(SpecConc[iBLMetal]))
#     } else {
#       CompConc[WhichMax] = CompConc[WhichMax] * (TotMoles[WhichMax] / CalcTotMoles[WhichMax])
#     }
#
#     SpecConc = RCalcSpecConc(CompConc = CompConc,
#                              SpecK = SpecK,
#                              SpecStoich = SpecStoich,
#                              SpecName = SpecName,
#                              NComp = NComp,
#                              NSpec = NSpec)
#     SpecConc[1:NComp] = CompConc
#
#     # Update Total Concentrations for Fixed Activity & Metal
#     for (iComp in which((CompType == "FixedConc") | (iTox & (CompName == MetalName)))){
#       TotMoles[iComp] = sum(SpecStoich[,iComp] * (SpecConc * SpecCtoM))
#       TotConc[iComp] = TotMoles[iComp] * SpecCtoM[iComp]
#     }
#
#     (RR = RCalcResidual(
#       NComp = NComp,
#       NSpec = NSpec,
#       SpecConc = SpecConc,
#       SpecStoich = SpecStoich,
#       TotMoles = TotMoles,
#       SpecCtoM = SpecCtoM,
#       CompName = CompName,
#       CompType = CompType,
#       iMetal = iMetal,
#       iBLMetal = iBLMetal,
#       CATarget = CATarget,
#       iTox = iTox
#     ))
#     Resid = RR$Resid
#     MaxError = RR$MaxError
#     WhichMax = RR$WhichMax
#     CalcTotConc = RR$CalcTotConc
#     CalcTotMoles = RR$CalcTotMoles
#   }
# }
