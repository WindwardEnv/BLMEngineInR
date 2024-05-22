rm(list = ls())
# devtools::clean_dll()
devtools::load_all()

DoTox = TRUE
iCA = 1L
QuietFlag ="Quiet"
ConvergenceCriteria = 0.0001
MaxIter = 100L
DoPartialStepsAlways = FALSE

ParamFile = "scrap/parameter file format/full_organic_WATER23dH.dat4"
InputFile = "scrap/input files/Cu_Cali_WQC.blm4"
iCA = 2L

ParamFile = "scrap/parameter file format/Zn_full_organic.dat4"
InputFile = "scrap/input files/Zn_Cali_WQC.blm4"
iCA = 2L

ThisProblem = DefineProblem(ParamFile, WriteLog = TRUE)

FunctionInputs = ThisProblem[
  which(names(ThisProblem) %in% formalArgs("GetData"))]
FunctionInputs$InputFile = InputFile
AllInput = do.call("GetData", args = FunctionInputs)

# test stuff
# capture.output(
ResultsTable <- BLM(
  ParamFile = ParamFile,
  InputFile = InputFile,
  DoTox = DoTox,
  iCA = iCA,
  QuietFlag = QuietFlag,
  ConvergenceCriteria = ConvergenceCriteria,
  MaxIter = MaxIter,
  DoPartialStepsAlways = DoPartialStepsAlways
)
# , file = "scrap/debug.txt")
beepr::beep()

save(ResultsTable, file = paste0(InputFile, "_SPEC.RData"))


load("scrap/input files/Zn_Cali_WQC.blm4_TOX.RData")
load("scrap/input files/Zn_Cali_WQC.blm4_SPEC.RData")
load("scrap/input files/Cu_Cali_WQC.blm4_TOX.RData")
load("scrap/input files/Cu_Cali_WQC.blm4_SPEC.RData")
aggregate(
  1:nrow(ResultsTable),
  by = list(
    Converged = !is.na(ResultsTable$FinalMaxError) &
      (ResultsTable$FinalMaxError < ConvergenceCriteria),
    ErraticBehavior = is.na(ResultsTable$FinalMaxError),
    MaxIterationsReached = (ResultsTable$FinalIter == MaxIter)
  ),
  FUN = length
)
# Zn Toxicity:
#   Converged ErraticBehavior MaxIterationsReached    x
# 1      TRUE           FALSE                FALSE 1043
# 2     FALSE            TRUE                FALSE   28
# 3     FALSE           FALSE                 TRUE  223
#
# Zn Speciation
#   Converged ErraticBehavior MaxIterationsReached    x
# 1      TRUE           FALSE                FALSE 1230
# 2     FALSE            TRUE                FALSE    3
# 3     FALSE           FALSE                 TRUE   61
#
# Cu Toxicity:
#   Converged ErraticBehavior MaxIterationsReached   x
# 1      TRUE           FALSE                FALSE 751
# 2     FALSE            TRUE                FALSE  72
# 3     FALSE           FALSE                 TRUE  91
#
# Cu Speciation
#   Converged ErraticBehavior MaxIterationsReached   x
# 1      TRUE           FALSE                FALSE 878
# 2     FALSE           FALSE                 TRUE  36

load("scrap/input files/Zn_Cali_WQC.blm4_TOX.RData")
openxlsx::write.xlsx(ResultsTable[, c(
  "ObsNum",
  "ID",
  "ID2",
  "Temp",
  "pH",
  "DOC",
  "HA",
  "Input.Zn",
  "Input.Ca",
  "Input.Mg",
  "Input.Na",
  "Input.K",
  "Input.SO4",
  "Input.Cl",
  "Input.CO3",
  "FinalIter",
  "FinalToxIter",
  "FinalMaxError",
  "IonicStrength",
  "T.Zn (mol/L)",
  intersect(colnames(ResultsTable),
            c(paste0(ThisProblem$SpecName, " (mol/L)"),
              paste0(ThisProblem$SpecName, " (mol/kg wet)"),
              paste0("TOrg.",ThisProblem$CompName, " (mol/L)")))
)], file = paste0(InputFile, "_TOX.xlsx"))
summary(ResultsTable$FinalMaxError > ConvergenceCriteria)
summary(ResultsTable$FinalIter == MaxIter)
summary(ResultsTable$FinalToxIter == MaxIter)
did_not_converge = which((ResultsTable$FinalIter == MaxIter))
AllInput.dncZn = lapply(
  AllInput,
  FUN = function(X) {
    if (is.null(dim(X))){
      if (length(X) == 1) {
        out = length(did_not_converge)
      } else {
        out = X[did_not_converge]
      }
    } else if (length(dim(X)) == 2) {
      out = X[did_not_converge,]
    } else {
      stop(dim(X))
    }
    out
  }
)
ResultsTable.dncZn <- BLM(
  ParamFile = ParamFile,
  AllInput = AllInput.dncZn,
  DoTox = TRUE,
  iCA = 2L,
  QuietFlag = "Quiet",
  ConvergenceCriteria = 0.0001,
  MaxIter = 1000L,
  DoPartialStepsAlways = DoPartialStepsAlways
)

load("scrap/input files/Cu_Cali_WQC.blm4_TOX.RData")
openxlsx::write.xlsx(ResultsTable[, c(
  "ObsNum",
  "ID",
  "ID2",
  "Temp",
  "pH",
  "DOC",
  "HA",
  "Input.Cu",
  "Input.Ca",
  "Input.Mg",
  "Input.Na",
  "Input.K",
  "Input.SO4",
  "Input.Cl",
  "Input.CO3",
  "FinalIter",
  "FinalToxIter",
  "FinalMaxError",
  "IonicStrength",
  "T.Cu (mol/L)",
  intersect(colnames(ResultsTable),
            c(paste0(ThisProblem$SpecName, " (mol/L)"),
              paste0(ThisProblem$SpecName, " (mol/kg wet)"),
              paste0("TOrg.",ThisProblem$CompName, " (mol/L)")))
)], file = paste0(InputFile, "_TOX.xlsx"))
summary(ResultsTable$FinalMaxError > ConvergenceCriteria)
summary(ResultsTable$FinalIter == MaxIter)
summary(ResultsTable$FinalToxIter == MaxIter)
did_not_converge = which((ResultsTable$FinalIter == MaxIter))
AllInput.dncCu = lapply(
  AllInput,
  FUN = function(X) {
    if (is.null(dim(X))){
      if (length(X) == 1) {
        out = length(did_not_converge)
      } else {
        out = X[did_not_converge]
      }
    } else if (length(dim(X)) == 2) {
      out = X[did_not_converge,]
    } else {
      stop(dim(X))
    }
    out
  }
)
ResultsTable.dncCu <- BLM(
  ParamFile = ParamFile,
  AllInput = AllInput.dncCu,
  DoTox = TRUE,
  iCA = 2L,
  QuietFlag = "Quiet",
  ConvergenceCriteria = 0.0001,
  MaxIter = 1000L,
  DoPartialStepsAlways = DoPartialStepsAlways
)

ResultsTable_3plus6_UpdateZEDalways_zeroNan = ResultsTable
save(ResultsTable_3plus6_UpdateZEDalways_zeroNan, file = "scrap/ResultsTable_3plus6_UpdateZEDalways_zeroNan.RData")


load("scrap/ResultsTable_3plus6_UpdateZEDalways_zeroNan.RData")
load("scrap/ResultsTable_3plus3_UpdateZEDalways_zeroNan.RData")
load("scrap/ResultsTable_3plus3_UpdateZEDalways_BreakNan.RData")
load("scrap/ResultsTable_3plus3_UpdateZEDalways.RData")
load("scrap/ResultsTable_3plus3_UpdateZEDin3.RData")
load("scrap/ResultsTable_3plus3.RData")
load("scrap/ResultsTable_3iters.RData")
load("scrap/ResultsTable_10iters.RData")

# save(ResultsTable_10iters, file = "scrap/ResultsTable_10iters.RData")
# save(ResultsTable_3iters, file = "scrap/ResultsTable_3iters.RData")
# save(ResultsTable_3plus3, file = "scrap/ResultsTable_3plus3.RData")
# save(ResultsTable_3plus3_UpdateZEDin3, file = "scrap/ResultsTable_3plus3_UpdateZEDin3.RData")
# save(ResultsTable_3plus3_UpdateZEDalways, file = "scrap/ResultsTable_3plus3_UpdateZEDalways.RData")
# save(ResultsTable_3plus3_UpdateZEDalways_zeroNan, file = "scrap/ResultsTable_3plus3_UpdateZEDalways_zeroNan.RData")
# save(ResultsTable_3plus3_UpdateZEDalways_BreakNan, file = "scrap/ResultsTable_3plus3_UpdateZEDalways_BreakNan.RData")

ResultsTable_10iters$FinalIter - ResultsTable_3iters$FinalIter
ResultsTable_10iters$FinalMaxError - ResultsTable_3iters$FinalMaxError
ResultsTable_3plus3$FinalIter - ResultsTable_3iters$FinalIter
ResultsTable_3plus3_UpdateZEDin3$FinalIter - ResultsTable_3iters$FinalIter
ResultsTable_3plus3_UpdateZEDalways$FinalIter - ResultsTable_3iters$FinalIter
ResultsTable_3plus3_UpdateZEDalways_BreakNan$FinalIter - ResultsTable_3plus3_UpdateZEDalways$FinalIter
ResultsTable_3plus3_UpdateZEDalways_BreakNan$FinalMaxError - ResultsTable_3plus3_UpdateZEDalways$FinalMaxError
ResultsTable_3plus3_UpdateZEDalways_zeroNan$FinalIter - ResultsTable_3plus3_UpdateZEDalways$FinalIter
ResultsTable_3plus3_UpdateZEDalways_zeroNan$FinalMaxError - ResultsTable_3plus3_UpdateZEDalways$FinalMaxError
ResultsTable_3plus6_UpdateZEDalways_zeroNan$FinalIter - ResultsTable_3plus3_UpdateZEDalways_zeroNan$FinalIter
ResultsTable_3plus6_UpdateZEDalways_zeroNan$FinalMaxError - ResultsTable_3plus3_UpdateZEDalways_zeroNan$FinalMaxError
