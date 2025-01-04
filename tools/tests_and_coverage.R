rm(list = ls())
devtools::load_all()
sapply(list.files("data-raw", full.names = TRUE), FUN = source)

files_with_tests = c(
  "BLM",
  "BlankProblem",
  "DefineProblem",
  "CheckBLMObject",
  "CHESSLog",
  "Components",
  "ConvertWHAMVThermoFile",
  "ConvertWindowsParamFile",
  "CriticalValues",
  "DefineProblem",
  "DefineWHAM",
  "ExpandWHAM",
  "GetData",
  "InLabs",
  "InVars",
  "ListCAT",
  "MassCompartments",
  "Phases",
  "ProblemConversionFunctions",
  "SpecialDefs",
  "Species",
  "StoichConversionFunctions",
  "WriteDetailedFile",
  "WriteInputFile",
  "WriteParamFile",
  "WriteWHAMFile"
)

for (i in 1:length(files_with_tests)) {
  print(files_with_tests[i])
  devtools::test_file(file = paste0("tests/testthat/test-", files_with_tests[i], ".R"))
  devtools::test_coverage_file(file = paste0("R/", files_with_tests[i], ".R"))
  tmp = scan()
}

devtools::test()
covr::file_coverage(
  source_files = paste0("R/", files_with_tests, ".R"),
  test_files = paste0("tests/testthat/test-", files_with_tests, ".R")
)

devtools::unload()
devtools::test_coverage()
devtools::load_all()


covr::file_coverage(
  source_files = list.files(".", pattern = "R", recursive = T),
  test_files = list.files("./tests/testthat", pattern="test-"),
)

devtools::test_coverage_active_file("src/CompUpdate.cpp")#88.89%
devtools::test_coverage_active_file("src/SimpleAdjustComp.cpp")#85.71%
devtools::test_coverage_active_file("src/UpdateFixedComps.cpp")#75.00%
devtools::test_coverage_active_file("src/CalcStep.cpp")#72.55%
devtools::test_coverage_active_file("src/AdjustForToxMode.cpp")#71.43%
devtools::test_coverage_active_file("src/CHESS.cpp")#68.79#
devtools::test_coverage_active_file("src/AdjustForWHAM.cpp")#56.90%
devtools::test_coverage_active_file("src/Jacobian.cpp")#41.11%
devtools::test_coverage_active_file("src/RcppArmaHelper.cpp")#4.30%
devtools::test_coverage_active_file("src/UpdateTotals.cpp")#0.00%
devtools::test_coverage_active_file("src/NumericalJacobian.cpp")#0.00%


lintr::lint("R/BlankProblem.R")
lintr::lint("R/BlankWHAM.R")
lintr::lint("R/BLM.R")
lintr::lint("R/CheckBLMObject.R")
lintr::lint("R/CHESSLog.R")
lintr::lint("R/CommonParameterDefinitions.R")
lintr::lint("R/Components.R")
lintr::lint("R/ConvertWHAMVThermoFile.R")
lintr::lint("R/ConvertWindowsParamFile.R")
lintr::lint("R/CriticalValues.R")
lintr::lint("R/data.R")
lintr::lint("R/DefineProblem.R")
lintr::lint("R/DefineWHAM.R")
lintr::lint("R/ExpandWHAM.R")
lintr::lint("R/GetData.R")
lintr::lint("R/InLabs.R")
lintr::lint("R/InVars.R")
lintr::lint("R/ListCAT.R")
lintr::lint("R/MassCompartments.R")
lintr::lint("R/Phases.R")
lintr::lint("R/ProblemConversionFunctions.R")
lintr::lint("R/SpecialDefs.R")
lintr::lint("R/Species.R")
lintr::lint("R/StoichConversionFunctions.R")
lintr::lint("R/WriteDetailedFile.R")
lintr::lint("R/WriteInputFile.R")
lintr::lint("R/WriteParamFile.R")
lintr::lint("R/WriteWHAMFile.R")
lintr::lint_package()

?styler
styler::style_pkg(style = styler::tidyverse_style(scope = c("indention", "line_breaks")), dry = "on")
