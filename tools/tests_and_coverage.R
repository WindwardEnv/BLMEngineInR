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

devtools::load_all()
devtools::unload()
devtools::test_coverage()


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


