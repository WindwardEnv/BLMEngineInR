test_that("GetData works", {
  expect_error(GetData)
  # NComp = as.integer(2)
  # CompNames = c("H","CO3")
  # TestResult = GetData("Test",NComp,CompNames)
  # expect_type(TestResult, "list")
  # expect_identical(names(TestResult),
  #                  c("NObs","obsLabels","totConcObs"))
  # expect_identical(dim(TestResult$obsLabels), c(TestResult$NObs,NComp))
  # expect_identical(dim(TestResult$totConcObs), c(TestResult$NObs,NComp))
})
