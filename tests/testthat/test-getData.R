test_that("getData works", {
  expect_error(getData())
  # NComp = as.integer(2)
  # CompNames = c("H","CO3")
  # TestResult = getData("Test",NComp,CompNames)
  # expect_type(TestResult, "list")
  # expect_identical(names(TestResult),
  #                  c("NObs","obsLabels","totConcObs"))
  # expect_identical(dim(TestResult$obsLabels), c(TestResult$NObs,NComp))
  # expect_identical(dim(TestResult$totConcObs), c(TestResult$NObs,NComp))
})
