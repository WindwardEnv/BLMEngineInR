test_that("DefineProblem works", {
  expect_error(DefineProblem())
  # TestResult = DefineProblem("Test")
  # expect_type(TestResult, "list")
  # expect_identical(names(TestResult),
  #                  c("NComp", "NSpec","CompNames","SpecNames",
  #                    "K", "logK", "Stoich", "CConc"))
  # expect_identical(TestResult$NComp, length(TestResult$CConc))
  # expect_identical(TestResult$NSpec, length(TestResult$K))
  # expect_identical(TestResult$NSpec, length(TestResult$logK))
  # expect_identical(c(TestResult$NSpec, TestResult$NComp),
  #                  dim(TestResult$Stoich))
})
