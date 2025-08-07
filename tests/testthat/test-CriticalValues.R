test_that("AddCriticalValues works", {
  expect_no_error(AddCriticalValues(
    ThisProblem = carbonate_system_problem,
    CA = 12345,
    Species = "A. species",
    TestType = "Acute",
    Duration = "24h",
    Lifestage = "adult",
    Endpoint = "survival",
    Quantifier = "LC50",
    References = "thin air",
    Miscellaneous = "individual data point"
  ))
  expect_no_error(AddCriticalValues(
    ThisProblem = carbonate_system_problem,
    CATab = data.frame(
      CA = 12345,
      Species = "A. species",
      TestType = "Acute",
      Duration = "24h",
      Lifestage = "adult",
      Endpoint = "survival",
      Quantifier = "LC50",
      References = "thin air",
      Miscellaneous = "individual data point"
    )))
  expect_no_error(AddCriticalValues(
    ThisProblem = carbonate_system_problem,
    CATab = data.frame(
      CA = 12345,
      Species = "A. species",
      TestType = "Acute",
      Endpoint = "survival",
      References = "thin air"
    )))
})
test_that("RemoveCriticalValues works", {

  mypfile = system.file("extdata","ParameterFiles","Cu_full_organic.dat4",
                        package = "BLMEngineInR",
                        mustWork = TRUE)
  myproblem = DefineProblem(ParamFile = mypfile)

  expect_no_error(RemoveCriticalValues(ThisProblem = myproblem, CAToRemove = 1))
  expect_error(RemoveCriticalValues(ThisProblem = myproblem, CAToRemove = -1))
  expect_error(RemoveCriticalValues(ThisProblem = myproblem, CAToRemove = 999))

})
