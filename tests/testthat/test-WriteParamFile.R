test_that("WriteParamFile works", {
  mypfile = system.file(
    file.path("extdata","ParameterFiles","carbonate_system_only.dat4"),
    package = "BLMEngineInR",
    mustWork = TRUE)
  myproblem = DefineProblem(ParamFile = mypfile)

  mytemppfile = tempfile()
  expect_no_error(WriteParamFile(ThisProblem = myproblem, ParamFile = mytemppfile))
  mtempproblem = DefineProblem(ParamFile = mytemppfile)

  expect_equal(mtempproblem, myproblem)
})
