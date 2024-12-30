test_that("WriteParamFile works", {

  mypfile = system.file(
    file.path("extdata","ParameterFiles","Cu_full_organic.dat4"),
    package = "BLMEngineInR",
    mustWork = TRUE)
  myproblem = DefineProblem(ParamFile = mypfile)
  myproblem = AddPhases(ThisProblem = myproblem,
                        PhaseEquation = "CO2(g) = 1 * CO3 + 2 * H",
                        PhaseLogK = -1.5,
                        PhaseDeltaH = 0,
                        PhaseTempKelvin = 298,
                        PhaseMoles = 10^-3.5)

  mytemppfile = withr::local_tempfile()
  expect_no_error(WriteParamFile(ThisProblem = myproblem, ParamFile = mytemppfile))

  mytempproblem = DefineProblem(ParamFile = mytemppfile)
  compare.names = setdiff(names(myproblem), c("ParamFile", "WHAM"))
  expect_equal(mytempproblem[compare.names], myproblem[compare.names])
  compare.names = setdiff(names(myproblem$WHAM), c("Ver", "File"))
  expect_equal(mytempproblem$WHAM[compare.names], myproblem$WHAM[compare.names])

})

