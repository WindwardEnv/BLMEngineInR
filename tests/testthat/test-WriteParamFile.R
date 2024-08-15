test_that("WriteParamFile works", {
  mypfile = system.file(
    file.path("extdata","ParameterFiles","Cu_full_organic_WATER23dH.dat4"),
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
  expect_equal(mytempproblem[names(mytempproblem) != "ParamFile"],
               myproblem[names(myproblem) != "ParamFile"])
})
