mypfile = system.file(file.path("extdata","ParameterFiles","carbonate_system_only.dat4"),
                      package = "BLMEngineInR",
                      mustWork = TRUE)
myinputfile = system.file(file.path("extdata","InputFiles","carbonate_system_test.blm4"),
                          package = "BLMEngineInR",
                          mustWork = TRUE)
myproblem = DefineProblem(ParamFile = mypfile)
myinputs = GetData(InputFile = myinputfile, ThisProblem = myproblem)


test_that("BLM function works", {

  tmp = BLM(ParamFile = mypfile,
            InputFile = myinputfile,
            DoTox = FALSE)
  tmp2 = BLM(ThisProblem = myproblem, AllInput = myinputs, DoTox = FALSE)

  expect_error(BLM())
  expect_identical(tmp$Miscellaneous$Status[1], "Okay")
  expect_equal(
    log10(tmp$Concentrations[1,"HCO3 (mol/L)"]),
    log10(8.171832e-05),
    tolerance = 0.0001
  )
  expect_equal(tmp, tmp2)

})

test_that("CHESS works", {

  testinputsDF = data.frame(ID = "",
                            Temp = 25,
                            pH = rep(c(7,5,9), each = 3),
                            CO3 = 10^rep(-5:-3, times = 3))
  K2 = carbonate_system_problem$Spec$K[match("HCO3", carbonate_system_problem$Spec$Name)]
  K1K2 = carbonate_system_problem$Spec$K[match("H2CO3", carbonate_system_problem$Spec$Name)]
  K1 = K1K2 / K2
  H = 10^-testinputsDF$pH
  testinputsDF$CO3.calcfree = testinputsDF$CO3 / (1 + K2 * H + K1K2 * H ^ 2)
  testinputs = MatchInputsToProblem(
    DFInputs = testinputsDF,
    ThisProblem = carbonate_system_problem
  )
  myproblem = carbonate_system_problem
  myproblem$Spec$ActCorr = "None"
  tmp = BLM(ThisProblem = myproblem,
            AllInput = testinputs, DoTox = FALSE)
  testinputsDF$CO3.CHESSfree = tmp$Concentrations$`CO3 (mol/L)`

  # with no activity correction, we expect these to match perfectly
  expect_equal(log10(testinputsDF$CO3.CHESSfree),
               log10(testinputsDF$CO3.calcfree))

})
