test_that("BLM function works - general functionality", {

  mypfile = withr::local_tempfile(fileext = ".dat4")
  myinputfile = withr::local_tempfile(fileext = ".blm4")

  myproblem = carbonate_system_problem
  WriteParamFile(ThisProblem = myproblem, ParamFile = mypfile)

  myinputs = MatchInputsToProblem(
    DFInputs = data.frame(
      ID = "Test",
      Temp = 25,
      pH = 7,
      CO3 = 1E-4
    ),
    ThisProblem = myproblem
  )
  WriteInputFile(AllInput = myinputs,
                 ThisProblem = myproblem,
                 InputFile = myinputfile)

  # for CRAN submission, we do not want any output by default
  expect_silent({
    tmp = BLM(ParamFile = mypfile, InputFile = myinputfile, DoTox = FALSE)
  })
  expect_silent({
    tmp2 = BLM(ThisProblem = myproblem, AllInput = myinputs, DoTox = FALSE)
  })
  tmp$TimeElapsed = "999 secs"
  tmp2$TimeElapsed = "999 secs"
  # tmp

  expect_error(BLM(), "Supply")
  expect_identical(tmp$Miscellaneous$Status[1], "Okay")
  expect_equal(
    log10(tmp$Concentrations[1, "HCO3 (mol/L)"]),
    log10(8.171832e-05),
    tolerance = 0.0001
  )
  expect_equal(tmp, tmp2)

  # expect_silent(BLM(ThisProblem = myproblem, AllInput = myinputs, DoTox = FALSE))

})

test_that("BLM function works - CHESS works", {

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
               log10(testinputsDF$CO3.calcfree),
               tolerance = 0.00001)

})

test_that("BLM function works - toxicity mode works", {

  # skip_on_cran()

  myproblem = Cu_full_organic_problem

  testinputsDF = data.frame(
    ObsNum = 1,
    ID = "Full_Organic",
    ID2 = "Hard ser 10",
    Temp = 15,
    pH = 7.57,
    DOC = 0.01,
    HA = 10,
    Cu = 4.2112e-7,
    Ca = 0.000037427,
    Mg = 6.274425e-05,
    Na = 0.0001375612,
    K = 6.713850e-06,
    SO4 = 9.993587e-05,
    Cl = 6.699012e-06,
    CO3 = 0.0001374837
  )
  testinputs = MatchInputsToProblem(
    DFInputs = testinputsDF,
    ThisProblem = myproblem
  )
  # for CRAN submission, we do not want any output by default
  expect_silent({
    tmp = BLM(ThisProblem = myproblem,
              AllInput = testinputs, DoTox = TRUE, CritAccumIndex = 1L)
  })

  expect_equal(log10(tmp$Concentrations$`T.Cu (mol/L)`),
               log10(7.62202e-10),
               tolerance = 0.05)

})

test_that("BLM function works - DOC binding works", {

  # skip_on_cran()

  myproblem = Cu_full_organic_problem

  testinputsDF = data.frame(
    ObsNum = c(1, 2),
    ID = "focused_input",
    ID2 = c("Obs9","Obs13"),
    Temp = 15,
    pH = 7,
    DOC = c(0.1, 20),
    HA = 10,
    Cu = 1.57366317313442e-09,
    Ca = 4.1141089251479e-05,
    Mg = 5.877298464497e-05,
    Na = 0.000132713191133803,
    K = 6.32939834638139e-06,
    SO4 = 0.000100344120125559,
    Cl = 6.32939834638139e-06,
    CO3 = 0.000132713191133803
  )

  testinputs = MatchInputsToProblem(
    DFInputs = testinputsDF,
    ThisProblem = myproblem
  )
  tmp = BLM(ThisProblem = myproblem,
            AllInput = testinputs, DoTox = FALSE)

  expect_equal(log10(tmp$Concentrations$`Cu (mol/L)`),
               log10(c(4.94586E-11, 1.92656E-13)),
               tolerance = 0.1)
  expect_equal(log10(tmp$Concentrations$`TOrg.Cu (mol/L)`),
               log10(c(1.43252E-09, 1.57324E-09)),
               tolerance = 0.05)

})
