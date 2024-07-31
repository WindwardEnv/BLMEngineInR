test_that("CheckBLMObject works", {

  # Test Problem
  mypfile = system.file(file.path("extdata","ParameterFiles","carbonate_system_only.dat4"),
                        package = "BLMEngineInR",
                        mustWork = TRUE)
  myinputfile = system.file(file.path("extdata","InputFiles","carbonate_system_test.blm4"),
                            package = "BLMEngineInR",
                            mustWork = TRUE)
  myproblem = DefineProblem(ParamFile = mypfile)
  myinputs = GetData(InputFile = myinputfile, ThisProblem = myproblem)


  # Reference problems
  refproblem = BlankProblem()
  refinputs = BlankInputList(ThisProblem = myproblem)


  # Nothing wrong
  expect_equal(CheckBLMObject(Object = myproblem, Reference = refproblem, BreakOnError = FALSE), character())
  expect_equal(CheckBLMObject(Object = myproblem, Reference = refproblem, BreakOnError = TRUE ), TRUE)
  expect_equal(CheckBLMObject(Object = myinputs, Reference = refinputs, BreakOnError = FALSE), character())
  expect_equal(CheckBLMObject(Object = myinputs, Reference = refinputs, BreakOnError = TRUE ), TRUE)


  # Remove something:
  mytestproblem = myproblem
  mytestproblem$N = NULL
  expect_error(CheckBLMObject(Object = mytestproblem, Reference = refproblem, BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestproblem, Reference = refproblem, BreakOnError = FALSE)), 1)

  mytestinputs = myinputs
  mytestinputs$SysTempKelvinObs = NULL
  expect_error(CheckBLMObject(Object = mytestinputs, Reference = refinputs, BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestinputs, Reference = refinputs, BreakOnError = FALSE)), 1)


  # Wrong type
  mytestproblem = myproblem
  mytestproblem$N[1] = as.numeric(mytestproblem$N[1])
  expect_error(CheckBLMObject(Object = mytestproblem, Reference = refproblem, BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestproblem, Reference = refproblem, BreakOnError = FALSE)), 1)

  mytestproblem = myproblem
  mytestproblem$Comp$Charge = as.numeric(mytestproblem$Comp$Charge)
  expect_error(CheckBLMObject(Object = mytestproblem, Reference = refproblem, BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestproblem, Reference = refproblem, BreakOnError = FALSE)), 1)

  mytestinputs = myinputs
  mytestinputs$SysTempKelvinObs = as.integer(mytestinputs$SysTempKelvinObs)
  expect_error(CheckBLMObject(Object = mytestinputs, Reference = refinputs, BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestinputs, Reference = refinputs, BreakOnError = FALSE)), 1)


  # Missing column or element
  mytestproblem = myproblem
  mytestproblem$Comp$SiteDens = NULL
  expect_error(CheckBLMObject(Object = mytestproblem, Reference = refproblem, BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestproblem, Reference = refproblem, BreakOnError = FALSE)), 1)

  mytestinputs = myinputs
  mytestinputs$TotConcObs = mytestinputs$TotConcObs[, -1]
  expect_error(CheckBLMObject(Object = mytestinputs, Reference = refinputs, BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestinputs, Reference = refinputs, BreakOnError = FALSE)), 1)


  # length mismatch
  mytestproblem = myproblem
  mytestproblem$Comp = rbind(
    mytestproblem$Comp,
    data.frame(Name = "Junk",
               Charge = 0L,
               MCName = "Water",
               MCR = 1L,
               Type = "MassBal",
               ActCorr = "Debye",
               SiteDens = 1.0)
  )
  expect_error(CheckBLMObject(Object = mytestproblem, Reference = refproblem, BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestproblem, Reference = refproblem, BreakOnError = FALSE)), 1)

  mytestinputs = myinputs
  mytestinputs$InLabObs = rbind(mytestinputs$InLabObs, matrix("Junk"))
  expect_error(CheckBLMObject(Object = mytestinputs, Reference = refinputs, BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestinputs, Reference = refinputs, BreakOnError = FALSE)), 1)

  mytestproblem = myproblem
  mytestproblem$Comp = mytestproblem$Comp[-1, ]
  expect_error(CheckBLMObject(Object = mytestproblem, Reference = refproblem, BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestproblem, Reference = refproblem, BreakOnError = FALSE)), 1)

  mytestinputs = myinputs
  mytestinputs$pHObs = mytestinputs$pHObs[-1]
  expect_error(CheckBLMObject(Object = mytestinputs, Reference = refinputs, BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestinputs, Reference = refinputs, BreakOnError = FALSE)), 1)

})
