test_that("CheckBLMObject works when nothing's wrong", {

  myproblem = carbonate_system_problem
  mylistproblem = ConvertToList(myproblem)
  myinputsdf = data.frame(
    ID = "Test",
    Temp = 25,
    pH = 7,
    CO3 = 1E-4
  )
  myinputs = MatchInputsToProblem(DFInputs = myinputsdf, ThisProblem = myproblem)

  # Nothing wrong
  expect_equal(CheckBLMObject(Object = myproblem, Reference = BlankProblem(), BreakOnError = FALSE), character())
  expect_equal(CheckBLMObject(Object = myproblem, Reference = BlankProblem(), BreakOnError = TRUE ), TRUE)
  expect_equal(CheckBLMObject(Object = mylistproblem, Reference = BlankProblemList(), BreakOnError = FALSE), character())
  expect_equal(CheckBLMObject(Object = mylistproblem, Reference = BlankProblemList(), BreakOnError = TRUE ), TRUE)
  expect_equal(CheckBLMObject(Object = myinputs, Reference = BlankInputList(ThisProblem = myproblem), BreakOnError = FALSE), character())
  expect_equal(CheckBLMObject(Object = myinputs, Reference = BlankInputList(ThisProblem = myproblem), BreakOnError = TRUE ), TRUE)

})
test_that("CheckBLMObject works when an element's been removed or added", {

  myproblem = carbonate_system_problem
  mylistproblem = ConvertToList(myproblem)
  myinputsdf = data.frame(
    ID = "Test",
    Temp = 25,
    pH = 7,
    CO3 = 1E-4
  )
  myinputs = MatchInputsToProblem(DFInputs = myinputsdf, ThisProblem = myproblem)

  # Remove something:
  mytestproblem = myproblem
  mytestproblem$N = NULL
  expect_error(CheckBLMObject(Object = mytestproblem, Reference = BlankProblem(), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestproblem, Reference = BlankProblem(), BreakOnError = FALSE)), 1)

  mytestlistproblem = mylistproblem
  mytestlistproblem$NMass = NULL
  expect_error(CheckBLMObject(Object = mytestlistproblem, Reference = BlankProblemList(), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestlistproblem, Reference = BlankProblemList(), BreakOnError = FALSE)), 1)

  mytestinputs = myinputs
  mytestinputs$SysTempKelvinObs = NULL
  expect_error(CheckBLMObject(Object = mytestinputs, Reference = BlankInputList(ThisProblem = myproblem), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestinputs, Reference = BlankInputList(ThisProblem = myproblem), BreakOnError = FALSE)), 1)


})
test_that("CheckBLMObject works when there's a wrong type", {

  myproblem = carbonate_system_problem
  mylistproblem = ConvertToList(myproblem)
  myinputsdf = data.frame(
    ID = "Test",
    Temp = 25,
    pH = 7,
    CO3 = 1E-4
  )
  myinputs = MatchInputsToProblem(DFInputs = myinputsdf, ThisProblem = myproblem)

  # Wrong type
  mytestproblem = myproblem
  mytestproblem$N[1] = as.numeric(mytestproblem$N[1])
  expect_error(CheckBLMObject(Object = mytestproblem, Reference = BlankProblem(), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestproblem, Reference = BlankProblem(), BreakOnError = FALSE)), 1)

  mytestproblem = myproblem
  mytestproblem$Comp$Charge = as.numeric(mytestproblem$Comp$Charge)
  expect_error(CheckBLMObject(Object = mytestproblem, Reference = BlankProblem(), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestproblem, Reference = BlankProblem(), BreakOnError = FALSE)), 1)

  mytestlistproblem = mylistproblem
  mytestlistproblem$NMass = as.numeric(mytestlistproblem$NMass)
  expect_error(CheckBLMObject(Object = mytestlistproblem, Reference = BlankProblemList(), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestlistproblem, Reference = BlankProblemList(), BreakOnError = FALSE)), 1)

  mytestinputs = myinputs
  mytestinputs$SysTempKelvinObs = as.integer(mytestinputs$SysTempKelvinObs)
  expect_error(CheckBLMObject(Object = mytestinputs, Reference = BlankInputList(ThisProblem = myproblem), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestinputs, Reference = BlankInputList(ThisProblem = myproblem), BreakOnError = FALSE)), 1)

})
test_that("CheckBLMObject works when a column in a matrix or data frame's been removed or added", {

  myproblem = carbonate_system_problem
  mylistproblem = ConvertToList(myproblem)
  myinputsdf = data.frame(
    ID = "Test",
    Temp = 25,
    pH = 7,
    CO3 = 1E-4
  )
  myinputs = MatchInputsToProblem(DFInputs = myinputsdf, ThisProblem = myproblem)

  # problem data.frame missing a column
  mytestproblem = myproblem
  mytestproblem$Comp$SiteDens = NULL
  expect_error(CheckBLMObject(Object = mytestproblem, Reference = BlankProblem(), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestproblem, Reference = BlankProblem(), BreakOnError = FALSE)), 1)

  # problem data.frame has extra column - we don't care about this because it
  # could be extra information that the user wants, and it will never be a
  # problem for our functions
  mytestproblem = myproblem
  mytestproblem$Comp$Junk = "Junk"
  expect_no_error(CheckBLMObject(Object = mytestproblem, Reference = BlankProblem(), BreakOnError = TRUE))

  # problem matrix missing a column
  mytestproblem = myproblem
  mytestproblem$SpecStoich = mytestproblem$SpecStoich[, -1, drop = FALSE]
  expect_error(CheckBLMObject(Object = mytestproblem, Reference = BlankProblem(), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestproblem, Reference = BlankProblem(), BreakOnError = FALSE)), 1)

  # listproblem matrix missing a column
  mytestlistproblem = mylistproblem
  mytestlistproblem$SpecStoich = mytestlistproblem$SpecStoich[, -1, drop = FALSE]
  expect_error(CheckBLMObject(Object = mytestlistproblem, Reference = BlankProblemList(), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestlistproblem, Reference = BlankProblemList(), BreakOnError = FALSE)), 1)

  # problem matrix has extra column
  mytestproblem = myproblem
  mytestproblem$SpecStoich = cbind(mytestproblem$SpecStoich, matrix(-999L, nrow = myproblem$N["Spec"], ncol = 1, dimnames = list(NULL, "Junk")))
  expect_error(CheckBLMObject(Object = mytestproblem, Reference = BlankProblem(), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestproblem, Reference = BlankProblem(), BreakOnError = FALSE)), 1)

  # listproblem matrix has extra column
  mytestlistproblem = mylistproblem
  mytestlistproblem$SpecStoich = cbind(mytestlistproblem$SpecStoich, matrix(-999L, nrow = mylistproblem$NSpec, ncol = 1, dimnames = list(NULL, "Junk")))
  expect_error(CheckBLMObject(Object = mytestlistproblem, Reference = BlankProblemList(), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestlistproblem, Reference = BlankProblemList(), BreakOnError = FALSE)), 1)

  # inputs matrix missing a column
  mytestinputs = myinputs
  mytestinputs$TotConcObs = mytestinputs$TotConcObs[, -1, drop = FALSE]
  expect_error(CheckBLMObject(Object = mytestinputs, Reference = BlankInputList(ThisProblem = myproblem), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestinputs, Reference = BlankInputList(ThisProblem = myproblem), BreakOnError = FALSE)), 1)

  # inputs matrix has extra column
  mytestinputs = myinputs
  mytestinputs$TotConcObs = cbind(mytestinputs$TotConcObs, matrix(-999, nrow = myinputs$NObs, ncol = 1, dimnames = list(NULL, "Junk")))
  expect_error(CheckBLMObject(Object = mytestinputs, Reference = BlankInputList(ThisProblem = myproblem), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestinputs, Reference = BlankInputList(ThisProblem = myproblem), BreakOnError = FALSE)), 1)

})
test_that("CheckBLMObject works when something's the wrong length", {

  myproblem = carbonate_system_problem
  mylistproblem = ConvertToList(myproblem)
  myinputsdf = data.frame(
    ID = "Test",
    Temp = 25,
    pH = 7,
    CO3 = 1E-4
  )
  myinputs = MatchInputsToProblem(DFInputs = myinputsdf, ThisProblem = myproblem)

  # problem data.frame has extra row
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
  expect_error(CheckBLMObject(Object = mytestproblem, Reference = BlankProblem(), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestproblem, Reference = BlankProblem(), BreakOnError = FALSE)), 1)

  # problem list array has extra element
  mytestlistproblem = mylistproblem
  mytestlistproblem$CompCharge = c(mytestlistproblem$CompCharge, -999L)
  expect_error(CheckBLMObject(Object = mytestlistproblem, Reference = BlankProblemList(), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestlistproblem, Reference = BlankProblemList(), BreakOnError = FALSE)), 1)

  # problem list matrix has extra row
  mytestlistproblem = mylistproblem
  mytestlistproblem$SpecStoich = rbind(mytestlistproblem$SpecStoich, matrix(-999L, nrow = 1, ncol = mylistproblem$NComp))
  expect_error(CheckBLMObject(Object = mytestlistproblem, Reference = BlankProblemList(), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestlistproblem, Reference = BlankProblemList(), BreakOnError = FALSE)), 1)

  # inputs array has extra row
  mytestinputs = myinputs
  mytestinputs$InLabObs = rbind(mytestinputs$InLabObs, matrix("Junk"))
  expect_error(CheckBLMObject(Object = mytestinputs, Reference = BlankInputList(ThisProblem = myproblem), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestinputs, Reference = BlankInputList(ThisProblem = myproblem), BreakOnError = FALSE)), 1)

  # problem data.frame missing a row
  mytestproblem = myproblem
  mytestproblem$Comp = mytestproblem$Comp[-1, ]
  expect_error(CheckBLMObject(Object = mytestproblem, Reference = BlankProblem(), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestproblem, Reference = BlankProblem(), BreakOnError = FALSE)), 1)

  # problem list array missing a row
  mytestlistproblem = mylistproblem
  mytestlistproblem$CompCharge = mytestlistproblem$CompCharge[-1]
  expect_error(CheckBLMObject(Object = mytestlistproblem, Reference = BlankProblemList(), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestlistproblem, Reference = BlankProblemList(), BreakOnError = FALSE)), 1)

  # problem list matrix missing a row
  mytestlistproblem = mylistproblem
  mytestlistproblem$SpecStoich = mytestlistproblem$SpecStoich[-1, , drop = FALSE]
  expect_error(CheckBLMObject(Object = mytestlistproblem, Reference = BlankProblemList(), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestlistproblem, Reference = BlankProblemList(), BreakOnError = FALSE)), 1)

  # inputs array missing row
  mytestinputs = myinputs
  mytestinputs$pHObs = mytestinputs$pHObs[-1]
  expect_error(CheckBLMObject(Object = mytestinputs, Reference = BlankInputList(ThisProblem = myproblem), BreakOnError = TRUE))
  expect_equal(length(CheckBLMObject(Object = mytestinputs, Reference = BlankInputList(ThisProblem = myproblem), BreakOnError = FALSE)), 1)

})
test_that("CheckBLMObject works when there's something extra", {

  myproblem = carbonate_system_problem
  mylistproblem = ConvertToList(myproblem)
  myinputsdf = data.frame(
    ID = "Test",
    Temp = 25,
    pH = 7,
    CO3 = 1E-4
  )
  myinputs = MatchInputsToProblem(DFInputs = myinputsdf, ThisProblem = myproblem)

  mymegaproblemandinputs = c(myproblem, mylistproblem, myinputs)
  mymegaDFproblemandinputs = c(myproblem, myinputs)
  mymegalistproblemandinputs = c(mylistproblem, myinputs)

  # Nothing wrong
  expect_no_error(CheckBLMObject(Object = mymegaDFproblemandinputs, Reference = BlankProblem(), BreakOnError = TRUE))
  expect_no_error(CheckBLMObject(Object = mymegalistproblemandinputs, Reference = BlankProblemList(), BreakOnError = TRUE))
  expect_no_error(CheckBLMObject(Object = mymegaDFproblemandinputs, Reference = BlankInputList(ThisProblem = myproblem), BreakOnError = TRUE))
  expect_no_error(CheckBLMObject(Object = mymegalistproblemandinputs, Reference = BlankInputList(ThisProblem = myproblem), BreakOnError = TRUE))

})
