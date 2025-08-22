test_that("ListCAT works", {

  # only one value
  mypfile1 = withr::local_tempfile(fileext = ".dat4")
  myCAT1 = data.frame(
    Num = 1,
    CA = 0.1234,
    Species = "Test species",
    TestType = "Test test type",
    Duration = "Test duration",
    Lifestage = "test ls",
    Endpoint = "test endpoint",
    Quantifier = "test quant",
    References = "test ref",
    Miscellaneous = "test misc"
  )
  colnames(myCAT1)[c(2, 4)] = c("CA (nmol/gw)", "Test Type")

  WriteParamFile(
    ThisProblem = AddCriticalValues(
      ThisProblem = carbonate_system_problem,
      CATab = myCAT1),
    ParamFile = mypfile1
  )
  expect_no_error({LC1 = ListCAT(ParamFile = mypfile1)})
  expect_equal(LC1, myCAT1)


  # Multiple values
  mypfile999 = withr::local_tempfile(fileext = ".dat4")
  myCAT999 = data.frame(
    Num = 1:999,
    CA = rlnorm(999),
    Species = paste0("Test species", 1:999),
    TestType = paste0("Test test type", 1:999),
    Duration = paste0("Test duration", 1:999),
    Lifestage = paste0("test ls", 1:999),
    Endpoint = paste0("test endpoint", 1:999),
    Quantifier = paste0("test quant", 1:999),
    References = paste0("test ref", 1:999),
    Miscellaneous = paste0("test misc", 1:999)
  )
  colnames(myCAT999)[c(2, 4)] = c("CA (nmol/gw)", "Test Type")
  WriteParamFile(
    ThisProblem = AddCriticalValues(
      ThisProblem = carbonate_system_problem,
      CATab = myCAT999),
    ParamFile = mypfile999
  )
  expect_no_error({LC999 = ListCAT(ParamFile = mypfile999)})
  expect_equal(LC999, myCAT999)

})
