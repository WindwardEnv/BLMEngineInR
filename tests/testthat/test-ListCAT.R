test_that("ListCAT works", {
  mypfile = system.file("extdata","ParameterFiles","Cu_full_organic.dat4",
                        package = "BLMEngineInR",
                        mustWork = TRUE)
  expect_no_error(ListCAT(ParamFile = mypfile))
  expect_snapshot(ListCAT(ParamFile = mypfile))
})
