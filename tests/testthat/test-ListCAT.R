test_that("ListCAT works", {
  mypfile = system.file(file.path("extdata","ParameterFiles","Cu_full_organic_WATER23dH.dat4"),
                        package = "BLMEngineInR",
                        mustWork = TRUE)
  expect_no_error(ListCAT(ParamFile = mypfile))
  expect_snapshot(ListCAT(ParamFile = mypfile))
})
