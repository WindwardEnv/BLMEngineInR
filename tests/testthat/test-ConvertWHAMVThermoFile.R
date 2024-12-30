test_that("ConvertWHAMVThermoFile works", {
  mydbsfile = "C:/Users/kellyc/Documents/BLM/thermo_databases/water23.dbs"
  testthat::skip_if_not(file.exists(mydbsfile))
  expect_no_error(ConvertWHAMVThermoFile(ThermoDBSName = mydbsfile))

  # testthat::skip_on_cran()
  converted.file = ConvertWHAMVThermoFile(ThermoDBSName = mydbsfile)
  previously.converted.file = DefineProblem(
    ParamFile = system.file(file.path("extdata","ParameterFiles", "All_Water23_reactions.dat4"),
                            package = "BLMEngineInR",
                            mustWork = TRUE))
  compare.names = setdiff(names(previously.converted.file), c("ParamFile","WHAM"))
  expect_equal(object = converted.file[compare.names],
               expected = previously.converted.file[compare.names])

  # testthat::skip_on_cran()
  compare.names = setdiff(names(previously.converted.file), c("Ver","File"))
  expect_equal(object = converted.file$WHAM[compare.names],
               expected = previously.converted.file$WHAM[compare.names])
})
