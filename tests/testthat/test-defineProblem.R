test_that("DefineProblem works", {

  expect_error(DefineProblem())

  mypfile = system.file(file.path("extdata","ParameterFiles","carbonate_system_only.dat4"),
                        package = "BLMEngineInR",
                        mustWork = TRUE)
  myproblem = DefineProblem(ParamFile = mypfile)
  expect_equal(CheckBLMObject(myproblem, BlankProblem(), FALSE),
               character())

})
