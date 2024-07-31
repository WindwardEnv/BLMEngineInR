test_that("GetData works", {

  expect_error(GetData())

  mypfile = system.file(file.path("extdata","ParameterFiles","carbonate_system_only.dat4"),
                        package = "BLMEngineInR",
                        mustWork = TRUE)
  myinputfile = system.file(file.path("extdata","InputFiles","carbonate_system_test.blm4"),
                            package = "BLMEngineInR",
                            mustWork = TRUE)
  myproblem = DefineProblem(ParamFile = mypfile)
  myinputs = GetData(InputFile = myinputfile, ThisProblem = myproblem)
  expect_equal(CheckBLMObject(myinputs, BlankInputList(myproblem), FALSE),
               character())

})
