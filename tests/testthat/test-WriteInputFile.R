test_that("WriteInputFile works", {

  mypfile = system.file("extdata","ParameterFiles","Cu_full_organic.dat4",
                        package = "BLMEngineInR",
                        mustWork = TRUE)
  myproblem = DefineProblem(ParamFile = mypfile)
  myinputfile = system.file("extdata","InputFiles","Cu_full_organic.blm4",
                            package = "BLMEngineInR",
                            mustWork = TRUE)
  myinputs = GetData(InputFile = myinputfile, ThisProblem = myproblem)

  mytestinputfile = withr::local_tempfile(fileext = ".blm4")

  expect_no_error(WriteInputFile(AllInput = myinputs, ThisProblem = myproblem, InputFile = mytestinputfile))

  expect_equal(scan(mytestinputfile, what = character(), quiet = TRUE, sep = ",", strip.white = TRUE)[1:13],
               scan(myinputfile, what = character(), quiet = TRUE, sep = ",", strip.white = TRUE)[1:13])



})
