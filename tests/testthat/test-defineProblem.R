test_that("DefineProblem works", {

  mypfile = system.file(file.path("extdata","ParameterFiles","Cu_full_organic_WATER23dH.dat4"),
                        package = "BLMEngineInR",
                        mustWork = TRUE)
  expect_no_error(DefineProblem(ParamFile = mypfile))

  myproblem = DefineProblem(ParamFile = mypfile)
  myproblem$WHAM$WdatFile = gsub(dirname(myproblem$WHAM$WdatFile), "", myproblem$WHAM$WdatFile)

  mytestproblem = Cu_full_organic_problem
  mytestproblem$WHAM$WdatFile = gsub(dirname(mytestproblem$WHAM$WdatFile), "", mytestproblem$WHAM$WdatFile)

  expect_equal(myproblem[names(myproblem) != "ParamFile"],
               mytestproblem[names(mytestproblem) != "ParamFile"])

})
