test_that("DefineProblem works", {

  mypfile = system.file("extdata","ParameterFiles","Cu_full_organic.dat4",
                        package = "BLMEngineInR",
                        mustWork = TRUE)
  expect_no_error(DefineProblem(ParamFile = mypfile))

  compare_names = names(Cu_full_organic_problem)[names(Cu_full_organic_problem) %in%
                                                   c("ParamFile","WHAM") == FALSE]
  compare_names_WHAM = names(Cu_full_organic_problem$WHAM)[
    names(Cu_full_organic_problem$WHAM) %in% "File" == FALSE]

  expect_equal(DefineProblem(ParamFile = mypfile)[compare_names],
               Cu_full_organic_problem[compare_names])
  expect_equal(DefineProblem(ParamFile = mypfile)$WHAM[compare_names_WHAM],
               Cu_full_organic_problem$WHAM[compare_names_WHAM])

  myproblem = DefineProblem(ParamFile = mypfile)
  myproblem$WHAM$File = basename(myproblem$WHAM$File)

  mytestproblem = Cu_full_organic_problem
  mytestproblem$WHAM$File = basename(mytestproblem$WHAM$File)

  expect_equal(myproblem[compare_names],
               mytestproblem[compare_names])
  expect_equal(myproblem$WHAM[compare_names_WHAM],
               mytestproblem$WHAM[compare_names_WHAM])

})



