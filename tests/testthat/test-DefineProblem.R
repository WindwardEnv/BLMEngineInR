test_that("DefineProblem works", {

  mypfile = withr::local_tempfile(fileext = ".dat4")
  WriteParamFile(ThisProblem = Cu_full_organic_problem, ParamFile = mypfile)
  expect_no_error({myproblem = DefineProblem(ParamFile = mypfile)})

  compare_names = setdiff(names(Cu_full_organic_problem), c("ParamFile","WHAM"))
  compare_names_WHAM = setdiff(names(Cu_full_organic_problem$WHAM),
                               c("File", "Ver", "Notes"))

  expect_equal(myproblem[compare_names],
               Cu_full_organic_problem[compare_names])
  expect_equal(myproblem$WHAM[compare_names_WHAM],
               Cu_full_organic_problem$WHAM[compare_names_WHAM])

  myproblem$WHAM$File = basename(myproblem$WHAM$File)

  mytestproblem = Cu_full_organic_problem
  mytestproblem$WHAM$File = basename(mytestproblem$WHAM$File)

  expect_equal(myproblem[compare_names],
               mytestproblem[compare_names])
  expect_equal(myproblem$WHAM[compare_names_WHAM],
               mytestproblem$WHAM[compare_names_WHAM])

})



