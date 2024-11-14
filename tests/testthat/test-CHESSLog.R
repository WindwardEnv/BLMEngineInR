test_that("CHESSLog works", {
  expect_no_error(CHESSLog(ThisProblem = carbonate_system_problem, LogFilename = withr::local_tempfile()))
  expect_true(CHESSLog(ThisProblem = carbonate_system_problem, LogFilename = withr::local_tempfile()))
  expect_no_error(CHESSLog(ThisProblem = Cu_full_organic_problem, LogFilename = withr::local_tempfile()))
})
