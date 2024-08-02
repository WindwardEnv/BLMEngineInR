test_that("AddInLabs works", {
  expect_no_error(AddInLabs(ThisProblem = carbonate_system_problem, InLabName = "Test"))
})
test_that("RemoveInLabs works", {
  expect_no_error(RemoveInLabs(ThisProblem = carbonate_system_problem, InLabToRemove = 1))
})
