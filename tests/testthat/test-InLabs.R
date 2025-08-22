test_that("AddInLabs works", {

  # All's well
  expect_no_error(AddInLabs(ThisProblem = carbonate_system_problem,
                            InLabName = "Test"))

  # NA names not allowed
  expect_error(AddInLabs(ThisProblem = carbonate_system_problem,
                         InLabName = NA))

  # names that already exist not allowed
  expect_error(AddInLabs(ThisProblem = carbonate_system_problem,
                         InLabName = carbonate_system_problem$InLabName[1]))

})
test_that("RemoveInLabs works", {

  # All's well
  expect_no_error(RemoveInLabs(ThisProblem = carbonate_system_problem,
                               InLabToRemove = 1))

  # removing name that doesn't exist
  expect_error(RemoveInLabs(ThisProblem = carbonate_system_problem,
                            InLabToRemove = "NotThere"),
               regexp = "does not exist")

  # removing number that doesn't exist
  expect_error(RemoveInLabs(ThisProblem = carbonate_system_problem,
                            InLabToRemove = 999),
               regexp = "trying to remove")

  # bad index
  expect_error(RemoveInLabs(ThisProblem = carbonate_system_problem,
                            InLabToRemove = -1),
               regexp = "Invalid index")

})
