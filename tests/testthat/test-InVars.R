test_that("AddInVars works", {
  expect_no_error(
    AddInVars(
      ThisProblem = carbonate_system_problem,
      InVarName = "DOC",
      InVarMCName = "Water",
      InVarType = "WHAM-HAFA"
    )
  )
  expect_error(
    AddInVars(
      ThisProblem = carbonate_system_problem,
      InVarName = "Junk",
      InVarMCName = "Water",
      InVarType = "Junk"
    )
  )
  expect_error(
    AddInVars(
      ThisProblem = carbonate_system_problem,
      InVarName = "DOC",
      InVarMCName = "Junk",
      InVarType = "WHAM-HAFA"
    )
  )
})
test_that("RemoveInVars works", {
  expect_no_error(
    RemoveInVars(
      ThisProblem = carbonate_system_problem,
      InVarToRemove = 1
    )
  )
  expect_no_error(
    RemoveInVars(
      ThisProblem = carbonate_system_problem,
      InVarToRemove = "pH"
    )
  )
})
