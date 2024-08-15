test_that("AddInVars works", {

  # all's well
  expect_no_error(AddInVars(ThisProblem = carbonate_system_problem, InVarName = "DOC",
                            InVarMCName = "Water",InVarType = "WHAM-HAFA"))
  expect_no_error(AddInVars(ThisProblem = water_MC_problem, InVarName = "pH",
                            InVarMCName = "Water", InVarType = "pH"))

  # duplicate InVarName
  expect_error(AddInVars(ThisProblem = carbonate_system_problem,
                         InVarName = carbonate_system_problem$InVar$Name[1],
                         InVarMCName = "Water", InVarType = "WHAM-HAFA"),
               regexp = "already exist")

  # NA arguments
  expect_error(AddInVars(ThisProblem = carbonate_system_problem, InVarName = NA,
                         InVarMCName = "Water", InVarType = "WHAM-HAFA"),
               regexp = "NA arguments not allowed")
  expect_error(AddInVars(ThisProblem = carbonate_system_problem, InVarName = "DOC",
                         InVarMCName = NA,  InVarType = "WHAM-HAFA"),
               regexp = "NA arguments not allowed")
  expect_error(AddInVars(ThisProblem = carbonate_system_problem, InVarName = "DOC",
                         InVarMCName = "Water", InVarType = NA),
               regexp = "NA arguments not allowed")

  # bad InVarType
  expect_error(AddInVars(ThisProblem = carbonate_system_problem, InVarName = "Junk",
                         InVarMCName = "Water", InVarType = "Junk"),
               regexp = "Invalid InVarType")

  # bad MCName
  expect_error(AddInVars(ThisProblem = carbonate_system_problem, InVarName = "DOC",
                         InVarMCName = "Junk", InVarType = "WHAM-HAFA"),
               regexp = "Mass compartment")

})
test_that("RemoveInVars works", {
  expect_no_error(RemoveInVars(ThisProblem = carbonate_system_problem, InVarToRemove = 1))
  expect_no_error(RemoveInVars(ThisProblem = carbonate_system_problem, InVarToRemove = "pH"))
  expect_error(RemoveInVars(ThisProblem = carbonate_system_problem, InVarToRemove = "Junk"),
               regexp = "do not exist")
  expect_error(RemoveInVars(ThisProblem = carbonate_system_problem, InVarToRemove = -999),
               regexp = "Invalid index")
  expect_error(RemoveInVars(ThisProblem = carbonate_system_problem, InVarToRemove = 999),
               regexp = "trying to remove")
})
