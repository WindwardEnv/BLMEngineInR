test_that("AddSpecialDefs works", {
  myproblem = AddInComps(ThisProblem = carbonate_system_problem,
                         InCompName = "Cu",
                         InCompCharge = 2,
                         InCompMCName = "Water",
                         InCompType = "MassBal",
                         InCompActCorr = "Debye")

  expect_no_error(AddSpecialDefs(ThisProblem = myproblem,
                                 Value = "Cu",
                                 SpecialDef = "Metal"))
  expect_error(AddSpecialDefs(ThisProblem = carbonate_system_problem,
                                 Value = "Cu",
                                 SpecialDef = "Metal"))

  myproblem = AddInVars(ThisProblem = myproblem,
                        InVarName = "Humics",
                        InVarMCName = "Water",
                        InVarType = "WHAM-HA")
  expect_no_error(AddSpecialDefs(ThisProblem = myproblem,
                                 Value = "V",
                                 SpecialDef = "WHAM"))
  expect_error(AddSpecialDefs(ThisProblem = myproblem,
                                 Value = "Junk",
                                 SpecialDef = "Junk"))
})
test_that("RemoveSpecialDefs works", {
  expect_no_error(RemoveSpecialDefs(ThisProblem = Cu_full_organic_problem,
                                    SpecialDefToRemove = "BL-Metal",
                                    Index = 2))
  expect_error(RemoveSpecialDefs(ThisProblem = Cu_full_organic_problem,
                                    SpecialDefToRemove = "Junk"))
})
