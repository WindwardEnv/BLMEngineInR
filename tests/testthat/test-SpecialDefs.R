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
                                 SpecialDef = "Metal"),
               regexp = "Unknown component")
  expect_error(AddSpecialDefs(ThisProblem = carbonate_system_problem,
                                 Value = "Cu",
                                 SpecialDef = "BL"),
               regexp = "Unknown component")
  expect_error(AddSpecialDefs(ThisProblem = carbonate_system_problem,
                                 Value = "BL1-Cu",
                                 SpecialDef = "BLMetal"),
               regexp = "Unknown species")
  expect_no_error(AddSpecialDefs(ThisProblem = Cu_full_organic_problem,
                                 Value = "BL1-Ca",
                                 SpecialDef = "BLMetal"))

  expect_error(AddSpecialDefs(ThisProblem = myproblem,
                                 Value = NA,
                                 SpecialDef = "Metal"),
               regexp = "NA inputs not allowed")
  expect_error(AddSpecialDefs(ThisProblem = myproblem,
                                 Value = "Cu",
                                 SpecialDef = NA),
               regexp = "NA inputs not allowed")

  myproblem = AddInVars(ThisProblem = myproblem,
                        InVarName = "Humics",
                        InVarMCName = "Water",
                        InVarType = "WHAM-HA")
  expect_no_error(AddSpecialDefs(ThisProblem = myproblem,
                                 Value = "V",
                                 SpecialDef = "WHAM"))
  expect_error(AddSpecialDefs(ThisProblem = carbonate_system_problem,
                              Value = "V",
                              SpecialDef = "WHAM"),
               regexp = "without a WHAM input variable")
  expect_error(AddSpecialDefs(ThisProblem = Cu_full_organic_problem,
                              Value = "VII",
                              SpecialDef = "WHAM"),
               regexp = "Only one WHAM version")
  expect_error(AddSpecialDefs(ThisProblem = myproblem,
                              Value = "Junk",
                              SpecialDef = "Junk"),
               regexp = "should be one of")
})
test_that("RemoveSpecialDefs works", {
  expect_no_error(RemoveSpecialDefs(ThisProblem = Cu_full_organic_problem,
                                    SpecialDefToRemove = "BL-Metal",
                                    Index = 2))
  expect_no_error(RemoveSpecialDefs(ThisProblem = Cu_full_organic_problem,
                                    SpecialDefToRemove = "WHAM"))
  expect_equal(RemoveSpecialDefs(ThisProblem = Cu_full_organic_problem,
                                 SpecialDefToRemove = "WHAM")$N,
               c(Mass=2L, InLab = 3L, InVar = 4L, InMass = 2L, InComp = 8L,
                 InDefComp = 2L, InSpec = 22L, DefComp = 2L, Comp = 10L,
                 Spec = 32L, Phase = 0L, BL = 1L, Metal = 1L, BLMetal = 2L,
                 CAT = 2L))
  expect_error(RemoveSpecialDefs(ThisProblem = Cu_full_organic_problem,
                                    SpecialDefToRemove = "Junk"),
               regexp = "should be one of")
})
