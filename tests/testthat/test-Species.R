test_that("AddSpecies works", {

  # all's well
  expect_no_error(AddSpecies(ThisProblem = carbonate_system_problem,
                             SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
                             SpecName = "H3CO3",
                             SpecCompNames = list(c("H", "CO3")),
                             SpecCompStoichs = list(c(3, 1)),
                             SpecStoich = matrix(c(3, 0, 1), nrow = 1, ncol = 3, dimnames = list("H3CO3", c("H","OH","CO3"))),
                             SpecMCName = "Water",
                             SpecActCorr = "Debye",
                             SpecLogK = 1.23,
                             SpecDeltaH = -1234,
                             SpecTempKelvin = 298))
  expect_no_error(AddSpecies(ThisProblem = carbonate_system_problem,
                             SpecName = "H3CO3",
                             SpecCompNames = list(c("H", "CO3")),
                             SpecCompStoichs = list(c(3, 1)),
                             SpecMCName = "Water",
                             SpecActCorr = "Debye",
                             SpecLogK = 1.23,
                             SpecDeltaH = -1234,
                             SpecTempKelvin = 298))
  expect_no_error(AddSpecies(ThisProblem = carbonate_system_problem,
                             SpecName = "H3CO3",
                             SpecStoich = matrix(c(3, 0, 1), nrow = 1, ncol = 3, dimnames = list("H3CO3", c("H","OH","CO3"))),
                             SpecMCName = "Water",
                             SpecActCorr = "Debye",
                             SpecLogK = 1.23,
                             SpecDeltaH = -1234,
                             SpecTempKelvin = 298))
  expect_no_error(AddSpecies(ThisProblem = carbonate_system_problem,
                             SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
                             SpecMCName = "Water",
                             SpecActCorr = "Debye",
                             SpecLogK = 1.23,
                             SpecDeltaH = -1234,
                             SpecTempKelvin = 298))
  expect_no_error(AddSpecies(ThisProblem = carbonate_system_problem,
                             SpecName = "H3CO3",
                             SpecEquation = "= 3 * H + 1 * CO3",
                             SpecMCName = "Water",
                             SpecActCorr = "Debye",
                             SpecLogK = 1.23,
                             SpecDeltaH = -1234,
                             SpecTempKelvin = 298))
  expect_no_error(AddSpecies(ThisProblem = carbonate_system_problem,
                             SpecName = "H3CO3",
                             SpecEquation = "3 * H + 1 * CO3",
                             SpecMCName = "Water",
                             SpecActCorr = "Debye",
                             SpecLogK = 1.23,
                             SpecDeltaH = -1234,
                             SpecTempKelvin = 298))
  expect_equal(AddSpecies(ThisProblem = carbonate_system_problem,
                             SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
                             SpecMCName = "Water",
                             SpecActCorr = "Debye",
                             SpecLogK = 1.23,
                             SpecDeltaH = -1234,
                             SpecTempKelvin = 298),
               AddSpecies(ThisProblem = carbonate_system_problem,
                          SpecName = "H3CO3",
                          SpecEquation = "= 3 * H + 1 * CO3",
                          SpecMCName = "Water",
                          SpecActCorr = "Debye",
                          SpecLogK = 1.23,
                          SpecDeltaH = -1234,
                          SpecTempKelvin = 298))
  expect_equal(AddSpecies(ThisProblem = carbonate_system_problem,
                          SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
                          SpecMCName = "Water",
                          SpecActCorr = "Debye",
                          SpecLogK = 1.23,
                          SpecDeltaH = -1234,
                          SpecTempKelvin = 298),
               AddSpecies(ThisProblem = carbonate_system_problem,
                          SpecName = "H3CO3",
                          SpecEquation = "3 * H + 1 * CO3",
                          SpecMCName = "Water",
                          SpecActCorr = "Debye",
                          SpecLogK = 1.23,
                          SpecDeltaH = -1234,
                          SpecTempKelvin = 298))
  expect_snapshot(AddSpecies(ThisProblem = carbonate_system_problem,
                             SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
                             SpecMCName = "Water",
                             SpecActCorr = "Debye",
                             SpecLogK = 1.23,
                             SpecDeltaH = -1234,
                             SpecTempKelvin = 298))

  # species exists
  expect_error(AddSpecies(ThisProblem = carbonate_system_problem,
                          SpecEquation = "H2CO3 = 2 * H + 1 * CO3",
                          SpecMCName = "Water",
                          SpecActCorr = "Debye",
                          SpecLogK = 1.23,
                          SpecDeltaH = -1234,
                          SpecTempKelvin = 298),
               regexp = "already exist")
  expect_error(AddSpecies(ThisProblem = carbonate_system_problem,
                          SpecEquation = "CO3 = 2 * H + 1 * CO3",
                          SpecMCName = "Water",
                          SpecActCorr = "Debye",
                          SpecLogK = 1.23,
                          SpecDeltaH = -1234,
                          SpecTempKelvin = 298),
               regexp = "already exist")

  # MC mismatch
  expect_error(AddSpecies(ThisProblem = Cu_full_organic_problem,
             SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
             SpecMCName = "Water",
             SpecActCorr = "Debye",
             SpecLogK = 1.23,
             SpecDeltaH = -1234,
             SpecTempKelvin = 298,
             SpecMCR = 2),
             regexp = "SpecMCName does not match SpecMCR")

  # MC does not exist
  expect_error(AddSpecies(ThisProblem = Cu_full_organic_problem,
             SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
             SpecMCName = "Junk",
             SpecActCorr = "Debye",
             SpecLogK = 1.23,
             SpecDeltaH = -1234,
             SpecTempKelvin = 298),
             regexp = "does not exist")

  # NA's
  expect_error(AddSpecies(ThisProblem = carbonate_system_problem,
                             SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
                             SpecMCName = "Water",
                             SpecActCorr = NA,
                             SpecLogK = 1.23,
                             SpecDeltaH = -1234,
                             SpecTempKelvin = 298),
               regexp = "NA arguments not allowed")
  expect_error(AddSpecies(ThisProblem = carbonate_system_problem,
                             SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
                             SpecMCName = "Water",
                             SpecActCorr = "Debye",
                             SpecLogK = NA,
                             SpecDeltaH = -1234,
                             SpecTempKelvin = 298),
               regexp = "NA arguments not allowed")
  expect_error(AddSpecies(ThisProblem = carbonate_system_problem,
                             SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
                             SpecMCName = "Water",
                             SpecActCorr = "Debye",
                             SpecLogK = 1.23,
                             SpecDeltaH = NA,
                             SpecTempKelvin = 298),
               regexp = "NA arguments not allowed")
  expect_error(AddSpecies(ThisProblem = carbonate_system_problem,
                             SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
                             SpecMCName = "Water",
                             SpecActCorr = "Debye",
                             SpecLogK = 1.23,
                             SpecDeltaH = -1234,
                             SpecTempKelvin = NA),
               regexp = "NA arguments not allowed")

  # Invalid SpecActCorr
  expect_error(AddSpecies(ThisProblem = carbonate_system_problem,
                          SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
                          SpecMCName = "Water",
                          SpecActCorr = "Junk",
                          SpecLogK = 1.23,
                          SpecDeltaH = -1234,
                          SpecTempKelvin = 298),
               regexp = "SpecActCorr values must be")

  # Invalid Type
  expect_error(AddSpecies(ThisProblem = carbonate_system_problem,
                          SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
                          SpecMCName = "Water",
                          SpecActCorr = "Debye",
                          SpecType = "Junk",
                          SpecLogK = 1.23,
                          SpecDeltaH = -1234,
                          SpecTempKelvin = 298),
               regexp = "SpecType values must be")

  # weird equations get fixed
  expect_equal(tail(AddSpecies(ThisProblem = carbonate_system_problem,
                               SpecEquation = "H3CO3 = 1 * H + 2 * H + 1 * CO3",
                               SpecMCName = "Water",
                               SpecActCorr = "Debye",
                               SpecLogK = 1.23,
                               SpecDeltaH = -1234,
                               SpecTempKelvin = 298)$Spec$Equation,1),
               "H3CO3 = 3 * H + 1 * CO3")
  expect_equal(tail(AddSpecies(ThisProblem = carbonate_system_problem,
                               SpecEquation = "H3CO3 = 1 * CO3 + 3 * H ",
                               SpecMCName = "Water",
                               SpecActCorr = "Debye",
                               SpecLogK = 1.23,
                               SpecDeltaH = -1234,
                               SpecTempKelvin = 298)$Spec$Equation,1),
               "H3CO3 = 3 * H + 1 * CO3")

  # missing stoichiometry and name
  expect_error(AddSpecies(ThisProblem = carbonate_system_problem,
                          SpecMCName = "Water",
                          SpecActCorr = "Debye",
                          SpecLogK = 1.23,
                          SpecDeltaH = -1234,
                          SpecTempKelvin = 298),
               regexp = "Must specify either")

  # missing stoichiometry
  expect_error(AddSpecies(ThisProblem = carbonate_system_problem,
                          SpecName = "H3CO3",
                          SpecMCName = "Water",
                          SpecActCorr = "Debye",
                          SpecLogK = 1.23,
                          SpecDeltaH = -1234,
                          SpecTempKelvin = 298),
               regexp = "Must specify either")

  # missing name
  expect_error(AddSpecies(ThisProblem = carbonate_system_problem,
             SpecCompNames = list(c("H", "CO3")),
             SpecCompStoichs = list(c(3, 1)),
             SpecMCName = "Water",
             SpecActCorr = "Debye",
             SpecLogK = 1.23,
             SpecDeltaH = -1234,
             SpecTempKelvin = 298),
             regexp = "Must specify either")

  # mismatched arguments
  expect_error(AddSpecies(ThisProblem = carbonate_system_problem,
                          SpecEquation = "H3CO3 = 999 * H + 1 * CO3",
                          SpecName = "H3CO3",
                          SpecCompNames = list(c("H", "CO3")),
                          SpecCompStoichs = list(c(3, 1)),
                          SpecStoich = matrix(c(3, 0, 1), nrow = 1, ncol = 3, dimnames = list("HCO3", c("H", "OH","CO3"))),
                          SpecMCName = "Water",
                          SpecActCorr = "Debye",
                          SpecLogK = 1.23,
                          SpecDeltaH = -1234,
                          SpecTempKelvin = 298))
  expect_error(AddSpecies(ThisProblem = carbonate_system_problem,
                             SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
                             SpecName = "H3CO3",
                             SpecCompNames = list(c("H", "CO3")),
                             SpecCompStoichs = list(c(999, 1)),
                             SpecStoich = matrix(c(3, 0, 1), nrow = 1, ncol = 3, dimnames = list("HCO3", c("H", "OH","CO3"))),
                             SpecMCName = "Water",
                             SpecActCorr = "Debye",
                             SpecLogK = 1.23,
                             SpecDeltaH = -1234,
                             SpecTempKelvin = 298))
  expect_error(AddSpecies(ThisProblem = carbonate_system_problem,
                             SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
                             SpecName = "H3CO3",
                             SpecCompNames = list(c("H", "CO3")),
                             SpecCompStoichs = list(c(3, 1)),
                             SpecStoich = matrix(c(999, 0, 1), nrow = 1, ncol = 3, dimnames = list("HCO3", c("H","OH","CO3"))),
                             SpecMCName = "Water",
                             SpecActCorr = "Debye",
                             SpecLogK = 1.23,
                             SpecDeltaH = -1234,
                             SpecTempKelvin = 298))
  expect_error(AddSpecies(ThisProblem = carbonate_system_problem,
                             SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
                             SpecName = "H999CO3",
                             SpecCompNames = list(c("H", "CO3")),
                             SpecCompStoichs = list(c(3, 1)),
                             SpecStoich = matrix(c(3, 0, 1), nrow = 1, ncol = 3, dimnames = list("HCO3", c("H","OH","CO3"))),
                             SpecMCName = "Water",
                             SpecActCorr = "Debye",
                             SpecLogK = 1.23,
                             SpecDeltaH = -1234,
                             SpecTempKelvin = 298))

})
test_that("RemoveSpecies works", {

  expect_no_error(RemoveSpecies(ThisProblem = carbonate_system_problem,
                                SpeciesToRemove = "HCO3"))
  expect_no_error(RemoveSpecies(ThisProblem = carbonate_system_problem,
                  SpeciesToRemove = 3))
  expect_error(RemoveSpecies(ThisProblem = carbonate_system_problem,
                  SpeciesToRemove = "Junk"),
               regexp = "does not exist")
  expect_error(RemoveSpecies(ThisProblem = carbonate_system_problem,
                  SpeciesToRemove = 999),
               regexp = "trying to remove the #")
  expect_error(RemoveSpecies(ThisProblem = carbonate_system_problem,
                  SpeciesToRemove = -999),
               regexp = "Invalid index")
  expect_snapshot(RemoveSpecies(ThisProblem = carbonate_system_problem,
                                SpeciesToRemove = "HCO3"))

})
