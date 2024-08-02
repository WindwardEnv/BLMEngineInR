test_that("AddSpecies works", {

  expect_no_error(AddSpecies(
    ThisProblem = carbonate_system_problem,
    SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
    SpecMCName = "Water",
    SpecActCorr = "Debye",
    SpecLogK = 1.23,
    SpecDeltaH = -1234,
    SpecTempKelvin = 298
  ))

  expect_snapshot(AddSpecies(
    ThisProblem = carbonate_system_problem,
    SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
    SpecMCName = "Water",
    SpecActCorr = "Debye",
    SpecLogK = 1.23,
    SpecDeltaH = -1234,
    SpecTempKelvin = 298
  ))

  expect_equal(tail(
    AddSpecies(
      ThisProblem = carbonate_system_problem,
      SpecEquation = "H3CO3 = 1 * H + 2 * H + 1 * CO3",
      SpecMCName = "Water",
      SpecActCorr = "Debye",
      SpecLogK = 1.23,
      SpecDeltaH = -1234,
      SpecTempKelvin = 298
    )$Spec$Equation,
    1
  ), "H3CO3 = 3 * H + 1 * CO3")

})
test_that("RemoveSpecies works", {

  expect_no_error(RemoveSpecies(ThisProblem = carbonate_system_problem,
                  SpeciesToRemove = "HCO3"))
  expect_no_error(RemoveSpecies(ThisProblem = carbonate_system_problem,
                  SpeciesToRemove = 3))
  expect_snapshot(RemoveSpecies(ThisProblem = carbonate_system_problem,
                                SpeciesToRemove = "HCO3"))

})
