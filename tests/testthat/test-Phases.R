test_that("AddPhases works", {
  # add with an equation
  expect_no_error(AddPhases(ThisProblem = carbonate_system_problem,
                            PhaseEquation = "CO2(g) = 1 * CO3 + 2 * H",
                            PhaseLogK = -1.5,
                            PhaseDeltaH = 0,
                            PhaseTempKelvin = 298,
                            PhaseMoles = 10^-3.5))

  # ERROR on add without equation or stoich
  expect_error(AddPhases(ThisProblem = carbonate_system_problem,
                            PhaseLogK = -1.5,
                            PhaseDeltaH = 0,
                            PhaseTempKelvin = 298,
                            PhaseMoles = 10^-3.5))

  # add with stoich lists
  expect_no_error(AddPhases(ThisProblem = carbonate_system_problem,
                            PhaseName = "CO2(g)",
                            PhaseCompNames = list(c("CO3","H")),
                            PhaseCompStoichs = list(c(1, 2)),
                            PhaseLogK = -1.5,
                            PhaseDeltaH = 0,
                            PhaseTempKelvin = 298,
                            PhaseMoles = 10^-3.5))

  # add with agreeing equation, stoich lists, and stoich matrix
  expect_no_error(AddPhases(ThisProblem = carbonate_system_problem,
                            PhaseEquation = "CO2(g) = 1 * CO3 + 2 * H",
                            PhaseName = "CO2(g)",
                            PhaseCompNames = list(c("CO3","H")),
                            PhaseCompStoichs = list(c(1, 2)),
                            PhaseStoich = matrix(c(1L, 2L, 0L), nrow = 1, ncol = 3, dimnames = list("CO2(g)",c("CO3","H","OH"))),
                            PhaseLogK = -1.5,
                            PhaseDeltaH = 0,
                            PhaseTempKelvin = 298,
                            PhaseMoles = 10^-3.5))

  # ERROR on disagreeing equation, stoich lists, and stoich matrix
  expect_error(AddPhases(ThisProblem = carbonate_system_problem,
                            PhaseEquation = "CO2(g) = 1 * CO3 + 2 * H",
                            PhaseName = "CO2(g)",
                            PhaseCompNames = list(c("CO3","H")),
                            PhaseCompStoichs = list(c(1, 2)),
                            PhaseStoich = matrix(c(1L, 2L, 1L), nrow = 1, ncol = 3,
                                                 dimnames = list("CO2(g)",c("CO3","H","OH"))),
                            PhaseLogK = -1.5,
                            PhaseDeltaH = 0,
                            PhaseTempKelvin = 298,
                            PhaseMoles = 10^-3.5))
})
test_that("RemovePhases works", {
  myproblem = carbonate_system_problem
  myproblem = AddPhases(ThisProblem = carbonate_system_problem,
                        PhaseEquation = "CO2(g) = 1 * CO3 + 2 * H",
                        PhaseLogK = -1.5,
                        PhaseDeltaH = 0,
                        PhaseTempKelvin = 298,
                        PhaseMoles = 10^-3.5)
  expect_no_error(RemovePhases(ThisProblem = myproblem,
                               PhasesToRemove = "CO2(g)"))
  expect_no_error(RemovePhases(ThisProblem = myproblem,
                               PhasesToRemove = 1))
  expect_error(RemovePhases(ThisProblem = carbonate_system_problem,
                            PhasesToRemove = 1))
})
