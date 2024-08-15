test_that("AddPhases works", {
  expect_no_error(AddPhases(ThisProblem = carbonate_system_problem,
                            PhaseEquation = "CO2(g) = 1 * CO3 + 2 * H",
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
