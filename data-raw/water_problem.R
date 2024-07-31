water_MC_problem = AddInLabs(
  ThisProblem = AddMassCompartments(ThisProblem = BlankProblem(),
                                     MassName = "Water",
                                     MassAmt = 1.0,
                                     MassUnit = "L"),
  InLabName = "ID")

water_problem = AddSpecies(
  ThisProblem = AddInVars(ThisProblem = water_MC_problem,
                          InVarName = c("Temperature", "pH"),
                          InVarMCName = "Water",
                          InVarType = c("Temperature", "pH")),
  SpecEquation = "OH = -1 * H",
  SpecMCName = "Water",
  SpecActCorr = "Debye",
  SpecLogK = -13.997,
  SpecDeltaH = -55810,
  SpecTempKelvin = 298.15)

usethis::use_data(water_MC_problem, overwrite = TRUE)
usethis::use_data(water_problem, overwrite = TRUE)
