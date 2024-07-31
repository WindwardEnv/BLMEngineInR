test_that("AddComponents works", {

  # nothing wrong
  expect_no_error(AddComponents(ThisProblem = carbonate_system_problem,
                                CompName = "test",
                                CompCharge = 0L,
                                CompMCName = "Water",
                                CompType = "MassBal",
                                CompActCorr = "None",
                                CompSiteDens = 1.0))

  # Invalid ThisProblem
  expect_error(AddComponents(ThisProblem = list(),
                             CompName = "test",
                             CompCharge = 0L,
                             CompMCName = "Water",
                             CompType = "MassBal",
                             CompActCorr = "None",
                             CompSiteDens = 1.0))

  # name already exists
  expect_error(AddComponents(ThisProblem = carbonate_system_problem,
                             CompName = "H",
                             CompCharge = 0L,
                             CompMCName = "Water",
                             CompType = "MassBal",
                             CompActCorr = "None",
                             CompSiteDens = 1.0))

  # Incorrect types
  expect_no_error(AddComponents(ThisProblem = carbonate_system_problem,
                                CompName = "test",
                                CompCharge = 0.5,#this will coerce to 0
                                CompMCName = "Water",
                                CompType = "MassBal",
                                CompActCorr = "None",
                                CompSiteDens = 1.0))
  expect_no_error(AddComponents(ThisProblem = carbonate_system_problem,
                                CompName = 123,#nothing technically wrong with this...coerces to a character
                                CompCharge = 0L,
                                CompMCName = "Water",
                                CompType = "MassBal",
                                CompActCorr = "None",
                                CompSiteDens = 1.0))

  # MC doesn't exist
  expect_error(AddComponents(ThisProblem = carbonate_system_problem,
                             CompName = "test",
                             CompCharge = 0L,
                             CompMCName = "Junk",#junk name
                             CompType = "MassBal",
                             CompActCorr = "None"))
  expect_error(AddComponents(ThisProblem = carbonate_system_problem,
                             CompName = "test",
                             CompCharge = 0L,
                             CompMCR = 999,#invalid index
                             CompType = "MassBal",
                             CompActCorr = "None"))

  # Invalid CompType
  expect_error(AddComponents(ThisProblem = carbonate_system_problem,
                             CompName = "test",
                             CompCharge = 0L,
                             CompMCName = "Water",
                             CompType = "Junk",
                             CompActCorr = "None"))

  # Invalid CompActCorr
  expect_error(AddComponents(ThisProblem = carbonate_system_problem,
                             CompName = "test",
                             CompCharge = 0L,
                             CompMCName = "Water",
                             CompType = "MassBal",
                             CompActCorr = "Junk",
                             CompSiteDens = 1.0))


})

test_that("AddInComp works", {

  # nothing wrong
  expect_no_error(AddInComps(ThisProblem = carbonate_system_problem,
                             InCompName = "test",
                             InCompCharge = 0L,
                             InCompMCName = "Water",
                             InCompType = "MassBal",
                             InCompActCorr = "None"))

  # name already exists - a different component
  expect_error(AddInComps(ThisProblem = carbonate_system_problem,
                          InCompName = "H",
                          InCompCharge = 0L,
                          InCompMCName = "Water",
                          InCompType = "MassBal",
                          InCompActCorr = "None"))

  # name already exists - another InComp
  expect_error(AddInComps(ThisProblem = carbonate_system_problem,
                          InCompName = "CO3",
                          InCompCharge = 0L,
                          InCompMCName = "Water",
                          InCompType = "MassBal",
                          InCompActCorr = "None"))

  # name already exists - species
  expect_error(AddInComps(ThisProblem = carbonate_system_problem,
                          InCompName = "OH",
                          InCompCharge = 0L,
                          InCompMCName = "Water",
                          InCompType = "MassBal",
                          InCompActCorr = "None"))

  # name already exists - neither component nor species
  expect_no_error(AddInComps(ThisProblem = carbonate_system_problem,
                             InCompName = "Water",
                             InCompCharge = 0L,
                             InCompMCName = "Water",
                             InCompType = "MassBal",
                             InCompActCorr = "None"))
  # This is silly and should probably be avoided, but since the component and MC
  # names should never really intermingle with each other, it's okay to do this,
  # and might be desired (e.g., for biotic ligand mass compartment and
  # components.)

})

test_that("AddDefComp works", {

  # nothing wrong
  expect_no_error(AddDefComps(ThisProblem = carbonate_system_problem,
                              DefCompName = "test",
                              DefCompFromNum = 1.0,
                              DefCompCharge = 0L,
                              DefCompMCName = "Water",
                              DefCompType = "MassBal",
                              DefCompActCorr = "None",
                              DefCompSiteDens = 1.0))
  expect_no_error(AddDefComps(ThisProblem = carbonate_system_problem,
                              DefCompName = "test",
                              DefCompFromVar = "CO3",
                              DefCompCharge = 0L,
                              DefCompMCName = "Water",
                              DefCompType = "MassBal",
                              DefCompActCorr = "None",
                              DefCompSiteDens = 1.0))
  expect_no_error(AddDefComps(ThisProblem = carbonate_system_problem,
                              DefCompName = "test",
                              DefCompFromVar = "Temp",
                              DefCompCharge = 0L,
                              DefCompMCName = "Water",
                              DefCompType = "MassBal",
                              DefCompActCorr = "None",
                              DefCompSiteDens = 1.0))

  # No DefCompFromNum or DefCompFromVar specified
  expect_error(AddDefComps(ThisProblem = carbonate_system_problem,
                           DefCompName = "test",
                           DefCompCharge = 0L,
                           DefCompMCName = "Water",
                           DefCompType = "MassBal",
                           DefCompActCorr = "None",
                           DefCompSiteDens = 1.0))

  # InValid DefCompFromVar
  expect_error(AddDefComps(ThisProblem = carbonate_system_problem,
                           DefCompName = "test",
                           DefCompFromVar = "Junk",
                           DefCompCharge = 0L,
                           DefCompMCName = "Water",
                           DefCompType = "MassBal",
                           DefCompActCorr = "None",
                           DefCompSiteDens = 1.0))

  # name already exists - a different component
  expect_error(AddDefComps(ThisProblem = carbonate_system_problem,
                           DefCompName = "CO3",
                           DefCompFromNum = 1.0,
                           DefCompCharge = 0L,
                           DefCompMCName = "Water",
                           DefCompType = "MassBal",
                           DefCompActCorr = "None",
                           DefCompSiteDens = 1.0))

  # name already exists - another DefComp
  expect_error(AddDefComps(ThisProblem = carbonate_system_problem,
                           DefCompName = "H",
                           DefCompFromNum = 1.0,
                           DefCompCharge = 0L,
                           DefCompMCName = "Water",
                           DefCompType = "MassBal",
                           DefCompActCorr = "None",
                           DefCompSiteDens = 1.0))

  # name already exists - species
  expect_error(AddDefComps(ThisProblem = carbonate_system_problem,
                           DefCompName = "OH",
                           DefCompFromNum = 1.0,
                           DefCompCharge = 0L,
                           DefCompMCName = "Water",
                           DefCompType = "MassBal",
                           DefCompActCorr = "None",
                           DefCompSiteDens = 1.0))

  # name already exists - neither component nor species
  expect_no_error(AddDefComps(ThisProblem = carbonate_system_problem,
                              DefCompName = "Water",
                              DefCompFromNum = 1.0,
                              DefCompCharge = 0L,
                              DefCompMCName = "Water",
                              DefCompType = "MassBal",
                              DefCompActCorr = "None",
                              DefCompSiteDens = 1.0))
  # This is silly and should probably be avoided, but since the component and MC
  # names should never really intermingle with each other, it's okay to do this,
  # and might be desired (e.g., for biotic ligand mass compartment and
  # components.)

  # check content
  expect_equal(AddDefComps(ThisProblem = carbonate_system_problem,
                           DefCompName = "test",
                           DefCompFromNum = 1.0,
                           DefCompCharge = 0L,
                           DefCompMCName = "Water",
                           DefCompType = "MassBal",
                           DefCompActCorr = "None",
                           DefCompSiteDens = 1.0)$N[
                               c("Mass","InLab","InVar","InMass","InComp",
                                 "InDefComp","InSpec","DefComp","Comp","Spec")],
               c(Mass = 1L, InLab = 1L, InVar = 2L, InMass = 1L, InComp = 1L,
                 InDefComp = 1L, InSpec = 3L, DefComp = 2L, Comp = 3L,
                 Spec = 6L))

})

test_that("RemoveInComps works", {

  # All good
  expect_no_error(RemoveInComps(ThisProblem = carbonate_system_problem,
                                InCompToRemove = 1))
  expect_no_error(RemoveInComps(ThisProblem = carbonate_system_problem,
                                InCompToRemove = "CO3"))

  # Try to remove non-input component
  expect_error(RemoveInComps(ThisProblem = carbonate_system_problem,
                                InCompToRemove = "H"))

  # check content
  expect_equal(RemoveInComps(ThisProblem = carbonate_system_problem,
                             InCompToRemove = "CO3")$N[
                               c("Mass","InLab","InVar","InMass","InComp",
                                 "InDefComp","InSpec","DefComp","Comp","Spec")],
               c(Mass = 1L, InLab = 1L, InVar = 2L, InMass = 1L, InComp = 0L,
                 InDefComp = 0L, InSpec = 1L, DefComp = 1L, Comp = 1L,
                 Spec = 2L))

})
