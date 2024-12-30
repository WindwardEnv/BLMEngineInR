test_that("ExpandWHAM works", {

  myproblem = Cu_full_organic_problem
  myproblem$WHAM$File = basename(myproblem$WHAM$File)
  myproblem_noDOC = RemoveSpecialDefs(ThisProblem = myproblem,
                                      SpecialDefToRemove = "WHAM")
  expect_no_error(ExpandWHAM(ThisProblem = myproblem_noDOC, ThisWHAM = DefineWHAM(WHAMVer = "V")))

  skip_on_cran()
  # expect_no_error(ExpandWHAM(ThisProblem = myproblem_noDOC))

  myproblem_WHAMV = ExpandWHAM(ThisProblem = myproblem_noDOC,
                               ThisWHAM = DefineWHAM(WHAMVer = "V"))
  myproblem_WHAMV$WHAM$File = basename(myproblem_WHAMV$WHAM$File)
  expect_equal(myproblem[names(myproblem) != "ParamFile"],
               myproblem_WHAMV[names(myproblem_WHAMV) != "ParamFile"])
  expect_equal(myproblem[c("N","Mass","DefComp","Comp","Spec","SpecStoich")],
               myproblem_WHAMV[c("N","Mass","DefComp","Comp","Spec","SpecStoich")])

  # WHAM V
  myproblem_humics = AddInComps(
    ThisProblem = AddInVars(
      ThisProblem = carbonate_system_problem,
      InVarName = "Humics",
      InVarMCName = "Water",
      InVarType = "WHAM-HA"
    ),
    InCompName = "Ca",
    InCompCharge = 2,
    InCompMCName = "Water",
    InCompType = "MassBal",
    InCompActCorr = "Debye"
  )
  myproblem_HA_V = ExpandWHAM(ThisProblem = myproblem_humics, ThisWHAM = DefineWHAM(WHAMVer = "V"))
  expect_equal(myproblem_HA_V$N,
               c(Mass = 2L, InLab = 1L, InVar = 3L, InMass = 1L, InComp = 2L,
                 InDefComp = 2L, InSpec = 2L, DefComp = 23L, Comp = 25L,
                 Spec = 96L, Phase = 0L, BL = 0L, Metal = 0L, BLMetal = 0L,
                 CAT = 0L))

  # WHAM VII
  myproblem_HA_VII = ExpandWHAM(ThisProblem = myproblem_humics, ThisWHAM = DefineWHAM(WHAMVer = "VII"))
  expect_equal(myproblem_HA_VII$N,
               c(Mass = 2L, InLab = 1L, InVar = 3L, InMass = 1L, InComp = 2L,
                 InDefComp = 2L, InSpec = 2L, DefComp = 25L, Comp = 27L,
                 Spec = 138L, Phase = 0L, BL = 0L, Metal = 0L, BLMetal = 0L,
                 CAT = 0L))

})
