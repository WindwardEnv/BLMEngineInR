test_that("ExpandWHAM works", {

  myproblem = Cu_full_organic_problem
  myproblem$WHAM$WdatFile = gsub(dirname(myproblem$WHAM$WdatFile), "", myproblem$WHAM$WdatFile)
  myproblem_noDOC = RemoveSpecialDefs(ThisProblem = myproblem,
                                      SpecialDefToRemove = "WHAM")
  expect_no_error(ExpandWHAM(ThisProblem = myproblem_noDOC, WHAMVer = "V"))

  skip_on_cran()
  expect_no_error(ExpandWHAM(ThisProblem = myproblem_noDOC))

  myproblem_WHAMV = ExpandWHAM(ThisProblem = myproblem_noDOC, WHAMVer = "V")
  myproblem_WHAMV$WHAM$WdatFile = gsub(dirname(myproblem_WHAMV$WHAM$WdatFile), "", myproblem_WHAMV$WHAM$WdatFile)
  expect_equal(myproblem[names(myproblem) != "ParamFile"],
               myproblem_WHAMV[names(myproblem_WHAMV) != "ParamFile"])
  expect_equal(myproblem[c("N","Mass","DefComp","Comp","Spec","SpecStoich")],
               ExpandWHAM(ThisProblem = myproblem_noDOC,
                          WdatFile = system.file(file.path("extdata","WHAM","WHAM_V.wdat"),
                                                 package = "BLMEngineInR", mustWork = TRUE))[c("N","Mass","DefComp","Comp","Spec","SpecStoich")])

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
  myproblem_HA_V = ExpandWHAM(ThisProblem = myproblem_humics, WHAMVer = "V")
  expect_equal(myproblem_HA_V$N,
               c(Mass = 2L, InLab = 1L, InVar = 3L, InMass = 1L, InComp = 2L,
                 InDefComp = 1L, InSpec = 3L, DefComp = 22L, Comp = 24L,
                 Spec = 96L, Phase = 0L, BL = 0L, Metal = 0L, BLMetal = 0L,
                 CAT = 0L))

  # WHAM VII
  myproblem_HA_VII = ExpandWHAM(ThisProblem = myproblem_humics, WHAMVer = "VII")
  expect_equal(myproblem_HA_VII$N,
               c(Mass = 2L, InLab = 1L, InVar = 3L, InMass = 1L, InComp = 2L,
                 InDefComp = 1L, InSpec = 3L, DefComp = 24L, Comp = 26L,
                 Spec = 138L, Phase = 0L, BL = 0L, Metal = 0L, BLMetal = 0L,
                 CAT = 0L))

})
