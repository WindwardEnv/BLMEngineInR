test_that("ReadInputsFromFile works", {

  myinputfile = withr::local_tempfile(
    lines = c(
      "1",
      "ID,    Temp,   pH,  CO3",
      "n/a,   deg C,  SU,  mol/L",
      "Test,  25,     7,   1e-04"
    ),
    fileext = ".dat4"
  )
  myproblem = carbonate_system_problem

  expect_no_error(GetData(InputFile = myinputfile, ThisProblem = myproblem))
  expect_no_error(GetData(InputFile = myinputfile,
                          NInLab = myproblem$N["InLab"],
                          InLabName = myproblem$InLabName,
                          NInVar = myproblem$N["InVar"],
                          InVarName = myproblem$InVar$Name,
                          InVarMCR = myproblem$InVar$MCR,
                          InVarType = myproblem$InVar$Type,
                          NInComp = myproblem$N["InComp"],
                          InCompName = myproblem$InCompName,
                          NComp = myproblem$N["Comp"],
                          CompName = myproblem$Comp$Name,
                          NDefComp = myproblem$N["DefComp"],
                          DefCompName = myproblem$DefComp$Name,
                          DefCompFromNum = myproblem$DefComp$FromNum,
                          DefCompFromVar = myproblem$DefComp$FromVar,
                          DefCompSiteDens = myproblem$DefComp$SiteDens))

})

test_that("MatchInputsToProblem works", {

  myproblem = Cu_full_organic_problem

  myinputsDF = data.frame(
    ObsNum = "1",
    ID = "1",
    ID2 = "test",
    Temp = 15,
    pH = 7,
    DOC = 1,
    HA = 10,
    AFA = 30,
    Cu = 1E-7,
    Ca = 1E-3,
    Mg = 1E-3,
    Na = 1E-4,
    K = 1E-5,
    SO4 = 1E-3,
    Cl = 1E-4,
    CO3 = 1E-3,
    Extra = "Junk"
  )

  expect_no_error(MatchInputsToProblem(DFInputs = myinputsDF, ThisProblem = myproblem))
  expect_error(MatchInputsToProblem(DFInputs = myinputsDF[, -1], ThisProblem = myproblem))
  expect_no_error(MatchInputsToProblem(NObs = nrow(myinputsDF),
                                       InLabObs = myinputsDF[, myproblem$InLabName, drop = FALSE],
                                       InVarObs = myinputsDF[, myproblem$InVar$Name, drop = FALSE],
                                       InCompObs = myinputsDF[, myproblem$InCompName, drop = FALSE],
                                       NInVar = myproblem$N["InVar"],
                                       InVarName = myproblem$InVar$Name,
                                       InVarMCR = myproblem$InVar$MCR,
                                       InVarType = myproblem$InVar$Type,
                                       NInComp = myproblem$N["InComp"],
                                       InCompName = myproblem$InCompName,
                                       NComp = myproblem$N["Comp"],
                                       CompName = myproblem$Comp$Name,
                                       NDefComp = myproblem$N["DefComp"],
                                       DefCompName = myproblem$DefComp$Name,
                                       DefCompFromNum = myproblem$DefComp$FromNum,
                                       DefCompFromVar = myproblem$DefComp$FromVar,
                                       DefCompSiteDens = myproblem$DefComp$SiteDens))
  expect_equal(MatchInputsToProblem(DFInputs = myinputsDF, ThisProblem = myproblem),
               MatchInputsToProblem(NObs = nrow(myinputsDF),
                                    InLabObs = myinputsDF[, myproblem$InLabName, drop = FALSE],
                                    InVarObs = myinputsDF[, myproblem$InVar$Name, drop = FALSE],
                                    InCompObs = myinputsDF[, myproblem$InCompName, drop = FALSE],
                                    NInVar = myproblem$N["InVar"],
                                    InVarName = myproblem$InVar$Name,
                                    InVarMCR = myproblem$InVar$MCR,
                                    InVarType = myproblem$InVar$Type,
                                    NInComp = myproblem$N["InComp"],
                                    InCompName = myproblem$InCompName,
                                    NComp = myproblem$N["Comp"],
                                    CompName = myproblem$Comp$Name,
                                    NDefComp = myproblem$N["DefComp"],
                                    DefCompName = myproblem$DefComp$Name,
                                    DefCompFromNum = myproblem$DefComp$FromNum,
                                    DefCompFromVar = myproblem$DefComp$FromVar,
                                    DefCompSiteDens = myproblem$DefComp$SiteDens))

  myproblem_DC = AddDefComps(ThisProblem = myproblem, DefCompName = c("Junk1", "Junk2"),
                             DefCompFromVar = c("Ca","Temp"), DefCompCharge = 3,
                             DefCompMCName = "Water", DefCompType = "FixedConc",
                             DefCompActCorr = "Debye", DefCompSiteDens = 2.0)
  expect_equal(MatchInputsToProblem(DFInputs = myinputsDF,
                                    ThisProblem = myproblem_DC)$TotConcObs[, c("Junk1", "Junk2"), drop = FALSE],
               matrix(as.numeric(myinputsDF[,c("Ca", "Temp")] * 2), nrow = 1,
                      dimnames = list(Obs = 1, Comp = c("Junk1", "Junk2"))))

  myproblem_AFA = AddInVars(ThisProblem = myproblem, InVarName = "AFA",
                            InVarMCName = "Water", InVarType = "PercAFA")
  expect_equal(MatchInputsToProblem(DFInputs = myinputsDF, ThisProblem = myproblem_AFA)$HumicSubstGramsPerLiterObs,
               c(1, 0.3) * MatchInputsToProblem(DFInputs = myinputsDF, ThisProblem = myproblem)$HumicSubstGramsPerLiterObs)

})

test_that("GetData works", {

  myproblem = carbonate_system_problem
  myinputfile = withr::local_tempfile(
    lines = c(
      "1",
      "ID,    Temp,   pH,  CO3",
      "n/a,   deg C,  SU,  mol/L",
      "Test,  25,     7,   1e-04"
    ),
    fileext = ".dat4"
  )

  expect_no_error(GetData(InputFile = myinputfile, ThisProblem = myproblem))

})

test_that("BlankInputList works", {

  myproblem = carbonate_system_problem
  expect_no_error(BlankInputList(ThisProblem = myproblem))
  expect_no_error(BlankInputList(ThisProblem = myproblem, NObs = 999L))
  expect_error(BlankInputList(ThisProblem = myproblem, NObs = -999L))
  expect_error(BlankInputList(list()))
  expect_true(all(c("NObs", "InLabObs","InVarObs", "InCompObs", "SysTempKelvinObs",
                    "pHObs", "TotConcObs", "HumicSubstGramsPerLiterObs") %in%
                    names(BlankInputList(ThisProblem = myproblem, NObs = 1L))))

})
