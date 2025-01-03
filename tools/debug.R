rm(list = ls())
# devtools::clean_dll()
devtools::load_all()

# ParamFile = "inst/extdata/ParameterFiles/Cu_full_organic.dat4"
# InputFile = "inst/extdata/InputFiles/Cu_full_organic.blm4"
ParamFile = "inst/extdata/ParameterFiles/Zn_full_organic_S.dat4"
InputFile = "inst/extdata/InputFiles/Zn_full_organic_S.blm4"

ThisProblem = DefineProblem(ParamFile, WriteLog = TRUE)
tmp = ReadInputsFromFile(InputFile, ThisProblem)
InputDF = data.frame(tmp$InLabObs, tmp$InVarObs, tmp$InCompObs)
# AllInput = GetData(InputFile, ThisProblem = ThisProblem)

InputDF = data.frame(
  Temp = 15,
  pH = 7,
  DOC = 20,
  HA = 0.01,
  Cu_ug.L = 0.1,
  Zn_ug.L = 0.1,
  Hard = 10,
  S = 1E-20
)
InputDF$ObsNum = InputDF$ID = InputDF$ID2 = as.character("#13")
InputDF$Cu = InputDF$Cu_ug.L / (10^6 * MW["Cu"])
InputDF$Zn = InputDF$Zn_ug.L / (10^6 * MW["Zn"])
InputDF$Ca = (InputDF$Hard / (1000 * MW["CaCO3"])) / (1 + 1 / 0.7)
InputDF$Mg = InputDF$Ca / 0.7
InputDF$Na = InputDF$Ca / 0.31
InputDF$K = InputDF$Ca / 6.50
InputDF$SO4 = InputDF$Ca / 0.41 # CaSO4 and MgSO4 used in EPA water
InputDF$Cl = InputDF$K #KCl used as salt in EPA water
InputDF$CO3 = InputDF$Na #NaHCO3 used as salt in EPA water

AllInput = MatchInputsToProblem(DFInputs = InputDF,
                                ThisProblem = ThisProblem)


ResultsTable <- BLM(
  ThisProblem = ThisProblem,
  AllInput = AllInput,
  DoTox = TRUE,
  QuietFlag ="Debug"#,
  # ConvergenceCriteria = 0.0001,
  # MaxIter = 100L,
  # DodVidCj = TRUE,
  # DodVidCjDonnan = FALSE,
  # DodKidCj = FALSE,
  # DoGammai = TRUE,
  # DoJacDonnan = FALSE,
  # DoJacWHAM = TRUE,
  # DoWHAMSimpleAdjust = TRUE,
  # DoDonnanSimpleAdjust = TRUE
)
ResultsTable$Miscellaneous

ThisProblem_toxmode = AddCriticalValues(
  ThisProblem = RemoveMassCompartments(
    ThisProblem = AddSpecialDefs(
      ThisProblem = AddSpecies(
        ThisProblem = AddDefComps(
          ThisProblem = AddInComps(
            ThisProblem = AddInVars(
              ThisProblem = AddMassCompartments(
                ThisProblem = water_problem,
                MassName = "BL",
                MassAmt = 1,
                MassUnit = "kg wet"
              ),
              InVarName = "DOC",
              InVarMCName = "Water",
              InVarType = "WHAM-FA"),
            InCompName = c("Cu","Ca","SO4"),
            InCompCharge = c(2, 2, -2),
            InCompMCName = "Water",
            InCompType = "MassBal",
            InCompActCorr = "None"
          ),
          DefCompName = "BL1",
          DefCompFromNum = 1.78e-5,
          DefCompCharge = -1,
          DefCompMCName = "BL",
          DefCompType = "MassBal",
          DefCompActCorr = "None",
          DefCompSiteDens = 3e-5
        ),
        SpecEquation = c("CaSO4 = 1 * Ca + 1 * SO4",
                         "BL1-Cu = 1 * BL1 + 1 * Cu"),
        SpecMCName = c("Water", "BL"),
        SpecLogK = c(2.3, 7.4),
        SpecActCorr = "None",
        SpecDeltaH = 0,
        SpecTempKelvin = 0
      ),
      Value = c("scrap/parameter file format/WHAM_TEST.wdat",
                "Cu",
                "BL1",
                "BL1-Cu"),
      SpecialDef = c("WHAM", "Metal", "BL", "BL-Metal")
    ),
    MCToRemove = "Water_DonnanFA"
  ),
  CA = 0.0001,
  Species = "Test",
  TestType = "test",
  Duration = "test",
  Lifestage = "test",
  Endpoint = "test",
  Quantifier = "test",
  References = "test",
  Miscellaneous = "test"
)
ThisProblem_toxmode$DefComp$ActCorr = "None"
ThisProblem_toxmode$Comp$ActCorr = "None"
ThisProblem_toxmode$Spec$ActCorr = "None"
ThisProblem_toxmode$InVar$Name[1] = "Temp"
AllInput_toxmode = MatchInputsToProblem(
  DFInputs = data.frame(
    ID = as.character(1:(3*2)),
    Temp = 15,
    pH = 10,
    DOC = rep(c(0.1, 5, 20), each = 2, times = 1),
    Cu = 1E-7,
    Ca = rep(c(3.78e-02, 3.78e-05), each = 1, times = 3),
    SO4 = rep(c(2.46e-02, 2.46e-05), each = 1, times = 3)
  ),
  ThisProblem = ThisProblem_toxmode
)


ThisProblem_scratch = RemoveMassCompartments(
  ThisProblem = AddSpecialDefs(
    ThisProblem = AddSpecies(
      ThisProblem = AddInComps(
        ThisProblem = AddInVars(
          ThisProblem = water_problem,
          InVarName = "DOC",
          InVarMCName = "Water",
          InVarType = "WHAM-FA"),
        InCompName = c("Ca","SO4"),
        InCompCharge = c(2, -2),
        InCompMCName = "Water",
        InCompType = "MassBal",
        InCompActCorr = "None"
      ),
      SpecEquation = "CaSO4 = 1 * Ca + 1 * SO4",
      SpecMCName = "Water",
      SpecLogK = 2.3,
      SpecActCorr = "None",
      SpecDeltaH = 0,
      SpecTempKelvin = 0
    ),
    Value = "scrap/parameter file format/WHAM_TEST.wdat",
    SpecialDef = "WHAM"
  ),
  MCToRemove = "Water_DonnanFA"
)
ThisProblem_scratch$DefComp$ActCorr = "None"
ThisProblem_scratch$Comp$ActCorr = "None"
ThisProblem_scratch$Spec$ActCorr = "None"
ThisProblem_scratch$InVar$Name[1] = "Temp"
ThisProblem_scratch$Spec$DeltaH = 0.0
ThisProblem_scratch$Spec$TempKelvin = 0.0
AllInput_scratch = MatchInputsToProblem(
  DFInputs = data.frame(
    ID = as.character(1:(3*2)),
    Temp = 15,
    pH = 10,
    DOC = rep(c(0.1, 5, 20), each = 2, times = 1),
    Ca = rep(c(3.78e-02, 3.78e-05), each = 1, times = 3),
    SO4 = rep(c(2.46e-02, 2.46e-05), each = 1, times = 3)
  ),
  ThisProblem = ThisProblem_scratch
)
ResultsTable <- BLM(
  ThisProblem = ThisProblem_scratch,
  AllInput = AllInput_scratch,
  DoTox = FALSE,
  QuietFlag = "Quiet",
  ConvergenceCriteria = 0.0001,
  MaxIter = 100
)
sum(ResultsTable$Miscellaneous$Status == "Not Converged") / nrow(ResultsTable$Miscellaneous)
ResultsTable$Inputs[ResultsTable$Miscellaneous$Status == "Okay",]


ThisProblem_noISFX = AddSpecies(
  ThisProblem = AddInComps(
    ThisProblem = water_problem,
    InCompName = c("Ca","SO4","FA1H","FA2H","FA12H"),
    InCompCharge = c(2, -2, 0, 0, 0),
    InCompMCName = "Water",
    InCompType = "MassBal",
    InCompActCorr = "None"
  ),
  SpecEquation = c(
    "CaSO4 = 1 * Ca + 1 * SO4",
    "FA1 = 1 * FA1H - 1 * H",
    "FA2 = 1 * FA2H - 1 * H",
    "FA1-Ca = 1 * FA1H - 1 * H + 1 * Ca",
    "FA2-Ca = 1 * FA2H - 1 * H + 1 * Ca",
    "FA1H2 = 1 * FA12H - 1 * H",
    "FA2H1 = 1 * FA12H - 1 * H",
    "FA12 = 1 * FA12H - 2 * H",
    "FA12-Ca = 1 * FA12H - 2 * H + 1 * Ca"
  ),
  SpecMCName = "Water",
  SpecLogK = c(2.3, -1.59, -2.703333, -2.2, -2.2, -1.59, -2.703333, -4.293333, -4.4),
  SpecActCorr = "None",
  SpecDeltaH = 0,
  SpecTempKelvin = 0
)
ThisProblem_noISFX$Spec$ActCorr = "None"
ThisProblem_noISFX$Comp$ActCorr = "None"
ThisProblem_noISFX$DefComp$ActCorr = "None"
ThisProblem_noISFX$InVar$Name[1] = "Temp"
ThisProblem_noISFX$Spec$DeltaH = 0.0
ThisProblem_noISFX$Spec$TempKelvin = 0.0

InputDF$FA1H = InputDF$DOC * 1.419e-06
InputDF$FA2H = InputDF$DOC * 1.419e-06
InputDF$FA12H = InputDF$DOC * 2.365e-07

AllInput_noISFX = MatchInputsToProblem(
  DFInputs = InputDF,
  # DFInputs = data.frame(
  #   ID = as.character(1:(3*2)),
  #   Temp = 15,
  #   pH = 10,
  #   Ca = rep(c(3.78e-02, 3.78e-05), each = 1, times = 3),
  #   SO4 = rep(c(2.46e-02, 2.46e-05), each = 1, times = 3),
  #   FA1H = rep(c(0.1, 5, 20), each = 2, times = 1) * 1.419e-06,
  #   FA2H = rep(c(0.1, 5, 20), each = 2, times = 1) * 1.419e-06,
  #   FA12H = rep(c(0.1, 5, 20), each = 2, times = 1) * 2.365e-07
  # ),
  ThisProblem = ThisProblem_noISFX
)
ResultsTable_noISFX <- BLM(
  ThisProblem = ThisProblem_noISFX,
  AllInput = AllInput_noISFX,
  DoTox = FALSE,
  QuietFlag = "Quiet",
  ConvergenceCriteria = 0.0001,
  MaxIter = 1000
)
sum(ResultsTable_noISFX$Miscellaneous$Status == "Not Converged") / nrow(ResultsTable_noISFX$Miscellaneous)
ResultsTable_noISFX$Inputs[ResultsTable_noISFX$Miscellaneous$Status == "Okay",]
summary(ResultsTable_noISFX$Miscellaneous[ResultsTable_noISFX$Miscellaneous$Status == "Okay",])
summary(ResultsTable_noISFX$Miscellaneous[ResultsTable_noISFX$Miscellaneous$Status == "Not Converged",])
ResultsTable_noISFX$Concentrations

plot(ResultsTable_noISFX$Miscellaneous[ResultsTable_noISFX$Miscellaneous$Status == "Okay",c("FinalIter","IonicStrength","ChargeBalance")])

# even when there's no IS effects on organic matter, we still end up with a
# problem that does not converge, but only at high IS...interesting.
