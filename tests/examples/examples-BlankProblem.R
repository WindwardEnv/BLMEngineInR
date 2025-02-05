# Make a blank problem:
ThisProblem = BlankProblem()
str(ThisProblem)

# Add Water as a mass compartment
ThisProblem = AddMassCompartments(
  ThisProblem,
  MassName = "Water",
  MassAmt = 1,
  MassUnit = "L"
)

# Add temperature and pH variables:
ThisProblem = AddInVars(ThisProblem,
                       InVarName = c("Temp", "pH"),
                       InVarMCName = rep("Water", 2),
                       InVarType = c("Temperature", "pH"))

# Add CO3 as a component:
ThisProblem = AddInComps(
  ThisProblem,
  InCompName = "CO3",
  InCompCharge = -2,
  InCompMCName = "Water",
  InCompType = "MassBal",
  InCompActCorr = "Debye"
)

# Add reactions (using SpecCompNames and SpecCompStoichs for arguments):
# HCO3 = H + CO3     logK = 10.329
# H2CO3 = 2*H + CO3  logK = 10.329 + 6.352 = 16.681
ThisProblem = AddSpecies(
  ThisProblem,
  SpecName = c("HCO3", "H2CO3"),
  SpecMCName = "Water",
  SpecActCorr = "Debye",
  SpecCompNames = list(c("H", "CO3"), c("H", "CO3")),
  SpecCompStoichs = list(c(1, 1), c(1, 2)),
  SpecLogK = c(10.329, 16.681),
  SpecDeltaH = c(-14997.55155, -24166.23162),
  SpecTempKelvin = 298.1514609
)
# ...ThisProblem now simulates carbonate reactions.

# Add major ions and copper as components
ThisProblem = AddInComps(
  ThisProblem,
  InCompName = c("Cu", "Ca", "Mg", "Na", "K", "SO4", "Cl", "S"),
  InCompCharge = c(2, 2, 2, 1, 1, -2, -1, -2),
  InCompMCName = "Water",
  InCompType = "MassBal",
  InCompActCorr = "Debye"
)

# Add reactions (using SpecEquation as an argument):
ThisProblem = AddSpecies(
  ThisProblem,
  SpecEquation = c(
    "CuOH = 1 * Cu + 1 * OH",
    "Cu(OH)2 = 1 * Cu + 2 * OH",
    "CuSO4 = 1 * Cu + 1 * SO4",
    "CuCl = 1 * Cu + 1 * Cl",
    "CuCO3 = 1 * Cu + 1 * CO3",
    "Cu(CO3)2 = 1 * Cu + 2 * CO3",
    "CuHCO3 = 1 * Cu + 1 * CO3 + 1 * H",
    "CaHCO3 = 1 * Ca + 1 * H + 1 * CO3",
    "CaCO3 = 1 * Ca + 1 * CO3",
    "CaSO4 = 1 * Ca + 1 * SO4",
    "MgHCO3 = 1 * Mg + 1 * H + 1 * CO3",
    "MgCO3 = 1 * Mg + 1 * CO3",
    "MgSO4 = 1 * Mg + 1 * SO4"
  ),
  SpecMCName = "Water",
  SpecActCorr = "Debye",
  SpecLogK = c(6.48, 11.78, 2.360, 0.400, 6.750, 9.920, 14.620,
               11.44, 3.22, 2.30, 11.4, 2.98, 2.37),
  SpecDeltaH = c(0, 0, 8844.385918, 6738.57975, 0, 0, 0,
                 -3664.102737, 14951.22381, 6949.160364, -11666.16619,
                 11413.46945, 19163.83616),
  SpecTempKelvin = 298.15
)

# Add BL mass compartment:
ThisProblem = AddMassCompartments(
  ThisProblem,
  MassName = "BL",
  MassAmt = 1,
  MassUnit = "kg wet"
)

# Add BL1 defined component:
ThisProblem = AddDefComps(ThisProblem,
                         DefCompName = "BL1",
                         DefCompFromNum = 1.78E-5,
                         DefCompCharge = -1,
                         DefCompMCName = "BL",
                         DefCompType = "MassBal",
                         DefCompActCorr = "None",
                         DefCompSiteDens = 3E-5)

# Add biotic ligand reactions (using SpecStoich):
spec_that_bind = c("Cu", "CuOH", "Ca", "Mg", "H", "Na")
temp_stoich_mat = ThisProblem$SpecStoich[spec_that_bind, ]
rownames(temp_stoich_mat) = paste0("BL1-", spec_that_bind)
temp_stoich_mat[, "BL1"] = 1L
temp_stoich_mat["BL1-CuOH", c("H","OH")] = c(-1L, 0L)
ThisProblem = AddSpecies(
  ThisProblem,
  SpecName = paste0("BL1-", spec_that_bind),
  SpecMCName = "BL",
  SpecActCorr = "None",
  SpecLogK = c(7.4, -1.3, 3.6, 3.6, 5.4, 3.0),
  SpecDeltaH = ThisProblem$Spec$DeltaH[match(spec_that_bind, ThisProblem$Spec$Name)],
  SpecTempKelvin = ThisProblem$Spec$TempKelvin[match(spec_that_bind, ThisProblem$Spec$Name)],
  SpecStoich = temp_stoich_mat
)

# Add special definitions for the toxic species:
ThisProblem = AddSpecialDefs(
  ThisProblem,
  Value = c("BL1","Cu","BL1-Cu","BL1-CuOH"),
  SpecialDef = c("BL","Metal","BLMetal","BLMetal")
)
# ...ThisProblem now simulates copper toxicity in the absence of organic matter.

# Add DOC: first we add DOC and HA input variables...
ThisProblem = AddInVars(
  ThisProblem,
  InVarName = c("DOC", "HA"),
  InVarMCName = "Water",
  InVarType = c("WHAM-HAFA", "PercHA")
)

# ...then we add a WHAM version as a special definition.
ThisProblem = AddSpecialDefs(
  ThisProblem,
  Value = "V",
  SpecialDef = "WHAM"
)

# As a finishing touch, we already know our critical values:
ThisProblem = AddCriticalValues(
  ThisProblem,
  CATab = data.frame(
    CA = c(0.05541, 0.03395),
    Species = c("Ceriodaphnia dubia","FAV"),
    TestType = "Acute",
    Duration = c("48h","DIV=2.00"),
    Lifestage = c("Neonate (<24h)","ACR=3.22"),
    Endpoint = c("Mortality","FAV"),
    Quantifier = c("EC50; LC50", "NA"),
    References = c("Gensemer et al. 2002; Hyne et al. 2005; 2002; 2003; Van Genderen et al. 2007",
                   "US EPA 2007"),
    Miscellaneous = c("SMEA calculated by median", NA)
  )
)
# ThisProblem can now calculate the Cu WQC

# Now what about CO2 dissolving from the atmosphere?
ThisProblem = AddPhases(
  ThisProblem,
  PhaseName = "CO2(g)",
  PhaseCompNames = list(c("CO3", "H")),
  PhaseCompStoichs = list(c(1, 2)),
  PhaseLogK = -1.5,
  PhaseDeltaH = 0,
  PhaseTempKelvin = 0,
  PhaseMoles = 10^-3.2
)

# Actually, scratch that - no CO2 dissolution
ThisProblem = RemovePhases(ThisProblem, "CO2(g)")

# Actually, I don't want C. dubia in this parameter file.
ThisProblem = RemoveCriticalValues(ThisProblem, 1)

# I know we usually have sulfide in there, but it's really not doing anything
# for us, so let's remove that.
ThisProblem = RemoveComponents(ThisProblem, "S")

# I kinda want to try this with WHAM VII instead of V...
ThisProblem = RemoveSpecialDefs(ThisProblem, SpecialDefToRemove = "WHAM")
ThisProblem = AddSpecialDefs(ThisProblem, Value = "VII", SpecialDef = "WHAM")

# Now what if I wanted to make this just a simulation of organic matter binding,
# sans biotic ligand?
ThisProblem = RemoveMassCompartments(ThisProblem, MCToRemove = "BL")
