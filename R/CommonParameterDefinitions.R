#' @name CommonParameterDefinitions
#'
#' @title Common Parameter Definitions
#'
#' @description These are parameters that are commonly used in the BLMEngineInR
#'   package. They will appear throughout the various internal functions, and
#'   this central repository of their definitions is helpful.
#'
#' @family BLMEngine Functions
#'
#' @param NMass integer, the number of mass compartments.
#' @param MassName character vector (NMass), The name of each mass compartment.
#' @param MassAmt numeric vector (NMass), The amount of each mass compartment.
#' @param MassUnit character vector (NMass), The units for each mass
#'   compartment.
#' @param AqueousMC integer, the (1-based) position of the water/aqueous mass
#'   compartment.
#' @param BioticLigMC integer, the (1-based) position of the biotic ligand mass
#'   compartment(s).
#' @param WHAMDonnanMC integer (2), the mass compartments corresponding to the
#'   humic acid (0) and fulvic acid (1) Donnan layers.
#'
#' @param NInLab integer, the number of input label fields
#' @param InLabName character vector (NInLab), The names of the input label
#'   fields.
#'
#' @param NInVar integer, the number of input variables
#' @param InVarName character vector (NInVar), The name of each input variable.
#' @param InVarMC integer vector (NInVar), The mass compartment of each input
#'   variable.
#' @param InVarType character vector (NInVar), The type of each input variable.
#'   Should be one of "Temperature" (the temperature in degrees C), "pH" (the
#'   -log[H]...you know, pH), "WHAM-HA", "WHAM-FA", "WHAM-HAFA" (Windemere Humic
#'   Aqueous Model organic matter (input mg C/L), as all humic acid, all fulvic
#'   acid, or a mix of humics and fulvics, respectively.), "PercHA" (optionally
#'   indicate the percent humic acid in a the WHAM-HAFA component for that
#'   compartment.), or "PercAFA" (optionally indicate the percent of active
#'   fulvic acid for the WHAM-FA or WHAM-HAFA component for that compartment)
#'
#' @param NInComp integer, the number of input components
#' @param InCompName character vector (NInComp), The names of the input
#'   components.
#'
#' @param NDefComp integer, the number of defined components
#' @param DefCompName character vector (NDefComp), the names of each defined
#'   component
#' @param DefCompFromNum numeric vector (NDefComp), the number used for deriving
#'   the concentration of each defined component
#' @param DefCompFromVar character vector (NDefComp), the variable used for
#'   deriving the concentration of each defined component
#' @param DefCompCharge signed integer vector (NDefComp), the charge of each
#'   defined component
#' @param DefCompMC integer vector (NDefComp), the mass compartment number each
#'   defined component
#' @param DefCompType character vector (NDefComp), the type of each defined
#'   component
#' @param DefCompActCorr character vector (NDefComp), the activity correction
#'   method to use with each defined component
#' @param DefCompSiteDens numeric vector (NDefComp), the site density of each
#'   defined component
#'
#' @param NComp integer, the combined number of components in the simulation,
#'   including the input components, defined components (and including the
#'   defined components that get added by ExpandWHAM)
#' @param CompName character vector (NComp), the name of each component in the
#'   simulation
#' @param CompCharge signed integer vector (NComp), the charge of each component
#'   in the simulation
#' @param CompMC integer vector (NComp), the mass compartment of each component
#'   in the simulation
#' @param CompCtoM numeric vector (NSpec), the concentration to mass conversion
#'   factor of the components
#' @param CompType character vector (NComp), the type of each component in the
#'   simulation
#' @param CompActCorr character vector (NComp), the activity correction method
#'   of each component in the simulation
#' @param CompSiteDens numeric vector (NComp), the site density of each
#'   component in the simulation
#' @param CompConc numeric vector (NComp), the free ion concentrations of each
#'   component in the simulation
#' @param TotConc numeric vector (NComp), the total concentrations of each
#'   component in the simulation (units of e.g., mol/L and mol/kg)
#' @param TotMoles numeric vector (NComp), the total moles of each component in
#'   the simulation (units of mol)
#'
#' @param NSpec integer, the number of chemical species for which we have
#'   formation reactions in the simulation
#' @param SpecName character vector (NSpec), the name of the chemical species
#'   for which we have formation reactions
#' @param SpecMC integer vector (NSpec), the mass compartment of the chemical
#'   species for which we have formation reactions
#' @param SpecActCorr character vector (NSpec), the activity correction method
#'   of the chemical species for which we have formation reactions
#' @param SpecNC integer vector (NSpec), the number of components for the
#'   formation reactions
#' @param SpecCompList integer matrix (NSpec x max(SpecNC)), the list of
#'   components for the formation reactions
#' @param SpecCtoM numeric vector (NSpec), the concentration to mass conversion
#'   factor of the chemical species for which we have formation reactions
#' @param SpecCharge signed integer vector (NSpec), the charge of the chemical
#'   species for which we have formation reactions
#' @param SpecK numeric vector (NSpec), the equilibrium coefficient of the
#'   formation reactions
#' @param SpecLogK numeric vector (NSpec), the log10-transformed equilibrium
#'   coefficient of the formation reactions
#' @param SpecDeltaH numeric vector (NSpec), the enthalpy change of the
#'   formation reactions
#' @param SpecTempKelvin numeric vector (NSpec), the temperature associated with
#'   K/logK and DeltaH of the formation reactions
#' @param SpecStoich signed integer matrix (NSpec x NComp), the reaction
#'   stoichiometry of the formation reactions
#' @param SpecConc numeric vector (NSpec), the concentrations of each species
#'   for which we have formation reactions
#' @param SpecMoles numeric vector (NSpec), the moles of each species for which
#'   we have formation reactions
#'
#' @param NPhase integer, the number of phases in the phase list
#' @param PhaseName character vector (NPhase), the name of the phases for which
#'   we have phase reactions
#' @param PhaseNC integer vector (NPhase), the number of components for the
#'   phase reactions
#' @param PhaseCompList integer matrix (NPhase x max(PhaseNC)), the list of
#'   components for the phase reactions
#' @param PhaseStoich signed integer matrix (NPhase x NComp), the reaction
#'   stoichiometry for the phase reactions
#' @param PhaseK numeric vector (NPhase), the equilibrium coefficient for the
#'   phase reactions
#' @param PhaseLogK numeric vector (NPhase), the log10-transformed equilibrium
#'   coefficient for the phase reactions
#' @param PhaseDeltaH numeric vector (NPhase), the enthalpy change for the phase
#'   reactions
#' @param PhaseTemp numeric vector (NPhase), the temperature associated with
#'   K/logK and DeltaH for the phase reactions
#' @param PhaseMoles numeric vector (NPhase), the number of moles of the phases
#'   for which we have phase reactions
#'
#' @param NSpecialDef integer, the number of special definitions in the
#'   parameter file, including biotic ligands, metals, WHAM versions, etc.
#' @param NBL integer, the number of biotic ligand components associated with
#'   toxic effects...typically one...and things might get messed up if it's not
#'   one.
#' @param NMetal integer, the number of metal components associated with toxic
#'   effects...typically one...and things might get messed up if it's not one.
#' @param NBLMetal integer, the number of biotic ligand-bound metal species that
#'   are associated with toxic effects.
#' @param BLName The name of the component that corresponds to the biotic ligand
#'   associated with toxic effects.
#' @param MetalName The name of the component that corresponds to the metal
#'   associated with toxic effects.
#' @param BLMetalName The names of the species that are the biotic ligand-bound
#'   metal associated with toxic effects.
#' @param BLComp integer vector (NBL), the (1-based) position of the biotic
#'   ligand component(s) in the component arrays
#' @param MetalComp integer vector (NMetal), the (1-based) position of the metal
#'   component(s) in the component arrays (i.e., which is the toxic metal
#'   component)
#' @param BLMetalSpecs integer vector (NBLMetal), the (1-based) positions of the
#'   species in the arrays which contribute to toxicity (i.e., which species are
#'   the toxic metal bound to the relevant biotic ligand)
#' @param DoWHAM logical, TRUE = there are WHAM species, FALSE = no WHAM species
#' @param wDLF numeric (2), WHAM's Double layer overlap factor
#' @param wKZED numeric (2), WHAM's Constant to control DDL at low ZED
#' @param SpecKsel numeric (NSpec, 2), WHAM's Selectivity coefficient Ksel for
#'   diffuse layer binding
#' @param wP numeric (2), WHAM's P parameter...
#' @param wRadius numeric (2), WHAM's molecular radius parameter for organic
#'   matter
#' @param wMolWt numeric (2), WHAM's molecular weight parameter for organic
#'   matter
#' @param HumicSubstGramsPerLiter numeric (2), grams per liter of each organic
#'   matter component (HA and FA) in solution
#'
#' @param CATab data frame, the critical accumulation table from the parameter
#'   file.
#' @param NCAT integer, the number of critical accumulations in the parameter
#'   file table.
#' @param CATarget numeric, the target critical accumulation in units of mol /
#'   kg (only used when DoTox == TRUE)
#'
#' @param NObs integer; the number of chemistry observations
#' @param InLabObs matrix with NObs rows and InLab columns; the input labels for
#'   each observation
#' @param InVarObs matrix with NObs rows and InVar columns; the input variables
#'   for each observation
#' @param InCompObs matrix with NObs rows and InComp columns; the input
#'   component concentrations for each observation
#' @param SysTempCelsiusObs numeric vector of length NObs; input temperatures,
#'   in Celsius
#' @param SysTempKelvinObs numeric vector of length NObs; input temperatures, in
#'   Kelvin
#' @param SysTempCelsius double; input temperature for the current observation,
#'   in Celsius
#' @param SysTempKelvin double; input temperature for the current observation,
#'   in Kelvin
#' @param TotConcObs numeric matrix with NObs rows and NComp columns; the total
#'   concentrations of each component, including derived components
#' @param pH numeric vector NObs; input pH for each observation
#'
#' @param FinalIter integer, the number of Newton-Raphson iterations that we
#'   needed to reach convergence
#' @param FinalMaxError numeric, the highest final absolute error fraction
#'   =max(abs(Resid / TotMoles))
#' @param MaxError numeric, the highest absolute error fraction in this
#'   iteration =max(abs(Resid / TotMoles))
#' @param CalcTotConc numeric vector (NComp), the calculated total
#'   concentrations of each component in the simulation (units of e.g., mol/L
#'   and mol/kg)
#' @param QuietFlag character, one of "Very Quiet" (only print out when run is
#'   done), "Quiet" (print out Obs=iObs), or "Debug" (print out lots of info)
#' @param DoTox logical, TRUE for toxicity mode where the MetalName component
#'   concentration is adjusted to try to match the CATarget with BLMetalSpecs
#' @param ConvergenceCriteria numeric, the maximum value of MaxError that counts
#'   as convergence by the Newton-Raphson root-finding algorithm
#' @param MaxIter integer, the maximum number of iterations the Newton-Raphson
#'   root-finding algorithm should do before giving up
#'
#' @param IonicStrength double, the ionic strength of the solution
#'
NULL
