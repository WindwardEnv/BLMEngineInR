% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CommonParameterDefinitions.R
\name{CommonParameterDefinitions}
\alias{CommonParameterDefinitions}
\title{Common Parameter Definitions}
\arguments{
\item{NMass}{integer, the number of mass compartments.}

\item{MassName}{character vector (NMass), The name of each mass compartment.}

\item{MassAmt}{numeric vector (NMass), The amount of each mass compartment.}

\item{MassUnit}{character vector (NMass), The units for each mass
compartment.}

\item{AqueousMCR}{integer, the (1-based) position of the water/aqueous mass
compartment.}

\item{BioticLigMCR}{integer, the (1-based) position of the biotic ligand mass
compartment(s).}

\item{WHAMDonnanMCR}{integer (2), the mass compartments corresponding to the
humic acid (1) and fulvic acid (2) Donnan layers.}

\item{WHAMDonnanMC}{integer (2), the mass compartments corresponding to the
humic acid (0) and fulvic acid (1) Donnan layers.}

\item{NInLab}{integer, the number of input label fields}

\item{InLabName}{character vector (NInLab), The names of the input label
fields.}

\item{NInVar}{integer, the number of input variables}

\item{InVarName}{character vector (NInVar), The name of each input variable.}

\item{InVarMCR}{integer vector (NInVar), The mass compartment of each input
variable. (1-based)}

\item{InVarMC}{integer vector (NInVar), The mass compartment of each input
variable. (0-based)}

\item{InVarType}{character vector (NInVar), The type of each input variable.
Should be one of "Temperature" (the temperature in degrees C), "pH" (the
-log[H]...you know, pH), "WHAM-HA", "WHAM-FA", "WHAM-HAFA" (Windemere Humic
Aqueous Model organic matter (input mg C/L), as all humic acid, all fulvic
acid, or a mix of humics and fulvics, respectively.), "PercHA" (optionally
indicate the percent humic acid in a the WHAM-HAFA component for that
compartment.), or "PercAFA" (optionally indicate the percent of active
fulvic acid for the WHAM-FA or WHAM-HAFA component for that compartment)}

\item{NInComp}{integer, the number of input components}

\item{InCompName}{character vector (NInComp), The names of the input
components.}

\item{NDefComp}{integer, the number of defined components}

\item{DefCompName}{character vector (NDefComp), the names of each defined
component}

\item{DefCompFromNum}{numeric vector (NDefComp), the number used for deriving
the concentration of each defined component}

\item{DefCompFromVar}{character vector (NDefComp), the variable used for
deriving the concentration of each defined component}

\item{DefCompCharge}{signed integer vector (NDefComp), the charge of each
defined component}

\item{DefCompMCR}{integer vector (NDefComp), the mass compartment number each
defined component (1-based)}

\item{DefCompMC}{integer vector (NDefComp), the mass compartment number each
defined component (0-based)}

\item{DefCompType}{character vector (NDefComp), the type of each defined
component}

\item{DefCompActCorr}{character vector (NDefComp), the activity correction
method to use with each defined component}

\item{DefCompSiteDens}{numeric vector (NDefComp), the site density of each
defined component}

\item{NComp}{integer, the combined number of components in the simulation,
including the input components, defined components (and including the
defined components that get added by ExpandWHAM)}

\item{CompName}{character vector (NComp), the name of each component in the
simulation}

\item{CompCharge}{signed integer vector (NComp), the charge of each component
in the simulation}

\item{CompMCR}{integer vector (NComp), the mass compartment of each component
in the simulation (1-based)}

\item{CompMC}{integer vector (NComp), the mass compartment of each component
in the simulation (0-based)}

\item{CompCtoM}{numeric vector (NSpec), the concentration to mass conversion
factor of the components}

\item{CompType}{character vector (NComp), the type of each component in the
simulation}

\item{CompActCorr}{character vector (NComp), the activity correction method
of each component in the simulation}

\item{CompSiteDens}{numeric vector (NComp), the site density of each
component in the simulation}

\item{CompConc}{numeric vector (NComp), the free ion concentrations of each
component in the simulation}

\item{TotConc}{numeric vector (NComp), the total concentrations of each
component in the simulation (units of e.g., mol/L and mol/kg)}

\item{TotMoles}{numeric vector (NComp), the total moles of each component in
the simulation (units of mol)}

\item{NSpec}{integer, the number of chemical species for which we have
formation reactions in the simulation}

\item{SpecName}{character vector (NSpec), the name of the chemical species
for which we have formation reactions}

\item{SpecMCR}{integer vector (NSpec), the mass compartment of the chemical
species for which we have formation reactions (1-based)}

\item{SpecMC}{integer vector (NSpec), the mass compartment of the chemical
species for which we have formation reactions (0-based)}

\item{SpecActCorr}{character vector (NSpec), the activity correction method
of the chemical species for which we have formation reactions}

\item{SpecNC}{integer vector (NSpec), the number of components for the
formation reactions}

\item{SpecCompList}{integer matrix (NSpec x max(SpecNC)), the list of
components for the formation reactions}

\item{SpecCtoM}{numeric vector (NSpec), the concentration to mass conversion
factor of the chemical species for which we have formation reactions}

\item{SpecCharge}{signed integer vector (NSpec), the charge of the chemical
species for which we have formation reactions}

\item{SpecK}{numeric vector (NSpec), the equilibrium coefficient of the
formation reactions}

\item{SpecLogK}{numeric vector (NSpec), the log10-transformed equilibrium
coefficient of the formation reactions}

\item{SpecDeltaH}{numeric vector (NSpec), the enthalpy change of the
formation reactions}

\item{SpecTempKelvin}{numeric vector (NSpec), the temperature associated with
K/logK and DeltaH of the formation reactions}

\item{SpecStoich}{signed integer matrix (NSpec x NComp), the reaction
stoichiometry of the formation reactions}

\item{SpecConc}{numeric vector (NSpec), the concentrations of each species
for which we have formation reactions}

\item{SpecMoles}{numeric vector (NSpec), the moles of each species for which
we have formation reactions}

\item{NPhase}{integer, the number of phases in the phase list}

\item{PhaseName}{character vector (NPhase), the name of the phases for which
we have phase reactions}

\item{PhaseNC}{integer vector (NPhase), the number of components for the
phase reactions}

\item{PhaseCompList}{integer matrix (NPhase x max(PhaseNC)), the list of
components for the phase reactions}

\item{PhaseStoich}{signed integer matrix (NPhase x NComp), the reaction
stoichiometry for the phase reactions}

\item{PhaseK}{numeric vector (NPhase), the equilibrium coefficient for the
phase reactions}

\item{PhaseLogK}{numeric vector (NPhase), the log10-transformed equilibrium
coefficient for the phase reactions}

\item{PhaseDeltaH}{numeric vector (NPhase), the enthalpy change for the phase
reactions}

\item{PhaseTemp}{numeric vector (NPhase), the temperature associated with
K/logK and DeltaH for the phase reactions}

\item{PhaseMoles}{numeric vector (NPhase), the number of moles of the phases
for which we have phase reactions}

\item{NSpecialDef}{integer, the number of special definitions in the
parameter file, including biotic ligands, metals, WHAM versions, etc.}

\item{NBL}{integer, the number of biotic ligand components associated with
toxic effects...typically one...and things might get messed up if it's not
one.}

\item{NMetal}{integer, the number of metal components associated with toxic
effects...typically one...and things might get messed up if it's not one.}

\item{NBLMetal}{integer, the number of biotic ligand-bound metal species that
are associated with toxic effects.}

\item{BLName}{The name of the component that corresponds to the biotic ligand
associated with toxic effects.}

\item{MetalName}{The name of the component that corresponds to the metal
associated with toxic effects.}

\item{BLMetalName}{The names of the species that are the biotic ligand-bound
metal associated with toxic effects.}

\item{BLCompR}{integer vector (NBL), the (1-based) position of the biotic
ligand component(s) in the component arrays}

\item{BLComp}{integer vector (NBL), the (0-based) position of the biotic
ligand component(s) in the component arrays}

\item{MetalCompR}{integer vector (NMetal), the (1-based) position of the
metal component(s) in the component arrays (i.e., which is the toxic metal
component)}

\item{MetalComp}{integer vector (NMetal), the (0-based) position of the metal
component(s) in the component arrays (i.e., which is the toxic metal
component)}

\item{BLMetalSpecsR}{integer vector (NBLMetal), the (1-based) positions of
the species in the arrays which contribute to toxicity (i.e., which species
are the toxic metal bound to the relevant biotic ligand)}

\item{BLMetalSpecs}{integer vector (NBLMetal), the (0-based) positions of the
species in the arrays which contribute to toxicity (i.e., which species are
the toxic metal bound to the relevant biotic ligand)}

\item{DoWHAM}{logical, TRUE = there are WHAM species, FALSE = no WHAM species}

\item{WHAMDLF}{numeric (2), WHAM's Double layer overlap factor}

\item{WHAMKZED}{numeric (2), WHAM's Constant to control DDL at low ZED}

\item{SpecKsel}{numeric (NSpec, 2), WHAM's Selectivity coefficient Ksel for
diffuse layer binding}

\item{WHAMP}{numeric (2), WHAM's P parameter...}

\item{WHAMRadius}{numeric (2), WHAM's molecular radius parameter for organic
matter}

\item{WHAMMolWt}{numeric (2), WHAM's molecular weight parameter for organic
matter}

\item{HumicSubstGramsPerLiter}{numeric (2), grams per liter of each organic
matter component (HA and FA) in solution}

\item{CATab}{data frame, the critical accumulation table from the parameter
file.}

\item{NCAT}{integer, the number of critical accumulations in the parameter
file table.}

\item{CATarget}{numeric, the target critical accumulation in units of mol /
kg (only used when DoTox == TRUE)}

\item{NObs}{integer; the number of chemistry observations}

\item{InLabObs}{character matrix with NObs rows and InLab columns; the input
labels for each observation}

\item{InVarObs}{matrix with NObs rows and InVar columns; the input variables
for each observation}

\item{InCompObs}{matrix with NObs rows and InComp columns; the input
component concentrations for each observation}

\item{SysTempCelsiusObs}{numeric vector of length NObs; input temperatures,
in Celsius}

\item{SysTempKelvinObs}{numeric vector of length NObs; input temperatures, in
Kelvin}

\item{SysTempCelsius}{double; input temperature for the current observation,
in Celsius}

\item{SysTempKelvin}{double; input temperature for the current observation,
in Kelvin}

\item{TotConcObs}{numeric matrix with NObs rows and NComp columns; the total
concentrations of each component, including derived components}

\item{pH}{numeric vector NObs; input pH for each observation}

\item{FinalIter}{integer, the number of Newton-Raphson iterations that we
needed to reach convergence}

\item{FinalMaxError}{numeric, the highest final absolute error fraction
=max(abs(Resid / TotMoles))}

\item{MaxError}{numeric, the highest absolute error fraction in this
iteration =max(abs(Resid / TotMoles))}

\item{CalcTotConc}{numeric vector (NComp), the calculated total
concentrations of each component in the simulation (units of e.g., mol/L
and mol/kg)}

\item{QuietFlag}{character, one of "Very Quiet" (only print out when run is
done), "Quiet" (print out Obs=iObs), or "Debug" (print out lots of info)}

\item{DoTox}{logical, TRUE for toxicity mode where the MetalName component
concentration is adjusted to try to match the CATarget with BLMetalSpecs}

\item{ConvergenceCriteria}{numeric, the maximum value of MaxError that counts
as convergence by the Newton-Raphson root-finding algorithm}

\item{MaxIter}{integer, the maximum number of iterations the Newton-Raphson
root-finding algorithm should do before giving up}

\item{IonicStrength}{double, the ionic strength of the solution}
}
\description{
These are parameters that are commonly used in the BLMEngineInR
  package. They will appear throughout the various internal functions, and
  this central repository of their definitions is helpful.
}
\seealso{
Other BLMEngine Functions: 
\code{\link{GetData}()},
\code{\link{MatchInputsToProblem}()},
\code{\link{ReadInputsFromFile}()}
}
\concept{BLMEngine Functions}
