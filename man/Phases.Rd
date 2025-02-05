% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Phases.R
\name{Phases}
\alias{Phases}
\alias{AddPhases}
\alias{RemovePhases}
\title{Add or remove phase reactions in a problem}
\usage{
AddPhases(
  ThisProblem,
  PhaseEquation = character(),
  PhaseName = character(),
  PhaseCompNames = list(),
  PhaseCompStoichs = list(),
  PhaseStoich = NULL,
  PhaseLogK,
  PhaseDeltaH,
  PhaseTempKelvin,
  PhaseMoles,
  DoCheck = TRUE
)

RemovePhases(ThisProblem, PhasesToRemove, DoCheck = TRUE)
}
\arguments{
\item{ThisProblem}{A list object with a structure like that returned by
`BlankProblem()`.}

\item{PhaseEquation}{A character vector giving the chemical equation for a
formation reaction. This must include the stoichiometric coefficients for
each reactant, even if it's 1. (e.g., the equation for the formation of
calcium chloride would be `"CaCl2 = 1 * Ca + 2 * Cl"`). If `PhaseName` is
also supplied, then a partial equation with just the right hand side
(reactants) can be supplied (i.e., `"= 1 * Ca + 2 * Cl"`). Can be omitted
if either `PhaseStoich` or both `PhaseCompNames` and `PhaseCompStoichs` are
supplied.}

\item{PhaseName}{A character vector with the name(s) of the species to add
formation reactions for. Can be omitted if `SpecEquation` indicates the
phase name.}

\item{PhaseCompNames}{A list where each element is a character vector of the
component names used to form each phase. See examples for clarification.
Can be omitted if `PhaseEquation` or `PhaseStoich` is supplied.}

\item{PhaseCompStoichs}{A list where each element is an integer vector of the
stoichiometric coefficients of each component used to form each phase. See
examples for clarification. Can be omitted if `PhaseEquation` or
`PhaseStoich` is supplied.}

\item{PhaseStoich}{A matrix of stoichiometric coefficients, where each row
corresponds to a phase reaction and each column corresponds to a component.
The columns should match `ThisProblem$Comp$Name` exactly. Can be omitted if
either `PhaseEquation` or both `PhaseCompNames` and `PhaseCompStoichs` are
supplied.}

\item{PhaseLogK}{A numeric vector with the log10-transformed equilibrium
coefficients of the phase formation reactions.}

\item{PhaseDeltaH}{A numeric vector with the change in enthalpy of the phase
formation reactions.}

\item{PhaseTempKelvin}{A numeric vector with the temperatures (in Kelvin)
corresponding to `PhaseDeltaH` values of the phase formation reactions.}

\item{PhaseMoles}{A numeric vector with the moles of the phase.}

\item{DoCheck}{A logical value indicating whether checks should be performed
on the incoming and outgoing problem objects. Defaults to `TRUE`, as you
usually want to make sure something isn't awry, but the value is often set
to `FALSE` when used internally (like in DefineProblem) so the problem is
only checked once at the end.}

\item{PhasesToRemove}{A character or integer vector indicating the names or
indices (respectively) of the phase formation reactions to remove.}
}
\value{
`ThisProblem`, with the phase reaction(s) changed.
}
\description{
PHASES ARE NOT CURRENTLY IMPLEMENTED. This function is here for
  as a placeholder since it will require much of the same support
  infrastructure once it is implemented, but no reactions are processed in
  CHESS.
}
\examples{
print(carbonate_system_problem$Phase)
my_new_problem = carbonate_system_problem
my_new_problem = AddPhases(ThisProblem = my_new_problem,
                           PhaseEquation = "CO2(g) = 1 * CO3 + 2 * H",
                           PhaseLogK = -1.5,
                           PhaseDeltaH = 0,
                           PhaseTempKelvin = 0,
                           PhaseMoles = 10^-3.5)
print(my_new_problem$Phase)

}
\seealso{
Other problem manipulation functions: 
\code{\link{BlankProblem}()},
\code{\link{Components}},
\code{\link{CriticalValues}},
\code{\link{InLabs}},
\code{\link{InVars}},
\code{\link{MassCompartments}},
\code{\link{SpecialDefs}},
\code{\link{Species}}
}
\concept{problem manipulation functions}
