# This script contains data documentation.

#' @title Molecular and atomic weights
#'
#' @description `MW` is a named list of molecular and atomic weights. The name
#'   of the list element is the symbol or formula of the chemical element or
#'   molecule (e.g. find hydrogen with "H", carbon dioxide as "CO2").
#'
#' @usage MW
#'
#' @source Prohaska, T., Irrgeher, J., Benefield, J., Böhlke, J., Chesson, L.,
#'   Coplen, T., Ding, T., Dunn, P., Gröning, M., Holden, N., Meijer, H.,
#'   Moossen, H., Possolo, A., Takahashi, Y., Vogl, J., Walczyk, T., Wang, J.,
#'   Wieser, M., Yoneda, S., Zhu, X. & Meija, J. (2022). Standard atomic weights
#'   of the elements 2021 (IUPAC Technical Report). Pure and Applied Chemistry,
#'   94(5), 573-600. https://doi.org/10.1515/pac-2019-0603
#'
#' @examples
#' # check that the molecular weight of CaCO3 is the same as Ca + C + O * 3
#' sum(MW[c("Ca", "C")], MW["O"] * 3) #=100.086
#' MW["CaCO3"] #=100.086
#'
"MW"

#' @title Carbonate system problem
#'
#' @description An example BLMEngineInR problem object, which describes a
#'   water-only system with only the (closed) carbonate system.
#'
#' @usage carbonate_system_problem
#'
#' @details This problem consists of two components (hydrogen "H" and carbonate
#'   "CO3") and three reactions (dissociation of water/formation of hydroxide
#'   "OH", formation of bicarbonate "HCO3" and formation of carbonic acid
#'   "H2CO3"). The pH and temperature are supplied as input variables, and the
#'   input label "ID" is supplied as well.
#'
#' @examples
#' print(carbonate_system_problem$Comp[, c("Name", "Charge", "Type")])
#' print(carbonate_system_problem$Spec[, c("Equation", "Charge", "ActCorr", "LogK", "DeltaH")])
#'
"carbonate_system_problem"

#' @title Water mass compartment only problem
#'
#' @description An example BLMEngineInR problem object, which describes a
#'   water-only system with no input variables or components yet, and the input
#'   label "ID".
"water_MC_problem"

#' @title Water-only problem
#'
#' @description An example BLMEngineInR problem object, which describes a
#'   water-only system with pH and temperature are supplied as input variables,
#'   and the input label "ID" is supplied as well. The only reaction is water
#'   dissociation (hydroxide "OH" formation reaction).
"water_problem"

#' @title Copper problem with WHAM V organic matter
#'
#' @description An example BLMEngineInR problem object, which describes a system
#'   with organic matter represented by WHAM V, and all of the common cations
#'   (Ca, Mg, Na, K) and anions (SO4, Cl, CO3) represented with their usual
#'   reactions. Copper is also represented as the toxic metal binding to a
#'   biotic ligand, and some example critical accumulations values are provided,
#'   including one for the United States Environmental Protection Agency's
#'   (USEPA) final acute value (FAV).
"Cu_full_organic_problem"

#' @title Cu problem with only inorganic components
#'
#' @description An example BLMEngineInR problem object, which describes a system
#'   with all of the common cations (Ca, Mg, Na, K) and anions (SO4, Cl, CO3)
#'   represented with their usual reactions. Copper is also represented as the
#'   toxic metal binding to a biotic ligand, and some example critical
#'   accumulations values are provided including one for the United States
#'   Environmental Protection Agency's (USEPA) final acute value (FAV). These
#'   critical accumulation values are the ones calibrated from the full organic
#'   model, as the DOC complexation should not affect the amount of organic
#'   matter required to induce a toxic effect, in theory. This will not give
#'   accurate predictions of toxicity when DOC is present in the water.
#'
"Cu_full_inorganic_problem"
