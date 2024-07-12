#' @name StoichConversionFunctions
#'
#' @title Stoichiometry conversion functions
#'
#' @description These functions help with converting between equations (the
#'   intuitive way people think about formation reactions), stoichiometric
#'   matrices (what CHESS needs to define the problem), and paired lists of
#'   component names and stoichiometries (how parameter files define formation
#'   reactions).
#'
#' @param CompName a character vector of the component names, in order. These
#'   are the columns of the `SpecStoich` matrix.
#' @param SpecName a character vector of the chemical species names, in order.
#'   These are the rows of the `SpecStoich` matrix.
#' @param SpecCompNames a list object where each element is a vector of
#'   component names in a formation reaction.
#' @param SpecCompStoichs a list object where each element is a vector of the
#'   stoichiometric coefficients in the formation reaction corresponding to the
#'   components in `SpecCompNames`.
#' @param SpecNC an integer vector giving the number of reactants in each
#'   formation reaction.
#' @param SpecStoich A matrix of stoichiometric coefficients, where each row
#'   corresponds to a chemical species and each column corresponds to a
#'   component.
#' @param SpecEquation A character vector giving the chemical equation for a
#'   formation reaction. This must include the stoichiometric coefficients for
#'   each reactant, even if it's 1. (e.g., the equation for the formation of
#'   calcium chloride would be `"CaCl2 = 1 * Ca + 2 * Cl"`).
#'
#' @details The naming scheme of these functions is simple
#'
#' @returns `XToStoichMatrix` returns `SpecStoich`, `XToEquation` returns
#'   `SpecEquation`, and `XToStoichComps` returns a list with at least the two
#'   items `SpecCompNames` and `SpecCompStoichs`. `EquationToStoich`
NULL


#' @rdname StoichConversionFunctions
#' @keywords internal
StoichCompsToStoichMatrix = function(
    SpecCompNames = list(),
    SpecCompStoichs = list(),
    CompName = unique(unlist(SpecCompNames)),
    SpecName = names(SpecCompNames),
    SpecNC = as.integer(sapply(SpecCompNames, function(X) {sum(X != "")}))) {

  SpecStoich = matrix(
    0L,
    nrow = length(SpecName),
    ncol = length(CompName),
    dimnames = list(SpecName, CompName))
  for (i in 1:length(SpecName)){
    if (SpecNC[i] > 0) {
      for (j in 1:SpecNC[i]) {
        SpecStoich[i, SpecCompNames[[i]][j]] =
          SpecStoich[i, SpecCompNames[[i]][j]] +
          as.integer(SpecCompStoichs[[i]][j])
      }
    }
  }

  return(SpecStoich)

}


#' @rdname StoichConversionFunctions
#' @keywords internal
StoichCompsToEquation = function(
    SpecCompNames = list(),
    SpecCompStoichs = list(),
    CompName = unique(unlist(SpecCompNames)),
    SpecName = names(SpecCompNames),
    SpecNC = as.integer(sapply(SpecCompNames, function(X) {sum(X != "")}))) {

  StoichMatrixToEquation(
    SpecStoich = StoichCompsToStoichMatrix(
      CompName = CompName,
      SpecName = SpecName,
      SpecCompNames = SpecCompNames,
      SpecCompStoichs = SpecCompStoichs,
      SpecNC = SpecNC
    ),
    CompName = CompName,
    SpecName = SpecName
  )

}


#' @rdname StoichConversionFunctions
#' @keywords internal
StoichMatrixToEquation = function(SpecStoich = matrix(),
                            SpecName = rownames(SpecStoich),
                            CompName = colnames(SpecStoich)) {
  paste(SpecName,
        apply(SpecStoich, MARGIN = 1, FUN = function(X){
          X.nonzero = X[X != 0]
          X.react = names(X.nonzero)
          return(gsub(" [+] -", " - ",
                      paste(paste(X.nonzero, X.react, sep = " * "),
                            collapse = " + ")))
        }),
        sep = " = "
  )
}


#' @rdname StoichConversionFunctions
#' @keywords internal
EquationToStoich = function(SpecEquation = character(), CompName) {

  Tmp = EquationToStoichComps(SpecEquation = SpecEquation,
                              CompName = CompName)

  FormedSpecName = Tmp$SpecName
  ReactantCompNames = Tmp$SpecCompNames
  ReactantCompStoichs = Tmp$SpecCompStoichs

  ReactionStoich = StoichCompsToStoichMatrix(
    SpecCompNames = ReactantCompNames,
    SpecCompStoichs = ReactantCompStoichs,
    CompName = CompName,
    SpecName = FormedSpecName
  )

  return(list(SpecStoich = ReactionStoich,
              SpecName = FormedSpecName,
              SpecCompNames = ReactantCompNames,
              SpecCompStoichs = ReactantCompStoichs))

}


#' @rdname StoichConversionFunctions
#' @keywords internal
EquationToStoichMatrix = function(SpecEquation = character(), CompName) {

  Tmp = EquationToStoich(SpecEquation = SpecEquation, CompName = CompName)
  return(Tmp$SpecStoich)

}


#' @rdname StoichConversionFunctions
#' @keywords internal
StoichMatrixToStoichComps = function(SpecStoich = matrix(), CompName) {



  SpecCompNames = apply(
    SpecStoich,
    MARGIN = 1,
    FUN = function(X) {
      CompName[X != 0]
    },
    simplify = FALSE
  )

  SpecCompStoichs = apply(
    SpecStoich,
    MARGIN = 1,
    FUN = function(X) {
      X[X != 0]
    },
    simplify = FALSE
  )

  return(list(SpecCompNames = SpecCompNames,
              SpecCompStoichs = SpecCompStoichs))

}


#' @rdname StoichConversionFunctions
#' @keywords internal
EquationToStoichComps = function(SpecEquation = character(), CompName) {

  SpecEquation = gsub(" -([[:digit:]])", " - \\1", SpecEquation)
  FormedSpecName = trimws(gsub("=.*", "", SpecEquation))
  ReactionText = gsub("([-+*]) ","\\1",
                      gsub("^([[:digit:]])","+\\1",
                           trimws(gsub("^.*=", "", SpecEquation))))

  ReactionList = strsplit(ReactionText, split = " ")

  ReactantCompStoichs = lapply(ReactionList, FUN = function(X) {
    if (length(X) > 0) {
      as.integer(X[seq(1, length(X), by = 2)])
    } else {
      0
    }
  })
  ReactantCompNames = lapply(ReactionList, FUN = function(X) {
    if (length(X) > 0) {
      gsub("[*]", "", X[seq(2, length(X), by = 2)])
    } else {
      ""
    }
  })

  AllReactantsAreComps =
    sapply(ReactantCompNames, FUN = function(X){all(X %in% c("",CompName))})
  if(!all(AllReactantsAreComps)) {
    stop(paste0(
      "Reactants from equation(s) [",
      paste(which(!AllReactantsAreComps), collapse = ", "),
      "] are not components."
    ))
  }

  # aggregate CompNames and CompStoichs lists to eliminate duplicates
  for (i in 1:length(ReactantCompNames)) {
    tmp = stats::aggregate(x = ReactantCompStoichs[[i]],
                    by = list(Names = ReactantCompNames[[i]]),
                    FUN = sum)
    ReactantCompStoichs[[i]] = as.integer(tmp$x)
    ReactantCompNames[[i]] = tmp$Names
  }

  return(list(
    SpecName = FormedSpecName,
    SpecCompNames = ReactantCompNames,
    SpecCompStoichs = ReactantCompStoichs
  ))

}
