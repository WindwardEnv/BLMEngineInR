#' @name ProblemConversionFunctions
#'
#' @title Problem Conversion functions
#'
#' @description These functions are for converting between a data.frame-based
#'   list object and a pure list-based object. The data.frame-based list object
#'   is easier to edit, generally, while the pure list-based object is more
#'   readily used by C++ functions such as CHESS.
#'
#' @param ThisProblemDF a list object with data.frames, lists, and vectors,
#'   following the template given by `BlankProblem()`.
#' @param ThisProblemList a list object where each element is an input for
#'   CHESS, following the template given by `BlankProblemList()`.
#'
#' @return a list object of the opposite type
#'


#' @rdname ProblemConversionFunctions
#' @keywords internal
ConvertToList = function(ThisProblemDF) {

  CheckBLMObject(ThisProblemDF, BlankProblem(), BreakOnError = TRUE)
  ThisProblemList = list()

  CompositeNames = c("N", "Mass", "InLab", "InVar", "InComp", "DefComp", "Comp",
                     "Spec", "Phase", "BL", "Metal", "BLMetal")
  OrganizedNames = c("Index", "WHAM")
  AsIsNames = setdiff(names(ThisProblemDF), c(CompositeNames, OrganizedNames))
  NoZeroLengthNames = c("MetalName","MetalCompR")

  for (i in CompositeNames) {
    if (is.data.frame(ThisProblemDF[[i]])) {
      for (j in names(ThisProblemDF[[i]])) {
        ThisProblemList[[paste0(i, j)]] = ThisProblemDF[[i]][, j]
        if (any(names(ThisProblemDF[[i]]) %in% "Name") & (j != "Name")) {
          names(ThisProblemList[[paste0(i, j)]]) = ThisProblemDF[[i]]$Name
        }
      }
    } else if (is.vector(ThisProblemDF[[i]]) || is.list(ThisProblemDF[[i]])) {
      for (j in names(ThisProblemDF[[i]])) {
        ThisProblemList[[paste0(i, j)]] = unname(ThisProblemDF[[i]][j])
      }
    }
  }

  for (i in OrganizedNames) {
    ThisProblemList = c(ThisProblemList, ThisProblemDF[[i]])
  }

  for (i in AsIsNames) {
    ThisProblemList[[i]] = ThisProblemDF[[i]]
  }

  for (i in NoZeroLengthNames) {
    if (typeof(ThisProblemList[[i]]) == "character"){
      ThisProblemList[[i]] = ""
    } else if (typeof(ThisProblemList[[i]]) == "integer") {
      ThisProblemList[[i]] = -1L
    } else {
      ThisProblemList[[i]] = NA
    }
  }

  # Re-organize list elements so they're grouped nicer
  BlankProblemNames = names(BlankProblemList())
  ThisProblemList = c(
    ThisProblemList[BlankProblemNames],
    ThisProblemList[setdiff(names(ThisProblemList), BlankProblemNames)]
  )

  CheckBLMObject(ThisProblemList, BlankProblemList(), BreakOnError = TRUE)
  return(ThisProblemList)

}


#' @rdname ProblemConversionFunctions
#' @keywords internal
ConvertToDF = function(ThisProblemList) {

  CheckBLMObject(ThisProblemList, BlankProblemList(), BreakOnError = TRUE)
  ThisProblemDF = list()

  ListNames = names(ThisProblemList)
  CompositeNames = c("N", "Mass", "InLab", "InVar", "InComp", "DefComp", "Comp",
                     "Spec", "Phase", "BL", "Metal", "BLMetal")
  OrganizedNames = list(
    Index = c("AqueousMCR", "BioticLigMCR", "WHAMDonnanMCR"),
    WHAM = c("DoWHAM", "WHAMVer", "WdatFile", "wDLF", "wKZED",
             "wP", "wRadius", "wMolWt")
  )
  NoZeroLengthNames = c("MetalName","MetalCompR")

  for (i in NoZeroLengthNames) {
    if (typeof(ThisProblemList[[i]]) == "character"){
      ThisProblemList[[i]] = character()
    } else if (typeof(ThisProblemList[[i]]) == "integer") {
      ThisProblemList[[i]] = integer()
    } else if (typeof(ThisProblemList[[i]]) == "double") {
      ThisProblemList[[i]] = double()
    } else {
      ThisProblemList[[i]] = NULL
    }
  }

  ConvertedNames = c()

  MatchedByMultipleCompositeNames = ListNames[rowSums(sapply(
    CompositeNames,
    FUN = function(X){
      grepl(paste0("^", X), ListNames)
    })) > 1]

  for (i in CompositeNames) {

    MatchedByI = ListNames[grepl(paste0("^", i), ListNames)]
    MatchedByINotMultiple = setdiff(MatchedByI, MatchedByMultipleCompositeNames)
    if (length(MatchedByINotMultiple) > 0) {
      MemberNames = MatchedByINotMultiple
    } else {
      MemberNames = MatchedByI
    }
    VectorNames = MemberNames[sapply(ThisProblemList[MemberNames], FUN = is.vector)]
    OtherNames = setdiff(MemberNames, VectorNames)

    if (length(VectorNames) == 1) {
      ThisProblemDF[[VectorNames]] = ThisProblemList[[VectorNames]]
    } else {
      # if (all(sapply(ThisProblemList[VectorNames], FUN = length) == 1)) {
      if (i %in% "N") {
        ThisProblemDF[[i]] = unlist(ThisProblemList[VectorNames])
      } else {
        ThisProblemDF[[i]] = as.data.frame(ThisProblemList[VectorNames])
        rownames(ThisProblemDF[[i]]) = NULL
      }
      names(ThisProblemDF[[i]]) = gsub(paste0("^", i), "", VectorNames)
    }

    for (j in OtherNames) {
      ThisProblemDF[[j]] = ThisProblemList[[j]]
    }

    ConvertedNames = c(ConvertedNames, MemberNames)

  }

  for (i in names(OrganizedNames)) {
    ThisProblemDF[[i]] = list()
    ThisProblemDF[[i]][OrganizedNames[[i]]] = ThisProblemList[OrganizedNames[[i]]]
    ConvertedNames = c(ConvertedNames, OrganizedNames[[i]])
  }


  AsIsNames = setdiff(ListNames, ConvertedNames)
  for (i in AsIsNames) {
    ThisProblemDF[[i]] = ThisProblemList[[i]]
  }

  ThisProblemDF$Spec$Equation = StoichMatrixToEquation(
    ThisProblemList$SpecStoich,
    ThisProblemList$SpecName,
    ThisProblemList$CompName)

  # Re-organize list elements so they're grouped nicer
  BlankProblemNames = names(BlankProblem())
  ThisProblemDF = c(
    ThisProblemDF[BlankProblemNames],
    ThisProblemDF[setdiff(names(ThisProblemDF), BlankProblemNames)]
  )

  CheckBLMObject(ThisProblemDF, BlankProblem(), BreakOnError = TRUE)
  return(ThisProblemDF)

}
