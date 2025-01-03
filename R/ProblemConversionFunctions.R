# Copyright 2024 Windward Environmental LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
ConvertToList = function(ThisProblemDF) { #nolint: cyclocomp_linter

  CheckBLMObject(ThisProblemDF, BlankProblem(), BreakOnError = TRUE)
  ThisProblemList = list()

  CompositeNames = c("N", "Mass", "InLab", "InVar", "InComp", "DefComp", "Comp",
                     "Spec", "Phase", "BL", "Metal", "BLMetal", "WHAM")
  OrganizedNames = c("Index")
  AsIsNames = setdiff(names(ThisProblemDF), c(CompositeNames, OrganizedNames))
  NoZeroLengthNames = c("MetalName", "MetalCompR", "BLName", "BLCompR")

  for (i in CompositeNames) {
    if (is.data.frame(ThisProblemDF[[i]])) {
      for (j in names(ThisProblemDF[[i]])) {
        ThisProblemList[[paste0(i, j)]] = ThisProblemDF[[i]][, j]
        if (any(names(ThisProblemDF[[i]]) %in% "Name") && (j != "Name")) {
          names(ThisProblemList[[paste0(i, j)]]) = ThisProblemDF[[i]]$Name
        }
      }
    } else if (is.list(ThisProblemDF[[i]])) {
      for (j in names(ThisProblemDF[[i]])) {
        ThisProblemList[[paste0(i, j)]] = ThisProblemDF[[i]][[j]]
      }
    } else if (is.vector(ThisProblemDF[[i]])) {
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
    if (length(ThisProblemList[[i]]) == 0) {
      if (typeof(ThisProblemList[[i]]) == "character") {
        ThisProblemList[[i]] = ""
      } else if (typeof(ThisProblemList[[i]]) == "integer") {
        ThisProblemList[[i]] = -1L
      } else {
        ThisProblemList[[i]] = NA
      }
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
                     "Spec", "Phase", "BL", "Metal", "BLMetal", "WHAM")
  OrganizedNames = list(
    Index = c("AqueousMCR", "BioticLigMCR", "WHAMDonnanMCR")
  )
  NoZeroLengthNames = c("MetalName", "MetalCompR", "BLName", "BLCompR")

  for (i in NoZeroLengthNames) {
    if (typeof(ThisProblemList[[i]]) == "character") {
      if (ThisProblemList[[i]] == "") {
        ThisProblemList[[i]] = character()
      }
    } else if (typeof(ThisProblemList[[i]]) == "integer") {
      if (ThisProblemList[[i]] <= 0L) {
        ThisProblemList[[i]] = integer()
      }
    } else if (typeof(ThisProblemList[[i]]) == "logical") {
      if (is.na(ThisProblemList[[i]])) {
        ThisProblemList[[i]] = logical()
      }
    } else {
      ThisProblemList[[i]] = NULL
    }
  }

  ConvertedNames = c()

  MatchedByMultCompositeNames = ListNames[
    rowSums(sapply(
      CompositeNames,
      FUN = function(X) {
        grepl(paste0("^", X), ListNames)
      }
    )) > 1
  ]

  for (i in CompositeNames) {

    MatchedByI = ListNames[grepl(paste0("^", i), ListNames)]
    MatchedByINotMultiple = setdiff(MatchedByI, MatchedByMultCompositeNames)
    if (length(MatchedByINotMultiple) > 0) {
      MemberNames = MatchedByINotMultiple
    } else {
      MemberNames = MatchedByI
    }
    VectorNames =
      MemberNames[sapply(ThisProblemList[MemberNames], FUN = is.vector)]
    OtherNames = setdiff(MemberNames, VectorNames)

    if (length(VectorNames) == 1) {
      ThisProblemDF[[VectorNames]] = ThisProblemList[[VectorNames]]
    } else {
      # if (all(sapply(ThisProblemList[VectorNames], FUN = length) == 1)) {
      if (i %in% "N") {
        ThisProblemDF[[i]] = unlist(ThisProblemList[VectorNames])
      } else if (i %in% "WHAM") {
        VectorNames = paste0("WHAM", names(BlankWHAM()))
        OtherNames =
          OtherNames[OtherNames %in% c(VectorNames, "WHAMDonnanMCR") == FALSE]
        ThisProblemDF[[i]] = as.list(ThisProblemList[VectorNames])
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
    ThisProblemDF[[i]][OrganizedNames[[i]]] =
      ThisProblemList[OrganizedNames[[i]]]
    ConvertedNames = c(ConvertedNames, OrganizedNames[[i]])
  }


  AsIsNames = setdiff(ListNames, ConvertedNames)
  for (i in AsIsNames) {
    ThisProblemDF[[i]] = ThisProblemList[[i]]
  }

  ThisProblemDF$Spec$Equation = StoichMatrixToEquation(
    ThisProblemList$SpecStoich,
    ThisProblemList$SpecName,
    ThisProblemList$CompName
  )

  # Re-organize list elements so they're grouped nicer
  BlankProblemNames = names(BlankProblem())
  ThisProblemDF = c(
    ThisProblemDF[BlankProblemNames],
    ThisProblemDF[setdiff(names(ThisProblemDF), BlankProblemNames)]
  )

  CheckBLMObject(ThisProblemDF, BlankProblem(), BreakOnError = TRUE)
  return(ThisProblemDF)

}
