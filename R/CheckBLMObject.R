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

#' @title Check an object for use in the BLMEngineInR package
#'
#' @description This function will compare an object to a reference object to
#'   make sure the required list elements are present and that they are the
#'   correct types.
#'
#' @param Object,Reference R objects that are to be compared. Assumed to be list
#'   objects.
#' @param BreakOnError A logical value indicating if an error should stop
#'   whatever function or code it might be embedded in (`TRUE`, the default) or
#'   allow it to proceed without stopping (`FALSE`).
#'
#' @return The returned value depends on the value of `BreakOnError`:
#'  \describe{
#'    \item{\code{TRUE}}{\code{TRUE} will be returned invisibly if all checks
#'      succeed, and an error with the error list as the text will be triggered
#'      if at least one check fails.}
#'    \item{\code{FALSE}}{The error list will be returned, which will be a
#'      zero-length vector if all checks succeed.}
#'  }
#'
#' @export
#'
#' @examples
#' # This one works:
#' myproblem = BlankProblem()
#' myproblem = AddMassCompartments(myproblem, data.frame(Name = "Water",Amt = 1.0, Unit = "L"))
#' CheckBLMObject(Object = myproblem, Reference = BlankProblem(), BreakOnError = FALSE)
#'
#' # This one fails:
#' myproblem$N = NULL
#' CheckBLMObject(Object = myproblem, Reference = BlankProblem(), BreakOnError = FALSE)
CheckBLMObject = function(Object, Reference, BreakOnError = TRUE) {

  ErrorList = character()

  advanced.typeof = function(X) {
    out = typeof(X)
    if (is.data.frame(X)) {out = "data.frame"}
    if (is.matrix(X)) {out = paste(out, "matrix")}
    if (is.array(X)) {out = paste(out, "array")}
    return(out)
  }

  # Level 1 compare
  Object.type = advanced.typeof(Object)
  Reference.type = advanced.typeof(Reference)
  if (Object.type != Reference.type) {
    tmp = paste0(Object.type, "-->", Reference.type)
    ErrorList = c(ErrorList, paste0("Invalid object - incorrect types (",tmp,"). "))
  }
  if (!all(names(Reference) %in% names(Object))) {
    ErrorList = c(ErrorList,
                  paste0("Invalid object - missing elements (",
                         paste(setdiff(names(Reference), names(Object)), collapse = ", "),
                         "). "))
  }

  # Level 2 compare
  for (i in names(Reference)){
    if (i %in% names(Object)) {
      if (Reference.type == "list") {
        tmp = CheckBLMObject(Object = Object[[i]],
                             Reference = Reference[[i]],
                             BreakOnError = FALSE)
      } else if (Reference.type == "data.frame") {
        tmp = CheckBLMObject(Object = Object[, i],
                             Reference = Reference[, i],
                             BreakOnError = FALSE)
      } else {tmp = character()}
      if (length(tmp) > 0) {
        tmp = gsub("[(]", paste0("(",i,":"), tmp)
        ErrorList = c(ErrorList, tmp)
      }

    }
  }

  # more specific comparisons
  if ("NObs" %in% intersect(names(Object), names(Reference))) {
    ObsVars = intersect(names(Object),
                        names(Reference)[grepl("Obs$", names(Reference)) &
                                           (names(Reference) != "NObs")])
    for (i in ObsVars) {

      # Each array and matrix need to have the correct NObs
      LenType = dim(Object[[i]])[1]
      if (is.null(LenType)) { LenType = length(Object[[i]]) }
      if (LenType != Object$NObs) {
        ErrorList = c(ErrorList, paste0("Observation counts mismatch in ", i, "."))
      }

      # Each matrix also needs to have the correct columns
      LenType = dim(Object[[i]])[2]
      if (!is.null(LenType)) {
        if (LenType != dim(Reference[[i]])[2]) {
          ErrorList = c(ErrorList, paste0("Column counts mismatch in ", i, "."))
        }
      }
    }
  }
  if ("N" %in% intersect(names(Object), names(Reference))) {
    # An array of "N"'s, each of the data.frames, matrices, and arrays
    # controlled by those "N"'s need to be the correct length
    for (iType in names(Reference$N)) {
      TypeVars = intersect(names(Object),
                           names(Reference)[grepl(paste0("^",iType), names(Reference))])
      if (iType == "BL") {
        TypeVars = setdiff(TypeVars, "BLMetal")
      }
      for (i in TypeVars) {
        LenType = dim(Object[[i]])[1]
        if (is.null(LenType)) { LenType = length(Object[[i]]) }
        if (LenType != Object$N[iType]) {
          ErrorList = c(ErrorList, paste0(iType, " counts mismatch in ", i,"."))
        }
      }

      if (iType == "Comp") {
        for (i in intersect(c("SpecStoich", "PhaseStoich"), names(Object))) {
          if (dim(Object[[i]])[2] != Object$N["Comp"]) {
            ErrorList = c(ErrorList, paste0(iType, " counts mismatch in ", i,"."))
          }
        }
      }
    }
  }

  # each of the arrays controlled by those "N"'s need to be the correct length
  NVarInBoth = intersect(c("NMass","NInLab","NInVar","NInMass","NInComp","NInDefComp","NInSpec",
                           "NDefComp","NComp","NSpec","NPhase","NBL","NMetal","NBLMetal","NCAT"),
                         intersect(names(Object), names(Reference)))
  TypeInBoth = gsub("^N", "", NVarInBoth)
  for (iType in TypeInBoth) {
    NVar = paste0("N", iType)
    TypeVars = intersect(names(Object),
                         names(Reference)[grepl(paste0("^",iType), names(Reference))])
    if (iType == "BL") {
      TypeVars = TypeVars[!grepl("BLMetal", TypeVars)]
    }
    for (i in TypeVars) {
      LenType = dim(Object[[i]])[1]
      if (is.null(LenType)) { LenType = length(Object[[i]]) }
      if ((iType %in% c("Metal", "BL")) && (Object[[NVar]] == 0)) {
        if ((typeof(Object[[i]]) == "character") & (Object[[i]] == "")){
          LenType = 0L
        }
        if ((typeof(Object[[i]]) == "integer") & (Object[[i]] == -1L)){
          LenType = 0L
        }
      }
      if (LenType != Object[[NVar]]) {
        ErrorList = c(ErrorList, paste0(iType, " counts mismatch in ", i,"."))
      }
    }
    if (iType == "Comp") {
      for (i in intersect(c("SpecStoich", "PhaseStoich"), names(Object))) {
        if (dim(Object[[i]])[2] != Object$NComp) {
          ErrorList = c(ErrorList, paste0(iType, " counts mismatch in ", i,"."))
        }
      }
    }

  }

  if (BreakOnError) {
    if (length(ErrorList) > 0) {
      stop(ErrorList)
    } else {
      return(invisible(TRUE))
    }
  } else {
    return(ErrorList)
  }

}

