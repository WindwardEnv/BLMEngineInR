#' @name InLabs
#'
#' @title Add or remove input labels in a problem
#'
#' @param ThisProblem A list object with a structure like that returned by
#'   `BlankProblem()`.
#' @param InLabName A character vector with the name(s) of the new input
#'   label(s).
#' @param InLabToRemove A character vector with names or indices of the input
#'   label(s) to remove from `ThisProblem`.
#'
#' @return `ThisProblem`, with the edited input labels.
#'
#' @inherit BlankProblem examples
#'
#' @family problem manipulation functions
NULL

#' @rdname InLabs
#' @export
AddInLabs = function(ThisProblem, InLabName) {

  NewProblem = ThisProblem

  # error checking
  if (any((InLabName %in% ThisProblem$InLabName))) {
    stop(paste0(
      "Input Label(s) (",
      paste(paste0("\"", InLabName[InLabName %in% ThisProblem$InLabName], "\""), collapse = ", "),
      ") already exist."
    ))
  }
  if (any(is.na(InLabName))) {
    stop("NA arguments not allowed.")
  }

  # add input labels
  NewProblem$InLabName = c(ThisProblem$InLabName,
                           trimws(as.character(InLabName)))
  NewProblem$N["InLab"] = ThisProblem$N["InLab"] + length(InLabName)

  return(NewProblem)

}

#' @rdname InLabs
#' @export
RemoveInLabs = function(ThisProblem, InLabToRemove) {

  NewProblem = ThisProblem

  InLabToRemoveOrig = InLabToRemove
  if (is.character(InLabToRemove)) {
    InLabToRemove = which(ThisProblem$InLabName %in% InLabToRemove)
  }
  if (length(InLabToRemove) < 1) {
    stop(paste0("Mass Compartment \"", InLabToRemoveOrig, "\" does not exist."))
  }
  if (any(InLabToRemove > ThisProblem$N["InLab"])) {
    stop(paste0("There are ", ThisProblem$N["InLab"], " input labels, ",
                "trying to remove the #(",
                InLabToRemove[InLabToRemove > ThisProblem$N["InLab"]],
                ") element(s)."))
  }

  # Remove input labels
  NewProblem$InLabName = NewProblem$InLabName[-InLabToRemove]
  NewProblem$N["InLab"] = ThisProblem$N["InLab"] - length(InLabToRemove)

  return(NewProblem)

}

