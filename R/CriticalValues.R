#' @name CriticalValues
#'
#' @title Edit Critical Values Table
#'
#' @description `AddCriticalValues` will add one or more rows to the critical
#'   accumulation table (CAT or CATab), while `RemoveCriticalValues` will remove
#'   one or more rows,
#'
#' @param ThisProblem A list object with a structure like that returned by
#'   `BlankProblem()`.
#' @param CATab a data.frame object with, at a minimum, columns named `CA`,
#'   `Species`, `Test.Type`, `Duration`, `Lifestage`, `Endpoint`, `Quantifier`,
#'   `References`, `Miscellaneous`. See optional parameter descriptions for
#'   further descriptions of each of those columns.
#' @param CA (optional) a numeric vector of the critical accumulation value(s)
#'   in nmol/gw.
#' @param Species (optional) a character vector of the species names to include
#'   for the corresponding `CA` value.
#' @param Test.Type (optional) a character vector of the test type (e.g.,
#'   "Acute" or "Chronic") to include for the corresponding `CA` value.
#' @param Duration (optional) a character vector of the Duration to include for
#'   the corresponding `CA` value. Can also be `"DIV=#.##"` for FAV, FCV, WQS,
#'   and HC5 critical values.
#' @param Lifestage (optional) a character vector of the organism Lifestage to
#'   include for the corresponding `CA` value. Can also be `"ACR=#.##"` for FAV,
#'   WQS, and HC5 critical values.
#' @param Endpoint (optional) a character vector of the Endpoint to include for
#'   the corresponding `CA` value. This can also be either `"FAV"`, `"FCV"`,
#'   `"HC5"`, `"WQS"`, `"CMC"`, or `"CCC"` to indicate this critical value
#'   calculates one of those water quality standards.
#' @param Quantifier (optional) a character vector of the Quantifier (e.g.,
#'   EC50, NOEC, ...) to include for the corresponding `CA` value. May also be
#'   `NA` if this is a WQS value.
#' @param References (optional) a character vector of the list of References to
#'   include for the corresponding `CA` value. Each `CA` value requires a single
#'   character string with no line breaks.
#' @param Miscellaneous (optional) a character vector of the miscellaneous
#'   information (e.g., how the value was calcualted, test conditions not
#'   covered by other columns, etc.) to include for the corresponding `CA`
#'   value.
#' @param CAToRemove an integer vector - the indices/row numbers of the critical
#'   values to remove from the table.
#'
#'
#' @return The edited version of `ThisProblem`.
#'
#' @inherit BlankProblem examples
#'
#' @family problem manipulation functions
NULL


#' @rdname CriticalValues
#' @export
AddCriticalValues = function(ThisProblem, CATab = data.frame(),
                             CA = CATab$CA,
                             Species = CATab$Species,
                             Test.Type = CATab$Test.Type,
                             Duration = CATab$Duration,
                             Lifestage = CATab$Lifestage,
                             Endpoint = CATab$Endpoint,
                             Quantifier = CATab$Quantifier,
                             References = CATab$References,
                             Miscellaneous = CATab$Miscellaneous) {

  NewProblem = ThisProblem

  if (is.null(Duration)) { Duration = NA }
  if (is.null(Lifestage)) { Lifestage = NA }
  if (is.null(Quantifier)) { Quantifier = NA }
  if (is.null(Miscellaneous)) { Miscellaneous = NA }

  NewCATab = data.frame(
    Num = NA,
    CA = CA,
    Species = Species,
    Test.Type = Test.Type,
    Duration = Duration,
    Lifestage = Lifestage,
    Endpoint = Endpoint,
    Quantifier = Quantifier,
    References = References,
    Miscellaneous = Miscellaneous
  )
  NumNewCA = nrow(NewCATab)
  NewCATab$Num = ThisProblem$N["CAT"] + seq(1, NumNewCA)
  NewProblem$N["CAT"] = NewProblem$N["CAT"] + NumNewCA

  NewProblem$CATab = rbind(
    ThisProblem$CATab,
    NewCATab
  )

  return(NewProblem)

}


#' @rdname CriticalValues
#' @export
RemoveCriticalValues = function(ThisProblem, CAToRemove) {

  NewProblem = ThisProblem

  CAToRemove = unique(CAToRemove)
  if (any(is.na(CAToRemove))) {
    stop(paste0("CA(s) (",
                paste(CAToRemoveOrig[is.na(CAToRemove)],
                      collapse = ", "),
                ") does not exist."))
  }
  if (any(CAToRemove > ThisProblem$N["CAT"])) {
    stop(paste0("There are ", ThisProblem$N["CAT"], " CAs, ",
                "trying to remove the #(",
                paste(CAToRemove[CAToRemove > ThisProblem$N["CAT"]],
                      collapse = ", "),
                ") element(s)."))
  }

  NewProblem$CATab = ThisProblem$CATab[-CAToRemove, ]
  NewProblem$N["CAT"] = ThisProblem$N["CAT"] - length(CAToRemove)
  if (NewProblem$N["CAT"] > 0) {
    NewProblem$CATab$Num = seq(1, NewProblem$N["CAT"])
  }

  return(NewProblem)

}
