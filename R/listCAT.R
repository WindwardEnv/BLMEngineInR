#' List Critical Accumulation Table
#'
#' List out the critical accumulation table for the user to allow them to pick
#' which CAT number they should specify for a toxicity run where the critical
#' value is coming from the table in the parameter file.
#'
#' @param paramFile the file name and path of the parameter file.
#'
#' @return A data frame with species and test information, and the CAT numbers
#'   that would correspond to each
#' @export
#'
#' @examples
#' ## Not Run
#' # listCAT
#' ## End Not Run
listCAT = function(paramFile){
  # Read in CAT table from parameter file
  # Return that table for the user
  out = data.frame(
    CAT = 1:10,
    Species = LETTERS[1:10],
    Lifestage = letters[1:10],
    Endpoint = NA,
    EffectLevel = NA,
    Duration = NA,
    References = NA,
    Misc = NA
  )

  return(out)
}
