#' @title Read a BLM Input File
#'
#' @description `ReadInputsFromFile` will read a BLM input file, assuming it
#' matches the problem as defined by the input arguments.
#'
#' @param InputFile character; the path and file name to a BLM input file
#' @param ThisProblem a list object following the template of BlankProblem
#' @param NInLab integer; number of input label columns
#' @param InLabName character vector of length `NInLab`; names of input columns
#' @param NInVar integer; Number of input variables
#' @param InVarName character vector of length `NInVar`; Names of input
#'   variables
#' @param NInComp integer; Number in input components
#' @param InCompName character vector of length `NInComp`; names of input
#'   components
#'
#' @return Returns a \code{list} object with the following components:
#' \describe{
#'  \item{\code{NObs}}{integer; the number of chemistry observations}
#'  \item{\code{InLabObs}}{matrix with \code{Obs} rows and \code{InLab}
#'    columns; the input labels for each observation}
#'  \item{\code{InVarObs}}{matrix with \code{Obs} rows and \code{InVar}
#'    columns; the input variables for each observation}
#'  \item{\code{InCompObs}}{matrix with \code{Obs} rows and \code{InComp}
#'    columns; the input component concentrations for each observation}
#' }
#'
#' @family BLMEngine Functions
#'
#' @export
#'
#' @examples
#' myinputfile = system.file(file.path("extdata","InputFiles","carbonate_system_test.blm4"),
#'                           package = "BLMEngineInR", mustWork = TRUE)
#' ReadInputsFromFile(InputFile = myinputfile,
#'                    ThisProblem = carbonate_system_problem)
#'
ReadInputsFromFile = function(
    InputFile,
    ThisProblem = NULL,
    NInLab = ThisProblem$N["InLab"],
    InLabName = ThisProblem$InLabName,
    NInVar = ThisProblem$N["InVar"],
    InVarName = ThisProblem$InVar$Name,
    NInComp = ThisProblem$N["InComp"],
    InCompName = ThisProblem$InCompName) {

  stopifnot(file.exists(InputFile))

  # read in input file
  # -get number of observations
  NObs = scan(file = InputFile, what = integer(), n = 1, quiet = TRUE)

  # -get data...skip the first two header rows
  Tmp = read.csv(file = InputFile, header = FALSE, strip.white = TRUE,
                 skip = 3, nrows = NObs)

  if (ncol(Tmp) != (NInLab + NInVar + NInComp)) {
    stop("Columns in input file do not match parameter file.")
  }
  rownames(Tmp) = 1:NObs
  colnames(Tmp) = c(InLabName, InVarName, InCompName)

  # -get labels
  InLabObs = array(NA, dim = c(NObs, NInLab),
                   dimnames = list(Obs = 1:NObs, InLabName))
  for (i in 1:NInLab) { InLabObs[, i] = as.character(Tmp[, InLabName[i]]) }

  # - get input variables
  InVarObs = array(NA, dim = c(NObs, NInVar),
                   dimnames = list(Obs = 1:NObs, InVarName))
  for (i in 1:NInVar) { InVarObs[, i] = as.numeric(Tmp[, InVarName[i]]) }

  # -get table of input concentrations
  InCompObs = array(NA, dim = c(NObs, NInComp),
                    dimnames = list(Obs = 1:NObs, InCompName))
  for (i in 1:NInComp) { InCompObs[, i] = as.numeric(Tmp[, InCompName[i]]) }

  Out = list(
    NObs = NObs,
    InLabObs = InLabObs,
    InVarObs = InVarObs,
    InCompObs = InCompObs
  )

  return(Out)
}



#' @title Match Inputs to Problem
#'
#' @description `MatchInputsToProblem` will take the input variables and
#'   component concentrations and match/transform them to the inputs for full
#'   list of components, including defined components and WHAM components.
#'
#' @param DFInputs A data.frame object with, at a minimum, columns named for
#'   `ThisProblem$InLabName`, `ThisProblem$InVar$Name` and
#'   `ThisProblem$InCompName`.
#' @param NObs integer; the number of chemistry observations
#' @param InLabObs character matrix with NObs rows and InLab columns; the input
#'   labels for each observation
#' @param InVarObs matrix with `NObs` rows and `NInVar` columns; the input
#'   variables for each observation
#' @param InCompObs matrix with `NObs` rows and `NInComp` columns; the input
#'   component concentrations for each observation
#' @param ThisProblem a list object following the template of BlankProblem
#' @param NInVar integer; Number of input variables
#' @param InVarName character vector of length `NInVar`; Names of input
#'   variables
#' @param InVarMCR integer vector of length `NInVar`;  Mass compartments of
#'   input variables
#' @param InVarType character vector of length `NInVar`; Types of input
#'   variables
#' @param NInComp integer; Number in input components
#' @param InCompName character vector of length `NInComp`; names of input
#'   components
#' @param NComp integer; Number of components
#' @param CompName character vector of length `NComp`; component names
#' @param NDefComp integer; Number of defined components
#' @param DefCompName character vector of length `NDefComp`; defined component
#'   names
#' @param DefCompFromNum numeric vector of length `NDefComp`; the number the
#'   defined component is formed from
#' @param DefCompFromVar character vector of length `NDefComp`; the column used
#'   to form the defined component
#' @param DefCompSiteDens numeric vector of length `NDefComp`; the binding site
#'   density of each defined component
#'
#' @return Returns a \code{list} object with the following components:
#' \describe{
#'   \item{NObs}{integer; the number of chemistry observations}
#'   \item{InLabObs}{matrix with NObs rows and InLab columns; the input labels
#'     for each observation}
#'   \item{InVarObs}{matrix with NObs rows and InVar columns; the input
#'     variables for each observation}
#'   \item{InCompObs}{matrix with NObs rows and InComp columns; the input
#'     component concentrations for each observation}
#'   \item{\code{SysTempCelsiusObs}}{numeric vector of length \code{NObs}; input
#'    temperatures, in Celsius}
#'   \item{\code{SysTempKelvinObs}}{numeric vector of length \code{NObs}; input
#'    temperatures, in Kelvin}
#'  \item{\code{pHObs}}{numeric vector (\code{NObs}); input pH for each
#'    observation}
#'  \item{\code{TotConcObs}}{numeric matrix with \code{NObs} rows and
#'    \code{NComp} columns; the total concentrations of each component,
#'    including derived components}
#'  \item{\code{HumicSubstGramsPerLiterObs}}{numeric matrix with \code{NObs}
#'    rows and 2 columns; the total concentration of humic substances (humic/HA
#'    and fulvic/FA) in grams per liter}
#' }
#'
#' @family BLMEngine Functions
#'
#' @export
#'
MatchInputsToProblem = function(
    #inputs from file
    DFInputs = data.frame(),
    NObs = nrow(DFInputs),
    InLabObs = DFInputs[, ThisProblem$InLabName, drop = FALSE],
    InVarObs = DFInputs[, ThisProblem$InVar$Name, drop = FALSE],
    InCompObs = DFInputs[, ThisProblem$InCompName, drop = FALSE],
    #information from DefineProblem:
    ThisProblem = NULL,
    NInVar = ThisProblem$N["InVar"],
    InVarName = ThisProblem$InVar$Name,
    InVarMCR = ThisProblem$InVar$MCR,
    InVarType = ThisProblem$InVar$Type,
    NInComp = ThisProblem$N["InComp"],
    InCompName = ThisProblem$InCompName,
    NComp = ThisProblem$N["Comp"],
    CompName = ThisProblem$Comp$Name,
    NDefComp = ThisProblem$N["DefComp"],
    DefCompName = ThisProblem$DefComp$Name,
    DefCompFromNum = ThisProblem$DefComp$FromNum,
    DefCompFromVar = ThisProblem$DefComp$FromVar,
    DefCompSiteDens = ThisProblem$DefComp$SiteDens) {

  # -get table of input concentrations
  TotConcObs = matrix(numeric(), nrow = NObs, ncol = NComp,
                      dimnames = list(Obs = 1:NObs, Comp = CompName))
  TotConcObs[, InCompName] = as.matrix(InCompObs)

  # -get temperatures
  SysTempCelsiusObs = as.numeric(InVarObs[, InVarName[InVarType == "Temperature"]])

  # -get pHObs - from InVarObs or InCompObs
  if (any(InVarType == "pH")) {
    pHObs = as.numeric(InVarObs[, InVarName[InVarType == "pH"]]) # nolint: object_name_linter, line_length_linter.
  } else {
    pHObs = -log10(InCompObs[, "H"]) # nolint: object_name_linter, line_length_linter.
  }

  # -get organic matter and parse into components
  HumicSubstGramsPerLiterObs = array(0.0, dim = c(NObs, 2),
                   dimnames = list(Obs = 1:NObs, c("HA", "FA")))
  if (any(grepl("WHAM", InVarType))) {
    OM = as.matrix(InVarObs[, InVarName[InVarType %in% c("Temperature", "pH") == FALSE], drop = FALSE]) # nolint: line_length_linter.
    for (i in which(grepl("WHAM", InVarType))) {
      OMColI = OM[, InVarName[i], drop = FALSE]#mg C / L

      # Initialize FracAFA and FracHA with the needed values
      FracAFACol = matrix(1, nrow = NObs, ncol = 1)
      FracHACol = matrix(NA, nrow = NObs, ncol = 1)
      if (any((InVarMCR %in% InVarMCR[i]) & (InVarType %in% "PercAFA"))) {
        FracAFACol = OM[, InVarName[(InVarMCR %in% InVarMCR[i]) &
                                      (InVarType %in% "PercAFA")],
                        drop = FALSE] / 100
      }
      if (any((InVarMCR %in% InVarMCR[i]) & (InVarType %in% "PercHA"))) {
        FracHACol = OM[, InVarName[(InVarMCR %in% InVarMCR[i]) &
                                     (InVarType %in% "PercHA")],
                       drop = FALSE] / 100
      }

      # Pull out relevant HA and FA components
      FAComps = which(DefCompFromVar %in% paste0(InVarName[i], "-FA_"))
      FACompName = DefCompName[FAComps]
      FACompSiteDens = DefCompSiteDens[FAComps] #mol / mg C
      HAComps = which(DefCompFromVar %in% paste0(InVarName[i], "-HA_"))
      HACompName = DefCompName[HAComps]
      HACompSiteDens = DefCompSiteDens[HAComps] #mol / mg C

      HumicSubstGramsPerLiterObs = array(0.0, dim = c(NObs, 2),
                       dimnames = list(Obs = 1:NObs, c("HA", "FA")))

      # Calculate the component concentrations
      if (InVarType[i] == "WHAM-FA") {
        # This is a FA-only input variable
        TotConcObs[, FACompName] = (OMColI * FracAFACol) %*% FACompSiteDens
        # TotConcObs[, "DonnanFA"] = 1E-4
        HumicSubstGramsPerLiterObs[, "FA"] = (OMColI * FracAFACol) / 1000 * 2
      } else if (InVarType[i] == "WHAM-HA") {
        # This is a HA-only input variable
        TotConcObs[, HACompName] = OMColI %*% HACompSiteDens
        # TotConcObs[, "DonnanHA"] = 1E-4
        HumicSubstGramsPerLiterObs[, "HA"] = OMColI / 1000 * 2
      } else if (InVarType[i]  == "WHAM-HAFA") {
        # This is a combined HA + FA component, and % HA is needed
        TotConcObs[, HACompName] = (OMColI * FracHACol) %*% HACompSiteDens
        TotConcObs[, FACompName] =
          (OMColI * (1 - FracHACol) * FracAFACol) %*% FACompSiteDens
        # TotConcObs[, c("DonnanHA", "DonnanFA")] = 1E-4
        HumicSubstGramsPerLiterObs[, "FA"] = (OMColI * (1 - FracHACol) * FracAFACol) / 1000 * 2
        HumicSubstGramsPerLiterObs[, "HA"] = (OMColI * FracHACol) / 1000 * 2
      }
      # TotConcObs[, HACompName] = SolHAObs %*% HACompSiteDens
      # TotConcObs[, FACompName] = SolFAObs %*% FACompSiteDens

      # # Calculate the moles of each OM component in solution...this doesn't
      # # # make sense...why would we multiply by the site density AGAIN?
      # # SolHAObs = TotConcObs[, HACompName] %*% HACompSiteDens
      # # SolFAObs = TotConcObs[, FACompName] %*% FACompSiteDens
      # HumicSubstGramsPerLiterObs = array(c(rowSums(TotConcObs[, HACompName, drop = FALSE]),
      #                    rowSums(TotConcObs[, FACompName, drop = FALSE])),
      #                  dim = c(NObs, 2),
      #                  dimnames = list(Obs = 1:NObs, c("HA", "FA")))

    }
  }

  # Other Defined Components - based on variables
  VarDefComps = which(DefCompFromVar %in% c(
    paste0(InVarName[InVarType %in% c("WHAM-HA","WHAM-FA","WHAM-HAFA")],
           "-", c("HA","FA"),"_"), NA) == FALSE)
  if (length(VarDefComps) > 0) {
    VarDefCompName = DefCompName[VarDefComps]
    for (i in which(DefCompName %in% VarDefCompName)){
      if (DefCompFromVar[i] == "pH") {
        TotConcObs[, DefCompName[i]] = 10^-pHObs
      } else if (DefCompFromVar[i] %in% CompName) {
        TotConcObs[, DefCompName[i]] =
          TotConcObs[, DefCompFromVar[i], drop = FALSE] *
          matrix(DefCompSiteDens[i], byrow = TRUE, nrow = NObs, ncol = 1,
                 dimnames = list(1:NObs, DefCompName[i]))
      } else if (DefCompFromVar[i] %in% InVarName) {
        TotConcObs[, DefCompName[i]] =
          as.matrix(InVarObs[, DefCompFromVar[i], drop = FALSE]) *
          matrix(DefCompSiteDens[i], byrow = TRUE, nrow = NObs, ncol = 1,
                 dimnames = list(1:NObs, DefCompName[i]))
      } else {
        stop("Unknown component or variable given in Defined Components 'From' column.") # nolint: line_length_linter.
      }
    }
  }

  # Other Defined Components - based on numbers
  NumDefComps = which(!is.na(DefCompFromNum))
  if (length(NumDefComps) > 0) {
    NumDefCompName = DefCompName[NumDefComps]
    TotConcObs[, NumDefCompName] =
      matrix(DefCompFromNum[NumDefComps] * DefCompSiteDens[NumDefComps],
             byrow = TRUE, nrow = NObs, ncol = length(NumDefCompName),
             dimnames = list(Obs = 1:NObs, Comp = NumDefCompName))
  }

  Out = list(
    NObs = as.integer(NObs),
    InLabObs = as.matrix(InLabObs),
    InVarObs = as.matrix(InVarObs),
    InCompObs = as.matrix(InCompObs),
    SysTempCelsiusObs = SysTempCelsiusObs,
    SysTempKelvinObs = SysTempCelsiusObs + 273.15,
    pHObs = pHObs,
    TotConcObs = TotConcObs,
    HumicSubstGramsPerLiterObs = HumicSubstGramsPerLiterObs
  )

  return(Out)

}


#' @title Get data from the input file
#'
#' @description `GetData` reads in the input file and prepares it for input to
#'   the BLM function.
#'
#' @param InputFile character(1); the path and file name to a BLM input file
#' @param ThisProblem a list object following the template of BlankProblem
#' @param NInLab integer; number of input label columns
#' @param InLabName character vector of length `NInLab`; names of input columns
#' @param NInVar integer; Number of input variables
#' @param InVarName character vector of length `NInVar`; Names of input
#'   variables
#' @param InVarMCR integer vector of length `NInVar`;  Mass compartments of
#'   input variables
#' @param InVarType character vector of length `NInVar`; Types of input
#'   variables
#' @param NInComp integer; Number in input components
#' @param InCompName character vector of length `NInComp`; names of input
#'   components
#' @param NComp integer; Number of components
#' @param CompName character vector of length `NComp`; component names
#' @param NDefComp integer; Number of defined components
#' @param DefCompName character vector of length `NDefComp`; defined component
#'   names
#' @param DefCompFromNum numeric vector of length `NDefComp`; the number the
#'   defined component is formed from
#' @param DefCompFromVar character vector of length `NDefComp`; the column used
#'   to form the defined component
#' @param DefCompSiteDens numeric vector of length `NDefComp`; the binding site
#'   density of each defined component
#'
#' @return Returns a `list` object with the following components:
#' \describe{
#'   \item{NObs}{integer; the number of chemistry observations}
#'   \item{InLabObs}{matrix with NObs rows and InLab columns; the input labels
#'     for each observation}
#'   \item{InVarObs}{matrix with NObs rows and InVar columns; the input
#'     variables for each observation}
#'   \item{InCompObs}{matrix with NObs rows and InComp columns; the input
#'     component concentrations for each observation}
#'   \item{\code{SysTempCelsiusObs}}{numeric vector of length \code{NObs}; input
#'    temperatures, in Celsius}
#'   \item{\code{SysTempKelvinObs}}{numeric vector of length \code{NObs}; input
#'    temperatures, in Kelvin}
#'  \item{\code{pHObs}}{numeric vector (\code{NObs}); input pH for each
#'    observation}
#'  \item{\code{TotConcObs}}{numeric matrix with \code{NObs} rows and
#'    \code{NComp} columns; the total concentrations of each component,
#'    including derived components}
#'  \item{\code{HumicSubstGramsPerLiterObs}}{numeric matrix with \code{NObs}
#'    rows and 2 columns; the total concentration of humic substances (humic/HA
#'    and fulvic/FA) in grams per liter}
#'  }
#'
#' @family BLMEngine Functions
#'
#' @export
#'
#' @examples
#' myinputfile = system.file(file.path("extdata","InputFiles","carbonate_system_test.blm4"),
#'                           package = "BLMEngineInR", mustWork = TRUE)
#' GetData(InputFile = myinputfile, ThisProblem = carbonate_system_problem)#'
#'
GetData = function(InputFile,
                   ThisProblem = NULL,
                   NInLab = ThisProblem$N["InLab"],
                   InLabName = ThisProblem$InLabName,
                   NInVar = ThisProblem$N["InVar"],
                   InVarName = ThisProblem$InVar$Name,
                   InVarMCR = ThisProblem$InVar$MCR,
                   InVarType = ThisProblem$InVar$Type,
                   NInComp = ThisProblem$N["InComp"],
                   InCompName = ThisProblem$InCompName,
                   NComp = ThisProblem$N["Comp"],
                   CompName = ThisProblem$Comp$Name,
                   NDefComp = ThisProblem$N["DefComp"],
                   DefCompName = ThisProblem$DefComp$Name,
                   DefCompFromNum = ThisProblem$DefComp$FromNum,
                   DefCompFromVar = ThisProblem$DefComp$FromVar,
                   DefCompSiteDens = ThisProblem$DefComp$SiteDens) {

  stopifnot(file.exists(InputFile))

  Out = ReadInputsFromFile(InputFile = InputFile,
                           NInLab = NInLab,
                           InLabName = InLabName,
                           NInVar = NInVar,
                           InVarName = InVarName,
                           NInComp = NInComp,
                           InCompName = InCompName)

  Out2 = MatchInputsToProblem(NObs = Out$NObs,
                              InLabObs = Out$InLabObs,
                              InVarObs = Out$InVarObs,
                              InCompObs = Out$InCompObs,
                              NInVar = NInVar,
                              InVarName = InVarName,
                              InVarMCR = InVarMCR,
                              InVarType = InVarType,
                              NInComp = NInComp,
                              InCompName = InCompName,
                              NComp = NComp,
                              CompName = CompName,
                              NDefComp = NDefComp,
                              DefCompName = DefCompName,
                              DefCompFromNum = DefCompFromNum,
                              DefCompFromVar = DefCompFromVar,
                              DefCompSiteDens = DefCompSiteDens)

  return(Out2)
}

#' @title Make a blank inputs list object
#'
#' @description This function is internal because the inputs can be specified in
#'   a data.frame object then matched to ThisProblem with
#'   `MatchInputsToProblem(DFInputs = , ThisProblem = )`. This function is useful
#'   as a comparison to make sure the outputs of `GetData` or
#'   `MatchInputsToProblem` and the inputs of `BLM` are acceptable.
#'
#' @param ThisProblem A list object with a structure like that returned by
#'   `BlankProblem()`.
#'
#' @return A list object with a template for defining the inputs for the
#'   chemical problem for the `BLMEngineInR` functions. Each element in the list
#'   is a vector, list, or data.frame object grouping related parameters
#'   together. See `str(BlankProblem())` for the structure and names of the list
#'   object.
#'
#' @keywords internal
BlankInputList = function(ThisProblem, NObs = 0L) {

  Out = list(
    NObs = as.integer(NObs),
    InLabObs = array("", dim = c(NObs, ThisProblem$N["InLab"]),
                     dimnames = list(NULL, ThisProblem$InLabName)),
    InVarObs = array(0.0, dim = c(NObs, ThisProblem$N["InVar"]),
                     dimnames = list(NULL, ThisProblem$InVar$Name)),
    InCompObs = array(0.0, dim = c(NObs, ThisProblem$N["InComp"]),
                      dimnames = list(NULL, ThisProblem$InCompName)),
    SysTempCelsiusObs = numeric(NObs),
    SysTempKelvinObs = numeric(NObs),
    pHObs = numeric(NObs),
    TotConcObs = array(0.0, dim = c(NObs, ThisProblem$N["Comp"]),
                       dimnames = list(NULL, ThisProblem$Comp$Name)),
    HumicSubstGramsPerLiterObs = array(0.0, dim = c(NObs,2),
                                       dimnames = list(NULL,c("HA","FA")))
  )

  return(Out)

}

