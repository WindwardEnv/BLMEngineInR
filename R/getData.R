#' Get data from the input file
#'
#' `getData` reads in the input file and prepares it for input to CHESS.
#'
#' @param inputFile character(1); the path and file name to a BLM input file
#' @param NInLab integer; number of input label columns
#' @param InLabName character vector of length `NInLab`; names of input columns
#' @param NInVar integer; Number of input variables
#' @param InVarName character vector of length `NInVar`; Names of input
#'   variables
#' @param InVarMC integer vector of length `NInVar`;  Mass compartments of input
#'   variables
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
#'  \item{\code{NObs}}{integer; the number of chemistry observations}
#'  \item{\code{Labels}}{matrix with \code{NObs} rows and \code{NInLab} columns; the label columns for each observation}
#'  \item{\code{InVarObs}}{matrix with \code{NObs} rows and \code{NInVar} columns; the input variables for each observation}
#'  \item{\code{SysTemp}}{numeric vector of length \code{NObs}; input temperatures, in Celsius}
#'  \item{\code{pH}}{numeric vector (\code{NObs}); input pH for each observation}
#'  \item{\cdoe{TotConcObs}}{numeric matrix with \code{NObs} rows and \code{NComp} columns; the total concentrations of each component, including derived components}
#' }
#'
#' @keywords internal
#'
#' @noRd
getData = function(inputFile,
                   NInLab, InLabName,
                   NInVar, InVarName, InVarMC, InVarType,
                   NInComp, InCompName,
                   NComp, CompName,
                   NDefComp, DefCompName, DefCompFromNum, DefCompFromVar, DefCompSiteDens){

  stopifnot(file.exists(inputFile))

  out = ReadInputsFromFile(inputFile = inputFile,
                           NInLab = NInLab,
                           InLabName = InLabName,
                           NInVar = NInVar,
                           InVarName = InVarName,
                           NInComp = NInComp,
                           InCompName = InCompName)

  out2 = MatchInputsToProblem(NObs = out$NObs,
                              InVarObs = out$InVarObs,
                              InCompObs = out$InCompObs,
                              NInVar = NInVar,
                              InVarName = InVarName,
                              InVarMC = InVarMC,
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

  out[names(out2)] = out2

  return(out)
}



#' Read a BLM Input File
#'
#' `ReadInputsFromFile` will read a BLM input file, assuming it matches the
#' problem as defined by the input arguments.
#'
#' @param inputFile character; the path and file name to a BLM input file
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
#'  \item{\code{InLabObs}}{matrix with \code{Obs} rows and \code{InLab} columns; the input labels for each observation}
#'  \item{\code{InVarObs}}{matrix with \code{Obs} rows and \code{InVar} columns; the input variables for each observation}
#'  \item{\code{InCompObs}}{matrix with \code{Obs} rows and \code{InComp} columns; the input component concentrations for each observation}
#' }
#'
#' @export
#'
#' @examples
#' ## Not Run
#' # thisProblem = defineProblem("my_parameter_file.dat")
#' # do.call(
#' #   "ReadInputsFromFile",
#' #   args = c(thisProblem[names(thisProblem) %in% formalArgs("ReadInputsFromFile")],
#' #            list(inputFile="my_input_file.dat"))
#' # )
#' ## End Not Run
ReadInputsFromFile = function(
    inputFile,
    NInLab, InLabName,
    NInVar, InVarName,
    NInComp, InCompName){

    stopifnot(file.exists(inputFile))

    # read in input file
    # -get number of observations
    NObs = scan(file = inputFile, what = integer(), n = 1, quiet = T)

    # -get data...skip the first two header rows
    tmp = read.csv(file = inputFile, header = F, strip.white = T,
                   skip = 3, nrows = NObs)

    if (ncol(tmp) != (NInLab + NInVar + NInComp)){
      stop("Columns in input file do not match parameter file.")
    }
    rownames(tmp) = 1:NObs
    colnames(tmp) = c(InLabName, InVarName, InCompName)

    # -get labels
    InLabObs = array(NA, dim = c(NObs, NInLab),
                     dimnames = list(Obs = 1:NObs, InLabName))
    for (i in 1:NInLab){ InLabObs[,i] = as.character(tmp[,InLabName[i]]) }

    # - get input variables
    InVarObs = array(NA, dim = c(NObs, NInVar),
                     dimnames = list(Obs = 1:NObs, InVarName))
    for (i in 1:NInVar){ InVarObs[,i] = as.numeric(tmp[,InVarName[i]]) }

    # -get table of input concentrations
    InCompObs = array(NA, dim = c(NObs, NInComp),
                          dimnames = list(Obs = 1:NObs, InCompName))
    for (i in 1:NInComp){ InCompObs[,i] = as.numeric(tmp[,InCompName[i]]) }

    out = list(
      NObs = NObs,
      InLabObs = InLabObs,
      InVarObs = InVarObs,
      InCompObs = InCompObs
    )

    return(out)
}



#' Match Inputs to Problem
#'
#' `MatchInputsToProblem` will take the input variables and component
#' concentrations and match/transform them to the inputs for full list of
#' components, including defined components and WHAM components.
#'
#' @param NObs integer; the number of chemistry observations
#' @param InVarObs matrix with `NObs` rows and `NInVar` columns; the input
#'   variables for each observation
#' @param InCompObs matrix with `NObs` rows and `NInComp` columns; the input
#'   component concentrations for each observation
#' @param NInVar integer; Number of input variables
#' @param InVarName character vector of length `NInVar`; Names of input
#'   variables
#' @param InVarMC integer vector of length `NInVar`;  Mass compartments of input
#'   variables
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
#'  \item{\code{SysTemp}}{numeric vector of length \code{NObs}; input temperatures, in Celsius}
#'  \item{\code{pH}}{numeric vector (\code{NObs}); input pH for each observation}
#'  \item{\cdoe{TotConcObs}}{numeric matrix with \code{NObs} rows and \code{NComp} columns; the total concentrations of each component, including derived components}
#' }
#'
#' @keywords internal
#'
#' @noRd
MatchInputsToProblem = function(
    NObs, InVarObs, InCompObs, #inputs from file
    #information from defineProblem:
    NInVar, InVarName, InVarMC, InVarType,
    NInComp, InCompName,
    NComp, CompName,
    NDefComp, DefCompName, DefCompFromNum, DefCompFromVar, DefCompSiteDens){

  # -get table of input concentrations
  TotConcObs = matrix(numeric(), nrow = NObs, ncol = NComp,
                      dimnames = list(Obs = 1:NObs, Comp = CompName))
  TotConcObs[,InCompName] = as.matrix(InCompObs)

  # -get temperatures
  SysTemp = as.numeric(InVarObs[,InVarName[InVarType == "Temperature"]])

  # -get pH - from InVarObs or InCompObs
  if (any(InVarType == "pH")){
    pH = as.numeric(InVarObs[,InVarName[InVarType == "pH"]])
    # convert to H
    TotConcObs[,"H"] = 10^-(pH)
  } else {
    pH = -log10(InCompObs[,"H"])
  }

  # -get organic matter and parse into components
  if (any(grepl("WHAM",InVarType))){
    OM = as.matrix(InVarObs[,InVarName[InVarType %in% c("Temperature","pH") == F], drop = F])
    for (i in which(grepl("WHAM",InVarType))){
      if (InVarType[i] == "WHAM-FA"){
        if (any((InVarMC %in% InVarMC[i]) & (InVarType %in% "PercAFA"))){
          TotConcObs[,DefCompName[DefCompFromVar %in% paste0(InVarName[i],"-FA_")]] =
            (OM[,InVarName[i], drop = F] * (OM[,InVarName[(InVarMC %in% InVarMC[i]) & (InVarType %in% "PercAFA")], drop=F]/100)) %*%
            DefCompSiteDens[DefCompFromVar %in% paste0(InVarName[i],"-FA_")]
        } else {
          TotConcObs[,DefCompName[DefCompFromVar %in% paste0(InVarName[i],"-FA_")]] =
            OM[,InVarName[i], drop = F] %*%
            DefCompSiteDens[DefCompFromVar %in% paste0(InVarName[i],"-FA_")]
        }
      } else if (InVarType[i] == "WHAM-HA"){
        TotConcObs[,DefCompName[DefCompFromVar %in% paste0(InVarName[i],"-HA_")]] =
          OM[,InVarName[i], drop = F] %*%
          DefCompSiteDens[DefCompFromVar %in% paste0(InVarName[i],"-HA_")]
      } else if (InVarType[i]  == "WHAM-HAFA"){
        DOC_HA = OM[,InVarName[i], drop = F] * (OM[,InVarName[(InVarMC %in% InVarMC[i]) & (InVarType %in% "PercHA")],drop=F]/100)
        DOC_FA = OM[,InVarName[i], drop = F] * (1-OM[,InVarName[(InVarMC %in% InVarMC[i]) & (InVarType %in% "PercHA")],drop=F]/100)
        if (any((InVarMC %in% InVarMC[i]) & (InVarType %in% "PercAFA"))){
          DOC_FA = DOC_FA * (OM[,InVarName[(InVarMC %in% InVarMC[i]) & (InVarType %in% "PercAFA")], drop=F]/100)
        }
        TotConcObs[,DefCompName[DefCompFromVar %in% paste0(InVarName[i],"-HA_")]] =
          DOC_HA %*%
          DefCompSiteDens[DefCompFromVar %in% paste0(InVarName[i],"-HA_")]
        TotConcObs[,DefCompName[DefCompFromVar %in% paste0(InVarName[i],"-FA_")]] =
          DOC_FA %*%
          DefCompSiteDens[DefCompFromVar %in% paste0(InVarName[i],"-FA_")]
      }
    }
  }

  # Other Defined Components
  OtherVarDefCompName = DefCompName[DefCompFromVar %in% c("DOC-HA_","DOC-FA_",NA) == F]
  if(length(OtherVarDefCompName) > 0){
    i = which(DefCompName %in% OtherVarDefCompName)
    if (DefCompName[i] %in% CompName){
      TotConcObs[,DefCompName[i]] =
        TotConcObs[,DefCompFromVar[i], drop = F] *
        matrix(DefCompSiteDens[i], byrow = T, nrow = NObs, ncol = length(i),
               dimnames = list(1:NObs, DefCompName[i]))
    } else if (DefCompName[i] %in% InVarName) {
      TotConcObs[,DefCompName[i]] =
        InVarObs[,DefCompFromVar[i], drop = F] *
        matrix(DefCompSiteDens[i], byrow = T, nrow = NObs, ncol = length(i),
               dimnames = list(1:NObs, DefCompName[i]))
    } else {
      stop("Unknown component or variable given in Defined Components 'From' column.")
    }

  }
  NumDefCompName = DefCompName[!is.na(DefCompFromNum)]
  if (length(NumDefCompName) > 0){
    TotConcObs[,NumDefCompName] =
      matrix(DefCompFromNum[!is.na(DefCompFromNum)], byrow = T,
             nrow = NObs, ncol = length(NumDefCompName),
             dimnames = list(Obs=1:NObs, Comp = NumDefCompName))
  }

  out = list(
    SysTemp = SysTemp,
    pH = pH,
    TotConcObs = TotConcObs
  )

  return(out)

}
