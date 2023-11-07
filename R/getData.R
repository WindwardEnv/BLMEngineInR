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
#'  \item{\code{NObs}}{integer; the number of species}
#'  \item{\code{Labels}}{matrix with \code{NObs} rows and \code{NInLab} columns; the label columns for each observation}
#'  \item{\code{InVar}}{matrix with \code{NObs} rows and \code{NInVar} columns; the input variables for each observation}
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

  # read in input file
  # -get number of observations
  NObs = scan(file = inputFile, what = integer(), n = 1, quiet = T)
  out = list(NObs = NObs)

  # -get data...skip the first two header rows
  tmp = read.csv(file = inputFile, header = F, strip.white = T, skip = 3,
                 nrows = NObs)
  if (ncol(tmp) != (NInLab + NInVar + NInComp)){
    stop("Columns in input file do not match parameter file.")
  }
  colnames(tmp) = c(InLabName, InVarName, InCompName)

  # -get labels
  out$Labels = array(tmp[,InLabName], dim = c(NObs, NInLab),
                     dimnames = list(Obs=1:NObs, Label=InLabName))

  # - get input variables
  out$InVar = tmp[,InVarName]

  # -get table of input concentrations
  TotConcObs = matrix(numeric(),nrow = NObs, ncol = NComp,
                      dimnames = list(Obs = 1:NObs, Comp = CompName))
  TotConcObs[,InCompName] = as.matrix(tmp[,InCompName])

  # -get temperatures
  out$SysTemp = as.numeric(tmp[,InVarName[InVarType == "Temperature"]])

  # -get pH
  if (any(InVarType == "pH")){
    out$pH = as.numeric(tmp[,InVarName[InVarType == "pH"]])
    # convert to H
    TotConcObs[,"H"] = 10^-(out$pH)
  } else {
    out$pH = -log10(TotConcObs[,"H"])
  }

  # -get organic matter and parse into components
  if (any(grepl("WHAM",InVarType))){
    OM = as.matrix(tmp[,InVarName[InVarType %in% c("Temperature","pH") == F], drop = F])
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
    TotConcObs[,DefCompName[i]] =
      TotConcObs[,DefCompFromVar[i], drop = F] *
      matrix(DefCompSiteDens[i], byrow = T, nrow = NObs, ncol = length(i),
             dimnames = list(1:NObs, DefCompName[i]))
  }
  NumDefCompName = DefCompName[!is.na(DefCompFromNum)]
  if (length(NumDefCompName) > 0){
    TotConcObs[,NumDefCompName] =
      matrix(DefCompFromNum[!is.na(DefCompFromNum)], byrow = T,
             nrow = NObs, ncol = length(NumDefCompName),
             dimnames = list(Obs=1:NObs, Comp = NumDefCompName))
  }

  out$TotConcObs = TotConcObs

  return(out)
}
