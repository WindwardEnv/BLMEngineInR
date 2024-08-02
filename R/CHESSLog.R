#' Create a file that shows the problem in a different, human-friendly format
#'
#' @param ThisProblem A problem list object, such as returned by
#'   `DefineProblem`.
#' @param LogFilename The path and file name of the created log file. By
#'   default, it is "CHESSLOG.txt", placed in the same directory as ParamFile.
#'
#' @return invisibly returns TRUE
#'
#' @keywords internal
CHESSLog = function(ThisProblem,
                    LogFilename = file.path(dirname(ThisProblem$ParamFile), "CHESSLOG.txt")) {

  # initialize
  CompName = character()
  SpecName = character()
  SpecStoich = matrix()
  SpecLogK = numeric()
  SpecDeltaH = numeric()
  CompType = character()
  SpecCharge = integer()

  # unpack the input list
  ThisProblemList = ConvertToList(ThisProblem)
  for (i in names(ThisProblemList)) {
    assign(i, ThisProblemList[[i]])
  }

  # initialize log file
  write(paste0("CHESS problem defined by '", ThisProblem$ParamFile, "' parameter file:\n",
               "(",Sys.time(),")"),
        file = LogFilename, append = FALSE)

  # Component List
  write(strrep("-", 80), file = LogFilename, append = TRUE)
  write("Component List:", file = LogFilename, append = TRUE)
  write(CompName, file = LogFilename, append = TRUE)

  # reactions list
  tmp = paste(SpecName,
        apply(SpecStoich, MARGIN = 1, FUN = function(X){
          X.nonzero = X[X != 0]
          X.react = names(X.nonzero)
          return(gsub(" [+] -", " - ",
                      paste(paste(X.nonzero, X.react, sep = " * "),
                            collapse = " + ")))
        }),
        sep = " = "
  )
  ColWidth = max(nchar(tmp))
  tmp = paste0(tmp, strrep(" ", ColWidth - nchar(tmp)))
  tmp = c(paste0("Reaction", strrep(" ", ColWidth - 8)),
          strrep("-", ColWidth),
          tmp)
  tmp = paste0(tmp, strrep(" ", 4))
  tmp = paste0(tmp, c(
    paste0(strrep(" ", 3), "LogK"),
    strrep("-", 7),
    formatC(SpecLogK, digits = 3, width = 7, format = "f", flag = " ")
  ))
  tmp = paste0(tmp, strrep(" ", 4))
  tmp = paste0(tmp, c(
    paste0(strrep(" ", 1), "DeltaH"),
    strrep("-", 7),
    formatC(SpecDeltaH, width = 7, format = "d", flag = " ")
  ))
  write(strrep("-", 80), file = LogFilename, append = TRUE)
  write(tmp, file = LogFilename, append = TRUE)

  # FixedAct Components
  write(strrep("-", 80), file = LogFilename, append = TRUE)
  write("FixedAct Components:", file = LogFilename, append = TRUE)
  write(CompName[CompType == "FixedAct"], file = LogFilename, append = TRUE)

  # FixedConc Components
  write(strrep("-", 80), file = LogFilename, append = TRUE)
  write("FixedConc Components:", file = LogFilename, append = TRUE)
  write(CompName[CompType == "FixedConc"], file = LogFilename, append = TRUE)

  # MassBal Totals
  tmp = paste0(
    "T.", CompName[CompType == "MassBal"],
    " = ",
    apply(SpecStoich[, CompType == "MassBal", drop = FALSE], MARGIN = 2, FUN = function(X){
      X.nonzero = X[X != 0]
      X.species = names(X.nonzero)
      return(gsub(" [+] -", " - ",
                  paste(paste(X.nonzero, X.species, sep = " * "),
                        collapse = " + ")))
    })
  )
  write(strrep("-", 80), file = LogFilename, append = TRUE)
  write("MassBal Totals:", file = LogFilename, append = TRUE)
  write(tmp, file = LogFilename, append = TRUE)

  # WHAM Totals
  tmp = paste0(
    "T.", CompName[CompType %in% c("WHAMHA","WHAMFA")],
    " = ",
    apply(SpecStoich[, CompType %in% c("WHAMHA","WHAMFA"), drop = FALSE], MARGIN = 2, FUN = function(X){
      X.nonzero = X[X != 0]
      X.species = names(X.nonzero)
      return(gsub(" [+] -", " - ",
                  paste(paste(X.nonzero, X.species, sep = " * "),
                        collapse = " + ")))
    })
  )
  write(strrep("-", 80), file = LogFilename, append = TRUE)
  write("WHAM Component Totals:", file = LogFilename, append = TRUE)
  write(tmp, file = LogFilename, append = TRUE)

  # DonnanChargeBal Totals
  write(strrep("-", 80), file = LogFilename, append = TRUE)
  write("Donnan Charge Balance Totals:", file = LogFilename, append = TRUE)
  for (iHS in c("HA", "FA")) {
    if (any(CompType == paste0("WHAM", iHS))) {
      DonnanComp = paste0("Donnan", iHS)
      tmp = paste0(
        "Z_Donnan", iHS,
        " = ",
        apply(SpecStoich[, DonnanComp, drop = FALSE], MARGIN = 2, FUN = function(X){
          X.nonzero = X[X != 0]
          X.species = names(X.nonzero)
          return(gsub(" [+] -", " - ",
                      paste(paste(X.nonzero, X.species, sep = " * "),
                            collapse = " + ")))
        })
      )
      write(tmp, file = LogFilename, append = TRUE)

      HSSpecName = SpecName[apply(SpecStoich[, CompType == paste0("WHAM", iHS), drop = FALSE], MARGIN = 1, FUN = function(X){any(X != 0)})]
      HSSpecCharge = SpecCharge[SpecName %in% HSSpecName]
      Nonzero = HSSpecCharge != 0
      tmp = paste0("Z_", iHS, " = ",
                   gsub(" [+] -", " - ",
                        paste0(
                          paste(HSSpecCharge[Nonzero],
                                HSSpecName[Nonzero], sep = " * "),
                          collapse = " + ")))
      write(tmp, file = LogFilename, append = TRUE)
    }
  }

  return(invisible(TRUE))

}
