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
                    LogFilename = file.path(dirname(ThisProblem$ParamFile),
                                            "CHESSLOG.txt")) {

  # initialize
  CompName = character()
  SpecName = character()
  SpecStoich = matrix()
  SpecLogK = numeric()
  SpecDeltaH = numeric()
  CompType = character()
  SpecCharge = integer()

  # initialize
  CompName = character()
  SpecName = character()
  SpecStoich = matrix()
  SpecLogK = numeric()
  SpecDeltaH = numeric()
  CompType = character()
  #CompActCorr = character()
  SpecCharge = integer()

  # unpack the input list
  ThisProblemList = ConvertToList(ThisProblem)
  for (i in names(ThisProblemList)) {
    assign(i, ThisProblemList[[i]])
  }

  # initialize log file
  write(paste0("CHESS problem defined by '", ThisProblem$ParamFile,
               "' parameter file:\n",
               "(", Sys.time(), ")"),
        file = LogFilename, append = FALSE)

  # Component List
  write(strrep("-", 80), file = LogFilename, append = TRUE)
  write("Component List:", file = LogFilename, append = TRUE)
  write(CompName, file = LogFilename, append = TRUE)

  # reactions list
  Tmp = paste(
    SpecName,
    apply(
      SpecStoich,
      MARGIN = 1,
      FUN = function(X) {
        XNonzero = X[X != 0]
        XReact = names(XNonzero)
        gsub(" [+] -", " - ",
             paste(paste(XNonzero, XReact, sep = " * "),
                   collapse = " + "))
      }
    ),
    sep = " = "
  )
  ColWidth = max(nchar(Tmp))
  Tmp = paste0(Tmp, strrep(" ", ColWidth - nchar(Tmp)))
  Tmp = c(paste0("Reaction", strrep(" ", ColWidth - 8)),
          strrep("-", ColWidth),
          Tmp)
  Tmp = paste0(Tmp, strrep(" ", 4))
  Tmp = paste0(Tmp, c(
    paste0(strrep(" ", 3), "LogK"),
    strrep("-", 7),
    formatC(SpecLogK, digits = 3, width = 7, format = "f", flag = " ")
  ))
  Tmp = paste0(Tmp, strrep(" ", 4))
  Tmp = paste0(Tmp, c(
    paste0(strrep(" ", 1), "DeltaH"),
    strrep("-", 7),
    formatC(SpecDeltaH, width = 7, format = "d", flag = " ")
  ))
  write(strrep("-", 80), file = LogFilename, append = TRUE)
  write(Tmp, file = LogFilename, append = TRUE)

  # FixedAct Components
  write(strrep("-", 80), file = LogFilename, append = TRUE)
  write("FixedAct Components:", file = LogFilename, append = TRUE)
  write(CompName[CompType == "FixedAct"], file = LogFilename, append = TRUE)

  # FixedConc Components
  write(strrep("-", 80), file = LogFilename, append = TRUE)
  write("FixedConc Components:", file = LogFilename, append = TRUE)
  write(CompName[CompType == "FixedConc"], file = LogFilename, append = TRUE)

  # MassBal Totals
  Tmp = paste0(
    "T.", CompName[CompType == "MassBal"],
    " = ",
    apply(
      SpecStoich[, CompType == "MassBal", drop = FALSE],
      MARGIN = 2,
      FUN = function(X) {
        XNonzero = X[X != 0]
        XSpecies = names(XNonzero)
        gsub(" [+] -", " - ",
             paste(paste(XNonzero, XSpecies, sep = " * "),
                   collapse = " + "))
      }
    )
  )
  write(strrep("-", 80), file = LogFilename, append = TRUE)
  write("MassBal Totals:", file = LogFilename, append = TRUE)
  write(Tmp, file = LogFilename, append = TRUE)

  # WHAM Totals
  Tmp = paste0(
    "T.", CompName[CompType %in% c("WHAMHA", "WHAMFA")],
    " = ",
    apply(
      SpecStoich[, CompType %in% c("WHAMHA", "WHAMFA"), drop = FALSE],
      MARGIN = 2,
      FUN = function(X) {
        XNonzero = X[X != 0]
        XSpecies = names(XNonzero)
        gsub(" [+] -", " - ",
             paste(paste(XNonzero, XSpecies, sep = " * "),
                   collapse = " + "))
      }
    )
  )
  write(strrep("-", 80), file = LogFilename, append = TRUE)
  write("WHAM Component Totals:", file = LogFilename, append = TRUE)
  write(Tmp, file = LogFilename, append = TRUE)

  # DonnanChargeBal Totals
  write(strrep("-", 80), file = LogFilename, append = TRUE)
  write("Donnan Charge Balance Totals:", file = LogFilename, append = TRUE)
  for (iHS in c("HA", "FA")) {
    if (any(CompType == paste0("WHAM", iHS))) {
      DonnanComp = paste0("Donnan", iHS)
      Tmp = paste0(
        "Z_Donnan", iHS,
        " = ",
        apply(
          SpecStoich[, DonnanComp, drop = FALSE],
          MARGIN = 2,
          FUN = function(X) {
            XNonzero = X[X != 0]
            XSpecies = names(XNonzero)
            gsub(" [+] -", " - ",
                 paste(paste(XNonzero, XSpecies, sep = " * "),
                       collapse = " + "))
          }
        )
      )
      write(Tmp, file = LogFilename, append = TRUE)

      HSSpecName = SpecName[
        apply(SpecStoich[, CompType == paste0("WHAM", iHS), drop = FALSE],
              MARGIN = 1, FUN = function(X) {any(X != 0)})
      ]
      HSSpecCharge = SpecCharge[SpecName %in% HSSpecName]
      Nonzero = HSSpecCharge != 0
      Tmp = paste0("Z_", iHS, " = ", gsub(" [+] -", " - ", paste0(
        paste(HSSpecCharge[Nonzero], HSSpecName[Nonzero], sep = " * "),
        collapse = " + "
      )))
      write(Tmp, file = LogFilename, append = TRUE)
    }
  }

  return(invisible(TRUE))

}
