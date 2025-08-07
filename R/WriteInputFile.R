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

#' @title Write a BLM input file
#'
#' @description This function will take a BLM inputs list object and turn it
#'   into an input file, effectively doing the opposite of `GetData`.
#'
#' @param AllInput A list object with a structure like that returned by
#'   `GetData()`.
#' @param ThisProblem A list object with a structure like that returned by
#' @param InputFile `BlankProblem()`.
#' @return TRUE (invisibly) if successful.
#'
#' @examples
#' tf = tempfile()
#' myinputfile = system.file(file.path("extdata", "InputFiles",
#'                                     "carbonate_system_test.blm4"),
#'                           package = "BLMEngineInR",
#'                           mustWork = TRUE)
#' myinputs = GetData(InputFile = myinputfile,
#'                    ThisProblem = carbonate_system_problem)
#' WriteInputFile(AllInput = myinputs, ThisProblem = carbonate_system_problem,
#'                InputFile = tf)
#' scan(tf, what = character(), sep = "\n")
#' scan(myinputfile, what = character(), sep = "\n")
#' file.remove(tf)
#'
#' @export
WriteInputFile = function(AllInput, ThisProblem, InputFile) {

  CheckBLMObject(AllInput, BlankInputList(ThisProblem), BreakOnError = TRUE)

  MakeUniformTXTColumn = function(X) {
    if (length(unique(nchar(X))) > 1) {
      X = paste0(X, strrep(" ", max(nchar(X)) - nchar(X)))
    }
    X
  }

  write(AllInput$NObs, file = InputFile, append = FALSE)

  # Initialize character vector
  TextToWrite = character(AllInput$NObs + 2)

  # Input labels
  if (ThisProblem$N["InLab"] > 0) {
    for (i in 1:ThisProblem$N["InLab"]) {
      TextToWrite = MakeUniformTXTColumn(
        paste0(TextToWrite,
               c(ThisProblem$InLabName[i], "n/a",
                 AllInput$InLabObs[, ThisProblem$InLabName[i]]),
               ",  ")
      )
    }
  }

  # Input variables
  for (i in 1:ThisProblem$N["InVar"]) {
    ThisColumn = c(ThisProblem$InVar$Name[i], "",
                   AllInput$InVarObs[, ThisProblem$InVar$Name[i]])
    if (ThisProblem$InVar$Type[i] == "pH") {
      ThisColumn[2] = "SU"
    } else if (ThisProblem$InVar$Type[i] == "Temperature") {
      ThisColumn[2] = "deg C"
    } else if (ThisProblem$InVar$Type[i] %in%
                 c("WHAM-HA", "WHAM-FA", "WHAM-HAFA")) {
      ThisColumn[2] = "mg C/L"
    } else if (ThisProblem$InVar$Type[i] %in% c("PercHA", "PercAFA")) {
      ThisColumn[2] = "%"
    }
    TextToWrite = MakeUniformTXTColumn(
      paste0(TextToWrite,
             ThisColumn,
             ",  ")
    )
  }

  # Input components
  for (i in 1:ThisProblem$N["InComp"]) {
    TextToWrite = MakeUniformTXTColumn(
      paste0(TextToWrite,
             c(ThisProblem$InCompName[i], "mol/L",
               AllInput$InCompObs[, ThisProblem$InCompName[i]]),
             ",  ")
    )
  }

  TextToWrite = gsub(",$", "", trimws(TextToWrite))

  write(TextToWrite, file = InputFile, append = TRUE)
  write("\nNotes:\n------------", file = InputFile, append = TRUE)
  write(paste0("ParameterFile = \"", ThisProblem$ParamFile, "\""),
        file = InputFile, append = TRUE)

  invisible(TRUE)
}
