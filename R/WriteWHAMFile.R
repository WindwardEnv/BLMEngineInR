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

#' @title Write a WHAM Parameter File
#'
#' @description This function will take a WHAM parameter list object and turn it
#'   into a WHAM parameter file, effectively doing the opposite of `DefineWHAM`.
#'
#' @param ThisWHAM A list object with a structure like that returned by
#'   `BlankWHAM()`.
#' @param WHAMFile a character value, indicating the file path and name of the
#'   WHAM parameter file to write.
#' @param Notes A character vector of additional notes to include at the bottom
#'   of the WHAM parameter file. The text "written by USERNAME from R:
#'   YYYY-MM-DD HH:MM:SS" will always be written, regardless of the value of
#'   this argument. Be default, this will be filled in with a "Notes" item in
#'   `ThisWHAM`, if available.
#'
#' @return ThisProblem, with the ParamFile element  changed to the ParamFile
#'   argument.
#'
#' @examples
#' tf = tempfile()
#' WriteWHAMFile(ThisWHAM = Cu_full_organic_problem$WHAM, WHAMFile = tf)
#' DefineWHAM(WHAMFile = tf)
#'
#' @export
WriteWHAMFile = function(ThisWHAM, WHAMFile, Notes = ThisWHAM$Notes) {

  SepLine = "--------------------------------------------------------"

  MakeUniformTXTColumn = function(X) {
    if (length(unique(nchar(X))) > 1) {
      X = paste0(X, strrep(" ", max(nchar(X)) - nchar(X)))
    }
    return(X)
  }

  # Header Info
  Tmp = "...WHAM"
  if (!is.na(ThisWHAM$Ver)) {
    Tmp = paste(Tmp, ThisWHAM$Ver)
  } else {
    Tmp = paste0(Tmp, " - CUSTOM")
  }
  Tmp = paste0(Tmp, "...")
  write(Tmp, file = WHAMFile)
  write(SepLine, file = WHAMFile, append = TRUE)

  # Counts
  Tmp = paste(
    c(
      "Number of monodentate sites",
      "Number of bidentate pairs",
      "Number of tridentate groups",
      "Number of metals-OM parameters",
      "Number of non-standard selectivity coefficients"
    ),
    c(
      nrow(ThisWHAM$MonodentTable),
      nrow(ThisWHAM$BidentTable),
      nrow(ThisWHAM$TridentTable),
      nrow(ThisWHAM$MetalsTable),
      nrow(ThisWHAM$SpecKselTable)
    ),
    sep = ","
  )
  write(Tmp, file = WHAMFile, append = TRUE)

  # Constants
  Tmp = paste(
    c("Double layer overlap factor", "Constant to control DDL at low ZED"),
    c(ThisWHAM$DLF, ThisWHAM$KZED),
    sep = ","
  )
  write(Tmp, file = WHAMFile, append = TRUE)
  write(SepLine, file = WHAMFile, append = TRUE)

  # Parameter table
  ParamNames = c("nA", "pKA", "pKB", "dpKA", "dpKB", "fprB", "fprT", "dLK1A",
                 "dLK1B", "P", "Radius", "MolWt")
  Tmp = MakeUniformTXTColumn(paste0(
    c("Param", ParamNames),
    ", "
  ))
  Tmp = MakeUniformTXTColumn(paste0(
    Tmp,
    c("units", "eq/g OM", "n/a", "n/a", "n/a", "n/a",  "fraction", "fraction",
      "n/a", "n/a", "n/a", "m", "g OM/mol"),
    ", "
  ))
  Tmp = MakeUniformTXTColumn(paste0(
    Tmp,
    c("HA", sapply(ThisWHAM[ParamNames], FUN = function(X) {X["HA"]})),
    ", "
  ))
  Tmp = paste0(
    Tmp,
    c("FA", sapply(ThisWHAM[ParamNames], FUN = function(X) {X["FA"]}))
  )
  write(Tmp, file = WHAMFile, append = TRUE)
  write(SepLine, file = WHAMFile, append = TRUE)

  # Monodentate Sites table
  write("Monodentate Sites", file = WHAMFile, append = TRUE)
  Tmp = MakeUniformTXTColumn(paste0(c("S", ThisWHAM$MonodentTable$S), ", "))
  Tmp = MakeUniformTXTColumn(
    paste0(Tmp, c("AbundDenom", ThisWHAM$MonodentTable$AbundDenom), ", ")
  )
  Tmp = paste0(Tmp, c("StrongWeak", ThisWHAM$MonodentTable$StrongWeak))
  write(Tmp, file = WHAMFile, append = TRUE)
  write(SepLine, file = WHAMFile, append = TRUE)

  # Bidentate Sites Table
  write("Bidentate Sites", file = WHAMFile, append = TRUE)
  Tmp = MakeUniformTXTColumn(paste0(c("S1", ThisWHAM$BidentTable$S1), ", "))
  Tmp =
    MakeUniformTXTColumn(paste0(Tmp, c("S2", ThisWHAM$BidentTable$S2), ", "))
  Tmp = paste0(Tmp, c("AbundDenom", ThisWHAM$BidentTable$AbundDenom))
  write(Tmp, file = WHAMFile, append = TRUE)
  write(SepLine, file = WHAMFile, append = TRUE)

  # Tridentate Sites Table
  write("Tridentate Sites", file = WHAMFile, append = TRUE)
  Tmp = MakeUniformTXTColumn(paste0(c("S1", ThisWHAM$TridentTable$S1), ", "))
  Tmp =
    MakeUniformTXTColumn(paste0(Tmp, c("S2", ThisWHAM$TridentTable$S2), ", "))
  Tmp =
    MakeUniformTXTColumn(paste0(Tmp, c("S3", ThisWHAM$TridentTable$S3), ", "))
  Tmp = paste0(Tmp, c("AbundDenom", ThisWHAM$TridentTable$AbundDenom))
  write(Tmp, file = WHAMFile, append = TRUE)
  write(SepLine, file = WHAMFile, append = TRUE)

  # Metals-OM Parameters Table
  write("Metals Parameters", file = WHAMFile, append = TRUE)
  Tmp = MakeUniformTXTColumn(
    paste0(c("Metal", ThisWHAM$MetalsTable$Metal), ", ")
  )
  Tmp = MakeUniformTXTColumn(
    paste0(Tmp, c("pKMAHA", ThisWHAM$MetalsTable$pKMAHA), ", ")
  )
  Tmp = MakeUniformTXTColumn(
    paste0(Tmp, c("pKMAFA", ThisWHAM$MetalsTable$pKMAFA), ", ")
  )
  Tmp = paste0(Tmp, c("dLK2", ThisWHAM$MetalsTable$dLK2))
  write(Tmp, file = WHAMFile, append = TRUE)
  write(SepLine, file = WHAMFile, append = TRUE)

  # Selectivity Coefficients
  write("Selectivity Coefficients for Non-specific binding",
        file = WHAMFile, append = TRUE)
  Tmp =
    MakeUniformTXTColumn(paste0(c("Spec", ThisWHAM$SpecKselTable$Spec), ", "))
  Tmp = MakeUniformTXTColumn(
    paste0(Tmp, c("KselHA", ThisWHAM$SpecKselTable$KselHA), ", ")
  )
  Tmp = paste0(Tmp, c("KselFA", ThisWHAM$SpecKselTable$KselFA))
  write(Tmp, file = WHAMFile, append = TRUE)
  write(SepLine, file = WHAMFile, append = TRUE)

  # Notes
  write("Notes", file = WHAMFile, append = TRUE)
  write(paste0("written by ", Sys.info()["user"], " from R: ", Sys.time()),
        file = WHAMFile, append = TRUE)
  if (!is.null(Notes)) { write(Notes, file = WHAMFile, append = TRUE) }

  ThisWHAM$File = WHAMFile
  return(invisible(ThisWHAM))

}
