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

#' @title Read a WHAM file and make a WHAM list
#'
#' @description A WHAM file is a text file (typically with the file extension
#'   ".wdat") that has all of the information necessary for defining organic
#'   matter binding, according to the Windemere Humic Aqueous Model (WHAM). Only
#'   constants relating directly to organic matter binding are in this object
#'   and file (i.e., nothing related to inorganic binding). This is the
#'   information needed by the `ExpandWHAM()` function and to do organic matter
#'   binding in the `CHESS` subroutine.
#'
#' @param WHAMVer a character string specifying the WHAM version to use, must be
#'   one of `"V"` (default), `"VI"`, or `"VII"`. Ignored if `WHAMFile` is not
#'   `NA`.
#' @param WHAMFile (optional) a character string specifying the file path of a
#'   WHAM parameter file
#'
#' @return A WHAM list in the format of `BlankWHAM()`.
#'
#' @export
#'
DefineWHAM = function(WHAMVer = "V", WHAMFile = NA) {


  # error catching and input cleanup
  WHAMVer = match.arg(WHAMVer, choices =  c("V", "VI", "VII", NA_character_))
  if (is.na(WHAMFile)) {
    if (WHAMVer == "V") {
      WHAMFile = system.file(file.path("extdata", "WHAM", "WHAM_V.wdat"),
                             package = "BLMEngineInR",
                             mustWork = TRUE)
    } else if (WHAMVer == "VI") {
      WHAMFile = system.file(file.path("extdata", "WHAM", "WHAM_VI.wdat"),
                             package = "BLMEngineInR",
                             mustWork = TRUE)
    } else if (WHAMVer  == "VII") {
      WHAMFile = system.file(file.path("extdata", "WHAM", "WHAM_VII.wdat"),
                             package = "BLMEngineInR",
                             mustWork = TRUE)
    }
  } else {
    WHAMVer = NA_character_
    WHAMFile = normalizePath(WHAMFile)
    stopifnot(file.exists(WHAMFile))
  }

  NewWHAM = BlankWHAM()
  NewWHAM$Ver = WHAMVer
  NewWHAM$File = WHAMFile

  # header info
  SkipRows = 2L
  Tmp = read.delim(
    file = WHAMFile,
    header = FALSE,
    sep = ",",
    skip = SkipRows,
    nrows = 7L
  )
  nMS = as.integer(Tmp[1, 2])#Number of monodentate sites
  nBP = as.integer(Tmp[2, 2])#Number of bidentate pairs
  nTG = as.integer(Tmp[3, 2])#Number of tridentate groups
  nMP = as.integer(Tmp[4, 2])#Number of metals-OM parameters
  nKsel = as.integer(Tmp[5, 2])#Number of non-standard selectivity coefficients
  NewWHAM$DLF = as.numeric(Tmp[6, 2])#Double layer overlap factor
  NewWHAM$KZED = as.numeric(Tmp[7, 2])#Constant to control DDL at low ZED

  # Parameters
  SkipRows = SkipRows + 7L + 1L
  ParamNames = c("nA", "pKA", "pKB", "dpKA", "dpKB", "fprB", "fprT", "dLK1A",
                 "dLK1B", "P", "Radius", "MolWt")
  Tmp = read.delim(
    file = WHAMFile,
    header = TRUE,
    sep = ",",
    skip = SkipRows,
    nrows = 12L,
    row.names = ParamNames
  )
  for (i in 1:12) {
    NewWHAM[[ParamNames[i]]] = as.numeric(Tmp[i, 3:4])
    names(NewWHAM[[ParamNames[i]]]) = c("HA", "FA")
  }

  # Monodentate Sites - these should always be the same, but we'll set things
  # up like this so we can add in this section if it's ever needed.
  # MonodentTable = data.frame(S=1:8, AbundDenom = c(rep(4,4),rep(8,4)))
  SkipRows = SkipRows + 12L + 3L
  if (nMS > 0) {
    NewWHAM$MonodentTable = read.delim(
      file = WHAMFile,
      header = TRUE,
      sep = ",",
      skip = SkipRows,
      nrows = nMS,
      col.names = colnames(NewWHAM$MonodentTable),
      strip.white = TRUE
    )
  }

  # Bidentate Pairs
  SkipRows = SkipRows + nMS + 3L
  if (nBP > 0) {
    NewWHAM$BidentTable = read.delim(
      file = WHAMFile,
      header = TRUE,
      sep = ",",
      skip = SkipRows,
      nrows = nBP,
      col.names = colnames(NewWHAM$BidentTable),
      strip.white = TRUE
    )
  }


  # Tridentate Groups
  SkipRows = SkipRows + nBP + 3L
  if (nTG > 0) {
    NewWHAM$TridentTable = read.delim(
      file = WHAMFile,
      header = TRUE,
      sep = ",",
      skip = SkipRows,
      nrows = nTG,
      col.names = colnames(NewWHAM$TridentTable),
      strip.white = TRUE
    )
  }

  # Metals Parameters Table
  SkipRows = SkipRows + nTG + 3L
  if (nMP > 0) {
    NewWHAM$MetalsTable = read.delim(
      file = WHAMFile,
      header = TRUE,
      sep = ",",
      skip = SkipRows,
      nrows = nMP,
      col.names = colnames(NewWHAM$MetalsTable),
      strip.white = TRUE
    )
    NewWHAM$MetalsTable$pKMAHA = as.numeric(NewWHAM$MetalsTable$pKMAHA)
    NewWHAM$MetalsTable$pKMAFA = as.numeric(NewWHAM$MetalsTable$pKMAFA)
    NewWHAM$MetalsTable$dLK2 = as.numeric(NewWHAM$MetalsTable$dLK2)
    # MetalsTable$pKMBHA = 3 * MetalsTable$pKMAHA - 3
    # MetalsTable$pKMBFA = 3.96 * MetalsTable$pKMAFA
  }

  # Non-standard selectivity coefficients
  SkipRows = SkipRows + nMP + 3L
  if (nKsel > 0) {
    NewWHAM$SpecKselTable = read.delim(
      file = WHAMFile,
      header = TRUE,
      sep = ",",
      skip = SkipRows,
      nrows = nKsel,
      col.names = colnames(NewWHAM$SpecKselTable),
      strip.white = TRUE
    )
  }

  CheckBLMObject(NewWHAM, BlankWHAM(), BreakOnError = TRUE)

  return(NewWHAM)

}
