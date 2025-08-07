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

#' @title Make a blank WHAM parameter list object
#'
#' @return A list object with a template for defining the organic matter binding
#'   in a chemical problem for the `BLMEngineInR` functions. Each element in the
#'   list is a vector, matrix, or data.frame object grouping related parameters
#'   together. See `str(BlankWHAM())` for the structure and names of the list
#'   object.
#'
#' @export
BlankWHAM = function() {

  NAVec = c(HA = NA_real_, FA = NA_real_)

  Out = list(
    Ver = NA_character_,
    File = NA_character_,
    DLF = NA_real_,
    KZED = NA_real_,
    nA = NAVec,
    pKA = NAVec,
    pKB = NAVec,
    dpKA = NAVec,
    dpKB = NAVec,
    fprB = NAVec,
    fprT = NAVec,
    dLK1A = NAVec,
    dLK1B = NAVec,
    P = NAVec,
    Radius = NAVec,
    MolWt = NAVec,
    # Constants = c(
    #   DLF = 0.0,
    #   KZED = 0.0
    # ),
    # Param = matrix(
    #   data = NA_real_,
    #   nrow = 12,
    #   ncol = 2,
    #   dimnames = list(c("nA", "pKA", "pKB", "dpKA", "dpKB", "fprB", "fprT",
    #                     "dLK1A", "dLK1B", "P", "radius", "MolWt"),
    #                   c("HA","FA"))
    # ),
    MonodentTable = data.frame(
      S = integer(),
      AbundDenom = integer(),
      StrongWeak = character()
    ),
    BidentTable = data.frame(
      S1 = integer(),
      S2 = integer(),
      AbundDenom = integer()
    ),
    TridentTable = data.frame(
      S1 = integer(),
      S2 = integer(),
      S3 = integer(),
      AbundDenom = integer()
    ),
    MetalsTable = data.frame(
      Metal = character(),
      pKMAHA = numeric(),
      pKMAFA = numeric(),
      dLK2 = numeric()
    ),
    SpecKselTable = data.frame(
      Spec = character(),
      KselHA = numeric(),
      KselFA = numeric()
    )
  )

  Out

}
