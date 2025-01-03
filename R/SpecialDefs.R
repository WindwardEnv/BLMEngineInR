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

#' @name SpecialDefs
#'
#' @title Add or remove species definitions
#'
#' @description The special definitions in a parameter file include indicating
#'   the biotic ligand species relevant to toxicity ("BL"), the toxic metal
#'   ("Metal"), the species responsible for the critical accumulation associated
#'   with toxicity at the biotic ligand ("BL-Metal"), and the model version of
#'   the Windemere Humic Aqueous Model to use to represent organic matter
#'   binding ("WHAM").
#'
#' @param ThisProblem A list object with a structure like that returned by
#'   `BlankProblem()`.
#' @param Value A character vector. When `SpecialDef` is either `"BL"` or
#'   `"Metal"`, this should be the name of a component in
#'   `ThisProblem$Chem$Name`. When `SpecialDef` is `"BL-Metal"`, this should be
#'   the name of a chemical species in `ThisProblem$Spec$Name`. When
#'   `SpecialDef` is `"WHAM"`, this should be either a supported WHAM version
#'   number (i.e., one of `"V"`, `"VI"`, or `"VII"`), or the file path to a WHAM
#'   parameters file (.wdat file) that follows the format of one of the standard
#'   versions supplied with this package (see
#'   `system.file("extdata/WHAM/WHAM_V.wdat", package = "BLMEngineInR")` for an
#'   example).
#' @param SpecialDef A character vector indicating which special definition to
#'   add a value for. Valid values are `"BL"`, `"Metal"`, `"BL-Metal"`,
#'   `"BLMetal"` (same as `"BL-Metal"`), and `"WHAM"`.
#' @param SpecialDefToRemove The name of the special definition to remove.
#' @param Index If applicable (such as if there are two BL-Metal species), the
#'   index of which to remove (i.e., the first one or second one).
#' @param DoCheck A logical value indicating whether checks should be performed
#'   on the incoming and outgoing problem objects. Defaults to `TRUE`, as you
#'   usually want to make sure something isn't awry, but the value is often set
#'   to `FALSE` when used internally (like in DefineProblem) so the problem is
#'   only checked once at the end.
#'
#' @return `ThisProblem`, with the special definitions changed.
#'
#' @family problem manipulation functions
#'
#' @examples
#' print(carbonate_system_problem[c("BL","Metal","BLMetal","WHAM")])
#' my_new_problem = carbonate_system_problem
#' my_new_problem = AddInComps(ThisProblem = my_new_problem, InCompName = "Cu",
#'                             InCompCharge = 2,
#'                             InCompMCName = "Water",
#'                             InCompType = "MassBal",
#'                             InCompActCorr = "Debye")
#' my_new_problem = AddSpecialDefs(ThisProblem = my_new_problem,
#'                                 Value = "Cu",
#'                                 SpecialDef = "Metal")
#' print(my_new_problem[c("BL","Metal","BLMetal","WHAM")])
#' my_new_problem = RemoveSpecialDefs(ThisProblem = my_new_problem,
#'                                    SpecialDefToRemove = "Metal")
#' print(my_new_problem[c("BL","Metal","BLMetal","WHAM")])
NULL


#' @rdname SpecialDefs
#' @export
AddSpecialDefs = function(ThisProblem, Value, SpecialDef, DoCheck = TRUE) { #nolint: cyclocomp_linter

  if (DoCheck) {
    CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  }
  NewProblem = ThisProblem

  if ((NewProblem$ParamFile != "") &&
        !grepl("[(]modified[)]$", NewProblem$ParamFile)) {
    NewProblem$ParamFile = paste0(NewProblem$ParamFile, " (modified)")
  }

  if (any(is.na(Value)) || any(is.na(SpecialDef))) {
    stop("NA inputs not allowed.")
  }
  NSpecialDef = length(Value)
  if (NSpecialDef > length(SpecialDef)) {
    SpecialDef = rep(SpecialDef, NSpecialDef)
  }
  Value = trimws(as.character(Value))

  for (i in 1:NSpecialDef) {
    SpecialDef[i] = match.arg(SpecialDef[i],
                              c("BL", "Metal", "BLMetal", "BL-Metal", "WHAM"))
    SpecialDef[i] = gsub("-", "", SpecialDef[i])
    if (SpecialDef[i] == "WHAM") {
      #--Name of WHAM file or WHAM version
      if (any(grepl("WHAM", ThisProblem$DefComp$Type)) ||
            !is.na(ThisProblem$WHAM$File) || !is.na(ThisProblem$WHAM$Ver)) {
        stop("Only one WHAM version or file can be specified.")
      }

      if (Value[i] %in% c("V", "VI", "VII")) {
        WHAMVer = Value[i]
        WHAMFile = NA_character_
      } else {
        WHAMVer = NA_character_
        WHAMFile = Value[i]
        if (!file.exists(WHAMFile)) {
          ParamFileDir = file.path(dirname(NewProblem$ParamFile), WHAMFile)
          SysFileDir = system.file(file.path("extdata", "WHAM", WHAMFile),
                                   package = "BLMEngineInR")
          if (file.exists(ParamFileDir)) {
            WHAMFile = ParamFileDir
          } else if (file.exists(SysFileDir)) {
            WHAMFile = SysFileDir
          } else {
            stop("Cannot find \"", WHAMFile, "\".")
          }
        }
      }
      NewProblem$WHAM = DefineWHAM(WHAMVer = WHAMVer, WHAMFile = WHAMFile)
      if (any(grepl("WHAM", ThisProblem$InVar$Type))) {
        NewProblem = ExpandWHAM(ThisProblem = NewProblem)
      }
    } else if (SpecialDef[i] %in% c("BL", "Metal")) {
      if (Value[i] %in% ThisProblem$Comp$Name) {
        NewProblem[[SpecialDef[i]]] = rbind(
          NewProblem[[SpecialDef[i]]],
          data.frame(
            Name = Value[i],
            CompR = match(Value[i], ThisProblem$Comp$Name)
          )
        )
        NewProblem$N[SpecialDef[i]] = NewProblem$N[SpecialDef[i]] + 1L
      } else {
        stop("Unknown component specifed for ", SpecialDef[i])
      }
    } else if (SpecialDef[i] %in% c("BL-Metal", "BLMetal")) {
      if (Value[i] %in% ThisProblem$Spec$Name) {
        NewProblem[[SpecialDef[i]]] = rbind(
          NewProblem[[SpecialDef[i]]],
          data.frame(
            Name = Value[i],
            SpecsR = match(Value[i], ThisProblem$Spec$Name)
          )
        )
        NewProblem$N[SpecialDef[i]] = NewProblem$N[SpecialDef[i]] + 1L
      } else {
        stop("Unknown species specifed for ", SpecialDef[i])
      }
    }
  }

  if (DoCheck) {
    CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  }
  return(NewProblem)

}


#' @rdname SpecialDefs
#' @export
RemoveSpecialDefs = function(ThisProblem, SpecialDefToRemove, Index = 1,
                             DoCheck = TRUE) {
  if (DoCheck) {
    CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  }
  NewProblem = ThisProblem

  if ((NewProblem$ParamFile != "") &&
        !grepl("[(]modified[)]$", NewProblem$ParamFile)) {
    NewProblem$ParamFile = paste0(NewProblem$ParamFile, " (modified)")
  }

  if (any(SpecialDefToRemove %in%
            c("BL", "Metal", "BLMetal", "BL-Metal", "WHAM") == FALSE)) {
    stop("SpecialDefToRemove should be one of 'BL', 'Metal', 'BLMetal', ",
         "'BL-Metal', or 'WHAM'.")
  }
  SpecialDefToRemove = gsub("[-]", "", SpecialDefToRemove)
  RemovalTable = data.frame(SpecialDef = SpecialDefToRemove, Index = Index)
  RemovalTable = RemovalTable[order(RemovalTable$SpecialDef,
                                    -RemovalTable$Index), ]

  for (i in unique(SpecialDefToRemove)) {
    if (i == "WHAM") {
      # Undo Expand WHAM

      # ...remove all WHAM defined components (this will remove components and
      # species along with it)
      NewProblem = RemoveDefComps(
        ThisProblem = NewProblem,
        DefCompToRemove =
          which(ThisProblem$DefComp$Type %in% c("WHAMHA", "WHAMFA")),
        DoCheck = DoCheck
      )

      # ...remove Donnan Mass compartments (this will remove components and
      # species along with it)
      NewProblem = RemoveMassCompartments(
        ThisProblem = NewProblem,
        MCToRemove = which(grepl("Donnan", ThisProblem$Mass$Name)),
        DoCheck = DoCheck
      )

      # ...reset WHAM values to NA
      NewProblem$WHAM = BlankProblem()$WHAM

    } else {
      NewProblem[[i]] =
        NewProblem[[i]][-RemovalTable$Index[RemovalTable$SpecialDef %in% i], ]
      NewProblem$N[i] = NewProblem$N[i] - sum(RemovalTable$SpecialDef %in% i)
    }
  }

  if (DoCheck) {
    CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  }
  return(NewProblem)
}
