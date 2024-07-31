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
AddSpecialDefs = function(ThisProblem, Value, SpecialDef) {

  CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  NewProblem = ThisProblem

  if (any(is.na(Value))) {
    stop("NA inputs not allowed.")
  }
  NSpecialDef = length(Value)
  if (NSpecialDef > length(SpecialDef)) {
    SpecialDef = rep(SpecialDef, NSpecialDef)
  }
  Value = trimws(as.character(Value))

  for (i in 1:NSpecialDef) {
    SpecialDef[i] = match.arg(SpecialDef[i], c("BL", "Metal", "BLMetal","BL-Metal", "WHAM"))
    SpecialDef[i] = gsub("-", "", SpecialDef[i])
    if (SpecialDef[i] == "WHAM") {
      #--Name of WHAM file or WHAM version
      if (any(grepl("WHAM", ThisProblem$InVar$Type))) {
        if (!is.na(ThisProblem$WHAM$WHAMVer) ||
            !is.na(ThisProblem$WHAM$WdatFile)) {
          stop("Only one WHAM version or file can be specified.")
        }


        if (Value[i] %in% c("V", "VI", "VII")) {
          WHAMVer = Value[i]
          WdatFile = NA
        } else {
          WHAMVer = NA
          WdatFile = Value[i]
        }
        NewProblem = ExpandWHAM(ThisProblem = NewProblem,
                                WHAMVer = WHAMVer,
                                WdatFile = WdatFile)
      } else {
        stop("WHAM version or file specified without a WHAM input variable.")
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
    } else if (SpecialDef[i] %in% c("BL-Metal","BLMetal")) {
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

  CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  return(NewProblem)

}


#' @rdname SpecialDefs
#' @export
RemoveSpecialDefs = function(ThisProblem, SpecialDefToRemove, Index = 1) {
  CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)
  NewProblem = ThisProblem

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
        DefCompToRemove = which(ThisProblem$DefComp$Type %in% c("WHAMHA","WHAMFA"))
      )

      # ...remove Donnan Mass compartments (this will remove components and
      # species along with it)
      NewProblem = RemoveMassCompartments(
        ThisProblem = NewProblem,
        MCToRemove = which(grepl("Donnan", ThisProblem$Mass$Name))
      )

      # ...reset WHAM values to NA
      NewProblem$WHAM = BlankProblem()$WHAM

    } else {
      NewProblem[[i]] = NewProblem[[i]][-RemovalTable$Index[RemovalTable$SpecialDef %in% i], ]
      NewProblem$N[i] = NewProblem$N[i] - sum(RemovalTable$SpecialDef %in% i)
    }
  }

  CheckBLMObject(NewProblem, BlankProblem(), BreakOnError = TRUE)
  return(NewProblem)
}
