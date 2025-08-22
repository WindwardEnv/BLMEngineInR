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

#' @title Write a BLM Parameter File
#'
#' @description This function will take a BLM chemical problem list object and
#'   turn it into a parameter file, effectively doing the opposite of
#'   `DefineProblem`.
#'
#' @param ThisProblem A list object with a structure like that returned by
#'   `BlankProblem()`.
#' @param ParamFile a character value, indicating the file path and name of the
#'   parameter file to write.
#' @param Notes A character vector of additional notes to include at the bottom
#'   of the parameter file. The text "written by USERNAME from R: YYYY-MM-DD
#'   HH:MM:SS" will always be written, regardless of the value of this argument.
#'   This will be filled in with a "Notes" item in `ThisProblem`, if available.
#'
#' @return ThisProblem, with the ParamFile element  changed to the ParamFile
#'   argument.
#'
#' @examples
#' tf = tempfile()
#' WriteParamFile(ThisProblem = carbonate_system_problem, ParamFile = tf)
#' DefineProblem(ParamFile = tf)
#'
#' @export
WriteParamFile = function(ThisProblem, ParamFile, Notes = ThisProblem$Notes) {


  CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)

  MakeUniformTXTColumn = function(X) {
    if (length(unique(nchar(X))) > 1) {
      X = paste0(X, strrep(" ", max(nchar(X)) - nchar(X)))
    }
    X
  }

  SectionBreak = strrep("-", 80)

  write("Column model parameter file, Ver 4.00",
        file = ParamFile, append = FALSE)
  write(SectionBreak, file = ParamFile, append = TRUE)

  Tmp = c(
    ThisProblem$N[c("InMass", "InLab", "InVar", "InComp", "InDefComp",
                    "InSpec", "Phase")],
    sum(ThisProblem$N[c("BL", "Metal", "BLMetal")]) +
      ((!is.na(ThisProblem$WHAM$Ver) | !is.na(ThisProblem$WHAM$File))),
    ThisProblem$N["CAT"]
  )
  Tmp = MakeUniformTXTColumn(paste0(Tmp, ", "))
  Tmp = paste0(Tmp, "Number of ",
               c("Mass Compartments",
                 "Input Labels",
                 "Input Variables",
                 "Input Components",
                 "Defined Components",
                 "Species",
                 "Phases",
                 "Special Definitions",
                 "Critical Values"))
  write(Tmp, file = ParamFile, append = TRUE)
  write(SectionBreak, file = ParamFile, append = TRUE)

  TmpTable =
    ThisProblem$Mass[ThisProblem$Mass$Name %in% ThisProblem$InMassName, ]
  write("Mass Compartment List", file = ParamFile, append = TRUE)
  Tmp = MakeUniformTXTColumn(paste0(c("Compartment", TmpTable$Name), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("C to M", TmpTable$Amt), ", "))
  Tmp = paste0(Tmp, c("Unit Label", TmpTable$Unit))
  write(Tmp, file = ParamFile, append = TRUE)
  write(SectionBreak, file = ParamFile, append = TRUE)

  write("Input Labels", file = ParamFile, append = TRUE)
  write(ThisProblem$InLabName, file = ParamFile, append = TRUE)
  write(SectionBreak, file = ParamFile, append = TRUE)

  write("Input Variables", file = ParamFile, append = TRUE)
  Tmp =
    MakeUniformTXTColumn(paste0(c("Variable", ThisProblem$InVar$Name), ", "))
  Tmp = MakeUniformTXTColumn(
    paste0(Tmp, c("Compartment", ThisProblem$InVar$MCName), ", ")
  )
  Tmp = paste0(Tmp, c("Type", ThisProblem$InVar$Type))
  write(Tmp, file = ParamFile, append = TRUE)
  write(SectionBreak, file = ParamFile, append = TRUE)

  TmpTable =
    ThisProblem$Comp[ThisProblem$Comp$Name %in% ThisProblem$InCompName, ]
  write("Input Components", file = ParamFile, append = TRUE)
  Tmp = MakeUniformTXTColumn(paste0(c("Component", TmpTable$Name), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Charge", TmpTable$Charge), ", "))
  Tmp =
    MakeUniformTXTColumn(paste0(Tmp, c("Compartment", TmpTable$MCName), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Type", TmpTable$Type), ", "))
  Tmp = paste0(Tmp, c("Activity", TmpTable$ActCorr))
  write(Tmp, file = ParamFile, append = TRUE)
  write(SectionBreak, file = ParamFile, append = TRUE)

  TmpTable = ThisProblem$DefComp[ThisProblem$DefComp$Name %in% ThisProblem$InDefCompName, ] #nolint: line_length_linter
  TmpTable$From = TmpTable$FromNum
  TmpTable$From[is.na(TmpTable$From)] = TmpTable$FromVar[is.na(TmpTable$From)]
  write("Defined Components", file = ParamFile, append = TRUE)
  Tmp = MakeUniformTXTColumn(paste0(c("Component", TmpTable$Name), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("From", TmpTable$From), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Charge", TmpTable$Charge), ", "))
  Tmp =
    MakeUniformTXTColumn(paste0(Tmp, c("Compartment", TmpTable$MCName), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Type", TmpTable$Type), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Activity", TmpTable$ActCorr), ", "))
  Tmp = paste0(Tmp, c("Site Den", TmpTable$SiteDens))
  write(Tmp, file = ParamFile, append = TRUE)
  write(SectionBreak, file = ParamFile, append = TRUE)

  TmpTable =
    ThisProblem$Spec[ThisProblem$Spec$Name %in% ThisProblem$InSpecName, ]
  TmpList = EquationToStoich(SpecEquation = TmpTable$Equation,
                             CompName = ThisProblem$Comp$Name)
  write("Formation Reactions", file = ParamFile, append = TRUE)
  Tmp = MakeUniformTXTColumn(paste0(c("Species", TmpTable$Name), ", "))
  Tmp =
    MakeUniformTXTColumn(paste0(Tmp, c("Compartment", TmpTable$MCName), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Activity", TmpTable$ActCorr), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("NC", TmpTable$NC), ", "))
  for (i in 1:max(TmpTable$NC)) {
    Tmp = MakeUniformTXTColumn(paste0(Tmp, c(
      paste0("C", i, ", "),
      sapply(TmpList$SpecCompNames, function(X) {
        if (length(X) < i) {
          ""
        } else{
          paste0(X[i], ", ")
        }
      })
    )))
    Tmp = MakeUniformTXTColumn(paste0(Tmp, c(
      paste0("S", i, ", "),
      sapply(TmpList$SpecCompStoichs, function(X) {
        if (length(X) < i) {
          ""
        } else{
          paste0(X[i], ", ")
        }
      })
    )))
  }

  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Log K", TmpTable$LogK), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Delta H", TmpTable$DeltaH), ", "))
  Tmp = paste0(Tmp, c("TempK", TmpTable$TempKelvin))
  write(Tmp, file = ParamFile, append = TRUE)
  write(SectionBreak, file = ParamFile, append = TRUE)


  TmpTable = ThisProblem$Phase
  write("Phase List", file = ParamFile, append = TRUE)
  if (ThisProblem$N["Phase"] > 0) {
    TmpList = EquationToStoich(SpecEquation = TmpTable$Equation,
                               CompName = ThisProblem$Comp$Name)
    Tmp = MakeUniformTXTColumn(paste0(c("Phases", TmpTable$Name), ", "))
    Tmp = MakeUniformTXTColumn(paste0(Tmp, c("NC", TmpTable$NC), ", "))
    for (i in 1:max(TmpTable$NC)) {
      Tmp = MakeUniformTXTColumn(paste0(Tmp, c(
        paste0("C", i, ", "),
        sapply(TmpList$SpecCompNames, function(X) {
          if (length(X) < i) {
            ""
          } else{
            paste0(X[i], ", ")
          }
        })
      )))
      Tmp = MakeUniformTXTColumn(paste0(Tmp, c(
        paste0("S", i, ", "),
        sapply(TmpList$SpecCompStoichs, function(X) {
          if (length(X) < i) {
            ""
          } else {
            paste0(X[i], ", ")
          }
        })
      )))
    }

    Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Log K", TmpTable$LogK), ", "))
    Tmp =
      MakeUniformTXTColumn(paste0(Tmp, c("Delta H", TmpTable$DeltaH), ", "))
    Tmp =
      MakeUniformTXTColumn(paste0(Tmp, c("TempK", TmpTable$TempKelvin), ", "))
    Tmp = paste0(Tmp, c("Moles", TmpTable$Moles))
  } else {
    Tmp = "Phase, NC, C1, S1, C2, S2, C3, S3, Log K, Delta H, Temp, Moles"
  }
  write(Tmp, file = ParamFile, append = TRUE)
  write(SectionBreak, file = ParamFile, append = TRUE)

  write("Special Definitions", file = ParamFile, append = TRUE)
  TmpTable = data.frame(A = "Definition", B = "Value")
  if (ThisProblem$N["BL"] > 0) {
    TmpTable = rbind(TmpTable,
                     data.frame(
                       A = rep("BL", ThisProblem$N["BL"]),
                       B = ThisProblem$BL$Name
                     ))
  }
  if (ThisProblem$N["Metal"] > 0) {
    TmpTable = rbind(TmpTable,
                     data.frame(
                       A = rep("Metal", ThisProblem$N["Metal"]),
                       B = ThisProblem$Metal$Name
                     ))
  }
  if (ThisProblem$N["BLMetal"] > 0) {
    TmpTable = rbind(TmpTable,
                     data.frame(
                       A = rep("BL-Metal", ThisProblem$N["BLMetal"]),
                       B = ThisProblem$BLMetal$Name
                     ))
  }
  # if (ThisProblem$DoWHAM) {
  if (!is.na(ThisProblem$WHAM$Ver) || !is.na(ThisProblem$WHAM$File)) {
    if (is.na(ThisProblem$WHAM$File)) {
      TmpTable = rbind(TmpTable,
                       data.frame(
                         A = "WHAM",
                         B = ThisProblem$WHAM$Ver
                       ))
    } else {
      Tmp = system.file(
        file.path("extdata", "WHAM", basename(ThisProblem$WHAM$File)),
        package = "BLMEngineInR"
      )
      if (file.exists(Tmp) &&
            (basename(ThisProblem$WHAM$File) %in%
               paste0("WHAM_", c("V", "VI", "VII")))) {
        TmpTable = rbind(TmpTable,
                         data.frame(
                           A = "WHAM",
                           B = ThisProblem$WHAM$Ver
                         ))
      } else {
        TmpTable = rbind(TmpTable,
                         data.frame(
                           A = "WHAM",
                           B = basename(ThisProblem$WHAM$File)
                         ))
      }
    }

  }
  Tmp = MakeUniformTXTColumn(paste0(TmpTable$A, ", "))
  Tmp = paste0(Tmp, TmpTable$B)
  write(Tmp, file = ParamFile, append = TRUE)
  write(SectionBreak, file = ParamFile, append = TRUE)

  TmpTable = ThisProblem$CATab
  write("Critical Accumulation Table", file = ParamFile, append = TRUE)
  Tmp = MakeUniformTXTColumn(paste0(c("Num", TmpTable$Num), ", "))
  for (i in c("CA", "Species", "TestType", "Duration", "Lifestage", "Endpoint",
              "Quantifier", "References")) {
    iBool = grepl(",", TmpTable[, i])
    TmpTable[iBool, i] = paste0("\"", TmpTable[iBool, i], "\"")
    Tmp = MakeUniformTXTColumn(
      paste0(Tmp, c(gsub("[.]", " ", i), TmpTable[, i]), ", ")
    )
  }
  Tmp = paste0(Tmp, c("Miscellaneous", TmpTable$Miscellaneous))
  write(Tmp, file = ParamFile, append = TRUE)
  write(SectionBreak, file = ParamFile, append = TRUE)

  write("-------Notes--------", file = ParamFile, append = TRUE)
  write(paste0("written by ", Sys.info()["user"], " from R: ", Sys.time()),
        file = ParamFile, append = TRUE)
  if (!is.null(Notes)) { write(Notes, file = ParamFile, append = TRUE) }

  ThisProblem$ParamFile = ParamFile
  return(invisible(ThisProblem))

}
