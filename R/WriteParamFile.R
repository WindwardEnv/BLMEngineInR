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
#'
#' @return (invisibly) TRUE, if successful.
#'
#' @examples
#' tf = tempfile()
#' WriteParamFile(ThisProblem = carbonate_system_problem, ParamFile = tf)
#' DefineProblem(tf)
#'
#' @export
WriteParamFile = function(ThisProblem, ParamFile) {

  CheckBLMObject(ThisProblem, BlankProblem(), BreakOnError = TRUE)

  MakeUniformTXTColumn = function(X) {
    if (length(unique(nchar(X))) > 1) {
      X = paste0(X, strrep(" ", max(nchar(X)) - nchar(X)))
    }
    return(X)
  }

  SectionBreak = "--------------------------------------------------------------------------------"

  write("Column model parameter file, Ver 4.00", file = ParamFile, append = FALSE)
  write(SectionBreak, file = ParamFile, append = TRUE)


  Tmp = c(
    ThisProblem$N[c("InMass","InLab","InVar","InComp","InDefComp","InSpec","Phase")],
    sum(ThisProblem$N[c("BL","Metal","BLMetal")]) + ThisProblem$WHAM$DoWHAM,
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

  TmpTable = ThisProblem$Mass[ThisProblem$Mass$Name %in% ThisProblem$InMassName, ]
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
  Tmp = MakeUniformTXTColumn(paste0(c("Variable", ThisProblem$InVar$Name), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Compartment", ThisProblem$InVar$MCName), ", "))
  Tmp = paste0(Tmp, c("Type", ThisProblem$InVar$Type))
  write(Tmp, file = ParamFile, append = TRUE)
  write(SectionBreak, file = ParamFile, append = TRUE)

  TmpTable = ThisProblem$Comp[ThisProblem$Comp$Name %in% ThisProblem$InCompName, ]
  write("Input Components", file = ParamFile, append = TRUE)
  Tmp = MakeUniformTXTColumn(paste0(c("Component", TmpTable$Name), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Charge", TmpTable$Charge), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Compartment", TmpTable$MCName), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Type", TmpTable$Type), ", "))
  Tmp = paste0(Tmp, c("Activity", TmpTable$ActCorr))
  write(Tmp, file = ParamFile, append = TRUE)
  write(SectionBreak, file = ParamFile, append = TRUE)

  TmpTable = ThisProblem$DefComp[ThisProblem$DefComp$Name %in% ThisProblem$InDefCompName, ]
  TmpTable$From = TmpTable$FromNum
  TmpTable$From[is.na(TmpTable$From)] = TmpTable$FromVar[is.na(TmpTable$From)]
  write("Defined Components", file = ParamFile, append = TRUE)
  Tmp = MakeUniformTXTColumn(paste0(c("Component", TmpTable$Name), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("From", TmpTable$From), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Charge", TmpTable$Charge), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Compartment", TmpTable$MCName), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Type", TmpTable$Type), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Activity", TmpTable$ActCorr), ", "))
  Tmp = paste0(Tmp, c("Site Den", TmpTable$SiteDens))
  write(Tmp, file = ParamFile, append = TRUE)
  write(SectionBreak, file = ParamFile, append = TRUE)

  TmpTable = ThisProblem$Spec[ThisProblem$Spec$Name %in% ThisProblem$InSpecName, ]
  TmpList = EquationToStoich(SpecEquation = TmpTable$Equation, CompName = ThisProblem$Comp$Name)
  write("Formation Reactions", file = ParamFile, append = TRUE)
  Tmp = MakeUniformTXTColumn(paste0(c("Species", TmpTable$Name), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Compartment", TmpTable$MCName), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Activity", TmpTable$ActCorr), ", "))
  Tmp = MakeUniformTXTColumn(paste0(Tmp, c("NC", TmpTable$NC), ", "))
  for (i in 1:max(TmpTable$NC)) {
    Tmp = MakeUniformTXTColumn(
      paste0(Tmp,
             c(paste0("C",i,", "),
               sapply(TmpList$SpecCompNames, function(X){if(length(X)<i){""}else{paste0(X[i],", ")}}))))
    Tmp = MakeUniformTXTColumn(
      paste0(Tmp,
             c(paste0("S", i,", "),
               sapply(TmpList$SpecCompStoichs, function(X){if(length(X)<i){""}else{paste0(X[i],", ")}}))))
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
      Tmp = MakeUniformTXTColumn(
        paste0(Tmp,
               c(paste0("C",i,", "),
                 sapply(TmpList$SpecCompNames, function(X){if(length(X)<i){""}else{paste0(X[i],", ")}}))))
      Tmp = MakeUniformTXTColumn(
        paste0(Tmp,
               c(paste0("S", i,", "),
                 sapply(TmpList$SpecCompStoichs, function(X){if(length(X)<i){""}else{paste0(X[i],", ")}}))))
    }

    Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Log K", TmpTable$LogK), ", "))
    Tmp = MakeUniformTXTColumn(paste0(Tmp, c("Delta H", TmpTable$DeltaH), ", "))
    Tmp = MakeUniformTXTColumn(paste0(Tmp, c("TempK", TmpTable$TempKelvin), ", "))
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
  if (ThisProblem$WHAM$DoWHAM) {
    TmpTable = rbind(TmpTable,
                     data.frame(
                       A = "WHAM",
                       B = ifelse(is.na(ThisProblem$WHAM$WHAMVer),
                                      ThisProblem$WHAM$WdatFile,
                                      ThisProblem$WHAM$WHAMVer)
                     ))
  }
  Tmp = MakeUniformTXTColumn(paste0(TmpTable$A, ", "))
  Tmp = paste0(Tmp, TmpTable$B)
  write(Tmp, file = ParamFile, append = TRUE)
  write(SectionBreak, file = ParamFile, append = TRUE)

  TmpTable = ThisProblem$CATab
  write("Critical Accumulation Table", file = ParamFile, append = TRUE)
  Tmp = MakeUniformTXTColumn(paste0(c("Num", TmpTable$Num), ", "))
  for (i in c("CA","Species","Test.Type","Duration","Lifestage","Endpoint",
              "Quantifier","References")) {
    i.bool = grepl(",", TmpTable[, i])
    TmpTable[i.bool, i] = paste0("\"", TmpTable[i.bool, i], "\"")
    Tmp = MakeUniformTXTColumn(paste0(Tmp, c(gsub("[.]"," ",i), TmpTable[, i]), ", "))
  }
  Tmp = paste0(Tmp, c("Miscellaneous", TmpTable$Miscellaneous))
  write(Tmp, file = ParamFile, append = TRUE)
  write(SectionBreak, file = ParamFile, append = TRUE)

  write("-------Notes--------", file = ParamFile, append = TRUE)
  write(paste0("written by ", Sys.info()["user"], " from R: ", Sys.time()), file = ParamFile, append = TRUE)

  invisible(TRUE)
}
