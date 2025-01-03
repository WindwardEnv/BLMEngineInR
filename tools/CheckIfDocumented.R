# Read in the CommonParameterDefinitions script, and parse out the @param
# arguments. If a BLM function has either an argument or an output value that is
# not documented in this code, it should be.
CommonDocumentation = function(Mode = c("check_outputs",
                                        "check_inputs",
                                        "generate_param_list",
                                        "generate_return_list"),
                               Func,
                               ...) {
  Tmp = scan(file = "R/CommonParameterDefinitions.R",
             what = character(),
             sep = "\n",
             quiet = TRUE)

  Tmp = gsub("#'", "", Tmp)
  Tmp = trimws(Tmp)
  Tmp = Tmp[Tmp != ""]

  i = 1
  while (i < length(Tmp)) {
    if (substr(Tmp[i], 1, 1) != "@") {
      Tmp[i - 1] = paste(Tmp[i - 1], Tmp[i])
      Tmp = Tmp[-i]
      next
    }

    i = i + 1
  }

  ParameterDocumentation = gsub("^@param ", "", Tmp[grepl("^@param ", Tmp)])

  ParameterList = data.frame(
    Name = gsub("^([[:alnum:]]+) .+", "\\1", ParameterDocumentation),
    Description = gsub("^[[:alnum:]]+ ", "", ParameterDocumentation)
  )

  if (any(grepl("check_", Mode))) {
    NotDocumented = c()
    if ("check_outputs" %in% Mode) {
      out = Func(...)
      NotDocumented = c(NotDocumented,
                        names(out)[names(out) %in% ParameterList$Name == FALSE])
    }
    if ("check_inputs" %in% Mode) {
      NotDocumented = c(NotDocumented,
                        formalArgs(Func)[
                          formalArgs(Func) %in% ParameterList$Name == FALSE])
    }
    if (length(NotDocumented) > 0) {
      print("These parameters were not documented:")
      return(NotDocumented)
    }
  }

  if ("generate_param_list" %in% Mode) {
    InputParameters = formalArgs(Func)
    ThisFunctionParams = ParameterList[
      ParameterList$Name %in% InputParameters, ]

    InputParametersNotDocumented = InputParameters[
      InputParameters %in% ThisFunctionParams$Name == FALSE]
    if (length(InputParametersNotDocumented) > 0) {
      ThisFunctionParams = rbind(
        ThisFunctionParams,
        data.frame(Name = InputParametersNotDocumented,
                   Description = "???")
      )
    }
    ThisFunctionParams = ThisFunctionParams[match(InputParameters,
                                                  ThisFunctionParams$Name), ]

    ParamList = paste("@param",
                       ThisFunctionParams$Name,
                       ThisFunctionParams$Description)
    return(paste0("#' ", strwrap(ParamList, width = 77, exdent = 2)))
  } else if ("generate_return_list" %in% Mode) {
    ReturnParameters = names(Func(...))
    ThisFunctionParams = ParameterList[
      ParameterList$Name %in% ReturnParameters, ]
    ReturnParametersNotDocumented = ReturnParameters[
      ReturnParameters %in% ThisFunctionParams$Name == FALSE]
    if (length(ReturnParametersNotDocumented) > 0) {
      ThisFunctionParams = rbind(
        ThisFunctionParams,
        data.frame(Name = ReturnParametersNotDocumented,
                   Description = "???")
      )
    }
    ThisFunctionParams = ThisFunctionParams[match(ReturnParameters,
                                                  ThisFunctionParams$Name), ]
    ParamList = paste0("\\item{", ThisFunctionParams$Name, "}",
                      "{", ThisFunctionParams$Description, "}")
    ParamList = strwrap(ParamList, width = 75, exdent = 2)
    ParamList = c("\\describe{", paste0("  ", ParamList), "}")
    ParamList = paste0("#' ", ParamList)
    return(ParamList)
  }

  return(TRUE)

}

GlobalTest = list(
  ParamFile = "scrap/parameter file format/full_organic.dat4",
  InputFile = "scrap/parameter file format/full_organic.blm4",
  DoTox = T,
  CritAccumIndex = 1L,
  QuietFlag ="Quiet",
  ConvergenceCriteria = 0.001,
  MaxIter = 30L
)
GlobalTest = c(GlobalTest, DefineProblem(ParamFile))
FunctionInputs = GlobalTest[
  which(names(GlobalTest) %in% formalArgs("GetData"))]
FunctionInputs$InputFile = InputFile
GlobalTest = c(GlobalTest, do.call("GetData", args = FunctionInputs))
GlobalTest$iObs = 1L
GlobalTest$InLab = GlobalTest$InLabObs[GlobalTest$iObs, ]
GlobalTest$SysTempKelvin = GlobalTest$SysTempKelvinObs[GlobalTest$iObs]
GlobalTest$TotConc = GlobalTest$TotConcObs[GlobalTest$iObs, ]
GlobalTest$CATargetDefault = GlobalTest$CATab$CA[GlobalTest$CritAccumIndex] * (10^-6) /
  GlobalTest$CompSiteDens[GlobalTest$BLComp]
GlobalTest$CATarget = CATargetDefault * GlobalTest$TotConc[GlobalTest$BLComp]

CommonDocumentation(Mode = "check_inputs", Func = DefineProblem)
CommonDocumentation(Mode = c("check_outputs", "generate_return_list"), Func = DefineProblem,
                    ParamFile = GlobalTest$ParamFile)

CommonDocumentation(Mode = "check_inputs", Func = GetData)
CommonDocumentation(Mode = c("check_outputs", "generate_return_list"),
                    Func = GetData,
                    InputFile = GlobalTest$InputFile,
                    NInLab = GlobalTest$NInLab,
                    InLabName  = GlobalTest$InLabName,
                    NInVar = GlobalTest$NInVar,
                    InVarName  = GlobalTest$InVarName,
                    InVarMCR  = GlobalTest$InVarMCR,
                    InVarType  = GlobalTest$InVarType,
                    NInComp  = GlobalTest$NInComp,
                    InCompName = GlobalTest$InCompName,
                    NComp  = GlobalTest$NComp,
                    CompName = GlobalTest$CompName,
                    NDefComp = GlobalTest$NDefComp,
                    DefCompName  = GlobalTest$DefCompName,
                    DefCompFromNum = GlobalTest$DefCompFromNum,
                    DefCompFromVar = GlobalTest$DefCompFromVar,
                    DefCompSiteDens  = GlobalTest$DefCompSiteDens)

paste0(formalArgs(MatchInputsToProblem), " = GlobalTest$", formalArgs(MatchInputsToProblem))
CommonDocumentation(Mode = "check_inputs", Func = MatchInputsToProblem)
CommonDocumentation(Mode = c("check_outputs", "generate_return_list"),
                    Func = MatchInputsToProblem,
                    NObs = GlobalTest$NObs,
                    InVarObs = GlobalTest$InVarObs,
                    InCompObs = GlobalTest$InCompObs,
                    NInVar = GlobalTest$NInVar,
                    InVarName = GlobalTest$InVarName,
                    InVarMCR = GlobalTest$InVarMCR,
                    InVarType = GlobalTest$InVarType,
                    NInComp = GlobalTest$NInComp,
                    InCompName = GlobalTest$InCompName,
                    NComp = GlobalTest$NComp,
                    CompName = GlobalTest$CompName,
                    NDefComp = GlobalTest$NDefComp,
                    DefCompName = GlobalTest$DefCompName,
                    DefCompFromNum = GlobalTest$DefCompFromNum,
                    DefCompFromVar = GlobalTest$DefCompFromVar,
                    DefCompSiteDens = GlobalTest$DefCompSiteDens)

paste0(formalArgs(ReadInputsFromFile), " = GlobalTest$", formalArgs(ReadInputsFromFile))
CommonDocumentation(Mode = "check_inputs", Func = ReadInputsFromFile)
CommonDocumentation(Mode = c("check_outputs", "generate_return_list"),
                    Func = ReadInputsFromFile,
                    InputFile = GlobalTest$InputFile,
                    NInLab = GlobalTest$NInLab,
                    InLabName = GlobalTest$InLabName,
                    NInVar = GlobalTest$NInVar,
                    InVarName = GlobalTest$InVarName,
                    NInComp = GlobalTest$NInComp,
                    InCompName = GlobalTest$InCompName)

paste0(formalArgs(CHESS), " = GlobalTest$", formalArgs(CHESS))
formalArgs(CHESS)[formalArgs(CHESS) %in% names(GlobalTest) == FALSE]
CommonDocumentation(Mode = c("check_inputs", "generate_param_list"), Func = CHESS)
CommonDocumentation(Mode = c("check_outputs", "generate_return_list"),
                    Func = CHESS,
                    QuietFlag = GlobalTest$QuietFlag,
                    ConvergenceCriteria = GlobalTest$ConvergenceCriteri,
                    MaxIter = GlobalTest$MaxIter,
                    NComp = GlobalTest$NComp,
                    NSpec = GlobalTest$NSpec,
                    NBLMetal = GlobalTest$NBLMetal,
                    SpecK = GlobalTest$SpecK,
                    SpecTempKelvin = GlobalTest$SpecTempKelvin,
                    SpecDeltaH = GlobalTest$SpecDeltaH,
                    SpecStoich = GlobalTest$SpecStoich,
                    SpecName = GlobalTest$SpecName,
                    CompType = GlobalTest$CompType,
                    CompName = GlobalTest$CompName,
                    TotConc = GlobalTest$TotConc,
                    SysTempKelvin = GlobalTest$SysTempKelvin,
                    DoTox = GlobalTest$DoTox,
                    MetalName = GlobalTest$MetalName,
                    MetalComp = GlobalTest$MetalComp,
                    BLMetalSpecs = GlobalTest$BLMetalSpecs,
                    CATarget = GlobalTest$CATarget)

