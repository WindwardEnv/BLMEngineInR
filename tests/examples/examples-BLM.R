# running the BLM function with a parameter file and input file:
mypfile = system.file("extdata","ParameterFiles","carbonate_system_only.dat4",
                      package = "BLMEngineInR",
                      mustWork = TRUE)
myinputfile = system.file("extdata","InputFiles","carbonate_system_test.blm4",
                          package = "BLMEngineInR",
                          mustWork = TRUE)
BLM(ParamFile = mypfile, InputFile = myinputfile, DoTox = FALSE)

# running the BLM with parameter and input objects
myinputs = GetData(InputFile = myinputfile, ThisProblem = carbonate_system_problem)
BLM(ThisProblem = carbonate_system_problem, AllInput = myinputs, DoTox = FALSE)

# here we only read in the same files, but the inputs could also be constructed

