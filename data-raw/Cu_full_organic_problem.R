mypfile = system.file(file.path("extdata","ParameterFiles","Cu_full_organic_WATER23dH.dat4"), package = "BLMEngineInR", mustWork = TRUE)
Cu_full_organic_problem = DefineProblem(ParamFile = mypfile, WriteLog = FALSE)

usethis::use_data(Cu_full_organic_problem, overwrite = TRUE)

