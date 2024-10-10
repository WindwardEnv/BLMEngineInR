mypfile = system.file(file.path("extdata","ParameterFiles","Cu_inorganic_only.dat4"), package = "BLMEngineInR", mustWork = TRUE)
Cu_full_inorganic_problem = DefineProblem(ParamFile = mypfile, WriteLog = FALSE)

usethis::use_data(Cu_full_inorganic_problem, overwrite = TRUE)
