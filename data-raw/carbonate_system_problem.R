mypfile = system.file(file.path("extdata","ParameterFiles","carbonate_system_only.dat4"),
                      package = "BLMEngineInR",
                      mustWork = TRUE)
carbonate_system_problem = DefineProblem(ParamFile = mypfile)

usethis::use_data(carbonate_system_problem, overwrite = TRUE)
