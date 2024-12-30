# The regular Ni file
Ni.w.pfile = normalizePath("C:/Users/kellyc/documents/BLM/parameter_files/fresh/Ni/Ni_all_studies_calibration_2024-12-18.dat")
Ni.r.pfile = file.path("inst","extdata","ParameterFiles","Ni_full_organic.dat4")

Ni_full_organic_problem = ConvertWindowsParamFile(
  WindowsParamFile = Ni.w.pfile,
  RParamFile = Ni.r.pfile)

usethis::use_data(Ni_full_organic_problem, overwrite = TRUE)

