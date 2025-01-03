# The Ni problem for organisms sensitive to bicarbonate
Ni.w.pfile = normalizePath("C:/Users/kellyc/documents/BLM/parameter_files/fresh/Ni/Ni_all_studies_calibration_toxic_carbonates_2024-04-10.dat")
Ni.r.pfile = file.path("inst", "extdata", "ParameterFiles", "Ni_HCO3_full_organic.dat4")

Ni_HCO3_full_organic_problem = ConvertWindowsParamFile(
  WindowsParamFile = Ni.w.pfile,
  RParamFile = Ni.r.pfile
)

usethis::use_data(Ni_HCO3_full_organic_problem, overwrite = TRUE)
