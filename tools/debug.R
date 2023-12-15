rm(list = ls())
# devtools::clean_dll()
devtools::load_all()

# test stuff
start.time = Sys.time()

ResultsTable = BLM(
  ParamFile = "scrap/parameter file format/full_organic.dat4",
  InputFile = "scrap/parameter file format/full_organic.blm4",
  DoTox = T,
  iCA = 1L,
  QuietFlag ="Very Quiet",
  ConvergenceCriteria = 0.001,
  MaxIter = 30L
)

end.time = Sys.time()
end.time - start.time


