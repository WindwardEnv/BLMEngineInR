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
  QuietFlag ="Quiet",
  ConvergenceCriteria = 0.001,
  MaxIter = 30L
)

end.time = Sys.time()
end.time - start.time

ResultsTable[, c("Obs","FinalIter","FinalMaxError")]
ResultsTable[, c("Obs","ID","ID2","T.Cu","Cu")]

# Old BLM
# Site Label    Sample Label  Dis. Cu      Free Cu
# Full_Organic  Hard ser 1    1.17628E-08  1.17741E-10
# Full_Organic  Hard ser 2    1.15568E-08  1.49896E-10
# Full_Organic  Hard ser 3    1.26915E-08  2.41360E-10
# Full_Organic  Hard ser 4    1.50443E-08  3.86284E-10
# Full_Organic  Hard ser 5    1.89049E-08  6.60671E-10
# Full_Organic  Hard ser 6    2.70029E-08  1.42593E-09
# Full_Organic  Hard ser 7    3.65246E-08  2.62630E-09


# No TempCorrection
#   Obs           ID        ID2         T.Cu           Cu
# 1   1 Full_Organic Hard ser 1 2.246596e-09 9.909165e-11
# 2   2 Full_Organic Hard ser 2 2.881979e-09 1.284319e-10
# 3   3 Full_Organic Hard ser 3 4.619184e-09 2.104946e-10
# 4   4 Full_Organic Hard ser 4 7.121356e-09 3.337472e-10
# 5   5 Full_Organic Hard ser 5 1.122678e-08 5.473762e-10
# 6   6 Full_Organic Hard ser 6 2.030386e-08 1.057243e-09
# 7   7 Full_Organic Hard ser 7 3.135992e-08 1.720029e-09

# With TempCorrection (15C)
#   Obs           ID        ID2         T.Cu           Cu
# 1   1 Full_Organic Hard ser 1 2.683681e-09 9.905403e-11
# 2   2 Full_Organic Hard ser 2 3.445260e-09 1.282940e-10
# 3   3 Full_Organic Hard ser 3 5.522465e-09 2.096111e-10
# 4   4 Full_Organic Hard ser 4 8.512210e-09 3.309170e-10
# 5   5 Full_Organic Hard ser 5 1.342391e-08 5.397773e-10
# 6   6 Full_Organic Hard ser 6 2.432579e-08 1.035269e-09
# 7   7 Full_Organic Hard ser 7 3.764053e-08 1.677717e-09
