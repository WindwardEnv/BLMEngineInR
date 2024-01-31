rm(list = ls())
# devtools::clean_dll()
devtools::load_all()

# test stuff
start.time = Sys.time()

ResultsTable = BLM(
  ParamFile = "scrap/parameter file format/full_organic_WATER23dH.dat4",
  InputFile = "scrap/parameter file format/full_organic.blm4",
  DoTox = F,
  # iCA = 1L,
  QuietFlag ="Quiet",
  ConvergenceCriteria = 0.001,
  MaxIter = 50L
)

end.time = Sys.time()
end.time - start.time

ResultsTable[, c("Obs","ID2","FinalIter","FinalMaxError")]
ResultsTable[, c("ID","ID2","T.Cu","Cu")]


OldBLMResultsTable = openxlsx::read.xlsx(xlsxFile = "scrap/old BLM/full_organic_SPEC.det.xlsx",
                                         sheet = 1, rows = c(5, 7:55))

(ResultsTable$DDL.Na = ResultsTable$T.Na - ResultsTable$Na)
ResultsTable$TOrg.Cu = rowSums(ResultsTable[, grepl("DOC", colnames(ResultsTable)) & grepl("Cu", colnames(ResultsTable))])
(OldBLMResultsTable$DDL.Na = OldBLMResultsTable$T.Na - OldBLMResultsTable$Na)
OldBLMResultsTable$ID = trimws(gsub("\"","", OldBLMResultsTable$Site.Label))
OldBLMResultsTable$ID2 = trimws(gsub("\"","", OldBLMResultsTable$Sample.Label))
ResultsTable[, c("ID","ID2","T.Cu","Cu","TOrg.Cu","DDL.Na")]
OldBLMResultsTable[, c("ID","ID2","T.Cu","Cu","TOrg.Cu","DDL.Na")]
OldBLMResultsTable$Act.OH = OldBLMResultsTable$OH


ResultsTable$SerLab = ""
ResultsTable$SerLab[grepl("Hard ser", ResultsTable$ID2)] =
  round((ResultsTable$T.Ca + ResultsTable$T.Mg) * 100086, 0)[grepl("Hard ser", ResultsTable$ID2)]
ResultsTable$SerLab[grepl("pH ser", ResultsTable$ID2)] =
  ResultsTable$pH[grepl("pH ser", ResultsTable$ID2)]
ResultsTable$SerLab[grepl("Temp ser", ResultsTable$ID2)] =
  ResultsTable$Temp[grepl("Temp ser", ResultsTable$ID2)]
ResultsTable$SerLab[grepl("DOC ser", ResultsTable$ID2)] =
  ResultsTable$DOC[grepl("DOC ser", ResultsTable$ID2)]

ResultsTable$Ser = gsub("ser .+","ser", ResultsTable$ID2)
leg.dat = data.frame(
  Ser = c("Hard ser", "DOC ser", "pH ser", "Temp ser",
          "hi DOC Hard ser", "hi DOC pH ser", "hi DOC Temp ser"),
  col = 2:8
)
ResultsTable$col = leg.dat$col[match(ResultsTable$Ser, leg.dat$Ser)]


if (FALSE) {

  log10(OldBLMResultsTable$Act.OH * ResultsTable$Act.H)[24:29]

  SysTempKelvin = (ResultsTable$Temp + 273.15)

  T1 = 1 / SysTempKelvin
  T0 = 1 / 298.15
  SpecK = 10^-14.03541
  Rcon = 8.314
  SpecDeltaH = 56190#-56197
  T2 = SpecDeltaH * (T0 - T1) / Rcon
  log10(SpecK * exp(T2))[24:29]

  OldBLMResultsTable = OldBLMResultsTable[!grepl("DOC ser", OldBLMResultsTable$ID2), ]
  ResultsTable = ResultsTable[!grepl("DOC ser", ResultsTable$ID2), ]


  pdf("scrap/speciation comparison.pdf")
  par(omi = rep(1, 4))
  # for (i.col in c("TOrg.Cu", "Cu", "CO3", "HCO3","Act.OH", "CuOH","Cu(OH)2")){
  for (i.col in colnames(ResultsTable)){
    if ((i.col %in% colnames(OldBLMResultsTable)) &&
        all(!is.na(is.numeric(ResultsTable[,i.col]))) &&
        !grepl("T[.]", i.col) &&
        (i.col %in% c("ID","ID2","#.Iter.", "DDL.Na", "TOrg.Cu") == FALSE)){
      ax.lim = range(OldBLMResultsTable[, i.col])
      ax.lim[1] = 10^(floor(log10(ax.lim[1])))
      ax.lim[2] = 10^(ceiling(log10(ax.lim[2])))
      plot(x = OldBLMResultsTable[,i.col],
           y = pmin(ax.lim[2], pmax(ax.lim[1], ResultsTable[,i.col])),
           col = c(rep(2, 7), rep(4, 7), rep(5, 6)),
           pch = ifelse(ResultsTable[,i.col] < ax.lim[1], 6,
                        ifelse(ResultsTable[,i.col] > ax.lim[2], 2, 16)),
           main = i.col, xlab = "BLM Version 3.41.2.45",
           ylab = "BLM in R",
           log = "xy", xlim = ax.lim, ylim = ax.lim, las = 2)
      text(x = OldBLMResultsTable[,i.col],
           y = pmin(ax.lim[2], pmax(ax.lim[1], ResultsTable[,i.col])),
           labels =ResultsTable$SerLab,
           col = c(rep(2, 7), rep(4, 7), rep(5, 6)),
           pos = 3, cex = 0.5)
      abline(a = 0, b = 1)
      legend("topleft", legend = c("Hard ser", "pH ser", "Temp ser"), pch = 16, col = c(2,4,5))
    }
  }
  dev.off()
}


for (i.col in c("TOrg.Cu", "Cu", "CO3", "HCO3","Act.OH", "CuOH","Cu(OH)2")){
# for (i.col in colnames(ResultsTable)){
    if ((i.col %in% colnames(OldBLMResultsTable)) &&
        all(!is.na(is.numeric(ResultsTable[,i.col]))) &&
        !grepl("T[.]", i.col) &&
        (i.col %in% c("ID","ID2","#.Iter.", "DDL.Na", "TOrg.Cu") == FALSE)){
    ax.lim = range(OldBLMResultsTable[, i.col])
    ax.lim[1] = 10^(floor(log10(ax.lim[1])))
    ax.lim[2] = 10^(ceiling(log10(ax.lim[2])))
    plot(x = OldBLMResultsTable[,i.col],
         y = pmin(ax.lim[2], pmax(ax.lim[1], ResultsTable[,i.col])),
         col = ResultsTable$col,
         pch = ifelse(ResultsTable[,i.col] < ax.lim[1], 6,
                      ifelse(ResultsTable[,i.col] > ax.lim[2], 2, 16)),
         main = i.col, xlab = "BLM Version 3.41.2.45",
         ylab = "BLM in R",
         log = "xy", xlim = ax.lim, ylim = ax.lim, las = 2)
    text(x = OldBLMResultsTable[,i.col],
         y = pmin(ax.lim[2], pmax(ax.lim[1], ResultsTable[,i.col])),
         labels = ResultsTable$SerLab, col = ResultsTable$col,
         pos = 3, cex = 0.5)
    abline(a = 0, b = 1)
    legend("topleft", legend = leg.dat$Ser, pch = 16, col = leg.dat$col)
  }
}


plot(
  x = OldBLMResultsTable$T.DOC,
  y = OldBLMResultsTable$TOrg.Cu / OldBLMResultsTable$T.Cu,
)


plot(x = OldBLMResultsTable$T.Cu, y = pmax(10^-9, ResultsTable$T.Cu),
     col = c(rep(2, 7), rep(3, 9), rep(4, 7)),
     pch = ifelse(ResultsTable$T.Cu < 10^-9, 6, 16),
     main = "Total Cu", xlab = "Old BLM", ylab = "Current BLM in R",
     log = "xy", xlim = c(1E-9, 1E-5), ylim = c(1E-9, 1E-5))
text(x = OldBLMResultsTable$T.Cu, y = pmax(10^-9, ResultsTable$T.Cu),
     labels = c(1:7, 1:9, 1:7), col = c(rep(2, 7), rep(3, 9), rep(4, 7)),
     pos = 3, cex = 0.5)
abline(a = 0, b = 1)
legend("topleft", legend = c("Hard ser", "DOC ser", "pH ser"), pch = 16, col = 2:4)
plot(x = OldBLMResultsTable$Cu, y = pmax(10^-11, ResultsTable$Cu),
     col = c(rep(2, 7), rep(3, 9), rep(4, 7)),
     pch = ifelse(ResultsTable$Cu < 10^-11, 6, 16),
     main = "Free Cu", xlab = "Old BLM", ylab = "Current BLM in R",
     log = "xy", xlim = c(1E-11, 1E-6), ylim = c(1E-11, 1E-6))
text(x = OldBLMResultsTable$Cu, y = pmax(10^-11, ResultsTable$Cu),
     labels = c(1:7, 1:9, 1:7), col = c(rep(2, 7), rep(3, 9), rep(4, 7)),
     pos = 3, cex = 0.5)
abline(a = 0, b = 1)
legend("topleft", legend = c("Hard ser", "DOC ser", "pH ser"), pch = 16, col = 2:4)

par(mfrow = c(1, 2))
plot(x = OldBLMResultsTable$Na, y = pmax(10^-5, ResultsTable$Na),
     col = c(rep(2, 7), rep(3, 9), rep(4, 7)),
     pch = ifelse(ResultsTable$Na < 10^-5, 6, 16),
     main = "Free Na", xlab = "Old BLM", ylab = "Current BLM in R",
     log = "xy", xlim = c(1E-5, 1E-2), ylim = c(1E-5, 1E-2))
text(x = OldBLMResultsTable$Na, y = pmax(10^-5, ResultsTable$Na),
     labels = c(1:7, 1:9, 1:7), col = c(rep(2, 7), rep(3, 9), rep(4, 7)),
     pos = 3, cex = 0.5)
abline(a = 0, b = 1)
legend("topleft", legend = c("Hard ser", "DOC ser", "pH ser"), pch = 16, col = 2:4)
plot(x = OldBLMResultsTable$DDL.Na, y = pmax(10^-8, ResultsTable$DDL.Na),
     col = c(rep(2, 7), rep(3, 9), rep(4, 7)),
     pch = ifelse(ResultsTable$DDL.Na < 10^-8, 6, 16),
     main = "Diffuse Bound Na", xlab = "Old BLM", ylab = "Current BLM in R",
     log = "xy", xlim = c(1E-8, 1E-6), ylim = c(1E-8, 1E-6))
text(x = OldBLMResultsTable$DDL.Na, y = pmax(10^-8, ResultsTable$DDL.Na),
     labels = c(1:7, 1:9, 1:7), col = c(rep(2, 7), rep(3, 9), rep(4, 7)),
     pos = 3, cex = 0.5)
abline(a = 0, b = 1)
legend("topleft", legend = c("Hard ser", "DOC ser", "pH ser"), pch = 16, col = 2:4)



# Old BLM (3 seconds)
# Site Label    Sample Label  Dis. Cu      Free Cu
# Full_Organic  Hard ser 1    1.17628E-08  1.17741E-10
# Full_Organic  Hard ser 2    1.15568E-08  1.49896E-10
# Full_Organic  Hard ser 3    1.26915E-08  2.41360E-10
# Full_Organic  Hard ser 4    1.50443E-08  3.86284E-10
# Full_Organic  Hard ser 5    1.89049E-08  6.60671E-10
# Full_Organic  Hard ser 6    2.70029E-08  1.42593E-09
# Full_Organic  Hard ser 7    3.65246E-08  2.62630E-09
# Full_Organic  DOC ser 1     3.04891E-08  3.93838E-10
# Full_Organic  DOC ser 2     1.13690E-07  3.94705E-10
# Full_Organic  DOC ser 3     2.18304E-07  3.96361E-10
# Full_Organic  DOC ser 4     4.28458E-07  3.98179E-10
# Full_Organic  DOC ser 5     1.06909E-06  4.03965E-10
# Full_Organic  DOC ser 6     1.61405E-06  4.08591E-10
# Full_Organic  DOC ser 7     2.16647E-06  4.12376E-10
# Full_Organic  DOC ser 8     3.28985E-06  4.18182E-10
# Full_Organic  DOC ser 9     4.43633E-06  4.22807E-10
# Full_Organic  pH ser 1      3.37413E-09  2.81654E-09
# Full_Organic  pH ser 2      1.30784E-09  6.52926E-10
# Full_Organic  pH ser 3      4.26151E-09  4.40352E-10
# Full_Organic  pH ser 4      1.64310E-08  4.15840E-10
# Full_Organic  pH ser 5      4.80817E-08  3.51482E-10
# Full_Organic  pH ser 6      1.15698E-07  1.37217E-10
# Full_Organic  pH ser 7      2.22631E-07  1.81455E-11

# With TempCorrection (15C), Activity Correction (Debye & Davies), and the start
# of Double layer corrections
# (7 seconds)
#              ID        ID2         T.Cu           Cu  FinalIter FinalMaxError
# 1  Full_Organic Hard ser 1 3.358391e-09 1.054186e-10          4  1.338656e-05
# 2  Full_Organic Hard ser 2 4.182023e-09 1.375034e-10          3  3.048735e-04
# 3  Full_Organic Hard ser 3 6.290904e-09 2.281302e-10          3  1.240858e-04
# 4  Full_Organic Hard ser 4 9.181675e-09 3.683168e-10          3  1.202428e-04
# 5  Full_Organic Hard ser 5 1.382160e-08 6.282487e-10          4  2.484206e-05
# 6  Full_Organic Hard ser 6 2.408016e-08 1.335511e-09          4  6.904417e-04
# 7  Full_Organic Hard ser 7 3.674944e-08 2.415648e-09          4  3.321452e-04
# 8  Full_Organic  DOC ser 1 1.406313e-08 3.773789e-10          4  5.387245e-05
# 9  Full_Organic  DOC ser 2 3.182102e-08 3.773075e-10          3  2.833375e-04
# 10 Full_Organic  DOC ser 3 5.401897e-08 3.772714e-10          3  2.983766e-04
# 11 Full_Organic  DOC ser 4 9.841585e-08 3.772380e-10          3  3.261347e-04
# 12 Full_Organic  DOC ser 5 2.316583e-07 3.771825e-10          3  5.502643e-04
# 13 Full_Organic  DOC ser 6 3.428767e-07 3.772230e-10          3  8.779916e-04
# 14 Full_Organic  DOC ser 7 4.721743e-07 3.848036e-10         31  1.049695e-03
# 15 Full_Organic  DOC ser 8 7.404287e-07 3.948084e-10         31  2.573287e-03
# 16 Full_Organic  DOC ser 9 1.065966e-06 4.112278e-10         31  5.511362e-03
# 17 Full_Organic   pH ser 1 3.137583e-09 2.632134e-09          6  4.508738e-04
# 18 Full_Organic   pH ser 2 1.223810e-09 6.283966e-10          6  3.160180e-04
# 19 Full_Organic   pH ser 3 3.060127e-09 4.281447e-10          5  9.768868e-04
# 20 Full_Organic   pH ser 4 8.032604e-09 4.006754e-10          5  4.333127e-04
# 21 Full_Organic   pH ser 5 2.477155e-08 3.348249e-10          6  4.865446e-04
# 22 Full_Organic   pH ser 6 1.418089e-35 7.777752e-43         31  4.887548e+04
# 23 Full_Organic   pH ser 7 8.283453e-12 1.000916e-31         31  3.533566e+15




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

# With TempCorrection (15C) and Activity Correction (Debye only)
#   Obs           ID        ID2         T.Cu           Cu
# 1   1 Full_Organic Hard ser 1 2.672889e-09 1.054182e-10
# 2   2 Full_Organic Hard ser 2 3.374221e-09 1.375028e-10
# 3   3 Full_Organic Hard ser 3 5.221387e-09 2.281308e-10
# 4   4 Full_Organic Hard ser 4 7.834061e-09 3.683141e-10
# 5   5 Full_Organic Hard ser 5 1.214413e-08 6.282100e-10
# 6   6 Full_Organic Hard ser 6 2.194183e-08 1.335361e-09
# 7   7 Full_Organic Hard ser 7 3.428998e-08 2.415803e-09


