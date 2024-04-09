rm(list = ls())
# devtools::clean_dll()
devtools::load_all()


ParamFile = "inst/extdata/ParameterFiles/full_inorg.dat4"
InputFile = "inst/extdata/InputFiles/full_inorg.blm4"

# ParamFile = "scrap/parameter file format/abbrev_organic.dat4"
# # ParamFile = "scrap/parameter file format/abbrev_organic (2).dat4"
# InputFile = "scrap/parameter file format/abbrev_organic.blm4"

# ParamFile = "scrap/parameter file format/full_organic_WATER23dH.dat4"
# InputFile = "scrap/parameter file format/full_organic.blm4"

# ParamFile = "scrap/parameter file format/full_organic_WATER23dH_FixedConcComps.dat4"
# InputFile = "scrap/parameter file format/full_organic_FixedConcComps.blm4"

DoTox = F
iCA = 1L
QuietFlag ="Debug"
ConvergenceCriteria = 0.0001
MaxIter = 100L
DoPartialStepsAlways = FALSE


ThisProblem = DefineProblem(ParamFile, WriteLog = TRUE)

FunctionInputs = ThisProblem[
  which(names(ThisProblem) %in% formalArgs("GetData"))]
FunctionInputs$InputFile = InputFile
AllInput = do.call("GetData", args = FunctionInputs)

# test stuff
start.time = Sys.time()
# capture.output(
ResultsTable <- BLM(
  ParamFile = ParamFile,
  InputFile = InputFile,
  DoTox = F,
  # iCA = 1L,
  QuietFlag ="Debug",
  ConvergenceCriteria = 0.0001,
  MaxIter = 100L,
  DoPartialStepsAlways = FALSE
)
# , file = "scrap/debug.txt")
end.time = Sys.time()
end.time - start.time
#
ResultsTable$Hard = (ResultsTable$Input.Ca + ResultsTable$Input.Mg) * 100086
#
# ResultsTable[, c("Obs","ID2","Hard","pH","DOC")]
ResultsTable[, c("Obs","ID2","FinalIter","FinalMaxError")]
# iObs = which((ResultsTable$FinalIter < 30) & !is.na(ResultsTable$FinalMaxError))[1]
# ResultsTable[, c("ID","ID2","T.Cu (mol/L)","Cu (mol/L)")]
# ResultsTable[, c("ID2","T.Cu (mol/L)","T.Cu (mol)","Water (L)", "Water_DonnanHA (L)","Water_DonnanFA (L)")]





OldBLMResultsTable = openxlsx::read.xlsx(xlsxFile = "scrap/old BLM/full_organic_SPEC.det.xlsx",
                                         sheet = 1, rows = c(5, 7:55))
OldBLMResultsTable$Obs = 1:nrow(OldBLMResultsTable)
OldBLMResultsTable$ID = trimws(gsub("\"","", OldBLMResultsTable$Site.Label))
OldBLMResultsTable$ID2 = trimws(gsub("\"","", OldBLMResultsTable$Sample.Label))
colnames(OldBLMResultsTable)[paste0(colnames(OldBLMResultsTable), " (mol/L)") %in% colnames(ResultsTable)] =
  paste0(colnames(OldBLMResultsTable)[paste0(colnames(OldBLMResultsTable), " (mol/L)") %in% colnames(ResultsTable)], " (mol/L)")
OldBLMResultsTable$`DonnanHA (mol/L)` = OldBLMResultsTable$Ratio_HS
OldBLMResultsTable$`DonnanFA (mol/L)` = OldBLMResultsTable$Ratio_FS
OldBLMResultsTable$Z_HA = OldBLMResultsTable$Z_HS
OldBLMResultsTable$Z_FA = OldBLMResultsTable$Z_FS

# ResultsTable[iObs, c("IonicStrength", "Water_DonnanHA (L)","Water_DonnanFA (L)",
#                      "DonnanHA (mol/L)","DonnanFA (mol/L)")]
# OldBLMResultsTable[iObs, c("Ionic.S.", "DDLMaxV_HS", "DDLMaxV_FS", "Ratio_HS", "Ratio_FS")]
#
# IS = 0
# for (iSpec in 1:ThisProblem$NSpec){
#   tmp = ResultsTable[iObs, paste0(ThisProblem$SpecName[iSpec],
#                                   " (mol/", ThisProblem$MassUnit[ThisProblem$SpecMCR[iSpec]], ")")] *
#     ThisProblem$SpecCharge[iSpec] ^ 2
#   # if (iSpec %in% c(11, 75:91)) {
#   #   tmp = tmp * 2.598275e-05
#   # } else if (iSpec %in% c(12, 92:108)) {
#   #   tmp = tmp * 0.001940053
#   # }
#   IS = IS + tmp
# }
# (IS = IS * 0.5)
#
# OldBLMResultsTable$CuHCO3[iObs]
# ResultsTable$`CuHCO3 (mol/L)`[iObs]
# ResultsTable$Act.CuHCO3[iObs]
#
# OldBLMResultsTable$CO3[iObs]
# ResultsTable$`CO3 (mol/L)`[iObs]
#
# OldBLMResultsTable$Cu[iObs]
# ResultsTable$`Cu (mol/L)`[iObs]
#
# ResultsTable[iObs,c("ID2","Input.Cu","Cu (mol/L)","CO3 (mol/L)","H (mol/L)","CuHCO3 (mol/L)")]
# OldBLMResultsTable[iObs,c("ID2","Cu","CO3","H","CuHCO3")]
#
# log10((ResultsTable$`Cu (mol/L)` * ResultsTable$`CO3 (mol/L)` * ResultsTable$`H (mol/L)`) / ResultsTable$`CuHCO3 (mol/L)`)[iObs]
# log10((OldBLMResultsTable$Cu * OldBLMResultsTable$CO3 * OldBLMResultsTable$H) / OldBLMResultsTable$CuHCO3)[iObs]
#
# OldBLMResultsTable$TOrg.Cu[iObs]
# ResultsTable$`TOrg.Cu (mol/L)`[iObs]

# par(mar = c(5,15,1,1))
# barplot(as.matrix(ResultsTable[iObs, Org.Cu.cols[grepl("Donnan", Org.Cu.cols)]]), horiz = TRUE, las = 2)
# abline(v = OldBLMResultsTable$TOrg.Cu[iObs], col = "red")
#
# (OldBLMResultsTable$Cu * OldBLMResultsTable$Act_z2)[iObs]
# ResultsTable$Act.Cu[iObs]
# (ResultsTable$Act.Cu / ResultsTable$`Cu (mol/L)`)[iObs]
# OldBLMResultsTable$Act_z2[iObs]
#
# ResultsTable$T.Cu = ResultsTable$`T.Cu (mol/L)`
# ResultsTable$Cu = ResultsTable$`Cu (mol/L)`
# (ResultsTable$DDL.Na = ResultsTable$`DonnanHA-Na (mol/L)` + ResultsTable$`DonnanFA-Na (mol/L)`)
# (OldBLMResultsTable$DDL.Na = OldBLMResultsTable$T.Na - OldBLMResultsTable$Na)
# ResultsTable[, c("ID","ID2","T.Cu","Cu","TOrg.Cu","DDL.Na")]
# OldBLMResultsTable[, c("ID","ID2","T.Cu","Cu","TOrg.Cu","DDL.Na")]
# OldBLMResultsTable$Act.OH = OldBLMResultsTable$OH


ResultsTable$SerLab = gsub(".+ ([0-9.]+)$", "\\1", ResultsTable$ID2)
ResultsTable$Ser = gsub("ser .+","ser", ResultsTable$ID2)
leg.dat = data.frame(
  Ser = c("Hard ser", "DOC ser", "pH ser", "Temp ser",
          "hi DOC Hard ser", "hi DOC pH ser", "hi DOC Temp ser"),
  col = 2:8
)
ResultsTable$col = leg.dat$col[match(ResultsTable$Ser, leg.dat$Ser)]

for (i.col in intersect(colnames(ResultsTable),
                        colnames(OldBLMResultsTable))) {
  if (all(!is.na(is.numeric(ResultsTable[,i.col]))) &&
      !grepl("T[.]", i.col) &&
      (i.col %in% c("Obs","ID","ID2","#.Iter.", "DDL.Na") == FALSE)) {

    ResultsTable[is.infinite(ResultsTable[, i.col]), i.col] = NA

    ax.lim = range(OldBLMResultsTable[, i.col], na.rm = T)
    if (all(ax.lim > 0)) {
      ax.lim[1] = 10^(floor(log10(ax.lim[1])))
      ax.lim[2] = 10^(ceiling(log10(ax.lim[2])))
    } else {
      ax.lim = c(min(pretty(ax.lim)),
                 max(pretty(ax.lim)))
    }

    layout(matrix(1:2, ncol = 2), widths = c(1, lcm(3)))
    par(mar = c(5,5,2,1))
    plot(x = OldBLMResultsTable[,i.col],
         y = pmin(ax.lim[2], pmax(ax.lim[1], ResultsTable[,i.col])),
         col = ResultsTable$col,
         pch = ifelse(ResultsTable[, i.col] < ax.lim[1], 6,
                      ifelse(ResultsTable[, i.col] > ax.lim[2], 2, 16)),
         main = i.col,
         xlab = "BLM Version 3.41.2.45",
         ylab = "BLM in R",
         log = ifelse(any(ax.lim<=0), "","xy"),
         xlim = ax.lim,
         ylim = ax.lim,
         las = 2)
    grid(nx = 5, ny = 5)
    text(x = OldBLMResultsTable[, i.col],
         y = pmin(ax.lim[2], pmax(ax.lim[1], ResultsTable[, i.col])),
         labels = ResultsTable$SerLab,
         col = ResultsTable$col,
         pos = 3,
         cex = 0.5)
    abline(a = 0, b = 1)

    par(mar = c(0,0,2,0))
    plot.new()
    legend("topleft", legend = leg.dat$Ser, pch = 16, col = leg.dat$col)
  }
}


if (FALSE) {

  # log10(OldBLMResultsTable$Act.OH * ResultsTable$Act.H)[24:29]
  #
  # SysTempKelvin = (ResultsTable$Temp + 273.15)
  #
  # T1 = 1 / SysTempKelvin
  # T0 = 1 / 298.15
  # SpecK = 10^-14.03541
  # Rcon = 8.314
  # SpecDeltaH = 56190#-56197
  # T2 = SpecDeltaH * (T0 - T1) / Rcon
  # log10(SpecK * exp(T2))[24:29]
  #
  # OldBLMResultsTable = OldBLMResultsTable[!grepl("DOC ser", OldBLMResultsTable$ID2), ]
  # ResultsTable = ResultsTable[!grepl("DOC ser", ResultsTable$ID2), ]


  pdf("scrap/speciation comparison.pdf")
  par(omi = rep(1, 4))
  # for (i.col in c("TOrg.Cu", "Cu", "CO3", "HCO3","Act.OH", "CuOH","Cu(OH)2")){
  for (i.col in intersect(colnames(ResultsTable), colnames(OldBLMResultsTable))){
    if (all(!is.na(is.numeric(ResultsTable[,i.col]))) &&
        !grepl("T[.]", i.col) &&
        (i.col %in% c("ID","ID2","#.Iter.", "DDL.Na", "TOrg.Cu") == FALSE)){
      ax.lim = range(OldBLMResultsTable[, i.col], na.rm = T)
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
#
# plot(x = ResultsTable$DOC, y = ResultsTable$Cu, type = "o", ylim = c(10^-10, 10^0), log = "y")
# lines(x = ResultsTable$DOC, y = ResultsTable$`DonnanFA-Cu`, type = "o", col = 2)
#
# plot(x = ResultsTable$DOC, y = ResultsTable$Na, type = "o", ylim = c(10^-10, 10^0), log = "y")
# lines(x = ResultsTable$DOC, y = ResultsTable$`DonnanFA-Na`, type = "o", col = 2)
#
# ResultsTable$`DonnanFA-Cu` * c(3.87787e-08, 1, 1, 1, 1, 0.00152631)
#
#
#
# for (i.col in c("TOrg.Cu", "Cu", "CO3", "HCO3","Act.OH", "CuOH","Cu(OH)2")){
# # for (i.col in colnames(ResultsTable)){
#     if ((i.col %in% colnames(OldBLMResultsTable)) &&
#         all(!is.na(is.numeric(ResultsTable[,i.col]))) &&
#         !grepl("T[.]", i.col) &&
#         (i.col %in% c("ID","ID2","#.Iter.", "DDL.Na", "TOrg.Cu") == FALSE)){
#     ax.lim = range(OldBLMResultsTable[, i.col])
#     ax.lim[1] = 10^(floor(log10(ax.lim[1])))
#     ax.lim[2] = 10^(ceiling(log10(ax.lim[2])))
#     plot(x = OldBLMResultsTable[,i.col],
#          y = pmin(ax.lim[2], pmax(ax.lim[1], ResultsTable[,i.col])),
#          col = ResultsTable$col,
#          pch = ifelse(ResultsTable[,i.col] < ax.lim[1], 6,
#                       ifelse(ResultsTable[,i.col] > ax.lim[2], 2, 16)),
#          main = i.col, xlab = "BLM Version 3.41.2.45",
#          ylab = "BLM in R",
#          log = "xy", xlim = ax.lim, ylim = ax.lim, las = 2)
#     text(x = OldBLMResultsTable[,i.col],
#          y = pmin(ax.lim[2], pmax(ax.lim[1], ResultsTable[,i.col])),
#          labels = ResultsTable$SerLab, col = ResultsTable$col,
#          pos = 3, cex = 0.5)
#     abline(a = 0, b = 1)
#     legend("topleft", legend = leg.dat$Ser, pch = 16, col = leg.dat$col)
#   }
# }
#
#
# plot(
#   x = OldBLMResultsTable$T.DOC,
#   y = OldBLMResultsTable$TOrg.Cu / OldBLMResultsTable$T.Cu,
# )
#
#
# plot(x = OldBLMResultsTable$T.Cu, y = pmax(10^-9, ResultsTable$T.Cu),
#      col = c(rep(2, 7), rep(3, 9), rep(4, 7)),
#      pch = ifelse(ResultsTable$T.Cu < 10^-9, 6, 16),
#      main = "Total Cu", xlab = "Old BLM", ylab = "Current BLM in R",
#      log = "xy", xlim = c(1E-9, 1E-5), ylim = c(1E-9, 1E-5))
# text(x = OldBLMResultsTable$T.Cu, y = pmax(10^-9, ResultsTable$T.Cu),
#      labels = c(1:7, 1:9, 1:7), col = c(rep(2, 7), rep(3, 9), rep(4, 7)),
#      pos = 3, cex = 0.5)
# abline(a = 0, b = 1)
# legend("topleft", legend = c("Hard ser", "DOC ser", "pH ser"), pch = 16, col = 2:4)
# plot(x = OldBLMResultsTable$Cu, y = pmax(10^-11, ResultsTable$Cu),
#      col = c(rep(2, 7), rep(3, 9), rep(4, 7)),
#      pch = ifelse(ResultsTable$Cu < 10^-11, 6, 16),
#      main = "Free Cu", xlab = "Old BLM", ylab = "Current BLM in R",
#      log = "xy", xlim = c(1E-11, 1E-6), ylim = c(1E-11, 1E-6))
# text(x = OldBLMResultsTable$Cu, y = pmax(10^-11, ResultsTable$Cu),
#      labels = c(1:7, 1:9, 1:7), col = c(rep(2, 7), rep(3, 9), rep(4, 7)),
#      pos = 3, cex = 0.5)
# abline(a = 0, b = 1)
# legend("topleft", legend = c("Hard ser", "DOC ser", "pH ser"), pch = 16, col = 2:4)
#
# par(mfrow = c(1, 2))
# plot(x = OldBLMResultsTable$Na, y = pmax(10^-5, ResultsTable$Na),
#      col = c(rep(2, 7), rep(3, 9), rep(4, 7)),
#      pch = ifelse(ResultsTable$Na < 10^-5, 6, 16),
#      main = "Free Na", xlab = "Old BLM", ylab = "Current BLM in R",
#      log = "xy", xlim = c(1E-5, 1E-2), ylim = c(1E-5, 1E-2))
# text(x = OldBLMResultsTable$Na, y = pmax(10^-5, ResultsTable$Na),
#      labels = c(1:7, 1:9, 1:7), col = c(rep(2, 7), rep(3, 9), rep(4, 7)),
#      pos = 3, cex = 0.5)
# abline(a = 0, b = 1)
# legend("topleft", legend = c("Hard ser", "DOC ser", "pH ser"), pch = 16, col = 2:4)
# plot(x = OldBLMResultsTable$DDL.Na, y = pmax(10^-8, ResultsTable$DDL.Na),
#      col = c(rep(2, 7), rep(3, 9), rep(4, 7)),
#      pch = ifelse(ResultsTable$DDL.Na < 10^-8, 6, 16),
#      main = "Diffuse Bound Na", xlab = "Old BLM", ylab = "Current BLM in R",
#      log = "xy", xlim = c(1E-8, 1E-6), ylim = c(1E-8, 1E-6))
# text(x = OldBLMResultsTable$DDL.Na, y = pmax(10^-8, ResultsTable$DDL.Na),
#      labels = c(1:7, 1:9, 1:7), col = c(rep(2, 7), rep(3, 9), rep(4, 7)),
#      pos = 3, cex = 0.5)
# abline(a = 0, b = 1)
# legend("topleft", legend = c("Hard ser", "DOC ser", "pH ser"), pch = 16, col = 2:4)
#
#
#
# # Old BLM (3 seconds)
# # Site Label    Sample Label  Dis. Cu      Free Cu
# # Full_Organic  Hard ser 1    1.17628E-08  1.17741E-10
# # Full_Organic  Hard ser 2    1.15568E-08  1.49896E-10
# # Full_Organic  Hard ser 3    1.26915E-08  2.41360E-10
# # Full_Organic  Hard ser 4    1.50443E-08  3.86284E-10
# # Full_Organic  Hard ser 5    1.89049E-08  6.60671E-10
# # Full_Organic  Hard ser 6    2.70029E-08  1.42593E-09
# # Full_Organic  Hard ser 7    3.65246E-08  2.62630E-09
# # Full_Organic  DOC ser 1     3.04891E-08  3.93838E-10
# # Full_Organic  DOC ser 2     1.13690E-07  3.94705E-10
# # Full_Organic  DOC ser 3     2.18304E-07  3.96361E-10
# # Full_Organic  DOC ser 4     4.28458E-07  3.98179E-10
# # Full_Organic  DOC ser 5     1.06909E-06  4.03965E-10
# # Full_Organic  DOC ser 6     1.61405E-06  4.08591E-10
# # Full_Organic  DOC ser 7     2.16647E-06  4.12376E-10
# # Full_Organic  DOC ser 8     3.28985E-06  4.18182E-10
# # Full_Organic  DOC ser 9     4.43633E-06  4.22807E-10
# # Full_Organic  pH ser 1      3.37413E-09  2.81654E-09
# # Full_Organic  pH ser 2      1.30784E-09  6.52926E-10
# # Full_Organic  pH ser 3      4.26151E-09  4.40352E-10
# # Full_Organic  pH ser 4      1.64310E-08  4.15840E-10
# # Full_Organic  pH ser 5      4.80817E-08  3.51482E-10
# # Full_Organic  pH ser 6      1.15698E-07  1.37217E-10
# # Full_Organic  pH ser 7      2.22631E-07  1.81455E-11
#
# # With TempCorrection (15C), Activity Correction (Debye & Davies), and the start
# # of Double layer corrections
# # (7 seconds)
# #              ID        ID2         T.Cu           Cu  FinalIter FinalMaxError
# # 1  Full_Organic Hard ser 1 3.358391e-09 1.054186e-10          4  1.338656e-05
# # 2  Full_Organic Hard ser 2 4.182023e-09 1.375034e-10          3  3.048735e-04
# # 3  Full_Organic Hard ser 3 6.290904e-09 2.281302e-10          3  1.240858e-04
# # 4  Full_Organic Hard ser 4 9.181675e-09 3.683168e-10          3  1.202428e-04
# # 5  Full_Organic Hard ser 5 1.382160e-08 6.282487e-10          4  2.484206e-05
# # 6  Full_Organic Hard ser 6 2.408016e-08 1.335511e-09          4  6.904417e-04
# # 7  Full_Organic Hard ser 7 3.674944e-08 2.415648e-09          4  3.321452e-04
# # 8  Full_Organic  DOC ser 1 1.406313e-08 3.773789e-10          4  5.387245e-05
# # 9  Full_Organic  DOC ser 2 3.182102e-08 3.773075e-10          3  2.833375e-04
# # 10 Full_Organic  DOC ser 3 5.401897e-08 3.772714e-10          3  2.983766e-04
# # 11 Full_Organic  DOC ser 4 9.841585e-08 3.772380e-10          3  3.261347e-04
# # 12 Full_Organic  DOC ser 5 2.316583e-07 3.771825e-10          3  5.502643e-04
# # 13 Full_Organic  DOC ser 6 3.428767e-07 3.772230e-10          3  8.779916e-04
# # 14 Full_Organic  DOC ser 7 4.721743e-07 3.848036e-10         31  1.049695e-03
# # 15 Full_Organic  DOC ser 8 7.404287e-07 3.948084e-10         31  2.573287e-03
# # 16 Full_Organic  DOC ser 9 1.065966e-06 4.112278e-10         31  5.511362e-03
# # 17 Full_Organic   pH ser 1 3.137583e-09 2.632134e-09          6  4.508738e-04
# # 18 Full_Organic   pH ser 2 1.223810e-09 6.283966e-10          6  3.160180e-04
# # 19 Full_Organic   pH ser 3 3.060127e-09 4.281447e-10          5  9.768868e-04
# # 20 Full_Organic   pH ser 4 8.032604e-09 4.006754e-10          5  4.333127e-04
# # 21 Full_Organic   pH ser 5 2.477155e-08 3.348249e-10          6  4.865446e-04
# # 22 Full_Organic   pH ser 6 1.418089e-35 7.777752e-43         31  4.887548e+04
# # 23 Full_Organic   pH ser 7 8.283453e-12 1.000916e-31         31  3.533566e+15
#
#
# # No TempCorrection
# #   Obs           ID        ID2         T.Cu           Cu
# # 1   1 Full_Organic Hard ser 1 2.246596e-09 9.909165e-11
# # 2   2 Full_Organic Hard ser 2 2.881979e-09 1.284319e-10
# # 3   3 Full_Organic Hard ser 3 4.619184e-09 2.104946e-10
# # 4   4 Full_Organic Hard ser 4 7.121356e-09 3.337472e-10
# # 5   5 Full_Organic Hard ser 5 1.122678e-08 5.473762e-10
# # 6   6 Full_Organic Hard ser 6 2.030386e-08 1.057243e-09
# # 7   7 Full_Organic Hard ser 7 3.135992e-08 1.720029e-09
#
# # With TempCorrection (15C)
# #   Obs           ID        ID2         T.Cu           Cu
# # 1   1 Full_Organic Hard ser 1 2.683681e-09 9.905403e-11
# # 2   2 Full_Organic Hard ser 2 3.445260e-09 1.282940e-10
# # 3   3 Full_Organic Hard ser 3 5.522465e-09 2.096111e-10
# # 4   4 Full_Organic Hard ser 4 8.512210e-09 3.309170e-10
# # 5   5 Full_Organic Hard ser 5 1.342391e-08 5.397773e-10
# # 6   6 Full_Organic Hard ser 6 2.432579e-08 1.035269e-09
# # 7   7 Full_Organic Hard ser 7 3.764053e-08 1.677717e-09
#
# # With TempCorrection (15C) and Activity Correction (Debye only)
# #   Obs           ID        ID2         T.Cu           Cu
# # 1   1 Full_Organic Hard ser 1 2.672889e-09 1.054182e-10
# # 2   2 Full_Organic Hard ser 2 3.374221e-09 1.375028e-10
# # 3   3 Full_Organic Hard ser 3 5.221387e-09 2.281308e-10
# # 4   4 Full_Organic Hard ser 4 7.834061e-09 3.683141e-10
# # 5   5 Full_Organic Hard ser 5 1.214413e-08 6.282100e-10
# # 6   6 Full_Organic Hard ser 6 2.194183e-08 1.335361e-09
# # 7   7 Full_Organic Hard ser 7 3.428998e-08 2.415803e-09
#
#
#
# start.time = Sys.time()
# ResultsTable_DPSTRUE <- BLM(
#   ParamFile = "scrap/parameter file format/full_organic_WATER23dH.dat4",
#   InputFile = "scrap/parameter file format/full_organic.blm4",
#   DoTox = F,
#   QuietFlag ="Quiet",
#   ConvergenceCriteria = 0.001,
#   MaxIter = 1000L,
#   DoPartialStepsAlways = TRUE
# )
# time.elapsed_DPSTRUE = Sys.time() - start.time
#
# start.time = Sys.time()
# ResultsTable_DPSFALSE <- BLM(
#   ParamFile = "scrap/parameter file format/full_organic_WATER23dH.dat4",
#   InputFile = "scrap/parameter file format/full_organic.blm4",
#   DoTox = F,
#   QuietFlag ="Quiet",
#   ConvergenceCriteria = 0.001,
#   MaxIter = 1000L,
#   DoPartialStepsAlways = FALSE
# )
# time.elapsed_DPSFALSE = Sys.time() - start.time
#
# ResultsTable_DPSTRUE$Hard = (ResultsTable_DPSTRUE$`T.Ca (mol/L)` + ResultsTable_DPSTRUE$`T.Mg (mol/L)`) * 100086
# ResultsTable_DPSFALSE$Hard = (ResultsTable_DPSFALSE$`T.Ca (mol/L)` + ResultsTable_DPSFALSE$`T.Mg (mol/L)`) * 100086
#
# ResultsTable_DPSTRUE[, c("Obs","ID2","FinalIter","FinalMaxError","Cu (mol/L)")]
# ResultsTable_DPSFALSE[, c("Obs","ID2","FinalIter","FinalMaxError","Cu (mol/L)")]
#
# # hardness series = 1 - 7
# # DOC series = 8 - 16
# # pH series = 17 - 23
# # Temp series = 24 - 29
# # hi DOC hard series = 30 - 36
# # hi DOC pH series = 37 - 43
# # hi DOC Temp series = 44 - 49
#
# mean(ResultsTable_DPSTRUE$FinalIter / ResultsTable_DPSFALSE$FinalIter)
# sum(ResultsTable_DPSFALSE$FinalIter == 1001L)#3 fail to converge when doing normal N-R
# sum(ResultsTable_DPSTRUE$FinalIter == 1001L)#5 fail to converge with DPS
# which((ResultsTable_DPSFALSE$FinalIter == 1001L)) #        9    24    41
# which((ResultsTable_DPSTRUE$FinalIter == 1001L))  #  2  3  9       25 41
# #                                                  Hard    DOC  temp  hdp
#
#
# sum(is.na(ResultsTable_DPSFALSE$`Cu (mol/L)`))#9 fail to converge when doing normal N-R
# sum(is.na(ResultsTable_DPSTRUE$`Cu (mol/L)`))#7 fail to converge with DPS
# which(is.na(ResultsTable_DPSFALSE$`Cu (mol/L)`)) #  1  2  3  4 21    23 26    35   43
# which(is.na(ResultsTable_DPSTRUE$`Cu (mol/L)`))  #  1          21 22 23 26 29      43
# #                                                   Hard       pH       temp  hdh  hdp
#
#
# plot(x = ResultsTable_DPSTRUE$`Cu (mol/L)`, y = ResultsTable_DPSFALSE$`Cu (mol/L)`, log = "xy"); abline(a = 0, b = 1)
# mean(ResultsTable_DPSTRUE$`Cu (mol/L)` / ResultsTable_DPSFALSE$`Cu (mol/L)`, na.rm = TRUE)
# # The answer's converging to the same point, when both converge, but it's a
# # toss-up of whether DPS actually helps...some fail to converge, some explode,
# # but there are some that do these in either case (with some overlap, but some
# # unique cases).
#
# sum(ResultsTable_DPSTRUE$FinalIter > ResultsTable_DPSFALSE$FinalIter)  #28
# sum(ResultsTable_DPSTRUE$FinalIter < ResultsTable_DPSFALSE$FinalIter)  #14
# sum(ResultsTable_DPSTRUE$FinalIter == ResultsTable_DPSFALSE$FinalIter) # 7
# which(ResultsTable_DPSTRUE$FinalIter > ResultsTable_DPSFALSE$FinalIter)
# which(ResultsTable_DPSTRUE$FinalIter < ResultsTable_DPSFALSE$FinalIter)
# which(ResultsTable_DPSTRUE$FinalIter == ResultsTable_DPSFALSE$FinalIter)
# #  best hard------------------  DOC-------------------  pH------------------  temp-------------  hi DOC hard---------   hi DOC pH-----------   hi DOC temp------
# #  N-R     2  3  4                10 11 12 13 14 15     17    19 20 21 22        25 26              31          35 36   37    39 40       43   44    46 47 48 49
# #  DPS  1           5  6  7  8                             18                 24       27 28 29  30       33 34                                   45
# #  same                         9                   16                    23                           32                  38       41 42
# #       hard------------------  DOC-------------------  pH------------------  temp-------------  hi DOC hard---------   hi DOC pH-----------   hi DOC temp------
#
# plot(x = ResultsTable_DPSTRUE$FinalIter, y = ResultsTable_DPSFALSE$FinalIter, log = "xy"); abline(a = 0, b = 1)
# plot(x = ResultsTable_DPSTRUE$Obs, y = ResultsTable_DPSTRUE$FinalIter, log = "y", type = "o"); points(x = ResultsTable_DPSFALSE$Obs, y = ResultsTable_DPSFALSE$FinalIter, col = 2, type = "o")
# plot(x = ResultsTable_DPSTRUE$DOC, y = ResultsTable_DPSTRUE$FinalIter, log = "xy", type = "o"); points(x = ResultsTable_DPSFALSE$DOC, y = ResultsTable_DPSFALSE$FinalIter, col = 2, type = "o")
# plot(x = ResultsTable_DPSTRUE$pH, y = ResultsTable_DPSTRUE$FinalIter, log = "y", type = "o"); points(x = ResultsTable_DPSFALSE$pH, y = ResultsTable_DPSFALSE$FinalIter, col = 2, type = "o")
# plot(x = ResultsTable_DPSTRUE$Hard, y = ResultsTable_DPSTRUE$FinalIter, log = "xy", type = "o"); points(x = ResultsTable_DPSFALSE$Hard, y = ResultsTable_DPSFALSE$FinalIter, col = 2, type = "o")
#
#
# # No partial steps: Time difference of 0.911727 secs
# #   Obs       ID2 FinalIter FinalMaxError
# # 1   1 DOC ser 1       451  6.060930e-04
# # 2   2 DOC ser 2         7  9.989368e-06
# # 3   3 DOC ser 3         7  2.465376e-04
# # 4   4 DOC ser 4         6  2.031681e-04
# # 5   5 DOC ser 5         6  2.070989e-04
# # 6   6 DOC ser 6         6  9.463851e-04
# #
# # Time difference of 1.696869 mins
# #> summary(ResultsTable$FinalIter)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# # 5      11      23     127      98    1001
# #
# # Partial steps: Time difference of 0.4131219 secs
# #   Obs       ID2 FinalIter FinalMaxError
# # 1   1 DOC ser 1        82  3.539385e-04
# # 2   2 DOC ser 2         6  9.035190e-04
# # 3   3 DOC ser 3         7  4.125571e-05
# # 4   4 DOC ser 4         6  2.031681e-04
# # 5   5 DOC ser 5         6  2.070989e-04
# # 6   6 DOC ser 6         6  9.463851e-04
# #
# # Time difference of 2.615645 mins
# # > summary(ResultsTable$FinalIter)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# # 5.0    14.0    44.0   175.9   145.0  1001.0
