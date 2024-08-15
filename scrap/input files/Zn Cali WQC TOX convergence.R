load("scrap/input files/Zn_Cali_WQC.blm4_TOX.RData")
dim(ResultsTable)
ResultsTable$Status = ifelse((ResultsTable$FinalMaxError > 0.0001) | is.na(ResultsTable$FinalMaxError), "Not Converged", "Okay")
summary(ResultsTable$Status)
aggregate(ResultsTable$Status, by = list(ResultsTable$Status), FUN = length)
ResultsTable$Hard = (ResultsTable$Input.Ca + ResultsTable$Input.Mg) * 100086

par(mfrow = c(2, 2))
boxplot(Temp ~ Status, data = ResultsTable)
boxplot(pH ~ Status, data = ResultsTable)
boxplot(DOC ~ Status, data = ResultsTable)
boxplot(Hard ~ Status, data = ResultsTable)

layout(mat = matrix(c(2,3,4,5,1,1), nrow = 3, ncol = 2, byrow = TRUE), heights = c(1,1,lcm(1)))
par(mar = rep(0,4))
plot.new()
legend("center", legend = c("Okay","Not Converged"), fill = c("green", "red"), cex = 1.5, horiz = TRUE)
par(mar = c(5,5,3,1))
for (icol in c("Temp","pH","DOC","Hard")) {

  hist.stat.all = hist(ResultsTable[, icol], plot = FALSE)
  hist.stat.okay = hist(ResultsTable[ResultsTable$Status == "Okay", icol], breaks = hist.stat.all$breaks, plot = FALSE)
  hist.stat.nc = hist(ResultsTable[ResultsTable$Status == "Not Converged", icol], breaks = hist.stat.all$breaks, plot = FALSE)

  barplot(matrix(data = c(hist.stat.okay$counts / hist.stat.all$counts * 100,
                          hist.stat.nc$counts / hist.stat.all$counts * 100),
                 nrow = 2, ncol = length(hist.stat.all$counts), byrow = TRUE,
                 dimnames = list(c("Okay","Not Converged"),
                                 bins = paste0(hist.stat.all$breaks[-length(hist.stat.all$breaks)],
                                               "-", hist.stat.all$breaks[-1]))),
          beside = FALSE, col = c("green", "red"),
          main = icol, ylab = "%", ylim = c(0,100), las = 2)
  box()

}
