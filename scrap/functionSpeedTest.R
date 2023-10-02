load("data/testData.RData")
source("R/KellyFunctions.R")
Rcpp::sourceCpp("scrap/KellyFunctions.cpp")
source("R/FKellyFunctions.R")

stoichV = as.vector(t(TestDataStoich[3:4,]))

time.rec = data.frame(KellyR = rep(NA, 10000), KellyCpp = NA, KellyF = NA)
for(i in 1:nrow(time.rec)){
  start.time = Sys.time()
  CalcSpecConc(CConc = TestDataFreeConc[1:2], K = TestDataK[3:4], Stoich = TestDataStoich[3:4,])
  time.1 = Sys.time()
  CppCalcSpecConc(CConc = TestDataFreeConc[1:2], K = TestDataK[3:4], Stoich = stoichV,
                  NComp = 2, NSpec = 2)
  time.2 = Sys.time()
  FCalcSpeciesConc(CConc = TestDataFreeConc[1:2], K = TestDataK[3:4],
                   Stoich = TestDataStoich[3:4,], NComp = 2, NSpec = 2)
  time.3 = Sys.time()

  time.rec$KellyR[i] = time.1 - start.time
  time.rec$KellyCpp[i] = time.2 - time.1
  time.rec$KellyF[i] = time.3 - time.2
}
summary(time.rec * 10^6)#microseconds
#     KellyR           KellyCpp          KellyF
# Min.   :  8.106   Min.   : 4.768   Min.   : 6.914
# 1st Qu.: 10.014   1st Qu.: 5.960   1st Qu.: 7.868
# Median : 10.967   Median : 5.960   Median : 8.106
# Mean   : 11.830   Mean   : 6.702   Mean   : 9.243
# 3rd Qu.: 11.206   3rd Qu.: 6.914   3rd Qu.: 9.060
# Max.   :450.134   Max.   :83.923   Max.   :50.068
apply(time.rec*10^6, MARGIN = 2, FUN = quantile, probs = c(0.05,0.25,0.5,0.75,0.95))
#        KellyR  KellyCpp    KellyF
# 5%   9.059906  5.006790  7.152557
# 25% 10.013580  5.960464  7.867813
# 50% 10.967255  5.960464  8.106232
# 75% 11.205673  6.914139  9.059906
# 95% 20.980835 10.967255 15.974045

# KellyR   0        .-|-.........                                            +++
# KellyCpp 0     |-....                                                        +
# KellyF   0      |-.......                                   0




