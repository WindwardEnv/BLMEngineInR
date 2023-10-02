# Create a set of test chemistry

# Total concentrations for components
TestDataTotalConc = c(
  T.H = 1.18332E-04,
  T.CO3 = 1.00000E-04
)
save(TestDataTotalConc, file = "data/TestDataTotalConc.RData")

# Free ion species concentrations
TestDataFreeConc = c(
  H = 1.00748E-07,
  CO3 = 3.90357E-08,
  HCO3 = 8.16805E-05,
  H2CO3 = 1.82705E-05
)
save(TestDataFreeConc, file = "data/TestDataFreeConc.RData")

# Reaction stoichiometry
TestDataStoich = matrix(c(
  # H  CO3
  1, 0,  #H = +1*H
  0, 1,  #CO3 = +1*CO3
  1, 1,  #HCO3 = +1*H +1*CO3
  2, 1   #H2CO3 = +2*H +1*CO3
), byrow = T, ncol = 2,
dimnames = list(Spec=c("H","CO3","HCO3","H2CO3"), Comp = c("H","CO3")))
save(TestDataStoich, file = "data/TestDataStoich.RData")

# Reaction Equilibrium constants
TestDataK = 10^c(
  H = 1,
  CO3 = 1,
  HCO3 = 10.329,
  H2CO3 = 16.686
)
save(TestDataK, file = "data/TestDataK.RData")
