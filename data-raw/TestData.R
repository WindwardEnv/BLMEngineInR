# Create a set of test chemistry

# Total concentrations for components
TestDataTotalConc = c(
  T.H = 1.18332E-04,
  T.CO3 = 1.00000E-04
)
usethis::use_data(TestDataTotalConc, overwrite = T)

# Free ion species concentrations
TestDataFreeConc = c(
  H = 1.00748E-07,
  CO3 = 3.90357E-08,
  HCO3 = 8.16805E-05,
  H2CO3 = 1.82705E-05
)
usethis::use_data(TestDataFreeConc, overwrite = T)

# Reaction stoichiometry
TestDataStoich = matrix(
  data = c(
    # H  CO3
    1, 0,  #H = +1*H
    0, 1,  #CO3 = +1*CO3
    1, 1,  #HCO3 = +1*H +1*CO3
    2, 1   #H2CO3 = +2*H +1*CO3
  ), byrow = T, nrow = 4, ncol = 2,
  dimnames = list(Spec = c("H","CO3","HCO3","H2CO3"),
                  Comp = c("H","CO3"))
)
usethis::use_data(TestDataStoich, overwrite = T)

# Reaction Equilibrium constants
TestDataK = 10^c(
  H = 0,
  CO3 = 0,
  HCO3 = 10.329,
  H2CO3 = 16.686
)
usethis::use_data(TestDataK, overwrite = T)
