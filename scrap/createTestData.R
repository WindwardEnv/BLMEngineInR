#create a set of test chemistry
labels = c(
  "T.H",
  "T.CO3",
  "H",
  "CO3",
  "HCO3",
  "H2CO3"
)
conc = c(
  1.18332E-04,  #"T.H",
  1.00000E-04,  #"T.CO3",
  1.00748E-07,  #"H",
  3.90357E-08,  #"CO3",
  8.16805E-05,  #"HCO3",
  1.82705E-05   #"H2CO3",
)

stoic = matrix(c(
 # H  CO3
   1, 0,  #T.H
   0, 1,  #T.CO3
   1, 1,  #HCO3
   2, 1   #H2CO3
), byrow = T, ncol = 2)

logK = c(
  1,      #T.H
  1,      #T.CO3
  10.329, #HCO3
  16.686  #H2CO3
)

save(conc, labels, stoic, logK, file = "data/testData.RData")
