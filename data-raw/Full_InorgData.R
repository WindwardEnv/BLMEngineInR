# Create a set of test chemistry

# Total concentrations for components
Full_InorgDataTotalConc = c(
  T.H   =  2.692990E-08,
  T.Cu  =  4.720990E-08,
  T.DOC =  2.000000E-09,
  T.Ca  =  2.994160E-04,
  T.Mg  =  5.019540E-04,
  T.Na  =  1.100490E-03,
  T.K   =  5.371080E-05,
  T.SO4 =  7.994870E-04,
  T.Cl  =  5.359210E-05,
  T.CO3 =  1.099870E-03,
  T.BL  =  1.780000E-05
)
usethis::use_data(Full_InorgDataTotalConc, overwrite = T)

# Free ion species concentrations
Full_InorgDataFreeConc = c(
  H        =  2.692990E-08,
  Cu       =  2.076250E-09,
  DOC      =  2.000000E-09,
  Ca       =  2.736991E-04,
  Mg       =  4.566249E-04,
  Na       =  1.100490E-03,
  K        =  5.371080E-05,
  SO4      =  7.362015E-04,
  Cl       =  5.359210E-05,
  CO3      =  2.100791E-06,
  BL       =  1.302789E-10,
  `BL-Cu`   =  4.779432E-12,
  `BL-CuOH` =  3.866660E-13,
  `BL-Ca`   =  9.985501E-11,
  `BL-Mg`   =  1.665927E-10,
  `BL-H`    =  8.070763E-13,
  `BL-Na`   =  1.313004E-10,
  HCO3     =  1.031935E-03,
  H2CO3    =  5.803086E-05,
  MgHCO3   =  4.166637E-06,
  MgCO3    =  4.994100E-07,
  MgSO4    =  4.066301E-05,
  CaHCO3   =  2.587318E-06,
  CaCO3    =  5.073102E-07,
  CaSO4    =  2.262226E-05,
  CuOH     =  1.381260E-09,
  `Cu(OH)2` =  6.870140E-11,
  CuSO4    =  1.944034E-10,
  CuCO3    =  1.449888E-08,
  `Cu(CO3)2`=  4.505236E-11,
  CuCl     =  2.048608E-13,
  CuHCO3   =  2.894471E-08
)
usethis::use_data(Full_InorgDataFreeConc, overwrite = T)

# Reaction stoichiometry
Full_InorgDataStoich = matrix(
  data = c(
    #H    Cu  DOC   Ca   Mg   Na   K    SO4  Cl  CO3   BL
    1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  #H        = H
    0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,  #Cu       = Cu
    0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,  #DOC      = DOC
    0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,  #Ca       = Ca
    0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,  #Mg       = Mg
    0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,  #Na       = Na
    0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,  #K        = K
    0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,  #SO4      = SO4
    0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,  #Cl       = Cl
    0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,  #CO3      = CO3
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,  #BL       = BL
    0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   1,  #BL-Cu    = BL + Cu
    -1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   1,  #BL-CuOH  = BL + Cu - H
    0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   1,  #BL-Ca    = BL + Ca
    0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   1,  #BL-Mg    = BL + Mg
    1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,  #BL-H     = BL + H
    0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   1,  #BL-Na    = BL + Na
    1,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,  #HCO3     = H + CO3
    2,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,  #H2CO3    = 2*H + CO3
    1,   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,  #MgHCO3   = H + Mg + CO3
    0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   0,  #MgCO3    = Mg + CO3
    0,   0,   0,   0,   1,   0,   0,   1,   0,   0,   0,  #MgSO4    = Mg + SO4
    1,   0,   0,   1,   0,   0,   0,   0,   0,   1,   0,  #CaHCO3   = H + Ca + CO3
    0,   0,   0,   1,   0,   0,   0,   0,   0,   1,   0,  #CaCO3    = Ca + CO3
    0,   0,   0,   1,   0,   0,   0,   1,   0,   0,   0,  #CaSO4    = Ca + SO4
    -1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,  #CuOH     = Cu - H
    -2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,  #Cu(OH)2  = Cu - 2*H
    0,   1,   0,   0,   0,   0,   0,   1,   0,   0,   0,  #CuSO4    = Cu + SO4
    0,   1,   0,   0,   0,   0,   0,   0,   0,   1,   0,  #CuCO3    = Cu + CO3
    0,   1,   0,   0,   0,   0,   0,   0,   0,   2,   0,  #Cu(CO3)2 = Cu + 2*CO3
    0,   1,   0,   0,   0,   0,   0,   0,   1,   0,   0,  #CuCl     = Cu + Cl
    1,   1,   0,   0,   0,   0,   0,   0,   0,   1,   0   #CuHCO3   = H + Cu + CO3
  ), byrow = T, nrow = 32, ncol = 11,
  dimnames = list(Spec = c("H","Cu","DOC","Ca","Mg","Na","K","SO4","Cl","CO3","BL",
                           "BL-Cu" ,"BL-CuOH" ,"BL-Ca" ,"BL-Mg" ,"BL-H" ,"BL-Na" ,
                           "HCO3" ,"H2CO3" ,"MgHCO3" ,"MgCO3" ,"MgSO4" ,"CaHCO3" ,
                           "CaCO3" ,"CaSO4" ,"CuOH" ,"Cu(OH)2" ,"CuSO4" ,"CuCO3" ,
                           "Cu(CO3)2" ,"CuCl" ,"CuHCO3"),
                  Comp = c("H","Cu","DOC","Ca","Mg","Na","K","SO4","Cl","CO3","BL"))
)
usethis::use_data(Full_InorgDataStoich, overwrite = T)

# Reaction Equilibrium constants
Full_InorgDataK = 10^c(
  H        =    0.0,
  Cu       =    0.0,
  DOC      =    0.0,
  Ca       =    0.0,
  Mg       =    0.0,
  Na       =    0.0,
  K        =    0.0,
  SO4      =    0.0,
  Cl       =    0.0,
  CO3      =    0.0,
  BL       =    0.0,
  `BL-Cu`   =    7.4,
  `BL-CuOH` =   -1.3,
  `BL-Ca`   =    3.6,
  `BL-Mg`   =    3.6,
  `BL-H`    =    5.4,
  `BL-Na`   =    3.0,
  HCO3     =   10.329,
  H2CO3    =   16.681,
  MgHCO3   =   11.4,
  MgCO3    =   2.98,
  MgSO4    =   2.37,
  CaHCO3   =  11.44,
  CaCO3    =   3.22,
  CaSO4    =   2.3,
  CuOH     =  -7.52,
  `Cu(OH)2` = -16.22,
  CuSO4    =   2.36,
  CuCO3    =   6.75,
  `Cu(CO3)2`=   9.92,
  CuCl     =   0.4,
  CuHCO3   =  14.62
)
usethis::use_data(Full_InorgDataK, overwrite = T)
