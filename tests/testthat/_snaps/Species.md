# AddSpecies works

    Code
      AddSpecies(ThisProblem = carbonate_system_problem, SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
        SpecMCName = "Water", SpecActCorr = "Debye", SpecLogK = 1.23, SpecDeltaH = -
        1234, SpecTempKelvin = 298)
    Output
      $ParamFile
      [1] "C:/Users/kellyc/Documents/BLM Development/engine/BLMEngineInR/inst/extdata/ParameterFiles/carbonate_system_only.dat4 (modified)"
      
      $N
           Mass     InLab     InVar    InMass    InComp InDefComp    InSpec   DefComp 
              1         1         2         1         1         1         4         1 
           Comp      Spec     Phase        BL     Metal   BLMetal       CAT 
              2         6         0         0         0         0         0 
      
      $Mass
         Name Amt Unit
      1 Water   1    L
      
      $InLabName
      [1] "ID"
      
      $InVar
        Name MCName MCR        Type
      1 Temp  Water   1 Temperature
      2   pH  Water   1          pH
      
      $InMassName
      [1] "Water"
      
      $InCompName
      [1] "CO3"
      
      $InDefCompName
      [1] "H"
      
      $InSpecName
      [1] "OH"    "HCO3"  "H2CO3" "H3CO3"
      
      $DefComp
        Name FromNum FromVar Charge MCName MCR     Type ActCorr SiteDens
      1    H      NA      pH      1  Water   1 FixedAct   Debye        1
      
      $Comp
        Name Charge MCName MCR     Type ActCorr SiteDens
      1  CO3     -2  Water   1  MassBal   Debye        1
      2    H      1  Water   1 FixedAct   Debye        1
      
      $Spec
         Name                Equation Charge MCName MCR   Type ActCorr    LogK
      1   CO3           CO3 = 1 * CO3     -2  Water   1 Normal   Debye   0.000
      2     H               H = 1 * H      1  Water   1 Normal   Debye   0.000
      3    OH             OH = -1 * H     -1  Water   1 Normal   Debye -13.997
      4  HCO3  HCO3 = 1 * CO3 + 1 * H     -1  Water   1 Normal   Debye  10.329
      5 H2CO3 H2CO3 = 1 * CO3 + 2 * H      0  Water   1 Normal   Debye  16.681
      6 H3CO3 H3CO3 = 1 * CO3 + 3 * H      1  Water   1 Normal   Debye   1.230
                   K DeltaH TempKelvin NC
      1 1.000000e+00      0       0.00  1
      2 1.000000e+00      0       0.00  1
      3 1.006932e-14 -55810     298.15  1
      4 2.133045e+10 -14600     298.15  2
      5 4.797334e+16 -23760     298.15  2
      6 1.698244e+01  -1234     298.00  2
      
      $SpecCompList
            [,1] [,2]
      CO3      1    0
      H        2    0
      OH       2    0
      HCO3     1    2
      H2CO3    1    2
      H3CO3    1    2
      
      $SpecStoich
            CO3  H
      CO3     1  0
      H       0  1
      OH      0 -1
      HCO3    1  1
      H2CO3   1  2
      H3CO3   1  3
      
      $Phase
      [1] Name       Equation   NC         LogK       DeltaH     TempKelvin Moles     
      <0 rows> (or 0-length row.names)
      
      $PhaseCompList
      <0 x 0 matrix>
      
      $PhaseStoich
           CO3 H
      
      $BL
      [1] Name  CompR
      <0 rows> (or 0-length row.names)
      
      $Metal
      [1] Name  CompR
      <0 rows> (or 0-length row.names)
      
      $BLMetal
      [1] Name   SpecsR
      <0 rows> (or 0-length row.names)
      
      $CATab
       [1] Num           CA            Species       Test.Type     Duration     
       [6] Lifestage     Endpoint      Quantifier    References    Miscellaneous
      <0 rows> (or 0-length row.names)
      
      $WHAM
      $WHAM$DoWHAM
      [1] FALSE
      
      $WHAM$WHAMVer
      [1] NA
      
      $WHAM$WdatFile
      [1] NA
      
      $WHAM$wDLF
      [1] NA
      
      $WHAM$wKZED
      [1] NA
      
      $WHAM$wP
      HA FA 
      NA NA 
      
      $WHAM$wRadius
      HA FA 
      NA NA 
      
      $WHAM$wMolWt
      HA FA 
      NA NA 
      
      
      $Index
      $Index$AqueousMCR
      [1] 1
      
      $Index$BioticLigMCR
      [1] NA
      
      $Index$WHAMDonnanMCR
      HA FA 
      -1 -1 
      
      

# RemoveSpecies works

    Code
      RemoveSpecies(ThisProblem = carbonate_system_problem, SpeciesToRemove = "HCO3")
    Output
      $ParamFile
      [1] "C:/Users/kellyc/Documents/BLM Development/engine/BLMEngineInR/inst/extdata/ParameterFiles/carbonate_system_only.dat4 (modified)"
      
      $N
           Mass     InLab     InVar    InMass    InComp InDefComp    InSpec   DefComp 
              1         1         2         1         1         1         2         1 
           Comp      Spec     Phase        BL     Metal   BLMetal       CAT 
              2         4         0         0         0         0         0 
      
      $Mass
         Name Amt Unit
      1 Water   1    L
      
      $InLabName
      [1] "ID"
      
      $InVar
        Name MCName MCR        Type
      1 Temp  Water   1 Temperature
      2   pH  Water   1          pH
      
      $InMassName
      [1] "Water"
      
      $InCompName
      [1] "CO3"
      
      $InDefCompName
      [1] "H"
      
      $InSpecName
      [1] "OH"    "H2CO3"
      
      $DefComp
        Name FromNum FromVar Charge MCName MCR     Type ActCorr SiteDens
      1    H      NA      pH      1  Water   1 FixedAct   Debye        1
      
      $Comp
        Name Charge MCName MCR     Type ActCorr SiteDens
      1  CO3     -2  Water   1  MassBal   Debye        1
      2    H      1  Water   1 FixedAct   Debye        1
      
      $Spec
         Name                Equation Charge MCName MCR   Type ActCorr    LogK
      1   CO3           CO3 = 1 * CO3     -2  Water   1 Normal   Debye   0.000
      2     H               H = 1 * H      1  Water   1 Normal   Debye   0.000
      3    OH             OH = -1 * H     -1  Water   1 Normal   Debye -13.997
      4 H2CO3 H2CO3 = 1 * CO3 + 2 * H      0  Water   1 Normal   Debye  16.681
                   K DeltaH TempKelvin NC
      1 1.000000e+00      0       0.00  1
      2 1.000000e+00      0       0.00  1
      3 1.006932e-14 -55810     298.15  1
      4 4.797334e+16 -23760     298.15  2
      
      $SpecCompList
            [,1] [,2]
      CO3      1    0
      H        2    0
      OH       2    0
      H2CO3    1    2
      
      $SpecStoich
            CO3  H
      CO3     1  0
      H       0  1
      OH      0 -1
      H2CO3   1  2
      
      $Phase
      [1] Name       Equation   NC         LogK       DeltaH     TempKelvin Moles     
      <0 rows> (or 0-length row.names)
      
      $PhaseCompList
      <0 x 0 matrix>
      
      $PhaseStoich
           CO3 H
      
      $BL
      [1] Name  CompR
      <0 rows> (or 0-length row.names)
      
      $Metal
      [1] Name  CompR
      <0 rows> (or 0-length row.names)
      
      $BLMetal
      [1] Name   SpecsR
      <0 rows> (or 0-length row.names)
      
      $CATab
       [1] Num           CA            Species       Test.Type     Duration     
       [6] Lifestage     Endpoint      Quantifier    References    Miscellaneous
      <0 rows> (or 0-length row.names)
      
      $WHAM
      $WHAM$DoWHAM
      [1] FALSE
      
      $WHAM$WHAMVer
      [1] NA
      
      $WHAM$WdatFile
      [1] NA
      
      $WHAM$wDLF
      [1] NA
      
      $WHAM$wKZED
      [1] NA
      
      $WHAM$wP
      HA FA 
      NA NA 
      
      $WHAM$wRadius
      HA FA 
      NA NA 
      
      $WHAM$wMolWt
      HA FA 
      NA NA 
      
      
      $Index
      $Index$AqueousMCR
      [1] 1
      
      $Index$BioticLigMCR
      [1] NA
      
      $Index$WHAMDonnanMCR
      HA FA 
      -1 -1 
      
      

