# AddSpecies works

    Code
      AddSpecies(ThisProblem = carbonate_system_problem, SpecEquation = "H3CO3 = 3 * H + 1 * CO3",
        SpecMCName = "Water", SpecActCorr = "Debye", SpecLogK = 1.23, SpecDeltaH = -
        1234, SpecTempKelvin = 298)
    Output
      $ParamFile
      [1] "C:\\Users\\kellyc\\Documents\\Projects\\BLM Development\\engine\\BLMEngineInR\\inst\\extdata\\ParameterFiles\\carbonate_system_only.dat4 (modified)"
      
      $N
           Mass     InLab     InVar    InMass    InComp InDefComp    InSpec   DefComp 
              1         1         2         1         1         2         3         2 
           Comp      Spec     Phase        BL     Metal   BLMetal       CAT 
              3         6         0         0         0         0         0 
      
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
      [1] "H"  "OH"
      
      $InSpecName
      [1] "HCO3"  "H2CO3" "H3CO3"
      
      $DefComp
        Name FromNum FromVar Charge MCName MCR     Type ActCorr SiteDens
      1    H      NA      pH      1  Water   1 FixedAct   Debye        1
      2   OH      NA    KW/H     -1  Water   1 FixedAct   Debye        1
      
      $Comp
        Name Charge MCName MCR     Type ActCorr SiteDens
      1    H      1  Water   1 FixedAct   Debye        1
      2   OH     -1  Water   1 FixedAct   Debye        1
      3  CO3     -2  Water   1  MassBal   Debye        1
      
      $Spec
         Name                Equation Charge MCName MCR   Type ActCorr   LogK
      1     H               H = 1 * H      1  Water   1 Normal   Debye  0.000
      2    OH             OH = 1 * OH     -1  Water   1 Normal   Debye  0.000
      3   CO3           CO3 = 1 * CO3     -2  Water   1 Normal   Debye  0.000
      4  HCO3  HCO3 = 1 * H + 1 * CO3     -1  Water   1 Normal   Debye 10.329
      5 H2CO3 H2CO3 = 2 * H + 1 * CO3      0  Water   1 Normal   Debye 16.681
      6 H3CO3 H3CO3 = 3 * H + 1 * CO3      1  Water   1 Normal   Debye  1.230
                   K DeltaH TempKelvin NC
      1 1.000000e+00      0       0.00  1
      2 1.000000e+00      0       0.00  1
      3 1.000000e+00      0       0.00  1
      4 2.133045e+10 -14600     298.15  2
      5 4.797334e+16 -23760     298.15  2
      6 1.698244e+01  -1234     298.00  2
      
      $SpecCompList
            [,1] [,2]
      H        1    0
      OH       2    0
      CO3      3    0
      HCO3     1    3
      H2CO3    1    3
      H3CO3    1    3
      
      $SpecStoich
            H OH CO3
      H     1  0   0
      OH    0  1   0
      CO3   0  0   1
      HCO3  1  0   1
      H2CO3 2  0   1
      H3CO3 3  0   1
      
      $Phase
      [1] Name       Equation   NC         LogK       DeltaH     TempKelvin Moles     
      <0 rows> (or 0-length row.names)
      
      $PhaseCompList
      <0 x 0 matrix>
      
      $PhaseStoich
           H OH CO3
      
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
       [1] Num           CA            Species       TestType      Duration     
       [6] Lifestage     Endpoint      Quantifier    References    Miscellaneous
      <0 rows> (or 0-length row.names)
      
      $DoWHAM
      [1] FALSE
      
      $WHAM
      $WHAM$Ver
      [1] NA
      
      $WHAM$File
      [1] NA
      
      $WHAM$DLF
      [1] NA
      
      $WHAM$KZED
      [1] NA
      
      $WHAM$nA
      HA FA 
      NA NA 
      
      $WHAM$pKA
      HA FA 
      NA NA 
      
      $WHAM$pKB
      HA FA 
      NA NA 
      
      $WHAM$dpKA
      HA FA 
      NA NA 
      
      $WHAM$dpKB
      HA FA 
      NA NA 
      
      $WHAM$fprB
      HA FA 
      NA NA 
      
      $WHAM$fprT
      HA FA 
      NA NA 
      
      $WHAM$dLK1A
      HA FA 
      NA NA 
      
      $WHAM$dLK1B
      HA FA 
      NA NA 
      
      $WHAM$P
      HA FA 
      NA NA 
      
      $WHAM$Radius
      HA FA 
      NA NA 
      
      $WHAM$MolWt
      HA FA 
      NA NA 
      
      $WHAM$MonodentTable
      [1] S          AbundDenom StrongWeak
      <0 rows> (or 0-length row.names)
      
      $WHAM$BidentTable
      [1] S1         S2         AbundDenom
      <0 rows> (or 0-length row.names)
      
      $WHAM$TridentTable
      [1] S1         S2         S3         AbundDenom
      <0 rows> (or 0-length row.names)
      
      $WHAM$MetalsTable
      [1] Metal  pKMAHA pKMAFA dLK2  
      <0 rows> (or 0-length row.names)
      
      $WHAM$SpecKselTable
      [1] Spec   KselHA KselFA
      <0 rows> (or 0-length row.names)
      
      
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
      [1] "C:\\Users\\kellyc\\Documents\\Projects\\BLM Development\\engine\\BLMEngineInR\\inst\\extdata\\ParameterFiles\\carbonate_system_only.dat4 (modified)"
      
      $N
           Mass     InLab     InVar    InMass    InComp InDefComp    InSpec   DefComp 
              1         1         2         1         1         2         1         2 
           Comp      Spec     Phase        BL     Metal   BLMetal       CAT 
              3         4         0         0         0         0         0 
      
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
      [1] "H"  "OH"
      
      $InSpecName
      [1] "H2CO3"
      
      $DefComp
        Name FromNum FromVar Charge MCName MCR     Type ActCorr SiteDens
      1    H      NA      pH      1  Water   1 FixedAct   Debye        1
      2   OH      NA    KW/H     -1  Water   1 FixedAct   Debye        1
      
      $Comp
        Name Charge MCName MCR     Type ActCorr SiteDens
      1    H      1  Water   1 FixedAct   Debye        1
      2   OH     -1  Water   1 FixedAct   Debye        1
      3  CO3     -2  Water   1  MassBal   Debye        1
      
      $Spec
         Name                Equation Charge MCName MCR   Type ActCorr   LogK
      1     H               H = 1 * H      1  Water   1 Normal   Debye  0.000
      2    OH             OH = 1 * OH     -1  Water   1 Normal   Debye  0.000
      3   CO3           CO3 = 1 * CO3     -2  Water   1 Normal   Debye  0.000
      4 H2CO3 H2CO3 = 2 * H + 1 * CO3      0  Water   1 Normal   Debye 16.681
                   K DeltaH TempKelvin NC
      1 1.000000e+00      0       0.00  1
      2 1.000000e+00      0       0.00  1
      3 1.000000e+00      0       0.00  1
      4 4.797334e+16 -23760     298.15  2
      
      $SpecCompList
            [,1] [,2]
      H        1    0
      OH       2    0
      CO3      3    0
      H2CO3    1    3
      
      $SpecStoich
            H OH CO3
      H     1  0   0
      OH    0  1   0
      CO3   0  0   1
      H2CO3 2  0   1
      
      $Phase
      [1] Name       Equation   NC         LogK       DeltaH     TempKelvin Moles     
      <0 rows> (or 0-length row.names)
      
      $PhaseCompList
      <0 x 0 matrix>
      
      $PhaseStoich
           H OH CO3
      
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
       [1] Num           CA            Species       TestType      Duration     
       [6] Lifestage     Endpoint      Quantifier    References    Miscellaneous
      <0 rows> (or 0-length row.names)
      
      $DoWHAM
      [1] FALSE
      
      $WHAM
      $WHAM$Ver
      [1] NA
      
      $WHAM$File
      [1] NA
      
      $WHAM$DLF
      [1] NA
      
      $WHAM$KZED
      [1] NA
      
      $WHAM$nA
      HA FA 
      NA NA 
      
      $WHAM$pKA
      HA FA 
      NA NA 
      
      $WHAM$pKB
      HA FA 
      NA NA 
      
      $WHAM$dpKA
      HA FA 
      NA NA 
      
      $WHAM$dpKB
      HA FA 
      NA NA 
      
      $WHAM$fprB
      HA FA 
      NA NA 
      
      $WHAM$fprT
      HA FA 
      NA NA 
      
      $WHAM$dLK1A
      HA FA 
      NA NA 
      
      $WHAM$dLK1B
      HA FA 
      NA NA 
      
      $WHAM$P
      HA FA 
      NA NA 
      
      $WHAM$Radius
      HA FA 
      NA NA 
      
      $WHAM$MolWt
      HA FA 
      NA NA 
      
      $WHAM$MonodentTable
      [1] S          AbundDenom StrongWeak
      <0 rows> (or 0-length row.names)
      
      $WHAM$BidentTable
      [1] S1         S2         AbundDenom
      <0 rows> (or 0-length row.names)
      
      $WHAM$TridentTable
      [1] S1         S2         S3         AbundDenom
      <0 rows> (or 0-length row.names)
      
      $WHAM$MetalsTable
      [1] Metal  pKMAHA pKMAFA dLK2  
      <0 rows> (or 0-length row.names)
      
      $WHAM$SpecKselTable
      [1] Spec   KselHA KselFA
      <0 rows> (or 0-length row.names)
      
      
      $Index
      $Index$AqueousMCR
      [1] 1
      
      $Index$BioticLigMCR
      [1] NA
      
      $Index$WHAMDonnanMCR
      HA FA 
      -1 -1 
      
      

