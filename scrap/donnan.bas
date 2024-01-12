SUB CalcDonnanLayer(ISTR#)
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  Calculates the effects from the Donnan layer/non-specific binding.   |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

'  *------------------------------------------*
'  |  Declare local variables and Initialize  |
'  *------------------------------------------*
    DIM N0 AS DOUBLE      'Avagadro's Number
    DIM IKappa AS DOUBLE  '1/kappa (Debye-Huckel characteristic distance)
    DIM VDP_HA AS DOUBLE  'Donnan Layer volume (L) per gram of HA
    DIM VDP_FA AS DOUBLE  '   "                     "          FA
    DIM iH_HADL AS INTEGER
    DIM iH_FADL AS INTEGER
    DIM iH_S AS INTEGER
    DIM Tmp AS DOUBLE
    N0 = 6.022E+23
    iH_HADL = FindSpecies("H_DL_HA")
    iH_FADL = FindSpecies("H_DL_FA")
    iH_S = FindSpecies("H")

'  *----------------------------------------------------------------*
'  |  Calculate the diffuse layer volume, VD_HA/FA, adjust CToM,    |
'  |  and set the equilibrium constant for DL                       |
'  *----------------------------------------------------------------*
    IKappa = 0.000000000304 / ISTR# ^ 0.5

    VDP_HA = 1000 * (N0 / wMOLWT(1)) * 4 * (3.14159 / 3) * ((wRADIUS(1) + IKappa) ^ 3 - wRADIUS(1) ^ 3)
    VDP_FA = 1000 * (N0 / wMOLWT(2)) * 4 * (3.14159 / 3) * ((wRADIUS(2) + IKappa) ^ 3 - wRADIUS(2) ^ 3)

    VDP_HA = VDP_HA * wKZED * ABS(Z_HA) / (1 + wKZED * ABS(Z_HA))
    VDP_FA = VDP_FA * wKZED * ABS(Z_FA) / (1 + wKZED * ABS(Z_FA))

    VDP_HA = VDP_HA * Sol_HA * wMOLWT(1)
    VDP_FA = VDP_FA * Sol_FA * wMOLWT(2)
    VD_HA_CToM = VDP_HA / (1 + VDP_HA / DLF + VDP_FA / DLF)
    VD_FA_CToM = VDP_FA / (1 + VDP_HA / DLF + VDP_FA / DLF)

    'Volumne should never actually be 0. We're setting a minimum value to avoid numerical issues.
    VD_FA_CToM = MAX(0#, VD_FA_CToM)
    VD_HA_CToM = MAX(0#, VD_HA_CToM)

    TempHA# = OC(iR_HA) 'OC(iH_HADL) / OC(iH_S)
    TempFA# = OC(iR_FA) 'OC(iH_FADL) / OC(iH_S)

    FOR S = 1 TO NSpecies
        SELECT CASE SType(S)
            CASE stAqueous
                CToM(S) = AqueousCToM - VD_FA_CToM - VD_HA_CToM
            CASE stDonnanHA
                CToM(S) = VD_HA_CToM
            CASE stDonnanFA
                CToM(S) = VD_FA_CToM
        END SELECT
    NEXT S

END SUB
