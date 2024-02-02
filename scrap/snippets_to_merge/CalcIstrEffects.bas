SUB CalcIStrEffects(ISTR#)
'  *-----------------------------------------------------------------------*
'  |  SUB CalcIStrEffects(ISTR#)                                           |
'  |                                                                       |
'  |  Calculates the net charge on humic and fulvic substances (units are  |
'  |  (equivalents / L) then adjusts the logK's accordingly.               |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

'  *------------------------------------------*
'  |  Declare local variables and Initialize  |
'  *------------------------------------------*
    DIM w_HA AS SINGLE
    DIM w_FA AS SINGLE
    DIM Min_Z_HA AS DOUBLE
    DIM Min_Z_FA AS DOUBLE
    DIM Max_Z_HA AS DOUBLE
    DIM Max_Z_FA AS DOUBLE
    DIM adj_HA AS DOUBLE
    DIM adj_FA AS DOUBLE
    DIM iP_HA AS INTEGER
    DIM iP_FA AS INTEGER
    DIM TmpD AS DOUBLE

'  *-------------------------------------------*
'  |  Check if WHAM database file gave values  |
'  |  for HS characteristics.                  |
'  *-------------------------------------------*
    IF (wMOLWT(1) = 0) THEN  'MW was not read from WHAM file, so know other values are missing too
        wMOLWT(1) = 15000
        wMOLWT(2) = 1500
        wP(1) = -196'-374
        wP(2) = -119'-103
        wRADIUS(1) = 1.72E-9
        wRADIUS(2) = 8E-10
        DLF = 0.25
        wKZED = 1000
        iR_HA = FindSpecies("R_HA")
        iR_FA = FindSpecies("R_FA")
    END IF

'  *------------------------------------------------*
'  |  Detemine the Net Charge on HA and FA species  |
'  *------------------------------------------------*
    iP_HA = FindSpecies("P_HA")
    iP_FA = FindSpecies("P_FA")
    w_HA! = wP(1) * LOG10(ISTR#)
    w_FA! = wP(2) * LOG10(ISTR#)

    Z_HA# = 0
    Z_FA# = 0
    FOR S = 1 TO NSpecies
       TmpD = OC(S) * V(S)
       IF (ActCorr(S) = acWHAMHA) THEN
          Z_HA# = Z_HA# + TmpD
       END IF
       IF (ActCorr(S) = acWHAMFA) THEN
          Z_FA# = Z_FA# + TmpD
       END IF
    NEXT S


    IF (iP_HA > 0) AND (Z_HA < 0) THEN
        KT(iR_HA) = -Z_HA / 100
        KT(iP_HA) = Z_HA / 100
        OC(iP_HA) = 1 / EXP(2 * w_HA * Z_HA)
    ELSE
        break
    END IF
	
    IF (iP_FA > 0) AND (Z_FA < 0) AND (Z_FA > -0.1) THEN
        KT(iR_FA) = -Z_FA / 100
        KT(iP_FA) = Z_FA / 100
    ELSE
        break
    END IF


END SUB                                                  