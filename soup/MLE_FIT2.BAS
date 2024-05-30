DECLARE FUNCTION ErfC(BYVAL X AS DOUBLE) AS DOUBLE
DECLARE FUNCTION FUNC (X AS EXTENDED, PCOM() AS EXTENDED, XICOM() AS EXTENDED, x_dat() AS EXTENDED, x_cf() AS LONG, NData AS LONG) AS EXTENDED
DECLARE FUNCTION FUNC2 (param() AS EXTENDED, x_dat() AS EXTENDED, x_cf() AS LONG, NData AS LONG) AS EXTENDED
DECLARE FUNCTION NormalCDF(X AS EXTENDED, x_mean AS EXTENDED, x_std AS EXTENDED)  AS EXTENDED
DECLARE FUNCTION NormalPDF(X AS EXTENDED, x_mean AS EXTENDED, x_std AS EXTENDED)  AS EXTENDED
DECLARE FUNCTION PlotPos(BYVAL i AS LONG, n AS LONG) AS DOUBLE
DECLARE FUNCTION SNormInv(p AS DOUBLE) AS EXTENDED
' --------------------------------------------------------------------------------------------
'
'   Routines for fitting an MLE distribution to a series of data that includes censored values
'   For instructions on use, see the subrouting MLE_Fit
'
'   Robert Santore
'   HydroQual, Inc
'   rsantore@hydroqual.com
'
' --------------------------------------------------------------------------------------------

' search for "TODO" or "EDIT" to find areas of code that need attention

' TODO : eliminate the need for PCOM() and XICOM() arrays to be declared AS GLOBAL
'           assigned in
'              LINMIN
'           used in
'              F1DIM, which is also used in FUNC
'


'#COMPILE EXE
'#DIM ALL

FUNCTION BRENT (AX AS EXTENDED, BX AS EXTENDED, CX AS EXTENDED, DUM AS EXTENDED, _
                TOL AS EXTENDED, XMIN AS EXTENDED, PCOM() AS EXTENDED, XICOM() AS EXTENDED, _
                x_dat() AS EXTENDED, x_cf() AS LONG, NData AS LONG) AS EXTENDED
    DIM ITMAX AS LONG
    DIM CGOLD AS DOUBLE
    DIM ZEPS AS DOUBLE
    DIM A AS EXTENDED
    DIM B AS EXTENDED
    DIM D AS EXTENDED
    DIM V AS EXTENDED
    DIM W AS EXTENDED
    DIM X AS EXTENDED
    DIM E AS EXTENDED
    DIM FX AS EXTENDED
    DIM FV AS EXTENDED
    DIM FW AS EXTENDED
    DIM ITER AS LONG
    DIM XM AS EXTENDED
    DIM TOL1 AS EXTENDED
    DIM TOL2 AS EXTENDED
    DIM ETEMP AS EXTENDED
    DIM DONE AS INTEGER
    DIM P AS EXTENDED
    DIM Q AS EXTENDED
    DIM R AS EXTENDED
    DIM U AS EXTENDED
    DIM FU AS EXTENDED


    ITMAX = 100
    CGOLD = .381966#
    ZEPS = 1E-10
    A = AX
    IF CX < AX THEN A = CX
    B = AX
    IF CX > AX THEN B = CX
    V = BX
    W = V
    X = V
    E = 0!
    FX = FUNC(X, PCOM(), XICOM(),x_dat(), x_cf(), NData)
    FV = FX
    FW = FX
    FOR ITER = 1 TO ITMAX
      XM = .5 * (A + B)
      TOL1 = TOL * ABS(X) + ZEPS
      TOL2 = 2! * TOL1
      IF ABS(X - XM) <= TOL2 - .5 * (B - A) THEN EXIT FOR
      DONE% = -1
      IF ABS(E) > TOL1 THEN
        R = (X - W) * (FX - FV)
        Q = (X - V) * (FX - FW)
        P = (X - V) * Q - (X - W) * R
        Q = 2! * (Q - R)
        IF Q > 0! THEN P = -P
        Q = ABS(Q)
        ETEMP = E
        E = D
        DUM = ABS(.5 * Q * ETEMP)
        IF ABS(P) < DUM AND P > Q * (A - X) AND P < Q * (B - X) THEN
          D = P / Q
          U = X + D
          IF U - A < TOL2 OR B - U < TOL2 THEN D = ABS(TOL1) * SGN(XM - X)
          DONE% = 0
        END IF
      END IF
      IF DONE% THEN
        IF X >= XM THEN
          E = A - X
        ELSE
          E = B - X
        END IF
        D = CGOLD * E
      END IF
      IF ABS(D) >= TOL1 THEN
        U = X + D
      ELSE
        U = X + ABS(TOL1) * SGN(D)
      END IF
      FU = FUNC(U, PCOM(), XICOM(), x_dat(), x_cf(), NData) '(X AS EXTENDED, PCOM() AS EXTENDED, XICOM() AS EXTENDED, x_dat() AS EXTENDED, x_cf() AS LONG, N AS LONG)
      IF FU <= FX THEN
        IF U >= X THEN
          A = X
        ELSE
          B = X
        END IF
        V = W
        FV = FW
        W = X
        FW = FX
        X = U
        FX = FU
      ELSE
        IF U < X THEN
          A = U
        ELSE
          B = U
        END IF
        IF FU <= FW OR W = X THEN
          V = W
          FV = FW
          W = U
          FW = FU
        ELSEIF FU <= FV OR V = X OR V = W THEN
          V = U
          FV = FU
        END IF
      END IF
    NEXT ITER

    ' EDIT : get rid of this print statement so that we can compile as dll
    IF ITER > ITMAX THEN PRINT "Brent exceed maximum iterations.": END
    XMIN = X
    BRENT = FX
END FUNCTION   ' BRENT

FUNCTION Erf_DONT_USE(BYVAL X AS DOUBLE) AS DOUBLE
'Error function
'
'The integral is
'
'                          x
'                           -
'                2         | |          2
'  erf(x)  =  --------     |    exp( - t  ) dt.
'             sqrt(pi)   | |
'                         -
'                          0
'
'For 0 <= |x| < 1, erf(x) = x * P4(x**2)/Q5(x**2); otherwise
'erf(x) = 1 - erfc(x).
'
'
'ACCURACY:
'
'                     Relative error:
'arithmetic   domain     # trials      peak         rms
'   IEEE      0,1         30000       3.7e-16     1.0e-16
'
'Cephes Math Library Release 2.8:  June, 2000
'Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
'

    DIM Result AS DOUBLE
    DIM XSq AS DOUBLE
    DIM S AS DOUBLE
    DIM P AS DOUBLE
    DIM Q AS DOUBLE

    S = SGN(X)
    X = ABS(X)
    IF X < 0.5 THEN
        XSq = X * X
        P = 7.54772803341863E-03
        P = 0.288805137207594 + XSq * P
        P = 14.3383842191748 + XSq * P
        P = 38.0140318123903 + XSq * P
        P = 3017.82788536508 + XSq * P
        P = 7404.07142710151 + XSq * P
        P = 80437.363096084 + XSq * P
        Q = 0#
        Q = 1# + XSq * Q
        Q = 38.0190713951939 + XSq * Q
        Q = 658.07015545924 + XSq * Q
        Q = 6379.60017324428 + XSq * Q
        Q = 34216.5257924629 + XSq * Q
        Q = 80437.363096084 + XSq * Q
        Result = S * 1.12837916709551 * X * P / Q
        Erf_DONT_USE = Result
        EXIT FUNCTION
    END IF
    IF X >= 10# THEN
        Result = S
        Erf_DONT_USE = Result
        EXIT FUNCTION
    END IF
    Result = S * (1# - ErfC(X))

    Erf_DONT_USE = Result
END FUNCTION

'
'  This seems to be a better ERF function
'
FUNCTION ERF(x AS DOUBLE) AS DOUBLE
' http://www.mrexcel.com/forum/showthread.php?t=234208
'
' This seems to be a decent ERF function
' tested against the internal function in Excel
' R. Santore Sep 2011
DIM f AS DOUBLE, c AS DOUBLE, pi AS DOUBLE
DIM j AS INTEGER
DIM S AS INTEGER

' note, following two lines were added (along with corresponding
' Dim S statement and addition of S to function result on last line)
' to allow negative input values
' R. Santore
'
S = SGN(x)
x = ABS(x)
c = 0
pi = 3.14159265358979
IF 1.5 < x THEN
    c = 2 - c
    j = 3 + INT(32 / x)
    f = 0
    DO WHILE j <> 0
        f = 1 / (f * j + x * SQR(2))
        j = j - 1
    LOOP
    f = f * c * (3 - c * c) * EXP(-x * x) / SQR(2 * pi) + (c - 1) * (3 * c * c + c - 8) / 6
ELSE
    j = 3 + INT(9 * x)
    f = 1
    DO WHILE j <> 0
        f = 1 + f * x * x * (0.5 - j) / j / (0.5 + j)
        j = j - 1
    LOOP
    f = c + f * x * (2 - 4 * c) / SQR(pi)
END IF
    ERF = f * S
END FUNCTION


'Complementary error function
'
' 1 - erf(x) =
'
'                          inf.
'                            -
'                 2         | |          2
'  erfc(x)  =  --------     |    exp( - t  ) dt
'              sqrt(pi)   | |
'                          -
'                           x
'
'
'For small x, erfc(x) = 1 - erf(x); otherwise rational
'approximations are computed.
'
'
'ACCURACY:
'
'                     Relative error:
'arithmetic   domain     # trials      peak         rms
'   IEEE      0,26.6417   30000       5.7e-14     1.5e-14
'
'Cephes Math Library Release 2.8:  June, 2000
'Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
'
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
FUNCTION ErfC(BYVAL X AS DOUBLE) AS DOUBLE
    DIM Result AS DOUBLE
    DIM P AS DOUBLE
    DIM Q AS DOUBLE

    IF X < 0# THEN
        Result = 2# - ErfC(-X)
        ErfC = Result
        EXIT FUNCTION
    END IF
    IF X < 0.5 THEN
        Result = 1# - Erf(X)
        ErfC = Result
        EXIT FUNCTION
    END IF
    IF X >= 10# THEN
        Result = 0#
        ErfC = Result
        EXIT FUNCTION
    END IF
    P = 0#
    P = 0.56418778255074 + X * P
    P = 9.67580788298727 + X * P
    P = 77.0816173036843 + X * P
    P = 368.519615471001 + X * P
    P = 1143.26207070389 + X * P
    P = 2320.43959025164 + X * P
    P = 2898.02932921677 + X * P
    P = 1826.33488422951 + X * P
    Q = 1#
    Q = 17.1498094362761 + X * Q
    Q = 137.125596050062 + X * Q
    Q = 661.736120710765 + X * Q
    Q = 2094.38436778954 + X * Q
    Q = 4429.61280388368 + X * Q
    Q = 6089.54242327244 + X * Q
    Q = 4958.82756472114 + X * Q
    Q = 1826.33488422951 + X * Q
    'Result = Exp(-Square(X)) * P / Q
    Result = EXP((X ^ 2)) * P / Q
    Result = EXP(-SQR(X)) * P / Q
    ErfC = Result
END FUNCTION

SUB EstimateProbabilities(x_dat() AS EXTENDED, x_mean AS EXTENDED, x_std AS EXTENDED, x_cf()AS LONG, x_prob()AS EXTENDED, prob_j AS EXTENDED, N AS LONG)
    DIM i AS LONG
    DIM prob_i AS EXTENDED
    prob_j = 1##
    FOR i=1 TO N
        IF (x_cf(i)=0) THEN
            '
            ' For this data point the value is an actual value
            prob_i = NormalPDF(x_dat(i), x_mean, x_std)
        ELSEIF (x_cf(i)=1) THEN
            '
            ' This is a less than value
            prob_i = NormalCDF(x_dat(i), x_mean, x_std)
        ELSEIF (x_cf(i)=2) THEN
            '
            ' This is a greater than value
            prob_i = 1 - NormalCDF(x_dat(i), x_mean, x_std)
        END IF
        x_prob(i) = prob_i
        prob_j = prob_j * prob_i
    NEXT i

END SUB

FUNCTION F1DIM (X AS EXTENDED, PCOM() AS EXTENDED, XICOM() AS EXTENDED, x_dat() AS EXTENDED, x_cf() AS LONG, N AS LONG) AS EXTENDED
    DIM NCOM AS LONG

    ' TODO Edit : why isn't this a passable parameter?
    NCOM = 2
    DIM XT(NCOM) AS EXTENDED
    DIM J AS LONG

    FOR J = 1 TO NCOM
      XT(J) = PCOM(J) + X * XICOM(J)
    NEXT J
    F1DIM = FUNC2(XT(), x_dat(), x_cf(), N)
    ERASE XT
END FUNCTION

FUNCTION FUNC (X AS EXTENDED, PCOM() AS EXTENDED, XICOM() AS EXTENDED, x_dat() AS EXTENDED, x_cf() AS LONG, N AS LONG) AS EXTENDED
   FUNC = F1DIM(X, PCOM(), XICOM(), x_dat(), x_cf(), N)
END FUNCTION

FUNCTION FUNC2 (param() AS EXTENDED, x_dat() AS EXTENDED, x_cf() AS LONG, NData AS LONG) AS EXTENDED
    'FUNC2 = .5 - BESSJ0((X(1) - 1!) ^ 2 + (X(2) - 2!) ^ 2 + (X(3) - 3!) ^ 2)
    'CALL EstimateProbabilities(x(), x_mean, x_std, x_cf(), x_prob(), prob_j, N)

    ' global vars
    '
    '   x_dat(i)
    '   x_cf(i)
    '   x_prob(i)

    DIM i AS LONG
    DIM prob_i AS EXTENDED
    DIM prob_j AS EXTENDED
    DIM x_mean AS EXTENDED
    DIM x_std AS EXTENDED

    x_mean = param(1)
    x_std = param(2)

    prob_j = 1##
    FOR i=1 TO NData
        IF (x_cf(i)=0) THEN
            '
            ' For this data point the value is an actual value
            prob_i = NormalPDF(x_dat(i), x_mean, x_std)
        ELSEIF (x_cf(i)=1) THEN
            '
            ' This is a less than value
            prob_i = NormalCDF(x_dat(i), x_mean, x_std)
        ELSEIF (x_cf(i)=2) THEN
            '
            ' This is a greater than value
            prob_i = 1 - NormalCDF(x_dat(i), x_mean, x_std)
        END IF
        'x_prob(i) = prob_i
        prob_j = prob_j * prob_i
    NEXT i

    FUNC2 = - prob_j

END FUNCTION

SUB IndSortExt (eArray() AS EXTENDED, Index() AS LONG, NEls AS LONG, iDir AS LONG)
    '  *-----------------------------------------------------------------------*
    '  |                                                                       |
    '  | SUB IndSortI                                                          |
    '  |                                                                       |
    '  | Used for performing an indexed sort on integer arrays.                |
    '  |                                                                       |
    '  | Input:                                                                |
    '  |    NEls      - integer -   Number of elements to be sorted            |
    '  |                                                                       |
    '  |    iDir      - integer -   Direction of sort.                         |
    '  |                            = 0, sort in ascending order               |
    '  |                            <> 0 , sort in descending order            |
    '  |                                                                       |
    '  |    eArray()  - integer -   Array dimensioned from at least 1 to NEls  |
    '  |                            contains the values to be sorted           |
    '  |                                                                       |
    '  | Output:                                                               |
    '  |    Index()   - integer -   Array that holds the sorted order of       |
    '  |                            the values in iArray.  Note that this array|
    '  |                            does not have to be initialized as it does |
    '  |                            in the assembly version.                   |
    '  |                                                                       |
    '  |  Developed for CHESS ver 3.00, May 1995                               |
    '  |                                                                       |
    '  |  Robert Santore                                                       |
    '  |  Dept. Civil and Environmental Eng, Syracuse University               |
    '  |  rsantore@mailbox.syr.edu                                             |
    '  |                                                                       |
    '  *-----------------------------------------------------------------------*
    DIM i AS INTEGER
    DIM j AS INTEGER
    DIM k AS INTEGER
    DIM InOrder AS INTEGER
    DIM TRUE AS LONG
    DIM FALSE AS LONG
    FALSE = 0
    TRUE = NOT FALSE

    ' *------------------------------*
    ' |  Initialize the Index array  |
    ' *------------------------------*
    FOR i = 1 TO NEls
       Index(i) = i
    NEXT i

    ' *-------------------*
    ' | Perform the sort  |
    ' *-------------------*
    DO
       InOrder = TRUE
       FOR i = 1 TO NEls - 1
          IF iDir THEN
             ' *----------------------------*
             ' |  Sort in descending order  |
             ' *----------------------------*
             j = i + 1
             k = i
          ELSE
             ' *---------------------------*
             ' |  Sort in ascending order  |
             ' *---------------------------*
             j = i
             k = i + 1
          END IF

          IF eArray(Index(j)) > eArray(Index(k)) THEN
             SWAP Index(j), Index(k)
             InOrder = FALSE
          END IF
       NEXT i
    LOOP UNTIL InOrder

END SUB

SUB LINMIN (P() AS EXTENDED, XI() AS EXTENDED, NDim AS LONG, NData AS LONG, FRET AS EXTENDED, x_dat() AS EXTENDED, x_cf() AS LONG)
'    CALL MNBRAK(AX, XX, BX, FA, FX, FB, DUM, PCOM(), XICOM(), x_dat(), x_cf(), NData)
    ' Given an N dimentional point P and an N dimensional direction
    ' XI, moves and resets P to where the function FUNC(P) takes on
    ' a minimum along the direction XI from P, and replaces XI by
    ' the actual vector displacement that P was moved.  Also
    ' returns as FRET the value of FUNC at the returned location P.
    ' This is actually all accomplished by calling the routines
    ' MNBRAK and BRENT.

    DIM TOL AS EXTENDED
    DIM NCOM AS LONG
    DIM J AS LONG
    DIM AX AS EXTENDED
    DIM XX AS EXTENDED
    DIM BX AS EXTENDED
    DIM FA AS EXTENDED
    DIM FX AS EXTENDED
    DIM FB AS EXTENDED
    DIM DUM AS EXTENDED
    DIM XMIN AS EXTENDED
    DIM PCOM(1 TO NDim) AS EXTENDED ' AS GLOBAL EXTENDED
    DIM XICOM(1 TO NDim) AS EXTENDED ' AS GLOBAL EXTENDED

    TOL = .0001
    NCOM = NDim
    FOR J = 1 TO NDim
      PCOM(J) = P(J)
      XICOM(J) = XI(J)
    NEXT J
    AX = 0!
    XX = 1!
    CALL MNBRAK(AX, XX, BX, FA, FX, FB, DUM, PCOM(), XICOM(), x_dat(), x_cf(), NData)
'SUB MNBRAK (AX AS EXTENDED, BX AS EXTENDED, CX AS EXTENDED, FA AS EXTENDED, FB AS EXTENDED, FC AS EXTENDED, DUM AS EXTENDED, _
'            PCOM() AS EXTENDED, XICOM() AS EXTENDED, x_dat() AS EXTENDED, x_cf() AS LONG, N AS LONG)

    FRET = BRENT(AX, XX, BX, DUM, TOL, XMIN, PCOM(), XICOM(), x_dat(), x_cf(), NData)
    FOR J = 1 TO NDim
      XI(J) = XMIN * XI(J)
      P(J) = P(J) + XI(J)
    NEXT J
END SUB

SUB Median1 (X() AS EXTENDED, X_ind() AS LONG, N AS LONG, XMED AS EXTENDED)
    '
    ' Returns the median of the values in an array
    '
    '    X()    : array of N elements for which a median value will be calculated
    '    X_ind(): array of N elements which contains the sorted order of array X
    '    N      : number of elements
    '    XMed   : will return the median on output

    'CALL SORT(N, X())
    DIM N2 AS LONG
    N2 = INT(N / 2)
    IF 2 * N2 = N THEN
      XMED = .5 * (X(X_ind(N2)) + X(X_ind(N2 + 1)))
    ELSE
      XMED = X(X_ind(N2 + 1))
    END IF
END SUB

SUB MLE_Fit(x_dat_in() AS EXTENDED, x_dl() AS LONG, NData AS LONG, iLog AS LONG, x_mle_mean AS EXTENDED, x_mle_std AS EXTENDED)
'        SUB MLE_Fit
'
'        This routine uses maximum likelihood estimation to determine the mean and standard deviation of
'        data that includes censored values.  Censored values may include observations where the true value
'        is either less than or greater than the reported value.
'
'        Input values:
'          x_dat : floating point vector containing values to be analyzed
'          x_dl  : integer vector with detection limit flags (0 = no flag, 1 = actual value is less than, 2 = actual value
'                  is greater than)
'          NData : integer variable with number of values
'          iLog  : boolean variable that is TRUE if data should be log transformed, or FALSE if not
'

'        Output values:
'        x_mle_mean  : a floating point variable with the mean of the MLE distribution
'        x_mle_std   : a floating point variable with the std of the MLE distribution
'
'        Robert Santore
'        HydroQual Inc
'        rsantore@hydroqual.com
'
    ' The number of parameters to optimize = 2
    DIM NDim AS LONG
    NDIM = 2
    DIM FTOL AS EXTENDED
    FTOL = .00000001
    DIM P(NDim) AS EXTENDED
    DIM XI(NDim, NDim)  AS EXTENDED
    DIM NP AS LONG
    DIM J AS LONG
    DIM ITER AS LONG
    DIM FRET AS EXTENDED
    DIM x_dat(1 TO NData) AS EXTENDED
    DIM x_AVE AS EXTENDED, x_ADEV AS EXTENDED, x_SDEV AS EXTENDED, x_VAR AS EXTENDED, x_SKEW AS EXTENDED, x_CURT AS EXTENDED

    IF (iLog) THEN
       FOR j=1 TO NData
          x_dat(j) = LOG10(x_dat_in(j))
       NEXT j
    ELSE
       FOR j=1 TO NData
          x_dat(j) = x_dat_in(j)
       NEXT j
    END IF
    '
    ' Set initial starting values for the parameters we will optimize
    '
    CALL Moment (x_dat(), NData, x_AVE, x_ADEV, x_SDEV, x_VAR, x_SKEW, x_CURT)
    x_mle_mean = x_AVE
    x_mle_std = x_SDEV
    P(1) = x_mle_mean
    P(2) = x_mle_std
    '
    ' Set the initial search directions
    '
    XI(1,1) = 1
    XI(1,2) = 0
    XI(2,1) = 0
    XI(2,2) = 1
'
'    Perform the optimization
     CALL POWELL(P(), XI(), NDIM, x_dat(), x_dl(), NData, FTOL, ITER, FRET)
'     PRINT "2nd"
'     CALL POWELL(P(), XI(), NDIM, x_dat(), x_dl(), NData, FTOL, ITER, FRET)
         'PRINT

    ' set the optimized parameter values back into the variables to return to the calling program
    x_mle_mean = P(1)
    x_mle_std  = P(2)

END SUB

SUB MNBRAK (AX AS EXTENDED, BX AS EXTENDED, CX AS EXTENDED, FA AS EXTENDED, FB AS EXTENDED, FC AS EXTENDED, DUM AS EXTENDED, _
            PCOM() AS EXTENDED, XICOM() AS EXTENDED, x_dat() AS EXTENDED, x_cf() AS LONG, N AS LONG)
    DIM GOLD  AS EXTENDED
    DIM GLIMIT AS EXTENDED
    DIM TINY  AS EXTENDED
    DIM DONE AS INTEGER
    DIM R  AS EXTENDED
    DIM Q  AS EXTENDED
    DIM U  AS EXTENDED
    DIM ULIM  AS EXTENDED
    DIM FU   AS EXTENDED


    GOLD = 1.618034
    GOLD = 1.618034 / 40

    GLIMIT = 100!
    TINY = 1E-20
    FA = FUNC(AX, PCOM(), XICOM(), x_dat(), x_cf(), N)
    FB = FUNC(BX, PCOM(), XICOM(), x_dat(), x_cf(), N)
    IF FB > FA THEN
      DUM = AX
      AX = BX
      BX = DUM
      DUM = FB
      FB = FA
      FA = DUM
    END IF
    CX = BX + GOLD * (BX - AX)
    FC = FUNC(CX, PCOM(), XICOM(), x_dat(), x_cf(), N)
    DO
      IF FB < FC THEN EXIT DO
      DONE% = -1
      R = (BX - AX) * (FB - FC)
      Q = (BX - CX) * (FB - FA)
      DUM = Q - R
      IF ABS(DUM) < TINY THEN DUM = TINY
      U = BX - ((BX - CX) * Q - (BX - AX) * R) / (2! * DUM)
      ULIM = BX + GLIMIT * (CX - BX)
      IF (BX - U) * (U - CX) > 0! THEN
        FU = FUNC(U, PCOM(), XICOM(), x_dat(), x_cf(), N)
        IF FU < FC THEN
          AX = BX
          FA = FB
          BX = U
          FB = FU
          EXIT SUB
        ELSEIF FU > FB THEN
          CX = U
          FC = FU
          EXIT SUB
        END IF
        U = CX + GOLD * (CX - BX)
        FU = FUNC(U, PCOM(), XICOM(), x_dat(), x_cf(), N)
      ELSEIF (CX - U) * (U - ULIM) > 0! THEN
        FU = FUNC(U, PCOM(), XICOM(), x_dat(), x_cf(), N)
        IF FU < FC THEN
          BX = CX
          CX = U
          U = CX + GOLD * (CX - BX)
          FB = FC
          FC = FU
          FU = FUNC(U, PCOM(), XICOM(), x_dat(), x_cf(), N)
        END IF
      ELSEIF (U - ULIM) * (ULIM - CX) >= 0! THEN
        U = ULIM
        FU = FUNC(U, PCOM(), XICOM(), x_dat(), x_cf(), N)
      ELSE
        U = CX + GOLD * (CX - BX)
        FU = FUNC(U, PCOM(), XICOM(), x_dat(), x_cf(), N)
      END IF
      IF DONE% THEN
        AX = BX
        BX = CX
        CX = U
        FA = FB
        FB = FC
        FC = FU
      ELSE
        DONE% = 0
      END IF
    LOOP WHILE NOT DONE%
END SUB

SUB MOMENT (DATQ() AS EXTENDED, N AS LONG, AVE AS EXTENDED, ADEV AS EXTENDED, SDEV AS EXTENDED, VAR AS EXTENDED, SKEW AS EXTENDED, CURT AS EXTENDED)
    DIM P AS EXTENDED
    DIM S  AS EXTENDED
    DIM J AS LONG

    ' EDIT : get rid of this print statement so that we can compile as dll
    IF N <= 1 THEN PRINT "N must be at least 2": EXIT SUB
    S = 0!
    FOR J = 1 TO N
      S = S + DATQ(J)
    NEXT J
    AVE = S / N
    ADEV = 0##
    VAR = 0##
    SKEW = 0##
    CURT = 0##
    FOR J = 1 TO N
      S = DATQ(J) - AVE
      ADEV = ADEV + ABS(S)
      P = S * S
      VAR = VAR + P
      P = P * S
      SKEW = SKEW + P
      P = P * S
      CURT = CURT + P
    NEXT J
    ADEV = ADEV / N
    VAR = VAR / (N - 1)
    SDEV = SQR(VAR)
    IF VAR <> 0! THEN
      SKEW = SKEW / (N * SDEV ^ 3)
      CURT = CURT / (N * VAR ^ 2) - 3!
    ELSE
      ' EDIT : get rid of this print statement so that we can compile as dll
       PRINT "no skew or kurtosis when zero variance"
    END IF
END SUB

FUNCTION NormalCDF(X AS EXTENDED, x_mean AS EXTENDED, x_std AS EXTENDED)  AS EXTENDED
   NormalCDF = 0.5 * (Erf((X - x_mean) / (x_std * 1.4142135623731)) + 1#)
END FUNCTION


FUNCTION NormalPDF(X AS EXTENDED, x_mean AS EXTENDED, x_std AS EXTENDED)  AS EXTENDED

'subroutine normal_pdf ( x, a, b, pdf )
'!
'!*******************************************************************************
'!
'!! NORMAL_PDF evaluates the Normal PDF.
'!
'!
'!  Formula:
'!
'!    PDF(X)(A,B)
'!      = EXP ( - 0.5E+00 * ( ( X - A ) / B )**2 )
'!      / ( B * SQRT ( 2 * PI ) )
'!
'!  Discussion:
'!
'!    The normal PDF is also known as the Gaussian PDF.
'!
'!
'!  Parameters:
'!
'!    Input, real X, the argument of the PDF.
'!
'!    Input, real A, B, the parameters of the PDF.
'!    0.0E+00 < B.
'!
'!    Output, real PDF, the value of the PDF.
'!


'Dim a, b, pdf, x, y As Double
DIM y  AS EXTENDED
DIM pdf  AS EXTENDED
DIM Pi AS EXTENDED
    Pi = 3.14159265

  y = (X - x_mean) / x_std
  pdf = EXP(-0.5 * y ^ 2) / (x_std * SQR(2# * Pi))

   NormalPDF = pdf

END FUNCTION

SUB POWELL (P() AS EXTENDED, XI() AS EXTENDED, NDim AS LONG, x_dat() AS EXTENDED, x_cf() AS LONG, NData AS LONG, FTOL AS EXTENDED, ITER AS LONG, FRET AS EXTENDED)
    DIM ITMAX AS LONG

    ITMAX = 200
    DIM PT(NDim)  AS EXTENDED
    DIM PTT(NDim) AS EXTENDED
    DIM XIT(NDim) AS EXTENDED
    DIM J AS LONG
    DIM FP AS EXTENDED
    DIM IBIG AS LONG
    DIM DEL AS EXTENDED
    DIM I AS LONG
    DIM FPTT AS EXTENDED
    DIM DUM AS EXTENDED
    DIM T AS EXTENDED
    'dim FRET as extended

 '         FUNC2 (param() AS EXTENDED, x_dat() AS EXTENDED, x_cf() AS LONG, NData AS LONG) AS EXTENDED
    FRET = FUNC2(P(), x_dat(), x_cf(), NData)
    FOR J = 1 TO NDim
      PT(J) = P(J)
    NEXT J
    ITER = 0
    DO
      DO
        DO
          ITER = ITER + 1
          FP = FRET
          IBIG = 0
          DEL = 0!
          FOR I = 1 TO NDim
            FOR J = 1 TO NDim
              XIT(J) = XI(J, I)
            NEXT J
            FPTT = FRET
            CALL LINMIN(P(), XIT(), NDim, NData, FRET, x_dat(), x_cf())

            ' EDIT : get rid of this print statement so that we can compile as dll
            PRINT USING$(" ##.######^^^^", P(1), P(2), -FRET)

            IF ABS(FPTT - FRET) > DEL THEN
              DEL = ABS(FPTT - FRET)
              IBIG = I
            END IF
          NEXT I
          IF 2! * ABS(FP - FRET) <= FTOL * (ABS(FP) + ABS(FRET)) THEN
            ERASE XIT, PTT, PT
            EXIT SUB
          END IF
          IF ITER = ITMAX THEN
            ' EDIT : get rid of this print statement so that we can compile as dll
            PRINT "Powell exceeding maximum iterations."
            EXIT SUB
          END IF
          FOR J = 1 TO NDim
            PTT(J) = 2! * P(J) - PT(J)
            XIT(J) = P(J) - PT(J)
            PT(J) = P(J)
          NEXT J
          FPTT = FUNC2(PTT(), x_dat(), x_cf(), NData)
        LOOP WHILE FPTT >= FP
        DUM = FP - 2! * FRET + FPTT
        T = 2! * DUM * (FP - FRET - DEL) ^ 2 - DEL * (FP - FPTT) ^ 2
      LOOP WHILE T >= 0!
      CALL LINMIN(P(), XIT(), NDim, NData, FRET, x_dat(), x_cf())
      FOR J = 1 TO NDim
        XI(J, IBIG) = XIT(J)
      NEXT J
    LOOP
END SUB

FUNCTION ZScore(X AS LONG, n AS LONG) AS EXTENDED
    '
    ' Calculate the ZScore for position x of n total items
    '
    ' This ZScore uses the same plotting position used by Helsel
    ' For the most part we should only use this when calculating ZScore
    ' for MLE optimization when consistency with Helsel cenmle is desired
    '
    ZScore = SNormInv(PlotPos(x, n))
END FUNCTION

FUNCTION ZScoreNorm(X AS LONG, n AS LONG) AS EXTENDED
    '
    ' Calculate the ZScore for position x of n total items
    '
    ZScoreNorm = SNormInv(x/n)
END FUNCTION

FUNCTION ZScore_gdp(X AS LONG, n AS LONG) AS EXTENDED
    ' Function ZScore returns a single precition zscore for
    ' position x, given number of points n
    DIM C0 AS EXTENDED
    DIM C1 AS EXTENDED
    DIM C2 AS EXTENDED
    DIM D1 AS EXTENDED
    DIM D2 AS EXTENDED
    DIM D3 AS EXTENDED

    DIM n2 AS LONG
    DIM YSign AS LONG
    DIM y  AS EXTENDED
    DIM t AS EXTENDED

    C0 = 2.515517
    C1 = 0.802853
    C2 = 0.010328
    D1 = 1.432788
    D2 = 0.189269
    D3 = 0.001308
    n2& = (n& + 1) / 2

    IF (X& > n2&) THEN
        'x = FLOAT(n - (J - 1)) / (FN + 1#)
        y = (n& - (X& - 1)) / (n& + 1#)
        YSign& = 1#
    ELSE
        'x = FLOAT(J) / (FN + 1#)
        y = (X&) / (n& + 1#)
        YSign& = -1#
    END IF
    t = SQR(LOG(1# / (y * y)))
    ZScore_gdp = YSign& * (t - (C0 + t * (C1 + C2 * t)) / (1# + t * (D1 + t * (D2 + D3 * t))))

END FUNCTION

   'FUNCTION SNormInv(p AS DOUBLE) AS DOUBLE
   FUNCTION SNormInv(p AS DOUBLE) AS EXTENDED

      DIM q AS DOUBLE
      DIM r AS DOUBLE

      'Coefficients in rational approximations.

      DIM A1 AS DOUBLE
      DIM A2 AS DOUBLE
      DIM A3 AS DOUBLE
      DIM A4 AS DOUBLE
      DIM A5 AS DOUBLE
      DIM A6 AS DOUBLE

      DIM B1 AS DOUBLE
      DIM B2 AS DOUBLE
      DIM B3 AS DOUBLE
      DIM B4 AS DOUBLE
      DIM B5 AS DOUBLE

      DIM C1 AS DOUBLE
      DIM C2 AS DOUBLE
      DIM C3 AS DOUBLE
      DIM C4 AS DOUBLE
      DIM C5 AS DOUBLE
      DIM C6 AS DOUBLE

      DIM D1 AS DOUBLE
      DIM D2 AS DOUBLE
      DIM D3 AS DOUBLE
      DIM D4 AS DOUBLE

      A1 = -39.6968302866538
      A2 = 220.946098424521
      A3 = -275.928510446969
      A4 = 138.357751867269
      A5 = -30.6647980661472
      A6 = 2.50662827745924

      B1 = -54.4760987982241
      B2 = 161.585836858041
      B3 = -155.698979859887
      B4 = 66.8013118877197
      B5 = -13.2806815528857

      C1 = -7.78489400243029E-03
      C2 = -0.322396458041136
      C3 = -2.40075827716184
      C4 = -2.54973253934373
      C5 = 4.37466414146497
      C6 = 2.93816398269878

      D1 = 7.78469570904146E-03
      D2 = 0.32246712907004
      D3 = 2.445134137143
      D4 = 3.75440866190742

      ' Define break-points.
      DIM p_LOW AS DOUBLE
      DIM P_HIGH AS DOUBLE
      P_LOW  = 0.02425
      P_HIGH = 1 - P_LOW

      IF (p > 0) AND (p < P_LOW) THEN

         'Rational approximation for lower region.
         q = SQR(-2 * LOG(p))

         SNormInv = (((((C1 * q + C2) * q + C3) * q + C4) * q + C5) * q + C6) / _
            ((((D1 * q + D2) * q + D3) * q + D4) * q + 1)
         EXIT FUNCTION

      ELSEIF (p >= P_LOW) AND (p <= P_HIGH) THEN

         'Rational approximation for central region.
         q = p - 0.5
         r = q * q

         SNormInv = (((((A1 * r + A2) * r + A3) * r + A4) * r + A5) * r + A6) * q / _
           (((((B1 * r + B2) * r + B3) * r + B4) * r + B5) * r + 1)
         EXIT FUNCTION

      ELSEIF (p > P_HIGH) AND (p < 1) THEN

         'Rational approximation for upper region.
         q = SQR(-2 * LOG(1 - p))

         SNormInv = -(((((C1 * q + C2) * q + C3) * q + C4) * q + C5) * q + C6) / _
                ((((D1 * q + D2) * q + D3) * q + D4) * q + 1)
         EXIT FUNCTION
      ELSE

         'Throw New ArgumentOutOfRangeException()

      END IF

   END FUNCTION

FUNCTION PlotPos(BYVAL i AS LONG, n AS LONG) AS DOUBLE
DIM a AS DOUBLE
IF (n <= 10) THEN
    a = 0.375
ELSE
    a = 0.5
END IF
PlotPos = (i - a) / (n + 1 - 2 * a)
END FUNCTION
