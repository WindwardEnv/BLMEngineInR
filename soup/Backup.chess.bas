

SUB CHESS (ErrorCode%, ErrorMsg$)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB CHESS                                                            |
'  |                                                                       |
'  |  Solves the equilibrium system defined by the parameter file opened   |
'  |  by sub DefineProblem.  Before any call to the CHESS subroutine is    |
'  |  made, a parameter file must be passed to SUB DefineProblem for       |
'  |  initialization.                                                      |
'  |                                                                       |
'  |  In addition to the information read into DefineProblem, it is up     |
'  |  to the calling program to store information in several variables     |
'  |  before a call to CHESS is made.                                      |
'  |                                                                       |
'  |     KT(1 to NComp)          Total Concentration (moles) of components |
'  |                                                                       |
'  |     OC(1 to NComp)          Free Conc. for ctFixed components, or     |
'  |                             first guess of the equil. conc. for       |
'  |                             ctFloat components (moles/liter)          |
'  |                                                                       |
'  |     PhaseMass(1 to NPhase)  Total mass of phases (moles).             |
'  |                                                                       |
'  |     AqueousCToM             Liters of water                           |
'  |                                                                       |
'  |     SurfaceCToM             Kilograms of soil                         |
'  |                                                                       |
'  |     SysTemp                 System temperature (Kelvin)               |
'  |                                                                       |
'  |                                                                       |
'  |  This version of CHESS uses Gaussian Elimination to solve the         |
'  |  linear approximation to the non-linear equations developed by CHESS  |
'  |                                                                       |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*
' *-------------------------------------------------------------*
' |  Type variables definitions and declarations                |
' *-------------------------------------------------------------*
' *---------------------------------------*
' |  Simple Variables Local to SUB CHESS  |
' *---------------------------------------*
DIM c AS INTEGER            '  Array/loop index variable
DIM CHARGE AS DOUBLE        '  Ionic Charge
DIM Converge AS INTEGER     '  Boolean flag to indicate convergence
DIM Iter AS INTEGER         '  Iteration count
DIM iSing AS INTEGER        '  Flag to indicate a singular matrix
DIM ISTR AS DOUBLE          '  Ionic Strength
DIM MaxError AS SINGLE      '  Value of the largest residual (relative units)
DIM OverFlag AS INTEGER     '  Flag to indicate an oversaturated phase
DIM Skip AS INTEGER         '  Flag to indicate if this is a full iteration
DIM UnderFlag AS INTEGER    '  Flag to indicate an undersaturated phase

' *--------------------------------------*
' |  Array Variables Local to SUB CHESS  |
' *--------------------------------------*
DIM X(1 TO NComp) AS DOUBLE'  Iterative improvement to component concen.

' *--------------------*
' |  Initialize Flags  |
' *--------------------*
Iter% = ZERO
UnderFlag% = FALSE
OverFlag% = FALSE
Skip% = FALSE

' *-----------------------------------------------*
' |  Set hEC() including temperature corrections  |
' *-----------------------------------------------*
CALL TempCorrection
IF (iDebug = 01 OR iDebug = 07) THEN CALL DebugLog ("CHES", "After TempCorrection")


' *-----------------------------------------------------------------*
' |  Set tSC(), tEC() including solid/gaseous phase transformation  |
' *-----------------------------------------------------------------*
CALL SSaturate
IF (iDebug = 01 OR iDebug = 07) THEN CALL DebugLog ("CHES", "After SSaturate")

' *-------------------------------------------------------*
' |  Generate first guesses for component concentrations  |
' |  if none were supplied by the calling program         |
' *-------------------------------------------------------*
CALL FirstGuess(ErrorCode%, ErrorMsg$, ISTR#)
IF (iDebug = 01 OR iDebug = 07) THEN CALL DebugLog ("CHES", "After FirstGuess")
IF ErrorCode% THEN
   EXIT SUB
END IF

' *---------------------------------------------------------*
' |  The species concentraions and residual values need to  |
' |  be calculated once to seed the iterative solution      |
' *---------------------------------------------------------*
CALL CalcSpeciesConc(ISTR#, ErrorCode%, ErrorMsg$)
IF (iDebug = 01 OR iDebug = 07) THEN CALL DebugLog ("CHES", "After CalcSpeciesConc (before loop)")

IF ErrorCode% THEN
   EXIT SUB
END IF
CALL CalcResidual(MaxError!, MaxComp%)
IF (iDebug = 01 OR iDebug = 07) THEN CALL DebugLog ("CHES", "After CalcResidual (before loop)")

Converge% = FALSE
DO
   Iter% = Iter% + 1
   iSing% = FALSE

   IF NOT Skip% THEN
      CALL JacobianAA
      CALL JacobianNum

      'call OutputJacobian

      IF (iDebug = 01 OR iDebug = 07) THEN CALL DebugLog ("CHES", "After JacobianAA")
      'CALL Gauss(Z(), R(), X(), NSolve, iSing)
      CALL Gauss(Zn(), R(), X(), NSolve, iSing)
      IF (iDebug = 01 OR iDebug = 07) THEN CALL DebugLog ("CHES", "After Gauss")
      IF iSing THEN
         ' Matrix was singular
         ErrorMsg$ = "singular matrix"
         ErrorCode% = TRUE
         EXIT SUB
      END IF
   END IF

   CALL SUpdate(X(), Iter%, MaxError!, Skip%)
   IF (iDebug = 01 OR iDebug = 07) THEN CALL DebugLog ("CHES", "After SUpdate")

   CALL CalcSpeciesConc(ISTR#, ErrorCode%, ErrorMsg$)
   IF (iDebug = 01 OR iDebug = 07) THEN CALL DebugLog ("CHES", "After CalcSpeciesConc")
   IF ErrorCode% THEN
      EXIT SUB
   END IF

   CALL CalcResidual(MaxError!, MaxComp%)
   IF (iDebug = 01 OR iDebug = 07) THEN CALL DebugLog ("CHES", "After CalcResidual")
   IF (iDebug = 01 OR iDebug = 07) THEN CALL DebugLog ("CHES", "Max Error " + Label$(MaxComp%))

   IF NOT iQuiet THEN
      CALL PrintResidual(Iter%, MaxError!, MaxComp%, iSing%)
   END IF

   IF ABS(MaxError) < ErrorTol! THEN Converge% = TRUE

   IF (PeekBuf = 27) THEN
      CALL FileClose
   END IF

   IF Converge% THEN
      CALL CheckSolidPhase(Converge%, UnderFlag%, OverFlag%)
      IF (iDebug = 01 OR iDebug = 07) THEN CALL DebugLog ("CHES", "After CheckSolidPhase")

      ' If Convergence flag is set back to FALSE, it must be due to
      '    change in the list of active solid phases.
      IF NOT Converge% THEN CALL SSaturate
   END IF

   IF (Iter% > MaxIter%) THEN
      ErrorCode% = TRUE
      ErrorMsg$ = "Iterations exceeded max."
      EXIT SUB
   END IF

LOOP UNTIL (Converge% AND NOT Skip%)

CALL ChargeBalance(CHARGE#, ISTR#)
IF (iDebug = 01 OR iDebug = 07) THEN CALL DebugLog ("CHES", "After ChargeBalance")

OC(0) = CDBL(Iter%)
FOR i = 2 TO NWHAM%
   ii = WHAMMAP%(i)
   IF (ii <> 0) THEN
      KT(i) = wT(ii) * CtoM!(i)
      'wA(ii) = wT(ii)
   END IF
NEXT i

END SUB

DEFINT A-Z
FUNCTION Davies# (ISTR AS DOUBLE, ZZ AS SINGLE)
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  FUNCTION Davies                                                      |
'  |                                                                       |
'  |  Calculates activity coefficients using the Davies equation.          |
'  |                                                                       |
'  |  Subroutine Parameters:                                               |
'  |                                                                       |
'  |     IStr   (double, input)   Ionic strength, (moles/liter)            |
'  |                                                                       |
'  |     ZZ     (single, input)   Ionic charge                             |
'  |                                                                       |
'  |     Davies (double, output)  Davies activity coefficient              |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

IF (ZZ = 0) OR (ISTR# <= 0) THEN
   ' *------------------------------------------------------------*
   ' | Davies coefficient is undefined for neutral species        |
   ' | For zero ionic strength systems, the activity coefficient  |
   ' | is one (1) by definition                                   |
   ' *------------------------------------------------------------*
   Davies# = 1#
   EXIT FUNCTION
END IF

' *---------------------------------------------*
' |  Simple Variables Local to FUNCTION Davies  |
' *---------------------------------------------*
DIM SqrIStr AS DOUBLE       '  Square root of ionic strength
DIM NewZ AS SINGLE          '  Absolute value of incoming ionic charge
DIM Temp AS DOUBLE          '  Temporary double precision variable

' *----------------------------------------*
' |  Here we begin our Davies calculation  |
' |  for a given IStr and ionic charge z   |
' *----------------------------------------*
   NewZ! = ABS(ZZ)
   SqrIStr# = SQR(ISTR#)
   Temp# = -deA * NewZ! ^ TWO * ((SqrIStr# / (ONE + SqrIStr#)) - deB * ISTR#)
   IF Temp# > THREE THEN
      ' *---------------------------------------*
      ' |  Avoid overflows                      |
      ' |  Just ignore numbers this big, they   |
      ' |  typically result from wild guesses   |
      ' |  for component concentrations far     |
      ' |  from the final equilibrium solution. |
      ' *---------------------------------------*
      Temp# = ONE
   ELSE
      Temp# = TEN ^ Temp#
   END IF

' *----------------*
' |  Return Value  |
' *----------------*
Davies# = Temp#

END FUNCTION

DEFSNG A-Z
FUNCTION DebyeHuckel# (ISTR#, ZZ!)
'
'  Extended Debye-Huckel calculation from WHAM
'

ACTA = .27 + (.0008 * TEMPK)                          ' const A in D-H
ACTB = .33                                            ' const B in D-H
ActC = -ACTA * SQR(ISTR#)
ACTD = ACTB * SQR(ISTR#)
ACTE = ABS(ZZ)
ACTF = ZZ * ZZ

SELECT CASE ACTE
   CASE 0
      DebyeHuckel# = 1
   CASE 1
      DebyeHuckel# = 10 ^ (ActC / (1 + (3 * ACTD)))             ' act coeff M+/-
   CASE 2
      DebyeHuckel# = 10 ^ (4 * ActC / (1 + (6 * ACTD)))         ' act coeff M2+/-
   CASE 3
      DebyeHuckel# = 10 ^ (9 * ActC / (1 + (9 * ACTD)))         ' act coeff M3+/-
   CASE 4
      DebyeHuckel# = 10 ^ (16 * ActC / (1 + (12 * ACTD)))       ' act coeff M4+/-
   CASE ELSE
      DebyeHuckel# = 10 ^ (ACTF * ActC / (1 + (3 * ACTE * ACTD)))
END SELECT

END FUNCTION

DEFINT A-Z

FUNCTION Delimit%(Work$, Delim$)
   '
   ' Count the number of delimited fields in Work$ while ignoring
   ' consecutive delimiters.  Also keep fields within quotes intact.
   '
   DIM Last AS INTEGER
   DIM BeginQuote AS INTEGER
   FALSE = 0
   TRUE = NOT FALSE

   Last = FALSE
   FOR X% = 1 TO LEN(Work$)
       T$ = MID$(Work$, X%, 1)
       IF(INSTR(T$, CHR$(34))) THEN
          ' Quotes are special, ignore anything between quotes
          IF BeginQuote = FALSE THEN
             BeginQuote = TRUE
          ELSE
             BeginQuote = FALSE
             Counter% = Counter% + 1
          END IF
       ELSE
          IF BeginQuote = FALSE THEN
             IF(INSTR(MID$(Work$, X%, 1), ANY Delim$)) THEN
                ' This character is a delimiter.
                ' Was the last character a delimiter as well?
                IF Last = TRUE THEN
                   ' last character was delimiter, skip this one
                ELSE
                   Counter% = Counter% + 1
                   Last = TRUE
                END IF
                IF (X% = LEN(Work$)) THEN
                   '
                   ' We are at the end of the string, so it dosn't
                   ' matter if the last character(s) is(are)
                   ' delimiters.
                   Counter% = Counter% - 1
                END IF
             ELSE
                ' This character is not a delimiter
                Last = FALSE
                'IF X% = LEN(Work$) THEN
                '   ' this is the last character, if the whole string
                '   ' doesn't end in a delimiter we still want to count
                '   ' the last field
                '   Counter% = Counter% + 1
                'END IF
             END IF
          END IF
       END IF
   NEXT X%
   Delimit% = Counter%
END FUNCTION

SUB DefineProblem300 (FName$, ErrorCode, ErrorMsg$)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB DefineProblem300                                                 |
'  |                                                                       |
'  |  Reads in a ver 3.00 file for use with the CHESS model.  Although a   |
'  |  large number of variables are initialized in this subroutine, most   |
'  |  are defined as global variables.  See the variable lists for         |
'  |  specific details.                                                    |
'  |                                                                       |
'  |  Subroutine Parameters:                                               |
'  |                                                                       |
'  |     FName$  (string, input)    Name of input file to read          |
'  |                                                                       |
'  |     ErrorCode% (integer, output)  Non-zero value indicates an error   |
'  |                                   occured                             |
'  |                                                                       |
'  |     ErrorMsg$  (string, output)   Text of an error msg for reporting  |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

' $INCLUDE: 'E:\Users\BLM_Devl\dp_dims.bi'

' *--------------------------------------*
' |  Extra dimension statements for BLM  |
' *--------------------------------------*
   'DIM sMetBL(10) AS STRING
   'DIM iBLM(10) AS INTEGER
   'DIM iBLMTox(10) AS INTEGER

' *-----------------------------------------------*
' |  Open input file                              |
' |                                               |
' |  Note that no attempt is made to verify that  |
' |  the file exists                              |
' *-----------------------------------------------*
iFileIn = FREEFILE
OPEN FName$ FOR INPUT AS #iFileIn

' *---------------------------------------*
' |  Check the file header to make sure   |
' |  this is the right type of data file  |
' *---------------------------------------*
LINE INPUT #iFileIn, TmpS$
TmpI = INSTR(TmpS$, "Column model parameter file")
IF TmpI = 0 THEN
   ErrorCode = TRUE
   ErrorMsg$ = FName$ + " is not a CHESS/COLUMN data file."
   EXIT SUB
END IF

' *---------------------------------------------*
' |  Read the version of the input file and     |
' |  split into major and minor verion numbers  |
' *---------------------------------------------*
TmpI = INSTR(TmpS$, "Ver ")
MajorVer = VAL(MID$(TmpS$, TmpI + 4, 1))
MinorVer = VAL(MID$(TmpS$, TmpI + 6, 2))

' *---------------------------------------------------------*
' |  Strip out some 'informative' text from the input file  |
' *---------------------------------------------------------*
LINE INPUT #iFileIn, TmpS$
LINE INPUT #iFileIn, TmpS$

' *----------------------------------*
' |  Read in the problem dimensions  |
' *----------------------------------*
INPUT #iFileIn, NComp, NSpecies, NPhase, NList
NSpecies = NSpecies + NComp
NFirstS = NComp + 1

' *------------------------*
' |  Check max dimensions  |
' *------------------------*
IF (NComp > MaxComp) THEN
   ErrorCode = TRUE
   ErrorMsg$ = "Number of components is greater than maximum."
END IF
IF (NSpecies > MaxSpecies) THEN
   ErrorCode = TRUE
   ErrorMsg$ = "Number of species is greater than maximum."
END IF
IF (NPhase > MaxPhases) THEN
   ErrorCode = TRUE
   ErrorMsg$ = "Number of phases is greater than maximum."
END IF
IF (NList > MaxLists) THEN
   ErrorCode = TRUE
   ErrorMsg$ = "Number of lists is greater than maximum."
END IF

' *---------------------------------------------*
' |  Avoid errors in sims where NSpecies=NComp  |
' *---------------------------------------------*
IF NFirstS > NSpecies THEN
   NFirstS = NSpecies
END IF

' *---------------------------------------------------------*
' |  Strip out some 'informative' text from the input file  |
' *---------------------------------------------------------*
LINE INPUT #iFileIn, TmpS$

' *----------------------------------------*
' |  Read in the list of Components, plus  |
' |  descriptive information about them    |
' *----------------------------------------*
FOR i = 1 TO NComp
   INPUT #iFileIn, TmpS$, Tmp!, TmpI2

   ' *------------------*
   ' |  Component name  |
   ' *------------------*
   Label$(i) = QPTrim$(TmpS$)

   ' *--------------------*
   ' |  Component charge  |
   ' *--------------------*
   V(i) = Tmp!

   ' *------------------------------*
   ' |  Species and Component Type  |
   ' *------------------------------*
   SType(i) = TmpI2 \ TEN
   CType(i) = TmpI2 MOD TEN

   ' *-----------------------*
   ' |  Activity correction  |
   ' *-----------------------*
   INPUT #iFileIn, ActCorr(i)

   ' *------------------------------*
   ' |  Report this error and exit  |
   ' *------------------------------*
   IF (ActCorr(i) = acChargeFraction) AND (V(i) = ZERO) THEN
      ErrorCode = TRUE
      TmpS$ = "A charge of zero was specified for a charge fraction component ("
      TmpS$ = TmpS$ + Label(i) + " )."
      ErrorMsg$ = TmpS$
      EXIT SUB
   END IF

   INPUT #iFileIn, SiteDen(i)
NEXT i

' *--------------------------------------------------------*
' |  For every component except P (CTSURFPOT), the       |
' |  stoichiometric coefficients for that component        |
' |  in it's own mole balance and mass actions equations=1 |
' *--------------------------------------------------------*
FOR c = 1 TO NComp
   IF (CType(c) <> ctSurfPot) THEN
      SCa(c, c) = ONE
      SCb(c, c) = ONE
      EC(c) = ONE
   ELSE
      SCa(c, c) = ZERO
      SCb(c, c) = ZERO
      EC(c) = ZERO
   END IF
NEXT c

' *---------------------------------------------------------*
' |  Strip out some 'informative' text from the input file  |
' *---------------------------------------------------------*
LINE INPUT #iFileIn, TmpS$
LINE INPUT #iFileIn, TmpS$

' *----------------------------------------------------*
' |  Load the species info into appropriate variables  |
' *----------------------------------------------------*
FOR i = NFirstS TO NSpecies
   INPUT #iFileIn, TmpS$, TmpI1, TmpI2

   ' *----------------*
   ' |  Species Name  |
   ' *----------------*
   Label$(i) = QPTrim$(TmpS$)

   ' *----------------*
   ' |  Species Type  |
   ' *----------------*
   SType(i) = TmpI1 \ TEN

   ' *------------------------------*
   ' |  Activity Correction Method  |
   ' *------------------------------*
   ActCorr(i) = TmpI2

   ' *--------------------------------------------*
   ' |  Flag to indicate use of temp. correction  |
   ' *--------------------------------------------*
   INPUT #iFileIn, TempCorr(i)

   INPUT #iFileIn, iMonteInfo
   IF iMonteInfo THEN
      ' *---------------------------------------------*
      ' |  Store Flag to tell whether to adjust this  |
      ' |  constant during Monte Carlo iterations     |
      ' |  and add to Monte Carlo list if appropriate |
      ' *---------------------------------------------*
      NMonteVars = NMonteVars + 1
      MonteList(NMonteVars) = i
   END IF

   ' *----------------------------------*
   ' |  Read in reaction stoichiometry  |
   ' *----------------------------------*
   FOR c = 1 TO NComp
      INPUT #iFileIn, TmpS$
      iPos = INSTR(TmpS$, ":")
      ' *------------------------------------------------*
      ' |  The presence of a colon in the SC field will  |
      ' |  indicate the use of different coefficients    |
      ' |  for mass action and mole balance              |
      ' *------------------------------------------------*
      IF iPos THEN
         ' Two types of coefficients. Split them apart
         SCa(i, c) = VAL(MID$(TmpS$, 1, iPos - 1))
         SCb(i, c) = VAL(MID$(TmpS$, iPos + 1, LEN(TmpS$) - iPos))
      ELSE
         ' No colon, both types of coefficients are the same
         Tmp! = VAL(TmpS$)
         SCa(i, c) = Tmp!
         SCb(i, c) = Tmp!
      END IF
   NEXT c

   ' *-----------------------------------------------*
   ' |  Input the log K and convert to non-Log form  |
   ' *-----------------------------------------------*
   INPUT #iFileIn, Tmp!

   ' *----------------------------------------------------*
   ' |  Input the variation on this Equilibrium constant  |
   ' *----------------------------------------------------*
   INPUT #iFileIn, Tmp2!

   ' *------------------------------------------------------*
   ' |  Store Eq. constant for species formation reactions  |
   ' *------------------------------------------------------*
   EC(i) = 10 ^ Tmp!

   IF iMonteInfo THEN
      ' *--------------------------*
      ' |  Store Monte Carlo Info  |
      ' *--------------------------*
      MeanLogK(i) = 10 ^ Tmp!
      VarLogK(i) = 10 ^ Tmp2!
   END IF

   ' *------------------------------------------------------*
   ' |  Delta enthalpy values and Temp that K was measured  |
   ' *------------------------------------------------------*
   INPUT #iFileIn, DeltaH(i), ECTemp(i)

   ' *----------------------------------------------------------*
   ' |  Finally, get the initial concentration of this species  |
   ' *----------------------------------------------------------*
   INPUT #iFileIn, InitialConc(i)
NEXT i

' *-----------------------------------------------------------------*
' |  Get linked species info for charge fraction and mole fraction  |
' |  activity corrections                                           |
' *-----------------------------------------------------------------*
    ' *---------------------------------------------------------*
    ' |  Strip out some 'informative' text from the input file  |
    ' *---------------------------------------------------------*
   LINE INPUT #iFileIn, TmpS$
   LINE INPUT #iFileIn, TmpS$

   FOR i = 1 TO NList
      INPUT #iFileIn, Nj

      FOR j = 1 TO Nj
         INPUT #iFileIn, TmpI
         LinkList(i, j) = TmpI

         IF LinkComp(TmpI) = 0 THEN
            LinkComp(TmpI) = LinkList(i, 1)
         ELSE
            TmpS$ = "The species " + Label$(TmpI) + " appears in more than one linked list."
            ErrorCode = TRUE
            EXIT SUB
         END IF

      NEXT j
      LinkList(i, 0) = Nj
      ListPos(LinkList(i, 1)) = i
   NEXT i


' *-------------------------------------------*
' |  Now calculate the charge on all species  |
' *-------------------------------------------*
FOR i = NFirstS TO NSpecies
   TmpD# = 0
   FOR c = 1 TO NComp
      TmpD# = TmpD# + SCa(i, c) * V(c)
   NEXT c
   V(i) = TmpD#
NEXT i

' *---------------------------------------------------------*
' |  Strip out some 'informative' text from the input file  |
' *---------------------------------------------------------*
LINE INPUT #iFileIn, TmpS$
LINE INPUT #iFileIn, TmpS$

' *----------------------------------------------------------------*
' |  Read in phase information and store in appropriate variables  |
' *----------------------------------------------------------------*
FOR i = 1 TO NPhase
   INPUT #iFileIn, TmpS$
   Label$(i + NSpecies) = QPTrim$(TmpS$)

   INPUT #iFileIn, STempCorr(i)

   INPUT #iFileIn, SP(i, 0)  ' this will be the Monte Carlo flag for phases

   FOR c = 1 TO NComp
      INPUT #iFileIn, SP(i, c)
   NEXT c

   INPUT #iFileIn, Tmp!, Tmp2!  ' Tmp2! is the MC var K for phases
   SEC(i) = 10 ^ Tmp!

   INPUT #iFileIn, Tmp!, Tmp2!
   SDeltaH(i) = Tmp!
   SECTemp(i) = Tmp2!


   INPUT #iFileIn, Tmp!
   PhaseMass(i) = Tmp!

   IF Tmp! = 0! THEN
      PType(i) = ptNotAvailable
   ELSE
      PType(i) = ptSaturated
   END IF
NEXT i

' *---------------------------------------------------------*
' |  Strip out some 'informative' text from the input file  |
' *---------------------------------------------------------*
LINE INPUT #iFileIn, TmpS$
LINE INPUT #iFileIn, TmpS$
LINE INPUT #iFileIn, TmpS$


' *-----------------------------------------------*
' |  Now for the file I/O and system description  |
' *-----------------------------------------------*
i = 0
NBLM = 0
iDone = FALSE
DO
   ExitFlag1 = FALSE
   DO   ' Keep grabbing lines until we get one with data
       ' *---------------------------------------------------------*
       ' |  Strip out some 'informative' text from the input file  |
       ' *---------------------------------------------------------*
      LINE INPUT #iFileIn, TmpS$
      OldTmpS$ = TmpS$
      TmpS$ = QPTrim$(LEFT$(TmpS$, 30))
      IF LEN(TmpS$) > 0 THEN ExitFlag1 = TRUE
   LOOP UNTIL ExitFlag1
   i = i + 1
   SELECT CASE i
      CASE 1         ' Input file name
         InputFile$ = QPTrim$(OldTmpS$)
      CASE 2         ' Input file type
         InputFileType = VAL(TmpS$)
      CASE 3         ' Output file name
         OutputFile$ = QPTrim$(OldTmpS$)
         iPos = INSTR(OutputFile$, "??")
         IF iPos <> 0 THEN
            OutputFile$ = MID$(OutputFile$, 1, iPos - 1)
         END IF
      CASE 4         ' Output file type
         OutputFileType = VAL(TmpS$)
      CASE 5         ' Simulation type
         SimType = VAL(TmpS$)
      CASE 6         ' Number of layers
         NLayers = VAL(TmpS$)
      CASE 7         ' Number of Monte Carlo Simulations
         NSims = VAL(TmpS$)
      CASE 8         ' Error Tolerance
         ErrorTol = VAL(TmpS$)
      CASE 9         ' Maximum number of iterations
         MaxIter = VAL(TmpS$)

      ' *----------------------------------*
      ' |  Determine Output file contents  |
      ' *----------------------------------*
      CASE 10        ' Concentrations?
         OFC.Conc = VAL(TmpS$)
      CASE 11        ' Totals?
         OFC.Tot = VAL(TmpS$)
      CASE 12        ' Equil. Constants?
         OFC.Keq = VAL(TmpS$)
      CASE 13        ' Residuals?
         OFC.Resid = VAL(TmpS$)
      CASE 14        ' Activity Coefficients?
         OFC.ActC = VAL(TmpS$)
      CASE 15        ' Phase Saturation Indices?
         OFC.PhaseSI = VAL(TmpS$)
      CASE 16        ' Phase Ion Activity Products?
         OFC.PhaseIAP = VAL(TmpS$)
      CASE ELSE
         iToxic = TRUE
         IF (INSTR(TmpS$, "[")) THEN
            p1 = INSTR(TmpS$, "[")
            p2 = INSTR(TmpS$, "]")
            TmpS2$ = MID$(TmpS$, p1 + 1, p2 - p1 - 1)
            TmpS3$ = RTRIM$(LTRIM$(RIGHT$(TmpS$, LEN(TmpS$) - p2 - 1)))
            IF (INSTR(TmpS$, "/N")) THEN
               p2 = INSTR(TmpS3$, "/N")
               TmpS3$ = RTRIM$(LTRIM$(LEFT$(TmpS3$, p2 - 1)))
               iToxic = FALSE
            END IF
         END IF
         SELECT CASE TmpS2$
            CASE "Metal"
               sMetal$ = TmpS3$
            CASE "Gill", "BL"
               sBL$ = TmpS3$
            CASE "Gill-Metal", "BL-Metal"
               NBLM = NBLM + 1
               sMetBL$(NBLM) = TmpS3$
               iBLMTox(NBLM) = iToxic
            CASE "DOC"
               sDOC$ = TmpS3$
            CASE "THERMO"
               sWHAM$ = TmpS3$
            CASE "CRITICAL", "LA50"
               CritBL = VAL(TmpS3$)
            CASE "LA50 ADJUST"
               iLA50Adjust = 1
               ' Most of these strings were truncated at 30 characters, to allow comments in
               ' the parameter file beyond 30 characters.  Since we need longer strings to
               ' store this info, let TmpS3$ revert back to the full length string in the
               ' parameter file.
               TmpS3$ = RTRIM$(LTRIM$(RIGHT$(OldTmpS$, LEN(OldTmpS$) - p2 - 1)))
               '
               ' Now parse it into individual parameters
               Delim$ = ","
               CALL NParse(TmpS3$, Delim$, LA50PARAMS$(), NLA50Params)
               IF LA50PARAMS$(1) = "IONICSTR" THEN
                  iLA50Variable = NSpecies + 2
               ELSE
                  iLA50Variable = FindSpecies(LA50PARAMS$(1))
               END IF
               LA50Slope = CDBL(VAL(LA50PARAMS$(2))) '* OKT(iBL) / 1000000!
               LA50Break = CDBL(VAL(LA50PARAMS$(3)))
               LA50MinFrac = CDBL(VAL(LA50PARAMS$(4)))
            CASE "LOG LA50 ADJUST"
               iLA50Adjust = 2
               ' Most of these strings were truncated at 30 characters, to allow comments in
               ' the parameter file beyond 30 characters.  Since we need longer strings to
               ' store this info, let TmpS3$ revert back to the full length string in the
               ' parameter file.
               TmpS3$ = RTRIM$(LTRIM$(RIGHT$(OldTmpS$, LEN(OldTmpS$) - p2 - 1)))
               '
               ' Now parse it into individual parameters
               Delim$ = ","
               CALL NParse(TmpS3$, Delim$, LA50PARAMS$(), NLA50Params)
               IF LA50PARAMS$(1) = "IONICSTR" THEN
                  iLA50Variable = NSpecies + 2
               ELSE
                  iLA50Variable = FindSpecies(LA50PARAMS$(1))
               END IF
               LA50Slope = CDBL(VAL(LA50PARAMS$(2))) '* OKT(iBL) / 1000000!
               LA50Break = CDBL(VAL(LA50PARAMS$(3)))
               'LA50MinFrac = CDBL(VAL(LA50PARAMS$(4)))
      '
      '        The following lines are for debugging to see if the line was
      '        parsed properly.
      '
      '        Call Pmsg(TmpS3$)
      '        CALL Pmsg(LA50PARAMS$(1))   ' the parameter that will be used to adjust
      '        CALL Pmsg(LA50PARAMS$(2))   ' the slope of LA50 to the parameter
      '        CALL Pmsg(LA50PARAMS$(3))   ' the minimum value of paramter below which LA50 is adjusted
      '        CALL Pmsg(LA50PARAMS$(4))   ' the minimum fraction of the default LA50 that is allowed after adjustment
      '        CALL Pmsg(LA50PARAMS$(5))
      '        Call Pmsg(Label(iLA50Variable))
            CASE "END"
               iDone = TRUE
            CASE ELSE
         END SELECT
   END SELECT

'LOOP UNTIL i = 15
'LOOP UNTIL i = 23
LOOP UNTIL iDone
CLOSE #iFileIn

' *-------------------------------------------------------------*
' |  Now that all the data has been read in, we can work on it  |
' *-------------------------------------------------------------*

' *------------------------------------*
' |  Layers only apply to column sims  |
' *------------------------------------*
IF NOT SoilCol THEN
   NLayers = 1
END IF

' *--------------------------------------------------------------*
' |  Set the output file name to the largest base we can afford  |
' *--------------------------------------------------------------*
MaxLen = 8
IF (NLayers > 1) THEN
   MaxLen = MaxLen - 1
END IF

IF (NSims > 0) THEN
   MaxLen = MaxLen - 3
END IF

IF (LEN(OutputFile$) > MaxLen) THEN
   OutputFile$ = LEFT$(OutputFile$, MaxLen)
END IF

IF NPhase > 0 THEN
   CALL SelectComp(SComp(), SPhase(), ErrorCode%, ErrorMsg$)
END IF

CALL SSaturate

CALL SimplifyProblem

END SUB


SUB DefineProblem300W (FName$, ErrorCode, ErrorMsg$)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB DefineProblem300W                                                 |
'  |                                                                       |
'  |  Reads in a ver 3.00 file for use with the CHESS model.  Although a   |
'  |  large number of variables are initialized in this subroutine, most   |
'  |  are defined as global variables.  See the variable lists for         |
'  |  specific details.                                                    |
'  |                                                                       |
'  |  Subroutine Parameters:                                               |
'  |                                                                       |
'  |     FName$  (string, input)    Name of input file to read          |
'  |                                                                       |
'  |     ErrorCode% (integer, output)  Non-zero value indicates an error   |
'  |                                   occured                             |
'  |                                                                       |
'  |     ErrorMsg$  (string, output)   Text of an error msg for reporting  |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

' *-----------------------------------------------*
' |  Simple Variables Local to SUB DefineProblem  |
' *-----------------------------------------------*
DIM c AS INTEGER            '  Array/loop index variable
DIM ExitFlag1 AS INTEGER    '  Boolean flag used to test a conditional loop
'DIM i AS INTEGER            '  Array/loop index variable
DIM iFileIn AS INTEGER      '  I/O handle for input file
DIM iMonteInfo AS INTEGER   '  Flag to indicate this species will Monte Carlo
DIM iPos AS INTEGER         '  Position locator of one string within another
DIM j AS INTEGER            '  Array/loop index variable
DIM Nj AS INTEGER           '  Number of elements in a linked/list
DIM Tmp AS SINGLE           '  Temporary single precision variable
DIM TmpD AS DOUBLE          '  Temporary double precision variable
DIM Tmp2 AS SINGLE          '  Temporary single precision variable
DIM TmpS AS STRING          '  Temporary string variable
DIM TmpS2 AS STRING         '  Temporary string variable
DIM TmpI AS INTEGER         '  Temporary integer variable
DIM TmpI1 AS INTEGER        '  Temporary integer variable
DIM TmpI2 AS INTEGER        '  Temporary integer variable
DIM NFields AS INTEGER          ' Number of data fields on input line
DIM MaxFields AS INTEGER        ' Maximum number of data fields expected
DIM DField$(1 TO 35)            ' input line parsed into delimited fields
DIM Delim AS STRING             ' string that holds delimiter characters

' $INCLUDE: 'E:\Users\BLM_Devl\dp_dims.bi'
IF (iWHAM) THEN
   'DIM wOrgBound!(MaxSpecies)

END IF

' *-----------------------------------------------*
' |  Open input file                              |
' |                                               |
' |  Note that no attempt is made to verify that  |
' |  the file exists                              |
' *-----------------------------------------------*
iFileIn = FREEFILE
OPEN FName$ FOR INPUT AS #iFileIn

' *---------------------------------------*
' |  Check the file header to make sure   |
' |  this is the right type of data file  |
' *---------------------------------------*
LINE INPUT #iFileIn, TmpS$
TmpI = INSTR(TmpS$, "Column model parameter file")
IF TmpI = 0 THEN
   ErrorCode = TRUE
   ErrorMsg$ = FName$ + " is not a CHESS/COLUMN data file."
   EXIT SUB
END IF

' *---------------------------------------------*
' |  Read the version of the input file and     |
' |  split into major and minor verion numbers  |
' *---------------------------------------------*
TmpI = INSTR(TmpS$, "Ver ")
MajorVer = VAL(MID$(TmpS$, TmpI + 4, 1))
MinorVer = VAL(MID$(TmpS$, TmpI + 6, 2))

' *---------------------------------------------------------*
' |  Strip out some 'informative' text from the input file  |
' *---------------------------------------------------------*
LINE INPUT #iFileIn, TmpS$
LINE INPUT #iFileIn, TmpS$

' *----------------------------------*
' |  Read in the problem dimensions  |
' *----------------------------------*
IF (MajorVer >=3) AND (MinorVer >= 10) THEN
    '  For vers 3.10 and above there are
    ' variable numbers of mass compartments
    INPUT #iFileIn, NMass, NComp, NSpecies, NPhase, NList
ELSE
    INPUT #iFileIn, NComp, NSpecies, NPhase, NList
END IF
NSpecies = NSpecies + NComp
NFirstS = NComp + 1

' *------------------------*
' |  Check max dimensions  |
' *------------------------*
IF (NMass > MaxMass) THEN
   ErrorCode = TRUE
   ErrorMsg$ = "Number of mass compartments is greater than maximum."
END IF
IF (NComp > MaxComp) THEN
   ErrorCode = TRUE
   ErrorMsg$ = "Number of components is greater than maximum."
END IF
IF (NSpecies > MaxSpecies) THEN
   ErrorCode = TRUE
   ErrorMsg$ = "Number of species is greater than maximum."
END IF
IF (NPhase > MaxPhases) THEN
   ErrorCode = TRUE
   ErrorMsg$ = "Number of phases is greater than maximum."
END IF
IF (NList > MaxLists) THEN
   ErrorCode = TRUE
   ErrorMsg$ = "Number of lists is greater than maximum."
END IF

' *---------------------------------------------*
' |  Avoid errors in sims where NSpecies=NComp  |
' *---------------------------------------------*
IF NFirstS > NSpecies THEN
   NFirstS = NSpecies
END IF

IF (MajorVer >=4) AND (MinorVer >= 10) THEN  '<--- Note I changed the major ver to 4 to allow 3.10 to have a critical accumulation table
    '  For vers 3.10 and above there are
    ' variable numbers of mass compartments
    ' *---------------------------------------------------------*
    ' |  Strip out some 'informative' text from the input file  |
    ' *---------------------------------------------------------*
    LINE INPUT #iFileIn, TmpS$
    LINE INPUT #iFileIn, TmpS$
    FOR i = 1 TO NMass
        INPUT #iFileIn, TmpS$, Tmp!, TmpS2$
        MassName(i) = QPTrim$(TmpS$)            ' The name of a mass compartment
        MassVal(i) = Tmp!                       ' The mass value
        MassUnitLabel(i) = QPTrim$(TmpS2$)       ' The units of a mass compartment (e.g. L or kg)
    NEXT i
    ' *---------------------------------------------------------*
    ' |  Strip out some 'informative' text from the input file  |
    ' *---------------------------------------------------------*
    LINE INPUT #iFileIn, TmpS$
    LINE INPUT #iFileIn, TmpS$
ELSE
    ' *---------------------------------------------------------*
    ' |  Strip out some 'informative' text from the input file  |
    ' *---------------------------------------------------------*
    LINE INPUT #iFileIn, TmpS$
    ' *-------------------------------------------------------------------------------------*
    ' |  For earlier versions, this info is not in the parameter file, so hardwire it here  |
    ' *-------------------------------------------------------------------------------------*
    NMass = 3
    MassName(1) = "Water"                   ' The compartment name
    MassVal(1) = 1.0                        ' The mass value
    MassUnitLabel(1) = "L"                  ' The units of a mass compartment (e.g. L or kg)

    MassName(2) = "BL"                      ' The compartment name
    MassVal(2) = 1.0                        ' The mass value
    MassUnitLabel(2) = "kg wet"             ' The units of a mass compartment (e.g. L or kg)

    MassName(3) = "Water"                   ' The compartment name
    MassVal(3) = 1.0                        ' The mass value
    MassUnitLabel(3) = "L"                  ' The units of a mass compartment (e.g. L or kg)

END IF

' *----------------------------------------*
' |  Read in the list of Components, plus  |
' |  descriptive information about them    |
' *----------------------------------------*
FOR i = 1 TO NComp
   INPUT #iFileIn, TmpS$, Tmp!, TmpI2

   ' *------------------*
   ' |  Component name  |
   ' *------------------*
   Label$(i) = QPTrim$(TmpS$)

   ' *--------------------*
   ' |  Component charge  |
   ' *--------------------*
   V(i) = Tmp!

   ' *------------------------------*
   ' |  Species and Component Type  |
   ' *------------------------------*
   SType(i) = TmpI2 \ TEN
   CType(i) = TmpI2 MOD TEN

   ' *-----------------------*
   ' |  Activity correction  |
   ' *-----------------------*
   INPUT #iFileIn, ActCorr(i)

   ' *------------------------------*
   ' |  Report this error and exit  |
   ' *------------------------------*
   IF (ActCorr(i) = acChargeFraction) AND (V(i) = ZERO) THEN
      ErrorCode = TRUE
      TmpS$ = "A charge of zero was specified for a charge fraction component ("
      TmpS$ = TmpS$ + Label(i) + " )."
      ErrorMsg$ = TmpS$
      EXIT SUB
   END IF

   INPUT #iFileIn, SiteDen(i)
NEXT i

' *--------------------------------------------------------*
' |  For every component except P (CTSURFPOT), the       |
' |  stoichiometric coefficients for that component        |
' |  in it's own mole balance and mass actions equations=1 |
' *--------------------------------------------------------*
FOR c = 1 TO NComp
   IF (CType(c) <> ctSurfPot) THEN
      SCa(c, c) = ONE
      SCb(c, c) = ONE
      EC(c) = ONE
   ELSE
      SCa(c, c) = ZERO
      SCb(c, c) = ZERO
      EC(c) = ZERO
   END IF
NEXT c

'iAl = FindSpecies%("Al")
'iAlp = FindSpecies%("Al(OH)3(s)")
'IF (iAl <> 0) AND (iAlp <> 0) THEN
'    'break
'    LKSOAL25 = 8.625
'    'adj up to slightly higher Al solubility as suggested by Niva data
'    LKSOAL25 = 9.76
'    IF (iLogKsAl = TRUE) THEN LKSOAL25 = LogKsAl     ' set the Log K for Al solubility to command line value
'
'    SCa(iAlp, FindSpecies%("H")) = -1*THREE
'    SCa(iAlp, iAl) = ONE
'    SCa(iAlp, iAlp) = ONE 'This should mathematically cancel the Al(OH3)(s) component out Kso = ([H]^3 * [Al(OH)3(s)]) / ([Al] * [Al(OH)3(s)]) = [H]^3 / [Al]
'    SCb(iAlp, FindSpecies%("H")) = -1*THREE
'    SCb(iAlp, iAl) = ONE
'    SCb(iAlp, iAlp) = ONE 'This should mathematically cancel the Al(OH3)(s) component out Kso = ([H]^3 * [Al(OH)3(s)]) / ([Al] * [Al(OH)3(s)]) = [H]^3 / [Al]
'    EC(iAlp) = LKSOAL25
'    DeltaH(iAlp) = -25000
'    ECTemp(iAlp) = 298
'    InitialConc(iAlp) = ZeroS
'
'end if

' *---------------------------------------------------------*
' |  Strip out some 'informative' text from the input file  |
' *---------------------------------------------------------*
LINE INPUT #iFileIn, TmpS$
LINE INPUT #iFileIn, TmpS$

' *----------------------------------------------------*
' |  Load the species info into appropriate variables  |
' *----------------------------------------------------*
FOR i = NFirstS TO NSpecies
   INPUT #iFileIn, TmpS$, TmpI1, TmpI2

   ' *----------------*
   ' |  Species Name  |
   ' *----------------*
   Label$(i) = QPTrim$(TmpS$)

   ' *----------------*
   ' |  Species Type  |
   ' *----------------*
   SType(i) = TmpI1 \ TEN

   ' *------------------------------*
   ' |  Activity Correction Method  |
   ' *------------------------------*
   ActCorr(i) = TmpI2

   ' *--------------------------------------------*
   ' |  Flag to indicate use of temp. correction  |
   ' *--------------------------------------------*
   INPUT #iFileIn, TempCorr(i)

   INPUT #iFileIn, iMonteInfo
   IF iMonteInfo THEN
      ' *---------------------------------------------*
      ' |  Store Flag to tell whether to adjust this  |
      ' |  constant during Monte Carlo iterations     |
      ' |  and add to Monte Carlo list if appropriate |
      ' *---------------------------------------------*
      NMonteVars = NMonteVars + 1
      MonteList(NMonteVars) = i
   END IF

   ' *----------------------------------*
   ' |  Read in reaction stoichiometry  |
   ' *----------------------------------*
   FOR c = 1 TO NComp
      INPUT #iFileIn, TmpS$
      iPos = INSTR(TmpS$, ":")
      ' *------------------------------------------------*
      ' |  The presence of a colon in the SC field will  |
      ' |  indicate the use of different coefficients    |
      ' |  for mass action and mole balance              |
      ' *------------------------------------------------*
      IF iPos THEN
         ' Two types of coefficients. Split them apart
         SCa(i, c) = VAL(MID$(TmpS$, 1, iPos - 1))
         SCb(i, c) = VAL(MID$(TmpS$, iPos + 1, LEN(TmpS$) - iPos))
      ELSE
         ' No colon, both types of coefficients are the same
         Tmp! = VAL(TmpS$)
         SCa(i, c) = Tmp!
         SCb(i, c) = Tmp!
      END IF
   NEXT c

   ' *-----------------------------------------------*
   ' |  Input the log K and convert to non-Log form  |
   ' *-----------------------------------------------*
   INPUT #iFileIn, Tmp!

   ' *----------------------------------------------------*
   ' |  Input the variation on this Equilibrium constant  |
   ' *----------------------------------------------------*
   INPUT #iFileIn, Tmp2!

   ' *------------------------------------------------------*
   ' |  Store Eq. constant for species formation reactions  |
   ' *------------------------------------------------------*
   EC(i) = 10 ^ Tmp!

   IF iMonteInfo THEN
      ' *--------------------------*
      ' |  Store Monte Carlo Info  |
      ' *--------------------------*
      MeanLogK(i) = 10 ^ Tmp!
      VarLogK(i) = 10 ^ Tmp2!
   END IF

   ' *------------------------------------------------------*
   ' |  Delta enthalpy values and Temp that K was measured  |
   ' *------------------------------------------------------*
   INPUT #iFileIn, DeltaH(i), ECTemp(i)

   ' *----------------------------------------------------------*
   ' |  Finally, get the initial concentration of this species  |
   ' *----------------------------------------------------------*
   INPUT #iFileIn, InitialConc(i)
NEXT i

' *-----------------------------------------------------------------*
' |  Get linked species info for charge fraction and mole fraction  |
' |  activity corrections                                           |
' *-----------------------------------------------------------------*
    ' *---------------------------------------------------------*
    ' |  Strip out some 'informative' text from the input file  |
    ' *---------------------------------------------------------*
   LINE INPUT #iFileIn, TmpS$
   LINE INPUT #iFileIn, TmpS$

   FOR i = 1 TO NList
      INPUT #iFileIn, Nj

      FOR j = 1 TO Nj
         INPUT #iFileIn, TmpI
         LinkList(i, j) = TmpI

         IF LinkComp(TmpI) = 0 THEN
            LinkComp(TmpI) = LinkList(i, 1)
         ELSE
            TmpS$ = "The species " + Label$(TmpI) + " appears in more than one linked list."
            ErrorCode = TRUE
            EXIT SUB
         END IF

      NEXT j
      LinkList(i, 0) = Nj
      ListPos(LinkList(i, 1)) = i
   NEXT i


' *-------------------------------------------*
' |  Now calculate the charge on all species  |
' *-------------------------------------------*
FOR i = NFirstS TO NSpecies
   TmpD# = 0
   FOR c = 1 TO NComp
      TmpD# = TmpD# + SCa(i, c) * V(c)
   NEXT c
   V(i) = TmpD#
NEXT i

' *---------------------------------------------------------*
' |  Strip out some 'informative' text from the input file  |
' *---------------------------------------------------------*
LINE INPUT #iFileIn, TmpS$
LINE INPUT #iFileIn, TmpS$

' *----------------------------------------------------------------*
' |  Read in phase information and store in appropriate variables  |
' *----------------------------------------------------------------*
FOR i = 1 TO NPhase
   INPUT #iFileIn, TmpS$
   Label$(i + NSpecies) = QPTrim$(TmpS$)

   INPUT #iFileIn, STempCorr(i)

   INPUT #iFileIn, SP(i, 0)  ' this will be the Monte Carlo flag for phases

   FOR c = 1 TO NComp
      INPUT #iFileIn, SP(i, c)
   NEXT c

   INPUT #iFileIn, Tmp!, Tmp2!  ' Tmp2! is the MC var K for phases
   SEC(i) = 10 ^ Tmp!

   INPUT #iFileIn, Tmp!, Tmp2!
   SDeltaH(i) = Tmp!
   SECTemp(i) = Tmp2!


   INPUT #iFileIn, Tmp!
   PhaseMass(i) = Tmp!

   IF Tmp! = 0! THEN
      PType(i) = ptNotAvailable
   ELSE
      PType(i) = ptSaturated
   END IF
NEXT i

' *---------------------------------------------------------*
' |  Strip out some 'informative' text from the input file  |
' *---------------------------------------------------------*
LINE INPUT #iFileIn, TmpS$
LINE INPUT #iFileIn, TmpS$
LINE INPUT #iFileIn, TmpS$


' *-----------------------------------------------*
' |  Now for the file I/O and system description  |
' *-----------------------------------------------*
i = 0
NBLM = 0
iDone = FALSE
DO
   ExitFlag1 = FALSE
   DO   ' Keep grabbing lines until we get one with data
       ' *---------------------------------------------------------*
       ' |  Strip out some 'informative' text from the input file  |
       ' *---------------------------------------------------------*
      LINE INPUT #iFileIn, TmpS$
      OldTmpS$ = TmpS$
      TmpS$ = QPTrim$(LEFT$(TmpS$, 30))
      IF LEN(TmpS$) > 0 THEN ExitFlag1 = TRUE
   LOOP UNTIL ExitFlag1
   i = i + 1
   SELECT CASE i
      CASE 1         ' Input file name
         InputFile$ = QPTrim$(OldTmpS$)
      CASE 2         ' Input file type
         InputFileType = VAL(TmpS$)
      CASE 3         ' Output file name
         OutputFile$ = QPTrim$(OldTmpS$)
         iPos = INSTR(OutputFile$, "??")
         IF iPos <> 0 THEN
            OutputFile$ = MID$(OutputFile$, 1, iPos - 1)
         END IF
      CASE 4         ' Output file type
         OutputFileType = VAL(TmpS$)
      CASE 5         ' Simulation type
         SimType = VAL(TmpS$)
      CASE 6         ' Number of layers
         NLayers = VAL(TmpS$)
      CASE 7         ' Number of Monte Carlo Simulations
         NSims = VAL(TmpS$)
      CASE 8         ' Error Tolerance
         ErrorTol = VAL(TmpS$)
      CASE 9         ' Maximum number of iterations
         MaxIter = VAL(TmpS$)

      ' *----------------------------------*
      ' |  Determine Output file contents  |
      ' *----------------------------------*
      CASE 10        ' Concentrations?
         OFC.Conc = VAL(TmpS$)
      CASE 11        ' Totals?
         OFC.Tot = VAL(TmpS$)
      CASE 12        ' Equil. Constants?
         OFC.Keq = VAL(TmpS$)
      CASE 13        ' Residuals?
         OFC.Resid = VAL(TmpS$)
      CASE 14        ' Activity Coefficients?
         OFC.ActC = VAL(TmpS$)
      CASE 15        ' Phase Saturation Indices?
         OFC.PhaseSI = VAL(TmpS$)
      CASE 16        ' Phase Ion Activity Products?
         OFC.PhaseIAP = VAL(TmpS$)
      CASE ELSE
         iToxic = TRUE
         IF (INSTR(TmpS$, "[")) THEN
            p1 = INSTR(TmpS$, "[")
            p2 = INSTR(TmpS$, "]")
            TmpS2$ = MID$(TmpS$, p1 + 1, p2 - p1 - 1)
            TmpS3$ = RTRIM$(LTRIM$(RIGHT$(TmpS$, LEN(TmpS$) - p2 - 1)))
            IF (INSTR(TmpS$, "/N")) THEN
               p2 = INSTR(TmpS3$, "/N")
               TmpS3$ = RTRIM$(LTRIM$(LEFT$(TmpS3$, p2 - 1)))
               iToxic = FALSE
            END IF
         END IF
         SELECT CASE TmpS2$
            CASE "Metal"
               sMetal$ = TmpS3$
            CASE "Gill", "BL"
               sBL$ = TmpS3$
            CASE "Gill-Metal", "BL-Metal"
               NBLM = NBLM + 1
               sMetBL$(NBLM) = TmpS3$
               iBLMTox(NBLM) = iToxic
            CASE "DOC"
               sDOC$ = TmpS3$
            CASE "THERMO"
               sWHAM$ = TmpS3$
            CASE "CRITICAL", "LA50"
               'CritBL = VAL(TmpS3$) 'commented as of 16/3/1
               DefaultCriticalValue = VAL(TmpS3$)
               iHaveCAT = FALSE
            CASE "LA50 ADJUST"
               iLA50Adjust = 1
               ' Most of these strings were truncated at 30 characters, to allow comments in
               ' the parameter file beyond 30 characters.  Since we need longer strings to
               ' store this info, let TmpS3$ revert back to the full length string in the
               ' parameter file.
               TmpS3$ = RTRIM$(LTRIM$(RIGHT$(OldTmpS$, LEN(OldTmpS$) - p2 - 1)))
               '
               ' Now parse it into individual parameters
               Delim$ = ","
               CALL NParse(TmpS3$, Delim$, LA50PARAMS$(), NLA50Params)
               IF LA50PARAMS$(1) = "IONICSTR" THEN
                  iLA50Variable = NSpecies + 2
               ELSE
                  iLA50Variable = FindSpecies(LA50PARAMS$(1))
               END IF
               LA50Slope = CDBL(VAL(LA50PARAMS$(2))) '* OKT(iBL) / 1000000!
               LA50Break = CDBL(VAL(LA50PARAMS$(3)))
               LA50MinFrac = CDBL(VAL(LA50PARAMS$(4)))
            CASE "LOG LA50 ADJUST"
               iLA50Adjust = 2
               ' Most of these strings were truncated at 30 characters, to allow comments in
               ' the parameter file beyond 30 characters.  Since we need longer strings to
               ' store this info, let TmpS3$ revert back to the full length string in the
               ' parameter file.
               TmpS3$ = RTRIM$(LTRIM$(RIGHT$(OldTmpS$, LEN(OldTmpS$) - p2 - 1)))
               '
               ' Now parse it into individual parameters
               Delim$ = ","
               CALL NParse(TmpS3$, Delim$, LA50PARAMS$(), NLA50Params)
               IF LA50PARAMS$(1) = "IONICSTR" THEN
                  iLA50Variable = NSpecies + 2
               ELSE
                  iLA50Variable = FindSpecies(LA50PARAMS$(1))
               END IF
               LA50Slope = CDBL(VAL(LA50PARAMS$(2))) '* OKT(iBL) / 1000000!
               LA50Break = CDBL(VAL(LA50PARAMS$(3)))
               'LA50MinFrac = CDBL(VAL(LA50PARAMS$(4)))
            CASE "DISS_A"
               DISS_A = VAL(TmpS3$)
            CASE "DISS_B"
               DISS_B = VAL(TmpS3$)
            CASE "PART_A"
               PART_A = VAL(TmpS3$)
            CASE "PART_B"
               PART_B = VAL(TmpS3$)
            CASE "PART_H"
               PART_H = VAL(TmpS3$)
            CASE "TARGET_R"
               TARGET_R = VAL(TmpS3$)
            CASE "PART_HI"
               PART_HI = VAL(TmpS3$)
            CASE "PART_HM"
               PART_HM = VAL(TmpS3$)

      '
      '        The following lines are for debugging to see if the line was
      '        parsed properly.
      '
      '        Call Pmsg(TmpS3$)
      '        CALL Pmsg(LA50PARAMS$(1))   ' the parameter that will be used to adjust
      '        CALL Pmsg(LA50PARAMS$(2))   ' the slope of LA50 to the parameter
      '        CALL Pmsg(LA50PARAMS$(3))   ' the minimum value of paramter below which LA50 is adjusted
      '        CALL Pmsg(LA50PARAMS$(4))   ' the minimum fraction of the default LA50 that is allowed after adjustment
      '        CALL Pmsg(LA50PARAMS$(5))
      '        Call Pmsg(Label(iLA50Variable))
            CASE "ACUTE_DIV"
               ACUTE_DIV = VAL(TmpS3$)
            CASE "CHRONIC_DIV"
               CHRONIC_DIV = VAL(TmpS3$)
            CASE "Z_EF"
               Z_EF_Def = VAL(TmpS3$)
            CASE "METAL_DL"
                Metal_DL$ = TmpS3$
            CASE "CRITICAL START"
                ' -----------------------------------------------------------------------------------------------------------
                ' new to ver 2.40: Read in a CAT table for setting critical accumulation to a list of organisms and endpoints
                ' -----------------------------------------------------------------------------------------------------------
               iHaveCAT = TRUE
               LINE INPUT #iFileIn, TmpS2$  ' first line has column headings
               Delim$ = ","
               MaxFields = 8
               ' Keep grabbing lines and parsing into variables until the table
               ' is read in and we get the "CRITICAL END" ending delimiter
               LINE INPUT #iFileIn, TmpS2$
               iCAT_NRow = 0
               DO
                  iCAT_NRow = iCAT_NRow + 1
                  NFields = Delimit%(TmpS2$, Delim$)
                  IF NFields > MaxFields THEN
                      '
                      '  Too many input fields, must be some problem
                      ErrorCode% = TRUE
                      ErrorMsg$ = "Too many data fields in critical accumulation table."
                      EXIT SUB
                  END IF
                  CALL NParse(Tmps2$, Delim$, DField$(), NFields)
                  CAT_CA(iCAT_NRow)            = VAL(DField$(1))
                  CAT_Species(iCAT_NRow)       = DField$(2)
                  CAT_TestType(iCAT_NRow)      = DField$(3)
                  CAT_Lifestage(iCAT_NRow)     = DField$(4)
                  CAT_Endpoint(iCAT_NRow)      = DField$(5)
                  CAT_Quantifier(iCAT_NRow)    = DField$(6)
                  CAT_References(iCAT_NRow)    = DField$(7)
                  CAT_Miscellaneous(iCAT_NRow) = DField$(8)
                  LINE INPUT #iFileIn, TmpS2$
                  TmpS2$=TRIM$(TmpS2$)
               LOOP UNTIL TmpS2$ = "[CRITICAL END]"


            CASE "END"
               iDone = TRUE
            CASE ELSE
         END SELECT
   END SELECT

'LOOP UNTIL i = 15
'LOOP UNTIL i = 23
LOOP UNTIL iDone
CLOSE #iFileIn

' *-------------------------------------------------------------*
' |  Now that all the data has been read in, we can work on it  |
' *-------------------------------------------------------------*

    ' *--------------------------------*
    ' |  Revert to default Z_EF if <0  |
    ' *--------------------------------*
    IF (Z_EF<=0) THEN
       Z_EF = Z_EF_Def
    END IF
    ' *------------------------------------*
    ' |  Layers only apply to column sims  |
    ' *------------------------------------*
    IF NOT SoilCol THEN
       NLayers = 1
    END IF

    ' *--------------------------------------------------------------*
    ' |  Set the output file name to the largest base we can afford  |
    ' *--------------------------------------------------------------*
    MaxLen = 8
    IF (NLayers > 1) THEN
       MaxLen = MaxLen - 1
    END IF

    IF (NSims > 0) THEN
       MaxLen = MaxLen - 3
    END IF

    IF (LEN(OutputFile$) > MaxLen) THEN
       OutputFile$ = LEFT$(OutputFile$, MaxLen)
    END IF


    TmpS$ = UCASE$(sWHAM$)
    IF (TmpS$ = "NONE") OR (TmpS$ = "") THEN
       iWHAM = FALSE
    ELSE
       iWHAM = TRUE
    END IF

'IF (FALSE) THEN
IF (iWHAM) THEN


   ' *-----------------------------------------------*
   ' |  If this simulation includes WHAM components  |
   ' |  increase component number                    |
   ' *-----------------------------------------------*

   iDebugOld = FALSE
   IF iDebugOld THEN
      ii = FREEFILE
      OPEN "debug.out" FOR OUTPUT AS #ii
      FOR i = 1 TO NSpecies
         PRINT #ii, Label$(i); EC(i)
      NEXT i
   END IF

   '
   '  Make room for WHAM components for FA (8) and HA(8) reactions
   '  by moving all species
   '
   WOffSet = 16
   iStartWHAM = NComp
   FOR i = 0 TO (NSpecies - NComp - 1)
      j = NSpecies - i
      k = NSpecies - i + WOffSet
      ActCorr(k) = ActCorr(j)
      ActCoef(k) = ActCoef(j)
      CtoM(k) = CtoM(j)
      DeltaH(k) = DeltaH(j)
      EC(k) = EC(j)
      ECTemp(k) = ECTemp(k)
      hEC(k) = hEC(j)
      InitialConc(k) = InitialConc(j)
      Label(k + NPhases) = Label(j + NPhases)
      LinkComp(k) = LinkComp(j)
      MeanLogK(k) = MeanLogK(j)
      MonteList(k) = MonteList(j)
      OC(k + FOUR) = OC(j + FOUR)
      SType(k) = SType(j)
      tEC(k) = tEC(j)
      TempCorr(k) = TempCorr(j)
      V(k) = V(j)
      VarLogK(k) = VarLogK(j)

      FOR ii = 1 TO NComp
         SCa(k, ii) = SCa(j, ii)
         SCb(k, ii) = SCb(j, ii)
         tSCa(k, ii) = tSCa(j, ii)
         tSCb(k, ii) = tSCb(j, ii)
         SCa(j, ii) = 0
         SCb(j, ii) = 0
         tSCa(j, ii) = 0
         tSCb(j, ii) = 0
      NEXT ii

      FOR ii = 1 TO NList
         LinkList(ii, k) = LinkList(ii, j)
      NEXT ii
   NEXT i

   '
   ' These variables shouln't be affected by the additional components
   '
   'CType(k)
   'CompOrder(k)
   'KT(k)
   'ListPos(k)
   'SiteDen(k)
   'SPhase(k)
   'ChangeMass (NPhases)
   'hSEC (NPhases)
   'tSEC (NPhases)
   'tSP(NPhases, NComp)
   'PhaseMass (NPhases)
   'PhaseOrder (NPhases)
   'PSubOrder (NPhases)
   'PType (NPhases)
   'SComp (NPhases)
   'SDeltaH (NPhases)
   'SEC (NPhases)
   'SECTemp (NPhases)
   'SP(NPhases, NComp)
   'STempCorr (NPhases)

   IF iDebugOld THEN
      FOR i = 1 TO NSpecies
         PRINT #ii, Label$(i); EC(i)
      NEXT i
      CLOSE
      'END
   END IF

   '
   '  Give the WHAM components names
   '
   Label(NComp + 1) = "HA1H"
   Label(NComp + 2) = "HA2H"
   Label(NComp + 3) = "HA3H"
   Label(NComp + 4) = "HA4H"
   Label(NComp + 5) = "HB1H"
   Label(NComp + 6) = "HB2H"
   Label(NComp + 7) = "HB3H"
   Label(NComp + 8) = "HB4H"
   Label(NComp + 9) = "FA1H"
   Label(NComp + 10) = "FA2H"
   Label(NComp + 11) = "FA3H"
   Label(NComp + 12) = "FA4H"
   Label(NComp + 13) = "FB1H"
   Label(NComp + 14) = "FB2H"
   Label(NComp + 15) = "FB3H"
   Label(NComp + 16) = "FB4H"

   '
   '  Save values of the indices of these components
   '  to make it easy to adjust component totals later
   '
   iHA1 = NComp + 1
   iHA2 = NComp + 2
   iHA3 = NComp + 3
   iHA4 = NComp + 4
   iHB1 = NComp + 5
   iHB2 = NComp + 6
   iHB3 = NComp + 7
   iHB4 = NComp + 8
   iFA1 = NComp + 9
   iFA2 = NComp + 10
   iFA3 = NComp + 11
   iFA4 = NComp + 12
   iFB1 = NComp + 13
   iFB2 = NComp + 14
   iFB3 = NComp + 15
   iFB4 = NComp + 16

   FOR i = NComp + 1 TO NComp + WOffSet
      '
      ' Component Charge
      '
      V(i) = ZERO

      '
      ' Species and Component Type
      '
      SType(i) = stAqueous
      CType(i) = ctFloat

      '
      ' Activity correction
      '
      ActCorr(i) = acWHAM

      '
      '  Moles of sites per mole C
      '
      SiteDen(i) = ONE

      '
      '  Equilibrium constant and Stoichiometry
      '  are trivial for components by definition
      '
      EC(i) = ONE
      SCa(i, i) = ONE
      SCb(i, i) = ONE
      KT(c) = .0001
   NEXT i


   '
   '  Adjust dimensions to reflect these additions
   '
   NComp = NComp + WOffSet
   NSpecies = NSpecies + WOffSet
   NFirstS = NComp + 1
   FOR i = 1 TO NList
      FOR j = 1 TO LinkList(i, 0)
         IF LinkList(i, J) > (NComp - WOffSet) THEN
             ' only offset species, with speciesnumbers more than the number of components before adding WHAM components
            LinkList(i, j) = LinkList(i, j) + WOffSet
         END IF
      NEXT j
   NEXT i
   '
   ' Open WHAM database and read in data
   '
   '
   ' Check to see if the file is there
   IF NOT Exist(sWHAM$) THEN
      '
      ' database does not exist in the current directory,
      ' check the directory where the parameter file was found
      TmpS$ = PathNameBLM$(FName$) + FirstName$(sWHAM$) + "." + ExtName$(sWHAM$)
      'CALL Pmsg(TmpS$)
      IF NOT Exist(TmpS$) THEN
         ErrorMsg$= "File " & sWHAM$ & " not found."
         ErrorCode% = TRUE
         EXIT SUB
      ELSE
         sWHAM$ = TmpS$
      END IF
   END IF
   OPEN sWHAM$ FOR INPUT AS #iWHAMDBS
   INPUT #iWHAMDBS, S$                                     ' database identifier
   INPUT #iWHAMDBS, S$, wNCOOH!(1), wPKHA!(1), wPKHB!(1), wDPKHA!(1), wDPKHB!(1), wP!(1), wFPR!(1), wRADIUS!(1), wMOLWT!(1)' HA properties
   INPUT #iWHAMDBS, S$, wNCOOH!(2), wPKHA!(2), wPKHB!(2), wDPKHA!(2), wDPKHB!(2), wP!(2), wFPR!(2), wRADIUS!(2), wMOLWT!(2)' FA properties
   INPUT #iWHAMDBS, S$, DLF!                                ' double layer overlap factor
   INPUT #iWHAMDBS, S$, wKZED!                               ' DL vol at low Z
   INPUT #iWHAMDBS, S$, wNODATA%                            ' number of species

   '
   '  Add proton binding to humic and fulvic sites
   '
   '  Humic acid A sites
   '
      ii = NSpecies + 1
      jj = iHA1
      kk = 1
      Label(ii) = LEFT$(Label(jj), 3)
      EC(ii) = 10 ^ (-wPKHA!(1) - (wDPKHA!(1) / 2))
      SCa(ii, jj) = ONE
      SCb(ii, jj) = ONE
      SCa(ii, kk) = -ONE
      SCb(ii, kk) = -ONE

      ii = NSpecies + 2
      jj = iHA2
      Label(ii) = LEFT$(Label(jj), 3)
      EC(ii) = 10 ^ (-wPKHA!(1) - (wDPKHA!(1) / 6))
      SCa(ii, jj) = ONE
      SCb(ii, jj) = ONE
      SCa(ii, kk) = -ONE
      SCb(ii, kk) = -ONE

      ii = NSpecies + 3
      jj = iHA3
      Label(ii) = LEFT$(Label(jj), 3)
      EC(ii) = 10 ^ (-wPKHA!(1) + (wDPKHA!(1) / 6))
      SCa(ii, jj) = ONE
      SCb(ii, jj) = ONE
      SCa(ii, kk) = -ONE
      SCb(ii, kk) = -ONE

      ii = NSpecies + 4
      jj = iHA4
      Label(ii) = LEFT$(Label(jj), 3)
      EC(ii) = 10 ^ (-wPKHA!(1) + (wDPKHA!(1) / 2))
      SCa(ii, jj) = ONE
      SCb(ii, jj) = ONE
      SCa(ii, kk) = -ONE
      SCb(ii, kk) = -ONE

   '
   '  Humic acid B sites
   '
      ii = NSpecies + 5
      jj = iHB1
      Label(ii) = LEFT$(Label(jj), 3)
      EC(ii) = 10 ^ (-wPKHB!(1) - (wDPKHB!(1) / 2))
      SCa(ii, jj) = ONE
      SCb(ii, jj) = ONE
      SCa(ii, kk) = -ONE
      SCb(ii, kk) = -ONE

      ii = NSpecies + 6
      jj = iHB2
      Label(ii) = LEFT$(Label(jj), 3)
      EC(ii) = 10 ^ (-wPKHB!(1) - (wDPKHB!(1) / 6))
      SCa(ii, jj) = ONE
      SCb(ii, jj) = ONE
      SCa(ii, kk) = -ONE
      SCb(ii, kk) = -ONE

      ii = NSpecies + 7
      jj = iHB3
      Label(ii) = LEFT$(Label(jj), 3)
      EC(ii) = 10 ^ (-wPKHB!(1) + (wDPKHB!(1) / 6))
      SCa(ii, jj) = ONE
      SCb(ii, jj) = ONE
      SCa(ii, kk) = -ONE
      SCb(ii, kk) = -ONE

      ii = NSpecies + 8
      jj = iHB4
      Label(ii) = LEFT$(Label(jj), 3)
      EC(ii) = 10 ^ (-wPKHB!(1) + (wDPKHB!(1) / 2))
      SCa(ii, jj) = ONE
      SCb(ii, jj) = ONE
      SCa(ii, kk) = -ONE
      SCb(ii, kk) = -ONE

   '
   '  Fulvic acid A sites
   '
      ii = NSpecies + 9
      jj = iFA1
      Label(ii) = LEFT$(Label(jj), 3)
      EC(ii) = 10 ^ (-wPKHA!(2) - (wDPKHA!(2) / 2))
      SCa(ii, jj) = ONE
      SCb(ii, jj) = ONE
      SCa(ii, kk) = -ONE
      SCb(ii, kk) = -ONE

      ii = NSpecies + 10
      jj = iFA2
      Label(ii) = LEFT$(Label(jj), 3)
      EC(ii) = 10 ^ (-wPKHA!(2) - (wDPKHA!(2) / 6))
      SCa(ii, jj) = ONE
      SCb(ii, jj) = ONE
      SCa(ii, kk) = -ONE
      SCb(ii, kk) = -ONE

      ii = NSpecies + 11
      jj = iFA3
      Label(ii) = LEFT$(Label(jj), 3)
      EC(ii) = 10 ^ (-wPKHA!(2) + (wDPKHA!(2) / 6))
      SCa(ii, jj) = ONE
      SCb(ii, jj) = ONE
      SCa(ii, kk) = -ONE
      SCb(ii, kk) = -ONE

      ii = NSpecies + 12
      jj = iFA4
      Label(ii) = LEFT$(Label(jj), 3)
      EC(ii) = 10 ^ (-wPKHA!(2) + (wDPKHA!(2) / 2))
      SCa(ii, jj) = ONE
      SCb(ii, jj) = ONE
      SCa(ii, kk) = -ONE
      SCb(ii, kk) = -ONE

   '
   '  Fulvic acid B sites
   '
      ii = NSpecies + 13
      jj = iFB1
      Label(ii) = LEFT$(Label(jj), 3)
      EC(ii) = 10 ^ (-wPKHB!(2) - (wDPKHB!(2) / 2))
      SCa(ii, jj) = ONE
      SCb(ii, jj) = ONE
      SCa(ii, kk) = -ONE
      SCb(ii, kk) = -ONE

      ii = NSpecies + 14
      jj = iFB2
      Label(ii) = LEFT$(Label(jj), 3)
      EC(ii) = 10 ^ (-wPKHB!(2) - (wDPKHB!(2) / 6))
      SCa(ii, jj) = ONE
      SCb(ii, jj) = ONE
      SCa(ii, kk) = -ONE
      SCb(ii, kk) = -ONE

      ii = NSpecies + 15
      jj = iFB3
      Label(ii) = LEFT$(Label(jj), 3)
      EC(ii) = 10 ^ (-wPKHB!(2) + (wDPKHB!(2) / 6))
      SCa(ii, jj) = ONE
      SCb(ii, jj) = ONE
      SCa(ii, kk) = -ONE
      SCb(ii, kk) = -ONE

      ii = NSpecies + 16
      jj = iFB4
      Label(ii) = LEFT$(Label(jj), 3)
      EC(ii) = 10 ^ (-wPKHB!(2) + (wDPKHB!(2) / 2))
      SCa(ii, jj) = ONE
      SCb(ii, jj) = ONE
      SCa(ii, kk) = -ONE
      SCb(ii, kk) = -ONE

   '
   '  This info.is the same for all species
   '
   FOR ii = NSpecies + 1 TO NSpecies + 16
      SType(ii) = stAqueous
      'ActCorr(ii) = acDebye
      ActCorr(ii) = acWHAM
      TempCorr(ii) = FALSE
      DeltaH(ii) = 0
      ECTemp(ii) = 273
      InitialConc(ii) = 0
      V(ii) = -ONE
      WSpec(ii) = TRUE
   NEXT ii

   '
   '  Adjust array sizes again
   '
   NSpecies = NSpecies + 16

'IF FALSE THEN
   NWHAM = 0
   FOR i% = 1 TO wNODATA%
      'INPUT #iWHAMDBS, wX%, wN$(wX%), wCH%(wX%), wM1%(wX%), wM2%(wX%), wM3%(wX%), wS1%(wX%), wS2%(wX%), wS3%(wX%), wLK!(wX%), wDH!(wX%), wPKMHA!(1, wX%), wPKMHA!(2, wX%)
      INPUT #iWHAMDBS, wtX%, wtN$, wtCH%, wtM1%, wtM2%, wtM3%, wtS1%, wtS2%, wtS3%, wtLK!, wtDH!, wtPKMHA1!, wtPKMHA2!
      IF wtPKMHA1! = 999 THEN
         wtBINDTEST% = FALSE     ' identifies species that don't bind to HS
      ELSE
         wtBINDTEST% = TRUE
      END IF
      wBINDTEST%(wtX%) = wtBINDTEST%
      '
      '  Check to see if all the components needed to make
      '  this species are included in this simulation
      '
      j = FindSpecies%(wtN$)
      IF (j <> 0) THEN
         NWHAM = NWHAM + 1
         'WHAMID%(NWHAM) = wtX%
         WHAMID(NWHAM) = wtX%
      ELSEIF (wtX% = 51) THEN
         '
         ' OH is a special case
         NWHAM = NWHAM + 1
         'WHAMID%(NWHAM) = wtX%
         WHAMID(NWHAM) = wtX%
      END IF

      IF AddSpecies%(wtM1%, wtM2%, wtM3%) THEN

         IF (wtX% > 100) THEN
            '
            '  This is a new species entirely and must be added
            '  even in an inorganic form
            '
            NSpecies = NSpecies + 1
            Label(NSpecies) = wtN$
            EC(NSpecies) = 10 ^ wtLK!
            DeltaH(NSpecies) = wtDH!
            SType(NSpecies) = stAqueous
            ActCorr(NSpecies) = acDebye
            TempCorr(NSpecies) = FALSE
            InitialConc(NSpecies) = 0
            V(NSpecies) = wtCH%
            WSpec(NSpecies) = TRUE
            CMap(wtX%) = NSpecies

            IF (wtM1% > 0) THEN
               jj = FindSpecies%(wNs$(wtM1%))
               SCa(NSpecies, jj) = wtS1%
               SCb(NSpecies, jj) = wtS1%
            END IF
            IF (wtM2% > 0) THEN
               jj = FindSpecies%(wNs$(wtM2%))
               SCa(NSpecies, jj) = wtS2%
               SCb(NSpecies, jj) = wtS2%
            END IF
            IF (wtM3% > 0) THEN
               jj = FindSpecies%(wNs$(wtM3%))
               SCa(NSpecies, jj) = wtS3%
               SCb(NSpecies, jj) = wtS3%
            END IF
         ELSE
            '
            '  This is a component and should already be in CHESS
            '
            CMap(wtX%) = FindSpecies%(wtN$)
         END IF
         IF wtBINDTEST% THEN
            '
            ' We need to add species to describe bound forms
            '
            '
            ' Start with mono-dentate forms
            '
            ii = iHA1 - 1
            jj = NSpecies
            kk = FindSpecies%(wtN$)
            FOR j = 1 TO 16
               NSpecies = NSpecies + 1
               Label(NSpecies) = LEFT$(Label(ii + j), 3) + "-" + wtN$
               IF (j <= 8) THEN
                  EC(NSpecies) = 10 ^ (-wtPKMHA1!)
               ELSE
                  EC(NSpecies) = 10 ^ (-wtPKMHA2!)
               END IF
               SCa(NSpecies, ii + j) = ONE
               SCb(NSpecies, ii + j) = ONE
               SCa(NSpecies, 1) = -ONE
               SCb(NSpecies, 1) = -ONE
               IF (kk < NComp) THEN
                  SCa(NSpecies, kk) = ONE
                  SCb(NSpecies, kk) = ONE
               ELSE
                  IF (wtM1% > 0) THEN
                     k = FindSpecies%(wNs$(wtM1%))
                     SCa(NSpecies, k) = ONE
                     SCb(NSpecies, k) = ONE
                  END IF
                  IF (wtM2% > 0) THEN
                     k = FindSpecies%(wNs$(wtM2%))
                     SCa(NSpecies, k) = ONE
                     SCb(NSpecies, k) = ONE
                  END IF
                  IF (wtM3% > 0) THEN
                     k = FindSpecies%(wNs$(wtM3%))
                     SCa(NSpecies, k) = ONE
                     SCb(NSpecies, k) = ONE
                  END IF
               END IF
               DeltaH(NSpecies) = wtDH!
               SType(NSpecies) = stAqueous
               'ActCorr(NSpecies) = acDebye
               ActCorr(NSpecies) = acWHAM
               TempCorr(NSpecies) = FALSE
               InitialConc(NSpecies) = 0
               V(NSpecies) = wtCH%
            NEXT j
            '
            ' Next, bi-dentate forms
            '
         END IF
      END IF
      '
      '  Add data to WHAM species list
      '
      wX% = wtX%
      wNs$(wX%) = wtN$
      wCH%(wX%) = wtCH%
      wM1%(wX%) = wtM1%
      wM2%(wX%) = wtM2%
      wM3%(wX%) = wtM3%
      wS1%(wX%) = wtS1%
      wS2%(wX%) = wtS2%
      wS3%(wX%) = wtS3%
      wLK!(wX%) = wtLK!
      wDH!(wX%) = wtDH!
      wPKMHA!(1, wX%) = wtPKMHA1!
      wPKMHA!(2, wX%) = wtPKMHA2!
   NEXT i%

   ' Add one to NWHAM to account for DOC
   NWHAM = NWHAM + 1

'END IF
   CLOSE #iWHAMDBS
   CALL WHAMDefine
   PHFIX$ = "YES"
   PCO2% = 999
   PRECISION! = .01

   'CALL PrintI(CMap%(), NSP%)
END IF

IF NPhase > 0 THEN
   CALL SelectComp(SComp(), SPhase(), ErrorCode%, ErrorMsg$)
END IF

CALL SSaturate

CALL SimplifyProblem


END SUB

DEFSNG A-Z
FUNCTION FindSpecies% (Spec$)
'
'  Returns the position of species Spec$
'  If no species is found, returns 0
'
FindSpecies% = 0
i% = 0
DO
   i% = i% + 1
   IF Label$(i%) = Spec$ THEN
      FindSpecies% = i%
      EXIT DO
   END IF
LOOP UNTIL (i% = NSpecies)
END FUNCTION

FUNCTION FindWHAMSpecies% (Spec$)
'
'  Returns the position of species Spec$
'  If no species is found, returns 0
'
FindWHAMSpecies% = 0
i% = 0
DO
   i% = i% + 1
   IF wNs$(i%) = Spec$ THEN
      FindWHAMSpecies% = i%
      EXIT DO
   END IF
LOOP UNTIL (i% = NSP%)

END FUNCTION

DEFINT A-Z
SUB FirstGuess (ErrorCode%, ErrorMsg$, ISTR#)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB FirstGuess                                                       |
'  |                                                                       |
'  |  This subroutine provides seed values (initial guesses) for the       |
'  |  values of component concentrations before the non-linear equation    |
'  |  solver iteration can begin.                                          |
'  |                                                                       |
'  |  If initial values were provided by the calling program (stored in    |
'  |  the OC() array) they will not be over-written.  However, if values   |
'  |  were not provided, a reasonable guess is provided depending on the   |
'  |  component type and reaction list.                                    |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*


' *--------------------------------------------*
' |  Simple Variables Local to SUB FirstGuess  |
' *--------------------------------------------*
DIM c AS INTEGER            '  Array/loop index variable
DIM cc AS INTEGER           '  Array/loop index variable
DIM IStr2 AS DOUBLE         '  First guess for ionic strength


FOR cc = 1 TO NComp
   c = CompOrder(cc)
   SELECT CASE cc
      CASE 1 TO NSolve
         ' *---------------------------------------*
         ' |  These components are not fixed, but  |
         ' |  are used in species reactions        |
         ' *---------------------------------------*
         IF CType(c) = ctSurfPot THEN
            ' *----------------------------------------------*
            ' |  First guess for the 'dummy' component used  |
            ' |  in the EDL model corresponds to a small     |
            ' |  surface potential of 1 mV                   |
            ' *----------------------------------------------*
            IF OC(c) = ZeroD THEN
               OC(c) = PsiToP#(.0000001)
            END IF
         ELSE
            ' All other non-fixed but used components that need a first guess
            IF OC(c) <= ZeroD THEN
               OC(c) = ABS(KT(c) / (1000 * CtoM!(c)))
            END IF

            ' If it's still equal to zero then just set it to something small
            IF OC(c) = ZeroD THEN
               OC(c) = .00001
            END IF
         END IF
      CASE NSolve + 1 TO NUsed
         ' *-----------------------------------------------------*
         ' |  These are fixed activity components.  If initial   |
         ' |  guesses haven't been specified then see if totals  |
         ' |  were specified.  If not, then there is a problem.  |
         ' *-----------------------------------------------------*
         IF OC(c) = ZeroD THEN
            IF KT(c) = ZeroD THEN
'              ErrorMsg$ = "A zero concentration was specified for a fixed activity component."
               ErrorMsg$ = "A zero concentration was specified for a fixed activity component ("
               ErrorMsg$ = ErrorMsg$ + Label(c) + ")."
               ErrorCode = TRUE
               EXIT SUB
            END IF
            OC(c) = KT(c) / CtoM!(c)
         END IF
      CASE NUsed TO NComp
         ' *--------------------------------------------------------------*
         ' |  These components are not used, so set equal to total conc.  |
         ' *--------------------------------------------------------------*
         OC(c) = KT(c) * SiteDen(c) / CtoM!(c)
   END SELECT

   ' *---------------------------------------*
   ' |  A first estimate for ionic strength  |
   ' *---------------------------------------*
   IF SType(c) = stAqueous THEN
      IStr2# = IStr2# + (KT(c) * SiteDen(c) / CtoM!(c)) * V(c) ^ 2
   END IF

   IF (iDebug = 01 OR iDebug = 03) THEN
      CALL DebugLog ("1stG", Label(cc) + STR$(OC(cc)))
   END IF
NEXT cc

' *----------------------------------------*
' |  Use the first guess for IStr only if  |
' |  there is currently no value.          |
' *----------------------------------------*
IF ISTR# = ZeroD THEN
   ISTR# = IStr2# * .5
END IF

iDoneWHAM = FALSE

IF (iDebug = 01 OR iDebug = 03) THEN
    FOR cc = 1 TO NComp
       CALL DebugLog ("1stG", Label(cc) + STR$(OC(cc)))
    NEXT cc
END IF

END SUB

FUNCTION IAP# (P AS INTEGER)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  FUNCTION IAP                                                         |
'  |                                                                       |
'  |  Calculates the ion activity product for phase p                      |
'  |                                                                       |
'  |  Function Parameters:                                                 |
'  |                                                                       |
'  |     p          (integer, input)   Array index of a phase              |
'  |                                                                       |
'  |     IAP        (double, output)   Ion activity product of phase P     |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

' *-----------------------------------------------*
' |  Simple Variables Local to SUB DefineProblem  |
' *-----------------------------------------------*
DIM Activity AS DOUBLE      '  Chemical activity for a given species
DIM c AS INTEGER            '  Array/loop index variable
DIM Temp AS DOUBLE          '  Temporary variable

Temp# = 1
FOR c = 1 TO NComp
   Activity# = OC(c) * ActCoef(c)
   Temp# = Temp# * (Activity# ^ SP(P, c))
NEXT c

IAP# = Temp#

END FUNCTION

DEFSNG A-H, O-Z
SUB IndSortI (iArray(), Index(), NEls, iDir)
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  | SUB IndSortI                                                            |
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
'  |    iArray()  - integer -   Array dimensioned from at least 1 to NEls  |
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

      IF iArray(Index(j)) > iArray(Index(k)) THEN
         SWAP Index(j), Index(k)
         InOrder = FALSE
      END IF
   NEXT i
LOOP UNTIL InOrder

END SUB


DEFINT A-H, O-Z
SUB JacobianAA
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB JacobianAA                                                       |
'  |                                                                       |
'  |  Calculates an analytical approximation to the Jacobian matrix.       |
'  |                                                                       |
'  |  This method is approximate because only the reaction stoichiometry   |
'  |  is considered.  All extra-thermodynamic calculations are ignored     |
'  |  (for example, activity coefficients).                                |
'  |                                                                       |
'  |  This analytical approximation is less accurate but much faster than  |
'  |  a numerical approximation to the exact residual equations and is     |
'  |  good enough for most problems.                                       |
'  |                                                                       |
'  |  The partial derivitives of all residuals with respect                |
'  |  to all components are stored in the Jacobian matrix Z().             |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

' *--------------------------------------------*
' |  Simple Variables Local to SUB JacobianAA  |
' *--------------------------------------------*
DIM Actv AS DOUBLE          '  Activity of a given species
DIM j AS INTEGER            '  Array/loop index variable
DIM jj AS INTEGER           '  Array/loop index variable
DIM k AS INTEGER            '  Array/loop index variable
DIM kk AS INTEGER           '  Array/loop index variable
DIM Sum AS DOUBLE           '  Array/loop index variable
DIM S  AS INTEGER           '  Array/loop index variable
DIM ss AS INTEGER           '  Array/loop index variable

FOR jj = 1 TO NSolve
   j = CompOrder(jj)
   IF (CType(j) = ctChargeBal) THEN
      FOR kk = 1 TO NSolve
         k = CompOrder(kk)
         Sum# = 0
         FOR S = 1 TO NSpecies
            IF (SType(S) = stAqueous) OR (SType(S) = stPrecipitate) THEN
               Sum# = Sum# + OC(S) * tSCb(S, k) * V(S)
            END IF
         NEXT S
         Z(jj, kk) = Sum# / OC(k)
      NEXT kk
   ELSE
      FOR kk = 1 TO NSolve
         k = CompOrder(kk)
         Sum# = 0
         FOR ss = 1 TO SpecList(jj, 0)
            S = SpecList(jj, ss)
            Sum# = Sum# + tSCa(S, k) * tSCb(S, j) * OC(S) * CtoM!(S)
         NEXT ss
         IF OC(k) <> 0 THEN
            Z(jj, kk) = Sum# / OC(k)
         END IF
      NEXT kk
   END IF
NEXT jj

END SUB

SUB JacobianNum
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB JacobianNum                                                      |
'  |                                                                       |
'  |  Calculates a numerical Jacobian matrix.                              |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

' *--------------------------------------------*
' |  Simple Variables Local to SUB JacobianAA  |
' *--------------------------------------------*
DIM Actv AS DOUBLE          '  Activity of a given species
DIM j AS INTEGER            '  Array/loop index variable
DIM jj AS INTEGER           '  Array/loop index variable
DIM k AS INTEGER            '  Array/loop index variable
DIM kk AS INTEGER           '  Array/loop index variable
DIM Sum AS DOUBLE           '  Array/loop index variable
DIM S  AS INTEGER           '  Array/loop index variable
DIM ss AS INTEGER           '  Array/loop index variable
DIM SaveR(1 TO NComp) AS DOUBLE
DIM SaveOC(1 TO NComp) AS DOUBLE
DIM Znum AS DOUBLE
DIM Zden AS DOUBLE
DIM ZnumMax AS DOUBLE
DIM ZnumMin AS DOUBLE
DIM OCScale(1 TO NComp) AS STATIC DOUBLE

CALL CalcSpeciesConc(ISTR#, ErrorCode%, ErrorMsg$)
CALL CalcResidual(MaxError!, MaxComp%)
FOR jj = 1 TO NComp
  'j = CompOrder(jj)
  SaveR(jj) = R(jj)
  SaveOC(jj) = OC(jj)
NEXT jj

IF OCScale(1) = 0 THEN
    FOR j=1 TO NComp
        OCScale(j) = 1.001
    NEXT j
END IF

FOR jj = 1 TO NSolve
  j = CompOrder(jj)

  OC(j) = OC(j) * OCScale(j)
  Zden = (OC(j) - SaveOC(j))
  CALL CalcSpeciesConc(ISTR#, ErrorCode%, ErrorMsg$)
  CALL CalcResidual(MaxError!, MaxComp%)

  FOR kk = 1 TO NSolve
     k = CompOrder(kk)
     Znum = R(kk) - SaveR(kk)
     Zn(kk,jj) = Znum/Zden
     IF (kk=1) THEN
         ZnumMax = Znum
         ZnumMin = Znum
     ELSE
         IF (Znum>ZnumMax) THEN ZnumMax = Znum
         IF (Znum<ZnumMax) THEN ZnumMin = Znum
     END IF
  NEXT kk
'  PRINT Label(j);" ";
'  PRINT FORMAT$(ZnumMin, "0.0E+##");" ";FORMAT$(ZnumMax, "0.0E+##")
  OC(j) = SaveOC(j)
  'OCScale(j) = OCScale(j) * (100000/ZnumMax)
  IF (ZnumMax > 1.0E-4) THEN
      OCScale(j) = OCScale(j)/1000
  ELSE
      OCScale(j) = OCScale(j)*10
  END IF
  'OCScale(j) = 1.001
NEXT jj

CALL CalcSpeciesConc(ISTR#, ErrorCode%, ErrorMsg$)
CALL CalcResidual(MaxError!, MaxComp%)

END SUB


DEFSNG A-Z
FUNCTION Log10D# (X#)

IF (X# = ZeroD) THEN X = ONE

Log10D# = LOG(X#) / LOG(10#)

END FUNCTION

FUNCTION LOG10S! (X!)

IF (X = ZeroS) THEN X = ONE

LOG10S! = LOG(X) / LOG(10)

END FUNCTION

DEFINT A-Z
FUNCTION PsiToP# (Psi#)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  FUNCTION PsiToP                                                      |
'  |                                                                       |
'  |  Calculate the EDL component 'P' for a given surface potential 'Psi'  |
'  |  where the value of Psi is in Volts.                                  |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  |                                                                       |
'  *-----------------------------------------------------------------------*


' *-----------------------------------------------*
' |  Simple Variables Local to SUB DefineProblem  |
' *-----------------------------------------------*
DIM T AS DOUBLE             '  Temporary single precision variable

T# = Fcon# * Psi# / (Rcon# * SysTemp)
PsiToP# = EXP(-T#)

END FUNCTION

SUB SSaturate()
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB SSaturate                                                        |
'  |                                                                       |
'  |  Transform stoichiometric coefficients and chemical equilibiurm       |
'  |  constants according to substitutions of selected components to       |
'  |  reflect saturation of solid or gaseous phases.                       |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

' *-----------------------------------------------*
' |  Simple Variables Local to SUB DefineProblem  |
' *-----------------------------------------------*
DIM A AS INTEGER            '  Array/loop index variable
DIM c AS INTEGER            '  Array/loop index variable
DIM d AS INTEGER            '  Array/loop index variable
DIM e AS INTEGER            '  Array/loop index variable
DIM FirstPhase AS INTEGER   '  Boolean flag
DIM j AS INTEGER            '  Array/loop index variable
DIM k AS INTEGER            '  Array/loop index variable
DIM LogMsg AS STRING        '  Variable used to write a msg to a log file
DIM P AS INTEGER            '  Array/loop index variable
DIM q AS INTEGER            '  Array/loop index variable
DIM S AS INTEGER            '  Array/loop index variable
DIM Temp AS SINGLE          '  Temporary variable

' *--------------------------------------------------*
' |  Set up proper component type for active phases  |
' *--------------------------------------------------*
FOR P = 1 TO NPhase
   IF PType(P) = ptSaturated THEN
      CType(SComp(P)) = ctSubstituted
   END IF
NEXT P

' *-----------------------------------------------------------------------*
' |  Transform saturated phase information by forward elimination scheme  |
' *-----------------------------------------------------------------------*
FirstPhase = FALSE
FOR d = 1 TO NPhase
   e = PSubOrder(d)
   IF PType(e) = ptSaturated THEN
      ' *----------------------------------------------------------------*
      ' |  Take the first ptSaturated phase in PSubOrder and keep as is  |
      ' *----------------------------------------------------------------*
      IF NOT FirstPhase THEN
         FOR c = 1 TO NComp
            tSP(e, c) = SP(e, c)
         NEXT c
         tSEC(e) = SEC(e)
         FirstPhase = TRUE
      END IF

      ' *----------------------------------------------*
      ' |  Transform all the other ptSaturated phases  |
      ' *----------------------------------------------*
      FOR q = d + 1 TO NPhase
         P = PSubOrder(q)
         IF PType(P) = ptSaturated THEN
            A = SComp(e)
            Temp! = -SP(P, A) / SP(e, A)
            FOR c = 1 TO NComp
               IF c = A THEN
                  tSP(P, c) = 0
               ELSE
                  tSP(P, c) = SP(P, c) + SP(e, c) * Temp!
               END IF
            NEXT c
            tSEC(P) = SEC(P) * SEC(e) ^ Temp!
         END IF
      NEXT q
   END IF
NEXT d

' *------------------------------------------------*
' |  Now use transformed phase info. to transform  |
' |  species stoichiometric coefficients           |
' *------------------------------------------------*

' work on copies to avoid destroying the original info.
FOR S = 1 TO NSpecies
   FOR c = 1 TO NComp
      ' Copy coefficients
      tSCa(S, c) = SCa(S, c)
      tSCb(S, c) = SCb(S, c)
   NEXT c
   tEC(S) = hEC(S)
NEXT S

FOR q = 1 TO NPhase
   P = PSubOrder(q)
   IF PType(P) = ptSaturated THEN
      A = SComp(P)
      FOR S = NFirstS TO NSpecies

         FOR c = 1 TO NComp
            ' *----------------------------------*
            ' |  Transform species coefficients  |
            ' *----------------------------------*
            tSCa(S, c) = tSCa(S, c) - ((tSCa(S, A) * tSP(P, c)) / tSP(P, A))
            tSCb(S, c) = tSCb(S, c) - ((tSCb(S, A) * tSP(P, c)) / tSP(P, A))

            IF ABS(tSCa(S, c)) < RoundOffLimit AND ABS(tSCa(S, c)) <> 0 THEN
               ' *------------------------------------------------------*
               ' |  WRITE TO LOG FILE WHEN THIS CONDITION IS TRUE!!!    |
               ' |  It may indicate an error associated with round-off  |
               ' |  due to the single precision nature of the SC array  |
               ' *------------------------------------------------------*
               LogMsg$ = "A transformed coefficient of " + STR$(tSCa(S, c))
               LogMsg$ = LogMsg$ + " was calculated for species " + Label(S)
               LogMsg$ = LogMsg$ + "and component " + Label(c) + "and set to zero."

               tSCa(S, c) = ZeroS
               tSCb(S, c) = ZeroS
            END IF
         NEXT c

         ' Transform equilibrium constants
         tEC(S) = hEC(S) * (tSEC(P) ^ (SCa(S, A) / tSP(P, A)))

      NEXT S
   END IF
NEXT q

' *------------------------------------------------*
' |  Insert revised stoichiometry for substituted  |
' |  components into tSC() array                   |
' *------------------------------------------------*
FOR j = 1 TO NComp
   IF CType(j) = 3 THEN
      tSCa(j, j) = 0
      tSCb(j, j) = 0

      P = SPhase(j)
      FOR k = 1 TO NComp
         ' *---------------------------------------------------*
         ' |  Substituted components are changed according to  |
         ' |  transformed solubility expression                |
         ' *---------------------------------------------------*
         IF k <> j THEN
            tSCa(j, k) = -tSP(P, k)
            tSCb(j, k) = -tSP(P, k)
         END IF
      NEXT k
      tEC(j) = tSEC(P) ^ (1 / tSP(P, j))
   ELSE
      tSCa(j, j) = SCa(j, j)
      tSCb(j, j) = SCb(j, j)
   END IF
NEXT j

' *------------------------------------------------------------------*
' |  Now that the transformed coefficients have changed, a new call  |
' |  to SUB SimplifyProblem is required to set up new species and    |
' |  component lists.                                                |
' *------------------------------------------------------------------*
CALL SimplifyProblem

END SUB

SUB SelectComp (SComp() AS INTEGER, SPhase() AS INTEGER, ErrorCode%, ErrorMsg$)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB SelectComp                                                       |
'  |                                                                       |
'  |  Select a component that will substituted out of the residual         |
'  |  equations for each phase                                             |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

' *--------------------------------------------*
' |  Simple Variables Local to SUB SelectComp  |
' *--------------------------------------------*
DIM c AS INTEGER            '  Array/loop index variable
DIM Direction AS INTEGER    '  Constant used to determine sort direction
DIM Num AS INTEGER          '  Counting variable
DIM MinSP AS INTEGER        '  Temporary variable
DIM MinC AS INTEGER         '  Temporary variable
DIM P AS INTEGER            '  Array/loop index variable
DIM q AS INTEGER            '  Array/loop index variable
DIM Start AS INTEGER        '  Counting variable

' *-------------------------------------------*
' |  Array Variables Local to SUB SelectComp  |
' *-------------------------------------------*
DIM Index(NPhase)               ' Table of sort info
DIM AvailComp(NComp) AS INTEGER ' Boolean array of component availability
DIM NumComp(NPhase) AS INTEGER  ' Number of comp w/non-zero coef

' *-------------------------------------------*
' |  Make a list of components available for  |
' |  substitution (no fixed components)       |
' *-------------------------------------------*
FOR c = 1 TO NComp
   IF CType(c) = ctFloat THEN
      AvailComp(c) = TRUE
   ELSE
      AvailComp(c) = FALSE
   END IF
NEXT c

FOR P = 1 TO NPhase
   Index(P) = P
NEXT P

Start = 0
DO
   ' *--------------------------------------------------------------*
   ' |  For each phase, make a list of components with non-zero     |
   ' |  coefficients AND that are still available for substitution  |
   ' *--------------------------------------------------------------*
   Start = Start + 1
   FOR q = Start TO NPhase
      P = Index(q)
      Num = 0
      FOR c = 1 TO NComp
         IF (SP(P, c) <> 0) AND (AvailComp(c) = TRUE) THEN
            Num = Num + 1
         END IF
      NEXT c

      IF Num = 0 THEN
         ErrorMsg$ = "Gibbs phase rule violation in list of phases."
         ErrorCode = TRUE
      ELSE
         NumComp(P) = Num
      END IF
   NEXT q

   ' *-------------------------------------------------------------*
   ' |  Now sort that list of the no. of components in each phase  |
   ' *-------------------------------------------------------------*
   Direction = 0                          ' Ascending

   CALL IndSortI(NumComp(), Index(), NPhase, Direction)

   ' *----------------------------------------*
   ' |  Select a component to substitute for  |
   ' |  the phase at the top of the list      |
   ' *----------------------------------------*
   MinSP = 0
   MinC = 0
   P = Index(Start)
   FOR c = 1 TO NComp
      IF ABS(SP(P, c)) <> 0 AND AvailComp(c) THEN
         IF MinSP = 0 OR ABS(SP(P, c)) <= MinSP THEN
            MinSP = ABS(SP(P, c))
            MinC = c
         END IF
      END IF
   NEXT c
   IF MinC <> 0 THEN
      ' *-----------------------------------------*
      ' |  Keep track of the component that       |
      ' |  this phase will substitute for . . .   |
      ' *-----------------------------------------*
      SPhase(MinC) = Index(Start)

      ' *-----------------------------------------*
      ' |  . . . and keep track of the phase      |
      ' | that can be substituted for this        |
      ' | component.                              |
      ' *-----------------------------------------*
      SComp(Index(Start)) = MinC

      IF PType(Start) = ptSaturated THEN
         CType(MinC) = ctSubstituted
      END IF

      ' *---------------------------------------*
      ' |  Take this component off the list of  |
      ' |  available comonents.                 |
      ' *---------------------------------------*
      AvailComp(MinC) = FALSE
      NumComp(Index(Start)) = -1

      ' *------------------------------------------------*
      ' |  Record the order that substitutions are made  |
      ' *------------------------------------------------*
      PhaseOrder(Start) = Index(Start)
   END IF
LOOP UNTIL Start = NPhase

' *----------------------------------*
' |  Now determine real phase order  |
' *----------------------------------*
SetPhaseOrder PhaseOrder(), PSubOrder()

' *-------------------------------------------*
' |  Free up memory assigned to local arrays  |
' *-------------------------------------------*
ERASE AvailComp, NumComp, Index
END SUB

SUB SetCtoM (AqueousCtoM!, SurfaceCtoM!)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB SetCtoM                                                          |
'  |                                                                       |
'  |  Stores the values of concentration to mass conversion for all        |
'  |  aqueous and surface species.                                         |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

' *-----------------------------------------------*
' |  Simple Variables Local to SUB ChargeBalance  |
' *-----------------------------------------------*
DIM S AS INTEGER            '  Array/loop index variable

FOR S = 1 TO NSpecies
   IF SType(S) = stAqueous THEN
      CtoM(S) = AqueousCtoM!
   ELSEIF SType(S) = stSurface THEN
      CtoM(S) = SurfaceCtoM!
   ELSEIF STYpe(S) = stPrecipitate THEN
      CtoM(S) = AqueousCtoM!
   END IF
NEXT S
END SUB

SUB SetCtoM310 ()
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB SetCtoM310                                                       |
'  |                                                                       |
'  |  Stores the values of concentration to mass conversion for all        |
'  |  aqueous and surface species.  Now with variable mass compartments.   |
'  |                                                                       |
'  |  Developed for CHESS ver 3.10, July 2007                              |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  HydroQual, Inc.                                                      |
'  |  rsantore@hydroqual.com                                               |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

' *-----------------------------------------------*
' |  Simple Variables Local to SUB ChargeBalance  |
' *-----------------------------------------------*
DIM S AS INTEGER            '  Array/loop index variable
DIM i AS INTEGER

FOR S = 1 TO NSpecies
   i = FIX(SType(S)/10) + 1
   i = SType(S) + 1
   CtoM(S) = MassVal(i)
NEXT S
END SUB

SUB SetPhaseOrder (PhaseOrder() AS INTEGER, PSubOrder() AS INTEGER)
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB SetPhaseOrder                                                    |
'  |                                                                       |
'  |  The order of operations involving phases are determined in this      |
'  |  subroutine.  Two orders are determined, the order for substitution   |
'  |  of phase solubility relationships for component activities are stored|
'  |  in PSubOrder() and the order for mass balance calculations on phases |
'  |  is stored in PhaseOrder().                                           |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

' *-----------------------------------------------*
' |  Simple Variables Local to SUB ChargeBalance  |
' *-----------------------------------------------*
DIM c AS INTEGER            '  Array/Loop index variable
DIM CompThere AS INTEGER
DIM ErrorCode AS INTEGER    '  Error flag
DIM ErrorMsg AS STRING      '  Error text
DIM GotOne AS INTEGER       '  Boolean flag
DIM NumberDone AS INTEGER   '  Array/Loop index variable
DIM Order AS INTEGER        '  Array/Loop index variable
DIM P AS INTEGER            '  Array/Loop index variable
DIM p2 AS INTEGER           '  Array/Loop index variable
DIM RefP AS INTEGER         '  Array/Loop index variable
DIM TestC AS INTEGER        '  Index variable
DIM TestP AS INTEGER        '  Index variable

' *----------------------------------------------*
' |  Array Variables Local to SUB ChargeBalance  |
' *----------------------------------------------*
DIM MarkPhase(MaxPhases) AS INTEGER

' *-----------------------------------------------------*
' |  Check that there are no two phases that each have  |
' |  the others substituted components                  |
' *-----------------------------------------------------*
ErrorCode = FALSE
FOR P = 1 TO NPhase
   TestC = SComp(P)
   FOR p2 = P + 1 TO NPhase
      IF SP(p2, TestC) <> 0 AND SP(P, SComp(p2)) <> 0 THEN
         ErrorCode = TRUE
      END IF
   NEXT p2
   IF ErrorCode THEN
      ErrorMsg$ = "The list of phases in the parameter file is causing a problem."
      ErrorMsg$ = ErrorMsg$ + "Although this may be a perfectly reasonable combination of phases,"
      ErrorMsg$ = ErrorMsg$ + "the algorithm used in CHESS cannot handle this combination."
   END IF
NEXT P

' Determine order for interdependant phase transformations
'     and place in PhaseOrder()
'
' Strategy: 1) Take a phase, note its associated component.
'           2) Examine the remaining list of phases to see if any
'                 other phase has a non-zero coef. for that comp.
'           3) If no, then mark the selected phase and the order
'                 it was marked.  Remove the phase from the list.
'           4) If yes, then skip this phase.
'           5) If all phases are marked then finish, otherwise
'                 get the next phase in list and begin at step (1)

RefP = 0
Order = 1
NumberDone = 0
DO
   GOSUB NextRefP

   ' Find the component for this phase
   TestC = SComp(RefP)

   ' Check all the other phases to see if they include this component
   CompThere = FALSE
   TestP = 1
   DO
      IF NOT MarkPhase(TestP) AND SP(TestP, TestC) <> 0 AND TestP <> RefP THEN CompThere = TRUE
      TestP = TestP + 1
   LOOP UNTIL TestP > NPhase

   IF NOT CompThere THEN
      ' Ok, take this phase off the list now and mark it's order
      MarkPhase(RefP) = TRUE
      PhaseOrder(Order) = RefP
      Order = Order + 1
      NumberDone = NumberDone + 1
   END IF
LOOP UNTIL NumberDone = NPhase

' Determine order for component substitions
'     and place in PSubOrder()
'
' Strategy: 1) Take a phase and examine its stoichiometric coefficients
'                 to find every component (other than it's associated
'                 component) that has a non-zero component.
'           2) For every component found, determine if this component
'                 is associated with another phase in list.
'           3) If no component for this phase is associated with another
'                 phase, then mark this phase as done and note it's order.
'                 Remove this phase from the list.
'           4) Otherwise, skip this phase for later.
'           5) If all phases are marked, then finish, otherwise get the
'                 next phase in the list and continue at (1)

' *---------------------------------------*
' |  Re-Zero the contents of MarkPhase()  |
' *---------------------------------------*
REDIM MarkPhase(MaxPhases) AS INTEGER

RefP = 0
Order = 1
NumberDone = 0
DO
   GOSUB NextRefP

   ' Check all the components in this phase with non-zero coef.
   CompThere = FALSE
   FOR c = 1 TO NComp
      IF SP(RefP, c) <> 0 THEN
         IF c <> SComp(RefP) THEN
             IF SPhase(c) <> 0 THEN
               IF NOT MarkPhase(SPhase(c)) THEN CompThere = TRUE
             END IF
         END IF
      END IF

   NEXT c

   IF NOT CompThere THEN
      ' Ok, take this phase off the list now and mark it's order
      MarkPhase(RefP) = TRUE
      PSubOrder(Order) = RefP
      Order = Order + 1
      NumberDone = NumberDone + 1
   END IF
LOOP UNTIL NumberDone = NPhase

' *-------------------------------------------*
' |  Free up memory assigned to local arrays  |
' *-------------------------------------------*
ERASE MarkPhase

EXIT SUB    ' END OF SUBPROGRAM ------------------------------------------

NextRefP:
   ' Find the next phase that hasn't been marked as done
   P = RefP + 1
   GotOne = FALSE
   DO
      IF P > NPhase THEN
         ' potential infinite loop here, but never a problem unless
         ' the above check for inappropriate phases was somehow avoided
         P = 1
      ELSEIF MarkPhase(P) THEN
         P = P + 1
      ELSEIF NOT MarkPhase(P) THEN
         GotOne = TRUE
      END IF
   LOOP UNTIL GotOne
   RefP = P
RETURN

END SUB

DEFSNG A-Z
FUNCTION Sign% (X#)
   IF (X# < 0) THEN
      Sign% = -1
   ELSE
      Sign% = 1
   END IF
END FUNCTION

DEFINT A-Z
SUB SimplifyProblem
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB SimplifyProblem                                                  |
'  |                                                                       |
'  |  The purpose of this subroutine is to simplify, as much as possible,  |
'  |  the development of the non-linear equations used to solve the        |
'  |  equilibrium problem.  In order to accomplish this task, components   |
'  |  that are not used to calculate species concentrations as well as     |
'  |  fixed acivity components are identified so they can be ignored,      |
'  |  while the root finding algorithm concentrations only on the          |
'  |  remaining components.  The information about active versus inactive  |
'  |  components and the species that depend on them is stored in the      |
'  |  CompList() and SpecList() arrays.                                    |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*


' *------------------------------------------------*
' |  Simple Variables Local to SUB SimplfyProblem  |
' *------------------------------------------------*
DIM c AS INTEGER            '  Array/loop index variable
DIM cc AS INTEGER           '  Array/loop index variable
DIM Count1 AS INTEGER       '  Temporary variable
DIM Count2 AS INTEGER       '  Temporary variable
DIM Used AS INTEGER         '  Boolean flag
DIM S AS INTEGER            '  Array/loop index variable

' *-----------------------------------------------*
' |  Array Variables Local to SUB SimplfyProblem  |
' *-----------------------------------------------*
DIM CpTmp(NComp) AS INTEGER  ' used just for temporary storage

' *----------------------------------------------------*
' |  Three types of components to be identified here:  |
' *----------------------------------------------------*
NSolve = 0     ' These are components needed to solve equil. problem
NFixed = 0     ' These are fixed activity componets
NotUsed = 0    ' the rest are neither fixed activity, nor used for species

' *--------------------------------------------------------------------*
' |  These three types of components are identified and thier numbers  |
' |  are stored in CompOrder() in the order NSolve, NFixed, NotUsed.   |
' *--------------------------------------------------------------------*
FOR c = 1 TO NComp
   'IF (CType(c) = ctFixed) OR (CType(c) = ctChargeBal) THEN
   IF (CType(c) = ctFixed) THEN
      NFixed = NFixed + 1
      CpTmp(NFixed) = c
   ELSE
      Used = FALSE
      FOR S = NFirstS TO NSpecies
         IF (tSCa(S, c) <> 0) OR (tSCb(S, c) <> 0) THEN Used = TRUE
      NEXT S
      IF Used THEN
         NSolve = NSolve + 1
         CompOrder(NSolve) = c
      ELSE
         NotUsed = NotUsed + 1
         CompOrder(NComp - NotUsed + 1) = c
      END IF
   END IF
NEXT c
' *--------------------------------------------------------------------------*
' |  Now CompOrder() has the list of components from 1 to NSolve that        |
' |  that need to be considered for the problem, and the list of components  |
' |  from NComp backward to NComp-NotUsed are those not used.  Now, just     |
' |  fill in the list of fixed activity components from CpTmp() into the     |
' |  middle of the CompOrder() array.                                        |
' *--------------------------------------------------------------------------*
FOR c = 1 TO NFixed
   CompOrder(NSolve + c) = CpTmp(c)
NEXT c

' *------------------------------------------*
' |  Determine the list of species relevant  |
' |  to the mole balance of components       |
' *------------------------------------------*
NUsed = NSolve + NFixed
REDIM CompList(0 TO NSpecies, 0 TO NUsed) AS INTEGER
REDIM SpecList(0 TO NUsed, 0 TO NSpecies)  AS INTEGER

FOR cc = 1 TO NUsed
   c = CompOrder(cc)
   Count1 = 0
   FOR S = 1 TO NSpecies
      IF tSCb(S, c) <> 0 THEN
         ' *-------------------------------------------------*
         ' |  Store a list of all species relevant to the    |
         ' |  mole balance of each component in SpeciesList  |
         ' *-------------------------------------------------*
         Count1 = Count1 + 1
         SpecList(cc, Count1) = S
      END IF
      IF tSCa(S, c) <> 0 THEN
         ' *---------------------------------------------------------*
         ' |  Store a list of all components relevant to the         |
         ' |  mass action expression for each species in CompList()  |
         ' *---------------------------------------------------------*
         Count2 = CompList(S, 0) + 1
         CompList(S, Count2) = c
         CompList(S, 0) = Count2
      END IF
   NEXT S
   SpecList(cc, 0) = Count1
NEXT cc

' *-------------------------------------------*
' |  Free up memory assigned to local arrays  |
' *-------------------------------------------*
ERASE CpTmp

END SUB

DEFSNG A-Z
SUB SubstituteComp
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB SubstitueComp                                                    |
'  |                                                                       |
'  |  Calculate the equilibrium concentration of components                |
'  |  substituted by phases.                                               |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

' *-----------------------------------------------*
' |  Simple Variables Local to SUB DefineProblem  |
' *-----------------------------------------------*
DIM A AS INTEGER            '  Array/loop index variable
DIM c AS INTEGER            '  Array/loop index variable
DIM q AS INTEGER            '  Array/loop index variable
DIM P AS INTEGER            '  Array/loop index variable
DIM Temp AS DOUBLE          '  Temporary variable

FOR q = 1 TO NPhase
   P = PSubOrder(q)
   IF PType(P) = 1 THEN
      A = SComp(P)
      Temp# = 1
      FOR c = 1 TO NComp
         IF c <> A THEN
            Temp# = Temp# * (OC(c) * ActCoef(c)) ^ (-SP(P, c) / SP(P, A))
         END IF
      NEXT c

      OC(A) = (SEC(P) ^ (1 / SP(P, A))) * Temp#
   END IF
NEXT q

END SUB

DEFINT A-Z
SUB SUpdate (X() AS DOUBLE, Iter%, MaxError!, Skip%)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB SUpdate                                                          |
'  |                                                                       |
'  |  The solution to the set of simutaneous equations derived from the    |
'  |  Jacobian matrix is used to update the free component concentrations. |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*


' *-----------------------------------------------*
' |  Static Variables Local to SUB DefineProblem  |
' *-----------------------------------------------*
STATIC BadCount      AS INTEGER  ' Count of iterations without improvement
STATIC DontBeHasty   AS INTEGER  ' Delay before reducing relaxation
STATIC GoodCount     AS INTEGER  ' Count of iterations with improvement
STATIC NothingCount  AS INTEGER  ' Count of iter. with little change
STATIC OldMaxError   AS SINGLE   ' Storage of MaxError from prev. iter
'STATIC oldOC()       AS DOUBLE   ' Store values of OC() to reset if nec.
STATIC Relax         AS SINGLE   ' Fraction of full NR step being used
STATIC RelaxOn       AS INTEGER  ' Boolean flag sets relaxation on/off
STATIC SkipCount     AS INTEGER  ' Count of disgarded iterations

' *-----------------------------------------------*
' |  Simple Variables Local to SUB DefineProblem  |
' *-----------------------------------------------*
DIM c AS INTEGER            '  Array/loop index variable
DIM cc AS INTEGER           '  Array/loop index variable
DIM Temp1 AS DOUBLE         '  Temporary single precision variable
DIM Temp2 AS DOUBLE         '  Temporary single precision variable

IF Iter% = 1 THEN
   'REDIM oldOC(0 TO NComp) AS DOUBLE
   FOR c = 1 TO NComp
      oldOC(c) = 0
   NEXT c
   oldOC(0) = MaxError!
END IF

IF Iter% > 20 THEN RelaxOn = TRUE

IF Skip THEN
   Skip = 0
   EXIT SUB
END IF

IF RelaxOn THEN
   ' This will provide a relaxed alternative to the Newton-Raphson
   '     iteration.  Sometimes the N-R routine will overcompensate
   '     when calculating adjusted component concentrations resulting
   '     in oscillations and preventing convergence.
   IF Iter% = 1 THEN
      IF Relax! < .05 THEN Relax! = 1
      DontBeHasty = 0
   ELSEIF MaxError! > 1.3 * OldMaxError! OR NothingCount > 10 THEN
      GoodCount = 0
      NothingCount = 0
      BadCount = BadCount + 1
      IF BadCount > 5 THEN
         SkipCount = SkipCount + 1
         IF SkipCount > 5 THEN BadCount = 1
      ELSE
         Relax! = Relax! / 1.5
         IF Relax! = 0! THEN Relax! = .05
         Skip = TRUE
      END IF
   ELSEIF MaxError! < OldMaxError! THEN
      BadCount = 0
      GoodCount = GoodCount + 1
      IF (GoodCount > (3 + Iter% / 3)) AND (GoodCount > DontBeHasty) AND (Relax! < 1!) THEN
         NothingCount = 0
         'DontBeHasty = GoodCount + GoodCount / 5
         DontBeHasty = GoodCount / 5
         Relax! = Relax! * 1.5
         IF Relax! > 1! THEN Relax! = 1!
      END IF
   ELSE
      NothingCount = NothingCount + 1
   END IF
ELSE
   Relax! = 1!
END IF

FOR cc = 1 TO NSolve
   c = CompOrder(cc)
   Temp2# = OC(c)
   IF NOT Skip THEN
      Temp1# = X(cc) * Relax!'/ CtoM!(c)
      oldOC(c) = Temp2#
      IF Temp1# >= Temp2# THEN
         OC(c) = Temp2# / 10!
      ELSE
         OC(c) = OC(c) - Temp1#
      END IF
   ELSE
      Temp1# = X(cc) * Relax! '/ CtoM!(c)
      Temp2# = oldOC(c)
      IF Temp1# >= Temp2# THEN
         OC(c) = Temp2# / 10!
      ELSE
         OC(c) = oldOC(c) - Temp1#
      END IF
   END IF
   IF OC(c) <= ZeroD THEN
      OC(c) = 1E-10
   END IF
NEXT cc

OldMaxError! = MaxError!
IF MaxError! < oldOC(0) THEN oldOC(0) = MaxError!

IF (iDebug = 01 OR iDebug = 05) THEN
'    FOR s = NFirstS TO NSpecies
    FOR c = 1 TO NComp
       CALL DebugLog ("Supd", Label(c) + STR$(OC(c)))
    NEXT c
END IF

END SUB

SUB TempCorrection
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB TempCorrection                                                   |
'  |                                                                       |
'  |  Adjust temperature correction according to van't Hoff equation       |
'  |  given the temperature of the current system, the temperature         |
'  |  that the thermodynamic constants were measured in, and the change    |
'  |  in enthalpy for each reaction.                                       |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

' *-----------------------------------------------*
' |  Simple Variables Local to SUB DefineProblem  |
' *-----------------------------------------------*
DIM i AS INTEGER            '  Array/loop index variable
DIM T0 AS DOUBLE            '  Temporary variable
DIM T1 AS DOUBLE            '  Temporary variable
DIM T2 AS DOUBLE            '  Temporary variable

T1# = 1 / SysTemp

FOR i = NFirstS TO NSpecies
   ' *--------------------------------*
   ' |  Check Temperature correction  |
   ' |  flag for this species         |
   ' *--------------------------------*
   IF TempCorr(i) THEN
      ' *---------------------------------*
      ' |  Flag is TRUE, make correction  |
      ' *---------------------------------*
      T0# = 1 / ECTemp(i)
      T2# = DeltaH(i) * (T0# - T1#) / Rcon#
      hEC(i) = EC(i) * EXP(T2#)
   ELSE
      ' *--------------------------------------------------*
      ' |  Flag is FALSE, just copy the uncorrected value  |
      ' *--------------------------------------------------*
      hEC(i) = EC(i)
   END IF
NEXT i

FOR i = 1 TO NPhase
   ' *--------------------------------*
   ' |  Check Temperature correction  |
   ' |  flag for this phase           |
   ' *--------------------------------*
   IF STempCorr(i) THEN
      ' *---------------------------------*
      ' |  Flag is TRUE, make correction  |
      ' *---------------------------------*
      T0# = 1 / SECTemp(i)
      T2# = SDeltaH(i) * (T0# - T1#) / Rcon#
      hSEC(i) = SEC(i) * EXP(T2#)
   ELSE
      ' *--------------------------------------------------*
      ' |  Flag is FALSE, just copy the uncorrected value  |
      ' *--------------------------------------------------*
      hSEC(i) = SEC(i)
   END IF
NEXT i

END SUB

DEFSNG A-Z
SUB WHAMCalcDiffuse (HS%)


'##############################################################################
WWSR7:     ' called by WWSR5                                                 ##
'##          calcs binding by DL accumulation, for FA and HA                 ##
'##############################################################################

' First clear the DDL array for this HS%

FOR j% = 1 TO NSP%
   wCDDL(HS%, j%) = 0
NEXT j%

'CALL PrintA2(wCDDL(), NSP%)

TOTCHN = ZED(HS%) * wTH(HS%)                      ' total charge to be
                                                 ' neutralised, per litre

TDCONC = -TOTCHN / DVOL(HS%)                     ' total conc of counterions
                                                 ' per litre of diffuse layer

' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF TDCONC < 0 THEN GOTO WWSR7L1                       ' anions attracted


'**********************************************************************
' Come to here if humics have a net negative charge (cations attracted)

WWSR7L2:

TDCONCCALC = 0

FOR X% = 1 TO NSP%
IF wCH%(X%) < 0 THEN wCDDL(HS%, X%) = 0
' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF wCH%(X%) < 0 THEN GOTO WWSR7L3
wCDDL(HS%, X%) = wC(X%) * (wRATIO(HS%) ^ (wCH%(X%)))
TDCONCCALC = TDCONCCALC + (wCDDL(HS%, X%) * wCH%(X%))
WWSR7L3:
NEXT X%

'CALL PrintA2(wCDDL(), NSP%)
'CALL PrintA(wC(), NSP%)

 TDCONCERR = (2 * (TDCONC - TDCONCCALC) / (TDCONC + TDCONCCALC)) ^ 2
' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF TDCONCERR < CRIT THEN GOTO WWSR7L5


' Adjust wRATIO and re-try

wRATIO(HS%) = ((wRATIO(HS%) * TDCONC / TDCONCCALC) + wRATIO(HS%)) / 2
GOTO WWSR7L2



'*********************************************************************
' Come to here if humics have a net positive charge (anions attracted)

WWSR7L1:

TDCONCCALC = 0
FOR X% = 1 TO NSP%
IF wCH%(X%) > 0 THEN wCDDL(HS%, X%) = 0
' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF wCH%(X%) > 0 THEN GOTO WWSR7L4
wCDDL(HS%, X%) = wC(X%) * (wRATIO(HS%) ^ (-wCH%(X%)))
TDCONCCALC = TDCONCCALC + (wCDDL(HS%, X%) * wCH%(X%))
WWSR7L4:
NEXT X%
'
' RCS EDIT:  The following line was rewritten to avoid division
'            by zero error when TCCONC = TDCONCCALC
'TDCONCERR = (2 * (TDCONC - TDCONCCALC) / (TDCONC + TDCONCCALC)) ^ 2
TDCONCERR = 2 * (TDCONC - TDCONCCALC) / (TDCONC + TDCONCCALC)
IF TDCONCERR <> 0 THEN TDCONCERR = TDCONCERR ^ 2
' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF TDCONCERR < CRIT THEN GOTO WWSR7L5


' Adjust wRATIO and re-try

wRATIO(HS%) = ((wRATIO(HS%) * TDCONC / TDCONCCALC) + wRATIO(HS%)) / 2
GOTO WWSR7L1

WWSR7L5:
'RETURN

END SUB

SUB WHAMCalcSpecies
'STATIC wISTR!

'##############################################################################
WWSR5:     ' called by WWSR4 - calls WWSR6 and WWSR7                         ##
'##          calcs OH- activity, act coeffs, CO3,2- if pCO2 fixed            ##
'##          calcs activities, concs of inorganic complexes                  ##
'##          every 2nd itertn, calls WWSR6 to get NU's and Z's for HS        ##
'##          calls WWSR7 to get binding by DL accumulation (FA,HA)           ##
'##############################################################################

' Calculate activity of OH- from A(1), temp, deltaH

LKW = -14 + (2935 * (.003354 - (1 / TEMPK)))
wA(51) = 10 ^ (LKW) / wA(1)


'************************************************************
' Calculate activity coefficients using extended Debye-Huckel

wGAMMA(0) = 1
ACTA = .27 + (.0008 * TEMPK)                          ' const A in D-H
ACTB = .33                                            ' const B in D-H
ActC = -ACTA * SQR(wISTR)
ACTD = ACTB * SQR(wISTR)
wGAMMA(1) = 10 ^ (ActC / (1 + (3 * ACTD)))             ' act coeff M+/-
wGAMMA(2) = 10 ^ (4 * ActC / (1 + (6 * ACTD)))         ' act coeff M2+/-
wGAMMA(3) = 10 ^ (9 * ActC / (1 + (9 * ACTD)))         ' act coeff M3+/-
wGAMMA(4) = 10 ^ (16 * ActC / (1 + (12 * ACTD)))       ' act coeff M4+/-

'CALL PrintA(wA(), NSP%)
'**********************************************
' Calculate concentrations of inorganic species

' First do the carbonate system

' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF PCO2% = 999 THEN GOTO WWSR5L1                            ' pCO2 not specified
'
ENTHTERM = 220 * (-.962) * (.003354 - (1 / TEMPK))
wA(55) = PCO2% / 10 ^ (18.149 + ENTHTERM) / wA(1) / wA(1) ' activity of CO32-


WWSR5L1:

FOR X% = 101 TO NSP%

' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF wS1%(X%) = 0 THEN GOTO WWSR5L2                           ' complex not defined

' RCS EDIT:
' The following three lines were rewritten with explicit GOTO
IF wS1%(X%) > 0 AND wA(wM1%(X%)) = 0 THEN GOTO WWSR5L2        ' no calculation if
IF wS2%(X%) > 0 AND wA(wM2%(X%)) = 0 THEN GOTO WWSR5L2        ' contributing species
IF wS3%(X%) > 0 AND wA(wM3%(X%)) = 0 THEN GOTO WWSR5L2        ' absent

ENTHTERM = 220 * wDH(X%) * (.003354 - (1 / TEMPK))

LOGACT = wLK(X%) + ENTHTERM
LOGACT = LOGACT + (wS1%(X%) * LOG10S(wA(wM1%(X%))))
LOGACT = LOGACT + (wS2%(X%) * LOG10S(wA(wM2%(X%))))

IF wS3%(X%) > 0 THEN LOGACT = LOGACT + (wS3%(X%) * LOG10S(wA(wM3%(X%))))

wA(X%) = 10 ^ (LOGACT)

WWSR5L2:
NEXT X%

'CALL PrintA(wA(), NSP%)

' Calculate concns from activities

FOR X% = 1 TO NSP%
CHARGE% = ABS(wCH%(X%))
wC(X%) = wA(X%) / wGAMMA(CHARGE%)
NEXT X%

'CALL PrintA(wC(), NSP%)

'***************************************************
' Calculate concentrations of humic-bound components

' Check if the iteration no. (NUMIT%) is a multiple of 2 ; if not, return
' Note that TESTITER, TESTITER% are also used in subroutine WWSR3

TESTITER = NUMIT% / 2
TESTITER% = INT(NUMIT% / 2)
' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF TESTITER > TESTITER% THEN GOTO WWSR5L5             ' no calc this iter; return


FOR HS% = 1 TO 2

' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF wTH(HS%) = 0 THEN GOTO WWSR5L3                      ' no calcn if HS absent

'GOSUB WWSR6                                      ' calculates wZCALC's and NU's
CALL WHAMCalcSpecific(HS%)


' Calculate concns of specifically bound species per litre from NU and [HS]

FOR X% = 1 TO NSP%
   wCHC(HS%, X%) = wNU(HS%, X%) * wTH(HS%)
NEXT X%

'CALL PrintA2(wCHC(), NSP%)

WWSR5L3:
NEXT HS%


' Calculate maximum volumes of HA and FA diffuse layers

FOR HS% = 1 TO 2

' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF wTH(HS%) = 0 THEN GOTO WWSR5L4                           ' no calcn if HS absent

VTERM1 = wRADIUS(HS%) + (3.04E-10 / SQR(wISTR))
VTERM2 = (VTERM1 ^ 3) - ((wRADIUS(HS%)) ^ 3)
VTERM3 = 4.19 * VTERM2
DDLVOL(HS%) = 6E+23 * VTERM3 * (1000 / wMOLWT(HS%))    ' litres/gHS
DVOLMAX(HS%) = DDLVOL(HS%) * wTH(HS%)


' Adjust diffuse layer volume for low ZED

ZTERM = wKZED * ABS(ZED(HS%))
DVOLMAX(HS%) = DVOLMAX(HS%) * ZTERM / (1 + ZTERM)     ' max vol, litres/litre

WWSR5L4:
NEXT HS%


' Calculate the actual diffuse layer volumes, using DLF

DENOM = 1 + ((DVOLMAX(1) + DVOLMAX(2)) / DLF)

FOR HS% = 1 TO 2
DVOL(HS%) = DVOLMAX(HS%) / DENOM
NEXT HS%

VOLSOL = 1 - DVOL(1) - DVOL(2)

FOR HS% = 1 TO 2
IF wTH(HS%) > 0 THEN
   'GOSUB WWSR7                       ' calc DDL concns
   CALL WHAMCalcDiffuse(HS%)
END IF
NEXT HS%


WWSR5L5:


END SUB

SUB WHAMCalcSpecific (HS%)

'##############################################################################
WWSR6:     ' called by WWSR5                                                 ##
'##          calcs NU (specific binding) and Z for FA, HA                    ##
'##############################################################################

W = wP(HS%) * LOG10S(wISTR)                           ' e'static interaction factor


FOR i% = 1 TO 8
TEMPVAL = wKH(HS%, i%) * EXP(2 * W * ZED(HS%)) / wA(1)
wFP(i%) = 1 / (1 + TEMPVAL)                       ' protonation factors
NEXT i%

'CALL PrintA(wFP(), 8)

'***************************
' Binding at bidentate sites

FOR j% = 1 TO 12                                 ' do each site in turn

   SYTE1% = wSITE%(j%, 1)                            ' identifies proton sites that
   SYTE2% = wSITE%(j%, 2)                            ' make up the bidentate sites


   SUMBITERM = 1

   FOR k% = 1 TO NSP%
      ' RCS EDIT:
      ' The following line was rewritten with explicit GOTO
      IF wBINDTEST%(k%) = 0 THEN GOTO WWSR6L1
         TEMPVAL = 2 * W * ZED(HS%) * (2 - wCH%(k%))
         TEMPVAL = EXP(TEMPVAL)
         TEMPVAL = wKMH(HS%, SYTE1%, k%) * wKMH(HS%, SYTE2%, k%) * wA(k%) * TEMPVAL
         wBIBTERM(k%) = TEMPVAL * wFP(SYTE1%) * wFP(SYTE2%) / (wA(1) ^ 2)
         SUMBITERM = SUMBITERM + wBIBTERM(k%)
WWSR6L1:
   NEXT k%
   'CALL PrintA(wBIBTERM(), NSP%)


   SUMBITHETA = 0                                   ' preparing to sum
   BIMETCH = 0                                      ' theta's and charges

   FOR k% = 1 TO NSP%
      ' RCS EDIT:
      ' The following line was rewritten with explicit GOTO
      IF wBINDTEST%(k%) = 0 THEN GOTO WWSR6L2
         BITHETA = wBIBTERM(k%) / SUMBITERM
         SUMBITHETA = SUMBITHETA + BITHETA
         wBISNU(j%, k%) = BITHETA * wSM(HS%, j%)
         BIMETCH = BIMETCH + (BITHETA * wCH%(k%))
WWSR6L2:
   NEXT k%


   PROT2CH = 2 * wFP(SYTE1%) * wFP(SYTE2%) * (1 - SUMBITHETA)   ' 2 bound H+
   PROT1CH1 = wFP(SYTE1%) * (1 - wFP(SYTE2%)) * (1 - SUMBITHETA)' H+ bound at site 1
   PROT1CH2 = wFP(SYTE2%) * (1 - wFP(SYTE1%)) * (1 - SUMBITHETA)' H+ bound at site 2


   BINETCH = BIMETCH + PROT2CH + PROT1CH1 + PROT1CH2 - 2      ' net charge
   wBIZ(j%) = BINETCH * wSM(HS%, j%)                            ' this site's charge

NEXT j%


' Now calculate total amounts bound, and total net charge (bidentate sites)

FOR k% = 1 TO NSP%
' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF wBINDTEST%(k%) = 0 THEN GOTO WWSR6L3
wBINU(k%) = 0

FOR j% = 1 TO 12
wBINU(k%) = wBINU(k%) + wBISNU(j%, k%)
NEXT j%

WWSR6L3:
NEXT k%


BIZCALC = 0
FOR j% = 1 TO 12
BIZCALC = BIZCALC + wBIZ(j%)
NEXT j%


'*****************************
' Binding at monodentate sites

FOR j% = 1 TO 8

SUMMONTERM = 1

FOR k% = 1 TO NSP%
' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF wBINDTEST%(k%) = 0 THEN GOTO WWSR6L4
TEMPVAL = 2 * W * ZED(HS%) * (1 - wCH%(k%))
TEMPVAL = EXP(TEMPVAL)
TEMPVAL = wKMH(HS%, j%, k%) * wA(k%) * TEMPVAL
wMONBTERM(k%) = TEMPVAL * wFP(j%) / wA(1)
SUMMONTERM = SUMMONTERM + wMONBTERM(k%)
WWSR6L4:
NEXT k%


SUMMONTHETA = 0                                  ' preparing to sum
MONMETCH = 0                                     ' theta's and charges

FOR k% = 1 TO NSP%
' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF wBINDTEST%(k%) = 0 THEN GOTO WWSR6L5
MONTHETA = wMONBTERM(k%) / SUMMONTERM
SUMMONTHETA = SUMMONTHETA + MONTHETA
wMONSNU(j%, k%) = MONTHETA * wSP(HS%, j%)
MONMETCH = MONMETCH + (MONTHETA * wCH%(k%))
WWSR6L5:
NEXT k%


PROT1CH = wFP(j%) * (1 - SUMMONTHETA)             ' H+ bound
MONNETCH = MONMETCH + PROT1CH - 1                ' net charge
wMONZ(j%) = MONNETCH * wSP(HS%, j%)

NEXT j%


' Now calculate total amounts bound, and total net charge (monodentate sites)

FOR k% = 1 TO NSP%
' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF wBINDTEST%(k%) = 0 THEN GOTO WWSR6L6
wMONNU(k%) = 0

FOR j% = 1 TO 8
wMONNU(k%) = wMONNU(k%) + wMONSNU(j%, k%)
NEXT j%

WWSR6L6:
NEXT k%


MONZCALC = 0

FOR j% = 1 TO 8
MONZCALC = MONZCALC + wMONZ(j%)
NEXT j%


'********************************************
' Overall summation ; bidentate + monodentate

FOR k% = 1 TO NSP%
' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF wBINDTEST%(k%) = 0 THEN GOTO WWSR6L7
wNU(HS%, k%) = wBINU(k%) + wMONNU(k%)
WWSR6L7:
NEXT k%


' Calculated value of Z, and Z error term

wZCALC(HS%) = BIZCALC + MONZCALC
' RCS EDIT: June 1, 1997
' The following line crashed when ZED = wZCALC (Division by zero).
' An equivalent line was substituted that overcomes this problem
'wZERR(HS%) = (2 * (ZED(HS%) - wZCALC(HS%)) / (ZED(HS%) + wZCALC(HS%))) ^ 2
TMP1! = ZED(HS%) - wZCALC(HS%)
Tmp2! = ZED(HS%) + wZCALC(HS%)
TMP3! = 2 * (TMP1!) / (Tmp2!)
IF TMP3! = 0 THEN
   wZERR(HS%) = 0
ELSE
   wZERR(HS%) = TMP3! ^ 2
END IF

'RETURN


END SUB

SUB WHAMDefine

'##############################################################################
WWSR1:     ' called by main program                                          ##
'##          DIMensions arrays                                               ##
'##          reads in data from data base                                    ##
'##          sets up Model V                                                 ##
'##############################################################################

' DIMension arrays ; inorganic chemistry

' RCS EDIT:
' The following DIMENSION statements must be near begining of program
'DIM N%(NSP%)                 ' numerical identification of species
'DIM wN$(NSP%)                 ' nominal identification of species
'DIM T(NSP%)                  ' total concentrations of components
'DIM TCALC(NSP%)              ' total calculated concentrations of components
'DIM A(NSP%)                  ' activities of solution species
'DIM wC(NSP%)                  ' concentations of solution species
'DIM CH%(NSP%)                ' charges on solution species
'DIM wM1%(NSP%)                ' master species identifiers
'DIM wM2%(NSP%)                ' master species identifiers
'DIM wM3%(NSP%)                ' master species identifiers
'DIM wS1%(NSP%)                ' stoichiometries
'DIM wS2%(NSP%)                ' stoichiometries
'DIM S3%(NSP%)                ' stoichiometries
'DIM LK(NSP%)                 ' equilibrium constants
'DIM DH(NSP%)                 ' enthalpies
'DIM GAMMA(5)                 ' activity coefficients
'
'
'' DIMension humic arrays ; 1 = HA, 2 = FA
'
'DIM wSP(2, 8)                 ' proton-binding sites not forming bidentate sites
'DIM wSM(2, 12)                ' bidentate metal-binding sites
'DIM wKH(2, 8)                 ' K's for proton-binding
'DIM PKH(2, 8)                ' PK's for proton-binding
'DIM wKMH(2, 8, NSP%)          ' K's for metal-proton exchange
'DIM wSITE%(12, 2)             ' proton sites making bidentate sites
'DIM wPKMHA(2, NSP%)           ' metal-proton exchange constants
'DIM PKMHB(2, NSP%)           ' metal-proton exchange constants
'DIM FP(8)                    ' dissociation factors for proton-binding
'DIM BIBTERM(NSP%)            ' terms calc'd in finding bidentate theta's
'DIM MONBTERM(NSP%)           ' terms calc'd in finding monodentate theta's
'DIM BISNU(12, NSP%)          ' values of NU at each bidentate site
'DIM MONSNU(8, NSP%)          ' values of NU at each monodentate site
'DIM BINU(NSP%)               ' overall NU for each metal/bidentate sites
'DIM MONNU(NSP%)              ' overall NU for each metal/monodentate sites
'DIM NU(2, NSP%)              ' overall NU for each complexed species
'DIM wCHC(2, NSP%)             ' conc of complexed species per litre
'DIM wCDDL(2, NSP%)            ' concs in HS DDL's
'DIM BIZ(12)                  ' net charge at each bidentate site
'DIM MONZ(8)                  ' net charge at each monodentate site
'DIM BINDTEST%(NSP%)          ' allows binding of a species to be tested


'****************************************************************
' Input of constants for inorganic speciation and HS complexation

'OPEN DBS$ + ".DBS" FOR INPUT AS #2
'
'INPUT #2, S$                                     ' database identifier
'INPUT #2, S$, wNCOOH(1), wPKHA(1), wPKHB(1), wDPKHA(1), wDPKHB(1), P(1), wFPR(1), RADIUS(1), MOLWT(1)' HA properties
'INPUT #2, S$, wNCOOH(2), wPKHA(2), wPKHB(2), wDPKHA(2), wDPKHB(2), P(2), wFPR(2), RADIUS(2), MOLWT(2)' FA properties
'INPUT #2, S$, DLF                                ' double layer overlap factor
'INPUT #2, S$, wKZED                               ' DL vol at low Z
'INPUT #2, S$, NODATA%                            ' number of species
'
'FOR i% = 1 TO NODATA%
'
'INPUT #2, X%, wN$(X%), CH%(X%), wM1%(X%), wM2%(X%), wM3%(X%), wS1%(X%), wS2%(X%), S3%(X%), LK(X%), DH(X%), wPKMHA(1, X%), wPKMHA(2, X%)
'
'IF wPKMHA(1, X%) < 999 THEN BINDTEST%(X%) = 1     ' identifies species that
'                                                 ' don't bind to HS
'NEXT i%
'
'CLOSE #2


'**********************************
' Set derived constants for Model V

FOR X% = 1 TO 2
wPKH(X%, 1) = wPKHA(X%) - (wDPKHA(X%) / 2)
wPKH(X%, 2) = wPKHA(X%) - (wDPKHA(X%) / 6)
wPKH(X%, 3) = wPKHA(X%) + (wDPKHA(X%) / 6)
wPKH(X%, 4) = wPKHA(X%) + (wDPKHA(X%) / 2)          ' individual site pK's
wPKH(X%, 5) = wPKHB(X%) - (wDPKHB(X%) / 2)
wPKH(X%, 6) = wPKHB(X%) - (wDPKHB(X%) / 6)
wPKH(X%, 7) = wPKHB(X%) + (wDPKHB(X%) / 6)
wPKH(X%, 8) = wPKHB(X%) + (wDPKHB(X%) / 2)
NEXT X%

FOR i% = 1 TO 2
FOR j% = 1 TO 8
wKH(i%, j%) = 10 ^ (-wPKH(i%, j%))                 ' convert pK's to K's
NEXT j%
NEXT i%


' Identify the proton-binding sites that combine to make the bidentate sites

wSITE%(1, 1) = 1:    wSITE%(1, 2) = 2
wSITE%(2, 1) = 1:    wSITE%(2, 2) = 4
wSITE%(3, 1) = 1:    wSITE%(3, 2) = 6
wSITE%(4, 1) = 1:    wSITE%(4, 2) = 8
wSITE%(5, 1) = 2:    wSITE%(5, 2) = 3
wSITE%(6, 1) = 2:    wSITE%(6, 2) = 5
wSITE%(7, 1) = 2:    wSITE%(7, 2) = 7
wSITE%(8, 1) = 3:    wSITE%(8, 2) = 4
wSITE%(9, 1) = 3:    wSITE%(9, 2) = 6
wSITE%(10, 1) = 3:   wSITE%(10, 2) = 8
wSITE%(11, 1) = 4:   wSITE%(11, 2) = 5
wSITE%(12, 1) = 4:   wSITE%(12, 2) = 7


' Set site concentrations

FOR i% = 1 TO 4
wSP(1, i%) = (1 - wFPR(1)) * wNCOOH(1) / 4
wSP(2, i%) = (1 - wFPR(2)) * wNCOOH(2) / 4
NEXT i%

FOR i% = 5 TO 8
wSP(1, i%) = (1 - wFPR(1)) * wNCOOH(1) / 8
wSP(2, i%) = (1 - wFPR(2)) * wNCOOH(2) / 8
NEXT i%

FOR i% = 1 TO 12
wSM(1, i%) = wFPR(1) * wNCOOH(1) / 16
wSM(2, i%) = wFPR(2) * wNCOOH(2) / 16
NEXT i%


' Set metal-proton exchange constants for individual sites

FOR i% = 1 TO NSP%

' RCS EDIT:
' The following two lines were rewritten with explicit GOTO's
'OLD LINE: IF wN$(I%) = "" THEN WWSR1L1                       ' species not defined
'OLD LINE: IF wPKMHA(1, I%) = 999 THEN WWSR1L1
IF wNs$(i%) = "" THEN GOTO WWSR1L1                       ' species not defined
IF wPKMHA(1, i%) = 999 THEN GOTO WWSR1L1

FOR j% = 1 TO 4
wPKMH = wPKMHA(1, i%)
wKMH(1, j%, i%) = 10 ^ (-wPKMH)
wPKMH = wPKMHA(2, i%)
wKMH(2, j%, i%) = 10 ^ (-wPKMH)
NEXT j%

FOR j% = 5 TO 8
wPKMH = (3 * wPKMHA(1, i%)) - 3                    ' 3*A - 3 conversion ; HA
wKMH(1, j%, i%) = 10 ^ (-wPKMH)
wPKMH = 3.96 * wPKMHA(2, i%)                       ' 3.96*A  conversion ; FA
wKMH(2, j%, i%) = 10 ^ (-wPKMH)
NEXT j%
GOTO WWSR1L2

' If come to here, then no binding of species I% can occur
WWSR1L1:

FOR j% = 1 TO 8
wKMH(1, j%, i%) = 0
wKMH(2, j%, i%) = 0
NEXT j%

WWSR1L2:
NEXT i%

'RETURN

END SUB

SUB WHAMIter
'##############################################################################
IF (iGlobalAl = 2) THEN
  ' Al is always controlled by solubility
  iALPPT = TRUE
  LKSOAL25 = 8.625
  'adj up to slightly higher Al solubility as suggested by Niva data
  LKSOAL25 = 9.76
  IF (iLogKsAl = TRUE) THEN LKSOAL25 = LogKsAl     ' set the Log K for Al solubility to command line value
  DHAL = -25000
  LKSOAL = LKSOAL25 + (0.219 * DHAL * ((1/298) - (1/SysTemp)))
ELSE
  ' assume no solubility, unless oversaturated and iGlobalAl = 1
  iALPPT = FALSE
END IF
IF (iGlobalFe = 2) THEN
  ' Fe is always controlled by solubility
  iFePPT = TRUE
  LKSOFE25 = 3.0
    LKSOFE25 = 1.0
    LKSOFE25 = 2.5
  DHFE = -25000
  LKSOFE = LKSOFE25 + (0.219 * DHFE * ((1/298) - (1/SysTemp)))
  LIAPFE = LOG10(wA(11)) + (3*PH)
ELSE
  ' assume no solubility, unless oversaturated and iGlobalFe = 1
  iFePPT = FALSE
END IF

IF (iGlobalPb = 2) THEN
  ' Pb is always controlled by solubility
  iPbPPT = TRUE
'  LKSOFE25 = 3.0
'    LKSOFE25 = 1.0
'    LKSOFE25 = 2.5
'  DHFE = -25000
'  LKSOFE = LKSOFE25 + (0.219 * DHFE * ((1/298) - (1/SysTemp)))
'  LIAPFE = LOG10(wA(11)) + (3*PH)
'  wA(55)
ELSE
  ' assume no solubility, unless oversaturated and iGlobalFe = 1
  iPbPPT = FALSE
END IF

IF (iGlobalCO2=TRUE) THEN
'    DIM iGlobalCO2 AS GLOBAL INTEGER
'    DIM iCO2PPT AS GLOBAL INTEGER
'    DIM LKH_CO2 AS GLOBAL DOUBLE         ' Log Henry's Law constant
'    DIM pPCO2 AS GLOBAL DOUBLE            ' negative log of the partial pressure of CO2
'    DIM Log_CO3 AS GLOBAL DOUBLE         ' Log of the CO3[2-] ion concentration
'    DIM Log_K_H2CO3 as global double     ' Log of the formation constant for H2CO3
    '
    ' Inorganic carbon is goverend by CO2 solubility
    '    Log_CO3 = pPCO2 + LKH_CO2 - Log(2* [H]) - Log_K_H2CO3
    '    Log_CO3 = pPCO2 + LKH_CO2 - 2 pH - Log_K_H2CO3
    '        where:
    '           pPCO3 is the negative log of the partial pressure for CO2
    '           LKH_CO2 is the log of the Henry's Law constant for CO2 ( -1.5)
    '           Log_K_H2CO3 is the log of the formation constant for H2CO3

    iCO2PPT = TRUE

ELSE
    iCO2PPT = FALSE
END IF


WWSR2:     ' called from main program - calls WWSR3 and WWSR4                ##
'##          sets initial trial values of master species                     ##
'##          initialises all activities and concns by calling WWSR4          ##
'##          controls pH improvement if pH not fixed                         ##
'##          controls level of precision                                     ##
'##          calls WWSR3 to calc speciation at a given pH                    ##
'##          tests for correct pH with CHRATIO if pH not fixed               ##
'##          sets up screen to report on progress of calculation             ##
'##############################################################################
' Set initial trial values

wISTR = .1                                          ' ionic strength

ZED(1) = -.0001                                  ' Z for HA
ZED(2) = -.0001                                  ' Z for FA

wRATIO(1) = 10                                    ' RATIO for HA
wRATIO(2) = 10                                    ' RATIO for FA

'PH = PHSTART                                     ' fixed or starting pH
wA(1) = OC(1) '* ActCoef(1)                        ' activity of H+
PH = -Log10D#(OC(1))

'FOR X% = 2 TO 50                                  ' species #1 is H+
''   wA(X%) = wT(X%)                                ' cationic master species
'   ii = CMap%(X%)
'   IF (ii <> 0) THEN
'      wA(X%) = OC(ii)
'      '
'      ' The following line is closer to WHAM method for debuggin
'      ' but is probably less efficient starting point
'      wA(X%) = KT(ii)
'   END IF
'NEXT X%
'CALL PrintA(wA(), 50)

TEMPK = SysTemp
FOR i = 2 TO NWHAM%
   ii = WHAMMAP%(i)
   IF (ii <> 0) THEN
      wA(ii) = OC(i) * ActCoef(i)
      wT(ii) = KT(i) / CtoM!(i)
      'wA(ii) = wT(ii)
   END IF
NEXT i
'CALL PrintA(wA(), NSP%)

wTH(1) = PercHA! * KT(iDOC)
wTH(2) = (1! - PercHA!) * KT(iDOC)

'*************************
' Calculate ionic strength

wISTR = 0                                           ' preparing to sum

FOR X% = 1 TO NSP%
IF wC(X%) > 0 THEN wISTR = wISTR + (.5 * (wC(X%) * wCH%(X%) * wCH%(X%)))
NEXT X%

IF (wISTR > 100) THEN wISTR = 100                        ' avoids initial high wISTR
IF (wISTR! = 0) THEN wISTR! = .1




'A(52) = T(52)                                    ' Cl
'A(53) = T(53)                                    ' NO3
'A(54) = T(54)                                    ' SO4
'A(56) = T(56) * .0001                            ' F
'A(57) = T(57) * 1E-15                            ' PO4

IF PCO2% = 999 THEN wA(55) = .000001 * wT(55)       ' CO3,2-



'**************************************
' Summary output to screen, and headers

'IF FALSE THEN
'LOCATE 13
'PRINT "                                                                "
'LOCATE 13
'PRINT "     SOURCE FILE : "; SF$
'PRINT "   PRECISION (%) : "; USING "#.##^^^^"; 1.0001 * PRECISION
'IF wTH(1) = 0 AND wTH(2) = 0 THEN PRINT "HUMIC SUBSTANCES : ABSENT"
'IF (wTH(1) + wTH(2)) > 0 THEN PRINT "HUMIC SUBSTANCES : PRESENT"
'IF PHFIX$ = "NO" THEN PRINT "              PH : VARIABLE"
'IF PHFIX$ = "YES" THEN PRINT "              PH : FIXED"
'IF PCO2% = 999 THEN PRINT "            PCO2 : VARIABLE"
'IF PCO2% < 999 THEN PRINT "            PCO2 : FIXED"
'
'
'PRINT
'
'PRINT "ITER"; TAB(15); "PH"; TAB(30); "IS"; TAB(45); "CHRATIO"; TAB(62); "DPH"
'END IF

'*************************************************************************
' Begin calculations

NUMIT% = 0                                       ' iteration counter

CRIT = .0001: GOSUB WWSR4                        ' initialization

' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF PHFIX$ = "NO" THEN GOTO WWSR2L1                    ' routine if pH variable

IF PHFIX$ = "YES" THEN CRTEST$ = "PHASE2"        ' routine for
IF PHFIX$ = "YES" THEN CRIT = (PRECISION / 100) ^ 2' fixed pH calcn
IF PHFIX$ = "YES" THEN GOSUB WWSR3
' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF PHFIX$ = "YES" THEN GOTO WWSR2L4                   ' fixed pH calcs done; return


'**********************************************
' Use pH adjusting routine, solving for each pH

WWSR2L1:

CRTEST$ = "PHASE1"

CRIT = .0001: GOSUB WWSR3                        ' find initial CHRATIO etc

DPHF = .2                                        ' initial pH adjust factor

' RCS EDIT:
' The following two lines were rewritten with explicit GOTO
IF CHRATIO < 1 THEN GOTO WWSR2L2
IF CHRATIO > 1 THEN GOTO WWSR2L3



'+++++++++++++++++++++++++++++++
' Come to here if CHRATIO is < 1

WWSR2L2:

DPH = DPHF / CHRATIO                             ' pH increment
IF DPH < .01 THEN CRIT = (PRECISION / 100) ^ 2   ' adjust precision
IF DPH < .0001 THEN CRTEST$ = "PHASE2"           ' refining
PH = PH - DPH                                    ' increment pH
GOSUB WWSR3                                      ' calculate speciation
' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF CRTEST$ = "PHASE2" AND ABS(CHRATIO - 1) < SQR(CRIT) THEN GOTO WWSR2L4       ' finished
IF CHRATIO > 1 THEN DPHF = DPHF / 3              ' change pH adjust factor
' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF CHRATIO > 1 THEN GOTO WWSR2L3                      ' go to > routine
GOTO WWSR2L2                                     ' continue to increment pH


'+++++++++++++++++++++++++++++++
' Come to here if CHRATIO is > 1

WWSR2L3:

DPH = DPHF * CHRATIO                             ' pH increment
IF DPH < .01 THEN CRIT = (PRECISION / 100) ^ 2   ' adjust precision
IF DPH < .0001 THEN CRTEST$ = "PHASE2"           ' refining
PH = PH + DPH                                    ' increment pH
GOSUB WWSR3                                      ' calculate speciation
' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF CRTEST$ = "PHASE2" AND ABS(CHRATIO - 1) < SQR(CRIT) THEN GOTO WWSR2L4       ' finished
IF CHRATIO < 1 THEN DPHF = DPHF / 3              ' change pH adjust factor
' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF CHRATIO < 1 THEN GOTO WWSR2L2                      ' go to < routine
GOTO WWSR2L3                                     ' continue to increment pH

WWSR2L4:


FOR i = 1 TO NSP%
   ii = CMap%(i)
   IF (ii <> 0) THEN
      OC(ii) = wC(i)
      wOrgBound!(ii) = DVOL(1) * wCDDL(1, i) + (DVOL(2) * wCDDL(2, i) + wCHC(1, i) + wCHC(2, i))
      HABound(ii) = DVOL(1) * wCDDL(1, i) + wCHC(1, i)
      FABound(ii) = DVOL(2) * wCDDL(2, i) + wCHC(2, i)
   END IF
NEXT i
'RETURN

''
'  Check if gibbsite solubility is exceeded
'
'
' No point in checking if there is no aluminum in the simulation
IF (wT(5) > 0) AND (iGlobalAl = 1) THEN
    LKSOAL25 = 8.625
    'adj up to slightly higher Al solubility as suggested by Niva data
    LKSOAL25 = 9.76
    IF (iLogKsAl = TRUE) THEN LKSOAL25 = LogKsAl     ' set the Log K for Al solubility to command line value
    DHAL = -25000
    LKSOAL = LKSOAL25 + (0.219 * DHAL * ((1/298) - (1/TEMPK)))
    LIAPAL = LOG10(wA(5)) + (3*PH)
    IF LIAPAL >= LKSOAL THEN
      iALPPT = TRUE
      GOSUB WWSR2             ' recalculate with Al(OH)3
    END IF
END IF
''
'  Check if Fe(OH)3 solubility is exceeded
'
'
' No point in checking if there is no iron in the simulation
IF (wT(11) > 0) AND (iGlobalFe = 1) THEN
    LKSOFE25 = 3.0
      LKSOFE25 =1
      LKSOFE25 = 2.5
    DHFE = -25000
    LKSOFE = LKSOFE25 + (0.219 * DHFE * ((1/298) - (1/TEMPK)))
    LIAPFE = LOG10(wA(11)) + (3*PH)
    IF LIAPFE >= LKSOFE THEN                        '<----------- shouldn't this be LIAPFE  ???  EDIT DEBUG ???
      iFePPT = TRUE
      GOSUB WWSR2             ' recalculate with Fe(OH)3
    END IF
END IF

' check Lead solubility
IF (wT(21) > 0) AND (iGlobalPb = 1) THEN
    LKSOPB = -13.419
    LIAPPB = LOG10(wA(21)) + LOG10(wA(55))
    IF LIAPPB > LKSOPB THEN
        iPbPPT = TRUE
        GOSUB WWSR2
    END IF
END IF

EXIT SUB



'##############################################################################
WWSR3:     ' called by WWSR3 - calls WWSR4                                   ##
'##          controls improvement of activities (not H+) and ZFA,ZHA         ##
'##          calls WWSR4 to do mass and charge calculations                  ##
'##          in PHASE 1 returns when CHRATIO has become nearly-constant      ##
'##          in PHASE 2 returns when mass and charge balances O.K.           ##
'##          reports progress of calculation to screen                       ##
'##############################################################################

'wA(1) = 10 ^ (-PH)                                ' activity of H+

LASTCHRATIO = 0: LAST1CHRATIO = 0                ' resetting values


'*************************************************************************
' Begin iterative cycle ; come back to here until convergence is achieved

WWSR3L1:

NUMIT% = NUMIT% + 1                              ' update iteration counter


'*************************
' Calculate ionic strength

wISTR = 0                                           ' preparing to sum

FOR X% = 1 TO NSP%
IF wC(X%) > 0 THEN wISTR = wISTR + (.5 * (wC(X%) * wCH%(X%) * wCH%(X%)))
NEXT X%

IF wISTR > 100 THEN wISTR = 100                        ' avoids initial high wISTR


'*********************
' Improve trial values
'CALL PrintA(wA(), NSP%)

FOR X% = 2 TO 50
' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF wT(X%) = 0 THEN GOTO WWSR3L2
CONCFACTOR = wT(X%) / wTCALC(X%)
wA(X%) = wA(X%) * SQR(CONCFACTOR)                  ' cation activities
WWSR3L2:
NEXT X%

' check for Al(OH)3 pptn
IF iALPPT THEN
   ' Al controlled by Al(OH)3
   wA(5) = (10^(LKSOAL)) * (wA(1)^3)
   ' SI = log (IAP / KS)
   ' SI = Log(IAP) - Log(KS)
END IF
' check for Fe(OH)3 pptn
IF iFePPT THEN
   ' Fe controlled by Fe(OH)3
   wA(11) = (10^(LKSOFE)) * (wA(1)^3)
END IF
IF iPbPPT THEN
    LKSOPB = -13.419
    wA(21) = 10^LKSOPB / wA(55)
END IF



FOR X% = 52 TO 100
' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF wT(X%) = 0 THEN GOTO WWSR3L3
CONCFACTOR = wT(X%) / wTCALC(X%)
wA(X%) = wA(X%) * CONCFACTOR                       ' anion activities
WWSR3L3:
NEXT X%

IF iCO2PPT THEN
    ' CO3 is controlled by CO2 solubility

    LKH_CO2 = -1.5
    'pPCO2 = -3.5
    Log_K_H2CO3 = 16.681
    Log_CO3 = pPCO2 + LKH_CO2 + 2*pH - Log_K_H2CO3
    wA(55) = 10^(Log_CO3)


END IF

'CALL PrintA(wA(), NSP%)

' Improve ZED values only when iteration no. is multiple of 2 (see WWSR5)

IF wTH(1) > 0 AND TESTITER = TESTITER% THEN ZED(1) = ZED(1) + ((wZCALC(1) - ZED(1)) / 5)
IF wTH(2) > 0 AND TESTITER = TESTITER% THEN ZED(2) = ZED(2) + ((wZCALC(2) - ZED(2)) / 5)


LAST1CHRATIO = LASTCHRATIO                       ' record values to avoid
LASTCHRATIO = CHRATIO                            ' finding complete soln. in
                                                 ' first phase - see below

' Calculate total concns of master species
GOSUB WWSR4


'********************************
' Output current status to screen
'IF FALSE THEN
'LOCATE 20
'PRINT NUMIT%; TAB(13); USING "##.###"; PH; TAB(28);
'PRINT USING "#.###^^^^"; wISTR; TAB(44);
'PRINT USING "##.####"; CHRATIO; TAB(60); DPH
'END IF

'***************************************************************************
' Test whether CHRATIO has been found to acceptable precision in first phase

' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF CRTEST$ = "PHASE1" AND (LASTCHRATIO / CHRATIO) > .9 AND (LASTCHRATIO / CHRATIO) < 1.1 AND (LAST1CHRATIO / CHRATIO) > .9 AND (LAST1CHRATIO / CHRATIO) < 1.1 THEN GOTO WWSR3L6      ' return

' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF CRTEST$ = "PHASE1" AND LAST1CHRATIO = CHRATIO THEN GOTO WWSR3L6             ' oscillating - return

'*****************************************************************
' Test cations, anions, Z's for convergence - note that final test
' is CHRATIO, done in pH-adjusting routine

FOR X% = 2 TO 50                                                ' cations
   ' RCS EDIT:
   ' The following line was rewritten with explicit GOTO
   IF wT(X%) = 0 THEN GOTO WWSR3L4
      IF (iALPPT AND X% = 5) THEN
          wT(X%) = wTCALC(X%)
          GOTO WWSR3L4  ' skip mass balance on total Al if we are precip Al(OH)3
      END IF
      IF (iFePPT AND X% = 11) THEN
          wT(X%) = wTCALC(X%)
          GOTO WWSR3L4  ' skip mass balance on total Fe if we are precip fe(OH)3
      END IF
      IF (iPbPPT AND X% = 21) THEN
          wT(X%) = wTCALC(X%)
          GOTO WWSR3L4  ' skip mass balance on total Pb if we are precip PbCO3
      END IF

      CONCERR = (2 * (wT(X%) - wTCALC(X%)) / (wT(X%) + wTCALC(X%))) ^ 2   ' error term
      ' RCS EDIT:
      ' The following line was rewritten with explicit GOTO
      IF CONCERR > CRIT THEN GOTO WWSR3L1                                  ' re-iterate
WWSR3L4:
NEXT X%

FOR X% = 52 TO 100                                              ' anions
   ' RCS EDIT:
   ' The following line was rewritten with explicit GOTO
   IF wT(X%) = 0 THEN GOTO WWSR3L5
      IF (iCO2PPT) AND (X% = 55) THEN
          wT(X%) = wTCALC(X%)
          GOTO WWSR3L5:  ' skip mass balance on total CO3 if we are precip CO2
      END IF
      CONCERR = (2 * (wT(X%) - wTCALC(X%)) / (wT(X%) + wTCALC(X%))) ^ 2   ' error term
      ' RCS EDIT:
      ' The following line was rewritten with explicit GOTO
      IF CONCERR > CRIT THEN GOTO WWSR3L1                                  ' re-iterate
WWSR3L5:
NEXT X%

' RCS EDIT:
' The following line was rewritten with explicit GOTO
IF (wTH(1) + wTH(2)) = 0 THEN GOTO WWSR3L6              ' no humics ; return

' RCS EDIT:
' The following two lines were rewritten with explicit GOTO
IF wZERR(1) > CRIT THEN GOTO WWSR3L1                   ' re-iterate
IF wZERR(2) > CRIT THEN GOTO WWSR3L1                   ' re-iterate


WWSR3L6:
RETURN



'##############################################################################
WWSR4:     ' called by WWSR2 and WWSR3 - calls WWSR5                         ##
'##          calls WWSR5 to calc activities, concs of complexes and          ##
'##               amounts bound by FA and HA                                 ##
'##          sums to get total calc'd concns                                 ##
'##          calcs +ve and -ve charge and CHRATIO                            ##
'##############################################################################

' Calc activities and concns of inorganic complexes, amnts bound to FA & HA

'GOSUB WWSR5
CALL WHAMCalcSpecies

'******************************************************************
' Do summations to obtain total calculated concns of master species

FOR X% = 1 TO 100

   IF wA(X%) = 0 THEN wTCALC(X%) = 0
   ' RCS EDIT:
   ' The following line was rewritten with explicit GOTO
   IF wA(X%) = 0 THEN GOTO WWSR8L1

      ' RCS EDIT:
      ' The following two lines were rewritten with explicit GOTO
      IF X% = 1 THEN GOTO WWSR8L1
      IF X% = 51 THEN GOTO WWSR8L1


         wTCALC(X%) = 0                                         ' preparing to sum

         FOR y% = 1 TO NSP%

            IF wM1%(y%) = X% THEN wTCALC(X%) = wTCALC(X%) + (VOLSOL * wC(y%) * wS1%(y%)) + (wCHC(1, y%) * wS1%(y%)) + (wCHC(2, y%) * wS1%(y%)) + (DVOL(1) * wCDDL(1, y%) * wS1%(y%)) + (DVOL(2) * wCDDL(2, y%) * wS1%(y%)) _
                                                                                         ' DDL (2)

            IF wM2%(y%) = X% THEN wTCALC(X%) = wTCALC(X%) + (VOLSOL * wC(y%) * wS2%(y%)) + (wCHC(1, y%) * wS2%(y%)) + (wCHC(2, y%) * wS2%(y%)) + (DVOL(1) * wCDDL(1, y%) * wS2%(y%)) + (DVOL(2) * wCDDL(2, y%) * wS2%(y%)) _
                                                                                         ' DDL (2)

            IF wM3%(y%) = X% THEN wTCALC(X%) = wTCALC(X%) + (VOLSOL * wC(y%) * wS3%(y%)) + (wCHC(1, y%) * wS3%(y%)) + (wCHC(2, y%) * wS3%(y%)) + (DVOL(1) * wCDDL(1, y%) * wS3%(y%)) + (DVOL(2) * wCDDL(2, y%) * wS3%(y%)) _
                                                                                         ' DDL (2)

         NEXT y%

WWSR8L1:
NEXT X%


'******************************************
' Calculate +ve and -ve charges and CHRATIO

POSCH = 0:     NEGCH = 0                              ' preparing to sum

FOR X% = 1 TO NSP%
IF wCH%(X%) > 0 THEN POSCH = POSCH + (wC(X%) * wCH%(X%))
IF wCH%(X%) < 0 THEN NEGCH = NEGCH - (wC(X%) * wCH%(X%))
NEXT X%

CHRATIO = (POSCH / NEGCH)

RETURN



END SUB

'ÕÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¸
'³ GAUSS.BAS                                                                ³
'³ Solves a set of simultaneous linear equations by Gaussian Elimination.   ³
'³ This entire module can be loaded into any program, and the solution can  ³
'³    be obtained by a single call to SUB GAUSS.  See the subroutine GAUSS  ³
'³    for further information.                                              ³
'³                                                                          ³
'ÔÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¾

DECLARE SUB PivotMatrix (A() AS DOUBLE, Order() AS INTEGER, P%, n%)

DEFINT A-Z

SUB Gauss (A() AS DOUBLE, B() AS DOUBLE, X() AS DOUBLE, n%, iSing%)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB Gauss                                                            |
'  |                                                                       |
'  |  This subroutine calls the related subroutines to solve a set of      |
'  |  simulataneous linear equations.  The matrix A() contains nxn         |
'  |  coefficients for n equations  expressed as A * X = B.                |
'  |                                                                       |
'  |  The coeficient matrix A() is destroyed on output.                    |
'  |                                                                       |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

REDIM Order(n) AS INTEGER, SCALE(n) AS DOUBLE

IF n > 0 THEN
   GOSUB ScaleOrder
   GOSUB ForwardElim
   GOSUB BackSub
END IF

ERASE Order, SCALE
EXIT SUB

ScaleOrder:
   FOR ROW = 1 TO n
      Order(ROW) = ROW
      SCALE(ROW) = ABS(A(ROW, 1))

      FOR COL = 1 TO n
         ' Get the largest element in A(Row, )
         IF ABS(A(ROW, COL)) > SCALE(ROW) THEN
            ' Copy the largest value of A(Row,Col) into Scale(Row)
            SCALE(ROW) = ABS(A(ROW, COL))
         END IF
      NEXT COL

      ' Now scale this row of the matrix so that the largest element is
      '    set to "1".  This helps to reduce roundoff error in subsequent steps.

      FOR COL = 1 TO n
         A(ROW, COL) = A(ROW, COL) / SCALE(ROW)
      NEXT COL
      B(ROW) = B(ROW) / SCALE(ROW)
   NEXT ROW
RETURN

ForwardElim:
   FOR P = 1 TO n - 1
      PivotMatrix A(), Order(), P, n

      FOR ROW = P + 1 TO n
         P1 = Order(ROW)
         P2 = Order(P)
         Factor# = A(P1, P) / A(P2, P)

         FOR COL = P + 1 TO n
            A(P1, COL) = A(P1, COL) - Factor# * A(P2, COL)
         NEXT COL
         B(P1) = B(P1) - Factor# * B(P2)
      NEXT ROW
   NEXT P

   ' Check for singularities
   Sum# = 1
   FOR ROW = 1 TO n
      ' A singular matrix will have a zero element in the diagonal
      Sum# = Sum# * A(Order(ROW), ROW)
   NEXT ROW

   IF Sum# = 0 THEN
      iSing = TRUE
      EXIT SUB
   END IF
RETURN

BackSub:
   ' First, solve for X(n)
   X(n) = B(Order(n)) / A(Order(n), n)

   ' Now backsubstitute through the rest of the rows

   FOR ROW = n - 1 TO 1 STEP -1
      ' Calculate summation term
      Sum# = 0
      FOR COL = ROW + 1 TO n
         Sum# = Sum# + A(Order(ROW), COL) * X(COL)
      NEXT COL

      X(ROW) = (B(Order(ROW)) - Sum#) / A(Order(ROW), ROW)
   NEXT ROW
RETURN
END SUB

SUB PivotMatrix (A() AS DOUBLE, Order() AS INTEGER, P%, n%)

Pivot = P
MaxVal# = ABS(A(Order(P), P))

FOR ROW = P + 1 TO n
   Check# = ABS(A(Order(ROW), P))
   IF Check# > MaxVal# THEN
      MaxVal# = Check#
      Pivot = ROW
   END IF
NEXT ROW

SWAP Order(Pivot), Order(P)

END SUB


' $INCLUDE: 'E:\Users\BLM_Devl\common1.bi'
' $INCLUDE: 'E:\Users\BLM_Devl\common2.bi'
' $INCLUDE: 'E:\Users\BLM_Devl\constant.bi'

DEFINT A-Z
SUB PrintResidual (Iter%, MaxError!, MaxComp%, iSing%)

LOCATE 9, 43
PRINT "Max Error"
LOCATE 10, 43
PRINT "Component"

LOCATE 9, 54
'PRINT USING "##.##^^^^"; MaxError!
'PRINT FORMAT$(MaxError!, "##.##^^^^")
'PRINT FORMAT$(MaxError!, "#.###E+##")
PRINT FORMAT$(MaxError!, "0.000E+##")

LOCATE 9, 54
'CALL QPrint0(M$, -1)

M$ = SPACE$(7)
MID$(M$, 1) = Label$(MaxComp%)

LOCATE 10, 54
PRINT M$

T$ = "Iteration " + STR$(Iter%)
IF iSing% THEN
   T$ = T$ + " S"
END IF

l$ = SPACE$(18)
MID$(l$, 1) = T$

LOCATE 11, 43
PRINT l$

END SUB

FUNCTION Exist(DirFileName AS STRING) AS LONG
   FUNCTION = GETATTR(DirFileName$): FUNCTION = (ERR = 0)
END FUNCTION

SUB FileClose

CLOSE
'
' We need to trap an error here that can go all the way
' back to PBMAIN and end the program
'END

END SUB

DEFINT I-N
SUB ISortI (iArray(), Index(), NEls, iDir)
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  | SUB ISortI                                                            |
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
'  |    iArray()  - integer -   Array dimensioned from at least 1 to NEls  |
'  |                            contains the values to be sorted           |
'  |                                                                       |
'  | Output:                                                               |
'  |    Index()   - integer -   Array that holds the sorted order of       |
'  |                            the values in iArray.  Note that this array|
'  |                            does not have to be initialized as it does |
'  |                            in the assembly version.                   |
'  |                                                                       |
'  |  Developed for CHESS ver 2.20, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

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

      IF iArray(Index(j)) > iArray(Index(k)) THEN
         SWAP Index(j), Index(k)
         InOrder = FALSE
      END IF
   NEXT i
LOOP UNTIL InOrder

END SUB


SUB Pmsg(Strng$)
   'DIM ibl AS STRING
   PRINT Strng$ + "   "
   'ibl$ = WAITKEY$
   DO UNTIL LEN(INKEY$)
   LOOP
END SUB



DEFSNG I-N
FUNCTION PeekBuf%
T$ = INKEY$
IF LEN(T$) THEN
   PeekBuf% = ASC(T$)
ELSE
   PeekBuf% = 0
END IF
END FUNCTION

SUB PrintMsg (TmpS$)
PRINT
PRINT TmpS$
PRINT
END SUB

FUNCTION QPTrim$ (TmpStr$)

QPTrim$ = RTRIM$(LTRIM$(TmpStr$))

END FUNCTION

DECLARE FUNCTION LOG10S! (X!)

'$INCLUDE: 'E:\Users\BLM_Devl\common1.bi'
'$INCLUDE: 'E:\Users\BLM_Devl\common2.bi'
'$INCLUDE: 'E:\Users\BLM_Devl\constant.bi'
'$INCLUDE: 'E:\Users\BLM_Devl\default.bi'

DEFINT A-Z
SUB CHESSLog (LogFileName$, iAppend, iClose, iLogHand, ParameterFile$)
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB CHESSLog                                                         |
'  |                                                                       |
'  |  Creates a log file summarizing the current simulation.               |
'  |                                                                       |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*
IF NOT (iLogHand) THEN
   iLogHand = FREEFILE
END IF

IF (iAppend) THEN
   OPEN LogFileName$ FOR APPEND AS #iLogHand
   PRINT #iLogHand,
ELSE
   OPEN LogFileName$ FOR OUTPUT AS #iLogHand
END IF
DivLine$ = STRING$(80, "-")
LabelMask$ = " \  \"


FOR i = 1 TO NSpecies
   PRINT #iLogHand, Label(i)
NEXT i
FOR j = 1 TO NSpecies
   FOR i = 1 TO NComp
      T$ = LTRIM$(RTRIM$(STR$(SCb(j, i))))
      'PRINT #iLogHand, USING LabelMask$; T$;
      PRINT #iLogHand, FormatStr$(T$, LabelMask$);
   NEXT i
   PRINT #iLogHand,
NEXT j


T$ = "Simulation Type: "
SELECT CASE SimType
   CASE stBatch
      T$ = "Batch"
   CASE stTitration
      T$ = "Titration"
   CASE stColumn
      T$ = "Column"
END SELECT

PRINT #iLogHand, "CHESS ver 3.00 Log for " + T$ + " simulation described in " + ParameterFile$
PRINT #iLogHand, DivLine$

PRINT #iLogHand, "Component";
PRINT #iLogHand, TAB(15); "Species Type";
PRINT #iLogHand, TAB(30); "Activity";
PRINT #iLogHand, TAB(45); "Site Density";
PRINT #iLogHand,
FOR i = 1 TO NComp
   PRINT #iLogHand, Label(i);
   SELECT CASE SType(i)
      CASE stAqueous
         T$ = "Aqueous"
      CASE stSurface
         T$ = "Surface"
      CASE stPrecipitate
         T$ = "Precipitate"
   END SELECT
   PRINT #iLogHand, TAB(15); T$;

   SELECT CASE ActCorr(i)
      CASE acNone
         T$ = "None"
      CASE acDavies
         T$ = "Davies"
      CASE acDebye
         T$ = "Debye-Huckel"
      CASE acMoleFraction
         T$ = "Mole Fraction"
      CASE acChargeFraction
         T$ = "Charge Fraction"
      CASE acWHAM
         T$ = "Donnan"
   END SELECT
   PRINT #iLogHand, TAB(30); T$;

   T$ = STR$(SiteDen!(i))
   PRINT #iLogHand, TAB(45); T$;

   PRINT #iLogHand,
NEXT i


PRINT #iLogHand,
PRINT #iLogHand, "Species Formation Reactions";
PRINT #iLogHand, TAB(65); "Log K"
PRINT #iLogHand, DivLine$
P$ = " + "
M$ = " - "
e$ = " = "
FOR i = NComp + 1 TO NSpecies
   T$ = Label$(i) + e$
   NN = CompList(i, 0)
   FOR jj = 1 TO NN
      j = CompList(i, jj)
      SELECT CASE SCa(i, j)
         CASE 0
            ' do nothing, but don't let these go to case else
         CASE 1
            T$ = T$ + P$ + Label$(j)
         CASE -1
            T$ = T$ + M$ + Label$(j)
         CASE ELSE
            IF (SCa(i, j) > 0) THEN
               T$ = T$ + P$ + STR$(SCa(i, j)) + "*" + Label(j)
            ELSE
               T$ = T$ + M$ + STR$(ABS(SCa(i, j))) + "*" + Label(j)
            END IF
      END SELECT
   NEXT jj
   PRINT #iLogHand, T$;
   T$ = STR$(LOG10S(EC(i)))
   PRINT #iLogHand, TAB(65); T$
NEXT i


PRINT #iLogHand,
PRINT #iLogHand, "Component Mole Balance"
PRINT #iLogHand, DivLine$

'FOR ii = 1 TO NComp
FOR ii = 1 TO NUsed
   i = CompOrder(ii)
   T$ = "T." + Label(i) + " ="
   Ns = SpecList(ii, 0)
   FOR j = 1 TO Ns
      s = SpecList(ii, j)
      IF (tSCb(s, i) > 0) THEN
        T2$ = " + "
      ELSEIF (tSCb(s, i) < 0) THEN
        T2$ = " - "
      END IF
      IF (ABS(tSCb(s, i)) > 1) THEN
         T2$ = T2$ + LTRIM$(RTRIM$(STR$(tSCb(s, i))))
      END IF
      T2$ = T2$ + Label(s)
      T$ = T$ + T2$
   NEXT j
   PRINT #iLogHand, T$
NEXT ii

IF (iClose) THEN
   CLOSE #iLogHand
END IF

END SUB

FUNCTION FormatStr$(Strng$, StrMask$)
   DIM TmpS AS STRING
   DIM i AS INTEGER
   i = LEN(StrMask$)
   IF i=0 THEN
      IF LEFT$(Strng$, 1) = Quote$ AND RIGHT$(Strng$, 1) = Quote$ THEN
         TmpS$ = Strng$ + Comma$
      ELSE
         TmpS$ = Quote$ + Strng$ + Quote$ + Comma$
      END IF
   ELSE
     TmpS$ = SPACE$(i)
     'LSET TmpS$ = Strng$
     RSET TmpS$ = Strng$
   END IF
   FUNCTION = TmpS$
END FUNCTION

FUNCTION AddSpecies% (wtM1%, wtM2%, wtM3%)
'
'  Function to test a species in the WHAM database to
'  determine whether all components required are
'  included in this simulation.
'
'  Input:
'    wtM1%, wtM2%, wtM3%   WHAM database identifiers for
'                          components needed for this species
'
'  Returns:
'    AddSpecies = TRUE     all components are present
'    AddSpecies = FALSE    one or more components are missing
'
TEST1% = FALSE
TEST2% = FALSE
TEST3% = FALSE
FOR i = 1 TO NWHAM%
'   IF (WHAMID%(i) = wtM1%) OR (wtM1% = 0) THEN TEST1% = TRUE
'   IF (WHAMID%(i) = wtM2%) OR (wtM2% = 0) THEN TEST2% = TRUE
'   IF (WHAMID%(i) = wtM3%) OR (wtM3% = 0) THEN TEST3% = TRUE
   IF (WHAMID(i) = wtM1%) OR (wtM1% = 0) THEN TEST1% = TRUE
   IF (WHAMID(i) = wtM2%) OR (wtM2% = 0) THEN TEST2% = TRUE
   IF (WHAMID(i) = wtM3%) OR (wtM3% = 0) THEN TEST3% = TRUE
NEXT i

IF (TEST1% AND TEST2% AND TEST3%) THEN
   AddSpecies = TRUE
ELSE
   AddSpecies = FALSE
END IF
'AddSpecies = FALSE
END FUNCTION

DEFINT A-Z
SUB CalcActivityCoef (ISTR AS DOUBLE, ErrorCode AS INTEGER, ErrorMsg AS STRING)
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB CalcActivityCoef                                                 |
'  |                                                                       |
'  |  Calculates activity coefficients for all species.                    |
'  |                                                                       |
'  |  Subroutine Parameters:                                               |
'  |                                                                       |
'  |     ErrorCode% (integer, output) Non-zero value indicates an error    |
'  |                                  occured                              |
'  |                                                                       |
'  |     ErrorMsg$  (string, output)  Text of an error msg for reporting   |
'  |                                                                       |
'  |     IStr       (double, input)   Ionic strength, (moles/liter)        |
'  |                                                                       |
'  |  COMMON SHARED Variables                                              |
'  |                                                                       |
'  |     ActCorr()  (integer, input)  Specifies method used to calculate   |
'  |                                  activity correction                  |
'  |                                                                       |
'  |     OC()       (double, input)   Species concentrations (moles/liter) |
'  |                                                                       |
'  |     V()        (double, input)   Ionic charge                         |
'  |                                                                       |
'  |     ActCoeff() (double, output)  Values of Activity Coefficients      |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*


' *------------------------------------------------------------*
' |  Simple and Array Variables Local to SUB CalcActivityCoef  |
' *------------------------------------------------------------*
DIM GAMMA AS DOUBLE         '  Temporary variable used to calculate coeff.
DIM i AS INTEGER            '  Array/loop index variable
DIM j AS INTEGER            '  Array/loop index variable
DIM k AS INTEGER            '  Array index variable
DIM LLDenom(NList) AS DOUBLE'  Denomenator of linked list calculations
DIM S AS INTEGER            '  Array/loop index variable
DIM TmpS AS STRING          '  Temporary string variable
DIM ACTA AS DOUBLE
DIM ACTB AS DOUBLE
DIM ActC AS DOUBLE
DIM ACTD AS DOUBLE

' *------------------------------*
' |  Initialize Error variables  |
' *------------------------------*
ErrorCode = FALSE
ErrorMsg$ = NullStr

' *-------------------------------------*
' | Calculate and store denomenator of  |
' | mole or charge fraction corrections |
' | in variable LLDenom()               |
' *-------------------------------------*
FOR i = 1 TO NList
   SELECT CASE ActCorr(LinkList(i, 1))
      CASE 4
         ' *-------------------*
         ' |  mole fraction    |
         ' *-------------------*
         LLDenom(i) = ZeroD
         FOR j = 1 TO LinkList(i, 0)
            LLDenom(i) = LLDenom(i) + OC(LinkList(i, j))
         NEXT j
      CASE 5
         ' *-------------------*
         ' |  charge fraction  |
         ' *-------------------*
         LLDenom(i) = ZeroD
         FOR j = 1 TO LinkList(i, 0)
            LLDenom(i) = LLDenom(i) + OC(LinkList(i, j)) * V(LinkList(i, j))
         NEXT j
      CASE ELSE
         ' *---------------------*
         ' |  Error notification |
         ' *---------------------*
         ErrorMsg$ = "The species " + Label$(LinkList(i, 1)) + " has an inappropriate activity" + CHR$(13)
         ErrorMsg$ = ErrorMsg$ + "correction method specified."
         ErrorCode = TRUE
         EXIT SUB

   END SELECT
   ' *-------------------------------------------------------------*
   ' |  Assume that activity coefficients of zero result from a    |
   ' |  first pass through the model when all concentrations are   |
   ' |  far from equilibrium.  Set equal to 1 to avoid errors.     |
   ' *-------------------------------------------------------------*
   IF LLDenom(i) = 0 THEN LLDenom(i) = 1#
NEXT i

' *------------------------------------------------*
' |  Store the activity coefficients for all ions  |
' *------------------------------------------------*
FOR S = 1 TO NSpecies
   SELECT CASE ActCorr(S)
      ' *--------------------------------*
      ' |  Key for ActCorr assignments:  |
      ' |  1 = none                      |
      ' |  2 = Davies                    |
      ' |  3 = Debye-Hu�ckel              |
      ' |  4 = mole fraction             |
      ' |  5 = charge fraction           |
      ' *--------------------------------*
      CASE acNone
         ' *---------------------------------*
         ' | No activity correction is used, |
         ' | activity coeff = 1              |
         ' *---------------------------------*
         GAMMA# = 1

      CASE acDavies
         ' *------------------------------------*
         ' |  Davies equation can be applied to |
         ' |  charged aqueous phase ions        |
         ' *------------------------------------*
         GAMMA# = Davies#(ISTR#, V(S))

      CASE acDebye
         ' *------------------------------------*
         ' |  Debye-Hu�ckel not yet implemented  |
         ' *------------------------------------*
         'ErrorMsg$ = "Debye-Hu�ckel was specified for species " + Label$(S) + CHR$(13)
         'ErrorMsg$ = ErrorMsg$ + "but has not been implemented in this version of CHESS"
         'ErrorCode = TRUE
         'EXIT SUB
         GAMMA# = DebyeHuckel#(ISTR#, V(S))

      CASE acMoleFraction
         ' *----------------*
         ' |  mole fraction |
         ' *----------------*
         GAMMA# = LLDenom(ListPos(LinkComp(S))) ^ -1

      CASE acChargeFraction
         ' *------------------*
         ' |  charge fraction |
         ' *------------------*
         k = LinkComp(S)
         GAMMA# = V(S) / LLDenom(ListPos(k))

      CASE acWHAM
         GAMMA# = 1#

      CASE ELSE
         ' *---------------------*
         ' |  Error notification |
         ' *---------------------*
         ErrorMsg$ = "Unknown activity correction method was specified for species " + Label$(S)
         ErrorCode = TRUE
         EXIT SUB
   END SELECT
   ActCoef(S) = GAMMA#
NEXT S

' *-------------------------------------------*
' |  Free up memory assigned to local arrays  |
' *-------------------------------------------*
ERASE LLDenom

END SUB

DEFSNG A-Z
FUNCTION CalcEN#

ok% = FALSE
count% = 0
DO
   count% = count% + 1
   IStr1# = ISTR#
   CALL ChargeBalance(CHARGE#, ISTR#)
   CALL CalcSpeciesConc(ISTR#, ErrorCode%, ErrorMsg$)
   IF (ABS((ISTR# - IStr1#) / ISTR#) < .001) THEN
      ok% = TRUE
   END IF
LOOP UNTIL ok% OR (count% > 5)

CalcEN# = CHARGE#

END FUNCTION

DEFINT A-Z
SUB CalcMassChange
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB CalcMassChange                                                   |
'  |                                                                       |
'  |  Calculates the mass (moles) of a solid or gaseous phase required     |
'  |  to precipitate or dissolve in order to reach saturation.             |
'  |                                                                       |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

' *-----------------------------------------------*
' |  Simple Variables Local to SUB DefineProblem  |
' *-----------------------------------------------*
DIM A1 AS INTEGER           '  Array index variable
DIM a2 AS INTEGER           '  Array index variable
DIM p1 AS INTEGER           '  Array/loop index variable
DIM p2 AS INTEGER           '  Array/loop index variable
DIM q1 AS INTEGER           '  Array/loop index variable
DIM q2 AS INTEGER           '  Array/loop index variable
DIM S1 AS INTEGER           '  Array/loop index variable
DIM Sum AS DOUBLE           '  Temporary variable
DIM Temp AS DOUBLE          '  Temporary variable

FOR q1 = 1 TO NPhase
   p1 = PhaseOrder(q1)
   IF PType(p1) = ptSaturated THEN
      ' *------------------------------------------------*
      ' | Estimate the solubility of this phase based on |
      ' | current species concentration only             |
      ' *------------------------------------------------*
      A1 = SComp(p1)
      Sum# = OC(A1) * CtoM!(A1)
      FOR S1 = NFirstS TO NSpecies
         Sum# = Sum# + OC(S1) * SCb(S1, A1) * CtoM!(S1)
      NEXT S1

      ' *----------------------------------------------------------------*
      ' |  Now add contributions from other dissolving, saturated phases |
      ' |  (according to design, they should come before this phase in   |
      ' |  PhaseOrder() since the component for this phase may have been |
      ' |  substituted out of all other expressions)                     |
      ' *----------------------------------------------------------------*
      FOR q2 = 1 TO q1 - 1
         p2 = PhaseOrder(q2)
         a2 = SComp(p2)
         Sum# = Sum# + ChangeMass(p2) * SP(p2, A1) / SP(p2, a2)
      NEXT q2

      ChangeMass(p1) = (KT(A1) * SiteDen(A1) - Sum#) / SP(p1, A1)

   ELSEIF PType(p1) = ptUnSaturated THEN
      ' *----------------------------------------*
      ' |  Everything that's here will dissolve  |
      ' *----------------------------------------*
      ChangeMass(p1) = -PhaseMass(p1)

   ELSEIF PType(p1) = ptNotAvailable THEN
      ' *---------------------------------*
      ' |  No amount of phase is present, |
      ' |  and the system is unsaturated  |
      ' *---------------------------------*
      ChangeMass(p1) = 0

   END IF
NEXT q1

' *---------------------------------------------------------------*
' |  Now go back and add the contribution of Type 2 phases to all |
' |  other phase mass changes.  This has to be done after the     |
' |  fact since Type 2 phases don't occur in any particular       |
' |  order in the list of phases.                                 |
' *---------------------------------------------------------------*
FOR p1 = 1 TO NPhase
   IF PType(p1) = ptUnsaturated THEN '2 THEN
      A1 = SComp(p1)
      FOR p2 = 1 TO NPhase
         IF PType(p2) = ptSaturated THEN '1 THEN
            Temp# = SP(p2, A1) * PhaseMass(p1) / SP(p1, A1)
            ChangeMass(p2) = ChangeMass(p2) + Temp#
         END IF
      NEXT p2
   END IF
NEXT p1

END SUB

SUB CalcResidual (MaxError!, MaxComp%)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB CalcResidual                                                     |
'  |                                                                       |
'  |  The residual function is evaluated here for all components.          |
'  |  This is calculated as the difference between the analytically        |
'  |  known total conc. and the total quantity calculated as a result      |
'  |  of mass balance relationships.                                       |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

' *-----------------------------------------------*
' |  Simple Variables Local to SUB DefineProblem  |
' *-----------------------------------------------*
DIM c AS INTEGER            '  Array/loop index variable
DIM c2 AS INTEGER           '  Array/loop index variable
DIM cc AS INTEGER           '  Array/loop index variable
DIM CT AS DOUBLE            '  Calculated Total conc.
DIM Ns AS INTEGER           '  Number of species in a given sum
DIM MaxElement AS DOUBLE    '  Temp. variable to store largest element of CT
DIM S AS INTEGER            '  Array/loop index variable
DIM ss AS INTEGER           '  Array/loop index variable
DIM Temp AS DOUBLE          '  Temporary variable
DIM Factors AS SINGLE
DIM EN AS DOUBLE
DIM BestEN AS DOUBLE
DIM BestKT AS DOUBLE
DIM Bracketed AS INTEGER

   ' *----------------------------------------------*
   ' |  First, calculate the change in mass due to  |
   ' |  dissolution or precipitation of phases.     |
   ' *----------------------------------------------*
   CalcMassChange

   MaxError! = 0
   FOR cc = 1 TO NUsed
      c = CompOrder(cc)

      ' *----------------------------------------*
      ' |  Add the contribution of all species   |
      ' |  to the mass balance of this component |
      ' *----------------------------------------*
      CT# = 0#
      MaxElement# = CT#
      Ns = SpecList(cc, 0)
      FOR ss = 1 TO Ns
         S = SpecList(cc, ss)
         Temp# = OC(S) * tSCb(S, c) * CtoM!(S)
         CT# = CT# + Temp#
         MaxElement# = MaxElement# + ABS(Temp#)
      NEXT ss

      IF (CType(c) = ctChargeBal) THEN
         ' *---------------------------------------------------*
         ' |  Find the value of OC(c) to result in a solution  |
         ' |  chemistry that is electrially neutral            |
         ' *---------------------------------------------------*

         '  Get the initial value of electroneutrality
         ENa# = CalcEN#

'        iDebugOld% = FREEFILE
'        OPEN "debug.prn" FOR OUTPUT AS #iDebugOld
'        PRINT #iDebugOld, CDBL(OC(c)), CDBL(ENa#)
'
'        '  Find a value of OC(c) that brackets EN = 0
'        Bracketed = FALSE
'        Dir% = 1
'        count% = 0
'        ChangedDir = FALSE
'        OCa# = OC(c)
'        oldEN# = ENa#
'        DO
'           count% = count% + 1
'           'OC(c) = OC(c) * (1 - .1 * Sign%(CDBL(V(c))) * Sign%(ENb#))
'           OC(c) = OC(c) * (1 - .1 * Dir%)
'           ENb# = CalcEN#
'           IF (ABS(ENb#) > ABS(oldEN#)) AND NOT (count% < oldc% + 3) THEN
'              oldc% = count%
'              oldEN# = ENb#
'              ChangedDir = TRUE
'              Dir% = Dir% * -1
'           END IF
'           PRINT #iDebugOld, USING "##.####^^^^ "; CSNG(OC(c)); CSNG(ENb#); Dir%
'           IF (Sign%(ENa#) <> Sign%(ENb#)) THEN Bracketed = TRUE
'        LOOP UNTIL Bracketed
'
'        CLOSE #iDebugOld
'        OCb# = OC(c)
'
'        ' *----------------------------------------------------*
'        ' |  Ok, now the root is bracketed, make sure point a  |
'        ' |  corresponds to positive EN                        |
'        ' *----------------------------------------------------*
'        IF (ENa# < 0) THEN
'           SWAP OCa#, OCb#
'           SWAP ENa#, ENb#
'        END IF
'
'        ' *------------------------------------------*
'        ' |  Now bisect the interval until we close  |
'        ' *------------------------------------------*
'        DO
'           OC(c) = (OCa# + OCb#) / 2
'           EN# = CalcEN#
'           IF (EN# > 0) THEN
'              OCa# = OC(c)
'           ELSE
'              OCb# = OC(c)
'           END IF
'
'        'LOOP UNTIL (ABS(EN# / OC(NSpecies + 2)) < .001)
'        LOOP UNTIL (((OCa# - OCb#) / OCa#) < .001)
      END IF

      ' *------------------------------------------*
      ' |  Sum up the component totals relevant to |
      ' |  this residual created from solubility   |
      ' |  product substitutions                   |
      ' *------------------------------------------*
      Temp# = 0#
      FOR c2 = 1 TO NComp
         S = SpecList(cc, c2)
         IF S > NComp THEN EXIT FOR
         Temp# = Temp# + KT(S) * SiteDen(S) * tSCb(S, c)
      NEXT c2

      SELECT CASE CType(c)
         CASE ctFloat
            R(cc) = CT# - Temp#
         CASE ctChargeBal
            KT(c) = CT#
            R(cc) = ENa#
         CASE ctFixed
            KT(c) = CT#
            R(cc) = ZeroD
         CASE ctSurfPot
            KT(c) = LOG(OC(c)) * Lambda# * SQR(OC(NSpecies + 1)) '/ CtoM!(c)
            R(cc) = CT# - KT(c)
         CASE ELSE
      END SELECT

      ' *--------------------------------------------*
      ' |  Store the maximum error so that it can be |
      ' |  checked against the convergence criteria. |
      ' *--------------------------------------------*
      IF MaxElement# <> 0 AND cc <= NSolve THEN
         IF (CType(c) = ctChargeBal) THEN
            Temp# = ABS(R(cc) / OC(NSpecies + TWO))
         ELSE
            Temp# = ABS(R(cc) / MaxElement#)
         END IF
         IF Temp# > MaxError! THEN
            ' *--------------------------------*
            ' | Avoid overflows for big errors |
            ' *--------------------------------*
            IF Temp# < 1E+33 THEN
               MaxError! = CSNG(Temp#)
               MaxComp% = c
            ELSE
               MaxError! = 1E+33
               MaxComp% = c
            END IF
         END IF
      END IF
   NEXT cc
    IF (iDebug = 01 OR iDebug = 06) THEN
    '    FOR s = NFirstS TO NSpecies
        FOR cc = 1 TO NUsed
            c = CompOrder(cc)
           CALL DebugLog ("CRes", Label(c) + STR$(OC(c)))
           CALL DebugLog ("CRes", "R." + Label(c) + STR$(R(cc)))
           CALL DebugLog ("CRes", "T." + Label(c) + STR$(KT(c)))
        NEXT cc
    END IF

END SUB

SUB CalcSpeciesConc (ISTR#, ErrorCode%, ErrorMsg$)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB CalcSpeciesConc                                                  |
'  |                                                                       |
'  |  The equilibrium concentration of all species is calculated from      |
'  |  the current values of component concentrations.                      |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

' *-----------------------------------------------*
' |  Simple Variables Local to SUB ChargeBalance  |
' *-----------------------------------------------*
DIM Activity AS DOUBLE      '  Temporary variable
DIM c AS INTEGER            '  Array/Loop index variable
DIM cc AS INTEGER           '  Array/Loop index variable
DIM CHARGE AS DOUBLE        '  Net charge balance
DIM N AS INTEGER            '  End of loop counter
DIM P AS INTEGER            '  Array/Loop index variable
DIM q AS INTEGER            '  Array/Loop index variable
DIM S AS INTEGER            '  Array/Loop index variable
DIM Temp AS DOUBLE          '  Temporary variable

CALL ChargeBalance(CHARGE#, ISTR#)
IF (iDebug = 01 OR iDebug = 04) THEN CALL DebugLog ("CSpe", "After Charge Balance, Charge " + STR$(CHARGE#) + ", ISTR " + STR$(ISTR#))

CALL CalcActivityCoef(ISTR#, ErrorCode%, ErrorMsg$)
IF (iDebug = 01 OR iDebug = 04) THEN CALL DebugLog ("CSpe", "After Charge Balance, ISTR " + STR$(ISTR#))

IF ErrorCode% THEN
   EXIT SUB
END IF

' *-----------------------------------------------------------*
' |  For any component replaced by a solubility expression,   |
' |  calculate it's concentration based on solubility control |
' *-----------------------------------------------------------*
FOR q = 1 TO NPhase
   P = PSubOrder(q)
   IF PType(P) = ptSaturated THEN
       S = SComp(P)
       Temp# = 1
       N = CompList(S, 0)
       FOR cc = 1 TO N
         c = CompList(S, cc)
         IF c <> S THEN
            Activity# = OC(c) * ActCoef(c)
            Temp# = Temp# * Activity# ^ tSCa(S, c)
         END IF
       NEXT cc
       OC(S) = Temp# * tEC(S) / ActCoef(S)
   END IF
NEXT q

' *-------------------------------------------------------------------*
' |  Calculate the concentrations of all species given the reactions  |
' |  and equilibrium constants that are: 1) corrected for temp, and   |
' |  2) transformed by solubility relationships                       |
' *-------------------------------------------------------------------*
FOR S = NFirstS TO NSpecies
 IF NOT (WSpec(S)) THEN
   Temp# = 1
   N = CompList(S, 0)
   FOR cc = 1 TO N
      c = CompList(S, cc)
      Activity# = OC(c) * ActCoef(c)
      Temp# = Temp# * Activity# ^ tSCa(S, c)
   NEXT cc
   OC(S) = Temp# * tEC(S) / ActCoef(S)
   IF (iDebug = 01 OR iDebug = 04) THEN CALL DebugLog ("CSpe", Label(s) + STR$(OC(s)))

 END IF
NEXT S

'CALL WHAMCalcSpecies
IF (NOT iDoneWHAM) AND (iWHAM) THEN
   CALL WHAMIter
   iDoneWHAM = TRUE
END IF

IF FALSE THEN
    IF (iAl <> 0) AND (iAlp <> 0) THEN
        IF (OKT(iAl) > KT(iAl)) THEN
            KT(iAlp) = OKT(iAl) - KT(iAl)
        ELSE
            KT(iAlp) = SmallNumber
        END IF
        'OC(iAlp) = KT(iAlp)
    END IF
END IF

IF (iDebug = 01 OR iDebug = 04) THEN
'    FOR s = NFirstS TO NSpecies
    FOR s = 1 TO NSpecies
       CALL DebugLog ("CSpe", Label(s) + STR$(OC(s)))
    NEXT s
END IF

END SUB

SUB ChargeBalance (CHARGE#, ISTR#)
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB ChargeBalance                                                    |
'  |                                                                       |
'  |  Calculates the sum of ion equivalents and net charge for             |
'  |  all dissolved  species.  Units are (equivalents / Liter)             |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

' *-----------------------------------------------*
' |  Simple Variables Local to SUB ChargeBalance  |
' *-----------------------------------------------*
DIM S AS INTEGER            '  Array/Loop index variable
DIM T1 AS SINGLE            '  Temporary single precision variable
DIM T2 AS DOUBLE            '  Temporary double precision variable
DIM SumCations AS DOUBLE    '  Sum of aqueous cation charge
DIM SumAnions AS DOUBLE     '  Sum of aqueous anion charge

CHARGE# = 0
ISTR# = 0
SumCations# = 0
SumAnions# = 0

' *---------------------------------------------*
' |  Sum ionic charge over all aqueous species  |
' *---------------------------------------------*
FOR S = 1 TO NSpecies
   IF (SType(S) = stAqueous) THEN
      T1! = V(S)
      T2# = OC(S) * T1!
      CHARGE# = CHARGE# + T2#
      IF (T1! < 0) THEN
         SumAnions# = SumAnions# + T2#
      ELSE
         SumCations# = SumCations# + T2#
      END IF
      ISTR# = ISTR# + T1! * T2#
   END IF
NEXT S
ISTR# = ISTR# * .5

' *-----------------------------*
' |  Store charge balance info  |
' *-----------------------------*
OC(NSpecies + ONE) = CHARGE#
OC(NSpecies + TWO) = ISTR#
OC(NSpecies + THREE) = SumCations#
OC(NSpecies + FOUR) = SumAnions#

END SUB

SUB CheckSolidPhase (Converge%, UnderFlag%, OverFlag%)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB CheckSolidPhase                                                  |
'  |                                                                       |
'  |  Checks for oversaturation or undersaturation in the current list     |
'  |  of phases that could change their definition and lead to a new       |
'  |  equilibrium problem.                                                 |
'  |                                                                       |
'  |  Type ptSaturated phases are present in sufficient quantities to      |
'  |  saturate the aqueous system at equilibrium.  All type ptSaturated    |
'  |  phases redefine a ctFloat type component as ctSubstituted.           |
'  |                                                                       |
'  |  Type ptUnSaturated phases are present, but in insufficient           |
'  |  quantities to saturate the aqueous phase.                            |
'  |                                                                       |
'  |  Type ptNotAvailable phases are not present in any quantity.          |
'  |                                                                       |
'  |  A phase defined as type ptSaturated with insufficient mass to        |
'  |  saturate the aqueous system is redefined as type ptUnSaturated.      |
'  |  A ptSaturated phase can be defined with infinite mass available      |
'  |  (i.e., PhaseMass() = ptUnlimited) and then it can never become       |
'  |  a type ptUnSaturated or ptNotAvailable phase.                        |
'  |                                                                       |
'  |  A phase defined as type NotAvailable may precipitate if the aqueous  |
'  |  system becomes oversaturated.  Then it will be redefined as a type   |
'  |  ptSaturated phase.                                                   |
'  |                                                                       |
'  |  Since a ptUnSaturated phase can dissolve and is not present in       |
'  |  sufficient quantities to saturate the aqueous system, it must become |
'  |  a type ptNotAvailable.  However, that redefinition will have no      |
'  |  effect on the current problem and must be handled by the calling     |
'  |  program.                                                             |
'  |                                                                       |
'  |                                                                       |
'  |  Subroutine Parameters:                                               |
'  |                                                                       |
'  |     Converge   (integer, output)  Boolean flag to indicate to calling |
'  |                                   program that the current system has |
'  |                                   not yet been solved.                |
'  |                                                                       |
'  |                                   A redefinition of any phase will    |
'  |                                   result in Converge = FALSE          |
'  |                                                                       |
'  |     UnderFlag  (integer, output)  Boolean flag to indicate that one   |
'  |                                   or more phases is undersaturated.   |
'  |                                                                       |
'  |     OverFlag   (integer, output)  Boolean flag to indicated that one  |
'  |                                   or more phases is oversaturated.    |
'  |                                                                       |
'  |  Developed for CHESS ver 3.00, May 1995                               |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |  Dept. Civil and Environmental Eng, Syracuse University               |
'  |  rsantore@mailbox.syr.edu                                             |
'  |                                                                       |
'  *-----------------------------------------------------------------------*

' *-----------------------------------------------*
' |  Simple Variables Local to SUB DefineProblem  |
' *-----------------------------------------------*
DIM P AS INTEGER            '  Array/loop index variable
DIM Temp AS DOUBLE          '  Temporary variable

FOR P = 1 TO NPhase

   SELECT CASE PType(P)

      CASE ptSaturated
         ' *----------------------------------*
         ' |  Phase is solubility controlled  |
         ' *----------------------------------*
         IF PhaseMass(P) <> ptUnlimited THEN
            ' *-----------------------------------------*
            ' |  Only consider phases with finite mass  |
            ' *-----------------------------------------*
            Temp# = IAP(P) / SEC(P)
            IF (Temp# - 1!) > .001 THEN
               ' *----------------------------------------------------*
               ' |  This phase is supposed to be saturated but isn't  |
               ' |  Don't change it's status, but make another call   |
               ' |  to Saturate to update solubility constants        |
               ' |  Another iteration can be accomplished by simply   |
               ' |  setting the Converge flag to FALSE.               |
               ' *----------------------------------------------------*
               Converge = FALSE

            ELSEIF PhaseMass(P) + ChangeMass(P) < 0 THEN
               ' *------------------------------------------*
               ' |  More phase is required than is present  |
               ' *------------------------------------------*
               Converge = FALSE
               UnderFlag = TRUE
               IF PhaseMass(P) > 0 THEN
                  ' *--------------------------------------------------------*
                  ' |  If some phase is present, move to type ptUnSaturated  |
                  ' *--------------------------------------------------------*
                  PType(P) = ptUnSaturated

               ELSE
                  ' *-------------------------------------------------------*
                  ' |  If no phase is present, move to type ptNotAvailable  |
                  ' *-------------------------------------------------------*
                  PType(P) = ptNotAvailable
                  PType(P) = ptUnSaturated
               END IF
               CType(SComp(P)) = ctFloat
            END IF
         END IF

      CASE ptUnSaturated OR ptNotAvailable
         ' *----- -----------------------------------------------------*
         ' |  Case ptUnSaturated:                                      |
         ' |     Phase should all dissolve and still be undersaturated |
         ' |                                                           |
         ' |  Case ptNotAvailable:                                     |
         ' |     Phase is not present and aqueous system should be     |
         ' |     undersaturated.                                       |
         ' *-----------------------------------------------------------*
         Temp# = IAP(P) / SEC(P)

         IF Temp# > 1 THEN
            ' *-----------------------------------------------------*
            ' |  Phase is oversaturated.  Move to Type ptSaturated  |
            ' *-----------------------------------------------------*
            Converge = FALSE
            OverFlag = TRUE
            PType(P) = ptSaturated
            CType(SComp(P)) = ctSubstituted
         END IF
   END SELECT

NEXT P

END SUB

DEFSNG A-Z

SUB CSFIT1 (N AS INTEGER, X() AS DOUBLE, Y() AS DOUBLE, XVAL AS DOUBLE, FVAL AS DOUBLE) STATIC

    ' Natural cubic spline interpolation subroutine

    ' Input

    '  N    = number of X and Y data points (N >= 2)
    '  X()  = vector of X data ( N rows )
    '  Y()  = vector of Y data ( N rows )
    '  XVAL = X argument

    ' Output

    '  FVAL = function value at XVAL

    DIM S(N) AS DOUBLE
    DIM G(N - 1) AS DOUBLE
    DIM WORK(N - 1) AS DOUBLE
    DIM I AS INTEGER
    DIM XI   AS DOUBLE
    DIM XIM1 AS DOUBLE
    DIM XIP1 AS DOUBLE
    DIM YI   AS DOUBLE
    DIM YIM1 AS DOUBLE
    DIM YIP1 AS DOUBLE
    DIM X    AS DOUBLE
    DIM H    AS DOUBLE
    DIM T    AS DOUBLE
    DIM W    AS DOUBLE
    DIM U    AS DOUBLE
    DIM Z    AS DOUBLE
    DIM SS   AS DOUBLE

    FOR I = 2 TO N - 1
        XI = X(I)
        XIM1 = X(I - 1)
        XIP1 = X(I + 1)
        YI = Y(I)
        YIM1 = Y(I - 1)
        YIP1 = Y(I + 1)
        X = XI - XIM1
        H = XIP1 - XIM1
        WORK(I) = .5# * X / H
        T = ((YIP1 - YI) / (XIP1 - XI) - (YI - YIM1) / X) / H
        S(I) = 2# * T
        G(I) = 3# * T
    NEXT I

    S(1) = 0#
    S(N) = 0#

    W = 8# - 4# * SQR(3#)

    DO
       U = 0#

       FOR I = 2 TO N - 1
           T = W * (-S(I) - WORK(I) * S(I - 1) - (.5# - WORK(I)) * S(I + 1) + G(I))
           H = ABS(T)
           IF (H > U) THEN U = H
           S(I) = S(I) + T
       NEXT I
    LOOP UNTIL (U < .00000001#)

    FOR I = 1 TO N - 1
        G(I) = (S(I + 1) - S(I)) / (X(I + 1) - X(I))
    NEXT I

    I = 1

    DO
       I = I + 1
    LOOP UNTIL (XVAL <= X(I))

    I = I - 1

    H = XVAL - X(I)
    T = XVAL - X(I + 1)
    X = H * T
    SS = S(I) + H * G(I)
    Z = 1# / 6#
    U = Z * (S(I) + S(I + 1) + SS)
    W = (Y(I + 1) - Y(I)) / (X(I + 1) - X(I))

    FVAL = W * H + Y(I) + X * U

    ERASE S, G, WORK

END SUB

SUB POLINT (XA() AS DOUBLE, YA() AS DOUBLE, N AS INTEGER, X AS DOUBLE, Y AS DOUBLE, DY AS DOUBLE)
DIM C(N) AS LOCAL DOUBLE
DIM D(N) AS LOCAL DOUBLE
DIM NS AS INTEGER
DIM DIFF AS DOUBLE
DIM I AS INTEGER
DIM M AS INTEGER
DIM W AS DOUBLE
NS = 1
DIF = ABS(X - XA(1))
FOR I = 1 TO N
  DIFT = ABS(X - XA(I))
  IF DIFT < DIF THEN
    NS = I
    DIF = DIFT
  END IF
  C(I) = YA(I)
  D(I) = YA(I)
NEXT I
Y = YA(NS)
NS = NS - 1
FOR M = 1 TO N - 1
  FOR I = 1 TO N - M
    HO = XA(I) - X
    HP = XA(I + M) - X
    W = C(I + 1) - D(I)
    DEN = HO - HP
    IF DEN = 0! THEN
       'PRINT "Abnormal exit"
       Y = -999
       EXIT SUB
    END IF
    DEN = W / DEN
    D(I) = HP * DEN
    C(I) = HO * DEN
  NEXT I
  IF 2 * NS < N - M THEN
    DY = C(NS + 1)
  ELSE
    DY = D(NS)
    NS = NS - 1
  END IF
  Y = Y + DY
NEXT M
ERASE D, C
END SUB

SUB OutputJacobian
    FOR j=1 TO NSolve
        PRINT Label(CompOrder(j));" ";
    NEXT j
    PRINT
    FOR j=1 TO NSolve
        PRINT Label(CompOrder(j));" ";
        FOR k=1 TO NSolve
            PRINT FORMAT$(Z(j,k), "0.0E+##");" ";
        NEXT k
        PRINT
    NEXT j

    PRINT
    FOR j=1 TO NSolve
        PRINT Label(CompOrder(j));" ";
    NEXT j
    PRINT
    FOR j=1 TO NSolve
        PRINT Label(CompOrder(j));" ";
        FOR k=1 TO NSolve
            PRINT FORMAT$(Zn(j,k), "0.0E+##");" ";
        NEXT k
        PRINT
    NEXT j
END SUB
