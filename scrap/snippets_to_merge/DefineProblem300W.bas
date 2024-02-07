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

IF (MajorVer >=3) AND (MinorVer >= 10) THEN
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
    NMass = 2
    MassName(1) = "Water"                   ' The compartment name
    MassVal(1) = 1.0                        ' The mass value
    MassUnitLabel(1) = "L"                  ' The units of a mass compartment (e.g. L or kg)

    MassName(2) = "BL"                      ' The compartment name
    MassVal(2) = 1.0                        ' The mass value
    MassUnitLabel(2) = "kg wet"             ' The units of a mass compartment (e.g. L or kg)

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


'Start of WHAM shenanigans
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
      TmpS$ = PathNameblm$(FName$) + FirstName$(sWHAM$) + "." + ExtName$(sWHAM$)
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
