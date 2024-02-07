SUB DefineProblem300B (FName$, ErrorCode, ErrorMsg$)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB DefineProblem300B                                                |
'  |                                                                       |
'  |  Reads in a ver 3.00 file for use with the CHESS model.  Although a   |
'  |  large number of variables are initialized in this subroutine, most   |
'  |  are defined as global variables.  See the variable lists for         |
'  |  specific details.                                                    |
'  |                                                                       |
'  |  Subroutine Parameters:                                               |
'  |                                                                       |
'  |     FName$  (string, input)       Name of input file to read          |
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
DIM i AS INTEGER            '  Array/loop index variable
DIM j AS INTEGER            '  Array/loop index variable
DIM k AS INTEGER            '  Array/loop index variable
DIM hs_i AS INTEGER         '  Array/loop index variable for humic substances
DIM Nj AS INTEGER           '  Number of elements in a linked/list
DIM Tmp AS SINGLE           '  Temporary single precision variable
DIM TmpD AS DOUBLE          '  Temporary double precision variable
DIM Tmp2 AS SINGLE          '  Temporary single precision variable
DIM TmpS AS STRING          '  Temporary string variable
DIM TmpS2 AS STRING         '  Temporary string variable
DIM TmpI AS INTEGER         '  Temporary integer variable
DIM TmpI1 AS INTEGER        '  Temporary integer variable
DIM TmpI2 AS INTEGER        '  Temporary integer variable
DIM wtPKMHA(2) AS SINGLE    '  metal binding constant from line i of WHAM database
DIM FirstHS AS INTEGER   '  starting index of each humic substance protonated component
DIM Num_DL AS INTEGER       '  number of DL species
DIM wtM(3) AS INTEGER
DIM wtS(3) AS INTEGER

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
   IF (Label$(i) = "H") THEN iH = i

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
   Label$(i + MaxSpecies) = QPTrim$(TmpS$)

   INPUT #iFileIn, STempCorr(i)     'iTemp

   INPUT #iFileIn, SP(i, 0)         'iMonte - this will be the Monte Carlo flag for phases

   FOR c = 1 TO NComp
      INPUT #iFileIn, SP(i, c)      'Stoic
   NEXT c

   INPUT #iFileIn, Tmp!, Tmp2!      'LogK; Tmp2! is the MC Var LogK for phases
   SEC(i) = 10 ^ Tmp!

   INPUT #iFileIn, Tmp!, Tmp2!
   SDeltaH(i) = Tmp!                'DeltaH
   SECTemp(i) = Tmp2!               'Temp


   INPUT #iFileIn, Tmp!
   PhaseMass(i) = Tmp!              'Moles

   IF Tmp! = 0! THEN
      PType(i) = ptNotAvailable     'Converts "0" values to "3" (ptNotAvailable)
   ELSE
      PType(i) = ptSaturated        'all other values get converted to "1" (ptSaturated)
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
               NBL = NBL + 1
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
   ' Open WHAM database and read in data
   '
   '
   ' Check to see if the file is there
   IF NOT Exist(sWHAM$) THEN
      '
      ' database does not exist in the current directory,
      ' check the directory where the parameter file was found
      TmpS$ = PATHNAMEblm$(FName$) + FirstName$(sWHAM$) + "." + ExtName$(sWHAM$)
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

    iStartWHAM = NComp
    WOffSet = 21

    IF ((NComp + 2*WOffSet) > MaxComp) THEN
           ErrorMsg$ = "Too many components."
           ErrorCode% = TRUE
           EXIT SUB
    END IF

    FOR hs_i = 1 TO 2
        Sol_HS(hs_i) = 0

       IF (hs_i = 1) THEN
           hs$ = "HA"
       ELSEIF (hs_i = 2) THEN
           hs$ = "FA"
       END IF

       ' *-----------------------------------------------*
       ' |  Make room for WHAM components for HS_i (20)  |
       ' |  reactions by moving all species              |
       ' *-----------------------------------------------*
       'Shift over numbers in species arrays
       FOR i = 0 TO (NSpecies - NComp + NBL + hs_i - 2) '-1
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
          Label(k) = Label(j)
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

       'Next, move components
       FOR i = 0 TO NBL + hs_i - 2 '- 1
           j = NComp - i
           k = NComp - i + WOffSet
           CType(k) = CType(j)
           CompOrder(k) = CompOrder(j)
           KT(k) = KT(j)
           OKT(k) = OKT(j)
           oldOC(k) = oldOC(j)
           ListPos(k) = ListPos(j)
           SiteDen(k) = SiteDen(j)
           SPhase(k) = SPhase(j)

           FOR ii = 1 TO NSpecies + WOffSet
               SCa(ii, k) = SCa(ii, j)
               SCb(ii, k) = SCb(ii, j)
               tSCa(ii, k) = tSCa(ii, j)
               tSCb(ii, k) = tSCb(ii, j)
               SCa(ii, j) = 0
               SCb(ii, j) = 0
               tSCa(ii, j) = 0
               tSCb(ii, j) = 0
           NEXT ii

           FOR ii = 1 TO NPhase
               SP(ii, k) = SP(ii, j)
               tSP(ii, k) = tSP(ii, j)
               SP(ii, j) = 0
               tSP(ii, j) = 0
           NEXT ii

       NEXT i
       '
       ' These variables shouln't be affected by the additional components
       '
       'ChangeMass (NPhases)
       'hSEC (NPhases)
       'tSEC (NPhases)
       'PhaseMass (NPhases)
       'PhaseOrder (NPhases)
       'PSubOrder (NPhases)
       'PType (NPhases)
       'SComp (NPhases)
       'SDeltaH (NPhases)
       'SEC (NPhases)
       'SECTemp (NPhases)
       'STempCorr (NPhases)

       IF iDebugOld THEN
          FOR i = 1 TO NSpecies
             PRINT #ii, Label$(i); EC(i)
          NEXT i
          CLOSE
          'END
       END IF

    ' *----------------------------------*
    ' |  Give the WHAM components names  |
    ' *----------------------------------*
       FOR i = 1 TO 8
           Label(NComp + i - NBL - hs_i + 1) = hs$ & "_H" & CHR$(i + 48)
       NEXT i
       FOR j = 1 TO 4
           FOR ii = (j+1) TO (7 + j MOD 2) STEP 2
               Label(NComp + i - NBL - hs_i + 1) = hs$ & "_H" & CHR$(j + 48) & CHR$(ii + 48)
               i = i + 1
           NEXT ii
       NEXT j
       Label(NComp + WOffSet - NBL - hs_i + 1) = "R_" & hs$

       FirstHS = NComp - NBL - hs_i + 2 '+1

       ' Add in the protonated components
       FOR i = FirstHS TO FirstHS + WOffSet - 2 '1
          ' Component Charge
          V(i) = ZERO

          ' Species and Component Type
          SType(i) = stAqueous
          CType(i) = ctFloat

          ' Activity correction
          IF (hs_i = 1) THEN
              ActCorr(i) = acWHAMHA
          ELSEIF (hs_i = 2) THEN
              ActCorr(i) = acWHAMFA
          END IF

          '  Moles of sites per mole C
          IF ((i - FirstHS + 1) <= 4) THEN '((i - NComp) <= 4) then
              SiteDen(i) = (1 - wFPR!(hs_i)) * wNCOOH!(hs_i) / 4000
          ELSEIF ((i - FirstHS + 1) <=8) THEN
              SiteDen(i) = (1 - wFPR!(hs_i)) * wNCOOH!(hs_i) / 8000
          ELSE
              SiteDen(i) = wFPR(hs_i) * wNCOOH!(hs_i) / 16000
          END IF
          Sol_HS(hs_i) = Sol_HS(hs_i) + SiteDen(i) ' at this stage, this is the amount of HS(hs_i) in solution per mgC/L of DOC

          '  Equilibrium constant and Stoichiometry
          '  are trivial for components by definition
          EC(i) = ONE
          SCa(i, i) = ONE
          SCb(i, i) = ONE
'          if (hs_i = 1) then
'              KT(i) = KT(iDOC) * PercHA / 100
'          else
'              KT(i) = KT(iDOC) * (1 - PercHA / 100)
'          end if
       NEXT i


       '  Adjust dimensions to reflect these additions
       NComp = NComp + WOffSet
       NSpecies = NSpecies + WOffSet
       NFirstS = NComp + 1
       FOR i = 1 TO NList
          FOR j = 1 TO LinkList(i, 0)
             IF LinkList(i, j) > (NComp - NBL - WOffSet - hs_i + 1) THEN
                 ' only offset species with numbers greater than the BL species (the ones that were moved)
                LinkList(i, j) = LinkList(i, j) + WOffSet
             END IF
          NEXT j
       NEXT i


       'Add proton binding to monodentate sites
       '  A sites
       kk% = FindSpecies("H")
       FOR i = 1 TO 4
          ii = NSpecies + i
          jj = FirstHS + i - 1
          Label(ii) = hs$ & "_" & CHR$(i+48)
          EC(ii) = 10 ^ (-1 * (wPKHA!(hs_i) + (2 * i - 5) / 6 * wDPKHA!(hs_i)))
          SCa(ii, jj) = ONE
          SCb(ii, jj) = ONE
          SCa(ii, kk) = -ONE
          SCb(ii, kk) = -ONE
          V(ii) = -ONE
          V_Me(ii) = ONE
       NEXT i

       '  B sites
       FOR i = 5 TO 8
          ii = NSpecies + i
          jj = FirstHS + i - 1
          Label(ii) = hs$ & "_" & CHR$(i+48)
          EC(ii) = 10 ^ (-1 * (wPKHB!(hs_i) + (2 * i - 13) / 6 * wDPKHB!(hs_i)))
          SCa(ii, jj) = ONE
          SCb(ii, jj) = ONE
          SCa(ii, kk) = -ONE
          SCb(ii, kk) = -ONE
          V(ii) = -ONE
          V_Me(ii) = ONE
       NEXT i

       'Add proton binding to bidentate sites
       '  losing both protons
       jj = FirstHS + 7
       FOR i = 1 TO 4
           FOR j = (i+1) TO (7+ i MOD 2) STEP 2
               ii = ii + 1
               jj = jj + 1
               Label(ii) = hs$ & "_" & CHR$(i+48) & CHR$(j+48)
               EC(ii) = EC(NSpecies + j) * EC(NSpecies + i)
               SCa(ii, jj) = ONE
               SCb(ii, jj) = ONE
               SCa(ii, kk) = -TWO
               SCb(ii, kk) = -TWO
               V(ii) = -TWO
               V_Me(ii) = TWO
           NEXT j
       NEXT i

       '  losing first proton
       jj = FirstHS + 7
       FOR i = 1 TO 4
           FOR j = (i+1) TO (7+ i MOD 2) STEP 2
               ii = ii + 1
               jj = jj + 1
               Label(ii) = hs$ & "_" & CHR$(i+48) & "-H" & CHR$(j+48)
               EC(ii) = EC(NSpecies + i)
               SCa(ii, jj) = ONE
               SCb(ii, jj) = ONE
               SCa(ii, kk) = -ONE
               SCb(ii, kk) = -ONE
               V(ii) = -ONE
               V_Me(ii) = ONE
           NEXT j
       NEXT i

       '  losing second proton
       jj = FirstHS + 7
       FOR i = 1 TO 4
           FOR j = (i+1) TO (7+ i MOD 2) STEP 2
               ii = ii + 1
               jj = jj + 1
               Label(ii) = hs$ & "_" & CHR$(j+48) & "-H" & CHR$(i+48)
               EC(ii) = EC(NSpecies + j)
               SCa(ii, jj) = ONE
               SCb(ii, jj) = ONE
               SCa(ii, kk) = -ONE
               SCb(ii, kk) = -ONE
               V(ii) = -ONE
               V_Me(ii) = ONE
           NEXT j
       NEXT i

       '  This info is the same for all species
       FOR ii = NSpecies + 1 TO NSpecies + 44
          SType(ii) = stAqueous
          IF (hs_i = 1) THEN
              ActCorr(ii) = acWHAMHA
          ELSEIF (hs_i = 2) THEN
              ActCorr(ii) = acWHAMFA
          END IF
          TempCorr(ii) = FALSE
          DeltaH(ii) = 0
          ECTemp(ii) = 273
          InitialConc(ii) = 0
          WSpec(ii) = FALSE
       NEXT ii

       '
       '  Adjust array sizes again
       NSpecies = NSpecies + 44

    NEXT hs_i

    iR_HA = NComp - NBL
    iR_FA = NComp - NBL - 1

    FirstHS = NComp - NBL - WOffset * 2 + 1
' *----------------------------------------*
' |  Read in reactions from WHAM database  |
' *----------------------------------------*
'IF FALSE THEN
   NWHAM = 0
   FOR i% = 1 TO wNODATA%
      INPUT #iWHAMDBS, wtX%, wtN$, wtCH%, wtM%(1), wtM%(2), wtM%(3), wtS%(1), wtS%(2), wtS%(3), wtLK!, wtDH!, wtPKMHA!(1), wtPKMHA!(2)
      IF wtPKMHA!(1) = 999 THEN
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
         WHAMID(NWHAM) = wtX%
      ELSEIF (wtX% = 51) THEN
         ' OH is a special case
         NWHAM = NWHAM + 1
         WHAMID(NWHAM) = wtX%
      END IF
      IF AddSpecies%(wtM(1), wtM(2), wtM(3)) THEN'wtM1%, wtM2%, wtM3%) THEN  'Checks if all required species were given at the beginning of the WHAM dbs
         IF (wtX% = 51) THEN
             NSpecies = NSpecies + 1
             Label(NSpecies) = wtN$
             EC(NSpecies) = 10 ^ -14
             DeltaH(NSpecies) = 2935
             ECTemp(NSpecies) = 298
             SType(NSpecies) = stAqueous
             k = FindSpecies("H")
             ActCorr(NSpecies) = ActCorr(k)
             TempCorr(NSpecies) = TRUE
             InitialConc(NSpecies) = 0
             V(NSpecies) = wtCH%
             WSpec(NSpecies) = FALSE
             CMap(wtX%) = NSpecies
             SCa(NSpecies, k) = -ONE
             SCb(NSpecies, k) = -ONE
             tSCa(NSpecies, k) = -ONE
             tSCb(NSpecies, k) = -ONE
         ELSEIF (wtX% > 100) THEN
            '
            '  This is a new species entirely and must be added
            '  even in an inorganic form
            '
            NSpecies = NSpecies + 1
            Label(NSpecies) = wtN$
            DeltaH(NSpecies) = wtDH!
            ECTemp(NSpecies) = 298
            SType(NSpecies) = stAqueous
            ActCorr(NSpecies) = acDebye
            TempCorr(NSpecies) = TRUE
            InitialConc(NSpecies) = 0
            V(NSpecies) = wtCH%
            WSpec(NSpecies) = FALSE
            CMap(wtX%) = NSpecies
            log_OH% = 0
            FOR m = 1 TO 3
                IF (wtM(m) = 51) THEN
                    jj = FindSpecies%("H")
                    SCa(NSpecies, jj) = -wtS(m)
                    SCb(NSpecies, jj) = -wtS(m)
                    tSCa(NSpecies, jj) = -wtS(m)
                    tSCb(NSpecies, jj) = -wtS(m)
                    log_OH% = m
                ELSEIF (wtM(m) > 0) THEN
                    jj = FindSpecies%(wNs$(wtM(m)))
                    SCa(NSpecies, jj) = wtS(m)
                    SCb(NSpecies, jj) = wtS(m)
                    tSCa(NSpecies, jj) = wtS(m)
                    tSCb(NSpecies, jj) = wtS(m)
                END IF
            NEXT m
            IF (log_OH% <> 0) THEN
                'EC(NSpecies) = 10 ^ (wtLK! - 13.997 * wtS(log_OH))
                EC(NSpecies) = 10 ^ (wtLK! - 14 * wtS(log_OH))
            ELSE
                EC(NSpecies) = 10 ^ wtLK!
            END IF
         ELSE
            '
            '  This is a component and should already be in CHESS
            '
            CMap(wtX%) = FindSpecies%(wtN$)
         END IF

         IF wtBINDTEST% THEN
            ' We need to add species to describe bound forms
            kk = FindSpecies%(wtN$)
            FOR j = 1 TO 40 '20 'TODO - Kelly - Change this to 40 when adding in both HA and FA
               ii = FirstHS - 1 + j
               NSpecies = NSpecies + 1
               'Label(NSpecies) = LEFT$(Label(ii + j), 3) + "-" + wtN$
               Label(Nspecies) = Label(ii) & "-" + wtN$
               REPLACE "_H" WITH "_" IN Label$(NSpecies)
               SELECT CASE j
                   CASE 1 TO 4
                       EC(NSpecies) = (10 ^ (-wtPKMHA!(1))) * EC(kk)                  'A site MON (HA)
                   CASE 5 TO 8
                       EC(NSpecies) = (10 ^ -(wtPKMHA!(1) * 3 - 3)) * EC(kk)          'B site MON (HA)
                   CASE 9,10,13,16
                       EC(NSpecies) = (10 ^ (-wtPKMHA!(1) * 2)) * EC(kk)              '2 x A sites BID (HA)
                   CASE 11,12,14,15,17,18,19,20
                       EC(NSpecies) = (10 ^ -(wtPKMHA!(1) * (1 + 3) - 3)) * EC(kk)    'A & B sites BID (HA)
                   CASE 21 TO 24 'TODO - Kelly - Change these case #'s to +20 when add in HA
                       EC(NSpecies) = (10 ^ (-wtPKMHA!(2))) * EC(kk)                  'A site MON (FA)
                   CASE 25 TO 28
                       EC(NSpecies) = (10 ^ (-wtPKMHA!(2) * 3.96)) * EC(kk)           'B site MON (FA)
                   CASE 29,30,33,36
                       EC(NSpecies) = (10 ^ (-wtPKMHA!(2) * 2)) * EC(kk)              '2 x A sites BID (FA)
                   CASE 31,32,34,35,37,38,39,40
                       EC(NSpecies) = (10 ^ (-wtPKMHA!(2) * (1 + 3.96))) * EC(kk)     'A & B sites BID (FA)
               END SELECT
               SCa(NSpecies, ii) = ONE
               SCb(NSpecies, ii) = ONE
               k = FindSpecies%("H")
               IF ((j<=8) OR (j>20 AND j<=28)) THEN
                   SCa(NSpecies, k) = -ONE
                   SCb(NSpecies, k) = -ONE
                   tSCa(NSpecies, k) = -ONE
                   tSCb(NSpecies, k) = -ONE
                   V(NSpecies) = wtCH% - ONE
               ELSEIF ((j<=20) OR (j>20 AND j<=40)) THEN
                   SCa(NSpecies, k) = -TWO
                   SCb(NSpecies, k) = -TWO
                   tSCa(NSpecies, k) = -TWO
                   tSCb(NSpecies, k) = -TWO
                   V(NSpecies) = wtCH% - TWO
               END IF
               IF (kk < NComp) THEN
                  SCa(NSpecies, kk) = ONE
                  SCb(NSpecies, kk) = ONE
                  tSCa(NSpecies, kk) = ONE
                  tSCb(NSpecies, kk) = ONE
               ELSE
                  FOR m = 1 TO 3
                     IF (wtM(m) = 51) THEN
                        k = FindSpecies%("H")
                        SCa(NSpecies, k) = SCa(NSpecies, k) - wtS(m)
                        SCb(NSpecies, k) = SCb(NSpecies, k) - wtS(m)
                        tSCa(NSpecies, k) = SCa(NSpecies, k) - wtS(m)
                        tSCb(NSpecies, k) = SCb(NSpecies, k) - wtS(m)
                     ELSEIF (wtM(m) > 0) THEN
                        k = FindSpecies%(wNs$(wtM(m)))
                        SCa(NSpecies, k) = wtS(m)
                        SCb(NSpecies, k) = wtS(m)
                        tSCa(NSpecies, k) = wtS(m)
                        tSCb(NSpecies, k) = wtS(m)
                     END IF
                     V_Me(NSpecies) = V_Me(NSpecies) + wtS(m) * V(k)
                  NEXT m
               END IF
               DeltaH(NSpecies) = wtDH!
               SType(NSpecies) = stAqueous
               IF (j <= 20) THEN
                   ActCorr(NSpecies) = acWHAMHA
               ELSEIF (j > 20) THEN
                   ActCorr(NSpecies) = acWHAMFA
               END IF
               TempCorr(NSpecies) = FALSE
               InitialConc(NSpecies) = 0
            NEXT j
         END IF 'if BINDTEST
      END IF 'if add species
      '
      '  Add data to WHAM species list
      '
      wX% = wtX%
      wNs$(wX%) = wtN$
      wCH%(wX%) = wtCH%
      wM1%(wX%) = wtM(1)
      wM2%(wX%) = wtM(2)
      wM3%(wX%) = wtM(3)
      wS1%(wX%) = wtS(1)
      wS2%(wX%) = wtS(2)
      wS3%(wX%) = wtS(3)
      wLK!(wX%) = wtLK!
      wDH!(wX%) = wtDH!
      wPKMHA!(1, wX%) = wtPKMHA!(1)
      wPKMHA!(2, wX%) = wtPKMHA!(2)
   NEXT i% 'for every line in WHAM database

'   ' *---------------------------*
'   ' |  Add in DL "R" component  |
'   ' *---------------------------*
'    WOffSet = 2
'    iStartWHAM = NComp
'    FOR i = 0 TO (NSpecies - NComp) '-1
'       j = NSpecies - i
'       k = NSpecies - i + WOffSet
'       ActCorr(k) = ActCorr(j)
'       ActCoef(k) = ActCoef(j)
'       CtoM(k) = CtoM(j)
'       DeltaH(k) = DeltaH(j)
'       EC(k) = EC(j)
'       ECTemp(k) = ECTemp(k)
'       hEC(k) = hEC(j)
'       InitialConc(k) = InitialConc(j)
'       Label(k) = Label(j)
'       LinkComp(k) = LinkComp(j)
'       MeanLogK(k) = MeanLogK(j)
'       MonteList(k) = MonteList(j)
'       OC(k + FOUR) = OC(j + FOUR)
'       SiteDen(k) = SiteDen(j)
'       SType(k) = SType(j)
'       tEC(k) = tEC(j)
'       TempCorr(k) = TempCorr(j)
'       V(k) = V(j)
'       VarLogK(k) = VarLogK(j)
'
'       FOR ii = 1 TO NComp
'          SCa(k, ii) = SCa(j, ii)
'          SCb(k, ii) = SCb(j, ii)
'          tSCa(k, ii) = tSCa(j, ii)
'          tSCb(k, ii) = tSCb(j, ii)
'          SCa(j, ii) = 0
'          SCb(j, ii) = 0
'          tSCa(j, ii) = 0
'          tSCb(j, ii) = 0
'       NEXT ii
'
'       FOR ii = 1 TO NList
'          LinkList(ii, k) = LinkList(ii, j)
'       NEXT ii
'    NEXT i
'
'    FOR ii = 1 TO NSpecies + wOffset
'        SCa(ii, NComp + wOffset) = SCa(ii, NComp)
'        SCb(ii, NComp + wOffset) = SCb(ii, NComp)
'        tSCa(ii, NComp + wOffset) = tSCa(ii, NComp)
'        tSCb(ii, NComp + wOffset) = tSCb(ii, NComp)
'        SCa(ii, NComp) = 0
'        SCb(ii, NComp) = 0
'        tSCa(ii, NComp) = 0
'        tSCb(ii, NComp) = 0
'    NEXT ii
'
'    iR_HA = NComp
'    Label(iR_HA) = "R_HA"
    V(iR_HA) = ZERO
    SType(iR_HA) = stAqueous
    CType(iR_HA) = ctFloat
    ActCorr(iR_HA) = acNone
    SiteDen(iR_HA) = ONE
    EC(iR_HA) = ONE
    SCa(iR_HA, iR_HA) = ONE
    SCb(iR_HA, iR_HA) = ONE
    tSCa(iR_HA, iR_HA) = ONE
    tSCb(iR_HA, iR_HA) = ONE
'
'    iR_FA = NComp + 1
'    Label(iR_FA) = "R_FA"
    V(iR_FA) = ZERO
    SType(iR_FA) = stAqueous
    CType(iR_FA) = ctFloat
    ActCorr(iR_FA) = acNone
    SiteDen(iR_FA) = ONE
    EC(iR_FA) = ONE
    SCa(iR_FA, iR_FA) = ONE
    SCb(iR_FA, iR_FA) = ONE
    tSCa(iR_FA, iR_FA) = ONE
    tSCb(iR_FA, iR_FA) = ONE
'
'    '  Adjust dimensions to reflect these additions
'    NComp = NComp + WOffSet
'    NSpecies = NSpecies + WOffSet
'    NFirstS = NComp + 1
'    FOR i = 1 TO NList
'       FOR j = 1 TO LinkList(i, 0)
'          IF LinkList(i, J) > (NComp - WOffSet) THEN
'              ' only offset species, with species numbers more than the number of components before adding WHAM components
'             LinkList(i, j) = LinkList(i, j) + WOffSet
'          END IF
'       NEXT j
'    NEXT i

   ' *---------------------*
   ' |  Add in DL species  |
   ' *---------------------*
   Num_DL = 0
   FOR i = 1 TO NSpecies
       IF ((V(i) > 0) AND (ActCorr(i) <> acWHAMHA) AND (ActCorr(i) <> acWHAMFA) AND (SType(i) = stAqueous)) THEN
           'HA Donnan Layer species
           Num_DL = Num_DL + 1
           k = NSpecies + Num_DL
           Label(k) = Label(i) & "_DL_HA"
           FOR c = 1 TO NComp
               SCa(k, c) = SCa(i, c)
               SCb(k, c) = SCb(i, c)
           NEXT c
           SCa(k, iR_HA) = V(i)
           SCb(k, iR_HA) = V(i)
           EC(k) = EC(i)
           V(k) = V(i)
           CType(k) = ctFloat
           SType(k) = stDonnanHA
           ActCorr(k) = acNone
           TempCorr(k) = FALSE
           InitialConc(k) = 0
           WSpec(k) = FALSE
           DeltaH(k) = 0
           ECTemp(k) = 273

           'FA Donnan Layer Species
           Num_DL = Num_DL + 1
           k = NSpecies + Num_DL
           Label(k) = Label(i) & "_DL_FA"
           FOR c = 1 TO NComp
               SCa(k, c) = SCa(i, c)
               SCb(k, c) = SCb(i, c)
           NEXT c
           SCa(k, iR_FA) = V(i)
           SCb(k, iR_FA) = V(i)
           EC(k) = EC(i)
           V(k) = V(i)
           CType(k) = ctFloat
           SType(k) = stDonnanFA
           ActCorr(k) = acNone
           TempCorr(k) = FALSE
           InitialConc(k) = 0
           WSpec(k) = FALSE
           DeltaH(k) = 0
           ECTemp(k) = 273

       END IF
   NEXT i
   NSpecies = NSpecies + Num_DL

   ' Add one to NWHAM to account for DOC
   NWHAM = NWHAM + 1

'END IF
   CLOSE #iWHAMDBS
'   CALL WHAMDefine
'   PHFIX$ = "YES"
'   PCO2% = 999
'   PRECISION! = .01

   'CALL PrintI(CMap%(), NSP%)

ELSE
    FOR i = NFirstS TO NSpecies
        IF ((ActCorr(i) = acWHAMHA) OR (ActCorr(i) = acWHAMFA)) THEN
            T$ = RIGHT$(Label(i), LEN(Label(i)) - INSTR(Label(i), "-"))
            IF (T$ <> Label(i)) THEN
                k = FindSpecies(T$)
                V_Me(i) = V(k)
            ELSE
                V_Me(i) = - V(k)
            END IF
'            V_Me(i) = 0
'            for j = 1 to NComp
'                V_Me(i) = V_Me(i) + SCa(i,j) * V(j)
'            next j
        END IF
    NEXT i
END IF

iWHAM = FALSE

iAppend = FALSE
iClose = TRUE
iLogHand = FREEFILE

IF NPhase > 0 THEN
   CALL SelectComp(SComp(), SPhase(), ErrorCode%, ErrorMsg$)
END IF

CALL SSaturate

CALL SimplifyProblem

END SUB

