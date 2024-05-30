
DEFINT I-N
DEFDBL A-H, O-Z
' ------------------------------------------------------------------------
'
'  Biotic Ligand Model
'
'  A metal toxicity simulation program based on the CHESS model
'
'  Originally developed by
'     Robert Santore (while at HydroQual, Inc.)
'     Currently at Windward Environmental, LLC
'     206-812-5450
'  for US EPA, Office of Water, Health and Ecological Criteria Division
'  Delivered September 1999
'
'  See SUB VersionInfo for more details regarding updates
'
'  old email addresses (no longer working)
'  bsantore@hydroqual.com
'  rsantore@hydroqual.com
'  robert.santore@hdrinc.com
'
'  current email
'  robert.santore@gmail.com
'  RobertS@windwardenv.com
'
' ------------------------------------------------------------------------
'
DECLARE FUNCTION AddSpecies% (wtM1%, wtM2%, wtM3%)
DECLARE FUNCTION CalcEN# ()
DECLARE FUNCTION CheckFile%(Strng$)
DECLARE FUNCTION Davies# (ISTR AS DOUBLE, ZZ AS SINGLE)
DECLARE FUNCTION DebyeHuckel# (ISTR#, ZZ!)
DECLARE FUNCTION Delimit%(Work$, Delim$)
DECLARE FUNCTION FirstName$ (WName$)
DECLARE FUNCTION FindSpecies% (Spec$)
DECLARE FUNCTION FindWHAMSpecies% (s$)
DECLARE FUNCTION FormatStr$(Strng$, StrMask$)
DECLARE FUNCTION FUsing$ (S$, T$)
DECLARE FUNCTION IAP# (P AS INTEGER)
DECLARE FUNCTION LOG10S! (X!)
DECLARE FUNCTION Log10D# (X#)
DECLARE FUNCTION QPTrim$ (TmpStr$)
DECLARE FUNCTION OpenEXCEL(mFileName$, xlsFileNumber&) AS LONG
DECLARE FUNCTION OpenEXCELWQC(mFileName$, xlsFileNumber&) AS LONG
DECLARE FUNCTION PATHNAMEBLM$ (FullName$)
DECLARE FUNCTION PeekBuf% ()
DECLARE FUNCTION PsiToP# (Psi#)
DECLARE FUNCTION QPTrim$ (TmpStr$)
DECLARE FUNCTION Sign% (X#)
DECLARE FUNCTION SignS% (X!)
'DECLARE FUNCTION ZScore(X AS LONG, N AS LONG) AS DOUBLE
'
' Subroutine declarations
'
DECLARE SUB AdjustCu (Diff#, Iter%)
DECLARE SUB CalcActivityCoef (ISTR AS DOUBLE, ErrorCode AS INTEGER, ErrorMsg AS STRING)
DECLARE SUB CalcMassChange ()
DECLARE SUB CalcResidual (MaxError!, MaxComp%)
DECLARE SUB CalcSpeciesConc (ISTR#, ErrorCode%, ErrorMsg$)
DECLARE SUB ChargeBalance (CHARGE#, ISTR#)
DECLARE SUB CheckSolidPhase (Converge%, UnderFlag%, OverFlag%)
DECLARE SUB CHESS (ErrorCode%, ErrorMsg$)
DECLARE SUB CHESSLog (LogFileName$, iAppend%, iClose%, iLogHand%, ParameterFile$)
DECLARE SUB DefineProblem300W (FName$, ErrorCode%, ErrorMsg$)
DECLARE SUB FileClose ()
DECLARE SUB FirstGuess (ErrorCode%, ErrorMsg$, ISTR#)
DECLARE SUB Gauss (A() AS DOUBLE, B() AS DOUBLE, X() AS DOUBLE, N%, iSing%)
DECLARE SUB GetData (iInFile%, LC50%, iExpert1%, CriticalValue#)
DECLARE SUB GetArgs (Arg1$, Arg2$, Arg3$, LC50ALL%, iDeriv%, iExpert1%, iExpert2%, iExpert3%, ErrorCode%, ErrorMsg$)
DECLARE SUB IndSortI (iArray%(), Index%(), NEls%, iDir%)
DECLARE SUB InitData (iInFile%, sInFile$, NData&, ErrorCode%, ErrorMsg$)
DECLARE SUB JacobianAA ()
DECLARE SUB NParse(Work$, Delim$, NArray$(), k%)
DECLARE SUB PrintResidual (Iter%, MaxError!, MaxComp%, iSing%)
DECLARE SUB SSaturate ()
DECLARE SUB SelectComp (SComp() AS INTEGER, SPhase() AS INTEGER, ErrorCode%, ErrorMsg$)
DECLARE SUB SetCtoM (AqueousCtoM!, SurfaceCtoM!)
DECLARE SUB SetPhaseOrder (PhaseOrder() AS INTEGER, PSubOrder() AS INTEGER)
DECLARE SUB SimplifyProblem ()
DECLARE SUB SUpdate (X() AS DOUBLE, Iter%, MaxError!, Skip%)
DECLARE SUB TempCorrection ()
DECLARE SUB WriteDeriv (iHand%, Iter%, a#(), DLC50#(), n%)
DECLARE SUB WriteDetailed (iHand%, Iter%)
DECLARE SUB WriteSimple (iHand%, Iter%, LC50%)
DECLARE SUB WriteDetailed2 (iHand&, Iter, iComplete)
DECLARE SUB WriteSimple2 (iHand&, Iter, LC50, iComplete)
DECLARE SUB WriteDetailedComma (iHand&, Iter, iComplete)
DECLARE SUB WriteSimpleComma (iHand&, Iter, LC50, iComplete)
'#OPTION VERSION3  Needed for NT 3x
'#OPTION VERSION4  Needed for 95/98/ME or NT4
'#OPTION VERSION5  Needed for Windows 2000 (NT5) or Windows XP (NT5.1).
'
' Inherit CHESS common variables
'
'$INCLUDE: 'E:\Users\BLM_Devl\common1.bi'
'$INCLUDE: 'E:\Users\BLM_Devl\common2.bi'
'$INCLUDE: 'E:\Users\BLM_Devl\blm_comm.bi'
'$INCLUDE: 'E:\Users\BLM_Devl\constant.bi'
'#OPTION VERSION5
'#INCLUDE "Excel.inc"

#INCLUDE "MLE_FIT2.BAS"  'include routines for MLE fitting
'
' Use #OPTION VERSION3 to make the compiled output file require
' a minimum of Windows 95 or NT 3.1. That includes Windows 95,
' 98, ME, Windows NT 3.1-4.0, Windows 2000, XP, Windows 2003
' (and later).
'
' Use #OPTION VERSION4  (default) to make the compiled output
' file require a minimum of Windows 95 or NT4. That includes
' Windows 95, 98, ME, Windows NT 4.0, Windows 2000, XP,
' Windows 2003 (and later).
'
' Use #OPTION VERSION5 to make the compiled output file
' require a minimum of Windows 2000. That includes Windows 2000,
' XP, Windows 2003 (and later).
'
#OPTION VERSION3
#INCLUDE "Excel.inc"

%XLSFALSE = 0
%XLSTRUE = NOT %XLSFALSE

' *-------------------------------------------------------------*
' |  Type variables definitions and declarations                |
' *-------------------------------------------------------------*
TYPE OFContents   ' Type variable for defining output file contents
   Conc AS INTEGER
   Tot AS INTEGER
   Keq AS INTEGER
   Resid AS INTEGER
   UserDF AS INTEGER
   ActC AS INTEGER
   PhaseSI AS INTEGER
   PhaseIAP AS INTEGER
END TYPE

'SUB Pmsg(Strng$)
'   'DIM ibl AS STRING
'   PRINT Strng$ + "   "
'   'ibl$ = WAITKEY$
'   DO UNTIL LEN(INKEY$)
'   LOOP
'END SUB
'
'FUNCTION Exist(DirFileName AS STRING) AS LONG
'   FUNCTION = GETATTR(DirFileName$): FUNCTION = (ERR = 0)
'END FUNCTION

'
'
' Global version string to ID model output
'
FUNCTION PBMAIN()
DIM StartTime AS LONG
DIM StopTime AS LONG
DIM TimeElapsed AS LONG
DIM Tmp1 AS DOUBLE
DIM Tmp2 AS DOUBLE
DIM Tmp3 AS DOUBLE
StartTime = TIMER
iDBIter = 0
'DIM SHARED VersionID AS STRING
DIM VersionID AS GLOBAL STRING
'VersionID$ = "ver 2.1.0, build 2004-03-14"
VersionID$ = "ver AP08, build 2003-05-22"
VersionID$ = "ver 2.2.6, build 2008-10-23"
VersionID$ = "ver 2.2.6, build 2009-04-08"     ' increased MaxSpecies to 300
VersionID$ = "ver 2.2.6, build 2010-03-29"     ' The higher Al sol Ksp (matching Niva) 'was prev only used
                                               ' for /Al1, and is now set to be used for /Al2 as well
VersionID$ = "ver 2.2.7, build 2011-09-11"
VersionID$ = "ver 2.2.8, build 2012-05-03"       ' added Pb solubility by Cerrusite
VersionID$ = "ver 2.2.9, build 2012-10-10"       ' added Al dissolved-particulate mixtures
VersionID$ = "ver 2.3.0, build 2014-07-03"       ' now uses last free ion solution as first guess for next run, added /B2 switch for output log BLM_LA50 for use in optimization
VersionID$ = "ver 2.3.0, build 2014-07-07"       ' changed total values in detailed xls
VersionID$ = "ver 2.3.1, build 2014-07-28"       ' changed variable type for NData from integer to longint to allow more than 32767 input observations
VersionID$ = "ver 2.3.1, build 2014-07-30"       ' fixed trailing '2' left from using /CO2 switch
VersionID$ = "ver 2.3.2, build 2014-08-26"       ' fixed trailing 'M' left from using /M switch
VersionID$ = "ver 2.3.3, build 2014-10-09"       ' revised particulate Al approach
VersionID$ = "ver 2.3.4, build 2014-10-16"       ' revised dissolved Al approach
VersionID$ = "ver 2.3.5, build 2014-10-22"       ' revised dissolved Al approach
VersionID$ = "ver 2.3.6, build 2014-11-07"       ' revised precipitated Al normalization approach
VersionID$ = "2.37, build 2015-10-22"           ' revised precipitated Al normalization approach, removed particulate species from output when not Al, added WQC, FMB, and /V switches (merge with blm_212g)
VersionID$ = "2.38, build 2016-02-28"           ' Fixed problem with solubility switches
VersionID$ = "2.39, build 2016-06-15"           ' Fixed the problem with using the low hardness adjustment in expert toxicity mode, made changes to prevent problems with "Iterations exceed max"
                                                ' and now output NA in incomplete cases or when iterations exceed max.
VersionID$ = "2.40, build 2016-08-18"           ' Fixed activity components will now print KT into the output file, rather than OKT
VersionID$ = "2.40, build 2016-12-21"           ' Added capability to use a critical accumulation table (CAT) in the parameter file.  This is a table that has the critical accumation along with organism & test information and references
VersionID$ = "2.40, build 2017-01-16"           ' Added error handling to opening the BLM input file.  If the file cannot be accessed, it will give an error code.
VersionID$ = "2.41, build 2017-09-22"           ' Added iteration cap to LC50 finding, and started implementing a message file for communication with the interface.
VersionID$ = "2.42, build 2017-11-27"           ' Corrected a bug in the /CA switch.  In 2.41, a number > 10 after this switch would be ignored and the first value would be used.
VersionID$ = "2.43, build 2018-03-29"           ' Updated the /CO2 switch so that it could be specified as /CO2x, where x is the pCO2 to use for the open system.  Will still default to 3.5 if not specified
VersionID$ = "2.43, build 2018-04-01"           ' Got the comma-separated output format working (/O2 switch).
VersionID$ = "2.44, build 2018-05-04"           ' Fixed the /CO2 switch - it was cutting off the last number when the switch was put at the end of the arguments list.
VersionID$ = "2.44, build 2018-09-12"           ' Fixed error with the FMB not working with files that are not the EPA Cu WQC file -  set defaults of 1 for ACUTE_DIV and CHRONIC_DIV, but have not enabled yet.
VersionID$ = "2.44, build 2018-10-30"           ' Fixed error where CMC/CCC and TU's are calculated for incomplete data as well, with just the last complete observation being subbed in, enabled defaults for ACUTE_DIC and CHRONIC_DIV
VersionID$ = "2.44, build 2018-11-05"           ' Fixed error when a non-3.01+ parameter file is run with a /CA switch
VersionID$ = "2.44, build 2018-11-08"           ' Fixed problem with CSV files writing weird quotations around site/sample labels
VersionID$ = "2.45, build 2018-11-20"           ' Switched to 2.45 version number, no other changes
VersionID$ = "2.45b, build 2019-02-28"           ' When using the /B3 switch, it will now give -999 for the two values if it's an incomplete observation.
VersionID$ = "2.45c, build 2021-04-05"           ' A minor update to version 2.45 to add activity correction info to the detailed output file
VersionID$ = "2.45d, build 2024-03-06"           ' Added additional WHAM related information to the output file
#INCLUDE "all_dims.bi"


'NOTES:
'build 2017-09-22:
'There may have been a version of 2.41 that was compiled but not documented.
'
'We have source code and compiled version of 2.41 with date/time of
'6/13/2017 7:53 AM.  We'll call this 2.41orig
'
'The BCMOE package 1.05 based on 3.15.2.41b has a compiled version of 2.41
'with a date/time stamp of 6/13/2017 11:16 AM. We'll call this 2.41BCMOE
'
'the 2.41orig version had a bug related to line 795 (search for "edit 2.41")
'
'The 2.41BCMOE compiled version does not seem to have this bug.  But we do not
'have code with the same time stamp to document what is different.
'
'The current version of 2.41 was compiled 9/22/2017 at 2:55 PM and fixes the
'bug at line 795.  We believe it is the same as 2.41BCMOE but we cannot confirm.



'build 2000-02-04:  Added line 222 "OC(j) = OKT(j)" to allow detailed
'                   output to have actual DOC in input file
'
'build 2000-06-12:  Added organic metal to the output file
'                   Now uses last successful LC50 as the initial guess
'
'build 2000-06-20:  Combined code to allow d[LC50]/d[input param]
'                   and "expert" input with critical BL numbers
'
'                   New line switches with this build:
'                   /d  -  model outputs derivitive information
'                          d[lc50] / d[C] for each input C
'                          Note that "/d" automatically generates LC50s
'                          so that "/l" is not needed.
'                   /e  -  allows "expert mode" with critical
'                          sensitivity values in the input file
'
'build 2000-07-10:  Changed the units in the derivitive file for DOC
'                   (to mg/L) and H (to pH)
'
'                   New line switch:
'                   /v  -  Report the version number and exit
'
'build 2001-03-16:  Added column to simple output TOrg which
'                   includes other metals species bound to organics
'
'build 2001-05-09:  Removed some WHAM species from detailed output
'
'New ver AP08, build 2001-09-25
'
'                   This version compiled with powerbasic compiler.  Many
'                   source code changes were required to adapt to new compiler.
'
'                   New line switches with this build:
'                   /s   use scripts
'                        usage:  blm_ap08 /s <script file name>
'                        The script file must be a text file with three lines:
'                        Line1 contains the parameter file name
'                        Line2 contains the input file name
'                        Line2 contains any other valid options
'
'                   /w   use input file created by windows interface
'
'                   /a1  use older method for LC50 determination (false position)
'                   /a2  use Newton-Raphson for LC50 determination (DEFAULT)
'                   /a3  use cubic spline interpolation for LC50 determination
'                   /a4  use cubic polynomial interpolation for LC50 determination
'
'                   /q   Quiet mode (report only obs number to screen)
'                   /qq  Very quiet mode (report nothing to screen)
'
'build 2001-10-02:  Fixed a bug with script files that use lower case for options
'
'                   Stopped printing time elapsed to screen when run in very
'                   quiet mode
'
'                   Now reports time elapsed to a file "elapsed" when run in
'                   very quiet mode
'
'build 2001-12-03:  Added new switch for debugging input files
'
'                   /I  Echo inputs to output without running anything
'
'build 2002-02-22:  Fixed ability to handle input files with incomplete records
'
'                   Can now also use input files that have added LA50 column when
'                   not run in expert mode
'
'build 2002-05-07:  Now labels columns in simple output as "total" concentrations
'
'build 2003-05-22:  Added variable convergence criteria and /C switch
'
'build 2003-06-05:  Increased number of significant digits in simple output file
'
' Simple local variables with initialization
'
DIM i AS INTEGER           ' Loop counter
    i = 0
DIM LC50 AS INTEGER        ' Boolean flag when TRUE we are predicting LC50
    LC50 = FALSE
DIM LC50ALL AS INTEGER     ' Boolean flag when TRUE we are predicting LC50
    LC50ALL = FALSE        '  no matter what the value of LC50
DIM iDeriv AS INTEGER       ' Boolean flag to determine if this run will
    iDeriv = FALSE          '  output derivitive information
DIM iExpert1 AS INTEGER     ' Flag to determine expert mode
    iExpert1 = FALSE
DIM iExpert2 AS INTEGER     ' Flag to determine expert mode
    iExpert2 = FALSE
DIM iExpert3 AS INTEGER     ' Flag to determine expert mode
    iExpert3 = FALSE
DIM LastLC50#              ' Used to store last successful LC50 value
DIM ConvCrit AS SINGLE

'iWHAM = TRUE
'iWHAM = FALSE


'
' Get input files and option switches from the command line
'
CALL GetArgs(Arg1$, Arg2$, Arg3$, LC50ALL, iDeriv, iExpert1, iExpert2, iExpert3, ErrorCode%, ErrorMsg$)
IF ErrorCode% THEN
   IF (LEN(ErrorMsg$)>0) THEN
      PRINT
      PRINT ErrorMsg$
   END IF
   EXIT FUNCTION
END IF
'sInFile$ = "copp_in.txt"               ' Input file with chemistry to solve
'sInFile$ = "b_test1.txt"               ' Input file with chemistry to solve
'sInFile$ = "b_test2.txt"               ' Input file with chemistry to solve
iInFile = 1
'sOutF1$ = "simple.out"                 ' Simple output, what most folks want
iOutF1& = 2
'sOutF2$ = "detail.out"                 ' Detailed output, the works
iOutF2& = 3
iOutF3 = 5
'sWHAM = "water10.dbs"                  ' WHAM thermo data base
iWHAMDBS = 4
sBLM$ = "fhminnow.dbs"                 ' Biotic ligand model parameters
iBLM = 5
'sCHESS$ = "copper.dat"                 ' CHESS model data base
sCHESS$ = Arg1$                        ' CHESS model data base
sInFile$ = Arg2$                       ' Input file with chemistry to solve
sOutF1$ = FirstName$(sInFile$) + ".sim"
sOutF2$ = FirstName$(sInFile$) + ".det"
IF (iOFormat = ONE) THEN
   sOutF1$ = sOutF1$ + ".txt"
   sOutF2$ = sOutF2$ + ".txt"
ELSEIF (iOFormat = TWO) THEN
   sOutF1$ = sOutF1$ + ".csv"
   sOutF2$ = sOutF2$ + ".csv"
ELSEIF (iOFormat = THREE) THEN
   sOutF1$ = sOutF1$ + ".xls"
   sOutF2$ = sOutF2$ + ".xls"
   IF (iCriteria = ONE) THEN
   END IF
END IF
IF (iCriteria <> ZERO) THEN
   sOutF4txt$ = FirstName$(sInFile$) + ".wqc.txt"
   sOutF4xls$ = FirstName$(sInFile$) + ".wqc.xls"
END IF

IF (iDebug > 0) THEN
   ' Open debug file
   'iDBOut = FREEFILE
   iDBOut = 6
   OPEN FirstName$(sInFile$) + ".debug" FOR OUTPUT AS iDBOut&
END IF

' Default convergence criterion when iConverge = ONE or is not specified
'
ConvCrit = .001
IF iConverge = TWO THEN
   ConvCrit = ConvCrit / 10
ELSEIF iConverge = THREE THEN
   ConvCrit = ConvCrit / 100
ELSEIF iConverge = FOUR THEN
   ConvCrit = ConvCrit / 1000
ELSEIF iConverge = FIVE THEN
   ConvCrit = ConvCrit / 10000
END IF
'T$ = FORMAT$(ConvCrit)
'CALL Pmsg(T$)

'
' Order of components in WHAM database
'
'DIM SHARED OKT(1 TO 40) AS DOUBLE
DIM OKT(1 TO 40) AS GLOBAL DOUBLE
'DIM WHAMID%(15)        ' WHAM database ID of needed components

IF (iDebug = 01 OR iDebug = 02) THEN CALL DebugLog ("main", "Before DefineProblem")

'
' Load appropriate model into CHESS
'
'IF iWHAM THEN
'   CALL DefineProblem300W(sCHESS$, ErrorCode%, ErrorMsg$)
'ELSE
'   CALL DefineProblem300(sCHESS$, ErrorCode%, ErrorMsg$)
'END IF
CALL DefineProblem300W(sCHESS$, ErrorCode%, ErrorMsg$)

' ---------------------------------------------------------------------
'    Need to check on whether we can do diss and part effects for Al
' ---------------------------------------------------------------------
'   these new parameters are specified in parameter file
'   and should now have values other than zero
'       DIM DISS_A AS GLOBAL DOUBLE        ' intercept of response curve based on BL-metal
'       DIM DISS_B AS GLOBAL DOUBLE        ' slope of response curve based on BL-metal
'       DIM PART_A AS GLOBAL DOUBLE        ' intercept of response curve based on particulate metal
'       DIM PART_B AS GLOBAL DOUBLE        ' slope of response curve based on particulate metal
'       DIM PART_H AS GLOBAL DOUBLE        ' hardness dependence of particulate response
'       DIM TARGET_R AS GLOBAL DOUBLE      ' target reduction in biological effect due to toxicity (as fraction of control)
'
'   These are calculated internally
'
'       DIM DISS_Effect AS GLOBAL DOUBLE   ' portion of target effect acheivable by dissolved metal
'       DIM Part_Conc AS GLOBAL DOUBLE     ' particulate metal concentration needed if Diss_Effect < TARGET_R
'
'   We also will shut off all the particulate toxicity calculations for speciation runs (i.e., LC50ALL = FALSE)
IF (DISS_A = 0) OR (DISS_B = 0) OR (PART_HI = 0) OR (PART_HM = 0) OR (PART_B = 0) OR (LC50ALL = FALSE) THEN
    iPartMe = FALSE
ELSE
    iPartMe = TRUE
END IF


IF (iDebug = 01 OR iDebug = 02) THEN CALL DebugLog ("main", "After DefineProblem")

IF ErrorCode% THEN
   PRINT
   PRINT ErrorMsg$
   EXIT FUNCTION
END IF
'
' Variables for Derivitive Information
'
DIM NDeriv AS DOUBLE
DIM DLC50(1 TO NComp + 1) AS DOUBLE
DIM KeepKT AS DOUBLE

'DIM WHAMMAP%(NComp)       ' Map of WHAM components to CHESS

IF iWHAM THEN
   FOR i = 1 TO NComp
      WHAMMAP%(i) = FindWHAMSpecies%(Label(i))
   NEXT i
END IF

'CALL PrintI(WHAMID%(), NWHAM)
'CALL PrintI(WHAMMAP%(), NWHAM)

'WHAMID%(1) = 1  ' H
'WHAMID%(2) = 3  ' Na
'WHAMID%(3) = 4  ' Mg
'WHAMID%(4) = 6  ' K
'WHAMID%(5) = 7  ' Ca
'WHAMID%(6) = 14 ' Cu
'WHAMID%(7) = 52 ' Cl
'WHAMID%(8) = 54 ' SO4
'WHAMID%(9) = 51 ' OH
'WHAMID%(10) = 55 ' CO3
'NWHAM = 10

'
' Order of components in CHESS
'
'WHAMMAP%(1) = 1  ' H
'WHAMMAP%(2) = 14 ' Cu
'WHAMMAP%(3) = 0  ' DOC
'WHAMMAP%(4) = 7  ' Ca
'WHAMMAP%(5) = 4  ' Mg
'WHAMMAP%(6) = 3  ' Na
'WHAMMAP%(7) = 6  ' K
'WHAMMAP%(8) = 54 ' SO4
'WHAMMAP%(9) = 52 ' Cl
'WHAMMAP%(10) = 55 ' CO3
'WHAMMAP%(11) = 0  ' S
'WHAMMAP%(12) = 0  ' Gill

' Assign array index for key species
iBL = FindSpecies%(sBL$)
FOR j = 1 TO NBLM
   iBLM(j) = FindSpecies%(sMetBL$(j))
NEXT j
iMetal = FindSpecies%(sMetal$)
iDOC = FindSpecies%(sDOC$)
iMetalDL = FindSpecies%(Metal_DL$)

'iGill = FindSpecies%("Gill")
'iGillM = FindSpecies%("Gill-Cu")
'iMetal = FindSpecies%("Cu")
'iDOC = FindSpecies%("DOC")

IF (iDebug = 01 OR iDebug = 02) THEN CALL DebugLog ("main", "Before CHESSLog")
iAppend = FALSE
iClose = TRUE
iLogHand = FREEFILE
CALL CHESSLog("BLM.LOG", iAppend, iClose, iLogHand, sCHESS$)
'
' Get first line of input file
'
IF (iDebug = 01 OR iDebug = 02) THEN CALL DebugLog ("main", "Before InitData")
CALL InitData(iInFile, sInFile$, NData&, ErrorCode%, ErrorMsg$)
IF ErrorCode% THEN
   IF (LEN(ErrorMsg$)>0) THEN
      PRINT
      PRINT ErrorMsg$
   END IF
   EXIT FUNCTION
END IF
IF (iDebug = 01 OR iDebug = 02) THEN CALL DebugLog ("main", "After InitData")
'
' Open output files
IF (iOFormat = ONE) OR (iOFormat = TWO) THEN
   OPEN sOutF1$ FOR OUTPUT AS iOutF1&
   OPEN sOutF2$ FOR OUTPUT AS iOutF2&
ELSEIF (iOFormat = THREE) THEN
   stat& = OpenEXCEL(sOutF1$, iOutF1&)
   stat& = OpenEXCEL(sOutF2$, iOutF2&)
END IF

IF (iCalcLA50 <> FALSE) THEN
   'T$ = "Opening blm_la50.txt"
   'CALL Pmsg(T$)
   iLAFile = FREEFILE
   OPEN "blm_la50.txt" FOR OUTPUT AS iLAFile
END IF
IF (iCriteria <> ZERO) THEN
    iOutF4txt = FREEFILE
   OPEN sOutF4txt$ FOR OUTPUT AS iOutF4txt&
   IF (iOFormat = THREE) OR (iCriteria <> ZERO) THEN   ' for now, let's always make a WQC file in xls format
      iOutF4xls = FREEFILE
      stat& = OpenEXCELWQC(sOutF4xls$, iOutF4xls&)
   END IF
END IF

'
' Put a version stamp on each output file
sTmp1$ = VersionID$
sTmp2$ = sCHESS$
sTmp3$ = sInFile$
sTmp4$ = UCASE$(LTRIM$(RTRIM$(COMMAND$))) + ", " + Arg3$
IF (iOFormat = ONE) THEN
   '  Space formated file
   PRINT #iOutF1&, sTmp1$
   PRINT #iOutF1&, sTmp2$
   PRINT #iOutF1&, sTmp3$
   PRINT #iOutF1&, sTmp4$
   PRINT #iOutF2&, sTmp1$
   PRINT #iOutF2&, sTmp2$
   PRINT #iOutF2&, sTmp3$
   PRINT #iOutF2&, sTmp4$
ELSEIF (iOFormat = TWO) THEN
   '  Comma delimited file
   PRINT #iOutF1&, Quote$ + sTmp1$ + Quote$
   PRINT #iOutF1&, Quote$ + sTmp2$ + Quote$
   PRINT #iOutF1&, Quote$ + sTmp3$ + Quote$
   PRINT #iOutF1&, Quote$ + sTmp4$ + Quote$
   PRINT #iOutF2&, Quote$ + sTmp1$ + Quote$
   PRINT #iOutF2&, Quote$ + sTmp2$ + Quote$
   PRINT #iOutF2&, Quote$ + sTmp3$ + Quote$
   PRINT #iOutF2&, Quote$ + sTmp4$ + Quote$
ELSEIF (iOFormat = THREE) THEN
   '  MS Excel file
    stat& = xlsWriteText(sTmp1$, 1, 1, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iOutF1&)
    stat& = xlsWriteText(sTmp2$, 2, 1, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iOutF1&)
    stat& = xlsWriteText(sTmp3$, 3, 1, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iOutF1&)
    stat& = xlsWriteText(sTmp4$, 4, 1, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iOutF1&)

    stat& = xlsWriteText(sTmp1$, 1, 1, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iOutF2&)
    stat& = xlsWriteText(sTmp2$, 2, 1, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iOutF2&)
    stat& = xlsWriteText(sTmp3$, 3, 1, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iOutF2&)
    stat& = xlsWriteText(sTmp4$, 4, 1, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iOutF2&)
END IF
IF (iCriteria <> ZERO) THEN
   PRINT #iOutF4txt&, sTmp1$
   PRINT #iOutF4txt&, sTmp2$
   PRINT #iOutF4txt&, sTmp3$
   PRINT #iOutF4txt&, sTmp4$

   IF (iOFormat = THREE) OR (iCriteria <> ZERO) THEN   ' for now, let's always make a WQC file in xls format
      stat& = xlsWriteText(sTmp1$, 1, 1, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls&)
      stat& = xlsWriteText(sTmp2$, 2, 1, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls&)
      stat& = xlsWriteText(sTmp3$, 3, 1, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls&)
      stat& = xlsWriteText(sTmp4$, 4, 1, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls&)
   END IF
END IF
'if TRUE THEN eXIT function

IF (iDeriv = TRUE) THEN
   '
   ' open a file for derivitive info. and version stamp it
   '
   sOutF3$ = FirstName$(sInFile$) + ".der"
   OPEN sOutF3$ FOR OUTPUT AS iOutF3
   PRINT #iOutF3, sTmp$
END IF
'
' Define critical gill concenration
'
' Now read in from CHESS file and units converted here
' from nmol/g to mol/kg
'CritGill = OKT(NComp) / (CtoM(NComp) * 2)
'CritGill = 5.3E-13 / 2
'CritGill = .00001

' Set BL Concentration
'
'  TODO - instead of hardwiring this, we should have a search function for the BL compartment
'
'OKT(iBL) = .0000178
OKT(iBL) = .0000178 * MassVal(NMass)
IF iLowMass THEN
   OKT(iBL) = OKT(iBL) / 1000000
END IF
IF (iBLM(1) > 10) THEN
   CritBL = DefaultCriticalValue * OKT(iBL) / 1000000!
   'CritBL = CritBL * OKT(iBL) / 1000000! 'commented as of 16/3/1
END IF
'DefaultCriticalValue = CritBL 'commented as of 16/3/1

' -----------------------------------------------------------------------------------
' new to ver 2.40: If appropriate, use a CAT table entry to set critical accumulation
' -----------------------------------------------------------------------------------
IF iHaveCAT THEN
   IF (iCAT_Row > iCAT_NRow) OR (iCAT_Row < 1) THEN
      IF LC50ALL THEN
         PRINT "Reverting to first row of CA Table."
      END IF
      iCAT_Row = 1
   END IF

   DefaultCriticalValue = CAT_CA(iCAT_Row)
   IF (iBLM(1) > 10) THEN
       CritBL = DefaultCriticalValue * OKT(iBL) / 1000000!
   ELSE
       CritBL = DefaultCriticalValue
   END IF
END IF

'T$ = "Outside Top of main loop"
'CALL Pmsg(T$)
'T$ = FORMAT$(NData)
'CALL Pmsg(T$)
'
'
' Begin loop solving chemical problem for each line of input data
'
FOR i = 1 TO NData&
    iObs = i
    IF (iDebug = 01 OR iDebug = 02) THEN CALL DebugLog ("main", "Begin new obs")

   'T$ = "Inside Top of main loop"
   'CALL Pmsg(T$)

   IF NOT iQuiet THEN
      IF (i=1) THEN
         CLS
      END IF
      LOCATE 8, 43
   END IF
   IF iVQuiet THEN
      IF Exist("elapsed") THEN
         KILL "elapsed"
      END IF
   ELSE
      PRINT "Obs "; i; " of "; NData&
   END IF
   IF iPartMe = TRUE THEN
       DISS_Effect =  0
       PART_Effect = 0
       Part_Conc = 0
   END IF


   '
   '
   ' Get next line of input data
   '
   'CALL GetData(iInFile, LC50, iExpert, CriticalValue)
   CALL GetData2 (iInFile, LC50, iExpert1, iExpert2, iExpert3, PFile$, CriticalValue, ErrorCode%, ErrorMsg$)
   IF (iDebug = 01 OR iDebug = 02) THEN CALL DebugLog ("main", "After GetData2")

   IF iCriteria <> ZERO THEN
      IF OKT(iMetal)>0 THEN
          ' metal conc was measured, don't set the censored flag
          AmbientMetal(i) = OKT(iMetal)
          CenFlag(i) = 0
      ELSEIF (OKT(iMetal)<0)AND (OKT(iMetal)>-990) THEN
          ' metal conc is below the reported detection limit,
          ' set the censored flag to 1 (indicating that the actual
          ' value is less than the reported value), and take away
          ' negative sign
          OKT(iMetal) = -OKT(iMetal)
          AmbientMetal(i) = OKT(iMetal)
          CenFlag(i) = 1
      END IF
    ' '
    '  IF iMetalDL = 0 THEN
    '     CenFlag(i) = 0
    '  ELSE
    '     CenFlag(i) = OKT(iMetalDL)
    '  END IF
   END IF 'If iCriteria <> ZERO

   IF (iExpert2 = TRUE) AND (PFile$ <> sCHESS$) THEN
      sCHESS$ = PFile$
      CALL DefineProblem300W(sCHESS$, ErrorCode%, ErrorMsg$)

      ' Set BL Concentration
      OKT(iBL) = .0000178
      IF (SType(iBL)=stSurface) THEN
         'CritBL = CritBL * OKT(iBL) / 1000000! 'commented as of 16/3/1
         CritBL = DefaultCriticalValue * OKT(iBL) / 1000000!
      END IF
      'DefaultCriticalValue = CritBL
   END IF 'if (iExpert2 = TRUE) AND (PFile$ <> sCHESS$)
   IF (iExpert3 = TRUE) THEN
      DISS_A = DISS_A_expert
      PART_HI = PART_HI_expert
   END IF

   'Pmsg(STR$(iComplete))
   IF EchoInput THEN
      ' Just echo all inputs to output files
      ' This is the equivalent of saying each line is not complete
      iComplete = FALSE
   END IF
   'Pmsg(STR$(iComplete))
   'Pmsg(STR$(EchoInput))
   IF (ErrorCode%) THEN
      PRINT ErrorMsg$
      EXIT FUNCTION
   END IF

   IF (iExpert1) THEN
     IF (iBLM(1) > 10) THEN
       CritBL = CriticalValue * OKT(iBL) / 1000000!
     ELSE
       CritBL = CriticalValue
     END IF
     'IF (CriticalValue < 0!) THEN
     '   ErrorCode% = TRUE
     '   ErrorMsg$ = "Invalid critical value in obs " + STR$(i)
     '   print ErrorMsg$
     '   exit function
     'END IF
   ELSE
      IF (iBLM(1) > 10) THEN
         CritBL = DefaultCriticalValue * OKT(iBL) / 1000000!
      ELSE
         CritBL = DefaultCriticalValue
      END IF
   END IF 'IF (iExpert1)

   IF (iPartMe = TRUE) THEN
     ' Calculate the critical BL amount needed to acheive
     ' the target effect
     ' where:
     '   R = Target effect = 1 - exp(a+b*log10(BL-Al))/(1+exp(a+b*log10(x)))
     '   R1 = 1 - R
     '   R2 = R1 / (1 - R1)
     '   Log Crit BL-Al = (ln(R2) - a)/b
'               DISS_A = VAL(TmpS3$)
'            CASE "DISS_B"
'               DISS_B = VAL(TmpS3$)
'            CASE "PART_A"
'               PART_A = VAL(TmpS3$)
'            CASE "PART_B"
'               PART_B = VAL(TmpS3$)
'            CASE "PART_H"
'               PART_H = VAL(TmpS3$)
'            CASE "TARGET_R"
'               TARGET_R = VAL(TmpS3$)
     DIM R1 AS DOUBLE
     DIM R2 AS DOUBLE
     R1 = 1 - TARGET_R
     R2 = R1 / (1 - R1)
     'CritBL = (LOG(R2) - DISS_A)/DISS_B         ' expression based on A', such that A' = A*B
     CritBL = (LOG(R2)+DISS_A*DISS_B )/DISS_B    ' expression based on A
     CritBL = 10^CritBL
     CritBL = CritBL * OKT(iBL) / 1000000!


   END IF 'IF (iPartMe = TRUE)
   IF LC50ALL THEN LC50 = TRUE
   '
   '  Begin loop for derivitive information
   '
   '  Value of NDeriv will determine behavior
   '     NDeriv =
   '       iBL             - LC50 determination with measured chem
   '
   '       j (j=1, iBL-1)  - Determine derivitive for component j
   '
   '       0                 - Done with derivitive info.
   '
   IF (iDeriv = TRUE) THEN
      NDeriv = iBL
   ELSE
      NDeriv = ZERO%
   END IF
   CritBL0 = CritBL
   DO 'LOOP UNTIL NDeriv = ZERO%
      '
      '  Begin loop for single LC50 determination
      '
      IF (LC50 = TRUE) AND (LastLC50# > 0) THEN
         OKT(iMetal) = LastLC50#
      END IF
      '
      '
      ' Set up problem for CHESS
      BLCrit = FALSE
      ii = 1
      ErrorNotConverged = FALSE
      DO 'LOOP UNTIL (NOT LC50) OR BLCrit OR NOT iComplete

         FOR j = 1 TO NComp
            KT(j) = OKT(j)
            'OC(j) = OKT(j) / MassVal(SType(i) + 1)
            OC(j) = OKT(j) / MassVal(SType(j) + 1)   ' edit 2.41 why do we do this?  This will break firstguess.  Should probably only do this on iter = 1 (DO NOT DELETE THIS COMMENT)
         NEXT j
       '  OC(1) = OKT(1)
         IF (iDebug = 01 OR iDebug = 02) THEN
             CALL DebugLog ("main", "After set totals")
             FOR j = 1 TO NComp
               CALL DebugLog ("main", Label(j) + ".OC" + STR$(OC(j)))
               CALL DebugLog ("main", Label(j) + ".KT" + STR$(KT(j)))
             NEXT j
         END IF 'IF (iDebug = 01 OR iDebug = 02)
         'IF TRUE THEN
         IF FALSE THEN
            LOCATE 1, 1
            FOR j = 1 TO NComp
               IF j = 15 THEN
                 DO UNTIL LEN(INKEY$)
                 LOOP
               END IF
               PRINT j;
               PRINT Label(j); " ";
               'PRINT USING "##.##^^^^"; TAB(5); OKT(j)
               'PRINT USING " ##.##^^^^"; OKT(j); KT(j); OC(j);
               PRINT FORMAT$(OKT(J), " 0.0E-##");
               PRINT FORMAT$(KT(j), " 0.0E-##");
               PRINT FORMAT$(OC(j), " 0.0E-##");
               PRINT CType(j)
            NEXT j
            'INPUT j
         END IF 'If FALSE

         'DIM PrnKT(1 TO NComp) AS SINGLE
         '
         'FOR j = 1 TO NComp
         '   KT(j) = OKT(j)
         '   PrnKT(j) = KT(j)
         'NEXT j
         'OC(1) = OKT(1)
         '
         'CALL PrintA(PrnKT(), NComp)
         IF (iComplete = TRUE) THEN
            IF iPartMe = TRUE THEN
               '
               ' for simulations with particulate metal effects,
               '    1 - make sure the toxic metal is Al
               '    2 - on the first iteration, set the Al solubility switch to saturated
               IF (sMetal$<>"Al") THEN
                   ' this isn't an Al simulation, so just ignore the particulate stuff
                   iPartMe = FALSE
               ELSE
                   ' this is Al so let's proceed
                   IF (ii = 1) THEN ' ii is the poorly named iteration count
                       iGlobalAl = 2
                   END IF
               END IF
            END IF 'IF iPartMe = TRUE
            IF (iDebug = 01 OR iDebug = 02) THEN CALL DebugLog ("main", "Before CHESS")
            CALL CHESS(ErrorCode%, ErrorMsg$)
            IF (iDebug = 01 OR iDebug = 02) THEN CALL DebugLog ("main", "After CHESS, Iter" + STR$(ii))

            IF (ErrorCode%) THEN
               PRINT ErrorMsg$
               IF (iDebug = 01 OR iDebug = 02) THEN CALL DebugLog ("main", ErrorMsg$)

               IF (ErrorMsg$ = "Iterations exceeded max.") THEN
                   ErrorNotConverged = TRUE
               ELSE
                 IF (iOFormat = THREE) THEN
                    ' close the spreadsheets
                    stat& = xlsCloseFile(iOutF1&)
                    stat& = xlsCloseFile(iOutF2&)
                 END IF
                 EXIT FUNCTION
               END IF
               ErrorCode% = FALSE
            END IF 'IF (ErrorCode%)

            '
            '
            IF (LC50 = TRUE) THEN
               ' Adjust the critical concentration if necessary, to allow for
               ' increased sensitivity in some environments
               IF (iLA50Adjust = 1) THEN
                  'LA50PARAMS$(1)   = the parameter that will be used to adjust
                  'LA50PARAMS$(2), LA50Slope   = the slope of LA50 to the parameter
                  'LA50PARAMS$(3), LA50Break   = the minimum value of parameter below which LA50 is adjusted
                  'LA50PARAMS$(4), LA50MinFrac   = the minimum fraction of the default LA50 that is allowed after adjustment
                  'iLA50Variable    = the position in OC()where the parameter is found
                  'OC(iLA50Variable)= 1e-5
                  IF LA50PARAMS$(1) = "IONICSTR" THEN
                     iLA50Variable = NSpecies + 2
                  ELSE
                     iLA50Variable = FindSpecies(LA50PARAMS$(1))
                  END IF

                  CritBL = LA50Slope * (OKT(iBL) / 1000000!) * (OC(iLA50Variable) - LA50Break) + CritBL0

                  IF (CritBL > CritBL0) THEN
                      CritBL = CritBL0
                  END IF
                  IF (CritBL < (CritBL0 * LA50MinFrac )) THEN
                      CritBL = CritBL0 * LA50MinFrac
                  END IF
               ELSEIF (iLA50Adjust = 2) THEN
                  'print "Reaching the LA50 adjustment..."
                  ' Adjust the LA50 down, using a relationship such as:
                  '
                  '  LA50PARAMS$(1) is the name of the parameter that will be used to adjust
                  '
                  '  iLA50Variable is the position in OC() where the parameter is found
                  '
                  '  LA50Break is a value in the parameter file that specifies the value
                  '            specified by iLA50Variable, below which an adjustment begins
                  '
                  '  LA50Slope is the slope of the adjustment when OC(iLA50Variable)
                  '            below values of LA50Break
                  '
                  '  The new LA50:
                  '
                  '  Adjustment factor
                  '  AdjFact = LA50Slope((LOG10(OC(iLA50Variable)) -  LOG10(LA50Break))
                  '  LA50_ADJ = LA50 * 10^(AF)
                  '

                  IF LA50PARAMS$(1) = "IONICSTR" THEN
                     iLA50Variable = NSpecies + 2
                  ELSE
                     iLA50Variable = FindSpecies(LA50PARAMS$(1))
                  END IF

                  ' Only adjust the LA50 if it's on the slope below LA50Break
                  IF (OC(iLA50Variable) < LA50Break) THEN
                     AdjFact = LA50Slope * (LOG10(OC(iLA50Variable)) - LOG10(LA50Break))
                  ELSE
                     AdjFact = 0
                  END IF
                  CritBL = CritBL0 * (10^AdjFact)

               END IF 'IF (iLA50Adjust = 1)

               ' Test to see if BL Metal conc equals the critical conc
               BLMConc! = 0!
               FOR j = 1 TO NBLM
                  IF (iBLMTox(j)) THEN
                     BLMConc! = BLMConc! + OC(iBLM(j))
                  END IF
               NEXT j
               Diff = BLMConc! - CritBL
               Fract = ABS(Diff) / CritBL
               IF (iDebug = 01 OR iDebug = 02) THEN CALL DebugLog ("main", "After Diff, OC" + Label(iMetal) + " " + STR$(OC(iMetal)) )
               IF (iDebug = 01 OR iDebug = 02) THEN CALL DebugLog ("main", "After Diff, BLMConc" + STR$(BLMConc!))
               IF (iDebug = 01 OR iDebug = 02) THEN CALL DebugLog ("main", "After Diff, Diff" + STR$(Diff))
               IF (iDebug = 01 OR iDebug = 02) THEN CALL DebugLog ("main", "After Diff, CritBL" + STR$(CritBL))
               IF (iDebug = 01 OR iDebug = 02) THEN CALL DebugLog ("main", "After Diff, Fract" + STR$(Fract))
               IF (iDebug = 01 OR iDebug = 02) THEN CALL DebugLog ("main", "After Diff, ConvCrit" + STR$(ConvCrit))

               IF NOT iQuiet THEN
                  LOCATE 12, 43
                  PRINT "LC50 % Error: ";
                  'PRINT USING "####.###"; FRAC * 100
                  PRINT FORMAT$(Fract * 100, "####.###")
                  PRINT "Diff"; Diff
                  PRINT "CritBL"; CritBL
               END IF 'IF NOT iQuiet
               IF (iPartMe = TRUE) AND (ii = 1) THEN
                  '
                  ' for simulations where particulate effects will be considered,
                  ' and on the first iteration only, the iGlobalAl switch will have
                  ' been set to perform a saturation Al calculation.  Now we can
                  ' check to see whether the amount of Al on the biotic ligand, is
                  ' is excess of the amount needed for all desired effects to be
                  ' attributed to dissolved metal.  If not, then a particulate effect
                  ' needs to be added.
                  '
                  ' does OKT mean anything at this point?
                  OKT(iMetal) = KT(iMetal)

                  IF Diff < 0 THEN
                      PRINT "Dissolved and precipitated Al required for target effect"
                      ' There is not sufficient solubility for all effects to be dissolved
                      '
                      ' estimate how much of an effect can be attributable to solubile metal
                      '
                      ' remaining effect from particulate metal
                    'DIM DISS_Effect AS GLOBAL DOUBLE   ' portion of target effect acheivable by dissolved metal
                    'DIM Part_Conc AS GLOBAL DOUBLE     ' particulate metal concentration needed if Diss_Effect < TARGET_R
                      'R_diss = 1 - exp(a+b*log10(BL-Al))/(1+exp(a+b*log10(x)))
                      R1 = BLMConc! * 1000000! / KT(iBL)
                      R2 = (KT(FindSpecies%("Ca")) + KT(FindSpecies%("Mg"))) * 1000 * 100 ' hardness in mg/L
                      'DISS_Effect =  1 - EXP(DISS_A+DISS_B*LOG10(R1))/(1+EXP(DISS_A+DISS_B*LOG10(R1)))
                      DISS_Effect = 1/(1+EXP(DISS_B*LOG10(R1) - DISS_A*DISS_B))

                      IF (iDebug = 01 OR iDebug = 09) THEN CALL DebugLog ("main", "After solubility calc, diss effect " + STR$(DISS_Effect))
                      '
                      ' R_part = 1-((1-0.1)/(1-R_diss))
                      PART_Effect = 1-((1-TARGET_R)/(1-DISS_Effect))

                      '
                      ' particulate Al in ug/L: 10^((ln(1/(1/(1-R_part) - 1)) - a - 0.025*hardness)/b)
                      PART_A = LOG10(R2 * PART_HM + PART_HI)
                      PART_A = (LOG10(R2) * PART_HM) + PART_HI
                      Part_Conc = 10^((LOG(1/(1/(1-PART_Effect) - 1)) - PART_A - PART_H*R2)/PART_B)
                      Part_Conc = 10^((LOG(1/(1/(1-PART_Effect) - 1)) + PART_A*PART_B)/PART_B)
                      '
                      ' signal to the convergence criteria that we're done
                      Fract = ConvCrit  / 10

'                      PRINT " Diff"; Diff
'                      PRINT "Crit BL"; CritBL
'                      PRINT "Diss effect"; DISS_Effect
'                      PRINT "Part effect"; PART_Effect
'                      PRINT "Part Conc"; Part_Conc
                  ELSE
                      ' There IS sufficient solubility for all effects to be dissolved
                      ' so remove the solubility contraint and solve normally
'                      PRINT "Diss Al sufficient"
'                      PRINT "Crit BL"; CritBL
'                      PRINT " Diff"; Diff
                      iGlobalAl = 0
                      DISS_Effect = TARGET_R
                  END IF 'IF Diff < 0
              END IF 'IF (iPartMe = TRUE) AND (ii = 1)

               IF (Fract < ConvCrit) OR (ErrorNotConverged = TRUE) THEN
                  '
                  ' We're done
                  BLCrit = TRUE
                  LastLC50# = OKT(iMetal)
               ELSE
                  '
                  ' Update the total copper concentration
                  SELECT CASE iAdjust
                     CASE ONE
                        CALL AdjustCu(Diff, ii)
                     CASE TWO
                        CALL AdjustCu2(Diff, ii)
                     CASE THREE
                        CALL AdjustCu3(Diff, ii)
                     CASE FOUR
                        CALL AdjustCu4(Diff, ii)
                     CASE ELSE
                        CALL AdjustCu4(Diff, ii)
                  END SELECT
'                      PRINT "OKT(Al)"; OKT(iMetal)
'                      PRINT "KT(Al)"; KT(iMetal)
'                      PRINT "OC(Al)"; OC(iMetal)

                  'PRINT "Iter ", ii
                  ii = ii + 1
                  IF (ii > MaxECxIter) THEN
                      CALL MessageFile(iEC_NonCritical_NoPost, iObs, "Toxicity endpoint not determined: Iterations exceeded max.")
                  END IF
               END IF'IF (Fract < .001)
            END IF 'IF (LC50 = TRUE)
         END IF 'IF (iComplete = TRUE)
      LOOP UNTIL (NOT LC50) OR BLCrit OR NOT iComplete OR (ii > MaxECxIter)

      IF (iDeriv = FALSE) OR (NDeriv = iBL) THEN
         '
         ' Write results with measured chem. to output files
         '
         'CALL WriteSimple(iOutF1, i, LC50)
         'CALL WriteDetailed(iOutF2&, i)
         IF (iOFormat = ONE) THEN
            CALL WriteSimple2(iOutF1&, i, LC50, iComplete)
            CALL WriteDetailed2(iOutF2&, i, iComplete)
         ELSEIF (iOFormat = TWO) THEN
            CALL WriteSimpleComma(iOutF1&, i, LC50, iComplete)
            CALL WriteDetailedComma(iOutF2&, i, iComplete)
         ELSEIF (iOFormat = THREE) THEN
            CALL WriteSimpleXLS(iOutF1&, i, LC50, iComplete)
            CALL WriteDetailedXLS(iOutF2&, i, iComplete)
         END IF
         IF iCriteria <> ZERO THEN
            HC5(i) = OKT(iMetal)
            'PRINT HC5(i)
            CALL WriteWQC_Acute_and_Chronic(iOutF4txt&, iOutF4xls&, i, NData&)
            'SELECT CASE iCriteria
            '    CASE ONE
            '        CALL WriteWQCFile(iOutF4txt&, iOutF4xls&, i, NData&)
            '    CASE TWO
            '        CALL WriteWQCFileAcute(iOutF4txt&, iOutF4xls&, i, NData&)
            '    CASE THREE
            '        CALL WriteWQCFileChronic(iOutF4txt&, iOutF4xls&, i, NData&)
            'END SELECT
         END IF 'IF iCriteria <> ZERO
      END IF 'if (iDeriv = FALSE) OR (NDeriv = iBL)

      IF (iCalcLA50) THEN
         '
         '  Write out either LA50s or LC50s to a file for access by other programs
         '
         IF LC50ALL THEN
             '
             ' In toxicity mode, write out the LC50
             '
             Tmp! = CSNG(OKT(iMetal))
             IF iPartMe THEN
                 Tmp! = Tmp! + Part_Conc /(1000000*26.981)
             END IF
         ELSE
             '
             ' In speciation mode, write out the LA50
             '
            'T$ = "writing to blm_la50.txt"
            'CALL Pmsg(T$)
            Tmp! = 0!
            FOR j = 1 TO NBLM
              IF (SType(iBLM(1))=stSurface) THEN
                Tmp! = Tmp! + OC(iBLM(j)) * 1000000! / KT(iBL)
              ELSE
                Tmp! = Tmp! + OC(iBLM(j))
              END IF
            NEXT j
        END IF 'if LC50ALL
        IF NOT iComplete THEN
          Tmp! = -999.
        END IF
        IF iCalcLA50Log THEN
            Tmp! = LOG10(Tmp!)
        END IF
        IF iCalcLA50Al AND NOT LC50ALL THEN
          '
          ' For precipitated Al runs in speciation mode, output DISS_A and PART_HI
          '
            ' Use already calculated BL_metal to calculate DISS_A
            Tmp1 =  LOG10(Tmp!)-(1/DISS_B)* LOG((1/TARGET_R)-1)
            IF NOT iComplete THEN
                Tmp1 = -999.
            END IF
            PRINT #iLAFile, FORMAT$(i, "* ####  ");
            PRINT #iLAFile, FORMAT$(Tmp1, "#.0000E+##  ");
            '
            ' Calc PART_A assuming all of ECx is due to precip Al
            '   PART_A = LOG10(ECx)-(1/B)*LN((1/(Rx))-1)
            Tmp1 = LOG10(OKT(iMetal)*27.1*1E6)-(1/PART_B)* LOG((1/TARGET_R)-1)
            'Tmp1 = 10^Tmp1
            ' Calculate hardness in mg/L
            Tmp2 = (KT(FindSpecies%("Ca")) + KT(FindSpecies%("Mg"))) * 1000 * 100
            ' Calculate PART_HI
            ' PART_HI = PART_A - Log10(H)*PART_HM
            Tmp3 = Tmp1 - LOG10(Tmp2)*PART_HM
            '
'            ' Calc PART_A assuming all of ECx is due to precip Al
'            '   PART_A = LOG10(ECx)-(1/B)*LN((1/(Rx))-1)
'            Tmp1 = LOG10(OKT(iMetal))-(1/PART_B)* LOG((1/TARGET_R)-1)
'            ' Calculate hardness in mg/L
'            Tmp2 = (KT(FindSpecies%("Ca")) + KT(FindSpecies%("Mg"))) * 1000 * 100
'            ' Calculate PART_HI
'            ' PART_HI = PART_A - Log10(H)*PART_HM
'            Tmp3 = Tmp1 - LOG10(Tmp2)*PART_HM
            IF NOT iComplete THEN
                Tmp3 = -999.
            END IF
            PRINT #iLAFile, FORMAT$(Tmp3, "#.0000E+##")
        ELSE
            '
            ' For non-Al simulations, just write out the BL_metal
            PRINT #iLAFile, FORMAT$(i, "* ####  ");
            PRINT #iLAFile, FORMAT$(Tmp!, "#.0000E+##")
        END IF 'if iCalcLA50Al and not LC50ALL
      END IF 'if iCalcLA50

      IF (iDeriv = TRUE) THEN
         '
         ' Store adjusted LC50 values for derivitive analysis
         '
         DLC50(NDeriv) = LastLC50#
         IF (NDeriv < iBL) THEN
            OKT(NDeriv) = KeepKT#
         END IF
         NDeriv = NDeriv - 1
         IF (NDeriv = iMetal) THEN
            NDeriv = NDeriv - 1
         END IF
         IF (NDeriv > 0) THEN
            IF (NDeriv = 1) THEN

            END IF
            '
            ' Save KT value to restore later
            KeepKT# = OKT(NDeriv)
            OKT(NDeriv) = KeepKT * 1.1#
            '
            ' Print current derivitive estimate on screen
            IF NOT iQuiet THEN
               LOCATE 13, 43
               TmpS$ = "d[LC50] / d[" + Label$(NDeriv) + "]"
               TmpS$ = TmpS$ + SPACE$(22 - LEN(TmpS$))
               PRINT TmpS$;
            END IF
         ELSE
            '
            ' Erase the last derivitive label
            IF NOT iQuiet THEN
               LOCATE 13, 43
               TmpS$ = SPACE$(22)
               PRINT TmpS$;
            END IF 'if not iQuiet
         END IF 'if NDeriv > 0
      END IF 'if iDeriv = TRUE
   LOOP UNTIL NDeriv = ZERO%
   '
   '  Output Derivitive info if calculated
   IF (iDeriv = TRUE) THEN
      CALL WriteDeriv(iOutF3, i, KT(), DLC50(), iBL)
   END IF
NEXT i

IF iCriteria <> ZERO THEN
    'CALL WriteWQCFile(iOutF4txt&, iOutF4xls&, Iter, NData)
END IF

StopTime = TIMER
TimeElapsed = StopTime - StartTime
IF (TimeElapsed < 0) THEN
   TimeElapsed = (60 * 60 * 24) - StartTime - StopTime
END IF

IF NOT iVQuiet THEN
   PRINT "Seconds elapsed: "; TimeElapsed
ELSE
   i = FREEFILE
   OPEN "elapsed" FOR OUTPUT AS #i
   PRINT #i, TimeElapsed
   CLOSE #i
END IF

IF (iOFormat = THREE) THEN
    ' close the spreadsheets
    stat& = xlsCloseFile(iOutF1&)
    stat& = xlsCloseFile(iOutF2&)
END IF
IF (iCriteria <> ZERO) THEN
    stat& = xlsCloseFile(iOutF4xls&)
END IF

END FUNCTION

SUB AdjustCu2 (Diff, Iter)
STATIC CuHigh AS DOUBLE
STATIC CuLow AS DOUBLE
STATIC BLHigh AS DOUBLE
STATIC BLLow AS DOUBLE
STATIC Bracket AS INTEGER
STATIC iDebugOld AS INTEGER

' Find a metal concentration that yields a Diff = 0
'
' Assume:
'
' Diff = a * b * c / (1 + a * c)   -  g
'
'    c = metal conc (dissolved or free, could try either)
'    b = max BL conc, constant for each sim
'    a = conditional stability constant
'    g = critical BL concentration
'
'        a = (Diff + g) / (c * (b - Diff - g))
'
' dd   = derivative of Diff with respect to c
'
'    dd = a * (b + a * c - a * b * c) / (1 + a * c * (2 + a * c))
'
' cnew = new guess for c
'
'    cnew = c - Diff / dd
'
DIM a AS DOUBLE
DIM b AS DOUBLE
DIM c AS DOUBLE
DIM d AS DOUBLE
DIM ddc AS DOUBLE
DIM g AS DOUBLE

g = CritBL * 1000000! / OKT(iBL)
b = 1000000! * OKT(iBL)
c = OKT(iMetal)
d = Diff * 1000000! / OKT(iBL)
a = (Diff + g) / (c * (b  - Diff  - g))
a = (d + g) / (c * (b  - d  - g))
ddc = a * (b + a * c - a * b * c) / (1 + a * c * (2 + a * c))
cnew = c - d / ddc

'DIM adT$
'adT$ = "Diff" & STR$(d)
'CALL Pmsg(adT$)
'adT$ = "a" & STR$(a)
'CALL Pmsg(adT$)
'adT$ = "b" & STR$(b)
'CALL Pmsg(adT$)
'adT$ = "c" & STR$(c)
'CALL Pmsg(adT$)
'adT$ = "ddc" & STR$(ddc)
'CALL Pmsg(adT$)
'adT$ = "g" & STR$(g)
'CALL Pmsg(adT$)
'adT$ = "cnew" & STR$(cnew)
'CALL Pmsg(adT$)

DIM Slope AS SINGLE
DIM SignDiff AS INTEGER

OldTCu = OKT(iMetal)
SignDiff = Sign(Diff)

'DIM adT$
'adT$ = STR$(KT(iBL)) + STR$(OC(iBLM(1)) * 1000000! / KT(iBL))
'CALL Pmsg(adT$)
'
'
' Check the current value of total metal and the BL-metal produced
' against the previous high and low values
IF Iter = 1 THEN
   '
   ' This is the first iteration for this sample,
   ' reset the value of all variables
   CuLow = OKT(iMetal)
   CuHigh = OKT(iMetal)
   BLLow = Diff
   BLHigh = Diff
   Bracket = FALSE

   ' BEGIN DEBUG
   'iDebugOld = FREEFILE
   'OPEN "debug.prn" FOR OUTPUT AS #iDebugOld
   'PRINT #iDebugOld, " TCu      "; " Next TCu    "; " Diff     "; " CuHigh   "; " BLHigh "; " CuLow    "; " BLLow"
   ' BEGIN DEBUG
   'INPUT j
END IF

' Do we bracket yet?
IF ((NOT Bracket) AND (SignDiff <> Sign(BLLow))) THEN
   '
   ' Great!, we bracket now set the flag so won't keep checking
   Bracket = TRUE
END IF

IF Bracket THEN
'
' We bracket, so use Newton and
' false position if Newton fails
' BUT - if iteration number gets above 10, assume the function
' has a bad shape for slope-based methods.  For these try bisection.
   IF (SignDiff < 0) THEN
      ' The last value replaces the low value
      BLLow = Diff
      CuLow = OKT(iMetal)
   ELSE
      ' The last value replaces the high value
      BLHigh = Diff
      CuHigh = OKT(iMetal)
   END IF
'      OKT(iMetal) = CuLow + (CuHigh - CuLow) * BLLow / (BLLow - BLHigh)
'      adT$ = "cfp" & STR$(OKT(iMetal))
'      CALL Pmsg(adT$)
   IF (Iter > 10) THEN
      OKT(iMetal) = (CuLow + CuHigh) / 2
   ELSE
      OKT(iMetal) = c - d / ddc
      IF (OKT(iMetal) < CuLow) OR (OKT(iMetal) > CuHigh) THEN
         ' New estimate is bad, revert to false position
         OKT(iMetal) = CuLow + (CuHigh - CuLow) * BLLow / (BLLow - BLHigh)
      END IF
   END IF
ELSE
'
' We don't bracket, keep expanding the range until we do
   IF (SignDiff < 0) THEN
      ' The last value replaces the high value
      BLLow = BLHigh
      CuLow = CuHigh
      BLHigh = Diff
      CuHigh = OKT(iMetal)
   ELSE
      ' The last value replaces the low value
      BLHigh = BLLow
      CuHigh = CuLow
      BLLow = Diff
      CuLow = OKT(iMetal)
   END IF
   '
   ' Use the value of Diff to determine the direction
   ' that we have to move to get the next best guess
   ' for Total metal.  When Diff < 0, increase the value of KT
   ' and when Diff > 0, decrease the value of KT
   '
   SignDiff = Sign(Diff)
   OKT(iMetal) = OKT(iMetal) * (10 ^ -SignDiff)
END IF

' BEGIN DEBUG
'PRINT #iDebugOld, USING "##.###^^^^"; OldTCu; OKT(iMetal); Diff; CuHigh; BLHigh; CuLow; BLLow
' END DEBUG

END SUB

SUB AdjustCu3 (Diff, Iter)
STATIC CuHigh AS DOUBLE
STATIC CuLow AS DOUBLE
STATIC BLHigh AS DOUBLE
STATIC BLLow AS DOUBLE
STATIC Bracket AS INTEGER
STATIC iDebugOld AS INTEGER
DIM Diffs(10) AS STATIC DOUBLE
DIM Mets(10) AS STATIC DOUBLE
'STATIC Diffss(1 TO 10) AS DOUBLE
'STATIC Metss(1 TO 10) AS DOUBLE

DIM id AS INTEGER
id = FALSE

' Find a metal concentration that yields a Diff = 0
'
' Assume:
'
' Diff = a * b * c / (1 + a * c)   -  g
'
'    c = metal conc (dissolved or free, could try either)
'    b = max BL conc, constant for each sim
'    a = conditional stability constant
'    g = critical BL concentration
'
'        a = (Diff + g) / (c * (b - Diff - g))
'
' ddc   = derivative of Diff with respect to c
'
'    ddc = a * (b + a * c - a * b * c) / (1 + a * c * (2 + a * c))
'
' cnew = new guess for c
'
'    cnew = c - Diff / ddc
'
DIM a AS DOUBLE
DIM b AS DOUBLE
DIM c AS DOUBLE
DIM d AS DOUBLE
DIM ddc AS DOUBLE
DIM g AS DOUBLE

g = CritBL * 1000000! / OKT(iBL)
b = 1000000! * OKT(iBL)
c = OKT(iMetal)
d = Diff * 1000000! / OKT(iBL)
a = (Diff + g) / (c * (b  - Diff  - g))
a = (d + g) / (c * (b  - d  - g))
ddc = a * (b + a * c - a * b * c) / (1 + a * c * (2 + a * c))
cnew = c - d / ddc

'DIM adT$
'adT$ = "Diff" & STR$(d)
'CALL Pmsg(adT$)
'adT$ = "a" & STR$(a)
'CALL Pmsg(adT$)
'adT$ = "b" & STR$(b)
'CALL Pmsg(adT$)
'adT$ = "c" & STR$(c)
'CALL Pmsg(adT$)
'adT$ = "ddc" & STR$(ddc)
'CALL Pmsg(adT$)
'adT$ = "g" & STR$(g)
'CALL Pmsg(adT$)
'adT$ = "cnew" & STR$(cnew)
'CALL Pmsg(adT$)
DIM i AS INTEGER
DIM j AS INTEGER
DIM n AS INTEGER
DIM Slope AS SINGLE
DIM SignDiff AS INTEGER
DIM mn AS INTEGER
mn = 4

'
' Insert this pair of metal conc. and Diff into the arrays
' such that the pairs are in ascending order (determined by
' the value of Diff).

IF Iter = 1 THEN
   i = 1
   FOR j = 1 TO mn
      Diffs(j) = 0#
      Mets(j) = 0#
   NEXT j
ELSE
   n = Iter
   IF n > mn THEN
      n = mn
      '
      ' The maximum number of elements have already been found,
      ' find the worst of the existing pairs to throw out
      j = 1
      DO
         d = ABS(Diffs(j))
         IF d > ddc THEN k = j
         j = j + 1
      LOOP UNTIL j = mn
      '
      ' Thow out element k
      FOR j = k TO n - 1
         Diffs(j) = Diffs(j + 1)
         Mets(j) = Mets(j + 1)
      NEXT j
   END IF
   '
   ' Find out where to insert the latest pair
   i = 1
   DO
      IF Diffs(i) > Diff THEN EXIT DO
      IF i > Iter THEN EXIT DO
      i = i + 1
   LOOP UNTIL i = mn
   '
   ' Clear out the insertion point
   FOR j = n TO i STEP -1
      Diffs(j) = Diffs(j - 1)
      Mets(j) = Mets(j - 1)
   NEXT j
END IF
'
' Do the insertion
Diffs(i) = Diff * 1000000! / KT(iBL)
Mets(i) = OKT(iMetal) * 1.0E6
'Diffs(i) = Diff
'Mets(i) = OKT(iMetal)

IF id THEN
LOCATE 1,1
FOR j = 1 TO mn
    PRINT Diffs(j), Mets(j)
NEXT j

PRINT "Iter = "+STR$(Iter)
PRINT "Diff = "+STR$(Diff)
PRINT "i = "+STR$(i)+ ", n = "+STR$(n)
END IF

IF id THEN
d = 0#
IF n > 3 THEN
   CALL CSFIT1(n, Diffs(), Mets(), d, cnew)
   'CALL Pmsg("Predicted c = "+STR$(cnew))
   PRINT "Predicted c = "+STR$(cnew)
ELSE
   PRINT "                     "
END IF
END IF

OldTCu = OKT(iMetal)
SignDiff = Sign(Diff)

'DIM adT$
'adT$ = STR$(KT(iBL)) + STR$(OC(iBLM(1)) * 1000000! / KT(iBL))
'CALL Pmsg(adT$)
'
'
' Check the current value of total metal and the BL-metal produced
' against the previous high and low values
IF Iter = 1 THEN
   '
   ' This is the first iteration for this sample,
   ' reset the value of all variables
   CuLow = OKT(iMetal)
   CuHigh = OKT(iMetal)
   BLLow = Diff
   BLHigh = Diff
   Bracket = FALSE

   ' BEGIN DEBUG
   'iDebugOld = FREEFILE
   'OPEN "debug.prn" FOR OUTPUT AS #iDebugOld
   'PRINT #iDebugOld, " TCu      "; " Next TCu    "; " Diff     "; " CuHigh   "; " BLHigh "; " CuLow    "; " BLLow"
   ' BEGIN DEBUG
   'INPUT j
END IF

' Do we bracket yet?
IF ((NOT Bracket) AND (SignDiff <> Sign(BLLow))) THEN
   '
   ' Great!, we bracket now set the flag so won't keep checking
   Bracket = TRUE
END IF

IF Bracket THEN
'
' We bracket, so use Newton and
' false position if Newton fails
' BUT - if iteration number gets above 10, assume the function
' has a bad shape for slope-based methods.  For these try bisection.
   IF (SignDiff < 0) THEN
      ' The last value replaces the low value
      BLLow = Diff
      CuLow = OKT(iMetal)
   ELSE
      ' The last value replaces the high value
      BLHigh = Diff
      CuHigh = OKT(iMetal)
   END IF
'      OKT(iMetal) = CuLow + (CuHigh - CuLow) * BLLow / (BLLow - BLHigh)
'      adT$ = "cfp" & STR$(OKT(iMetal))
'      CALL Pmsg(adT$)
   IF (Iter > 10) THEN
      OKT(iMetal) = (CuLow + CuHigh) / 2
   ELSE
      d = 0#
      IF n > 3 THEN
         CALL CSFIT1(n, Diffs(), Mets(), d, cnew)
      ELSE
         g = CritBL * 1000000! / OKT(iBL)
         b = 1000000! * OKT(iBL)
         c = OKT(iMetal)
         d = Diff * 1000000! / OKT(iBL)
         a = (Diff + g) / (c * (b  - Diff  - g))
         a = (d + g) / (c * (b  - d  - g))
         ddc = a * (b + a * c - a * b * c) / (1 + a * c * (2 + a * c))
         cnew = c - d / ddc
      END IF
      OKT(iMetal) = cnew
      IF (OKT(iMetal) < CuLow) OR (OKT(iMetal) > CuHigh) THEN
         ' New estimate is bad, revert to false position
         OKT(iMetal) = CuLow + (CuHigh - CuLow) * BLLow / (BLLow - BLHigh)
      END IF
   END IF
ELSE
'
' We don't bracket, keep expanding the range until we do
   IF (SignDiff < 0) THEN
      ' The last value replaces the high value
      BLLow = BLHigh
      CuLow = CuHigh
      BLHigh = Diff
      CuHigh = OKT(iMetal)
   ELSE
      ' The last value replaces the low value
      BLHigh = BLLow
      CuHigh = CuLow
      BLLow = Diff
      CuLow = OKT(iMetal)
   END IF
   '
   ' Use the value of Diff to determine the direction
   ' that we have to move to get the next best guess
   ' for Total metal.  When Diff < 0, increase the value of KT
   ' and when Diff > 0, decrease the value of KT
   '
   SignDiff = Sign(Diff)
   OKT(iMetal) = OKT(iMetal) * (10 ^ -SignDiff)
END IF

' BEGIN DEBUG
'PRINT #iDebugOld, USING "##.###^^^^"; OldTCu; OKT(iMetal); Diff; CuHigh; BLHigh; CuLow; BLLow
' END DEBUG

END SUB

SUB AdjustCu4 (Diff, Iter)
STATIC CuHigh AS DOUBLE
STATIC CuLow AS DOUBLE
STATIC BLHigh AS DOUBLE
STATIC BLLow AS DOUBLE
STATIC Bracket AS INTEGER
STATIC iDebugOld AS INTEGER
DIM Diffs(10) AS STATIC DOUBLE
DIM Mets(10) AS STATIC DOUBLE
'STATIC Diffss(1 TO 10) AS DOUBLE
'STATIC Metss(1 TO 10) AS DOUBLE

DIM id AS INTEGER
id = FALSE

' Find a metal concentration that yields a Diff = 0
'
' Assume:
'
' Diff = a * b * c / (1 + a * c)   -  g
'
'    c = metal conc (dissolved or free, could try either)
'    b = max BL conc, constant for each sim
'    a = conditional stability constant
'    g = critical BL concentration
'
'        a = (Diff + g) / (c * (b - Diff - g))
'
' ddc   = derivative of Diff with respect to c
'
'    ddc = a * (b + a * c - a * b * c) / (1 + a * c * (2 + a * c))
'
' cnew = new guess for c
'
'    cnew = c - Diff / ddc
'
DIM a AS DOUBLE
DIM b AS DOUBLE
DIM c AS DOUBLE
DIM d AS DOUBLE
DIM ddc AS DOUBLE
DIM g AS DOUBLE
DIM cnew AS DOUBLE
DIM dcnew AS DOUBLE

g = CritBL * 1000000! / OKT(iBL)
b = 1000000! * OKT(iBL)
c = OKT(iMetal)
d = Diff * 1000000! / OKT(iBL)
a = (Diff + g) / (c * (b  - Diff  - g))
a = (d + g) / (c * (b  - d  - g))
ddc = a * (b + a * c - a * b * c) / (1 + a * c * (2 + a * c))
cnew = c - d / ddc

'DIM adT$
'adT$ = "Diff" & STR$(d)
'CALL Pmsg(adT$)
'adT$ = "a" & STR$(a)
'CALL Pmsg(adT$)
'adT$ = "b" & STR$(b)
'CALL Pmsg(adT$)
'adT$ = "c" & STR$(c)
'CALL Pmsg(adT$)
'adT$ = "ddc" & STR$(ddc)
'CALL Pmsg(adT$)
'adT$ = "g" & STR$(g)
'CALL Pmsg(adT$)
'adT$ = "cnew" & STR$(cnew)
'CALL Pmsg(adT$)
DIM i AS INTEGER
DIM j AS INTEGER
DIM n AS INTEGER
DIM Slope AS SINGLE
DIM SignDiff AS INTEGER
DIM mn AS INTEGER
mn = 4

'
' Insert this pair of metal conc. and Diff into the arrays
' such that the pairs are in ascending order (determined by
' the value of Diff).

IF Iter = 1 THEN
   i = 1
   FOR j = 1 TO mn
      Diffs(j) = 0#
      Mets(j) = 0#
   NEXT j
ELSE
   n = Iter
   IF n > mn THEN
      n = mn
      '
      ' The maximum number of elements have already been found,
      ' find the worst of the existing pairs to throw out
      j = 1
      DO
         d = ABS(Diffs(j))
         IF d > ddc THEN k = j
         j = j + 1
      LOOP UNTIL j = mn
      '
      ' Thow out element k
      FOR j = k TO n - 1
         Diffs(j) = Diffs(j + 1)
         Mets(j) = Mets(j + 1)
      NEXT j
   END IF
   '
   ' Find out where to insert the latest pair
   i = 1
   DO
      IF Diffs(i) > Diff THEN EXIT DO
      IF i > Iter THEN EXIT DO
      i = i + 1
   LOOP UNTIL i = mn
   '
   ' Clear out the insertion point
   FOR j = n TO i STEP -1
      Diffs(j) = Diffs(j - 1)
      Mets(j) = Mets(j - 1)
   NEXT j
END IF
'
' Do the insertion
Diffs(i) = Diff * 1000000! / KT(iBL)
Mets(i) = OKT(iMetal) * 1.0E6
'Diffs(i) = Diff
'Mets(i) = OKT(iMetal)

IF id THEN
LOCATE 1,1
FOR j = 1 TO mn
    PRINT Diffs(j), Mets(j)
NEXT j

PRINT "Iter = "+STR$(Iter)
PRINT "Diff = "+STR$(Diff)
PRINT "i = "+STR$(i)+ ", n = "+STR$(n)

d = 0#
IF n > 3 THEN
   'CALL CSFIT1(n, Diffs(), Mets(), d, cnew)
   'CALL Pmsg("Predicted c = "+STR$(cnew))
   CALL POLINT(Diffs(), Mets(), n, d, cnew, dcnew)
   PRINT "Predicted c = "+STR$(cnew)
ELSE
   PRINT "                     "
END IF

END IF

OldTCu = OKT(iMetal)
SignDiff = Sign(Diff)

'DIM adT$
'adT$ = STR$(KT(iBL)) + STR$(OC(iBLM(1)) * 1000000! / KT(iBL))
'CALL Pmsg(adT$)
'
'
' Check the current value of total metal and the BL-metal produced
' against the previous high and low values
IF Iter = 1 THEN
   '
   ' This is the first iteration for this sample,
   ' reset the value of all variables
   CuLow = OKT(iMetal)
   CuHigh = OKT(iMetal)
   BLLow = Diff
   BLHigh = Diff
   Bracket = FALSE

   ' BEGIN DEBUG
   'iDebugOld = FREEFILE
   'OPEN "debug.prn" FOR OUTPUT AS #iDebugOld
   'PRINT #iDebugOld, " TCu      "; " Next TCu    "; " Diff     "; " CuHigh   "; " BLHigh "; " CuLow    "; " BLLow"
   ' BEGIN DEBUG
   'INPUT j
END IF

' Do we bracket yet?
IF ((NOT Bracket) AND (SignDiff <> Sign(BLLow))) THEN
   '
   ' Great!, we bracket now set the flag so won't keep checking
   Bracket = TRUE
END IF

IF Bracket THEN
'
' We bracket, so use Newton and
' false position if Newton fails
' BUT - if iteration number gets above 10, assume the function
' has a bad shape for slope-based methods.  For these try bisection.
   IF (SignDiff < 0) THEN
      ' The last value replaces the low value
      BLLow = Diff
      CuLow = OKT(iMetal)
   ELSE
      ' The last value replaces the high value
      BLHigh = Diff
      CuHigh = OKT(iMetal)
   END IF
'      OKT(iMetal) = CuLow + (CuHigh - CuLow) * BLLow / (BLLow - BLHigh)
'      adT$ = "cfp" & STR$(OKT(iMetal))
'      CALL Pmsg(adT$)
   IF (Iter > 10) THEN
      OKT(iMetal) = (CuLow + CuHigh) / 2
   ELSE
      d = 0#
      IF n > 3 THEN
         CALL POLINT(Diffs(), Mets(), n, d, cnew, dcnew)
      ELSE
         g = CritBL * 1000000! / OKT(iBL)
         b = 1000000! * OKT(iBL)
         c = OKT(iMetal)
         d = Diff * 1000000! / OKT(iBL)
         a = (Diff + g) / (c * (b  - Diff  - g))
         a = (d + g) / (c * (b  - d  - g))
         ddc = a * (b + a * c - a * b * c) / (1 + a * c * (2 + a * c))
         cnew = c - d / ddc
      END IF
      OKT(iMetal) = cnew
      OKT(iMetal) = cnew * 1e-6
      IF (OKT(iMetal) < CuLow) OR (OKT(iMetal) > CuHigh) THEN
         ' New estimate is bad, revert to false position
         OKT(iMetal) = CuLow + (CuHigh - CuLow) * BLLow / (BLLow - BLHigh)
      END IF
   END IF
ELSE
'
' We don't bracket, keep expanding the range until we do
   IF (SignDiff < 0) THEN
      ' The last value replaces the high value
      BLLow = BLHigh
      CuLow = CuHigh
      BLHigh = Diff
      CuHigh = OKT(iMetal)
   ELSE
      ' The last value replaces the low value
      BLHigh = BLLow
      CuHigh = CuLow
      BLLow = Diff
      CuLow = OKT(iMetal)
   END IF
   '
   ' Use the value of Diff to determine the direction
   ' that we have to move to get the next best guess
   ' for Total metal.  When Diff < 0, increase the value of KT
   ' and when Diff > 0, decrease the value of KT
   '
   SignDiff = Sign(Diff)
   OKT(iMetal) = OKT(iMetal) * (10 ^ -SignDiff)
END IF

' BEGIN DEBUG
'PRINT #iDebugOld, USING "##.###^^^^"; OldTCu; OKT(iMetal); Diff; CuHigh; BLHigh; CuLow; BLLow
' END DEBUG

END SUB


SUB AdjustCu (Diff, Iter)
STATIC CuHigh AS DOUBLE
STATIC CuLow AS DOUBLE
STATIC BLHigh AS DOUBLE
STATIC BLLow AS DOUBLE
STATIC Bracket AS INTEGER
STATIC iDebugOld AS INTEGER

DIM Slope AS SINGLE
DIM SignDiff AS INTEGER

OldTCu = OKT(iMetal)
SignDiff = Sign(Diff)

'DIM adT$
'adT$ = STR$(KT(iBL)) + STR$(OC(iBLM(1)) * 1000000! / KT(iBL))
'CALL Pmsg(adT$)
'
'
' Check the current value of total metal and the BL-metal produced
' against the previous high and low values
IF (Iter = 1 OR Iter = 50 OR Iter = 100 OR Iter = 150 OR Iter = 200) THEN
   '
   ' This is the first iteration for this sample,
   ' reset the value of all variables
   CuLow = OKT(iMetal)
   CuHigh = OKT(iMetal)
   BLLow = Diff
   BLHigh = Diff
   Bracket = FALSE

   ' BEGIN DEBUG
   'iDebugOld = FREEFILE
   'OPEN "debug.prn" FOR OUTPUT AS #iDebugOld
   'PRINT #iDebugOld, " TCu      "; " Next TCu    "; " Diff     "; " CuHigh   "; " BLHigh "; " CuLow    "; " BLLow"
   ' BEGIN DEBUG
   'INPUT j
END IF

' Do we bracket yet?
IF ((NOT Bracket) AND (SignDiff <> Sign(BLLow))) THEN
   '
   ' Great!, we bracket now set the flag so won't keep checking
   Bracket = TRUE
END IF

IF Bracket THEN
'
' We bracket, so use false position
   IF (SignDiff < 0) THEN
      ' The last value replaces the low value
      BLLow = Diff
      CuLow = OKT(iMetal)
   ELSE
      ' The last value replaces the high value
      BLHigh = Diff
      CuHigh = OKT(iMetal)
   END IF
   OKT(iMetal) = CuLow + (CuHigh - CuLow) * BLLow / (BLLow - BLHigh)
ELSE
'
' We don't bracket, keep expanding the range until we do
   IF (SignDiff < 0) THEN
      ' The last value replaces the high value
      BLLow = BLHigh
      CuLow = CuHigh
      BLHigh = Diff
      CuHigh = OKT(iMetal)
   ELSE
      ' The last value replaces the low value
      BLHigh = BLLow
      CuHigh = CuLow
      BLLow = Diff
      CuLow = OKT(iMetal)
   END IF
   '
   ' Use the value of Diff to determine the direction
   ' that we have to move to get the next best guess
   ' for Total metal.  When Diff < 0, increase the value of KT
   ' and when Diff > 0, decrease the value of KT
   '
   SignDiff = Sign(Diff)
   OKT(iMetal) = OKT(iMetal) * (10 ^ -SignDiff)
END IF

IF (iDebug = 01) THEN CALL DebugLog ("A1Cu", "Iter:  "+STR$(Iter))
IF (iDebug = 01) THEN CALL DebugLog ("A1Cu", "CuLow:  "+STR$(CuLow))
IF (iDebug = 01) THEN CALL DebugLog ("A1Cu", "CuHigh: "+STR$(CuHigh))
IF (iDebug = 01) THEN CALL DebugLog ("A1Cu", "BLLow:  "+STR$(BLLow))
IF (iDebug = 01) THEN CALL DebugLog ("A1Cu", "BLHigh: "+STR$(BLHigh))
IF (iDebug = 01) THEN CALL DebugLog ("A1Cu", "NextMe: "+STR$(OKT(iMetal)))

' BEGIN DEBUG
'PRINT #iDebugOld, USING "##.###^^^^"; OldTCu; OKT(iMetal); Diff; CuHigh; BLHigh; CuLow; BLLow
' END DEBUG

END SUB

DEFDBL P-Z
FUNCTION FirstName$ (WName$)

n = INSTR(-1, WName$, ".")
IF n = 0 THEN
  FirstName$ = WName$
ELSE
  FirstName$ = LEFT$(WName$, n - 1)
END IF

END FUNCTION

FUNCTION ExtName$ (WName$)

n = INSTR(WName$, ".")
IF n = 0 THEN
  ExtName$ = ""
ELSE
  ExtName$ = RIGHT$(WName$, LEN(WName$) - n)
END IF

END FUNCTION
SUB GetArgs (Arg1$, Arg2$, Arg3$, LC50ALL, iDeriv, iExpert1, iExpert2, iExpert3, ErrorCode%, ErrorMsg$)
   '
   ' SUB GetArgs
   '
   ' Grabs and processes command line arguments to determine how the BLM will be executed
   '
   ' /Ax - select LC50 determination method (x can be 1,2,3, or 4)
   '       1 is just false position
   '       2 is a bounded Newton with false position used for failed Newton iterations
   '       3 is a cubic spline interpolation
   '       4 is a bounded polynominal interpolation
   '
   ' /Alx - Aluminum solubility
   '      0 - no solubility
   '      1 - precip only if over saturated
   '      2 - precip or dissolve all the time
   '
   ' /Bx - Calculate LA50 from LC50 data
   '     1 Calcuate a separate LA50 for each observation
   '     2 Calculate as log value
   '     3 For Al toxicity output values for both dissolved and precipitated mechanisms
   '
   ' /Cx - Convergence criteria adjustment, x can be 1, 2, or 3
   '     1 uses a default value
   '         Fract = ABS(BLMConc! - CritBL) / CritBL
   '         where BLMConc is the predicted value and CritBL is the desired value
   '         Default value of Fract is 0.001
   '     2 uses 10x the default value
   '     3 uses 100x the default value
   '     4 uses 1000x the default value
   '     5 uses 10000x the default value
   '
   ' /CAxx - Obtain a critical accumulation value from from row xx of a table of values from the
   '         parameter file. When this option is selected, the critical value and other useful info
   '         will be read from a table in the parameter file.  The table is delimited by the fields
   '         [CRITICAL START] and [CRITICAL END]. The columns in the table are comma delimited
   '         (csv format) and include:
   '            CA:              floating point numerical, the critical accumulation
   '            Species:         string, the species name (latin or common)
   '            Test Type:       string, typically will be acute or chronic
   '            Lifestage        string,
   '            Endpoint         string,
   '            Quantifier       string,
   '            References       string,
   '            Miscellaneous    string,
   '         Any cell with string info needs to be surrounded by quotes if it contains a comma.
   '
   '         As an example, here is a two row table entry:
   '
   '         [CRITICAL START]:
   '         CA,      Species,        Test Type, Lifestage, Endpoint, Quantifier, References,               Miscellaneous
   '         134.267, Mytilus edulis, Acute,     Embryo,    Survival, EC50,       "Martin et al. 1981; ...",
   '         103.88,  Mytilus edulis, Chronic,   Embryo,    Survival, EC20,       "Martin et al. 1981; ...",Estimated from EC50 SMAV by dividing by 1.15 before speciation
   '         [CRITICAL END]
   '
   '         The spacing to get columns lined up in the above example is not necessary, i.e., a
   '         typical csv format is all that is required.
   '
   ' /CO2(x+) - Sets the carbonate system to an open system. The pCO2 is set to 3.5 by default,
   '            or whatever is specified in (x+), which can go out several significant digits
   '            (looks for a space, the next switch, or the end of the string.
   '
   ' /D - Derivative mode
   '
   ' /E or /Ex - Expert mode with a couple of options (x can be 1, 2, or 3)
   '      1 (default) requires an extra column with LA50 values in the input file
   '      2 requires two additional columns, first with LA50 second with parameter file
   '      3 requires two additional columns for Al toxicity - first with DISS_A, second with PART_HI
   '
   ' /F - Write out separate information for humic and fulvic binding
   '
   ' /Fex - Iron solubility
   '      0 - no solubility
   '      1 - precip only if over saturated
   '      2 - preip or dissolve all the time
   '
   ' /Gxx - Debug mode - can turn all debug switches on or only for individual code sections
   '      00 - turn all debug statements off (default)
   '      01 - turn all debug statements on
   '      02 - only debug main
   '      03 - 0nly debug FirstGuess
   '      04 - only debug CalcSpeciesConc
   '      05 - only debug SUpdate
   '      06 - only debug CalcResidual
   '      07 - only debug CHESS
   '      08 - only debug output files
   '
   ' /H - display help and end
   '
   ' /I - Echo inputs - don't run the model, just read the input file and echo to output file
   '      Used for de-bugging
   '
   ' /Kx - Calculate a criteria report file including FMB (if minimum data are available)
   '       1 traditional, default format, both acute and chronic, used for Cu WQC
   '       2 acute only, could be any metal
   '       3 chronic only, could be any metal
   '
   ' /KMEANx - as /Kx above, but base FMB on means of log distributions (default if not specified,
   '           so same as /Kx)
   '
   ' /KMEDIANx - as /Kx above, but base the FMB on medians of the log distributions
   '
   ' /KSAlxxxxx - Aluminum solubility constant, where "xxxxx" is the log of the solubility constant.
   '
   ' /L - ignore all metals input data for toxic metal species, run in toxicity mode
   '
   ' /M - use a lower mass BL site density (useful in marine)
   '
   ' /Ox - Select output file format
   '     1 space delimited (traditional, default format)
   '     2 comma delimited (not implemented)
   '     3 MS Excel
   '
   ' /Pbx - Lead solubility
   '      0 - no solubility
   '      1 - precip only if over saturated
   '      2 - preip or dissolve all the time
   '
   ' /Q - Run in quiet mode (no iterations printed to screen)
   '
   ' /QQ - Run in very quiet mode (nothing printed to screen)
   '
   ' /S - use script file for command line input
   '
   ' /VER(x+) - Specify version text that will be part of the reported version in output files where
   '            (x+) is a string of just about any length.
   '
   ' /V - display version number and end
   '
   ' /W - Read a windows interface format file
   '
   ' /Zxxxxx - Use the value xxxxx as the Z-score for the target exceedance frequency for FMB
   '
   ' /ZDxxxxx - Use the value xxxxx as the 1 in xxxx day target exceedance frequency for FMB
   '            Default is 1095 (if no value is specified with the /Z or /ZD switch)




DIM iScript AS INTEGER
DIM iFile AS INTEGER
DIM sFile AS STRING
DIM TmpS AS STRING
DIM Ver_Prefix AS STRING

T$ = LTRIM$(RTRIM$(COMMAND$))
T$ = UCASE$(T$)

p1% = INSTR(T$, "/S")
IF (p1% <> 0) THEN
  '
  '  Command line arguments will come from a script file
  '
  MID$(T$, p1%, 2) = "  "
  iScript = TRUE
  sFile$ = LTRIM$(RTRIM$(T$))
  '
  ' Check to see if the file is there
  IF NOT Exist(sFile$) THEN
     ErrorMsg$= "File " & sFile$ & " not found."
     ErrorCode% = TRUE
     EXIT SUB
  END IF
  '
  ' Open the file
  iFile = FREEFILE
  OPEN sFile$ FOR INPUT AS iFile
  '
  ' Parameter file
  LINE INPUT #iFile, Arg1$
  Arg1$ = LTRIM$(RTRIM$(Arg1$))
  '
  ' Input File
  LINE INPUT #iFile, Arg2$
  Arg2$ = LTRIM$(RTRIM$(Arg2$))
  '
  '  Any arguments
  LINE INPUT #iFile, T$
  T$ = UCASE$(T$)
  Arg3$ = T$
  '
  '
  CLOSE iFile
ELSE
  iScript = FALSE
END IF

IF (LEN(T$) = 0) OR (INSTR(T$, "/H") <> 0) THEN
   CALL PrintHelp
   ErrorCode% = TRUE
   EXIT SUB
END IF

' /VERxxxxx
p1% = INSTR(T$, "/VER")
IF (p1% <> 0) THEN
  ' the value past the VER switch is the version info
  p2% =  p1% + 1
  DO UNTIL (MID$(T$,p2%,1)=" ") OR  (MID$(T$,p2%,1)="/") OR p2%=>LEN(T$)
    p2% =  p2% + 1
  LOOP
  'D_EF = VAL(LEFT$(MID$(T$,p1%+2,p2%-p1%-2),5))
  Ver_Prefix = MID$(T$,p1%+4,p2%-p1%-4)
  Ver_Prefix = "Ver "+Ver_Prefix +"."
  MID$(T$, p1%, p2%-p1%+1) = SPACE$(p2%-p1%+1)
  VersionID$ = Ver_Prefix+VersionID$
END IF

'CALL Pmsg(T$)
'
' Check to see if the users just wants version info
'
IF (INSTR(T$, "/V") <> 0) THEN
  '
  '  Report the version and exit
  '
  ErrorMsg$ = "BLM " + VersionID$
  ErrorCode% = TRUE
  EXIT SUB
END IF


'CALL Pmsg(T$)

p1% = INSTR(T$, "/KSAL")
iLogKsAl = FALSE
IF (p1% <> 0) THEN
  '
  '  Change the Log K for the Aluminum solubility constant
  '
  ' DIM LogKsAl AS GLOBAL DOUBLE
  ' DIM iLogKsAl AS GLOBAL INTEGER
  TmpS = RIGHT$(T$, LEN(T$)-p1%)
  p2% = INSTR(TmpS, ANY " /")
  IF p2% = 0 THEN p2% = LEN(TmpS)
  TmpS = MID$(TmpS, 5, p2%-4)
  iLogKsAl = TRUE
  LogKsAl = VAL(TmpS)
  IF (LogKsAl < 7) OR (LogKsAl > 11) THEN
      ErrorCode% = TRUE
      ErrorMsg$ = "Log KsAl should be between 7 and 11. Current value is "+ STR$( LogKsAl)
      EXIT SUB
  END IF
  'MID$(T$, p1%, 3) = "    "
END IF

iGlobalAl = 0
p1% = INSTR(T$, "/AL2")
IF (p1% <> 0) THEN
  '
  '  Aluminum will always be controled by solubility
  '
  iGlobalAl = 2
  MID$(T$, p1%, 4) = "    "
END IF

p1% = INSTR(T$, "/AL1")
IF (p1% <> 0) THEN
  '
  '  Aluminum will be controlled by solubility only when over saturated
  '
  iGlobalAl = 1
  MID$(T$, p1%, 4) = "    "
END IF
p1% = INSTR(T$, "/AL0")
IF (p1% <> 0) THEN
'
END IF
' /VERxxxxx
iGlobalFe = 0
p1% = INSTR(T$, "/FE2")
IF (p1% <> 0) THEN
  '
  '  Iron will always be controled by solubility
  '
  iGlobalFe = 2
  MID$(T$, p1%, 4) = "    "
END IF
p1% = INSTR(T$, "/FE1")
IF (p1% <> 0) THEN
  '
  '  Iron will be controled by solubility only when over saturated
  '
  iGlobalFe = 1
  MID$(T$, p1%, 4) = "    "
END IF
p1% = INSTR(T$, "/FE0")
IF (p1% <> 0) THEN
  '
  '  Iron will never be controled by solubility
  '
  iGlobalFe = 0
  MID$(T$, p1%, 4) = "    "
END IF

'------------- lead solubility ---------
iGlobalPb = 0
p1% = INSTR(T$, "/PB2")
IF (p1% <> 0) THEN
  '
  '  Lead will always be controled by solubility
  '
  iGlobalPb = 2
  MID$(T$, p1%, 3) = "    "
END IF
p1% = INSTR(T$, "/PB1")
IF (p1% <> 0) THEN
  '
  '  Lead will be controled by solubility only when over saturated
  '
  iGlobalPb = 1
  MID$(T$, p1%, 4) = "    "
END IF
p1% = INSTR(T$, "/PB0")
IF (p1% <> 0) THEN
  '
  '  Lead will never be controled by solubility
  '
  iGlobalPb = 0
  MID$(T$, p1%, 4) = "    "
END IF

'--------------end lead solublity ------


'
'  Determine if carbonate system is open or closed
'
p1% = INSTR(T$, "/CO2")
IF (p1% <> 0) THEN
   '
   '  Open system
   '  Inorganic carbon will always be controled by CO2 solubility
   '
   iGlobalCO2 = TRUE
   ' the value past the VER switch is the version info
   p2% =  p1% + 4
   TmpS = MID$(T$,p2%,1)
   IF (TmpS$ = " ") OR (TmpS$ = "/") OR (p2% >= LEN(T$)) THEN
      p2% = p2% - 1
      pPCO2 = -3.5
   ELSE
      DO UNTIL (MID$(T$,p2%,1)=" ") OR  (MID$(T$,p2%,1)="/") OR p2%>LEN(T$)
         p2% =  p2% + 1
      LOOP
      TmpS = MID$(T$, p1%+4, p2%-p1%-4)
      pPCO2 = -1 * VAL(TmpS)
   END IF
   MID$(T$, p1%, p2%-p1%+1) = SPACE$(p2%-p1%+1)

   IF (-pPCO2 < 0) THEN'this seems a little backwards just because we're asking for the (+) number then saving it as the (-) number
      ErrorCode% = TRUE
      ErrorMsg$ = "pCO2 must be greater than 0. Current value is "+ STR$(-pPCO2)
      EXIT SUB
   END IF

ELSE
  '
  '  Close system
  '  Inorganic carbon will always be controled by mass balance
  '
  iGlobalCO2 = FALSE
END IF

p1% = INSTR(T$, "/I")
IF (p1% <> 0) THEN
  '
  '  Don't run anything, just echo inputs to the output file
  '
  EchoInput = TRUE
  MID$(T$, p1%, 2) = "  "
ELSE
  EchoInput = FALSE
END IF

p1% = INSTR(T$, "/L")
IF (p1% <> 0) THEN
  '
  '  Run LC50 determination for all observations regardless of
  '  whether there is dissolved metal or -999 input to model
  '
  LC50ALL = TRUE
  MID$(T$, p1%, 2) = "  "
ELSE
  LC50ALL = FALSE
END IF
'
' Expert mode checks (/E2 and /E3 need to be done before /E)
'
'  /E  : Additional column in input file will be used to specify critical accumulation (LA50)
'  /E1 : same as /E
'  /E2 : two additional columns, first with LA50 second with parameter file
'  /E3 : two additional columns for Al toxicity - first with DISS_A second with PART_HI
p1% = INSTR(T$, "/E2")
IF (p1% <> 0) THEN
  '
  '  Second type of "expert mode" with parameter file in extra
  '  column
  '
  iExpert2 = TRUE
  MID$(T$, p1%, 3) = "   "
ELSE
  iExpert2 = FALSE
END IF

p1% = INSTR(T$, "/E3")
IF (p1% <> 0) THEN
  '
  '  Third type of "expert mode" with Al toxicity parameters in
  '  two extra columns
  '
  iExpert3 = TRUE
  MID$(T$, p1%, 3) = "   "
ELSE
  iExpert3 = FALSE
END IF

p1% = INSTR(T$, "/E")
p2% = INSTR(T$, "/E1")
IF (p1% <> 0) OR (p2% <> 0) THEN
  '
  '  Assume "expert mode" with critical senstivity values in extra
  '  column
  '
  iExpert1 = TRUE
  IF (p1% <> 0) THEN
    MID$(T$, p1%, 2) = "  "
  ELSEIF (p2% <> 0) THEN
    MID$(T$, p2%, 3) = "   "
  END IF
ELSE
  iExpert1 = FALSE
END IF

p1% = INSTR(T$, "/W")
IF (p1% <> 0) THEN
  '  The /W switch specifies that the windows format will be used.
  '  Note, the windows format has additional info. used for monte carlo.
  '
  iOrig = FALSE
  MID$(T$, p1%, 2) = "  "
ELSE
  '
  '  Assume that the input file has the original (i.e., pre windows version)
  '  format.
  iOrig = TRUE
END IF

p1% = INSTR(T$, "/Q")
IF (p1% <> 0) THEN
  '  The /Q switch means run in quiet mode (Report no. obs only on screen)
  '  The /QQ switch means run in very quiet mode (Report nothing to screen)
  '
  iQuiet = TRUE
  IF (MID$(T$, p1% + 2, 1) = "Q") THEN
     iVQuiet = TRUE
     MID$(T$, p1%, 3) = "   "
  ELSE
     iVQuiet = FALSE
     MID$(T$, p1%, 2) = "  "
  END IF
ELSE
  '
  '  Verbose ouput (report everything to the screen)
  iQuiet = FALSE
  iVQuiet = FALSE
END IF

'
' Check to see if a smaller mass balance for BL sites is needed
'
iLowMass = FALSE
p1% = INSTR(T$, "/M")
IF (p1% <> 0) THEN
  '
  '  Use Low Mass BL site density
  '
  MID$(T$, p1%, 2) = "  "
  iLowMass = TRUE
END IF

p1% = INSTR(T$, "/A1")
IF (p1% <> 0) THEN
  '
  '  Adjust method 1 is just false position
  '
  iAdjust = ONE
  MID$(T$, p1%, 3) = "   "
END IF

p1% = INSTR(T$, "/A2")
IF (p1% <> 0) THEN
  '
  '  Adjust method 2 is a bounded Newton with
  '  false position used for failed Newton iterations
  '
  iAdjust = TWO
  'MID$(T$, p1%, 2) = "  "
  MID$(T$, p1%, 3) = "   "
END IF

p1% = INSTR(T$, "/A3")
IF (p1% <> 0) THEN
  '
  '  Adjust method 3 is a cubic spline interpolation
  '
  iAdjust = THREE
  MID$(T$, p1%, 3) = "   "
END IF

p1% = INSTR(T$, "/A4")
IF (p1% <> 0) THEN
  '
  '  Adjust method 4 is a bounded polynominal interpolation
  '
  iAdjust = FOUR
  MID$(T$, p1%, 3) = "   "
END IF

'CALL Pmsg(T$)
iCalcLA50 = FALSE
iCalcLA50Log = FALSE
p1% = INSTR(T$, "/B1")
IF (p1% <> 0) THEN
  '
  '  Calculate LA50 values from input LC50
  '
  iCalcLA50 = TRUE
  MID$(T$, p1%, 3) = "   "
END IF
p1% = INSTR(T$, "/B2")
IF (p1% <> 0) THEN
  '
  '  Calculate LA50 values from input LC50
  '
  iCalcLA50 = TRUE
  iCalcLA50Log = TRUE
  MID$(T$, p1%, 3) = "   "
END IF
p1% = INSTR(T$, "/B3")
IF (p1% <> 0) THEN
  '
  '  Calculate LA50 values from input LC50
  '
  iCalcLA50 = TRUE
  iCalcLA50Al = TRUE
  MID$(T$, p1%, 3) = "   "
END IF

'
' Critical Accumulation Table (CAT) /CAxx
'   (note, this must be checked before the /Cx switch)
'
p1% = INSTR(T$, "/CA")
'p2% = len(T$)
'TmpS = MID$(T$,p1%,p2%-p1%)
IF (p1% <> 0) THEN
  iHaveCAT = TRUE
  ' the value past the /CA switch is the CAT row that will be used for this simulation
  p1% = p1% + 3
  p2% =  p1% + 1
  DO UNTIL (MID$(T$,p2%,1)=" ") OR  (MID$(T$,p2%,1)="/") OR p2%>LEN(T$)
    p2% =  p2% + 1
  LOOP
  iCAT_Row = VAL(MID$(T$,p1%,p2%-p1%))
  PRINT "iCAT_Row = " & STR$(iCAT_Row)
  IF (iCAT_Row <=0) OR (iCAT_Row >=100) THEN
      PRINT "Reverting to first row of CA Table."
      iCAT_Row = 1
  END IF
  MID$(T$, p1%-3, p2%-p1%+3) = "   " & REPEAT$(p2%-p1%, " ")
END IF


'
' Convergence criteria switch /C
'
p1% = INSTR(T$, "/C1")
IF (p1% <> 0) THEN
  '
  '  Default convergence criterion
  '
  iConverge = ONE
  MID$(T$, p1%, 3) = "   "
END IF
'
p1% = INSTR(T$, "/C2")
IF (p1% <> 0) THEN
  '
  '  Default convergence criterion * 10
  '
  iConverge = TWO
  MID$(T$, p1%, 3) = "   "
END IF
'
p1% = INSTR(T$, "/C3")
IF (p1% <> 0) THEN
  '
  '  Default convergence criterion * 100
  '
  iConverge = THREE
  MID$(T$, p1%, 3) = "   "
END IF
p1% = INSTR(T$, "/C4")
IF (p1% <> 0) THEN
  '
  '  Default convergence criterion * 1000
  '
  iConverge = FOUR
  MID$(T$, p1%, 3) = "   "
END IF
p1% = INSTR(T$, "/C5")
IF (p1% <> 0) THEN
  '
  '  Default convergence criterion * 10000
  '
  iConverge = FIVE
  MID$(T$, p1%, 3) = "   "
END IF

'
'
' Output file format switch /O
'
iOFormat = ONE   ' set default value
p1% = INSTR(T$, "/O1")
IF (p1% <> 0) THEN
  '
  '  Space delimited (Default) output file
  '
  iOFormat = ONE
  MID$(T$, p1%, 3) = "   "
END IF
'
p1% = INSTR(T$, "/O2")
IF (p1% <> 0) THEN
  '
  '  Comma delimited output file
  '
  iOFormat = TWO
  MID$(T$, p1%, 3) = "   "
END IF
'
p1% = INSTR(T$, "/O3")
IF (p1% <> 0) THEN
  '
  '  MS Excel output file
  '
  iOFormat = THREE
  MID$(T$, p1%, 3) = "   "
END IF
'
p1% = INSTR(T$, "/O4")
IF (p1% <> 0) THEN
  '
  '  MS Excel output file
  '
  iOFormat = FOUR
  MID$(T$, p1%, 3) = "   "
END IF
'
' Comma delimited files are not yet implemented
' so default to space delimited
'IF (iOFormat = TWO) THEN iOFormat = ONE

'
' FMB calculation based on median values
'
FMB_option = FMBMean    ' default to calculation based on mean values
p1% = INSTR(T$, "/KMEDIAN")
IF (p1% <> 0) THEN
  '
  '  Calculate FMB based on median values
  '
  FMB_option = FMBMedian
  IF (p1% <> 0) THEN  MID$(T$, p1%, 8) = "      /K"
END IF
'
' FMB calculation based on mean values
'
p1% = INSTR(T$, "/KMEAN")
IF (p1% <> 0) THEN
  '
  '  Calculate FMB based on median values
  '
  FMB_option = FMBMean
  IF (p1% <> 0) THEN  MID$(T$, p1%, 8) = "    /K"
END IF

'
'
' Criteria calculation and output file creation /K1
'
iCriteria = ZERO   ' set default value
p1% = INSTR(T$, "/K1")
p2% = INSTR(T$, "CuOH5%le")
p3% = INSTR(Arg1$, "CuOH5%le")
IF (p1% <> 0) OR (p2% <> 0) OR (p3% <> 0) THEN
  '
  '  Calculate criteria and FMB, and create criteria report file
  '
  iCriteria = ONE
  IF (p1% <> 0) THEN  MID$(T$, p1%, 3) = "   "
  '
  '  A criteria only makes sense in toxicity mode so set that flag as well
  LC50ALL = TRUE
END IF
'
'
' Criteria calculation and output file creation /K2
'
p1% = INSTR(T$, "/K2")
IF (p1% <> 0) THEN
  '
  '  Calculate Acute only criteria and FMB, and create criteria report file
  '
  iCriteria = TWO
  IF (p1% <> 0) THEN  MID$(T$, p1%, 3) = "   "
  '
  '  A criteria only makes sense in toxicity mode so set that flag as well
  LC50ALL = TRUE
END IF
'
' Criteria calculation and output file creation /K3
'
p1% = INSTR(T$, "/K3")
IF (p1% <> 0) THEN
  '
  '  Calculate Acute only criteria and FMB, and create criteria report file
  '
  iCriteria = THREE
  IF (p1% <> 0) THEN  MID$(T$, p1%, 3) = "   "
  '
  '  A criteria only makes sense in toxicity mode so set that flag as well
  LC50ALL = TRUE
END IF

' /ZDxxxxx
' exceedance frequency specified as once in xxxxx days
p1% = INSTR(T$, "/ZD")
IF (p1% <> 0) THEN
  iCalc_FMB = TRUE
  ' the value past the Z switch is the target exceedance frequency in days
  p2% =  p1% + 1
  DO UNTIL (MID$(T$,p2%,1)=" ") OR  (MID$(T$,p2%,1)="/") OR p2%=>LEN(T$)
    p2% =  p2% + 1
  LOOP

  D_EF = VAL(MID$(T$,p1%+3,p2%-p1%-2))
  Z_EF = ZScoreNorm(D_EF-1, D_EF)

  PRINT "D_EF = " & STR$(D_EF)
  PRINT "Z_EF = " & STR$(Z_EF)
  IF Z_EF <=0 THEN
      PRINT "Reverting to default Z_EF = " & STR$(Z_EF_DEF)
      Z_EF = Z_EF_DEF
  END IF

  MID$(T$, p1%, p2%-p1%+1) = SPACE$(p2%-p1%+1)

END IF

' /Zxxxxx
' Exceedance frequency specified as a z-score
p1% = INSTR(T$, "/Z")
IF (p1% <> 0) THEN
  iCalc_FMB = TRUE
  ' the value past the Z switch is the target z score
  p2% =  p1% + 1
  DO UNTIL (MID$(T$,p2%,1)=" ") OR  (MID$(T$,p2%,1)="/") OR p2%=>LEN(T$)
    p2% =  p2% + 1
  LOOP
  Z_EF = VAL(LEFT$(MID$(T$,p1%+2,p2%-p1%-2),5))
  PRINT "Z_EF = " & STR$(Z_EF)
  IF Z_EF <=0 THEN
      PRINT "Reverting to default Z_EF = " & STR$(Z_EF_DEF)
      Z_EF = Z_EF_DEF
  END IF
END IF

' If neither a z-score or days for exceedance frequency was specified
' don't do the FMB calculation (previously: use the default z-score)
IF (iCriteria > 0) AND (Z_EF <= 0) THEN
   PRINT "Not performing FMB calculation." '"Using default Z_EF = " & STR$(Z_EF_DEF)
   Z_EF = Z_EF_DEF
END IF

' if days for exceedance frequency was not already specified, calculate it
IF Z_EF > 0 THEN
    IF D_EF = 0 THEN
  '
        ' Calculate D_EF from Z_EF
  '
        DIM tmp AS EXTENDED
        tmp = NormalCDF((Z_EF), 0, 1)
        D_EF = 1 / (1 - tmp)
    END IF
END IF

   ' /Gxx - Debug mode - can turn all debug switches on or only for individual code sections
   '      00 - turn all debug statements off (default)
   '      01 - turn all debug statements on
   '      02 - only debug main
   '      03 - 0nly debug FirstGuess
   '      04 - only debug CalcSpeciesConc
   '      05 - only debug SUpdate
   '      06 - only debug CalcResidual

p1% = INSTR(T$, "/G")
IF (p1% <> 0) THEN
  '
  '  Write out debug information
  '  The specific debug info will depend on the values after "G"
  '
  IF INSTR(p1%, T$, "/G00") THEN        '      00 - turn all debug statements off (default)
      iDebug = 00
      MID$(T$, p1%, 4) = "    "
  ELSEIF INSTR(p1%, T$, "/G01") THEN    '      01 - turn all debug statements on
      iDebug = 01
      MID$(T$, p1%, 4) = "    "
  ELSEIF INSTR(p1%, T$, "/G02") THEN    '      02 - only debug main
      iDebug = 02
      MID$(T$, p1%, 4) = "    "
  ELSEIF INSTR(p1%, T$, "/G03") THEN    '      03 - 0nly debug FirstGuess
      iDebug = 03
      MID$(T$, p1%, 4) = "    "
  ELSEIF INSTR(p1%, T$, "/G04") THEN    '      04 - only debug CalcSpeciesConc
      iDebug = 04
      MID$(T$, p1%, 4) = "    "
  ELSEIF INSTR(p1%, T$, "/G05") THEN    '      05 - only debug SUpdate
      iDebug = 05
      MID$(T$, p1%, 4) = "    "
  ELSEIF INSTR(p1%, T$, "/G06") THEN    '      06 - only debug CalcResidual
      iDebug = 06
      MID$(T$, p1%, 4) = "    "
  ELSEIF INSTR(p1%, T$, "/G07") THEN    '      07 - only debug CHESS main sub
      iDebug = 07
      MID$(T$, p1%, 4) = "    "
  ELSEIF INSTR(p1%, T$, "/G08") THEN    '      08 - only debug output files
      iDebug = 08
      MID$(T$, p1%, 4) = "    "
  ELSEIF INSTR(p1%, T$, "/G09") THEN    '      09 - only debug Al dose response info
      iDebug = 09
      MID$(T$, p1%, 4) = "    "
  ELSE                                  '      01 - turn all debug statements on (alternate)
      '
      ' For this condition, we know that the switch /G was used with no number,
      ' and this is the equivalent of turning on all debug code
      iDebug = 01
      MID$(T$, p1%, 2) = "  "
  END IF

ELSE                                    '      00 - turn all debug statements off (default)
  ' Debugging should stay off
  iDebug = 00
END IF

p1% = INSTR(T$, "/F")
IF NOT iScript THEN
  p1% = INSTR(T$, " ")
  Arg1$ = LTRIM$(RTRIM$(LEFT$(T$, p1% - 1)))
  Arg2$ = LTRIM$(RTRIM$(RIGHT$(T$, LEN(T$) - p1%)))
END IF

'
' Check to see if the file is there
IF NOT Exist(Arg1$) THEN
  ErrorMsg$= "File " & Arg1$ & " not found.  Use /H for Help."
  ErrorCode% = TRUE
  EXIT SUB
END IF
'
' Check to see if the file is there
IF NOT Exist(Arg2$) THEN
  ErrorMsg$= "File " & Arg2$ & " not found.  Use /H for Help."
  ErrorCode% = TRUE
  EXIT SUB
END IF

EXIT SUB

  ErrorTrap:
    SELECT CASE ERRCLEAR
      CASE 52   : PRINT "Bad file name or number."
      CASE 53   : PRINT "File not found."
      CASE 70   : PRINT "Permission denied (network access error)."
      CASE 71   : PRINT "Drive not ready."
      CASE 76   : PRINT "That path doesn't exist."
      CASE ELSE : PRINT "Unknown error!"
    END SELECT
END SUB

SUB GetData (iInFile, LC50, iExpert1, CriticalValue)
'
'
' Read in a line of input data and determine if we are
' predicting LC50 (copper = -999)
'
' Local variables
'DIM PercHA AS SINGLE            ' Percent of DOC that is humic acid
DIM i AS INTEGER                ' Loop counter
DIM Tmp AS STRING                 ' used for garbage string collection

'Temp %HumicA pH Cu DOC Ca Mg Na K SO4 Cl CO3 S
IF NOT iOrig THEN
   INPUT #iInFile, Tmp$
   IF (Tmp$ = "T") THEN
      iComplete = TRUE
   ELSE
      iComplete = FALSE
   END IF

   INPUT #iInFile, Tmp$
   SiteL$ = Tmp$
END IF

INPUT #iInFile, Tmp$
DataL$ = Tmp$

INPUT #iInFile, SysTemp

SysTemp = SysTemp + 273 'CHECK C to K conversion

INPUT #iInFile, PercHA!
PercHA! = PercHA! / 100!

FOR i = 1 TO iBL - 1
   INPUT #iInFile, OKT(i)
NEXT i

IF (iExpert1) THEN
   ' Read in critical BL number
   INPUT #iInFile, CriticalValue
   '
   '  Convert units from nmol/g wet
   IF (CriticalValue > 0!) THEN
      CriticalValue = CriticalValue * OKT(iBL) / 1000000!
   END IF
END IF

' Convert pH to H
OKT(1) = 10 ^ -OKT(1)
OC(1) = OKT(1)

' Check copper to see if we need to predict LC50
IF (OKT(iMetal) <= -990) THEN
   '
   ' We are predicting LC50, set flags accordingly
   LC50 = TRUE
   '
   ' Use something small for initial Cu conc
   OKT(iMetal) = .00000001#
END IF
'LC50 = FALSE
'LC50 = TRUE

' Convert Cu from ug/L to mole/L
'OKT(iMetal) = OKT(iMetal) * .000001# / 64.5 ' CHECK

IF (iWHAM) THEN
   ' Convert DOC from mg C/L to g OM/liter
   OKT(iDOC) = OKT(iDOC) * .001 * 2
   OKT(iHA1) = .0001
   OKT(iHA2) = .0001
   OKT(iHA3) = .0001
   OKT(iHA4) = .0001
   OKT(iHB1) = .0001
   OKT(iHB2) = .0001
   OKT(iHB3) = .0001
   OKT(iHB4) = .0001
   OKT(iFA1) = .0001
   OKT(iFA2) = .0001
   OKT(iFA3) = .0001
   OKT(iFA4) = .0001
   OKT(iFB1) = .0001
   OKT(iFB2) = .0001
   OKT(iFB3) = .0001
   OKT(iFB4) = .0001
END IF

END SUB

SUB GetData2 (iInFile, LC50, iExpert1, iExpert2, iExpert3, PFile$, CriticalValue, ErrorCode%, ErrorMsg$)
'
'
' Read in a line of input data and determine if we are
' predicting LC50 (copper = -999)
'
' Local variables
'DIM PercHA AS SINGLE           ' Percent of DOC that is humic acid
DIM i AS INTEGER                ' Loop counter
DIM j AS INTEGER
DIM Tmp AS STRING               ' Temporary string
DIM NFields AS INTEGER          ' Number of data fields on input line
DIM MaxFields AS INTEGER        ' Maximum number of data fields expected
DIM DField$(1 TO 35)            ' input line parsed into delimited fields
DIM Delim AS STRING             ' string that holds delimiter characters
DIM TmpD AS DOUBLE
Delim$ = " ,"

LINE INPUT #iInFile, Tmp$
NFields = Delimit%(Tmp$, Delim$)
MaxFields = 57
IF iExpert1 THEN MaxFields = MaxFields + 1
IF iExpert2 THEN MaxFields = MaxFields + 1
IF iExpert3 THEN MaxFields = MaxFields + 2
IF NFields > MaxFields THEN
   '
   '  Too many input fields, must be some problem
   ErrorCode% = TRUE
   ErrorMsg$ = "Too many data fields in input file."
   EXIT SUB
END IF

CALL NParse(Tmp$, Delim$, DField$(), NFields)
IF NFields < 14 THEN
   '
   '  Too few input fields, must be some problem
   ErrorCode% = TRUE
   ErrorMsg$ = "Too few data fields in input file (" & STR$(NFields) & ")."
   EXIT SUB
END IF

' debug code
'PRINT Tmp$
'FOR i=1 TO NFields
'    PRINT DField$(i)
'NEXT i
'   ErrorCode% = TRUE
'   ErrorMsg$ = "end of test."
'   EXIT SUB

'Temp %HumicA pH Cu DOC Ca Mg Na K SO4 Cl CO3 S
j  = 1
IF iOrig THEN
   '
   ' This is the original file format (i.e., pre-windows)
   ' we must assume that non-windows input files are complete
   iComplete = TRUE
ELSE
   ' These fields only exist on windows version input files
   '
   'Pmsg(DField$(j))
   IF (INSTR(DField$(j),"T")) THEN
      iComplete = TRUE
   ELSE
      iComplete = FALSE
   END IF
   'Pmsg(STR$(iComplete))
   j = j + 1

   SiteL$ = DField$(j)
   j = j + 1
END IF
'CALL PMsg(STR$(iComplete))

DataL$ = DField$(j)
j = j + 1

SysTemp = VAL(DField$(j)) + 273
j = j + 1

PercHA! = VAL(DField$(j)) / 100!
j = j + 1

FOR i = 1 TO iBL - 1
   TmpD = VAL(DField$(j))
   OKT(i) = TmpD * MassVal(SType(i) + 1)
   j = j + 1
NEXT i

IF (iExpert1) THEN
   ' Check to see if there is an LA50 column
   IF (j > NFields) THEN
      ErrorCode% = TRUE
      ErrorMsg$ = "Must include a valid LA50 column in the input file to run in expert mode."
      EXIT SUB
   END IF
'   IF (CriticalValue < 0) THEN
'      '  There is a missing value flag in LA50 column
'      ErrorCode% = TRUE
'      ErrorMsg$ = "Must input valid LA50 to run in expert mode."
'      EXIT SUB
'   END IF
   ' Read in critical BL number
   CriticalValue = VAL(DField$(j))
   j = j + 1
   IF (CriticalValue < 0) AND iComplete THEN
      '  There is a missing value flag in LA50 column
      ErrorCode% = TRUE
      ErrorMsg$ = "Must input valid LA50 to run in expert mode."
      EXIT SUB
   END IF
   '
   '  Convert units from nmol/g wet
'
'  We don't need this, since it is being done within the main loop
'
'   IF (CriticalValue > 0!) THEN
'      IF (iBLM(1) > 10) THEN
'         CriticalValue = CriticalValue * OKT(iBL) / 1000000!
'      END IF
'   END IF
END IF

IF (iExpert2) THEN
   ' Check to see if there is a column with a parameter file name
   IF (j > NFields) THEN
      ErrorCode% = TRUE
      ErrorMsg$ = "You must include a valid parameter file in the input file to run in expert2 mode."
      EXIT SUB
   END IF
   ' Read in parameter file
   PFile$ = DField$(j)
   j = j + 1
   IF (LEN(PFile$) = 0) THEN
      '  There is a missing parameter file name
      ErrorCode% = TRUE
      ErrorMsg$ = "You must include a valid parameter file in the input file to run in expert2 mode."
      EXIT SUB
   END IF
END IF

IF (iExpert3) THEN
   ' Check to see if there is are two additional columns
   IF (j+1 > NFields) THEN
      ErrorCode% = TRUE
      ErrorMsg$ = "Must include columns for DISS_A and PART_HI in the input file to run in expert mode 3."
      EXIT SUB
   END IF

   DISS_A_expert = VAL(DField$(j))
   j = j + 1

   PART_HI_expert = VAL(DField$(j))
   j = j + 1
END IF

' Convert pH to H
OKT(1) = 10 ^ -(OKT(1) / MassVal(SType(1) + 1))
OC(1) = OKT(1)

' Check copper to see if we need to predict LC50
LC50 = FALSE
IF (OKT(iMetal) <= -990) THEN
   '
   ' We are predicting LC50, set flags accordingly
   LC50 = TRUE
   '
   ' Use something small for initial Cu conc
   OKT(iMetal) = .00000001#
END IF
'LC50 = FALSE
'LC50 = TRUE

' Convert Cu from ug/L to mole/L
'OKT(iMetal) = OKT(iMetal) * .000001# / 64.5 ' CHECK

IF (iWHAM) THEN
   ' Convert DOC from mg C/L to g OM/liter
   OKT(iDOC) = OKT(iDOC) * .001 * 2
   OKT(iHA1) = .0001
   OKT(iHA2) = .0001
   OKT(iHA3) = .0001
   OKT(iHA4) = .0001
   OKT(iHB1) = .0001
   OKT(iHB2) = .0001
   OKT(iHB3) = .0001
   OKT(iHB4) = .0001
   OKT(iFA1) = .0001
   OKT(iFA2) = .0001
   OKT(iFA3) = .0001
   OKT(iFA4) = .0001
   OKT(iFB1) = .0001
   OKT(iFB2) = .0001
   OKT(iFB3) = .0001
   OKT(iFB4) = .0001
END IF

END SUB

SUB InitData (iInFile, sInFile$, NData&, ErrorCode%, ErrorMsg$)
'
'
' Initialize input data file
'   - open the file with the correct file handle
'   - read in the number of rows of data
'   - get rid of header information
'
' Local variables
DIM Tmp AS STRING
DIM i AS INTEGER
DIM n AS INTEGER
'T$ = "Top of InitData"
'CALL Pmsg(T$)
'
' Open the BLM input file
'OPEN sInFile$ FOR INPUT AS iInFile
TRY
   OPEN sInFile$ FOR INPUT AS iInFile
CATCH
   ErrorCode% = TRUE
   ErrorMsg$ = "Cannot open or access " + sInFile$ + "."
   EXIT SUB
END TRY
'
' RCS:  Note, we really should read in the version number
' at this point, instead of skipping it
'
' Skip the monte carlo info in the input file
' unless the "original file format" flag is set
IF NOT iOrig THEN
'
' Edit for allowing Rooni's extra lines, original line
' follows (as comment), modified line after that
'  FOR i = 1 TO 16
   FOR i = 1 TO 19
      LINE INPUT #iInFile, Tmp$
   NEXT i
END IF
'
' Get the number of data points
INPUT #iInFile, NData&
'
' Strip out the comments and column headings
IF iOrig THEN
   n = 2      ' column labels and unit labels in original file format
ELSE
   n = 3      ' plus another line to describe the dataset
END IF
FOR i = 1 TO n
   LINE INPUT #iInFile, Tmp$
NEXT i


IF (MajorVer >=3) AND (MinorVer >= 10) THEN
    CALL SetCtoM310
ELSE
    AqueousCtoM! = 1!
    SurfaceCtoM! = 1!
    CALL SetCtoM(AqueousCtoM!, SurfaceCtoM!)
END IF

'T$ = "End of InitData"
'CALL Pmsg(T$)
'T$ = FORMAT$(NData)
'CALL Pmsg(T$)
END SUB

SUB MessageFile(iErrCode&, iObs%, MsgString$)
    DIM iFile AS INTEGER
    DIM T AS STRING

    iFile = FREEFILE
    OPEN "blm.message" FOR OUTPUT AS iFile
    PRINT #iFile, iErrCode&; Comma$;iObs%; Comma$; Quote$+MsgString$+Quote$

END SUB

SUB NParse(Work$, Delim$, NArray$(), k%)
   '
   ' Parse delimited fields in Work$ into individual array elements in
   ' NArray$() using delimiter characters specified in Delim$.  Ignore
   ' consecutive delimiters and keep fields within quotes intact.
   '
   DIM Last AS INTEGER
   DIM BeginPtr AS INTEGER
   DIM BeginQuote AS INTEGER
   DIM j AS INTEGER
   FALSE = 0
   TRUE = NOT FALSE

   Last = FALSE
   BeginPtr = 0
   j = 1
   FOR X% = 1 TO LEN(Work$)
       T$ = MID$(Work$, X%, 1)
       IF(INSTR(T$, CHR$(34))) THEN
          ' Quotes are special, ignore anything between quotes
          IF BeginQuote = FALSE THEN
             '  This is the first quote of a quoted fieled
             BeginQuote = TRUE
             BeginPtr = X%
          ELSE
             '  This is the end quote of a quoted fieled
             BeginQuote = FALSE
             NArray$(j) = MID$(Work$,BeginPtr,X%-BeginPtr+1)
             j = j + 1
             Last = TRUE
          END IF
       ELSE
          IF BeginQuote = FALSE THEN
             IF(INSTR(MID$(Work$, X%, 1), ANY Delim$)) THEN
                ' This character is a delimiter.
                ' Was the last character a delimiter as well?
                IF Last = TRUE THEN
                   ' last character was delimiter, skip this one
                ELSE
                   ' This character is a delimiter and last was not
                   ' so a data field must have just ended
                   IF BeginPtr > 0 THEN
                      ' unless BeginPtr = 0 then these delimiters are
                      ' at the begining of the string
                      NArray$(j) = MID$(Work$,BeginPtr,X%-BeginPtr)
                      j = j + 1
                   END IF
                   Last = TRUE
                END IF
                IF (X% = LEN(Work$)) THEN
                   '
                   ' We are at the end of the string, so it dosn't
                   ' matter if the last character(s) is(are)
                   ' delimiters.
                   j = j - 1
                END IF
             ELSE
                ' This character is not a delimiter
                IF (Last = TRUE) OR (X% = 1) THEN
                   ' The last character was a delimiter and this
                   ' one in not, so must be the start of a field
                   ' Or, this is the first character of the string
                   ' and not a delimiter so also the start of a field
                   BeginPtr = X%
                END IF
                Last = FALSE
                IF X% = LEN(Work$) THEN
                   ' this is the last character, if the whole string
                   ' doesn't end in a delimiter we still want to save
                   ' the last field
                   NArray$(j) = MID$(Work$,BeginPtr,X%-BeginPtr+1)
                END IF
             END IF
          END IF
       END IF
   NEXT X%
   k% = j
END SUB

FUNCTION OpenEXCEL(mFileName$, xlsFileNumber&) AS LONG
    'enable the buffer. This is optional, but enabling a buffer will speed
    'up file access.
    'stat& = xlsBuffer(%XLSTRUE, (512 * 1024))  'a 512K buffer

    'Create the new spreadsheet
    'mFileName$ = ".\vbtest.xls"  'create spreadsheet in the current directory
    stat& = xlsCreateFile(mFileName$, xlsFileNumber&)
    'stat& returns non-zero if an error occured.

    'set a Password for the file. If set, the rest of the spreadsheet will
    'be encrypted. If a password is used it must immediately follow the
    'xlsCreateFile function call.

    'This is different then protecting the spreadsheet (see below).
    'NOTE: For some reason this function does not work. Excel will
    'recognize that the file is password protected, but entering the password
    'will not work. Also, the file is not encrypted. Therefore, do not use
    'this function until I can figure out why it doesn't work. There is not
    'much documentation on this function available.
    'stat& = xlsSetFilePassword("PAUL_SQUIRES")


    'specify whether to print the gridlines or not
    'this should come before the setting of fonts and margins
    stat& = xlsPrintGridLines(%XLSTRUE, xlsFileNumber&)


    'it is a good idea to set margins, fonts and column widths
    'prior to writing any text/numerics to the spreadsheet. These
    'should come before setting the fonts.

    stat& = xlsSetMargin(%xlsTopMargin, 1.5, xlsFileNumber&)   'set to 1.5 inches
    stat& = xlsSetMargin(%xlsLeftMargin, 1.5, xlsFileNumber&)
    stat& = xlsSetMargin(%xlsRightMargin, 1.5, xlsFileNumber&)
    stat& = xlsSetMargin(%xlsBottomMargin, 1.5, xlsFileNumber&)


    'Up to 4 fonts can be specified for the spreadsheet. This is a
    'limitation of the Excel 2.1 format. For each value written to the
    'spreadsheet you can specify which font to use.

    stat& = xlsSetFont("Arial", 10, %xlsNoFormat, xlsFileNumber&)             'font0
    stat& = xlsSetFont("Arial", 10, %xlsBold, xlsFileNumber&)                 'font1
    stat& = xlsSetFont("Arial", 10, %xlsBold + %xlsUnderline, xlsFileNumber&) 'font2
    stat& = xlsSetFont("Courier", 18, %xlsItalic, xlsFileNumber&)             'font3

    'Column widths are specified in Excel as 1/256th of a character.
    stat& = xlsSetColumnWidth(1, 2, 25, xlsFileNumber&)
    stat& = xlsSetColumnWidth(3, 3, 12, xlsFileNumber&)

    'set the global row height for the entire spreadsheet
    stat& = xlsSetDefaultRowHeight(16, xlsFileNumber&)

    'set the height of the first two rows a little bigger to allow for the
    'title of the spreadsheet.
    'stat& = xlsSetRowHeight(1, 24)
    'stat& = xlsSetRowHeight(2, 24)

    'set any header or footer that you want to print on
    'every page. This text will be centered at the top and/or
    'bottom of each page. The font will always be the font that
    'is specified as font0, therefore you should only set the
    'header/footer after specifying the fonts through SetFont.
    T$ = "BLM version " + VersionID$
    stat& = xlsSetHeader(T$, xlsFileNumber&)
    T$ = "BLM Output File, " + mFileName$ + ", created " + DATE$
    stat& = xlsSetFooter(T$, xlsFileNumber&)
END FUNCTION


FUNCTION OpenEXCELWQC(mFileName$, xlsFileNumber&) AS LONG
    'enable the buffer. This is optional, but enabling a buffer will speed
    'up file access.
    'stat& = xlsBuffer(%XLSTRUE, (512 * 1024))  'a 512K buffer

    'Create the new spreadsheet
    'mFileName$ = ".\vbtest.xls"  'create spreadsheet in the current directory
    stat& = xlsCreateFile(mFileName$, xlsFileNumber&)
    'stat& returns non-zero if an error occured.

    'set a Password for the file. If set, the rest of the spreadsheet will
    'be encrypted. If a password is used it must immediately follow the
    'xlsCreateFile function call.

    'This is different then protecting the spreadsheet (see below).
    'NOTE: For some reason this function does not work. Excel will
    'recognize that the file is password protected, but entering the password
    'will not work. Also, the file is not encrypted. Therefore, do not use
    'this function until I can figure out why it doesn't work. There is not
    'much documentation on this function available.
    'stat& = xlsSetFilePassword("PAUL_SQUIRES")


    'specify whether to print the gridlines or not
    'this should come before the setting of fonts and margins
    stat& = xlsPrintGridLines(%XLSTRUE, xlsFileNumber&)


    'it is a good idea to set margins, fonts and column widths
    'prior to writing any text/numerics to the spreadsheet. These
    'should come before setting the fonts.

    stat& = xlsSetMargin(%xlsTopMargin, 1.5, xlsFileNumber&)   'set to 1.5 inches
    stat& = xlsSetMargin(%xlsLeftMargin, 1.5, xlsFileNumber&)
    stat& = xlsSetMargin(%xlsRightMargin, 1.5, xlsFileNumber&)
    stat& = xlsSetMargin(%xlsBottomMargin, 1.5, xlsFileNumber&)


    'Up to 4 fonts can be specified for the spreadsheet. This is a
    'limitation of the Excel 2.1 format. For each value written to the
    'spreadsheet you can specify which font to use.

    stat& = xlsSetFont("Arial", 10, %xlsNoFormat, xlsFileNumber&)             'font0
    stat& = xlsSetFont("Arial", 10, %xlsBold, xlsFileNumber&)                 'font1
    stat& = xlsSetFont("Arial", 10, %xlsBold + %xlsUnderline, xlsFileNumber&) 'font2
    stat& = xlsSetFont("Courier", 18, %xlsItalic, xlsFileNumber&)             'font3

    'Column widths are specified in Excel as 1/256th of a character.
    'xlsSetColumnWidth(BYVAL FirstColumn&, BYVAL LastColumn&, BYVAL WidthValue&, BYVAL xlsFileNumber&) AS LONG
    stat& = xlsSetColumnWidth(1, 7, 20, xlsFileNumber&)
    stat& = xlsSetColumnWidth(8, 8, 12, xlsFileNumber&)

    'set the global row height for the entire spreadsheet
    stat& = xlsSetDefaultRowHeight(16, xlsFileNumber&)

    'set the height of the first two rows a little bigger to allow for the
    'title of the spreadsheet.
    'stat& = xlsSetRowHeight(1, 24)
    'stat& = xlsSetRowHeight(2, 24)

    'set any header or footer that you want to print on
    'every page. This text will be centered at the top and/or
    'bottom of each page. The font will always be the font that
    'is specified as font0, therefore you should only set the
    'header/footer after specifying the fonts through SetFont.
    T$ = "BLM version " + VersionID$
    stat& = xlsSetHeader(T$, xlsFileNumber&)
    T$ = "BLM Output File, " + mFileName$ + ", created " + DATE$
    stat& = xlsSetFooter(T$, xlsFileNumber&)
END FUNCTION

FUNCTION PATHNAMEBLM$ (FullName$)
'
' Strips the file name off the end of a full path and returns just the path
DIM i AS INTEGER
DIM n AS INTEGER
n = LEN(FullName$)
DO
   IF MID$(FullName$, n, 1) = "\" THEN
      EXIT DO
   END IF
   n = n - 1
LOOP UNTIL n = 0

IF n > 0 THEN
   '
   ' Return the path
   PATHNAMEBLM$ = LEFT$(FullName$, n)
ELSE
   '
   ' No path was found, return an empty string
   PATHNAMEBLM$ = ""
END IF
END FUNCTION

FUNCTION SignS% (X!)
   IF (X! < 0) THEN
      SignS% = -1
   ELSE
      SignS% = 1
   END IF
END FUNCTION

SUB WriteDeriv (iHand, Iter, a(), DLC50(), n)

LabelMask$ = " \        \"
OutputMask$ = " ##.###^^^^"

'
' On the first line of output, write column labels
IF (Iter = 1) THEN
   'PRINT #iHand, USING "\                   \"; "Data Label";
   PRINT #iHand, FormatStr$("Data Label", "\                   \");

   FOR c% = 1 TO n - 1
      'IF ((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL)) THEN
      IF (c% = 1) THEN
         TmpS$ = "pH"
      ELSE
         TmpS$ = Label(c%)
      END IF
      IF (c% <> iMetal) THEN
         'PRINT #iHand, USING LabelMask$; TmpS$;
         PRINT #iHand, FormatStr$(TmpS$, LabelMask$);
      END IF
   NEXT c%
   FOR c% = 1 TO n - 1
      'IF ((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL)) THEN
      IF (c% <> iMetal) THEN
         'PRINT #iHand, USING LabelMask$; Label(c%);
         PRINT #iHand, FormatStr$(Label(c%), LabelMask$);
      END IF
   NEXT c%
   PRINT #iHand,

   IF FALSE THEN
   '
   ' Second line, print units
   '
   'PRINT #iHand, USING "\                   \"; "          ";
   PRINT #iHand, FormatStr$("          ", "\                   \");
   'PRINT #iHand, USING LabelMask$; "          ";
   PRINT #iHand, FormatStr$("          ", LabelMask$);
   'PRINT #iHand, USING LabelMask$; "          ";
   PRINT #iHand, FormatStr$("          ", LabelMask$);
   'PRINT #iHand, USING LabelMask$; "mol/L     ";
   PRINT #iHand, FormatStr$("mol/L     ", LabelMask$);
   'PRINT #iHand, USING LabelMask$; "mol/L     ";
   PRINT #iHand, FormatStr$("mol/L     ", LabelMask$);
   'PRINT #iHand, USING LabelMask$; "mol/L     ";
   PRINT #iHand, FormatStr$("mol/L     ", LabelMask$);
   'PRINT #iHand, USING LabelMask$; "mol/L     ";
   PRINT #iHand, FormatStr$("mol/L     ", LabelMask$);

   FOR j = 1 TO NBLM
      'PRINT #iHand, USING LabelMask$; "nmol/g wet";
      PRINT #iHand, FormatStr$("nmol/g wet", LabelMask$);
   NEXT j
   'PRINT #iHand, USING LabelMask$; "mg/L      ";
   PRINT #iHand, FormatStr$("mg/L   ", LabelMask$);
   'PRINT #iHand, USING LabelMask$; "          ";
   PRINT #iHand, FormatStr$("     ", LabelMask$);
   FOR c% = 1 TO 12
      IF ((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL)) THEN
'         PRINT #iHand, USING LabelMask$; "mol/L";
          PRINT #iHand, FormatStr$("mol/L", LabelMask$);
     END IF
   NEXT c%
   PRINT #iHand,
   END IF
END IF


'
'
' Write the actual data
'
'PRINT #iHand, USING "\                   \"; DataL$;
PRINT #iHand, FormatStr$(DataL$, "\                   \");

LastLC50# = DLC50(n)
FOR j = 1 TO n - 1
   IF (j <> iMetal) THEN
      Tmp! = (DLC50(j) - LastLC50#) / (a(j) * .1#)

      IF (j = 1) THEN
         ' Convert H to pH
         dpH! = -LOG10S!(CSNG(OC(1))) * .1#
         dpH! = -LOG10S!(CSNG(OC(1))) - -LOG10S!(CSNG(OC(1) * .1#))
         dpH! = -LOG10S!(CSNG(OC(1))) - -LOG10S!(CSNG(OC(1) * (1! - .1#)))
         Tmp! = (DLC50(j) - LastLC50#) / dpH!
      ELSEIF (j = iDOC) THEN
         ' Convert DOC from g OM/liter to mg C/L
         Tmp! = Tmp! * .001 * 2
      END IF
      'PRINT #iHand, USING OutputMask$; Tmp!;
      PRINT #iHand, FORMAT$(Tmp!, OutputMask$);
   END IF
NEXT j
FOR j = 1 TO n - 1
   IF (j <> iMetal) THEN
      Tmp! = DLC50(j)
      'PRINT #iHand, USING OutputMask$; Tmp!;
      PRINT #iHand, FORMAT$(Tmp!, OutputMask$);
   END IF
NEXT j
PRINT #iHand,


END SUB



SUB WriteDetailed (iHand, Iter)

LabelMask$ = " \       \"
OutputMask$ = "##.###^^^^"


'
'
' On the first line of output, write column labels
IF (Iter = 1) THEN
   'PRINT #iHand, USING "\                   \"; "Site Label";
   PRINT #iHand, FormatStr$("Site Label", "\                   \");
   'PRINT #iHand, USING "\                   \"; "Data Label";
   PRINT #iHand, FormatStr$("Data Label", "\                   \");

   FOR s% = 1 TO NSpecies
      IF (ActCorr(s%) <> acWHAM) THEN
         L$ = Label$(s%)
         'PRINT #iHand, USING LabelMask$; L$;
         PRINT #iHand, FormatStr$(L$, LabelMask$);
      END IF
   NEXT s%

   FOR p% = 1 TO NPhase
      L$ = "MC." + Label$(NSpecies + p%)
      'PRINT #iHand, USING LabelMask$; L$;
      PRINT #iHand, FormatStr$(L$, LabelMask$);
   NEXT p%

   ' Write label for charge column
   L$ = "Charge"
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);

   ' Write label for ionic strength
   L$ = "Ionic S."
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);

   ' Write label for temperature
   L$ = "Temp (K)"
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);

   ' Write label for water volume
   L$ = "Water (L)"
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);

   ' Write label for BL mass
   L$ = "BL (kg wet)"
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);

   ' Write label for number of iterations
   L$ = "# Iter."
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);

   FOR c% = 1 TO NComp
      IF (ActCorr(c%) <> acWHAM) THEN
         L$ = "T." + Label$(c%)
         'PRINT #iHand, USING LabelMask$; L$;
         PRINT #iHand, FormatStr$(L$, LabelMask$);
      END IF
   NEXT c%

   FOR p% = 1 TO NPhase
      L$ = "T." + Label$(NSpecies + p%)
      'PRINT #iHand, USING LabelMask$; L$;
      PRINT #iHand, FormatStr$(L$, LabelMask$);
   NEXT p%

   PRINT #iHand,
END IF


'
'
' Write the actual data
'
   'PRINT #iHand, USING "\                   \"; SiteL$;
   'PRINT #iHand, USING "\                   \"; DataL$;
   PRINT #iHand, FormatStr$(SiteL$, "\                   \");
   PRINT #iHand, FormatStr$(DataL$, "\                   \");
   FOR s% = 1 TO NSpecies
      IF (ActCorr(s%) <> acWHAM) THEN
         'PRINT #iHand, USING OutputMask$; CSNG(OC(s%));
         PRINT #iHand, FORMAT$(CSNG(OC(s%)), OutputMask$);
      END IF
   NEXT s%

   FOR p% = 1 TO NPhase
      'PRINT #iHand, USING OutputMask$; CSNG(ChangeMass(p));
      PRINT #iHand, FORMAT$(CSNG(ChangeMass(p)), OutputMask$);
   NEXT p%

   ' Write charge balance to output file
   'PRINT #iHand, USING OutputMask$; CSNG(OC(NSpecies + 1));
   PRINT #iHand, FORMAT$(CSNG(OC(NSpecies + 1)), OutputMask$);

   ' Write ionic strength to output file
   'PRINT #iHand, USING OutputMask$; CSNG(OC(NSpecies + 2));
   PRINT #iHand, FORMAT$(CSNG(OC(NSpecies + 2)), OutputMask$);

   ' Write temperature to output file
   'PRINT #iHand, USING OutputMask$; CSNG(SysTemp);
   PRINT #iHand, FORMAT$(CSNG(SysTemp), OutputMask$);

   ' Write water volume to output file
   'PRINT #iHand, USING OutputMask$; CSNG(AqueousCtoM!);
   PRINT #iHand, FORMAT$(CSNG(AqueousCtoM!), OutputMask$);

   ' Write BL mass to output file
   'PRINT #iHand, USING OutputMask$; CSNG(SurfaceCtoM!);
   PRINT #iHand, FORMAT$(CSNG(SurfaceCtoM!), OutputMask$);

   ' Write the # iterations to output file
   'PRINT #iHand, USING OutputMask$; CINT(OC(0));
   PRINT #iHand, FORMAT$(CINT(OC(0)), OutputMask$);

   FOR c% = 1 TO NComp
      IF (ActCorr(c%) <> acWHAM) THEN
         Temp# = OKT(c%) '/ CtoM!(c)
         'PRINT #iHand, USING OutputMask$; CSNG(Temp#);
         PRINT #iHand, FORMAT$(CSNG(Temp#), OutputMask$);
      END IF
   NEXT c%

   FOR p% = 1 TO NPhase
      'PRINT #iHand, USING OutputMask$; CSNG(PhaseMass(p%));
      PRINT #iHand, FORMAT$(CSNG(PhaseMass(p%)), OutputMask$);
   NEXT p%

IF FALSE THEN
   FOR c% = 1 TO NSpecies
'      PRINT #iHand, USING OutputMask$; ActCoef(c);
   NEXT c%

   FOR s = NFirstS TO NSpecies
      Temp# = CDBL(hEC(s))
'      PRINT #iHand, USING OutputMask$; Temp#;
   NEXT s

   FOR c% = 1 TO NComp
'      PRINT #iHand, USING OutputMask$; R(c%);
   NEXT c%

   FOR c% = 1 TO NPhase
      Temp# = IAP#(c%)
'      PRINT #iHand, USING OutputMask$; Temp#;
   NEXT c%

   FOR c% = 1 TO NPhase
      Temp# = LOG10S!(CSNG(IAP#(c%) / SEC(c%)))
'      PRINT #iHand, USING OutputMask$; Temp#;
   NEXT c%
END IF
PRINT #iHand,

END SUB
SUB WriteDetailed2 (iHand&, Iter, iComplete)

DIM TmpS AS STRING
DIM Tmp#

LabelMask$ =  "\               \"
OutputMask$ = "       ##.###^^^^"
'
'
' On the first line of output, write column labels
IF (Iter = 1) THEN
   PRINT #iHand&, FormatStr$("Site Label", "\                    \");
   PRINT #iHand&, FormatStr$("Sample Label", "\                    \");

   FOR s% = 1 TO NSpecies
      IF (ActCorr(s%) <> acWHAM) THEN
         L$ = Label$(s%)
         PRINT #iHand&, FormatStr$(L$, LabelMask$);
      END IF
   NEXT s%

   IF (iWHAM = TRUE) THEN
       FOR s% = 1 TO NComp
          IF (V(s%)>0) THEN
             IF (SepNOM = TRUE) THEN
                L$ = "FA." + Label$(s%)
                PRINT #iHand&, FormatStr$(L$, LabelMask$);

                L$ = "HA." + Label$(s%)
                PRINT #iHand&, FormatStr$(L$, LabelMask$);
             END IF
             L$ = "TOrg." + Label$(s%)
             PRINT #iHand&, FormatStr$(L$, LabelMask$);
          END IF
       NEXT s%
   END IF

   FOR p% = 1 TO NPhase
      L$ = "MC." + Label$(NSpecies + p%)
      PRINT #iHand&, FormatStr$(L$, LabelMask$);
   NEXT p%

   ' Write label for charge column
   L$ = "Charge"
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   ' Write label for ionic strength
   L$ = "Ionic S."
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   ' Write labels for activity corrections  (added to ver 2.45, build 2021-04-05)
   FOR s% = 1 TO 4
     iCol = iCol + 1
     L$ = "Act_z"+TRIM$(STR$(s%))
     PRINT #iHand&, FormatStr$(L$, LabelMask$);
   NEXT s%

   ' Write label for temperature
   L$ = "Temp (K)"
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   ' Write label for water volume
   L$ = "Water (L)"
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   ' Write label for BL mass
   L$ = "BL (kg wet)"
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   ' Write label for number of iterations
   L$ = "# Iter."
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   FOR c% = 1 TO NComp
      IF (ActCorr(c%) <> acWHAM) THEN
         L$ = "T." + Label$(c%)
         PRINT #iHand&, FormatStr$(L$, LabelMask$);
      END IF
   NEXT c%

   FOR p% = 1 TO NPhase
      L$ = "T." + Label$(NSpecies + p%)
      PRINT #iHand&, FormatStr$(L$, LabelMask$);
   NEXT p%

   PRINT #iHand&,
END IF


'
'
' Write the actual data
'
   PRINT #iHand&, FormatStr$(SiteL$, "\                    \");
   PRINT #iHand&, FormatStr$(DataL$, "\                    \");
   FOR s% = 1 TO NSpecies
      IF (ActCorr(s%) <> acWHAM) THEN
         IF iComplete AND NOT ErrorNotConverged THEN
            Tmp# = OC(s%)
            PRINT #iHand&, USING$(OutputMask$, Tmp#);
         ELSE
            PRINT #iHand&, "               NA";
         END IF
      END IF
   NEXT s%

   IF (iWHAM = TRUE) THEN
       FOR s% = 1 TO NComp
          IF (V(s%)>0) THEN
              IF (SepNOM = TRUE) THEN
                IF iComplete AND NOT ErrorNotConverged THEN
                    ' FA Fraction
                    Tmp# = 0!
                    FOR j = 1 TO NSpecies
                      Tmp# = Tmp# + FABound(j) * SCb(j, s%)
                    NEXT j
                    PRINT #iHand&, USING$(OutputMask$, Tmp#);
                    ' HA Fraction
                    Tmp# = 0!
                    FOR j = 1 TO NSpecies
                      Tmp# = Tmp# + HABound(j) * SCb(j, s%)
                    NEXT j
                    PRINT #iHand&, USING$(OutputMask$, Tmp#);
                ELSE  ' if not complete
                    Tmp# = -999.#
                    PRINT #iHand&, "               NA";
                    PRINT #iHand&, "               NA";
                END IF
              END IF
              IF iComplete AND NOT ErrorNotConverged THEN
                  ' Combined FA and HA as Total Organic
                  Tmp# = 0!
                  FOR j = 1 TO NSpecies
                    Tmp# = Tmp# + wOrgBound!(j) * SCb(j, s%)
                  NEXT j
                  PRINT #iHand&, USING$(OutputMask$, Tmp#);
              ELSE  ' if not complete
                  Tmp# = -999.#
                  PRINT #iHand&, , "               NA";
              END IF
          END IF
       NEXT s%
   END IF


   FOR p% = 1 TO NPhase
      IF iComplete AND NOT ErrorNotConverged THEN
         Tmp# = ChangeMass(p%)
         PRINT #iHand&, USING$(OutputMask$, Tmp#);
      ELSE
         PRINT #iHand&, , "       NA";
      END IF
   NEXT p%

   ' Write charge balance to output file
   IF iComplete AND NOT ErrorNotConverged THEN
      Tmp# = CSNG(OC(NSpecies + 1))
      PRINT #iHand&, USING$(OutputMask$, Tmp#);
   ELSE
      PRINT #iHand&, "               NA";
   END IF

   ' Write ionic strength to output file
   IF iComplete AND NOT ErrorNotConverged THEN
      Tmp# = OC(NSpecies + 2)
      PRINT #iHand&, USING$(OutputMask$, Tmp#);
   ELSE
      PRINT #iHand&, "               NA";
   END IF

   ' Write activity corrections (added to ver 2.45, build 2021-04-05)
   FOR s% = 1 TO 4
       IF iComplete AND NOT ErrorNotConverged THEN
         iCol = iCol + 1
         Tmp# = wGAMMA(s%)
         PRINT #iHand&, USING$(OutputMask$, Tmp#);
       ELSE  ' if not complete
          PRINT #iHand&, "               NA";
       END IF
   NEXT s%

   ' Write temperature to output file

   IF iComplete AND NOT ErrorNotConverged THEN
      Tmp# = SysTemp
      PRINT #iHand&, USING$(OutputMask$, Tmp#);
   ELSE
      PRINT #iHand&, "               NA";
   END IF

   ' Write water volume to output file
   IF iComplete AND NOT ErrorNotConverged THEN
      Tmp# = AqueousCtoM!
      PRINT #iHand&, USING$(OutputMask$, Tmp#);
   ELSE
      PRINT #iHand&, "               NA";
   END IF

   ' Write BL mass to output file
   IF iComplete AND NOT ErrorNotConverged THEN
      Tmp# = SurfaceCtoM!
      PRINT #iHand&, USING$(OutputMask$, Tmp#);
   ELSE
      PRINT #iHand&, "               NA";
   END IF


   ' Write the # iterations to output file
   IF iComplete AND NOT ErrorNotConverged THEN
      Tmp# = OC(0)
      PRINT #iHand&, USING$(OutputMask$, Tmp#);
   ELSE
      PRINT #iHand&, "               NA";
   END IF

   FOR c% = 1 TO NComp
      IF (ActCorr(c%) <> acWHAM) THEN
         IF iComplete AND NOT ErrorNotConverged THEN
            Tmp# = OKT(c%)
            IF CType(c%)= ctFixed THEN
                Tmp# = KT(c%)
            END IF
            PRINT #iHand&, USING$(OutputMask$, Tmp#);
         ELSE
            PRINT #iHand&, "               NA";
         END IF
      END IF
   NEXT c%

   FOR p% = 1 TO NPhase
       IF iComplete AND NOT ErrorNotConverged THEN
          Tmp# = PhaseMass(p%)
          PRINT #iHand&, USING$(OutputMask$, Tmp#);
       ELSE
          PRINT #iHand&, "               NA";
       END IF
   NEXT p%

IF FALSE THEN
   FOR c% = 1 TO NSpecies
'      PRINT #iHand&, USING OutputMask$; ActCoef(c);
   NEXT c%

   FOR s = NFirstS TO NSpecies
      Temp# = CDBL(hEC(s))
'      PRINT #iHand&, USING OutputMask$; Temp#;
   NEXT s

   FOR c% = 1 TO NComp
'      PRINT #iHand&, USING OutputMask$; R(c%);
   NEXT c%

   FOR c% = 1 TO NPhase
      Temp# = IAP#(c%)
'      PRINT #iHand&, USING OutputMask$; Temp#;
   NEXT c%

   FOR c% = 1 TO NPhase
      Temp# = LOG10S!(CSNG(IAP#(c%) / SEC(c%)))
'      PRINT #iHand&, USING OutputMask$; Temp#;
   NEXT c%
END IF
PRINT #iHand&,

END SUB

SUB WriteDetailedComma (iHand&, Iter, iComplete)

DIM TmpS AS STRING
DIM Tmp#

LabelMask$ =  ""
OutputMask1$ = " ##.##############^^^^_,"
OutputMask2$ = " ##.##############_,"
'
'
' On the first line of output, write column labels
IF (Iter = 1) THEN
   PRINT #iHand&, FormatStr$("Site Label", LabelMask$);
   PRINT #iHand&, FormatStr$("Sample Label", LabelMask$);

   FOR s% = 1 TO NSpecies
      IF (ActCorr(s%) <> acWHAM) THEN
         L$ = Label$(s%)
         PRINT #iHand&, FormatStr$(L$, LabelMask$);
      END IF
   NEXT s%

   IF (iWHAM = TRUE) THEN
       FOR s% = 1 TO NComp
          IF (V(s%)>0) THEN
             IF (SepNOM = TRUE) THEN
                L$ = "FA." + Label$(s%)
                PRINT #iHand&, FormatStr$(L$, LabelMask$);

                L$ = "HA." + Label$(s%)
                PRINT #iHand&, FormatStr$(L$, LabelMask$);
             END IF
             L$ = "TOrg." + Label$(s%)
             PRINT #iHand&, FormatStr$(L$, LabelMask$);
          END IF
       NEXT s%
   END IF

   FOR p% = 1 TO NPhase
      L$ = "MC." + Label$(NSpecies + p%)
      PRINT #iHand&, FormatStr$(L$, LabelMask$);
   NEXT p%

   ' Write label for charge column
   L$ = "Charge"
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   ' Write label for ionic strength
   L$ = "Ionic S."
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   ' Write labels for activity corrections  (added to ver 2.45, build 2021-04-05)
   FOR s% = 1 TO 4
     iCol = iCol + 1
     L$ = "Act_z"+TRIM$(STR$(s%))
     PRINT #iHand&, FormatStr$(L$, LabelMask$);
   NEXT s%

   ' Write label for temperature
   L$ = "Temp (K)"
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   ' Write label for water volume
   L$ = "Water (L)"
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   ' Write label for BL mass
   L$ = "BL (kg wet)"
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   ' Write label for number of iterations
   L$ = "# Iter."
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   FOR c% = 1 TO NComp
      IF (ActCorr(c%) <> acWHAM) THEN
         L$ = "T." + Label$(c%)
         PRINT #iHand&, FormatStr$(L$, LabelMask$);
      END IF
   NEXT c%

   FOR p% = 1 TO NPhase
      L$ = "T." + Label$(NSpecies + p%)
      PRINT #iHand&, FormatStr$(L$, LabelMask$);
   NEXT p%

   PRINT #iHand&,
END IF


'
'
' Write the actual data
'
   PRINT #iHand&, FormatStr$(SiteL$, LabelMask$);
   PRINT #iHand&, FormatStr$(DataL$, LabelMask$);
   FOR s% = 1 TO NSpecies
      IF (ActCorr(s%) <> acWHAM) THEN
         IF iComplete AND NOT ErrorNotConverged THEN
            Tmp# = OC(s%)
            IF (Tmp# > 9999) OR (Tmp# < 0.01) THEN
               OutputMask$ = OutputMask1$
            ELSE
               OutputMask$ = OutputMask2$
            END IF
            PRINT #iHand&, USING$(OutputMask$, Tmp#);
         ELSE
            PRINT #iHand&, FormatStr$("NA", LabelMask$);
         END IF
      END IF
   NEXT s%

   IF (iWHAM = TRUE) THEN
       FOR s% = 1 TO NComp
          IF (V(s%)>0) THEN
              IF (SepNOM = TRUE) THEN
                IF iComplete AND NOT ErrorNotConverged THEN
                    ' FA Fraction
                    Tmp# = 0!
                    FOR j = 1 TO NSpecies
                      Tmp# = Tmp# + FABound(j) * SCb(j, s%)
                    NEXT j
                    IF (Tmp# > 9999) OR (Tmp# < 0.01) THEN
                       OutputMask$ = OutputMask1$
                    ELSE
                       OutputMask$ = OutputMask2$
                    END IF
                    PRINT #iHand&, USING$(OutputMask$, Tmp#);
                    ' HA Fraction
                    Tmp# = 0!
                    FOR j = 1 TO NSpecies
                      Tmp# = Tmp# + HABound(j) * SCb(j, s%)
                    NEXT j
                    IF (Tmp# > 9999) OR (Tmp# < 0.01) THEN
                       OutputMask$ = OutputMask1$
                    ELSE
                       OutputMask$ = OutputMask2$
                    END IF
                    PRINT #iHand&, USING$(OutputMask$, Tmp#);
                ELSE  ' if not complete
                    Tmp# = -999.#
                    PRINT #iHand&, FormatStr$("NA", LabelMask$);
                    PRINT #iHand&, FormatStr$("NA", LabelMask$);
                END IF
              END IF
              IF iComplete AND NOT ErrorNotConverged THEN
                  ' Combined FA and HA as Total Organic
                  Tmp# = 0!
                  FOR j = 1 TO NSpecies
                    Tmp# = Tmp# + wOrgBound!(j) * SCb(j, s%)
                  NEXT j
                  IF (Tmp# > 9999) OR (Tmp# < 0.01) THEN
                     OutputMask$ = OutputMask1$
                  ELSE
                     OutputMask$ = OutputMask2$
                  END IF
                  PRINT #iHand&, USING$(OutputMask$, Tmp#);
              ELSE  ' if not complete
                  Tmp# = -999.#
                  PRINT #iHand&, , FormatStr$("NA", LabelMask$);
              END IF
          END IF
       NEXT s%
   END IF


   FOR p% = 1 TO NPhase
      IF iComplete AND NOT ErrorNotConverged THEN
         Tmp# = ChangeMass(p%)
         IF (Tmp# > 9999) OR (Tmp# < 0.01) THEN
            OutputMask$ = OutputMask1$
         ELSE
            OutputMask$ = OutputMask2$
         END IF
         PRINT #iHand&, USING$(OutputMask$, Tmp#);
      ELSE
         PRINT #iHand&, FormatStr$("NA", LabelMask$);
      END IF
   NEXT p%

   ' Write charge balance to output file
   IF iComplete AND NOT ErrorNotConverged THEN
      Tmp# = CSNG(OC(NSpecies + 1))
      IF (Tmp# > 9999) OR (Tmp# < 0.01) THEN
         OutputMask$ = OutputMask1$
      ELSE
         OutputMask$ = OutputMask2$
      END IF
      PRINT #iHand&, USING$(OutputMask$, Tmp#);
   ELSE
      PRINT #iHand&, FormatStr$("NA", LabelMask$);
   END IF

   ' Write ionic strength to output file
   IF iComplete AND NOT ErrorNotConverged THEN
      Tmp# = OC(NSpecies + 2)
      IF (Tmp# > 9999) OR (Tmp# < 0.01) THEN
         OutputMask$ = OutputMask1$
      ELSE
         OutputMask$ = OutputMask2$
      END IF
      PRINT #iHand&, USING$(OutputMask$, Tmp#);
   ELSE
      PRINT #iHand&, FormatStr$("NA", LabelMask$);
   END IF


   ' Write activity corrections (added to ver 2.45, build 2021-04-05)
   FOR s% = 1 TO 4
       IF iComplete AND NOT ErrorNotConverged THEN
          Tmp# = wGAMMA(s%)
          IF (Tmp# > 9999) OR (Tmp# < 0.01) THEN
             OutputMask$ = OutputMask1$
          ELSE
             OutputMask$ = OutputMask2$
          END IF
          PRINT #iHand&, USING$(OutputMask$, Tmp#);
       ELSE
          PRINT #iHand&, FormatStr$("NA", LabelMask$);
       END IF
   NEXT s%

   ' Write temperature to output file

   IF iComplete AND NOT ErrorNotConverged THEN
      Tmp# = SysTemp
      IF (Tmp# > 9999) OR (Tmp# < 0.01) THEN
         OutputMask$ = OutputMask1$
      ELSE
         OutputMask$ = OutputMask2$
      END IF
      PRINT #iHand&, USING$(OutputMask$, Tmp#);
   ELSE
      PRINT #iHand&, FormatStr$("NA", LabelMask$);
   END IF

   ' Write water volume to output file
   IF iComplete AND NOT ErrorNotConverged THEN
      Tmp# = AqueousCtoM!
      IF (Tmp# > 9999) OR (Tmp# < 0.01) THEN
         OutputMask$ = OutputMask1$
      ELSE
         OutputMask$ = OutputMask2$
      END IF
      PRINT #iHand&, USING$(OutputMask$, Tmp#);
   ELSE
      PRINT #iHand&, FormatStr$("NA", LabelMask$);
   END IF

   ' Write BL mass to output file
   IF iComplete AND NOT ErrorNotConverged THEN
      Tmp# = SurfaceCtoM!
      IF (Tmp# > 9999) OR (Tmp# < 0.01) THEN
         OutputMask$ = OutputMask1$
      ELSE
         OutputMask$ = OutputMask2$
      END IF
      PRINT #iHand&, USING$(OutputMask$, Tmp#);
   ELSE
      PRINT #iHand&, FormatStr$("NA", LabelMask$);
   END IF


   ' Write the # iterations to output file
   IF iComplete AND NOT ErrorNotConverged THEN
      Tmp# = OC(0)
      IF (Tmp# > 9999) OR (Tmp# < 0.01) THEN
         OutputMask$ = OutputMask1$
      ELSE
         OutputMask$ = OutputMask2$
      END IF
      PRINT #iHand&, USING$(OutputMask$, Tmp#);
   ELSE
      PRINT #iHand&, FormatStr$("NA", LabelMask$);
   END IF

   FOR c% = 1 TO NComp
      IF (ActCorr(c%) <> acWHAM) THEN
         IF iComplete AND NOT ErrorNotConverged THEN
            Tmp# = OKT(c%)
            IF CType(c%)= ctFixed THEN
                Tmp# = KT(c%)
            END IF
            IF (Tmp# > 9999) OR (Tmp# < 0.01) THEN
               OutputMask$ = OutputMask1$
            ELSE
               OutputMask$ = OutputMask2$
            END IF
            PRINT #iHand&, USING$(OutputMask$, Tmp#);
         ELSE
            PRINT #iHand&, FormatStr$("NA", LabelMask$);
         END IF
      END IF
   NEXT c%

   FOR p% = 1 TO NPhase
       IF iComplete AND NOT ErrorNotConverged THEN
          Tmp# = PhaseMass(p%)
          IF (Tmp# > 9999) OR (Tmp# < 0.01) THEN
             OutputMask$ = OutputMask1$
          ELSE
             OutputMask$ = OutputMask2$
          END IF
          PRINT #iHand&, USING$(OutputMask$, Tmp#);
       ELSE
          PRINT #iHand&, FormatStr$("NA", LabelMask$);
       END IF
   NEXT p%

IF FALSE THEN
   FOR c% = 1 TO NSpecies
'      PRINT #iHand&, USING OutputMask$; ActCoef(c);
   NEXT c%

   FOR s = NFirstS TO NSpecies
      Temp# = CDBL(hEC(s))
'      PRINT #iHand&, USING OutputMask$; Temp#;
   NEXT s

   FOR c% = 1 TO NComp
'      PRINT #iHand&, USING OutputMask$; R(c%);
   NEXT c%

   FOR c% = 1 TO NPhase
      Temp# = IAP#(c%)
'      PRINT #iHand&, USING OutputMask$; Temp#;
   NEXT c%

   FOR c% = 1 TO NPhase
      Temp# = LOG10S!(CSNG(IAP#(c%) / SEC(c%)))
'      PRINT #iHand&, USING OutputMask$; Temp#;
   NEXT c%
END IF
PRINT #iHand&,

END SUB

SUB WriteDetailedOld (iHand&, Iter, iComplete)

LabelMask$ = " \       \"
OutputMask$ = "##.###^^^^"


'
'
' On the first line of output, write column labels
IF (Iter = 1) THEN
   'PRINT #iHand, USING "\                   \"; "Site Label";
   PRINT #iHand, FormatStr$("Site Label", "\                   \");
   'PRINT #iHand, USING "\                   \"; "Data Label";
   PRINT #iHand, FormatStr$("Data Label", "\                   \");

   FOR s% = 1 TO NSpecies
      IF (ActCorr(s%) <> acWHAM) THEN
         L$ = Label$(s%)
         'PRINT #iHand, USING LabelMask$; L$;
         PRINT #iHand, FormatStr$(L$, LabelMask$);
      END IF
   NEXT s%

   FOR p% = 1 TO NPhase
      L$ = "MC." + Label$(NSpecies + p%)
      'PRINT #iHand, USING LabelMask$; L$;
      PRINT #iHand, FormatStr$(L$, LabelMask$);
   NEXT p%

   ' Write label for charge column
   L$ = "Charge"
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);

   ' Write label for ionic strength
   L$ = "Ionic S."
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);

   ' Write label for temperature
   L$ = "Temp (K)"
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);

   ' Write label for water volume
   L$ = "Water (L)"
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);

   ' Write label for BL mass
   L$ = "BL (kg wet)"
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);

   ' Write label for number of iterations
   L$ = "# Iter."
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);

   FOR c% = 1 TO NComp
      IF (ActCorr(c%) <> acWHAM) THEN
         L$ = "T." + Label$(c%)
         'PRINT #iHand, USING LabelMask$; L$;
         PRINT #iHand, FormatStr$(L$, LabelMask$);
      END IF
   NEXT c%

   FOR p% = 1 TO NPhase
      L$ = "T." + Label$(NSpecies + p%)
      'PRINT #iHand, USING LabelMask$; L$;
      PRINT #iHand, FormatStr$(L$, LabelMask$);
   NEXT p%

   PRINT #iHand,
END IF


'
'
' Write the actual data
'
   'PRINT #iHand, USING "\                   \"; SiteL$;
   'PRINT #iHand, USING "\                   \"; DataL$;
   PRINT #iHand, FormatStr$(SiteL$, "\                   \");
   PRINT #iHand, FormatStr$(DataL$, "\                   \");
   IF NOT iComplete THEN
      ' Only write the site and data label for incomplete obs
      PRINT #iHand, FormatStr$("Incomplete", "\                   \");
      PRINT #iHand,
      EXIT SUB
   END IF
   FOR s% = 1 TO NSpecies
      IF (ActCorr(s%) <> acWHAM) THEN
         'PRINT #iHand, USING OutputMask$; CSNG(OC(s%));
         PRINT #iHand, FORMAT$(CSNG(OC(s%)), OutputMask$);
      END IF
   NEXT s%

   FOR p% = 1 TO NPhase
      'PRINT #iHand, USING OutputMask$; CSNG(ChangeMass(p));
      PRINT #iHand, FORMAT$(CSNG(ChangeMass(p)), OutputMask$);
   NEXT p%

   ' Write charge balance to output file
   'PRINT #iHand, USING OutputMask$; CSNG(OC(NSpecies + 1));
   PRINT #iHand, FORMAT$(CSNG(OC(NSpecies + 1)), OutputMask$);

   ' Write ionic strength to output file
   'PRINT #iHand, USING OutputMask$; CSNG(OC(NSpecies + 2));
   PRINT #iHand, FORMAT$(CSNG(OC(NSpecies + 2)), OutputMask$);

   ' Write temperature to output file
   'PRINT #iHand, USING OutputMask$; CSNG(SysTemp);
   PRINT #iHand, FORMAT$(CSNG(SysTemp), OutputMask$);

   ' Write water volume to output file
   'PRINT #iHand, USING OutputMask$; CSNG(AqueousCtoM!);
   PRINT #iHand, FORMAT$(CSNG(AqueousCtoM!), OutputMask$);

   ' Write BL mass to output file
   'PRINT #iHand, USING OutputMask$; CSNG(SurfaceCtoM!);
   PRINT #iHand, FORMAT$(CSNG(SurfaceCtoM!), OutputMask$);

   ' Write the # iterations to output file
   'PRINT #iHand, USING OutputMask$; CINT(OC(0));
   PRINT #iHand, FORMAT$(CINT(OC(0)), OutputMask$);

   FOR c% = 1 TO NComp
      IF (ActCorr(c%) <> acWHAM) THEN
         Temp# = OKT(c%) '/ CtoM!(c)
         'PRINT #iHand, USING OutputMask$; CSNG(Temp#);
         PRINT #iHand, FORMAT$(CSNG(Temp#), OutputMask$);
      END IF
   NEXT c%

   FOR p% = 1 TO NPhase
      'PRINT #iHand, USING OutputMask$; CSNG(PhaseMass(p%));
      PRINT #iHand, FORMAT$(CSNG(PhaseMass(p%)), OutputMask$);
   NEXT p%

IF FALSE THEN
   FOR c% = 1 TO NSpecies
'      PRINT #iHand, USING OutputMask$; ActCoef(c);
   NEXT c%

   FOR s = NFirstS TO NSpecies
      Temp# = CDBL(hEC(s))
'      PRINT #iHand, USING OutputMask$; Temp#;
   NEXT s

   FOR c% = 1 TO NComp
'      PRINT #iHand, USING OutputMask$; R(c%);
   NEXT c%

   FOR c% = 1 TO NPhase
      Temp# = IAP#(c%)
'      PRINT #iHand, USING OutputMask$; Temp#;
   NEXT c%

   FOR c% = 1 TO NPhase
      Temp# = LOG10S!(CSNG(IAP#(c%) / SEC(c%)))
'      PRINT #iHand, USING OutputMask$; Temp#;
   NEXT c%
END IF
PRINT #iHand,

END SUB

SUB WriteDetailedXLS (iHand&, Iter, iComplete)

LabelMask$ = " \       \"
OutputMask$ = "##.###^^^^"
DIM iRow AS LONG
DIM iRowOS AS LONG

IF (iDebug = 01 OR iDebug = 08) THEN CALL DebugLog ("wDET", "Activities")
IF (iDebug = 01 OR iDebug = 08) THEN
  FOR iRow = 1 TO 4
    sTmp1$ = "GAMMA" + STR$(iRow)
    sTmp1$ = sTmp1$ + STR$(wGAMMA(iRow))
    CALL DebugLog ("wDET", sTmp1$)
  NEXT X%

  FOR iRow = 1 TO NSP%
    IF (wA(iRow) <> 0) THEN
       sTmp1$ = "A." + wNs$(iRow)
       sTmp1$ = sTmp1$ + STR$(wA(iRow))
       CALL DebugLog ("wDET", sTmp1$)
    END IF
  NEXT X%

  FOR iRow = 1 TO NSP%
    IF (wC(iRow) <> 0) THEN
       sTmp1$ = "C." + wNs$(iRow)
       sTmp1$ = sTmp1$ + STR$(wC(iRow))
       CALL DebugLog ("wDET", sTmp1$)
    END IF
  NEXT X%
END IF

iRowOS = 6

'
'
' On the first line of output, write column labels
IF (Iter = 1) THEN
   iRow = 5
   iCol = 1
   sTmp1$ = "Site Label"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   iCol = 2
   sTmp1$ = "Sample Label"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   FOR s% = 1 TO NSpecies
      IF (ActCorr(s%) <> acWHAM) THEN
         iCol = iCol + 1
         sTmp1$ = Label$(s%)
         stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      END IF
   NEXT s%

   IF (iWHAM = TRUE) THEN
       FOR s% = 1 TO NComp
          IF (V(s%)>0) THEN
             IF (SepNOM = TRUE) THEN
                iCol = iCol + 1
                sTmp1$ = "FA." + Label$(s%)
                stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
                iCol = iCol + 1
                sTmp1$ = "HA." + Label$(s%)
                stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
             END IF
             iCol = iCol + 1
             sTmp1$ = "TOrg." + Label$(s%)
             stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
          END IF
       NEXT s%
   END IF

   FOR p% = 1 TO NPhase
      iCol = iCol + 1
      sTmp1$ = "MC." + Label$(NSpecies + p%)
      stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   NEXT p%

   ' Write label for charge column

   iCol = iCol + 1
   sTmp1$ = "Charge"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write label for ionic strength
   iCol = iCol + 1
   sTmp1$ = "Ionic S."
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write labels for activity corrections  (added to ver 2.45, build 2021-04-05)
   FOR s% = 1 TO 4
     iCol = iCol + 1
     sTmp1$ = "Act_z"+TRIM$(STR$(s%))
     stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   NEXT s%

   ' Write charge for Humics
   iCol = iCol + 1
   sTmp1$ = "Z_HS"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write charge for fulvics
   iCol = iCol + 1
   sTmp1$ = "Z_FS"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write diffuse layer vol for humics
   iCol = iCol + 1
   sTmp1$ = "DDLVol_HS"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write diffuse layer vol for humics
   iCol = iCol + 1
   sTmp1$ = "DDLVol_FS"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write max diffuse layer vol for humics
   iCol = iCol + 1
   sTmp1$ = "DDLMaxV_HS"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write max diffuse layer vol for fulvics
   iCol = iCol + 1
   sTmp1$ = "DDLMaxV_FS"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write mass humics
   iCol = iCol + 1
   sTmp1$ = "TH_HS"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write mass fulvics
   iCol = iCol + 1
   sTmp1$ = "TH_FS"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write Ratio humics
   iCol = iCol + 1
   sTmp1$ = "Ratio_HS"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write Ratio fulvics
   iCol = iCol + 1
   sTmp1$ = "Ratio_FS"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' write DVol
   iCol = iCol + 1
   sTmp1$ = "DVol(1)"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   iCol = iCol + 1
   sTmp1$ = "DVol(2)"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' write individual species concentrations for organically bound substances
    FOR i = 1 TO NSP%
       ii = CMap%(i)
       IF (ii <> 0) THEN
           ' humic bound specific
           iCol = iCol + 1
           sTmp1$ = "HS-" + Label(ii)
           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

           ' humic bound diffuse
           iCol = iCol + 1
           sTmp1$ = "HD-" + Label(ii)
           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

           ' total humic bound
           iCol = iCol + 1
           sTmp1$ = "HT-" + Label(ii)
           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

           ' fulvic bound specific
           iCol = iCol + 1
           sTmp1$ = "FS-" + Label(ii)
           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

           ' fulvic bound diffuse
           iCol = iCol + 1
           sTmp1$ = "FD-" + Label(ii)
           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

          ' total fulvic bound
           iCol = iCol + 1
           sTmp1$ = "FT-" + Label(ii)
           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
       END IF
    NEXT i


   ' Write label for temperature
   iCol = iCol + 1
   sTmp1$ = "Temp (K)"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write label for solution vol (total vol - ddl vol)
   iCol = iCol + 1
   sTmp1$ = "VOLSOL"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write label for BL mass
   FOR i = 1 TO NMass
      ' Write label for the name of each mass compartment
      iCol = iCol + 1
      sTmp1$ = MassName(i)
      stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   NEXT i

   ' Write label for number of iterations
   iCol = iCol + 1
   sTmp1$ = "# Iter."
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   FOR c% = 1 TO NComp
      IF (ActCorr(c%) <> acWHAM) THEN
         iCol = iCol + 1
         sTmp1$ = "T." + Label$(c%)
         stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      END IF
   NEXT c%

   FOR p% = 1 TO NPhase
      iCol = iCol + 1
      sTmp1$ = "T." + Label$(NSpecies + p%)
      stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   NEXT p%


   '   *-------------------------------*
   '   |   Now write the unit labels   |
   '   *-------------------------------*
   iRow = iRow + 1
   iCol = 1              ' The "Site Label" column has no need for units
   iCol = iCol + 1       ' The "Sample Label" column has no need for units

   '   *----------------------------------------------------------------------------------------------------*
   '   |  The species labels have to be constructed depending on which mass compartment this species is in  |
   '   *----------------------------------------------------------------------------------------------------*
      FOR s% = 1 TO NSpecies
      IF (ActCorr(s%) <> acWHAM) THEN
         iCol = iCol + 1
         sTmp1$ = "mol / " + MassUnitLabel(SType(s%) + 1)
         stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      END IF
   NEXT s%

   IF (iWHAM = TRUE) THEN
       FOR s% = 1 TO NComp
          IF (V(s%)>0) THEN
             IF (SepNOM = TRUE) THEN
                iCol = iCol + 1
                sTmp1$ = "mol / " + MassUnitLabel(SType(S) + 1)
                stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
                iCol = iCol + 1
                sTmp1$ = "mol / " + MassUnitLabel(SType(S) + 1)
                stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
             END IF
             iCol = iCol + 1
             sTmp1$ = "mol / " + MassUnitLabel(SType(S) + 1)
             stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
          END IF
       NEXT s%
   END IF

   FOR p% = 1 TO NPhase
      iCol = iCol + 1
      sTmp1$ = "mol"
      stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   NEXT p%

   ' Write units label for charge column
   iCol = iCol + 1
   sTmp1$ = "eq"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write units label for ionic strength
   iCol = iCol + 1
   sTmp1$ = "mol / L"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Skip unit labels for activity corrections  (added to ver 2.45, build 2021-04-05)
   FOR s% = 1 TO 4
     iCol = iCol + 1
   NEXT s%

   ' Write charge for Humics
   iCol = iCol + 1
   sTmp1$ = "charge/(g HS)"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write charge for fulvics
   iCol = iCol + 1
   sTmp1$ = "charge/(g FS)"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write diffuse layer vol for humics
   iCol = iCol + 1
   sTmp1$ = "L/(g HS)"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write diffuse layer vol for humics
   iCol = iCol + 1
   sTmp1$ = "L/(g FS)"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write max diffuse layer vol for humics
   iCol = iCol + 1
   sTmp1$ = "L"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write max diffuse layer vol for fulvics
   iCol = iCol + 1
   sTmp1$ = "L"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write mass humics
   iCol = iCol + 1
   sTmp1$ = "g HS"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write mass fulvics
   iCol = iCol + 1
   sTmp1$ = "g FS"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write Ratio humics
   iCol = iCol + 1
   sTmp1$ = "Ratio_HS"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write Ratio fulvics
   iCol = iCol + 1
   sTmp1$ = "Ratio_FS"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' write DVol
   iCol = iCol + 1
   sTmp1$ = "DVol(1)"
   stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   iCol = iCol + 1
   sTmp1$ = "DVol(2)"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' write individual species concentrations for organically bound substances
    FOR i = 1 TO NSP%
       ii = CMap%(i)
       IF (ii <> 0) THEN
           ' humic bound specific
           iCol = iCol + 1
           sTmp1$ = "HS-" + Label(ii)
           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

           ' humic bound diffuse
           iCol = iCol + 1
           sTmp1$ = "HD-" + Label(ii)
           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

           ' total humic bound
           iCol = iCol + 1
           sTmp1$ = "HT-" + Label(ii)
           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

           ' fulvic bound specific
           iCol = iCol + 1
           sTmp1$ = "FS-" + Label(ii)
           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

           ' fulvic bound diffuse
           iCol = iCol + 1
           sTmp1$ = "FD-" + Label(ii)
           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

          ' total fulvic bound
           iCol = iCol + 1
           sTmp1$ = "FT-" + Label(ii)
           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
       END IF
    NEXT i


   ' Write units label for temperature
   iCol = iCol + 1
   sTmp1$ = "K"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write label for solution vol (total vol - ddl vol)
   iCol = iCol + 1
   sTmp1$ = "L"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   FOR i = 1 TO NMass
      ' Write units label for the name of each mass compartment
      iCol = iCol + 1
      sTmp1$ = MassUnitLabel(i)
      stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   NEXT i

   ' no units are needed for the number of iterations
   iCol = iCol + 1

   '   *---------------------*
   '   |  The total labels   |
   '   *---------------------*
   FOR c% = 1 TO NComp
      IF (ActCorr(c%) <> acWHAM) THEN
         iCol = iCol + 1
         sTmp1$ = "mol"
         stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      END IF
   NEXT c%

   FOR p% = 1 TO NPhase
      iCol = iCol + 1
      sTmp1$ = "mol"
      stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   NEXT p%

   'PRINT #iHand,
END IF


'
'
' Write the actual data
'
   iRow = iRowOS + Iter
   iCol = 1
   sTmp1$ = SiteL$
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   iCol = iCol + 1
   sTmp1$ = DataL$
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   FOR s% = 1 TO NSpecies
      IF (ActCorr(s%) <> acWHAM) THEN
         iCol = iCol + 1
         IF iComplete AND NOT ErrorNotConverged THEN
            value# = OC(s%)
            stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
         ELSE
            stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
         END IF
      END IF
   NEXT s%
   IF (iWHAM = TRUE) THEN
       FOR s% = 1 TO NComp
          IF (V(s%)>0) THEN
              IF (SepNOM = TRUE) THEN
                IF iComplete AND NOT ErrorNotConverged THEN
                    ' FA Fraction
                    iCol = iCol + 1
                    value# = 0!
                    FOR j = 1 TO NSpecies
                      value# = value# + FABound(j) * SCb(j, s%)
                    NEXT j
                    stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
                    ' HA Fraction
                    iCol = iCol + 1
                    value# = 0!
                    FOR j = 1 TO NSpecies
                      value# = value# + HABound(j) * SCb(j, s%)
                    NEXT j
                    stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
                ELSE  ' if not complete
                    iCol = iCol + 1
                    stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
                    iCol = iCol + 1
                    stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
                END IF
              END IF
              IF iComplete AND NOT ErrorNotConverged THEN
                  ' Combined FA and HA as Total Organic
                  iCol = iCol + 1
                  value# = 0!
                  FOR j = 1 TO NSpecies
                    value# = value# + wOrgBound!(j) * SCb(j, s%)
                  NEXT j
                  stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
              ELSE  ' if not complete
                  iCol = iCol + 1
                  stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
              END IF
          END IF
       NEXT s%
   END IF

   FOR p% = 1 TO NPhase
      iCol = iCol + 1
      IF iComplete AND NOT ErrorNotConverged THEN
         value# = ChangeMass(p%)
         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      ELSE  ' if not complete
         stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      END IF
   NEXT p%

   ' Write charge balance to output file
   iCol = iCol + 1
   IF iComplete AND NOT ErrorNotConverged THEN
      value# = OC(NSpecies + 1)
      stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   ELSE  ' if not complete
      stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   END IF


   ' Write ionic strength to output file
   iCol = iCol + 1
   IF iComplete AND NOT ErrorNotConverged THEN
      value# = OC(NSpecies + 2)
      stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   ELSE  ' if not complete
      stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   END IF

   ' Write activity corrections (added to ver 2.45, build 2021-04-05)
   FOR s% = 1 TO 4
       IF iComplete AND NOT ErrorNotConverged THEN
         iCol = iCol + 1
         value# = wGAMMA(s%)
         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
       ELSE  ' if not complete
          stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
       END IF
   NEXT s%

   ' Write charge for Humics
   iCol = iCol + 1
   value# = ZED(1)
   stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write charge for fulvics
   iCol = iCol + 1
   value# = ZED(2)
   stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write diffuse layer vol for humics
   iCol = iCol + 1
   value# = DDLVOL(1)
   stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write diffuse layer vol for humics
   iCol = iCol + 1
   value# = DDLVOL(2)
   stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write max diffuse layer vol for humics
   iCol = iCol + 1
   value# = DVOLMAX(1)
   stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write max diffuse layer vol for fulvics
   iCol = iCol + 1
   value# = DVOLMAX(2)
   stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write mass humics
   iCol = iCol + 1
   value# = wTH(1)
   stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write mass fulvics
   iCol = iCol + 1
   value# = wTH(2)
   stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write Ratio humics
   iCol = iCol + 1
   value# = wRATIO(1)
   stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write Ratio fulvics
   iCol = iCol + 1
   value# = wRATIO(2)
   stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' write DVol
   iCol = iCol + 1
   value# = DVol(1)
   stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   iCol = iCol + 1
   value# = DVol(2)
   stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' write individual species concentrations for organically bound substances
    FOR i = 1 TO NSP%
       ii = CMap%(i)
       IF (ii <> 0) THEN
           ' humic bound specific
           iCol = iCol + 1
           value# = wCHC(1, i)
           stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

           ' humic bound diffuse
           iCol = iCol + 1
           value# = wCDDL(1, i)
           stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

           ' total humic bound
           iCol = iCol + 1
           value# = HABound(ii)
           stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

           ' fulvic bound specific
           iCol = iCol + 1
           value# = wCHC(2, i)
           stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

           ' fulvic bound diffuse
           iCol = iCol + 1
           value# = wCDDL(2, i)
           stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

          ' total fulvic bound
           iCol = iCol + 1
           value# = FABound(ii)
           stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
       END IF
    NEXT i


   ' Write temperature to output file
   iCol = iCol + 1
   value# = SysTemp
   stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write solution vol (total vol - ddl vol)
   iCol = iCol + 1
   value# = VOLSOL
   stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   ' Write BL mass to output file
   FOR i = 1 TO NMass
      ' Write the mass of each mass compartment
      iCol = iCol + 1
      value# = MassVal(i)
      stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   NEXT i

   ' Write the # iterations to output file
   iCol = iCol + 1
   value# = OC(0)
   stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   FOR c% = 1 TO NComp
      IF (ActCorr(c%) <> acWHAM) THEN
         iCol = iCol + 1
         IF iComplete THEN
            value# = OKT(c%)
            '
            ' 2014-07-07 - RCS: added IF block so that reported totals will be consistent with
            '              the sum of species for fixed concentration components
            '
            IF CType(c%)= ctFixed THEN
                value# = KT(c%)
            END IF

            stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
         ELSE
            stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
         END IF

      END IF
   NEXT c%

   FOR p% = 1 TO NPhase
      iCol = iCol + 1
      value# = PhaseMass(p%)
      stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   NEXT p%

IF FALSE THEN
   FOR c% = 1 TO NSpecies
'      PRINT #iHand, USING OutputMask$; ActCoef(c);
   NEXT c%

   FOR s = NFirstS TO NSpecies
      Temp# = CDBL(hEC(s))
'      PRINT #iHand, USING OutputMask$; Temp#;
   NEXT s

   FOR c% = 1 TO NComp
'      PRINT #iHand, USING OutputMask$; R(c%);
   NEXT c%

   FOR c% = 1 TO NPhase
      Temp# = IAP#(c%)
'      PRINT #iHand, USING OutputMask$; Temp#;
   NEXT c%

   FOR c% = 1 TO NPhase
      Temp# = LOG10S!(CSNG(IAP#(c%) / SEC(c%)))
'      PRINT #iHand, USING OutputMask$; Temp#;
   NEXT c%
END IF
'PRINT #iHand,

END SUB


'SUB WriteLabel (iHand, LabelMask$, L$)
'   'PRINT #iHand, USING LabelMask$; L$
'   PRINT #iHand, FORMAT$(L$, LabelMask$)
'END SUB

SUB WriteSimple (iHand, Iter, LC50)

LabelMask$ = " \        \"
OutputMask$ = " ##.###^^^^"

'
' On the first line of output, write column labels
IF (Iter = 1) THEN
   'PRINT #iHand, USING "\                   \"; "Data Label";
   'PRINT #iHand, USING "\                   \"; "Site Label";
   PRINT #iHand, FormatStr$("Site Label", "\                   \");
   PRINT #iHand, FormatStr$("Data Label", "\                   \");

   'PRINT #iHand, USING LabelMask$; "Mode";
   'PRINT #iHand, USING LabelMask$; "pH";
   PRINT #iHand, FormatStr$("Mode", LabelMask$);
   PRINT #iHand, FormatStr$("pH", LabelMask$);

   L$ = "Dis. " + Label(iMetal)
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);
   L$ = "Free " + Label(iMetal)
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);

   L$ = "Actv " + Label(iMetal)
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);

   L$ = "Org " + Label(iMetal)
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);

   L$ = "TOrg " + Label(iMetal)
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);

   FOR j = 1 TO NBLM
      L$ = Label(iBLM(j))
      'PRINT #iHand, USING LabelMask$; L$;
      PRINT #iHand, FormatStr$(L$, LabelMask$);
   NEXT j
   L$ = "DOC"
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);
   L$ = "HA%"
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);
   FOR c% = 1 TO 12
      IF ((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL)) THEN
         'PRINT #iHand, USING LabelMask$; Label(c%);
         PRINT #iHand, FormatStr$(Label(c%), LabelMask$);
      END IF
   NEXT c%
   PRINT #iHand,

   '
   ' Second line, print units
   '
   'PRINT #iHand, USING "\                   \"; "          ";
   'PRINT #iHand, USING "\                   \"; "          ";
   'PRINT #iHand, USING LabelMask$; "          ";
   'PRINT #iHand, USING LabelMask$; "          ";
   'PRINT #iHand, USING LabelMask$; "mol/L     ";
   'PRINT #iHand, USING LabelMask$; "mol/L     ";
   'PRINT #iHand, USING LabelMask$; "mol/L     ";
   'PRINT #iHand, USING LabelMask$; "mol/L     ";
   'PRINT #iHand, USING LabelMask$; "mol/L     ";
   PRINT #iHand, FormatStr$(" ", "\                   \");
   PRINT #iHand, FormatStr$(" ", "\                   \");
   PRINT #iHand, FormatStr$(" ", LabelMask$);
   PRINT #iHand, FormatStr$(" ", LabelMask$);
   PRINT #iHand, FormatStr$("mol/L", LabelMask$);
   PRINT #iHand, FormatStr$("mol/L", LabelMask$);
   PRINT #iHand, FormatStr$("mol/L", LabelMask$);
   PRINT #iHand, FormatStr$("mol/L", LabelMask$);
   PRINT #iHand, FormatStr$("mol/L", LabelMask$);

   FOR j = 1 TO NBLM
      'PRINT #iHand, USING LabelMask$; "nmol/g wet";
      PRINT #iHand, FormatStr$("nmol/g wet", LabelMask$);
   NEXT j
'  PRINT #iHand, USING LabelMask$; "mg/L      ";
'  PRINT #iHand, USING LabelMask$; "          ";
   PRINT #iHand, FormatStr$("mg/L", LabelMask$);
   PRINT #iHand, FormatStr$(" ", LabelMask$);
  FOR c% = 1 TO 12
      IF ((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL)) THEN
         'PRINT #iHand, USING LabelMask$; "mol/L";
         PRINT #iHand, FormatStr$("mol/L", LabelMask$);
      END IF
   NEXT c%
   PRINT #iHand,

END IF


'
'
' Write the actual data
'
   'PRINT #iHand, USING "\                   \"; DataL$;
   'PRINT #iHand, USING "\                   \"; SiteL$;
   PRINT #iHand, FormatStr$(DataL$, "\                   \");
   PRINT #iHand, FormatStr$(SiteL$, "\                   \");

   IF LC50 THEN
      L$ = "LC50"
   ELSE
      L$ = "Dissolved"
   END IF
   'PRINT #iHand, USING LabelMask$; L$;
   PRINT #iHand, FormatStr$(L$, LabelMask$);

   IF NOT (iWHAM) THEN
      PH! = -LOG10S!(CSNG(OC(1)))
   END IF
'   PRINT #iHand, USING "##.####   "; PH!;
   PRINT #iHand, FORMAT$(PH!, "##.####");

   Tmp! = OKT(iMetal)
'   PRINT #iHand, USING OutputMask$; Tmp!;
   PRINT #iHand, FORMAT$(Tmp!, OutputMask$);

   Tmp! = OC(iMetal)
'   PRINT #iHand, USING OutputMask$; Tmp!;
   PRINT #iHand, FORMAT$(Tmp!, OutputMask$);

   Tmp! = OC(iMetal) * ActCoef(iMetal)
'   PRINT #iHand, USING OutputMask$; Tmp!;
   PRINT #iHand, FORMAT$(Tmp!, OutputMask$);

   Tmp! = wOrgBound!(iMetal)
'   PRINT #iHand, USING OutputMask$; Tmp!;
   PRINT #iHand, FORMAT$(Tmp!, OutputMask$);

   Tmp! = 0!
   FOR j = 1 TO NSpecies
      Tmp! = Tmp! + wOrgBound!(j) * SCb(j, iMetal)
   NEXT j
'   PRINT #iHand, USING OutputMask$; Tmp!;
   PRINT #iHand, FORMAT$(Tmp!, OutputMask$);

   FOR j = 1 TO NBLM
      Tmp! = OC(iBLM(j)) * 1000000! / KT(iBL)
'   PRINT #iHand, USING OutputMask$; Tmp!;
      PRINT #iHand, FORMAT$(Tmp!, OutputMask$);
   NEXT j

   Tmp! = OKT(iDOC) * 1000! / 2
   'PRINT #iHand, USING OutputMask$; Tmp!;
   PRINT #iHand, FORMAT$(Tmp!, OutputMask$);

   Tmp! = PercHA! * 100!
   'PRINT #iHand, USING OutputMask$; Tmp!;
   PRINT #iHand, FORMAT$(Tmp!, OutputMask$);

   FOR c% = 1 TO 12
      IF ((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL)) THEN
         Tmp! = OKT(c%)
         'PRINT #iHand, USING OutputMask$; Tmp!;
         PRINT #iHand, FORMAT$(Tmp!, OutputMask$);
      END IF
   NEXT c%
   PRINT #iHand,

END SUB

SUB WriteSimple2 (iHand&, Iter, LC50, iComplete)

LabelMask$ = " \        \"
OutputMask$ = " ##.###^^^^"

LabelMask$ =  "\               \"
OutputMask$ = "       ##.###^^^^"

LabelMask$ =  "\             \"
OutputMask$ = "     ##.###^^^^"

'
' On the first line of output, write column labels
IF (Iter = 1) THEN
   PRINT #iHand&, FormatStr$("Site Label", SPACE$(22));
   PRINT #iHand&, FormatStr$("Sample Label", SPACE$(22));

   PRINT #iHand&, FormatStr$("Mode", LabelMask$);
   PRINT #iHand&, FormatStr$("pH", LabelMask$);

   L$ = "Dis. " + Label(iMetal)
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   IF iLogKsAl THEN
      L$ = "Dis. " + Label(iMetal) + " Eff"
      PRINT #iHand&, FormatStr$(L$, LabelMask$);

      L$ = "Part " + Label(iMetal)
      PRINT #iHand&, FormatStr$(L$, LabelMask$);

      L$ = "Part " + Label(iMetal) + " Eff"
      PRINT #iHand&, FormatStr$(L$, LabelMask$);
   END IF

   L$ = "Free " + Label(iMetal)
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   L$ = "TOrg " + Label(iMetal)
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   IF (SepNOM = TRUE) THEN
     L$ = "HA " + Label(iMetal)
     PRINT #iHand&, FormatStr$(L$, LabelMask$);

     L$ = "FA " + Label(iMetal)
     PRINT #iHand&, FormatStr$(L$, LabelMask$);
   END IF

   FOR j = 1 TO NBLM
      L$ = Label(iBLM(j))
      PRINT #iHand&, FormatStr$(L$, LabelMask$);
   NEXT j

   L$ = "DOC"
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   L$ = "HA%"
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   FOR c% = 1 TO iBL - 1
      IF ((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL)) THEN
         L$ = "T." + Label(c%)
         PRINT #iHand&, FormatStr$(L$, LabelMask$);
      END IF
   NEXT c%
   PRINT #iHand&,

   '
   ' Second line, print units
   '
   PRINT #iHand&, FormatStr$(" ", SPACE$(22));
   PRINT #iHand&, FormatStr$(" ", SPACE$(22));
   PRINT #iHand&, FormatStr$(" ", LabelMask$);
   PRINT #iHand&, FormatStr$(" ", LabelMask$);
   PRINT #iHand&, FormatStr$("mol/L", LabelMask$);

   IF iLogKsAl THEN
      PRINT #iHand&, FormatStr$("frac cont", LabelMask$);
      PRINT #iHand&, FormatStr$("mol/L", LabelMask$);
      PRINT #iHand&, FormatStr$("frac cont", LabelMask$);
   END IF

   PRINT #iHand&, FormatStr$("mol/L", LabelMask$);
   PRINT #iHand&, FormatStr$("mol/L", LabelMask$);

   FOR j = 1 TO NBLM
      PRINT #iHand&, FormatStr$("nmol/gw", LabelMask$);
   NEXT j
   PRINT #iHand&, FormatStr$("mg/L", LabelMask$);
   PRINT #iHand&, FormatStr$(" ", LabelMask$);
   FOR c% = 1 TO iBL - 1
      IF ((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL)) THEN
         PRINT #iHand&, FormatStr$("mol/L", LabelMask$);
      END IF
   NEXT c%
   IF (SepNOM = TRUE) THEN
         PRINT #iHand&, FormatStr$("mol/L", LabelMask$);
         PRINT #iHand&, FormatStr$("mol/L", LabelMask$);
   END IF
   PRINT #iHand&,

END IF


'
'
' Write the actual data
'
   'Site Label
   PRINT #iHand&, FormatStr$(SiteL$, SPACE$(22));

   'Sample Label
   PRINT #iHand&, FormatStr$(DataL$, SPACE$(22));

   'Mode Flag
   IF LC50 THEN
      L$ = "Toxicity"
   ELSE
      L$ = "Speciation"
   END IF
   IF ErrorNotConverged THEN
      L$ = "Not Converged"
   END IF
   IF NOT iComplete THEN
      L$ = "Incomplete"
   END IF
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   'pH
   IF NOT (iWHAM) THEN
      PH! = -LOG10S!(CSNG(OC(1)))
   END IF
   IF NOT (iComplete) THEN
      PH! = -LOG10S!(CSNG(KT(1)))
   END IF
   L$ = "##.####"
   j = LEN(OutputMask$) - LEN(L$)
   IF (j<0) THEN j=0
   L$ = SPACE$(j) + L$
   PRINT #iHand&, USING$(L$, PH!);

   'Diss. Metal
   IF iComplete AND NOT ErrorNotConverged THEN
      Tmp! = OKT(iMetal)
      IF CType(iMetal)= ctFixed THEN
          Tmp! = KT(iMetal)
      END IF
      PRINT #iHand&, USING$(OutputMask$, Tmp!);
   ELSE
      PRINT #iHand&, "             NA";
   END IF

   'Al Dissolved/Precipitate information
   IF iLogKsAl THEN
      IF iComplete AND NOT ErrorNotConverged THEN
         Tmp! = DISS_Effect
         PRINT #iHand&, USING$(OutputMask$, Tmp!);

         Tmp! = PART_Conc /(1000000*26.981)
         PRINT #iHand&, USING$(OutputMask$, Tmp!);

         Tmp! = PART_Effect
         PRINT #iHand&, USING$(OutputMask$, Tmp!);
      ELSE
         PRINT #iHand&, "             NA";
         PRINT #iHand&, "             NA";
         PRINT #iHand&, "             NA";
      END IF
   END IF

   'Free Metal
   IF iComplete AND NOT ErrorNotConverged THEN
      Tmp! = OC(iMetal)
      PRINT #iHand&, USING$(OutputMask$, Tmp!);
   ELSE
      PRINT #iHand&, "             NA";
   END IF

   IF iComplete AND NOT ErrorNotConverged THEN
      Tmp! = 0!
      FOR j = 1 TO NSpecies
         Tmp! = Tmp! + wOrgBound!(j) * SCb(j, iMetal)
      NEXT j
      PRINT #iHand&, USING$(OutputMask$, Tmp!);
   ELSE
      PRINT #iHand&, "             NA";
   END IF

   'DOC fractions
   IF (SepNOM = TRUE) THEN
     ' HA Fraction
     IF NOT iComplete THEN
        Tmp! = -999.#
        PRINT #iHand&, "             NA";
     ELSE
        Tmp! = 0!
        FOR j = 1 TO NSpecies
          Tmp! = Tmp! + HABound(j) * SCb(j, iMetal)
        NEXT j
        PRINT #iHand&, USING$(OutputMask$, Tmp!);
     END IF


     ' FA Fraction
     IF NOT iComplete THEN
        Tmp! = -999.#
        PRINT #iHand&, "             NA";
     ELSE
        Tmp! = 0!
        FOR j = 1 TO NSpecies
          Tmp! = Tmp! + FABound(j) * SCb(j, iMetal)
        NEXT j
        PRINT #iHand&, USING$(OutputMask$, Tmp!);
     END IF
   END IF

   FOR j = 1 TO NBLM
      IF iComplete AND NOT ErrorNotConverged THEN
         Tmp! = OC(iBLM(j)) * 1000000! / KT(iBL)
         PRINT #iHand&, USING$(OutputMask$, Tmp!);
      ELSE
         PRINT #iHand&, "             NA";
      END IF
   NEXT j

   Tmp! = OKT(iDOC) * 1000! / 2
   PRINT #iHand&, USING$(OutputMask$, Tmp!);

   Tmp! = PercHA! * 100!
   PRINT #iHand&, USING$(OutputMask$, Tmp!);

   FOR c% = 1 TO iBL - 1
      IF ((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL)) THEN
         Tmp! = OKT(c%)
         IF CType(c%)= ctFixed THEN
            Tmp! = KT(c%)
         END IF
         PRINT #iHand&, USING$(OutputMask$, Tmp!);
      END IF
   NEXT c%
   PRINT #iHand&,

END SUB


SUB WriteSimpleold (iHand&, Iter, LC50, iComplete)

LabelMask$ = " \        \"
OutputMask$ = " ##.###^^^^"

'x& = FILEATTR(#iHand&, 0)
'L$ = str$(x&)
'Call PMsg(L$)

'
' On the first line of output, write column labels
IF (Iter = 1) THEN
   'PRINT #iHand&, USING "\                   \"; "Data Label";
   'PRINT #iHand&, USING "\                   \"; "Site Label";
   PRINT #iHand&, FormatStr$("Site Label", "\                   \");
   PRINT #iHand&, FormatStr$("Data Label", "\                   \");

   'PRINT #iHand&, USING LabelMask$; "Mode";
   'PRINT #iHand&, USING LabelMask$; "pH";
   PRINT #iHand&, FormatStr$("Mode", LabelMask$);
   PRINT #iHand&, FormatStr$("pH", LabelMask$);

   L$ = "Dis. " + Label(iMetal)
   'PRINT #iHand&, USING LabelMask$; L$;
   PRINT #iHand&, FormatStr$(L$, LabelMask$);
   L$ = "Free " + Label(iMetal)
   'PRINT #iHand&, USING LabelMask$; L$;
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

'   L$ = "Actv " + Label(iMetal)
'   'PRINT #iHand&, USING LabelMask$; L$;
'   PRINT #iHand&, FormatStr$(L$, LabelMask$);
'
'   L$ = "Org " + Label(iMetal)
'   'PRINT #iHand&, USING LabelMask$; L$;
'   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   L$ = "TOrg " + Label(iMetal)
   'PRINT #iHand&, USING LabelMask$; L$;
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   FOR j = 1 TO NBLM
      L$ = Label(iBLM(j))
      'PRINT #iHand&, USING LabelMask$; L$;
      PRINT #iHand&, FormatStr$(L$, LabelMask$);
   NEXT j
   L$ = "DOC"
   'PRINT #iHand&, USING LabelMask$; L$;
   PRINT #iHand&, FormatStr$(L$, LabelMask$);
   L$ = "HA%"
   'PRINT #iHand&, USING LabelMask$; L$;
   PRINT #iHand&, FormatStr$(L$, LabelMask$);
   FOR c% = 1 TO 12
      IF ((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL)) THEN
         'PRINT #iHand&, USING LabelMask$; Label(c%);
         L$ = "T." + Label(c%)
         'IF (L$ = "CO3") THEN
         '   L$ = "DIC"
         'END IF
         PRINT #iHand&, FormatStr$(L$, LabelMask$);
      END IF
   NEXT c%
   PRINT #iHand&,

   '
   ' Second line, print units
   '
   'PRINT #iHand&, USING "\                   \"; "          ";
   'PRINT #iHand&, USING "\                   \"; "          ";
   'PRINT #iHand&, USING LabelMask$; "          ";
   'PRINT #iHand&, USING LabelMask$; "          ";
   'PRINT #iHand&, USING LabelMask$; "mol/L     ";
   'PRINT #iHand&, USING LabelMask$; "mol/L     ";
   'PRINT #iHand&, USING LabelMask$; "mol/L     ";
   'PRINT #iHand&, USING LabelMask$; "mol/L     ";
   'PRINT #iHand&, USING LabelMask$; "mol/L     ";
   PRINT #iHand&, FormatStr$(" ", "\                   \");
   PRINT #iHand&, FormatStr$(" ", "\                   \");
   PRINT #iHand&, FormatStr$(" ", LabelMask$);
   PRINT #iHand&, FormatStr$(" ", LabelMask$);
   PRINT #iHand&, FormatStr$("mol/L", LabelMask$);
   PRINT #iHand&, FormatStr$("mol/L", LabelMask$);
'   PRINT #iHand&, FormatStr$("mol/L", LabelMask$);
'   PRINT #iHand&, FormatStr$("mol/L", LabelMask$);
   PRINT #iHand&, FormatStr$("mol/L", LabelMask$);

   FOR j = 1 TO NBLM
      'PRINT #iHand&, USING LabelMask$; "nmol/g wet";
      PRINT #iHand&, FormatStr$("nmol/gw", LabelMask$);
   NEXT j
'  PRINT #iHand&, USING LabelMask$; "mg/L      ";
'  PRINT #iHand&, USING LabelMask$; "          ";
   PRINT #iHand&, FormatStr$("mg/L", LabelMask$);
   PRINT #iHand&, FormatStr$(" ", LabelMask$);
  FOR c% = 1 TO 12
      IF ((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL)) THEN
         'PRINT #iHand&, USING LabelMask$; "mol/L";
         PRINT #iHand&, FormatStr$("mol/L", LabelMask$);
      END IF
   NEXT c%
   PRINT #iHand&,

END IF


'
'
' Write the actual data
'
   'PRINT #iHand&, USING "\                   \"; DataL$;
   'PRINT #iHand&, USING "\                   \"; SiteL$;
   PRINT #iHand&, FormatStr$(DataL$, "\                   \");
   PRINT #iHand&, FormatStr$(SiteL$, "\                   \");

   IF LC50 THEN
      L$ = "LC50"
      L$ = "Toxicity"
   ELSE
      L$ = "Dissolved"
      L$ = "Speciation"
   END IF
   IF NOT iComplete THEN
      L$ = "Incomplete"
   END IF
   'PRINT #iHand&, USING LabelMask$; L$;
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   IF NOT (iWHAM) THEN
      PH! = -LOG10S!(CSNG(OC(1)))
   END IF
   IF NOT (iComplete) THEN
      PH! = -LOG10S!(CSNG(KT(1)))
   END IF
'   PRINT #iHand&, USING "##.####   "; PH!;
   PRINT #iHand&, FORMAT$(PH!, "##.####");

   Tmp! = OKT(iMetal)
'   PRINT #iHand&, USING OutputMask$; Tmp!;
   PRINT #iHand&, FORMAT$(Tmp!, OutputMask$);

   Tmp! = OC(iMetal)
'   PRINT #iHand&, USING OutputMask$; Tmp!;
   IF NOT iComplete THEN
      Tmp! = -999
   END IF
   PRINT #iHand&, FORMAT$(Tmp!, OutputMask$);

'   Tmp! = OC(iMetal) * ActCoef(iMetal)
'   PRINT #iHand&, USING OutputMask$; Tmp!;
'   IF NOT iComplete THEN
'      Tmp! = -999
'   END IF
'   PRINT #iHand&, FORMAT$(Tmp!, OutputMask$);
'
'   Tmp! = wOrgBound!(iMetal)
''   PRINT #iHand&, USING OutputMask$; Tmp!;
'   IF NOT iComplete THEN
'      Tmp! = -999
'   END IF
'   PRINT #iHand&, FORMAT$(Tmp!, OutputMask$);

   Tmp! = 0!
   FOR j = 1 TO NSpecies
      Tmp! = Tmp! + wOrgBound!(j) * SCb(j, iMetal)
   NEXT j
'   PRINT #iHand&, USING OutputMask$; Tmp!;
   IF NOT iComplete THEN
      Tmp! = -999
   END IF
   PRINT #iHand&, FORMAT$(Tmp!, OutputMask$);

   FOR j = 1 TO NBLM
      Tmp! = OC(iBLM(j)) * 1000000! / KT(iBL)
'   PRINT #iHand&, USING OutputMask$; Tmp!;
      IF NOT iComplete THEN
         Tmp! = -999
      END IF
      PRINT #iHand&, FORMAT$(Tmp!, OutputMask$);
   NEXT j

   Tmp! = OKT(iDOC) * 1000! / 2
   'PRINT #iHand&, USING OutputMask$; Tmp!;
   PRINT #iHand&, FORMAT$(Tmp!, OutputMask$);

   Tmp! = PercHA! * 100!
   'PRINT #iHand&, USING OutputMask$; Tmp!;
   PRINT #iHand&, FORMAT$(Tmp!, OutputMask$);

   FOR c% = 1 TO 12
      IF ((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL)) THEN
         Tmp! = OKT(c%)
         'PRINT #iHand&, USING OutputMask$; Tmp!;
         PRINT #iHand&, FORMAT$(Tmp!, OutputMask$);
      END IF
   NEXT c%
   PRINT #iHand&,
END SUB



SUB WriteSimpleXLS (iHand&, Iter, LC50, iComplete)

LabelMask$ = " \        \"
OutputMask$ = " ##.###^^^^"
DIM iRow AS LONG
DIM iRowOS AS LONG

iRowOS = 6
'
' On the first line of output, write column labels
IF (Iter = 1) THEN
   iRow = 5

   iCol = 1
   sTmp1$ = "Site Label"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   iCol = 2
   sTmp1$ = "Sample Label"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   iCol = 3
   sTmp1$ = "Mode"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   iCol = 4
   sTmp1$ = "pH"
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   iCol = 5
   sTmp1$ = "Dis. " + Label(iMetal)
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   IF iLogKsAl THEN
      iCol = 6
      sTmp1$ = "Dis. " + Label(iMetal) + " Eff"
      stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      iCol = 7
      sTmp1$ = "Part " + Label(iMetal)
      stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      iCol = 8
      sTmp1$ = "Part " + Label(iMetal) + " Eff"
      stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   END IF

   sTmp1$ = "Free " + Label(iMetal)
   iCol = iCol + 1
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   L$ = "TOrg " + Label(iMetal)
   sTmp1$ = L$
   iCol = iCol + 1
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   IF (SepNOM = TRUE) THEN
     L$ = "HA " + Label(iMetal)
     sTmp1$ = L$
     iCol = iCol + 1
     stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

     L$ = "FA " + Label(iMetal)
     sTmp1$ = L$
     iCol = iCol + 1
     stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   END IF 'SepNOM = TRUE

   FOR j = 1 TO NBLM
      sTmp1$ = Label(iBLM(j))
      iCol = iCol + 1
      stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   NEXT j '= 1 TO NBLM

   L$ = "DOC"
   sTmp1$ = L$
   iCol = iCol + 1
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   L$ = "HA%"
   sTmp1$ = L$
   iCol = iCol + 1
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   FOR c% = 1 TO 12
      IF ((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL)) THEN
         L$ = "T." + Label(c%)
         sTmp1$ = L$
         iCol = iCol + 1
         stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      END IF '((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL))
   NEXT c% '= 1 TO 12

   '
   ' Second line, print units
   '
   iRow = 6
   iCol = 5  'the first few columns have no units, so start at 5 (Diss. Me)
   stat& = xlsWriteText("mol/L", iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   IF iLogKsAL THEN
      iCol = iCol + 1 '=6
      stat& = xlsWriteText("frac cont", iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

      iCol = iCol + 1 '=7
      stat& = xlsWriteText("mol/L", iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

      iCol = iCol + 1 '=8
      stat& = xlsWriteText("frac cont", iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   END IF 'iLogKsAl

   iCol = iCol + 1    '=6 or 9
   stat& = xlsWriteText("mol/L", iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   iCol = iCol + 1
   stat& = xlsWriteText("mol/L", iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   IF (SepNOM = TRUE) THEN
     iCol = iCol + 1
     stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

     iCol = iCol + 1
     stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   END IF 'SepNOM = TRUE

   sTmp1$ = "nmol/gw"
   FOR j = 1 TO NBLM
      iCol = iCol + 1
      stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   NEXT j '= 1 TO NBLM

   sTmp1$ = "mg/L"
   iCol = iCol + 1
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   sTmp1$ = ""
   iCol = iCol + 1
   stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   FOR c% = 1 TO 12
      IF ((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL)) THEN
         sTmp1$ = "mol/L"
         iCol = iCol + 1
         stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      END IF '((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL))
   NEXT c% '= 1 TO 12

END IF

'
' Write the actual data
'
   'Site Label
   iCol = 1
   iRow = iRowOS + Iter
   stat& = xlsWriteText(SiteL$, iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   'Sample Label
   iCol = iCol + 1
   stat& = xlsWriteText(DataL$, iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   'Mode Flag
   iCol = iCol + 1
   IF LC50 THEN
      L$ = "Toxicity"
   ELSE
      L$ = "Speciation"
   END IF
   IF ErrorNotConverged THEN
      L$ = "Not converged"
      'ErrorNotConverged = FALSE  'moved to right before the loop starts
   END IF
   IF NOT iComplete THEN
      L$ = "Incomplete"
   END IF
   stat& = xlsWriteText(L$, iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   'pH
   iCol = iCol + 1
   IF NOT (iWHAM) THEN
      PH! = -LOG10S!(CSNG(OC(1)))
   END IF
   IF NOT (iComplete) THEN
      PH! = -LOG10S!(CSNG(KT(1)))
   END IF
   value# = pH!
   stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

   'Diss. Metal
   iCol = iCol + 1
   IF iComplete AND NOT ErrorNotConverged THEN
      value# = OKT(iMetal)
      IF CType(iMetal)= ctFixed THEN
        value# = KT(iMetal)
      END IF
      stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   ELSE
      stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   END IF

   'Aluminum dissolved/precipitate info
   IF iLogKsAl THEN
      iCol = iCol + 1
      IF iComplete AND NOT ErrorNotConverged THEN
         value# = DISS_Effect
         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      ELSE
         stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      END IF

      iCol = iCol + 1
      IF iComplete AND NOT ErrorNotConverged THEN
         value# = PART_Conc /(1000000*26.981)
         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      ELSE
         stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      END IF

      iCol = iCol + 1
      IF iComplete AND NOT ErrorNotConverged THEN
         value# = PART_Effect
         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      ELSE
         stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      END IF
   END IF

   'Free Metal
   iCol = iCol + 1
   IF iComplete AND NOT ErrorNotConverged THEN
      value# = OC(iMetal)
      stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   ELSE
      stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   END IF

   'Organically-bound metal
   iCol = iCol + 1
   IF iComplete AND NOT ErrorNotConverged THEN
      value# = 0!
      FOR j = 1 TO NSpecies
         value# = value# + wOrgBound!(j) * SCb(j, iMetal)
      NEXT j
      stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   ELSE
      stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   END IF

   'Fraction bound to separate NOM fractions
   IF (SepNOM = TRUE) THEN
     ' HA Fraction
     iCol = iCol + 1
     IF iComplete AND NOT ErrorNotConverged THEN
        value# = 0!
        FOR j = 1 TO NSpecies
          value# = value# + HABound(j) * SCb(j, iMetal)
        NEXT j
        stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
     ELSE
        stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
     END IF


     ' FA Fraction
     iCol = iCol + 1
     IF iComplete AND NOT ErrorNotConverged THEN
        value# = 0!
        FOR j = 1 TO NSpecies
          value# = value# + FABound(j) * SCb(j, iMetal)
        NEXT j
        stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
     ELSE
        stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
     END IF

   END IF


   'Biotic ligand species
   FOR j = 1 TO NBLM
      iCol = iCol + 1
      IF iComplete AND NOT ErrorNotConverged THEN
         value# = OC(iBLM(j)) * 1000000! / KT(iBL)
         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      ELSE
         stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      END IF
   NEXT j

   'DOC
   iCol = iCol + 1
   IF iComplete AND NOT ErrorNotConverged THEN
      IF iWHAM THEN
          value# = OKT(iDOC) * 1000! / 2
      ELSEIF (iDOC > 0 ) THEN
          Value# = OKT(iDOC)
      END IF
      stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   ELSE
      stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   END IF

   'HA%
   iCol = iCol + 1
   IF iComplete AND NOT ErrorNotConverged THEN
      value# = PercHA! * 100!
      stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   ELSE
      stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
   END IF

   'other chemical species
   FOR c% = 1 TO 12
      IF ((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL)) THEN
         iCol = iCol + 1
         value# = OKT(c%)
         IF CType(c%)= ctFixed THEN
             value# = KT(c%)
         END IF
         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
      END IF
   NEXT c%

END SUB

SUB WriteSimpleComma (iHand&, Iter, LC50, iComplete)

LabelMask$ = ""
'OutputMask$ = " ##.###^^^^,"
OutputMask1$ = " ##.#####^^^^_,"
OutputMask2$ = " ##.#####_,"

'
' On the first line of output, write column labels
IF (Iter = 1) THEN
   PRINT #iHand&, FormatStr$("Site Label", LabelMask$);
   PRINT #iHand&, FormatStr$("Sample Label", LabelMask$);

   PRINT #iHand&, FormatStr$("Mode", LabelMask$);
   PRINT #iHand&, FormatStr$("pH", LabelMask$);

   L$ = "Dis. " + Label(iMetal)
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   IF iLogKsAl THEN
      L$ = "Dis. " + Label(iMetal) + " Eff"
      PRINT #iHand&, FormatStr$(L$, LabelMask$);

      L$ = "Part " + Label(iMetal)
      PRINT #iHand&, FormatStr$(L$, LabelMask$);

      L$ = "Part " + Label(iMetal) + " Eff"
      PRINT #iHand&, FormatStr$(L$, LabelMask$);
   END IF

   L$ = "Free " + Label(iMetal)
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   L$ = "TOrg " + Label(iMetal)
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   IF (SepNOM = TRUE) THEN
     L$ = "HA " + Label(iMetal)
     PRINT #iHand&, FormatStr$(L$, LabelMask$);

     L$ = "FA " + Label(iMetal)
     PRINT #iHand&, FormatStr$(L$, LabelMask$);
   END IF

   FOR j = 1 TO NBLM
      L$ = Label(iBLM(j))
      PRINT #iHand&, FormatStr$(L$, LabelMask$);
   NEXT j

   L$ = "DOC"
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   L$ = "HA%"
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   FOR c% = 1 TO iBL - 1
      IF ((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL)) THEN
         L$ = "T." + Label(c%)
         PRINT #iHand&, FormatStr$(L$, LabelMask$);
      END IF
   NEXT c%
   PRINT #iHand&,

   '
   ' Second line, print units
   '
   PRINT #iHand&, FormatStr$(" ", LabelMask$);
   PRINT #iHand&, FormatStr$(" ", LabelMask$);
   PRINT #iHand&, FormatStr$(" ", LabelMask$);
   PRINT #iHand&,  FormatStr$(" ", LabelMask$);
   PRINT #iHand&, FormatStr$("mol/L", LabelMask$);

   IF iLogKsAl THEN
      PRINT #iHand&, FormatStr$("frac cont", LabelMask$);
      PRINT #iHand&, FormatStr$("mol/L", LabelMask$);
      PRINT #iHand&, FormatStr$("frac cont", LabelMask$);
   END IF

   PRINT #iHand&, FormatStr$("mol/L", LabelMask$);
   PRINT #iHand&, FormatStr$("mol/L", LabelMask$);

   FOR j = 1 TO NBLM
      PRINT #iHand&, FormatStr$("nmol/gw", LabelMask$);
   NEXT j
   PRINT #iHand&, FormatStr$("mg/L", LabelMask$);
   PRINT #iHand&, FormatStr$(" ", LabelMask$);
   FOR c% = 1 TO iBL - 1
      IF ((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL)) THEN
         PRINT #iHand&, FormatStr$("mol/L", LabelMask$);
      END IF
   NEXT c%
   IF (SepNOM = TRUE) THEN
         PRINT #iHand&, FormatStr$("mol/L", LabelMask$);
         PRINT #iHand&, FormatStr$("mol/L", LabelMask$);
   END IF
   PRINT #iHand&,

END IF


'
'
' Write the actual data
'
   'Site Label
   PRINT #iHand&, FormatStr$(SiteL$, LabelMask$);

   'Sample Label
   PRINT #iHand&, FormatStr$(DataL$, LabelMask$);

   'Mode Flag
   IF LC50 THEN
      L$ = "Toxicity"
   ELSE
      L$ = "Speciation"
   END IF
   IF ErrorNotConverged THEN
      L$ = "Not Converged"
   END IF
   IF NOT iComplete THEN
      L$ = "Incomplete"
   END IF
   PRINT #iHand&, FormatStr$(L$, LabelMask$);

   'pH
   IF NOT (iWHAM) THEN
      PH! = -LOG10S!(CSNG(OC(1)))
   END IF
   IF NOT (iComplete) THEN
      PH! = -LOG10S!(CSNG(KT(1)))
   END IF
   PRINT #iHand&, USING$(OutputMask2$, PH!);

   'Diss. Metal
   IF iComplete AND NOT ErrorNotConverged THEN
      Tmp! = OKT(iMetal)
      IF CType(iMetal)= ctFixed THEN
          Tmp! = KT(iMetal)
      END IF
      IF (Tmp! > 9999) OR (Tmp! < 0.01) THEN
         OutputMask$ = OutputMask1$
      ELSE
         OutputMask$ = OutputMask2$
      END IF
      PRINT #iHand&, USING$(OutputMask$, Tmp!);
   ELSE
      PRINT #iHand&, FormatStr$("NA", LabelMask$);
   END IF

   'Al Dissolved/Precipitate information
   IF iLogKsAl THEN
      IF iComplete AND NOT ErrorNotConverged THEN
         Tmp! = DISS_Effect
         IF (Tmp! > 9999) OR (Tmp! < 0.01) THEN
            OutputMask$ = OutputMask1$
         ELSE
            OutputMask$ = OutputMask2$
         END IF
         PRINT #iHand&, USING$(OutputMask$, Tmp!);

         Tmp! = PART_Conc /(1000000*26.981)
         IF (Tmp! > 9999) OR (Tmp! < 0.01) THEN
            OutputMask$ = OutputMask1$
         ELSE
            OutputMask$ = OutputMask2$
         END IF
         PRINT #iHand&, USING$(OutputMask$, Tmp!);

         Tmp! = PART_Effect
         IF (Tmp! > 9999) OR (Tmp! < 0.01) THEN
            OutputMask$ = OutputMask1$
         ELSE
            OutputMask$ = OutputMask2$
         END IF
         PRINT #iHand&, USING$(OutputMask$, Tmp!);
      ELSE
         PRINT #iHand&, FormatStr$("NA", LabelMask$);
         PRINT #iHand&, FormatStr$("NA", LabelMask$);
         PRINT #iHand&, FormatStr$("NA", LabelMask$);
      END IF
   END IF

   'Free Metal
   IF iComplete AND NOT ErrorNotConverged THEN
      Tmp! = OC(iMetal)
      IF (Tmp! > 9999) OR (Tmp! < 0.01) THEN
         OutputMask$ = OutputMask1$
      ELSE
         OutputMask$ = OutputMask2$
      END IF
      PRINT #iHand&, USING$(OutputMask$, Tmp!);
   ELSE
      PRINT #iHand&, FormatStr$("NA", LabelMask$);
   END IF

   IF iComplete AND NOT ErrorNotConverged THEN
      Tmp! = 0!
      FOR j = 1 TO NSpecies
         Tmp! = Tmp! + wOrgBound!(j) * SCb(j, iMetal)
      NEXT j
      IF (Tmp! > 9999) OR (Tmp! < 0.01) THEN
         OutputMask$ = OutputMask1$
      ELSE
         OutputMask$ = OutputMask2$
      END IF
      PRINT #iHand&, USING$(OutputMask$, Tmp!);
   ELSE
      PRINT #iHand&, FormatStr$("NA", LabelMask$);
   END IF

   'DOC fractions
   IF (SepNOM = TRUE) THEN
     ' HA Fraction
     IF NOT iComplete THEN
        Tmp! = -999.#
        PRINT #iHand&, FormatStr$("NA", LabelMask$);
     ELSE
        Tmp! = 0!
        FOR j = 1 TO NSpecies
          Tmp! = Tmp! + HABound(j) * SCb(j, iMetal)
        NEXT j
        IF (Tmp! > 9999) OR (Tmp! < 0.01) THEN
           OutputMask$ = OutputMask1$
        ELSE
           OutputMask$ = OutputMask2$
        END IF
        PRINT #iHand&, USING$(OutputMask$, Tmp!);
     END IF


     ' FA Fraction
     IF NOT iComplete THEN
        Tmp! = -999.#
        PRINT #iHand&, "NA,";
     ELSE
        Tmp! = 0!
        FOR j = 1 TO NSpecies
          Tmp! = Tmp! + FABound(j) * SCb(j, iMetal)
        NEXT j
        IF (Tmp! > 9999) OR (Tmp! < 0.01) THEN
           OutputMask$ = OutputMask1$
        ELSE
           OutputMask$ = OutputMask2$
        END IF
        PRINT #iHand&, USING$(OutputMask$, Tmp!);
     END IF
   END IF

   FOR j = 1 TO NBLM
      IF iComplete AND NOT ErrorNotConverged THEN
         Tmp! = OC(iBLM(j)) * 1000000! / KT(iBL)
         IF (Tmp! > 9999) OR (Tmp! < 0.01) THEN
            OutputMask$ = OutputMask1$
         ELSE
            OutputMask$ = OutputMask2$
         END IF
         PRINT #iHand&, USING$(OutputMask$, Tmp!);
      ELSE
         PRINT #iHand&, FormatStr$("NA", LabelMask$);
      END IF
   NEXT j

   Tmp! = OKT(iDOC) * 1000! / 2
   PRINT #iHand&, USING$(OutputMask$, Tmp!);

   Tmp! = PercHA! * 100!
   PRINT #iHand&, USING$(OutputMask$, Tmp!);

   FOR c% = 1 TO iBL - 1
      IF ((c% <> 1) AND (c% <> iDOC) AND (c% <> iMetal) AND (c% <> iBL)) THEN
         Tmp! = OKT(c%)
         IF CType(c%)= ctFixed THEN
            Tmp! = KT(c%)
         END IF
         IF (Tmp! > 9999) OR (Tmp! < 0.01) THEN
            OutputMask$ = OutputMask1$
         ELSE
            OutputMask$ = OutputMask2$
         END IF
         PRINT #iHand&, USING$(OutputMask$, Tmp!);
       END IF
   NEXT c%
   PRINT #iHand&,

END SUB


SUB WriteWQCFile(iOutF4txt AS LONG, iOutF4xls AS LONG, Iter AS INTEGER, NData AS LONG)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB WriteWQCFile                                                     |
'  |                                                                       |
'  |  This version does acute and chronic WQC and FMB                      |
'  |                                                                       |
'  |  Analyzes the output from an entire BLM run and creates an output     |
'  |  file with with acute and chronic instantaneous WQC information,      |
'  |  toxic units, and if metal concentrations are available, it will      |
'  |  also calculate a fixed monitoring benchmark                          |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |                                                                       |
'  *-----------------------------------------------------------------------*
DIM value_int AS INTEGER
DIM MolWt AS DOUBLE

SELECT CASE Label(iMetal)
    CASE "Cu"
        MolWt = 63.546
    CASE "Zn"
        MolWt = 65.38
    CASE "Pb"
        MolWt = 207.2
    CASE ELSE
        ' shouldn't come here
        PRINT "Need to define another metal WQC"
END SELECT

  ' This is the first row that will contain the column labels for the bulk of the data
  iRowOS  = 6
  iRowOS  = 10
  iRowOS  = 11
  IF (Iter = 1) THEN
    '
    '   Write column headings
    '
    IF (iOFormat = THREE) OR (iCriteria <> ZERO) THEN   ' for now, let's always make a WQC file in xls format

'      xlsWriteDefaultFormats(iOutF4xls)

       iRow = iRowOS - 1
       iCol = 1
       sTmp1$ = "Site Label"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       iCol = 2
       sTmp1$ = "Sample Label"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       iCol = 3
       ' 1st line has the label
       sTmp1$ = "Final Acute Value"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       ' 2nd line has the units
       sTmp1$ = "(FAV), ug/L"
       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       iCol = 4
       sTmp1$ = "CMC"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       ' 2nd line has the units
       sTmp1$ = "(CMC=FAV/2), ug/L"
       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       iCol = 5
       sTmp1$ = "CCC"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       ' 2nd line has the units
       sTmp1$ = "(CCC=FAV/ACR), ug/L"
       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       iCol = 6
       sTmp1$ = sMetal$
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       ' 2nd line has the units
       sTmp1$ = "ug/L"
       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       iCol = 7
       sTmp1$ = "Acute Toxic Units"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       ' 2nd line has the units
       sTmp1$ = "(Acute TU=Cu/CMC)"
       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       iCol = 8
       sTmp1$ = "Chronic Toxic Units"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       ' 2nd line has the units
       sTmp1$ = "(Chronic TU=Cu/CCC)"
       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       iCol = 9
       sTmp1$ = "Censored Flag"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
        ' 2nd line has the units
       sTmp1$ = "(0 = quantified, 1 = BDL)"
       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       IF (iCalc_FMB = TRUE) THEN
           ' FMB Labels
           iCol = 4
           iRow = 5
           IF (FMB_option = FMBMedian) THEN
              sTmp1$ = "Acute FMB (median)"
           ELSE
              sTmp1$ = "Acute FMB (geo mean)"
           END IF
           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
           sTmp1$ = "ug/L"
           stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

           iCol = 5
           iRow = 5
           IF (FMB_option = FMBMedian) THEN
              sTmp1$ = "Chronic (median)"
           ELSE
              sTmp1$ = "Chronic (geo mean)"
           END IF
           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
           sTmp1$ = "ug/L"
           stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

           iCol = 4
           iRow = 8
           sTmp1$ = TRIM$(STR$(D_EF))
           sTmp1$ = "Based on an exceedance frequency of 1/"+sTmp1$+" days (once in "+TRIM$(FORMAT$(D_EF/365, "###.0#"))+" years)."

           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
        END IF 'iCalc_FMB = TRUE
    END IF 'IF (iOFormat = THREE) OR (iCriteria <> ZERO)
 END IF 'If (Iter = 1)
  '
  '   Write out data
  '
  'FOR j = 1 TO NData
      j = Iter
      IF (iOFormat = THREE) OR (iCriteria <> ZERO) THEN   ' for now, let's always make a WQC file in xls format
         iCol = 1
         iRow = iRowOS + j
         stat& = xlsWriteText(SiteL$, iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
         iCol = iCol + 1
         stat& = xlsWriteText(DataL$, iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)

         IF iComplete AND NOT ErrorNotConverged THEN
            '  Final acute value
            value# = HC5(j) * MolWt * 1.E6
            iCol = iCol + 1
            stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)

            '  CMC
            value# = HC5(j) * MolWt * 1.E6 / ACUTE_DIV
            WQC_Acute(j) = value#
            iCol = iCol + 1
            stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)

            '  CCC
            value# = HC5(j) * MolWt * 1.E6 / CHRONIC_DIV
            WQC_Chron(j) = value#
            iCol = iCol + 1
            stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
         ELSE
            'FAV
            iCol = iCol + 1
            stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

            'CMC
            iCol = iCol + 1
            stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

            'CCC
            iCol = iCol + 1
            stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
         END IF
         '
         '  Ambient metal concentration
         value# = AmbientMetal(j) * MolWt * 1.E6
         MetalConc(j) = value#
         iCol = iCol + 1
         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)

         IF iComplete AND NOT ErrorNotConverged THEN
            '  Acute toxic units
            value# = AmbientMetal(j) / (HC5(j) / ACUTE_DIV)
            TU_Acute(j) = value#
            iCol = iCol + 1
            stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
            '        '
            '  Chronic toxic units
            value# = AmbientMetal(j) / (HC5(j) / CHRONIC_DIV)
            TU_Chron(j) = value#
            iCol = iCol + 1
            stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
         ELSE
            'Acute TUs
            iCol = iCol + 1
            stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)

            'Chronic TUs
            iCol = iCol + 1
            stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 0, iHand&)
         END IF
         '
         '  Censored flag
         value_int = CenFlag(j)
         iCol = iCol + 1
         stat& = xlsWriteInteger(value_int, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)

     END IF
  'NEXT j
  IF (Iter = NData) AND (iCalc_FMB = TRUE) THEN
      ' Calculate FMB
      'MLE_Fit(x_dat_in() AS EXTENDED, x_dl() AS LONG, NData AS LONG, iLog AS LONG, x_mle_mean AS EXTENDED, x_mle_std AS EXTENDED)
      'AmbientMetal(i) = OKT(iMetal)
      'CenFlag(i) = OKT(iMetalDL)
      'HC5(j)
      'ToxicUnit(1 TO MaxObs)
      DIM TU_Mean AS EXTENDED
      DIM TU_median AS EXTENDED
      DIM TU_STD AS EXTENDED
      DIM TU_STD_4D AS EXTENDED
      DIM TU_CHR_Mean AS EXTENDED
      DIM TU_CHR_median AS EXTENDED
      DIM TU_CHR_STD AS EXTENDED
      DIM TU_CHR_STD_4D AS EXTENDED
      DIM Me_Mean AS EXTENDED
      DIM Me_median AS EXTENDED
      DIM Me_STD AS EXTENDED
      DIM Me_STD_4D AS EXTENDED
      DIM Me_Compliance AS EXTENDED
      DIM FMB_Acute AS EXTENDED
      DIM FMB_Chron AS EXTENDED
      DIM me_conc_log(1 TO NData) AS EXTENDED
      DIM TU_acut_log(1 TO NData) AS EXTENDED
      DIM TU_chro_log(1 TO NData) AS EXTENDED
      DIM Me_Conc_est(1 TO NData) AS EXTENDED
      DIM TU_acute_est(1 TO NData) AS EXTENDED
      DIM TU_chron_est(1 TO NData) AS EXTENDED
      DIM Me_Conc_est_log(1 TO NData) AS EXTENDED
      DIM TU_acute_est_log(1 TO NData) AS EXTENDED
      DIM TU_chron_est_log(1 TO NData) AS EXTENDED
      DIM x_AVE AS EXTENDED, x_ADEV AS EXTENDED, x_VAR AS EXTENDED, x_SKEW AS EXTENDED, x_CURT AS EXTENDED


      DIM N_Eff AS EXTENDED
      DIM N_Avg AS EXTENDED
      DIM Rho AS EXTENDED
      DIM NCen AS LONG

      DIM iLog AS LONG
      DIM NDataL AS LONG
      NDataL = NData
      iLog = TRUE

      '  ------------------------------------------------------------------
      '
      '  Calculate an FMB
      '
      '  Note - we have duplicate variables defined for the ambient metal conc
      '         MetalConc(i) and MetalFlag(i) and also AmbientMetal(i)
      '         Lets decide which we want to use.
      '
      '  This entire process should be moved into an independant subroutine
      '
      '  The exceedance frequency Z_EF needs to be read in from the parameter
      '         file, or from some other source (command line switch?)
      '         Maybe it could be read in from either location, with the
      '         command line switch over-riding the value in the parameter
      '         file, so that the interface will have a convenient way to
      '         provide a new value. I don't think we have a /Zhh switch yet
      '
      '  1 - use the metal concentrations, considering censored flags with MLE,
      '      and estimate the mean, median, and std dev
      '      using variables Me_Mean, Me_STD, Me_median
      '           - use actual data whenever possible, and only use MLE
      '             estimates for data when original data are BDL
      '  2 - do the same with the TU distribution using the same censored flags
      '      using variables TU_Mean, TU_STD, TU_median
      '  3 - Calculate the TU value at the exceedance frequency Z_EF
      '          TU_EF = 10^(Z_EF * TU_STD + LOG10(TU_median))
      '  4 - Calculate the adjustment factor that will result in a TU value
      '      of one at the exceedance frequency
      '          Adj_Fact = 1 / TU_EF
      '  5 - Calculate the median metal concentration that would correspond
      '      to a TU distribution that is in compliance
      '      Me_Compliance = Me_median * Adj_Fact
      '  6 - Calculate the FMB
      '      FMB = 10^(Z_EF * Me_STD + LOG10(Me_Compliance))
      '
      '  ------------------------------------------------------------------
      '  Perform an MLE fit with the TU values to calculate the degree of
      '  adjustment allowed for a TU distribution exactly in compliance
      '  with the target exceedance frequency
      '  ------------------------------------------------------------------
      NCen = 0
      FOR i=1 TO NData
          IF (CenFlag(i)=1) OR (CenFlag(i)=2) THEN
              NCen = NCen + 1
          END IF 'CenFlag(i)=1 OR CenFlag(i)=2
      NEXT i '=1 to NData

      IF (NCen > 0) THEN
        CALL MLE_Fit(TU_Acute(), CenFlag(), NDataL, iLog, TU_Mean, TU_STD)
        CALL MLE_Fit(TU_Chron(), CenFlag(), NDataL, iLog, TU_CHR_Mean, TU_CHR_STD)
        CALL MLE_Fit(MetalConc(), CenFlag(), NDataL, iLog, Me_Mean, Me_STD)
      END IF 'NCen>0

      DIM iDir AS LONG
      DIM Me_index(1 TO NData) AS LONG
      DIM TU_index(1 TO NData) AS LONG
      'dim Me_Zscore as extended
      'dim TU_Zscore as extended
      DIM x_index(1 TO NData) AS LONG
      'DIM x_zscore(1 TO NData) AS EXTENDED
      'DIM x_dat_mle(1 TO NData) AS EXTENDED
      'DIM x_median AS EXTENDED
      iDir = 0      ' Direction of sort - 0 = ascending, non-zero = descending
      'DIM i AS LONG
      '
      ' Sort once for toxic units
      CALL IndSortExt (TU_Acute(), TU_index(), NDataL, iDir)
      '
      ' and sort again for metal concentrations
      CALL IndSortExt (MetalConc(), Me_index(), NDataL, iDir)

      '
      '  Create an estimated dataset with substituted values from observations
      '  that were originally below detection limit (as determined by MLE fit),
      '  and with actual values whenever possible
      '
      DIM NDataL1 AS LONG      ' copy NData into a LONG variable to pass to ZScore() function
      NDataL1 = NData
      PRINT "------------------  debug  -----------------"
      PRINT "NDataL:"; NDataL
      PRINT "NData:"; NData
      PRINT "CenFlag"; TAB(20); "Me_Conc_est"; TAB(45); "TU_acute_est"; TAB(60); "TU_chron_est"
      FOR i=1 TO NData
          IF (CenFlag(i)=1) OR (CenFlag(i)=2) THEN
              value# = Me_Mean + ZScore(Me_index(i), NDataL1) * Me_STD
              Me_Conc_est(i) = 10^value#

              value# = TU_Mean + ZScore(TU_index(i), NDataL1) * TU_STD
              TU_acute_est(i) = 10^value#

              value# = TU_CHR_Mean + ZScore(TU_index(i), NDataL1) * TU_CHR_STD
              TU_chron_est(i) = 10^value#
          ELSE
              '        tmp1 = x_zscore(i)*x_std_mle + x_median_mle
              Me_Conc_est(i) = MetalConc(i)
              TU_Acute_est(i) = TU_Acute(i)
              TU_Chron_est(i) = TU_Chron(i)
          END IF 'CenFlag(i)=1 OR CenFlag(i)=2
          Me_Conc_est_log(i)  = LOG10(Me_Conc_est(i))
          TU_acute_est_log(i) = LOG10(TU_acute_est(i))
          TU_chron_est_log(i) = LOG10(TU_chron_est(i))
          PRINT CenFlag(i); TAB(15); USING$("##.###^^^^",Me_Conc_est(i)); TAB(30); USING$("##.###^^^^",TU_acute_est(i)); TAB(45); USING$("##.###^^^^",TU_chron_est(i))

      NEXT i '=1 to NData
      PRINT "---------------- end debug  -----------------"


      ' print out the original and adjusted metal conc
      DIM OMaskT$
      OMaskT$ =  "##.#######^^^^"
      PRINT
      PRINT "Orig_Me"; TAB(20); "Est_Me"; TAB(45); "Cen_Flag"
      FOR i=1 TO NData
          PRINT USING$(OMaskT$, MetalConc(i)); TAB(20); USING$(OMaskT$, Me_Conc_est(i)); TAB(45); CenFlag(i)
      NEXT i


      ' print out the original and adjusted acute toxic units
      PRINT
      PRINT "Orig_ATU"; TAB(20); "Est_ATU"; TAB(45); "Cen_Flag"
      FOR i=1 TO NData
          PRINT USING$(OMaskT$, TU_Acute(i)); TAB(20); USING$(OMaskT$, TU_acute_est(i)); TAB(45); CenFlag(i)
      NEXT i


      ' print out the original and adjusted chronic toxic units
      PRINT
      PRINT "Orig_CTU"; TAB(20); "Est_CTU"; TAB(45); "Cen_Flag"
      FOR i=1 TO NData
          PRINT USING$(OMaskT$, TU_chron(i)); TAB(20); USING$(OMaskT$, TU_chron_est(i)); TAB(45); CenFlag(i)
      NEXT i

      ' Sort again for toxic units using estimates, and this time calc a median
      ' For this - use log transformed values
      CALL IndSortExt (TU_Acute_est_log(), TU_index(), NDataL, iDir)
      CALL Median1(TU_Acute_est_log(), TU_index(), NDataL, TU_median)
      CALL Median1(TU_Chron_est_log(), TU_index(), NDataL, TU_CHR_median)
      TU_median = 10^TU_median
      TU_CHR_median = 10^TU_CHR_median
      '
      ' and do the same for metal concentrations
      ' again, using log transformed values
      CALL IndSortExt (Me_Conc_est_log(), Me_index(), NDataL, iDir)
      CALL Median1(Me_Conc_est_log(), Me_index(), NDataL, Me_median)
      Me_median = 10^Me_median
      '
      '  Now adjust the STD using Log transformed values
      CALL Moment (Me_Conc_est_log(), NDataL, x_AVE, x_ADEV, Me_STD, x_VAR, x_SKEW, x_CURT)
          Me_mean = 10^x_AVE
      CALL Moment (TU_Acute_est_log(), NDataL, x_AVE, x_ADEV, TU_STD, x_VAR, x_SKEW, x_CURT)
          TU_mean = 10^x_AVE
      CALL Moment (TU_Chron_est_log(), NDataL, x_AVE, x_ADEV, TU_CHR_STD  , x_VAR, x_SKEW, x_CURT)
          TU_CHR_mean = 10^x_AVE

      '
      '  Use this latest estimate of the median and STD to
      '  calulate acute FMB
      '
      IF (FMB_option = FMBMedian) THEN
        TU_EF = 10^(Z_EF * TU_STD + LOG10(TU_median))
        Adj_Fact = 1 / TU_EF
        Me_Compliance = Me_median * Adj_Fact
        FMB_Acute = 10^(Z_EF * Me_STD + LOG10(Me_Compliance))
      ELSE
        TU_EF = 10^(Z_EF * TU_STD + LOG10(TU_mean))
        Adj_Fact = 1 / TU_EF
        Me_Compliance = Me_mean * Adj_Fact
        FMB_Acute = 10^(Z_EF * Me_STD + LOG10(Me_Compliance))
      END IF 'FMB_Option = FMBMedian


      PRINT
      PRINT "metal_median" ; Me_median
      PRINT "metal_mean"   ; Me_mean
      PRINT "metal_std"    ; Me_STD
      PRINT
      PRINT "TU_median" ; TU_median
      PRINT "TU_mean"   ; TU_mean
      PRINT "TU_std"    ; TU_STD

      PRINT
      PRINT "Z_EF"           ;  Z_EF
      PRINT "TU_EFF"         ;  TU_EF
      PRINT "ADJ_F"          ;  Adj_Fact
      PRINT "A_Me_Comp"      ;  Me_Compliance

      PRINT
      PRINT "Chron_TU_median" ; TU_CHR_median
      PRINT "Chron_TU_mean"   ; TU_CHR_mean
      PRINT "Chron_TU_std"    ; TU_CHR_STD

      '
      '  Estimate the 4-day avg toxic unit STD
      '
      N_Avg = 4.
      Rho = 0.8
      N_Eff = N_Avg^2*(1-Rho)^2/(N_Avg*(1-Rho^2)-(2*Rho*(1-Rho^N_Avg)))
      Me_STD_4D = (LOG(1+(EXP(Me_STD^2)-1)/N_Eff))^0.5
      TU_CHR_STD_4D = (LOG(1+(EXP(TU_CHR_STD^2)-1)/N_Eff))^0.5
      PRINT
      PRINT "N_Eff"           ;  N_Eff
      PRINT "metal_4D_std"    ;  Me_STD_4D
      PRINT "Chron_TU_4D_std" ;  TU_CHR_STD_4D
      PRINT       '
      '  and use that to calculate a chronic FMB
      '
      IF (FMB_option = FMBMedian) THEN
        TU_EF = 10^(Z_EF * TU_CHR_STD_4D + LOG10(TU_CHR_median))
        Adj_Fact = 1 / TU_EF
        Me_Compliance = Me_median * Adj_Fact
        FMB_Chron = 10^(Z_EF * Me_STD_4D + LOG10(Me_Compliance))
      ELSE
        TU_EF = 10^(Z_EF * TU_CHR_STD_4D + LOG10(TU_CHR_mean))
        Adj_Fact = 1 / TU_EF
        Me_Compliance = Me_mean * Adj_Fact
        FMB_Chron = 10^(Z_EF * Me_STD_4D + LOG10(Me_Compliance))
      END IF 'FMB_option = FMBMedian

     value# = FMB_Acute
     iRow = 7
     iCol = 4
     'stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
     stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
     value# = FMB_Chron
     iRow = 7
     iCol = 5
     stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)

      PRINT "TU_EFF_C"        ;   TU_EF
      PRINT "ADJ_F_C"         ;   Adj_Fact
      PRINT "C_Me_Comp"       ;   Me_Compliance
      PRINT
      PRINT "FMB_Acute"    ; FMB_Acute
      PRINT "FMB_Chron"    ; FMB_Chron
      PRINT

'     FOR i=1 TO NData
'        me_conc_log(i) =  LOG10(MetalConc(i))
'        TU_acut_log(i) =  LOG10(TU_Acute(i))
'        TU_chro_log(i) =  LOG10(TU_Chron(i))
'     NEXT i
'     CALL IndSortExt (TU_acut_log(), x_index(), NDataL, iDir)
'     CALL Median1(TU_acut_log(), x_index(), NDataL, TU_median)
'
'     CALL IndSortExt (TU_chro_log(), x_index(), NDataL, iDir)
'     CALL Median1(TU_chro_log(), x_index(), NDataL, TU_CHR_median)
'
'     CALL IndSortExt (me_conc_log(), x_index(), NDataL, iDir)
'     CALL Median1(me_conc_log(), x_index(), NDataL, Me_median)
'
'     PRINT "log metal median" ; Me_median
'     PRINT "log TU median" ; TU_median
'     PRINT "log Chron TU median" ; TU_CHR_median
'
'     TU_EF = 10^(Z_EF * TU_CHR_STD_4D + LOG10(TU_CHR_median))
'     Adj_Fact = 1 / TU_EF
'     Me_Compliance = Me_median * Adj_Fact
'     FMB_Chron = 10^(Z_EF * Me_STD_4D + LOG10(Me_Compliance))
'
'     TU_EF = 10^(Z_EF * TU_STD + LOG10(TU_median))
'     Adj_Fact = 1 / TU_EF
'     Me_Compliance = Me_median * Adj_Fact
'     FMB_Acute = 10^(Z_EF * Me_STD + LOG10(Me_Compliance))
'     print
'     print "using log median:"
'     PRINT "FMB Acute"    ; FMB_Acute
'     PRINT "FMB Chron"    ; FMB_Chron

  END IF 'Iter = NData AND iCalc_FMB = TRUE
END SUB 'WriteWQCFile

FUNCTION FMB_AcuteFunc(NData AS LONG) AS EXTENDED
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  FUNCTION FMB_AcuteFunc                                               |
'  |                                                                       |
'  |  Calculate an FMB
'  |
'  |  Note - we have duplicate variables defined for the ambient metal conc
'  |      MetalConc(i) and MetalFlag(i) and also AmbientMetal(i)
'  |      Lets decide which we want to use.
'  |
'  |  This entire process should be moved into an independant subroutine
'  |
'  |  The exceedance frequency Z_EF needs to be read in from the parameter
'  |        file, or from some other source (command line switch?)
'  |        Maybe it could be read in from either location, with the
'  |        command line switch over-riding the value in the parameter
'  |        file, so that the interface will have a convenient way to
'  |        provide a new value. I don't think we have a /Zhh switch yet
'  |
'  |  1 - use the metal concentrations, considering censored flags with MLE,
'  |      and estimate the mean, median, and std dev
'  |      using variables Me_Mean, Me_STD, Me_median
'  |         - use actual data whenever possible, and only use MLE
'  |           estimates for data when original data are BDL
'  |  2 - do the same with the TU distribution using the same censored flags
'  |      using variables TU_Mean, TU_STD, TU_median
'  |  3 - Calculate the TU value at the exceedance frequency Z_EF
'  |         TU_EF = 10^(Z_EF * TU_STD + LOG10(TU_median))
'  |  4 - Calculate the adjustment factor that will result in a TU value
'  |      of one at the exceedance frequency
'  |         Adj_Fact = 1 / TU_EF
'  |  5 - Calculate the median metal concentration that would correspond
'  |      to a TU distribution that is in compliance
'  |      Me_Compliance = Me_median * Adj_Fact
'  |  6 - Calculate the FMB
'  |      FMB = 10^(Z_EF * Me_STD + LOG10(Me_Compliance))
  DIM TU_Mean AS EXTENDED
  DIM TU_median AS EXTENDED
  DIM TU_STD AS EXTENDED
  DIM TU_STD_4D AS EXTENDED
  DIM TU_CHR_Mean AS EXTENDED
  DIM TU_CHR_median AS EXTENDED
  DIM TU_CHR_STD AS EXTENDED
  DIM TU_CHR_STD_4D AS EXTENDED
  DIM Me_Mean AS EXTENDED
  DIM Me_median AS EXTENDED
  DIM Me_STD AS EXTENDED
  DIM Me_STD_4D AS EXTENDED
  DIM Me_Compliance AS EXTENDED
  DIM FMB_Acute AS EXTENDED
  DIM FMB_Chron AS EXTENDED
  DIM me_conc_log(1 TO NData) AS EXTENDED
  DIM TU_acut_log(1 TO NData) AS EXTENDED
  DIM TU_chro_log(1 TO NData) AS EXTENDED
  DIM Me_Conc_est(1 TO NData) AS EXTENDED
  DIM TU_acute_est(1 TO NData) AS EXTENDED
  DIM TU_chron_est(1 TO NData) AS EXTENDED
  DIM Me_Conc_est_log(1 TO NData) AS EXTENDED
  DIM TU_acute_est_log(1 TO NData) AS EXTENDED
  DIM TU_chron_est_log(1 TO NData) AS EXTENDED
  DIM x_AVE AS EXTENDED, x_ADEV AS EXTENDED, x_VAR AS EXTENDED, x_SKEW AS EXTENDED, x_CURT AS EXTENDED

  DIM N_Eff AS EXTENDED
  DIM N_Avg AS EXTENDED
  DIM Rho AS EXTENDED
  DIM NCen AS LONG

  DIM iLog AS LONG
    iLog = TRUE
  DIM NDataL AS LONG
    NDataL = NData

  NCen = 0
  FOR i=1 TO NData
    IF (CenFlag(i)=1) OR (CenFlag(i)=2) THEN
      NCen = NCen + 1
    END IF
  NEXT i

  '  *-----------------------------------------------------------------------------*
  '  | If there are any censored values, perform an MLE fit with the TU values     |
  '  | to calculate the most likely distribution that corresponds to the measured  |
  '  | and censored data.                                                          |
  '  *-----------------------------------------------------------------------------*
  IF (NCen > 0) THEN
    CALL MLE_Fit(TU_Acute(), CenFlag(), NDataL, iLog, TU_Mean, TU_STD)
    CALL MLE_Fit(TU_Chron(), CenFlag(), NDataL, iLog, TU_CHR_Mean, TU_CHR_STD)
    CALL MLE_Fit(MetalConc(), CenFlag(), NDataL, iLog, Me_Mean, Me_STD)
  END IF

  '  *-----------------------------------------------------------------------------*
  '  | Sort the metal concentrations and toxic unit values and store the sorted    |
  '  | order in an index variable.  We'll need the sorted order in order to find   |
  '  | the position in the distribution for calculating the zscore for values we   |
  '  | need to estimate from the MLE distribution.                                 |
  '  *-----------------------------------------------------------------------------*
      DIM iDir AS LONG
      DIM Me_index(1 TO NData) AS LONG
      DIM TU_index(1 TO NData) AS LONG
      'dim Me_Zscore as extended
      'dim TU_Zscore as extended
      'DIM x_index(1 TO NData) AS LONG
      'DIM x_zscore(1 TO NData) AS EXTENDED
      'DIM x_dat_mle(1 TO NData) AS EXTENDED
      'DIM x_median AS EXTENDED
      iDir = 0      ' Direction of sort - 0 = ascending, non-zero = descending
      'DIM i AS LONG
      '
      ' Sort once for toxic units
      CALL IndSortExt (TU_Acute(), TU_index(), NDataL, iDir)
      '
      ' and sort again for metal concentrations
      CALL IndSortExt (MetalConc(), Me_index(), NDataL, iDir)

      '
      '  Create an estimated dataset with substituted values from observations
      '  that were originally below detection limit (as determined by MLE fit),
      '  and with actual values whenever possible
      '
      DIM NDataL1 AS LONG      ' copy NData into a LONG variable to pass to ZScore() function
      NDataL1 = NData
      PRINT "------------------  debug  -----------------"
      PRINT "NDataL:"; NDataL
      PRINT "NData:"; NData
      PRINT "CenFlag"; TAB(20); "Me_Conc_est"; TAB(45); "TU_acute_est"; TAB(60); "TU_chron_est"
      FOR i=1 TO NData
          IF (CenFlag(i)=1) OR (CenFlag(i)=2) THEN
              value# = Me_Mean + ZScore(Me_index(i), NDataL1) * Me_STD
              Me_Conc_est(i) = 10^value#

              value# = TU_Mean + ZScore(TU_index(i), NDataL1) * TU_STD
              TU_acute_est(i) = 10^value#

              value# = TU_CHR_Mean + ZScore(TU_index(i), NDataL1) * TU_CHR_STD
              TU_chron_est(i) = 10^value#
          ELSE
              '        tmp1 = x_zscore(i)*x_std_mle + x_median_mle
              Me_Conc_est(i) = MetalConc(i)
              TU_Acute_est(i) = TU_Acute(i)
              TU_Chron_est(i) = TU_Chron(i)
          END IF
          Me_Conc_est_log(i)  = LOG10(Me_Conc_est(i))
          TU_acute_est_log(i) = LOG10(TU_acute_est(i))
          TU_chron_est_log(i) = LOG10(TU_chron_est(i))
          PRINT CenFlag(i); TAB(15); USING$("##.###^^^^",Me_Conc_est(i)); TAB(30); USING$("##.###^^^^",TU_acute_est(i)); TAB(45); USING$("##.###^^^^",TU_chron_est(i))

      NEXT i
      PRINT "---------------- end debug  -----------------"


      ' print out the original and adjusted metal conc
      DIM OMaskT$
      OMaskT$ =  "##.#######^^^^"
      PRINT
      PRINT "Orig_Me"; TAB(20); "Est_Me"; TAB(45); "Cen_Flag"
      FOR i=1 TO NData
          PRINT USING$(OMaskT$, MetalConc(i)); TAB(20); USING$(OMaskT$, Me_Conc_est(i)); TAB(45); CenFlag(i)
      NEXT i


      ' print out the original and adjusted acute toxic units
      PRINT
      PRINT "Orig_ATU"; TAB(20); "Est_ATU"; TAB(45); "Cen_Flag"
      FOR i=1 TO NData
          PRINT USING$(OMaskT$, TU_Acute(i)); TAB(20); USING$(OMaskT$, TU_acute_est(i)); TAB(45); CenFlag(i)
      NEXT i


      ' print out the original and adjusted chronic toxic units
      PRINT
      PRINT "Orig_CTU"; TAB(20); "Est_CTU"; TAB(45); "Cen_Flag"
      FOR i=1 TO NData
          PRINT USING$(OMaskT$, TU_chron(i)); TAB(20); USING$(OMaskT$, TU_chron_est(i)); TAB(45); CenFlag(i)
      NEXT i

      ' Sort again for toxic units using estimates, and this time calc a median
      ' For this - use log transformed values
      CALL IndSortExt (TU_Acute_est_log(), TU_index(), NDataL, iDir)
      CALL Median1(TU_Acute_est_log(), TU_index(), NDataL, TU_median)
      CALL Median1(TU_Chron_est_log(), TU_index(), NDataL, TU_CHR_median)
      TU_median = 10^TU_median
      TU_CHR_median = 10^TU_CHR_median
      '
      ' and do the same for metal concentrations
      ' again, using log transformed values
      CALL IndSortExt (Me_Conc_est_log(), Me_index(), NDataL, iDir)
      CALL Median1(Me_Conc_est_log(), Me_index(), NDataL, Me_median)
      Me_median = 10^Me_median
      '
      '  Now adjust the STD using Log transformed values
      CALL Moment (Me_Conc_est_log(), NDataL, x_AVE, x_ADEV, Me_STD, x_VAR, x_SKEW, x_CURT)
          Me_mean = 10^x_AVE
      CALL Moment (TU_Acute_est_log(), NDataL, x_AVE, x_ADEV, TU_STD, x_VAR, x_SKEW, x_CURT)
          TU_mean = 10^x_AVE
      CALL Moment (TU_Chron_est_log(), NDataL, x_AVE, x_ADEV, TU_CHR_STD  , x_VAR, x_SKEW, x_CURT)
          TU_CHR_mean = 10^x_AVE

      '
      '  Use this latest estimate of the median and STD to
      '  calulate acute FMB
      '
      IF (FMB_option = FMBMedian) THEN
        TU_EF = 10^(Z_EF * TU_STD + LOG10(TU_median))
        Adj_Fact = 1 / TU_EF
        Me_Compliance = Me_median * Adj_Fact
        FMB_Acute = 10^(Z_EF * Me_STD + LOG10(Me_Compliance))
      ELSE
        TU_EF = 10^(Z_EF * TU_STD + LOG10(TU_mean))
        Adj_Fact = 1 / TU_EF
        Me_Compliance = Me_mean * Adj_Fact
        FMB_Acute = 10^(Z_EF * Me_STD + LOG10(Me_Compliance))
      END IF


      PRINT
      PRINT "metal_median" ; Me_median
      PRINT "metal_mean"   ; Me_mean
      PRINT "metal_std"    ; Me_STD
      PRINT
      PRINT "TU_median" ; TU_median
      PRINT "TU_mean"   ; TU_mean
      PRINT "TU_std"    ; TU_STD

      PRINT
      PRINT "Z_EF"           ;  Z_EF
      PRINT "TU_EFF"         ;  TU_EF
      PRINT "ADJ_F"          ;  Adj_Fact
      PRINT "A_Me_Comp"      ;  Me_Compliance

      PRINT
      PRINT "Chron_TU_median" ; TU_CHR_median
      PRINT "Chron_TU_mean"   ; TU_CHR_mean
      PRINT "Chron_TU_std"    ; TU_CHR_STD

      '
      '  Estimate the 4-day avg toxic unit STD
      '
      N_Avg = 4.
      Rho = 0.8
      N_Eff = N_Avg^2*(1-Rho)^2/(N_Avg*(1-Rho^2)-(2*Rho*(1-Rho^N_Avg)))
      Me_STD_4D = (LOG(1+(EXP(Me_STD^2)-1)/N_Eff))^0.5
      TU_CHR_STD_4D = (LOG(1+(EXP(TU_CHR_STD^2)-1)/N_Eff))^0.5
      PRINT
      PRINT "N_Eff"           ;  N_Eff
      PRINT "metal_4D_std"    ;  Me_STD_4D
      PRINT "Chron_TU_4D_std" ;  TU_CHR_STD_4D
      PRINT       '
      '  and use that to calculate a chronic FMB
      '
      IF (FMB_option = FMBMedian) THEN
        TU_EF = 10^(Z_EF * TU_CHR_STD_4D + LOG10(TU_CHR_median))
        Adj_Fact = 1 / TU_EF
        Me_Compliance = Me_median * Adj_Fact
        FMB_Chron = 10^(Z_EF * Me_STD_4D + LOG10(Me_Compliance))
      ELSE
        TU_EF = 10^(Z_EF * TU_CHR_STD_4D + LOG10(TU_CHR_mean))
        Adj_Fact = 1 / TU_EF
        Me_Compliance = Me_mean * Adj_Fact
        FMB_Chron = 10^(Z_EF * Me_STD_4D + LOG10(Me_Compliance))
      END IF
  ' *---------------------------------------------------------------------------------
  ' |  Return FMB values and end the function
  ' *---------------------------------------------------------------------------------
  FMB_AcuteFunc = FMB_Acute
END FUNCTION 'FMB_AcuteFunc

SUB FMBCalc(NData AS LONG, FMB_Acute AS DOUBLE, FMB_Chron AS DOUBLE)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  FUNCTION FMB_AcuteFunc                                               |
'  |                                                                       |
'  |  Calculate an FMB
'  |
'  |  Note - we have duplicate variables defined for the ambient metal conc
'  |      MetalConc(i) and MetalFlag(i) and also AmbientMetal(i)
'  |      Lets decide which we want to use.
'  |
'  |  This entire process should be moved into an independant subroutine
'  |
'  |  The exceedance frequency Z_EF needs to be read in from the parameter
'  |        file, or from some other source (command line switch?)
'  |        Maybe it could be read in from either location, with the
'  |        command line switch over-riding the value in the parameter
'  |        file, so that the interface will have a convenient way to
'  |        provide a new value. I don't think we have a /Zhh switch yet
'  |
'  |  1 - use the metal concentrations, considering censored flags with MLE,
'  |      and estimate the mean, median, and std dev
'  |      using variables Me_Mean, Me_STD, Me_median
'  |         - use actual data whenever possible, and only use MLE
'  |           estimates for data when original data are BDL
'  |  2 - do the same with the TU distribution using the same censored flags
'  |      using variables TU_Mean, TU_STD, TU_median
'  |  3 - Calculate the TU value at the exceedance frequency Z_EF
'  |         TU_EF = 10^(Z_EF * TU_STD + LOG10(TU_median))
'  |  4 - Calculate the adjustment factor that will result in a TU value
'  |      of one at the exceedance frequency
'  |         Adj_Fact = 1 / TU_EF
'  |  5 - Calculate the median metal concentration that would correspond
'  |      to a TU distribution that is in compliance
'  |      Me_Compliance = Me_median * Adj_Fact
'  |  6 - Calculate the FMB
'  |      FMB = 10^(Z_EF * Me_STD + LOG10(Me_Compliance))
  DIM TU_Mean AS EXTENDED
  DIM TU_median AS EXTENDED
  DIM TU_STD AS EXTENDED
  DIM TU_STD_4D AS EXTENDED
  DIM TU_CHR_Mean AS EXTENDED
  DIM TU_CHR_median AS EXTENDED
  DIM TU_CHR_STD AS EXTENDED
  DIM TU_CHR_STD_4D AS EXTENDED
  DIM Me_Mean AS EXTENDED
  DIM Me_median AS EXTENDED
  DIM Me_STD AS EXTENDED
  DIM Me_STD_4D AS EXTENDED
  DIM Me_Compliance AS EXTENDED
  'DIM FMB_Acute AS EXTENDED
  'DIM FMB_Chron AS EXTENDED
  DIM me_conc_log(1 TO NData) AS EXTENDED
  DIM TU_acut_log(1 TO NData) AS EXTENDED
  DIM TU_chro_log(1 TO NData) AS EXTENDED
  DIM Me_Conc_est(1 TO NData) AS EXTENDED
  DIM TU_acute_est(1 TO NData) AS EXTENDED
  DIM TU_chron_est(1 TO NData) AS EXTENDED
  DIM Me_Conc_est_log(1 TO NData) AS EXTENDED
  DIM TU_acute_est_log(1 TO NData) AS EXTENDED
  DIM TU_chron_est_log(1 TO NData) AS EXTENDED
  DIM x_AVE AS EXTENDED, x_ADEV AS EXTENDED, x_VAR AS EXTENDED, x_SKEW AS EXTENDED, x_CURT AS EXTENDED

  DIM N_Eff AS EXTENDED
  DIM N_Avg AS EXTENDED
  DIM Rho AS EXTENDED
  DIM NCen AS LONG

  DIM iLog AS LONG
    iLog = TRUE
  DIM NDataL AS LONG
    NDataL = NData

  NCen = 0
  FOR i=1 TO NData
    IF (CenFlag(i)=1) OR (CenFlag(i)=2) THEN
      NCen = NCen + 1
    END IF
  NEXT i



  '  *-----------------------------------------------------------------------------*
  '  | If there are any censored values, perform an MLE fit with the TU values     |
  '  | to calculate the most likely distribution that corresponds to the measured  |
  '  | and censored data.                                                          |
  '  *-----------------------------------------------------------------------------*
  IF (NCen > 0) THEN
    IF (iCriteria <> THREE) THEN
       CALL MLE_Fit(TU_Acute(), CenFlag(), NDataL, iLog, TU_Mean, TU_STD)
    END IF

    IF (iCriteria <> TWO) THEN
       CALL MLE_Fit(TU_Chron(), CenFlag(), NDataL, iLog, TU_CHR_Mean, TU_CHR_STD)
    END IF

    CALL MLE_Fit(MetalConc(), CenFlag(), NDataL, iLog, Me_Mean, Me_STD)
  END IF

  '  *-----------------------------------------------------------------------------*
  '  | Sort the metal concentrations and toxic unit values and store the sorted    |
  '  | order in an index variable.  We'll need the sorted order in order to find   |
  '  | the position in the distribution for calculating the zscore for values we   |
  '  | need to estimate from the MLE distribution.                                 |
  '  *-----------------------------------------------------------------------------*
      DIM iDir AS LONG
      DIM Me_index(1 TO NData) AS LONG
      DIM TU_index(1 TO NData) AS LONG
      'dim Me_Zscore as extended
      'dim TU_Zscore as extended
      'DIM x_index(1 TO NData) AS LONG
      'DIM x_zscore(1 TO NData) AS EXTENDED
      'DIM x_dat_mle(1 TO NData) AS EXTENDED
      'DIM x_median AS EXTENDED
      iDir = 0      ' Direction of sort - 0 = ascending, non-zero = descending
      'DIM i AS LONG
      '
      ' Sort once for toxic units
      CALL IndSortExt (TU_Acute(), TU_index(), NDataL, iDir)
      '
      ' and sort again for metal concentrations
      CALL IndSortExt (MetalConc(), Me_index(), NDataL, iDir)

      '
      '  Create an estimated dataset with substituted values from observations
      '  that were originally below detection limit (as determined by MLE fit),
      '  and with actual values whenever possible
      '
      DIM NDataL1 AS LONG      ' copy NData into a LONG variable to pass to ZScore() function
      NDataL1 = NData
      PRINT "------------------  debug  -----------------"
      PRINT "NDataL:"; NDataL
      PRINT "NData:"; NData
      PRINT "CenFlag"; TAB(20); "Me_Conc_est";
      IF (icriteria <> THREE) THEN
         PRINT TAB(45); "TU_acute_est";
      END IF
      IF (iCriteria <> TWO) THEN
         PRINT TAB(60); "TU_chron_est"
      END IF
      FOR i=1 TO NData
          IF (CenFlag(i)=1) OR (CenFlag(i)=2) THEN
              value# = Me_Mean + ZScore(Me_index(i), NDataL1) * Me_STD
              Me_Conc_est(i) = 10^value#

              IF (iCriteria <> THREE) THEN
                  value# = TU_Mean + ZScore(TU_index(i), NDataL1) * TU_STD
                  TU_acute_est(i) = 10^value#
              END IF

              IF (iCriteria <> TWO) THEN
                 value# = TU_CHR_Mean + ZScore(TU_index(i), NDataL1) * TU_CHR_STD
                 TU_chron_est(i) = 10^value#
              END IF

          ELSE
              '        tmp1 = x_zscore(i)*x_std_mle + x_median_mle
              Me_Conc_est(i) = MetalConc(i)
              IF (iCriteria <> THREE) THEN
                 TU_Acute_est(i) = TU_Acute(i)
              END IF
              IF (iCriteria <> TWO) THEN
                 TU_Chron_est(i) = TU_Chron(i)
              END IF
          END IF
          Me_Conc_est_log(i)  = LOG10(Me_Conc_est(i))
          IF (iCriteria <> THREE) THEN
             TU_acute_est_log(i) = LOG10(TU_acute_est(i))
          END IF
          IF (iCriteria <> TWO) THEN
             TU_chron_est_log(i) = LOG10(TU_chron_est(i))
          END IF
          PRINT CenFlag(i); TAB(15); USING$("##.###^^^^",Me_Conc_est(i));
          IF (iCriteria <> THREE) THEN
             PRINT TAB(30); USING$("##.###^^^^",TU_acute_est(i));
          END IF
          IF (iCriteria <> TWO) THEN
             PRINT TAB(45); USING$("##.###^^^^",TU_chron_est(i))
          END IF

      NEXT i
      PRINT "---------------- end debug  -----------------"


      ' print out the original and adjusted metal conc
      DIM OMaskT$
      OMaskT$ =  "##.#######^^^^"
      PRINT
      PRINT "Orig_Me"; TAB(20); "Est_Me"; TAB(45); "Cen_Flag"
      FOR i=1 TO NData
          PRINT USING$(OMaskT$, MetalConc(i)); TAB(20); USING$(OMaskT$, Me_Conc_est(i)); TAB(45); CenFlag(i)
      NEXT i


      ' print out the original and adjusted acute toxic units
      IF (iCriteria <> THREE) THEN
         PRINT
         PRINT "Orig_ATU"; TAB(20); "Est_ATU"; TAB(45); "Cen_Flag"
         FOR i=1 TO NData
             PRINT USING$(OMaskT$, TU_Acute(i)); TAB(20); USING$(OMaskT$, TU_acute_est(i)); TAB(45); CenFlag(i)
         NEXT i
      END IF


      ' print out the original and adjusted chronic toxic units
      IF (iCriteria <> TWO) THEN
         PRINT
         PRINT "Orig_CTU"; TAB(20); "Est_CTU"; TAB(45); "Cen_Flag"
         FOR i=1 TO NData
             PRINT USING$(OMaskT$, TU_chron(i)); TAB(20); USING$(OMaskT$, TU_chron_est(i)); TAB(45); CenFlag(i)
         NEXT i
      END IF

      ' Sort again for toxic units using estimates, and this time calc a median
      ' For this - use log transformed values
      IF (iCriteria <> THREE) THEN
         CALL IndSortExt (TU_Acute_est_log(), TU_index(), NDataL, iDir)
      ELSE
         CALL IndSortExt (TU_Chron_est_log(), TU_index(), NDataL, iDir)
      END IF
      IF (iCriteria <> THREE) THEN
         CALL Median1(TU_Acute_est_log(), TU_index(), NDataL, TU_median)
         TU_median = 10^TU_median
      END IF
      IF (iCriteria <> TWO) THEN
         CALL Median1(TU_Chron_est_log(), TU_index(), NDataL, TU_CHR_median)
         TU_CHR_median = 10^TU_CHR_median
      END IF

      '
      ' and do the same for metal concentrations
      ' again, using log transformed values
      CALL IndSortExt (Me_Conc_est_log(), Me_index(), NDataL, iDir)
      CALL Median1(Me_Conc_est_log(), Me_index(), NDataL, Me_median)
      Me_median = 10^Me_median
      '
      '  Now adjust the STD using Log transformed values
      CALL Moment (Me_Conc_est_log(), NDataL, x_AVE, x_ADEV, Me_STD, x_VAR, x_SKEW, x_CURT)
          Me_mean = 10^x_AVE
      IF (iCriteria <> THREE) THEN
         CALL Moment (TU_Acute_est_log(), NDataL, x_AVE, x_ADEV, TU_STD, x_VAR, x_SKEW, x_CURT)
         TU_mean = 10^x_AVE
      END IF
      IF (iCriteria <> TWO) THEN
         CALL Moment (TU_Chron_est_log(), NDataL, x_AVE, x_ADEV, TU_CHR_STD  , x_VAR, x_SKEW, x_CURT)
         TU_CHR_mean = 10^x_AVE
      END IF

      PRINT
      PRINT "metal_median" ; Me_median
      PRINT "metal_mean"   ; Me_mean
      PRINT "metal_std"    ; Me_STD

      PRINT
      PRINT "Z_EF"           ;  Z_EF

      '
      '  Use this latest estimate of the median and STD to
      '  calulate acute FMB
      '
      FMB_Acute = -999
      IF (iCriteria <> THREE) THEN
         IF (FMB_option = FMBMedian) THEN
           TU_EF = 10^(Z_EF * TU_STD + LOG10(TU_median))
           Adj_Fact = 1 / TU_EF
           Me_Compliance = Me_median * Adj_Fact
           FMB_Acute = 10^(Z_EF * Me_STD + LOG10(Me_Compliance))
         ELSE
           TU_EF = 10^(Z_EF * TU_STD + LOG10(TU_mean))
           Adj_Fact = 1 / TU_EF
           Me_Compliance = Me_mean * Adj_Fact
           FMB_Acute = 10^(Z_EF * Me_STD + LOG10(Me_Compliance))
         END IF

         PRINT
         PRINT "TU_median" ; TU_median
         PRINT "TU_mean"   ; TU_mean
         PRINT "TU_std"    ; TU_STD
         PRINT
         PRINT "TU_EFF"         ;  TU_EF
         PRINT "ADJ_F"          ;  Adj_Fact
         PRINT "A_Me_Comp"      ;  Me_Compliance
      END IF


      FMB_Chron = -999
      IF (iCriteria <> TWO) THEN
         PRINT
         PRINT "Chron_TU_median" ; TU_CHR_median
         PRINT "Chron_TU_mean"   ; TU_CHR_mean
         PRINT "Chron_TU_std"    ; TU_CHR_STD

         '
         '  Estimate the 4-day avg toxic unit STD
         '
         N_Avg = 4.
         Rho = 0.8
         N_Eff = N_Avg^2*(1-Rho)^2/(N_Avg*(1-Rho^2)-(2*Rho*(1-Rho^N_Avg)))
         Me_STD_4D = (LOG(1+(EXP(Me_STD^2)-1)/N_Eff))^0.5
         TU_CHR_STD_4D = (LOG(1+(EXP(TU_CHR_STD^2)-1)/N_Eff))^0.5
         PRINT
         PRINT "N_Eff"           ;  N_Eff
         PRINT "metal_4D_std"    ;  Me_STD_4D
         PRINT "Chron_TU_4D_std" ;  TU_CHR_STD_4D
         PRINT       '
         '  and use that to calculate a chronic FMB
         '
         IF (FMB_option = FMBMedian) THEN
           TU_EF = 10^(Z_EF * TU_CHR_STD_4D + LOG10(TU_CHR_median))
           Adj_Fact = 1 / TU_EF
           Me_Compliance = Me_median * Adj_Fact
           FMB_Chron = 10^(Z_EF * Me_STD_4D + LOG10(Me_Compliance))
         ELSE
           TU_EF = 10^(Z_EF * TU_CHR_STD_4D + LOG10(TU_CHR_mean))
           Adj_Fact = 1 / TU_EF
           Me_Compliance = Me_mean * Adj_Fact
           FMB_Chron = 10^(Z_EF * Me_STD_4D + LOG10(Me_Compliance))
         END IF
      END IF
  ' *---------------------------------------------------------------------------------
  ' |  Return FMB values and end the function
  ' *---------------------------------------------------------------------------------



END SUB 'FMBCalc

SUB WriteWQCFileAcute(iOutF4txt AS LONG, iOutF4xls AS LONG, Iter AS INTEGER, NData AS LONG)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB WriteWQCFileAcute                                                |
'  |                                                                       |
'  |  This version does acute WQC and FMB (no chronic)                     |
'  |                                                                       |
'  |  Analyzes the output from an entire BLM run and creates an output     |
'  |  file with with acute instantaneous WQC information, toxic units,     |
'  |  and if metal concentrations are available, it will also calculate    |
'  |  a fixed monitoring benchmark                                         |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |                                                                       |
'  *-----------------------------------------------------------------------*
DIM value_int AS INTEGER
DIM MolWt AS DOUBLE

SELECT CASE Label(iMetal)
    CASE "Cu"
        MolWt = 63.546
    CASE "Zn"
        MolWt = 65.38
    CASE "Pb"
        MolWt = 207.2
    CASE ELSE
        ' shouldn't come here
        PRINT "Need to define another metal WQC"
END SELECT

  ' This is the first row that will contain the column labels for the bulk of the data
  iRowOS  = 6
  iRowOS  = 10
  iRowOS  = 11
  IF (Iter = 1) THEN
    '
    '   Write column headings
    '
    IF (iOFormat = THREE) OR (iCriteria <> ZERO) THEN   ' for now, let's always make a WQC file in xls format
       '
       ' Site label and sample label headings
       '
       iRow = iRowOS - 1
       iCol = 1
       sTmp1$ = "Site Label"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       iCol = 2
       sTmp1$ = "Sample Label"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       '
       ' The final acute value (for EPA approved calculations) or the acute HC5 (for non EPA)
       '
       iCol = 3
       ' 1st line has the label
       sTmp1$ = "Acute HC5"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       '
       ' 2nd line has the units
       sTmp1$ = "(AHC5), ug/L"
       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       '
       ' The acute water quality criteria (for EPA approved calculations) or the acute safety factor (i.e., FAV/2) (for non EPA)
       '
       iCol = 4
       sTmp1$ = "Acute Safety Factor"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       ' 2nd line has the units
       sTmp1$ = "(ASF=AHC5/2), ug/L"
       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       '
       ' The chronic water quality criteria (for EPA approved calculations) or leave blank
       '
'       iCol = 5
'       sTmp1$ = "CCC"
'       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
'       ' 2nd line has the units
'       sTmp1$ = "(CCC=FAV/ACR), ug/L"
'       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       '
       ' Metal concentration in ug/L
       '
       iCol = 5
       sTmp1$ = sMetal$
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       ' 2nd line has the units
       sTmp1$ = "ug/L"
       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       '
       ' Acute toxic units = metal/(acute wqc)
       '
       iCol = 6
       sTmp1$ = "Acute Toxic Units"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       ' 2nd line has the units
       'sTmp1$ = "(Acute TU=Cu/CMC)"
       sTmp1$ = "(Acute TU="+sMetal$+"/ASF)"
       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       '
       ' Chronic toxic units = metal/(chronic wqc)
       '
'       iCol = 8
'       sTmp1$ = "Chronic Toxic Units"
'       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
'       ' 2nd line has the units
'       sTmp1$ = "(Chronic TU=Cu/CCC)"
'       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       '
       ' Censored flag for metal concentrations below detection limits
       '
       iCol = 7
       sTmp1$ = "Censored Flag"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
        ' 2nd line has the units
       sTmp1$ = "(0 = quantified, 1 = BDL)"
       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       IF (iCalc_FMB = TRUE) THEN
           '
           ' Acute FMB Labels
           '
           iCol = 4
           iRow = 5
           IF (FMB_option = FMBMedian) THEN
              sTmp1$ = "Acute FMB (median)"
           ELSE
              sTmp1$ = "Acute FMB (geo mean)"
           END IF
           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
           sTmp1$ = "ug/L"
           stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

           '
           '  Chronic FMB Labels
           '
    '       iCol = 5
    '       iRow = 5
    '       IF (FMB_option = FMBMedian) THEN
    '          sTmp1$ = "Chronic (median)"
    '       ELSE
    '          sTmp1$ = "Chronic (geo mean)"
    '       END IF
    '       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
    '       sTmp1$ = "ug/L"
    '       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

           '
           ' Report the exceedance frequency
           '
           iCol = 4
           iRow = 8
           sTmp1$ = TRIM$(STR$(D_EF))
           sTmp1$ = "Based on an exceedance frequency of 1/"+sTmp1$+" days (once in "+TRIM$(FORMAT$(D_EF/365, "###.0#"))+" years)."

           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       END IF 'iCalc_FMB = TRUE
    END IF
 END IF
  '
  '   Write out data
  '
  'FOR j = 1 TO NData
      j = Iter
      IF (iOFormat = THREE) OR (iCriteria <> ZERO) THEN   ' for now, let's always make a WQC file in xls format
         iCol = 1
         iRow = iRowOS + j
         stat& = xlsWriteText(SiteL$, iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
         iCol = iCol + 1
         stat& = xlsWriteText(DataL$, iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)

         '
         '  Final acute value
         value# = HC5(j) * MolWt * 1.E6
         iCol = iCol + 1
         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)

         '
         '  CMC
         value# = HC5(j) * MolWt * 1.E6 / ACUTE_DIV
         WQC_Acute(j) = value#
         iCol = iCol + 1
         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
'         '
'         '  CCC
         value# = HC5(j) * MolWt * 1.E6 / CHRONIC_DIV
         WQC_Chron(j) = value#
'         iCol = iCol + 1
'         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
         '
         '  Ambient metal concentration
         value# = AmbientMetal(j) * MolWt * 1.E6
         MetalConc(j) = value#
         iCol = iCol + 1
         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
         '        '
         '  Acute toxic units
         value# = AmbientMetal(j) / (HC5(j) / ACUTE_DIV)
         TU_Acute(j) = value#
         iCol = iCol + 1
         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
'         '        '
'         '  Chronic toxic units
         value# = AmbientMetal(j) / (HC5(j) / CHRONIC_DIV)
         TU_Chron(j) = value#
'         iCol = iCol + 1
'         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
         '
         '  Censored flag
         value_int = CenFlag(j)
         iCol = iCol + 1
         stat& = xlsWriteInteger(value_int, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)

     END IF
  'NEXT j
  IF (Iter = NData) AND (iCalc_FMB = TRUE) THEN
     '
     ' If this is the last observation, we can now Calculate the FMB
     ' Call the FMB_Acute function
     '
     FMB_Acute = FMB_AcuteFunc(NData&)
     value# = FMB_Acute
     iRow = 7
     iCol = 4
     'stat&c = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
     stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
'     value# = FMB_Chron
'     iRow = 7
'     iCol = 5
'     stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)

      PRINT "TU_EFF_C"        ;   TU_EF
      PRINT "ADJ_F_C"         ;   Adj_Fact
      PRINT "C_Me_Comp"       ;   Me_Compliance
      PRINT
      PRINT "FMB_Acute"    ; FMB_Acute
      PRINT "FMB_Chron"    ; FMB_Chron
      PRINT

  END IF 'Iter = NData AND iCalc_FMB = TRUE
END SUB 'WriteWQCFileAcute

SUB WriteWQCFileChronic(iOutF4txt AS LONG, iOutF4xls AS LONG, Iter AS INTEGER, NData AS LONG)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB WriteWQCFileChronic                                              |
'  |                                                                       |
'  |  This version does chronic WQC and FMB (no acute)                     |
'  |                                                                       |
'  |  Analyzes the output from an entire BLM run and creates an output     |
'  |  file with with chronic instantaneous WQC information, toxic units,   |
'  |  and if metal concentrations are available, it will also calculate    |
'  |  a fixed monitoring benchmark                                         |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |                                                                       |
'  *-----------------------------------------------------------------------*
DIM value_int AS INTEGER
DIM MolWt AS DOUBLE

SELECT CASE Label(iMetal)
    CASE "Cu"
        MolWt = 63.546
    CASE "Zn"
        MolWt = 65.38
    CASE "Pb"
        MolWt = 207.2
    CASE ELSE
        ' shouldn't come here
        PRINT "Need to define another metal WQC"
END SELECT

  ' This is the first row that will contain the column labels for the bulk of the data
  iRowOS  = 6
  iRowOS  = 10
  iRowOS  = 11
  IF (Iter = 1) THEN
    '
    '   Write column headings
    '
    IF (iOFormat = THREE) OR (iCriteria <> ZERO) THEN   ' for now, let's always make a WQC file in xls format

'      xlsWriteDefaultFormats(iOutF4xls)

       iRow = iRowOS - 1
       iCol = 1
       sTmp1$ = "Site Label"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       iCol = 2
       sTmp1$ = "Sample Label"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       iCol = 3
       ' 1st line has the label
       'sTmp1$ = "Final Chronic Value"
       sTmp1$ = "Chronic HC5"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       ' 2nd line has the units
       'sTmp1$ = "(FCV), ug/L"
       sTmp1$ = "(CHC5), ug/L"
       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

'       iCol = 4
'       sTmp1$ = "CMC"
'       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
'       ' 2nd line has the units
'       sTmp1$ = "(CMC=FAV/2), ug/L"
'       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

'       iCol = 4
'       sTmp1$ = "CCC"
'       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
'       ' 2nd line has the units
'       sTmp1$ = "(CCC=FAV/ACR), ug/L"
'       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       iCol = 4
       sTmp1$ = sMetal$
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       ' 2nd line has the units
       sTmp1$ = "ug/L"
       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

'       iCol = 7
'       sTmp1$ = "Acute Toxic Units"
'       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
'       ' 2nd line has the units
'       sTmp1$ = "(Acute TU=Cu/CMC)"
'       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       iCol = 5
       sTmp1$ = "Chronic Toxic Units"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       ' 2nd line has the units
       sTmp1$ = "(Chronic TU="+sMetal$+"/CHC5)"
       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       iCol = 6
       sTmp1$ = "Censored Flag"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
        ' 2nd line has the units
       sTmp1$ = "(0 = quantified, 1 = BDL)"
       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       IF (iCalc_FMB = TRUE) THEN

           ' FMB Labels
    '       iCol = 4
    '       iRow = 5
    '       IF (FMB_option = FMBMedian) THEN
    '          sTmp1$ = "Acute FMB (median)"
    '       ELSE
    '          sTmp1$ = "Acute FMB (geo mean)"
    '       END IF
    '       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
    '       sTmp1$ = "ug/L"
    '       stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

           iCol = 4
           iRow = 5
           IF (FMB_option = FMBMedian) THEN
              sTmp1$ = "Chronic (median)"
           ELSE
              sTmp1$ = "Chronic (geo mean)"
           END IF
           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
           sTmp1$ = "ug/L"
           stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

           iCol = 4
           iRow = 8
           sTmp1$ = TRIM$(STR$(D_EF))
           sTmp1$ = "Based on an exceedance frequency of 1/"+sTmp1$+" days (once in "+TRIM$(FORMAT$(D_EF/365, "###.0#"))+" years)."

           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       END IF 'iCalc_FMB = TRUE
    END IF '(iOFormat = THREE) OR (iCriteria <> ZERO)
 END IF
  '
  '   Write out data
  '
  'FOR j = 1 TO NData
      j = Iter
      IF (iOFormat = THREE) OR (iCriteria <> ZERO) THEN   ' for now, let's always make a WQC file in xls format
         iCol = 1
         iRow = iRowOS + j
         stat& = xlsWriteText(SiteL$, iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
         iCol = iCol + 1
         stat& = xlsWriteText(DataL$, iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)

         '
         '  Final acute value
         value# = HC5(j) * MolWt * 1.E6
         iCol = iCol + 1
         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)

'         '
'         '  CMC
'         value# = HC5(j) * MolWt * 1.E6 / ACUTE_DIV
'         WQC_Acute(j) = value#
'         iCol = iCol + 1
'         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
'         '
'         '  CCC
'         value# = HC5(j) * MolWt * 1.E6 / CHRONIC_DIV
'         WQC_Chron(j) = value#
'         iCol = iCol + 1
'         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
         '
         '  Ambient metal concentration
         value# = AmbientMetal(j) * MolWt * 1.E6
         MetalConc(j) = value#
         iCol = iCol + 1
         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
         '        '
'         '  Acute toxic units
'         value# = AmbientMetal(j) / (HC5(j) / ACUTE_DIV)
'         TU_Acute(j) = value#
'         iCol = iCol + 1
'         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
         '        '
         '  Chronic toxic units
         'value# = AmbientMetal(j) / (HC5(j) / CHRONIC_DIV)
         value# = AmbientMetal(j) / HC5(j)
         TU_Chron(j) = value#
         iCol = iCol + 1
         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
         '
         '  Censored flag
         value_int = CenFlag(j)
         iCol = iCol + 1
         stat& = xlsWriteInteger(value_int, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)

     END IF
  'NEXT j
  IF (Iter = NData) AND (iCalc_FMB = TRUE) THEN
      ' Calculate FMB
      'MLE_Fit(x_dat_in() AS EXTENDED, x_dl() AS LONG, NData AS LONG, iLog AS LONG, x_mle_mean AS EXTENDED, x_mle_std AS EXTENDED)
      'AmbientMetal(i) = OKT(iMetal)
      'CenFlag(i) = OKT(iMetalDL)
      'HC5(j)
      'ToxicUnit(1 TO MaxObs)
      DIM TU_Mean AS EXTENDED
      DIM TU_median AS EXTENDED
      DIM TU_STD AS EXTENDED
      DIM TU_STD_4D AS EXTENDED
      DIM TU_CHR_Mean AS EXTENDED
      DIM TU_CHR_median AS EXTENDED
      DIM TU_CHR_STD AS EXTENDED
      DIM TU_CHR_STD_4D AS EXTENDED
      DIM Me_Mean AS EXTENDED
      DIM Me_median AS EXTENDED
      DIM Me_STD AS EXTENDED
      DIM Me_STD_4D AS EXTENDED
      DIM Me_Compliance AS EXTENDED
      DIM FMB_Acute AS EXTENDED
      DIM FMB_Chron AS EXTENDED
      DIM me_conc_log(1 TO NData) AS EXTENDED
      DIM TU_acut_log(1 TO NData) AS EXTENDED
      DIM TU_chro_log(1 TO NData) AS EXTENDED
      DIM Me_Conc_est(1 TO NData) AS EXTENDED
      DIM TU_acute_est(1 TO NData) AS EXTENDED
      DIM TU_chron_est(1 TO NData) AS EXTENDED
      DIM Me_Conc_est_log(1 TO NData) AS EXTENDED
      DIM TU_acute_est_log(1 TO NData) AS EXTENDED
      DIM TU_chron_est_log(1 TO NData) AS EXTENDED
      DIM x_AVE AS EXTENDED, x_ADEV AS EXTENDED, x_VAR AS EXTENDED, x_SKEW AS EXTENDED, x_CURT AS EXTENDED


      DIM N_Eff AS EXTENDED
      DIM N_Avg AS EXTENDED
      DIM Rho AS EXTENDED
      DIM NCen AS LONG

      DIM iLog AS LONG
      DIM NDataL AS LONG
      NDataL = NData
      iLog = TRUE


      '  ------------------------------------------------------------------
      '
      '  Calculate an FMB
      '
      '  Note - we have duplicate variables defined for the ambient metal conc
      '         MetalConc(i) and MetalFlag(i) and also AmbientMetal(i)
      '         Lets decide which we want to use.
      '
      '  This entire process should be moved into an independant subroutine
      '
      '  The exceedance frequency Z_EF needs to be read in from the parameter
      '         file, or from some other source (command line switch?)
      '         Maybe it could be read in from either location, with the
      '         command line switch over-riding the value in the parameter
      '         file, so that the interface will have a convenient way to
      '         provide a new value. I don't think we have a /Zhh switch yet
      '
      '  1 - use the metal concentrations, considering censored flags with MLE,
      '      and estimate the mean, median, and std dev
      '      using variables Me_Mean, Me_STD, Me_median
      '           - use actual data whenever possible, and only use MLE
      '             estimates for data when original data are BDL
      '  2 - do the same with the TU distribution using the same censored flags
      '      using variables TU_Mean, TU_STD, TU_median
      '  3 - Calculate the TU value at the exceedance frequency Z_EF
      '          TU_EF = 10^(Z_EF * TU_STD + LOG10(TU_median))
      '  4 - Calculate the adjustment factor that will result in a TU value
      '      of one at the exceedance frequency
      '          Adj_Fact = 1 / TU_EF
      '  5 - Calculate the median metal concentration that would correspond
      '      to a TU distribution that is in compliance
      '      Me_Compliance = Me_median * Adj_Fact
      '  6 - Calculate the FMB
      '      FMB = 10^(Z_EF * Me_STD + LOG10(Me_Compliance))
      '
      '  ------------------------------------------------------------------
      '  Perform an MLE fit with the TU values to calculate the degree of
      '  adjustment allowed for a TU distribution exactly in compliance
      '  with the target exceedance frequency
      '  ------------------------------------------------------------------
      NCen = 0
      FOR i=1 TO NData
          IF (CenFlag(i)=1) OR (CenFlag(i)=2) THEN
              NCen = NCen + 1
          END IF
      NEXT i

      IF (NCen > 0) THEN
        CALL MLE_Fit(TU_Acute(), CenFlag(), NDataL, iLog, TU_Mean, TU_STD)
        CALL MLE_Fit(TU_Chron(), CenFlag(), NDataL, iLog, TU_CHR_Mean, TU_CHR_STD)
        CALL MLE_Fit(MetalConc(), CenFlag(), NDataL, iLog, Me_Mean, Me_STD)
      END IF

      DIM iDir AS LONG
      DIM Me_index(1 TO NData) AS LONG
      DIM TU_index(1 TO NData) AS LONG
      'dim Me_Zscore as extended
      'dim TU_Zscore as extended
      DIM x_index(1 TO NData) AS LONG
      'DIM x_zscore(1 TO NData) AS EXTENDED
      'DIM x_dat_mle(1 TO NData) AS EXTENDED
      'DIM x_median AS EXTENDED
      iDir = 0      ' Direction of sort - 0 = ascending, non-zero = descending
      'DIM i AS LONG
      '
      ' Sort once for toxic units
      CALL IndSortExt (TU_Acute(), TU_index(), NDataL, iDir)
      '
      ' and sort again for metal concentrations
      CALL IndSortExt (MetalConc(), Me_index(), NDataL, iDir)

      '
      '  Create an estimated dataset with substituted values from observations
      '  that were originally below detection limit (as determined by MLE fit),
      '  and with actual values whenever possible
      '
      DIM NDataL1 AS LONG      ' copy NData into a LONG variable to pass to ZScore() function
      NDataL1 = NData
      PRINT "------------------  debug  -----------------"
      PRINT "NDataL:"; NDataL
      PRINT "NData:"; NData
      PRINT "CenFlag"; TAB(20); "Me_Conc_est"; TAB(45); "TU_acute_est"; TAB(60); "TU_chron_est"
      FOR i=1 TO NData
          IF (CenFlag(i)=1) OR (CenFlag(i)=2) THEN
              value# = Me_Mean + ZScore(Me_index(i), NDataL1) * Me_STD
              Me_Conc_est(i) = 10^value#

              value# = TU_Mean + ZScore(TU_index(i), NDataL1) * TU_STD
              TU_acute_est(i) = 10^value#

              value# = TU_CHR_Mean + ZScore(TU_index(i), NDataL1) * TU_CHR_STD
              TU_chron_est(i) = 10^value#
          ELSE
              '        tmp1 = x_zscore(i)*x_std_mle + x_median_mle
              Me_Conc_est(i) = MetalConc(i)
              TU_Acute_est(i) = TU_Acute(i)
              TU_Chron_est(i) = TU_Chron(i)
          END IF
          Me_Conc_est_log(i)  = LOG10(Me_Conc_est(i))
          TU_acute_est_log(i) = LOG10(TU_acute_est(i))
          TU_chron_est_log(i) = LOG10(TU_chron_est(i))
          PRINT CenFlag(i); TAB(15); USING$("##.###^^^^",Me_Conc_est(i)); TAB(30); USING$("##.###^^^^",TU_acute_est(i)); TAB(45); USING$("##.###^^^^",TU_chron_est(i))

      NEXT i
      PRINT "---------------- end debug  -----------------"


      ' print out the original and adjusted metal conc
      DIM OMaskT$
      OMaskT$ =  "##.#######^^^^"
      PRINT
      PRINT "Orig_Me"; TAB(20); "Est_Me"; TAB(45); "Cen_Flag"
      FOR i=1 TO NData
          PRINT USING$(OMaskT$, MetalConc(i)); TAB(20); USING$(OMaskT$, Me_Conc_est(i)); TAB(45); CenFlag(i)
      NEXT i


      ' print out the original and adjusted acute toxic units
      PRINT
      PRINT "Orig_ATU"; TAB(20); "Est_ATU"; TAB(45); "Cen_Flag"
      FOR i=1 TO NData
          PRINT USING$(OMaskT$, TU_Acute(i)); TAB(20); USING$(OMaskT$, TU_acute_est(i)); TAB(45); CenFlag(i)
      NEXT i


      ' print out the original and adjusted chronic toxic units
      PRINT
      PRINT "Orig_CTU"; TAB(20); "Est_CTU"; TAB(45); "Cen_Flag"
      FOR i=1 TO NData
          PRINT USING$(OMaskT$, TU_chron(i)); TAB(20); USING$(OMaskT$, TU_chron_est(i)); TAB(45); CenFlag(i)
      NEXT i

      ' Sort again for toxic units using estimates, and this time calc a median
      ' For this - use log transformed values
      CALL IndSortExt (TU_Acute_est_log(), TU_index(), NDataL, iDir)
      CALL Median1(TU_Acute_est_log(), TU_index(), NDataL, TU_median)
      CALL Median1(TU_Chron_est_log(), TU_index(), NDataL, TU_CHR_median)
      TU_median = 10^TU_median
      TU_CHR_median = 10^TU_CHR_median
      '
      ' and do the same for metal concentrations
      ' again, using log transformed values
      CALL IndSortExt (Me_Conc_est_log(), Me_index(), NDataL, iDir)
      CALL Median1(Me_Conc_est_log(), Me_index(), NDataL, Me_median)
      Me_median = 10^Me_median
      '
      '  Now adjust the STD using Log transformed values
      CALL Moment (Me_Conc_est_log(), NDataL, x_AVE, x_ADEV, Me_STD, x_VAR, x_SKEW, x_CURT)
          Me_mean = 10^x_AVE
      CALL Moment (TU_Acute_est_log(), NDataL, x_AVE, x_ADEV, TU_STD, x_VAR, x_SKEW, x_CURT)
          TU_mean = 10^x_AVE
      CALL Moment (TU_Chron_est_log(), NDataL, x_AVE, x_ADEV, TU_CHR_STD  , x_VAR, x_SKEW, x_CURT)
          TU_CHR_mean = 10^x_AVE

      '
      '  Use this latest estimate of the median and STD to
      '  calulate acute FMB
      '
      IF (FMB_option = FMBMedian) THEN
        TU_EF = 10^(Z_EF * TU_STD + LOG10(TU_median))
        Adj_Fact = 1 / TU_EF
        Me_Compliance = Me_median * Adj_Fact
        FMB_Acute = 10^(Z_EF * Me_STD + LOG10(Me_Compliance))
      ELSE
        TU_EF = 10^(Z_EF * TU_STD + LOG10(TU_mean))
        Adj_Fact = 1 / TU_EF
        Me_Compliance = Me_mean * Adj_Fact
        FMB_Acute = 10^(Z_EF * Me_STD + LOG10(Me_Compliance))
      END IF


      PRINT
      PRINT "metal_median" ; Me_median
      PRINT "metal_mean"   ; Me_mean
      PRINT "metal_std"    ; Me_STD
      PRINT
      PRINT "TU_median" ; TU_median
      PRINT "TU_mean"   ; TU_mean
      PRINT "TU_std"    ; TU_STD

      PRINT
      PRINT "Z_EF"           ;  Z_EF
      PRINT "TU_EFF"         ;  TU_EF
      PRINT "ADJ_F"          ;  Adj_Fact
      PRINT "A_Me_Comp"      ;  Me_Compliance

      PRINT
      PRINT "Chron_TU_median" ; TU_CHR_median
      PRINT "Chron_TU_mean"   ; TU_CHR_mean
      PRINT "Chron_TU_std"    ; TU_CHR_STD

      '
      '  Estimate the 4-day avg toxic unit STD
      '
      N_Avg = 4.
      Rho = 0.8
      N_Eff = N_Avg^2*(1-Rho)^2/(N_Avg*(1-Rho^2)-(2*Rho*(1-Rho^N_Avg)))
      Me_STD_4D = (LOG(1+(EXP(Me_STD^2)-1)/N_Eff))^0.5
      TU_CHR_STD_4D = (LOG(1+(EXP(TU_CHR_STD^2)-1)/N_Eff))^0.5
      PRINT
      PRINT "N_Eff"           ;  N_Eff
      PRINT "metal_4D_std"    ;  Me_STD_4D
      PRINT "Chron_TU_4D_std" ;  TU_CHR_STD_4D
      PRINT       '
      '  and use that to calculate a chronic FMB
      '
      IF (FMB_option = FMBMedian) THEN
        TU_EF = 10^(Z_EF * TU_CHR_STD_4D + LOG10(TU_CHR_median))
        Adj_Fact = 1 / TU_EF
        Me_Compliance = Me_median * Adj_Fact
        FMB_Chron = 10^(Z_EF * Me_STD_4D + LOG10(Me_Compliance))
      ELSE
        TU_EF = 10^(Z_EF * TU_CHR_STD_4D + LOG10(TU_CHR_mean))
        Adj_Fact = 1 / TU_EF
        Me_Compliance = Me_mean * Adj_Fact
        FMB_Chron = 10^(Z_EF * Me_STD_4D + LOG10(Me_Compliance))
      END IF

'     value# = FMB_Acute
'     iRow = 7
'     iCol = 4
'     'stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
'     stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
     value# = FMB_Chron
     iRow = 7
     iCol = 4
     stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)

      PRINT "TU_EFF_C"        ;   TU_EF
      PRINT "ADJ_F_C"         ;   Adj_Fact
      PRINT "C_Me_Comp"       ;   Me_Compliance
      PRINT
      PRINT "FMB_Acute"    ; FMB_Acute
      PRINT "FMB_Chron"    ; FMB_Chron
      PRINT


  END IF 'Iter = NData AND iCalc_FMB = TRUE
END SUB 'WriteWQCFileChronic

SUB WriteWQC_Acute_and_Chronic(iOutF4txt AS LONG, iOutF4xls AS LONG, Iter AS INTEGER, NData AS LONG)
'
'  *-----------------------------------------------------------------------*
'  |                                                                       |
'  |  SUB WriteWQCFile                                                     |
'  |                                                                       |
'  |  This version does acute and chronic WQC and FMB                      |
'  |                                                                       |
'  |  Analyzes the output from an entire BLM run and creates an output     |
'  |  file with with acute and chronic instantaneous WQC information,      |
'  |  toxic units, and if metal concentrations are available, it will      |
'  |  also calculate a fixed monitoring benchmark                          |
'  |                                                                       |
'  |  Robert Santore                                                       |
'  |                                                                       |
'  *-----------------------------------------------------------------------*
DIM value_int AS INTEGER
DIM MolWt AS DOUBLE

SELECT CASE Label(iMetal)
    CASE "Cu"
        MolWt = 63.546
    CASE "Zn"
        MolWt = 65.38
    CASE "Pb"
        MolWt = 207.2
    CASE "Ni"
       MolWt = 58.693
    CASE "Cd"
       MolWt = 112.414
    CASE "Al"
       MolWt = 26.982
    CASE "Ag"
       MolWt = 107.868
    CASE ELSE
        ' shouldn't come here
        PRINT "Need to define another metal WQC"
END SELECT

' This is the first row that will contain the column labels for the bulk of the data
iRowOS  = 6
iRowOS  = 10
iRowOS  = 11

  IF (Iter = 1) THEN
    '
    '   Write column headings
    '
    IF (iOFormat = THREE) OR (iCriteria <> ZERO) THEN   ' for now, let's always make a WQC file in xls format

       iRow = iRowOS - 1
       iCol = 1
       sTmp1$ = "Site Label"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       iCol = 2
       sTmp1$ = "Sample Label"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       ' Column three is the BLM-predicted metal, label depends on the mode
       iCol = 3
       SELECT CASE iCriteria
          CASE ONE
             sTmp1$ = "Final Acute Value"
             sTmp2$ = "(FAV), ug/L"
          CASE TWO
             sTmp1$ = "Acute HC5"
             sTMp2$ = "(AHC5), ug/L"
          CASE THREE
             sTmp1$ = "Chronic HC5"
             sTMp2$ = "(CHC5), ug/L"
       END SELECT
       ' 1st line has the label
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       ' 2nd line has the units
       stat& = xlsWriteText(sTmp2$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       ' Next column has CMC / ASF, nothing for chronic only
       IF (iCriteria <> THREE) THEN
          iCol = iCol + 1
          SELECT CASE iCriteria
             CASE ONE
                sTmp1$ = "CMC"
                sTmp2$ = "(CMC=FAV/" + STR$(ACUTE_DIV) + "), ug/L"
             CASE TWO
                sTmp1$ = "Acute Safety Factor"
                sTmp2$ = "(ASF=AHC5/" + STR$(ACUTE_DIV) + "), ug/L"
          END SELECT
          ' 1st line has the label
          stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
          ' 2nd line has the units
          stat& = xlsWriteText(sTmp2$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       END IF

       ' next column, only for A+C run
       IF (iCriteria = ONE) THEN
          iCol = iCol + 1
          sTmp1$ = "CCC"
          stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
          ' 2nd line has the units
          sTmp2$ = "(CCC=FAV/ACR), ug/L"
          stat& = xlsWriteText(sTmp2$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       END IF

       ' Next col always has the metal
       iCol = iCol + 1
       sTmp1$ = sMetal$
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       ' 2nd line has the units
       sTmp2$ = "ug/L"
       stat& = xlsWriteText(sTmp2$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)

       ' Acute TU's when it's run
       IF (iCriteria <> THREE) THEN
          iCol = iCol + 1
          sTmp1$ = "Acute Toxic Units"
          stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
          ' 2nd line has the units
          SELECT CASE iCriteria
             CASE ONE
                sTmp2$ = "(Acute TU=" + sMetal$ + "/CMC)"
             CASE TWO
                sTmp2$ = "(Acute TU=" + sMetal$ + "/ASF)"
          END SELECT
          stat& = xlsWriteText(sTmp2$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       END IF

       ' Chronic TU's when it's run
       IF (iCriteria <> TWO) THEN
          iCol = iCol + 1
          sTmp1$ = "Chronic Toxic Units"
          stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
          ' 2nd line has the units
          SELECT CASE iCriteria
             CASE ONE
                sTmp2$ = "(Chronic TU=" + sMetal$ + "/CCC)"
             CASE THREE
                sTmp2$ = "(Chronic TU=" + sMetal$ + "/CHC5)"
          END SELECT
          stat& = xlsWriteText(sTmp2$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
       END IF

       ' Censored flag always
       iCol = iCol + 1
       sTmp1$ = "Censored Flag"
       stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
        ' 2nd line has the units
       sTmp2$ = "(0 = quantified, 1 = BDL)"
       stat& = xlsWriteText(sTmp2$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)


       IF (iCalc_FMB = TRUE) THEN
           '
           ' FMB Labels
           '
           iCol = 3
           IF (iCriteria <> THREE) THEN
              iCol = iCol + 1
              iRow = 5
              IF (FMB_option = FMBMedian) THEN
                 sTmp1$ = "Acute FMB (median)"
              ELSE
                 sTmp1$ = "Acute FMB (geo mean)"
              END IF
              stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
              sTmp1$ = "ug/L"
              stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
           END IF

           IF (iCriteria <> TWO) THEN
              iCol = iCol + 1
              iRow = 5
              IF (FMB_option = FMBMedian) THEN
                 sTmp1$ = "Chronic (median)"
              ELSE
                 sTmp1$ = "Chronic (geo mean)"
              END IF
              stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
              sTmp1$ = "ug/L"
              stat& = xlsWriteText(sTmp1$, iRow+1, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
           END IF

           iCol = 4
           iRow = 8
           sTmp1$ = TRIM$(STR$(D_EF))
           sTmp1$ = "Based on an exceedance frequency of 1/"+sTmp1$+" days (once in "+TRIM$(FORMAT$(D_EF/365, "###.0#"))+" years)."

           stat& = xlsWriteText(sTmp1$, iRow, iCol, %xlsFont1, %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
        END IF 'iCalc_FMB = TRUE
    END IF 'IF (iOFormat = THREE) OR (iCriteria <> ZERO)
 END IF 'IF (ITER = 1)
  '
  '   Write out data
  '
  'FOR j = 1 TO NData
      j = Iter
      IF (iOFormat = THREE) OR (iCriteria <> ZERO) THEN   ' for now, let's always make a WQC file in xls format
         iCol = 1
         iRow = iRowOS + j
         stat& = xlsWriteText(SiteL$, iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)

         iCol = iCol + 1
         stat& = xlsWriteText(DataL$, iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)

         '
         '  Final acute value
         iCol = iCol + 1
         value# = HC5(j) * MolWt * 1.E6
         IF iComplete AND NOT ErrorNotConverged THEN
            stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
         ELSE
            stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
         END IF

         '
         '  CMC
         IF (iCriteria <> THREE) THEN
            iCol = iCol + 1
            IF iComplete AND NOT ErrorNotConverged THEN
               value# = HC5(j) * MolWt * 1.E6 / ACUTE_DIV
               WQC_Acute(j) = value#
               stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
            ELSE
               WQC_Acute(j) = NULL
               stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
            END IF
         END IF
         '
         '  CCC
         IF (iCriteria = ONE) THEN
            iCol = iCol + 1
            IF iComplete AND NOT ErrorNotConverged THEN
               value# = HC5(j) * MolWt * 1.E6 / CHRONIC_DIV
               WQC_Chron(j) = value#
               stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
            ELSE
               WQC_Chron(j) = NULL
               stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
            END IF
         ELSEIF (iCriteria = THREE) THEN
            WQC_Chron(j) = HC5(j) * MolWt * 1.E6
         END IF
         '
         '  Ambient metal concentration
         value# = AmbientMetal(j) * MolWt * 1.E6
         MetalConc(j) = value#
         iCol = iCol + 1
         stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
         '
         '  Acute toxic units
         IF (iCriteria <> THREE) THEN
            iCol = iCol + 1
            IF iComplete AND NOT ErrorNotConverged THEN
               value# = MetalConc(j) / WQC_Acute(j)
               TU_Acute(j) = value#
               stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
            ELSE
               TU_Acute(j) = NULL
               stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
            END IF
         END IF
         '
         '  Chronic toxic units
         IF  (iCriteria <> TWO) THEN
            iCol = iCol + 1
            IF iComplete AND NOT ErrorNotConverged THEN
               value# = MetalConc(j) / WQC_Chron(j)
               TU_Chron(j) = value#
               stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
            ELSE
               TU_Chron(j) = NULL
               stat& = xlsWriteText("NA", iRow, iCol, %xlsFont0, %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
            END IF
         END IF
         '
         '  Censored flag
         value_int = CenFlag(j)
         iCol = iCol + 1
         stat& = xlsWriteInteger(value_int, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)

     END IF 'IF (iOFormat = THREE) OR (iCriteria <> ZERO)
  'NEXT j


     IF (Iter = NData) AND (iCalc_FMB = TRUE) THEN

       'Remove NULL values
       DIM AmbientMetal_tmp (1 TO NData) AS DOUBLE
       DIM HC5_tmp (1 TO NData) AS DOUBLE
       DIM CenFlag_tmp (1 TO NData) AS LONG
       DIM TU_Acute_tmp (1 TO NData) AS EXTENDED
       DIM TU_Chron_tmp (1 TO NData) AS EXTENDED
       DIM WQC_Acute_tmp (1 TO NData) AS EXTENDED
       DIM WQC_Chron_tmp (1 TO NData) AS EXTENDED
       DIM MetalConc_tmp (1 TO NData) AS EXTENDED
       FOR i = 1 TO NData
          AmbientMetal_tmp(i) = AmbientMetal(i)
          HC5_tmp(i) = HC5(i)
          CenFlag_tmp(i) = CenFlag(i)
          TU_Acute_tmp(i) = TU_Acute(i)
          TU_Chron_tmp(i) = TU_Chron(i)
          WQC_Acute_tmp(i) = WQC_Acute(i)
          WQC_Chron_tmp(i) = WQC_Chron(i)
          MetalConc_tmp(i) = MetalConc(i)
       NEXT i
       ii = 0
       FOR i = 1 TO NData
          AmbientMetal(i) = NULL
          HC5(i) = NULL
          CenFlag(i) = NULL
          TU_Acute(i) = NULL
          TU_Chron(i) = NULL
          WQC_Acute(i) = NULL
          WQC_Chron(i) = NULL
          MetalConc(i) = NULL
          IF TU_Acute_tmp(i) = NULL THEN
             PRINT STR$(i)
          ELSE
             ii = ii + 1
             AmbientMetal(ii) = AmbientMetal_tmp(i)
             HC5(ii) = HC5_tmp(i)
             CenFlag(ii) = CenFlag_tmp(i)
             TU_Acute(ii) = TU_Acute_tmp(i)
             TU_Chron(ii) = TU_Chron_tmp(i)
             WQC_Acute(ii) = WQC_Acute_tmp(i)
             WQC_Chron(ii) = WQC_Chron_tmp(i)
             MetalConc(ii) = MetalConc_tmp(i)
          END IF
       NEXT i
       NData = ii

       CALL FMBCalc(NData, FMB_Acute, FMB_Chron)

       iCol = 3
       IF(iCriteria <> THREE) THEN
           value# = FMB_Acute
           iRow = 7
           iCol = iCol + 1
           'stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 0, iOutF4xls)
           stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
       END IF

       IF (iCriteria <> TWO) THEN
        value# = FMB_Chron
        iRow = 7
        iCol = iCol + 1
        stat& = xlsWriteNumber(value#, iRow, iCol, %xlsFont0,  %xlsLeftAlign, %xlsCellNormal, 2, iOutF4xls)
       END IF

     END IF '(Iter = NData) AND (iCalc_FMB = TRUE)

END SUB


SUB PrintHelp
   TmpS$ = MID$(VersionID$, 1, 1) + MID$(VersionID$, 3, 2)
   PRINT
   PRINT "Usage: BLM_" + TmpS$ + "  <parameter file name>  <input file name>  [Options]"
   PRINT "Alternate Usage: BLM_" + TmpS$ + " /S <script file name>"
   PRINT "     where the script file has <parameter file name>,  <input file name>, and
   PRINT "     and [Options] specified in the first 3 lines."
   PRINT
   PRINT "Options:"
   PRINT "/D - Create a file with derivative information"
   PRINT "/L - Run all observations in metal toxicity mode"
   PRINT "/E - Read LA50 values from an extra column in the input file"
   PRINT "/E2 - Read LA50 and parameter file info from 2 extra columns in the input file"
   'PRINT "/E3 - (for Al toxicity) read DISS_A and PART_HI from 2 extra columns in the input file"
   PRINT "/Q - Quiet mode (no diagnostic information provided) - somewhat faster!"
   PRINT "/QQ - Very quiet mode"
   PRINT "/V - Print the version and build information for this copy of the BLM"
   PRINT "/VER(x+) - Specify the version information to appear in output files with (x+)"
   PRINT "/W - Use the input file format created by the MS Windows BLM version"
   PRINT "/F - Write out separate information for humic and fulvic binding."
   PRINT "/M - Use a lower mass BL site density (useful for marine)."
   PRINT "/I - Echo inputs (don't run the model, just read the input file and "
   PRINT "     return what was read - used for de-bugging)."
   PRINT "/CAxx - Obtain critical accumulation value from row xx of a comma-delimited"
   PRINT "        table in the parameter file."
   PRINT "/CO2(x+) -  Make the carbonate system an open system, with a pCO2 of whatever "
   PRINT "            is specified in (x+), or 3.5 if not specified."
   PRINT
   PRINT "Options to select convergence algorithm:"
   PRINT " /A1 - uses older algorithm used by BLM_A008, robust but slow."
   PRINT " /A2 - should be much faster than A1 for fish parameter files."
   PRINT " /A3 - (default), and /A4 should be faster than A1 or A2."
   PRINT
   PRINT "Options to select convergence criteria:"
   PRINT " /C1 - (default) specifies a convergence criteria for the function:
   PRINT "          'ABS(BLMConc! - CritBL) / CritBL' equal to 0.001."
   PRINT " /C2 - reduces the convergence criteria to 0.0001."
   PRINT " /C3 - reduces the convergence criteria to 0.00001."
   PRINT " /C4 - reduces the convergence criteria to 0.000001."
   PRINT " /C5 - reduces the convergence criteria to 0.0000001."
   PRINT
   PRINT "Options to select output format:"
   PRINT " /O1 - Specifies a space-delimited text file as the output file format."
   PRINT "       (default)"
   PRINT " /O2 - Specifies a comma-delimited text file as the output file format."
   'PRINT "       (NOTE /O2 is not yet implemented and currently reverts to /O1)."
   PRINT " /O3 - Specifies an MS Excel file as the output file format."
   PRINT
   PRINT "Options for calculating LA50/LC50 output:"
   PRINT " /B1 - Calcuate a separate LA50 for each observation"
   PRINT " /B2 - Calculate as log value"
   PRINT " /B3 - For Al toxicity, output values for both dissolved and precipitated"
   PRINT "       mechanisms"
   PRINT
   PRINT "Options for Solubility of Aluminum (Al), Iron (Fe), or Lead (Pb):"
   PRINT " /Al0, /Fe0, or /Pb0 - solubility is not considered."
   PRINT " /Al1, /Fe1, or /Pb1 - precipitate only if over-saturated."
   PRINT " /Al2, /Fe2, or /Pb2 - precip or dissolve all the time."
   PRINT " /KSAlxxxxx - Aluminum solubility constant, where 'xxxxx' is the log of the"
   PRINT "              solubility constant."
   PRINT
   PRINT "Robert Santore, RobertS@WindwardEnv.com, BLM " + VersionID$
END SUB

SUB DebugLog (Location$, DebugInfo$)
   PRINT #iDBOut&, FormatStr$(Location$, "/  /");
   PRINT #iDBOut&, USING$("######  ",iObs);
   PRINT #iDBOut&, DebugInfo$
END SUB

#INCLUDE "CHESS.bas"
