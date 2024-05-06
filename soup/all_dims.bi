'
'  START CONSTANT LIST
'
' *-------------------------------------------------*
' |  Expressions used for Boolean type comparisons  |
' *-------------------------------------------------*
DIM FALSE AS GLOBAL INTEGER
DIM TRUE AS GLOBAL INTEGER
 FALSE = 0
 TRUE = NOT FALSE

' *----------------------------*
' |  Simple Numeric constants  |
' *----------------------------*
DIM ZERO  AS GLOBAL INTEGER
DIM ZeroS AS GLOBAL SINGLE
DIM ZeroD  AS GLOBAL DOUBLE
DIM ONE  AS GLOBAL INTEGER
DIM TWO  AS GLOBAL INTEGER
DIM THREE  AS GLOBAL INTEGER
DIM FOUR  AS GLOBAL INTEGER
DIM FIVE  AS GLOBAL INTEGER
DIM TEN  AS GLOBAL INTEGER
 ZERO = 0
 ZeroS = 0!
 ZeroD = 0#
 ONE = 1
 TWO = 2
 THREE = 3
 FOUR = 4
 FIVE = 5
 TEN = 10

' *-------------------------*
' |  String Initialization  |
' *-------------------------*
DIM NullStr AS GLOBAL STRING
DIM Comma AS GLOBAL STRING
DIM Quote AS GLOBAL STRING
 Comma$ = ","
 Quote$ = CHR$(34)
 NullStr = ""


'
' Common Integers
'
' File Handles
DIM iInFile AS GLOBAL INTEGER
DIM iOutF1 AS GLOBAL LONG 'INTEGER
DIM iOutF2 AS GLOBAL LONG ' INTEGER
DIM iOutF3 AS GLOBAL INTEGER
DIM iOutF4 AS GLOBAL LONG
DIM iWHAMDBS AS GLOBAL INTEGER
DIM iBLM AS GLOBAL INTEGER
DIM iMsgF AS GLOBAL INTEGER                 ' file handle for blm.message file

DIM iHA1 AS GLOBAL INTEGER
DIM iHA2 AS GLOBAL INTEGER
DIM iHA3 AS GLOBAL INTEGER
DIM iHA4 AS GLOBAL INTEGER
DIM iHB1 AS GLOBAL INTEGER
DIM iHB2 AS GLOBAL INTEGER
DIM iHB3 AS GLOBAL INTEGER
DIM iHB4 AS GLOBAL INTEGER
DIM iFA1 AS GLOBAL INTEGER
DIM iFA2 AS GLOBAL INTEGER
DIM iFA3 AS GLOBAL INTEGER
DIM iFA4 AS GLOBAL INTEGER
DIM iFB1 AS GLOBAL INTEGER
DIM iFB2 AS GLOBAL INTEGER
DIM iFB3 AS GLOBAL INTEGER
DIM iFB4 AS GLOBAL INTEGER
DIM MajorVer AS GLOBAL INTEGER     '  Variable for major portion of verion number
DIM MinorVer AS GLOBAL INTEGER     '  Variable for minor portion of verion number

'
' Common Strings
'
' File Name
DIM sInFile AS GLOBAL STRING
DIM sOutF1 AS GLOBAL STRING
DIM sOutF2 AS GLOBAL STRING
DIM sBLM AS GLOBAL STRING
DIM sCHESS AS GLOBAL STRING
DIM DataL AS GLOBAL STRING
DIM SiteL AS GLOBAL STRING
'
'  Debugging info
'
DIM iDebug AS GLOBAL INTEGER
DIM iDBOut AS GLOBAL LONG
DIM iObs AS GLOBAL INTEGER
'

' *-------------------------------------------------------------*
' |  GlobalVars common block                                    |
' |  These are variables that other programs will need access   |
' |  if they wish to call the soil chemistry module             |
' *-------------------------------------------------------------*

' *-------------------------------------------------------------*
' |  Simple variables                                           |
' *-------------------------------------------------------------*
DIM NComp AS GLOBAL INTEGER
DIM NPhase AS GLOBAL INTEGER
DIM NSpecies AS GLOBAL INTEGER
DIM NList AS GLOBAL INTEGER
DIM MaxComp AS GLOBAL INTEGER
DIM MaxPhases AS GLOBAL INTEGER
DIM MaxSpecies AS GLOBAL INTEGER
DIM MaxLists AS GLOBAL INTEGER
DIM SysTemp AS GLOBAL SINGLE
DIM SimType AS GLOBAL INTEGER
DIM MaxIter AS GLOBAL INTEGER
DIM MaxECxIter AS GLOBAL INTEGER
DIM SepNOM AS GLOBAL INTEGER
DIM NMass AS GLOBAL INTEGER   ' new for ver 3.10

' *-------------------------------------------------------------*
' |  Error codes for use with message files                     |
' *-------------------------------------------------------------*
DIM iEC_Critical            AS GLOBAL LONG
DIM iEC_NonCritical_NoPost  AS GLOBAL LONG
DIM iEC_NonCritical_PostOK  AS GLOBAL LONG
    iEC_Critical           = 1
    iEC_NonCritical_NoPost = 2
    iEC_NonCritical_PostOK = 3

' *-------------------------------------------------------------*
' |  Array variables                                            |
' *-------------------------------------------------------------*
'dim CtoM() AS GLOBAL SINGLE
'dim CType() AS INTEGER
'dim EC() AS SINGLE
'dim IC() AS SINGLE
'dim KT() AS DOUBLE
'dim Label() AS STRING
'dim OC() AS DOUBLE
'dim OldKT() AS DOUBLE
'dim PType() AS INTEGER
'dim R() AS DOUBLE
'dim SCa() AS SINGLE
'dim SCb() AS SINGLE
'dim SEC() AS SINGLE
'dim SiteDen() AS SINGLE
'dim SType() AS INTEGER
'dim tEC() AS SINGLE
'dim tSCa() AS SINGLE
'dim tSCb() AS SINGLE
'dim V() AS SINGLE



' *-------------------------------------------------------------*
' |  SoilVars common block                                      |
' |                                                             |
' |  These variables are not likely to be needed by other       |
' |  programs that need to call the soil chemistry module       |
' *-------------------------------------------------------------*

' *-------------------------------------------------------------*
' |  Simple variables                                           |
' *-------------------------------------------------------------*
DIM AqueousCtoM AS GLOBAL SINGLE
DIM BulkDensity AS GLOBAL SINGLE
DIM CatchArea AS GLOBAL SINGLE
DIM ErrorTol AS GLOBAL SINGLE
DIM FileType AS GLOBAL INTEGER
DIM InputFile AS GLOBAL STRING
DIM InputFileType AS GLOBAL INTEGER
DIM Mass AS GLOBAL SINGLE
DIM NFirstS AS GLOBAL INTEGER
DIM NLayers AS GLOBAL INTEGER
DIM OutputFile AS GLOBAL STRING
DIM OutputFileType AS GLOBAL INTEGER
DIM PrecipRate AS GLOBAL SINGLE
DIM SoilCol AS GLOBAL INTEGER
DIM Surface AS GLOBAL INTEGER
DIM SurfaceCtoM AS GLOBAL SINGLE
DIM TimeStep AS GLOBAL SINGLE
DIM WaterVol AS GLOBAL SINGLE
DIM ErrorNotConverged AS GLOBAL INTEGER

' *-------------------------------------------------------------*
' |  Simple variables related to aluminum solubility            |
' *-------------------------------------------------------------*
DIM LKSOAL25 AS GLOBAL DOUBLE
DIM DHAL AS GLOBAL DOUBLE
DIM LKSOAL AS GLOBAL DOUBLE
DIM LIAPAL  AS GLOBAL DOUBLE
DIM iALPPT AS GLOBAL INTEGER
DIM iGlobalAl AS GLOBAL INTEGER
DIM LogKsAl AS GLOBAL DOUBLE
DIM iLogKsAl AS GLOBAL INTEGER
DIM SatInd AS GLOBAL DOUBLE
DIM iFGDone AS GLOBAL INTEGER
DIM iAl AS GLOBAL INTEGER
DIM iAlp AS GLOBAL INTEGER
DIM SmallNumber AS GLOBAL DOUBLE
    SmallNumber# = 1E-20

' *---------------------------------------------------------------------------------------*
' |  Simple variables related to dissolved and particulate aluminum toxicity              |
' *---------------------------------------------------------------------------------------*
'   these new parameters specified in parameter file
'       [DISS_A]:  3.08882
'       [DISS_B]: -1.09887
'
'       [PART_A]:  8.5368
'       [PART_B]: -3.0622
'       [PART_H]:  0.025
'
'       [TARGET_R]: 0.1
DIM DISS_A AS GLOBAL DOUBLE        ' intercept of response curve based on BL-metal
DIM DISS_B AS GLOBAL DOUBLE        ' slope of response curve based on BL-metal
DIM PART_A AS GLOBAL DOUBLE        ' intercept of response curve based on particulate metal
DIM PART_B AS GLOBAL DOUBLE        ' slope of response curve based on particulate metal
DIM PART_H AS GLOBAL DOUBLE        ' hardness dependence of particulate response
DIM PART_HI AS GLOBAL DOUBLE       ' hardness intercept for the relation PART_A = H * PART_HM +  PART_HI
DIM PART_HM AS GLOBAL DOUBLE       ' hardness slope for the relation PART_A = H * PART_HM +  PART_HI
DIM TARGET_R AS GLOBAL DOUBLE      ' target reduction in biological effect due to toxicity (as fraction of control)
DIM DISS_A_expert AS GLOBAL DOUBLE   ' value of DISS_A specified in expert mode
DIM PART_HI_expert AS GLOBAL DOUBLE  ' value of PART_HI specified in expert mode
DIM TARGET_R_expert AS GLOBAL DOUBLE ' value of TARGET_R specified in expert mode
DIM LogKsAl_expert AS GLOBAL DOUBLE  ' value of LogKsAl specified in expert mode

'
'   These are calculated internally only if iPartMe = True
'
DIM iPartMe AS GLOBAL INTEGER
DIM DISS_Effect AS GLOBAL DOUBLE   ' portion of target effect acheivable by dissolved metal
DIM PART_Effect AS GLOBAL DOUBLE   ' portion of target effect acheivable by particulate meta
DIM Part_Conc AS GLOBAL DOUBLE     ' particulate metal concentration needed if Diss_Effect < TARGET_R
DIM Surf_Charge AS GLOBAL DOUBLE   ' surface charge on the precipitated Al and adsorbed DOC
DIM Inorg_Surf_Charge AS GLOBAL DOUBLE
DIM Surf_Charge_Al AS GLOBAL DOUBLE
DIM Surf_Charge_Al_i AS GLOBAL DOUBLE
DIM Hardness_mgL AS GLOBAL DOUBLE
DIM DOC_mgL AS GLOBAL DOUBLE
DIM DOC_Charge_eqL AS GLOBAL DOUBLE
DIM DOC_Surf_Charge AS GLOBAL DOUBLE



' *-------------------------------------------------------------*
' |  Simple variables related to iron solubility                |
' *-------------------------------------------------------------*
DIM LKSOFE25 AS GLOBAL DOUBLE
DIM DHFE AS GLOBAL DOUBLE
DIM LKSOFE AS GLOBAL DOUBLE
DIM LIAPFE AS GLOBAL DOUBLE
DIM iFePPT AS GLOBAL INTEGER
DIM iGlobalFe AS GLOBAL INTEGER


' *---------------------------------------------------------*
' |  Simple variables related to lead solubility            |
' *---------------------------------------------------------*
DIM LKSOPB25 AS GLOBAL DOUBLE
DIM DHPB AS GLOBAL DOUBLE
DIM LKSOPB AS GLOBAL DOUBLE
DIM LIAPPB  AS GLOBAL DOUBLE
DIM iPbPPT AS GLOBAL INTEGER
DIM iGlobalPb AS GLOBAL INTEGER


' *-------------------------------------------------------------*
' |  Simple variables related to CO2(g) solubility              |
' *-------------------------------------------------------------*
DIM iGlobalCO2 AS GLOBAL INTEGER
DIM iCO2PPT AS GLOBAL INTEGER
DIM LKH_CO2 AS GLOBAL DOUBLE         ' Log Henry's Law constant
DIM pPCO2 AS GLOBAL DOUBLE            ' negative log of the partial pressure of CO2
DIM Log_CO3 AS GLOBAL DOUBLE         ' Log of the CO3[2-] ion concentration
DIM Log_K_HCO3 AS GLOBAL DOUBLE     ' Log of the formation constant for H2CO3

' *-------------------------------------------------------------*
' |  Array variables                                            |
' *-------------------------------------------------------------*
'DIM ActCorr() AS INTEGER
'DIM ActCoef() AS DOUBLE
'DIM Activity() AS STRING
'DIM ChangeMass() AS DOUBLE
'DIM DeltaH() AS SINGLE
'DIM ECTemp() AS SINGLE
'DIM hEC() AS SINGLE
'DIM hSEC() AS SINGLE
'DIM InitialConc() AS SINGLE
'DIM LinkComp() AS INTEGER
'DIM LinkList() AS INTEGER
'DIM ListPos() AS INTEGER
'DIM PhaseMass()  AS SINGLE
'DIM PhaseOrder() AS INTEGER
'DIM PSubOrder() AS INTEGER
'DIM SComp() AS INTEGER
'DIM SDeltaH() AS SINGLE
'DIM SECTemp() AS SINGLE
'DIM SPhase() AS INTEGER
'DIM SP() AS SINGLE
'DIM STempCorr() AS INTEGER
'DIM TempCorr() AS INTEGER
'DIM tSEC() AS SINGLE
'DIM tSP() AS SINGLE
'DIM Z() AS DOUBLE


DIM OFC AS GLOBAL OFContents

' *--------------------------------------------------------------------*
' |  These next set of variables are to be set in SUB SimplifyProblem  |
' *--------------------------------------------------------------------*
DIM NSolve AS GLOBAL INTEGER
DIM NFixed AS GLOBAL INTEGER
DIM NUsed AS GLOBAL INTEGER
DIM NotUsed AS GLOBAL INTEGER
'DIM CompOrder() AS INTEGER
'DIM CompList() AS INTEGER
'DIM SpecList() AS INTEGER

' *---------------------------------------------------------------*
' |  The following variables are useful in Monte Carlo scenarios  |
' *---------------------------------------------------------------*
'DIM MonteList() AS INTEGER  ' List of constants to be used in Monte Carlo
'DIM MeanLogK() AS SINGLE       ' Mean of Monte Carlo constants
'DIM VarLogK() AS SINGLE        ' Var of Monte Carlo constants
DIM NSims   AS GLOBAL INTEGER      ' Number of Mon. Carlo iter
DIM NMonteVars  AS GLOBAL INTEGER  ' Number of variables to adj.


' *--------------------------*
' |  WHAM Related Variables  |
' *--------------------------*
DIM iWHAM AS GLOBAL INTEGER ' Flag if true then incl WHAM
DIM iBL AS GLOBAL INTEGER ' position of BL in comp list
'DIM iBLM() AS INTEGER ' position of BL-Metal in spec list
'DIM iBLMTox() AS INTEGER ' toxicity of BL-Metal in spec list (true or false)
DIM iMetal AS GLOBAL INTEGER ' position of metal in comp list
DIM iDOC   AS GLOBAL INTEGER ' position of DOC in comp list
DIM PercHA AS GLOBAL SINGLE  ' Percent of DOC that is HA
DIM SHA AS GLOBAL SINGLE     ' Moles strong sites per g HA carbon
DIM SFA AS GLOBAL SINGLE     ' Moles strong sites per g FA carbon
DIM iSHA AS GLOBAL INTEGER   ' Position of HA strong sites
DIM iSFA AS GLOBAL INTEGER   ' Position of FA strong sites
DIM NSP    AS GLOBAL INTEGER ' Number of WHAM species
DIM sMetBL(10) AS GLOBAL STRING
DIM iBLM(10) AS GLOBAL INTEGER
DIM iBLMTox(10) AS GLOBAL INTEGER
DIM WHAMID(25) AS GLOBAL INTEGER       ' WHAM database ID of needed components

'DIM wOrgBound!()      ' Organically bound concentrations

'DIM wN%()            ' numerical identification of species
'DIM wN$()            ' nominal identification of species
'DIM wT!()            ' total concentrations of components
'DIM wTH!()           ' total concentrations of humics, fulvics
'DIM wA!()            ' activities of solution species
'DIM wC!()            ' concentations of solution species
'DIM wCH%()           ' charges on solution species
'DIM wM1%()           ' master species identifiers
'DIM wM2%()           ' master species identifiers
'DIM wM3%()           ' master species identifiers
'DIM wS1%()           ' stoichiometries
'DIM wS2%()           ' stoichiometries
'DIM wS3%()           ' stoichiometries
'DIM wLK!()           ' equilibrium constants
'DIM wDH!()           ' enthalpies
'DIM wGAMMA!()        ' activity coeficients
'
'DIM wSP!()           ' proton-binding sites not forming bidentate site
'DIM wSM!()           ' bidentate metal-binding sites
'DIM wKH!()           ' K's for proton-binding
'DIM wPKH!()          ' PK's for proton-binding
'DIM wKMH!()          ' K's for metal-proton exchange
'DIM wSITE%()         ' proton sites making bidentate sites
'DIM wPKMHA!()        ' metal-proton exchange constants
'DIM wPKMHB!()        ' metal-proton exchange constants
'DIM wFP!()           ' disociation factors for proton-binding
'DIM wBIBTERM!()      ' terms calc'd in finding bidentate theta's
'DIM wMONBTERM!()     ' terms calc'd in finding monodentate theta's
'DIM wBISNU!()        ' values of NU at each bidentate site
'DIM wMONSNU!()       ' values of NU at each monodentate site
'DIM wBINU!()         ' overal NU for each metal/bidentate sites
'DIM wMONNU!()        ' overal NU for each metal/monodentate sites
'DIM wNU!()           ' overal NU for each complexed species
'DIM wCHC!()          ' conc of complexed species per litre
'DIM wCDDL!()         ' concs in HS DDL's
'DIM wBIZ!()          ' net charge at each bidentate site
'DIM wMONZ!()         ' net charge at each monodentate site
'DIM wBINDTEST%()     ' alows binding of a species to be tested
'
'DIM wTCALC!()        ' total calculted concentrations of components
'DIM wRADIUS!()       ' hydrated radius of humics, fulvics
'DIM wMOLWT!()        ' molecular weight of humics, fulvics
'DIM DDLVOL!()        ' Diffuse layer volume of humics, fulvics
'DIM DVOLMAX!()       ' Max. diffuse layer volume of humics, fulvics
'DIM ZED!()           ' Charge of humics, fulvics
'DIM DVOL!()          ' Diffuse layer volume of humics, fulvics
DIM DLF AS GLOBAL SINGLE
DIM wKZED AS GLOBAL SINGLE

'DIM WHAMID%()        ' WHAM database ID of needed components
DIM NWHAM AS GLOBAL INTEGER          ' Number of WHAM components
'DIM WSpec%()         ' Flag (when true, use WHAM calc spec)
'DIM CMap%()          ' Map of CHESS species to WHAM
DIM TEMPK AS GLOBAL SINGLE  ' Duplicate of SysTemp variable
DIM wISTR AS GLOBAL SINGLE           ' Dupl of ISTR#              0
DIM CRIT AS GLOBAL SINGLE
DIM PCO2 AS GLOBAL INTEGER
DIM VOLSOL AS GLOBAL SINGLE
DIM PHFIX  AS GLOBAL STRING
DIM TDCONCERR AS GLOBAL SINGLE
DIM TDCONC AS GLOBAL SINGLE
DIM TDCONCCALC AS GLOBAL SINGLE
DIM PRECISION AS GLOBAL SINGLE
DIM iDoneWHAM AS GLOBAL INTEGER
DIM PH AS GLOBAL SINGLE

'dim wPKH!()
'DIM wPKHA!()
'DIM wDPKHA!()
'DIM wPKHB!()
'DIM wDPKHB!()
'DIM wFPR!()
'DIM wNCOOH!()
'DIM wRATIO!()
'DIM wP!()
'DIM wZCALC!()
'DIM wZERR!()

DIM sWHAM AS GLOBAL STRING
DIM sBL AS GLOBAL STRING
'DIM sMetBL() AS STRING
DIM sMetal AS GLOBAL STRING
DIM sDOC AS GLOBAL STRING
DIM CritBL AS GLOBAL DOUBLE
DIM DefaultCriticalValue AS GLOBAL DOUBLE
DIM NBLM AS GLOBAL INTEGER
DIM iLA50Adjust AS GLOBAL INTEGER
DIM LA50PARAMS(1 TO 5) AS GLOBAL STRING
DIM NLA50Params AS GLOBAL INTEGER
DIM iLA50Variable AS GLOBAL INTEGER
DIM LA50Slope AS GLOBAL DOUBLE
DIM LA50Break AS GLOBAL DOUBLE
DIM LA50MinFrac AS GLOBAL DOUBLE
'
'  New to version 2.40: The following variables are for storing info related to a critical accumulation table (CAT)
'
DIM CAT_CA(1 TO 100) AS GLOBAL DOUBLE ' the critical accumulation
DIM CAT_Species(1 TO 100) AS GLOBAL STRING        ' the species name (latin or common)
DIM CAT_TestType(1 TO 100) AS GLOBAL STRING       ' the test type (e.g. acute or chronic)
DIM CAT_Lifestage(1 TO 100) AS GLOBAL STRING      ' lifestage
DIM CAT_Endpoint(1 TO 100) AS GLOBAL STRING       ' endoint (e.g., survival, reproduction, growth)
DIM CAT_Quantifier(1 TO 100) AS GLOBAL STRING     ' quantifier (e.g., EC20, EC10, etc)
DIM CAT_References(1 TO 100) AS GLOBAL STRING     ' references
DIM CAT_Miscellaneous(1 TO 100) AS GLOBAL STRING  ' any other comments relevant to this entry
DIM iHaveCAT AS GLOBAL INTEGER                    ' True/False indicates that a CAT will be used
DIM iCAT_Row AS GLOBAL INTEGER                    ' which row of the CAT will be used for this simulation
DIM iCAT_NRow AS GLOBAL INTEGER                   ' the number of rows in the table (must be >0 and <100)
 iHaveCAT = FALSE                                 ' Default is false

DIM ACUTE_DIV AS GLOBAL DOUBLE              '<- New for FMB version
DIM CHRONIC_DIV  AS GLOBAL DOUBLE           '<- New for FMB version
   ACUTE_DIV = 1                            '<- Set default to 1 so that the file doesn't need to specify it if it's only an HC5
   CHRONIC_DIV = 1                          '<- Set default to 1 so that the file doesn't need to specify it if it's only an HC5
DIM Metal_DL AS GLOBAL STRING               '<- New for FMB version
DIM Z_EF AS GLOBAL DOUBLE                   '<- New for FMB version
DIM D_EF AS GLOBAL LONG                     '<- New for FMB version, days of EF
DIM Z_EF_Def AS GLOBAL DOUBLE               '<- New for FMB version
 Z_EF_Def = 3.117085807
DIM FMB_option AS GLOBAL INTEGER            '<- New for FMB version
DIM FMBMean AS GLOBAL INTEGER               '<- New for FMB version, option to specify FMB based on means
DIM FMBMedian AS GLOBAL INTEGER             '<- New for FMB version, option to specify FMB based on medians
 FMBMean = 1
 FMBMedian = 2
 FMB_option = FMBMean                       'default to FMB calculation based on means
DIM iCalc_FMB AS GLOBAL INTEGER                    '<- New for FMB version, option to calculate FMB or just WQC
 iCalc_FMB = FALSE

'COMMON SHARED Quote$


' *--------------------*
' |  Simulation Types  |
' *--------------------*
DIM stBatch AS GLOBAL INTEGER
DIM stTitration AS GLOBAL INTEGER
DIM stColumn AS GLOBAL INTEGER
 stBatch = 1
 stTitration = 2
 stColumn = 3

' *-----------------------------*
' |  Davies Equation constants  |
' *-----------------------------*
DIM deA AS GLOBAL SINGLE
DIM deB AS GLOBAL SINGLE
 deA = .5
 deB = .2

' *-------------------------------*
' |  Faraday constant (Coul/mol)  |
' *-------------------------------*
DIM Fcon AS GLOBAL DOUBLE
 Fcon# = 96480!

' *----------------------------*
' |  Gas Constant (J/[Kúmol])  |
' *----------------------------*
DIM Rcon AS GLOBAL DOUBLE
 Rcon# = 8.314

' *-------------------------------*
' |  Default temperature assumed  |
' |  for version 1 data files     |
' *-------------------------------*
DIM DefaultTemp AS GLOBAL SINGLE
 DefaultTemp = 298

' *-------------------------------------------------------*
' |  Lambda used for EDL calculations, From Morel, p 416  |
' *-------------------------------------------------------*
DIM Lambda AS GLOBAL DOUBLE
 Lambda# = -.00063

' *---------------------------------------------------------*
' |  Smallest allowable tranformed species stoichiometric   |
' |  coefficient.  Coefficients with values smaller than    |
' |  this will be assumed due to round off and set to zero  |
' *---------------------------------------------------------*
DIM RoundOffLimit AS GLOBAL SINGLE
 RoundOffLimit = .1

' *-----------------*
' |  Species Types  |
' *-----------------*
DIM stAqueous AS GLOBAL INTEGER
DIM stSurface AS GLOBAL INTEGER
DIM stPrecipitate AS GLOBAL INTEGER
 stAqueous = 0
 stSurface = 1
 stPrecipitate = 2

' *-------------------*
' |  Component Types  |
' *-------------------*
DIM ctFloat AS GLOBAL INTEGER
DIM ctFixed AS GLOBAL INTEGER
DIM ctSubstituted AS GLOBAL INTEGER
DIM ctChargeBal AS GLOBAL INTEGER
DIM ctSurfPot AS GLOBAL INTEGER
 ctFloat = 1
 ctFixed = 2
 ctSubstituted = 3
 ctChargeBal = 4
 ctSurfPot = 15

' *---------------*
' |  Phase Types  |
' *---------------*
DIM ptUnlimited AS GLOBAL INTEGER
DIM ptSaturated AS GLOBAL INTEGER
DIM ptUnSaturated AS GLOBAL INTEGER
DIM ptNotAvailable AS GLOBAL INTEGER
 ptUnlimited = -1
 ptSaturated = 1
 ptUnSaturated = 2
 ptNotAvailable = 3

' *-------------------------------*
' |  Activity Correction Methods  |
' *-------------------------------*
DIM acNone AS GLOBAL INTEGER
DIM acDavies AS GLOBAL INTEGER
DIM acDebye AS GLOBAL INTEGER
DIM acMoleFraction AS GLOBAL INTEGER
DIM acChargeFraction AS GLOBAL INTEGER
DIM acWHAM AS GLOBAL INTEGER
 acNone = 1
 acDavies = 2
 acDebye = 3
 acMoleFraction = 4
 acChargeFraction = 5
 acWHAM = 6

' *--------------*
' |  File Types  |
' *--------------*
DIM ftASCII AS GLOBAL INTEGER
DIM ftLotus AS GLOBAL INTEGER
 ftASCII = 1
 ftLotus = 2
'
' END CONSTANT LIST
'
'
' *--------------------------*
' |  WHAM Related Constants  |
' *--------------------------*

'COMMON SHARED /DefaultVals/ BuffSize AS INTEGER
'COMMON SHARED /DefaultVals/ OutputMask AS STRING

'
' These WHAM related variables used to be in DefineProblem
'
   ' Dimension statements were moved to here
   NSP% = 415
   DIM wN(NSP%) AS GLOBAL INTEGER          ' numerical identification of species
   DIM wNs(NSP%) AS GLOBAL STRING         ' nominal identification of species
   DIM CMap(NSP%) AS GLOBAL INTEGER
   DIM wT(NSP%) AS GLOBAL SINGLE          ' total concentrations of components
   DIM wTH(2) AS GLOBAL SINGLE
   DIM wTCALC(NSP%) AS GLOBAL SINGLE               ' total calculated concentrations of components
   DIM wA(NSP%) AS GLOBAL SINGLE                 ' activities of solution species
   DIM wC(NSP%) AS GLOBAL SINGLE                 ' concentations of solution species
   DIM wCH(NSP%) AS GLOBAL INTEGER               ' charges on solution species
   DIM wM1(NSP%) AS GLOBAL INTEGER               ' master species identifiers
   DIM wM2(NSP%) AS GLOBAL INTEGER               ' master species identifiers
   DIM wM3(NSP%) AS GLOBAL INTEGER              ' master species identifiers
   DIM wS1(NSP%) AS GLOBAL INTEGER              ' stoichiometries
   DIM wS2(NSP%) AS GLOBAL INTEGER              ' stoichiometries
   DIM wS3(NSP%) AS GLOBAL INTEGER              ' stoichiometries
   DIM wLK(NSP%) AS GLOBAL SINGLE               ' equilibrium constants
   DIM wDH(NSP%) AS GLOBAL SINGLE                ' enthalpies
   DIM wGAMMA(5) AS GLOBAL SINGLE                ' activity coefficients                                                    w

   DIM wNCOOH(2) AS GLOBAL SINGLE
   DIM wPKHA(2) AS GLOBAL SINGLE
   DIM wPKHB(2) AS GLOBAL SINGLE
   DIM wDPKHA(2) AS GLOBAL SINGLE
   DIM wDPKHB(2) AS GLOBAL SINGLE
   DIM wP(2) AS GLOBAL SINGLE
   DIM wFPR(2) AS GLOBAL SINGLE
   DIM wDLF AS GLOBAL SINGLE
   'DIM wKZED!
   DIM wNODATA AS GLOBAL INTEGER
   DIM wRADIUS(2) AS GLOBAL SINGLE
   DIM wRATIO(2) AS GLOBAL SINGLE
   DIM wMOLWT(2) AS GLOBAL SINGLE
   DIM DVOLMAX(2) AS GLOBAL SINGLE
   DIM DVOL(2) AS GLOBAL SINGLE
   DIM DDLVOL(2) AS GLOBAL SINGLE
   DIM ZED(2) AS GLOBAL SINGLE

   DIM wX AS GLOBAL INTEGER
   'DIM wN$(wX%)
   'DIM wCH%(wX%)
   'DIM wM1%(wX%)
   'DIM wM2%(wX%)
   'DIM wM3%(wX%)
   'DIM wS1%(wX%)
   'DIM wS2%(wX%)
   'DIM wS3%(wX%)
   'DIM wLK(wX%)
   'DIM wDH(wX%)
   'wDH (wX%)
   'DIM wLK(wX%)
   DIM wLK(NSP%) AS GLOBAL SINGLE
   'DIM wPKMHA!(1, wX%)
   'DIM wPKMHA(2, wX%)
   DIM wZCALC(2) AS GLOBAL SINGLE
   DIM wZERR(2) AS GLOBAL SINGLE
   ' DIMension humic arrays ; 1 = HA, 2 = FA

   DIM wSP(2, 8) AS GLOBAL SINGLE               ' proton-binding sites not forming bidentate sites
   DIM wSM(2, 12) AS GLOBAL SINGLE              ' bidentate metal-binding sites
   DIM wKH(2, 8) AS GLOBAL SINGLE               ' K's for proton-binding
   DIM wPKH(2, 8) AS GLOBAL SINGLE              ' PK's for proton-binding
   DIM wKMH(2, 8, NSP%) AS GLOBAL SINGLE        ' K's for metal-proton exchange
   DIM wSITE(12, 2) AS GLOBAL INTEGER           ' proton sites making bidentate sites
   DIM wPKMHA(2, NSP%) AS GLOBAL SINGLE         ' metal-proton exchange constants
   DIM wPKMHB(2, NSP%) AS GLOBAL SINGLE         ' metal-proton exchange constants
   DIM wFP(8) AS GLOBAL SINGLE                  ' dissociation factors for proton-binding
   DIM wBIBTERM(NSP%) AS GLOBAL SINGLE          ' terms calc'd in finding bidentate theta's
   DIM wMONBTERM(NSP%) AS GLOBAL SINGLE         ' terms calc'd in finding monodentate theta's
   DIM wBISNU(12, NSP%) AS GLOBAL SINGLE        ' values of NU at each bidentate site
   DIM wMONSNU(8, NSP%) AS GLOBAL SINGLE        ' values of NU at each monodentate site
   DIM wBINU(NSP%) AS GLOBAL SINGLE             ' overall NU for each metal/bidentate sites
   DIM wMONNU(NSP%) AS GLOBAL SINGLE            ' overall NU for each metal/monodentate sites
   DIM wNU(2, NSP%) AS GLOBAL SINGLE            ' overall NU for each complexed species
   DIM wCHC(2, NSP%) AS GLOBAL SINGLE           ' conc of complexed species per litre
   DIM wCDDL(2, NSP%) AS GLOBAL SINGLE          ' concs in HS DDL's
   DIM wBIZ(12) AS GLOBAL SINGLE                ' net charge at each bidentate site
   DIM wMONZ(8) AS GLOBAL SINGLE                ' net charge at each monodentate site
   DIM wBINDTEST(NSP%) AS GLOBAL INTEGER          ' allows binding of a species to be tested

   DIM sMetBL(10) AS GLOBAL STRING
   DIM iBLM(10) AS GLOBAL INTEGER
   DIM iBLMTox(10) AS GLOBAL INTEGER
   DIM iOrig AS GLOBAL INTEGER
   DIM iAdjust AS GLOBAL INTEGER
   DIM iQuiet AS GLOBAL INTEGER
   DIM iVQuiet AS GLOBAL INTEGER
   DIM EchoInput AS GLOBAL INTEGER
   DIM iComplete AS GLOBAL INTEGER
   DIM iConverge AS GLOBAL INTEGER
   DIM iOFormat AS GLOBAL INTEGER
   DIM iCalcLA50 AS GLOBAL INTEGER
   DIM iCalcLA50Log AS GLOBAL INTEGER
   DIM iCalcLA50Al AS GLOBAL INTEGER
   DIM iLAFile AS GLOBAL INTEGER
   DIM iLowMass  AS GLOBAL INTEGER
   DIM iCriteria AS GLOBAL INTEGER
   DIM iNumZ AS GLOBAL INTEGER

BufSize = 800
OutputMask$ = "##.######^^^^ "

' *------------------*
' |  Default values  |
' *------------------*
MaxIter% = 1500
MaxECxIter = 200

' *-----------------------------------------------*
' |  Simple Variables Local to SUB DefineProblem  |
' *-----------------------------------------------*
'DIM c AS INTEGER            '  Array/loop index variable
'DIM ExitFlag1 AS INTEGER    '  Boolean flag used to test a conditional loop
'DIM i AS INTEGER            '  Array/loop index variable
'DIM iFileIn AS INTEGER      '  I/O handle for input file
'DIM iMonteInfo AS INTEGER   '  Flag to indicate this species will Monte Carlo
'DIM iPos AS INTEGER         '  Position locator of one string within another
'DIM j AS INTEGER            '  Array/loop index variable
'DIM MajorVer AS INTEGER     '  Variable for major portion of verion number
'DIM MinorVer AS INTEGER     '  Variable for minor portion of verion number
'DIM Nj AS INTEGER           '  Number of elements in a linked/list
'DIM Tmp AS SINGLE           '  Temporary single precision variable
'DIM TmpD AS DOUBLE          '  Temporary double precision variable
'DIM Tmp2 AS SINGLE          '  Temporary single precision variable
'DIM TmpS AS STRING          '  Temporary string variable
'DIM TmpI AS INTEGER         '  Temporary integer variable
'DIM TmpI1 AS INTEGER        '  Temporary integer variable
'DIM TmpI2 AS INTEGER        '  Temporary integer variable

' *-------------------------------------*
' |  Default Sizes for Array Variables  |
' *-------------------------------------*
DIM MaxComp  AS GLOBAL INTEGER
DIM MaxSpecies  AS GLOBAL INTEGER
DIM MaxPhases  AS GLOBAL INTEGER
DIM MaxLists  AS GLOBAL INTEGER
DIM MaxMass  AS GLOBAL INTEGER
 MaxComp = 70
 MaxSpecies = 600
 MaxPhases = 12
 MaxLists = 2
 MaxMass = 5
 MaxObs = 1000
DIM AmbientMetal(1 TO MaxObs) AS GLOBAL DOUBLE            '<- New for FMB version
DIM HC5(1 TO MaxObs) AS GLOBAL DOUBLE                     '<- New for FMB version
DIM CenFlag(1 TO MaxObs) AS GLOBAL LONG                   '<- New for FMB version
DIM TU_Acute(1 TO MaxObs) AS GLOBAL EXTENDED              '<- New for FMB version
DIM TU_Chron(1 TO MaxObs) AS GLOBAL EXTENDED              '<- New for FMB version
DIM WQC_Acute(1 TO MaxObs) AS GLOBAL EXTENDED              '<- New for FMB version
DIM WQC_Chron(1 TO MaxObs) AS GLOBAL EXTENDED              '<- New for FMB version
DIM MetalConc(1 TO MaxObs) AS GLOBAL EXTENDED          '<- New for FMB version
'DIM MetalFlag(1 TO MaxObs) AS LONG              '<- New for FMB version

' *-----------------------------------------------------------------*
' |  Array Variables Defined Globally and used throughout the model |
' *-----------------------------------------------------------------*
DIM ActCorr(MaxSpecies) AS GLOBAL  INTEGER
DIM ActCoef(MaxSpecies) AS GLOBAL  DOUBLE
DIM CtoM(MaxSpecies) AS GLOBAL  SINGLE
DIM CType(MaxComp) AS GLOBAL  INTEGER
DIM ChangeMass(MaxPhases) AS GLOBAL  DOUBLE
DIM CompOrder(MaxComp) AS GLOBAL  INTEGER
DIM DeltaH(MaxSpecies) AS GLOBAL  SINGLE
DIM EC(MaxSpecies) AS GLOBAL  SINGLE
DIM ECTemp(MaxSpecies) AS GLOBAL  SINGLE
DIM hEC(MaxSpecies) AS GLOBAL  SINGLE
DIM hSEC(MaxPhases)  AS GLOBAL  SINGLE
DIM InitialConc(MaxSpecies) AS GLOBAL  SINGLE
DIM KT(MaxComp) AS GLOBAL  DOUBLE
DIM OKT(MaxComp) AS GLOBAL DOUBLE
DIM Label(MaxSpecies + MaxPhases) AS GLOBAL  STRING
DIM LinkComp(MaxSpecies) AS GLOBAL  INTEGER
DIM LinkList(MaxLists, MaxSpecies) AS GLOBAL  INTEGER
DIM ListPos(MaxComp) AS GLOBAL  INTEGER
DIM MeanLogK(MaxSpecies) AS GLOBAL  SINGLE
DIM MonteList(MaxSpecies) AS GLOBAL  INTEGER
DIM OC(MaxSpecies + FOUR) AS GLOBAL  DOUBLE
DIM oldOC(0 TO MaxComp) AS GLOBAL DOUBLE
DIM PhaseMass(MaxPhases) AS GLOBAL  SINGLE
DIM PhaseOrder(MaxPhases) AS GLOBAL  INTEGER
DIM PSubOrder(MaxPhases) AS GLOBAL  INTEGER
DIM PType(MaxPhases) AS GLOBAL  INTEGER
DIM SCa(MaxSpecies, MaxComp) AS GLOBAL  SINGLE
DIM SCb(MaxSpecies, MaxComp) AS GLOBAL  SINGLE
DIM SComp(MaxPhases) AS GLOBAL  INTEGER
DIM SDeltaH(MaxPhases) AS GLOBAL  SINGLE
DIM SEC(MaxPhases)  AS GLOBAL  SINGLE
DIM SECTemp(MaxPhases) AS GLOBAL  SINGLE
DIM SiteDen(MaxComp) AS GLOBAL  SINGLE
DIM SP(MaxPhases, MaxComp) AS GLOBAL  SINGLE
DIM SPhase(MaxComp) AS GLOBAL  INTEGER
DIM STempCorr(MaxPhases) AS GLOBAL  INTEGER
DIM SType(MaxSpecies) AS GLOBAL  INTEGER
DIM tEC(MaxSpecies) AS GLOBAL  SINGLE
DIM TempCorr(MaxSpecies) AS GLOBAL  INTEGER
DIM tSCa(MaxSpecies, MaxComp) AS GLOBAL  SINGLE
DIM tSCb(MaxSpecies, MaxComp) AS GLOBAL  SINGLE
DIM tSEC(MaxPhases)  AS GLOBAL  SINGLE
DIM tSP(MaxPhases, MaxComp) AS GLOBAL  SINGLE
DIM V(MaxSpecies) AS GLOBAL  SINGLE
DIM VarLogK(MaxSpecies) AS GLOBAL  SINGLE

DIM MassName(MaxMass) AS GLOBAL  STRING           ' The name of a mass compartment
DIM MassVal(MaxMass) AS GLOBAL  DOUBLE            ' The mass value
DIM MassUnitLabel(MaxMass) AS GLOBAL  STRING      ' The units of a mass compartment (e.g. L or kg)

DIM WSpec(MaxSpecies) AS GLOBAL INTEGER
DIM WHAMMAP(0 TO MaxComp) AS GLOBAL INTEGER       ' Map of WHAM components to CHESS
DIM wOrgBound(MaxSpecies) AS GLOBAL SINGLE
DIM HABound(MaxSpecies) AS GLOBAL DOUBLE
DIM FABound(MaxSpecies) AS GLOBAL DOUBLE

DIM CompList(0 TO MaxSpecies, 0 TO MaxComp) AS GLOBAL INTEGER
DIM SpecList(0 TO MaxComp, 0 TO MaxSpecies)  AS GLOBAL INTEGER


' *------------------------------------------------------*
' |  Dimension variables needed for N-R iteration        |
' *------------------------------------------------------*
DIM Z(MaxComp, MaxComp) AS GLOBAL  DOUBLE
DIM Zn(MaxComp, MaxComp) AS GLOBAL  DOUBLE
DIM R(MaxComp) AS GLOBAL  DOUBLE
DIM X(MaxComp) AS GLOBAL  DOUBLE
