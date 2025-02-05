test_that("ConvertWindowsParamFile 3.00 works", {

  # Make test data
  Y_3.00_text = c(
    "Column model parameter file, Ver 3.00",
    "--------------------------------------------------------------------------------",
    "Number of Components, Species, Phases, Linked Lists",
    "12,  4,  0,  0",
    "Component Charge  Type  Activity  Site Den",
    "H,         1,        2,        2,        1    ",
    "Y,         2,        1,        2,        1",
    "DOC,      -1,        1,        1,        .0006",
    "Ca,        2,        1,        2,        1",
    "Mg,        2,        1,        2,        1",
    "Na,        1,        1,        2,        1",
    "K,         1,        1,        2,        1",
    "SO4,      -2,        1,        2,        1",
    "Cl,       -1,        1,        2,        1",
    "CO3,      -2,        1,        2,        1",
    "S,        -2,        1,        2,        1",
    "Gill,     -1,       11,        1,        1000.E-06",
    "--------------------------------------------------------------------------------",
    "Species  Type Act. iTemp iMonte pH  Y   DOC Ca   Mg  Na  K   SO4 Cl  CO3 S   BL   Log K      Var Log K  Delta H  Temp   Conc",
    "Gill-Y,    11,  1,   0,    0,     0,  1,  0,  0,   0,  0,  0,  0,  0,  0,  0,  1,   4.000      0.000      0.000    0.000  0.0000E+00",
    "Gill-Ca,   11,  1,   0,    0,     0,  0,  0,  1,   0,  0,  0,  0,  0,  0,  0,  1,   4.250      0.000      0.000    0.000  0.0000E+00",
    "Gill-Mg,   11,  1,   0,    0,     0,  0,  0,  0,   1,  0,  0,  0,  0,  0,  0,  1,   3.600      0.000      0.000    0.000  0.0000E+00",
    "Gill-H,    11,  1,   0,    0,     1,  0,  0,  0,   0,  0,  0,  0,  0,  0,  0,  1,   4.700      0.000      0.000    0.000  0.0000E+00",
    "--------------------------------------------------------------------------------",
    "Linked Lists",
    "--------------------------------------------------------------------------------",
    "Phase    iTemp iMonte  H   Ni DOC Ca   Mg   Na  K SO4  Cl CO3   S BL  Log Ks  Var Log K   Delta H    Temp        Moles",
    "--------------------------------------------------------------------------------",
    "File I/O and System Description",
    "                              |  File with sequence of input chemistry (next line)",
    "D:\\MODELS\\CHESS315\\CU\\LC50.WK1                                              ",
    " 2                            |  Type of file (1 = ASCII, 2 = LOTUS)",
    "                              |",
    "                              |  Name of output file (next line)",
    "lc50??.WKS                                                                  ",
    " 2                            |  Type of file (1 = ASCII, 2 = LOTUS)",
    "                              |",
    " 1                            |  Simulation Type",
    "                              |  (1 = Batch, 2 = Titration, 3 = Column)",
    "                              |",
    " 1                            |  Number of Layers",
    "                              |",
    " 0                            |  Number of Monte Carlo Simulations",
    "                              |",
    " 1.0000E-06                   |  Error Tolerance",
    "                              |",
    " 500                          |  Maximum Iterations",
    "                              |",
    "                              |  Output file contents, do you want:",
    "-1                            |     Concentrations?",
    "-1                            |     Totals?",
    " 0                            |     Equilibrium constants?",
    " 0                            |     Residuals?",
    " 0                            |     Activity Coefficients?",
    "-1                            |     Phase Sat. Indices?",
    " 0                            |     Phase Ion Activity Products?",
    "",
    "User Notes:",
    "",
    "[Metal]: Y",
    "",
    "[Gill]: Gill",
    "",
    "[Gill-Metal]: Gill-Y",
    "",
    "[DOC]: DOC",
    "",
    "[THERMO]: NIST_20170203.dbs",
    "",
    "[CRITICAL]: 0.03395",
    "",
    "[ACUTE_DIV]: 2.0",
    "[CHRONIC_DIV]: 3.22",
    "[Z_EF]: 3.117",
    "",
    "[END]:",
    "",
    "",
    "",
    ""
  )
  Y_3.00_file = withr::local_tempfile(fileext = ".dat")
  write(x = Y_3.00_text, file = Y_3.00_file)
  ex_Rfile = withr::local_tempfile(fileext = ".dat4")

  expect_no_error(ConvertWindowsParamFile(WindowsParamFile = Y_3.00_file))
  expect_no_error(ConvertWindowsParamFile(WindowsParamFile = Y_3.00_file,
                                          RParamFile = ex_Rfile))
  expect_no_error(DefineProblem(ParamFile = ex_Rfile))

})
test_that("ConvertWindowsParamFile 3.01 works", {

  # Make test data
  Cu_3.01_text = c(
    "Column model parameter file, Ver 3.01",
    "--------------------------------------------------------------------------------",
    "Number of Components, Species, Phases, Linked Lists",
    "12,  6,  0,  0",
    "Component Charge  Type  Activity  Site Den",
    "H,         1,        2,        2,        1    ",
    "Cu,        2,        2,        2,        1",
    "DOC,      -1,        2,        1,        .0006",
    "Ca,        2,        2,        2,        1",
    "Mg,        2,        2,        2,        1",
    "Na,        1,        2,        2,        1",
    "K,         1,        2,        2,        1",
    "SO4,      -2,        2,        2,        1",
    "Cl,       -1,        2,        2,        1",
    "CO3,      -2,        2,        2,        1",
    "S,        -2,        2,        2,        1",
    "BL,       -1,       11,        1,        3E-05",
    "--------------------------------------------------------------------------------",
    "Species  Type Act. iTemp iMonte pH  Cu DOC   Ca  Mg  Na   K SO4  Cl CO3  S BL   Log K  Var Log K   Delta H    Temp        Conc",
    "BL-Cu,   11,  1,   0,    0,     0,  1,  0,  0,   0,  0,  0,  0,  0,  0,  0,  1,   7.400      0.000     0.000   0.000  0.0000E+00",
    "BL-CuOH, 11,  1,   0,    0,    -1,  1,  0,  0,   0,  0,  0,  0,  0,  0,  0,  1,  -1.300      0.000     0.000   0.000  0.0000E+00",
    "BL-Ca,   11,  1,   0,    0,     0,  0,  0,  1,   0,  0,  0,  0,  0,  0,  0,  1,   3.600      0.000     0.000   0.000  0.0000E+00",
    "BL-Mg,   11,  1,   0,    0,     0,  0,  0,  0,   1,  0,  0,  0,  0,  0,  0,  1,   3.600      0.000     0.000   0.000  0.0000E+00",
    "BL-H,    11,  1,   0,    0,     1,  0,  0,  0,   0,  0,  0,  0,  0,  0,  0,  1,   5.400      0.000     0.000   0.000  0.0000E+00",
    "BL-Na,   11,  1,   0,    0,     0,  0,  0,  0,   0,  1,  0,  0,  0,  0,  0,  1,   3.000      0.000     0.000   0.000  0.0000E+00",
    "--------------------------------------------------------------------------------",
    "Linked Lists",
    "--------------------------------------------------------------------------------",
    "Phase    iTemp iMonte  H   Cu DOC Ca   Mg   Na  K SO4  Cl CO3   S BL  Log Ks  Var Log K   Delta H    Temp        Moles",
    "--------------------------------------------------------------------------------",
    "File I/O and System Description",
    "                              |  File with sequence of input chemistry (next line)",
    "D:\\MODELS\\CHESS315\\CU\\LC50.WK1                                              ",
    " 2                            |  Type of file (1 = ASCII, 2 = LOTUS)",
    "                              |",
    "                              |  Name of output file (next line)",
    "lc50??.WKS                                                                  ",
    " 2                            |  Type of file (1 = ASCII, 2 = LOTUS)",
    "                              |",
    " 1                            |  Simulation Type",
    "                              |  (1 = Batch, 2 = Titration, 3 = Column)",
    "                              |",
    " 1                            |  Number of Layers",
    "                              |",
    " 0                            |  Number of Monte Carlo Simulations",
    "                              |",
    " 1.0000E-06                   |  Error Tolerance",
    "                              |",
    " 500                          |  Maximum Iterations",
    "                              |",
    "                              |  Output file contents, do you want:",
    "-1                            |     Concentrations?",
    "-1                            |     Totals?",
    " 0                            |     Equilibrium constants?",
    " 0                            |     Residuals?",
    " 0                            |     Activity Coefficients?",
    "-1                            |     Phase Sat. Indices?",
    " 0                            |     Phase Ion Activity Products?",
    "",
    "User Notes:",
    "",
    "[Metal]: Cu",
    "",
    "[BL]: BL",
    "",
    "[BL-Metal]: BL-Cu",
    "[BL-Metal]: BL-CuOH",
    "",
    "[DOC]: DOC",
    "",
    "[THERMO]: WATER23.DBS",
    "",
    "[CRITICAL START]:",
    "CA,Species,Test Type,Lifestage,Endpoint,Quantifier,References,Miscellaneous",
    "0.05541,Ceriodaphnia dubia,Acute,Neonate (<24 hr),,EC50; LC50,Gensemer et al. 2002; Hyne et al. 2005; Naddy et al. 2002; Naddy et al. 2003; Van Genderen et al. 2007,SMEA calculated by median",
    "0.0099515,Ceriodaphnia dubia,Chronic,Neonate (<24 hr),reproduction - # of young,EC20,Wang et al. 2011a,SMEA calculated by geomean",
    "1.0223,Chironomus tentans,Acute,\"Larva, 1st instar\",death and immobility,EC50,Gauss et al. 1985,SMEA calculated by geomean",
    "0.013496,Acipenser transmontanus,Acute,\"Columbia River, 26 dph, 0.08 g, 2.5 cm; Kootenai River, 38 dph, 0.07 g, 2.4 cm; Juvenile, 40 dph\",,LC50,Little et al. 2012; Vardy et al. 2013,SMEA calculated by geomean",
    "3.70,Oncorhynchus mykiss,Acute,,,,From the \"Cu_Rainbow_Trout_06-10-07.DAT\" parameter file,",
    "0.119,Daphnia magna,Acute,,,,From the \"Cu_Daphnia_Magna_06-10-07.DAT\" parameter file,",
    "0.0447,Daphnia pulex,Acute,,,,From the \"Cu_Daphnia_Pulex_06-10-07.DAT\" parameter file,",
    "5.48,Pimephales promelas,Acute,,,,From the \"Cu_Fathead_Minnow_06-10-07.DAT\" parameter file,",
    "0.021704,Ceriodaphnia dubia,Acute,Neonate (<24h),survival,EC50,\"Larry Walker Associates, 2013\",",
    "[CRITICAL END]",
    "",
    "[END]:",
    "",
    "",
    "",
    "#==================================================================================================",
    "# NOTES:",
    "",
    "2017-01-16: (KEC) This file is based on \"CuOH5%le_10-11-07.DAT\", the 2007 EPA FW Cu criteria file.  ",
    "It was converted to a version 3.01 parameter file by adding in a table of critical accumulation ",
    "data in place of the [CRITICAL] line.  The critical table is based on data in the 2016 ",
    "copper freshwater toxicity database.  This is part of a project for ICA.  The new file name is ",
    "\"Cu_freshwater_acute_and_chronic_2017-01-16.dat\".",
    "",
    "2017-01-17: (KEC) Renamed to \"Cu_freshwater_acute_and_chronic_2017-01-17.dat\".  Added lines in the ",
    "critical table so that the parameter files typically included in the BLM UI have their own lines.",
    "This is mostly to consolidate the number of files, and to allow for more precise labeling of ",
    "these organisms in the organism selection table. Also added a C. dubia calibration from the ",
    "Pacific EcoRisk LA River dataset:",
    "  ",
    "Larry Walker Associates, 2013. Comparison of Biotic Ligand Model (BLM) Results with Toxicity ",
    "Testing Data for the Copper Water-Effect Ratio (WER) Study for the Los Angeles River and ",
    "its Tributaries.  Memo to Jenny Newman, LA River Water Quality Control Board, May 2013)",
    "",
    "This dataset is the same which is used in the \"Cu_Ceriodaphnia_dubia_2013-06-21.dat\" parameter ",
    "file, but with speciation performed with a parameter file that has BL-Mg as a chemical species.",
    "",
    "#=================================================================================================="
  )
  Cu_3.01_file = withr::local_tempfile(fileext = ".dat")
  write(x = Cu_3.01_text, file = Cu_3.01_file)
  ex_Rfile = withr::local_tempfile(fileext = ".dat4")

  expect_no_error(ConvertWindowsParamFile(WindowsParamFile = Cu_3.01_file))
  expect_no_error(ConvertWindowsParamFile(WindowsParamFile = Cu_3.01_file,
                                          RParamFile = ex_Rfile))
  expect_no_error(DefineProblem(ParamFile = ex_Rfile))

})
test_that("ConvertWindowsParamFile marine works", {

  # Make test data
  Marine_text = c(
    "Column model parameter file, Ver 3.00. Stability constants and enthalpies from NIST Database 46 v. 8, except for Cu+2 + OH- complexes, which are from Paulson & Kester (1980) (stability constants) and Byrne et al. (1988) (enthalpies), entha",
    "--------------------------------------------------------------------------------",
    "Number of Components, Species, Phases, Linked Lists",
    "16,       58,     1,               0",
    "Component Charge  Type    Activity  Site Den",
    "H,        1,      2,      2,               1",
    "Cu,       2,      1,      2,               1",
    "L1,       -1,     1,      1,        3.2351474E-8",
    "L2,       -1,     1,      1,        1.2436842E-7",
    "L3,       -1,     1,      1,        6.8075071E-7",
    "L4,       -1,     1,      1,        9.0381860E-7",
    "Ca,       2,      1,      2,               1",
    "Mg,       2,      1,      2,               1",
    "Na,       1,      1,      2,               1",
    "K,        1,      1,      2,               1",
    "SO4,      -2,     1,      2,               1",
    "Cl,       -1,     1,      2,               1",
    "CO3,      -2,     1,      2,               1",
    "PO4,      -3,     1,      2,               1",
    "BL1,      -1,     11,     1,        2.71E-05",
    "BL2,      -1,     11,     1,         0.09113",
    "--------------------------------------------------------------------------------",
    "Species   Type    Act.    iTemp     iMonte  pH      Cu      L1      L2      L3      L4      Ca      Mg      Na      K       SO4     Cl      CO3     PO4     BL1     BL2     Log K    Var Log Delta H Temp    Conc",
    "BL1-Cu,   11,     1,      -1,       0,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      16.796017       0       0  298.15  0.00E+00",
    "BL1-CuOH, 11,     1,      -1,       0,      -1,     1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      8.0960167       0       0  298.15  0.00E+00",
    "BL1-Ca,   11,     1,      -1,       0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0,      0,      1,      0,            4.3       0       0  298.15  0.00E+00",
    "BL1-H,    11,     1,      -1,       0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,            5.4       0       0  298.15  0.00E+00",
    "BL1-Na,   11,     1,      -1,       0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0,      1,      0,            4.0       0       0  298.15  0.00E+00",
    "BL2-Cu,   11,     1,      -1,       0,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,           8.75       0       0  298.15  0.00E+00",
    "BL2-CuOH, 11,     1,      -1,       0,      -1,     1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,           0.06       0       0  298.15  0.00E+00",
    "BL2-Ca,   11,     1,      -1,       0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      1,            1.0       0       0  298.15  0.00E+00",
    "BL2-H,    11,     1,      -1,       0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,          -2.25       0       0  298.15  0.00E+00",
    "BL2-Na,   11,     1,      -1,       0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0,      1,            0.5       0       0  298.15  0.00E+00",
    "OH,       1,      2,      -1,       0,      -1,     0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,        -13.997       0  -55810  298.15  0.00E+00",
    "HCO3,     1,      2,      -1,       0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,         10.329       0  -14600  298.15  0.00E+00",
    "H2CO3,    1,      2,      -1,       0,      2,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,         16.681       0  -23760  298.15  0.00E+00",
    "HPO4,     1,      2,      -1,       0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,         12.375       0  -15000  298.15  0.00E+00",
    "H2PO4,    1,      2,      -1,       0,      2,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,         19.573       0  -18000  298.15  0.00E+00",
    "H3PO4,    1,      2,      -1,       0,      3,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,         21.727       0  -10500  298.15  0.00E+00",
    "MgHCO3,   1,      2,      -1,       0,      1,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      1,      0,      0,      0,         11.339       0   -9600  298.15  0.00E+00",
    "MgCO3,    1,      2,      -1,       0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      1,      0,      0,      0,           2.92       0   10000  298.15  0.00E+00",
    "MgSO4,    1,      2,      -1,       0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      1,      0,      0,      0,      0,      0,           2.26       0    5800  298.15  0.00E+00",
    "MgHPO4,   1,      2,      -1,       0,      1,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0,      1,      0,      0,         15.175       0   -3000  298.15  0.00E+00",
    "CaHCO3,   1,      2,      -1,       0,      1,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0,      1,      0,      0,      0,         11.529       0    4400  298.15  0.00E+00",
    "CaCO3,    1,      2,      -1,       0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0,      1,      0,      0,      0,           3.22       0   15000  298.15  0.00E+00",
    "CaSO4,    1,      2,      -1,       0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      1,      0,      0,      0,      0,      0,           2.36       0    7100  298.15  0.00E+00",
    "CaHPO4,   1,      2,      -1,       0,      1,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0,      1,      0,      0,         15.035       0   -3000  298.15  0.00E+00",
    "NaCO3,    1,      2,      -1,       0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      1,      0,      0,      0,           1.27       0       0  298.15         0",
    "CuOH,     1,      2,      -1,       0,      -1,     1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,          -7.96       0   50200  298.15  0.00E+00",
    "Cu(OH)2,  1,      2,      -1,       0,      -2,     1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,          -16.2       0   92100  298.15  0.00E+00",
    "CuHCO3,   1,      2,      -1,       0,      1,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,         12.129       0       0  298.15  0.00E+00",
    "CuCO3,    1,      2,      -1,       0,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,           6.77       0   13000  298.15  0.00E+00",
    "Cu(CO3)2, 1,      2,      -1,       0,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      2,      0,      0,      0,           10.2       0       0  298.15  0.00E+00",
    "CuSO4,    1,      2,      -1,       0,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0,           2.36       0    8700  298.15  0.00E+00",
    "CuCl,     1,      2,      -1,       0,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,            0.4       0    8300  298.15  0.00E+00",
    "CuHPO4,   1,      2,      -1,       0,      1,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,          15.68       0       0  298.15  0.00E+00",
    "CuH2PO4,  1,      2,      -1,       0,      2,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,         22.473       0   11000  298.15  0.00E+00",
    "CuL1,     1,      1,      -1,       0,      0,      1,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,  1.4456988E+01       0       0  298.15      0.00",
    "CuL2,     1,      1,      -1,       0,      0,      1,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,  1.1117723E+01       0       0  298.15      1.91",
    "CuL3,     1,      1,      -1,       0,      0,      1,      0,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,  9.3865796E+00       0       0  298.15      3.78",
    "CuL4,     1,      1,      -1,       0,      0,      1,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,  8.2358274E+00       0       0  298.15      5.78",
    "HL1,      1,      1,      -1,       0,      1,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,  7.5000000E+00       0       0  298.15      6.41",
    "HL2,      1,      1,      -1,       0,      1,      0,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,  7.5000000E+00       0       0  298.15      6.41",
    "HL3,      1,      1,      -1,       0,      1,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,  7.5000000E+00       0       0  298.15      6.41",
    "HL4,      1,      1,      -1,       0,      1,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,  7.5000000E+00       0       0  298.15      6.41",
    "CaL1,     1,      1,      -1,       0,      0,      0,      1,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,  3.5000000E+00       0       0  298.15      2.07",
    "CaL2,     1,      1,      -1,       0,      0,      0,      0,      1,      0,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,  7.1000000E-01       0       0  298.15      2.07",
    "CaL3,     1,      1,      -1,       0,      0,      0,      0,      0,      1,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,  -1.5600000000       0       0  298.15      2.07",
    "CaL4,     1,      1,      -1,       0,      0,      0,      0,      0,      0,      1,      1,      0,      0,      0,      0,      0,      0,      0,      0,      0,  -3.0900000000       0       0  298.15      2.07",
    "MgL1,     1,      1,      -1,       0,      0,      0,      1,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,  3.3000000E+00       0       0  298.15      2.24",
    "MgL2,     1,      1,      -1,       0,      0,      0,      0,      1,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,  5.1000000E-01       0       0  298.15      2.24",
    "MgL3,     1,      1,      -1,       0,      0,      0,      0,      0,      1,      0,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,  -1.7600000000       0       0  298.15      2.24",
    "MgL4,     1,      1,      -1,       0,      0,      0,      0,      0,      0,      1,      0,      1,      0,      0,      0,      0,      0,      0,      0,      0,  -3.2900000000       0       0  298.15      2.24",
    "NaL1,     1,      1,      -1,       0,      0,      0,      1,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0,      0,  1.1000000E+00       0       0  298.15      5.30",
    "NaL2,     1,      1,      -1,       0,      0,      0,      0,      1,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0,      0,  -1.6900000000       0       0  298.15      5.30",
    "NaL3,     1,      1,      -1,       0,      0,      0,      0,      0,      1,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0,      0,  -3.9600000000       0       0  298.15      5.30",
    "NaL4,     1,      1,      -1,       0,      0,      0,      0,      0,      0,      1,      0,      0,      1,      0,      0,      0,      0,      0,      0,      0,  -5.4900000000       0       0  298.15      5.30",
    "KL1,      1,      1,      -1,       0,      0,      0,      1,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0,  4.0000000E-01       0       0  298.15       5.5",
    "KL2,      1,      1,      -1,       0,      0,      0,      0,      1,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0,  -2.3900000000       0       0  298.15       5.5",
    "KL3,      1,      1,      -1,       0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0,  -4.6600000000       0       0  298.15       5.5",
    "KL4,      1,      1,      -1,       0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0,  -6.1900000000       0       0  298.15       5.5",
    "--------------------------------------------------------------------------------",
    "Linked Lists",
    "--------------------------------------------------------------------------------",
    "Phase    iTemp iMonte pH  Cu L1  L2  L3  L4  Ca  Mg   Na  K  SO4  Cl CO3 PO4 BL1 BL2   Log K    Var Log K   Delta H    Temp        Moles",
    "CO2(g),      0,    0, 2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  -21.567     0.000       0.000    0.000        -1",
    "--------------------------------------------------------------------------------",
    "File I/O and System Description",
    "                              |  File with sequence of input chemistry (next line)",
    "D:\\MODELS\\CHESS315\\CU\\LC50.WK1",
    " 2                            |  Type of file (1 = ASCII, 2 = LOTUS)",
    "                              |",
    "                              |  Name of output file (next line)",
    "lc50??.WKS",
    " 2                            |  Type of file (1 = ASCII, 2 = LOTUS)",
    "                              |",
    " 1                            |  Simulation Type",
    "                              |  (1 = Batch, 2 = Titration, 3 = Column)",
    "                              |",
    " 1                            |  Number of Layers",
    "                              |",
    " 0                            |  Number of Monte Carlo Simulations",
    "                              |",
    " 1.0000E-06                   |  Error Tolerance",
    "                              |",
    " 500                          |  Maximum Iterations",
    "                              |",
    "                              |  Output file contents, do you want:",
    "-1                            |     Concentrations?",
    "-1                            |     Totals?",
    " 0                            |     Equilibrium constants?",
    " 0                            |     Residuals?",
    " 0                            |     Activity Coefficients?",
    "-1                            |     Phase Sat. Indices?",
    " 0                            |     Phase Ion Activity Products?",
    " ",
    "User Notes:",
    " ",
    "[Metal]: Cu",
    " ",
    "[BL]: BL2",
    " ",
    "[BL-Metal]: BL1-Cu",
    "[BL-Metal]: BL1-CuOH",
    "[BL-Metal]: BL2-Cu",
    "[BL-Metal]: BL2-CuOH",
    " ",
    "[DOC]: DOC",
    " ",
    "[THERMO]:",
    " ",
    "[CRITICAL]: 303.6464738",
    "",
    " ",
    " ",
    "[END]:",
    " ",
    "          From Oxalate                              From NOM",
    "          Log K                                     Target dMin     Max",
    " ",
    "Cu:           4.85(25, 0.1)                0",
    "Zn:           3.58(25, 0.7)            -1.27           -0.49   -0.67   -0.36           -0.78",
    "Ca:           3.19(25, 0.0)            -1.66           -1.07   -1.35   -0.90           -0.59",
    "Mg:           2.76(20, 0.1)            -2.09           -1.24   -1.56   -1.06           -0.85",
    "Na:            0.5(25, 0.1)            -4.35",
    "K:            0.33(25, 0.1)            -4.52",
    "lpy for CuCO3 is from Byrne et al. (1988), also CuHPO4 is from a previous version, source is unknown",
    " ",
    " ",
    " ",
    " ",
    " ",
    " ",
    " ",
    " ",
    " ",
    " 09-NOV-2012:  This is the latest calibrated parameter file; calibrated to both the chemistry simulations",
    "               and the accumulation analysis.",
    " 02-FEB-2012:  Entered EPA-delivered WQC database calibrated critical value into parameter file.",
    "               This file was renamed 'BLM_4ligand_M-edulis_SMAV_2012-02-01.dat' on 02-FEB-2012; it",
    "               was derived from the file 'BLM_parameters_4ligand_09nov2011_12.dat'.",
    "               NOTE:  Critical value was then re-calibrated from '323.8' (median BL-accum for spec output with",
    "               this file) to '303.6464738'.  This number was the BL-accumulation of copper predicted by the 4-site",
    "               model in speciation mode given an input file with normalized water chemistry and the 3-site toxicity",
    "               prediction as the copper input (given the same, normalized chemistry.)",
    ""
  )
  Marine_file = withr::local_tempfile(fileext = ".dat")
  write(x = Marine_text, file = Marine_file)
  ex_Rfile = withr::local_tempfile(fileext = ".dat4")

  expect_no_error(ConvertWindowsParamFile(WindowsParamFile = Marine_file,
                                          RParamFile = ex_Rfile,
                                          MarineFile = TRUE))

  expect_no_error(DefineProblem(ParamFile = ex_Rfile))

})
test_that("ConvertWindowsParamFile unknown dbs file works", {

  Zn_file = "C:/Users/kellyc/Documents/BLM/parameter_files/fresh/Zn/Zn_algae_tox_chronic_09-12-24.dat"
  testthat::skip_if_not(file.exists(Zn_file))
  ex_Rfile = withr::local_tempfile(fileext = ".dat4")

  expect_no_error(ConvertWindowsParamFile(WindowsParamFile = Zn_file,
                                          RParamFile = ex_Rfile))

expect_no_error(DefineProblem(ParamFile = ex_Rfile))

})
