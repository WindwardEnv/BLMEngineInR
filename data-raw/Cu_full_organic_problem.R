# Copyright 2024 Windward Environmental LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# mypfile = system.file("extdata", "ParameterFiles", "Cu_full_organic.dat4",
#                       package = "BLMEngineInR", mustWork = TRUE)

Cu.r.pfile = file.path("inst", "extdata", "ParameterFiles", "Cu_full_organic.dat4")
Cu.w.pfile = withr::local_tempfile(fileext = ".dat")

write(c(
  # Cu_freshwater_acute_and_chronic_2017-01-17.dat with
  # CuOH5%le_10-11-07.DAT critical value in CA Table:
  "Column model parameter file, Ver 3.00",
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
  # "[CRITICAL]: 0.03395",
  # "[ACUTE_DIV]: 2.0",
  # "[CHRONIC_DIV]: 3.22",
  # "[Z_EF]: 3.117",
  #
  "[CRITICAL START]:",
  "CA,Species,Test Type,Lifestage,Endpoint,Quantifier,References,Miscellaneous",
  "0.03395,FAV,Acute,ACR=3.22,FAV,NA,US EPA 2007,test duration: DIV=2.0",
  "0.05541,Ceriodaphnia dubia,Acute,Neonate (<24 hr),,EC50; LC50,Gensemer et al. 2002; Hyne et al. 2005; Naddy et al. 2002; Naddy et al. 2003; Van Genderen et al. 2007,SMEA calculated by median",
  "0.0099515,Ceriodaphnia dubia,Chronic,Neonate (<24 hr),reproduction - # of young,EC20,Wang et al. 2011a,SMEA calculated by geomean",
  "1.0223,Chironomus tentans,Acute,\"Larva, 1st instar\",death and immobility,EC50,Gauss et al. 1985,SMEA calculated by geomean",
  "0.035846,Daphnia magna,Acute,Neonate (<24 hr),death/immobility,LC50; EC50,Al-Reasi et al. 2012; Fulton and Meyer 2014; Ryan et al. 2009; Villavicencio et al. 2011,SMEA calculated by geomean",
  "0.076833,Daphnia magna,Chronic,Neonate (<24 hr),reproduction,MATC,De Schamphelaere and Janssen 2004,SMEA calculated by median",
  "0.062603,Daphnia pulex,Acute,Neonate (<24 hr),,LC50,Van Genderen et al. 2007,SMEA calculated by geomean",
  "0.010919,Daphnia pulex,Chronic,Neonate (<24 hr),survival,EC20,Winner 1985,SMEA calculated by geomean",
  "0.0097598,Daphnia pulicaria,Acute,,,LC50,Lind et al. Manuscript 1978,SMEA calculated by geomean",
  "2.2092,Lampsilis fasciola,Acute,Glochidia,,EC50,Gillis et al. 2008,SMEA calculated by median",
  "0.25703,Lampsilis siliquoidea,Acute,Glochidia,,LC50,Wang et al. 2007a; 2007c; Wang et al. 2007c,SMEA calculated by geomean",
  "0.037046,Lampsilis siliquoidea,Acute,Juvenile,,LC50,Wang et al. 2009,SMEA calculated by median",
  "15.656,Lepomis macrochirus,Acute,\"4.2 cm, 0.59g\",,LC50,Inglis and Davis 1972,SMEA calculated by geomean",
  "0.21719,Oncorhynchus mykiss,Acute,\"swim-up, 0.25 g\",,LC50,Cacela et al. 1996,SMEA calculated by geomean",
  "0.99652,Oncorhynchus tshawytscha,Acute,\"Swim-up, 0.36-0.45 g\",,LC50,Welsh et al. 2000,SMEA calculated by geomean",
  "0.40192,Pimephales promelas,Acute,\"Larva, 1-7 d; <24 h, 0.68 mg\",,LC50,Lind et al. Manuscript 1978; Van Genderen et al. 2007; Welsh et al. 1993,SMEA calculated by median",
  "0.074349,Pimephales promelas,Chronic,\"Larva, <24 hr\",biomass,EC20,Besser et al. 2001; 2005,SMEA calculated by geomean",
  "1.5533,Utterbackia imbecillis,Acute,Juvenile,,LC50,Keller unpublished 2000 memo,SMEA calculated by geomean",
  "0.0054227,Villosa iris,Chronic,\"3 mo. old, 2.0 mm\",growth-weight,EC20,Wang et al. 2011a,SMEA calculated by geomean",
  "0.013496,Acipenser transmontanus,Acute,\"Columbia River, 26 dph, 0.08 g, 2.5 cm; Kootenai River, 38 dph, 0.07 g, 2.4 cm; Juvenile, 40 dph\",,LC50,Little et al. 2012; Vardy et al. 2013,SMEA calculated by geomean",
  "0.005297,Acipenser transmontanus,Chronic,\"Larva, 1 dph, 12.7 mm, 8.57 mg dry weight\",growth-dry weight,EC20,Wang et al. 2014,SMEA calculated by geomean",
  "3.70,Oncorhynchus mykiss,Acute,,,,From the \"Cu_Rainbow_Trout_06-10-07.DAT\" parameter file,",
  "0.119,Daphnia magna,Acute,,,,From the \"Cu_Daphnia_Magna_06-10-07.DAT\" parameter file,",
  "0.0447,Daphnia pulex,Acute,,,,From the \"Cu_Daphnia_Pulex_06-10-07.DAT\" parameter file,",
  "5.48,Pimephales promelas,Acute,,,,From the \"Cu_Fathead_Minnow_06-10-07.DAT\" parameter file,",
  "0.021704,Ceriodaphnia dubia,Acute,Neonate (<24h),survival,EC50,\"Larry Walker Associates, 2013\",",
  "[CRITICAL END]",
  "",
  "",
  "[END]:",
  "",
  "",
  "",
  ""
), file = Cu.w.pfile)

Cu_full_organic_problem = RemoveInComps(
  ThisProblem = AddInLabs(
    ThisProblem = RemoveInLabs(
      ThisProblem = ConvertWindowsParamFile(WindowsParamFile = Cu.w.pfile),
      InLabToRemove = c("Site Label", "Sample Label")
    ),
    InLabName = c("ObsNum", "ID", "ID2")
  ),
  InCompToRemove = "S"
)
WriteParamFile(ThisProblem = Cu_full_organic_problem, ParamFile = Cu.r.pfile)

# Cu_full_organic_problem = DefineProblem(ParamFile = Cu.r.pfile, WriteLog = FALSE)
Cu_full_organic_problem$ParamFile = basename(Cu.r.pfile)
Cu_full_organic_problem$WHAM$File = basename(Cu_full_organic_problem$WHAM$File)

usethis::use_data(Cu_full_organic_problem, overwrite = TRUE)
