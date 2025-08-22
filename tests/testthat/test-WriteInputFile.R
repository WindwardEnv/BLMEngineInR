test_that("WriteInputFile works", {

  myproblem = Cu_full_organic_problem
  myinputfile = withr::local_tempfile(fileext = ".blm4")
  write(c(
    "5",
    "ObsNum, ID,            ID2,                   Temp,  pH,   DOC,   HA,   Na,           Mg,           K,            Ca,           Cu,         SO4,          Cl,           CO3",
    "n/a,    n/a,           n/a,                   deg C, SU,   mgC/L, %,    mol/L,        mol/L,        mol/L,        mol/L,        mol/L,      mol/L,        mol/L,        mol/L",
    "1,      Full_Organic,  Hard ser 10,           15,    7.57, 0.01,  10,   1.375612e-04, 6.274425e-05, 6.713850e-06, 3.742700e-05, 4.2112E-07, 9.993587e-05, 6.699012e-06, 1.374837e-04",
    "2,      Full_Organic,  Hard ser 20,           15,    7.57, 0.01,  10,   2.751225e-04, 1.254885e-04, 1.342770e-05, 7.485400e-05, 4.2112E-07, 1.998717e-04, 1.339802e-05, 2.749675e-04",
    "3,      Full_Organic,  Hard ser 50,           15,    7.57, 0.01,  10,   6.878062e-04, 3.137212e-04, 3.356925e-05, 1.871350e-04, 4.2112E-07, 4.996794e-04, 3.349506e-05, 6.874188e-04",
    "4,      Full_Organic,  Hard ser 100,          15,    7.57, 0.01,  10,   1.375612e-03, 6.274425e-04, 6.713850e-05, 3.742700e-04, 4.2112E-07, 9.993587e-04, 6.699012e-05, 1.374838e-03",
    "5,      Full_Organic,  Hard ser 200,          15,    7.57, 0.01,  10,   0.0027512250, 0.0012548850, 0.0001342770, 0.0007485400, 4.2112E-07, 0.0019987175, 0.0001339802, 0.0027496750"
  ), file = myinputfile)


  myinputs = GetData(InputFile = myinputfile, ThisProblem = myproblem)

  mytestinputfile = withr::local_tempfile(fileext = ".blm4")

  expect_no_error(WriteInputFile(AllInput = myinputs, ThisProblem = myproblem, InputFile = mytestinputfile))

  expect_equal(scan(mytestinputfile, what = character(), quiet = TRUE, sep = ",", strip.white = TRUE)[1:13],
               scan(myinputfile, what = character(), quiet = TRUE, sep = ",", strip.white = TRUE)[1:13])

})
