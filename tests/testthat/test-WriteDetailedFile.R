test_that("WriteDetailedFile works", {

  myinputsdf = data.frame(
    ID = "Test",
    Temp = 25,
    pH = 7,
    CO3 = 1E-4
  )
  myinputs = MatchInputsToProblem(DFInputs = myinputsdf, ThisProblem = carbonate_system_problem)

  myoutput = BLM(ThisProblem = carbonate_system_problem,
                 AllInput = myinputs, DoTox = FALSE)
  expect_no_error(WriteDetailedFile(OutList = myoutput,
                                    FileName = withr::local_tempfile(fileext = "xlsx")))
})
