test_that("WriteDetailedFile works", {
  myoutput = BLM(ThisProblem = carbonate_system_problem,
                 InputFile = system.file("extdata","InputFiles","carbonate_system_test.blm4",
                                         package = "BLMEngineInR",
                                         mustWork = TRUE), DoTox = FALSE)
  expect_no_error(WriteDetailedFile(OutList = myoutput,
                                    FileName = withr::local_tempfile(fileext = "xlsx")))
})
