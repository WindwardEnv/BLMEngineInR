test_that("WriteWHAM works", {
  mywfile = system.file(file.path("extdata", "WHAM", "WHAM_V.wdat"),
                        package = "BLMEngineInR", mustWork = TRUE)
  mywlist = DefineWHAM(WHAMFile = mywfile)

  mytmpwfile = withr::local_tempfile()
  expect_no_error(WriteWHAMFile(ThisWHAM = mywlist, WHAMFile = mytmpwfile))

  mytmpwlist = DefineWHAM(WHAMFile = mytmpwfile)
  compare.names = setdiff(names(mywlist), c("File"))
  expect_equal(mywlist[compare.names], mytmpwlist[compare.names])

})
