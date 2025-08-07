test_that("DefineWHAM works", {

  mywfile = system.file("extdata","WHAM","WHAM_V.wdat",
                        package = "BLMEngineInR", mustWork = TRUE)
  expect_no_error(DefineWHAM(WHAMFile = mywfile))
  expect_no_error(DefineWHAM(WHAMVer = "VII"))

  compare_names = names(Cu_full_organic_problem$WHAM)[
    names(Cu_full_organic_problem$WHAM) %in% c("File","Ver") == FALSE]
  expect_equal(DefineWHAM(WHAMFile = mywfile)[compare_names],
               Cu_full_organic_problem$WHAM[compare_names])

})



