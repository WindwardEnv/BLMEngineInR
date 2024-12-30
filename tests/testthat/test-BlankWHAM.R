test_that("BlankWHAM works", {
  expect_error(BlankWHAM(1), "unused argument")
  expect_no_error(BlankWHAM())
})
