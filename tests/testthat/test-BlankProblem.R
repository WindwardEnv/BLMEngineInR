test_that("BlankProblem works", {
  expect_error(BlankProblem(1), "unused argument")
  expect_equal(BlankProblemList(), ConvertToList(BlankProblem()))
  expect_equal(BlankProblem(), ConvertToDF(BlankProblemList()))
})
