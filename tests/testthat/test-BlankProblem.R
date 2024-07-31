test_that("BlankProblem works", {
  expect_equal(BlankProblemList(), ConvertToList(BlankProblem()))
  expect_equal(BlankProblem(), ConvertToDF(BlankProblemList()))
})
