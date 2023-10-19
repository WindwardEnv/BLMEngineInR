test_that("defineProblem works", {
  expect_error(defineProblem())
  expect_type(defineProblem("Test"), "list")
  expect_equal(names(defineProblem("Test")), c("NComp","NSpec","K","logK","Stoich","CConc"))
  expect_equal(defineProblem("Test")$NComp, length(defineProblem("Test")$CConc))
})
