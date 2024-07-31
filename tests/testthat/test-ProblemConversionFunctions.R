SomethingExtra = list(ExtraIntValue = 1L,
                      ExtraUnnamedVector = 1:3,
                      ExtraNamedVector = c(A = 1, B = 2, C = 3),
                      ExtraDFTable = data.frame(A = 1:5, B = 6:10),
                      ExtraMatrix = matrix(runif(12), 3, 4))

test_that("ConvertToDF works", {

  BlankDF = BlankProblem()
  ConvertedBlankDF = ConvertToDF(BlankProblemList())
  expect_equal(BlankDF, ConvertedBlankDF)

  BlankDFWithExtra = c(BlankDF, SomethingExtra)
  ConvertedWithExtraDF = ConvertToDF(c(BlankProblemList(), SomethingExtra))
  expect_equal(BlankDFWithExtra, ConvertedWithExtraDF)


})

test_that("ConvertToList works", {
  BlankList = BlankProblemList()
  ConvertedBlankList = ConvertToList(BlankProblem())
  expect_equal(BlankList, ConvertedBlankList)

  BlankListWithExtra = c(BlankList, SomethingExtra)
  ConvertedWithExtraList = ConvertToList(c(BlankProblem(), SomethingExtra))
  expect_equal(BlankListWithExtra, ConvertedWithExtraList)
})

test_that("Convert Back and Forth works", {
  expect_equal(carbonate_system_problem,
               ConvertToDF(ConvertToList(carbonate_system_problem)))
})
