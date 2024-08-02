test_that("StoichCompsToStoichMatrix works", {

  expect_equal(
    StoichCompsToStoichMatrix(
      SpecCompNames = TestSpecCompNames,
      SpecCompStoichs = TestSpecCompStoichs,
      CompName = TestCompName,
      SpecName = TestSpecName
    ),
    TestSpecStoichMatrix
  )
  expect_equal(
    StoichCompsToStoichMatrix(
      SpecCompNames = TestSpecCompNames,
      SpecCompStoichs = TestSpecCompStoichs
    ),
    TestSpecStoichMatrix
  )
  expect_equal(
    StoichCompsToStoichMatrix(
      SpecCompNames = TestSpecCompNames_diff_order,
      SpecCompStoichs = TestSpecCompStoichs_diff_order
    ),
    TestSpecStoichMatrix
  )
})
test_that("StoichCompsToEquation works", {

  expect_equal(
    StoichCompsToEquation(
      SpecCompNames = TestSpecCompNames,
      SpecCompStoichs = TestSpecCompStoichs,
      CompName = TestCompName,
      SpecName = TestSpecName
    ),
    TestSpecEquation
  )
  expect_equal(
    StoichCompsToEquation(
      SpecCompNames = TestSpecCompNames,
      SpecCompStoichs = TestSpecCompStoichs
    ),
    TestSpecEquation
  )
  expect_equal(
    StoichCompsToEquation(
      SpecCompNames = TestSpecCompNames_diff_order,
      SpecCompStoichs = TestSpecCompStoichs_diff_order
    ),
    TestSpecEquation
  )
})
test_that("StoichMatrixToEquation works", {

  expect_equal(
    StoichMatrixToEquation(
      SpecStoich = TestSpecStoichMatrix,
      SpecName = TestSpecName
    ),
    TestSpecEquation
  )
  expect_equal(
    StoichMatrixToEquation(
      SpecStoich = TestSpecStoichMatrix
    ),
    TestSpecEquation
  )
})
test_that("EquationToStoichMatrix works", {

  expect_equal(
    EquationToStoichMatrix(
      SpecEquation = TestSpecEquation,
      CompName = TestCompName
    ),
    TestSpecStoichMatrix
  )
  expect_equal(
    EquationToStoichMatrix(
      SpecEquation = TestSpecEquation_diff_order,
      CompName = TestCompName
    ),
    TestSpecStoichMatrix
  )
  expect_error(
    EquationToStoichMatrix(
      SpecEquation = TestSpecEquation_not_comps,
      CompName = TestCompName
    )
  )
  expect_equal(
    EquationToStoichMatrix(
      SpecEquation = TestSpecEquation_extra_terms,
      CompName = TestCompName
    ),
    TestSpecStoichMatrix
  )
})
test_that("StoichMatrixToStoichComps works", {

  expect_equal(
    StoichMatrixToStoichComps(
      SpecStoich = TestSpecStoichMatrix,
      CompName = TestCompName
    ),
    list(
      SpecCompNames = TestSpecCompNames,
      SpecCompStoichs = TestSpecCompStoichs
    )
  )
})
test_that("EquationToStoichComps works", {

  expect_equal(
    EquationToStoichComps(
      SpecEquation = TestSpecEquation,
      CompName = TestCompName
    ),
    list(
      SpecName = TestSpecName,
      SpecCompNames = TestSpecCompNames,
      SpecCompStoichs = TestSpecCompStoichs
    )
  )
  expect_equal(
    EquationToStoichComps(
      SpecEquation = TestSpecEquation_diff_order,
      CompName = TestCompName
    ),
    list(
      SpecName = TestSpecName,
      SpecCompNames = TestSpecCompNames,
      SpecCompStoichs = TestSpecCompStoichs
    )
  )
  expect_error(
    EquationToStoichComps(
      SpecEquation = TestSpecEquation_not_comps,
      CompName = TestCompName
    )
  )
  expect_equal(
    EquationToStoichComps(
      SpecEquation = TestSpecEquation_extra_terms,
      CompName = TestCompName
    ),
    list(
      SpecName = TestSpecName,
      SpecCompNames = TestSpecCompNames,
      SpecCompStoichs = TestSpecCompStoichs
    )
  )
})
test_that("EquationToStoich works", {

  expect_equal(
    EquationToStoich(
      SpecEquation = TestSpecEquation,
      CompName = TestCompName
    ),
    list(
      SpecStoich = TestSpecStoichMatrix,
      SpecName = TestSpecName,
      SpecCompNames = TestSpecCompNames,
      SpecCompStoichs = TestSpecCompStoichs
    )
  )
  expect_equal(
    EquationToStoich(
      SpecEquation = TestSpecEquation_diff_order,
      CompName = TestCompName
    ),
    list(
      SpecStoich = TestSpecStoichMatrix,
      SpecName = TestSpecName,
      SpecCompNames = TestSpecCompNames,
      SpecCompStoichs = TestSpecCompStoichs
    )
  )
  expect_error(
    EquationToStoich(
      SpecEquation = TestSpecEquation_not_comps,
      CompName = TestCompName
    )
  )
  expect_equal(
    EquationToStoich(
      SpecEquation = TestSpecEquation_extra_terms,
      CompName = TestCompName
    ),
    list(
      SpecStoich = TestSpecStoichMatrix,
      SpecName = TestSpecName,
      SpecCompNames = TestSpecCompNames,
      SpecCompStoichs = TestSpecCompStoichs
    )
  )
})
