TestCompName = c("H", "CO3", "Cu")
TestSpecName = c("H", "OH", "HCO3", "H2CO3", "CuOH")
TestSpecCompNames = list(
  H = "H",
  OH = "H",
  HCO3 = c("H", "CO3"),
  H2CO3 = c("H", "CO3"),
  CuOH = c("H", "Cu")
)
TestSpecCompStoichs = list(
  H = 1L,
  OH = -1L,
  HCO3 = c(1L, 1L),
  H2CO3 = c(2L, 1L),
  CuOH = c(-1L, 1L)
)
TestSpecStoichMatrix = matrix(c(
#  H    CO3   Cu
  1L,   0L,   0L,
  -1L,   0L,  0L,
  1L,    1L,  0L,
  2L,    1L,  0L,
  -1L,   0L,  1L
), nrow = 5, ncol = 3, byrow = TRUE,
dimnames = list(TestSpecName, TestCompName))
TestSpecEquation = c(
  "H = 1 * H" ,
  "OH = -1 * H",
  "HCO3 = 1 * H + 1 * CO3",
  "H2CO3 = 2 * H + 1 * CO3",
  "CuOH = -1 * H + 1 * Cu"
)
TestSpecCompNames_diff_order = list(
  H = "H",
  OH = "H",
  HCO3 = c("CO3", "H"),#here
  H2CO3 = c("CO3", "H"),#here
  CuOH = c("Cu", "H")#here
)
TestSpecCompStoichs_diff_order = list(
  H = 1L,
  OH = -1L,
  HCO3 = c(1L, 1L),#here
  H2CO3 = c(1L, 2L),#here
  CuOH = c(1L, -1L)#here
)
TestSpecEquation_diff_order = c(
  "H = 1 * H" ,
  "OH = -1 * H",
  "HCO3 = 1 * CO3 + 1 * H",#here
  "H2CO3 = 1 * CO3 + 2 * H",#here
  "CuOH = 1 * Cu - 1 * H"#here
)
TestSpecEquation_not_comps = c(
  "H = 1 * H" ,
  "OH = 1 * OH",#here
  "HCO3 = 1 * H + 1 * CO3",
  "H2CO3 = 2 * H + 1 * CO3",
  "CuOH = 1 * Cu + 1 * OH"#here
)
TestSpecEquation_extra_terms = c(
  "H = 1 * H" ,
  "OH = -1 * H",
  "HCO3 = 1 * H + 1 * CO3 + 0 * Cu",#here
  "H2CO3 = 1 * H + 1 * H + 1 * CO3",#here
  "CuOH = -1 * H + 1 * Cu"
)


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
