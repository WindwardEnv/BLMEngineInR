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
