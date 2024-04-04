test_that("BLM function works", {
  expect_equal(
    log10(
      unlist(BLM(
        ParamFile = system.file(
          file.path("extdata","ParameterFiles","full_inorg.dat4"),
          package = "BLMEngineInR",
          mustWork = TRUE),
        InputFile = system.file(
          file.path("extdata","InputFiles","full_inorg.blm4"),
          package = "BLMEngineInR",
          mustWork = TRUE),
        DoTox = FALSE
      )[1, c("IonicStrength", "Cu (mol/L)", "Ca (mol/L)", "Mg (mol/L)",
             "Na (mol/L)", "K (mol/L)", "SO4 (mol/L)",
             "Cl (mol/L)", "CO3 (mol/L)", "H (mol/L)",
             "BL1 (mol/g wt)", "BL1-Cu (mol/g wt)", "BL1-CuOH (mol/g wt)",
             "BL1-Ca (mol/g wt)", "BL1-Mg (mol/g wt)", "BL1-H (mol/g wt)",
             "BL1-Na (mol/g wt)", "OH (mol/L)", "HCO3 (mol/L)",
             "H2CO3 (mol/L)", "MgHCO3 (mol/L)", "MgCO3 (mol/L)",
             "MgSO4 (mol/L)", "CaHCO3 (mol/L)", "CaCO3 (mol/L)",
             "CaSO4 (mol/L)", "CuOH (mol/L)", "Cu(OH)2 (mol/L)",
             "CuSO4 (mol/L)", "CuCl (mol/L)", "CuCO3 (mol/L)",
             "Cu(CO3)2 (mol/L)", "CuHCO3 (mol/L)")])
    ),
    c(`IonicStrength` = -2.39401042386244,
      `Cu (mol/L)` = -8.85026986596705,
      `Ca (mol/L)` = -3.56493916654534,
      `Mg (mol/L)` = -3.34574119840722,
      `Na (mol/L)` = -2.9585825633873,
      `K (mol/L)` = -4.26993837894087,
      `SO4 (mol/L)` = -3.13634624376631,
      `Cl (mol/L)` = -4.27089922602816,
      `CO3 (mol/L)` = -5.66130216464019,
      `H (mol/L)` = -7.53960433551007,
      `BL1 (mol/g wt)` = -5.38028314032168,
      `BL1-Cu (mol/g wt)` = -6.94534296039693,
      `BL1-CuOH (mol/g wt)` = -8.07534194021004,
      `BL1-Ca (mol/g wt)` = -5.46001226097523,
      `BL1-Mg (mol/g wt)` = -5.24081429283711,
      `BL1-H (mol/g wt)` = -7.55028416050858,
      `BL1-Na (mol/g wt)` = -5.36926136819891,
      `OH (mol/L)` = -6.39660331532318,
      `HCO3 (mol/L)` = -2.98669747444535,
      `H2CO3 (mol/L)` = -4.23509415912218,
      `MgHCO3 (mol/L)` = -5.37622862696078,
      `MgCO3 (mol/L)` = -6.25662327126381,
      `MgSO4 (mol/L)` = -4.34166735038993,
      `CaHCO3 (mol/L)` = -5.5554265950989,
      `CaCO3 (mol/L)` = -6.23582123940193,
      `CaSO4 (mol/L)` = -4.63086531852805,
      `CuOH (mol/L)` = -8.88466313539843,
      `Cu(OH)2 (mol/L)` = -10.0450577797015,
      `CuSO4 (mol/L)` = -9.85619601794976,
      `CuCl (mol/L)` = -12.8359590461034,
      `CuCO3 (mol/L)` = -7.99115193882364,
      `Cu(CO3)2 (mol/L)` = -10.4824541034638,
      `CuHCO3 (mol/L)` = -7.66075729452061),
    tolerance = 0.00001
  )

  expect_error(BLM())
})
