test_that("BLM function works", {
  # expect_equal(
  #   log10(
  #     unlist(BLM(
  #       ParamFile = system.file(
  #         file.path("extdata", "ParameterFiles", "full_inorg.dat4"),
  #         package = "BLMEngineInR",
  #         mustWork = TRUE),
  #       InputFile = system.file(
  #         file.path("extdata", "InputFiles", "full_inorg.blm4"),
  #         package = "BLMEngineInR",
  #         mustWork = TRUE),
  #       DoTox = FALSE
  #          )[1, c("IonicStrength", "Cu (mol/L)", "Ca (mol/L)", "Mg (mol/L)",
  #            "Na (mol/L)", "K (mol/L)", "SO4 (mol/L)",
  #            "Cl (mol/L)", "CO3 (mol/L)", "H (mol/L)",
  #            "BL1 (mol/kg wet)", "BL1-Cu (mol/kg wet)", "BL1-CuOH (mol/kg wet)",
  #            "BL1-Ca (mol/kg wet)", "BL1-Mg (mol/kg wet)", "BL1-H (mol/kg wet)",
  #            "BL1-Na (mol/kg wet)", "OH (mol/L)", "HCO3 (mol/L)",
  #            "H2CO3 (mol/L)", "MgHCO3 (mol/L)", "MgCO3 (mol/L)",
  #            "MgSO4 (mol/L)", "CaHCO3 (mol/L)", "CaCO3 (mol/L)",
  #            "CaSO4 (mol/L)", "CuOH (mol/L)", "Cu(OH)2 (mol/L)",
  #            "CuSO4 (mol/L)", "CuCl (mol/L)", "CuCO3 (mol/L)",
  #            "Cu(CO3)2 (mol/L)", "CuHCO3 (mol/L)")])
  #   ),
  #   c(
  #     `IonicStrength` = -2.39381947120872,
  #     `Cu (mol/L)` = -8.72065380062928,
  #     `Ca (mol/L)` = -3.56442861566559,
  #     `Mg (mol/L)` = -3.34523470072935,
  #     `Na (mol/L)` = -2.95841394992089,
  #     `K (mol/L)` = -4.26993837894087,
  #     `SO4 (mol/L)` = -3.13638512256423,
  #     `Cl (mol/L)` = -4.27089922643923,
  #     `CO3 (mol/L)` = -5.66130724198874,
  #     `H (mol/L)` = -7.53959804806886,
  #     `BL1 (mol/kg wet)` = -9.90448378839183,
  #     `BL1-Cu (mol/kg wet)` = -11.3399499609509,
  #     `BL1-CuOH (mol/kg wet)` = -12.4699489029081,
  #     `BL1-Ca (mol/kg wet)` = -9.98372477598718,
  #     `BL1-Mg (mol/kg wet)` = -9.76453086105093,
  #     `BL1-H (mol/kg wet)` = -12.0744848464346,
  #     `BL1-Na (mol/kg wet)` = -9.89329969024387,
  #     `OH (mol/L)` = -6.3996452244215,
  #     `HCO3 (mol/L)` = -2.98670584514424,
  #     `H2CO3 (mol/L)` = -4.2351009841194,
  #     `MgHCO3 (mol/L)` = -5.37575577768391,
  #     `MgCO3 (mol/L)` = -6.25617648465137,
  #     `MgSO4 (mol/L)` = -4.34126101779357,
  #     `CaHCO3 (mol/L)` = -5.55495656211838,
  #     `CaCO3 (mol/L)` = -6.23537343662893,
  #     `CaSO4 (mol/L)` = -4.63044444770621,
  #     `CuOH (mol/L)` = -8.75511139698053,
  #     `Cu(OH)2 (mol/L)` = -9.91556052526432,
  #     `CuSO4 (mol/L)` = -9.72667125965631,
  #     `CuCl (mol/L)` = -12.7063711838389,
  #     `CuCO3 (mol/L)` = -7.86158578647752,
  #     `Cu(CO3)2 (mol/L)` = -10.3528930284663,
  #     `CuHCO3 (mol/L)` = -7.53118489258915
  #   ),
  #
  #   tolerance = 0.0001
  # )
  expect_error(BLM())

})
