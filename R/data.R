# This script contains data documentation.

#' @rdname TestData
#'
#' @title Simple Test Data with Two Components
#'
#' @docType data
#'
#' @description A collection of chemistry inputs and outputs from an existing
#'   BLM run. The only components are H and CO3. The only species are H(+),
#'   CO3(-2), HCO3(-), and H2CO3(0).
#'
#' @details Windward BLM (Ver 3.57.2.49, build 2022-03-30)\cr Parameter File:
#'   Cu_freshwater_acute_and_chronic_2017-01-17.dat\cr Switches: /W /Q /VER3.57
#'   /O3\cr
#' \itemize{
#'    \item Temp. = 25
#'    \item HA = 10
#'    \item pH = 7
#'    \item Cu = 1.57366E-028
#'    \item DOC = 1E-008
#'    \item Ca = 2.49513E-025
#'    \item Mg = 4.11438E-025
#'    \item Na = 4.34976E-025
#'    \item K = 2.55766E-025
#'    \item SO4 = 1.041E-025
#'    \item Cl = 2.82064E-025
#'    \item DIC = 0.0001
#'    \item S = 3.11876E-025
#'  }
#' @source Windward BLM (Ver 3.57.2.49, build 2022-03-30)
#'
#' @format \code{TestDataTotalConc} is an object of class \code{numeric} of
#'   length 2 containing total concentrations.
#'   \enumerate{
#'     \item{H}
#'     \item{CO3}
#'   }
#'
#' @keywords datasets
"TestDataTotalConc"


#' @rdname TestData
#' @format \code{TestDataFreeConc} is an object of class \code{numeric} of
#'   length 4 containing free concentrations.
#'   \enumerate{
#'     \item{H}
#'     \item{CO3}
#'     \item{HCO3}
#'     \item{H2CO3}
#'   }
"TestDataFreeConc"


#' @rdname TestData
#' @format \code{TestDataK} is an object of class \code{numeric} of
#'   length 4 containing equilibrium coefficients.
#'   \enumerate{
#'     \item{H = H}
#'     \item{CO3 = CO3}
#'     \item{HCO3 = H + CO3}
#'     \item{H2CO3 = 2*H + CO3}
#'   }
"TestDataK"


#' @rdname TestData
#' @format \code{TestDataStoich} is An object of class \code{matrix} (inherits
#'   from \code{array}) with 4 rows and 2 columns containing stoichiometry
#'   information for each reaction.
#'   \describe{
#'     \item{\code{rows}}{Species reactions:
#'       \enumerate{
#'         \item{H = H}
#'         \item{CO3 = CO3}
#'         \item{HCO3 = H + CO3}
#'         \item{H2CO2 = 2*H + CO3}
#'       }
#'     }
#'     \item{\code{cols}}{Components:
#'       \enumerate{
#'         \item{H}
#'         \item{CO3}
#'       }
#'     }
#'   }
"TestDataStoich"


#' @rdname Full_InorgData
#'
#' @docType data
#'
#' @title Simple Test DAta with Full Inorganic Components
#'
#' @description A collection of chemistry inputs and outputs from an existing
#'   BLM run with all inorganic components (H, Cu, DOC, Ca, Mg, Na, K, SO4, Cl,
#'   CO3, BL).
#'
#' @details Windward BLM (Ver 3.57.2.49, build 2022-03-30)\cr Parameter File:
#'   Cu_freshwater_acute_and_chronic_2017-01-17.dat\cr Switches: /W /Q /VER3.57
#'   /O3\cr
#' \itemize{
#'    \item Temp. = 25
#'    \item HA = 10
#'    \item pH = 7
#'    \item Cu = 1.57366E-028
#'    \item DOC = 1E-008
#'    \item Ca = 2.49513E-025
#'    \item Mg = 4.11438E-025
#'    \item Na = 4.34976E-025
#'    \item K = 2.55766E-025
#'    \item SO4 = 1.041E-025
#'    \item Cl = 2.82064E-025
#'    \item DIC = 0.0001
#'    \item S = 3.11876E-025
#'  }
#' @source Windward BLM (Ver 3.57.2.49, build 2022-03-30)
#'
#' @format \code{Full_InorgDataTotalConc} is an object of class \code{numeric} of
#'   length 11 containing total concentrations.
#'
#' @keywords datasets
"Full_InorgDataTotalConc"

#' @rdname Full_InorgData
#' @format \code{Full_InorgDataFreeConc} is an object of class \code{numeric} of
#'   length 32 containing free concentrations.
"Full_InorgDataFreeConc"

#' @rdname Full_InorgData
#' @format \code{Full_InorgDataK} is an object of class \code{numeric} of
#'   length 31 containing equilibrium coefficients.
"Full_InorgDataK"

#' @rdname Full_InorgData
#' @format \code{Full_InorgDataStoich} is An object of class \code{matrix} (inherits
#'   from \code{array}) with 32 rows and 11 columns containing stoichiometry
#'   information for each reaction.
#'   \describe{
#'     \item{\code{rows}}{Species reactions.}
#'     \item{\code{cols}}{Components.}
#'   }
"Full_InorgDataStoich"
