\name{TestDataStoich}
\alias{TestDataStoich}
\docType{data}
\title{Reaction Stoiciometry for Test Data}
\description{
  Reaction stoichiometry for species reactions in a simple carbonate equilibrium
  system. The only components are hydrogen (H) and carbonate (CO3). The only
  species are H(+), CO3(-2), HCO3(-), and H2CO3(0).
}
\usage{TestDataStoich}
\format{
  A matrix with 4 rows and 2 columns.
  \describe{
    \item{\code{rows}}{Species reactions:\cr
        \enumerate{
          \item{H = H}
          \item{CO3 = CO3}
          \item{HCO3 = H + CO3}
          \item{H2CO2 = 2*H + CO3}
        }
      }
    \item{\code{cols}}{Components: \cr
        \enumerate{
            \item{H}
            \item{CO3}
        }
      }
  }
}
\details{
  Windward BLM (Ver 3.57.2.49, build 2022-03-30)\cr
  PFile: Cu_freshwater_acute_and_chronic_2017-01-17.dat\cr
  Switches: /W /Q /VER3.57 /O3\cr
  Input chemistry:
  \itemize{
    \item Temp. = 25
    \item HA = 10
    \item pH = 7
    \item Cu = 1.57366E-028
    \item DOC = 1E-008
    \item Ca = 2.49513E-025
    \item Mg = 4.11438E-025
    \item Na = 4.34976E-025
    \item K = 2.55766E-025
    \item SO4 = 1.041E-025
    \item Cl = 2.82064E-025
    \item DIC = 0.0001
    \item S = 3.11876E-025
  }
}
\source{Windward BLM (Ver 3.57.2.49, build 2022-03-30)}
\examples{data(TestDataStoich)}
\keyword{datasets}
