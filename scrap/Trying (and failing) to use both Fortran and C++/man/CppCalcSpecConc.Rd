\name{CppCalcSpecConc}
\alias{CppCalcSpecConc}
\title{Calculate Species Concentrations with C++}
\description{
  Use this function to calculate species concentrations from a known set of free
  component ion concentrations. Version coded in C++.
}
\usage{
  CppCalcSpecConc(CConc, K, Stoich, NComp, NSpec)
}
\arguments{
  \item{CConc}{
    A vector of component concentrations for each of \code{NComp} components.
  }
  \item{K}{
    A vector of reaction equilibrium constants for each of \code{NSpec}
    reactions.
  }
  \item{Stoich}{
    A matrix of reaction stoichiometry, with \code{NSpec} rows and \code{NComp}
    columns.
  }
  \item{NComp}{
    The number of components in the equilibrium system.
  }
  \item{NSpec}{
    The number of species (reactions) in the equilibrium system.
  }
}
\details{
  SConc[i] = K[i] * prod(CConc^Stoich[i,])
}
\value{
  A vector of \code{NSpec} species concentrations.
}
\author{
  Kelly Croteau \email{kellyc@windwardenv.com}
}
\examples{
## Not Run
# data(TestDataFreeConc, TestDataK, TestDataStoich)
# CppCalcSpecConc(
#   CConc = TestDataFreeConc[1:2],
#   K = TestDataK[3:4],
#   Stoich = TestDataStoich[3:4,],
#   NComp = 2,
#   NSpec = 2
# )
# # [1] 8.388772e-05 1.922803e-05
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
