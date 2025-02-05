% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WriteWHAMFile.R
\name{WriteWHAMFile}
\alias{WriteWHAMFile}
\title{Write a WHAM Parameter File}
\usage{
WriteWHAMFile(ThisWHAM, WHAMFile, Notes = ThisWHAM$Notes)
}
\arguments{
\item{ThisWHAM}{A list object with a structure like that returned by
`BlankWHAM()`.}

\item{WHAMFile}{a character value, indicating the file path and name of the
WHAM parameter file to write.}

\item{Notes}{A character vector of additional notes to include at the bottom
of the WHAM parameter file. The text "written by USERNAME from R:
YYYY-MM-DD HH:MM:SS" will always be written, regardless of the value of
this argument. Be default, this will be filled in with a "Notes" item in
`ThisWHAM`, if available.}
}
\value{
ThisProblem, with the ParamFile element  changed to the ParamFile
  argument.
}
\description{
This function will take a WHAM parameter list object and turn it
  into a WHAM parameter file, effectively doing the opposite of `DefineWHAM`.
}
\examples{
tf = tempfile()
WriteWHAMFile(ThisWHAM = Cu_full_organic_problem$WHAM, WHAMFile = tf)
DefineWHAM(WHAMFile = tf)

}
