% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CHESSLog.R
\name{CHESSLog}
\alias{CHESSLog}
\title{Create a file that shows the problem in a different, human-friendly format}
\usage{
CHESSLog(
  ThisProblem,
  LogFilename = file.path(dirname(ThisProblem$ParamFile), "CHESSLOG.txt")
)
}
\arguments{
\item{ThisProblem}{A problem list object, such as returned by
`DefineProblem`.}

\item{LogFilename}{The path and file name of the created log file. By
default, it is "CHESSLOG.txt", placed in the same directory as ParamFile.}
}
\value{
invisibly returns TRUE
}
\description{
Create a file that shows the problem in a different, human-friendly format
}
\keyword{internal}
