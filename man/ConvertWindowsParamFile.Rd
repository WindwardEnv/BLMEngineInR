% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConvertWindowsParamFile.R
\name{ConvertWindowsParamFile}
\alias{ConvertWindowsParamFile}
\title{Convert From a Windows BLM Parameter File}
\usage{
ConvertWindowsParamFile(
  WindowsParamFile,
  RParamFile = NULL,
  RWHAMFile = NULL,
  MarineFile = FALSE
)
}
\arguments{
\item{WindowsParamFile}{Character string with the file path of the
Windows-format BLM parameter file. Typically will have the extension
".dat".}

\item{RParamFile}{(optional) Character string with the file path of the
R-format BLM parameter file to save.}

\item{RWHAMFile}{(optional) Character string with the file path of the
R-format WHAM parameter file to save.}

\item{MarineFile}{Boolean value - is this a marine file? If so, it uses a
lower mass value. In the Windows BLM, this is equivalent to using the "/M"
switch. Defaults to `FALSE`.}
}
\value{
The BLMEngineInR-compatible chemistry problem object. If RParamFile
  is provided, this will return invisibly.
}
\description{
Convert From a Windows BLM Parameter File
}
