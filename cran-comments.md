## Submission

Changes made:

-   DARMA_USE_LEGACY was added to PKG_CXXFLAGS in Makevars and Makevars.win. I attempted to use DARMA_USE_CURRENT, but was getting the message "#pragma message: NOTE: option ARMA_CRIPPLED_LAPACK is not supported". I have not been able to diagnose this problem, and legacy arma should work for now, just to get this thing on CRAN.

## R CMD check environment

Windows x86_64-w64-mingw32 (64-bit), R version 4.5.1 (2025-06-13 ucrt)

## R CMD check results

0 errors ✔ \| 0 warnings ✔ \| 0 notes ✔

-   This is a new release.
