#' Calculate the Newton-Raphson step
#'
#' @param JacobianMatrix numeric matrix (NComp x NComp), the Jacobian matrix
#'   ("Z")
#' @param Resid numeric vector (NComp), the residuals = calculated totals -
#'   known totals
#' @param NComp integer, the number of components
#' @param CompType character vector (NComp), the type of component. It should be
#'   a fixed set of values (MassBal, FixedConc, Substituted, ChargeBal, SurfPot)
#' @param CompName character vector (NComp), the names of the components
#'
#' @return numeric vector (NComp), the N-R step to take for each component ("X"
#'   in C(i+1) = C(i) - X)
#' @export
RCalcStep = function(JacobianMatrix, Resid, NComp, CompType, CompName){

  i.solve = which(CompType != "FixedConc")
  n.solve = length(i.solve)

  (Resid.solve = matrix(Resid[i.solve], nrow = n.solve, ncol = 1,
                        dimnames = list(CompName[i.solve],NULL)))
  (JacobianMatrix.solve = JacobianMatrix[i.solve, i.solve, drop = F])

  # find the matrix inverse of JacobianMatrix by SVD
  JacobianMatrixsvd = svd(JacobianMatrix.solve)
  dinv = diag(1/JacobianMatrixsvd$d, nrow = n.solve, ncol = n.solve)
  U = JacobianMatrixsvd$u
  V = JacobianMatrixsvd$v
  UT = t(U)
  (JacobianMatrixinv = V %*% dinv %*% UT)

  # # If we wanted to not do SVD...
  # JacobianMatrixinv = solve(JacobianMatrix.solve)

  (CompConcStep.solve = as.numeric(JacobianMatrixinv %*% Resid.solve))

  CompConcStep = array(0, dim = NComp, dimnames = list(CompName))
  CompConcStep[i.solve] = CompConcStep.solve

  # CompConcStep[CompType == "FixedConc"] = 0

  return(CompConcStep)

}
