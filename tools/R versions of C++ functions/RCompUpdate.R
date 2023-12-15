#' @title Iterative step improvement in component concentrations
#'
#' @description RCompUpdate calculates an iterative improvement on the component
#' concentrations based on the Newton-Raphson solution from the current
#' iteration.
#'
#' @details If the iteration would cause the adjusted component concentrations
#' to be less than zero, then the component concentration is simply divided by
#' 10 for this iteration.
#'
#'
#' @param CompConcStep : Vector (NComp) of adjustments to the component
#'   concentrations
#' @param CompConc : Vector (NComp) of component concentrations, input values
#'   are from this iteration
#' @param CompName : Vector (NComp) with the names of the components
#'
#' @return  Vector CompConc (NComp) modified for the next iteration
RCompUpdate = function(CompConcStep, CompConc, CompName){

  (oldCompConc = CompConc)
  # CompConc[CompConcStep < oldCompConc] = (oldCompConc - CompConcStep)[CompConcStep < oldCompConc]
  # CompConc = (oldCompConc - CompConcStep / CompCtoM)
  CompConc = (oldCompConc - CompConcStep)
  # CompConc[CompConcStep >= oldCompConc] = oldCompConc[CompConcStep >= oldCompConc] / 10
  CompConc[CompConc <= 0] = oldCompConc[CompConc <= 0] / 10


  # for (iComp in 1:NComp){
  #   if (CompConcStep[iComp] >= oldCompConc[iComp]) {
  #     CompConc[iComp] = oldCompConc[iComp] / 10
  #   } else {
  #     CompConc[iComp] = CompConc[iComp] - CompConcStep[iComp]
  #   }
  #   if (CompConc[iComp] <= 0.0) {
  #     CompConc[iComp] = 1E-10
  #   }
  # }#NEXT cc

  names(CompConc) = CompName
  return(CompConc)

}
