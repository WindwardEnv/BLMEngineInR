#' @title Write a VERY Detailed Output File
#'
#' @description
#' This will write an output XLSX file with everything that is returned by the
#' `BLM` function. This includes inputs, concentrations, activities, etc.
#'
#' @param OutList The list object returned by the BLM function.
#' @param FileName The name of the file you'd like to write.
#' @param AdditionalInfo This vector will be included in the "Additional Info".
#'   By default, it will give the date/time the file was saved.
#'
#' @return Returns TRUE (invisibly) if successful.
#'
#' @export
WriteDetailedFile = function(
    OutList,
    FileName,
    AdditionalInfo = paste0("Saved on: ", Sys.time())
) {

  # wb = openxlsx::createWorkbook()
  # for (i in names(OutList)){
  #   openxlsx::addWorksheet(wb, sheetName = i)
  #   openxlsx::writeDataTable(wb, sheet = i, x = OutList[[i]])
  # }
  AdditionalInfo = c(OutList$TimeElapsed, AdditionalInfo)
  OutList$TimeElapsed = NULL
  wb = openxlsx::buildWorkbook(OutList, asTable = TRUE, colWidths = "auto")
  openxlsx::addWorksheet(wb, sheetName = "Additional Info")
  openxlsx::writeData(wb, sheet = "Additional Info", x = AdditionalInfo)
  openxlsx::saveWorkbook(wb, file = FileName, overwrite = TRUE)

  return(invisible(TRUE))
}
