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
  wb = openxlsx::buildWorkbook(OutList, asTable = TRUE, colWidths = "auto")
  openxlsx::addWorksheet(wb, sheetName = "Additional Info")
  openxlsx::writeData(wb, sheet = "Additional Info", x = AdditionalInfo)
  openxlsx::saveWorkbook(wb, file = FileName, overwrite = TRUE)

  return(TRUE)
}
