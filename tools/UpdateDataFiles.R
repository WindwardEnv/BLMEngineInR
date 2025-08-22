# Data-creation scripts in /data-raw should be sourced in a particular sequence
# to ensure dependencies are all updated.

source.order = c(
  "MW.R",
  "INTERNAL_DATA.R",
  "water_problem.R",
  "All_WATER23_reactions.R",
  "All_NIST20170203_reactions.R",
  "carbonate_system_problem.R",
  "Cu_full_organic_problem.R",
  "Cu_full_inorganic_problem.R",
  "Ni_full_organic_problem.R",
  "Ni_HCO3_full_organic_problem.R"
)

if (!all(source.order %in% list.files(path = "data-raw"))) {
  stop(
    "Some data-creation do not exist anymore:",
    paste(setdiff(source.order, list.files(path = "data-raw")),
          collapse = ", ")
  )
}

for (i in source.order) {
  devtools::load_all()
  source(file.path("data-raw", i))
}

if (!all(list.files(path = "data-raw") %in% source.order)) {
  warning(
    "Some data-creation scripts were not updated:",
    paste(setdiff(list.files(path = "data-raw"), source.order),
          collapse = ", ")
  )
}



