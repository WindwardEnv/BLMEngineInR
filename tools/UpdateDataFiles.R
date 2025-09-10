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
    "Some data-creation files do not exist anymore:",
    paste(setdiff(source.order, list.files(path = "data-raw")),
          collapse = ", ")
  )
}


# This will stop data files from being re-sourced if nothing has changed
load(file.path("tools", "data-raw_mtime_previous.RData"))
data.raw.mtime = merge(
  data.frame(source.file = source.order,
             mtime.new = file.info(file.path("data-raw", source.order))$mtime),
  data.raw.mtime,
  all.x = TRUE,
  sort = FALSE
)
stopifnot(all(data.raw.mtime$source.file == source.order))
data.raw.mtime$update =
  as.logical(cumsum((data.raw.mtime$mtime.new > data.raw.mtime$mtime.prev)))

data.raw.mtime$update = data.raw.mtime$update |
  switch(menu(
    choices = c("Yes", "No"),
    title = paste0("These files have been updated since the most recent data ",
                   "file update. Do any of these warrant a data update?\n",
                   paste(list.files(path = "R")[
                     file.info(list.files(path = "R", full.names = TRUE))$mtime >
                       max(file.info(list.files(path = "data", full.names = TRUE))$mtime)
                   ], collapse = "\n"))
  ), Yes = TRUE, No = FALSE)

for (i in 1:length(source.order)) {
  if (data.raw.mtime$update[i]) {
    devtools::load_all()
    source(file.path("data-raw", source.order[i]))
  }
}

if (!all(list.files(path = "data-raw") %in% source.order)) {
  warning(
    "Some data-creation scripts were not updated:",
    paste(setdiff(list.files(path = "data-raw"), source.order),
          collapse = ", ")
  )
}

if (any(data.raw.mtime$update)) {
  data.raw.mtime$mtime.prev = data.raw.mtime$mtime.new
  data.raw.mtime$mtime.new = NULL
  data.raw.mtime$update = NULL
  save(data.raw.mtime, file = file.path("tools", "data-raw_mtime_previous.RData"))
}
