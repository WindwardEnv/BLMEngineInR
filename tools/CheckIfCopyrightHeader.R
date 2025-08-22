copyright.boilerplate.R = trimws(c(
  '# Copyright 2024 Windward Environmental LLC',
  '#',
  '# Licensed under the Apache License, Version 2.0 (the "License");',
  '# you may not use this file except in compliance with the License.',
  '# You may obtain a copy of the License at',
  '#',
  '# http://www.apache.org/licenses/LICENSE-2.0',
  '#',
  '# Unless required by applicable law or agreed to in writing, software',
  '# distributed under the License is distributed on an "AS IS" BASIS,',
  '# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.',
  '# See the License for the specific language governing permissions and',
  '# limitations under the License.'
))
generated.files = c("R/BLMEngineInR-package.R", "R/RcppExports.R", "src/RcppExports.cpp")
for (i in setdiff(list.files(path = "R", pattern = "[.]R", full.names = TRUE), generated.files)) {
  tmp = trimws(scan(file = i, what = character(), nlines = length(copyright.boilerplate.R), sep = "\n", quiet = TRUE))
  if (!all(tmp == copyright.boilerplate.R)) {
    print(i)
  }
}
for (i in setdiff(list.files(path = "data-raw", pattern = "[.]R", full.names = TRUE), generated.files)) {
  tmp = trimws(scan(file = i, what = character(), nlines = length(copyright.boilerplate.R), sep = "\n", quiet = TRUE))
  if (!all(tmp == copyright.boilerplate.R)) {
    print(i)
  }
}

copyright.boilerplate.cpp = gsub("^#", "//", copyright.boilerplate.R)
for (i in setdiff(list.files(path = "src", pattern = "[.]cpp", full.names = TRUE), generated.files)) {
  tmp = trimws(scan(file = i, what = character(), nlines = length(copyright.boilerplate.R), sep = "\n", quiet = TRUE))
  if (!all(tmp == copyright.boilerplate.cpp)) {
    print(i)
  }
}
for (i in setdiff(list.files(path = "src", pattern = "[.]h", full.names = TRUE), generated.files)) {
  tmp = trimws(scan(file = i, what = character(), nlines = length(copyright.boilerplate.R), sep = "\n", quiet = TRUE))
  if (!all(tmp == copyright.boilerplate.cpp)) {
    print(i)
  }
}
