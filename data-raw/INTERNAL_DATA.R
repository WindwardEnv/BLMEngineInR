# Copyright 2024 Windward Environmental LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

WHAM_V_LIST = DefineWHAM(
  WHAMFile = system.file("extdata", "WHAM", "WHAM_V.wdat",
                         package = "BLMEngineInR",
                         mustWork = TRUE)
)
WHAM_V_LIST$File = basename(WHAM_V_LIST$File)
WHAM_VI_LIST = DefineWHAM(
  WHAMFile = system.file("extdata", "WHAM", "WHAM_VI.wdat",
                         package = "BLMEngineInR",
                         mustWork = TRUE)
)
WHAM_VI_LIST$File = basename(WHAM_VI_LIST$File)
WHAM_VII_LIST = DefineWHAM(
  WHAMFile = system.file("extdata", "WHAM", "WHAM_VII.wdat",
                         package = "BLMEngineInR",
                         mustWork = TRUE)
)
WHAM_VII_LIST$File = basename(WHAM_VII_LIST$File)

usethis::use_data(WHAM_V_LIST, WHAM_VI_LIST, WHAM_VII_LIST,
                  internal = TRUE, overwrite = TRUE)
