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

mypfile = file.path("inst", "extdata", "ParameterFiles", "carbonate_system_only.dat4")

carbonate_system_problem = AddSpecies(
  ThisProblem = AddInComps(
    ThisProblem = water_problem,
    InCompName = "CO3",
    InCompCharge = -2,
    InCompMCName = "Water",
    InCompType = "MassBal",
    InCompActCorr = "Debye"
  ),
  SpecEquation = c(
    "HCO3 = 1 * H + 1 * CO3",
    "H2CO3 = 2 * H + 1 * CO3"
  ),
  SpecMCName = "Water",
  SpecActCorr = "Debye",
  SpecLogK = c(
    10.329, #HCO3
    16.681  #H2CO3
  ),
  SpecDeltaH = c(
    -14997.55155,  #HCO3
    -24166.23162  #H2CO3
  ),
  SpecTempKelvin = 298.1514609
)
carbonate_system_problem$ParamFile = basename(mypfile)
WriteParamFile(ThisProblem = carbonate_system_problem, ParamFile = mypfile)

usethis::use_data(carbonate_system_problem, overwrite = TRUE)
