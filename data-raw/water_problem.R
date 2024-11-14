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

water_MC_problem = AddInLabs(
  ThisProblem = AddMassCompartments(ThisProblem = BlankProblem(),
                                     MassName = "Water",
                                     MassAmt = 1.0,
                                     MassUnit = "L"),
  InLabName = "ID")

water_problem = AddSpecies(
  ThisProblem = AddInVars(ThisProblem = water_MC_problem,
                          InVarName = c("Temperature", "pH"),
                          InVarMCName = "Water",
                          InVarType = c("Temperature", "pH")),
  SpecEquation = "OH = -1 * H",
  SpecMCName = "Water",
  SpecActCorr = "Debye",
  SpecLogK = -13.997,
  SpecDeltaH = -55810,
  SpecTempKelvin = 298.15)

usethis::use_data(water_MC_problem, overwrite = TRUE)
usethis::use_data(water_problem, overwrite = TRUE)
