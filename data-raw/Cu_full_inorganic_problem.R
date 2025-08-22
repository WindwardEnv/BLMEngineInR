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


mypfile = file.path("inst", "extdata", "ParameterFiles", "Cu_inorganic_only.dat4")

Cu_full_inorganic_problem = RemoveInVars(
  ThisProblem = Cu_full_organic_problem,
  InVarToRemove = c("DOC", "HA")
)
Cu_full_inorganic_problem$DoWHAM = FALSE
Cu_full_inorganic_problem$WHAM = BlankWHAM()
Cu_full_inorganic_problem$ParamFile = basename(mypfile)
WriteParamFile(ThisProblem = Cu_full_inorganic_problem, ParamFile = mypfile)

Cu_full_inorganic_noBL_problem = RemoveCriticalValues(
  RemoveMassCompartments(
    ThisProblem = Cu_full_inorganic_problem,
    MCToRemove = "Biotic"
  ),
  CAToRemove = 1:Cu_full_inorganic_problem$N["CAT"]
)
Cu_full_inorganic_noBL_problem$ParamFile = "Cu_inorganic_only_noBL.dat4"
WriteParamFile(ThisProblem = Cu_full_inorganic_noBL_problem,
               ParamFile = file.path("inst", "extdata", "ParameterFiles",
                                     "Cu_inorganic_only_noBL.dat4"))

usethis::use_data(Cu_full_inorganic_problem, overwrite = TRUE)
